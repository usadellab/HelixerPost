use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use threadpool::ThreadPool;
use std::path::Path;
use std::sync::mpsc::channel;

use crate::analysis::extractor::{BasePredictionExtractor, ComparisonExtractor};
use crate::analysis::gff_conv::hmm_solution_to_gff;
use crate::analysis::hmm::{HmmConfig, HmmStateRegion, PredictionHmm};
use crate::analysis::rater::{RatingWriter, SequenceRater, SequenceRating};
use crate::analysis::window::BasePredictionWindowThresholdIterator;
use crate::gff::GffWriter;
use crate::results::conv::{ArrayConvInto, Bases, ClassPrediction, PhasePrediction};
use crate::results::{Sequence, Species};
use std::io::Write;

pub mod extractor;
pub mod gff_conv;
pub mod hmm;
pub mod rater;
pub mod window;

/// Top-level config file shape: one section per pipeline pass. Each section's
/// missing fields fall back to the section's `Default` (via `#[serde(default)]`),
/// so a config file can override only the knobs it cares about.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FileConfig {
    pub window: WindowConfig,
    pub hmm: HmmConfig,
    pub filter: FilterConfig,
}

impl FileConfig {
    pub fn load_yaml(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading config file {}", path.display()))?;
        serde_yaml::from_str(&text)
            .with_context(|| format!("parsing config file {}", path.display()))
    }

    pub fn to_yaml(&self) -> Result<String> {
        serde_yaml::to_string(self).context("serialising config to YAML")
    }
}

/// Pass 1: sliding-window candidate gene-region detection.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct WindowConfig {
    pub window_size: usize,
    pub edge_threshold: f32,
    pub peak_threshold: f32,
}

impl Default for WindowConfig {
    fn default() -> Self {
        Self {
            window_size: 100,
            edge_threshold: 0.1,
            peak_threshold: 0.8,
        }
    }
}

impl WindowConfig {
    pub fn with_overrides(
        self,
        window_size: Option<usize>,
        edge_threshold: Option<f32>,
        peak_threshold: Option<f32>,
    ) -> Self {
        Self {
            window_size: window_size.unwrap_or(self.window_size),
            edge_threshold: edge_threshold.unwrap_or(self.edge_threshold),
            peak_threshold: peak_threshold.unwrap_or(self.peak_threshold),
        }
    }
}

/// Pass 3: post-HMM filtering of decoded gene models. More knobs to come.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FilterConfig {
    pub min_coding_length: usize,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self { min_coding_length: 60 }
    }
}

impl FilterConfig {
    pub fn with_overrides(self, min_coding_length: Option<usize>) -> Self {
        Self {
            min_coding_length: min_coding_length.unwrap_or(self.min_coding_length),
        }
    }
}

pub struct Analyzer<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>> {
    bp_extractor: BasePredictionExtractor<'a, TC, TP>,
    comp_extractor: ComparisonExtractor<'a>,
    window_cfg: WindowConfig,
    filter_cfg: FilterConfig,
    hmm_cfg: HmmConfig,
    thread_pool: ThreadPool,
}

impl<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>>
    Analyzer<'a, TC, TP>
{
    pub fn new(
        bp_extractor: BasePredictionExtractor<'a, TC, TP>,
        comp_extractor: ComparisonExtractor<'a>,
        window_cfg: WindowConfig,
        filter_cfg: FilterConfig,
        hmm_cfg: HmmConfig,
        thread_count: usize,
    ) -> Analyzer<'a, TC, TP> {
        let thread_pool = ThreadPool::new(thread_count);

        Analyzer {
            bp_extractor,
            comp_extractor,
            window_cfg,
            filter_cfg,
            hmm_cfg,
            thread_pool,
        }
    }

    pub fn has_ref(&self) -> bool {
        self.comp_extractor.has_ref()
    }

    fn process_sequence_1d<W, BPI>(
        &self,
        species: &Species,
        seq: &Sequence,
        rev: bool,
        bp_iter: BasePredictionWindowThresholdIterator<BPI>,
        gene_idx: &mut usize,
        rater: &mut SequenceRater,
        gff_writer: &mut GffWriter<W>,
    ) -> Result<(usize, usize)>
    where
        W: Write,
        BPI: Iterator<Item = (Bases, ClassPrediction, PhasePrediction)>,
    {
        let min_coding_length = self.filter_cfg.min_coding_length;
        let mut window_count = 0;
        let mut window_length_total = 0;

        let (tx, rx) = channel();

        for (bp_vec, _total_vec, start_pos, _peak) in bp_iter {
            window_length_total += bp_vec.len();

            let end_pos = start_pos + bp_vec.len();
            println!("Queuing a window from {} to {} (length: {})", start_pos, end_pos, bp_vec.len());

            let tx = tx.clone();
            let hmm_cfg = self.hmm_cfg.clone();
            self.thread_pool.execute(move || {
                let len = bp_vec.len();
                let hmm = PredictionHmm::new(bp_vec, hmm_cfg);
                let maybe_solution = hmm.solve();
                // Send can only fail if the receiver was dropped; rx lives on the
                // current thread until the for-loop below drains the channel, so
                // a failure here is a genuine invariant break.
                tx.send((window_count, start_pos, end_pos, maybe_solution))
                    .expect("HMM worker lost its result channel");

                println!("Solved a window from {} to {} (length: {})", start_pos, end_pos, len);
            });

            window_count += 1;
        }

        let mut results = Vec::with_capacity(window_count);
        for _ in 0..window_count {
            results.push((0, 0, None));
        }

        for _ in 0..window_count {
            let (index, start_pos, end_pos, maybe_solution) =
                rx.recv().expect("HMM worker dropped a window result");
            results[index] = (start_pos, end_pos, maybe_solution);
        }

        for (start_pos, end_pos, maybe_solution) in results.into_iter() {
            let Some(solution) = maybe_solution else {
                bail!(
                    "HMM produced no solution for {} window {}-{}",
                    seq.get_name(),
                    start_pos,
                    end_pos
                );
            };

            let solution_regions = solution.trace_regions();
            let genes = HmmStateRegion::split_genes(solution_regions);

            for (gene_regions, coding_length) in genes.iter() {
                rater.rate_window_regions(
                    start_pos,
                    gene_regions,
                    (*coding_length < min_coding_length) && (*coding_length > 0),
                );
            }

            let gff_records = hmm_solution_to_gff(
                genes,
                species.get_name(),
                seq.get_name(),
                "Helixer",
                rev,
                start_pos,
                seq.get_length(),
                min_coding_length,
                gene_idx,
            );
            gff_writer
                .write_records(&gff_records)
                .with_context(|| format!("writing GFF records for {}", seq.get_name()))?;
        }

        Ok((window_count, window_length_total))
    }

    pub fn process_sequence<GW: Write, RW: Write>(
        &self,
        species: &Species,
        seq: &Sequence,
        fwd_rating: &mut SequenceRating,
        rev_rating: &mut SequenceRating,
        gff_writer: &mut GffWriter<GW>,
        rating_writer: &mut RatingWriter<RW>,
    ) -> Result<(usize, usize)> {
        let id = seq.get_id();
        println!(
            "  BP_Extractor for Sequence {} - ID {}",
            seq.get_name(),
            id.inner()
        );

        let mut gene_idx = 1;

        let fwd_bp_iter = BasePredictionWindowThresholdIterator::new(
            self.bp_extractor.fwd_iterator(id),
            self.window_cfg.window_size,
            self.window_cfg.edge_threshold,
            self.window_cfg.peak_threshold,
        )
        .context("forward base-prediction iterator failed to initialise")?;

        let mut fwd_comp_rater = SequenceRater::new(
            self.comp_extractor.fwd_iterator(id),
            seq.get_length() as usize,
        );

        let (fwd_windows, fwd_window_len) = self.process_sequence_1d(
            species,
            seq,
            false,
            fwd_bp_iter,
            &mut gene_idx,
            &mut fwd_comp_rater,
            gff_writer,
        )?;
        let fwd_seq_rating = fwd_comp_rater
            .calculate_stats(species, seq, false, rating_writer)
            .with_context(|| format!("writing forward rating for {}", seq.get_name()))?;
        println!(
            "Rating for Sequence {} - ID {} - Forward",
            seq.get_name(),
            id.inner()
        );
        fwd_seq_rating.dump(self.comp_extractor.has_ref());

        fwd_rating.accumulate(&fwd_seq_rating);

        let rev_bp_iter = BasePredictionWindowThresholdIterator::new(
            self.bp_extractor.rev_iterator(id),
            self.window_cfg.window_size,
            self.window_cfg.edge_threshold,
            self.window_cfg.peak_threshold,
        )
        .context("reverse base-prediction iterator failed to initialise")?;
        let mut rev_comp_rater = SequenceRater::new(
            self.comp_extractor.rev_iterator(id),
            seq.get_length() as usize,
        );

        let (rev_windows, rev_window_len) = self.process_sequence_1d(
            species,
            seq,
            true,
            rev_bp_iter,
            &mut gene_idx,
            &mut rev_comp_rater,
            gff_writer,
        )?;
        let rev_seq_rating = rev_comp_rater
            .calculate_stats(species, seq, false, rating_writer)
            .with_context(|| format!("writing reverse rating for {}", seq.get_name()))?;
        println!(
            "Rating for Sequence {} - ID {} - Reverse",
            seq.get_name(),
            id.inner()
        );
        rev_seq_rating.dump(self.comp_extractor.has_ref());

        rev_rating.accumulate(&rev_seq_rating);

        Ok((fwd_windows + rev_windows, fwd_window_len + rev_window_len))
    }
}
