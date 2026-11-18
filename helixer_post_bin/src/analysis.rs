use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use threadpool::ThreadPool;
use std::path::Path;
use std::collections::VecDeque;
use std::sync::mpsc::{channel, Receiver};

use crate::analysis::extractor::{BasePredictionExtractor, ComparisonExtractor};
use crate::analysis::gff_conv::hmm_solution_to_gff;
use crate::analysis::hmm::{HmmConfig, HmmStateRegion, PredictionHmm};
use crate::analysis::rater::{RatingWriter, SequenceRater, SequenceRating};
use crate::analysis::window::BasePredictionWindowThresholdIterator;
use crate::gff::GffWriter;
use crate::progress::Progress;
use crate::results::conv::{ArrayConvInto, ClassPrediction, PhasePrediction};
use crate::results::{Sequence, Species};
use crate::verbose::is_verbose;
use crate::{log_info, vprintln};
use std::io::Write;

pub mod extractor;
pub mod gff_conv;
pub mod hmm;
pub mod rater;
pub mod window;

/// Cap on the windows queued or being decoded at a time, per worker thread. Each one holds the base-wise predictions
/// of its window in memory; a few per thread keep the workers busy while the main thread finalises strands.
const WINDOWS_IN_FLIGHT_PER_THREAD: usize = 4;

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

    /// Decodes the windows of all `seqs` on the worker pool and finalises the strands into `sink` in input order,
    /// forward before reverse strand. Strands are not decoded one after the other: while the windows of one strand
    /// are being decoded, the windows of the following strands are queued as well, at most
    /// `WINDOWS_IN_FLIGHT_PER_THREAD` per thread. A strand is finalised once all of its windows and those of all
    /// earlier strands are decoded, so the output does not depend on the order in which the workers finish.
    pub fn process_sequences<GW: Write, RW: Write>(
        &self,
        species: &Species,
        seqs: &[&Sequence],
        sink: &mut Sink<'_, GW, RW>,
    ) -> Result<()> {
        let max_in_flight = WINDOWS_IN_FLIGHT_PER_THREAD * self.thread_pool.max_count();
        let (tx, rx) = channel();
        let mut pending = PendingStrands::default();
        let mut in_flight = 0;

        for seq in seqs {
            vprintln!(
                "Processing Sequence {} - ID {} ({}bp)",
                seq.get_name(),
                seq.get_id().inner(),
                seq.get_length()
            );
            for rev in [false, true] {
                let strand_idx = pending.push(seq, rev);
                let id = seq.get_id();
                let bp_iter = BasePredictionWindowThresholdIterator::new(
                    if rev { self.bp_extractor.rev_iterator(id) } else { self.bp_extractor.fwd_iterator(id) },
                    self.window_cfg.window_size,
                    self.window_cfg.edge_threshold,
                    self.window_cfg.peak_threshold,
                )
                .with_context(|| format!("{} base-prediction iterator failed to initialise", strand_name(rev)))?;

                for (bp_vec, _total_vec, start_pos, _peak) in bp_iter {
                    while in_flight >= max_in_flight {
                        self.receive_window(&rx, species, &mut pending, sink)?;
                        in_flight -= 1;
                    }

                    let end_pos = start_pos + bp_vec.len();
                    vprintln!("Queuing a window from {} to {} (length: {})", start_pos, end_pos, bp_vec.len());
                    let window_idx = pending.add_window(strand_idx);
                    in_flight += 1;

                    let tx = tx.clone();
                    let hmm_cfg = self.hmm_cfg.clone();
                    self.thread_pool.execute(move || {
                        let genes = PredictionHmm::new(bp_vec, hmm_cfg)
                            .solve()
                            .map(|solution| HmmStateRegion::split_genes(solution.trace_regions()));
                        vprintln!("Solved a window from {} to {} (length: {})", start_pos, end_pos, end_pos - start_pos);
                        // fails only after the main thread stopped on an error and dropped the receiver;
                        // the result is not needed then
                        let _ = tx.send((strand_idx, window_idx, DecodedWindow { start_pos, end_pos, genes }));
                    });
                }
                pending.close(strand_idx);
                self.finalize_ready(species, &mut pending, sink)?;
            }
        }

        while in_flight > 0 {
            self.receive_window(&rx, species, &mut pending, sink)?;
            in_flight -= 1;
        }
        Ok(())
    }

    /// Waits for the next decoded window, stores it and finalises the strands it completes.
    fn receive_window<GW: Write, RW: Write>(
        &self,
        rx: &Receiver<WindowMessage>,
        species: &Species,
        pending: &mut PendingStrands<'_>,
        sink: &mut Sink<'_, GW, RW>,
    ) -> Result<()> {
        let (strand_idx, window_idx, window) = rx.recv().expect("HMM worker dropped a window result");
        pending.store(strand_idx, window_idx, window);
        self.finalize_ready(species, pending, sink)
    }

    /// Finalises the leading strands whose windows are all queued and decoded.
    fn finalize_ready<GW: Write, RW: Write>(
        &self,
        species: &Species,
        pending: &mut PendingStrands<'_>,
        sink: &mut Sink<'_, GW, RW>,
    ) -> Result<()> {
        while let Some(strand) = pending.pop_ready() {
            let windows = strand.windows.into_iter().map(|window| window.expect("window result missing")).collect();
            self.finalize_strand(species, strand.seq, strand.rev, windows, sink)?;
        }
        Ok(())
    }

    /// Writes the GFF records and the rating of one decoded strand and adds it to the totals of
    /// `sink`. The forward strand opens the sequence (region header, gene numbering from 1), the
    /// reverse strand closes it (gene count, progress).
    fn finalize_strand<GW: Write, RW: Write>(
        &self,
        species: &Species,
        seq: &Sequence,
        rev: bool,
        windows: Vec<DecodedWindow>,
        sink: &mut Sink<'_, GW, RW>,
    ) -> Result<()> {
        let min_coding_length = self.filter_cfg.min_coding_length;
        let id = seq.get_id();
        if !rev {
            sink.gff_writer
                .write_region_header(seq.get_name(), seq.get_length())
                .with_context(|| format!("writing sequence-region header for {}", seq.get_name()))?;
            sink.gene_idx = 1;
        }

        let comp_iter = if rev { self.comp_extractor.rev_iterator(id) } else { self.comp_extractor.fwd_iterator(id) };
        let mut rater = SequenceRater::new(comp_iter, seq.get_length() as usize);

        for DecodedWindow { start_pos, end_pos, genes } in windows {
            let Some(genes) = genes else {
                bail!(
                    "HMM produced no solution for {} window {}-{}",
                    seq.get_name(),
                    start_pos,
                    end_pos
                );
            };
            sink.windows += 1;
            sink.bp += end_pos - start_pos;

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
                &mut sink.gene_idx,
            );
            sink.gff_writer
                .write_records(&gff_records)
                .with_context(|| format!("writing GFF records for {}", seq.get_name()))?;
        }

        let rating = rater
            .calculate_stats(species, seq, rev, sink.rating_writer)
            .with_context(|| format!("writing {} rating for {}", strand_name(rev), seq.get_name()))?;
        if is_verbose() {
            log_info!(
                "Rating for Sequence {} - ID {} - {}",
                seq.get_name(),
                id.inner(),
                if rev { "Reverse" } else { "Forward" }
            );
            rating.dump(self.comp_extractor.has_ref());
        }

        if rev {
            sink.rev_rating.accumulate(&rating);
            // gene_idx starts at 1 and is incremented once per written gene
            let genes = sink.gene_idx - 1;
            sink.genes += genes;
            sink.progress.sequence_done(seq.get_length(), sink.block_counts[id.inner()], genes);
        } else {
            sink.fwd_rating.accumulate(&rating);
        }
        Ok(())
    }
}

fn strand_name(rev: bool) -> &'static str {
    if rev { "reverse" } else { "forward" }
}

/// Genes of one window, each as its HMM state regions plus its coding length.
type Genes = Vec<(Vec<HmmStateRegion>, usize)>;

/// Strand index, window index within the strand and HMM result of one window, sent from a worker.
type WindowMessage = (usize, usize, DecodedWindow);

/// HMM result of one candidate window; `genes` is `None` when the HMM found no solution.
struct DecodedWindow {
    start_pos: usize,
    end_pos: usize,
    genes: Option<Genes>,
}

/// Strands whose windows are queued or decoded but not yet finalised, in input order. A strand is addressed by its
/// index among all strands of the run, `first_idx` being the index of the front strand.
#[derive(Default)]
struct PendingStrands<'s> {
    strands: VecDeque<PendingStrand<'s>>,
    first_idx: usize,
}

struct PendingStrand<'s> {
    seq: &'s Sequence,
    rev: bool,
    /// one entry per queued window, filled in when its result arrives
    windows: Vec<Option<DecodedWindow>>,
    received: usize,
    /// set once all windows of the strand are queued
    closed: bool,
}

impl<'s> PendingStrands<'s> {
    /// Appends a strand and returns its index.
    fn push(&mut self, seq: &'s Sequence, rev: bool) -> usize {
        self.strands.push_back(PendingStrand { seq, rev, windows: Vec::new(), received: 0, closed: false });
        self.first_idx + self.strands.len() - 1
    }

    /// Reserves the slot of the next window of a strand and returns its index within the strand.
    fn add_window(&mut self, strand_idx: usize) -> usize {
        let windows = &mut self.strands[strand_idx - self.first_idx].windows;
        windows.push(None);
        windows.len() - 1
    }

    /// Marks a strand as complete: no further windows are queued for it.
    fn close(&mut self, strand_idx: usize) {
        self.strands[strand_idx - self.first_idx].closed = true;
    }

    fn store(&mut self, strand_idx: usize, window_idx: usize, window: DecodedWindow) {
        let strand = &mut self.strands[strand_idx - self.first_idx];
        strand.windows[window_idx] = Some(window);
        strand.received += 1;
    }

    /// Removes and returns the front strand if all of its windows are queued and decoded.
    fn pop_ready(&mut self) -> Option<PendingStrand<'s>> {
        let front = self.strands.front()?;
        if !front.closed || front.received < front.windows.len() {
            return None;
        }
        self.first_idx += 1;
        self.strands.pop_front()
    }
}

/// Output side of the pipeline: the GFF and rating writers, the species ratings, the progress
/// lines and the totals of all finalised strands.
pub struct Sink<'w, GW: Write, RW: Write> {
    gff_writer: &'w mut GffWriter<GW>,
    rating_writer: &'w mut RatingWriter<RW>,
    fwd_rating: &'w mut SequenceRating,
    rev_rating: &'w mut SequenceRating,
    progress: &'w mut Progress,
    /// HDF5 blocks of both strands per sequence, indexed by sequence ID, for the progress ETA
    block_counts: &'w [u64],
    /// Number of the next gene within the current sequence (GFF IDs)
    gene_idx: usize,
    pub windows: usize,
    pub bp: usize,
    pub genes: usize,
}

impl<'w, GW: Write, RW: Write> Sink<'w, GW, RW> {
    pub fn new(
        gff_writer: &'w mut GffWriter<GW>,
        rating_writer: &'w mut RatingWriter<RW>,
        fwd_rating: &'w mut SequenceRating,
        rev_rating: &'w mut SequenceRating,
        progress: &'w mut Progress,
        block_counts: &'w [u64],
    ) -> Self {
        Self {
            gff_writer,
            rating_writer,
            fwd_rating,
            rev_rating,
            progress,
            block_counts,
            gene_idx: 1,
            windows: 0,
            bp: 0,
            genes: 0,
        }
    }
}
