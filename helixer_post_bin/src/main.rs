use anyhow::{bail, Context, Result};
use clap::Parser;
use helixer_post_bin::analysis::extractor::{BasePredictionExtractor, ComparisonExtractor};
use helixer_post_bin::analysis::hmm::show_hmm_config;
use helixer_post_bin::analysis::rater::{RatingWriter, SequenceRating};
use helixer_post_bin::analysis::{Analyzer, FileConfig};
use helixer_post_bin::gff::GffWriter;
use helixer_post_bin::results::raw::RawHelixerPredictions;
use helixer_post_bin::results::HelixerResults;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "helixer_post_bin",
    about = "Post-process Helixer base-level predictions into a GFF gene annotation.",
    long_about = "Reads the HDF5 genome and prediction files produced by the Helixer pipeline, \
                  runs sliding-window candidate detection followed by a per-window HMM decode, \
                  and writes a GFF3 annotation plus a precision/recall rating file."
)]
struct Cli {
    /// HDF5 genome data input (as produced by fasta2h5py). Not required with --print-default-config.
    #[arg(required_unless_present = "print_default_config")]
    genome: Option<PathBuf>,

    /// HDF5 base-level predictions input. Not required with --print-default-config.
    #[arg(required_unless_present = "print_default_config")]
    predictions: Option<PathBuf>,

    /// Output GFF3 path. Not required with --print-default-config.
    #[arg(required_unless_present = "print_default_config")]
    gff: Option<PathBuf>,

    /// YAML config file. CLI flags override values from this file; missing
    /// fields fall back to built-in defaults. Use --print-default-config to
    /// dump a starting template.
    #[arg(long, short = 'c')]
    config: Option<PathBuf>,

    /// Print the default config as YAML and exit. Useful for generating a
    /// starting template (e.g. `helixer_post_bin --print-default-config > tune.yml`).
    #[arg(long, exclusive = true)]
    print_default_config: bool,

    /// Output rating / stats path. Defaults to <gff>.rating.
    #[arg(long)]
    rating: Option<PathBuf>,

    /// Worker threads for HMM evaluation (defaults to available parallelism).
    #[arg(long, short = 'j')]
    threads: Option<usize>,

    /// Sliding-window width in bases (pass 1).
    #[arg(long)]
    window_size: Option<usize>,

    /// Mean genic score required to enter/leave a candidate window (pass 1).
    #[arg(long)]
    edge_threshold: Option<f32>,

    /// Peak genic score required to accept a candidate window (pass 1).
    #[arg(long)]
    peak_threshold: Option<f32>,

    /// Drop genes whose total coding length is below this (pass 3).
    #[arg(long)]
    min_coding_length: Option<usize>,

    // ---- HMM tuning overrides (pass 2). Full splice/penalty tables come from the config file once Phase C lands. ----
    /// Floor applied to raw model probabilities before negative-log conversion.
    #[arg(long)]
    hmm_prob_floor: Option<f64>,

    /// Weight given to the predicted phase signal (0.0 = ignore, 1.0 = trust).
    #[arg(long)]
    hmm_phase_retain: Option<f64>,

    /// Multiplicative weight on the start codon base-match penalty.
    #[arg(long)]
    hmm_start_weight: Option<f64>,

    /// Multiplicative weight on the stop codon base-match penalty.
    #[arg(long)]
    hmm_stop_weight: Option<f64>,

    /// Multiplicative weight on the splice donor base-match penalty.
    #[arg(long)]
    hmm_donor_weight: Option<f64>,

    /// Multiplicative weight on the splice acceptor base-match penalty.
    #[arg(long)]
    hmm_acceptor_weight: Option<f64>,
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn default_rating_path(gff: &PathBuf) -> PathBuf {
    let mut p = gff.clone();
    let new_name = match p.file_name() {
        Some(name) => {
            let mut s = name.to_os_string();
            s.push(".rating");
            s
        }
        None => "helixer.rating".into(),
    };
    p.set_file_name(new_name);
    p
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.print_default_config {
        print!("{}", FileConfig::default().to_yaml()?);
        return Ok(());
    }

    // required_unless_present above guarantees these three are Some when we reach this point.
    let genome_path = cli.genome.expect("genome path required");
    let predictions_path = cli.predictions.expect("predictions path required");
    let gff_path = cli.gff.expect("gff path required");

    let file_cfg = match &cli.config {
        Some(path) => FileConfig::load_yaml(path)?,
        None => FileConfig::default(),
    };

    let window_cfg = file_cfg.window.with_overrides(
        cli.window_size,
        cli.edge_threshold,
        cli.peak_threshold,
    );
    let filter_cfg = file_cfg.filter.with_overrides(cli.min_coding_length);
    let hmm_cfg = file_cfg.hmm.with_overrides(
        cli.hmm_prob_floor,
        cli.hmm_phase_retain,
        cli.hmm_start_weight,
        cli.hmm_stop_weight,
        cli.hmm_donor_weight,
        cli.hmm_acceptor_weight,
    );
    let thread_count = cli.threads.unwrap_or_else(default_threads);
    let rating_path = cli.rating.unwrap_or_else(|| default_rating_path(&gff_path));

    let helixer_res = HelixerResults::new(&predictions_path, &genome_path)
        .context("opening Helixer prediction / genome HDF5 inputs")?;

    let bp_extractor = BasePredictionExtractor::new_from_prediction(&helixer_res)
        .context("opening Base / ClassPrediction / PhasePrediction datasets")?;

    let comp_extractor = ComparisonExtractor::new(&helixer_res).context(
        "opening ClassReference / PhaseReference / ClassPrediction / PhasePrediction datasets",
    )?;

    let analyzer = Analyzer::new(
        bp_extractor,
        comp_extractor,
        window_cfg,
        filter_cfg,
        hmm_cfg,
        thread_count,
    );

    show_hmm_config(&hmm_cfg);
    println!(
        "Pipeline config: window_size={} edge_threshold={} peak_threshold={} min_coding_length={} threads={}",
        window_cfg.window_size,
        window_cfg.edge_threshold,
        window_cfg.peak_threshold,
        filter_cfg.min_coding_length,
        thread_count,
    );

    let gff_file = File::create(&gff_path)
        .with_context(|| format!("creating GFF output {}", gff_path.display()))?;
    let mut gff_writer = GffWriter::new(BufWriter::new(gff_file));

    let all_species = helixer_res.get_all_species();
    if all_species.len() != 1 {
        bail!(
            "Expected exactly one species for GFF output, found {}",
            all_species.len()
        );
    }

    let rhg = RawHelixerPredictions::new(&predictions_path)
        .context("re-opening predictions HDF5 for model_md5sum lookup")?;
    let model_md5sum = rhg.get_model_md5sum().ok();
    let species_name = all_species.first().map(|x| x.get_name());

    gff_writer
        .write_global_header(species_name, model_md5sum)
        .with_context(|| format!("writing GFF header to {}", gff_path.display()))?;

    let rating_file = File::create(&rating_path)
        .with_context(|| format!("creating rating output {}", rating_path.display()))?;
    let mut rating_writer = RatingWriter::new(BufWriter::new(rating_file));

    let mut total_count = 0;
    let mut total_length = 0;

    for species in helixer_res.get_all_species() {
        let mut fwd_species_rating = SequenceRating::new();
        let mut rev_species_rating = SequenceRating::new();

        let id = species.get_id();
        println!(
            "Sequences for Species {} - {}",
            species.get_name(),
            id.inner()
        );
        for seq_id in helixer_res.get_sequences_for_species(id) {
            let seq = helixer_res.get_sequence_by_id(*seq_id);
            gff_writer
                .write_region_header(seq.get_name(), seq.get_length())
                .with_context(|| {
                    format!("writing sequence-region header for {}", seq.get_name())
                })?;

            let (count, length) = analyzer.process_sequence(
                species,
                seq,
                &mut fwd_species_rating,
                &mut rev_species_rating,
                &mut gff_writer,
                &mut rating_writer,
            )?;

            total_count += count;
            total_length += length;
        }

        println!(
            "Forward for Species {} - {}",
            species.get_name(),
            id.inner()
        );
        fwd_species_rating.dump(analyzer.has_ref());

        println!(
            "Reverse for Species {} - {}",
            species.get_name(),
            id.inner()
        );
        rev_species_rating.dump(analyzer.has_ref());

        let mut species_rating = SequenceRating::new();
        species_rating.accumulate(&fwd_species_rating);
        species_rating.accumulate(&rev_species_rating);

        println!("Total for Species {} - {}", species.get_name(), id.inner());
        species_rating.dump(analyzer.has_ref());
    }

    println!("Total: {}bp across {} windows", total_length, total_count);
    Ok(())
}
