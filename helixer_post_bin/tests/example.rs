//! End-to-end tests running the binary on the bundled `example/` data.
//!
//! The example data is not curated enough to serve as a reference annotation, so these tests
//! compare runs against each other (config round trip, CLI precedence, thread count and in-flight cap)
//! or check properties of the output (minimum intron length) instead of comparing against committed
//! expected output.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn example_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("example")
}

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_helixer_post_bin")
}

/// Runs helixer_post_bin on the example data with the given extra arguments (`--threads 1` unless
/// given) and returns the GFF, the per-base rating file (only written with `rating`) and stdout.
fn run_all(args: &[&str], rating: bool) -> (String, Option<String>, String) {
    let tmp = tempdir().expect("tempdir");
    let gff = tmp.path().join("out.gff");
    let rating_path = tmp.path().join("rating.tsv");
    let example = example_dir();

    let mut cmd = Command::new(bin());
    if !args.contains(&"--threads") {
        cmd.args(["--threads", "1"]);
    }
    if rating {
        cmd.arg("--rating").arg(&rating_path);
    }
    let output = cmd
        .args(args)
        .arg(example.join("genome_data.h5"))
        .arg(example.join("predictions.h5"))
        .arg(&gff)
        .output()
        .expect("failed to spawn helixer_post_bin");
    assert!(
        output.status.success(),
        "helixer_post_bin exited non-zero: {}\nstderr: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr),
    );
    (
        std::fs::read_to_string(&gff).expect("reading GFF output"),
        rating.then(|| std::fs::read_to_string(&rating_path).expect("reading rating output")),
        String::from_utf8_lossy(&output.stdout).into_owned(),
    )
}

/// Runs helixer_post_bin on the example data with the given extra arguments and returns the GFF.
fn run(args: &[&str]) -> String {
    run_all(args, false).0
}

/// Settings that split the example into 20 candidate windows on both strands instead of 1, so that
/// several windows are decoded at the same time and can finish out of order, and that loosen the HMM
/// (no start/stop codon penalty, high probability floor, no coding length filter) so that several of
/// them yield a gene in the GFF instead of 1 in total.
const MANY_WINDOWS: [&str; 14] = [
    "--window-size", "1",
    "--edge-threshold", "0.1",
    "--peak-threshold", "0.2",
    "--hmm-start-weight", "0",
    "--hmm-stop-weight", "0",
    "--hmm-prob-floor", "0.3",
    "--min-coding-length", "0",
];

fn gene_count(gff: &str) -> usize {
    gff.lines().filter(|l| l.split('\t').nth(2) == Some("gene")).count()
}

/// Returns the decoded window bp and window count from the `Total:` line of stdout
/// (`Total: <genes> genes, <bp>bp across <windows> windows in <time>`).
fn window_totals(stdout: &str) -> (usize, usize) {
    let total = stdout
        .lines()
        .find_map(|l| l.split_once("Total: ").map(|(_, rest)| rest))
        .expect("no Total: line in stdout");
    let words: Vec<&str> = total.split_whitespace().collect();
    let bp = words[2].trim_end_matches("bp").parse().expect("window bp");
    let windows = words[4].parse().expect("window count");
    (bp, windows)
}

/// Returns the lengths of all introns in a GFF, i.e. the gaps between consecutive exons of each
/// transcript.
fn intron_lengths(gff: &str) -> Vec<u64> {
    let mut exons: HashMap<&str, Vec<(u64, u64)>> = HashMap::new();
    for line in gff.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() == 9 && cols[2] == "exon" {
            let parent = cols[8]
                .split(';')
                .find_map(|attr| attr.strip_prefix("Parent="))
                .expect("exon without Parent attribute");
            let start = cols[3].parse().expect("exon start");
            let end = cols[4].parse().expect("exon end");
            exons.entry(parent).or_default().push((start, end));
        }
    }
    exons
        .into_values()
        .flat_map(|mut e| {
            e.sort_unstable();
            e.windows(2).map(|w| w[1].0 - w[0].1 - 1).collect::<Vec<_>>()
        })
        .collect()
}

#[test]
fn default_config_yaml_matches_default_invocation() {
    let cfg = example_dir().join("default_config.yml");
    assert_eq!(
        run(&["--config", cfg.to_str().unwrap()]),
        run(&[]),
        "--config default_config.yml diverged from the built-in defaults"
    );
}

#[test]
fn cli_override_beats_config_file() {
    let tmp = tempdir().expect("tempdir");
    let cfg_path = tmp.path().join("tweak.yml");
    std::fs::write(&cfg_path, "hmm:\n  prob_floor: 0.5\n").unwrap();
    let cfg = cfg_path.to_str().unwrap();
    let default = run(&[]);

    assert_ne!(
        run(&["--config", cfg]),
        default,
        "precondition: prob_floor 0.5 in the config file must change the output"
    );
    assert_eq!(
        run(&["--config", cfg, "--hmm-prob-floor", "1e-9"]),
        default,
        "CLI --hmm-prob-floor should override the config-file value"
    );
}

#[test]
fn min_intron_length_is_enforced() {
    let min_length = 351;

    assert!(
        intron_lengths(&run(&[])).iter().any(|&l| l < min_length),
        "precondition: the default run must contain an intron shorter than {min_length} bp"
    );
    let lengths = intron_lengths(&run(&["--min-intron-length", &min_length.to_string()]));
    assert!(
        lengths.iter().all(|&l| l >= min_length),
        "introns shorter than --min-intron-length {min_length}: {lengths:?}"
    );
}

#[test]
fn thread_count_does_not_change_output() {
    let (gff, rating, stdout) = run_all(&MANY_WINDOWS, true);
    let (_, windows) = window_totals(&stdout);
    assert!(windows > 4, "precondition: more windows ({windows}) than threads (4)");
    assert!(gene_count(&gff) > 1, "precondition: several genes in the GFF");

    let (gff_4, rating_4, _) = run_all(&[&MANY_WINDOWS[..], &["--threads", "4"]].concat(), true);
    assert_eq!(gff_4, gff, "--threads 4 changed the GFF");
    assert_eq!(rating_4, rating, "--threads 4 changed the rating");
}

#[test]
fn max_in_flight_mbp_does_not_change_output() {
    let args = [&MANY_WINDOWS[..], &["--threads", "4"]].concat();
    let (gff, rating, stdout) = run_all(&args, true);
    let (bp, _) = window_totals(&stdout);
    assert!(bp > 2_000, "precondition: the windows ({bp} bp) exceed the 0.002 Mbp cap");
    assert!(gene_count(&gff) > 1, "precondition: several genes in the GFF");

    // 1 bp: one window at a time; 0.002 Mbp: some windows wait for others to finish
    for cap in ["0.000001", "0.002"] {
        let (gff_cap, rating_cap, _) = run_all(&[&args[..], &["--max-in-flight-mbp", cap]].concat(), true);
        assert_eq!(gff_cap, gff, "--max-in-flight-mbp {cap} changed the GFF");
        assert_eq!(rating_cap, rating, "--max-in-flight-mbp {cap} changed the rating");
    }
}
