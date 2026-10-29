//! End-to-end tests running the binary on the bundled `example/` data.
//!
//! The example data is not curated enough to serve as a reference annotation, so these tests
//! compare runs against each other (config round trip, CLI precedence) or check properties of the
//! output (minimum intron length) instead of comparing against committed expected output.

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

/// Runs helixer_post_bin on the example data with the given extra arguments and returns the GFF.
fn run(args: &[&str]) -> String {
    let tmp = tempdir().expect("tempdir");
    let gff = tmp.path().join("out.gff");
    let example = example_dir();

    let output = Command::new(bin())
        .args(["--threads", "1"])
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
    std::fs::read_to_string(&gff).expect("reading GFF output")
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
