//! End-to-end smoke tests using the bundled `example/` data.
//!
//! These exercise the full pipeline against the committed `example/output.gff3`
//! and `example/output.txt`. A diff lands here whenever any part of the pipeline
//! changes the gene model or rating numbers — this is the project's safety net
//! for the dozens of behaviours that aren't covered by per-module unit tests.

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

fn run(args: &[&Path]) -> (PathBuf, PathBuf) {
    let tmp = tempdir().expect("tempdir");
    let gff = tmp.path().join("out.gff");
    let rating = tmp.path().join("out.rating");

    let mut cmd = Command::new(bin());
    cmd.arg("--threads")
        .arg("1")
        .arg("--rating")
        .arg(&rating);
    for a in args {
        cmd.arg(a);
    }
    let example = example_dir();
    cmd.arg(example.join("genome_data.h5"))
        .arg(example.join("predictions.h5"))
        .arg(&gff);

    let output = cmd.output().expect("failed to spawn helixer_post_bin");
    assert!(
        output.status.success(),
        "helixer_post_bin exited non-zero: {}\nstderr: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr),
    );

    // Move the produced files out of the tempdir so they survive teardown.
    let gff_persist = tmp.path().join("out.gff.keep");
    let rating_persist = tmp.path().join("out.rating.keep");
    std::fs::rename(&gff, &gff_persist).unwrap();
    std::fs::rename(&rating, &rating_persist).unwrap();
    // Leak the tempdir so the kept paths stay valid for the assertions below.
    let leaked = tmp.keep();
    (leaked.join("out.gff.keep"), leaked.join("out.rating.keep"))
}

#[test]
fn default_invocation_matches_committed_output() {
    let (gff, rating) = run(&[]);
    let expected_gff = std::fs::read_to_string(example_dir().join("output.gff3")).unwrap();
    let expected_rating = std::fs::read_to_string(example_dir().join("output.txt")).unwrap();
    let actual_gff = std::fs::read_to_string(&gff).unwrap();
    let actual_rating = std::fs::read_to_string(&rating).unwrap();
    assert_eq!(actual_gff, expected_gff, "GFF differs from example/output.gff3");
    assert_eq!(
        actual_rating, expected_rating,
        "rating differs from example/output.txt"
    );
    let _ = std::fs::remove_file(&gff);
    let _ = std::fs::remove_file(&rating);
}

#[test]
fn default_config_yaml_matches_default_invocation() {
    let cfg = example_dir().join("default_config.yml");
    let (gff, _rating) = run(&[Path::new("--config"), &cfg]);
    let expected = std::fs::read_to_string(example_dir().join("output.gff3")).unwrap();
    let actual = std::fs::read_to_string(&gff).unwrap();
    assert_eq!(actual, expected, "--config default_config.yml diverged from default");
    let _ = std::fs::remove_file(&gff);
}

#[test]
fn cli_override_beats_config_file() {
    // YAML sets a deliberately disruptive prob_floor; the CLI flag below puts
    // it back to the built-in default, so the GFF must match the committed one.
    let tmp = tempdir().expect("tempdir");
    let cfg_path = tmp.path().join("tweak.yml");
    std::fs::write(&cfg_path, "hmm:\n  prob_floor: 0.5\n").unwrap();

    let (gff, _rating) = run(&[
        Path::new("--config"),
        &cfg_path,
        Path::new("--hmm-prob-floor"),
        Path::new("1e-9"),
    ]);
    let expected = std::fs::read_to_string(example_dir().join("output.gff3")).unwrap();
    let actual = std::fs::read_to_string(&gff).unwrap();
    assert_eq!(
        actual, expected,
        "CLI --hmm-prob-floor should override config-file value"
    );
    let _ = std::fs::remove_file(&gff);
}
