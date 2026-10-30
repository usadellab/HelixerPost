# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

HelixerPost is the Rust post-processing stage for the Helixer gene-prediction pipeline (https://github.com/weberlab-hhu/Helixer). It reads HDF5 input (the preprocessed genome and base-wise predictions from Helixer's neural network) and produces a GFF gene annotation plus a rating/comparison output.

## Build & run

Single Cargo workspace with one member crate (`helixer_post_bin`).

```
cargo build --release          # release binary at target/release/helixer_post_bin
cargo build                    # debug
cargo test                     # unit tests in src/ plus end-to-end tests in tests/example.rs
cargo test <name>              # run a single test by substring match
```

`tests/example.rs` runs the binary on `example/` and compares runs against each other or checks output properties (e.g. minimum intron length), each with an explicit precondition so a test cannot pass vacuously. There are deliberately no committed expected-output files: the example data is not curated enough to serve as a reference annotation.

System dependency: HDF5 development headers (`hdf5-devel` on Fedora, `libhdf5-dev` on Ubuntu). The crate uses `hdf5-metno` (a fork of the abandoned `hdf5` crate) for compatibility with current native HDF5 versions.

### Current CLI signature

The binary uses `clap` derive (`Cli` struct in `src/main.rs`). Three positional args are required; everything else is a flag:

```
helixer_post_bin [OPTIONS] <genome.h5> <predictions.h5> <output.gff>
```

Key flags:

- `-c, --config <path>` — YAML config file (full schema in `example/default_config.yml`).
- `--print-default-config` — dump built-in defaults to stdout as YAML, then exit. Positional args are not required when this flag is present.
- `--rating <path>` — per-base rating file, only written when given (`RatingWriter::disabled()` otherwise skips all per-base formatting); roughly 230 bytes per bp.
- `-v, --verbose` — enables the per-sequence/per-strand rating dumps, the per-sequence and per-window progress lines and the per-region HMM traces. Without it, stdout holds the HMM/pipeline config, the progress lines, the species total rating and the `Total:` line (genes, bp, windows, run time).- `-j, --threads <N>` — defaults to `std::thread::available_parallelism()`.
- `--window-size`, `--edge-threshold`, `--peak-threshold`, `--min-coding-length` — pass-1/3 overrides.
- `--hmm-prob-floor`, `--hmm-phase-retain`, `--hmm-{start,stop,donor,acceptor}-weight` — pass-2 HMM tuning overrides.
- `--min-intron-length` (sets `u2_gt_ag` and `u2_gc_ag`; U12-type GT-AG introns are decoded as U2 and follow it too) and `--min-intron-length-u12` (sets `u12_at_ac`) — pass-2 minimum intron lengths, enforced inside the decoder (floor 2).

Precedence: **CLI flag > `--config` file value > built-in default**. The splice-flag booleans and the three per-donor fixed penalties are file-only (no CLI flags); set them via the YAML config.

The Helixer Python wrapper used to call this binary with 9 positional args; that contract is broken and the upstream invocation needs to be updated to the flag form. A worked example with input/output files lives in `example/example.md`, and the example folder also contains `default_config.yml` which round-trips byte-identical output via `--config`.

## Architecture

The pipeline is roughly: HDF5 → blocked datasets → per-sequence iterators → sliding-window candidate detection → per-window HMM → GFF records + rating stats. The crate is split into three top-level modules: `results` (input), `analysis` (the pipeline), and `gff` (output formatting).

### `results/` — HDF5 input layer

- `results/raw/{genome,predictions}.rs`: thin wrappers around the two HDF5 files. `RawHelixerPredictions` reads `predictions.h5` (class + phase probabilities, `model_md5sum` and `helixer_version` attributes; the latter is absent in files from Helixer <= 0.3.7); `RawHelixerGenome` reads `genome_data.h5` (one-hot bases X, optional reference Y/phases, sample weights, species/seqid metadata).
- `results/index.rs`: `HelixerIndex` walks the genome metadata once at startup, builds `Species`/`Sequence`/`BlockID` tables, and produces forward/reverse block-ID lists per sequence. Both block count and blocksize come from the predictions file (`get_blocks_and_blocksize`), and the genome is opened against those dimensions — predictions and genome must agree.
- `results/conv.rs`: `ArrayConvInto` and the typed wrappers (`Bases`, `ClassPrediction`, `PhasePrediction`, `ClassReference`, `PhaseReference`, `Transitions`). The same physical dataset can be viewed as predictions (`f32`) or pseudo-predictions derived from references (`i8`) — this is why extractors and most generics are parameterized on `TC: ArrayConvInto<ClassPrediction>` and `TP: ArrayConvInto<PhasePrediction>`.
- `results/iter.rs`: `BlockedDataset1D`/`BlockedDataset2D` provide block-by-block streaming over the typed datasets without loading the full HDF5 file into memory. Blocks are read through the `BlockElement` trait (a supertrait of `ArrayConvInto`): `f32` datasets stored as float16 (the Helixer genome X and both prediction datasets) are read as `half::f16` and widened in Rust, because the HDF5 1.10 float16 to float32 conversion dominated the runtime (about 240 vs 36 ns per bp and dataset).
- `results::HelixerResults` is the façade that ties raw files + index together and is what `main` and the extractors see.

### `analysis/` — the pipeline

The pipeline is split into three passes, each with a corresponding config struct (all defined in `analysis.rs` or `analysis/hmm.rs`):

- **Pass 1: sliding-window candidate detection** — driven by `WindowConfig` (`window_size`, `edge_threshold`, `peak_threshold`).
- **Pass 2: per-window HMM decode** — driven by `HmmConfig` (`prob_floor`, `phase_retain`, plus `SpliceFlags`, `HmmWeights`, `DonorFixedPenalties` sub-structs).
- **Pass 3: post-HMM gene filtering** — driven by `FilterConfig` (`min_coding_length`; expected to grow).

`FileConfig { window, filter, hmm }` is the top-level YAML shape. All four structs and their sub-structs use `#[serde(default, deny_unknown_fields)]` so partial YAML overlays the built-in `Default` impl and typos surface as parse errors. Each struct also has a `with_overrides(Option<T>, ...)` that maps the CLI's `Option<T>` flags onto the value. `main.rs` builds the final configs as `file_section.with_overrides(cli_flags)`, giving CLI > file > default precedence in one line per section.

The dataflow per sequence/strand is built up as nested iterators:

1. `extractor::BasePredictionExtractor::{fwd,rev}_iterator(seq_id)` yields per-base `(Bases, ClassPrediction, PhasePrediction)` triples for a given sequence and strand by zipping the bases / class / phase blocked datasets. `ComparisonExtractor` does the same but also includes reference labels when present, for rating.
2. `window::BasePredictionWindowThresholdIterator` wraps the per-base iterator into a sliding-window scanner: it accumulates the genic score over `window_size` bases, finds runs where the mean exceeds `edge_threshold`, and emits only the runs whose peak exceeds `peak_threshold`. Each emitted item is one candidate gene-containing region.
3. `Analyzer::process_sequence_1d` (in `analysis.rs`) dispatches each candidate window to a `threadpool::ThreadPool` worker. The analyzer holds an `Arc<HmmConfig>`; each worker gets a clone before `thread_pool.execute(move || ...)`, and the worker calls `PredictionHmm::new(bp_vec, hmm_cfg).solve()`. Results come back via an `mpsc` channel and are re-sorted into emission order before downstream processing.
4. `hmm::PredictionHmm` runs a Viterbi-style decode using per-base penalties (`ClassPredPenalty`, `PhasePredPenalty`, `PredPenalty`, `BasesPenalty`). These are now built via explicit `::new(pred, prob_floor[, phase_retain])` constructors (formerly `impl From`) so the floor and phase-blending parameters flow from `HmmConfig`. `TransitionContext` carries `cfg: &HmmConfig`; the donor/acceptor methods and `HmmState::get_common_state_entrance_penalty` read weights and fixed penalties off it, and `populate_successor_states_and_transition_penalties` reads splice booleans from `cfg.splice.*` at each `consider_transition` call site. Output is a sequence of `HmmStateRegion`s, which `HmmStateRegion::split_genes` groups into individual genes.
5. `gff_conv::hmm_solution_to_gff` converts each gene into GFF records honouring the strand and the active `min_coding_length`; records are streamed out via `gff::GffWriter`.
6. In parallel with GFF emission, `rater::SequenceRater` consumes the `ComparisonExtractor` iterator and tallies a confusion matrix per base class, producing precision/recall/F1 per sequence and per species; these are dumped to stdout (per sequence and per strand only with `--verbose`) and to the rating output file via `rater::RatingWriter`.

Things to keep in mind when editing the pipeline:

- Anything user-tunable in pass 1, pass 2, or pass 3 belongs in `WindowConfig` / `HmmConfig` / `FilterConfig` respectively. Plain `const` is reserved for structural numbers (`HMM_STATES`, `PENALTY_SCALE`, `MAX_EVALS`) and the deliberately-narrow `APPROX_ZERO` heuristic in `BasesPenalty::as_str`. When adding a new tuning knob, add it to the relevant config struct's `Default` impl and `with_overrides` signature, then expose a `--…` flag in `main.rs::Cli` if you want CLI access; YAML access is automatic via the serde derive.
- The HMM is the expensive step and runs on a worker pool sized by the CLI `--threads` flag. Reordering or re-buffering must preserve candidate-window emission order (`process_sequence_1d` does this via a pre-allocated `results` vec keyed by `window_count` index).
- Forward and reverse strands are processed independently with separate ratings; both feed the same `GffWriter`. The `gene_idx` counter is shared so gene IDs stay unique within a sequence.
- `main.rs` rejects multi-species inputs with a clear `bail!` for GFF output (was an `assert!` pre-Phase A).
- Errors propagate through `process_sequence` / `process_sequence_1d` as `anyhow::Result`; the previous `panic!("No solution at …")` is now a `bail!`. Note: the workspace still has `panic = "abort"`, so a worker-thread panic (e.g. `MAX_EVALS` exceeded inside `PredictionHmm::solve`) will still abort the whole process — those paths haven't been migrated to `Result`.
- Verbosity is a process-wide flag in `verbose.rs` (`set_verbose` in `main`, `is_verbose()`), so worker threads and the HMM trace can read it without parameter threading. Log lines use the macros in `verbose.rs`, which prefix the local time and a level in the Helixer log format (`%Y-%m-%d %H:%M:%S, I: ...`, via `chrono`): `log_info!` (stdout), `log_warn!` (stderr, regardless of verbosity) and `vprintln!` (stdout, `D`, only with `--verbose`). Multi-line dumps (rating tables, HMM config details) stay unprefixed and are wrapped in `if is_verbose()` where verbose-only. `progress::Progress` is fed by the sequence loop in `main` and prints a progress line per `PROGRESS_STEP` (5) percent of the genome in bp; there are no flags for it. Its ETA is the average of a bp-based and an HDF5-block-based extrapolation (block counts per sequence from the index), since short scaffolds cost more than their bp share and less than their block share.
- `analysis::hmm::show_hmm_config(&HmmConfig)` is called once at startup so the active HMM tuning is visible in any run log; `main` also prints the pipeline (window/filter) config on the next line.

### `gff.rs` — output

Stateless GFF3 writer. `write_global_header` emits the `##gff-version`, optional `##species`, the `model_md5sum` comment from the predictions HDF5, and `# Helixer version <v>` (from the `helixer_version` attribute, `unknown` if absent) plus `# HelixerPost version <crate version>` comments; `write_region_header` emits one `##sequence-region` per input sequence; `write_records` streams a `Vec<GffRecord>` for a single gene group.
