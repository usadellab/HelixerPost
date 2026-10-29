//! Progress lines for the per-sequence loop in `main`, measured in processed base pairs.
//!
//! Base pairs rather than sequences keep the number of lines independent of assembly
//! fragmentation: a genome in 50,000 scaffolds gets as many lines as one in 5 chromosomes.

use crate::log_info;
use std::time::{Duration, Instant};

/// Percent of the total bp between two progress lines (at most 20 lines per run).
const PROGRESS_STEP: u64 = 5;

/// Counts processed sequences, base pairs and genes, and prints a progress line each time
/// another `PROGRESS_STEP` percent of the genome is done.
pub struct Progress {
    total_bp: u64,
    total_seqs: usize,
    start: Instant,
    next_percent: u64,
    done_bp: u64,
    done_seqs: usize,
    genes: usize,
}

impl Progress {
    pub fn new(total_bp: u64, total_seqs: usize) -> Self {
        Self {
            total_bp,
            total_seqs,
            start: Instant::now(),
            next_percent: PROGRESS_STEP,
            done_bp: 0,
            done_seqs: 0,
            genes: 0,
        }
    }

    /// Records one finished sequence and prints a progress line if it crossed a milestone.
    pub fn sequence_done(&mut self, length: u64, genes: usize) {
        self.done_bp += length;
        self.done_seqs += 1;
        self.genes += genes;

        let percent = self.done_bp * 100 / self.total_bp.max(1);
        if percent >= self.next_percent {
            self.next_percent = (percent / PROGRESS_STEP + 1) * PROGRESS_STEP;
            self.print(percent);
        }
    }

    fn print(&self, percent: u64) {
        let elapsed = self.start.elapsed();
        let eta = if self.done_bp < self.total_bp {
            let remaining = (self.total_bp - self.done_bp) as f64 / self.done_bp as f64;
            format!(", ETA {}", hms(elapsed.mul_f64(remaining)))
        } else {
            String::new()
        };
        log_info!(
            "Progress: {}% ({:.1}/{:.1} Mbp, {}/{} sequences, {} genes), elapsed {}{}",
            percent,
            self.done_bp as f64 / 1e6,
            self.total_bp as f64 / 1e6,
            self.done_seqs,
            self.total_seqs,
            self.genes,
            hms(elapsed),
            eta
        );
    }

    pub fn elapsed(&self) -> Duration {
        self.start.elapsed()
    }
}

/// Formats a duration as `HH:MM:SS` (hours are not wrapped at 24).
pub fn hms(duration: Duration) -> String {
    let secs = duration.as_secs();
    format!("{:02}:{:02}:{:02}", secs / 3600, secs / 60 % 60, secs % 60)
}
