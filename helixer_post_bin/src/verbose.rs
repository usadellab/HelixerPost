//! Process-wide verbosity switch and timestamped log macros.
//!
//! The verbosity flag is set once in `main` from the `--verbose` flag and read from the
//! analysis code, including the HMM worker threads, without passing a flag through every
//! call. Log lines follow the Helixer log format (`2026-01-31 12:00:00, I: message`, local
//! time), so HelixerPost output reads consistently inside a Helixer.py log.

use std::fmt::Display;
use std::sync::atomic::{AtomicBool, Ordering};

static VERBOSE: AtomicBool = AtomicBool::new(false);

pub fn set_verbose(verbose: bool) {
    VERBOSE.store(verbose, Ordering::Relaxed);
}

pub fn is_verbose() -> bool {
    VERBOSE.load(Ordering::Relaxed)
}

/// Current local time in the Helixer log format.
pub fn timestamp() -> impl Display {
    chrono::Local::now().format("%Y-%m-%d %H:%M:%S")
}

/// `println!` prefixed with a timestamp and the info level.
#[macro_export]
macro_rules! log_info {
    ($($arg:tt)*) => {
        println!("{}, I: {}", $crate::verbose::timestamp(), format_args!($($arg)*))
    };
}

/// `eprintln!` prefixed with a timestamp and the warning level.
#[macro_export]
macro_rules! log_warn {
    ($($arg:tt)*) => {
        eprintln!("{}, W: {}", $crate::verbose::timestamp(), format_args!($($arg)*))
    };
}

/// `println!` prefixed with a timestamp and the debug level, only printed when verbose
/// output is enabled.
#[macro_export]
macro_rules! vprintln {
    ($($arg:tt)*) => {
        if $crate::verbose::is_verbose() {
            println!("{}, D: {}", $crate::verbose::timestamp(), format_args!($($arg)*));
        }
    };
}
