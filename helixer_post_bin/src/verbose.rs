//! Process-wide verbosity switch for detailed progress and debug output.
//!
//! Set once in `main` from the `--verbose` flag; read from the analysis code, including
//! the HMM worker threads, without passing a flag through every call.

use std::sync::atomic::{AtomicBool, Ordering};

static VERBOSE: AtomicBool = AtomicBool::new(false);

pub fn set_verbose(verbose: bool) {
    VERBOSE.store(verbose, Ordering::Relaxed);
}

pub fn is_verbose() -> bool {
    VERBOSE.load(Ordering::Relaxed)
}

/// `println!` that only prints when verbose output is enabled.
#[macro_export]
macro_rules! vprintln {
    ($($arg:tt)*) => {
        if $crate::verbose::is_verbose() {
            println!($($arg)*);
        }
    };
}
