use std::fmt::{self, Debug, Display};

pub enum Error {
    Hdf5(hdf5_metno::Error),
    MismatchedDimensions(usize, usize),
    MismatchedBlockCount(usize, usize),
    MismatchedBlockSize(usize, usize),
    MismatchedDataSize(usize, usize),
    InvalidValue(String),
    DuplicateValue(String),
}

impl Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)
    }
}

impl Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)
    }
}

impl From<hdf5_metno::Error> for Error {
    fn from(hdf5_error: hdf5_metno::Error) -> Self {
        Self::Hdf5(hdf5_error)
    }
}

impl From<&str> for Error {
    fn from(msg: &str) -> Self {
        Self::InvalidValue(msg.to_string())
    }
}

impl Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Hdf5(err) => Display::fmt(err, f),
            Error::MismatchedDimensions(found, expected) => write!(
                f,
                "Mismatched Dimensions Found {}, expected {}",
                found, expected
            ),
            Error::MismatchedBlockCount(found, expected) => write!(
                f,
                "Mismatched Block Count: Found {}, expected {}",
                found, expected
            ),
            Error::MismatchedBlockSize(found, expected) => write!(
                f,
                "Mismatched Block Size: Found {}, expected {}",
                found, expected
            ),
            Error::MismatchedDataSize(found, expected) => write!(
                f,
                "Mismatched Data Size: Found {}, expected {}",
                found, expected
            ),
            Error::InvalidValue(msg) => write!(f, "Invalid Value: {}", msg),
            Error::DuplicateValue(msg) => write!(f, "Duplicate Value: {}", msg),
        }
    }

}
