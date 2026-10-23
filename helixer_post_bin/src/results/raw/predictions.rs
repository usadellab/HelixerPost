use super::super::{Error, Result};
use hdf5_metno::{types::VarLenUnicode, Dataset, File};
use std::path::Path;
use crate::results::error::FromCtx;

pub struct RawHelixerPredictions {
    predictions_file: File,
}

const CLASS_DATASIZE: usize = 4;
const PHASE_DATASIZE: usize = 4;

impl RawHelixerPredictions {
    pub fn new(predictions_file_path: &Path) -> Result<RawHelixerPredictions> {
        let predictions_file = match File::open(predictions_file_path) {
            Ok(f) => f,
            Err(e) => {
                let ctx = format!("opening {predictions_file_path:?} as predictions HDF5");
                return Err(Error::from_ctx(&ctx, e))
            }
        };
        Ok(RawHelixerPredictions { predictions_file })
    }

    pub fn get_model_md5sum(&self) -> Result<String> {
        self.get_string_attr("model_md5sum")
    }

    /// Version of Helixer that wrote the predictions; absent in files from Helixer <= 0.3.7.
    pub fn get_helixer_version(&self) -> Result<String> {
        self.get_string_attr("helixer_version")
    }

    fn get_string_attr(&self, name: &str) -> Result<String> {
        let attr = match self.predictions_file.attr(name) {
            Ok(a) => a,
            Err(e) => {
                let ctx = format!("reading '{name}' attribute from predictions HDF5");
                return Err(Error::from_ctx(&ctx, e))
            }
        };

        let r: VarLenUnicode = match attr.as_reader().read_scalar() {
            Ok(r) => r,
            Err(e) => {
                let ctx = format!("reading '{name}' attribute value from predictions HDF5");
                return Err(Error::from_ctx(&ctx, e))
            }
        };

        Ok(r.to_string())
    }

    fn get_dataset_raw(&self, name: &str) -> Result<Dataset> {
        let dataset = match self.predictions_file.dataset(name) {
            Ok(f) => f,
            Err(e) => {
                let ctx = format!("reading '{name}' dataset from predictions HDF5");
                return Err(Error::from_ctx(&ctx, e))
            }
        };

        Ok(dataset)
    }

    pub fn get_class_raw(&self) -> Result<Dataset> {
        let pred_dataset = self.get_dataset_raw("predictions")?;

        let shape = pred_dataset.shape();
        if shape.len() != 3 {
            return Err(Error::MismatchedDimensions(shape.len(), 3));
        }

        if shape[2] != CLASS_DATASIZE {
            return Err(Error::MismatchedDataSize(shape[2], CLASS_DATASIZE));
        }

        Ok(pred_dataset)
    }

    pub fn get_phase_raw(&self) -> Result<Dataset> {
        let phase_dataset = self.get_dataset_raw("predictions_phase")?;

        let shape = phase_dataset.shape();
        if shape.len() != 3 {
            return Err(Error::MismatchedDimensions(shape.len(), 3));
        }

        let (blocks, blocksize) = self.get_blocks_and_blocksize()?;

        if shape[0] != blocks {
            return Err(Error::MismatchedBlockCount(shape[0], blocks));
        }

        if shape[1] != blocksize {
            return Err(Error::MismatchedDataSize(shape[1], blocksize));
        }

        if shape[2] != PHASE_DATASIZE {
            return Err(Error::MismatchedDataSize(shape[2], PHASE_DATASIZE));
        }

        Ok(phase_dataset)
    }

    pub fn get_blocks_and_blocksize(&self) -> Result<(usize, usize)> {
        let pred_dataset = self.get_class_raw()?;
        let shape = pred_dataset.shape();

        let blocks = shape[0];
        let blocksize = shape[1];

        Ok((blocks, blocksize))
    }
}
