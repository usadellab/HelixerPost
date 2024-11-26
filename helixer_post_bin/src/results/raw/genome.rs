use super::super::{Error, Result};
use hdf5_metno::types::FixedAscii;
use hdf5_metno::{Dataset, File};
use std::path::Path;

//use ndarray::iter::Iter;

pub struct RawHelixerGenome {
    genome_file: File,

    blocks: usize,
    blocksize: usize,
}

const X_DATASIZE: usize = 4;
const PHASE_DATASIZE: usize = 4;
const STARTENDS_DATASIZE: usize = 2;
const TRANSITIONS_DATASIZE: usize = 6;
const Y_DATASIZE: usize = 4;

/*
#[derive(hdf5::H5Type, Clone, PartialEq)]
#[repr(i8)]
pub enum ErrSamples
{
    FALSE = 0,
    TRUE = 1
}
*/

type SeqidsType = FixedAscii<50>;
type SpeciesType = FixedAscii<25>;

impl RawHelixerGenome {
    pub fn new(
        genome_file_path: &Path,
        blocks: usize,
        blocksize: usize,
    ) -> Result<RawHelixerGenome> {
        let genome_file = File::open(genome_file_path)?;
        Ok(RawHelixerGenome {
            genome_file,
            blocks,
            blocksize,
        })
    }

    // Should be 1D - [Blocks]
    fn validate_dataset_shape_scalar(&self, dataset: &Dataset) -> Result<()> {
        let shape = dataset.shape();
        if shape.len() != 1 {
            return Err(Error::MismatchedDimensions(shape.len(), 1));
        }

        if shape[0] != self.blocks {
            return Err(Error::MismatchedBlockCount(shape[0], self.blocks));
        }

        Ok(())
    }

    // Should be 2D - [Blocks][size]
    fn validate_dataset_shape_array(&self, dataset: &Dataset, data_size: usize) -> Result<()> {
        let shape = dataset.shape();
        if shape.len() != 2 {
            return Err(Error::MismatchedDimensions(shape.len(), 2));
        }

        if shape[0] != self.blocks {
            return Err(Error::MismatchedBlockCount(shape[0], self.blocks));
        }

        if shape[1] != data_size {
            return Err(Error::MismatchedDataSize(shape[2], data_size));
        }

        Ok(())
    }

    // Should be 2D - [Blocks][Blocksize]
    fn validate_dataset_shape_blocksize_scalar(&self, dataset: &Dataset) -> Result<()> {
        let shape = dataset.shape();
        if shape.len() != 2 {
            return Err(Error::MismatchedDimensions(shape.len(), 2));
        }

        if shape[0] != self.blocks {
            return Err(Error::MismatchedBlockCount(shape[0], self.blocks));
        }

        if shape[1] != self.blocksize {
            return Err(Error::MismatchedBlockSize(shape[1], self.blocksize));
        }

        Ok(())
    }

    // Should be 3D - [Blocks][Blocksize][size]
    fn validate_dataset_shape_blocksize_array(
        &self,
        dataset: &Dataset,
        data_size: usize,
    ) -> Result<()> {
        let shape = dataset.shape();
        if shape.len() != 3 {
            return Err(Error::MismatchedDimensions(shape.len(), 3));
        }

        if shape[0] != self.blocks {
            return Err(Error::MismatchedBlockCount(shape[0], self.blocks));
        }

        if shape[1] != self.blocksize {
            return Err(Error::MismatchedBlockSize(shape[1], self.blocksize));
        }

        if shape[2] != data_size {
            return Err(Error::MismatchedDataSize(shape[2], data_size));
        }

        Ok(())
    }

    pub fn get_x_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/X")?;
        self.validate_dataset_shape_blocksize_array(&dataset, X_DATASIZE)?;
        Ok(dataset)
    }

    pub fn get_err_samples_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/err_samples")?;
        self.validate_dataset_shape_scalar(&dataset)?;
        Ok(dataset)
    }

    pub fn get_err_samples(&self) -> Result<Vec<bool>> {
        let err_samples_ds = self.get_err_samples_raw()?;
        let err_samples_array = err_samples_ds.read_1d::<bool>()?;
        Ok(err_samples_array.to_vec())

        //      let err_samples_array = err_samples_ds.read_1d::<ErrSamples>()?;
        //      Ok(err_samples_array.iter().map(|x| *x == ErrSamples::TRUE).collect())
    }

    pub fn get_fully_intergenic_samples_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/fully_intergenic_samples")?;
        self.validate_dataset_shape_scalar(&dataset)?;
        Ok(dataset)
    }

    pub fn get_fully_intergenic_samples(&self) -> Result<Vec<bool>> {
        let fully_ig_samples_ds = self.get_fully_intergenic_samples_raw()?;
        let fully_ig_samples_array = fully_ig_samples_ds.read_1d::<bool>()?;
        Ok(fully_ig_samples_array.to_vec())
    }

    // For pre-annotated
    pub fn get_gene_lengths_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/gene_lengths")?;
        self.validate_dataset_shape_blocksize_scalar(&dataset)?;
        Ok(dataset)
    }

    // For pre-annotated
    pub fn get_is_annotated_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/is_annotated")?;
        self.validate_dataset_shape_scalar(&dataset)?;
        Ok(dataset)
    }

    pub fn get_is_annotated(&self) -> Result<Vec<bool>> {
        let is_annotated_ds = self.get_is_annotated_raw()?;
        let is_annotated_array = is_annotated_ds.read_1d::<bool>()?;
        Ok(is_annotated_array.to_vec())
    }

    // For pre-annotated
    pub fn get_phases_raw(&self) -> Result<Option<Dataset>> {
        if !self.genome_file.link_exists("data/phases") {
            return Ok(None);
        }

        let dataset = self.genome_file.dataset("data/phases")?;
        self.validate_dataset_shape_blocksize_array(&dataset, PHASE_DATASIZE)?;
        Ok(Some(dataset))
    }

    pub fn get_sample_weights_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/sample_weights")?;
        self.validate_dataset_shape_blocksize_scalar(&dataset)?;
        Ok(dataset)
    }

    pub fn get_seqids_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/seqids")?;
        self.validate_dataset_shape_scalar(&dataset)?;
        Ok(dataset)
    }

    pub fn get_seqids(&self) -> Result<Vec<String>> {
        let seqids_ds = self.get_seqids_raw()?;
        let seqids_array = seqids_ds.read_1d::<SeqidsType>()?;

        Ok(seqids_array.iter().map(|x| x.to_string()).collect())
    }

    pub fn get_species_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/species")?;
        self.validate_dataset_shape_scalar(&dataset)?;
        Ok(dataset)
    }

    pub fn get_species(&self) -> Result<Vec<String>> {
        let species_ds = self.get_species_raw()?;
        let species_array = species_ds.read_1d::<SpeciesType>()?;

        Ok(species_array.iter().map(|x| x.to_string()).collect())
    }

    pub fn get_start_ends_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/start_ends")?;
        self.validate_dataset_shape_array(&dataset, STARTENDS_DATASIZE)?;
        Ok(dataset)
    }

    pub fn get_start_ends(&self) -> Result<Vec<(u64, u64)>> {
        let startends_ds = self.get_start_ends_raw()?;
        let startends_array = startends_ds.read_2d::<i64>()?;

        Ok(startends_array
            .outer_iter()
            .map(|x| (x[0] as u64, x[1] as u64))
            .collect())
    }

    // For pre-annotated
    pub fn get_transitions_raw(&self) -> Result<Dataset> {
        let dataset = self.genome_file.dataset("data/transitions")?;
        self.validate_dataset_shape_blocksize_array(&dataset, TRANSITIONS_DATASIZE)?;
        Ok(dataset)
    }

    // For pre-annotated
    pub fn get_y_raw(&self) -> Result<Option<Dataset>> {
        if !self.genome_file.link_exists("data/y") {
            return Ok(None);
        }

        let dataset = self.genome_file.dataset("data/y")?;
        self.validate_dataset_shape_blocksize_array(&dataset, Y_DATASIZE)?;
        Ok(Some(dataset))
    }
}
