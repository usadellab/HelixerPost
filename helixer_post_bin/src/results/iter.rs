use super::BlockID;
use hdf5_metno::{Dataset, H5Type};

use ndarray::{s, Array1, Array2, ArrayView1, Axis};

use super::HelixerIndex;
use crate::results::conv::ArrayConvInto;
use crate::results::SequenceID;
use std::marker::PhantomData;
use std::ops::Range;

// Should be possible to merge the 1D / 2D versions at some stage with a non-trivial amount of generic magic

pub struct BlockedDataset1D<'a, T: H5Type + Clone + Copy> {
    index: &'a HelixerIndex,
    dataset: Dataset,
    phantom: PhantomData<T>,
}

impl<'a, T: H5Type + Clone + Copy> BlockedDataset1D<'a, T> {
    pub fn new(index: &'a HelixerIndex, dataset: Dataset) -> BlockedDataset1D<'a, T> {
        BlockedDataset1D {
            index,
            dataset,
            phantom: PhantomData,
        }
    }

    pub fn get_index(&self) -> &HelixerIndex {
        self.index
    }

    pub fn get_dataset(&self) -> &Dataset {
        &self.dataset
    }

    pub fn fwd_iter(&'a self, id: SequenceID) -> BlockedDataset1DIter<'a, T> {
        let (fwd, _rev) = self.index.get_block_ids_for_sequence(id);
        let block_offsets = self.index.get_all_block_offsets();

        BlockedDataset1DIter::new(&self, block_offsets, fwd)
    }

    pub fn rev_iter(&'a self, id: SequenceID) -> BlockedDataset1DIter<'a, T> {
        let (_fwd, rev) = self.index.get_block_ids_for_sequence(id);
        let block_offsets = self.index.get_all_block_offsets();

        BlockedDataset1DIter::new(&self, block_offsets, rev)
    }

    fn get_data_for_block(&self, block_id: BlockID) -> hdf5_metno::Result<Array1<T>> {
        let slice = s![block_id.inner(), ..];
        self.dataset.read_slice_1d(slice)
    }
}

pub struct BlockedDataset1DIter<'a, T: H5Type + Clone + Copy> {
    blocked_dataset: &'a BlockedDataset1D<'a, T>,
    block_offsets: &'a [(u64, u64)],
    block_iter: std::slice::Iter<'a, BlockID>,

    block: Option<(BlockID, Array1<T>, Range<usize>)>,
}

impl<'a, T: H5Type + Clone + Copy> BlockedDataset1DIter<'a, T> {
    fn new(
        blocked_dataset: &'a BlockedDataset1D<'a, T>,
        block_offsets: &'a [(u64, u64)],
        blocks: &'a [BlockID],
    ) -> BlockedDataset1DIter<'a, T> {
        let block_iter = blocks.iter();
        let mut iter = BlockedDataset1DIter {
            blocked_dataset,
            block_offsets,
            block_iter,
            block: None,
        };
        iter.block = iter.next_block();
        iter
    }

    fn next_block(&mut self) -> Option<(BlockID, Array1<T>, Range<usize>)> {
        self.block_iter.next().map(|id| {
            let data = self
                .blocked_dataset
                .get_data_for_block(*id)
                .expect("Failed to read block");
            let (start, end) = self.block_offsets[id.inner()];
            let length = if start < end {
                end - start
            } else {
                start - end
            } as usize;

            (*id, data, 0..length)
        })
    }

    fn next_position(&mut self) -> (bool, Option<usize>) {
        if let Some((_, _, ref mut block_range)) = self.block {
            (true, block_range.next())
        } else {
            (false, None)
        }
    }
}

impl<'a, T: H5Type + Clone + Copy> Iterator for BlockedDataset1DIter<'a, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        let (_, maybe_offset) = match self.next_position() {
            (true, None) => {
                self.block = self.next_block();
                self.next_position()
            }
            keep => keep,
        };

        if let (Some(offset), Some((_, block_array, _))) = (maybe_offset, self.block.as_ref()) {
            Some(block_array[offset])
        } else {
            None
        }
    }
}

pub struct BlockedDataset2D<'a, T: ArrayConvInto<O>, O> {
    index: &'a HelixerIndex,
    dataset: Dataset,
    phantom: PhantomData<(T, O)>,
}

impl<'a, T: ArrayConvInto<O>, O> BlockedDataset2D<'a, T, O> {
    pub fn new(index: &'a HelixerIndex, dataset: Dataset) -> BlockedDataset2D<'a, T, O> {
        BlockedDataset2D {
            index,
            dataset,
            phantom: PhantomData,
        }
    }

    pub fn get_index(&self) -> &HelixerIndex {
        self.index
    }

    pub fn get_dataset(&self) -> &Dataset {
        &self.dataset
    }

    pub fn fwd_iter(&'a self, id: SequenceID) -> BlockedDataset2DIter<'a, T, O> {
        let (fwd, _rev) = self.index.get_block_ids_for_sequence(id);
        let block_offsets = self.index.get_all_block_offsets();

        BlockedDataset2DIter::new(&self, block_offsets, fwd)
    }

    pub fn rev_iter(&'a self, id: SequenceID) -> BlockedDataset2DIter<'a, T, O> {
        let (_fwd, rev) = self.index.get_block_ids_for_sequence(id);
        let block_offsets = self.index.get_all_block_offsets();

        BlockedDataset2DIter::new(&self, block_offsets, rev)
    }

    fn get_data_for_block(&self, block_id: BlockID) -> hdf5_metno::Result<Array2<T>> {
        let slice = s![block_id.inner(), .., ..];
        self.dataset.read_slice_2d(slice)
    }
}

pub struct BlockedDataset2DIter<'a, T: ArrayConvInto<O>, O> {
    blocked_dataset: &'a BlockedDataset2D<'a, T, O>,
    block_offsets: &'a [(u64, u64)],
    block_iter: std::slice::Iter<'a, BlockID>,

    block: Option<(BlockID, Array2<T>, Range<usize>)>,
}

impl<'a, T: ArrayConvInto<O>, O> BlockedDataset2DIter<'a, T, O> {
    fn new(
        blocked_dataset: &'a BlockedDataset2D<'a, T, O>,
        block_offsets: &'a [(u64, u64)],
        blocks: &'a [BlockID],
    ) -> BlockedDataset2DIter<'a, T, O> {
        let block_iter = blocks.iter();
        let mut iter = BlockedDataset2DIter {
            blocked_dataset,
            block_offsets,
            block_iter,
            block: None,
        };
        iter.block = iter.next_block();
        iter
    }

    fn next_block(&mut self) -> Option<(BlockID, Array2<T>, Range<usize>)> {
        self.block_iter.next().map(|id| {
            let data = self
                .blocked_dataset
                .get_data_for_block(*id)
                .expect("Failed to read block");
            let (start, end) = self.block_offsets[id.inner()];
            let length = if start < end {
                end - start
            } else {
                start - end
            } as usize;

            (*id, data, 0..length)
        })
    }

    fn next_position(&mut self) -> (bool, Option<usize>) {
        if let Some((_, _, ref mut block_range)) = self.block {
            (true, block_range.next())
        } else {
            (false, None)
        }
    }
}

impl<'a, T: ArrayConvInto<O>, O> Iterator for BlockedDataset2DIter<'a, T, O> {
    type Item = O;

    fn next(&mut self) -> Option<Self::Item> {
        let (_, maybe_offset) = match self.next_position() {
            (true, None) => {
                self.block = self.next_block();
                self.next_position()
            }
            keep => keep,
        };

        if let (Some(offset), Some((_, block_array, _))) = (maybe_offset, self.block.as_ref()) {
            let slice: ArrayView1<T> = block_array.index_axis(Axis(0), offset); //.to_owned();
            Some(T::into(slice))
        } else {
            None
        }
    }
}
