use hdf5_metno::H5Type;
use ndarray::ArrayView1;

pub trait ArrayConvFrom<T: Sized + H5Type> {
    fn from(array: ArrayView1<'_, T>) -> Self;
}

pub trait ArrayConvInto<T>: Sized + H5Type {
    fn into(array: ArrayView1<'_, Self>) -> T;
}

impl<T: Sized + H5Type, U> ArrayConvInto<U> for T
where
    U: ArrayConvFrom<T>,
{
    fn into(array: ArrayView1<'_, Self>) -> U {
        U::from(array)
    }
}

#[derive(Clone, Copy)]
pub struct ClassPrediction {
    values: [f32; 4], // Ordering is intergenic, utr, coding, intron
}

impl ClassPrediction {
    pub fn get(&self) -> &[f32; 4] {
        &self.values
    }

    pub fn get_intergenic(&self) -> f32 {
        self.values[0]
    }

    pub fn get_utr(&self) -> f32 {
        self.values[1]
    }

    pub fn get_coding(&self) -> f32 {
        self.values[2]
    }

    pub fn get_intron(&self) -> f32 {
        self.values[3]
    }

    pub fn get_genic(&self) -> f32 {
        1.0 - self.values[0]
    }

    pub fn get_max_idx(&self) -> usize {
        let mut max = self.values[0];
        let mut max_idx = 0;

        for i in 1..4 {
            if self.values[i] > max {
                max = self.values[i];
                max_idx = i;
            }
        }

        max_idx
    }
}

impl ArrayConvFrom<f32> for ClassPrediction {
    fn from(array: ArrayView1<'_, f32>) -> Self {
        let values: [f32; 4] = [array[0], array[1], array[2], array[3]];
        ClassPrediction { values }
    }
}

impl ArrayConvFrom<i8> for ClassPrediction {
    fn from(array: ArrayView1<'_, i8>) -> Self {
        let values: [f32; 4] = [
            array[0] as f32,
            array[1] as f32,
            array[2] as f32,
            array[3] as f32,
        ];
        ClassPrediction { values }
    }
}

#[derive(Clone, Copy)]
pub struct PhasePrediction {
    values: [f32; 4], // Ordering is non_coding, phase 0, phase 1, phase 2
}

impl PhasePrediction {
    pub fn get(&self) -> &[f32; 4] {
        &self.values
    }

    pub fn get_non_coding(&self) -> f32 {
        self.values[0]
    }

    pub fn get_phase0(&self) -> f32 {
        self.values[1]
    }

    pub fn get_phase1(&self) -> f32 {
        self.values[2]
    }

    pub fn get_phase2(&self) -> f32 {
        self.values[3]
    }

    pub fn get_max_idx(&self) -> usize {
        let mut max = self.values[0];
        let mut max_idx = 0;

        for i in 1..4 {
            if self.values[i] > max {
                max = self.values[i];
                max_idx = i;
            }
        }

        max_idx
    }
}

impl ArrayConvFrom<f32> for PhasePrediction {
    fn from(array: ArrayView1<'_, f32>) -> Self {
        let values: [f32; 4] = [array[0], array[1], array[2], array[3]];
        PhasePrediction { values }
    }
}

impl ArrayConvFrom<i8> for PhasePrediction {
    fn from(array: ArrayView1<'_, i8>) -> Self {
        let values: [f32; 4] = [
            array[0] as f32,
            array[1] as f32,
            array[2] as f32,
            array[3] as f32,
        ];
        PhasePrediction { values }
    }
}

#[derive(Clone, Copy)]
pub struct Bases {
    values: [f32; 4], // Ordering is C, A, T, G
}

impl Bases {
    pub fn get(&self) -> &[f32; 4] {
        &self.values
    }
}

impl ArrayConvFrom<f32> for Bases {
    fn from(array: ArrayView1<'_, f32>) -> Self {
        let values: [f32; 4] = [array[0], array[1], array[2], array[3]];
        Bases { values }
    }
}

pub struct Transitions {
    values: [i8; 6], // Ordering is +TR, +CDS, +In, -TR, -CDS, -In
                     // [transcription start site, start codon, donor splice site, transcription stop site, stop codon, acceptor splice site]
}

impl Transitions {
    pub fn get(&self) -> &[i8; 6] {
        &self.values
    }
}

impl ArrayConvFrom<i8> for Transitions {
    fn from(array: ArrayView1<'_, i8>) -> Self {
        let values: [i8; 6] = [array[0], array[1], array[2], array[3], array[4], array[5]];
        Transitions { values }
    }
}

pub struct PhaseReference {
    values: [i8; 4], // Ordering is Non-Coding, Phase 0, Phase 1, Phase2
}

impl PhaseReference {
    pub fn get(&self) -> &[i8; 4] {
        &self.values
    }

    pub fn get_max_idx(&self) -> usize {
        let mut max = self.values[0];
        let mut max_idx = 0;

        for i in 1..4 {
            if self.values[i] > max {
                max = self.values[i];
                max_idx = i;
            }
        }

        max_idx
    }
}

impl Default for PhaseReference {
    fn default() -> Self {
        let values: [i8; 4] = [i8::MAX, 0, 0, 0];
        PhaseReference { values }
    }
}

impl ArrayConvFrom<i8> for PhaseReference {
    fn from(array: ArrayView1<'_, i8>) -> Self {
        let values: [i8; 4] = [array[0], array[1], array[2], array[3]];
        PhaseReference { values }
    }
}

pub struct ClassReference {
    values: [i8; 4], // Ordering is intergenic, utr, coding, intron
}

impl ClassReference {
    pub fn get(&self) -> &[i8; 4] {
        &self.values
    }

    pub fn get_max_idx(&self) -> usize {
        let mut max = self.values[0];
        let mut max_idx = 0;

        for i in 1..4 {
            if self.values[i] > max {
                max = self.values[i];
                max_idx = i;
            }
        }

        max_idx
    }
}

impl Default for ClassReference {
    fn default() -> Self {
        let values: [i8; 4] = [i8::MAX, 0, 0, 0];
        ClassReference { values }
    }
}

impl ArrayConvFrom<i8> for ClassReference {
    fn from(array: ArrayView1<'_, i8>) -> Self {
        let values: [i8; 4] = [array[0], array[1], array[2], array[3]];
        ClassReference { values }
    }
}
