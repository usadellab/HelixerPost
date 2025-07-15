use crate::results::conv::{
    ArrayConvInto, Bases, ClassPrediction, ClassReference, PhasePrediction, PhaseReference,
};
use crate::results::iter::{BlockedDataset2D, BlockedDataset2DIter};
use crate::results::{HelixerResults, Result, SequenceID, SpeciesID};

// Hard coded as Bases / Predictions for now. Would make more sense to have all base-level datasets optionally available (and seekable)
pub struct BasePredictionExtractor<
    'a,
    TC: ArrayConvInto<ClassPrediction>,
    TP: ArrayConvInto<PhasePrediction>,
> {
    helixer_res: &'a HelixerResults,

    bases_blocked_dataset: BlockedDataset2D<'a, f32, Bases>,

    class_pred_blocked_dataset: BlockedDataset2D<'a, TC, ClassPrediction>,
    phase_pred_blocked_dataset: BlockedDataset2D<'a, TP, PhasePrediction>,
}

impl<'a> BasePredictionExtractor<'a, f32, f32> {
    pub fn new_from_prediction(
        helixer_res: &'a HelixerResults,
    ) -> Result<BasePredictionExtractor<'a, f32, f32>> {
        let bases_blocked_dataset = helixer_res.get_x()?;
        let class_pred_blocked_dataset = helixer_res.get_class_predictions()?;
        let phase_pred_blocked_dataset = helixer_res.get_phase_predictions()?;

        Ok(BasePredictionExtractor {
            helixer_res,
            bases_blocked_dataset,
            class_pred_blocked_dataset,
            phase_pred_blocked_dataset,
        })
    }
}

impl<'a> BasePredictionExtractor<'a, i8, i8> {
    pub fn new_from_pseudo_predictions(
        helixer_res: &'a HelixerResults,
    ) -> Result<Option<BasePredictionExtractor<'a, i8, i8>>> {
        let bases_blocked_dataset = helixer_res.get_x()?;

        match (
            helixer_res.get_class_reference_as_pseudo_predictions()?,
            helixer_res.get_phase_reference_as_pseudo_predictions()?,
        ) {
            (Some(class_pred_blocked_dataset), Some(phase_pred_blocked_dataset)) => {
                Ok(Some(BasePredictionExtractor {
                    helixer_res,
                    bases_blocked_dataset,
                    class_pred_blocked_dataset,
                    phase_pred_blocked_dataset,
                }))
            }
            (Some(_), None) => {
                eprintln!("Warning: Ref Class predictions without Ref Phase predictions");
                Ok(None)
            }
            (None, Some(_)) => {
                eprintln!("Warning: Ref Phase predictions without Ref Class predictions");
                Ok(None)
            }
            (None, None) => Ok(None),
        }
    }
}

impl<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>>
    BasePredictionExtractor<'a, TC, TP>
{
    pub fn fwd_iterator(&'a self, sequence_id: SequenceID) -> BasePredictionIterator<'a, TC, TP> {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let base_iter = self.bases_blocked_dataset.fwd_iter(sequence_id);
        let class_pred_iter = self.class_pred_blocked_dataset.fwd_iter(sequence_id);
        let phase_pred_iter = self.phase_pred_blocked_dataset.fwd_iter(sequence_id);

        BasePredictionIterator::new(
            self,
            sequence.get_species_id(),
            sequence_id,
            false,
            base_iter,
            class_pred_iter,
            phase_pred_iter,
        )
    }

    pub fn rev_iterator(&'a self, sequence_id: SequenceID) -> BasePredictionIterator<'a, TC, TP> {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let base_iter = self.bases_blocked_dataset.rev_iter(sequence_id);
        let class_pred_iter = self.class_pred_blocked_dataset.rev_iter(sequence_id);
        let phase_pred_iter = self.phase_pred_blocked_dataset.rev_iter(sequence_id);

        BasePredictionIterator::new(
            self,
            sequence.get_species_id(),
            sequence_id,
            true,
            base_iter,
            class_pred_iter,
            phase_pred_iter,
        )
    }
}

pub struct BasePredictionIterator<
    'a,
    TC: ArrayConvInto<ClassPrediction>,
    TP: ArrayConvInto<PhasePrediction>,
> {
    extractor: &'a BasePredictionExtractor<'a, TC, TP>,

    species_id: SpeciesID,
    sequence_id: SequenceID,
    rc: bool,

    base_iter: BlockedDataset2DIter<'a, f32, Bases>,

    class_pred_iter: BlockedDataset2DIter<'a, TC, ClassPrediction>,
    phase_pred_iter: BlockedDataset2DIter<'a, TP, PhasePrediction>,
}

impl<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>>
    BasePredictionIterator<'a, TC, TP>
{
    fn new(
        extractor: &'a BasePredictionExtractor<'a, TC, TP>,
        species_id: SpeciesID,
        sequence_id: SequenceID,
        rc: bool,
        base_iter: BlockedDataset2DIter<'a, f32, Bases>,
        class_pred_iter: BlockedDataset2DIter<'a, TC, ClassPrediction>,
        phase_pred_iter: BlockedDataset2DIter<'a, TP, PhasePrediction>,
    ) -> BasePredictionIterator<'a, TC, TP> {
        BasePredictionIterator {
            extractor,
            species_id,
            sequence_id,
            rc,
            base_iter,
            class_pred_iter,
            phase_pred_iter,
        }
    }

    pub fn get_extractor(&self) -> &BasePredictionExtractor<TC, TP> {
        self.extractor
    }

    pub fn get_species_id(&self) -> SpeciesID {
        self.species_id
    }

    pub fn get_sequence_id(&self) -> SequenceID {
        self.sequence_id
    }

    pub fn get_rc(&self) -> bool {
        self.rc
    }
}

impl<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>> Iterator
    for BasePredictionIterator<'a, TC, TP>
{
    type Item = (Bases, ClassPrediction, PhasePrediction);

    fn next(&mut self) -> Option<Self::Item> {
        let bases = self.base_iter.next();
        let class_pred = self.class_pred_iter.next();
        let phase_pred = self.phase_pred_iter.next();

        if bases.is_none() && class_pred.is_none() && phase_pred.is_none() {
            return None;
        }

        if bases.is_none() || class_pred.is_none() || phase_pred.is_none() {
            panic!("Different lengths in base/class_pred/phase_pred iterators")
        }

        let bases = bases.expect("Unexpected end of base iter");
        let class_pred = class_pred.expect("Unexpected end of class_pred iter");
        let phase_pred = phase_pred.expect("Unexpected end of phase_pred iter");

        return Some((bases, class_pred, phase_pred));
    }
}

// Hard coded as Class/Phase reference/reference for now. Would make more sense to have all base-level datasets optionally available (and seekable)
pub struct ComparisonExtractor<'a> {
    helixer_res: &'a HelixerResults,

    class_ref_blocked_dataset: Option<BlockedDataset2D<'a, i8, ClassReference>>,
    phase_ref_blocked_dataset: Option<BlockedDataset2D<'a, i8, PhaseReference>>,

    class_pred_blocked_dataset: BlockedDataset2D<'a, f32, ClassPrediction>,
    phase_pred_blocked_dataset: BlockedDataset2D<'a, f32, PhasePrediction>,
}

impl<'a> ComparisonExtractor<'a> {
    pub fn new(helixer_res: &'a HelixerResults) -> Result<ComparisonExtractor<'a>> {
        let class_ref_blocked_dataset = helixer_res.get_class_reference()?;
        let phase_ref_blocked_dataset = helixer_res.get_phase_reference()?;

        let class_pred_blocked_dataset = helixer_res.get_class_predictions()?;
        let phase_pred_blocked_dataset = helixer_res.get_phase_predictions()?;

        Ok(ComparisonExtractor {
            helixer_res,
            class_ref_blocked_dataset,
            phase_ref_blocked_dataset,
            class_pred_blocked_dataset,
            phase_pred_blocked_dataset,
        })
    }

    pub fn fwd_iterator(&'a self, sequence_id: SequenceID) -> ComparisonIterator<'a> {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let class_ref_iter = self
            .class_ref_blocked_dataset
            .as_ref()
            .map(|d| d.fwd_iter(sequence_id));
        let phase_ref_iter = self
            .phase_ref_blocked_dataset
            .as_ref()
            .map(|d| d.fwd_iter(sequence_id));

        let class_pred_iter = self.class_pred_blocked_dataset.fwd_iter(sequence_id);
        let phase_pred_iter = self.phase_pred_blocked_dataset.fwd_iter(sequence_id);

        ComparisonIterator::new(
            self,
            sequence.get_species_id(),
            sequence_id,
            false,
            class_ref_iter,
            phase_ref_iter,
            class_pred_iter,
            phase_pred_iter,
        )
    }

    pub fn rev_iterator(&'a self, sequence_id: SequenceID) -> ComparisonIterator<'a> {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let class_ref_iter = self
            .class_ref_blocked_dataset
            .as_ref()
            .map(|d| d.rev_iter(sequence_id));
        let phase_ref_iter = self
            .phase_ref_blocked_dataset
            .as_ref()
            .map(|d| d.rev_iter(sequence_id));

        let class_pred_iter = self.class_pred_blocked_dataset.rev_iter(sequence_id);
        let phase_pred_iter = self.phase_pred_blocked_dataset.rev_iter(sequence_id);

        ComparisonIterator::new(
            self,
            sequence.get_species_id(),
            sequence_id,
            true,
            class_ref_iter,
            phase_ref_iter,
            class_pred_iter,
            phase_pred_iter,
        )
    }

    pub fn has_ref(&self) -> bool {
        self.class_ref_blocked_dataset.is_some() && self.phase_ref_blocked_dataset.is_some()
    }
}

pub struct ComparisonIterator<'a> {
    extractor: &'a ComparisonExtractor<'a>,

    species_id: SpeciesID,
    sequence_id: SequenceID,
    rc: bool,

    has_ref: bool,

    class_ref_iter: Option<BlockedDataset2DIter<'a, i8, ClassReference>>,
    phase_ref_iter: Option<BlockedDataset2DIter<'a, i8, PhaseReference>>,

    class_pred_iter: BlockedDataset2DIter<'a, f32, ClassPrediction>,
    phase_pred_iter: BlockedDataset2DIter<'a, f32, PhasePrediction>,
}

impl<'a> ComparisonIterator<'a> {
    fn new(
        extractor: &'a ComparisonExtractor<'a>,
        species_id: SpeciesID,
        sequence_id: SequenceID,
        rc: bool,
        class_ref_iter: Option<BlockedDataset2DIter<'a, i8, ClassReference>>,
        phase_ref_iter: Option<BlockedDataset2DIter<'a, i8, PhaseReference>>,
        class_pred_iter: BlockedDataset2DIter<'a, f32, ClassPrediction>,
        phase_pred_iter: BlockedDataset2DIter<'a, f32, PhasePrediction>,
    ) -> ComparisonIterator<'a> {
        let has_ref = extractor.has_ref();
        ComparisonIterator {
            extractor,
            species_id,
            sequence_id,
            rc,
            has_ref,
            class_ref_iter,
            phase_ref_iter,
            class_pred_iter,
            phase_pred_iter,
        }
    }

    pub fn get_extractor(&self) -> &ComparisonExtractor {
        self.extractor
    }

    pub fn get_species_id(&self) -> SpeciesID {
        self.species_id
    }

    pub fn get_sequence_id(&self) -> SequenceID {
        self.sequence_id
    }

    pub fn has_ref(&self) -> bool { self.has_ref }

    pub fn get_rc(&self) -> bool {
        self.rc
    }
}

impl<'a> Iterator for ComparisonIterator<'a> {
    type Item = (
        ClassReference,
        PhaseReference,
        ClassPrediction,
        PhasePrediction,
    );

    fn next(&mut self) -> Option<Self::Item> {
        let class_pred = self.class_pred_iter.next();
        let phase_pred = self.phase_pred_iter.next();

        let (class_ref, phase_ref) = if self.has_ref {
            let class_ref = self.class_ref_iter.as_mut().unwrap().next();
            let phase_ref = self.phase_ref_iter.as_mut().unwrap().next();

            if class_ref.is_none()
                && phase_ref.is_none()
                && class_pred.is_none()
                && phase_pred.is_none()
            {
                return None;
            }

            if class_ref.is_none()
                || phase_ref.is_none()
                || class_pred.is_none()
                || phase_pred.is_none()
            {
                panic!("Different lengths in class_ref/phase_ref/class_pred/phase_pred iterators")
            }

            (class_ref.unwrap(), phase_ref.unwrap())
        } else {
            if class_pred.is_none() && phase_pred.is_none() {
                return None;
            }

            if class_pred.is_none() || phase_pred.is_none() {
                panic!("Different lengths in class_pred/phase_pred iterators")
            }

            (ClassReference::default(), PhaseReference::default())
        };

        let class_pred = class_pred.unwrap();
        let phase_pred = phase_pred.unwrap();

        return Some((class_ref, phase_ref, class_pred, phase_pred));
    }
}
