use std::io::{BufWriter, Write};
use crate::analysis::extractor::ComparisonIterator;
use crate::analysis::hmm::{HmmAnnotationLabel, HmmStateRegion};
use crate::results::conv::{ClassPrediction, ClassReference, PhasePrediction, PhaseReference};
use crate::results::{Sequence, Species};
/*
    Rows are based on 'reference', Columns are based on 'prediction'

    Precision = TP / ( TP + FP )
    Recall = TP / ( TP + FN )

    F1 = (2 * Precision * Recall) / (Precision + Recall)
    F1 = TP / ( TP + ( FP + FN ) / 2 )
    F1 = (2 * TP) / (2 * TP + FP + FN)
*/

fn calc_precision_recall_f1(true_pos: u64, false_pos: u64, false_neg: u64) -> (f64, f64, f64) {
    let true_pos = true_pos as f64;
    let false_pos = false_pos as f64;
    let false_neg = false_neg as f64;

    let precision = true_pos / (true_pos + false_pos);
    let recall = true_pos / (true_pos + false_neg);
    let f1 = (2.0 * true_pos) / (2.0 * true_pos + false_pos + false_neg);

    (precision, recall, f1)
}

#[derive(Copy, Clone)]
pub struct ConfusionMatrix<const N: usize> {
    count: [[u64; N]; N],
}

impl<const N: usize> ConfusionMatrix<N> {
    //const SMALLER: usize = N - 1;
    //const BIGGER: usize = N + 1;

    pub fn new() -> ConfusionMatrix<N> {
        let count = [[0; N]; N];

        ConfusionMatrix { count }
    }

    pub fn size(&self) -> usize {
        N
    }

    pub fn accumulate(&mut self, other: &Self) {
        for p in 0..N {
            for r in 0..N {
                self.count[r][p] += other.count[r][p];
            }
        }
    }

    /*
        fn set(&mut self, ref_idx: usize, pred_idx: usize, count: u64)
        {
            self.count[ref_idx][pred_idx]=count;
        }
    */

    pub fn increment(&mut self, ref_idx: usize, pred_idx: usize) {
        self.count[ref_idx][pred_idx] += 1;
    }

    pub fn get_tp(&self, idx: usize) -> u64 {
        self.count[idx][idx]
    }

    pub fn get_fp(&self, idx: usize) -> u64 {
        let mut false_pos = 0;
        for r in 0..idx {
            false_pos += self.count[r][idx];
        }

        for r in idx + 1..N {
            false_pos += self.count[r][idx];
        }

        false_pos
    }

    pub fn get_fn(&self, idx: usize) -> u64 {
        let mut false_neg = 0;
        for p in 0..idx {
            false_neg += self.count[idx][p];
        }

        for p in idx + 1..N {
            false_neg += self.count[idx][p];
        }

        false_neg
    }

    pub fn get_precision_recall_f1(&self, idx: usize) -> (f64, f64, f64) {
        calc_precision_recall_f1(self.get_tp(idx), self.get_fp(idx), self.get_fn(idx))
    }
}

#[derive(Clone, Copy)]
enum Annotation {
    OutsideWindow,
    Filtered,
    Intergenic,
    UTR,
    CodingPhase0,
    CodingPhase1,
    CodingPhase2,
    Intron,
}

impl Annotation {
    pub fn to_str(&self) -> &str {
        match self {
            Annotation::OutsideWindow => "OutsideWindow",
            Annotation::Filtered => "Filtered",
            Annotation::Intergenic => "Intergenic",
            Annotation::UTR => "UTR",
            Annotation::CodingPhase0 => "CodingPhase0",
            Annotation::CodingPhase1 => "CodingPhase1",
            Annotation::CodingPhase2 => "CodingPhase2",
            Annotation::Intron => "Intron",
        }
    }

    fn get_class_idx(self) -> usize {
        match self {
            Annotation::OutsideWindow => 0,
            Annotation::Filtered => 0,
            Annotation::Intergenic => 0,
            Annotation::UTR => 1,
            Annotation::CodingPhase0 => 2,
            Annotation::CodingPhase1 => 2,
            Annotation::CodingPhase2 => 2,
            Annotation::Intron => 3,
        }
    }

    fn get_phase_idx(self) -> usize {
        match self {
            Annotation::OutsideWindow => 0,
            Annotation::Filtered => 0,
            Annotation::Intergenic => 0,
            Annotation::UTR => 0,
            Annotation::CodingPhase0 => 1,
            Annotation::CodingPhase1 => 2,
            Annotation::CodingPhase2 => 3,
            Annotation::Intron => 0,
        }
    }
}

impl ClassReference {
    fn get_class_idx(&self) -> usize {
        self.get_max_idx()
    }
}

impl PhaseReference {
    fn get_phase_idx(&self) -> usize {
        self.get_max_idx()
    }
}

impl ClassPrediction {
    fn get_class_idx(&self) -> usize {
        self.get_max_idx()
    }
}

impl PhasePrediction {
    fn get_phase_idx(&self) -> usize {
        self.get_max_idx()
    }
}

pub struct SequenceRater<'a> {
    comp_iterator: ComparisonIterator<'a>,
    annotation: Vec<Annotation>,
}

impl<'a> SequenceRater<'a> {
    pub fn new(comp_iterator: ComparisonIterator<'a>, seq_length: usize) -> SequenceRater<'a> {
        let mut annotation = Vec::with_capacity(seq_length);
        annotation.resize(seq_length, Annotation::OutsideWindow);

        SequenceRater {
            comp_iterator,
            annotation,
        }
    }

    pub fn rate_window_regions(
        &mut self,
        start_offset: usize,
        gene_regions: &[HmmStateRegion],
        filtered: bool,
    ) {
        if filtered {
            for region in gene_regions {
                for pos in
                    region.get_start_pos() + start_offset..region.get_end_pos() + start_offset
                {
                    self.annotation[pos] = Annotation::Filtered;
                }
            }
        } else {
            let mut coding_annotation = Annotation::CodingPhase0;

            for region in gene_regions {
                println!("Rating region {} {} {}", region.get_start_pos(), region.get_end_pos(), region.get_annotation_label().to_str());

                if region.get_annotation_label() != HmmAnnotationLabel::Coding {
                    let annotation = match region.get_annotation_label() {
                        HmmAnnotationLabel::Intergenic => Annotation::Intergenic,
                        HmmAnnotationLabel::UTR5 => Annotation::UTR,
                        HmmAnnotationLabel::Intron => Annotation::Intron,
                        HmmAnnotationLabel::UTR3 => Annotation::UTR,
                        _ => panic!(
                            "Unexpected Hmm Annotation Label {}",
                            region.get_annotation_label().to_str()
                        ),
                    };

                    for pos in
                        region.get_start_pos() + start_offset..region.get_end_pos() + start_offset
                    {
                        self.annotation[pos] = annotation;
                    }
                } else {
                    for pos in
                        region.get_start_pos() + start_offset..region.get_end_pos() + start_offset
                    {
                        self.annotation[pos] = coding_annotation;

                        coding_annotation = match coding_annotation {
                            Annotation::CodingPhase0 => Annotation::CodingPhase2,
                            Annotation::CodingPhase1 => Annotation::CodingPhase0,
                            Annotation::CodingPhase2 => Annotation::CodingPhase1,
                            _ => panic!(
                                "Unexpected Coding Annotation Label {}",
                                coding_annotation.to_str()
                            ),
                        }
                    }
                }
            }
        }
    }

    pub fn calculate_stats<RW: Write>(self, species: &Species, seq: &Sequence, rev: bool, rating_writer: &mut RatingWriter<RW>) -> std::io::Result<SequenceRating> {
        let mut rating = SequenceRating::new();
        rating.rate(self.comp_iterator, self.annotation, species, seq, rev, rating_writer)?;

        Ok(rating)
    }
}

// Window States = Outside / Inside * Intergenic / Genic
// Class = Intergenic / UTR / Coding / Intron
// Phase = NonCoding / Ph0 / Ph1 / Ph2

pub struct SequenceRating {
    ref_ml_class_confusion: ConfusionMatrix<4>,
    ref_ml_phase_confusion: ConfusionMatrix<4>,

    ref_hp_class_confusion: ConfusionMatrix<4>,
    ref_hp_phase_confusion: ConfusionMatrix<4>,

    ml_hp_class_confusion: ConfusionMatrix<4>,
    ml_hp_phase_confusion: ConfusionMatrix<4>,

    outside_window_count: u64,
    filtered_count: u64,
}

impl SequenceRating {
    pub fn new() -> SequenceRating {
        SequenceRating {
            ref_ml_class_confusion: ConfusionMatrix::<4>::new(),
            ref_ml_phase_confusion: ConfusionMatrix::<4>::new(),
            ref_hp_class_confusion: ConfusionMatrix::<4>::new(),
            ref_hp_phase_confusion: ConfusionMatrix::<4>::new(),
            ml_hp_class_confusion: ConfusionMatrix::<4>::new(),
            ml_hp_phase_confusion: ConfusionMatrix::<4>::new(),
            outside_window_count: 0,
            filtered_count: 0,
        }
    }

    pub fn accumulate(&mut self, other: &Self) {
        self.ref_ml_class_confusion
            .accumulate(&other.ref_ml_class_confusion);
        self.ref_ml_phase_confusion
            .accumulate(&other.ref_ml_phase_confusion);

        self.ref_hp_class_confusion
            .accumulate(&other.ref_hp_class_confusion);
        self.ref_hp_phase_confusion
            .accumulate(&other.ref_hp_phase_confusion);

        self.ml_hp_class_confusion
            .accumulate(&other.ml_hp_class_confusion);
        self.ml_hp_phase_confusion
            .accumulate(&other.ml_hp_phase_confusion);

        self.outside_window_count += other.outside_window_count;
        self.filtered_count += other.filtered_count;
    }

    fn rate<RW: Write>(&mut self, comp_iterator: ComparisonIterator, annotation: Vec<Annotation>,
                       species: &Species, seq: &Sequence, rev: bool, rating_writer: &mut RatingWriter<RW>) -> std::io::Result<()> {
        let has_ref = comp_iterator.has_ref();

        for (idx, ((class_ref, phase_ref, class_ml, phase_ml), annotation)) in
            comp_iterator.zip(annotation.into_iter()).enumerate()
        {
            let ref_class_idx = class_ref.get_class_idx();
            let ref_phase_idx = phase_ref.get_phase_idx();

            let ml_class_idx = class_ml.get_class_idx();
            let ml_phase_idx = phase_ml.get_phase_idx();

            let hp_class_idx = annotation.get_class_idx();
            let hp_phase_idx = annotation.get_phase_idx();

            self.ref_ml_class_confusion
                .increment(ref_class_idx, ml_class_idx);
            self.ref_ml_phase_confusion
                .increment(ref_phase_idx, ml_phase_idx);

            self.ref_hp_class_confusion
                .increment(ref_class_idx, hp_class_idx);
            self.ref_hp_phase_confusion
                .increment(ref_phase_idx, hp_phase_idx);

            self.ml_hp_class_confusion
                .increment(ml_class_idx, hp_class_idx);
            self.ml_hp_phase_confusion
                .increment(ml_phase_idx, hp_phase_idx);

            if ref_class_idx != 0 {
                match annotation {
                    Annotation::OutsideWindow => self.outside_window_count += 1,
                    Annotation::Filtered => self.filtered_count += 1,
                    _ => (),
                }
            }

            if has_ref {
                rating_writer.write(species, seq, rev, idx, Some(&class_ref), Some(&phase_ref), &class_ml, &phase_ml, annotation)?;
            }
            else {
                rating_writer.write(species, seq, rev, idx, None, None, &class_ml, &phase_ml, annotation)?;
            }
        }

        Ok(())
    }

    const CLASS_NAMES: [&'static str; 4] = ["Intergenic", "UTR", "Coding", "Intron"];
    const PHASE_NAMES: [&'static str; 4] = ["Non Coding", "Phase 0", "Phase 1", "Phase 2"];

    pub fn show_confusion_matrices(
        class_confusion: &ConfusionMatrix<4>,
        phase_confusion: &ConfusionMatrix<4>,
        corner_label: &str,
    ) {
        /* Confusion Matrices, two column layout */
        print!("{:>10} Class\t", corner_label);
        for j in 0..4 {
            print!("{:>12}\t", Self::CLASS_NAMES[j]);
        }

        print!("\t");

        print!("{:>10} Phase\t", corner_label);
        for j in 0..4 {
            print!("{:>12}\t", Self::PHASE_NAMES[j]);
        }

        println!();

        for i in 0..4 {
            print!("{:>16}\t", Self::CLASS_NAMES[i]);
            for j in 0..4 {
                print!("{:12}\t", class_confusion.count[i][j]);
            }

            print!("\t");

            print!("{:>16}\t", Self::PHASE_NAMES[i]);
            for j in 0..4 {
                print!("{:12}\t", phase_confusion.count[i][j]);
            }

            println!();
        }

        println!();

        /* Summaries, two column layout */

        print!(
            "{:>16}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\t",
            "", "Precision", "Recall", "F1", ""
        );
        print!("\t");
        println!(
            "{:>16}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\t",
            "", "Precision", "Recall", "F1", ""
        );

        for i in 0..4 {
            let (prec, rec, f1) = class_confusion.get_precision_recall_f1(i);
            print!(
                "{:>16}\t{:12.5}\t{:12.5}\t{:12.5}\t{:>12}\t",
                Self::CLASS_NAMES[i],
                prec,
                rec,
                f1,
                ""
            );

            print!("\t");

            let (prec, rec, f1) = phase_confusion.get_precision_recall_f1(i);
            print!(
                "{:>16}\t{:12.5}\t{:12.5}\t{:12.5}\t{:>12}\t",
                Self::PHASE_NAMES[i],
                prec,
                rec,
                f1,
                ""
            );

            println!();
        }

        println!();

        let subg_true_pos = class_confusion.get_tp(2) + class_confusion.get_tp(3);
        let subg_false_pos = class_confusion.get_fp(2) + class_confusion.get_fp(3);
        let subg_false_neg = class_confusion.get_fn(2) + class_confusion.get_fn(3);

        let (subg_prec, subg_rec, subg_f1) =
            calc_precision_recall_f1(subg_true_pos, subg_false_pos, subg_false_neg);

        let gen_true_pos = subg_true_pos + class_confusion.get_tp(1);
        let gen_false_pos = subg_false_pos + class_confusion.get_fp(1);
        let gen_false_neg = subg_false_neg + class_confusion.get_fn(1);

        let (gen_prec, gen_rec, gen_f1) =
            calc_precision_recall_f1(gen_true_pos, gen_false_pos, gen_false_neg);

        let coding_true_pos =
            phase_confusion.get_tp(1) + phase_confusion.get_tp(2) + phase_confusion.get_tp(3);
        let coding_false_pos =
            phase_confusion.get_fp(1) + phase_confusion.get_fp(2) + phase_confusion.get_fp(3);
        let coding_false_neg =
            phase_confusion.get_fn(1) + phase_confusion.get_fn(2) + phase_confusion.get_fn(3);

        let (coding_prec, coding_rec, coding_f1) =
            calc_precision_recall_f1(coding_true_pos, coding_false_pos, coding_false_neg);

        print!(
            "{:>16}\t{:12.5}\t{:12.5}\t{:12.5}\t{:>12}\t",
            "Subgenic", subg_prec, subg_rec, subg_f1, ""
        );
        print!("\t");
        println!(
            "{:>16}\t{:12.5}\t{:12.5}\t{:12.5}\t{:>12}\t",
            "Coding", coding_prec, coding_rec, coding_f1, ""
        );

        println!(
            "{:>16}\t{:12.5}\t{:12.5}\t{:12.5}\t",
            "Genic", gen_prec, gen_rec, gen_f1
        );

        println!();
        println!();
    }

    pub fn dump(&self, has_ref: bool) {
        if has_ref {
            println!(
                "Lost Ref Genic: Outside Window {}, Filtered {}",
                self.outside_window_count, self.filtered_count
            );
            println!();

            Self::show_confusion_matrices(
                &self.ref_ml_class_confusion,
                &self.ref_ml_phase_confusion,
                "Ref v ML",
            );
            Self::show_confusion_matrices(
                &self.ref_hp_class_confusion,
                &self.ref_hp_phase_confusion,
                "Ref v HP",
            );
        }
        Self::show_confusion_matrices(
            &self.ml_hp_class_confusion,
            &self.ml_hp_phase_confusion,
            "ML v HP",
        );
    }
}


pub struct RatingWriter<W: Write> {
    writer: BufWriter<W>,
}

const CLASS_NAMES: [&str; 4] = ["IG", "UTR", "CDS", "INT"];
const PHASE_NAMES: [&str; 4] = ["NP", "P0", "P1", "P2"];

impl<W: Write> RatingWriter<W> {
    pub fn new(writer: BufWriter<W>) -> RatingWriter<W> {
        RatingWriter { writer }
    }

    fn write(&mut self, _species: &Species, seq: &Sequence, rev: bool, pos: usize,
                 class_ref: Option<&ClassReference>, phase_ref: Option<&PhaseReference>,
                 class_ml: &ClassPrediction, phase_ml: &PhasePrediction,
                 annotation: Annotation) -> std::io::Result<()> {
//write!(self.writer, "# {}\n", helixer_model_md5sum)?; println!("{:.3}", x)

        let rev = if rev {"-"} else {"+"};

        let (class_ref_str, class_ref_scores) = if let Some(class_ref) = class_ref {
            let class_ref_str = CLASS_NAMES[class_ref.get_class_idx()];
            let class_ref_vals = class_ref.get();
            let class_ref_scores = format!("{},{},{},{}", class_ref_vals[0], class_ref_vals[1], class_ref_vals[2], class_ref_vals[3]);
            (class_ref_str, class_ref_scores)
        }
        else {
            ("-", String::from("-,-,-,-"))
        };

        let class_ml_str = CLASS_NAMES[class_ml.get_class_idx()];
        let class_ml_vals = class_ml.get();
        let class_ml_scores = format!("{:.3},{:.3},{:.3},{:.3}", class_ml_vals[0], class_ml_vals[1], class_ml_vals[2], class_ml_vals[3]);

        let (phase_ref_str, phase_ref_scores) = if let Some(phase_ref) = phase_ref {
            let phase_ref_str = PHASE_NAMES[phase_ref.get_phase_idx()];
            let phase_ref_vals = phase_ref.get();
            let phase_ref_scores = format!("{},{},{},{}", phase_ref_vals[0], phase_ref_vals[1], phase_ref_vals[2], phase_ref_vals[3]);
            (phase_ref_str, phase_ref_scores)
        } else {
            ("-", String::from("-,-,-,-"))
        };

        let phase_ml_str = PHASE_NAMES[phase_ml.get_phase_idx()];
        let phase_ml_vals = phase_ml.get();
        let phase_ml_scores = format!("{:.3},{:.3},{:.3},{:.3}", phase_ml_vals[0], phase_ml_vals[1], phase_ml_vals[2], phase_ml_vals[3]);

        writeln!(self.writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
               seq.get_name(), rev, pos,
               class_ref_str, class_ref_scores, class_ml_str, class_ml_scores,
               phase_ref_str, phase_ref_scores, phase_ml_str, phase_ml_scores,
               annotation.to_str())?;

        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use crate::analysis::rater::ConfusionMatrix;

    #[test]
    fn test_confusion_empty_matrix() {
        let empty_matrix = ConfusionMatrix::<3>::new();
        assert_eq!(empty_matrix.size(), 3);

        for i in 0..3 {
            assert_eq!(empty_matrix.get_tp(i), 0);
            assert_eq!(empty_matrix.get_fp(i), 0);
            assert_eq!(empty_matrix.get_fn(i), 0);
        }
    }

    #[test]
    fn test_confusion_perfect_matrix() {
        let mut matrix = ConfusionMatrix::<3>::new();

        for i in 0..3 {
            matrix.increment(i, i);
        }

        for i in 0..3 {
            assert_eq!(matrix.get_tp(i), 1);
            assert_eq!(matrix.get_fp(i), 0);
            assert_eq!(matrix.get_fn(i), 0);
        }
    }

    #[test]
    fn test_confusion_fp_matrix() {
        for i in 0..3 {
            let mut matrix = ConfusionMatrix::<3>::new();

            for j in 0..3 {
                matrix.increment(j, i);
            } // Same pred for all refs

            assert_eq!(matrix.get_tp(i), 1);
            assert_eq!(matrix.get_fp(i), 2);
            assert_eq!(matrix.get_fn(i), 0);
        }
    }

    #[test]
    fn test_confusion_fn_matrix() {
        for i in 0..3 {
            let mut matrix = ConfusionMatrix::<3>::new();

            for j in 0..3 {
                matrix.increment(i, j);
            } // Same pred for all refs

            assert_eq!(matrix.get_tp(i), 1);
            assert_eq!(matrix.get_fp(i), 0);
            assert_eq!(matrix.get_fn(i), 2);
        }
    }
}
