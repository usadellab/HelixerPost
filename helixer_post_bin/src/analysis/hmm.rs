use crate::results::conv::{Bases, ClassPrediction, PhasePrediction};
use crate::{log_info, vprintln};
use serde::{Deserialize, Serialize};
use std::cmp::{max, Ordering};
use std::collections::BinaryHeap;


/// User-tunable HMM parameters. Built with `HmmConfig::default()` and
/// optionally overridden by a YAML config file and/or CLI flags.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct HmmConfig {
    /// Lower bound applied to raw class/phase/base probabilities before taking
    /// the negative log, to keep penalties finite when the model emits ~0.
    pub prob_floor: f64,
    /// Weight given to the predicted phase signal when blended with the
    /// uniform coding-phase dilution target. 1.0 = trust phase fully, 0.0 =
    /// ignore phase predictions entirely.
    pub phase_retain: f64,
    /// Per-transition booleans gating which splice-junction types the decoder
    /// will consider in each surrounding context.
    pub splice: SpliceFlags,
    /// Minimum intron length in bases, for each supported intron class.
    /// Enforced in the decoder: the donor state of each class spans
    /// `length - 1` bases, so shorter introns cannot be decoded (floor: 2).
    pub minimum_intron_lengths: MinimalIntronLengths,
    /// Multiplicative weights for the start / stop / donor / acceptor signals
    /// when they enter the per-base penalty sum.
    pub weights: HmmWeights,
    /// Constant penalty added when accepting each donor variant, independent
    /// of the local base signal. Use to bias against rare splice types.
    pub donor_fixed_penalty: DonorFixedPenalties,
}

impl Default for HmmConfig {
    fn default() -> Self {
        Self {
            prob_floor: 0.000_000_001,
            phase_retain: 0.20,
            splice: SpliceFlags::default(),
            minimum_intron_lengths: MinimalIntronLengths::default(),
            weights: HmmWeights::default(),
            donor_fixed_penalty: DonorFixedPenalties::default(),
        }
    }
}

impl HmmConfig {
    pub fn phase_dilute(&self) -> f64 {
        1.0 - self.phase_retain
    }

    pub fn with_overrides(
        mut self,
        prob_floor: Option<f64>,
        phase_retain: Option<f64>,
        start_weight: Option<f64>,
        stop_weight: Option<f64>,
        donor_weight: Option<f64>,
        acceptor_weight: Option<f64>,
        min_intron_length: Option<usize>,
        min_intron_length_u12: Option<usize>,
    ) -> Self {
        if let Some(v) = prob_floor {
            self.prob_floor = v;
        }
        if let Some(v) = phase_retain {
            self.phase_retain = v;
        }
        if let Some(v) = start_weight {
            self.weights.start = v;
        }
        if let Some(v) = stop_weight {
            self.weights.stop = v;
        }
        if let Some(v) = donor_weight {
            self.weights.donor = v;
        }
        if let Some(v) = acceptor_weight {
            self.weights.acceptor = v;
        }
        if let Some(v) = min_intron_length {
            self.minimum_intron_lengths.u2_gt_ag = v;
            self.minimum_intron_lengths.u2_gc_ag = v;
        }
        if let Some(v) = min_intron_length_u12 {
            self.minimum_intron_lengths.u12_at_ac = v;
        }
        self
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SpliceFlags {
    pub utr5: bool,
    pub utr5_start: bool,
    pub start: bool,
    pub start_coding: bool,
    pub coding: bool,
    pub coding_stop: bool,
    pub stop: bool,
    pub stop_utr3: bool,
    pub utr3: bool,
}

impl Default for SpliceFlags {
    fn default() -> Self {
        Self {
            utr5: true,
            utr5_start: true,
            start: true,
            start_coding: true,
            coding: true,
            coding_stop: true,
            stop: true,
            stop_utr3: true,
            utr3: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct MinimalIntronLengths {
    pub u2_gt_ag: usize,
    pub u2_gc_ag: usize,
    //    pub u12_gt_ag: usize, // Not currently considered a unique intron class
    pub u12_at_ac: usize,
}

impl Default for MinimalIntronLengths {
    fn default() -> Self {
        Self {
            u2_gt_ag: 50,
            u2_gc_ag: 50,
            //            u12_gt_ag: 30 // Not currently considered a unique intron class
            u12_at_ac: 30,
        }
    }
}


#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct HmmWeights {
    pub start: f64,
    pub stop: f64,
    pub donor: f64,
    pub acceptor: f64,
}

impl Default for HmmWeights {
    fn default() -> Self {
        Self {
            start: 1_000.0,
            stop: 1_000.0,
            donor: 1.0,
            acceptor: 1.0,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct DonorFixedPenalties {
    pub u2_gt_ag: f64,
    pub u2_gc_ag: f64,
//    pub u12_gt_ag: f64, // Not currently considered a unique intron class
    pub u12_at_ac: f64,
}



pub fn show_hmm_config(cfg: &HmmConfig) {
    let s = &cfg.splice;
    let l = &cfg.minimum_intron_lengths;
    let w = &cfg.weights;
    let p = &cfg.donor_fixed_penalty;

    log_info!("HMM Config");
    println!(
        "  Splicing Flags: U:{} US:{} S:{} SC:{} C:{} CS:{} S:{} SU:{} U:{}",
        s.utr5, s.utr5_start, s.start, s.start_coding, s.coding,
        s.coding_stop, s.stop, s.stop_utr3, s.utr3,
    );
    println!(
        "  Splicing - Minimum Intron Lengths: U2-GT-AG {}, U2-GC-AG {} U12-AT-AC {}",
        l.u2_gt_ag, l.u2_gc_ag, /*l.u12_gt_ag,*/ l.u12_at_ac,  // Removed unused U12-GT-AG
    );
    println!(
        "  Splicing - Weights: Donor {}, Acceptor {}",
        w.donor, w.acceptor,
    );
    println!(
        "  Splicing - Fixed Penalties: U2-GT-AG {}, U2-GC-AG {} U12-AT-AC {}",
        p.u2_gt_ag, p.u2_gc_ag, /*p.u12_gt_ag,*/ p.u12_at_ac, // Removed unused U12-GT-AG
    );

    println!("  Coding - Weights: Start {}, Stop {}", w.start, w.stop);
    println!(
        "  Phase Mode: Implementation 1, Dilute to Total, Retention: {}",
        cfg.phase_retain,
    );
    println!("  Prob Floor: {}", cfg.prob_floor);
    println!();
}


fn convert_raw_pred<const N: usize>(raw_pred: &[f32; N]) -> [f64; N] {
    let mut pred: [f64; N] = [0.0; N];

    for i in 0..N {
        pred[i] = raw_pred[i] as f64
    }

    pred
}

fn raw_pred_to_neg_log_prob<const N: usize>(raw_pred: &[f64; N], floor: f64) -> [f64; N] {

    let mut neg_log_prob = [0.0; N];
    for i in 0..N {
        let adjusted_pred = raw_pred[i] as f64;

        let adjusted_pred = if adjusted_pred > floor {
            adjusted_pred
        } else {
            floor
        };
        neg_log_prob[i] = -f64::log2(adjusted_pred);
    }

    neg_log_prob
}

fn neg_log_prob_to_penalty<const N: usize>(neg_log_prob: &[f64; N]) -> [f64; N] {
    let mut min_penalty = neg_log_prob[0];
    for i in 1..N {
        if neg_log_prob[i] < min_penalty {
            min_penalty = neg_log_prob[i]
        }
    }

    let mut penalty = [0.0; N];

    for i in 0..N {
        penalty[i] = neg_log_prob[i] - min_penalty;
    }

    penalty
}

#[derive(Clone, Copy)]
struct ClassPredPenalty { // Ordering is intergenic, utr, coding, intron
    //neg_log_prob: [f64; 4], // Negated log probability, lower value is more likely
    penalty: [f64; 4], // Penalty is adjusted negated log probability with min prob (most likely) subtracted from all
}

#[allow(dead_code)]
impl ClassPredPenalty {
    pub fn get_intergenic_penalty(&self) -> f64 {
        self.penalty[0]
    }

    pub fn get_utr_penalty(&self) -> f64 {
        self.penalty[1]
    }

    pub fn get_coding_penalty(&self) -> f64 {
        self.penalty[2]
    }

    pub fn get_intron_penalty(&self) -> f64 {
        self.penalty[3]
    }
}

impl ClassPredPenalty {
    fn new(pred: &ClassPrediction, prob_floor: f64) -> Self {
        let raw_pred = pred.get();

        let converted_pred = convert_raw_pred(raw_pred);
        let neg_log_prob = raw_pred_to_neg_log_prob(&converted_pred, prob_floor);
        let penalty = neg_log_prob_to_penalty(&neg_log_prob);

        ClassPredPenalty { /*neg_log_prob,*/ penalty }
    }
}

#[derive(Clone, Copy)]
struct PhasePredPenalty { // Ordering is intergenic, utr, coding, intron
    //neg_log_prob: [f64; 4], // Negated log probability, lower value is more likely
    penalty: [f64; 4], // Penalty is adjusted negated log probability with min prob (most likely) subtracted from all
}

#[allow(dead_code)]
impl PhasePredPenalty {
    pub fn get_non_coding_penalty(&self) -> f64 {
        self.penalty[0]
    }

    pub fn get_phase0_penalty(&self) -> f64 {
        self.penalty[1]
    }

    pub fn get_phase1_penalty(&self) -> f64 {
        self.penalty[2]
    }

    pub fn get_phase2_penalty(&self) -> f64 {
        self.penalty[3]
    }
}

impl PhasePredPenalty {
    fn new(pred: &PhasePrediction, prob_floor: f64) -> Self {
        let raw_pred = pred.get();

        let converted_pred = convert_raw_pred(raw_pred);
        let neg_log_prob = raw_pred_to_neg_log_prob(&converted_pred, prob_floor);
        let penalty = neg_log_prob_to_penalty(&neg_log_prob);

        PhasePredPenalty { /* neg_log_prob, */ penalty }
    }
}

#[derive(Clone, Copy)]
struct PredPenalty { // Ordering is intergenic, utr, coding_phase0, coding_phase1, coding_phase2, intron
    neg_log_prob: [f64; 6], // Negated log probability, lower value is more likely
    penalty: [f64; 6], // Penalty is adjusted negated log probability with min prob (most likely) subtracted from all
}

#[allow(dead_code)]
impl PredPenalty {
    pub fn get_intergenic_neg_log_prob(&self) -> f64 {
        self.neg_log_prob[0]
    }
    pub fn get_intergenic_penalty(&self) -> f64 {
        self.penalty[0]
    }

    pub fn get_utr_neg_log_prob(&self) -> f64 {
        self.neg_log_prob[1]
    }
    pub fn get_utr_penalty(&self) -> f64 {
        self.penalty[1]
    }

    pub fn get_coding_phase0_neg_log_prob(&self) -> f64 {
        self.neg_log_prob[2]
    }
    pub fn get_coding_phase0_penalty(&self) -> f64 {
        self.penalty[2]
    }

    pub fn get_coding_phase1_neg_log_prob(&self) -> f64 {
        self.neg_log_prob[3]
    }
    pub fn get_coding_phase1_penalty(&self) -> f64 {
        self.penalty[3]
    }

    pub fn get_coding_phase2_neg_log_prob(&self) -> f64 {
        self.neg_log_prob[4]
    }
    pub fn get_coding_phase2_penalty(&self) -> f64 {
        self.penalty[4]
    }

    pub fn get_intron_neg_log_prob(&self) -> f64 {
        self.neg_log_prob[5]
    }
    pub fn get_intron_penalty(&self) -> f64 {
        self.penalty[5]
    }

}

impl PredPenalty {
    fn new(
        class_pred: &ClassPrediction,
        phase_pred: &PhasePrediction,
        prob_floor: f64,
        phase_retain: f64,
    ) -> Self {
        let phase_dilute = 1.0 - phase_retain;

        let phase0 = phase_pred.get_phase0() as f64;
        let phase1 = phase_pred.get_phase1() as f64;
        let phase2 = phase_pred.get_phase2() as f64;

        // Approach 1: rescale total phase to match coding and blend to dilution target (mean coding or total coding)
        let coding = class_pred.get_coding() as f64;
        let total_coding_phase = phase0 + phase1 + phase2;
        let (phase0, phase1, phase2) = if total_coding_phase > 0.0 {
            let phase_scale = coding / total_coding_phase;
            (
                phase0 * phase_scale,
                phase1 * phase_scale,
                phase2 * phase_scale,
            )
        } else {
            (coding / 3.0, coding / 3.0, coding / 3.0)
        };

        let dilution_target = coding;

        let phase0 = phase0 * phase_retain + dilution_target * phase_dilute;
        let phase1 = phase1 * phase_retain + dilution_target * phase_dilute;
        let phase2 = phase2 * phase_retain + dilution_target * phase_dilute;

        let raw_probs = [
            class_pred.get_intergenic() as f64,
            class_pred.get_utr() as f64,
            phase0,
            phase1,
            phase2,
            class_pred.get_intron() as f64,
        ];

        let neg_log_prob = raw_pred_to_neg_log_prob(&raw_probs, prob_floor);
        let penalty = neg_log_prob_to_penalty(&neg_log_prob);

        PredPenalty { neg_log_prob, penalty }
    }
}

#[derive(Clone, Copy)]
pub struct BasesPenalty { // Ordering is C, A, T, G
    //neg_log_prob: [f64; 4], // Negated log probability, lower value is more likely
    penalty: [f64; 4], // Penalty is adjusted negated log probability with min prob (most likely) subtracted from all
}

fn min2(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

fn min3(a: f64, b: f64, c: f64) -> f64 {
    let ab = if a < b { a } else { b };
    if ab < c {
        ab
    } else {
        c
    }
}

impl BasesPenalty {
    pub fn get_c(&self) -> f64 {
        self.penalty[0]
    }

    pub fn get_a(&self) -> f64 {
        self.penalty[1]
    }

    pub fn get_t(&self) -> f64 {
        self.penalty[2]
    }

    pub fn get_g(&self) -> f64 {
        self.penalty[3]
    }

    pub fn as_str(&self) -> char {
        // Pick the base whose normalised penalty is essentially zero (the
        // most-likely class). Threshold is loose because the penalty array is
        // already min-subtracted, so the winner sits at 0.0 in normal data.
        const APPROX_ZERO: f64 = 1e-9;
        if self.penalty[0] < APPROX_ZERO {
            'C'
        } else if self.penalty[1] < APPROX_ZERO {
            'A'
        } else if self.penalty[2] < APPROX_ZERO {
            'T'
        } else {
            'G'
        }
    }
}

impl BasesPenalty {
    fn new(bases: &Bases, prob_floor: f64) -> Self {
        let raw_bases = bases.get();

        let converted_pred = convert_raw_pred(raw_bases);
        let neg_log_prob = raw_pred_to_neg_log_prob(&converted_pred, prob_floor);
        let penalty = neg_log_prob_to_penalty(&neg_log_prob);

        BasesPenalty { /*neg_log_prob, */ penalty }
    }
}

struct TransitionContext<'a> {
    cfg: &'a HmmConfig,
    class_pred_pen: &'a [ClassPredPenalty],
    phase_pred_pen: &'a [PhasePredPenalty],
    pred_pen: &'a [PredPenalty],

    base_pen: &'a [BasesPenalty],
    offset: usize,
}

#[allow(dead_code)]
impl<'a> TransitionContext<'a> {
    fn new(
        cfg: &'a HmmConfig,
        class_pred_pen: &'a [ClassPredPenalty],
        phase_pred_pen: &'a [PhasePredPenalty],
        pred_pen: &'a [PredPenalty],
        base_pen: &'a [BasesPenalty],
        offset: usize,
    ) -> TransitionContext<'a> {
        TransitionContext {
            cfg,
            class_pred_pen,
            phase_pred_pen,
            pred_pen,
            base_pen,
            offset,
        }
    }

    fn get_class_pred(&self, position: usize) -> Option<&ClassPredPenalty> {
        self.class_pred_pen.get(self.offset + position)
    }

    fn get_phase_pred(&self, position: usize) -> Option<&PhasePredPenalty> {
        self.phase_pred_pen.get(self.offset + position)
    }

    fn get_pred(&self, position: usize) -> Option<&PredPenalty> {
        self.pred_pen.get(self.offset + position)
    }

    /*
        fn get_upstream(&self, len: usize) -> Option<&[BasesPenalty]>
        {
            if len <= self.offset
                { Some(&self.base_pen[self.offset-len .. self.offset]) }
            else { None }
        }
    */

    fn get_downstream(&self, len: usize) -> Option<&[BasesPenalty]> {
        if self.offset + len <= self.base_pen.len() {
            Some(&self.base_pen[self.offset..self.offset + len])
        } else {
            None
        }
    }

    fn get_ctx(&self, ulen: usize, dlen: usize) -> Option<&[BasesPenalty]> {
        if ulen <= self.offset && self.offset + dlen <= self.base_pen.len() {
            Some(&self.base_pen[self.offset - ulen..self.offset + dlen])
        } else {
            None
        }
    }

    fn get_donor_penalty_u2_gt_ag(&self, can_splice: bool) -> Option<f64> {
        if let (Some(ds), true) = (self.get_ctx(0, 2), can_splice) {
            let pen = ds[0].get_g() + ds[1].get_t();
            Some(pen * self.cfg.weights.donor + self.cfg.donor_fixed_penalty.u2_gt_ag)
        } else {
            None
        }
    }

    fn get_acceptor_penalty_u2_gt_ag(&self) -> Option<f64> {
        if let Some(us) = self.get_ctx(2, 0) {
            let pen = us[0].get_a() + us[1].get_g();
            Some(pen * self.cfg.weights.acceptor)
        } else {
            None
        }
    }

    fn get_donor_penalty_u2_gc_ag(&self, can_splice: bool) -> Option<f64> {
        if let (Some(ds), true) = (self.get_ctx(2, 2), can_splice) {
            let pen =
                //ds[0].get_a() +
                ds[1].get_g() +
                ds[2].get_g() +
                ds[3].get_c();
            Some(pen * self.cfg.weights.donor + self.cfg.donor_fixed_penalty.u2_gc_ag)
        } else {
            None
        }
    }

    fn get_acceptor_penalty_u2_gc_ag(&self) -> Option<f64> {
        if let Some(us) = self.get_ctx(2, 0) {
            let pen = us[0].get_a() + us[1].get_g();
            Some(pen * self.cfg.weights.acceptor)
        } else {
            None
        }
    }

    fn get_donor_penalty_u12_at_ac(&self, can_splice: bool) -> Option<f64> {
        if let (Some(ds), true) = (self.get_ctx(0, 7), can_splice) {
            let pen = // ATATCCT
                ds[0].get_a() +
                ds[1].get_t() +
                ds[2].get_a() +
                ds[3].get_t() +
                ds[4].get_c() +
                ds[5].get_c() +
                ds[6].get_t();

            Some(pen * self.cfg.weights.donor + self.cfg.donor_fixed_penalty.u12_at_ac)
        } else {
            None
        }
    }

    fn get_acceptor_penalty_u12_at_ac(&self) -> Option<f64> {
        if let Some(us) = self.get_ctx(2, 0) {
            let pen = us[0].get_a() + us[1].get_c();
            Some(pen * self.cfg.weights.acceptor)
        } else {
            None
        }
    }
}

#[allow(dead_code)]
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum HmmAnnotationLabel {
    Intergenic,
    UTR5,
    Start,
    Coding,
    Intron,
    Stop,
    UTR3,
}

impl HmmAnnotationLabel {
    pub fn to_str(&self) -> &str {
        match self {
            HmmAnnotationLabel::Intergenic => "Intergenic",
            HmmAnnotationLabel::UTR5 => "UTR5",
            HmmAnnotationLabel::Start => "Start",
            HmmAnnotationLabel::Coding => "Coding",
            HmmAnnotationLabel::Intron => "Intron",
            HmmAnnotationLabel::Stop => "Stop",
            HmmAnnotationLabel::UTR3 => "UTR3",
        }
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum HmmPrimaryState {
    Intergenic,
    UTR5,
    Start0, // Possible Start - After A
    Start1, // Possible Start - After AT
    Start2, // Possible Start - After ATG
    Coding0,
    Coding1,
    Coding2,
    Stop0T,  // Possible Stop - After T
    Stop1TA, // Possible Stop - After TA
    Stop1TG, // Possible Stop - After TG
    Stop2,   // Possible Stop - After TAA / TAG / TGA
    UTR3,
}

impl HmmPrimaryState {
    pub fn to_str(&self) -> &str {
        match self {
            HmmPrimaryState::Intergenic => "Intergenic",
            HmmPrimaryState::UTR5 => "UTR5",
            HmmPrimaryState::Start0 => "Start0",
            HmmPrimaryState::Start1 => "Start1",
            HmmPrimaryState::Start2 => "Start2",
            HmmPrimaryState::Coding0 => "Coding0",
            HmmPrimaryState::Coding1 => "Coding1",
            HmmPrimaryState::Coding2 => "Coding2",
            HmmPrimaryState::Stop0T => "Stop0T",
            HmmPrimaryState::Stop1TA => "Stop1TA",
            HmmPrimaryState::Stop1TG => "Stop1TG",
            HmmPrimaryState::Stop2 => "Stop2",
            HmmPrimaryState::UTR3 => "UTR3",
        }
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum HmmIntronState {
    None = 0,
    U2GtAgDSS = 1,
    U2GtAg = 2,
    U2GcAgDSS = 3,
    U2GcAg = 4,
    U12AtAcDSS = 5,
    U12AtAc = 6,
}

//const HMM_INTRON_STATES: usize = 7;

impl HmmIntronState {
    pub fn to_str(&self) -> &str {
        match self {
            HmmIntronState::None => "None",
            HmmIntronState::U2GtAgDSS => "U2GtAgDSS",
            HmmIntronState::U2GtAg => "U2GtAg",
            HmmIntronState::U2GcAgDSS => "U2GcAgDSS",
            HmmIntronState::U2GcAg => "U2GcAg",
            HmmIntronState::U12AtAcDSS => "U12AtAcDSS",
            HmmIntronState::U12AtAc => "U12AtAc",
        }
    }
}

const HMM_STATES: usize = 73;

#[derive(Clone, Copy, Eq, PartialEq, Ord)]
enum HmmState {
    Intergenic = 0,

    UTR5 = 1,
    UTR5IntronU2GtAgDSS = 2,
    UTR5IntronU2GtAg = 3,
    UTR5IntronU2GcAgDSS = 4,
    UTR5IntronU2GcAg = 5,
    UTR5IntronU12AtAcDSS = 6,
    UTR5IntronU12AtAc = 7,

    Start0 = 8, // After A
    Start0IntronU2GtAgDSS = 9,
    Start0IntronU2GtAg = 10,
    Start0IntronU2GcAgDSS = 11,
    Start0IntronU2GcAg = 12,
    Start0IntronU12AtAcDSS = 13,
    Start0IntronU12AtAc = 14,

    Start1 = 15, // After AT
    Start1IntronU2GtAgDSS = 16,
    Start1IntronU2GtAg = 17,
    Start1IntronU2GcAgDSS = 18,
    Start1IntronU2GcAg = 19,
    Start1IntronU12AtAcDSS = 20,
    Start1IntronU12AtAc = 21,

    Start2 = 22, // After ATG

    Coding0 = 23,
    Coding0IntronU2GtAgDSS = 24,
    Coding0IntronU2GtAg = 25,
    Coding0IntronU2GcAgDSS = 26,
    Coding0IntronU2GcAg = 27,
    Coding0IntronU12AtAcDSS = 28,
    Coding0IntronU12AtAc = 29,

    Coding1 = 30,
    Coding1IntronU2GtAgDSS = 31,
    Coding1IntronU2GtAg = 32,
    Coding1IntronU2GcAgDSS = 33,
    Coding1IntronU2GcAg = 34,
    Coding1IntronU12AtAcDSS = 35,
    Coding1IntronU12AtAc = 36,

    Coding2 = 37,
    Coding2IntronU2GtAgDSS = 38,
    Coding2IntronU2GtAg = 39,
    Coding2IntronU2GcAgDSS = 40,
    Coding2IntronU2GcAg = 41,
    Coding2IntronU12AtAcDSS = 42,
    Coding2IntronU12AtAc = 43,

    Stop0T = 44,
    Stop0TIntronU2GtAgDSS = 45,
    Stop0TIntronU2GtAg = 46,
    Stop0TIntronU2GcAgDSS = 47,
    Stop0TIntronU2GcAg = 48,
    Stop0TIntronU12AtAcDSS = 49,
    Stop0TIntronU12AtAc = 50,

    Stop1TA = 51,
    Stop1TAIntronU2GtAgDSS = 52,
    Stop1TAIntronU2GtAg = 53,
    Stop1TAIntronU2GcAgDSS = 54,
    Stop1TAIntronU2GcAg = 55,
    Stop1TAIntronU12AtAcDSS = 56,
    Stop1TAIntronU12AtAc = 57,

    Stop1TG = 58,
    Stop1TGIntronU2GtAgDSS = 59,
    Stop1TGIntronU2GtAg = 60,
    Stop1TGIntronU2GcAgDSS = 61,
    Stop1TGIntronU2GcAg = 62,
    Stop1TGIntronU12AtAcDSS = 63,
    Stop1TGIntronU12AtAc = 64,

    Stop2 = 65,

    UTR3 = 66,
    UTR3IntronU2GtAgDSS = 67,
    UTR3IntronU2GtAg = 68,
    UTR3IntronU2GcAgDSS = 69,
    UTR3IntronU2GcAg = 70,
    UTR3IntronU12AtAcDSS = 71,
    UTR3IntronU12AtAc = 72,
}

impl std::cmp::PartialOrd for HmmState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let s = *self as u8;
        let o = *other as u8;

        Some(s.cmp(&o))
    }
}

impl HmmState {
    fn get_component_states(self) -> (HmmPrimaryState, HmmIntronState) {
        match self {
            HmmState::Intergenic => (HmmPrimaryState::Intergenic, HmmIntronState::None),

            HmmState::UTR5 => (HmmPrimaryState::UTR5, HmmIntronState::None),
            HmmState::UTR5IntronU2GtAgDSS => (HmmPrimaryState::UTR5, HmmIntronState::U2GtAgDSS),
            HmmState::UTR5IntronU2GtAg => (HmmPrimaryState::UTR5, HmmIntronState::U2GtAg),
            HmmState::UTR5IntronU2GcAgDSS => (HmmPrimaryState::UTR5, HmmIntronState::U2GcAgDSS),
            HmmState::UTR5IntronU2GcAg => (HmmPrimaryState::UTR5, HmmIntronState::U2GcAg),
            HmmState::UTR5IntronU12AtAcDSS => (HmmPrimaryState::UTR5, HmmIntronState::U12AtAcDSS),
            HmmState::UTR5IntronU12AtAc => (HmmPrimaryState::UTR5, HmmIntronState::U12AtAc),

            HmmState::Start0 => (HmmPrimaryState::Start0, HmmIntronState::None),
            HmmState::Start0IntronU2GtAgDSS => (HmmPrimaryState::Start0, HmmIntronState::U2GtAgDSS),
            HmmState::Start0IntronU2GtAg => (HmmPrimaryState::Start0, HmmIntronState::U2GtAg),
            HmmState::Start0IntronU2GcAgDSS => (HmmPrimaryState::Start0, HmmIntronState::U2GcAgDSS),
            HmmState::Start0IntronU2GcAg => (HmmPrimaryState::Start0, HmmIntronState::U2GcAg),
            HmmState::Start0IntronU12AtAcDSS => {
                (HmmPrimaryState::Start0, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Start0IntronU12AtAc => (HmmPrimaryState::Start0, HmmIntronState::U12AtAc),

            HmmState::Start1 => (HmmPrimaryState::Start1, HmmIntronState::None),
            HmmState::Start1IntronU2GtAgDSS => (HmmPrimaryState::Start1, HmmIntronState::U2GtAgDSS),
            HmmState::Start1IntronU2GtAg => (HmmPrimaryState::Start1, HmmIntronState::U2GtAg),
            HmmState::Start1IntronU2GcAgDSS => (HmmPrimaryState::Start1, HmmIntronState::U2GcAgDSS),
            HmmState::Start1IntronU2GcAg => (HmmPrimaryState::Start1, HmmIntronState::U2GcAg),
            HmmState::Start1IntronU12AtAcDSS => {
                (HmmPrimaryState::Start1, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Start1IntronU12AtAc => (HmmPrimaryState::Start1, HmmIntronState::U12AtAc),

            HmmState::Start2 => (HmmPrimaryState::Start2, HmmIntronState::None),

            HmmState::Coding0 => (HmmPrimaryState::Coding0, HmmIntronState::None),
            HmmState::Coding0IntronU2GtAgDSS => {
                (HmmPrimaryState::Coding0, HmmIntronState::U2GtAgDSS)
            }
            HmmState::Coding0IntronU2GtAg => (HmmPrimaryState::Coding0, HmmIntronState::U2GtAg),
            HmmState::Coding0IntronU2GcAgDSS => {
                (HmmPrimaryState::Coding0, HmmIntronState::U2GcAgDSS)
            }
            HmmState::Coding0IntronU2GcAg => (HmmPrimaryState::Coding0, HmmIntronState::U2GcAg),
            HmmState::Coding0IntronU12AtAcDSS => {
                (HmmPrimaryState::Coding0, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Coding0IntronU12AtAc => (HmmPrimaryState::Coding0, HmmIntronState::U12AtAc),

            HmmState::Coding1 => (HmmPrimaryState::Coding1, HmmIntronState::None),
            HmmState::Coding1IntronU2GtAgDSS => {
                (HmmPrimaryState::Coding1, HmmIntronState::U2GtAgDSS)
            }
            HmmState::Coding1IntronU2GtAg => (HmmPrimaryState::Coding1, HmmIntronState::U2GtAg),
            HmmState::Coding1IntronU2GcAgDSS => {
                (HmmPrimaryState::Coding1, HmmIntronState::U2GcAgDSS)
            }
            HmmState::Coding1IntronU2GcAg => (HmmPrimaryState::Coding1, HmmIntronState::U2GcAg),
            HmmState::Coding1IntronU12AtAcDSS => {
                (HmmPrimaryState::Coding1, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Coding1IntronU12AtAc => (HmmPrimaryState::Coding1, HmmIntronState::U12AtAc),

            HmmState::Coding2 => (HmmPrimaryState::Coding2, HmmIntronState::None),
            HmmState::Coding2IntronU2GtAgDSS => {
                (HmmPrimaryState::Coding2, HmmIntronState::U2GtAgDSS)
            }
            HmmState::Coding2IntronU2GtAg => (HmmPrimaryState::Coding2, HmmIntronState::U2GtAg),
            HmmState::Coding2IntronU2GcAgDSS => {
                (HmmPrimaryState::Coding2, HmmIntronState::U2GcAgDSS)
            }
            HmmState::Coding2IntronU2GcAg => (HmmPrimaryState::Coding2, HmmIntronState::U2GcAg),
            HmmState::Coding2IntronU12AtAcDSS => {
                (HmmPrimaryState::Coding2, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Coding2IntronU12AtAc => (HmmPrimaryState::Coding2, HmmIntronState::U12AtAc),

            HmmState::Stop0T => (HmmPrimaryState::Stop0T, HmmIntronState::None),
            HmmState::Stop0TIntronU2GtAgDSS => (HmmPrimaryState::Stop0T, HmmIntronState::U2GtAgDSS),
            HmmState::Stop0TIntronU2GtAg => (HmmPrimaryState::Stop0T, HmmIntronState::U2GtAg),
            HmmState::Stop0TIntronU2GcAgDSS => (HmmPrimaryState::Stop0T, HmmIntronState::U2GcAgDSS),
            HmmState::Stop0TIntronU2GcAg => (HmmPrimaryState::Stop0T, HmmIntronState::U2GcAg),
            HmmState::Stop0TIntronU12AtAcDSS => {
                (HmmPrimaryState::Stop0T, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Stop0TIntronU12AtAc => (HmmPrimaryState::Stop0T, HmmIntronState::U12AtAc),

            HmmState::Stop1TA => (HmmPrimaryState::Stop1TA, HmmIntronState::None),
            HmmState::Stop1TAIntronU2GtAgDSS => {
                (HmmPrimaryState::Stop1TA, HmmIntronState::U2GtAgDSS)
            }
            HmmState::Stop1TAIntronU2GtAg => (HmmPrimaryState::Stop1TA, HmmIntronState::U2GtAg),
            HmmState::Stop1TAIntronU2GcAgDSS => {
                (HmmPrimaryState::Stop1TA, HmmIntronState::U2GcAgDSS)
            }
            HmmState::Stop1TAIntronU2GcAg => (HmmPrimaryState::Stop1TA, HmmIntronState::U2GcAg),
            HmmState::Stop1TAIntronU12AtAcDSS => {
                (HmmPrimaryState::Stop1TA, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Stop1TAIntronU12AtAc => (HmmPrimaryState::Stop1TA, HmmIntronState::U12AtAc),

            HmmState::Stop1TG => (HmmPrimaryState::Stop1TG, HmmIntronState::None),
            HmmState::Stop1TGIntronU2GtAgDSS => {
                (HmmPrimaryState::Stop1TG, HmmIntronState::U2GtAgDSS)
            }
            HmmState::Stop1TGIntronU2GtAg => (HmmPrimaryState::Stop1TG, HmmIntronState::U2GtAg),
            HmmState::Stop1TGIntronU2GcAgDSS => {
                (HmmPrimaryState::Stop1TG, HmmIntronState::U2GcAgDSS)
            }
            HmmState::Stop1TGIntronU2GcAg => (HmmPrimaryState::Stop1TG, HmmIntronState::U2GcAg),
            HmmState::Stop1TGIntronU12AtAcDSS => {
                (HmmPrimaryState::Stop1TG, HmmIntronState::U12AtAcDSS)
            }
            HmmState::Stop1TGIntronU12AtAc => (HmmPrimaryState::Stop1TG, HmmIntronState::U12AtAc),

            HmmState::Stop2 => (HmmPrimaryState::Stop2, HmmIntronState::None),

            HmmState::UTR3 => (HmmPrimaryState::UTR3, HmmIntronState::None),
            HmmState::UTR3IntronU2GtAgDSS => (HmmPrimaryState::UTR3, HmmIntronState::U2GtAgDSS),
            HmmState::UTR3IntronU2GtAg => (HmmPrimaryState::UTR3, HmmIntronState::U2GtAg),
            HmmState::UTR3IntronU2GcAgDSS => (HmmPrimaryState::UTR3, HmmIntronState::U2GcAgDSS),
            HmmState::UTR3IntronU2GcAg => (HmmPrimaryState::UTR3, HmmIntronState::U2GcAg),
            HmmState::UTR3IntronU12AtAcDSS => (HmmPrimaryState::UTR3, HmmIntronState::U12AtAcDSS),
            HmmState::UTR3IntronU12AtAc => (HmmPrimaryState::UTR3, HmmIntronState::U12AtAc),
        }
    }

    fn get_annotation_label(self) -> HmmAnnotationLabel {
        let (primary, intron) = self.get_component_states();

        if intron != HmmIntronState::None {
            return HmmAnnotationLabel::Intron;
        }

        match primary {
            HmmPrimaryState::Intergenic => HmmAnnotationLabel::Intergenic,

            HmmPrimaryState::UTR5 => HmmAnnotationLabel::UTR5,

            HmmPrimaryState::Start0 | HmmPrimaryState::Start1 | HmmPrimaryState::Start2 => {
                HmmAnnotationLabel::Coding
            } //Start,

            HmmPrimaryState::Coding0 | HmmPrimaryState::Coding1 | HmmPrimaryState::Coding2 => {
                HmmAnnotationLabel::Coding
            }

            HmmPrimaryState::Stop0T
            | HmmPrimaryState::Stop1TA
            | HmmPrimaryState::Stop1TG
            | HmmPrimaryState::Stop2 => HmmAnnotationLabel::Coding, //Stop,

            HmmPrimaryState::UTR3 => HmmAnnotationLabel::UTR3,
        }
    }

    #[allow(unused_variables)]
    fn get_state_penalty(
        self,
        class_pred: &ClassPredPenalty,
        phase_pred: &PhasePredPenalty,
        pred: &PredPenalty,
    ) -> (f64, f64) {
        let (primary, intron) = self.get_component_states();

        if intron != HmmIntronState::None {
            return (pred.get_intron_neg_log_prob(), pred.get_intron_penalty());
        }

        match primary {
            HmmPrimaryState::Intergenic => (pred.get_intergenic_neg_log_prob(), pred.get_intergenic_penalty()),

            HmmPrimaryState::UTR5 | HmmPrimaryState::UTR3 =>
                (pred.get_utr_neg_log_prob(), pred.get_utr_penalty()),

            HmmPrimaryState::Coding0 => (pred.get_coding_phase0_neg_log_prob(), pred.get_coding_phase0_penalty()),
            HmmPrimaryState::Start0 | HmmPrimaryState::Stop0T =>
                (pred.get_coding_phase0_neg_log_prob(), pred.get_coding_phase0_penalty()),

            //HmmPrimaryState::Coding1 => (pred.get_coding_phase2_penalty(), pred.get_coding_phase2_neg_log_prob()),
            HmmPrimaryState::Coding1 => (pred.get_coding_phase2_neg_log_prob(), pred.get_coding_phase2_penalty()),
            HmmPrimaryState::Start1 | HmmPrimaryState::Stop1TA | HmmPrimaryState::Stop1TG =>
                (pred.get_coding_phase2_neg_log_prob(), pred.get_coding_phase2_penalty()),

            HmmPrimaryState::Coding2 => (pred.get_coding_phase1_neg_log_prob(), pred.get_coding_phase1_penalty()),
            HmmPrimaryState::Start2 | HmmPrimaryState::Stop2 =>
                (pred.get_coding_phase1_neg_log_prob(), pred.get_coding_phase1_penalty()),
        }
    }

    fn get_base_count(self, mil: &MinimalIntronLengths) -> usize {
        let (_, intron) = self.get_component_states();

        match intron {
            HmmIntronState::U2GtAgDSS => max(mil.u2_gt_ag, 2) - 1,
            HmmIntronState::U2GcAgDSS => max(mil.u2_gc_ag, 2) - 1,
            HmmIntronState::U12AtAcDSS => max(mil.u12_at_ac, 2) - 1,
            _ => 1,
        }
    }

    // Calculate the common (minimum) penalty for each 'destination' state - based on either DSS (intron start) or primary states with base matches (start/stop)
    // Valid for Intron DSS and all non-intron states
    fn get_common_state_entrance_penalty(
        self: HmmState,
        trans_ctx: &TransitionContext,
    ) -> Option<f64> {
        let (primary, intron) = self.get_component_states();

        if intron != HmmIntronState::None {
            match intron {
                HmmIntronState::U2GtAgDSS => return trans_ctx.get_donor_penalty_u2_gt_ag(true),
                HmmIntronState::U2GcAgDSS => return trans_ctx.get_donor_penalty_u2_gc_ag(true),
                HmmIntronState::U12AtAcDSS => return trans_ctx.get_donor_penalty_u12_at_ac(true),

                _ => panic!(
                    "Called get_state_entrance_penalty with unexpected state {} - {} {}",
                    self as u8,
                    primary.to_str(),
                    intron.to_str()
                ),
            }
        }

        match primary {
            HmmPrimaryState::Intergenic => Some(0.0),

            HmmPrimaryState::UTR5 => Some(0.0),

            HmmPrimaryState::Start0 => trans_ctx
                .get_downstream(1)
                .map(|ds| ds[0].get_a() * trans_ctx.cfg.weights.start),
            HmmPrimaryState::Start1 => trans_ctx
                .get_downstream(1)
                .map(|ds| ds[0].get_t() * trans_ctx.cfg.weights.start),
            HmmPrimaryState::Start2 => trans_ctx
                .get_downstream(1)
                .map(|ds| ds[0].get_g() * trans_ctx.cfg.weights.start),

            HmmPrimaryState::Coding0 => trans_ctx
                .get_downstream(1)
                .map(|ds| {
                    min3(ds[0].get_a(), ds[0].get_c(), ds[0].get_g())
                        * trans_ctx.cfg.weights.stop
                }),
            HmmPrimaryState::Coding1 => Some(0.0),
            HmmPrimaryState::Coding2 => Some(0.0),

            HmmPrimaryState::Stop0T => trans_ctx
                .get_downstream(1)
                .map(|ds| ds[0].get_t() * trans_ctx.cfg.weights.stop),
            HmmPrimaryState::Stop1TA => trans_ctx
                .get_downstream(1)
                .map(|ds| ds[0].get_a() * trans_ctx.cfg.weights.stop),
            HmmPrimaryState::Stop1TG => trans_ctx
                .get_downstream(1)
                .map(|ds| ds[0].get_g() * trans_ctx.cfg.weights.stop),
            HmmPrimaryState::Stop2 => Some(0.0),

            HmmPrimaryState::UTR3 => Some(0.0),
        }
    }

    fn populate_successor_states_and_transition_penalties(
        self,
        trans_ctx: &TransitionContext,
        successors: &mut Vec<(HmmState, f64)>,
    ) {
        let consider_transition = |successors: &mut Vec<(HmmState, f64)>,
                                   new_state: HmmState,
                                   other_pen: Option<f64>,
                                   allow: bool| {
            if let (true, Some(ex), Some(en)) = (
                allow,
                other_pen,
                new_state.get_common_state_entrance_penalty(trans_ctx),
            ) {
                successors.push((new_state, ex + en));
            }
        };

        let (_, intron) = self.get_component_states();

        let acceptor_penalty = match intron {
            HmmIntronState::None => Some(0.0),
            HmmIntronState::U2GtAg => trans_ctx.get_acceptor_penalty_u2_gt_ag(),
            HmmIntronState::U2GcAg => trans_ctx.get_acceptor_penalty_u2_gc_ag(),
            HmmIntronState::U12AtAc => trans_ctx.get_acceptor_penalty_u12_at_ac(),
            _ => None,
        };

        match self {
            HmmState::Intergenic => {
                successors.push((self, 0.0));
                consider_transition(successors, HmmState::UTR5, Some(0.0), true);
            }

            HmmState::UTR5 => {
                successors.push((self, 0.0));
                consider_transition(successors, HmmState::Start0, Some(0.0), true);
                consider_transition(
                    successors,
                    HmmState::UTR5IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.utr5,
                );
                consider_transition(
                    successors,
                    HmmState::UTR5IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.utr5,
                );
                consider_transition(
                    successors,
                    HmmState::UTR5IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.utr5,
                );
            }

            HmmState::UTR5IntronU2GtAgDSS => successors.push((HmmState::UTR5IntronU2GtAg, 0.0)),
            HmmState::UTR5IntronU2GcAgDSS => successors.push((HmmState::UTR5IntronU2GcAg, 0.0)),
            HmmState::UTR5IntronU12AtAcDSS => successors.push((HmmState::UTR5IntronU12AtAc, 0.0)),

            HmmState::UTR5IntronU2GtAg
            | HmmState::UTR5IntronU2GcAg
            | HmmState::UTR5IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::UTR5,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.utr5,
                );
                consider_transition(
                    successors,
                    HmmState::Start0,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.utr5_start,
                );
            }

            HmmState::Start0 => {
                consider_transition(successors, HmmState::Start1, Some(0.0), true);
                consider_transition(
                    successors,
                    HmmState::Start0IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.start,
                );
                consider_transition(
                    successors,
                    HmmState::Start0IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.start,
                );
                consider_transition(
                    successors,
                    HmmState::Start0IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.start,
                );
            }

            HmmState::Start0IntronU2GtAgDSS => successors.push((HmmState::Start0IntronU2GtAg, 0.0)),
            HmmState::Start0IntronU2GcAgDSS => successors.push((HmmState::Start0IntronU2GcAg, 0.0)),
            HmmState::Start0IntronU12AtAcDSS => {
                successors.push((HmmState::Start0IntronU12AtAc, 0.0))
            }

            HmmState::Start0IntronU2GtAg
            | HmmState::Start0IntronU2GcAg
            | HmmState::Start0IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::Start1,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.start,
                );
            }

            HmmState::Start1 => {
                consider_transition(successors, HmmState::Start2, Some(0.0), true);
                consider_transition(
                    successors,
                    HmmState::Start1IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.start,
                );
                consider_transition(
                    successors,
                    HmmState::Start1IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.start,
                );
                consider_transition(
                    successors,
                    HmmState::Start1IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.start,
                );
            }

            HmmState::Start1IntronU2GtAgDSS => successors.push((HmmState::Start1IntronU2GtAg, 0.0)),
            HmmState::Start1IntronU2GcAgDSS => successors.push((HmmState::Start1IntronU2GcAg, 0.0)),
            HmmState::Start1IntronU12AtAcDSS => {
                successors.push((HmmState::Start1IntronU12AtAc, 0.0))
            }

            HmmState::Start1IntronU2GtAg
            | HmmState::Start1IntronU2GcAg
            | HmmState::Start1IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::Start2,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.start,
                );
            }

            HmmState::Coding0 => {
                consider_transition(successors, HmmState::Coding1, Some(0.0), true);

                consider_transition(
                    successors,
                    HmmState::Coding0IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Coding0IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Coding0IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
            }

            HmmState::Coding0IntronU2GtAgDSS => {
                successors.push((HmmState::Coding0IntronU2GtAg, 0.0))
            }
            HmmState::Coding0IntronU2GcAgDSS => {
                successors.push((HmmState::Coding0IntronU2GcAg, 0.0))
            }
            HmmState::Coding0IntronU12AtAcDSS => {
                successors.push((HmmState::Coding0IntronU12AtAc, 0.0))
            }

            HmmState::Coding0IntronU2GtAg
            | HmmState::Coding0IntronU2GcAg
            | HmmState::Coding0IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::Coding1,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.coding,
                );
            }

            HmmState::Coding1 => {
                consider_transition(successors, HmmState::Coding2, Some(0.0), true);

                consider_transition(
                    successors,
                    HmmState::Coding1IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Coding1IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Coding1IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
            }

            HmmState::Coding1IntronU2GtAgDSS => {
                successors.push((HmmState::Coding1IntronU2GtAg, 0.0))
            }
            HmmState::Coding1IntronU2GcAgDSS => {
                successors.push((HmmState::Coding1IntronU2GcAg, 0.0))
            }
            HmmState::Coding1IntronU12AtAcDSS => {
                successors.push((HmmState::Coding1IntronU12AtAc, 0.0))
            }

            HmmState::Coding1IntronU2GtAg
            | HmmState::Coding1IntronU2GcAg
            | HmmState::Coding1IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::Coding2,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.coding,
                );
            }

            HmmState::Start2 | HmmState::Coding2 => {
                consider_transition(successors, HmmState::Coding0, Some(0.0), true);
                consider_transition(successors, HmmState::Stop0T, Some(0.0), true);

                consider_transition(
                    successors,
                    HmmState::Coding2IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Coding2IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Coding2IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.coding,
                );
            }

            HmmState::Coding2IntronU2GtAgDSS => {
                successors.push((HmmState::Coding2IntronU2GtAg, 0.0))
            }
            HmmState::Coding2IntronU2GcAgDSS => {
                successors.push((HmmState::Coding2IntronU2GcAg, 0.0))
            }
            HmmState::Coding2IntronU12AtAcDSS => {
                successors.push((HmmState::Coding2IntronU12AtAc, 0.0))
            }

            HmmState::Coding2IntronU2GtAg
            | HmmState::Coding2IntronU2GcAg
            | HmmState::Coding2IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::Coding0,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.coding,
                );
                consider_transition(
                    successors,
                    HmmState::Stop0T,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.coding_stop,
                );
            }

            HmmState::Stop0T =>
            // Equivalent to Coding0, but potentially a stop codon (Txx)
            {
                consider_transition(successors, HmmState::Stop1TA, Some(0.0), true);
                consider_transition(successors, HmmState::Stop1TG, Some(0.0), true);

                if let Some(ds) = trans_ctx.get_downstream(1) {
                    consider_transition(
                        successors,
                        HmmState::Coding1,
                        Some(min2(ds[0].get_c(), ds[0].get_t()) * trans_ctx.cfg.weights.stop),
                        true,
                    );
                }

                consider_transition(
                    successors,
                    HmmState::Stop0TIntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
                consider_transition(
                    successors,
                    HmmState::Stop0TIntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
                consider_transition(
                    successors,
                    HmmState::Stop0TIntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
            }

            HmmState::Stop0TIntronU2GtAgDSS => successors.push((HmmState::Stop0TIntronU2GtAg, 0.0)),
            HmmState::Stop0TIntronU2GcAgDSS => successors.push((HmmState::Stop0TIntronU2GcAg, 0.0)),
            HmmState::Stop0TIntronU12AtAcDSS => {
                successors.push((HmmState::Stop0TIntronU12AtAc, 0.0))
            }

            HmmState::Stop0TIntronU2GtAg
            | HmmState::Stop0TIntronU2GcAg
            | HmmState::Stop0TIntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(successors, HmmState::Stop1TA, acceptor_penalty, true);
                consider_transition(successors, HmmState::Stop1TG, acceptor_penalty, true);

                if let (Some(acceptor_penalty), Some(ds)) =
                    (acceptor_penalty, trans_ctx.get_downstream(1))
                {
                    consider_transition(
                        successors,
                        HmmState::Coding1,
                        Some(acceptor_penalty + min2(ds[0].get_c(), ds[0].get_t()) * trans_ctx.cfg.weights.stop),
                        trans_ctx.cfg.splice.stop,
                    );
                }
            }

            HmmState::Stop1TA =>
            // Equivalent to Coding1, but potentially a stop codon (TAx)
            {
                if let Some(ds) = trans_ctx.get_downstream(1) {
                    consider_transition(
                        successors,
                        HmmState::Stop2,
                        Some(min2(ds[0].get_a(), ds[0].get_g()) * trans_ctx.cfg.weights.stop),
                        true,
                    );
                    consider_transition(
                        successors,
                        HmmState::Coding2,
                        Some(min2(ds[0].get_c(), ds[0].get_t()) * trans_ctx.cfg.weights.stop),
                        true,
                    );
                }

                consider_transition(
                    successors,
                    HmmState::Stop1TAIntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
                consider_transition(
                    successors,
                    HmmState::Stop1TAIntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
                consider_transition(
                    successors,
                    HmmState::Stop1TAIntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
            }

            HmmState::Stop1TAIntronU2GtAgDSS => {
                successors.push((HmmState::Stop1TAIntronU2GtAg, 0.0))
            }
            HmmState::Stop1TAIntronU2GcAgDSS => {
                successors.push((HmmState::Stop1TAIntronU2GcAg, 0.0))
            }
            HmmState::Stop1TAIntronU12AtAcDSS => {
                successors.push((HmmState::Stop1TAIntronU12AtAc, 0.0))
            }

            HmmState::Stop1TAIntronU2GtAg
            | HmmState::Stop1TAIntronU2GcAg
            | HmmState::Stop1TAIntronU12AtAc => {
                successors.push((self, 0.0));

                if let (Some(acceptor_penalty), Some(ds)) =
                    (acceptor_penalty, trans_ctx.get_downstream(1))
                {
                    consider_transition(
                        successors,
                        HmmState::Stop2,
                        Some(acceptor_penalty + min2(ds[0].get_a(), ds[0].get_g()) * trans_ctx.cfg.weights.stop),
                        true,
                    );
                    consider_transition(
                        successors,
                        HmmState::Coding2,
                        Some(acceptor_penalty + min2(ds[0].get_c(), ds[0].get_t()) * trans_ctx.cfg.weights.stop),
                        true,
                    );
                }
            }

            HmmState::Stop1TG =>
            // Equivalent to Coding1, but potentially a stop codon (TGx)
            {
                if let Some(ds) = trans_ctx.get_downstream(1) {
                    consider_transition(
                        successors,
                        HmmState::Stop2,
                        Some(ds[0].get_a() * trans_ctx.cfg.weights.stop),
                        true,
                    );
                    consider_transition(
                        successors,
                        HmmState::Coding2,
                        Some(min3(ds[0].get_c(), ds[0].get_g(), ds[0].get_t()) * trans_ctx.cfg.weights.stop),
                        true,
                    );
                }

                consider_transition(
                    successors,
                    HmmState::Stop1TGIntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
                consider_transition(
                    successors,
                    HmmState::Stop1TGIntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
                consider_transition(
                    successors,
                    HmmState::Stop1TGIntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop,
                );
            }

            HmmState::Stop1TGIntronU2GtAgDSS => {
                successors.push((HmmState::Stop1TGIntronU2GtAg, 0.0))
            }
            HmmState::Stop1TGIntronU2GcAgDSS => {
                successors.push((HmmState::Stop1TGIntronU2GcAg, 0.0))
            }
            HmmState::Stop1TGIntronU12AtAcDSS => {
                successors.push((HmmState::Stop1TGIntronU12AtAc, 0.0))
            }

            HmmState::Stop1TGIntronU2GtAg
            | HmmState::Stop1TGIntronU2GcAg
            | HmmState::Stop1TGIntronU12AtAc => {
                successors.push((self, 0.0));

                if let (Some(acceptor_penalty), Some(ds)) =
                    (acceptor_penalty, trans_ctx.get_downstream(1))
                {
                    consider_transition(
                        successors,
                        HmmState::Stop2,
                        Some(acceptor_penalty + ds[0].get_a() * trans_ctx.cfg.weights.stop),
                        true,
                    );
                    consider_transition(
                        successors,
                        HmmState::Coding2,
                        Some(
                            acceptor_penalty
                                + min3(ds[0].get_c(), ds[0].get_g(), ds[0].get_t()) * trans_ctx.cfg.weights.stop,
                        ),
                        true,
                    );
                }
            }

            HmmState::Stop2 => {
                successors.push((HmmState::UTR3, 0.0));
                consider_transition(
                    successors,
                    HmmState::UTR3IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop_utr3,
                );
                consider_transition(
                    successors,
                    HmmState::UTR3IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop_utr3,
                );
                consider_transition(
                    successors,
                    HmmState::UTR3IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.stop_utr3,
                );
            }

            HmmState::UTR3 => {
                successors.push((self, 0.0));
                successors.push((HmmState::Intergenic, 0.0));
                consider_transition(
                    successors,
                    HmmState::UTR3IntronU2GtAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.utr3,
                );
                consider_transition(
                    successors,
                    HmmState::UTR3IntronU2GcAgDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.utr3,
                );
                consider_transition(
                    successors,
                    HmmState::UTR3IntronU12AtAcDSS,
                    Some(0.0),
                    trans_ctx.cfg.splice.utr3,
                );
            }

            HmmState::UTR3IntronU2GtAgDSS => successors.push((HmmState::UTR3IntronU2GtAg, 0.0)),
            HmmState::UTR3IntronU2GcAgDSS => successors.push((HmmState::UTR3IntronU2GcAg, 0.0)),
            HmmState::UTR3IntronU12AtAcDSS => successors.push((HmmState::UTR3IntronU12AtAc, 0.0)),

            HmmState::UTR3IntronU2GtAg
            | HmmState::UTR3IntronU2GcAg
            | HmmState::UTR3IntronU12AtAc => {
                successors.push((self, 0.0));
                consider_transition(
                    successors,
                    HmmState::UTR3,
                    acceptor_penalty,
                    trans_ctx.cfg.splice.utr3,
                );
            }
        }
    }
}

const PENALTY_SCALE: f64 = 1_000_000.0; // Convert to u64 to avoid FP annoyances

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct HmmEval {
    start_position: usize, // Number of bases produced before this state
    end_position: usize,   // Number of bases produced including this state
    state: HmmState,
    previous_state: HmmState,

    accum_penalty: u64, // Accumulated penalty, target to minimise
    trans_penalty: u64, // Transition penalty into this state
    neg_log_prob: u64   // Negative Log probability of all bases within this state
}

impl HmmEval {
    fn new_root() -> HmmEval {
        HmmEval {
            start_position: 0,
            end_position: 0,
            state: HmmState::Intergenic,
            previous_state: HmmState::Intergenic,
            accum_penalty: 0,
            trans_penalty: 0,
            neg_log_prob: 0
        }
    }

    fn new_successor(
        start_position: usize,
        end_position: usize,
        state: HmmState,
        previous_state: HmmState,
        accum_penalty: u64,
        trans_penalty: u64,
        neg_log_prob: u64
    ) -> HmmEval {
        HmmEval {
            start_position,
            end_position,
            state,
            previous_state,
            accum_penalty,
            //trans_penalty: trans_penalty + 1, // FAKE
            trans_penalty,
            neg_log_prob
        }
    }
}

impl Ord for HmmEval {
    fn cmp(&self, other: &Self) -> Ordering {
        // Order on penalty (lowest), then position (highest)
        other
            .accum_penalty
            .cmp(&self.accum_penalty)
            .then_with(|| self.end_position.cmp(&other.end_position))
        //.then_with(|| self.state.cmp(&other.state))
    }
}

impl PartialOrd for HmmEval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

const MAX_EVALS: u64 = 100_000_000_000;

pub struct PredictionHmm {
    cfg: HmmConfig,
    class_pred_pen: Vec<ClassPredPenalty>,
    phase_pred_pen: Vec<PhasePredPenalty>,

    pred_pen: Vec<PredPenalty>,
    bases_pen: Vec<BasesPenalty>,

    best_eval: Vec<Option<HmmEval>>,

    eval_heap: BinaryHeap<HmmEval>,
}

impl PredictionHmm {
    pub fn new(
        bp_vector: Vec<(Bases, ClassPrediction, PhasePrediction)>,
        cfg: HmmConfig,
    ) -> PredictionHmm {
        let prob_floor = cfg.prob_floor;
        let phase_retain = cfg.phase_retain;

        let mut class_pred_pen = Vec::with_capacity(bp_vector.len());
        let mut phase_pred_pen = Vec::with_capacity(bp_vector.len());
        let mut pred_pen = Vec::with_capacity(bp_vector.len());

        let mut bases_pen = Vec::with_capacity(bp_vector.len());

        for (bases, class_pred, phase_pred) in bp_vector.iter() {
            class_pred_pen.push(ClassPredPenalty::new(class_pred, prob_floor));
            phase_pred_pen.push(PhasePredPenalty::new(phase_pred, prob_floor));
            pred_pen.push(PredPenalty::new(class_pred, phase_pred, prob_floor, phase_retain));

            bases_pen.push(BasesPenalty::new(bases, prob_floor));
        }

        let total_states = (bp_vector.len() + 1) * HMM_STATES;
        let best_eval = vec![None; total_states];

        let eval_heap = BinaryHeap::new();
        PredictionHmm {
            cfg,
            class_pred_pen,
            phase_pred_pen,
            pred_pen,
            bases_pen,
            best_eval,
            eval_heap,
        }
    }

    fn consider_eval(&mut self, eval: HmmEval) {
        let idx = eval.end_position * HMM_STATES + (eval.state as usize);

        let maybe_old_eval = &self.best_eval[idx];

        if let Some(old_eval) = maybe_old_eval {
            if eval.accum_penalty >= old_eval.accum_penalty {
                return;
            }
        }

        self.best_eval[idx] = Some(eval);
        self.eval_heap.push(eval);
    }

    fn is_eval_current(&self, eval: &HmmEval) -> bool {
        let idx = eval.end_position * HMM_STATES + (eval.state as usize);
        let best_eval = &self.best_eval[idx].expect("Eval from heap not in best_eval");

        best_eval == eval
    }

    fn process_eval(&mut self, eval: &HmmEval) {
        if !self.is_eval_current(eval) {
            return;
        }

        let trans_ctx = TransitionContext::new(
            &self.cfg,
            &self.class_pred_pen,
            &self.phase_pred_pen,
            &self.pred_pen,
            &self.bases_pen,
            eval.end_position,
        );

        let mil = self.cfg.minimum_intron_lengths.clone();

        let mut successors = Vec::with_capacity(HMM_STATES);
        eval.state
            .populate_successor_states_and_transition_penalties(&trans_ctx, &mut successors);

        for (next_state, trans_penalty) in successors.into_iter() {
            let mut local_neg_log_prob = trans_penalty;
            let mut local_penalty = trans_penalty;

            let start_position = eval.end_position;
            let end_position = start_position + next_state.get_base_count(&mil);

            if end_position <= self.class_pred_pen.len()
            // Drop 'long' state picked near end
            {
                for pos in start_position..end_position {
                    let (nlg, pen) = next_state.get_state_penalty(
                        &self.class_pred_pen[pos],
                        &self.phase_pred_pen[pos],
                        &self.pred_pen[pos],
                    );
                    local_neg_log_prob += nlg;
                    local_penalty += pen;
                }

                let accum_penalty = eval.accum_penalty + ((local_penalty * PENALTY_SCALE) as u64);
                let scaled_trans_penalty = (trans_penalty * PENALTY_SCALE) as u64;
                let scaled_neg_log_prob = (local_neg_log_prob * PENALTY_SCALE) as u64;

                let next_eval = HmmEval::new_successor(
                    start_position,
                    end_position,
                    next_state,
                    eval.state,
                    accum_penalty,
                    scaled_trans_penalty,
                    scaled_neg_log_prob
                );

                self.consider_eval(next_eval);
            }
        }
    }

    pub fn solve(mut self) -> Option<PredictionHmmSolution> {
        let initial_eval = HmmEval::new_root();
        self.consider_eval(initial_eval);

        let mut evals = 0;

        while evals < MAX_EVALS {
            if let Some(eval) = self.eval_heap.pop() {
                if eval.end_position == self.class_pred_pen.len() {
                    return Some(PredictionHmmSolution::new(self, eval));
                }

                self.process_eval(&eval);
            } else {
                break;
            } // Nothing left to do

            evals += 1;
        }

        panic!("MAX_EVALS exceeded - raise limit or window thresholds");
        //        None
    }
}

pub struct HmmStateRegion {
    start_pos: usize,
    end_pos: usize,
    annotation_label: HmmAnnotationLabel,

    entry_trans_penalty: u64, // Transition penalty into first state in the region
    mid_trans_penalty: u64, // Total transition penalty between each state pair inside the region
    exit_trans_penalty: u64, // Transition penalty out of the final state in the region

    neg_log_prob: u64   // Negative Log probability of all bases within this region
}

impl HmmStateRegion {
    fn new(start_pos: usize, end_pos: usize, annotation_label: HmmAnnotationLabel,
           entry_trans_penalty: u64, mid_trans_penalty: u64, exit_trans_penalty: u64,
           neg_log_prob: u64) -> HmmStateRegion {
        HmmStateRegion {
            start_pos,
            end_pos,
            annotation_label,
            entry_trans_penalty,
            mid_trans_penalty,
            exit_trans_penalty,
            neg_log_prob
        }
    }

    pub fn get_start_pos(&self) -> usize {
        self.start_pos
    }

    pub fn get_end_pos(&self) -> usize {
        self.end_pos
    }

    pub fn get_annotation_label(&self) -> HmmAnnotationLabel {
        self.annotation_label
    }

    pub fn len(&self) -> usize {
        self.end_pos - self.start_pos
    }

    pub fn split_genes(regions: Vec<HmmStateRegion>) -> Vec<(Vec<HmmStateRegion>, usize)> {
        let mut vec_of_vecs = Vec::new();

        let mut current_vec = Vec::new();
        let mut coding_length = 0;

        for region in regions {

            if region.annotation_label == HmmAnnotationLabel::Intergenic {
                if current_vec.len() > 0 {
                    vec_of_vecs.push((current_vec, coding_length));
                    current_vec = Vec::new();
                    coding_length = 0;
                }
            } else {
                if region.annotation_label == HmmAnnotationLabel::Coding {
                    coding_length += region.end_pos - region.start_pos;
                }
            }

            current_vec.push(region);
        }

        if current_vec.len() > 0 {
            vec_of_vecs.push((current_vec, coding_length));
        }

        vec_of_vecs
    }
}

pub struct PredictionHmmSolution {
    hmm: PredictionHmm,
    eval: HmmEval,
}

impl PredictionHmmSolution {
    fn new(hmm: PredictionHmm, eval: HmmEval) -> PredictionHmmSolution {
        PredictionHmmSolution { hmm, eval }
    }

    pub fn trace_regions(&self) -> Vec<HmmStateRegion> {
        let mut regions = Vec::new();

        let mut eval = &self.eval;
        let mut region_end_pos = eval.end_position;

        let mut boundary_trans_penalty = 0u64;
        let mut accum_mid_trans_penalty = 0u64;
        let mut accum_neg_log_prob = 0u64;

        // State positions are inclusive start, exclusive end: state.start_pos = 0, state.end_pos = 1 is the first real state (almost always)
        // A 'dummy root' state exists at the start with start_pos = 0, end_pos = 0, with intergenic state

        // Extract all but the 'start dummy' HmmEval entry
        while eval.end_position > 0 {

            accum_neg_log_prob += eval.neg_log_prob;

            if eval.state.get_annotation_label() == eval.previous_state.get_annotation_label() { // If prev state matches annotation label
                accum_mid_trans_penalty += eval.trans_penalty;

            } else {
                let neg_log_rate = accum_neg_log_prob / ((region_end_pos - eval.start_position) as u64);

                let entry_trans_penalty = eval.trans_penalty;
                let exit_trans_penalty = boundary_trans_penalty;

                vprintln!("Generate a standard HmmStateRegion S: {} E: {} Label: {} Trans: {} {} {} NegLog: {} {}",
                         eval.start_position, region_end_pos, eval.state.get_annotation_label().to_str(),
                         entry_trans_penalty, accum_mid_trans_penalty, exit_trans_penalty,
                         accum_neg_log_prob, neg_log_rate);

                regions.push(HmmStateRegion::new(
                    eval.start_position,
                    region_end_pos,
                    eval.state.get_annotation_label(),
                    entry_trans_penalty,
                    accum_mid_trans_penalty,
                    exit_trans_penalty,
                    accum_neg_log_prob
                ));
                region_end_pos = eval.start_position; // Equivalent to previous_state.end_position

                boundary_trans_penalty = entry_trans_penalty;
                accum_mid_trans_penalty = 0;
                accum_neg_log_prob = 0;
            }

            let prev_position = eval.start_position;
            let idx = prev_position * HMM_STATES + (eval.previous_state as usize);

            eval = self.hmm.best_eval.get(idx).unwrap().as_ref().unwrap();
        }

        // Drain accumulated intergenic region if non-zero length (almost always)
        if region_end_pos > 0 {
            let neg_log_rate = accum_neg_log_prob / (region_end_pos as u64);

            let entry_trans_penalty = eval.trans_penalty;
            let exit_trans_penalty = boundary_trans_penalty;

            vprintln!("Generate a starting Intergenic HmmStateRegion S: {} E: {} Label: {} Trans: {} {} {} NegLog: {} {}",
                     0, region_end_pos, eval.state.get_annotation_label().to_str(),
                     entry_trans_penalty, accum_mid_trans_penalty, exit_trans_penalty,
                     accum_neg_log_prob, neg_log_rate);

            regions.push(HmmStateRegion::new(
                eval.start_position,
                region_end_pos,
                eval.state.get_annotation_label(),
                entry_trans_penalty,
                accum_mid_trans_penalty,
                exit_trans_penalty,
                accum_neg_log_prob
            ));
        }

        regions.reverse();

        regions
    }

    pub fn dump(&self, position: usize) {
        println!(
            "Solution Penalty: {} over {} bp starting at {}",
            self.eval.accum_penalty,
            self.hmm.class_pred_pen.len(),
            position
        );

        let regions = self.trace_regions();

        for region in regions.iter() {
            //            println!("{} to {} is {}", region.start_pos+position+1, region.end_pos+position, region.state.to_str()); // Biologist coordinates

            let mut seq = String::with_capacity(region.end_pos - region.start_pos);
            for idx in region.start_pos..region.end_pos {
                seq.push(self.hmm.bases_pen[idx].as_str());
            }

            println!(
                "{} to {} aka {} to {} is {} - {}",
                region.start_pos,
                region.end_pos,
                region.start_pos + position,
                region.end_pos + position,
                region.annotation_label.to_str(),
                seq
            );
        }

        //        if position > 30000
        //            { panic!("First 30k"); }
        panic!("Show only one region");
    }
}

/*

    base[0]=A
    base[1]=T
    base[2]=G


0   UTR5    - : A
1   Start1  A : T
    Start2  AT : G
    Coding0 ATG : -
    Coding1 ATG - : -
    Coding2 ATG -- : -
    Coding0 ATG --- : T
    StopT   ATG --- T : A
    StopTA  ATG --- TA : A
    Stop3   ATG --- TAA : -
    UTR3    ATG --- TAA -

                           0120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120
    191 to 345 is Coding - ATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCT

                           12012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
    427 to 708 is Coding - TCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAG

                           012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
   917 to 1037 is Coding - AGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTG

                           012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
  1137 to 1527 is Coding - GTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATG

                           012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
  1605 to 1758 is Coding - ATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAG

                         012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
1870 to 2059 is Coding - GTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGT

2059 to 2062 is Stop - TAA

3760 to 3762 is Start
3763 to 3913 is Coding
3996 to 4276 is Coding
4486 to 4605 is Coding
4706 to 5095 is Coding
5174 to 5326 is Coding
5439 to 5627 is Coding
5628 to 5630 is Stop

Chr1	phytozomev10	CDS	3760	3913
Chr1	phytozomev10	CDS	3996	4276
Chr1	phytozomev10	CDS	4486	4605
Chr1	phytozomev10	CDS	4706	5095
Chr1	phytozomev10	CDS	5174	5326
Chr1	phytozomev10	CDS	5439	5630

23136 to 23349 is UTR5
23350 to 23514 is Intron
23515 to 23524 is UTR5
23525 to 24451 is Coding
24542 to 24655 is Coding
24752 to 24962 is Coding
25041 to 25435 is Coding
25524 to 25743 is Coding
25825 to 25997 is Coding
26081 to 26203 is Coding
26292 to 26452 is Coding
26543 to 26776 is Coding
26862 to 27012 is Coding
27099 to 27281 is Coding
27372 to 27533 is Coding
27618 to 27713 is Coding
27803 to 28431 is Coding
28708 to 28805 is Coding
28890 to 29080 is Coding
29193 to 30065 is Coding
30147 to 30311 is Coding
30410 to 30816 is Coding
30902 to 31079 is Coding
31080 to 31200 is UTR3


                                                                    YYYYYYYYNCAG                      YYYYYYYYNCAG
 1234567890123456789012345678901234567890123456789012345678901234567890123456789 012345678901234567890123456789012 34567
 GTATATATATATATATTATGCTTAGTGTCTTTTTTTTTTTTGTTGAAACTATCTAATCATATTTGGTATATATATGTAG ATTCTTGAAGCCTTGACTGCCGCCTCGTGCCAG GAAAC

                                                                                            YTRAC




Chr1    HelixerPost     CDS     8594    8646    .       -       0       ID=Athaliana_Chr1_006990.1.CDS.1;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     8417    8464    .       -       2       ID=Athaliana_Chr1_006990.1.CDS.2;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     8236    8325    .       -       2       ID=Athaliana_Chr1_006990.1.CDS.3;Parent=Athaliana_Chr1_006990.1

Chr1    HelixerPost     CDS     7762    7886    .       -       2       ID=Athaliana_Chr1_006990.1.CDS.4;Parent=Athaliana_Chr1_006990.1

Chr1    HelixerPost     CDS     7564    7649    .       -       1       ID=Athaliana_Chr1_006990.1.CDS.5;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     7384    7450    .       -       0       ID=Athaliana_Chr1_006990.1.CDS.6;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     7159    7232    .       -       1       ID=Athaliana_Chr1_006990.1.CDS.7;Parent=Athaliana_Chr1_006990.1

Chr1    HelixerPost     CDS     6428    6655    .       -       0       ID=Athaliana_Chr1_006991.1.CDS.1;Parent=Athaliana_Chr1_006991.1


Chr1    phytozomev10    CDS     8571    8666    .       -       0       ID=AT1G01020.1.TAIR10.CDS.1;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     8417    8464    .       -       0       ID=AT1G01020.1.TAIR10.CDS.2;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     8236    8325    .       -       0       ID=AT1G01020.1.TAIR10.CDS.3;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7942    7987    .       -       0       ID=AT1G01020.1.TAIR10.CDS.4;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7762    7835    .       -       2       ID=AT1G01020.1.TAIR10.CDS.5;Parent=AT1G01020.1.TAIR10;pacid=19655142

Chr1    phytozomev10    CDS     7564    7649    .       -       0       ID=AT1G01020.1.TAIR10.CDS.6;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7384    7450    .       -       1       ID=AT1G01020.1.TAIR10.CDS.7;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7157    7232    .       -       0       ID=AT1G01020.1.TAIR10.CDS.8;Parent=AT1G01020.1.TAIR10;pacid=19655142

Chr1    phytozomev10    CDS     6915    7069    .       -       2       ID=AT1G01020.1.TAIR10.CDS.9;Parent=AT1G01020.1.TAIR10;pacid=19655142











Chr1    HelixerPost     CDS     32547   32670   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.1;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     32431   32477   .       -       1       ID=Athaliana_Chr1_006988.1.CDS.2;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     32282   32347   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.3;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     32088   32195   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.4;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31933   31998   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.5;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31693   31813   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.6;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31521   31602   .       -       1       ID=Athaliana_Chr1_006988.1.CDS.7;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31382   31424   .       -       2       ID=Athaliana_Chr1_006988.1.CDS.8;Parent=Athaliana_Chr1_006988.1

Chr1    phytozomev10    CDS     32547   32670   .       -       0       ID=AT1G01050.1.TAIR10.CDS.1;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     32431   32459   .       -       2       ID=AT1G01050.1.TAIR10.CDS.2;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     32282   32347   .       -       0       ID=AT1G01050.1.TAIR10.CDS.3;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     32088   32195   .       -       0       ID=AT1G01050.1.TAIR10.CDS.4;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31933   31998   .       -       0       ID=AT1G01050.1.TAIR10.CDS.5;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31693   31813   .       -       0       ID=AT1G01050.1.TAIR10.CDS.6;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31521   31602   .       -       2       ID=AT1G01050.1.TAIR10.CDS.7;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31382   31424   .       -       1       ID=AT1G01050.1.TAIR10.CDS.8;Parent=AT1G01050.1.TAIR10;pacid=19652974




Chr1    HelixerPost     gene    23136   31200   .       +       .       ID=Athaliana_Chr1_000002
Chr1    HelixerPost     mRNA    23136   31200   .       +       .       ID=Athaliana_Chr1_000002.1;Parent=Athaliana_Chr1_000002

Chr1    HelixerPost     exon    23136   23265   .       +       .       ID=Athaliana_Chr1_000002.1.exon.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     five_prime_UTR  23136   23264   .       +       .       ID=Athaliana_Chr1_000002.1.five_prime_UTR.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23265   23265   .       +       0       ID=Athaliana_Chr1_000002.1.CDS.1;Parent=Athaliana_Chr1_000002.1

Chr1    HelixerPost     exon    23267   23339   .       +       .       ID=Athaliana_Chr1_000002.1.exon.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23267   23339   .       +       2       ID=Athaliana_Chr1_000002.1.CDS.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     exon    23515   24451   .       +       .       ID=Athaliana_Chr1_000002.1.exon.3;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23515   24451   .       +       1       ID=Athaliana_Chr1_000002.1.CDS.3;Parent=Athaliana_Chr1_000002.1





Chr1    HelixerPost     mRNA    23136   31200   .       +       .       ID=Athaliana_Chr1_000002.1;Parent=Athaliana_Chr1_000002
Chr1    HelixerPost     exon    23136   23344   .       +       .       ID=Athaliana_Chr1_000002.1.exon.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     five_prime_UTR  23136   23344   .       +       .       ID=Athaliana_Chr1_000002.1.five_prime_UTR.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     exon    23515   24451   .       +       .       ID=Athaliana_Chr1_000002.1.exon.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     five_prime_UTR  23515   23524   .       +       .       ID=Athaliana_Chr1_000002.1.five_prime_UTR.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23525   24451   .       +       0       ID=Athaliana_Chr1_000002.1.CDS.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     exon    24546   24655   .       +       .       ID=Athaliana_Chr1_000002.1.exon.3;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     24546   24655   .       +       0       ID=Athaliana_Chr1_000002.1.CDS.2;Parent=Athaliana_Chr1_000002.1






 */

#[cfg(test)]
mod tests {
    use super::*;

    // --- numeric primitives ---

    /// Inputs below the floor are clamped before the log so penalties stay
    /// finite. `1.0` collapses to a penalty of `0.0`. Tweak these constants
    /// at your peril.
    #[test]
    fn raw_pred_to_neg_log_prob_clamps_to_floor_and_logs() {
        let floor = 1e-9_f64;
        let out = raw_pred_to_neg_log_prob(&[0.0, 1e-12, 0.5, 1.0], floor);
        let floored = -(floor.log2());
        assert!((out[0] - floored).abs() < 1e-12, "0.0 should be clamped to floor");
        assert!((out[1] - floored).abs() < 1e-12, "1e-12 should be clamped to floor");
        assert!((out[2] - 1.0).abs() < 1e-12, "-log2(0.5) == 1.0");
        assert!(out[3].abs() < 1e-12, "-log2(1.0) == 0.0");
    }

    #[test]
    fn neg_log_prob_to_penalty_subtracts_minimum() {
        let pen = neg_log_prob_to_penalty(&[5.0_f64, 2.0, 3.0, 7.0]);
        assert_eq!(pen, [3.0, 0.0, 1.0, 5.0]);
    }

    /// The "winning" base in a one-hot Bases has penalty 0, so `as_str`
    /// returns it. Two bases pre-normalisation -> after penalty conversion
    /// the wrong one still differs from zero by more than APPROX_ZERO.
    #[test]
    fn bases_penalty_as_str_returns_winner() {
        // Pure-A
        let pen = BasesPenalty::new(&Bases::new([0.0, 1.0, 0.0, 0.0]), 1e-9);
        assert_eq!(pen.as_str(), 'A');
        // Pure-T
        let pen = BasesPenalty::new(&Bases::new([0.0, 0.0, 1.0, 0.0]), 1e-9);
        assert_eq!(pen.as_str(), 'T');
        // Pure-G
        let pen = BasesPenalty::new(&Bases::new([0.0, 0.0, 0.0, 1.0]), 1e-9);
        assert_eq!(pen.as_str(), 'G');
        // Pure-C
        let pen = BasesPenalty::new(&Bases::new([1.0, 0.0, 0.0, 0.0]), 1e-9);
        assert_eq!(pen.as_str(), 'C');
    }

    // --- PredPenalty phase-blending math ---

    /// When the model emits zero for all three coding-phase channels the
    /// constructor falls back to an even coding/3 split per phase, then
    /// blends with the dilution target. With phase_retain == 0 the result
    /// must equal pure-coding probability for each phase.
    /// (Tolerance is loose because raw probs go through f32->f64.)
    #[test]
    fn pred_penalty_handles_zero_phase_with_phase_retain_zero() {
        let class = ClassPrediction::new([0.7, 0.0, 0.3, 0.0]); // 0.3 coding
        let phase = PhasePrediction::new([1.0, 0.0, 0.0, 0.0]); // all in non-coding
        let pp = PredPenalty::new(&class, &phase, 1e-9, 0.0);

        // raw_probs layout (post-blend): [intergenic, utr, phase0, phase1, phase2, intron]
        // With phase_retain=0 the dilution_target (= coding) wins; phase0/1/2 = coding = 0.3.
        let intergenic_nlp = -(0.7_f32 as f64).log2();
        let coding_nlp = -(0.3_f32 as f64).log2();
        assert!((pp.neg_log_prob[0] - intergenic_nlp).abs() < 1e-6);
        assert!((pp.neg_log_prob[2] - coding_nlp).abs() < 1e-6);
        assert!((pp.neg_log_prob[3] - coding_nlp).abs() < 1e-6);
        assert!((pp.neg_log_prob[4] - coding_nlp).abs() < 1e-6);
    }

    /// With phase_retain == 1 (trust the predicted phase fully) and a perfect
    /// phase-0 prediction, the three phase channels should split coding
    /// fully into phase 0.
    #[test]
    fn pred_penalty_full_phase_retain_routes_coding_into_predicted_phase() {
        let class = ClassPrediction::new([0.1, 0.0, 0.9, 0.0]); // 0.9 coding
        let phase = PhasePrediction::new([0.0, 1.0, 0.0, 0.0]); // phase 0
        let pp = PredPenalty::new(&class, &phase, 1e-9, 1.0);

        let phase0_nlp = -(0.9_f32 as f64).log2();
        // Phase 1 and 2 receive ~0 probability — penalty bottoms out at the floor.
        let floor_nlp = -(1e-9_f64).log2();
        assert!((pp.neg_log_prob[2] - phase0_nlp).abs() < 1e-6);
        assert!((pp.neg_log_prob[3] - floor_nlp).abs() < 1e-6);
        assert!((pp.neg_log_prob[4] - floor_nlp).abs() < 1e-6);
    }

    // --- HmmStateRegion::split_genes ---

    fn region(start: usize, end: usize, label: HmmAnnotationLabel) -> HmmStateRegion {
        HmmStateRegion::new(start, end, label, 0, 0, 0, 0)
    }

    /// `split_genes` flushes the current bucket each time it hits an
    /// intergenic *after* the bucket already has content. The intergenic
    /// itself is pushed into the *next* bucket. So a leading intergenic
    /// stays inside the first gene, and a trailing intergenic becomes a
    /// final bucket of its own.
    #[test]
    fn split_genes_partitions_on_intergenic_and_sums_coding() {
        let regions = vec![
            region(0, 100, HmmAnnotationLabel::Intergenic),
            region(100, 110, HmmAnnotationLabel::UTR5),
            region(110, 113, HmmAnnotationLabel::Start),
            region(113, 200, HmmAnnotationLabel::Coding),
            region(200, 203, HmmAnnotationLabel::Stop),
            region(203, 210, HmmAnnotationLabel::UTR3),
            region(210, 500, HmmAnnotationLabel::Intergenic),
            region(500, 600, HmmAnnotationLabel::Coding),
            region(600, 700, HmmAnnotationLabel::Intergenic),
        ];
        let genes = HmmStateRegion::split_genes(regions);
        assert_eq!(genes.len(), 3);

        // Bucket 0: leading Intergenic + first full gene structure.
        // Coding length = 200 - 113 = 87.
        assert_eq!(genes[0].1, 87);
        assert_eq!(genes[0].0.len(), 6);
        assert_eq!(
            genes[0].0[0].get_annotation_label(),
            HmmAnnotationLabel::Intergenic
        );

        // Bucket 1: middle Intergenic + Coding(500,600). Coding length = 100.
        assert_eq!(genes[1].1, 100);
        assert_eq!(genes[1].0.len(), 2);
        assert_eq!(
            genes[1].0[0].get_annotation_label(),
            HmmAnnotationLabel::Intergenic
        );

        // Bucket 2: trailing Intergenic alone, no coding.
        assert_eq!(genes[2].1, 0);
        assert_eq!(genes[2].0.len(), 1);
        assert_eq!(
            genes[2].0[0].get_annotation_label(),
            HmmAnnotationLabel::Intergenic
        );
    }

    /// Regions with no intergenic at all stay in a single bucket. Coding
    /// length still totals only Coding regions, not Start/Stop/UTR.
    #[test]
    fn split_genes_no_intergenic_single_bucket() {
        let regions = vec![
            region(0, 5, HmmAnnotationLabel::UTR5),
            region(5, 8, HmmAnnotationLabel::Start),
            region(8, 50, HmmAnnotationLabel::Coding),
            region(50, 53, HmmAnnotationLabel::Stop),
        ];
        let genes = HmmStateRegion::split_genes(regions);
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].1, 42); // only the 8..50 Coding region counts
        assert_eq!(genes[0].0.len(), 4);
    }

    #[test]
    fn split_genes_empty_input_yields_no_genes() {
        let genes = HmmStateRegion::split_genes(Vec::new());
        assert!(genes.is_empty());
    }
}
