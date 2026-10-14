use crate::results::conv::{Bases, ClassPrediction, PhasePrediction};
use std::collections::vec_deque::Iter;
use std::collections::VecDeque;

/// The per-base item type the sliding window consumes. Any iterator yielding
/// this triple can drive the window — the HDF5-backed `BasePredictionIterator`
/// in production, or a `Vec<...>::into_iter()` in tests.
type BpItem = (Bases, ClassPrediction, PhasePrediction);

pub struct BasePredictionWindow<I: Iterator<Item = BpItem>> {
    bp_iter: I,

    window_size: usize,
    scale: f32,
    window_total: u64,

    window: VecDeque<BpItem>,
    position: usize,
}

impl<I: Iterator<Item = BpItem>> BasePredictionWindow<I> {
    pub fn new(bp_iter: I, window_size: usize, scale: f32) -> Option<BasePredictionWindow<I>> {
        let window = VecDeque::with_capacity(window_size);

        let mut bp_window = BasePredictionWindow {
            bp_iter,
            window_size,
            scale,
            window_total: 0,
            window,
            position: 0,
        };
        bp_window.fill_window();

        Some(bp_window)
    }

    fn fill_window(&mut self) -> bool {
        while self.window.len() < self.window_size {
            if !self.push() {
                return false;
            }
        }
        true
    }

    fn push(&mut self) -> bool {
        let maybe_next = self.bp_iter.next();

        if let Some((bases, class_pred, phase_pred)) = maybe_next {
            self.window_total += (class_pred.get_genic() * self.scale) as u64;
            self.window.push_back((bases, class_pred, phase_pred));
            true
        } else {
            false
        }
    }

    fn pop(&mut self) -> Option<BpItem> {
        let maybe_next = self.window.pop_front();

        if let Some((_bases, class_pred, _phase_pred)) = &maybe_next {
            self.window_total -= (class_pred.get_genic() * self.scale) as u64;
            self.position += 1;
        }

        maybe_next
    }

    pub fn get_window_total(&self) -> u64 {
        self.window_total
    }

    pub fn is_window_full(&self) -> bool {
        self.window.len() == self.window_size
    }

    pub fn get_window_iter(&self) -> Iter<'_, BpItem> { self.window.iter() }
}

pub struct BasePredictionWindowThresholdScanner<I: Iterator<Item = BpItem>> {
    bp_window: BasePredictionWindow<I>,
    edge_threshold: u64,
}

impl<I: Iterator<Item = BpItem>> BasePredictionWindowThresholdScanner<I> {
    pub fn new(
        bp_window: BasePredictionWindow<I>,
        edge_threshold: f32,
    ) -> BasePredictionWindowThresholdScanner<I> {
        let threshold = (edge_threshold * bp_window.scale * (bp_window.window_size as f32)) as u64;

        BasePredictionWindowThresholdScanner {
            bp_window,
            edge_threshold: threshold,
        }
    }

    fn scan_for_start(&mut self) -> bool {
        while self.bp_window.is_window_full()
            && self.bp_window.get_window_total() < self.edge_threshold
        {
            self.bp_window.pop(); // Pop and discard
            self.bp_window.push();
        }

        self.bp_window.is_window_full()
    }

    fn accumulate_above_threshold(&mut self) -> (Vec<BpItem>, Vec<u64>, usize, u64) {
        let mut accum = Vec::new();
        let mut total_accum = Vec::new();

        if !self.bp_window.is_window_full()
            || self.bp_window.get_window_total() < self.edge_threshold
        {
            panic!("Accumulate called with window not past threshold");
        }

        let mut peak = 0;
        let position = self.bp_window.position;

        while self.bp_window.is_window_full()
            && self.bp_window.get_window_total() >= self.edge_threshold
        {
            let total = self.bp_window.get_window_total();
            total_accum.push(total);
            peak = std::cmp::max(peak, total);

            let (bases, class_pred, phase_pred) = self.bp_window.pop().unwrap();
            accum.push((bases, class_pred, phase_pred));

            self.bp_window.push();
        }

        for (b, cp, pp) in self.bp_window.get_window_iter() {
            accum.push((*b, *cp, *pp))
        }

        if self.bp_window.is_window_full()
        // If window is still full, the last element crossed below the threshold, remove it
        {
            accum.pop();
        }

        (accum, total_accum, position, peak)
    }
}

const THRESHOLD_SCALE: f32 = 1_000_000.0;

pub struct BasePredictionWindowThresholdIterator<I: Iterator<Item = BpItem>> {
    bp_scanner: BasePredictionWindowThresholdScanner<I>,
    peak_threshold: u64,
    peak_scale: f32,
}

impl<I: Iterator<Item = BpItem>> BasePredictionWindowThresholdIterator<I> {
    pub fn new(
        bp_iter: I,
        window_size: usize,
        edge_threshold: f32,
        peak_threshold: f32,
    ) -> Option<BasePredictionWindowThresholdIterator<I>> {
        let bp_window = BasePredictionWindow::new(bp_iter, window_size, THRESHOLD_SCALE)?;
        let bp_scanner = BasePredictionWindowThresholdScanner::new(bp_window, edge_threshold);

        let window_size = window_size as f32;

        let peak_threshold = (peak_threshold * THRESHOLD_SCALE * window_size) as u64;
        let peak_scale = 1.0 / (THRESHOLD_SCALE * window_size);

        Some(BasePredictionWindowThresholdIterator {
            bp_scanner,
            peak_threshold,
            peak_scale,
        })
    }
}

impl<I: Iterator<Item = BpItem>> Iterator for BasePredictionWindowThresholdIterator<I> {
    type Item = (Vec<BpItem>, Vec<u64>, usize, f32);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if !self.bp_scanner.scan_for_start() {
                return None;
            }

            let (bp_accum, total_accum, position, peak) =
                self.bp_scanner.accumulate_above_threshold();

            if peak > self.peak_threshold {
                let peak = (peak as f32) * self.peak_scale;
                return Some((bp_accum, total_accum, position, peak));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic per-base trace from a list of genic scores.
    /// `intergenic = 1 - genic`; the remaining channels and phase are zero,
    /// which is enough for the windowing logic (only `get_genic()` is read).
    fn trace(scores: &[f32]) -> Vec<BpItem> {
        scores
            .iter()
            .map(|&g| {
                let class = ClassPrediction::new([1.0 - g, 0.0, g, 0.0]);
                let phase = PhasePrediction::new([1.0, 0.0, 0.0, 0.0]);
                let bases = Bases::new([0.0, 1.0, 0.0, 0.0]); // 'A', irrelevant here
                (bases, class, phase)
            })
            .collect()
    }

    /// A flat all-intergenic trace (score 0) never crosses the edge threshold
    /// so no candidate windows are emitted.
    #[test]
    fn no_candidates_when_below_edge_threshold() {
        let bp = trace(&vec![0.0; 500]);
        let it = BasePredictionWindowThresholdIterator::new(bp.into_iter(), 100, 0.1, 0.8).unwrap();
        assert_eq!(it.count(), 0);
    }

    /// A region whose mean exceeds `edge_threshold` but whose peak never
    /// reaches `peak_threshold` is filtered out by the peak gate.
    #[test]
    fn peak_threshold_filters_low_peak_regions() {
        // 200 bases at 0.2 (above edge_threshold=0.1, below peak_threshold=0.8)
        // sandwiched in intergenic.
        let mut scores = vec![0.0_f32; 100];
        scores.extend(vec![0.2_f32; 200]);
        scores.extend(vec![0.0_f32; 100]);

        let it = BasePredictionWindowThresholdIterator::new(
            trace(&scores).into_iter(),
            100,
            0.1,
            0.8,
        )
        .unwrap();
        assert_eq!(it.count(), 0, "peak below threshold must be rejected");
    }

    /// A clear genic block — long enough above edge_threshold and with a peak
    /// above peak_threshold — must be emitted exactly once.
    #[test]
    fn single_candidate_for_clear_genic_block() {
        // 200 intergenic, 300 high-genic, 200 intergenic.
        let mut scores = vec![0.0_f32; 200];
        scores.extend(vec![0.95_f32; 300]);
        scores.extend(vec![0.0_f32; 200]);

        let candidates: Vec<_> = BasePredictionWindowThresholdIterator::new(
            trace(&scores).into_iter(),
            100,
            0.1,
            0.8,
        )
        .unwrap()
        .collect();

        assert_eq!(candidates.len(), 1);
        let (_bp_vec, _total_vec, start_pos, peak) = &candidates[0];
        // Start lands somewhere inside the leading intergenic — the window
        // reports `position` from `BasePredictionWindow::position`, which
        // advances as the window slides forward. We don't pin an exact
        // value (depends on the exact crossing point) but it must precede
        // the high-genic block.
        assert!(*start_pos < 300, "start should be before the genic block end");
        assert!(*peak >= 0.8, "peak must clear the peak threshold");
    }

    /// Two well-separated genic blocks produce two candidate windows.
    #[test]
    fn two_separated_candidates() {
        let mut scores = vec![0.0_f32; 200];
        scores.extend(vec![0.95_f32; 300]);
        scores.extend(vec![0.0_f32; 400]); // wide intergenic gap
        scores.extend(vec![0.95_f32; 300]);
        scores.extend(vec![0.0_f32; 200]);

        let candidates: Vec<_> = BasePredictionWindowThresholdIterator::new(
            trace(&scores).into_iter(),
            100,
            0.1,
            0.8,
        )
        .unwrap()
        .collect();

        assert_eq!(candidates.len(), 2);
        // Each candidate must clear the peak gate.
        for (_, _, _, peak) in &candidates {
            assert!(*peak >= 0.8);
        }
        // The two candidates must be in start-position order.
        assert!(candidates[0].2 < candidates[1].2);
    }

    /// An input shorter than window_size cannot fill the window, so the
    /// iterator yields nothing without panicking.
    #[test]
    fn short_input_yields_no_candidates() {
        let it = BasePredictionWindowThresholdIterator::new(
            trace(&vec![0.95; 50]).into_iter(),
            100,
            0.1,
            0.8,
        );
        // `new` returns Some even with a too-short iter (window simply never
        // fills); the iterator then immediately yields None.
        let it = it.expect("constructor returns Some regardless");
        assert_eq!(it.count(), 0);
    }
}
