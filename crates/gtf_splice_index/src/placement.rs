use serde::{Deserialize, Serialize};

/// Stable integer identifier used by the placement engine.
/// Biological meaning is supplied by the owning index.
pub type PlacementId = usize;

/// Per-chromosome bucket index: bin -> placed model ids.
///
/// This is a pre-filter only: it returns candidate placed model IDs that overlap bins.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChrBuckets {
    pub bin_width: u32,
    pub bins: Vec<Vec<PlacementId>>,
    pub max_end: u32,
}

impl ChrBuckets {
    pub fn new(bin_width: u32) -> Self {
        Self {
            bin_width,
            bins: Vec::new(),
            max_end: 0,
        }
    }

    fn ensure_len_for_end(&mut self, end0: u32) {
        self.max_end = self.max_end.max(end0);

        let need_bins =
            ((self.max_end as u64 + self.bin_width as u64 - 1) / self.bin_width as u64) as usize;
        if self.bins.len() < need_bins {
            self.bins.resize_with(need_bins, Vec::new);
        }
    }

    pub fn add_span(&mut self, model_id: PlacementId, start0: u32, end0: u32) {
        if end0 <= start0 {
            return;
        }

        self.ensure_len_for_end(end0);

        let b0 = (start0 / self.bin_width) as usize;
        let b1 = ((end0.saturating_sub(1)) / self.bin_width) as usize;

        for b in b0..=b1 {
            self.bins[b].push(model_id);
        }
    }

    pub fn finalize_by_start(&mut self, span_start: &[u32], span_end: &[u32]) {
        for bin in &mut self.bins {
            // Sort by (start, end, id) to make it deterministic.
            bin.sort_unstable_by(|a, b| {
                let sa = span_start[*a];
                let sb = span_start[*b];
                match sa.cmp(&sb) {
                    std::cmp::Ordering::Equal => {
                        let ea = span_end[*a];
                        let eb = span_end[*b];
                        match ea.cmp(&eb) {
                            std::cmp::Ordering::Equal => a.cmp(b),
                            other => other,
                        }
                    }
                    other => other,
                }
            });

            // Now adjacent duplicates are guaranteed adjacent (because equal ids compare equal)
            bin.dedup();
        }
    }
}


#[inline]
pub(crate) fn partition_point<T, F>(slice: &[T], mut pred: F) -> usize
where
    F: FnMut(&T) -> bool,
{
    let mut left = 0usize;
    let mut right = slice.len();
    while left < right {
        let mid = left + (right - left) / 2;
        if pred(&slice[mid]) { left = mid + 1; } else { right = mid; }
    }
    left
}
