use onehot_dna::OneHotSequence;

#[derive(Debug, Clone)]
pub struct AnchorSearch {
    fixed: Vec<u8>,
    anchor: Vec<u8>,
    anchor_onehot: OneHotSequence,
    anchor_offset: usize,
    max_mismatches: usize,
}

impl AnchorSearch {
    pub fn new(fixed: &[u8], max_mismatches: usize) -> Option<Self> {
        if fixed.len() < 8 {
            return None;
        }

        // skip noisy read-start bases if possible
        let anchor_offset = if fixed.len() >= 12 { 3 } else { 0 };

        let anchor = fixed[anchor_offset..].to_vec();

        let anchor_onehot = OneHotSequence::from_iupac_bytes(&anchor);
        Some(Self { fixed: fixed.to_vec(), anchor, anchor_onehot, anchor_offset, max_mismatches })
    }

    pub fn identify_cell_start(&self, read: &[u8]) -> Option<usize> {
        self.find_next_cell_start(read, 0)
    }

    /// Find the next candidate primer start at or after `from`.
    ///
    /// This is intentionally streaming: callers that only need the next
    /// candidate do not have to allocate a `Vec` containing every anchor hit
    /// in the read.
    pub fn find_next_cell_start(&self, read: &[u8], from: usize) -> Option<usize> {
        let observed = OneHotSequence::from_iupac_bytes(read);
        self.find_next_cell_start_packed(&observed, from)
    }

    /// Packed variant used by the detector so each read is encoded exactly once.
    pub fn find_next_cell_start_packed(
        &self,
        observed: &OneHotSequence,
        from: usize,
    ) -> Option<usize> {
        if observed.len() < self.anchor.len() || from >= observed.len() {
            return None;
        }

        let first_anchor = if from == 0 {
            0
        } else {
            from.checked_add(self.anchor_offset)?
        };
        let last_anchor = observed.len().checked_sub(self.anchor.len())?;
        if first_anchor > last_anchor { return None; }

        // Eight compatible bases are the fast gate. Four packed bytes give
        // canonical anchors ~1/65,536 random selectivity before the complete
        // anchor/mismatch policy is evaluated.
        const FAST_SEED: usize = 8;
        let seed_len = FAST_SEED.min(self.anchor.len());
        let mut anchor_start = first_anchor;
        while anchor_start <= last_anchor {
            let candidate = observed.find_next_compatible_seed_with_mismatches(
                &self.anchor_onehot, anchor_start, seed_len, self.max_mismatches,
            )?;
            if candidate > last_anchor { return None; }
            let Some((informative, compatible)) = observed.compatibility_counts(
                candidate, &self.anchor_onehot, 0, self.anchor.len(),
            ) else { return None; };
            let mismatches = self.anchor.len().saturating_sub(compatible);
            if informative == self.anchor.len() && mismatches <= self.max_mismatches {
                let primer_start = candidate.saturating_sub(self.anchor_offset);
                if primer_start >= from { return Some(primer_start); }
            }
            anchor_start = candidate.saturating_add(1);
        }
        None
    }

    /// Find a FIXED anchor only within an explicit primer-start window.
    ///
    /// `min_start` and `max_start` are primer-start coordinates, not anchor
    /// coordinates. This is the fast path for `SEARCH:a..b + FIXED:...`: the
    /// SEARCH range bounds one OneHot anchor scan instead of retrying the full
    /// primer at every offset.
    pub fn find_cell_start_in_range_packed(
        &self,
        observed: &OneHotSequence,
        min_start: usize,
        max_start: usize,
    ) -> Option<usize> {
        if min_start > max_start || observed.len() < self.anchor.len() {
            return None;
        }

        let first_anchor = min_start.checked_add(self.anchor_offset)?;
        let max_anchor_from_read = observed.len().checked_sub(self.anchor.len())?;
        let last_anchor = max_start
            .checked_add(self.anchor_offset)
            .unwrap_or(usize::MAX)
            .min(max_anchor_from_read);
        if first_anchor > last_anchor {
            return None;
        }

        const FAST_SEED: usize = 8;
        let seed_len = FAST_SEED.min(self.anchor.len());
        let mut anchor_start = first_anchor;
        while anchor_start <= last_anchor {
            let candidate = observed.find_next_compatible_seed_with_mismatches(
                &self.anchor_onehot,
                anchor_start,
                seed_len,
                self.max_mismatches,
            )?;
            if candidate > last_anchor {
                return None;
            }

            let Some((informative, compatible)) = observed.compatibility_counts(
                candidate,
                &self.anchor_onehot,
                0,
                self.anchor.len(),
            ) else {
                return None;
            };
            let mismatches = self.anchor.len().saturating_sub(compatible);
            if informative == self.anchor.len() && mismatches <= self.max_mismatches {
                let primer_start = candidate.checked_sub(self.anchor_offset)?;
                if primer_start >= min_start && primer_start <= max_start {
                    return Some(primer_start);
                }
            }
            anchor_start = candidate.saturating_add(1);
        }
        None
    }

    pub fn identify_all_cell_starts(&self, read: &[u8]) -> Vec<usize> {
        let mut starts = Vec::new();
        let mut from = 0usize;

        while let Some(start) = self.find_next_cell_start(read, from) {
            starts.push(start);
            let Some(next) = start.checked_add(1) else {
                break;
            };
            from = next;
        }

        starts
    }

    pub fn anchor_len(&self) -> usize {
        self.anchor.len()
    }

    pub fn anchor_offset(&self) -> usize {
        self.anchor_offset
    }

    pub fn fixed(&self) -> &[u8] {
        &self.fixed
    }
}
