#[derive(Debug, Clone)]
pub struct AnchorSearch {
    fixed: Vec<u8>,
    anchor: Vec<u8>,
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

        Some(Self {
            fixed: fixed.to_vec(),
            anchor,
            anchor_offset,
            max_mismatches,
        })
    }

    #[inline]
    fn mismatches(a: &[u8], b: &[u8]) -> usize {
        debug_assert_eq!(a.len(), b.len());
        a.iter().zip(b).filter(|(x, y)| x != y).count()
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
        if read.len() < self.anchor.len() || from >= read.len() {
            return None;
        }

        // primer_start = anchor_start - anchor_offset.  Starting the anchor
        // scan here guarantees that returned primer starts are >= `from`.
        // For from == 0 keep the leading, saturating candidates accepted by
        // the old implementation.
        let first_anchor = if from == 0 {
            0
        } else {
            from.checked_add(self.anchor_offset)?
        };

        let last_anchor = read.len().checked_sub(self.anchor.len())?;
        if first_anchor > last_anchor {
            return None;
        }

        for anchor_start in first_anchor..=last_anchor {
            let obs = &read[anchor_start..anchor_start + self.anchor.len()];

            if Self::mismatches(obs, &self.anchor) <= self.max_mismatches {
                let primer_start = anchor_start.saturating_sub(self.anchor_offset);
                if primer_start >= from {
                    return Some(primer_start);
                }
            }
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
