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
        self.identify_all_cell_starts(read).into_iter().next()
    }

    pub fn identify_all_cell_starts(&self, read: &[u8]) -> Vec<usize> {
        let mut starts = Vec::new();

        if read.len() < self.anchor.len() {
            return starts;
        }

        for anchor_start in 0..=read.len() - self.anchor.len() {
            let obs = &read[anchor_start..anchor_start + self.anchor.len()];

            if Self::mismatches(obs, &self.anchor) <= self.max_mismatches {
                let primer_start = anchor_start.saturating_sub(self.anchor_offset);
                starts.push(primer_start);
            }
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
