//alignement.rs

use crate::IntToStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_penalty: i32,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_score: -1,
            gap_penalty: -2,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentResult {
    pub score: i32,
    pub matches: usize,
    pub mismatches: usize,
    pub gaps: usize,
    pub aligned_self: String,
    pub aligned_other: String,
    pub normalized_distance: f32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExactOverlap {
    /// Offset of `other` relative to `self`.
    ///
    /// offset = 0   => same start
    /// offset > 0   => other starts to the right within/after self
    /// offset < 0   => other starts to the left of self
    pub offset: isize,
    pub overlap: usize,
    pub mismatches: usize,
    pub self_start: usize,
    pub other_start: usize,
}

impl IntToStr {
    #[inline]
    pub fn len_bases(&self) -> usize {
        self.size
    }

    #[inline]
    pub fn base_code_at(&self, pos: usize) -> Option<u8> {
        if pos >= self.size {
            return None;
        }
        let byte_idx = pos / 4;
        let within = pos % 4;
        let shift = within * 2;
        self.u8_encoded.get(byte_idx).map(|b| (b >> shift) & 0b11)
    }

    #[inline]
    pub fn dna_char_at(&self, pos: usize) -> Option<char> {
        match self.base_code_at(pos)? {
            0 => Some('A'),
            1 => Some('C'),
            2 => Some('G'),
            3 => Some('T'),
            _ => None,
        }
    }

    #[inline]
    fn push_base_char(dst: &mut String, base: u8) {
        dst.push(match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        });
    }

    #[inline]
    fn push_gap(dst: &mut String) {
        dst.push('-');
    }

    /// Count mismatches across a fixed ungapped overlap window.
    pub fn overlap_mismatches(
        &self,
        self_start: usize,
        other: &Self,
        other_start: usize,
        overlap_len: usize,
    ) -> Option<usize> {
        if self_start + overlap_len > self.size || other_start + overlap_len > other.size {
            return None;
        }

        let mut mismatches = 0usize;
        for i in 0..overlap_len {
            if self.base_code_at(self_start + i)? != other.base_code_at(other_start + i)? {
                mismatches += 1;
            }
        }
        Some(mismatches)
    }

    /// Return overlap details for one concrete offset if the overlap is long enough
    /// and below the mismatch-fraction cutoff.
    pub fn exact_overlap_with_offset(
        &self,
        other: &Self,
        offset: isize,
        min_overlap: usize,
        max_mismatch_fraction: f32,
    ) -> Option<ExactOverlap> {
        let self_len = self.size as isize;
        let other_len = other.size as isize;

        let ov_start = 0isize.max(offset);
        let ov_end = self_len.min(offset + other_len);

        if ov_end <= ov_start {
            return None;
        }

        let overlap = (ov_end - ov_start) as usize;
        if overlap < min_overlap {
            return None;
        }

        let self_start = ov_start as usize;
        let other_start = if offset < 0 {
            (-offset) as usize
        } else {
            0usize
        };

        let mismatches = self.overlap_mismatches(self_start, other, other_start, overlap)?;
        let frac = mismatches as f32 / overlap as f32;

        if frac <= max_mismatch_fraction {
            Some(ExactOverlap {
                offset,
                overlap,
                mismatches,
                self_start,
                other_start,
            })
        } else {
            None
        }
    }

    /// Search all ungapped relative offsets and return the best overlap.
    ///
    /// Preference order:
    /// 1. longer overlap
    /// 2. fewer mismatches
    /// 3. smaller absolute offset
    pub fn best_exact_overlap(
        &self,
        other: &Self,
        min_overlap: usize,
        max_mismatch_fraction: f32,
    ) -> Option<ExactOverlap> {
        if self.size == 0 || other.size == 0 {
            return None;
        }

        let min_offset = -(other.size as isize) + min_overlap as isize;
        let max_offset = self.size as isize - min_overlap as isize;

        let mut best: Option<ExactOverlap> = None;

        for offset in min_offset..=max_offset {
            if let Some(hit) =
                self.exact_overlap_with_offset(other, offset, min_overlap, max_mismatch_fraction)
            {
                let replace = match best {
                    None => true,
                    Some(prev) => {
                        hit.overlap > prev.overlap
                            || (hit.overlap == prev.overlap && hit.mismatches < prev.mismatches)
                            || (hit.overlap == prev.overlap
                                && hit.mismatches == prev.mismatches
                                && hit.offset.abs() < prev.offset.abs())
                    }
                };

                if replace {
                    best = Some(hit);
                }
            }
        }

        best
    }

    pub fn needleman_wunsch(&self, other: &Self) -> AlignmentResult {
        self.needleman_wunsch_with(other, AlignmentConfig::default())
    }

    pub fn needleman_wunsch_distance(&self, other: &Self) -> f32 {
        self.needleman_wunsch(other).normalized_distance
    }

    pub fn needleman_wunsch_with(&self, other: &Self, cfg: AlignmentConfig) -> AlignmentResult {
        let rows = self.size + 1;
        let cols = other.size + 1;

        if self.size == 0 && other.size == 0 {
            return AlignmentResult {
                score: 0,
                matches: 0,
                mismatches: 0,
                gaps: 0,
                aligned_self: String::new(),
                aligned_other: String::new(),
                normalized_distance: 0.0,
            };
        }

        let mut dp = vec![vec![0i32; cols]; rows];
        let mut trace = vec![vec![0u8; cols]; rows];
        // 0 = diag, 1 = up, 2 = left

        for i in 1..rows {
            dp[i][0] = dp[i - 1][0] + cfg.gap_penalty;
            trace[i][0] = 1;
        }

        for j in 1..cols {
            dp[0][j] = dp[0][j - 1] + cfg.gap_penalty;
            trace[0][j] = 2;
        }

        for i in 1..rows {
            for j in 1..cols {
                let a = self.base_code_at(i - 1).unwrap();
                let b = other.base_code_at(j - 1).unwrap();

                let diag = dp[i - 1][j - 1]
                    + if a == b {
                        cfg.match_score
                    } else {
                        cfg.mismatch_score
                    };
                let up = dp[i - 1][j] + cfg.gap_penalty;
                let left = dp[i][j - 1] + cfg.gap_penalty;

                if diag >= up && diag >= left {
                    dp[i][j] = diag;
                    trace[i][j] = 0;
                } else if up >= left {
                    dp[i][j] = up;
                    trace[i][j] = 1;
                } else {
                    dp[i][j] = left;
                    trace[i][j] = 2;
                }
            }
        }

        let mut i = self.size;
        let mut j = other.size;

        let mut aligned_self_rev = String::new();
        let mut aligned_other_rev = String::new();

        let mut matches = 0usize;
        let mut mismatches = 0usize;
        let mut gaps = 0usize;

        while i > 0 || j > 0 {
            let dir = if i == 0 {
                2
            } else if j == 0 {
                1
            } else {
                trace[i][j]
            };

            match dir {
                0 => {
                    let a = self.base_code_at(i - 1).unwrap();
                    let b = other.base_code_at(j - 1).unwrap();

                    Self::push_base_char(&mut aligned_self_rev, a);
                    Self::push_base_char(&mut aligned_other_rev, b);

                    if a == b {
                        matches += 1;
                    } else {
                        mismatches += 1;
                    }

                    i -= 1;
                    j -= 1;
                }
                1 => {
                    let a = self.base_code_at(i - 1).unwrap();
                    Self::push_base_char(&mut aligned_self_rev, a);
                    Self::push_gap(&mut aligned_other_rev);
                    gaps += 1;
                    i -= 1;
                }
                2 => {
                    let b = other.base_code_at(j - 1).unwrap();
                    Self::push_gap(&mut aligned_self_rev);
                    Self::push_base_char(&mut aligned_other_rev, b);
                    gaps += 1;
                    j -= 1;
                }
                _ => unreachable!(),
            }
        }

        let aligned_self: String = aligned_self_rev.chars().rev().collect();
        let aligned_other: String = aligned_other_rev.chars().rev().collect();

        let aln_len = matches + mismatches + gaps;
        let normalized_distance = if aln_len == 0 {
            0.0
        } else {
            (mismatches + gaps) as f32 / aln_len as f32
        };

        AlignmentResult {
            score: dp[self.size][other.size],
            matches,
            mismatches,
            gaps,
            aligned_self,
            aligned_other,
            normalized_distance,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::IntToStr;

    #[test]
    fn test_base_access_roundtrip_is_sane() {
        let seq = IntToStr::new(b"ACGTTCAG");
        let recovered: String = (0..seq.len_bases())
            .map(|i| seq.dna_char_at(i).unwrap())
            .collect();

        assert_eq!(recovered, "ACGTTCAG");
    }

    #[test]
    fn test_needleman_wunsch_prefers_single_gap_for_internal_insertion() {
        let a = IntToStr::new(b"ACGTACGT");
        let b = IntToStr::new(b"ACGTTACGT");

        let aln = a.needleman_wunsch(&b);

        assert_eq!(aln.matches, 8);
        assert_eq!(aln.mismatches, 0);
        assert_eq!(aln.gaps, 1);
        assert_eq!(aln.aligned_self, "ACG-TACGT");
        assert_eq!(aln.aligned_other, "ACGTTACGT");
        assert!((aln.normalized_distance - (1.0 / 9.0)).abs() < 1e-6);
    }

    #[test]
    fn test_needleman_wunsch_detects_two_distributed_mismatches_without_gaps() {
        let a = IntToStr::new(b"ACGTCGTAAC");
        let b = IntToStr::new(b"ACGTTGTAAT");

        let aln = a.needleman_wunsch(&b);

        assert_eq!(aln.matches, 8);
        assert_eq!(aln.mismatches, 2);
        assert_eq!(aln.gaps, 0);
        assert_eq!(aln.aligned_self.len(), 10);
        assert_eq!(aln.aligned_other.len(), 10);
        assert!((aln.normalized_distance - 0.2).abs() < 1e-6);
    }

    #[test]
    fn test_best_exact_overlap_finds_realistic_shifted_primer_overlap() {
        let a = IntToStr::new(b"AAGCAGTGGTATCAACGC");
        let b = IntToStr::new(b"TGGTATCAACGCAGAGTAA");

        let hit = a.best_exact_overlap(&b, 10, 0.0).unwrap();

        assert_eq!(hit.offset, 6);
        assert_eq!(hit.overlap, 12);
        assert_eq!(hit.mismatches, 0);
        assert_eq!(hit.self_start, 6);
        assert_eq!(hit.other_start, 0);
    }

    #[test]
    fn test_best_exact_overlap_rejects_short_spurious_matches() {
        let a = IntToStr::new(b"AAGCAGTGGTATCAACGC");
        let b = IntToStr::new(b"TTTTTGGTAAAAACCCC");

        let hit = a.best_exact_overlap(&b, 10, 0.0);
        assert!(hit.is_none());
    }

    #[test]
    fn test_best_exact_overlap_respects_mismatch_fraction_cutoff() {
        let a = IntToStr::new(b"AAGCAGTGGTATCAACGC");
        let b = IntToStr::new(b"TGGTATCAATGCAGAGTAA");

        let strict = a.best_exact_overlap(&b, 12, 0.0);
        assert!(strict.is_none());

        let permissive = a.best_exact_overlap(&b, 12, 0.10).unwrap();
        assert_eq!(permissive.offset, 6);
        assert_eq!(permissive.overlap, 12);
        assert_eq!(permissive.mismatches, 1);
    }

    #[test]
    fn test_exact_overlap_with_negative_offset_is_detected() {
        let a = IntToStr::new(b"TGGTATCAACGCAGAGTAA");
        let b = IntToStr::new(b"AAGCAGTGGTATCAACGC");

        let hit = a.best_exact_overlap(&b, 10, 0.0).unwrap();

        assert_eq!(hit.offset, -6);
        assert_eq!(hit.overlap, 12);
        assert_eq!(hit.mismatches, 0);
        assert_eq!(hit.self_start, 0);
        assert_eq!(hit.other_start, 6);
    }

    #[test]
    fn test_needleman_wunsch_is_symmetric_for_distance() {
        let a = IntToStr::new(b"TACATGCAACTCAGCAGC");
        let b = IntToStr::new(b"TACATGCAACTCAACAGC");

        let d1 = a.needleman_wunsch_distance(&b);
        let d2 = b.needleman_wunsch_distance(&a);

        assert!((d1 - d2).abs() < 1e-6);
        assert!(d1 > 0.0);
        assert!(d1 < 0.2);
    }
}
