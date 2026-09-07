use anyhow::{bail, Result};

/// Expand constrained IUPAC ambiguity into concrete reference variants.
///
/// `N` remains unknown. Constrained symbols are expanded while the variant
/// count remains within `max_variants`. If another ambiguity would exceed the
/// cap, that position becomes N in all variants rather than inventing an allele.
pub fn expand_iupac(seq: &[u8], max_variants: usize) -> Result<Vec<Vec<u8>>> {
    if max_variants == 0 {
        bail!("max_variants must be greater than zero");
    }

    let mut variants = vec![Vec::with_capacity(seq.len())];

    for &raw in seq {
        let base = raw.to_ascii_uppercase();
        let choices: &[u8] = match base {
            b'A' => b"A",
            b'C' => b"C",
            b'G' => b"G",
            b'T' => b"T",
            b'N' => b"N",
            b'R' => b"AG",
            b'Y' => b"CT",
            b'S' => b"CG",
            b'W' => b"AT",
            b'K' => b"GT",
            b'M' => b"AC",
            b'B' => b"CGT",
            b'D' => b"AGT",
            b'H' => b"ACT",
            b'V' => b"ACG",
            other => bail!("unsupported FASTA base '{}'", other as char),
        };

        if choices.len() == 1 {
            for variant in &mut variants {
                variant.push(choices[0]);
            }
            continue;
        }

        let expanded_len = variants.len().saturating_mul(choices.len());
        if expanded_len > max_variants {
            for variant in &mut variants {
                variant.push(b'N');
            }
            continue;
        }

        let previous = std::mem::take(&mut variants);
        variants = Vec::with_capacity(expanded_len);
        for variant in previous {
            for &choice in choices {
                let mut next = variant.clone();
                next.push(choice);
                variants.push(next);
            }
        }
    }

    Ok(variants)
}

#[cfg(test)]
mod iupac_tests {
    use super::*;

    #[test]
    fn reference_normalization_preserves_iupac() {
        assert_eq!(
            normalize_reference_dna(b"acgtnryswkmbdhv").unwrap(),
            b"ACGTNRYSWKMBDHV"
        );
    }

    #[test]
    fn constrained_iupac_is_expanded() {
        let mut variants = expand_iupac(b"AYR", 8).unwrap();
        variants.sort();
        assert_eq!(
            variants,
            vec![
                b"ACA".to_vec(),
                b"ACG".to_vec(),
                b"ATA".to_vec(),
                b"ATG".to_vec(),
            ]
        );
    }

    #[test]
    fn expansion_cap_falls_back_to_unknown_not_fake_allele() {
        let variants = expand_iupac(b"YYY", 4).unwrap();
        assert_eq!(variants.len(), 4);
        assert!(variants.iter().all(|variant| variant[2] == b'N'));
    }
}

pub fn normalize_dna(seq: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(seq.len());
    for &base in seq {
        let b = base.to_ascii_uppercase();
        match b {
            b'A' | b'C' | b'G' | b'T' | b'N' => out.push(b),
            _ if b.is_ascii_whitespace() => {}
            _ => bail!("unsupported FASTA base {:?}", b as char),
        }
    }
    Ok(out)
}

/// Normalize a reference FASTA sequence while preserving standard IUPAC
/// ambiguity symbols.  Unlike read normalization, reference sequence may
/// legitimately contain constrained ambiguity such as Y=C/T or R=A/G.
pub fn normalize_reference_dna(seq: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(seq.len());
    for &base in seq {
        let b = base.to_ascii_uppercase();
        match b {
            b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B'
            | b'D' | b'H' | b'V' => out.push(b),
            _ if b.is_ascii_whitespace() => {}
            _ => bail!("unsupported FASTA base {:?}", b as char),
        }
    }
    Ok(out)
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            _ => b'N',
        })
        .collect()
}

pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    if seq.len() > 31 {
        return None;
    }

    let mut value = 0u64;
    for &base in seq {
        let bits = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        value = (value << 2) | bits;
    }
    Some(value)
}

/// Expand one reference k-mer containing constrained IUPAC ambiguity into
/// concrete A/C/G/T encodings.  This deliberately happens per k-mer, not per
/// biological V/D/J/C segment, so reference counts and locus geometry remain
/// unchanged.
///
/// N is genuinely unknown and therefore makes the k-mer unseedable.  If the
/// number of concrete combinations would exceed `max_variants`, the k-mer is
/// skipped rather than selecting an arbitrary allele.
pub fn encode_reference_kmers(seq: &[u8], max_variants: usize) -> Vec<u64> {
    if seq.len() > 31 || max_variants == 0 {
        return Vec::new();
    }

    let mut variants = vec![0u64];
    for &base in seq {
        let choices: &[u64] = match base.to_ascii_uppercase() {
            b'A' => &[0],
            b'C' => &[1],
            b'G' => &[2],
            b'T' => &[3],
            b'R' => &[0, 2],    // A/G
            b'Y' => &[1, 3],    // C/T
            b'S' => &[1, 2],    // C/G
            b'W' => &[0, 3],    // A/T
            b'K' => &[2, 3],    // G/T
            b'M' => &[0, 1],    // A/C
            b'B' => &[1, 2, 3], // C/G/T
            b'D' => &[0, 2, 3], // A/G/T
            b'H' => &[0, 1, 3], // A/C/T
            b'V' => &[0, 1, 2], // A/C/G
            b'N' => return Vec::new(),
            _ => return Vec::new(),
        };

        let needed = variants.len().saturating_mul(choices.len());
        if needed > max_variants {
            return Vec::new();
        }

        let old = std::mem::take(&mut variants);
        variants = Vec::with_capacity(needed);
        for prefix in old {
            for &bits in choices {
                variants.push((prefix << 2) | bits);
            }
        }
    }

    variants
}

/// True when the query base is one of the bases represented by the reference
/// IUPAC symbol. Query sequence itself is expected to be concrete A/C/G/T/N.
pub fn reference_base_matches(reference: u8, query: u8) -> bool {
    let q = query.to_ascii_uppercase();
    match reference.to_ascii_uppercase() {
        b'A' => q == b'A',
        b'C' => q == b'C',
        b'G' => q == b'G',
        b'T' => q == b'T',
        b'R' => matches!(q, b'A' | b'G'),
        b'Y' => matches!(q, b'C' | b'T'),
        b'S' => matches!(q, b'C' | b'G'),
        b'W' => matches!(q, b'A' | b'T'),
        b'K' => matches!(q, b'G' | b'T'),
        b'M' => matches!(q, b'A' | b'C'),
        b'B' => matches!(q, b'C' | b'G' | b'T'),
        b'D' => matches!(q, b'A' | b'G' | b'T'),
        b'H' => matches!(q, b'A' | b'C' | b'T'),
        b'V' => matches!(q, b'A' | b'C' | b'G'),
        b'N' => matches!(q, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_is_correct() {
        assert_eq!(reverse_complement(b"ACGTTN"), b"NAACGT");
    }

    #[test]
    fn reverse_complement_preserves_iupac_semantics() {
        assert_eq!(reverse_complement(b"ARYKBDHV"), b"BDHVMRYT");
    }

    #[test]
    fn kmer_rejects_n() {
        assert!(encode_kmer(b"ACNT").is_none());
        assert!(encode_kmer(b"ACGT").is_some());
    }

    #[test]
    fn reference_kmer_expands_y() {
        let mut got = encode_reference_kmers(b"AY", 8);
        got.sort_unstable();
        let mut expected = vec![encode_kmer(b"AC").unwrap(), encode_kmer(b"AT").unwrap()];
        expected.sort_unstable();
        assert_eq!(got, expected);
    }

    #[test]
    fn reference_kmer_cap_skips_excessive_ambiguity() {
        assert!(encode_reference_kmers(b"YYY", 4).is_empty());
        assert_eq!(encode_reference_kmers(b"YYY", 8).len(), 8);
    }

    #[test]
    fn iupac_reference_base_matches_allowed_bases() {
        assert!(reference_base_matches(b'Y', b'C'));
        assert!(reference_base_matches(b'Y', b'T'));
        assert!(!reference_base_matches(b'Y', b'A'));
        assert!(reference_base_matches(b'N', b'G'));
    }
}


#[derive(Debug, Clone)]
pub(crate) struct SequenceMerge {
    pub(crate) sequence: Vec<u8>,
    pub(crate) overlap: usize,
    pub(crate) mismatches: usize,
    pub(crate) extension: usize,
    /// Start of the right sequence relative to the left sequence.
    pub(crate) right_offset: isize,
    /// True when the reverse-complemented right sequence produced this merge.
    pub(crate) right_reverse: bool,
}

impl SequenceMerge {
    pub(crate) fn better_than(&self, other: &Self) -> bool {
        self.overlap > other.overlap
            || (self.overlap == other.overlap
                && (self.mismatches < other.mismatches
                    || (self.mismatches == other.mismatches
                        && (self.extension > other.extension
                            || (self.extension == other.extension
                                && self.sequence < other.sequence)))))
    }
}

pub(crate) fn merge_sequences_any_orientation(
    left: &[u8],
    right: &[u8],
    min_overlap: usize,
) -> Option<SequenceMerge> {
    let reverse = reverse_complement(right);
    let forward = merge_oriented_sequences(left, right, min_overlap, false);
    let reverse = merge_oriented_sequences(left, &reverse, min_overlap, true);
    match (forward, reverse) {
        (Some(a), Some(b)) => Some(if a.better_than(&b) { a } else { b }),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn merge_oriented_sequences(
    left: &[u8],
    right: &[u8],
    min_overlap: usize,
    right_reverse: bool,
) -> Option<SequenceMerge> {
    if left.is_empty() || right.is_empty() {
        return None;
    }
    let min_overlap = min_overlap.min(left.len()).min(right.len());
    if min_overlap == 0 {
        return None;
    }

    let min_offset = -(right.len() as isize) + min_overlap as isize;
    let max_offset = left.len() as isize - min_overlap as isize;
    let mut best: Option<SequenceMerge> = None;

    for offset in min_offset..=max_offset {
        let left_start = offset.max(0) as usize;
        let right_start = (-offset).max(0) as usize;
        let overlap = (left.len() - left_start).min(right.len() - right_start);
        if overlap < min_overlap {
            continue;
        }

        let mut mismatches = 0usize;
        for i in 0..overlap {
            let a = left[left_start + i].to_ascii_uppercase();
            let b = right[right_start + i].to_ascii_uppercase();
            if a != b && a != b'N' && b != b'N' {
                mismatches += 1;
            }
        }
        let max_mismatches = (overlap / 20).max(1);
        if mismatches > max_mismatches {
            continue;
        }

        let start = offset.min(0);
        let end = (left.len() as isize).max(offset + right.len() as isize);
        let mut sequence = Vec::with_capacity((end - start) as usize);
        for pos in start..end {
            let a = if pos >= 0 && pos < left.len() as isize {
                Some(left[pos as usize].to_ascii_uppercase())
            } else {
                None
            };
            let right_pos = pos - offset;
            let b = if right_pos >= 0 && right_pos < right.len() as isize {
                Some(right[right_pos as usize].to_ascii_uppercase())
            } else {
                None
            };
            sequence.push(match (a, b) {
                (Some(x), None) | (None, Some(x)) => x,
                (Some(x), Some(y)) if x == y => x,
                (Some(b'N'), Some(y)) => y,
                (Some(x), Some(b'N')) => x,
                (Some(_), Some(_)) => b'N',
                (None, None) => unreachable!(),
            });
        }

        let candidate = SequenceMerge {
            extension: sequence.len().saturating_sub(left.len().max(right.len())),
            sequence,
            overlap,
            mismatches,
            right_offset: offset,
            right_reverse,
        };
        if best.as_ref().map_or(true, |old| candidate.better_than(old)) {
            best = Some(candidate);
        }
    }
    best
}

/// Translate one DNA reading frame using the standard genetic code.
/// Incomplete trailing codons are ignored; non-ACGT codons become `X`.
pub fn translate_frame(sequence: &[u8], frame: usize) -> Vec<u8> {
    sequence
        .get(frame..)
        .unwrap_or_default()
        .chunks_exact(3)
        .map(translate_codon)
        .collect()
}

/// Translate a single DNA codon with the standard genetic code.
pub fn translate_codon(codon: &[u8]) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    let a = codon[0].to_ascii_uppercase();
    let b = codon[1].to_ascii_uppercase();
    let c = codon[2].to_ascii_uppercase();
    match (a, b, c) {
        (b'T', b'T', b'T' | b'C') => b'F',
        (b'T', b'T', b'A' | b'G') => b'L',
        (b'T', b'C', _) => b'S',
        (b'T', b'A', b'T' | b'C') => b'Y',
        (b'T', b'A', b'A' | b'G') => b'*',
        (b'T', b'G', b'T' | b'C') => b'C',
        (b'T', b'G', b'A') => b'*',
        (b'T', b'G', b'G') => b'W',
        (b'C', b'T', _) => b'L',
        (b'C', b'C', _) => b'P',
        (b'C', b'A', b'T' | b'C') => b'H',
        (b'C', b'A', b'A' | b'G') => b'Q',
        (b'C', b'G', _) => b'R',
        (b'A', b'T', b'T' | b'C' | b'A') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'A', b'C', _) => b'T',
        (b'A', b'A', b'T' | b'C') => b'N',
        (b'A', b'A', b'A' | b'G') => b'K',
        (b'A', b'G', b'T' | b'C') => b'S',
        (b'A', b'G', b'A' | b'G') => b'R',
        (b'G', b'T', _) => b'V',
        (b'G', b'C', _) => b'A',
        (b'G', b'A', b'T' | b'C') => b'D',
        (b'G', b'A', b'A' | b'G') => b'E',
        (b'G', b'G', _) => b'G',
        _ => b'X',
    }
}
