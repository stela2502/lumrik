use crate::sequence::{reference_base_matches, reverse_complement};

/// Coordinates for the best Smith-Waterman local alignment.
/// Query coordinates refer to the sequence passed to `local_alignment`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalAlignment {
    pub score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub reference_start: usize,
    pub reference_end: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrientedLocalAlignment {
    pub alignment: LocalAlignment,
    /// True when query coordinates refer to the reverse complement of the input query.
    pub reverse_complement: bool,
}

#[derive(Debug, Clone, Copy, Default)]
struct Cell {
    score: i32,
    query_start: usize,
    reference_start: usize,
}

/// Small Smith-Waterman implementation used only after seed pre-selection.
/// Match=+2, mismatch=-2, gap=-3. In addition to score it keeps the query and
/// reference interval of the winning local alignment so V/J geometry can be
/// validated rather than inferred from two independent scores.
pub fn local_alignment(query: &[u8], reference: &[u8]) -> LocalAlignment {
    if query.is_empty() || reference.is_empty() {
        return LocalAlignment {
            score: 0,
            query_start: 0,
            query_end: 0,
            reference_start: 0,
            reference_end: 0,
        };
    }

    let mut prev = vec![Cell::default(); reference.len() + 1];
    let mut curr = vec![Cell::default(); reference.len() + 1];
    let mut best = LocalAlignment {
        score: 0,
        query_start: 0,
        query_end: 0,
        reference_start: 0,
        reference_end: 0,
    };

    for (i, &q) in query.iter().enumerate() {
        curr[0] = Cell::default();
        for (j, &r) in reference.iter().enumerate() {
            let match_score = if reference_base_matches(r, q) { 2 } else { -2 };
            let diag_score = prev[j].score + match_score;
            let up_score = prev[j + 1].score - 3;
            let left_score = curr[j].score - 3;

            let mut cell = Cell::default();
            if diag_score > 0 && diag_score >= up_score && diag_score >= left_score {
                cell.score = diag_score;
                if prev[j].score > 0 {
                    cell.query_start = prev[j].query_start;
                    cell.reference_start = prev[j].reference_start;
                } else {
                    cell.query_start = i;
                    cell.reference_start = j;
                }
            } else if up_score > 0 && up_score >= left_score {
                cell.score = up_score;
                cell.query_start = prev[j + 1].query_start;
                cell.reference_start = prev[j + 1].reference_start;
            } else if left_score > 0 {
                cell.score = left_score;
                cell.query_start = curr[j].query_start;
                cell.reference_start = curr[j].reference_start;
            }
            curr[j + 1] = cell;

            let candidate = LocalAlignment {
                score: cell.score,
                query_start: cell.query_start,
                query_end: i + 1,
                reference_start: cell.reference_start,
                reference_end: j + 1,
            };
            if candidate.score > best.score
                || (candidate.score == best.score
                    && candidate.score > 0
                    && (candidate.query_start, candidate.query_end, candidate.reference_start, candidate.reference_end)
                        < (best.query_start, best.query_end, best.reference_start, best.reference_end))
            {
                best = candidate;
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        curr.fill(Cell::default());
    }
    best
}

/// Best local alignment in either read orientation. Forward wins exact ties.
pub fn best_oriented_local_alignment(query: &[u8], reference: &[u8]) -> OrientedLocalAlignment {
    let reverse = reverse_complement(query);
    best_oriented_local_alignment_with_reverse(query, &reverse, reference)
}

/// Best local alignment when the caller already has the reverse complement.
/// This avoids rebuilding the same reverse-complemented read for every germline
/// candidate in the posterior hot path.
pub fn best_oriented_local_alignment_with_reverse(
    query: &[u8],
    reverse: &[u8],
    reference: &[u8],
) -> OrientedLocalAlignment {
    let forward = local_alignment(query, reference);
    let reverse_alignment = local_alignment(reverse, reference);
    if reverse_alignment.score > forward.score {
        OrientedLocalAlignment {
            alignment: reverse_alignment,
            reverse_complement: true,
        }
    } else {
        OrientedLocalAlignment {
            alignment: forward,
            reverse_complement: false,
        }
    }
}

pub fn local_score(query: &[u8], reference: &[u8]) -> i32 {
    local_alignment(query, reference).score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_beats_mutant() {
        assert!(local_score(b"AACCGGTT", b"AACCGGTT") > local_score(b"AACCGGTT", b"AATTTGTT"));
    }

    #[test]
    fn iupac_reference_matches_allowed_read_base() {
        assert_eq!(local_score(b"T", b"Y"), 2);
        assert_eq!(local_score(b"C", b"Y"), 2);
        assert_eq!(local_score(b"A", b"Y"), 0);
    }

    #[test]
    fn reports_query_geometry() {
        let hit = local_alignment(b"TTTAACCGGTTGGG", b"AACCGGTT");
        assert_eq!(hit.score, 16);
        assert_eq!((hit.query_start, hit.query_end), (3, 11));
        assert_eq!((hit.reference_start, hit.reference_end), (0, 8));
    }

    #[test]
    fn oriented_alignment_can_choose_reverse_complement() {
        let hit = best_oriented_local_alignment(b"AACCGGTT", b"AACCGGTT");
        assert!(!hit.reverse_complement);
        let hit = best_oriented_local_alignment(b"AAAACCCC", b"GGGGTTTT");
        assert!(hit.reverse_complement);
        assert_eq!(hit.alignment.score, 16);
    }
}
