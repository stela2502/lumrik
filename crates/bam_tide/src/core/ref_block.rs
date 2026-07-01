//ref_block.rs

use rust_htslib::bam::Record;
use rust_htslib::bam::record::Cigar;

use gtf_splice_index::RefBlock;

/// Convert a BAM record into reference-covered blocks using its CIGAR string.
///
/// - Only reference-consuming *and* coverage-producing ops (M, =, X) create blocks.
/// - Ops that consume reference but do not represent coverage (D, N) advance the reference
///   without emitting blocks.
/// - Other ops (I, S, H, P) do not consume reference and do not emit blocks.
/// - Output blocks are 0-based, half-open and adjacent blocks are merged.
///
/// This function assumes the record is mapped; callers should filter unmapped reads beforehand.
pub fn record_to_blocks(rec: &Record) -> Vec<RefBlock> {
    let pos0 = rec.pos();
    if pos0 < 0 {
        return Vec::new();
    }

    let mut ref_pos: u64 = pos0 as u64;
    let mut blocks: Vec<RefBlock> = Vec::new();

    for op in rec.cigar().iter() {
        match *op {
            // Coverage on reference
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_u64 = len as u64;
                if len_u64 == 0 {
                    continue;
                }

                let start = ref_pos;
                let end = ref_pos + len_u64;

                if let Some(last) = blocks.last_mut() {
                    // Merge if overlapping or directly adjacent
                    if last.end as u64 >= start {
                        last.end = end.min(u32::MAX as u64) as u32;
                    } else {
                        blocks.push(RefBlock::new(
                            start.min(u32::MAX as u64) as u32,
                            end.min(u32::MAX as u64) as u32,
                        ));
                    }
                } else {
                    blocks.push(RefBlock::new(
                        start.min(u32::MAX as u64) as u32,
                        end.min(u32::MAX as u64) as u32,
                    ));
                }

                ref_pos = end;
            }

            // Deletion: advance reference but DO NOT break the block
            Cigar::Del(len) => {
                ref_pos += len as u64;
                // do not end the block
            }

            // Splice: break block
            Cigar::RefSkip(len) => {
                ref_pos += len as u64;
                // forces new block next time
            }

            // Does not consume reference
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    blocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn fake_record(pos: i64, cigar: CigarString) -> Record {
        let mut rec = Record::new();
        rec.set_pos(pos);
        rec.set_cigar(Some(&cigar));
        rec
    }

    #[test]
    fn test_simple_match() {
        let rec = fake_record(100, CigarString(vec![Cigar::Match(10)]));
        let blocks = record_to_blocks(&rec);
        assert_eq!(
            blocks,
            vec![RefBlock {
                start: 100,
                end: 110
            }]
        );
    }

    #[test]
    fn test_splice() {
        let rec = fake_record(
            100,
            CigarString(vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)]),
        );
        let blocks = record_to_blocks(&rec);
        assert_eq!(
            blocks,
            vec![
                RefBlock {
                    start: 100,
                    end: 105
                },
                RefBlock {
                    start: 115,
                    end: 120
                }
            ]
        );
    }

    /// Helper: build a BAM record with given position and CIGAR
    fn make_record(pos: i64, cigar: Vec<Cigar>) -> Record {
        let mut rec = Record::new();
        rec.set_pos(pos);
        rec.set_cigar(Some(&CigarString(cigar)));
        rec
    }

    #[test]
    fn test_deletion_splits_blocks() {
        // CIGAR: 56M1D35M
        let rec = make_record(
            108_996,
            vec![Cigar::Match(56), Cigar::Del(1), Cigar::Match(35)],
        );

        let blocks = record_to_blocks(&rec);

        // Deletion creates a gap in coverage
        assert_eq!(blocks.len(), 2);

        assert_eq!(blocks[0].start, 108_996);
        assert_eq!(blocks[0].end, 109_052);

        assert_eq!(blocks[1].start, 109_053);
        assert_eq!(blocks[1].end, 109_088);
    }
    #[test]
    fn test_split_blocks_across_splice() {
        // CIGAR: 50M1000N50M
        // Expected: two blocks
        let rec = make_record(
            1_000,
            vec![Cigar::Match(50), Cigar::RefSkip(1000), Cigar::Match(50)],
        );

        let blocks = record_to_blocks(&rec);

        assert_eq!(blocks.len(), 2, "Splice should split blocks");

        assert_eq!(blocks[0].start, 1_000);
        assert_eq!(blocks[0].end, 1_050);

        assert_eq!(blocks[1].start, 2_050);
        assert_eq!(blocks[1].end, 2_100);
    }
}
