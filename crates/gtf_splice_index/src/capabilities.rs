use crate::RefBlock;
use serde::{Deserialize, Serialize};

/// Biological coordinate systems understood by Ommverse-capable objects.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Axis {
    Genomic,
    Transcriptomic,
    Proteomic,
}

/// An object that retains every useful name/identifier by which it is known.
/// Fast string -> internal-ID resolution belongs to the owning index HashMap;
/// this trait exposes the aliases retained on the resolved object.
pub trait Identifiable {
    fn aliases(&self) -> &[String];
}

/// An object that has drawable interval geometry on one biological axis.
/// `RefBlock` remains the common, cheap 0-based half-open geometry type.
pub trait Plottable {
    fn axis(&self) -> Axis;
    fn blocks(&self) -> Vec<RefBlock>;
}

/// An object that knows how to translate interval geometry between biological
/// coordinate systems. Unsupported transformations return `None`; a supported
/// transformation with no overlap returns `Some(Vec::new())`.
pub trait CoordinateMapper {
    fn axes(&self) -> &[Axis];
    fn project(&self, from: Axis, to: Axis, blocks: &[RefBlock]) -> Option<Vec<RefBlock>>;
}

/// One typed edge exposed by a graph/tree-like biological object.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Connection<'a> {
    pub target: &'a str,
    pub relation: &'a str,
}

/// Capability for objects whose useful structure is relational rather than
/// spatial: STRING interactions, medical/disease graphs, InterPro hierarchy,
/// pathways, ontologies, and similar future providers.
pub trait Connected {
    fn connections(&self) -> Vec<Connection<'_>>;
}

pub(crate) fn compact_positions(mut positions: Vec<u32>) -> Vec<RefBlock> {
    positions.sort_unstable();
    positions.dedup();
    let Some(&first) = positions.first() else { return Vec::new(); };
    let mut out = Vec::new();
    let mut start = first;
    let mut prev = first;
    for &pos in positions.iter().skip(1) {
        if pos == prev + 1 {
            prev = pos;
        } else {
            out.push(RefBlock { start, end: prev + 1 });
            start = pos;
            prev = pos;
        }
    }
    out.push(RefBlock { start, end: prev + 1 });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Strand, Transcript};

    #[test]
    fn transcript_projects_protein_across_splice_junction() {
        let mut tx = Transcript::new(0, 0, "TX", 0, Strand::Plus);
        tx.add_exon(RefBlock::new(100, 106));
        tx.add_exon(RefBlock::new(200, 209));
        tx.add_cds(RefBlock::new(100, 106));
        tx.add_cds(RefBlock::new(200, 209));
        tx.finalize();

        let got = tx.project(Axis::Proteomic, Axis::Genomic, &[RefBlock::new(1, 4)]).unwrap();
        assert_eq!(got, vec![RefBlock::new(103, 106), RefBlock::new(200, 206)]);
    }

    #[test]
    fn transcript_advertises_all_three_axes() {
        let tx = Transcript::new(0, 0, "TX", 0, Strand::Plus);
        assert_eq!(tx.axes(), &[Axis::Genomic, Axis::Transcriptomic, Axis::Proteomic]);
        assert_eq!(tx.aliases(), &["TX".to_string()]);
    }
}
