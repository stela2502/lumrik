//! Direct fragment-level linkage between an assembled V(D)J call and a
//! transcript constant region.
//!
//! Constant exons are spliced onto receptor transcripts, so they need not be
//! contiguous with the compact sequence used to reconstruct V(D)J.  This
//! module therefore assigns a missing constant region only when the same
//! physical fragment (`EvidenceId`) directly supports the called V or J and a
//! constant segment.  Paired mates and same-QNAME multimapping BAM records
//! share an `EvidenceId` during BAM ingestion.

use super::{ConstantRegionEvidence, Recombination};
use crate::cellrep::{CellEvidence, EvidenceId};
use crate::index::{Chain, SegmentId, SegmentKind, VdjIndex};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Default)]
struct FragmentSegmentLinks {
    receptor_segments: HashSet<SegmentId>,
    constant_segments: HashSet<SegmentId>,
}

pub(super) fn rescue_missing_constant_regions(
    calls: &mut [Recombination],
    cell: &CellEvidence,
    index: &VdjIndex,
    chain: Chain,
) {
    let fragments = collect_fragment_segment_links(cell, index, chain);

    for call in calls.iter_mut().filter(|call| call.constant.is_none()) {
        let Some((constant_id, supporting_fragments)) =
            best_direct_constant_link(call, &fragments)
        else {
            continue;
        };
        let Some(segment) = index.segment(constant_id) else {
            continue;
        };

        call.constant = Some(ConstantRegionEvidence {
            segment: constant_id,
            supporting_features: supporting_fragments,
            sequence: segment.sequence.clone(),
        });
    }
}

fn collect_fragment_segment_links(
    cell: &CellEvidence,
    index: &VdjIndex,
    chain: Chain,
) -> HashMap<EvidenceId, FragmentSegmentLinks> {
    let mut fragments = HashMap::<EvidenceId, FragmentSegmentLinks>::new();

    for feature in cell.features_for_chain(index, chain) {
        let links = fragments.entry(feature.id).or_default();
        for mapping in &feature.mappings {
            let Some(segment) = index.segment(mapping.segment_id) else {
                continue;
            };
            if segment.chain != chain {
                continue;
            }

            match segment.kind {
                SegmentKind::V | SegmentKind::J => {
                    links.receptor_segments.insert(mapping.segment_id);
                }
                SegmentKind::C => {
                    links.constant_segments.insert(mapping.segment_id);
                }
                SegmentKind::D => {}
            }
        }
    }

    fragments
}

fn best_direct_constant_link(
    call: &Recombination,
    fragments: &HashMap<EvidenceId, FragmentSegmentLinks>,
) -> Option<(SegmentId, u32)> {
    let mut support = HashMap::<SegmentId, u32>::new();

    for links in fragments.values() {
        let supports_call =
            links.receptor_segments.contains(&call.v) || links.receptor_segments.contains(&call.j);
        if !supports_call {
            continue;
        }

        for &constant_id in &links.constant_segments {
            *support.entry(constant_id).or_default() += 1;
        }
    }

    support
        .into_iter()
        .max_by_key(|(segment_id, count)| (*count, std::cmp::Reverse(*segment_id)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragment_links_default_empty() {
        let links = FragmentSegmentLinks::default();
        assert!(links.receptor_segments.is_empty());
        assert!(links.constant_segments.is_empty());
    }
}
