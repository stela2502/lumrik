//! Direct fragment-level linkage between an assembled V(D)J call and a
//! transcript constant region.
//!
//! Raw BAM records are discarded in bounded batches.  During each batch we
//! collapse physical-fragment V/J -> C linkage into compact signature counters;
//! this module consumes those counters without retaining read-level evidence.

use super::{ConstantRegionEvidence, Recombination};
use crate::cellrep::CellEvidence;
use crate::index::{Chain, SegmentId, VdjIndex};
use std::collections::HashMap;

pub(super) fn rescue_missing_constant_regions(
    calls: &mut [Recombination],
    cell: &CellEvidence,
    index: &VdjIndex,
    chain: Chain,
) {
    for call in calls.iter_mut().filter(|call| call.constant.is_none()) {
        let Some((constant_id, supporting_fragments)) =
            best_direct_constant_link(call, cell, chain)
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

fn best_direct_constant_link(
    call: &Recombination,
    cell: &CellEvidence,
    chain: Chain,
) -> Option<(SegmentId, u32)> {
    let mut support = HashMap::<SegmentId, u32>::new();

    for (signature, &count) in cell.fragment_link_support() {
        if signature.chain != chain || count == 0 {
            continue;
        }
        // A constant call must be physically linked to the reconstructed J.
        // V-only -> C linkage is too permissive for RNA data because abundant
        // constant/sterile transcription can otherwise dominate assignment.
        if !signature.receptor_segments.contains(&call.j) {
            continue;
        }
        for &constant_id in &signature.constant_segments {
            let n = support.entry(constant_id).or_default();
            *n = n.saturating_add(count);
        }
    }

    support
        .into_iter()
        .max_by_key(|(segment_id, count)| (*count, std::cmp::Reverse(*segment_id)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cellrep::FragmentLinkSignature;

    #[test]
    fn fragment_link_signature_can_be_empty() {
        let signature = FragmentLinkSignature {
            chain: Chain::Igh,
            receptor_segments: Vec::new(),
            constant_segments: Vec::new(),
        };
        assert!(signature.receptor_segments.is_empty());
        assert!(signature.constant_segments.is_empty());
    }
}
