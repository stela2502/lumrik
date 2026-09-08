use super::CellEvidence;
use crate::index::{Chain, VdjIndex};
use std::fmt;

/// Human-readable view of the compact receptor summaries retained after
/// batched BAM ingestion.
pub struct SummarizedEvidenceDisplay<'a> {
    evidence: &'a CellEvidence,
    index: &'a VdjIndex,
    chain: Chain,
}

/// Raw read-level evidence is intentionally not retained in the production
/// runner.  Keep this display type so the diagnostic CLI remains source/API
/// compatible and can explain the limitation instead of silently showing
/// fabricated data.
pub struct RawEvidenceDisplay<'a> {
    chain: Chain,
    _evidence: &'a CellEvidence,
}

impl CellEvidence {
    pub fn display_raw_chain<'a>(
        &'a self,
        _index: &'a VdjIndex,
        chain: Chain,
        _min_overlap: usize,
    ) -> RawEvidenceDisplay<'a> {
        RawEvidenceDisplay {
            chain,
            _evidence: self,
        }
    }

    pub fn display_summarized_chain<'a>(
        &'a self,
        index: &'a VdjIndex,
        chain: Chain,
        _min_overlap: usize,
    ) -> SummarizedEvidenceDisplay<'a> {
        SummarizedEvidenceDisplay {
            evidence: self,
            index,
            chain,
        }
    }
}

impl fmt::Display for SummarizedEvidenceDisplay<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let summaries = self.evidence.summaries_for_chain(self.chain);

        writeln!(f, "{} POST-MERGE SUMMARIES", self.chain)?;
        if summaries.is_empty() {
            writeln!(f, "  <no summaries>")?;
            return Ok(());
        }

        for (summary_no, summary) in summaries.iter().enumerate() {
            writeln!(f)?;
            writeln!(
                f,
                "summary {}  len={}  supporting_features={}",
                summary_no + 1,
                summary.len(),
                summary.support_features
            )?;

            let consensus = summary.consensus(self.index);
            writeln!(f, "consensus {}", String::from_utf8_lossy(&consensus))?;

            writeln!(f, "segment_support:")?;
            let mut support = summary.segment_support.clone();
            support
                .sort_by_key(|(id, _)| self.index.segment(*id).map(|s| (s.kind, s.name.clone())));
            for (id, count) in support {
                if let Some(segment) = self.index.segment(id) {
                    writeln!(
                        f,
                        "  {:?} {}  support={}",
                        segment.kind, segment.name, count
                    )?;
                } else {
                    writeln!(f, "  {:?}  support={}", id, count)?;
                }
            }

            writeln!(f, "germline_anchors:")?;
            if summary.germline_anchors.is_empty() {
                writeln!(f, "  <none>")?;
            } else {
                let mut anchors = summary.germline_anchors.clone();
                anchors.sort_by_key(|anchor| {
                    self.index
                        .segment(anchor.segment_id)
                        .map(|segment| (segment.kind, segment.name.clone(), anchor.summary_start))
                });
                for anchor in anchors {
                    if let Some(segment) = self.index.segment(anchor.segment_id) {
                        writeln!(
                            f,
                            "  {:?} {}  summary_start={}  germline_len={}",
                            segment.kind,
                            segment.name,
                            anchor.summary_start,
                            segment.sequence.len()
                        )?;
                    } else {
                        writeln!(
                            f,
                            "  {:?}  summary_start={}",
                            anchor.segment_id, anchor.summary_start
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}

impl fmt::Display for RawEvidenceDisplay<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{} PRE-MERGE RAW EVIDENCE", self.chain)?;
        writeln!(
            f,
            "  <not retained: production BAM ingestion compacts raw evidence in 20,000-record batches>"
        )
    }
}
