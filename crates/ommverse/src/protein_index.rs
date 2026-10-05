use std::collections::HashMap;

use gtf_splice_index::placement::ChrBuckets;
use gtf_splice_index::{RefBlock, Strand, Transcript, TranscriptId};
use serde::{Deserialize, Serialize};

use crate::ProteinFeatureKind;

pub type ProteinFeatureId = usize;

/// One protein annotation after projection from amino-acid coordinates onto
/// the reference genome.  `blocks` contain coding genomic bases only; a
/// feature crossing splice junctions is therefore represented by multiple
/// blocks rather than by an intron-spanning interval.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProjectedProteinFeature {
    pub id: ProteinFeatureId,
    pub label: String,
    pub kind: ProteinFeatureKind,
    pub description: String,
    pub signature: String,
    pub chr_id: usize,
    pub strand: Strand,
    pub blocks: Vec<RefBlock>,

    /// External protein on which the source annotation was made (for example
    /// a UniProt accession from protein2ipr).  This is provenance, not
    /// structural ownership.
    pub source_protein: String,
    /// Transcript used to perform the protein -> genome projection.  Again,
    /// provenance only: other transcripts are matched from genomic geometry.
    pub source_transcript: TranscriptId,
    /// Original 0-based half-open amino-acid range.
    pub source_protein_range: (u32, u32),
}

impl ProjectedProteinFeature {
    pub fn span(&self) -> Option<(u32, u32)> {
        Some((self.blocks.first()?.start, self.blocks.last()?.end))
    }
}

/// Spatial index of protein features in genomic coordinates.
///
/// This deliberately mirrors only the lean placement part of
/// `gtf_splice_index`: InterPro/UniProt semantics stay in Ommverse, while the
/// chromosome/bin acceleration is shared with transcript placement.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinFeatureIndex {
    pub bin_width: u32,
    pub features: Vec<ProjectedProteinFeature>,
    pub chr_buckets: Vec<ChrBuckets>,
    span_start: Vec<u32>,
    span_end: Vec<u32>,
    by_label: HashMap<String, Vec<ProteinFeatureId>>,
}

impl Default for ProteinFeatureIndex {
    fn default() -> Self {
        Self::new(64, 0)
    }
}

impl ProteinFeatureIndex {
    pub fn new(bin_width: u32, chromosome_count: usize) -> Self {
        Self {
            bin_width,
            features: Vec::new(),
            chr_buckets: (0..chromosome_count).map(|_| ChrBuckets::new(bin_width)).collect(),
            span_start: Vec::new(),
            span_end: Vec::new(),
            by_label: HashMap::new(),
        }
    }

    /// Project one amino-acid feature through an annotated transcript onto the
    /// genome and add the resulting splice model to the index.
    pub fn add_from_transcript(
        &mut self,
        transcript: &Transcript,
        source_protein: impl Into<String>,
        label: impl Into<String>,
        kind: ProteinFeatureKind,
        description: impl Into<String>,
        signature: impl Into<String>,
        aa_start: u32,
        aa_end: u32,
    ) -> Option<ProteinFeatureId> {
        let blocks = transcript.protein_range_to_genomic_blocks(aa_start, aa_end);
        if blocks.is_empty() {
            return None;
        }
        let start = blocks.first()?.start;
        let end = blocks.last()?.end;
        if start >= end {
            return None;
        }

        let id = self.features.len();
        let label = label.into();
        let feature = ProjectedProteinFeature {
            id,
            label: label.clone(),
            kind,
            description: description.into(),
            signature: signature.into(),
            chr_id: transcript.chr_id,
            strand: transcript.strand,
            blocks,
            source_protein: source_protein.into(),
            source_transcript: transcript.id,
            source_protein_range: (aa_start, aa_end),
        };

        while self.chr_buckets.len() <= feature.chr_id {
            self.chr_buckets.push(ChrBuckets::new(self.bin_width));
        }
        self.chr_buckets[feature.chr_id].add_span(id, start, end);
        self.span_start.push(start);
        self.span_end.push(end);
        self.by_label.entry(label).or_default().push(id);
        self.features.push(feature);
        Some(id)
    }

    /// Collapse repeated observations of selected canonical feature labels at the
    /// same genomic locus.  A very large annotation is not allowed to absorb a
    /// much smaller one: candidates must overlap most of the shorter feature and
    /// have comparable genomic sizes.  Equal labels at distinct loci remain
    /// separate features.
    pub fn collapse_same_label_loci(&mut self, labels: &std::collections::HashSet<String>) {
        fn bases(blocks: &[RefBlock]) -> u64 {
            blocks.iter().map(|block| (block.end - block.start) as u64).sum()
        }
        fn overlap(a: &[RefBlock], b: &[RefBlock]) -> u64 {
            a.iter()
                .flat_map(|left| b.iter().map(move |right| (left, right)))
                .map(|(left, right)| {
                    left.end
                        .min(right.end)
                        .saturating_sub(left.start.max(right.start)) as u64
                })
                .sum()
        }
        fn same_locus(a: &ProjectedProteinFeature, b: &ProjectedProteinFeature) -> bool {
            if a.chr_id != b.chr_id || a.strand != b.strand || a.label != b.label {
                return false;
            }
            let a_len = bases(&a.blocks);
            let b_len = bases(&b.blocks);
            if a_len == 0 || b_len == 0 {
                return false;
            }
            let shared = overlap(&a.blocks, &b.blocks);
            let shorter = a_len.min(b_len);
            let longer = a_len.max(b_len);
            shared * 100 >= shorter * 80 && shorter * 2 >= longer
        }

        let mut order: Vec<usize> = (0..self.features.len()).collect();
        order.sort_by_key(|&id| {
            let feature = &self.features[id];
            (
                feature.chr_id,
                feature.label.clone(),
                bases(&feature.blocks),
                feature.span().map(|span| span.0).unwrap_or(u32::MAX),
            )
        });

        let mut kept = Vec::<ProjectedProteinFeature>::new();
        for id in order {
            let feature = &self.features[id];
            if labels.contains(&feature.label)
                && kept.iter().any(|existing| same_locus(existing, feature))
            {
                continue;
            }
            kept.push(feature.clone());
        }

        self.features = kept;
        self.chr_buckets = (0..self.chr_buckets.len())
            .map(|_| ChrBuckets::new(self.bin_width))
            .collect();
        self.span_start.clear();
        self.span_end.clear();
        self.by_label.clear();
        for (id, feature) in self.features.iter_mut().enumerate() {
            feature.id = id;
            let Some((start, end)) = feature.span() else { continue; };
            self.chr_buckets[feature.chr_id].add_span(id, start, end);
            self.span_start.push(start);
            self.span_end.push(end);
            self.by_label.entry(feature.label.clone()).or_default().push(id);
        }
    }

    /// Sort/deduplicate placement buckets after bulk construction.
    pub fn finalize(&mut self) {
        for chr in &mut self.chr_buckets {
            chr.finalize_by_start(&self.span_start, &self.span_end);
        }
    }

    pub fn ids_for_label(&self, label: &str) -> &[ProteinFeatureId] {
        self.by_label.get(label).map(Vec::as_slice).unwrap_or(&[])
    }

    /// Candidate features whose genomic spans overlap `[start, end)`.
    /// Exact block/splice compatibility is intentionally a separate operation.
    pub fn candidates_for_region(&self, chr_id: usize, start: u32, end: u32) -> Vec<ProteinFeatureId> {
        if start >= end {
            return Vec::new();
        }
        let Some(chr) = self.chr_buckets.get(chr_id) else {
            return Vec::new();
        };
        if chr.bins.is_empty() {
            return Vec::new();
        }
        let first_bin = (start / chr.bin_width) as usize;
        let last_bin = ((end - 1) / chr.bin_width) as usize;
        if first_bin >= chr.bins.len() {
            return Vec::new();
        }

        let mut ids = Vec::new();
        for bin_idx in first_bin..=last_bin.min(chr.bins.len() - 1) {
            for &id in &chr.bins[bin_idx] {
                if self.span_start[id] < end && self.span_end[id] > start {
                    ids.push(id);
                }
            }
        }
        ids.sort_unstable();
        ids.dedup();
        ids
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_feature_projects_to_spliced_genome_and_is_region_queryable() {
        let mut tx = Transcript::new(0, 0, "TX1", 0, Strand::Plus);
        tx.add_exon(RefBlock::new(100, 109));
        tx.add_exon(RefBlock::new(200, 212));
        tx.add_cds(RefBlock::new(100, 109));
        tx.add_cds(RefBlock::new(200, 212));
        tx.finalize();

        let mut idx = ProteinFeatureIndex::new(64, 1);
        let id = idx
            .add_from_transcript(&tx, "P1", "IPR_TEST", ProteinFeatureKind::Domain, "test domain", "SIG_TEST", 2, 5)
            .expect("feature should project");
        idx.finalize();

        assert_eq!(idx.features[id].blocks, vec![RefBlock::new(106, 109), RefBlock::new(200, 206)]);
        assert_eq!(idx.ids_for_label("IPR_TEST"), &[id]);
        assert_eq!(idx.candidates_for_region(0, 108, 202), vec![id]);
        assert!(idx.candidates_for_region(0, 300, 400).is_empty());
    }
}
