use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result};
use gtf_splice_index::{GeneId, SpliceIndex};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{HeaderView, Record};
use scdata::FeatureIndex;

/// TE feature view over Lumrik's existing splice annotation index.
///
/// The splice index remains the genomic lookup structure. `TeIndex` only maps
/// `(chromosome, splice-index bin, TE gene)` to compact `u64` feature ids used
/// by Scdata. Feature ids are created lazily as overlaps are observed.
#[derive(Debug, Clone)]
pub struct TeIndex {
    splice: SpliceIndex,
    feature_ids: HashMap<(usize, usize, GeneId), u64>,
    feature_keys: Vec<(usize, usize, GeneId)>,
    feature_names: Vec<String>,
    feature_by_name: HashMap<String, u64>,
}

impl TeIndex {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let splice = SpliceIndex::load(path)
            .with_context(|| format!("failed to load TE splice index {}", path.display()))?;
        Ok(Self::from_splice_index(splice))
    }

    pub fn from_splice_index(splice: SpliceIndex) -> Self {
        Self {
            splice,
            feature_ids: HashMap::new(),
            feature_keys: Vec::new(),
            feature_names: Vec::new(),
            feature_by_name: HashMap::new(),
        }
    }

    pub fn splice_index(&self) -> &SpliceIndex {
        &self.splice
    }

    /// Number of spatial TE features actually observed so far.
    pub fn len(&self) -> usize {
        self.feature_keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.feature_keys.is_empty()
    }

    pub fn annotation_gene_count(&self) -> usize {
        self.splice.genes.len()
    }

    pub fn annotation_transcript_count(&self) -> usize {
        self.splice.transcripts.len()
    }

    pub fn feature_key(&self, feature_id: u64) -> Option<(usize, usize, GeneId)> {
        self.feature_keys.get(feature_id as usize).copied()
    }

    pub fn feature_coordinates(&self, feature_id: u64) -> Option<(&str, u64, u64, &str)> {
        let (chr_id, bin_id, gene_id) = self.feature_key(feature_id)?;
        let chrom = self.splice.chr_names.get(chr_id)?.as_str();
        let gene = self.splice.genes.get(gene_id)?.primary_name().unwrap_or("TE");
        let start = (bin_id as u64).saturating_mul(self.splice.bin_width as u64);
        let end = start.saturating_add(self.splice.bin_width as u64);
        Some((chrom, start, end, gene))
    }

    fn feature_for_key(&mut self, key: (usize, usize, GeneId)) -> u64 {
        if let Some(&feature_id) = self.feature_ids.get(&key) {
            return feature_id;
        }

        let (chr_id, bin_id, gene_id) = key;
        let feature_id = self.feature_keys.len() as u64;
        let chrom = self
            .splice
            .chr_names
            .get(chr_id)
            .map(String::as_str)
            .unwrap_or("NA");
        let gene = self
            .splice
            .genes
            .get(gene_id)
            .and_then(|g| g.primary_name())
            .unwrap_or("TE");
        let start = (bin_id as u64).saturating_mul(self.splice.bin_width as u64);
        let end = start.saturating_add(self.splice.bin_width as u64);
        let name = format!("{chrom}:{start}-{end}:{gene}");

        self.feature_ids.insert(key, feature_id);
        self.feature_keys.push(key);
        self.feature_by_name.insert(name.clone(), feature_id);
        self.feature_names.push(name);
        feature_id
    }

    /// Return spatial TE feature ids touched by aligned read bases.
    ///
    /// Deletions and skipped introns consume reference sequence but do not
    /// represent aligned read bases, so they do not contribute TE overlaps.
    pub fn record_overlaps(&mut self, rec: &Record, header: &HeaderView) -> Result<Vec<u64>> {
        if rec.is_unmapped() || rec.tid() < 0 {
            return Ok(Vec::new());
        }

        let chrom = std::str::from_utf8(header.tid2name(rec.tid() as u32))
            .context("non-UTF8 reference name in BAM header")?;
        let Some(chr_id) = self.splice.chr_id(chrom) else {
            return Ok(Vec::new());
        };

        let mut pos = rec.pos().max(0) as u32;
        let mut keys = HashSet::<(usize, usize, GeneId)>::new();

        for op in rec.cigar().iter() {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let end = pos.saturating_add(n);
                    for (bin_id, gene_id) in
                        self.splice.overlapping_genes_with_bins(chr_id, pos, end)
                    {
                        keys.insert((chr_id, bin_id, gene_id));
                    }
                    pos = end;
                }
                Cigar::Del(n) | Cigar::RefSkip(n) => pos = pos.saturating_add(n),
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }

        let mut keys: Vec<_> = keys.into_iter().collect();
        keys.sort_unstable();
        let mut ids = Vec::with_capacity(keys.len());
        for key in keys {
            ids.push(self.feature_for_key(key));
        }
        ids.sort_unstable();
        Ok(ids)
    }
}

impl FeatureIndex for TeIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        self.feature_names
            .get(feature_id as usize)
            .map(String::as_str)
            .unwrap_or("NA")
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.feature_by_name.get(name).copied()
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.feature_names.len() as u64).collect()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let feature = self.feature_name(feature_id);
        let gene = self
            .feature_key(feature_id)
            .and_then(|(_, _, gene_id)| self.splice.genes.get(gene_id))
            .and_then(|gene| gene.primary_name())
            .unwrap_or("TE");
        format!("{feature}\t{gene}\tTransposable Element")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtf_splice_index::AnnotationBuilder;
    use std::io::Cursor;

    const TEST_BIN_SIZE: u32 = 1_000_000;

    fn test_index(gtf: &[u8]) -> TeIndex {
        let splice = AnnotationBuilder::new(TEST_BIN_SIZE)
            .build_from_reader(Cursor::new(gtf))
            .unwrap();
        TeIndex::from_splice_index(splice)
    }

    #[test]
    fn creates_spatial_gene_features_lazily() {
        let gtf = b"chr1\tx\texon\t101\t200\t.\t+\t.\tgene_id \"L1HS\"; transcript_id \"L1HS_1\";\n";
        let index = test_index(gtf);

        assert_eq!(index.annotation_gene_count(), 1);
        assert_eq!(index.annotation_transcript_count(), 1);
        assert_eq!(index.len(), 0);
        assert_eq!(index.splice_index().transcripts[0].span(), Some((100, 200)));
    }

    #[test]
    fn reports_position_and_gene() {
        let gtf = b"chr1\tx\texon\t101\t200\t.\t+\t.\tgene_id \"L1HS\"; transcript_id \"L1HS_1\";\n";
        let mut index = test_index(gtf);
        let feature = index.feature_for_key((0, 0, 0));

        assert_eq!(feature, 0);
        assert_eq!(index.feature_name(0), "chr1:0-1000000:L1HS");
        assert_eq!(index.feature_id("chr1:0-1000000:L1HS"), Some(0));
        assert_eq!(
            index.to_10x_feature_line(0),
            "chr1:0-1000000:L1HS\tL1HS\tTransposable Element"
        );
    }
}
