//! Gene exonic-union length versus aligned read coverage, without per-read storage.
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use anyhow::{Context, Result};
use gtf_splice_index::SpliceIndex;

const THRESHOLDS: [u64; 6] = [10, 25, 50, 75, 90, 100];
const LABELS: [&str; 8] = ["<250 bp", "250-500 bp", "500-1000 bp", "1000-2000 bp", "2000-3000 bp", "3000-4000 bp", "4000-5000 bp", ">=5000 bp"];

#[derive(Clone, Copy)]
struct LengthBin { min: u64, max: u64 }
impl LengthBin {
    #[inline]
    fn fits(self, length: u64) -> bool { self.min <= length && length < self.max }
}
const LENGTH_BINS: [LengthBin; 8] = [
    LengthBin { min: 0, max: 250 }, LengthBin { min: 250, max: 500 },
    LengthBin { min: 500, max: 1000 }, LengthBin { min: 1000, max: 2000 },
    LengthBin { min: 2000, max: 3000 }, LengthBin { min: 3000, max: 4000 },
    LengthBin { min: 4000, max: 5000 }, LengthBin { min: 5000, max: u64::MAX },
];

#[derive(Clone, Default)]
struct LengthBinStats {
    reads: u64,
    coverage: [u64; 6],
}

/// Small per-worker counters; the gene lengths are stored only in SummaryTable.
#[derive(Clone)]
pub struct CoverageCounts {
    bins: [LengthBinStats; 8],
    observed_genes: Vec<bool>,
}
impl CoverageCounts {
    pub fn new(n_genes: usize) -> Self {
        Self { bins: std::array::from_fn(|_| LengthBinStats::default()), observed_genes: vec![false; n_genes] }
    }
    pub fn add(&mut self, gene_id: usize, gene_length: u64, aligned_bases: u64) {
        if gene_length == 0 { return; }
        let Some(bin_id) = LENGTH_BINS.iter().position(|b| b.fits(gene_length)) else { return; };
        let Some(seen) = self.observed_genes.get_mut(gene_id) else { return; };
        *seen = true;
        let bin = &mut self.bins[bin_id];
        bin.reads += 1;
        for (i, threshold) in THRESHOLDS.iter().enumerate() {
            if u128::from(aligned_bases) * 100 >= u128::from(gene_length) * u128::from(*threshold) {
                bin.coverage[i] += 1;
            }
        }
    }
    pub fn merge(&mut self, other: &Self) {
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            dst.reads += src.reads;
            for (a, b) in dst.coverage.iter_mut().zip(src.coverage) { *a += b; }
        }
        for (dst, src) in self.observed_genes.iter_mut().zip(&other.observed_genes) { *dst |= *src; }
    }
}

fn union_length(exons: &[gtf_splice_index::RefBlock]) -> u64 {
    let mut blocks = exons.to_vec();
    blocks.sort_unstable_by_key(|b| (b.start, b.end));
    let (mut total, mut end) = (0u64, 0u64);
    for exon in blocks {
        let start = u64::from(exon.start);
        let next_end = u64::from(exon.end);
        if next_end > end {
            total += next_end - start.max(end);
            end = next_end;
        }
    }
    total
}

/// Initialized once from the splice index; never stores per-read or per-UMI lengths.
pub struct SummaryTable {
    gene_lengths: Vec<u64>,
    counts: CoverageCounts,
}
impl SummaryTable {
    pub fn new(index: &SpliceIndex) -> Self {
        let gene_lengths = (0..index.genes.len()).map(|gene_id| {
            index.gene_exon_catalogue(gene_id).map(union_length).unwrap_or(0)
        }).collect::<Vec<_>>();
        let counts = CoverageCounts::new(gene_lengths.len());
        Self { gene_lengths, counts }
    }
    pub fn lengths(&self) -> &[u64] { &self.gene_lengths }
    pub fn empty_counts(&self) -> CoverageCounts { CoverageCounts::new(self.gene_lengths.len()) }
    pub fn merge(&mut self, counts: &CoverageCounts) { self.counts.merge(counts); }
    pub fn write_tsv(&self, path: &Path) -> Result<()> {
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let mut out = BufWriter::new(file);
        writeln!(out, "Statistic\t{}", LABELS.join("\t"))?;
        let lengths: [Vec<u64>; 8] = std::array::from_fn(|bin_id| {
            let mut values: Vec<u64> = self.gene_lengths.iter().enumerate()
                .filter(|(gene, len)| self.counts.observed_genes[*gene] && LENGTH_BINS[bin_id].fits(**len))
                .map(|(_, len)| *len).collect();
            values.sort_unstable(); values
        });
        write!(out, "Total reads")?;
        for bin in &self.counts.bins { write!(out, "\t{}", bin.reads)?; }
        writeln!(out)?;
        for label in ["Min gene exonic length", "Median gene exonic length", "Max gene exonic length"] {
            write!(out, "{label}")?;
            for values in &lengths {
                if values.is_empty() { write!(out, "\tNA")?; continue; }
                let value = match label {
                    "Min gene exonic length" => values[0] as f64,
                    "Max gene exonic length" => values[values.len()-1] as f64,
                    _ => { let m = values.len()/2; if values.len()%2 == 1 { values[m] as f64 } else { (values[m-1] as f64 + values[m] as f64)/2.0 } },
                };
                write!(out, "\t{value}")?;
            }
            writeln!(out)?;
        }
        for (i, threshold) in THRESHOLDS.iter().enumerate() {
            write!(out, "Reads covering >= {threshold}%")?;
            for bin in &self.counts.bins {
                if bin.reads == 0 { write!(out, "\tNA")?; }
                else { write!(out, "\t{:.2}", 100.0 * bin.coverage[i] as f64 / bin.reads as f64)?; }
            }
            writeln!(out)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bins_and_cumulative_thresholds() {
        assert!(LENGTH_BINS[0].fits(249));
        assert!(!LENGTH_BINS[0].fits(250));
        assert!(LENGTH_BINS[7].fits(5000));
        let mut counts = CoverageCounts::new(1);
        counts.add(0, 1000, 500);
        assert_eq!(counts.bins[3].coverage, [1, 1, 1, 0, 0, 0]);
    }
    #[test]
    fn union_merges_overlaps() {
        use gtf_splice_index::RefBlock;
        assert_eq!(union_length(&[RefBlock::new(100, 200), RefBlock::new(150, 250)]), 150);
    }
}
