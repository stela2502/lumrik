use crate::reference::VdjReference;
use crate::types::Chain;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct SterileBin {
    pub bin: usize,
    pub start_fraction: f64,
    pub end_fraction: f64,
    pub unique_umis: usize,
    pub reads: usize,
}

#[derive(Debug, Clone)]
pub struct SupportedInterval {
    pub start: u32,
    pub end: u32,
    pub unique_umis: usize,
    pub reads: usize,
}

#[derive(Debug, Clone)]
pub struct SterileProfile {
    pub chain: Chain,
    pub bins: Vec<SterileBin>,
    pub supported_intervals: Vec<SupportedInterval>,
    pub breadth: f64,
    pub centroid: f64,
    pub proximal_fraction: f64,
    pub distal_fraction: f64,
    pub total_unique_umis: usize,
    pub total_reads: usize,
}

#[derive(Debug)]
pub struct SterileAccumulator {
    chain: Chain,
    chr: String,
    locus_start: u32,
    locus_end: u32,
    bins: usize,
    bin_umis: Vec<HashSet<String>>,
    bin_reads: Vec<usize>,
    intervals: HashMap<(u32, u32), (HashSet<String>, usize)>,
}

impl SterileAccumulator {
    pub fn new(reference: &VdjReference, chain: Chain, bins: usize) -> Option<Self> {
        let (chr, start, end) = reference.locus_bounds(chain)?;
        Some(Self {
            chain,
            chr: chr.to_string(),
            locus_start: start,
            locus_end: end,
            bins: bins.max(8),
            bin_umis: (0..bins.max(8)).map(|_| HashSet::new()).collect(),
            bin_reads: vec![0; bins.max(8)],
            intervals: HashMap::new(),
        })
    }

    pub fn observe(&mut self, chr: &str, start: u32, end: u32, umi: &str) {
        self.observe_n(chr, start, end, umi, 1);
    }

    pub fn observe_n(&mut self, chr: &str, start: u32, end: u32, umi: &str, reads: usize) {
        if reads == 0 {
            return;
        }
        if chr != self.chr
            || end <= self.locus_start
            || start >= self.locus_end
            || self.locus_end <= self.locus_start
        {
            return;
        }
        let s = start.max(self.locus_start);
        let e = end.min(self.locus_end);
        if e <= s {
            return;
        }
        let span = (self.locus_end - self.locus_start) as f64;
        let b0 = (((s - self.locus_start) as f64 / span) * self.bins as f64).floor() as usize;
        let b1 = ((((e - 1 - self.locus_start) as f64 / span) * self.bins as f64).floor() as usize)
            .min(self.bins - 1);
        for b in b0.min(self.bins - 1)..=b1 {
            self.bin_umis[b].insert(umi.to_string());
            self.bin_reads[b] += reads;
        }
        let entry = self
            .intervals
            .entry((s, e))
            .or_insert_with(|| (HashSet::new(), 0));
        entry.0.insert(umi.to_string());
        entry.1 += reads;
    }

    pub fn finish(self) -> SterileProfile {
        let mut bins = Vec::with_capacity(self.bins);
        let mut occupied = 0usize;
        let mut weighted = 0.0;
        let mut mass = 0.0;
        let mut all_umis = HashSet::new();
        let mut total_reads = 0;
        for i in 0..self.bins {
            let n = self.bin_umis[i].len();
            if n > 0 {
                occupied += 1
            }
            let mid = (i as f64 + 0.5) / self.bins as f64;
            weighted += mid * n as f64;
            mass += n as f64;
            all_umis.extend(self.bin_umis[i].iter().cloned());
            total_reads += self.bin_reads[i];
            bins.push(SterileBin {
                bin: i,
                start_fraction: i as f64 / self.bins as f64,
                end_fraction: (i + 1) as f64 / self.bins as f64,
                unique_umis: n,
                reads: self.bin_reads[i],
            });
        }
        let mut supported_intervals: Vec<_> = self
            .intervals
            .into_iter()
            .map(|((start, end), (u, r))| SupportedInterval {
                start,
                end,
                unique_umis: u.len(),
                reads: r,
            })
            .collect();
        supported_intervals.sort_by_key(|x| x.start);
        let proximal_mass: usize = bins
            .iter()
            .filter(|b| b.end_fraction <= 0.25)
            .map(|b| b.unique_umis)
            .sum();
        let distal_mass: usize = bins
            .iter()
            .filter(|b| b.start_fraction >= 0.75)
            .map(|b| b.unique_umis)
            .sum();
        SterileProfile {
            chain: self.chain,
            bins,
            supported_intervals,
            breadth: occupied as f64 / self.bins as f64,
            centroid: if mass > 0.0 { weighted / mass } else { 0.0 },
            proximal_fraction: if mass > 0.0 {
                proximal_mass as f64 / mass
            } else {
                0.0
            },
            distal_fraction: if mass > 0.0 {
                distal_mass as f64 / mass
            } else {
                0.0
            },
            total_unique_umis: all_umis.len(),
            total_reads,
        }
    }
}
