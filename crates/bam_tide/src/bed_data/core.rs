use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use bigtools::BigWigWrite;
use bigtools::beddata::BedParserStreamingIterator;
use rust_htslib::bam::{self, Read, Reader};

use rayon::prelude::*;

use crate::core::ref_block::record_to_blocks;
use crate::data_iter::DataIter;
use gtf_splice_index::types::RefBlock; // your new CIGAR-derived blocks

use crate::cli::CoverageCli;
use crate::core::alignment_policy::AlignmentPolicy;

use clap::ValueEnum;

/// Represents a single value in a bedGraph / bigWig stream
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl Value {
    pub fn flat(&self) -> (u32, u32, f32) {
        (self.start, self.end, self.value)
    }
}

/// Supported normalization options
#[derive(ValueEnum, Clone, Debug)]
pub enum Normalize {
    Not,
    Rpkm,
    Cpm,
    Bpm,
    Rpgc,
}

#[derive(Debug)]
pub struct BedData {
    pub genome_info: Vec<(String, usize, usize)>, // (chr, length, bin_offset)
    pub search: HashMap<String, usize>,           // chr -> genome_info index
    pub coverage_data: Vec<f32>,                  // bins
    pub bin_width: usize,
    pub threads: usize,
    pub nreads: usize,
}

impl fmt::Display for BedData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "BedData Report:")?;
        writeln!(f, "  Bin width: {}", self.bin_width)?;
        writeln!(f, "  Processed reads: {}", self.nreads)?;
        writeln!(f, "  Genome Info:")?;
        for (chr, len, offset) in &self.genome_info {
            writeln!(
                f,
                "    - Chr: {}, Length: {}, Bin offset: {}",
                chr, len, offset
            )?;
        }
        Ok(())
    }
}

impl BedData {
    // ============================
    // Initialization (NEW DESIGN)
    // ============================

    pub fn from_bam_with_policy(opts: &CoverageCli) -> Result<BedData, String> {
        let mut reader = Reader::from_path(&opts.bam)
            .map_err(|e| format!("bam file could not be read: {e:?}"))?;
        let header = reader.header().clone();

        let mut bed = BedData::init_from_header_view(&header, opts.width as usize, 1, None);

        let policy = AlignmentPolicy::from_cli(opts);

        for rec in reader.records() {
            let rec = rec.map_err(|e| format!("BAM read error: {e:?}"))?;
            if !policy.passes_filter(&rec) {
                continue;
            }

            let chr = std::str::from_utf8(header.tid2name(rec.tid() as u32))
                .map_err(|e| format!("Invalid chromosome name in BAM header: {e:?}"))?;
            let blocks = record_to_blocks(&rec);
            bed.add_ref_blocks(chr, &blocks);
        }

        bed.normalize(&opts.normalize);
        Ok(bed)
    }

    pub fn init_from_header_view(
        header: &bam::HeaderView,
        bin_width: usize,
        threads: usize,
        limit_to: Option<&[String]>,
    ) -> Self {
        let genome_info = Self::create_ref_id_to_name_vec(header, bin_width, limit_to);
        let search = Self::genome_info_to_search(&genome_info);

        let num_bins = genome_info
            .iter()
            .map(|(_, length, _)| length.div_ceil(bin_width))
            .sum::<usize>();

        let coverage_data = vec![0.0_f32; num_bins];

        Self {
            genome_info,
            search,
            coverage_data,
            bin_width,
            threads,
            nreads: 0,
        }
    }

    // ============================
    // Core accumulation
    // ============================
    pub fn add_ref_blocks(&mut self, chr: &str, blocks: &[RefBlock]) {
        let chr_id = match self.search.get(chr) {
            Some(id) => *id,
            None => return,
        };

        let (_chrom_name, chrom_length, chrom_offset) = &self.genome_info[chr_id];

        self.nreads += 1;

        let mut hit_bins: HashSet<usize> = HashSet::new();

        for block in blocks {
            let start = block.start as usize;
            let end = block.end.min(*chrom_length as u32) as usize;

            if end <= start {
                continue;
            }

            let start_window = start / self.bin_width;
            let end_window = (end - 1) / self.bin_width; // correct half-open handling

            for id in start_window..=end_window {
                let bin_start = id * self.bin_width;
                let bin_end = bin_start + self.bin_width;

                let overlap = (end.min(bin_end)).saturating_sub(start.max(bin_start));

                if overlap > 0 {
                    hit_bins.insert(id);
                }
            }
        }

        for id in hit_bins {
            let index = *chrom_offset + id;
            self.coverage_data[index] += 1.0;
        }
    }

    // ============================
    // Normalization
    // ============================

    pub fn normalize(&mut self, by: &Normalize) {
        // Step 0: convert from "covered bases per bin" to "mean depth per base in bin"
        // This matches deepTools default (normalizeUsing None).
        // HERE AFFECTS THE WAY WE COUNT
        //let bw = self.bin_width as f32;
        //self.coverage_data.par_iter_mut().for_each(|x| *x /= bw);

        match by {
            Normalize::Not => {
                // already mean depth
            }

            Normalize::Cpm => {
                // CPM: scale by mapped reads per million
                // mean_depth / (nreads / 1e6)
                let denom = (self.nreads as f32) / 1_000_000.0;
                if denom > 0.0 {
                    self.coverage_data.par_iter_mut().for_each(|x| *x /= denom);
                }
            }

            Normalize::Rpkm => {
                // RPKM for binned coverage:
                // RPKM = (reads_in_bin) / (reads_in_millions * bin_length_in_kb)
                //
                // We currently have mean depth (bases/bin_width).
                // For a fixed bin width, mean depth is proportional to reads_in_bin/bin_width.
                //
                // The simplest compatible definition used in coverage tools is:
                // mean_depth / (reads_in_millions) * (1 / bin_length_in_kb)
                //
                // Since bin_length_in_kb = bin_width / 1000:
                // divide by (reads_in_millions * (bin_width/1000)) == divide by reads_in_millions, then multiply by 1000/bin_width.
                //
                // But because we already divided by bin_width to get mean depth, we should NOT divide by bin_width again.
                // So: RPKM = mean_depth / reads_in_millions * 1000
                //
                // (This matches: (bases/bin_width) / (reads/1e6) * 1000
                //  = bases * 1e6 * 1000 / (reads * bin_width), which is the classic scaling.)
                let rpm = (self.nreads as f32) / 1_000_000.0;
                if rpm > 0.0 {
                    self.coverage_data
                        .par_iter_mut()
                        .for_each(|x| *x = (*x / rpm) * 1000.0);
                }
            }

            Normalize::Bpm => {
                // BPM: bins per million mapped reads (deepTools calls this BPM / RPGC variants differently)
                // Your previous implementation divided by sum(signal). That is more like "fraction of total signal".
                // If you want "per million bins" style, you'd typically scale by total signal per million.
                //
                // Keeping your original intent (normalize by total signal):
                let total: f32 = self.coverage_data.par_iter().sum();
                if total > 0.0 {
                    self.coverage_data.par_iter_mut().for_each(|x| *x /= total);
                }
            }

            Normalize::Rpgc => {
                panic!("Rpgc not implemented");
            }
        }
    }

    // ============================
    // Genome layout helpers
    // ============================

    pub fn genome_info_to_search(genome_info: &[(String, usize, usize)]) -> HashMap<String, usize> {
        genome_info
            .iter()
            .enumerate()
            .map(|(index, (name, _, _))| (name.clone(), index))
            .collect()
    }

    pub fn create_ref_id_to_name_vec(
        header: &bam::HeaderView,
        bin_width: usize,
        limit_to: Option<&[String]>,
    ) -> Vec<(String, usize, usize)> {
        let mut result = Vec::new();
        let mut total_bins = 0;

        for (rid, name) in header.target_names().iter().enumerate() {
            let chr = String::from_utf8_lossy(name).to_string();
            let len = header.target_len(rid as u32).unwrap() as usize;

            if let Some(allowed) = limit_to
                && !allowed.iter().any(|x| x == &chr)
            {
                continue;
            }

            let bins = len.div_ceil(bin_width);
            result.push((chr, len, total_bins));
            total_bins += bins;
        }

        result
    }

    pub fn id_for_chr_start(&self, chr: &str, start: usize) -> Option<usize> {
        self.search
            .get(chr)
            .map(|id| self.genome_info[*id].2 + start / self.bin_width)
    }

    pub fn current_chr_for_id(&self, id: usize) -> Option<(String, usize, usize)> {
        for (chr, length, offset) in &self.genome_info {
            if id >= *offset && (id - offset) * self.bin_width < *length {
                return Some((chr.clone(), *length, *offset));
            }
        }
        None
    }

    pub fn idx_to_ucsc_pos(&self, idx: usize) -> Option<(String, usize, usize)> {
        let (chr, chr_len, offset) = self.current_chr_for_id(idx)?;
        let rel = idx - offset;
        let start = rel * self.bin_width;
        if start >= chr_len {
            return None;
        }
        let end = (start + self.bin_width).min(chr_len);
        Some((chr, start, end))
    }

    // ============================
    // Writers (unchanged)
    // ============================

    pub fn write_bedgraph(&self, file_path: &str) -> std::io::Result<()> {
        let mut file = File::create(file_path)?;
        let iter = DataIter::new(self);

        for values in iter {
            writeln!(
                file,
                "{}\t{}\t{}\t{}",
                values.0, values.1.start, values.1.end, values.1.value
            )?;
        }
        Ok(())
    }

    pub fn write_bigwig(&self, file: &str) -> Result<(), String> {
        let outfile = Path::new(file);

        let chrom_map: HashMap<String, u32> = self
            .genome_info
            .iter()
            .map(|(chrom, len, _)| (chrom.clone(), *len as u32))
            .collect();

        let mut outb = BigWigWrite::create_file(outfile, chrom_map)
            .map_err(|e| format!("Failed to create BigWig file: {}", e))?;

        outb.options.channel_size = 0;
        outb.options.max_zooms = 1;
        outb.options.compress = true;
        outb.options.inmemory = false;

        let runtime = tokio::runtime::Builder::new_current_thread()
            .build()
            .unwrap();

        let iter = DataIter::new(self);
        let data = BedParserStreamingIterator::wrap_infallible_iter(iter, true);

        outb.write(data, runtime)
            .map_err(|e| format!("Failed to write BigWig file: {}", e))?;

        Ok(())
    }
}

#[cfg(test)]
mod binning_tests {
    use super::*;
    use gtf_splice_index::types::RefBlock; // adjust path to where RefBlock lives

    const EPS: f32 = 1e-6;

    /// Helper: make a minimal BedData with 1 chromosome and predictable offsets.
    fn bed_one_chr(chr: &str, chr_len: usize, bin_width: usize) -> BedData {
        let bins = (chr_len + bin_width - 1) / bin_width;

        let genome_info = vec![(chr.to_string(), chr_len, 0usize)];
        let search = BedData::genome_info_to_search(&genome_info);
        let coverage_data = vec![0.0_f32; bins];

        BedData {
            genome_info,
            search,
            coverage_data,
            bin_width,
            threads: 1,
            nreads: 0,
        }
    }

    fn get_bin(bed: &BedData, chr: &str, bin_id: usize) -> f32 {
        let idx = bed.search.get(chr).unwrap();
        let offset = bed.genome_info[*idx].2;
        bed.coverage_data[offset + bin_id]
    }

    #[test]
    fn test_add_ref_blocks_single_bin_full_overlap() {
        // chr length 100, bin=10 -> bins [0..10),[10..20),...
        let mut bed = bed_one_chr("chr1", 100, 10);

        // block exactly equals bin1: [10,20)
        bed.add_ref_blocks("chr1", &[RefBlock { start: 10, end: 20 }]);

        assert_eq!(bed.nreads, 1);
        assert!((get_bin(&bed, "chr1", 0) - 0.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 1) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 2) - 0.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_crosses_two_bins_partial() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        // [5,15) overlaps bin0 by 5 bp and bin1 by 5 bp
        bed.add_ref_blocks("chr1", &[RefBlock { start: 5, end: 15 }]);

        assert_eq!(bed.nreads, 1);
        assert!((get_bin(&bed, "chr1", 0) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 1) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_exact_bin_boundary_no_spill() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        // [0,10) should contribute only to bin0
        bed.add_ref_blocks("chr1", &[RefBlock { start: 0, end: 10 }]);

        assert_eq!(bed.nreads, 1);
        assert!((get_bin(&bed, "chr1", 0) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 1) - 0.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_multi_bins() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        // [0,25) => bin0:10, bin1:10, bin2:5
        bed.add_ref_blocks("chr1", &[RefBlock { start: 0, end: 25 }]);

        assert_eq!(bed.nreads, 1);
        assert!((get_bin(&bed, "chr1", 0) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 1) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 2) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 3) - 0.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_clips_to_chr_end() {
        let mut bed = bed_one_chr("chr1", 23, 10);
        // chr length 23: bins are [0..10),[10..20),[20..23)

        // block extends beyond chr end, should clip at 23: [15,30) -> [15,23)
        // overlaps bin1 [10..20): 5 bp (15..20)
        // overlaps bin2 [20..23): 3 bp (20..23)
        bed.add_ref_blocks("chr1", &[RefBlock { start: 15, end: 30 }]);

        assert_eq!(bed.nreads, 1);
        assert!((get_bin(&bed, "chr1", 1) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 2) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_ignores_empty_and_outside() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        // empty block should be ignored
        // block outside chr entirely should contribute nothing (after clipping)
        bed.add_ref_blocks(
            "chr1",
            &[
                RefBlock { start: 10, end: 10 },
                RefBlock {
                    start: 200,
                    end: 210,
                },
            ],
        );

        // Decide your intended behavior for nreads:
        // I recommend: if after clipping there is no contribution, still counts as a read OR not?
        //
        // Your spec said: "increments nreads once per read".
        // So even if the blocks are empty/outside, the call represents a read -> nreads++.
        assert_eq!(bed.nreads, 1);

        for b in 0..10 {
            assert!((get_bin(&bed, "chr1", b) - 0.0).abs() < EPS);
        }
    }

    #[test]
    fn test_add_ref_blocks_multiple_blocks_one_read_counts_once() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        // two blocks from same read
        bed.add_ref_blocks(
            "chr1",
            &[
                RefBlock { start: 0, end: 10 },  // bin0 += 10
                RefBlock { start: 20, end: 25 }, // bin2 += 5
            ],
        );

        assert_eq!(bed.nreads, 1);
        assert!((get_bin(&bed, "chr1", 0) - 1.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 2) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_accumulates_across_reads() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        bed.add_ref_blocks("chr1", &[RefBlock { start: 0, end: 10 }]); // +10
        bed.add_ref_blocks("chr1", &[RefBlock { start: 5, end: 15 }]); // +5 to bin0, +5 to bin1

        assert_eq!(bed.nreads, 2);
        assert!((get_bin(&bed, "chr1", 0) - 2.0).abs() < EPS);
        assert!((get_bin(&bed, "chr1", 1) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_add_ref_blocks_unknown_chr_is_error_or_noop() {
        let mut bed = bed_one_chr("chr1", 100, 10);

        // Decide intended behavior:
        // - either panic (strict)
        // - or ignore unknown chr (lenient)
        //
        // If your implementation returns Result, test that.
        // If it's a noop, test that nreads increments or not accordingly.

        bed.add_ref_blocks("chrX", &[RefBlock { start: 0, end: 10 }]);

        // Most robust: ignore unknown chr and DO NOT count read.
        // If you currently count it anyway, change expected here.
        assert_eq!(bed.nreads, 0);
        assert!((get_bin(&bed, "chr1", 0) - 0.0).abs() < EPS);
    }

    #[test]
    fn add_ref_blocks_counts_read_once_per_bin_even_with_two_blocks() {
        use std::collections::HashMap;

        // Minimal BedData with one chr length 200, bin width 50 => 4 bins
        let mut bed = BedData {
            genome_info: vec![("MT".to_string(), 200usize, 0usize)],
            search: {
                let mut h = HashMap::new();
                h.insert("MT".to_string(), 0usize);
                h
            },
            coverage_data: vec![0.0_f32; 4],
            bin_width: 50,
            threads: 1,
            nreads: 0,
        };

        // Two blocks, both inside bin 0 [0,50)
        let blocks = vec![RefBlock::new(10, 20), RefBlock::new(30, 40)];

        bed.add_ref_blocks("MT", &blocks);

        assert_eq!(bed.nreads, 1);
        assert_eq!(
            bed.coverage_data[0], 1.0,
            "bin 0 should be incremented once"
        );
        assert_eq!(bed.coverage_data[1], 0.0);
        assert_eq!(bed.coverage_data[2], 0.0);
        assert_eq!(bed.coverage_data[3], 0.0);
    }
}
