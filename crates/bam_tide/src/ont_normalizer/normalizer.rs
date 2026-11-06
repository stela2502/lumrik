use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::ngs_normalizer::{
    CHUNK_SIZE, NgsNormalizerSupport, NormalizedMolecule,
};
use crate::ont_normalizer::cli::Cli;
use crate::{AdditionalFeatureSource, FeatureTagCounts};
use read_tag_table::ReadTagTable;

use fast_tag_mapper::{BuiltinTagSet, FastTagMapper};

use anyhow::{Context, Result};
use int_to_dna::IntToDna;
use rayon::prelude::*;
use run_stats::RunStats;
use rust_htslib::bam::{Read, Reader};
use sc_primer::{Orientation, PrimerDetector};
use scdata::{GeneUmiHash, ScdataInsertState};

use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::File;
use std::collections::HashSet;
use crate::illumina_normalizer::normalizer::DedupKey;

#[derive(Debug, Clone)]
pub struct OntNormalizerConfig {
    pub bam: PathBuf,
    /// Top N completely rejected reads per CHUNK_SIZE input records (0 disables).
    pub dump_failed_reads: usize,
    /// Enable full, expensive per-batch complexity diagnostics.
    pub ont_diagnostics: bool,
    pub out: PathBuf,
    pub read_tags: PathBuf,
    pub min_transcript_len: usize,

    /// Maximum number of raw FASTQ pairs to process from one input pair.
    pub max_reads: Option<usize>,

    pub primer: PrimerDetector,
    pub additional_features: Vec<AdditionalFeatureSource>,
    pub additional_feature_min_hits: u32,

    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

impl OntNormalizerConfig {
    fn process_read_into(
        &self,
        read: &FastqRecord,
        feature_tag_mapper: &FastTagMapper,
        out: &mut OntPartial,
    ) {
        out.stats.report("total_records");
        out.stats.report("reads_processed");

        // ONT-only: reject extreme low-complexity reads before cassette detection.
        if is_low_complexity(&read.seq) {
            out.stats.report("low_complexity");
            out.stats.report("low_complexity_filtered");
            return;
        }

        let matches = match self.primer.detect_all(&read.seq, &read.qual) {
            Ok(matches) => matches,
            Err(_) => {
                out.stats.report("zero_cassette");
                out.stats.report("no_primer_match");
                out.stats.report("no_cell_umi");
                return;
            }
        };

        for diagnosis in self.primer.cassette_failure_diagnostics(&read.seq, &matches) {
            out.stats.report(diagnosis);
        }

        for diagnosis in self.primer.crosscheck_rejected(&read.seq, &read.qual, &matches) {
            out.stats.report(diagnosis);
        }

        match matches.len() {
            0 => {
                out.stats.report("zero_cassette");
                out.stats.report("no_cell_umi");
                return;
            }
            1 => out.stats.report("one_cassette"),
            _ => out.stats.report("multi_cassette"),
        }

        for (match_index, primer_match) in matches.iter().enumerate() {
            let insert = match primer_match.get_insert(&read.seq, &read.qual) {
                Ok(x) => x,
                Err(_) => {
                    out.stats.report("bad_insert_slice");
                    continue;
                }
            };

            if insert.seq.len() < self.min_transcript_len {
                out.stats.report("short_insert");
                continue;
            }

            let cell = match primer_match.get_cell(&read.seq, &read.qual) {
                Ok(x) => x,
                Err(_) => {
                    out.stats.report("bad_cell_slice");
                    continue;
                }
            };

            let umi = match primer_match.get_umi(&read.seq, &read.qual) {
                Ok(x) => x,
                Err(_) => {
                    out.stats.report("bad_umi_slice");
                    continue;
                }
            };

            NgsNormalizerSupport::report_orientation(&mut out.stats, primer_match.orientation);

            let insert_id = NgsNormalizerSupport::normalized_molecule_id(&read.id, match_index);
            let mut insert_record = FastqRecord {
                id: insert_id,
                seq: insert.seq.to_vec(),
                qual: insert.qual.to_vec(),
            };

            if primer_match.orientation == Orientation::ReverseComplement {
                insert_record = insert_record.revcomp();
            }

            let cell_id = IntToDna::new(&cell.seq).into_u64();
            let umi_id = IntToDna::new(&umi.seq).into_u64();

            if let Some(id) = feature_tag_mapper.map_feature_id(&insert_record.seq, &mut out.stats)
            {
                let state =
                    out.feature_tag_table
                        .try_insert(&cell_id, GeneUmiHash(id, umi_id), 1.0);
                out.stats.report(state.as_str());
                if state == ScdataInsertState::Inserted {
                    out.stats.report("unique_feature");
                    out.stats.report("feature_tag_match");
                }

                continue;
            }

            // Reuse the exact Illumina hard-UMI identity, using this cassette's
            // oriented insert rather than the entire concatenated ONT read.
            let normalized_cell = primer_match.cell_seq.as_deref().unwrap_or(&cell.seq);
            let identity = match self.primer.grammar_for_match(primer_match)
                .molecule_identity_if_exact(
                    Some(normalized_cell), Some(&umi.seq), &read.seq, &insert_record.seq,
                ) {
                Ok(Some(identity)) => identity,
                Ok(None) => {
                    out.stats.report("ambiguous_molecule_identity");
                    continue;
                }
                Err(_) => {
                    out.stats.report("invalid_molecule_identity");
                    continue;
                }
            };
            out.candidates.push((DedupKey {
                cell_id: identity.cell_id,
                hard_umi: identity.molecule_id,
            }, NormalizedMolecule {
                fastq: insert_record,
                original_read_id: Some(read.id.clone()),
                orientation: primer_match.orientation,
                cell_seq: cell.seq.to_vec(),
                cell_qual: cell.qual.to_vec(),
                umi_seq: umi.seq.to_vec(),
                umi_qual: umi.qual.to_vec(),
            }, matches.len() > 1));
        }
    }
}

// Parallel workers collect candidates; deduplication happens during the
// sequential merge so identities are shared across threads and BAM chunks.
struct OntPartial {
    candidates: Vec<(DedupKey, NormalizedMolecule, bool)>,
    stats: RunStats,
    feature_tag_table: scdata::Scdata,
}

impl OntPartial {
    fn new() -> Self {
        Self {
            candidates: Vec::new(),
            stats: RunStats::new(None, 0.0, 0),
            feature_tag_table: NgsNormalizerSupport::new_feature_tag_table(),
        }
    }
}

pub struct OntNormalizer {
    config: OntNormalizerConfig,
    stats: RunStats,
    read_tags: ReadTagTable,
    feature_tag_counts: FeatureTagCounts,
    seen: HashSet<DedupKey>,
}

impl OntNormalizer {
    pub fn new(config: OntNormalizerConfig) -> Result<Self> {
        let feature_tag_counts = FeatureTagCounts::from_sources(
            &config.additional_features,
            config.additional_feature_min_hits,
        )?;

        Ok(Self {
            config,
            stats: NgsNormalizerSupport::new_stats(),
            read_tags: ReadTagTable::new(),
            feature_tag_counts,
            seen: HashSet::new(),
        })
    }

    pub fn take_feature_tag_counts(&mut self) -> FeatureTagCounts {
        std::mem::take(&mut self.feature_tag_counts)
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        cli.out.prepare()?;
        cli.read_tags.prepare()?;
        let mut additional_features = Vec::new();
        if let Some(species) = cli.feature_tags.species {
            additional_features.push(match species {
                BuiltinTagSet::Human => AdditionalFeatureSource::BdSampleHuman,
                BuiltinTagSet::Mouse => AdditionalFeatureSource::BdSampleMouse,
            });
        }
        if let Some(path) = cli.feature_tags.tags {
            additional_features.push(AdditionalFeatureSource::Fasta(path.into_path_buf()));
        }

        Self::new(OntNormalizerConfig {
            bam: cli.bam.into_path_buf(),
            dump_failed_reads: 0,
            ont_diagnostics: false,
            out: cli.out.into_path_buf(),
            read_tags: cli.read_tags.into_path_buf(),
            min_transcript_len: cli.min_transcript_len,
            max_reads: cli.max_reads,
            primer: cli
                .primer
                .detector()
                .map_err(anyhow::Error::msg)?
                .with_read_wide_primer_search(true),
            additional_features,
            additional_feature_min_hits: cli.feature_tags.min_hits,
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        })
    }

    /// Normalizes input reads and streams mapper-ready FASTQ records to `emit`.
    ///
    /// Reads are processed in chunks using the normalizer's existing parallel
    /// normalization logic. Each accepted molecule is annotated with its
    /// `ReadTagRecord` metadata in the FASTQ QNAME before being passed to `emit`.
    ///
    /// Each callback receives one completed chunk as a borrowed slice of
    /// `(optional_r1, r2)` records. `r1` is `Some` when normalization produced
    /// a usable paired read and `None` for single-read output. The batch is
    /// mapper-ready and can be submitted without another per-read callback.
    ///
    /// This function does not write normalized FASTQ files. Instead, each
    /// normalized molecule is handed directly to the callback, allowing callers
    /// such as Nelrune to stream the reads into an external mapper.
    ///
    /// Returning `Ok(true)` from `emit` continues processing. Returning
    /// `Ok(false)` stops processing further input without treating this as an
    /// error.
    ///
    /// Normalization statistics, read tags, and feature-tag data are updated in
    /// the same way as during the normal file-based normalization path.    
    pub fn nelrune_run<F, P>(&mut self, mut emit: F, mut report_progress: P) -> Result<()>
    where
        F: FnMut(&[(Option<FastqRecord>, FastqRecord)]) -> Result<bool>,
        P: FnMut(&RunStats),
    {
        let mut bam = Reader::from_path(&self.config.bam)
            .with_context(|| format!("failed to open BAM: {}", self.config.bam.display()))?;

        if self.config.threads > 1 {
            bam.set_threads(self.config.threads.saturating_sub(1).max(1))
                .context("failed to set BAM reader threads")?;
        }

        let mut chunk = Vec::<FastqRecord>::with_capacity(CHUNK_SIZE);

        let mut output = Vec::<(Option<FastqRecord>, FastqRecord)>::new();

        let mut records = bam.records();
        let mut processed_from_input = 0usize;
        let mut failed_writer = if self.config.dump_failed_reads > 0 {
            let path = self.config.out.parent().unwrap_or_else(|| std::path::Path::new(".")).join("failed_reads.tsv");
            let mut writer = BufWriter::new(File::create(&path).with_context(|| format!("creating {}", path.display()))?);
            writeln!(writer, "batch\tread_id\tread_length\treason\tlow_entropy_window_fraction\tmean_4mer_entropy\tsequence")?;
            Some(writer)
        } else { None };
        let mut batch_number = 0usize;
        // Explicit opt-in: avoid creating files or recomputing metrics overnight.
        let mut diagnostic_writers = if self.config.ont_diagnostics {
            let dir = self.config.out.parent().unwrap_or_else(|| std::path::Path::new("."));
            let quantile_path = dir.join("read_complexity_by_batch.tsv");
            let length_path = dir.join("read_length_by_batch.tsv");
            let mut quantiles = BufWriter::new(File::create(&quantile_path)
                .with_context(|| format!("creating {}", quantile_path.display()))?);
            let mut lengths = BufWriter::new(File::create(&length_path)
                .with_context(|| format!("creating {}", length_path.display()))?);
            writeln!(quantiles, "batch\tpopulation\tmetric\tn\tmean\tsd\tmin\tp01\tp05\tp10\tp25\tp50\tp75\tp90\tp95\tp99\tmax")?;
            writeln!(lengths, "batch\tlength_bin\tall_reads\tcell_id\tno_cell_id\tno_cassette\tlow_complexity")?;
            Some((quantiles, lengths))
        } else { None };

        while self
            .config
            .max_reads
            .map_or(true, |max| processed_from_input < max)
        {
            let Some(rec_result) = records.next() else {
                break;
            };

            let rec = rec_result.context("failed to read BAM record")?;

            processed_from_input += 1;
            chunk.push(FastqRecord::from_bam_record(&rec));

            if chunk.len() >= CHUNK_SIZE {
                output.clear();

                self.process_chunk(&chunk, &mut output)?;
                batch_number += 1;
                if let Some((quantiles, lengths)) = diagnostic_writers.as_mut() {
                    self.write_batch_diagnostics(&chunk, batch_number, quantiles, lengths)?;
                }
                if let Some(writer) = failed_writer.as_mut() {
                    self.dump_failed_batch(&chunk, batch_number, writer)?;
                }

                // Publish normalization progress before handing the batch to the
                // downstream mapper. Mapper submission/draining may block while STAR
                // warms up or applies backpressure; the health server must still show
                // that this chunk has already been processed.
                report_progress(&self.stats);

                NgsNormalizerSupport::prepare_emit_batch(&mut output, &mut self.read_tags)?;

                if !output.is_empty() && !emit(&output)? {
                    return Ok(());
                }

                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            output.clear();

            self.process_chunk(&chunk, &mut output)?;
            batch_number += 1;
            if let Some((quantiles, lengths)) = diagnostic_writers.as_mut() {
                    self.write_batch_diagnostics(&chunk, batch_number, quantiles, lengths)?;
                }
            if let Some(writer) = failed_writer.as_mut() {
                self.dump_failed_batch(&chunk, batch_number, writer)?;
            }

            report_progress(&self.stats);

            NgsNormalizerSupport::prepare_emit_batch(&mut output, &mut self.read_tags)?;

            if !output.is_empty() && !emit(&output)? {
                return Ok(());
            }
        }

        if let Some(mut writer) = failed_writer { writer.flush()?; }
        if let Some((mut quantiles, mut lengths)) = diagnostic_writers {
            quantiles.flush()?;
            lengths.flush()?;
        }
        Ok(())
    }

    /// Complete, unsampled per-batch diagnostics on raw input sequences.
    /// Populations overlap: all_reads includes every read; low_complexity is
    /// a separate subset, and cell_id is based on successful cassette slicing.
    fn write_batch_diagnostics(
        &self, reads: &[FastqRecord], batch: usize,
        quantiles: &mut impl Write, lengths: &mut impl Write,
    ) -> Result<()> {
        const NAMES: [&str; 5] = ["all_reads", "cell_id", "no_cell_id", "no_cassette", "low_complexity"];
        const METRICS: [&str; 5] = ["read_length", "mean_4mer_entropy", "low_entropy_fraction", "longest_homopolymer", "tandem_repeat_coverage"];
        const BINS: [&str; 7] = ["lt250", "250_499", "500_999", "1000_1999", "2000_4999", "5000_9999", "ge10000"];
        let mut values: Vec<Vec<Vec<f64>>> = (0..5).map(|_| (0..5).map(|_| Vec::new()).collect()).collect();
        let mut bin_counts = [[0usize; 5]; 7];
        for read in reads {
            let seq = &read.seq;
            let scores = complexity_scores(seq);
            let low_complexity = scores.is_some_and(|(fraction, _)| fraction > 0.6);
            let matches = self.config.primer.detect_all(seq, &read.qual).ok();
            let no_cassette = matches.as_ref().is_none_or(|m| m.is_empty());
            let has_cell_id = !low_complexity && matches.as_ref().is_some_and(|ms| ms.iter().any(|m| {
                m.get_insert(seq, &read.qual).is_ok_and(|insert| insert.seq.len() >= self.config.min_transcript_len)
                    && m.get_cell(seq, &read.qual).is_ok()
                    && m.get_umi(seq, &read.qual).is_ok()
            }));
            let populations = [true, has_cell_id, !has_cell_id && !low_complexity, no_cassette && !low_complexity, low_complexity];
            let len = seq.len();
            let bin = match len { 0..=249 => 0, 250..=499 => 1, 500..=999 => 2, 1000..=1999 => 3, 2000..=4999 => 4, 5000..=9999 => 5, _ => 6 };
            let homopolymer = longest_homopolymer(seq) as f64;
            let repeat_coverage = tandem_repeat_coverage(seq) as f64;
            let metrics = [Some(len as f64), scores.map(|(_, h)| h), scores.map(|(f, _)| f), Some(homopolymer), Some(repeat_coverage)];
            for (pop, included) in populations.into_iter().enumerate() {
                if !included { continue; }
                bin_counts[bin][pop] += 1;
                for (metric, value) in metrics.iter().enumerate() {
                    if let Some(value) = value { values[pop][metric].push(*value); }
                }
            }
        }
        println!("[nelrune][batch {batch}] Complexity quantiles (all reads, no sampling):");
        println!("population\tmetric\tn\tmean\tsd\tmin\tp01\tp05\tp10\tp25\tp50\tp75\tp90\tp95\tp99\tmax");
        for (pop, name) in NAMES.iter().enumerate() {
            for (metric, metric_name) in METRICS.iter().enumerate() {
                let data = &mut values[pop][metric];
                if data.is_empty() { continue; }
                data.sort_unstable_by(f64::total_cmp);
                let n = data.len();
                let mean = data.iter().sum::<f64>() / n as f64;
                let sd = if n > 1 { (data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64).sqrt() } else { 0.0 };
                let mut row = format!("{batch}\t{name}\t{metric_name}\t{n}\t{mean:.6}\t{sd:.6}");
                for q in [0.0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0] {
                    let index = q * (n - 1) as f64;
                    let lo = index.floor() as usize;
                    let hi = index.ceil() as usize;
                    let value = data[lo] + (data[hi] - data[lo]) * (index - lo as f64);
                    row.push_str(&format!("\t{value:.6}"));
                }
                writeln!(quantiles, "{row}")?;
                println!("{}", row.split_once('\t').map_or(row.as_str(), |(_, rest)| rest));
            }
        }
        println!("[nelrune][batch {batch}] Read length bins:");
        println!("length_bin\tall_reads\tcell_id\tno_cell_id\tno_cassette\tlow_complexity");
        for (bin, name) in BINS.iter().enumerate() {
            let mut row = format!("{batch}\t{name}");
            for count in bin_counts[bin] { row.push_str(&format!("\t{count}")); }
            writeln!(lengths, "{row}")?;
            println!("{}", row.split_once('\t').map_or(row.as_str(), |(_, rest)| rest));
        }
        quantiles.flush()?;
        lengths.flush()?;
        Ok(())
    }

    /// Export the highest-entropy failed reads per batch, restricted to
    /// typical transcript lengths and without substantial repetitive sequence.
    fn dump_failed_batch(&self, reads: &[FastqRecord], batch: usize, writer: &mut impl Write) -> Result<()> {
        let limit = self.config.dump_failed_reads;
        if limit == 0 { return Ok(()); }
        let mut best: Vec<(&FastqRecord, &'static str, f64, f64)> = Vec::with_capacity(limit.min(reads.len()));
        for read in reads {
            if !(500..=2000).contains(&read.len()) { continue; }
            let Some((low_fraction, entropy)) = complexity_scores(&read.seq) else { continue; };
            if !entropy.is_finite() || low_fraction > 0.05 || tandem_repeat_coverage(&read.seq) > 0.05 { continue; }
            let matches = match self.config.primer.detect_all(&read.seq, &read.qual) {
                Ok(m) => m,
                Err(_) => continue,
            };
            // Focus on the cassette detection failure, not barcode/UMI slicing.
            if !matches.is_empty() { continue; }
            let reason = "no_cassette";
            if best.len() < limit {
                best.push((read, reason, low_fraction, entropy));
            } else if let Some((index, _)) = best.iter().enumerate().min_by(|(_, a), (_, b)| a.3.total_cmp(&b.3)) {
                if entropy > best[index].3 { best[index] = (read, reason, low_fraction, entropy); }
            }
        }
        best.sort_unstable_by(|a, b| b.3.total_cmp(&a.3).then_with(|| a.0.id.cmp(&b.0.id)));
        for (read, reason, low_fraction, entropy) in best {
            writeln!(writer, "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}", batch, read.id, read.len(), reason, low_fraction, entropy, String::from_utf8_lossy(&read.seq))?;
        }
        writer.flush()?;
        Ok(())
    }

    pub fn stats(&self) -> &RunStats {
        &self.stats
    }

    pub fn config(&self) -> &OntNormalizerConfig {
        &self.config
    }

    pub fn collect_fastqs(&mut self) -> Result<Vec<(Option<FastqRecord>, FastqRecord)>> {
        let mut bam = Reader::from_path(&self.config.bam)
            .with_context(|| format!("failed to open BAM: {}", self.config.bam.display()))?;

        if self.config.threads > 1 {
            bam.set_threads(self.config.threads.saturating_sub(1).max(1))
                .context("failed to set BAM reader threads")?;
        }

        let mut output: Vec<(Option<FastqRecord>, FastqRecord)> = Vec::new();

        let mut chunk: Vec<FastqRecord> = Vec::with_capacity(CHUNK_SIZE);

        let mut records = bam.records();
        let mut processed_from_input = 0usize;

        while self
            .config
            .max_reads
            .map_or(true, |max| processed_from_input < max)
        {
            let Some(rec_result) = records.next() else {
                break;
            };

            let rec = rec_result.context("failed to read BAM record")?;

            processed_from_input += 1;
            chunk.push(FastqRecord::from_bam_record(&rec));

            if chunk.len() >= CHUNK_SIZE {
                self.process_chunk(&chunk, &mut output)?;
                chunk.clear();

                eprintln!(
                    "processed {} reads; emitted {} molecules; sample-tag reads {}",
                    self.stats.get_issue_count("total_records"),
                    self.stats.get_issue_count("emitted_molecules"),
                    self.stats.get_issue_count("feature_tag_match"),
                );
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &mut output)?;
        }

        Ok(output)
    }

    pub fn run(&mut self) -> Result<()> {
        let fastqs = self.collect_fastqs()?;

        let mut fastq =
            FastqWriter::new(&self.config.out, self.config.gzip, self.config.gzip_level)
                .with_context(|| {
                    format!("failed to create FASTQ: {}", self.config.out.display())
                })?;

        for (_, read) in fastqs {
            fastq.write(&read)?;
            self.stats.report("fastq_reads_written");
        }

        fastq.finish()?;

        self.read_tags
            .save(&self.config.read_tags)
            .with_context(|| {
                format!(
                    "failed to write read-tag table {}",
                    self.config.read_tags.display()
                )
            })?;

        let out_root = self
            .config
            .out
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));

        let cell_barcode_len = self.config.primer.cell_len();

        self.feature_tag_counts
            .finalize_and_write_all(cell_barcode_len, out_root)?;

        Ok(())
    }

    fn process_chunk(
        &mut self,
        input: &[FastqRecord],
        output: &mut Vec<(Option<FastqRecord>, FastqRecord)>,
    ) -> Result<()> {
        let config = self.config.clone();

        self.stats.start_counter();
        self.stats
            .start_timer("bam_tide/multi_cpu/ont_normalize_chunk");

        let threads = rayon::current_num_threads().max(1);
        let chunk_size = (input.len() / threads).max(10_000);
        let partials: Vec<OntPartial> = input
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut out = OntPartial::new();
                for read in chunk {
                    config.process_read_into(read, self.feature_tag_counts.mapper(), &mut out);
                }
                out
            })
            .collect();

        self.stats
            .stop_timer("bam_tide/multi_cpu/ont_normalize_chunk");
        self.stats.stop_multi_processor_time();

        for partial in partials {
            self.stats.merge(&partial.stats);
            self.feature_tag_counts.merge_table(&partial.feature_tag_table);
            for (key, molecule, is_multi) in partial.candidates {
                if self.seen.insert(key) {
                    output.push((None, molecule.insert_read_tag(&mut self.read_tags)));
                    self.stats.report("emitted_molecules");
                    self.stats.report("unique_genomic");
                    self.stats.report(if is_multi {
                        "multi_cassette_emitted_molecules"
                    } else {
                        "single_cassette_emitted_molecules"
                    });
                } else {
                    self.stats.report("duplicate_molecules");
                    self.stats.report("duplicate");
                }
            }
        }

        Ok(())
    }

    pub fn stats_report(&self) -> String {
        let info = &self.stats;

        let total = info.get_issue_count("total_records");
        let zero = info.get_issue_count("zero_cassette");
        let one = info.get_issue_count("one_cassette");
        let multi = info.get_issue_count("multi_cassette");

        let emitted = info.get_issue_count("emitted_molecules");
        let written = info.get_issue_count("fastq_reads_written");

        let fwd = info.get_issue_count("forward_molecules");
        let rev = info.get_issue_count("reverse_molecules");

        let feature_tags = info.get_issue_count("feature_tag_match");
        let bad_insert = info.get_issue_count("bad_insert_slice");
        let bad_cell = info.get_issue_count("bad_cell_slice");
        let bad_umi = info.get_issue_count("bad_umi_slice");
        let too_short = info.get_issue_count("short_insert");
        let low_complexity = info.get_issue_count("low_complexity");
        let low_complexity_filtered = info.get_issue_count("low_complexity_filtered");

        let with_cassette = one + multi;

        let pct = |n: usize, d: usize| {
            if d == 0 {
                0.0
            } else {
                (n as f64 / d as f64) * 100.0
            }
        };

        let mean = if total == 0 {
            0.0
        } else {
            emitted as f64 / total as f64
        };

        format!(
            r#"bam-ont-normalizer summary
============================

Input
-----
reads processed        : {total}
low complexity         : {low_complexity} ({low_complexity_pct:.2}%)
low complexity filtered: {low_complexity_filtered}
reads w/o cassette     : {zero} ({zero_pct:.2}%)
reads w/ cassette      : {with_cassette} ({cassette_pct:.2}%)
  ├─ one cassette      : {one}
  └─ multi cassette    : {multi} ({multi_pct:.2}%)

Output
------
FASTQ reads written    : {written}
molecules emitted      : {emitted}
mean molecules/read    : {mean:.3}
feature-tag molecules  : {feature_tags}

Orientation
-----------
forward                : {fwd} ({fwd_pct:.2}%)
reverse                : {rev} ({rev_pct:.2}%)

Rejected
--------
bad insert slice       : {bad_insert}
bad cell slice         : {bad_cell}
bad UMI slice          : {bad_umi}
too short insert       : {too_short}
"#,
            total = total,
            low_complexity = low_complexity,
            low_complexity_pct = pct(low_complexity, total),
            low_complexity_filtered = low_complexity_filtered,
            zero = zero,
            zero_pct = pct(zero, total),
            with_cassette = with_cassette,
            cassette_pct = pct(with_cassette, total),
            one = one,
            multi = multi,
            multi_pct = pct(multi, total),
            written = written,
            emitted = emitted,
            mean = mean,
            feature_tags = feature_tags,
            fwd = fwd,
            rev = rev,
            fwd_pct = pct(fwd, emitted),
            rev_pct = pct(rev, emitted),
            bad_insert = bad_insert,
            bad_cell = bad_cell,
            bad_umi = bad_umi,
            too_short = too_short,
        )
    }
}

/// Reject reads whose bodies are dominated by low 4-mer Shannon entropy.
/// Uses IntToDna's packed 2-bit encoding (four bases per byte), with a
/// fixed 256-bin frequency table for each sampled 200-base window.
fn is_low_complexity(seq: &[u8]) -> bool {
    complexity_scores(seq).is_some_and(|(fraction, _)| fraction > 0.6)
}

/// Returns (fraction of windows below 4.0 bits, mean 4-mer entropy).
/// Uses precisely the same windows and entropy calculation as the filter.
fn complexity_scores(seq: &[u8]) -> Option<(f64, f64)> {
    const WINDOW: usize = 200;
    const STEP: usize = 100;
    const ENTROPY_LIMIT: f64 = 4.0;
    if seq.len() < WINDOW {
        return None;
    }

    // IntToDna encodes N as A. Avoid false low-entropy calls on ambiguous
    // windows by excluding them from the denominator altogether.
    let Ok(packed) = IntToDna::enc_bytes(seq) else {
        return None;
    };
    let mut low = 0usize;
    let mut total = 0usize;
    let mut entropy_sum = 0.0;
    const KMERS: usize = WINDOW - 3;
    // c * log2(c), precomputed once, no logarithms in the window loop.
    static C_LOG_C: std::sync::OnceLock<[f64; KMERS + 1]> = std::sync::OnceLock::new();
    let lookup = C_LOG_C.get_or_init(|| {
        let mut table = [0.0; KMERS + 1];
        for (c, value) in table.iter_mut().enumerate().skip(1) {
            *value = c as f64 * (c as f64).log2();
        }
        table
    });
    let max_weighted_sum = KMERS as f64 * ((KMERS as f64).log2() - ENTROPY_LIMIT);

    for start in (0..=seq.len() - WINDOW).step_by(STEP) {
        if seq[start..start + WINDOW]
            .iter()
            .any(|&base| !matches!(base, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
        {
            continue;
        }
        let mut counts = [0u8; 256];
        let mut kmer = 0u8;
        for pos in start..start + WINDOW {
            let base = (packed[pos / 4] >> ((pos % 4) * 2)) & 3;
            kmer = (kmer << 2) | base;
            if pos >= start + 3 {
                counts[kmer as usize] += 1;
            }
        }
        let weighted_sum: f64 = counts.iter().map(|&count| lookup[count as usize]).sum();
        total += 1;
        entropy_sum += (KMERS as f64).log2() - weighted_sum / KMERS as f64;
        if weighted_sum > max_weighted_sum {
            low += 1;
        }
    }
    if total == 0 {
        None
    } else {
        Some((low as f64 / total as f64, entropy_sum / total as f64))
    }
}

/// Longest identical-base run in the untrimmed read.
fn longest_homopolymer(seq: &[u8]) -> usize {
    let mut best = 0;
    let mut current = 0;
    let mut previous = 0;
    for &base in seq {
        current = if base == previous && matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') { current + 1 } else { 1 };
        previous = base;
        best = best.max(current);
    }
    best
}

/// Fraction of bases in runs of >= 24 bp with a repeated motif of 1..=6 bp.
/// A run is counted once even when multiple motif periods match.
fn tandem_repeat_coverage(seq: &[u8]) -> f64 {
    if seq.is_empty() { return 0.0; }
    let mut covered = vec![false; seq.len()];
    for period in 1..=6 {
        if seq.len() < period + 24 { continue; }
        let mut start = 0;
        while start + period < seq.len() {
            let mut end = start + period;
            while end < seq.len() && seq[end] == seq[end - period] && matches!(seq[end], b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
                end += 1;
            }
            if end - start >= 24 {
                covered[start..end].fill(true);
                start = end;
            } else {
                start += 1;
            }
        }
    }
    covered.iter().filter(|&&x| x).count() as f64 / seq.len() as f64
}

#[cfg(test)]
mod low_complexity_tests {
    use super::is_low_complexity;

    #[test]
    fn rejects_simple_repeats() {
        assert!(is_low_complexity(&b"AC".repeat(1000)));
        assert!(is_low_complexity(&b"A".repeat(2000)));
        assert!(is_low_complexity(&b"ATGGA".repeat(2000)));
    }

    #[test]
    fn keeps_short_and_diverse_sequences() {
        assert!(!is_low_complexity(&b"A".repeat(199)));
        assert!(!is_low_complexity(&b"N".repeat(2000)));
        let mut seq = Vec::new();
        for i in 0..2000u32 {
            let x = i.wrapping_mul(747796405).wrapping_add(2891336453);
            seq.push(b"ACGT"[((x ^ (x >> 13) ^ (x >> 19)) & 3) as usize]);
        }
        assert!(!is_low_complexity(&seq));
    }
}
