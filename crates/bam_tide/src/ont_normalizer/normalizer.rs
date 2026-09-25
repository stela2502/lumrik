use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::ngs_normalizer::{
    CHUNK_SIZE, NgsNormalizerSupport, NormalizedMolecule, NormalizerPartial,
};
use crate::ont_normalizer::cli::Cli;
use crate::{AdditionalFeatureSource, FeatureTagCounts};
use read_tag_table::ReadTagTable;

use fast_tag_mapper::{BuiltinTagSet, FastTagMapper};

use anyhow::{Context, Result};
use int_to_dna::IntToDna;
use mapping_info::MappingInfo;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Reader};
use sc_primer::{Orientation, PrimerDetector};
use scdata::{GeneUmiHash, ScdataInsertState};

use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct OntNormalizerConfig {
    pub bam: PathBuf,
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
        out: &mut NormalizerPartial,
    ) {
        out.stats.report("total_records");
        out.stats.report("reads_processed");

        let matches = match self.primer.detect_all(&read.seq, &read.qual) {
            Ok(matches) => matches,
            Err(_) => {
                out.stats.report("zero_cassette");
                out.stats.report("no_primer_match");
                out.stats.report("no_cell_umi");
                return;
            }
        };

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
                let state = out
                    .feature_tag_table
                    .try_insert(&cell_id, GeneUmiHash(id, umi_id), 1.0);
                out.stats.report(state.as_str());
                if state == ScdataInsertState::Inserted {
                    out.stats.report("unique_feature");
                    out.stats.report("feature_tag_match");
                }

                continue;
            }

            out.push_molecule(NormalizedMolecule {
                fastq: insert_record,
                original_read_id: Some(read.id.clone()),
                orientation: primer_match.orientation,
                cell_seq: cell.seq.to_vec(),
                cell_qual: cell.qual.to_vec(),
                umi_seq: umi.seq.to_vec(),
                umi_qual: umi.qual.to_vec(),
            });
            out.stats.report("emitted_molecules");
            out.stats.report("unique_genomic");
        }
    }
}

pub struct OntNormalizer {
    config: OntNormalizerConfig,
    stats: MappingInfo,
    read_tags: ReadTagTable,
    feature_tag_counts: FeatureTagCounts,
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
        })
    }

    pub fn take_feature_tag_counts(&mut self) -> FeatureTagCounts {
        std::mem::take(&mut self.feature_tag_counts)
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        let mut additional_features = Vec::new();
        if let Some(species) = cli.feature_tags.species {
            additional_features.push(match species {
                BuiltinTagSet::Human => AdditionalFeatureSource::BdSampleHuman,
                BuiltinTagSet::Mouse => AdditionalFeatureSource::BdSampleMouse,
            });
        }
        if let Some(path) = cli.feature_tags.tags {
            additional_features.push(AdditionalFeatureSource::Fasta(path));
        }

        Self::new(OntNormalizerConfig {
            bam: cli.bam,
            out: cli.out,
            read_tags: cli.read_tags,
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
        P: FnMut(&MappingInfo),
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

            report_progress(&self.stats);

            NgsNormalizerSupport::prepare_emit_batch(&mut output, &mut self.read_tags)?;

            if !output.is_empty() && !emit(&output)? {
                return Ok(());
            }
        }

        Ok(())
    }

    pub fn stats(&self) -> &MappingInfo {
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
        let partials: Vec<NormalizerPartial> = input
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut out = NormalizerPartial::new();
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
            let fastqs = partial.merge_into(
                &mut self.stats,
                &mut self.read_tags,
                self.feature_tag_counts.data_mut(),
            );

            output.extend(fastqs.into_iter().map(|record| (None, record)));
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
