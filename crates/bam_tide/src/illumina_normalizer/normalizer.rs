use crate::fastq::{FastqPairReader, FastqRecord, FastqWriter};
use crate::illumina_normalizer::cli::{Cli, InsertRead, PrimerRead};
use crate::ngs_normalizer::{CHUNK_SIZE, NgsNormalizerSupport};
use crate::{AdditionalFeatureSource, FeatureTagCounts};
use read_tag_table::{ReadTagRecord, ReadTagTable};

use anyhow::{Context, Result, bail};
use int_to_str::IntToStr;
use mapping_info::MappingInfo;
use rayon::prelude::*;
use sc_primer::PrimerDetector;
use scdata::{GeneUmiHash, Scdata};

use fast_tag_mapper::{BuiltinTagSet, FastTagMapper};
use std::collections::HashSet;
use std::path::Path;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct IlluminaNormalizerConfig {
    pub out: PathBuf,
    pub read_tags: PathBuf,
    pub primer_read: PrimerRead,
    pub insert_read: InsertRead,
    pub primer: PrimerDetector,
    pub additional_features: Vec<AdditionalFeatureSource>,
    pub additional_feature_min_hits: u32,
    pub min_insert_len: usize,

    /// Maximum number of raw FASTQ pairs to process from one input pair.
    pub max_reads: Option<usize>,

    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

pub struct IlluminaPartial {
    pub candidates: Vec<IlluminaCandidate>,
    pub feature_tag_table: Scdata,
    pub stats: MappingInfo,
}

pub struct IlluminaCandidate {
    pub dedup_key: DedupKey,
    pub fastq_record: FastqRecord, // exported R2 / insert
    pub paired_r1_record: Option<FastqRecord>,
    pub read_tag: ReadTagRecord,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DedupKey {
    pub cell_id: u64,
    pub hard_umi: u64,
}

impl IlluminaPartial {
    pub fn new() -> Self {
        Self {
            candidates: Vec::new(),
            feature_tag_table: NgsNormalizerSupport::new_feature_tag_table(),
            stats: NgsNormalizerSupport::new_stats(),
        }
    }

    pub fn normalize_pair(
        &mut self,
        r1: &FastqRecord,
        r2: &FastqRecord,
        config: &IlluminaNormalizerConfig,
        feature_tag_mapper: &FastTagMapper,
    ) -> Result<()> {
        self.stats.report("total_pairs");
        self.stats.report("reads_processed");

        let primer_match = match config.primer.detect_first(&r1.seq, &r1.qual) {
            Ok(Some(x)) => x,
            Ok(None) => {
                self.stats.report("no_primer_match");
                self.stats.report("no_cell_umi");
                bail!("no primer match");
            }
            Err(err) => {
                self.stats.report("no_primer_match");
                self.stats.report("no_cell_umi");
                bail!("primer detection failed: {err}");
            }
        };

        let cell = primer_match.get_cell(&r1.seq, &r1.qual).map_err(|err| {
            self.stats.report("bad_cell_slice");
            self.stats.report("no_cell_umi");
            anyhow::anyhow!(err)
        })?;

        let umi = primer_match.get_umi(&r1.seq, &r1.qual).map_err(|err| {
            self.stats.report("bad_umi_slice");
            self.stats.report("no_cell_umi");
            anyhow::anyhow!(err)
        })?;

        let cell_id = IntToStr::new(&cell.seq).into_u64();
        let umi_id = IntToStr::new(&umi.seq).into_u64();

        if let Some(id) = feature_tag_mapper.map_feature_id(&r2.seq, &mut self.stats) {
            if self.feature_tag_table.try_insert(
                &cell_id,
                GeneUmiHash(id, umi_id),
                1.0,
                &mut self.stats,
            ) {
                self.stats.report("unique_feature");
                self.stats.report("feature_tag_match");
            } else {
                self.stats.report("duplicate");
                self.stats.report("duplicate_molecules");
            }
            return Ok(());
        }

        let mut emitted_r2 = r2.clone();
        emitted_r2.id = NgsNormalizerSupport::normalized_molecule_id(&r2.id, 0);

        let paired_r1_record = match primer_match.get_insert(&r1.seq, &r1.qual) {
            Ok(insert) if usable_insert(&insert.seq, 30, 0.5) => {
                self.stats.report("paired_r1_insert_found");
                Some(FastqRecord::new(&emitted_r2.id, &insert.seq, &insert.qual))
            }
            _ => {
                self.stats.report("no_usable_paired_r1_insert");
                None
            }
        };

        let cell_str = std::str::from_utf8(&cell.seq).unwrap();
        let cell_id = IntToStr::str_to_u64(cell_str)
            .expect("cell barcode should be <=32 bp ACGT after primer extraction");

        let mut hard_key = Vec::with_capacity(32);
        hard_key.extend_from_slice(&umi.seq);

        let remaining = 32usize.saturating_sub(umi.seq.len());
        hard_key.extend_from_slice(&r2.seq[..remaining.min(r2.seq.len())]);

        let hard_umi_str = std::str::from_utf8(&hard_key).unwrap();
        let hard_umi = IntToStr::str_to_u64(hard_umi_str)
            .expect("hard UMI should be <=32 bp ACGT after construction");

        let dedup_key = DedupKey { cell_id, hard_umi };

        let read_tag = ReadTagRecord::new(
            emitted_r2.id.clone(),
            Some(r2.id.clone()),
            &cell.seq,
            &cell.qual,
            &umi.seq,
            &umi.qual,
        );

        NgsNormalizerSupport::report_orientation(&mut self.stats, primer_match.orientation);

        self.candidates.push(IlluminaCandidate {
            dedup_key,
            fastq_record: emitted_r2,
            paired_r1_record,
            read_tag,
        });

        self.stats.report("candidate_pairs");

        Ok(())
    }
}

fn usable_insert(seq: &[u8], min_len: usize, max_single_base_fraction: f64) -> bool {
    if seq.len() < min_len {
        return false;
    }

    let mut counts = [0usize; 4];
    let mut acgt = 0usize;

    for b in seq {
        match b.to_ascii_uppercase() {
            b'A' => {
                counts[0] += 1;
                acgt += 1;
            }
            b'C' => {
                counts[1] += 1;
                acgt += 1;
            }
            b'G' => {
                counts[2] += 1;
                acgt += 1;
            }
            b'T' => {
                counts[3] += 1;
                acgt += 1;
            }
            _ => {}
        }
    }

    if acgt < min_len {
        return false;
    }

    let max_count = counts.into_iter().max().unwrap_or(0);

    (max_count as f64 / acgt as f64) <= max_single_base_fraction
}

pub struct IlluminaNormalizer {
    config: IlluminaNormalizerConfig,
    stats: MappingInfo,
    read_tags: ReadTagTable,
    feature_tag_counts: FeatureTagCounts,
    seen: HashSet<DedupKey>,
}

impl IlluminaNormalizer {
    pub fn new(config: IlluminaNormalizerConfig) -> Result<Self> {
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

        Self::new(IlluminaNormalizerConfig {
            out: cli.out,
            read_tags: cli.read_tags,
            primer_read: cli.primer_read,
            insert_read: cli.insert_read,
            primer: cli.primer.detector().map_err(anyhow::Error::msg)?,
            additional_features,
            additional_feature_min_hits: cli.feature_tags.min_hits,
            min_insert_len: cli.min_insert_len,
            max_reads: cli.max_reads,
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
    pub fn nelrune_run<F, P>(
        &mut self,
        r1_path: &PathBuf,
        r2_path: &PathBuf,
        mut emit: F,
        mut report_progress: P,
    ) -> Result<()>
    where
        F: FnMut(&[(Option<FastqRecord>, FastqRecord)]) -> Result<bool>,
        P: FnMut(&MappingInfo),
    {
        NgsNormalizerSupport::configure_rayon_threads(self.config.threads);

        let mut reader = FastqPairReader::from_paths(r1_path, r2_path).with_context(|| {
            format!(
                "failed to open FASTQ pair: {} and {}",
                r1_path.display(),
                r2_path.display(),
            )
        })?;

        let mut chunk = Vec::<(FastqRecord, FastqRecord)>::with_capacity(CHUNK_SIZE);

        let mut output = Vec::<(Option<FastqRecord>, FastqRecord)>::new();

        let mut processed_from_input = 0usize;

        while self
            .config
            .max_reads
            .map_or(true, |max| processed_from_input < max)
        {
            let Some(pair) = reader.next_pair()? else {
                break;
            };

            processed_from_input += 1;
            chunk.push(pair);

            if chunk.len() >= CHUNK_SIZE {
                output.clear();

                self.process_chunk(&chunk, &mut output)?;

                NgsNormalizerSupport::prepare_emit_batch(&mut output, &mut self.read_tags)?;

                if !output.is_empty() && !emit(&output)? {
                    return Ok(());
                }

                report_progress(&self.stats);
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            output.clear();

            self.process_chunk(&chunk, &mut output)?;

            NgsNormalizerSupport::prepare_emit_batch(&mut output, &mut self.read_tags)?;

            if !output.is_empty() && !emit(&output)? {
                return Ok(());
            }

            report_progress(&self.stats);
        }

        Ok(())
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &IlluminaNormalizerConfig {
        &self.config
    }

    /// Process the configured input FASTQs and return all accepted reads.
    ///
    /// The first element is the optional R1 insert read.
    /// The second element is always the emitted R2 read.
    ///
    /// This function:
    /// - reads the input FASTQs
    /// - performs primer/cell/UMI processing
    /// - performs feature-tag detection
    /// - performs deduplication
    /// - fills self.read_tags
    /// - fills self.feature_tag_table
    /// - updates self.stats
    ///
    /// It does NOT write output files.
    pub fn collect_fastqs(
        &mut self,
        r1_path: &Path,
        r2_path: &Path,
    ) -> Result<Vec<(Option<FastqRecord>, FastqRecord)>> {
        NgsNormalizerSupport::configure_rayon_threads(self.config.threads);

        let mut reader = FastqPairReader::from_paths(r1_path, r2_path).with_context(|| {
            format!(
                "failed to open FASTQ pair: {} and {}",
                r1_path.display(),
                r2_path.display()
            )
        })?;

        let mut output: Vec<(Option<FastqRecord>, FastqRecord)> = Vec::new();

        let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(CHUNK_SIZE);

        let mut processed_from_input = 0usize;

        while self
            .config
            .max_reads
            .map_or(true, |max| processed_from_input < max)
        {
            let Some(pair) = reader.next_pair()? else {
                break;
            };

            processed_from_input += 1;
            chunk.push(pair);

            if chunk.len() >= CHUNK_SIZE {
                self.process_chunk(&chunk, &mut output)?;

                chunk.clear();

                let total = self.stats.get_issue_count("total_pairs");

                let accepted = self.stats.get_issue_count("accepted_pairs");

                let pct = if total > 0 {
                    100.0 * accepted as f64 / total as f64
                } else {
                    0.0
                };

                eprintln!(
                    "processed {} read pairs; accepted {} ({:.2}%); failed {} pairs; {} duplicates detected; sample-tag reads {}",
                    total,
                    accepted,
                    pct,
                    self.stats.get_issue_count("failed_pairs"),
                    self.stats.get_issue_count("duplicate_molecules"),
                    self.stats.get_issue_count("feature_tag_match"),
                );
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &mut output)?;
        }

        Ok(output)
    }

    /// Original file-producing entry point.
    ///
    /// Collection/normalization happens in collect_fastqs().
    /// This function only adds the output-file side effects.
    pub fn run(&mut self, r1_path: &Path, r2_path: &Path) -> Result<()> {
        let fastqs = self.collect_fastqs(r1_path, r2_path)?;

        let (paired_r1_path, paired_r2_path) = self.paired_out_paths(&self.config.out);

        std::fs::create_dir_all(paired_r1_path.parent().unwrap_or_else(|| Path::new(".")))?;

        let mut fastq =
            FastqWriter::new(&self.config.out, self.config.gzip, self.config.gzip_level)
                .with_context(|| {
                    format!("failed to create FASTQ: {}", self.config.out.display())
                })?;

        let mut paired_r1 =
            FastqWriter::new(&paired_r1_path, self.config.gzip, self.config.gzip_level)
                .with_context(|| {
                    format!(
                        "failed to create paired R1 FASTQ: {}",
                        paired_r1_path.display()
                    )
                })?;

        let mut paired_r2 =
            FastqWriter::new(&paired_r2_path, self.config.gzip, self.config.gzip_level)
                .with_context(|| {
                    format!(
                        "failed to create paired R2 FASTQ: {}",
                        paired_r2_path.display()
                    )
                })?;

        for (r1, r2) in fastqs {
            if let Some(r1) = r1 {
                paired_r1.write(&r1)?;
                paired_r2.write(&r2)?;

                self.stats.report("paired_fastq_pairs_written");
            }

            fastq.write(&r2)?;
            self.stats.report("fastq_reads_written");
        }

        fastq.finish()?;
        paired_r1.finish()?;
        paired_r2.finish()?;

        self.read_tags
            .save(&self.config.read_tags)
            .with_context(|| {
                format!(
                    "failed to write read-tag table {}",
                    self.config.read_tags.display()
                )
            })?;

        let out_root = self.config.out.parent().unwrap_or_else(|| Path::new("."));

        let cell_barcode_len = self.config.primer.cell_len();

        self.feature_tag_counts
            .finalize_and_write_all(cell_barcode_len, out_root)?;

        Ok(())
    }

    fn paired_out_paths(&self, out: &Path) -> (PathBuf, PathBuf) {
        let parent = out.parent().unwrap_or_else(|| Path::new("."));
        let stem = self.fastq_stem(out);

        let paired_dir = parent.join("paired");

        (
            paired_dir.join(format!("{stem}.R1.fastq.gz")),
            paired_dir.join(format!("{stem}.R2.fastq.gz")),
        )
    }

    fn fastq_stem(&self, path: &Path) -> String {
        let name = path
            .file_name()
            .and_then(|x| x.to_str())
            .unwrap_or("normalized");

        name.strip_suffix(".fastq.gz")
            .or_else(|| name.strip_suffix(".fq.gz"))
            .or_else(|| name.strip_suffix(".fastq"))
            .or_else(|| name.strip_suffix(".fq"))
            .unwrap_or(name)
            .to_string()
    }

    fn process_chunk(
        &mut self,
        input: &[(FastqRecord, FastqRecord)],
        output: &mut Vec<(Option<FastqRecord>, FastqRecord)>,
    ) -> Result<()> {
        self.stats.start_counter();
        self.stats
            .start_timer("bam_tide/multi_cpu/illumina_normalize_chunk");

        let partials: Vec<IlluminaPartial> = input
            .par_iter()
            .map(|(r1, r2)| {
                let mut out = IlluminaPartial::new();

                if out
                    .normalize_pair(r1, r2, &self.config, self.feature_tag_counts.mapper())
                    .is_err()
                {
                    out.stats.report("failed_pairs");
                }

                out
            })
            .collect();

        self.stats
            .stop_timer("bam_tide/multi_cpu/illumina_normalize_chunk");
        self.stats.stop_multi_processor_time();

        for partial in partials {
            self.stats.merge(&partial.stats);
            self.feature_tag_counts
                .merge_table(&partial.feature_tag_table);

            for candidate in partial.candidates {
                if self.seen.insert(candidate.dedup_key) {
                    self.read_tags.insert(candidate.read_tag);

                    output.push((candidate.paired_r1_record, candidate.fastq_record));

                    self.stats.report("accepted_pairs");
                    self.stats.report("unique_genomic");
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

        let total = info.get_issue_count("total_pairs");
        let written = info.get_issue_count("fastq_reads_written");
        let failed = info.get_issue_count("failed_pairs");
        let accepted = info.get_issue_count("accepted_pairs");
        let no_primer = info.get_issue_count("no_primer_match");
        let short_insert = info.get_issue_count("short_insert");
        let bad_insert = info.get_issue_count("bad_insert_slice");
        let bad_cell = info.get_issue_count("bad_cell_slice");
        let bad_umi = info.get_issue_count("bad_umi_slice");
        let feature_tags = info.get_issue_count("feature_tag_match");
        let fwd = info.get_issue_count("forward_molecules");
        let rev = info.get_issue_count("reverse_molecules");

        let pct = |n: usize, d: usize| {
            if d == 0 {
                0.0
            } else {
                (n as f64 / d as f64) * 100.0
            }
        };

        format!(
            r#"bam-illumina-normalizer summary
=================================

Input
-----
FASTQ pairs processed : {total}
primer read           : {primer_read}
insert read           : {insert_read}

Output
------
accepted pairs        : {accepted} ({accepted_pct:.2}%)
FASTQ reads written   : {written} ({written_pct:.2}%)
feature-tag molecules : {feature_tags}
failed pairs          : {failed} ({failed_pct:.2}%)

Orientation
-----------
forward               : {fwd}
reverse               : {rev}

Rejected
--------
no primer match       : {no_primer}
short insert/read     : {short_insert}
bad insert slice      : {bad_insert}
bad cell slice        : {bad_cell}
bad UMI slice         : {bad_umi}
"#,
            total = total,
            primer_read = self.config.primer_read.as_str(),
            insert_read = self.config.insert_read.as_str(),
            accepted = accepted,
            accepted_pct = pct(accepted, total),
            written = written,
            written_pct = pct(written, total),
            feature_tags = feature_tags,
            failed = failed,
            failed_pct = pct(failed, total),
            fwd = fwd,
            rev = rev,
            no_primer = no_primer,
            short_insert = short_insert,
            bad_insert = bad_insert,
            bad_cell = bad_cell,
            bad_umi = bad_umi,
        )
    }
}
