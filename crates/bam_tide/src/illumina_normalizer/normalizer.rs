use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::illumina_normalizer::cli::{Cli, InsertRead, PrimerRead};
use crate::ngs_normalizer::{
    FastqPairReader, NgsNormalizerSupport, NormalizedMolecule, NormalizerPartial, CHUNK_SIZE,
};
use crate::read_tag_table::{ReadTagTable,ReadTagRecord};

use anyhow::{bail, Context, Result};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use sc_primer::{Orientation, PrimerDetector};
use scdata::{Scdata, GeneUmiHash};
use int_to_str::IntToStr;

use std::path::Path;
use std::path::PathBuf;
use std::collections::HashSet;
use fast_tag_mapper::FastTagMapper;


#[derive(Debug, Clone)]
pub struct IlluminaNormalizerConfig {
    pub r1: PathBuf,
    pub r2: PathBuf,
    pub out: PathBuf,
    pub read_tags: PathBuf,
    pub primer_read: PrimerRead,
    pub insert_read: InsertRead,
    pub primer: PrimerDetector,
    pub feature_tag_mapper: Option<FastTagMapper>,
    pub min_insert_len: usize,
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
    pub fastq_record: FastqRecord,              // exported R2 / insert
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
    ) -> Result<()> {
        self.stats.report("total_pairs");

        let primer_match = match config.primer.detect_first(&r1.seq, &r1.qual) {
            Ok(Some(x)) => x,
            Ok(None) => {
                self.stats.report("no_primer_match");
                bail!("no primer match");
            }
            Err(err) => {
                self.stats.report("no_primer_match");
                bail!("primer detection failed: {err}");
            }
        };

        let cell = primer_match.get_cell(&r1.seq, &r1.qual).map_err(|err| {
            self.stats.report("bad_cell_slice");
            anyhow::anyhow!(err)
        })?;

        let umi = primer_match.get_umi(&r1.seq, &r1.qual).map_err(|err| {
            self.stats.report("bad_umi_slice");
            anyhow::anyhow!(err)
        })?;

        let cell_id = IntToStr::new(&cell.seq).into_u64();
        let umi_id = IntToStr::new(&umi.seq).into_u64();

        if let Some(mapper) = &config.feature_tag_mapper {
            if let Some( id ) =  mapper.map_feature_id(
            &r2.seq,
            &mut self.stats,
            ) {

            self.feature_tag_table.try_insert(
                &cell_id,
                GeneUmiHash( id as u64, umi_id),
                1.0,
                &mut self.stats,
            );
            return Ok(())
            }
        };

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
    feature_tag_table: Scdata,
}

impl IlluminaNormalizer {
    pub fn new(config: IlluminaNormalizerConfig) -> Self {
        Self {
            config,
            stats: NgsNormalizerSupport::new_stats(),
            read_tags: ReadTagTable::new(),
            feature_tag_table: NgsNormalizerSupport::new_feature_tag_table(),
        }
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        Ok(Self::new(IlluminaNormalizerConfig {
            r1: cli.r1,
            r2: cli.r2,
            out: cli.out,
            read_tags: cli.read_tags,
            primer_read: cli.primer_read,
            insert_read: cli.insert_read,
            primer: cli.primer.detector().map_err(anyhow::Error::msg)?,
            feature_tag_mapper: Some(cli.feature_tags.mapper().map_err(anyhow::Error::msg)?),
            min_insert_len: cli.min_insert_len,
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        }))
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &IlluminaNormalizerConfig {
        &self.config
    }

    pub fn run(&mut self) -> Result<()> {
        NgsNormalizerSupport::configure_rayon_threads(self.config.threads);

        let mut seen = HashSet::<DedupKey>::new();

        let mut reader = FastqPairReader::from_paths(&self.config.r1, &self.config.r2)
            .with_context(|| {
                format!(
                    "failed to open FASTQ pair: {} and {}",
                    self.config.r1.display(),
                    self.config.r2.display()
                )
            })?;

        let (paired_r1_path, paired_r2_path) = self.paired_out_paths(&self.config.out);

        std::fs::create_dir_all(
            paired_r1_path
                .parent()
                .unwrap_or_else(|| std::path::Path::new(".")),
        )?;

        let mut fastq = FastqWriter::new(
            &self.config.out,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!( "failed to create FASTQ: {}", self.config.out.display()))?;

        let mut paired_r1 = FastqWriter::new(
            &paired_r1_path,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!("failed to create paired R1 FASTQ: {}", paired_r1_path.display()))?;

        let mut paired_r2 = FastqWriter::new(
            &paired_r2_path,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!("failed to create paired R2 FASTQ: {}", paired_r2_path.display()))?;

        let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(CHUNK_SIZE);
        let mut processed_pairs = 0;

        while let Some(pair) = reader.next_pair()? {
            chunk.push(pair);

            if chunk.len() >= CHUNK_SIZE {
                processed_pairs += CHUNK_SIZE;
                self.process_chunk(&chunk, &mut fastq, &mut paired_r1, &mut paired_r2, &mut seen)?;
                chunk.clear();

                let total = self.stats.get_issue_count("total_pairs") as f64;
                let written = self.stats.get_issue_count("fastq_reads_written") as f64;

                let pct = if total > 0.0 {
                    100.0 * written / total
                } else {
                    0.0
                };

                eprintln!(
                    "processed {} read pairs; written {} ({:.2}%) FASTQ reads; failed {} pairs; {} duplicates detected; sample-tag reads {}",
                    self.stats.get_issue_count("total_pairs"),
                    self.stats.get_issue_count("fastq_reads_written"),
                    pct,
                    self.stats.get_issue_count("failed_pairs"), 
                    self.stats.get_issue_count("duplicate_molecules"),
                    self.stats.get_issue_count("feature_tag_match"),
                );

                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(
                &chunk,
                &mut fastq,
                &mut paired_r1,
                &mut paired_r2,
                &mut seen,
            )?;
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

        NgsNormalizerSupport::write_feature_tag_table_if_present(
            &mut self.feature_tag_table,
            self.config.feature_tag_mapper.as_ref(),
            &self.config.out,
        )?;

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
        fastq: &mut FastqWriter,
        paired_r1: &mut FastqWriter,
        paired_r2: &mut FastqWriter,
        seen: &mut HashSet<DedupKey>,
    ) -> Result<()> {


        let partials: Vec<IlluminaPartial> = input
            .par_iter()
            .map(|(r1, r2)| {
                let mut out = IlluminaPartial::new();

                if out.normalize_pair(r1, r2, &self.config).is_err() {
                    out.stats.report("failed_pairs");
                }

                out
            })
            .collect();

        for partial in partials {
            self.stats.merge(&partial.stats);
            self.feature_tag_table.merge(&partial.feature_tag_table);

            self.stats.stop_multi_processor_time();

            for candidate in partial.candidates {
                if seen.insert(candidate.dedup_key) {
                    if let Some(r1) = candidate.paired_r1_record {
                        paired_r1.write(&r1)?;
                        paired_r2.write(&candidate.fastq_record)?;
                        self.stats.report("paired_fastq_pairs_written");
                    } 
                    fastq.write(&candidate.fastq_record)?;
                    self.stats.report("fastq_reads_written");
                    self.read_tags.insert(candidate.read_tag);
                    self.stats.report("accepted_pairs");
                } else {
                    self.stats.report("duplicate_molecules");
                }
            }
            self.stats.stop_file_io_time();

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
