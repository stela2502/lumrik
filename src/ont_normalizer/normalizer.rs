use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::ngs_normalizer::{
    NgsNormalizerSupport, NormalizedMolecule, NormalizerPartial, CHUNK_SIZE,
};
use crate::ont_normalizer::cli::Cli;
use crate::read_tag_table::ReadTagTable;

use fast_tag_mapper::FastTagMapper;

use anyhow::{Context, Result};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Reader};
use sc_primer::{Orientation, PrimerDetector};
use scdata::{Scdata, GeneUmiHash};
use int_to_str::IntToStr;

use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct OntNormalizerConfig {
    pub bam: PathBuf,
    pub out: PathBuf,
    pub read_tags: PathBuf,
    pub min_transcript_len: usize,

    pub primer: PrimerDetector,
    pub feature_tag_mapper: Option<FastTagMapper>,

    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

impl OntNormalizerConfig {
    fn process_read(&self, read: &FastqRecord ) -> NormalizerPartial {
        let mut out = NormalizerPartial::new();
        out.stats.report("total_records");

        let matches = match self.primer.detect_all(&read.seq, &read.qual) {
            Ok(matches) => matches,
            Err(_) => {
                out.stats.report("zero_cassette");
                out.stats.report("no_primer_match");
                return out;
            }
        };

        match matches.len() {
            0 => {
                out.stats.report("zero_cassette");
                return out;
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
            
            let cell_id = IntToStr::new(&cell.seq).into_u64();
            let umi_id = IntToStr::new(&umi.seq).into_u64();


            if let Some(mapper) = &self.feature_tag_mapper
                && let Some(id) = mapper.map_feature_id(&insert_record.seq, &mut out.stats)
            {
                out.feature_tag_table.try_insert(
                    &cell_id,
                    GeneUmiHash(id, umi_id),
                    1.0,
                    &mut out.stats,
                );

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
        }

        out
    }
}

pub struct OntNormalizer {
    config: OntNormalizerConfig,
    stats: MappingInfo,
    read_tags: ReadTagTable,
    feature_tag_table: Scdata,
}

impl OntNormalizer {
    pub fn new(config: OntNormalizerConfig) -> Self {
        Self {
            config,
            stats: NgsNormalizerSupport::new_stats(),
            read_tags: ReadTagTable::new(),
            feature_tag_table: NgsNormalizerSupport::new_feature_tag_table(),
        }
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        Ok(Self::new(OntNormalizerConfig {
            bam: cli.bam,
            out: cli.out,
            read_tags: cli.read_tags,
            min_transcript_len: cli.min_transcript_len,
            primer: cli.primer.detector().map_err(anyhow::Error::msg)?,
            feature_tag_mapper: Some(cli.feature_tags.mapper()?),
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        }))
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &OntNormalizerConfig {
        &self.config
    }

    pub fn write_feature_tag_table_if_present(&mut self) -> Result<()> {
        NgsNormalizerSupport::write_feature_tag_table_if_present(
            &mut self.feature_tag_table,
            self.config.feature_tag_mapper.as_ref(),
            &self.config.out,
        )
    }

    pub fn run(&mut self) -> Result<()> {
        let mut bam = Reader::from_path(&self.config.bam)
            .with_context(|| format!("failed to open BAM: {}", self.config.bam.display()))?;

        if self.config.threads > 1 {
            bam.set_threads(self.config.threads.saturating_sub(1).max(1))
                .context("failed to set BAM reader threads")?;
        }

        let mut fastq = FastqWriter::new(
            &self.config.out,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!("failed to create FASTQ: {}", self.config.out.display()))?;

        let mut chunk: Vec<FastqRecord> = Vec::with_capacity(CHUNK_SIZE);
        let mut processed_pairs = 0;
        for rec_result in bam.records() {
            let rec = rec_result.context("failed to read BAM record")?;
            chunk.push(FastqRecord::from_bam_record(&rec));

            if chunk.len() >= CHUNK_SIZE {
                processed_pairs += CHUNK_SIZE;
                self.process_chunk(&chunk, &mut fastq )?;
                chunk.clear();

                eprintln!(
                    "processed {} read pairs; written {} FASTQ reads; failed {} pairs; sample-tag reads {}",
                    self.stats.get_issue_count("total_pairs"),
                    self.stats.get_issue_count("fastq_reads_written"),
                    self.stats.get_issue_count("failed_pairs"),
                    self.stats.get_issue_count("feature_tag_match"),
                );

                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &mut fastq )?;
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

        self.write_feature_tag_table_if_present()?;

        Ok(())
    }

    fn process_chunk(&mut self, input: &[FastqRecord], fastq: &mut FastqWriter) -> Result<()> {
        let config = self.config.clone();

        let partials: Vec<NormalizerPartial> = input
            .par_iter()
            .map(|read| config.process_read(read ))
            .collect();

        for partial in partials {
            partial.merge_into(
                &mut self.stats,
                &mut self.read_tags,
                &mut self.feature_tag_table,
                fastq,
            )?;
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
