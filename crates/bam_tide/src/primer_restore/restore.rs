use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use mapping_info::MappingInfo;
use sc_primer::{Grammar, PrimerDetector};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use crate::fastq::record::FastqRecord;
use crate::fastq::writer::FastqWriter;
use crate::primer_restore::cli::Cli;

use read_tag_table::{
    ReadTagTable,
    ReadTagTableConfig,
    ReadTagRecord,
};


#[derive(Debug, Clone)]
pub struct PrimerRestoreConfig {
    pub fastq: PathBuf,
    pub out: PathBuf,
    pub read_tags: ReadTagTable,
    pub primer: PrimerDetector,
    pub gzip: bool,
    pub gzip_level: u32,
    pub die_on_error: bool,
}

pub struct PrimerRestore {
    config: PrimerRestoreConfig,
    stats: MappingInfo,
}

impl PrimerRestore {
    pub fn new(config: PrimerRestoreConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
        }
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        Ok(Self::new(PrimerRestoreConfig {
            fastq: cli.fastq,
            out: cli.out,
            read_tags: cli.read_tags.load()?,
            primer: cli.primer.detector().map_err(anyhow::Error::msg)?,
            gzip: !cli.no_gzip,
            gzip_level: cli.gzip_level,
            die_on_error: !cli.die_on_error,
        }))
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn run(&mut self) -> Result<()> {
        let mut reader = FastqReader::from_path(&self.config.fastq)
            .with_context(|| format!("failed to open FASTQ: {}", self.config.fastq.display()))?;

        let mut writer = FastqWriter::new(
            &self.config.out,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| {
            format!(
                "failed to create restored FASTQ: {}",
                self.config.out.display()
            )
        })?;

        while let Some(record) = reader.next_record()? {
            self.stats.report("total_reads");

            let Some(tag) = self.config.read_tags.get(record.id.as_str()).cloned() else {
                self.stats.report("missing_read_tag");
                continue;
            };

            let (_target_cell, primer_seq) = match self.config.primer.generate(&tag) {
                Ok(x) => x,
                Err(err) => {
                    self.stats.report("primer_restore_failed");
                    eprintln!(
                        "primer-restore: failed to generate primer for read '{}': {err}",
                        record.id
                    );
                    continue;
                }
            };

            let primer_qual = vec![40u8; primer_seq.len()];

            let mut seq = Vec::with_capacity(primer_seq.len() + record.seq.len());
            seq.extend_from_slice(&primer_seq);
            seq.extend_from_slice(&record.seq);

            let mut qual = Vec::with_capacity(primer_qual.len() + record.qual.len());
            qual.extend_from_slice(&primer_qual);
            qual.extend_from_slice(&record.qual);

            let restored = FastqRecord::new(record.id.clone(), &seq, &qual);

            writer.write(&restored)?;
            self.stats.report("restored_reads");
        }

        writer.finish()?;
        self.finish_or_fail()
    }

    fn finish_or_fail(&self) -> Result<()> {
        let total = self.stats.get_issue_count("total_reads");
        let restored = self.stats.get_issue_count("restored_reads");
        let missing = self.stats.get_issue_count("missing_read_tag");
        let failed = self.stats.get_issue_count("primer_restore_failed");
        let errors = missing + failed;

        eprintln!(
            "\
primer-restore summary
======================
FASTQ reads processed      : {total}
FASTQ reads restored       : {restored}
reads without tag record   : {missing}
primer generation failures : {failed}
"
        );

        if errors > 0 {
            eprintln!(
                "WARNING: primer-restore did not restore {errors} reads. \
                 Output FASTQ is incomplete."
            );
        }

        Ok(())
    }
}

struct FastqReader {
    reader: Box<dyn BufRead>,
    line: String,
}

impl FastqReader {
    fn from_path(path: &Path) -> Result<Self> {
        let reader = open_maybe_gz(path)?;
        Ok(Self {
            reader: Box::new(BufReader::new(reader)),
            line: String::new(),
        })
    }

    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        let id = match self.read_line()? {
            None => return Ok(None),
            Some(line) => line,
        };

        if !id.starts_with('@') {
            bail!("invalid FASTQ record: id line does not start with '@': {id}");
        }

        let seq = self
            .read_line()?
            .context("truncated FASTQ record: missing sequence line")?;

        let plus = self
            .read_line()?
            .context("truncated FASTQ record: missing plus line")?;

        if !plus.starts_with('+') {
            bail!("invalid FASTQ record for {id}: plus line does not start with '+'");
        }

        let qual = self
            .read_line()?
            .context("truncated FASTQ record: missing quality line")?;

        if seq.len() != qual.len() {
            bail!(
                "invalid FASTQ record for {id}: sequence length {} != quality length {}",
                seq.len(),
                qual.len()
            );
        }

        let clean_id = id[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        let qual_phred: Vec<u8> = qual
            .as_bytes()
            .iter()
            .map(|q| q.saturating_sub(33))
            .collect();

        Ok(Some(FastqRecord::new(
            clean_id,
            seq.as_bytes(),
            &qual_phred,
        )))
    }

    fn read_line(&mut self) -> Result<Option<String>> {
        self.line.clear();
        let n = self.reader.read_line(&mut self.line)?;

        if n == 0 {
            return Ok(None);
        }

        while self.line.ends_with('\n') || self.line.ends_with('\r') {
            self.line.pop();
        }

        Ok(Some(self.line.clone()))
    }
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}
