use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;

use crate::fastq::FastqRecord;

pub struct FastqPairReader {
    r1: Box<dyn FastqRead>,
    r2: Box<dyn FastqRead>,
}

impl FastqPairReader {
    pub fn from_paths(r1: &PathBuf, r2: &PathBuf) -> Result<Self> {
        Ok(Self {
            r1: Box::new(SimpleFastqReader::new(r1)?),
            r2: Box::new(SimpleFastqReader::new(r2)?),
        })
    }

    pub fn next_pair(&mut self) -> Result<Option<(FastqRecord, FastqRecord)>> {
        let r1 = self.r1.next_record()?;
        let r2 = self.r2.next_record()?;

        match (r1, r2) {
            (Some(r1), Some(r2)) => Ok(Some((r1, r2))),
            (None, None) => Ok(None),
            (Some(_), None) => bail!("R1 has more records than R2"),
            (None, Some(_)) => bail!("R2 has more records than R1"),
        }
    }
}

pub trait FastqRead: Send {
    fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}

pub struct SimpleFastqReader {
    reader: Box<dyn BufRead + Send>,
    line: String,
}

impl SimpleFastqReader {
    pub fn new(path: &PathBuf) -> Result<Self> {
        Self::from_path(path.as_path())
    }

    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;

        let reader: Box<dyn BufRead + Send> = if path.extension().is_some_and(|ext| ext == "gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(Self {
            reader,
            line: String::new(),
        })
    }
}

impl FastqRead for SimpleFastqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            return Ok(None);
        }

        let id_line = self.line.trim_end().to_string();
        if !id_line.starts_with('@') {
            bail!("invalid FASTQ record: expected @ header, got {id_line}");
        }

        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            bail!("invalid FASTQ record {id_line}: missing sequence line");
        }
        let seq = self.line.trim_end().as_bytes().to_vec();

        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            bail!("invalid FASTQ record {id_line}: missing + line");
        }
        let plus = self.line.trim_end().to_string();
        if !plus.starts_with('+') {
            bail!("invalid FASTQ record: expected + line, got {plus}");
        }

        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            bail!("invalid FASTQ record {id_line}: missing quality line");
        }
        let qual_ascii = self.line.trim_end().as_bytes().to_vec();

        if seq.len() != qual_ascii.len() {
            bail!(
                "invalid FASTQ record {}: sequence length {} != quality length {}",
                id_line,
                seq.len(),
                qual_ascii.len()
            );
        }

        let qual: Vec<u8> = qual_ascii.iter().map(|q| q.saturating_sub(33)).collect();

        Ok(Some(FastqRecord::new(
            id_line.trim_start_matches('@').to_string(),
            &seq,
            &qual,
        )))
    }
}