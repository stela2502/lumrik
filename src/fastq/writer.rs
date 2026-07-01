use anyhow::Result;
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use crate::fastq::record::FastqRecord;

const BUFFER_SIZE: usize = 4 * 1024 * 1024;

pub struct FastqWriter {
    writer: Box<dyn Write>,
    qual_buf: Vec<u8>,
}

impl FastqWriter {
    pub fn new<P: AsRef<Path>>(path: P, gzip: bool, gzip_level: u32) -> Result<Self> {
        let path = path.as_ref();

        let writer: Box<dyn Write> = if path == Path::new("-") {
            // stdout is always plain FASTQ.
            // Use: bam-ont-normalizer --out - --no-tags | pigz -p 8 -1 > out.fastq.gz
            Box::new(BufWriter::with_capacity(64 * 1024, io::stdout()))
        } else {
            let file = File::create(path)?;
            let buffered = BufWriter::with_capacity(BUFFER_SIZE, file);

            if gzip {
                Box::new(GzEncoder::new(
                    buffered,
                    Compression::new(gzip_level.min(9)),
                ))
            } else {
                Box::new(buffered)
            }
        };

        Ok(Self {
            writer,
            qual_buf: Vec::new(),
        })
    }

    pub fn write(&mut self, record: &FastqRecord) -> Result<()> {
        record.write_to(&mut self.writer, &mut self.qual_buf)?;
        Ok(())
    }

    pub fn write_record(&mut self, id: &str, seq: &[u8], qual: &[u8]) -> Result<()> {
        self.writer.write_all(b"@")?;
        self.writer.write_all(id.as_bytes())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(seq)?;
        self.writer.write_all(b"\n+\n")?;

        self.qual_buf.clear();
        self.qual_buf.reserve(qual.len());
        self.qual_buf
            .extend(qual.iter().map(|q| q.saturating_add(33)));

        self.writer.write_all(&self.qual_buf)?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}
