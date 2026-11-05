//! An independent annotation + UMI matrix for newly observed splice isoforms.
use anyhow::{Result, Context};
use std::io::{BufWriter, Write};
use std::path::Path;
use gtf_splice_index::{SpliceIndex, SpliceMatchMode};
use scdata::{MatrixValueType, Scdata, GeneUmiHash, SnpUmiEvidence};

pub struct NewTranscriptData {
    pub scdata: Scdata<SnpUmiEvidence>,
    pub gtf_splice_index: SpliceIndex,
}

impl NewTranscriptData {
    pub fn new(bin_width: u32) -> Self {
        Self {
            scdata: Scdata::new_with_evidence(1, MatrixValueType::Integer),
            gtf_splice_index: SpliceIndex::new(bin_width).with_match_mode(SpliceMatchMode::Transcript),
        }
    }

    pub fn add_molecule(&mut self, transcript_id: usize, cell: u64, umi: u64, evidence: SnpUmiEvidence) {
        self.scdata.try_insert_evidence(&cell, GeneUmiHash(transcript_id as u64, umi), evidence);
    }

    /// Export novel isoforms and per-molecule SNP support without changing normal SNP output.
    pub fn write_reports(&mut self, output: impl AsRef<Path>) -> Result<()> {
        let output = output.as_ref();
        if self.gtf_splice_index.transcripts.is_empty() { return Ok(()); }
        std::fs::create_dir_all(output)?;
        let gtf_path = output.join("new_transcripts.gtf");
        let mut gtf = BufWriter::new(std::fs::File::create(&gtf_path)
            .with_context(|| format!("creating {}", gtf_path.display()))?);
        for tx in &self.gtf_splice_index.transcripts {
            let gene = &self.gtf_splice_index.genes[tx.gene_id];
            let gene_name = gene.primary_name().unwrap_or("NA");
            let name = tx.primary_name().unwrap_or("NA");
            let chr = &self.gtf_splice_index.chr_names[tx.chr_id];
            let strand = match tx.strand {
                gtf_splice_index::Strand::Plus => "+",
                gtf_splice_index::Strand::Minus => "-",
                gtf_splice_index::Strand::Unknown => ".",
            };
            for exon in tx.exons() {
                writeln!(gtf, "{chr}\tlumrik\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{name}\";", exon.start + 1, exon.end)?;
            }
        }
        let detail_path = output.join("new_transcript_snp_evidence.tsv");
        let mut detail = BufWriter::new(std::fs::File::create(&detail_path)?);
        writeln!(detail, "cell_id\ttranscript_id\tumi\tsnp_id\tref_reads\talt_reads")?;
        let mut write_error = None;
        self.scdata.for_each_umi(|cell, feature_umi, evidence| {
            if write_error.is_some() { return; }
            let name = self.gtf_splice_index.transcripts
                .get(feature_umi.0 as usize)
                .and_then(|tx| tx.primary_name()).unwrap_or("NA");
            if evidence.snps.is_empty() {
                if let Err(e) = writeln!(detail, "{cell}\t{name}\t{}\t.\t0\t0", feature_umi.1) { write_error = Some(e); }
            } else {
                for snp in &evidence.snps {
                    if let Err(e) = writeln!(detail, "{cell}\t{name}\t{}\t{}\t{}\t{}", feature_umi.1, snp.snp_id, snp.ref_reads, snp.alt_reads) {
                        write_error = Some(e); break;
                    }
                }
            }
        });
        if let Some(err) = write_error { return Err(err.into()); }
        let mut conflicts = 0usize;
        self.scdata.for_each_umi(|cell, key, evidence| {
            for snp in evidence.conflicts() {
                conflicts += 1;
                if conflicts <= 20 {
                    eprintln!("[bam_tide] WARNING: conflicting novel-transcript SNP cell={cell} transcript={} umi={} snp={} REF={} ALT={}",
                        key.0, key.1, snp.snp_id, snp.ref_reads, snp.alt_reads);
                }
            }
        });
        if conflicts > 0 {
            eprintln!("[bam_tide] WARNING: {conflicts} conflicting novel-transcript molecule/SNP observations (first 20 shown)");
        }

        // Keep the per-UMI SNP report above independent from the ordinary
        // transcript-by-cell UMI matrix. Export every observed novel-transcript
        // cell, without applying the normal GEX cell-calling threshold.
        let feature_index = self.gtf_splice_index.feature_index();
        let cells = self.scdata.cell_ids();
        self.scdata.finalize_for_cells(&cells, &feature_index);
        let matrix_dir = output.join("new_transcripts_mtx");
        let summary = self.scdata.write_sparse(&matrix_dir, &feature_index)
            .map_err(anyhow::Error::msg)
            .with_context(|| format!("writing novel transcript matrix to {}", matrix_dir.display()))?;
        eprintln!("[bam_tide] {summary}");
        Ok(())
    }

    /// Index first, then rewrite the source matrix in place, then deduplicate UMIs.
    pub fn merge(&mut self, mut source: Self) -> Result<()> {
        let mapping = self.gtf_splice_index.merge_novel_transcripts(&source.gtf_splice_index)?;
        source.scdata.remap_feature_ids(&mapping).map_err(anyhow::Error::msg)?;
        self.scdata.merge_single_thread(source.scdata);
        Ok(())
    }
}
