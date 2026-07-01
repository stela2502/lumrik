// src/encoder/sam.rs

//use std::fs;
//use std::path::Path;

use crate::QuantData;
use crate::encoder::cli::TestDataCli;
use crate::encoder::cli::TruthFeatureMode;
use crate::encoder::model::{EncoderConfig, SamRead};
use crate::index::{GeneFeatureIndex, TranscriptFeatureIndex};

use gtf_splice_index::{
    GeneId, MatchClass, RefBlock, SpliceIndex, Strand, Transcript, TranscriptId,
};
use int_to_str::IntToStr;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use scdata::cell_data::GeneUmiHash;
use snp_index::{Genome, SnpIndex, VcfReadOptions};

pub struct SamEncoder {
    genome: Genome,
    splice_index: SpliceIndex,
    snp_index: SnpIndex,
    config: EncoderConfig,
    read_index: usize,
}

//clippy thinks this is unused...
#[allow(dead_code)]
struct GeneModel<'a> {
    gene_id: GeneId,
    gene_name: String,
    transcript_ids: Vec<TranscriptId>,
    exons: Vec<(usize, &'a RefBlock, Strand)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SnpTruth {
    snp_id: u64,
    is_alt: bool,
}

impl SamEncoder {
    pub fn new(cli: TestDataCli) -> Result<Self, String> {
        let genome = Genome::from_fasta(&cli.fasta)
            .map_err(|e| format!("failed to load FASTA {:?}: {e}", cli.fasta))?;

        let splice_index = SpliceIndex::load(&cli.gtf)
            .map_err(|e| format!("failed to load splice index {:?}: {e}", cli.gtf))?;

        let snp_index = match &cli.vcf {
            Some(vcf) => SnpIndex::from_vcf_path(
                vcf,
                genome.chr_names(),
                genome.chr_lengths(),
                10_000,
                &VcfReadOptions::default(),
            )
            .map_err(|e| format!("failed to load VCF {:?}: {e}", vcf))?,

            None => SnpIndex::new(
                &genome.chr_names(),
                &genome.chr_lengths(),
                Vec::new(),
                10_000,
            )
            .map_err(|e| format!("failed to create empty SNP index: {e}"))?,
        };

        Ok(Self {
            genome,
            splice_index,
            snp_index,
            config: cli.config(),
            read_index: 0,
        })
    }

    /// Generate synthetic SAM records and generator-side truth data.
    ///
    /// Truth data is always collected in memory while reads are generated. If
    /// `truth_path` is `Some(path)`, the collected truth matrices are written by
    /// `QuantData::write()` after all SAM records have been generated.
    ///
    /// Current hidden truth classes:
    ///
    /// - exonic / spliced sense reads -> exonic truth matrix
    /// - exon-intron boundary reads -> intronic truth matrix
    /// - VCF-derived reference / alternate observations -> SNP truth matrices
    /// - antisense reads -> emitted as negative-control SAM records, not counted
    pub fn generate(&mut self) -> Result<Vec<String>, String> {
        self.config.validate()?;
        let truth_path = self.config.truth_path.clone();
        let mut rng = SmallRng::seed_from_u64(self.config.seed);
        let mut lines = Vec::new();
        let mut truth = QuantData::new();

        self.write_header(&mut lines)?;

        let gene_ids: Vec<GeneId> = self.splice_index.genes.iter().map(|gene| gene.id).collect();
        let cells = self.config.cell_barcodes.clone();

        let mut umi = 0_u64;

        if gene_ids.is_empty() {
            return Err("GTF contains no genes".into());
        }

        for cell in &cells {
            for _ in 0..self.config.genes_per_cell {
                let transcript = loop {
                    let gene =
                        &self.splice_index.genes[rng.gen_range(0..self.splice_index.genes.len())];

                    let tx_ids = gene.transcript_ids();

                    if !tx_ids.is_empty() {
                        let tx_id = tx_ids[rng.gen_range(0..tx_ids.len())];
                        break self.splice_index.transcripts[tx_id].clone();
                    }
                };
                self.generate_cell_gene_reads(
                    cell,
                    &transcript,
                    &mut lines,
                    &mut rng,
                    &mut truth,
                    &mut umi,
                )?;
            }
        }

        if let Some(path) = truth_path {
            match self.config.truth_feature_mode {
                TruthFeatureMode::Gene => {
                    let gtf_index = GeneFeatureIndex::new(&self.splice_index);

                    truth.write(
                        path,
                        self.config.min_cell_counts,
                        &gtf_index,
                        Some(&self.snp_index),
                    )?;
                }

                TruthFeatureMode::Transcript => {
                    let gtf_index = TranscriptFeatureIndex::new(&self.splice_index);

                    truth.write(
                        path,
                        self.config.min_cell_counts,
                        &gtf_index,
                        Some(&self.snp_index),
                    )?;
                }
            }
        }

        Ok(lines)
    }

    fn write_header(&self, out: &mut Vec<String>) -> Result<(), String> {
        out.push("@HD\tVN:1.6\tSO:unsorted".into());

        for (chr_id, name) in self.genome.chr_names.iter().enumerate() {
            let len = self
                .genome
                .chr_len(chr_id)
                .ok_or_else(|| format!("missing chromosome length for {name}"))?;

            out.push(format!("@SQ\tSN:{name}\tLN:{len}"));
        }

        out.push("@PG\tID:bam-quant-testdata\tPN:bam-quant-testdata\tVN:0.1".into());
        Ok(())
    }

    /// Generate all synthetic reads for one `(cell, gene)` pair.
    ///
    /// UMI identity is stored as `u64` in truth matrices. The SAM record receives
    /// the nucleotide-string encoding produced by `int_to_str`.
    fn generate_cell_gene_reads<R: Rng>(
        &mut self,
        cell: &str,
        transcript: &Transcript,
        out: &mut Vec<String>,
        rng: &mut R,
        truth: &mut QuantData,
        umi: &mut u64,
    ) -> Result<(), String> {
        if transcript.exons().is_empty() {
            return Err(format!(
                "transcript {:?} cannot generate reads: no exons",
                transcript.id
            ));
        }

        let n_total = self.config.reads_per_cell;
        let n_antisense = (n_total as f64 * self.config.antisense_fraction) as usize;
        let n_intronic = (n_total as f64 * self.config.unspliced_fraction) as usize;
        let n_spliced = n_total.saturating_sub(n_antisense + n_intronic);

        let clean = cell.split_once('-').map_or(cell, |(barcode, _)| barcode);
        let cell_id = IntToStr::new(clean.as_bytes()).into_u64();

        let feature_id = transcript.gene_id as u64;

        for _ in 0..n_spliced {
            let mut read = self.make_spliced_read(transcript, cell, rng)?;
            read.umi = Self::umi_to_string(*umi);

            let snp_truth = self.write_finalized_read(read, out, rng)?;

            truth.gene.try_insert(
                &cell_id,
                GeneUmiHash(feature_id, *umi),
                1.0,
                &mut truth.report,
            );

            Self::record_snp_truth(&cell_id, *umi, snp_truth, truth);
            *umi += 1;
        }

        for _ in 0..n_intronic {
            let mut read = self.make_intronic_read(transcript, cell, rng)?;
            read.umi = Self::umi_to_string(*umi);

            let snp_truth = self.write_finalized_read(read, out, rng)?;

            truth.intron.try_insert(
                &cell_id,
                GeneUmiHash(feature_id, *umi),
                1.0,
                &mut truth.report,
            );

            Self::record_snp_truth(&cell_id, *umi, snp_truth, truth);
            *umi += 1;
        }

        for _ in 0..n_antisense {
            let mut read = self.make_exonic_read(transcript, cell, true, rng)?;
            read.umi = Self::umi_to_string(*umi);

            let snp_truth = self.write_finalized_read(read, out, rng)?;

            Self::record_snp_truth(&cell_id, *umi, snp_truth, truth);
            *umi += 1;
        }

        Ok(())
    }

    fn record_snp_truth(cell: &u64, umi: u64, snp: Option<SnpTruth>, truth: &mut QuantData) {
        if let Some(snp) = snp {
            let target = if snp.is_alt {
                &mut truth.snp_alt
            } else {
                &mut truth.snp_ref
            };

            target.try_insert(cell, GeneUmiHash(snp.snp_id, umi), 1.0, &mut truth.report);
        }
    }

    fn write_finalized_read<R: Rng>(
        &mut self,
        mut read: SamRead,
        out: &mut Vec<String>,
        rng: &mut R,
    ) -> Result<Option<SnpTruth>, String> {
        let snp_truth = self.inject_vcf_snp_if_covered(&mut read, rng)?;
        self.inject_end_weighted_errors(&mut read, rng);

        out.push(read.to_sam_line());
        self.read_index += 1;

        Ok(snp_truth)
    }

    fn make_spliced_read<R: Rng>(
        &self,
        transcript: &Transcript,
        cell: &str,
        rng: &mut R,
    ) -> Result<SamRead, String> {
        let exons = transcript.exons();

        if exons.len() < 2 {
            return self.make_exonic_read(transcript, cell, false, rng);
        }

        let j = rng.gen_range(0..exons.len() - 1);
        let left = &exons[j];
        let right = &exons[j + 1];

        let left_anchor = self.config.seq_len / 2;
        let right_anchor = self.config.seq_len - left_anchor;
        let right_end = right.start + right_anchor as u32;

        if left.end < left_anchor as u32
            || right.start <= left.end
            || right_end > self.genome.chr_len(transcript.chr_id).unwrap_or(0)
        {
            return self.make_exonic_read(transcript, cell, false, rng);
        }

        let start0 = left.end - left_anchor as u32;

        let mut seq = self.fetch(transcript.chr_id, start0, left.end)?;
        seq.push_str(&self.fetch(transcript.chr_id, right.start, right_end)?);

        let mut flag = 0_u16;
        if matches!(transcript.strand, Strand::Minus) {
            seq = Self::revcomp(&seq);
            flag |= 16;
        }

        self.base(
            transcript,
            cell,
            transcript.chr_id,
            start0,
            flag,
            format!("{left_anchor}M{}N{right_anchor}M", right.start - left.end),
            seq,
            MatchClass::ExactJunctionChain,
        )
    }
    /// Create a contiguous genomic read that crosses an exon-intron boundary.
    ///
    /// This intentionally hides pre-mRNA / intronic signal in a normal `M` CIGAR.
    /// It should later be classified as intronic by the splice/gene matcher.
    fn make_intronic_read<R: Rng>(
        &self,
        transcript: &Transcript,
        cell: &str,
        rng: &mut R,
    ) -> Result<SamRead, String> {
        let exons = transcript.exons();

        if exons.len() < 2 {
            return self.make_exonic_read(transcript, cell, false, rng);
        }

        let j = rng.gen_range(0..exons.len() - 1);
        let left = &exons[j];
        let right = &exons[j + 1];

        let read_len = self.config.seq_len as u32;
        let intron_len = right.start.saturating_sub(left.end);
        let exon_anchor = (read_len / 2).min(left.end.saturating_sub(left.start));
        let intron_anchor = read_len.saturating_sub(exon_anchor);

        if right.start <= left.end
            || exon_anchor == 0
            || intron_anchor == 0
            || intron_anchor > intron_len
        {
            return self.make_exonic_read(transcript, cell, false, rng);
        }

        let start0 = left.end - exon_anchor;
        let mut seq = self.fetch(transcript.chr_id, start0, start0 + read_len)?;

        let mut flag = 0_u16;
        if matches!(transcript.strand, Strand::Minus) {
            seq = Self::revcomp(&seq);
            flag |= 16;
        }

        self.base(
            transcript,
            cell,
            transcript.chr_id,
            start0,
            flag,
            format!("{}M", self.config.seq_len),
            seq,
            MatchClass::Intronic,
        )
    }

    fn make_exonic_read<R: Rng>(
        &self,
        transcript: &Transcript,
        cell: &str,
        antisense: bool,
        rng: &mut R,
    ) -> Result<SamRead, String> {
        let exons = transcript.exons();

        if exons.is_empty() {
            return Err(format!(
                "transcript {:?} cannot generate exonic read: no exons",
                transcript.id
            ));
        }

        let exon = &exons[rng.gen_range(0..exons.len())];

        let exon_len = exon.end.saturating_sub(exon.start);
        let read_len = self.config.seq_len as u32;

        if exon.end <= exon.start || exon_len < read_len {
            return Err(format!(
                "transcript {:?} has invalid/short exon interval {}-{} for read length {}",
                transcript.id, exon.start, exon.end, read_len
            ));
        }

        let start0 = exon.start + rng.gen_range(0..=(exon_len - read_len));
        let mut seq = self.fetch(transcript.chr_id, start0, start0 + read_len)?;

        let sense_is_reverse = matches!(transcript.strand, Strand::Minus);
        let read_is_reverse = (antisense && !sense_is_reverse) || (!antisense && sense_is_reverse);

        let mut flag = 0_u16;
        if read_is_reverse {
            seq = Self::revcomp(&seq);
            flag |= 16;
        }

        self.base(
            transcript,
            cell,
            transcript.chr_id,
            start0,
            flag,
            format!("{}M", self.config.seq_len),
            seq,
            MatchClass::Compatible,
        )
    }

    fn base(
        &self,
        transcript: &Transcript,
        cell: &str,
        chr_id: usize,
        pos0: u32,
        flag: u16,
        cigar: String,
        seq: String,
        expected: MatchClass,
    ) -> Result<SamRead, String> {
        let gene_name = self
            .splice_index
            .genes
            .get(transcript.gene_id)
            .and_then(|gene| gene.names.first())
            .cloned()
            .unwrap_or_else(|| format!("{:?}", transcript.gene_id));

        let transcript_name = transcript
            .names
            .first()
            .cloned()
            .unwrap_or_else(|| format!("{:?}", transcript.id));

        Ok(SamRead {
            qname: format!("r{:08}", self.read_index),
            flag,
            chrom: self.genome.chr_names[chr_id].clone(),
            pos_1based: pos0 as u64 + 1,
            mapq: 60,
            cigar,
            seq,
            qual: "I".repeat(self.config.seq_len),
            cell_barcode: cell.into(),
            umi: String::new(),
            gene_id: Some(gene_name),
            gene_name: None,
            transcript_id: Some(transcript_name),
            nm: 0,
            expected,
        })
    }

    /// Optionally record one VCF SNP opportunity for a contiguous read.
    ///
    /// If a covered SNP is selected and `snp_fraction` triggers, the ALT allele
    /// is injected into the read and ALT truth is returned. If a covered SNP is
    /// selected but no ALT is injected, REF truth is returned.
    fn inject_vcf_snp_if_covered<R: Rng>(
        &self,
        read: &mut SamRead,
        rng: &mut R,
    ) -> Result<Option<SnpTruth>, String> {
        if read.cigar.contains('N') || read.cigar.contains('I') || read.cigar.contains('D') {
            return Ok(None);
        }

        let chr_id = self
            .genome
            .chr_id(&read.chrom)
            .ok_or_else(|| format!("read chromosome not in genome: {}", read.chrom))?;

        let start0 = read.pos_1based as u32 - 1;
        let end0 = start0 + read.seq.len() as u32;

        let covered: Vec<_> = self
            .snp_index
            .loci
            .iter()
            .filter(|locus| {
                locus.chr_id == chr_id
                    && locus.pos0 >= start0
                    && locus.pos0 < end0
                    && !locus.alternates.is_empty()
            })
            .collect();

        if covered.is_empty() {
            return Ok(None);
        }

        let snp = covered[rng.gen_range(0..covered.len())];
        let valid_alts: Vec<u8> = snp
            .alternates
            .iter()
            .copied()
            .filter(|alt| *alt != snp.reference)
            .collect();

        if valid_alts.is_empty() {
            return Ok(Some(SnpTruth {
                snp_id: snp.id as u64,
                is_alt: false,
            }));
        }

        if rng.gen_range(0.0..1.0) >= self.config.snp_fraction {
            return Ok(Some(SnpTruth {
                snp_id: snp.id as u64,
                is_alt: false,
            }));
        }

        let mut seq = read.seq.as_bytes().to_vec();
        let genomic_offset = (snp.pos0 - start0) as usize;
        let alt_base = valid_alts[rng.gen_range(0..valid_alts.len())].to_ascii_uppercase();

        let (read_offset, alt_base) = if read.flag & 16 != 0 {
            (
                seq.len() - 1 - genomic_offset,
                Self::complement_base(alt_base),
            )
        } else {
            (genomic_offset, alt_base)
        };

        if read_offset < seq.len() {
            seq[read_offset] = alt_base;
            read.seq = String::from_utf8(seq).map_err(|e| e.to_string())?;
            read.nm += 1;

            Ok(Some(SnpTruth {
                snp_id: snp.id as u64,
                is_alt: true,
            }))
        } else {
            Ok(None)
        }
    }

    fn inject_end_weighted_errors<R: Rng>(&self, read: &mut SamRead, rng: &mut R) {
        let mut seq = read.seq.as_bytes().to_vec();
        let mut qual = vec![b'I'; seq.len()];

        for i in 0..seq.len() {
            let dist_to_right = seq.len() - 1 - i;
            let dist_to_nearest_end = i.min(dist_to_right);

            if dist_to_right < self.config.bad_end_bases {
                qual[i] = b'#';
            }

            let error_probability = if dist_to_nearest_end >= self.config.bad_end_bases {
                self.config.body_error_rate
            } else {
                let t = 1.0 - dist_to_nearest_end as f64 / self.config.bad_end_bases as f64;
                self.config.body_error_rate
                    + t * (self.config.end_error_rate - self.config.body_error_rate)
            };

            if rng.gen_range(0.0..1.0) < error_probability {
                seq[i] = Self::random_different_base(seq[i], rng);
                qual[i] = b'!';
                read.nm += 1;
            }
        }

        read.seq = String::from_utf8(seq).expect("generated sequence should be UTF-8");
        read.qual = String::from_utf8(qual).expect("generated quality should be UTF-8");
    }

    /*
    fn end_weighted_error_probability(&self, pos: usize, len: usize) -> f64 {
        if len == 0 {
            return self.config.body_error_rate;
        }

        let dist_to_end = pos.min(len - 1 - pos);

        if dist_to_end >= self.config.bad_end_bases {
            self.config.body_error_rate
        } else {
            let t = 1.0 - dist_to_end as f64 / self.config.bad_end_bases as f64;
            self.config.body_error_rate
                + t * (self.config.end_error_rate - self.config.body_error_rate)
        }
    }


    fn gene_model(&self, gene_id: GeneId) -> Result<GeneModel<'_>, String> {
        let gene = self.gene(gene_id)?;
        let gene_name = gene
            .names
            .first()
            .cloned()
            .ok_or_else(|| format!("gene {:?} has no names", gene_id))?;

        let transcript_ids = gene.transcript_ids().to_vec();
        let mut exons = Vec::new();

        for tx_id in &transcript_ids {
            let tx = self.transcript(*tx_id)?;
            for exon in tx.exons() {
                exons.push((tx.chr_id, exon, tx.strand));
            }
        }

        Ok(GeneModel {
            gene_id,
            gene_name,
            transcript_ids,
            exons,
        })
    }
    */

    fn fetch(&self, chr_id: usize, start0: u32, end0: u32) -> Result<String, String> {
        let chr_name = self.genome.chr_name(chr_id).unwrap_or("<unknown>");
        let chr_len = self
            .genome
            .chr_len(chr_id)
            .ok_or_else(|| format!("unknown chromosome id {chr_id}"))?;

        if start0 > end0 || end0 > chr_len {
            return Err(format!(
                "invalid genome slice {chr_name}:{start0}-{end0} with chromosome length {chr_len}"
            ));
        }

        let seq = self
            .genome
            .slice(chr_id, start0, end0)
            .ok_or_else(|| format!("failed to fetch genome slice {chr_name}:{start0}-{end0}"))?;

        std::str::from_utf8(seq)
            .map(str::to_string)
            .map_err(|e| format!("reference slice {chr_name}:{start0}-{end0} is not UTF-8: {e}"))
    }

    /*
    fn gene(&self, gene_id: GeneId) -> Result<&gtf_splice_index::Gene, String> {
        self.splice_index
            .genes
            .get(gene_id)
            .ok_or_else(|| format!("bad gene id {:?}", gene_id))
    }

    fn transcript(
        &self,
        transcript_id: TranscriptId,
    ) -> Result<&gtf_splice_index::Transcript, String> {
        self.splice_index
            .transcripts
            .get(transcript_id)
            .ok_or_else(|| format!("bad transcript id {:?}", transcript_id))
    }
    */

    fn umi_to_string(umi: u64) -> String {
        IntToStr::from_u64(umi).to_string(12)
    }

    fn complement_base(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        }
    }

    fn random_different_base<R: Rng>(base: u8, rng: &mut R) -> u8 {
        let choices = match base.to_ascii_uppercase() {
            b'A' => [b'C', b'G', b'T'],
            b'C' => [b'A', b'G', b'T'],
            b'G' => [b'A', b'C', b'T'],
            b'T' => [b'A', b'C', b'G'],
            _ => [b'A', b'C', b'G'],
        };

        choices[rng.gen_range(0..choices.len())]
    }

    fn revcomp(seq: &str) -> String {
        seq.bytes()
            .rev()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                _ => 'N',
            })
            .collect()
    }
}
