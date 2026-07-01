//! Translate transcriptome-coordinate BAM records back to genome coordinates.
//!
//! This module intentionally keeps the design small:
//!
//! - `BamTranscriptomeMapper` owns the `SpliceIndex`, chromosome lengths, and stats.
//! - `make_header()` creates a genome-coordinate BAM header.
//! - `map_record_in_place()` rewrites one BAM record from transcriptome coordinates to
//!   genome coordinates.
//! - CIGAR translation is implemented as private methods on `BamTranscriptomeMapper`.
//!
//! Assumed external API from `gtf_splice_index` / local crates:
//!
//! - `SpliceIndex::from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<SpliceIndex>`
//! - `SpliceIndex::transcript_name(&self, name: &str) -> anyhow::Result<&Transcript>`
//! - `SpliceIndex::gene_name(&self, gene_id: usize) -> anyhow::Result<&Gene>`
//! - `Transcript::exons(&self) -> &[RefBlock]`
//! - `Transcript::primary_name(&self) -> &str`
//! - `Gene::primary_name(&self) -> &str`
//!
//! If these names differ, only the tiny lookup/accessor calls below need adjustment.

use anyhow::{anyhow, bail, Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Cigar, CigarString};
use std::fmt;
use std::path::Path;

use gtf_splice_index::{RefBlock, SpliceIndex, Strand, Transcript, IdNameKeys};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MappedCigar {
    pub chr_id: usize, //local might differ from genome!
    pub start0: u32,
    pub cigar: Vec<Cigar>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TranscriptomeToGenomeStats {
    pub seen: u64,
    pub written: u64,
    pub unmapped: u64,
    pub unknown_transcript: u64,
    pub failed: u64,
    pub reverse_complemented: u64,
}

impl fmt::Display for TranscriptomeToGenomeStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Transcriptome-to-genome mapping stats")?;
        writeln!(f, "  seen                    : {}", self.seen)?;
        writeln!(f, "  written                 : {}", self.written)?;
        writeln!(f, "  unmapped                : {}", self.unmapped)?;
        writeln!(f, "  unknown_transcript      : {}", self.unknown_transcript)?;
        writeln!(f, "  failed                  : {}", self.failed)?;
        writeln!(f, "  reverse_complemented    : {}", self.reverse_complemented)?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct BamTranscriptomeMapper {
    index: SpliceIndex,
    chr_lengths: Vec<u32>,
    genome_chr_ids: Vec<usize>,
    stats: TranscriptomeToGenomeStats,
    genome: snp_index::Genome,
}

impl BamTranscriptomeMapper {

    pub fn new<P: AsRef<Path>, Q: AsRef<Path>>(gtf: P, fasta: Q) -> Result<Self> {
        let gtf = gtf.as_ref();
        let fasta = fasta.as_ref();

        let index = match gtf.extension().and_then(|s| s.to_str()) {
            Some("dat") => {
                SpliceIndex::load(gtf)
                    .with_context(|| {
                        format!("loading binary SpliceIndex from {}", gtf.display())
                    })?
            }

            _ => {
                SpliceIndex::from_path(gtf, 1_000_000, IdNameKeys::default() )
                    .with_context(|| {
                        format!("building SpliceIndex from {}", gtf.display())
                    })?
            }
        };

        let genome = snp_index::Genome::from_fasta(fasta)
            .with_context(|| format!("reading genome FASTA {}", fasta.display()))?;

        Self::from_index_and_genome(index, genome)
    }

    pub fn from_index_and_genome(
        index: SpliceIndex,
        genome: snp_index::Genome,
    ) -> Result<Self> {

        /*
        if genome.chr_names.len() != index.chr_names.len() {
            bail!(
                "genome FASTA has {} chromosomes but SpliceIndex has {}",
                genome.chr_names.len(),
                index.chr_names.len()
            );
        }*/

        let mut genome_chr_ids = Vec::with_capacity(index.chr_names.len());
        let mut chr_lengths = Vec::with_capacity(index.chr_names.len());

        
        for chr_name in &index.chr_names {
            let genome_chr_id = genome
                .chr_id(chr_name)
                .ok_or_else(|| {
                    anyhow!(
                        "chromosome {chr_name:?} from GTF not present in genome FASTA"
                    )
                })?;

            let len = genome
                .chr_len(genome_chr_id)
                .ok_or_else(|| {
                    anyhow!(
                        "failed to get length for chromosome {chr_name:?}"
                    )
                })?;

            genome_chr_ids.push(genome_chr_id);
            chr_lengths.push(len);
        }

        Ok(Self {
            chr_lengths,
            index,
            genome,
            genome_chr_ids,
            stats: TranscriptomeToGenomeStats::default(),
        })
    }


    pub fn from_path_and_fasta<P: AsRef<Path>, Q: AsRef<Path>>(
        gtf: P,
        fasta: Q,
    ) -> Result<Self> {
        Self::new(gtf, fasta)
    }

    pub fn index(&self) -> &SpliceIndex {
        &self.index
    }

    pub fn stats(&self) -> &TranscriptomeToGenomeStats {
        &self.stats
    }

    pub fn make_header(&self, old_header: &bam::HeaderView) -> Result<Header> {
        let mut header = Header::new();

        self.copy_non_sq_header_records(old_header, &mut header)?;
        self.push_genome_sq_records(&mut header)?;

        Ok(header)
    }

    /// Map one BAM record in place.
    ///
    /// Returns:
    /// - `Ok(true)` if the record should be written.
    /// - `Ok(false)` if the record was intentionally skipped.
    ///
    /// Unmapped input records are currently passed through unchanged and should be written.
    /// Unknown transcript names are skipped because we cannot place them on the genome.
    pub fn map_record_in_place(
        &mut self,
        record: &mut bam::Record,
        old_header: &bam::HeaderView,
    ) -> Result<bool> {
        self.stats.seen += 1;

        if record.is_unmapped() || record.tid() < 0 {
            self.stats.unmapped += 1;
            self.stats.written += 1;
            return Ok(true);
        }

        let transcript_name = self.transcript_name_for_record(record, old_header)?;

        let tx = match self.index.transcript_by_name(transcript_name) {
            Ok(tx) => tx,
            Err(_) => {
                self.stats.unknown_transcript += 1;
                return Ok(false);
            }
        };

        let tx_start0 = u32::try_from(record.pos()).with_context(|| {
            format!(
                "negative transcriptome position {} for read {:?}",
                record.pos(),
                String::from_utf8_lossy(record.qname())
            )
        })?;

        let old_cigar: Vec<Cigar> = record.cigar().iter().copied().collect();

        let strand = if tx.strand == Strand::Unknown {
            self.infer_record_strand(&tx, record, tx_start0, &old_cigar)?
        } else {
            tx.strand
        };

        let mut mapped_tx = tx.clone();
        mapped_tx.strand = strand;

        let mapped = match self.map_cigar(&mapped_tx, tx_start0, &old_cigar) {
            Ok(mapped) => mapped,
            Err(err) => {
                self.stats.failed += 1;
                eprintln!(
                    "WARNING: failed to map read {:?} on transcript {transcript_name:?}: {err}",
                    String::from_utf8_lossy(record.qname())
                );

                return Ok(false);
            }
        };

        if mapped_tx.strand == Strand::Minus {
            self.reverse_complement_record(record);
            self.stats.reverse_complemented += 1;
        }

        self.add_original_alignment_tags(record, transcript_name, tx_start0, &old_cigar)?;

        record.set_tid(i32::try_from(mapped_tx.chr_id).context("chr_id does not fit i32")?);
        record.set_pos(i64::from(mapped.start0));
        record.set_cigar(Some(&CigarString(mapped.cigar)));

        self.stats.written += 1;
        Ok(true)
    }

    fn copy_non_sq_header_records(
        &self,
        old_header: &bam::HeaderView,
        new_header: &mut Header,
    ) -> Result<()> {
        let text = std::str::from_utf8(old_header.as_bytes())
            .context("BAM header is not valid UTF-8")?;

        for line in text.lines() {
            if line.starts_with("@SQ") {
                continue;
            }

            if line.trim().is_empty() {
                continue;
            }

            self.push_raw_header_line(new_header, line)?;
        }

        Ok(())
    }

    fn push_raw_header_line(&self, header: &mut Header, line: &str) -> Result<()> {
        let mut fields = line.split('\t');
        let kind = fields
            .next()
            .ok_or_else(|| anyhow!("empty BAM header line"))?;

        if !kind.starts_with('@') || kind.len() != 3 {
            bail!("invalid BAM header record kind in line {line:?}");
        }

        let tag = &kind.as_bytes()[1..3];
        let mut rec = HeaderRecord::new(tag);

        for field in fields {
            let Some((key, value)) = field.split_once(':') else {
                // @CO records are free text after @CO. Preserve them as a CO tag.
                if tag == b"CO" {
                    rec.push_tag(b"CO", &field);
                    continue;
                }
                bail!("invalid BAM header field {field:?} in line {line:?}");
            };

            if key.len() != 2 {
                bail!("invalid BAM header tag {key:?} in line {line:?}");
            }

            rec.push_tag(key.as_bytes(), &value);
        }

        header.push_record(&rec);
        Ok(())
    }

    fn push_genome_sq_records(&self, header: &mut Header) -> Result<()> {
        for (chr_id, chr_name) in self.index.chr_names.iter().enumerate() {
            let len = self
                .chr_lengths
                .get(chr_id)
                .copied()
                .ok_or_else(|| anyhow!("missing length for chromosome id {chr_id}"))?;

            let mut rec = HeaderRecord::new(b"SQ");
            rec.push_tag(b"SN", chr_name);
            rec.push_tag(b"LN", &len);
            header.push_record(&rec);
        }

        Ok(())
    }

    fn transcript_name_for_record<'a>(
        &self,
        record: &bam::Record,
        old_header: &'a bam::HeaderView,
    ) -> Result<&'a str> {
        let tid = record.tid();
        if tid < 0 {
            bail!("cannot get transcript name for negative tid {tid}");
        }

        let name = old_header.tid2name(tid as u32);
        std::str::from_utf8(name).with_context(|| {
            format!(
                "transcript name for tid {tid} is not valid UTF-8: {:?}",
                name
            )
        })
    }

    fn add_original_alignment_tags(
        &self,
        record: &mut bam::Record,
        transcript_name: &str,
        tx_start0: u32,
        old_cigar: &[Cigar],
    ) -> Result<()> {
        let old_cigar_string = CigarString(old_cigar.to_vec()).to_string();

        record
            .push_aux(b"TX", bam::record::Aux::String(transcript_name))
            .context("adding TX tag")?;
        record
            .push_aux(b"XP", bam::record::Aux::I32(i32::try_from(tx_start0)?))
            .context("adding XP tag")?;
        record
            .push_aux(b"XC", bam::record::Aux::String(&old_cigar_string))
            .context("adding XC tag")?;

        Ok(())
    }

    fn infer_record_strand(
        &self,
        tx: &Transcript,
        record: &bam::Record,
        tx_start0: u32,
        cigar: &[Cigar],
    ) -> Result<Strand> {

        let mut tx = tx.clone();
        tx.strand = Strand::Plus;
        let plus = self.map_cigar(&tx, tx_start0, cigar)?;
        tx.strand = Strand::Minus;
        let minus = self.map_cigar(&tx, tx_start0, cigar)?;

        let read_seq = record.seq().as_bytes();

        
        let plus_score = self.score_read_against_genome(&read_seq, &plus, Strand::Plus)?;
        let minus_score = self.score_read_against_genome(&read_seq, &minus, Strand::Minus)?;
    

        if plus_score > minus_score {
            Ok(Strand::Plus)
        } else if minus_score > plus_score {
            Ok(Strand::Minus)
        } else {
            bail!(
                "cannot infer unknown strand for read {:?} transcript {:?}: plus/minus score tie {plus_score}",
                String::from_utf8_lossy(record.qname()),
                tx.names
            )
        }
    }

    /*
    fn debug_score_read_against_genome(
        &self,
        label: &str,
        read_seq: &[u8],
        mapped: &MappedCigar,
        strand: Strand,
    ) -> Result<i64> {
        let mut read_pos = 0usize;
        let mut genome_pos = mapped.start0;
        let mut score = 0i64;

        let genome_chr_id = *self
            .genome_chr_ids
            .get(mapped.chr_id)
            .ok_or_else(|| anyhow!("missing genome chromosome id for mapped index chr {}", mapped.chr_id))?;

        eprintln!();
        eprintln!("DEBUG score {label}");
        eprintln!("  mapped.chr_id={}", mapped.chr_id);
        eprintln!("  genome_chr_id={genome_chr_id}");
        eprintln!("  start0={}", mapped.start0);
        eprintln!("  strand={strand:?}");
        eprintln!("  cigar={}", CigarString(mapped.cigar.clone()));

        for op in &mapped.cigar {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for _ in 0..n {
                        let raw_read_base = read_seq[read_pos];

                        let read_base = match strand {
                            Strand::Plus | Strand::Unknown => raw_read_base,
                            Strand::Minus => {
                                self.complement_base(read_seq[read_seq.len() - 1 - read_pos])
                            }
                        };

                        let genome_base = self
                            .genome
                            .base(genome_chr_id, genome_pos)
                            .ok_or_else(|| anyhow!("no genome base at FASTA chr id {genome_chr_id}, position {genome_pos}"))?;

                        let s = self.base_match_score(read_base, genome_base);
                        score += s;

                        eprintln!(
                            "  read_pos={read_pos:02} genome_pos={genome_pos:03} raw={} scored_read={} genome={} score={s:+} total={score:+}",
                            raw_read_base as char,
                            read_base as char,
                            genome_base as char,
                        );

                        read_pos += 1;
                        genome_pos += 1;
                    }
                }
                Cigar::Ins(n) | Cigar::SoftClip(n) => {
                    eprintln!("  {op:?}: consume read {n}");
                    read_pos += n as usize;
                }
                Cigar::Del(n) | Cigar::RefSkip(n) => {
                    eprintln!("  {op:?}: consume genome {n}");
                    genome_pos += n;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }

        eprintln!("  FINAL {label}: {score}");

        Ok(score)
    }*/

    fn score_read_against_genome(
        &self,
        read_seq: &[u8],
        mapped: &MappedCigar,
        strand: Strand,
    ) -> Result<i64> {
        let mut read_pos = 0usize;
        let mut genome_pos = mapped.start0;
        let mut score = 0i64;

        for op in &mapped.cigar {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for _ in 0..n {
                        if read_pos >= read_seq.len() {
                            bail!("read sequence shorter than CIGAR");
                        }

                        let read_base = match strand {
                            Strand::Plus | Strand::Unknown => read_seq[read_pos],
                            Strand::Minus => {
                                self.complement_base(read_seq[read_seq.len() - 1 - read_pos])
                            }
                        };

                        let genome_chr_id = *self
                            .genome_chr_ids
                            .get(mapped.chr_id)
                            .ok_or_else(|| {
                                anyhow!(
                                    "missing genome chromosome id for mapped index chr {}",
                                    mapped.chr_id
                                )
                            })?;

                        let genome_base = self
                            .genome
                            .base(genome_chr_id, genome_pos)
                            .ok_or_else(|| {
                                anyhow!(
                                    "no genome base at FASTA chr id {genome_chr_id}, position {genome_pos}"
                                )
                            })?;

                        score += self.base_match_score(read_base, genome_base);

                        read_pos += 1;
                        genome_pos += 1;
                    }
                }
                Cigar::Ins(n) | Cigar::SoftClip(n) => {
                    read_pos += n as usize;
                }
                Cigar::Del(n) | Cigar::RefSkip(n) => {
                    genome_pos += n;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }

        Ok(score)
    }

    fn base_match_score(&self, a: u8, b: u8) -> i64 {
        let a = a.to_ascii_uppercase();
        let b = b.to_ascii_uppercase();

        if a == b'N' || b == b'N' {
            0
        } else if a == b {
            2
        } else {
            -1
        }
    }


    fn map_cigar(
        &self,
        tx: &Transcript,
        tx_start0: u32,
        cigar: &[Cigar],
    ) -> Result<MappedCigar> {
        let exons = self.exons_in_transcript_order(tx);

        if exons.is_empty() {
            bail!("transcript {:?} has no exons", tx.names);
        }

        let mut exon_idx = self.exon_index_for_tx_pos(&exons, tx_start0)
            .with_context(|| format!("transcript start {tx_start0} outside transcript {:?}", tx.names))?;

        let mut tx_pos0 = tx_start0;
        let mut genome_pos0 = self.tx_pos_to_genome_in_ordered_exons(tx, &exons, tx_start0)?;
        let mut leftmost: Option<u32> = None;
        let mut out = Vec::with_capacity(cigar.len() + 8);

        for op in cigar {
            match *op {
                Cigar::Match(n) => {
                    self.consume_transcript_ref_op(
                        tx,
                        &exons,
                        &mut exon_idx,
                        &mut tx_pos0,
                        &mut genome_pos0,
                        &mut leftmost,
                        &mut out,
                        n,
                        Cigar::Match,
                    )?;
                }
                Cigar::Equal(n) => {
                    self.consume_transcript_ref_op(
                        tx,
                        &exons,
                        &mut exon_idx,
                        &mut tx_pos0,
                        &mut genome_pos0,
                        &mut leftmost,
                        &mut out,
                        n,
                        Cigar::Equal,
                    )?;
                }
                Cigar::Diff(n) => {
                    self.consume_transcript_ref_op(
                        tx,
                        &exons,
                        &mut exon_idx,
                        &mut tx_pos0,
                        &mut genome_pos0,
                        &mut leftmost,
                        &mut out,
                        n,
                        Cigar::Diff,
                    )?;
                }
                Cigar::Del(n) => {
                    self.consume_transcript_ref_op(
                        tx,
                        &exons,
                        &mut exon_idx,
                        &mut tx_pos0,
                        &mut genome_pos0,
                        &mut leftmost,
                        &mut out,
                        n,
                        Cigar::Del,
                    )?;
                }
                Cigar::Ins(n) => out.push(Cigar::Ins(n)),
                Cigar::SoftClip(n) => out.push(Cigar::SoftClip(n)),
                Cigar::HardClip(n) => out.push(Cigar::HardClip(n)),
                Cigar::Pad(n) => out.push(Cigar::Pad(n)),
                Cigar::RefSkip(n) => {
                    bail!("unexpected RefSkip/N in transcriptome-space CIGAR: {n}N")
                }
            }
        }

        let start0 = leftmost.ok_or_else(|| {
            anyhow!(
                "CIGAR for transcript {:?} produced no genome-consuming aligned bases",
                tx.names
            )
        })?;

        if tx.strand == Strand::Minus {
            out.reverse();
        }

        self.merge_cigar(&mut out);

        Ok(MappedCigar {
            chr_id: tx.chr_id,
            start0,
            cigar: out,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn consume_transcript_ref_op<F>(
        &self,
        tx: &Transcript,
        exons: &[RefBlock],
        exon_idx: &mut usize,
        tx_pos0: &mut u32,
        genome_pos0: &mut u32,
        leftmost: &mut Option<u32>,
        out: &mut Vec<Cigar>,
        mut len: u32,
        make_op: F,
    ) -> Result<()>
    where
        F: Fn(u32) -> Cigar + Copy,
    {
        while len > 0 {
            if self.cursor_is_past_current_exon(tx.strand, exons, *exon_idx, *genome_pos0)? {
                let intron_len = self.move_to_next_exon(tx, exons, exon_idx, genome_pos0)?;
                out.push(Cigar::RefSkip(intron_len));
            }
            let exon = exons
                .get(*exon_idx)
                .ok_or_else(|| anyhow!("CIGAR extends beyond transcript {:?}", tx.names))?;

            let available = self.remaining_in_current_exon(tx.strand, exon, *genome_pos0)?;
            let take = available.min(len);

            self.mark_leftmost(leftmost, tx.strand, *genome_pos0, take);
            out.push(make_op(take));

            self.advance_position(tx.strand, genome_pos0, take)?;
            *tx_pos0 = tx_pos0
                .checked_add(take)
                .ok_or_else(|| anyhow!("transcript coordinate overflow"))?;
            len -= take;

            if len > 0 {
                let intron_len = self.move_to_next_exon(
                    tx,
                    exons,
                    exon_idx,
                    genome_pos0,
                )?;
                out.push(Cigar::RefSkip(intron_len));
            }
        }

        Ok(())
    }

    fn cursor_is_past_current_exon(
        &self,
        strand: Strand,
        exons: &[RefBlock],
        exon_idx: usize,
        genome_pos0: u32,
    ) -> Result<bool> {
        let exon = exons
            .get(exon_idx)
            .ok_or_else(|| anyhow!("CIGAR extends beyond transcript"))?;

        Ok(match strand {
            Strand::Plus => genome_pos0 >= exon.end,
            Strand::Minus => genome_pos0 < exon.start,
            _ => bail!("cannot map unknown-strand transcript"),
        })
    }

    fn exons_in_transcript_order(&self, tx: &Transcript) -> Vec<RefBlock> {
        let mut exons = tx.exons().to_vec();
        if tx.strand == Strand::Minus {
            exons.reverse();
        }
        exons
    }

    fn exon_index_for_tx_pos(&self, exons_in_tx_order: &[RefBlock], tx_pos0: u32) -> Result<usize> {
        let mut cursor = 0u32;

        for (idx, exon) in exons_in_tx_order.iter().enumerate() {
            let len = self.ref_block_len(*exon)?;
            let next = cursor
                .checked_add(len)
                .ok_or_else(|| anyhow!("transcript length overflow"))?;

            if tx_pos0 >= cursor && tx_pos0 < next {
                return Ok(idx);
            }

            cursor = next;
        }

        bail!("position {tx_pos0} is outside transcript length {cursor}")
    }

    fn tx_pos_to_genome_in_ordered_exons(
        &self,
        tx: &Transcript,
        exons_in_tx_order: &[RefBlock],
        tx_pos0: u32,
    ) -> Result<u32> {
        let mut cursor = 0u32;

        for exon in exons_in_tx_order {
            let len = self.ref_block_len(*exon)?;
            let next = cursor
                .checked_add(len)
                .ok_or_else(|| anyhow!("transcript length overflow"))?;

            if tx_pos0 >= cursor && tx_pos0 < next {
                let offset = tx_pos0 - cursor;
                return match tx.strand {
                    Strand::Plus => Ok(exon.start + offset),
                    Strand::Minus => Ok(exon.end - 1 - offset),
                    Strand::Unknown => bail!("cannot map unknown-strand transcript"),
                };
            }

            cursor = next;
        }

        bail!("position {tx_pos0} is outside transcript length {cursor}")
    }

    fn remaining_in_current_exon(
        &self,
        strand: Strand,
        exon: &RefBlock,
        genome_pos0: u32,
    ) -> Result<u32> {
        match strand {
            Strand::Plus => exon
                .end
                .checked_sub(genome_pos0)
                .filter(|n| *n > 0)
                .ok_or_else(|| {
                    anyhow!(
                        "genome position {genome_pos0} not inside plus-strand exon [{}, {})",
                        exon.start,
                        exon.end
                    )
                }),
            Strand::Minus => genome_pos0
                .checked_sub(exon.start)
                .map(|n| n + 1)
                .filter(|n| *n > 0)
                .ok_or_else(|| {
                    anyhow!(
                        "genome position {genome_pos0} not inside minus-strand exon [{}, {})",
                        exon.start,
                        exon.end
                    )
                }),
            _ => bail!("cannot map unknown-strand transcript"),
        }
    }

    fn mark_leftmost(
        &self,
        leftmost: &mut Option<u32>,
        strand: Strand,
        genome_pos0: u32,
        len: u32,
    ) {
        if len == 0 {
            return;
        }

        let candidate = match strand {
            Strand::Plus => genome_pos0,
            Strand::Minus => genome_pos0 + 1 - len,
            _ => unreachable!(),
        };

        *leftmost = Some(leftmost.map_or(candidate, |old| old.min(candidate)));
    }

    fn advance_position(
        &self,
        strand: Strand,
        genome_pos0: &mut u32,
        len: u32,
    ) -> Result<()> {
        match strand {
            Strand::Plus => {
                *genome_pos0 = genome_pos0
                    .checked_add(len)
                    .ok_or_else(|| anyhow!("genomic coordinate overflow"))?;
            }
            Strand::Minus => {
                *genome_pos0 = genome_pos0
                    .checked_sub(len)
                    .ok_or_else(|| anyhow!("genomic coordinate underflow"))?;
            }
            _ => bail!("cannot map unknown-strand transcript"),
        }

        Ok(())
    }

    fn move_to_next_exon(
        &self,
        tx: &Transcript,
        exons: &[RefBlock],
        exon_idx: &mut usize,
        genome_pos0: &mut u32,
    ) -> Result<u32> {
        let current = *exons
            .get(*exon_idx)
            .ok_or_else(|| anyhow!("missing current exon"))?;

        *exon_idx += 1;

        let next = *exons
            .get(*exon_idx)
            .ok_or_else(|| anyhow!("CIGAR extends beyond transcript"))?;

        let intron_len = match tx.strand {
            Strand::Plus => next
                .start
                .checked_sub(current.end)
                .ok_or_else(|| anyhow!("overlapping or unsorted plus-strand exons"))?,
            Strand::Minus => current
                .start
                .checked_sub(next.end)
                .ok_or_else(|| anyhow!("overlapping or unsorted minus-strand exons"))?,
            gtf_splice_index::Strand::Unknown => bail!("cannot map unknown-strand transcript" ),
        };

        *genome_pos0 = match tx.strand {
            Strand::Plus => next.start,
            Strand::Minus => next.end - 1,
            _ => bail!("cannot map unknown-strand transcript"),
        };

        Ok(intron_len)
    }

    fn ref_block_len(&self, block: RefBlock) -> Result<u32> {
        block
            .end
            .checked_sub(block.start)
            .filter(|len| *len > 0)
            .ok_or_else(|| anyhow!("invalid RefBlock [{}, {})", block.start, block.end))
    }

    fn reverse_complement_record(&self, record: &mut bam::Record) {
        let mut seq = record.seq().as_bytes();
        let mut qual = record.qual().to_vec();
        let qname = record.qname().to_vec();
        let cigar = record.cigar().take();

        seq.reverse();
        for base in &mut seq {
            *base = self.complement_base(*base);
        }
        qual.reverse();

        record.set(&qname, Some(&cigar), &seq, &qual);
        record.set_reverse();
    }


    fn complement_base(&self, base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'U' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            b'N' => b'N',
            _ => b'N',
        }
    }

    fn merge_cigar(&self, cigar: &mut Vec<Cigar>) {
        let mut out: Vec<Cigar> = Vec::with_capacity(cigar.len());

        for op in cigar.drain(..) {
            match (out.last_mut(), op) {
                (Some(Cigar::Match(a)), Cigar::Match(b)) => *a += b,
                (Some(Cigar::Ins(a)), Cigar::Ins(b)) => *a += b,
                (Some(Cigar::Del(a)), Cigar::Del(b)) => *a += b,
                (Some(Cigar::RefSkip(a)), Cigar::RefSkip(b)) => *a += b,
                (Some(Cigar::SoftClip(a)), Cigar::SoftClip(b)) => *a += b,
                (Some(Cigar::HardClip(a)), Cigar::HardClip(b)) => *a += b,
                (Some(Cigar::Pad(a)), Cigar::Pad(b)) => *a += b,
                (Some(Cigar::Equal(a)), Cigar::Equal(b)) => *a += b,
                (Some(Cigar::Diff(a)), Cigar::Diff(b)) => *a += b,
                _ => out.push(op),
            }
        }

        *cigar = out;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use std::path::PathBuf;
    use gtf_splice_index::IdNameKeys;

    fn cigar_string(cigar: &[Cigar]) -> String {
        CigarString(cigar.to_vec()).to_string()
    }

    fn assert_cigar_eq(got: &[Cigar], expected: &[Cigar]) {
        assert_eq!(
            cigar_string(got),
            cigar_string(expected),
            "CIGAR mismatch: got {}, expected {}",
            cigar_string(got),
            cigar_string(expected)
        );
    }

    fn tmp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        let pid = std::process::id();
        let t = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("bam_tide_tx2genome_{name}_{pid}_{t}.gtf"));
        p
    }

    fn write_test_gtf() -> PathBuf {
        let path = tmp_path("projection");
        let mut f = fs::File::create(&path).unwrap();

        let gtf = "\
chrA\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"gene_plus_single\"; transcript_id \"tx_plus_single\";\n\
chrA\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"gene_plus_two\"; transcript_id \"tx_plus_two\";\n\
chrA\ttest\texon\t201\t250\t.\t+\t.\tgene_id \"gene_plus_two\"; transcript_id \"tx_plus_two\";\n\
chrA\ttest\texon\t101\t150\t.\t-\t.\tgene_id \"gene_minus_two\"; transcript_id \"tx_minus_two\";\n\
chrA\ttest\texon\t201\t250\t.\t-\t.\tgene_id \"gene_minus_two\"; transcript_id \"tx_minus_two\";\n\
chrA\ttest\texon\t101\t110\t.\t+\t.\tgene_id \"gene_plus_three\"; transcript_id \"tx_plus_three\";\n\
chrA\ttest\texon\t201\t210\t.\t+\t.\tgene_id \"gene_plus_three\"; transcript_id \"tx_plus_three\";\n\
chrA\ttest\texon\t301\t310\t.\t+\t.\tgene_id \"gene_plus_three\"; transcript_id \"tx_plus_three\";\n\
chrA\ttest\texon\t101\t110\t.\t.\t.\tgene_id \"gene_unknown_two\"; transcript_id \"tx_unknown_two\";\n\
chrA\ttest\texon\t201\t210\t.\t.\t.\tgene_id \"gene_unknown_two\"; transcript_id \"tx_unknown_two\";\n";

        f.write_all(gtf.as_bytes()).unwrap();

        path
    }
    fn write_test_fasta() -> PathBuf {
        let path = tmp_path("genome_fasta");
        let mut seq = vec![b'A'; 1_000];

        seq[100..110].copy_from_slice(b"ACGTACGTAA");
        seq[200..210].copy_from_slice(b"GGGTTTCCCA");

        let mut f = fs::File::create(&path).unwrap();
        writeln!(f, ">chrA").unwrap();
        writeln!(f, "{}", String::from_utf8(seq).unwrap()).unwrap();

        path
    }
    fn dummy_mapper() -> BamTranscriptomeMapper {
        let gtf = write_test_gtf();
        let fasta = write_test_fasta();

        let index = SpliceIndex::from_path(&gtf, 10_000, IdNameKeys::default() ).unwrap();
        let genome = snp_index::Genome::from_fasta(&fasta).unwrap();

        let mapper = BamTranscriptomeMapper::from_index_and_genome(index, genome ).unwrap();
        let _ = fs::remove_file(gtf);
        let _ = fs::remove_file(fasta);

        mapper
    }

    fn tx<'a>(mapper: &'a BamTranscriptomeMapper, name: &str) -> &'a Transcript {
        mapper.index.transcript_by_name(name).unwrap()
    }

    #[test]
    fn plus_single_exon_maps_start_and_cigar_unchanged() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_single");

        let mapped = mapper
            .map_cigar(tx, 10, &[Cigar::Match(20)])
            .unwrap();

        assert_eq!(mapped.start0, 110);
        assert_cigar_eq(&mapped.cigar, &[Cigar::Match(20)]);
    }

    #[test]
    fn plus_two_exon_read_crossing_boundary_gets_refskip() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_two");

        let mapped = mapper
            .map_cigar(tx, 40, &[Cigar::Match(20)])
            .unwrap();

        assert_eq!(mapped.start0, 140);
        assert_cigar_eq(
            &mapped.cigar,
            &[Cigar::Match(10), Cigar::RefSkip(50), Cigar::Match(10)],
        );
    }

    #[test]
    fn plus_two_exon_start_at_second_exon() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_two");

        let mapped = mapper
            .map_cigar(tx, 55, &[Cigar::Match(10)])
            .unwrap();

        assert_eq!(mapped.start0, 205);
        assert_cigar_eq(&mapped.cigar, &[Cigar::Match(10)]);
    }

    #[test]
    fn plus_two_exon_insertion_does_not_consume_transcript() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_two");

        let mapped = mapper
            .map_cigar(
                tx,
                40,
                &[Cigar::Match(5), Cigar::Ins(3), Cigar::Match(10)],
            )
            .unwrap();

        assert_eq!(mapped.start0, 140);
        assert_cigar_eq(
            &mapped.cigar,
            &[
                Cigar::Match(5),
                Cigar::Ins(3),
                Cigar::Match(5),
                Cigar::RefSkip(50),
                Cigar::Match(5),
            ],
        );
    }

    #[test]
    fn plus_two_exon_deletion_consumes_transcript_and_can_cross_boundary() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_two");

        let mapped = mapper
            .map_cigar(tx, 45, &[Cigar::Match(3), Cigar::Del(5), Cigar::Match(4)])
            .unwrap();

        assert_eq!(mapped.start0, 145);
        assert_cigar_eq(
            &mapped.cigar,
            &[
                Cigar::Match(3),
                Cigar::Del(2),
                Cigar::RefSkip(50),
                Cigar::Del(3),
                Cigar::Match(4),
            ],
        );
    }

    #[test]
    fn plus_three_exon_read_crosses_two_introns() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_three");

        let mapped = mapper
            .map_cigar(tx, 5, &[Cigar::Match(25)])
            .unwrap();

        assert_eq!(mapped.start0, 105);
        assert_cigar_eq(
            &mapped.cigar,
            &[
                Cigar::Match(5),
                Cigar::RefSkip(90),
                Cigar::Match(10),
                Cigar::RefSkip(90),
                Cigar::Match(10),
            ],
        );
    }

    #[test]
    fn minus_tx_position_zero_maps_to_last_base_of_last_genomic_exon() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_minus_two");
        let exons = mapper.exons_in_transcript_order(tx);

        assert_eq!(
            mapper
                .tx_pos_to_genome_in_ordered_exons(tx, &exons, 0)
                .unwrap(),
            249
        );
        assert_eq!(
            mapper
                .tx_pos_to_genome_in_ordered_exons(tx, &exons, 49)
                .unwrap(),
            200
        );
        assert_eq!(
            mapper
                .tx_pos_to_genome_in_ordered_exons(tx, &exons, 50)
                .unwrap(),
            149
        );
        assert_eq!(
            mapper
                .tx_pos_to_genome_in_ordered_exons(tx, &exons, 99)
                .unwrap(),
            100
        );
    }

    #[test]
    fn minus_two_exon_single_exon_read_has_leftmost_start_and_reversed_cigar() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_minus_two");

        let mapped = mapper
            .map_cigar(tx, 10, &[Cigar::Match(20)])
            .unwrap();

        assert_eq!(mapped.start0, 220);
        assert_cigar_eq(&mapped.cigar, &[Cigar::Match(20)]);
    }

    #[test]
    fn minus_two_exon_read_crossing_boundary_gets_refskip_and_leftmost_start() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_minus_two");

        let mapped = mapper
            .map_cigar(tx, 40, &[Cigar::Match(20)])
            .unwrap();

        assert_eq!(mapped.start0, 140);
        assert_cigar_eq(
            &mapped.cigar,
            &[Cigar::Match(10), Cigar::RefSkip(50), Cigar::Match(10)],
        );
    }

    #[test]
    fn minus_two_exon_mixed_cigar_is_reversed_for_genomic_bam_order() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_minus_two");

        let mapped = mapper
            .map_cigar(
                tx,
                40,
                &[Cigar::Match(5), Cigar::Ins(3), Cigar::Match(10)],
            )
            .unwrap();

        assert_eq!(mapped.start0, 145);
        assert_cigar_eq(
            &mapped.cigar,
            &[
                Cigar::Match(5),
                Cigar::RefSkip(50),
                Cigar::Match(5),
                Cigar::Ins(3),
                Cigar::Match(5),
            ],
        );
    }

    #[test]
    fn equal_and_diff_are_preserved_and_split() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_two");

        let mapped = mapper
            .map_cigar(tx, 45, &[Cigar::Equal(3), Cigar::Diff(6)])
            .unwrap();

        assert_eq!(mapped.start0, 145);
        assert_cigar_eq(
            &mapped.cigar,
            &[
                Cigar::Equal(3),
                Cigar::Diff(2),
                Cigar::RefSkip(50),
                Cigar::Diff(4),
            ],
        );
    }

    #[test]
    fn adjacent_ops_are_merged_after_mapping() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_single");

        let mapped = mapper
            .map_cigar(tx, 10, &[Cigar::Match(5), Cigar::Match(7)])
            .unwrap();

        assert_eq!(mapped.start0, 110);
        assert_cigar_eq(&mapped.cigar, &[Cigar::Match(12)]);
    }

    #[test]
    fn refskip_in_transcriptome_cigar_is_rejected() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_single");

        let err = mapper
            .map_cigar(tx, 10, &[Cigar::Match(5), Cigar::RefSkip(3)])
            .unwrap_err();

        assert!(
            err.to_string().contains("unexpected RefSkip"),
            "unexpected error: {err:?}"
        );
    }

    #[test]
    fn cigar_extending_beyond_transcript_errors() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_single");

        let err = mapper
            .map_cigar(tx, 95, &[Cigar::Match(10)])
            .unwrap_err();

        assert!(
            err.to_string().contains("extends beyond transcript")
                || err.to_string().contains("outside transcript"),
            "unexpected error: {err:?}"
        );
    }

    #[test]
    fn invalid_start_outside_transcript_errors() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_single");

        let err = mapper
            .map_cigar(tx, 100, &[Cigar::Match(1)])
            .unwrap_err();

        assert!(
            err.to_string().contains("outside transcript"),
            "unexpected error: {err:?}"
        );
    }

    #[test]
    fn soft_and_hard_clips_are_preserved() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_plus_two");

        let mapped = mapper
            .map_cigar(
                tx,
                45,
                &[Cigar::SoftClip(4), Cigar::Match(10), Cigar::HardClip(2)],
            )
            .unwrap();

        assert_eq!(mapped.start0, 145);
        assert_cigar_eq(
            &mapped.cigar,
            &[
                Cigar::SoftClip(4),
                Cigar::Match(5),
                Cigar::RefSkip(50),
                Cigar::Match(5),
                Cigar::HardClip(2),
            ],
        );
    }

    #[test]
    fn complement_base_handles_iupac() {
        let mapper = dummy_mapper();

        assert_eq!(mapper.complement_base(b'A'), b'T');
        assert_eq!(mapper.complement_base(b'C'), b'G');
        assert_eq!(mapper.complement_base(b'G'), b'C');
        assert_eq!(mapper.complement_base(b'T'), b'A');
        assert_eq!(mapper.complement_base(b'R'), b'Y');
        assert_eq!(mapper.complement_base(b'Y'), b'R');
        assert_eq!(mapper.complement_base(b'N'), b'N');
        assert_eq!(mapper.complement_base(b'?'), b'N');
    }


    fn make_record(seq: &[u8], cigar: &[Cigar]) -> bam::Record {
        let mut rec = bam::Record::new();
        let qual = vec![30u8; seq.len()];
        rec.set(b"read1", Some(&CigarString(cigar.to_vec())), seq, &qual);
        rec
    }

    #[test]
    fn unknown_strand_is_inferred_as_plus_from_genome_score() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_unknown_two");

        // Plus transcript order:
        // exon1 chrA:100..110 = ACGTACGTAA
        // exon2 chrA:200..210 = GGGTTTCCCA
        let rec = make_record(b"ACGTACGTAAGGGTTTCCCA", &[Cigar::Match(20)]);

        let strand = mapper
            .infer_record_strand(tx, &rec, 0, &[Cigar::Match(20)])
            .unwrap();

        assert_eq!(strand, Strand::Plus);
    }

    #[test]
    fn unknown_strand_is_inferred_as_minus_from_genome_score() {
        let mapper = dummy_mapper();
        let tx = tx(&mapper, "tx_unknown_two");

        // Minus transcript order is reverse-complement of:
        // exon2 + exon1 = GGGTTTCCCA ACGTACGTAA
        //
        // revcomp("GGGTTTCCCAACGTACGTAA") = TTACGTACGTTGGGAAACCC
        let rec = make_record(b"TTACGTACGTTGGGAAACCC", &[Cigar::Match(20)]);

        let strand = mapper
            .infer_record_strand(tx, &rec, 0, &[Cigar::Match(20)])
            .unwrap();

        assert_eq!(strand, Strand::Minus);
    }
}
