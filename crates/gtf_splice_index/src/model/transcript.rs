use crate::model::types::{MatchClass, MatchHit, MatchOptions, TranscriptId};
use crate::types::{RefBlock, SplicedRead, Strand};
use serde::{Deserialize, Serialize};

const MIN_TRANSCRIPT_END_OVERHANG_BP: u32 = 100;
use int_to_dna::IntToDna;
use int_to_prot::IntToProt;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Transcript {
    pub id: TranscriptId,
    pub gene_id: usize,
    pub names: Vec<String>,
    pub chr_id: usize,
    pub strand: Strand,
    exons: Vec<RefBlock>,
    /// Genomic 0-based half-open CDS span from the annotation.  This stays in
    /// the same coordinate system as the exon blocks; consumers can project it
    /// onto a spliced cDNA only when they actually materialize one.
    cds_start: Option<u32>,
    cds_end: Option<u32>,
    finalized: bool,
}

impl Transcript {
    pub fn new(
        id: TranscriptId,
        gene_id: usize,
        primary_name: impl Into<String>,
        chr_id: usize,
        strand: Strand,
    ) -> Self {
        Self {
            id,
            gene_id,
            names: vec![primary_name.into()],
            chr_id,
            strand,
            exons: Vec::new(),
            cds_start: None,
            cds_end: None,
            finalized: false,
        }
    }

    pub fn add_name(&mut self, name: &str) {
        let name = name.trim();
        if name.is_empty() {
            return;
        }
        if !self.names.iter().any(|n| n == name) {
            self.names.push(name.to_string());
        }
    }

    pub fn primary_name(&self) -> Option<&str> {
        self.names.first().map(|s| s.as_str())
    }

    pub fn add_exon(&mut self, block: RefBlock) {
        self.exons.push(block);
        self.finalized = false;
    }

    pub fn exons(&self) -> &[RefBlock] {
        &self.exons
    }

    pub fn add_cds(&mut self, block: RefBlock) {
        self.cds_start = Some(self.cds_start.map_or(block.start, |x| x.min(block.start)));
        self.cds_end = Some(self.cds_end.map_or(block.end, |x| x.max(block.end)));
    }

    /// Genomic 0-based half-open coding span.  The bounds deliberately remain
    /// genomic; strand-aware projection onto spliced transcript coordinates is
    /// a downstream operation.
    pub fn cds_span(&self) -> Option<(u32, u32)> {
        Some((self.cds_start?, self.cds_end?))
    }

    /// Length of the spliced transcript in bases.
    pub fn transcript_len(&self) -> usize {
        self.exons.iter().map(|b| b.len() as usize).sum()
    }

    /// Map a genomic base coordinate onto the spliced transcript.
    /// Transcript coordinates are always 5' -> 3', so minus-strand
    /// transcripts run opposite to genomic coordinate order.
    pub fn transcript_position_of_genomic(&self, genomic_pos: u32) -> Option<usize> {
        let mut offset = 0usize;
        match self.strand {
            Strand::Plus | Strand::Unknown => {
                for exon in &self.exons {
                    if genomic_pos >= exon.start && genomic_pos < exon.end {
                        return Some(offset + (genomic_pos - exon.start) as usize);
                    }
                    offset += exon.len() as usize;
                }
            }
            Strand::Minus => {
                for exon in self.exons.iter().rev() {
                    if genomic_pos >= exon.start && genomic_pos < exon.end {
                        return Some(offset + (exon.end - 1 - genomic_pos) as usize);
                    }
                    offset += exon.len() as usize;
                }
            }
        }
        None
    }

    /// Map a spliced transcript base coordinate back to the genome.
    pub fn genomic_position_of_transcript(&self, transcript_pos: usize) -> Option<u32> {
        let mut offset = 0usize;
        match self.strand {
            Strand::Plus | Strand::Unknown => {
                for exon in &self.exons {
                    let len = exon.len() as usize;
                    if transcript_pos < offset + len {
                        return Some(exon.start + (transcript_pos - offset) as u32);
                    }
                    offset += len;
                }
            }
            Strand::Minus => {
                for exon in self.exons.iter().rev() {
                    let len = exon.len() as usize;
                    if transcript_pos < offset + len {
                        return Some(exon.end - 1 - (transcript_pos - offset) as u32);
                    }
                    offset += len;
                }
            }
        }
        None
    }

    /// CDS bounds in spliced transcript coordinates, 0-based and half-open.
    pub fn cds_transcript_span(&self) -> Option<(usize, usize)> {
        let (start, end) = self.cds_span()?;
        if start >= end {
            return None;
        }
        let a = self.transcript_position_of_genomic(start)?;
        let b = self.transcript_position_of_genomic(end - 1)?;
        Some((a.min(b), a.max(b) + 1))
    }

    /// Materialize the spliced cDNA from one chromosome/reference sequence.
    /// `chromosome` must use the same 0-based genomic coordinate system as
    /// this transcript's exon blocks.
    pub fn cdna(&self, chromosome: &[u8]) -> Result<IntToDna, String> {
        let mut seq = Vec::with_capacity(self.transcript_len());
        match self.strand {
            Strand::Plus | Strand::Unknown => {
                for exon in &self.exons {
                    append_reference_block(&mut seq, chromosome, *exon)?;
                }
            }
            Strand::Minus => {
                for exon in self.exons.iter().rev() {
                    let start = exon.start as usize;
                    let end = exon.end as usize;
                    let bases = chromosome.get(start..end).ok_or_else(|| {
                        format!(
                            "exon {}..{} exceeds reference length {}",
                            start,
                            end,
                            chromosome.len()
                        )
                    })?;
                    seq.extend(bases.iter().rev().map(|&b| complement(b)));
                }
            }
        }
        IntToDna::try_new(seq)
    }

    /// Materialize only the protein-coding portion of the spliced transcript.
    pub fn coding_dna(&self, chromosome: &[u8]) -> Result<Option<IntToDna>, String> {
        let Some((start, end)) = self.cds_transcript_span() else {
            return Ok(None);
        };
        let cdna = self.cdna(chromosome)?;
        let seq = cdna.to_string(cdna.size);
        Ok(Some(IntToDna::try_new(&seq.as_bytes()[start..end])?))
    }

    /// Materialize and translate this transcript's annotated CDS.
    pub fn protein(&self, chromosome: &[u8]) -> Result<Option<IntToProt>, String> {
        Ok(self.coding_dna(chromosome)?.map(|cds| cds.translate()))
    }

    /// sorts the transcripts exons and returns (total start: u32, total end: u32)
    pub fn finalize(&mut self) -> (u32, u32) {
        if self.exons.is_empty() {
            self.finalized = true;
            return (0, 0);
        }

        self.exons.sort_by_key(|b| (b.start, b.end));

        let mut merged: Vec<RefBlock> = Vec::with_capacity(self.exons.len());
        let mut cur = self.exons[0];

        for &b in &self.exons[1..] {
            if b.start <= cur.end {
                cur.end = cur.end.max(b.end);
            } else {
                merged.push(cur);
                cur = b;
            }
        }
        merged.push(cur);

        self.exons = merged;
        self.finalized = true;

        let start = self.exons.first().unwrap().start;
        let end = self.exons.last().unwrap().end;

        (start, end)
    }

    pub fn span(&self) -> Option<(u32, u32)> {
        if self.exons.is_empty() {
            return None;
        }
        Some((
            self.exons.first().unwrap().start,
            self.exons.last().unwrap().end,
        ))
    }

    pub fn junctions(&self) -> Vec<(u32, u32)> {
        RefBlock::junctions_from_blocks(&self.exons, 0)
    }

    /// For every observed read junction that is not an exact annotated junction,
    /// return the signed donor and acceptor displacement from the nearest
    /// annotated junction of this transcript. Coordinates are the native
    /// 0-based half-open block boundaries used by both the BAM-derived read and
    /// the splice index. Positive means the observed boundary is at a larger
    /// genomic coordinate than the annotation.
    pub fn junction_mismatch_offsets(
        &self,
        read: &SplicedRead,
        allowed_gap_size: u32,
    ) -> Vec<(i32, i32)> {
        read.assert_finalized();
        let read_junctions = RefBlock::junctions_from_blocks(&read.blocks, allowed_gap_size);
        let tx_junctions = self.junctions();
        if tx_junctions.is_empty() {
            return Vec::new();
        }

        read_junctions
            .into_iter()
            .filter(|junction| !tx_junctions.contains(junction))
            .filter_map(|(read_donor, read_acceptor)| {
                tx_junctions
                    .iter()
                    .min_by_key(|&&(tx_donor, tx_acceptor)| {
                        read_donor.abs_diff(tx_donor) as u64
                            + read_acceptor.abs_diff(tx_acceptor) as u64
                    })
                    .map(|&(tx_donor, tx_acceptor)| {
                        (
                            read_donor as i64 - tx_donor as i64,
                            read_acceptor as i64 - tx_acceptor as i64,
                        )
                    })
                    .and_then(|(donor, acceptor)| {
                        Some((i32::try_from(donor).ok()?, i32::try_from(acceptor).ok()?))
                    })
            })
            .collect()
    }

    pub fn match_spliced_read(&self, read: &SplicedRead, opts: MatchOptions) -> MatchHit {
        read.assert_finalized();
        self.match_read_blocks(read.chr_id, read.strand, &read.blocks, opts)
    }
    /*
    fn overlaps(a0: u32, a1: u32, b0: u32, b1: u32) -> bool {
        a0 < b1 && b0 < a1
    }

    fn introns(&self) -> Vec<RefBlock> {
        let mut out = Vec::new();
        if self.exons.len() < 2 {
            return out;
        }
        for w in self.exons.windows(2) {
            let a = w[0];
            let b = w[1];
            // intron is [a.end, b.start)
            if a.end < b.start {
                out.push(RefBlock::new(a.end, b.start));
            }
        }
        out
    }
    */

    fn match_read_blocks(
        &self,
        read_chr_id: usize,
        read_strand: Strand,
        read_blocks: &[RefBlock],
        opts: MatchOptions,
    ) -> MatchHit {
        if !self.finalized {
            panic!("Transcript::match_read_blocks called before finalize()");
        }

        // Overhangs are reported for diagnostics and ranking/debugging.
        // They are initialized to zero and filled once chromosome/span compatibility
        // has been established.
        let mut over5 = 0u32;
        let mut over3 = 0u32;

        // ------------------------------------------------------------
        // 1) Basic validity / chromosome checks
        // ------------------------------------------------------------
        if read_blocks.is_empty() {
            return MatchHit::new(MatchClass::NoOverlap, over5, over3);
        }

        if self.chr_id != read_chr_id {
            return MatchHit::new(MatchClass::NoOverlap, over5, over3);
        }

        // Read span is the union span from first block start to last block end.
        // Assumption: read_blocks are sorted by genomic coordinate.
        let r0 = read_blocks[0].start;
        let r1 = read_blocks[read_blocks.len() - 1].end;

        let Some((t0, t1)) = self.span() else {
            return MatchHit::new(MatchClass::NoOverlap, over5, over3);
        };

        let read_span = RefBlock { start: r0, end: r1 };
        let tx_span = RefBlock { start: t0, end: t1 };

        if !read_span.overlaps(tx_span) {
            return MatchHit::new(MatchClass::NoOverlap, over5, over3);
        }

        // From here on the read overlaps the transcript span, so overhangs are
        // meaningful to report.
        if let Some((o5, o3)) = self.compute_overhangs_strand_aware(read_blocks) {
            over5 = o5;
            over3 = o3;
        }

        // ------------------------------------------------------------
        // 2) Optional strand check
        // ------------------------------------------------------------
        if opts.require_strand && !self.strand.is_compatible_with(read_strand) {
            return MatchHit::new(MatchClass::StrandMismatch, over5, over3);
        }

        // ------------------------------------------------------------
        // 3) Transcript boundary overhang check
        // ------------------------------------------------------------
        // Transcript ends, especially annotated 3-prime/poly(A) ends, are not
        // exact biological boundaries. Never make the matcher stricter than
        // 100 bp even when an older caller explicitly passes zero. Larger user
        // tolerances remain honored.
        let max_5p_overhang_bp = opts.max_5p_overhang_bp.max(MIN_TRANSCRIPT_END_OVERHANG_BP);
        let max_3p_overhang_bp = opts.max_3p_overhang_bp.max(MIN_TRANSCRIPT_END_OVERHANG_BP);
        if over5 > max_5p_overhang_bp || over3 > max_3p_overhang_bp {
            return MatchHit::new(MatchClass::OverhangTooLarge, over5, over3);
        }

        // ------------------------------------------------------------
        // 4) Exon-fit guard
        // ------------------------------------------------------------
        // This is the important biological routing guard.
        //
        // Before we even look at junction compatibility, every read block must be
        // compatible with transcript exons:
        //
        // - a purely intronic single-block read must fail here
        // - a read crossing exon/intron sequence incorrectly must fail here
        // - a valid exon-contained single-block read may pass here
        // - a valid spliced read with blocks fitting exons may pass here
        //
        // Without this guard, an unspliced read with no junctions can be incorrectly
        // classified as Compatible because an empty iterator satisfies `.all(...)`.
        if !self.blocks_fit_exons_allowing_end_overhang(read_blocks, 10) {
            // Failing exon compatibility is NOT sufficient evidence for an
            // intronic molecule. Only route to Intronic when aligned sequence
            // positively overlaps an annotated intron of this transcript.
            let class = if self.blocks_overlap_introns(read_blocks) {
                MatchClass::Intronic
            } else {
                MatchClass::Incompatible
            };
            return MatchHit::new(class, over5, over3);
        }

        // ------------------------------------------------------------
        // 5) Junction-chain classification
        // ------------------------------------------------------------
        let read_junctions =
            RefBlock::junctions_from_blocks(read_blocks, opts.allowed_intronic_gap_size);
        let tx_junctions = self.junctions();

        // If exact mode is requested, a read must reproduce the complete transcript
        // junction chain exactly.
        //
        // Important bug fix:
        // Do NOT allow a single-block read with an empty junction list to be
        // ExactJunctionChain for a single-exon transcript merely because
        // `read_junctions == tx_junctions == []`.
        //
        // A single-block exon-compatible read is Compatible, not ExactJunctionChain.
        if opts.require_exact_junction_chain {
            let class = if !read_junctions.is_empty() && read_junctions == tx_junctions {
                MatchClass::ExactJunctionChain
            } else if read_junctions.is_empty() {
                MatchClass::Compatible
            } else {
                MatchClass::JunctionMismatch
            };

            return MatchHit::new(class, over5, over3);
        }

        // Non-exact mode:
        //
        // - Empty read_junctions:
        //   The read has no splice junctions. Since it passed exon-fit above, it is
        //   exon-compatible, but it does not prove the transcript junction chain.
        //   Therefore it is Compatible, not ExactJunctionChain.
        //
        // - Non-empty read_junctions:
        //   Every read junction must exist in the transcript junction chain.
        //   If the full chain is identical, it is ExactJunctionChain.
        //   If it is a proper subset, it is Compatible.
        //   Otherwise it is JunctionMismatch.
        let class = if read_junctions.is_empty() {
            MatchClass::Compatible
        } else if read_junctions.iter().all(|j| tx_junctions.contains(j)) {
            if read_junctions == tx_junctions {
                MatchClass::ExactJunctionChain
            } else {
                MatchClass::Compatible
            }
        } else {
            MatchClass::JunctionMismatch
        };

        MatchHit::new(class, over5, over3)
    }

    fn compute_overhangs_strand_aware(&self, read_blocks: &[RefBlock]) -> Option<(u32, u32)> {
        let (t0, t1) = self.span()?;
        let r0 = read_blocks.first()?.start;
        let r1 = read_blocks.last()?.end;

        let left_over = t0.saturating_sub(r0);
        let right_over = r1.saturating_sub(t1);

        match self.strand {
            Strand::Plus => Some((left_over, right_over)),
            Strand::Minus => Some((right_over, left_over)),
            Strand::Unknown => Some((left_over, right_over)),
        }
    }

    fn blocks_overlap_introns(&self, read_blocks: &[RefBlock]) -> bool {
        if self.exons.len() < 2 {
            return false;
        }

        self.exons.windows(2).any(|pair| {
            let left = pair[0];
            let right = pair[1];
            if left.end >= right.start {
                return false;
            }
            let intron = RefBlock::new(left.end, right.start);
            read_blocks.iter().any(|block| block.overlaps(intron))
        })
    }

    fn blocks_fit_exons_allowing_end_overhang(
        &self,
        read_blocks: &[RefBlock],
        max_in_exon_gap_bp: u32,
    ) -> bool {
        let mut exon_idx = 0usize;
        let num_blocks = read_blocks.len();

        // Track previous block and which exon it matched,
        // so we can validate gaps inside an exon.
        let mut prev_block_end: Option<u32> = None;
        let mut prev_exon_idx: Option<usize> = None;

        for (block_idx, &read_block) in read_blocks.iter().enumerate() {
            // Move exon pointer until we reach an exon whose end is beyond the block start.
            while exon_idx < self.exons.len() && self.exons[exon_idx].end <= read_block.start {
                exon_idx += 1;
            }
            if exon_idx == self.exons.len() {
                return false;
            }

            let exon = self.exons[exon_idx];

            // Block must overlap the exon it is assigned to.
            if !read_block.overlaps(exon) {
                return false;
            }

            let is_first_block = block_idx == 0;
            let is_last_block = block_idx + 1 == num_blocks;

            // Sequence outside an exon is allowed only at the OUTER transcript
            // boundaries. This is genomic-coordinate based, so it works for both
            // strands and also for single-exon transcripts. The actual tolerance
            // was already enforced by compute_overhangs_strand_aware() above.
            let is_first_exon = exon_idx == 0;
            let is_last_exon = exon_idx + 1 == self.exons.len();
            let allow_left_overhang = is_first_block && is_first_exon;
            let allow_right_overhang = is_last_block && is_last_exon;

            if read_block.start < exon.start && !allow_left_overhang {
                return false;
            }
            if read_block.end > exon.end && !allow_right_overhang {
                return false;
            }

            // -----------------------------
            // NEW: allow small gaps inside the same exon
            // -----------------------------
            if let (Some(prev_end), Some(prev_ei)) = (prev_block_end, prev_exon_idx) {
                if read_block.start > prev_end {
                    let gap = read_block.start - prev_end;

                    // If two consecutive blocks map to the SAME exon, the gap must be small
                    // and fully internal to that exon. Otherwise it's a structural mismatch.
                    if prev_ei == exon_idx {
                        if gap > max_in_exon_gap_bp {
                            return false;
                        }

                        // Gap endpoints must lie within the exon bounds.
                        // (We allow equality at exon end/start since coordinates are half-open.)
                        let prev_end_in_exon = prev_end >= exon.start && prev_end <= exon.end;
                        let next_start_in_exon =
                            read_block.start >= exon.start && read_block.start <= exon.end;

                        if !prev_end_in_exon || !next_start_in_exon {
                            return false;
                        }
                    }
                }
            }

            prev_block_end = Some(read_block.end);
            prev_exon_idx = Some(exon_idx);
        }

        true
    }
}

fn append_reference_block(
    out: &mut Vec<u8>,
    chromosome: &[u8],
    block: RefBlock,
) -> Result<(), String> {
    let start = block.start as usize;
    let end = block.end as usize;
    let bases = chromosome.get(start..end).ok_or_else(|| {
        format!(
            "exon {}..{} exceeds reference length {}",
            start,
            end,
            chromosome.len()
        )
    })?;
    out.extend_from_slice(bases);
    Ok(())
}

#[inline]
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::MatchOptions;

    // ---------- small helpers ----------
    fn opts() -> MatchOptions {
        MatchOptions {
            max_5p_overhang_bp: 15,
            max_3p_overhang_bp: 15,
            require_strand: false,
            require_exact_junction_chain: false,
            ..Default::default()
        }
    }

    fn tx_two_exons(chr_id: usize, strand: Strand) -> Transcript {
        let mut tx = Transcript::new(0, 0, "T", chr_id, strand);
        tx.add_exon(RefBlock::new(100, 150));
        tx.add_exon(RefBlock::new(200, 250));
        tx.finalize();
        tx
    }

    fn tx_three_exons(chr_id: usize, strand: Strand) -> Transcript {
        let mut tx = Transcript::new(0, 0, "T", chr_id, strand);
        tx.add_exon(RefBlock::new(100, 150));
        tx.add_exon(RefBlock::new(200, 250));
        tx.add_exon(RefBlock::new(300, 350));
        tx.finalize();
        tx
    }

    fn read(chr_id: usize, strand: Strand, blocks: Vec<RefBlock>) -> SplicedRead {
        let mut r = SplicedRead::new(chr_id, strand, blocks);
        r.finalize();
        r
    }

    // ---------- your existing tests ----------
    #[test]
    fn transcript_names_dedup() {
        let mut t = Transcript::new(0, 0, "T1", 0, Strand::Plus);
        t.add_name("T1");
        t.add_name("ENST0000");
        t.add_name("ENST0000");
        assert_eq!(t.names, vec!["T1".to_string(), "ENST0000".to_string()]);
        assert_eq!(t.primary_name(), Some("T1"));
    }

    #[test]
    fn exact_chain_and_overhang_reported() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(90, 150), RefBlock::new(200, 260)],
        );

        let hit = tx.match_spliced_read(
            &read,
            MatchOptions {
                max_5p_overhang_bp: 15,
                max_3p_overhang_bp: 15,
                ..Default::default()
            },
        );

        assert_eq!(hit.class, MatchClass::ExactJunctionChain);
        assert_eq!(hit.overhang_5p_bp, 10);
        assert_eq!(hit.overhang_3p_bp, 10);
    }

    #[test]
    fn exact_chain_with_sequenceing_errors() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(
            1,
            Strand::Plus,
            vec![
                RefBlock::new(90, 150),
                RefBlock::new(200, 229),
                RefBlock::new(230, 260),
            ],
        );

        let hit = tx.match_spliced_read(
            &read,
            MatchOptions {
                max_5p_overhang_bp: 15,
                max_3p_overhang_bp: 15,
                allowed_intronic_gap_size: 10,
                ..Default::default()
            },
        );

        assert_eq!(hit.class, MatchClass::ExactJunctionChain);
        assert_eq!(hit.overhang_5p_bp, 10);
        assert_eq!(hit.overhang_3p_bp, 10);
    }

    // ---------- new comprehensive coverage tests ----------

    #[test]
    fn matchclass_nooverlap_same_chr_far_away() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(1, Strand::Plus, vec![RefBlock::new(1000, 1050)]);

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::NoOverlap);
        assert_eq!(hit.overhang_5p_bp, 0);
        assert_eq!(hit.overhang_3p_bp, 0);
    }

    #[test]
    fn matchclass_strand_mismatch_when_required() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(
            1,
            Strand::Minus,
            vec![RefBlock::new(100, 150), RefBlock::new(200, 250)],
        );

        let mut o = opts();
        o.require_strand = true;

        let hit = tx.match_spliced_read(&read, o);

        assert_eq!(hit.class, MatchClass::StrandMismatch);
        // overhangs should still be reported (typically 0 here)
        assert_eq!(hit.overhang_5p_bp, 0);
        assert_eq!(hit.overhang_3p_bp, 0);
    }

    #[test]
    fn matchclass_overhang_too_large() {
        let tx = tx_two_exons(1, Strand::Plus);

        // 5' overhang = 30 (start 70 vs tx start 100)
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(70, 150), RefBlock::new(200, 250)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::OverhangTooLarge);
        assert_eq!(hit.overhang_5p_bp, 30);
        assert_eq!(hit.overhang_3p_bp, 0);
    }

    #[test]
    fn matchclass_intronic_due_to_intron_overlap() {
        let tx = tx_two_exons(1, Strand::Plus);

        // This second block overlaps the intron region [150,200) (and also a bit into exon2 start),
        // which should fail blocks_fit_exons_allowing_end_overhang(...) and be classified as Intronic.
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(170, 210)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::Intronic);
    }

    #[test]
    fn matchclass_junction_mismatch_skipping_middle_exon() {
        let tx = tx_three_exons(1, Strand::Plus);

        // Read uses exon1 + exon3, skipping exon2.
        // Blocks still fit exons, but junction (150 -> 300) is not in tx junction set.
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(300, 350)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::JunctionMismatch);
    }

    #[test]
    fn matchclass_compatible_when_junctions_are_subset() {
        let tx = tx_three_exons(1, Strand::Plus);

        // Read uses exon1 + exon2 only (junction subset, not equal to tx's full chain).
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(200, 250)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::Compatible);
    }

    #[test]
    fn require_exact_junction_chain_rejects_subset_as_junction_mismatch() {
        let tx = tx_three_exons(1, Strand::Plus);

        // Same read as Compatible case, but now exact is required.
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(200, 250)],
        );

        let mut o = opts();
        o.require_exact_junction_chain = true;

        let hit = tx.match_spliced_read(&read, o);

        assert_eq!(hit.class, MatchClass::JunctionMismatch);
    }

    #[test]
    fn single_block_exonic_and_intronic_reads_are_not_both_compatible() {
        let mut tx = Transcript {
            id: 0,
            gene_id: 0,
            names: vec!["tx1".to_string()],
            chr_id: 0,
            strand: Strand::Plus,
            exons: vec![
                RefBlock {
                    start: 100,
                    end: 150,
                },
                RefBlock {
                    start: 250,
                    end: 300,
                },
            ],
            cds_start: None,
            cds_end: None,
            finalized: false,
        };

        tx.finalize();

        let opts = MatchOptions::default();

        let exon_internal = vec![RefBlock {
            start: 110,
            end: 130,
        }];
        let intron_internal = vec![RefBlock {
            start: 180,
            end: 220,
        }];

        let exon_hit = tx.match_read_blocks(0, Strand::Plus, &exon_internal, opts);
        let intron_hit = tx.match_read_blocks(0, Strand::Plus, &intron_internal, opts);

        assert_eq!(
            exon_hit.class,
            MatchClass::Compatible,
            "single-block read fully inside exon should be exon-compatible"
        );

        assert_eq!(
            intron_hit.class,
            MatchClass::Intronic,
            "single-block read fully inside intron must not become Compatible"
        );
    }

    #[test]
    fn plus_strand_materializes_cdna_cds_and_protein() {
        let mut tx = Transcript::new(1, 1, "tx", 0, Strand::Plus);
        tx.add_exon(RefBlock::new(0, 6));
        tx.add_exon(RefBlock::new(10, 19));
        tx.add_cds(RefBlock::new(0, 6));
        tx.add_cds(RefBlock::new(10, 19));
        tx.finalize();

        let chromosome = b"ATGGCTNNNNGAATTTTAA";
        assert_eq!(
            tx.cdna(chromosome).unwrap().to_string(15),
            "ATGGCTGAATTTTAA"
        );
        assert_eq!(tx.cds_transcript_span(), Some((0, 15)));
        assert_eq!(
            tx.coding_dna(chromosome).unwrap().unwrap().to_string(15),
            "ATGGCTGAATTTTAA"
        );
        assert_eq!(
            tx.protein(chromosome).unwrap().unwrap().to_string(),
            "MAEF*"
        );
        assert_eq!(tx.transcript_position_of_genomic(10), Some(6));
        assert_eq!(tx.genomic_position_of_transcript(6), Some(10));
    }

    #[test]
    fn minus_strand_materializes_cdna_cds_and_protein_in_transcript_orientation() {
        let mut tx = Transcript::new(1, 1, "tx", 0, Strand::Minus);
        tx.add_exon(RefBlock::new(0, 9));
        tx.add_exon(RefBlock::new(12, 18));
        tx.add_cds(RefBlock::new(0, 9));
        tx.add_cds(RefBlock::new(12, 18));
        tx.finalize();

        // reverse-complement(exon 12..18) + reverse-complement(exon 0..9)
        // = ATGGCT + GAATTTTAA
        let chromosome = b"TTAAAATTCNNNAGCCAT";
        assert_eq!(
            tx.cdna(chromosome).unwrap().to_string(15),
            "ATGGCTGAATTTTAA"
        );
        assert_eq!(tx.cds_transcript_span(), Some((0, 15)));
        assert_eq!(
            tx.protein(chromosome).unwrap().unwrap().to_string(),
            "MAEF*"
        );
        assert_eq!(tx.transcript_position_of_genomic(17), Some(0));
        assert_eq!(tx.genomic_position_of_transcript(0), Some(17));
    }

    #[test]
    fn utr_is_removed_before_translation() {
        let mut tx = Transcript::new(1, 1, "tx", 0, Strand::Plus);
        tx.add_exon(RefBlock::new(0, 9));
        tx.add_exon(RefBlock::new(12, 24));
        tx.add_cds(RefBlock::new(3, 9));
        tx.add_cds(RefBlock::new(12, 21));
        tx.finalize();

        let chromosome = b"CCCATGGCTNNNGAATTTTAAGGG";
        assert_eq!(
            tx.cdna(chromosome).unwrap().to_string(21),
            "CCCATGGCTGAATTTTAAGGG"
        );
        assert_eq!(tx.cds_transcript_span(), Some((3, 18)));
        assert_eq!(
            tx.protein(chromosome).unwrap().unwrap().to_string(),
            "MAEF*"
        );
    }

    #[test]
    fn exon_fit_failure_without_intron_overlap_is_not_intronic() {
        let tx = tx_two_exons(1, Strand::Plus);
        let read = read(1, Strand::Plus, vec![RefBlock::new(90, 160)]);
        let mut o = opts();
        o.max_5p_overhang_bp = 100;
        o.max_3p_overhang_bp = 100;
        let hit = tx.match_spliced_read(&read, o);
        assert_eq!(
            hit.class,
            MatchClass::Intronic,
            "block entering the annotated intron is positive intronic evidence"
        );

        let mut single = Transcript::new(0, 0, "single", 1, Strand::Plus);
        single.add_exon(RefBlock::new(100, 150));
        single.finalize();
        let read = read(1, Strand::Plus, vec![RefBlock::new(90, 160)]);
        let hit = single.match_spliced_read(&read, o);
        assert_eq!(
            hit.class,
            MatchClass::Compatible,
            "terminal overhang within tolerance remains exonic-compatible"
        );
    }

    #[test]
    fn default_match_options_allow_100bp_transcript_end_overhang() {
        let mut tx = Transcript::new(0, 0, "single", 1, Strand::Plus);
        tx.add_exon(RefBlock::new(100, 150));
        tx.finalize();

        let left = read(1, Strand::Plus, vec![RefBlock::new(1, 150)]);
        let right = read(1, Strand::Plus, vec![RefBlock::new(100, 249)]);
        assert_eq!(
            tx.match_spliced_read(&left, MatchOptions::default()).class,
            MatchClass::Compatible
        );
        assert_eq!(
            tx.match_spliced_read(&right, MatchOptions::default()).class,
            MatchClass::Compatible
        );

        let too_far = read(1, Strand::Plus, vec![RefBlock::new(100, 251)]);
        assert_eq!(
            tx.match_spliced_read(&too_far, MatchOptions::default())
                .class,
            MatchClass::OverhangTooLarge
        );

        let legacy_zero = MatchOptions {
            max_5p_overhang_bp: 0,
            max_3p_overhang_bp: 0,
            ..MatchOptions::default()
        };
        assert_eq!(
            tx.match_spliced_read(&right, legacy_zero).class,
            MatchClass::Compatible
        );
    }

    #[test]
    fn junction_mismatch_offsets_report_signed_one_base_shifts() {
        let mut tx = Transcript::new(0, 0, "T", 0, Strand::Plus);
        tx.add_exon(RefBlock::new(100, 150));
        tx.add_exon(RefBlock::new(200, 250));
        tx.finalize();

        let mut exact = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(200, 240)],
        );
        exact.finalize();
        assert!(tx.junction_mismatch_offsets(&exact, 0).is_empty());

        let mut donor_plus_one = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(110, 151), RefBlock::new(200, 240)],
        );
        donor_plus_one.finalize();
        assert_eq!(
            tx.junction_mismatch_offsets(&donor_plus_one, 0),
            vec![(1, 0)]
        );

        let mut acceptor_minus_one = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(199, 240)],
        );
        acceptor_minus_one.finalize();
        assert_eq!(
            tx.junction_mismatch_offsets(&acceptor_minus_one, 0),
            vec![(0, -1)]
        );
    }

    #[test]
    fn exact_junction_coordinates_are_zero_based_half_open() {
        let mut tx = Transcript::new(0, 0, "T", 0, Strand::Plus);
        tx.add_exon(RefBlock::new(100, 150));
        tx.add_exon(RefBlock::new(200, 250));
        tx.finalize();

        assert_eq!(tx.junctions(), vec![(150, 200)]);

        let mut read = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(120, 150), RefBlock::new(200, 230)],
        );
        read.finalize();
        assert_eq!(read.junctions(), vec![(150, 200)]);
        assert_eq!(
            tx.match_spliced_read(&read, MatchOptions::default()).class,
            MatchClass::ExactJunctionChain
        );
    }
}
