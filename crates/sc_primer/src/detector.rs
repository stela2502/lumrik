use crate::chemistry::Chemistry;
use crate::error::{PrimerError, PrimerResult};
use crate::grammar::{Grammar, GrammarOp};
use crate::model::{
    BdPrimerDiagnostics, Orientation, PrimerAttempt, PrimerMatch, PrimerMatchDiagnostics,
    PrimerSegmentAttempt,
};

use crate::single_cell_systems::*;

use int_to_str::IntToStr;
use onehot_dna::OneHotSequence;
use read_tag_table::ReadTagRecord;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct PrimerDetector {
    grammar: Grammar,
    pub single_cell_system: Option<SingleCellSystem>,
    detect_reverse_complement: bool,
    // source_cell -> (target_cell, count)
    primer_translation: HashMap<u64, (u64, usize)>,
    umi_translation: HashMap<u64, usize>,
    alternates: Vec<PrimerDetector>,
}

impl PrimerDetector {
    pub fn from_chemistry(chemistry: Chemistry) -> PrimerResult<Self> {
        Self::from_grammar(chemistry.grammar()?)
    }

    /// Compile multiple independent chemistry grammars into one detector.
    ///
    /// Detection order is optimized independently of CLI order. Grammars with
    /// a leading FIXED anchor are preferred because they can reject unrelated
    /// reads before invoking more expensive barcode/cassette matching.
    pub fn from_chemistries<I>(chemistries: I) -> PrimerResult<Self>
    where
        I: IntoIterator<Item = Chemistry>,
    {
        let grammars = chemistries
            .into_iter()
            .map(Chemistry::grammar)
            .collect::<PrimerResult<Vec<_>>>()?;
        Self::from_grammars(grammars)
    }

    pub fn from_grammars(mut grammars: Vec<Grammar>) -> PrimerResult<Self> {
        if grammars.is_empty() {
            return Err(PrimerError::invalid_grammar(
                "at least one primer grammar is required",
            ));
        }

        grammars.sort_by_key(Self::grammar_rejection_cost);
        let mut compiled = grammars
            .into_iter()
            .map(Self::from_grammar)
            .collect::<PrimerResult<Vec<_>>>()?;
        let mut primary = compiled.remove(0);
        primary.alternates = compiled;
        Ok(primary)
    }

    fn grammar_rejection_cost(grammar: &Grammar) -> (u8, usize, usize) {
        match grammar.ops.first() {
            Some(GrammarOp::Fixed { seq, mismatches }) => {
                (0, *mismatches, usize::MAX - seq.len())
            }
            _ if grammar.anchor_search().is_some() => (1, 0, 0),
            Some(GrammarOp::Search { .. }) => (3, 0, 0),
            _ => (2, 0, 0),
        }
    }

    pub fn from_grammar(grammar: Grammar) -> PrimerResult<Self> {
        let single_cell_system = grammar.system();

        Ok(Self {
            grammar,
            single_cell_system,
            detect_reverse_complement: true,
            primer_translation: HashMap::new(),
            umi_translation: HashMap::new(),
            alternates: Vec::new(),
        })
    }

    pub fn with_reverse_complement_detection(mut self, enabled: bool) -> Self {
        self.detect_reverse_complement = enabled;
        for detector in &mut self.alternates {
            detector.detect_reverse_complement = enabled;
        }
        self
    }

    /// Replace/install the cell whitelist from a line-delimited text file.
    ///
    /// TENX_CELL grammars retain their chemistry coordinates but use the
    /// supplied 16 bp whitelist. A generic grammar must contain exactly one
    /// CELL:N operation; its N (1..=32) selects the const-generic matcher
    /// behind the runtime CLI facade. BD uses three independent 9 bp block
    /// whitelists and therefore requires its built-in component tables.
    pub fn with_whitelist_path<P: AsRef<Path>>(
        mut self,
        path: P,
        max_mismatches: u32,
    ) -> PrimerResult<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| {
            PrimerError::invalid_grammar(format!(
                "failed to read whitelist '{}': {e}",
                path.display()
            ))
        })?;

        match &self.single_cell_system {
            Some(SingleCellSystem::Tenx(tenx)) => {
                let version = tenx.version();
                let replacement =
                    TenxWhitelist::from_text_with_mismatches(version, &text, max_mismatches)
                        .map_err(PrimerError::invalid_grammar)?;
                self.single_cell_system = Some(SingleCellSystem::Tenx(replacement));
            }
            Some(SingleCellSystem::Rhapsody(_)) => {
                return Err(PrimerError::invalid_grammar(
                    "--whitelist cannot replace BD's three independent C1/C2/C3 block whitelists",
                ));
            }
            Some(SingleCellSystem::Whitelist(_)) => unreachable!("detector is newly constructed"),
            None => {
                let cell_ops = self
                    .grammar
                    .ops
                    .iter()
                    .filter_map(|op| match op {
                        GrammarOp::Cell { len } => Some(*len),
                        _ => None,
                    })
                    .collect::<Vec<_>>();
                if cell_ops.len() != 1 {
                    return Err(PrimerError::invalid_grammar(
                        "--whitelist with a custom grammar requires exactly one CELL:N operation",
                    ));
                }
                let whitelist = crate::whitelist_hash::RuntimeWhitelistHash::from_text(
                    cell_ops[0],
                    &text,
                    max_mismatches,
                )
                .map_err(PrimerError::invalid_grammar)?;
                self.single_cell_system = Some(SingleCellSystem::Whitelist(whitelist));
            }
        }

        Ok(self)
    }

    /// Length of the normalized cell barcode returned in `PrimerMatch::cell_seq`.
    ///
    /// BD Rhapsody cassettes contain linker/spacing sequence in the read, but the
    /// normalized cell barcode is the three corrected 9 bp blocks (27 bp).
    /// `grammar.cell_len()` describes the full BD cassette representation used by
    /// synthesis and must therefore not be used for matrix barcode export.
    pub fn cell_len(&self) -> usize {
        match &self.single_cell_system {
            Some(SingleCellSystem::Rhapsody(_)) => 27,
            Some(SingleCellSystem::Tenx(system)) => system.version().cell_len(),
            Some(SingleCellSystem::Whitelist(_)) => self.grammar.cell_len(),
            None => self.grammar.cell_len(),
        }
    }
    /// 	Convert a normalized cell-barcode sequence into the chemistry-specific
    /// 	numeric cell ID.
    /// 	Returns `Some(cell_id)` when the configured single-cell system supports
    /// 	sequence-to-ID translation and the barcode is recognized. Returns `None`
    /// 	otherwise.
    pub fn cell_id_for_seq(&self, seq: &[u8]) -> Option<u64> {
        match &self.single_cell_system {
            Some(SingleCellSystem::Rhapsody(rhapsody)) => rhapsody.cell_id_for_seq(seq),
            _ => None,
        }
    }

    pub fn detect(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerMatch>> {
        self.detect_all(seq, qual)
    }

    /// Generate a new primer from the internal Grammar and an external ReadTagRecord.
    ///
    /// Returns:
    /// - `(new_cell_seq, synthetic_primer_seq)`
    ///
    /// If the incoming cell/UMI lengths match the grammar, reuse them.
    /// Otherwise translate the source cell to a stable generated target cell
    /// and generate a fresh/random UMI of the grammar UMI length.
    pub fn generate(&mut self, data: &ReadTagRecord) -> PrimerResult<(Vec<u8>, Vec<u8>)> {
        let target_cell = self.get_or_create_export_cell(&data.cell_seq)?;
        let target_umi = self.get_or_create_export_umi(&data.umi_seq)?;

        let primer = self.grammar.synthesize(&target_cell, &target_umi)?;

        Ok((target_cell, primer))
    }

    fn get_or_create_export_umi(&mut self, source_umi: &[u8]) -> PrimerResult<Vec<u8>> {
        if source_umi.len() == self.grammar.umi_len() {
            return Ok(source_umi.to_vec());
        }

        let source_id = IntToStr::new(source_umi).into_u64();

        let target_index = if let Some(index) = self.umi_translation.get(&source_id) {
            *index
        } else {
            let index = self.umi_translation.len();
            self.umi_translation.insert(source_id, index);
            index
        };

        Ok(IntToStr::from_u64(target_index as u64)
            .to_string(self.grammar.umi_len())
            .into_bytes())
    }

    pub fn primer_translation(&self) -> &HashMap<u64, (u64, usize)> {
        &self.primer_translation
    }

    /// Creates a tab separated file:
    /// old_cell_id  new_cell_id  reads_detected
    pub fn save_cell_translation_table<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut out = BufWriter::new(file);

        writeln!(out, "old_cell_id\tnew_cell_id\treads_detected")?;

        for (old_cell, (new_cell, count)) in &self.primer_translation {
            writeln!(
                out,
                "{}\t{}\t{}",
                IntToStr::from_u64(*old_cell).to_string(self.grammar.cell_len()),
                IntToStr::from_u64(*new_cell).to_string(self.grammar.umi_len()),
                count,
            )?;
        }

        Ok(())
    }

    fn get_or_create_export_cell(&mut self, source_cell: &[u8]) -> PrimerResult<Vec<u8>> {
        // Already valid in the target system? Use it as-is.
        if source_cell.len() == self.grammar.cell_len() {
            if let Some(system) = &self.single_cell_system {
                if system.cell_index_for_seq(source_cell).is_some() {
                    return Ok(source_cell.to_vec());
                }
            } else {
                return Ok(source_cell.to_vec());
            }
        }

        // Otherwise translate source seq -> source id -> target allocation index.
        let source_id = IntToStr::new(source_cell).into_u64();

        let target_index =
            if let Some((target_index, count)) = self.primer_translation.get_mut(&source_id) {
                *count += 1;
                *target_index
            } else {
                let target_index = self.primer_translation.len() as u64;

                self.primer_translation.insert(source_id, (target_index, 1));

                target_index
            };

        self.cell_seq_for_index(target_index)
    }

    pub fn cell_seq_for_index(&self, target_index: u64) -> PrimerResult<Vec<u8>> {
        if let Some(system) = &self.single_cell_system {
            return system.cell_seq_for_index(target_index).ok_or_else(|| {
                PrimerError::invalid_coordinates(format!(
                    "single-cell whitelist exhausted at index {target_index}"
                ))
            });
        }

        Ok(IntToStr::from_u64(target_index)
            .to_string(self.grammar.cell_len())
            .into_bytes())
    }

    pub fn detect_first(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Option<PrimerMatch>> {
        for detector in std::iter::once(self).chain(self.alternates.iter()) {
            if let Some(hit) = detector.detect_first_single(seq, qual)? {
                return Ok(Some(hit));
            }
        }
        Ok(None)
    }

    fn detect_first_single(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Option<PrimerMatch>> {
        Self::validate_read(seq, qual)?;

        if let Some(mut hit) = self.detect_first_forward(seq, qual, Orientation::Forward)? {
            hit.insert_end = seq.len();
            return Ok(Some(hit));
        }

        if self.detect_reverse_complement {
            let rc_seq = Self::reverse_complement(seq);
            let rc_qual = Self::reverse(qual);

            if let Some(mut hit) =
                self.detect_first_forward(&rc_seq, &rc_qual, Orientation::ReverseComplement)?
            {
                hit.insert_end = rc_seq.len();
                hit.remap_reverse_coordinates(seq.len());
                return Ok(Some(hit));
            }
        }

        Ok(None)
    }

    /// Detect the first primer match and return the compact diagnostics that
    /// were collected while producing that match. Callers uninterested in
    /// diagnostics can continue using `detect_first`, or destructure this as
    /// `(hit, _)` without re-running any matching work.
    pub fn detect_first_with_diagnostics(
        &self,
        seq: &[u8],
        qual: &[u8],
    ) -> PrimerResult<(Option<PrimerMatch>, Option<PrimerMatchDiagnostics>)> {
        let hit = self.detect_first(seq, qual)?;
        let diagnostics = hit.as_ref().map(|hit| hit.diagnostics);
        Ok((hit, diagnostics))
    }

    pub fn detect_all(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerMatch>> {
        for detector in std::iter::once(self).chain(self.alternates.iter()) {
            let hits = detector.detect_all_single(seq, qual)?;
            if !hits.is_empty() {
                return Ok(hits);
            }
        }
        Ok(Vec::new())
    }

    fn detect_all_single(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerMatch>> {
        Self::validate_read(seq, qual)?;

        let hits = self.detect_all_forward(seq, qual, Orientation::Forward)?;
        if !hits.is_empty() || !self.detect_reverse_complement {
            return Ok(hits);
        }

        let rc_seq = Self::reverse_complement(seq);
        let rc_qual = Self::reverse(qual);
        let mut reverse_hits =
            self.detect_all_forward(&rc_seq, &rc_qual, Orientation::ReverseComplement)?;

        for hit in &mut reverse_hits {
            hit.remap_reverse_coordinates(seq.len());
        }

        Ok(reverse_hits)
    }

    fn detect_all_forward(
        &self,
        seq: &[u8],
        qual: &[u8],
        orientation: Orientation,
    ) -> PrimerResult<Vec<PrimerMatch>> {
        let mut hits = Vec::new();
        let mut cursor = 0usize;
        let packed = self
            .leading_bd_search()
            .map(|_| OneHotSequence::from_bytes(seq));

        while let Some(offset) = self.next_candidate_start(seq, packed.as_ref(), cursor) {
            match self.try_from_start_with_packed(
                seq,
                qual,
                offset,
                orientation,
                packed.as_ref(),
            )? {
                Some(hit) => {
                    // A successful match owns the whole primer span. Continue
                    // at its end instead of re-testing SEARCH shifts or other
                    // candidate starts inside the primer we just accepted.
                    cursor = hit.primer_end.max(offset.saturating_add(1));
                    hits.push(hit);
                }
                None => {
                    cursor = offset.saturating_add(1);
                }
            }
        }

        self.finish_insert_ends(&mut hits, seq.len());
        Ok(hits)
    }

    fn detect_first_forward(
        &self,
        seq: &[u8],
        qual: &[u8],
        orientation: Orientation,
    ) -> PrimerResult<Option<PrimerMatch>> {
        let mut cursor = 0usize;
        let packed = self
            .leading_bd_search()
            .map(|_| OneHotSequence::from_bytes(seq));

        while let Some(offset) = self.next_candidate_start(seq, packed.as_ref(), cursor) {
            if let Some(hit) =
                self.try_from_start_with_packed(seq, qual, offset, orientation, packed.as_ref())?
            {
                return Ok(Some(hit));
            }
            cursor = offset.saturating_add(1);
        }

        Ok(None)
    }

    #[inline]
    fn validate_read(seq: &[u8], qual: &[u8]) -> PrimerResult<()> {
        if seq.len() != qual.len() {
            return Err(PrimerError::invalid_coordinates(
                "sequence and quality have different lengths",
            ));
        }
        Ok(())
    }

    /// Return the next plausible grammar start without first materializing all
    /// offsets in the read.
    fn next_candidate_start(
        &self,
        seq: &[u8],
        _packed: Option<&OneHotSequence>,
        from: usize,
    ) -> Option<usize> {
        if let Some(anchor) = self.grammar.anchor_search() {
            return anchor.find_next_cell_start(seq, from);
        }

        // Built-in 10x read structures are positional: the cell barcode begins
        // at the start of R1.  Do not turn the millions-entry whitelist into a
        // read-wide search index by retrying it at every base.  Long-read / ONT
        // cassette discovery must provide an explicit anchor/search grammar
        // instead of relying on TENX_CELL itself to discover coordinates.
        if self.leading_tenx_cell() {
            return (from == 0 && !seq.is_empty()).then_some(0);
        }

        // Built-in BD read structures are positional, just like 10x. The
        // leading SEARCH operation describes the small tolerated shift of the
        // cassette itself (normally 0..=4); it must not turn the linker into a
        // read-wide search index. This rule is applied identically to forward
        // and reverse-complement detection.
        if self.leading_bd_search().is_some()
            && matches!(
                &self.single_cell_system,
                Some(SingleCellSystem::Rhapsody(_))
            )
        {
            return (from == 0 && !seq.is_empty()).then_some(0);
        }

        (from < seq.len()).then_some(from)
    }

    /// A built-in TENX_CELL without an explicit leading anchor/search is a
    /// positional R1 grammar.  Its barcode starts at offset zero.
    fn leading_tenx_cell(&self) -> bool {
        matches!(
            self.grammar.ops.as_slice(),
            [GrammarOp::TenxCell { .. }, ..]
        ) && matches!(&self.single_cell_system, Some(SingleCellSystem::Tenx(_)))
    }

    /// SEARCH directly before BD_CELL is the only flexible prefix used by the
    /// built-in BD v2 chemistries.  Keep this deliberately narrow so arbitrary
    /// custom grammars retain the generic detector semantics.
    fn leading_bd_search(&self) -> Option<(usize, usize)> {
        match self.grammar.ops.as_slice() {
            [GrammarOp::Search { start, end }, GrammarOp::BdCell { .. }, ..] => {
                Some((*start, *end))
            }
            [GrammarOp::BdCell { .. }, ..] => Some((0, 0)),
            _ => None,
        }
    }

    pub fn explain_all(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<Vec<PrimerAttempt>> {
        if seq.len() != qual.len() {
            return Err(PrimerError::invalid_coordinates(
                "sequence and quality have different lengths",
            ));
        }

        let mut attempts = Vec::new();

        for offset in 0..seq.len() {
            //let sub_seq = &seq[offset..];
            //let sub_qual = &qual[offset..];
            match self.try_from_start(seq, qual, offset, Orientation::Forward)? {
                Some(hit) => {
                    attempts.push(PrimerAttempt {
                        offset,
                        orientation: Orientation::Forward,
                        ok: true,
                        reason: "matched".to_string(),
                        cell_seq: hit
                            .cell_seq
                            .as_ref()
                            .map(|seq| String::from_utf8_lossy(seq).to_string()),
                        segments: hit
                            .segments
                            .iter()
                            .flat_map(|segment| {
                                segment.ranges.iter().map(|range| PrimerSegmentAttempt {
                                    name: segment.name.clone(),
                                    range: range.clone(),
                                    dna: String::from_utf8_lossy(&seq[range.start..range.end])
                                        .to_string(),
                                    ok: true,
                                    reason: "matched".to_string(),
                                })
                            })
                            .collect(),
                    });

                    break;
                }

                None => {
                    attempts.push(PrimerAttempt {
                        offset,
                        orientation: Orientation::Forward,
                        ok: false,
                        reason: "no complete primer match".to_string(),
                        cell_seq: None,
                        segments: Vec::new(),
                    });
                }
            }
        }

        Ok(attempts)
    }

    pub fn explain_failure(&self, seq: &[u8], qual: &[u8]) -> PrimerResult<String> {
        Self::validate_read(seq, qual)?;

        if let Some((search_start, search_end)) = self.leading_bd_search() {
            if let Some(SingleCellSystem::Rhapsody(rhapsody)) = &self.single_cell_system {
                return Ok(rhapsody.explain_call_failure(seq, 0, search_start, search_end));
            }
        }

        Ok("no complete primer match".to_string())
    }

    /// Primary grammar after rejection-cost ordering. Kept for legacy callers
    /// that configure one chemistry or need one representative grammar.
    pub fn grammar(&self) -> &Grammar {
        &self.grammar
    }

    pub fn grammars(&self) -> impl Iterator<Item = &Grammar> {
        std::iter::once(&self.grammar).chain(self.alternates.iter().map(|d| &d.grammar))
    }

    /// Resolve the grammar that produced a PrimerMatch.
    pub fn grammar_for_match(&self, hit: &PrimerMatch) -> &Grammar {
        self.grammars()
            .find(|grammar| grammar.name == hit.chemistry_name)
            .unwrap_or(&self.grammar)
    }

    pub fn detect_one_orientation(
        &self,
        seq: &[u8],
        qual: &[u8],
        orientation: Orientation,
    ) -> PrimerResult<Option<PrimerMatch>> {
        self.try_from_start(seq, qual, 0, orientation)
    }

    pub fn try_from_start(
        &self,
        seq: &[u8],
        qual: &[u8],
        start: usize,
        orientation: Orientation,
    ) -> PrimerResult<Option<PrimerMatch>> {
        self.try_from_start_with_packed(seq, qual, start, orientation, None)
    }

    fn try_from_start_with_packed(
        &self,
        seq: &[u8],
        qual: &[u8],
        start: usize,
        orientation: Orientation,
        packed: Option<&OneHotSequence>,
    ) -> PrimerResult<Option<PrimerMatch>> {
        // Most candidate starts fail.  Do not clone the chemistry name for
        // every failed attempt; attach it only when the whole grammar matches.
        let mut primer_match = PrimerMatch::new(String::new(), orientation);
        primer_match.primer_start = start;
        let mut pos = start;
        let mut search = (0usize, 0usize);

        //let mut saw_insert_constraint = false;
        //let matched_insert = false;

        for op in &self.grammar.ops {
            match op {
                GrammarOp::Fixed {
                    seq: fixed,
                    mismatches,
                } => {
                    let chosen = Self::find_fixed(seq, pos, fixed, *mismatches, search)?;
                    let Some(next_pos) = chosen else {
                        return Ok(None);
                    };
                    pos = next_pos + fixed.len();
                    search = (0, 0);
                }
                GrammarOp::Cell { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) {
                        return Ok(None);
                    }

                    if let Some(SingleCellSystem::Whitelist(whitelist)) = &self.single_cell_system {
                        let observed = &seq[pos..pos + *len];
                        let Some((cell_index, _distance)) = whitelist.best_match(observed) else {
                            return Ok(None);
                        };
                        let Some(corrected) = whitelist.sequence(cell_index) else {
                            return Ok(None);
                        };
                        primer_match.add_cell_seq(&corrected);
                    }

                    primer_match.add_segment("CELL", pos..pos + *len);
                    pos += *len;
                }
                GrammarOp::Umi { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) {
                        return Ok(None);
                    }
                    primer_match.add_segment("UMI", pos..pos + *len);
                    pos += *len;
                }
                GrammarOp::PolyT { min } => {
                    let count = Self::count_base(seq, pos, b'T');
                    if count < *min {
                        return Ok(None);
                    }
                    pos += count;
                }
                GrammarOp::Insert {
                    seq: fixed,
                    mismatches,
                } => {
                    //saw_insert_constraint = true;

                    let Ok(Some(_chosen)) = Self::find_fixed(seq, pos, fixed, *mismatches, search)
                    else {
                        /*eprintln!(
                            "I have not found the INSERT {} in sequence {} at position {}",
                            String::from_utf8_lossy(fixed),
                            String::from_utf8_lossy(seq),
                            pos
                        );*/
                        return Ok(None);
                    };

                    primer_match.insert_start = pos;
                    primer_match.insert_end = seq.len();
                    primer_match.primer_end = pos;
                    primer_match.chemistry_name = self.grammar.name.clone();
                    return Ok(Some(primer_match));
                }
                GrammarOp::Skip { len } => {
                    if !Self::has_range(seq, pos, *len) {
                        return Ok(None);
                    }
                    pos += *len;
                }
                GrammarOp::Search { start, end } => {
                    search = (*start, *end);
                }
                GrammarOp::BdCell { version } => {
                    let Some(SingleCellSystem::Rhapsody(rhapsody)) = &self.single_cell_system
                    else {
                        return Ok(None);
                    };

                    if rhapsody.version() != *version {
                        return Ok(None);
                    }

                    let call = match packed {
                        Some(packed)
                            if matches!(
                                rhapsody.version(),
                                BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj
                            ) =>
                        {
                            rhapsody.call_with_packed(seq, qual, packed, pos, search.0, search.1)
                        }
                        _ => rhapsody.call(seq, qual, pos, search.0, search.1),
                    };
                    let Some(call) = call else {
                        return Ok(None);
                    };

                    primer_match.bd_cell_id = Some(call.cell_id);
                    if let Some(diagnostics) = call.diagnostics {
                        primer_match.diagnostics.bd = Some(BdPrimerDiagnostics {
                            shift: call.shift,
                            linker_signature: diagnostics.linker_signature,
                            linker_mismatches: diagnostics.linker_mismatches,
                            c1_mismatches: diagnostics.c1_mismatches,
                            c2_mismatches: diagnostics.c2_mismatches,
                            c3_mismatches: diagnostics.c3_mismatches,
                        });
                    }

                    // probably use full cassette if this is later used for synthesize()
                    primer_match.add_cell_seq(&call.cell_seq);

                    primer_match.add_segment_ranges(
                        "BD_CELL",
                        vec![
                            call.c1.0..call.c1.1,
                            call.c2.0..call.c2.1,
                            call.c3.0..call.c3.1,
                        ],
                    );

                    primer_match.add_segment("UMI", call.umi.0..call.umi.1);

                    pos = pos.checked_add(call.consumed).ok_or_else(|| {
                        PrimerError::invalid_coordinates("BD consumed coordinate overflow")
                    })?;
                    search = (0, 0);
                }

                GrammarOp::TenxCell { version } => {
                    let len = version.cell_len();

                    if seq.len() < pos + len {
                        return Ok(None);
                    }

                    let Some(SingleCellSystem::Tenx(tenx)) = &self.single_cell_system else {
                        return Err(PrimerError::invalid_grammar(
                            "TENX_CELL grammar has no 10x whitelist",
                        ));
                    };

                    let observed = &seq[pos..pos + len];
                    let Some(cell_index) = tenx.index_cell(observed) else {
                        return Ok(None);
                    };
                    let Some(corrected) = tenx.cell_id_to_seq(cell_index + 1) else {
                        return Ok(None);
                    };

                    primer_match.add_segment("CELL", pos..pos + len);
                    primer_match.add_cell_seq(&corrected);

                    pos += len;
                }
                GrammarOp::Tag { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) {
                        return Ok(None);
                    }
                    primer_match.add_segment("TAG", pos..pos + *len);
                    pos += *len;
                }
                GrammarOp::Feature { len } => {
                    if !Self::has_range(seq, pos, *len) || !Self::has_range(qual, pos, *len) {
                        return Ok(None);
                    }
                    primer_match.add_segment("FEATURE", pos..pos + *len);
                    pos += *len;
                }
            }
        }
        primer_match.primer_end = pos;
        primer_match.insert_start = pos;
        primer_match.insert_end = seq.len();
        primer_match.chemistry_name = self.grammar.name.clone();
        Ok(Some(primer_match))
    }

    pub fn finish_insert_ends(&self, hits: &mut [PrimerMatch], seq_len: usize) {
        if hits.is_empty() {
            return;
        }
        for idx in 0..hits.len() {
            hits[idx].insert_end = if idx + 1 < hits.len() {
                hits[idx + 1].primer_start
            } else {
                seq_len
            };
        }
    }

    pub fn has_range(seq: &[u8], start: usize, len: usize) -> bool {
        start.checked_add(len).is_some_and(|end| end <= seq.len())
    }

    pub fn count_base(seq: &[u8], start: usize, base: u8) -> usize {
        seq.iter().skip(start).take_while(|b| **b == base).count()
    }

    pub fn find_fixed(
        seq: &[u8],
        pos: usize,
        fixed: &[u8],
        mismatches: usize,
        search: (usize, usize),
    ) -> PrimerResult<Option<usize>> {
        let start = pos + search.0;
        let end = pos + search.1;
        for candidate in start..=end {
            if !Self::has_range(seq, candidate, fixed.len()) {
                continue;
            }
            if Self::hamming(&seq[candidate..candidate + fixed.len()], fixed)? <= mismatches {
                return Ok(Some(candidate));
            }
        }
        Ok(None)
    }

    pub fn hamming(left: &[u8], right: &[u8]) -> PrimerResult<usize> {
        if left.len() != right.len() {
            return Err(PrimerError::invalid_coordinates(
                "hamming requires equal lengths",
            ));
        }
        Ok(left
            .iter()
            .zip(right.iter())
            .filter(|(a, b)| a != b)
            .count())
    }

    pub fn reverse(seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().copied().collect()
    }

    pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().map(|b| Self::complement(*b)).collect()
    }

    pub fn complement(base: u8) -> u8 {
        match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other,
        }
    }

    pub fn version_from_bd_op(version: BdCellVersion) -> BdCellVersion {
        version
    }
}
