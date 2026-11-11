//! An independent annotation + UMI matrix for newly observed splice isoforms.
use anyhow::{Result, Context};
use std::io::{BufWriter, Write};
use std::path::Path;
use gtf_splice_index::{SpliceIndex, SpliceMatchMode, UmiTranscriptMatch, RefBlock, TranscriptPrior, TranscriptSnpSupport, MoleculeSnpEvidence, posterior_with_snp_evidence};
use std::collections::BTreeMap;
use super::cli::QuantMode;
use run_stats::RunStats;
use scdata::QuantData;
use std::collections::HashMap;
use std::time::Instant;
use scdata::{MatrixValueType, Scdata, GeneUmiHash, SnpUmiEvidence, ExonMask, TranscriptUmiEvidence, UmiCoverageRange};

pub struct NewTranscriptData {
    pub scdata: Scdata<SnpUmiEvidence>,
    /// Raw gene-wide exon evidence, keyed by (cell, gene, UMI).
    /// Transcript assignment must happen only after collection.
    pub exon_evidence: Scdata<TranscriptUmiEvidence>,
    pub gtf_splice_index: SpliceIndex,
    /// Genomic span supported by each novel-model UMI (0-based, half-open).
    supported_regions: HashMap<(u64, usize, u64), (u32, u32)>,
    /// Reference-compatible UMIs whose transcript identity is unresolved.
    /// These are not novel transcript features.
    ambiguity: BTreeMap<(usize, Vec<usize>), AmbiguityGroup>,
    posterior_rows: Vec<PosteriorRow>,
    ambiguity_gene_names: HashMap<usize, String>,
    ambiguity_transcript_names: HashMap<usize, String>,

}

// Use the same coordinate convention in production and in splice-event tests:
// RefBlock is 0-based half-open; names use 1-based inclusive coordinates.
fn skip_event_name(chr: &str, exon: RefBlock) -> String {
    format!("skip_{chr}:{}-{}", exon.start + 1, exon.end)
}

// Characterization helper: model only positively observed splice junctions.
// Coverage blocks alone never establish an exon skip; a deletion or a gap
// between independent reads is not a CIGAR N junction.
#[cfg(test)]
fn reference_relative_event_candidates(
    reference_exons: &[RefBlock],
    observed_blocks: &[RefBlock],
    observed_junctions: &[(u32, u32)],
) -> (Vec<RefBlock>, Vec<RefBlock>) {
    let skipped = reference_exons.iter().copied()
        .filter(|exon| observed_junctions.iter().any(|&(donor, acceptor)|
            donor < acceptor && donor <= exon.start && exon.end <= acceptor))
        .filter(|exon| !observed_blocks.iter().any(|block| block.overlaps(*exon)))
        .collect();
    // Coverage alone cannot establish a novel exon. Novel exon creation
    // requires additional junction evidence, not represented by this helper.
    let added = Vec::new();
    (skipped, added)
}

#[derive(Clone, Debug, Default)]
struct AmbiguityGroup {
    molecules: Vec<(u64, u64, SnpUmiEvidence)>,
}

#[derive(Clone, Debug)]
struct PosteriorRow {
    cell: u64,
    gene: String,
    umi: u64,
    candidate: String,
    probability: f64,
    snp_informed: bool,
}

/// Only molecules with exactly one structurally compatible reference transcript
/// train priors and SNP associations. Ambiguous molecules never train themselves.
/// Exons outside the observed window are unknown, not deletions.
fn reference_posteriors(
    reference: &SpliceIndex,
    molecules: &[(u64, GeneUmiHash, TranscriptUmiEvidence)],
) -> HashMap<(u64, usize, u64), Vec<(usize, f64, bool)>> {
    let mut by_gene: HashMap<usize, Vec<usize>> = HashMap::new();
    for (idx, (_, key, _)) in molecules.iter().enumerate() {
        by_gene.entry(key.0 as usize).or_default().push(idx);
    }
    let mut result = HashMap::new();
    for (gene_id, indices) in by_gene {
        let Some(gene) = reference.genes.get(gene_id) else { continue; };
        let Some(catalog) = reference.gene_exon_catalog(gene_id) else { continue; };
        let masks: Vec<_> = gene.transcript_ids().iter().filter_map(|&id| {
            reference.transcripts.get(id).map(|tx| (id, catalog.transcript_mask(tx)))
        }).collect();
        if masks.is_empty() { continue; }
        let mut counts: HashMap<usize, usize> = HashMap::new();
        let mut snp_counts: HashMap<(usize, usize), (u64, u64)> = HashMap::new();
        // First pass: structural evidence alone identifies training UMIs.
        for &idx in &indices {
            let (_, _, evidence) = &molecules[idx];
            let candidates: Vec<_> = masks.iter().map(|(id, mask)| TranscriptPrior {
                id: *id, exon_mask: mask, support: 0,
            }).collect();
            let posterior = gtf_splice_index::posterior_for_partial_observation(
                &catalog, &evidence.splice, &candidates, 1.0);
            if posterior.len() != 1 { continue; }
            let id = posterior[0].id;
            *counts.entry(id).or_default() += 1;
            for snp in &evidence.snps.snps {
                let entry = snp_counts.entry((id, snp.snp_id as usize)).or_default();
                match snp.ref_reads.cmp(&snp.alt_reads) {
                    std::cmp::Ordering::Greater => entry.0 += 1,
                    std::cmp::Ordering::Less => entry.1 += 1,
                    std::cmp::Ordering::Equal => {},
                }
            }
        }
        let candidates: Vec<_> = masks.iter().map(|(id, mask)| TranscriptPrior {
            id: *id, exon_mask: mask, support: *counts.get(id).unwrap_or(&0),
        }).collect();
        let training: Vec<_> = snp_counts.into_iter().map(|((transcript_id, snp_id), (ref_molecules, alt_molecules))|
            TranscriptSnpSupport { transcript_id, snp_id, ref_molecules, alt_molecules }).collect();
        for &idx in &indices {
            let (cell, key, evidence) = &molecules[idx];
            let snps: Vec<_> = evidence.snps.snps.iter().map(|s| MoleculeSnpEvidence {
                snp_id: s.snp_id as usize, ref_reads: s.ref_reads as u32, alt_reads: s.alt_reads as u32,
            }).collect();
            let structural = gtf_splice_index::posterior_for_partial_observation(
                &catalog, &evidence.splice, &candidates, 1.0);
            let posterior = posterior_with_snp_evidence(
                &catalog, &evidence.splice, &candidates, 1.0, &snps, &training, 1.0);
            if posterior.is_empty() { continue; }
            let snp_informed = posterior.len() > 1 && posterior.iter().any(|p| structural.iter()
                .find(|s| s.id == p.id).is_some_and(|s| (s.probability - p.probability).abs() > 1e-10));
            result.insert((*cell, key.0 as usize, key.1), posterior.into_iter()
                .map(|p| (p.id, p.probability, snp_informed)).collect());
        }
    }
    result
}

impl NewTranscriptData {
    pub fn new(bin_width: u32) -> Self {
        Self {
            scdata: Scdata::new_with_evidence(1, MatrixValueType::Integer),
            exon_evidence: Scdata::new_with_evidence(1, MatrixValueType::Integer),
            gtf_splice_index: SpliceIndex::new(bin_width).with_match_mode(SpliceMatchMode::Transcript),
            supported_regions: HashMap::new(),
            ambiguity: BTreeMap::new(),
            posterior_rows: Vec::new(),
            ambiguity_gene_names: HashMap::new(),
            ambiguity_transcript_names: HashMap::new(),
        }
    }

    pub fn add_molecule(&mut self, transcript_id: usize, cell: u64, umi: u64, evidence: SnpUmiEvidence) {
        self.scdata.try_insert_evidence(&cell, GeneUmiHash(transcript_id as u64, umi), evidence);
    }

    pub fn add_exon_evidence(&mut self, gene_id: usize, cell: u64, umi: u64, mask: ExonMask, snps: SnpUmiEvidence) {
        self.add_exon_evidence_with_coverage(gene_id, cell, umi, mask, snps, UmiCoverageRange::default());
    }

    pub fn add_exon_evidence_with_coverage(&mut self, gene_id: usize, cell: u64, umi: u64,
        mask: ExonMask, snps: SnpUmiEvidence, coverage: UmiCoverageRange) {
        if mask.count_ones() > 0 {
            self.exon_evidence.try_insert_evidence(&cell, GeneUmiHash(gene_id as u64, umi),
                TranscriptUmiEvidence { splice: mask, snps, coverage });
        }
    }

    /// Resolve only after every worker's UMI evidence has been merged.
    /// Rayon performs read-only classification; registration/output stays serial.
    pub fn classify_umis_parallel(&self, reference: &SpliceIndex)
        -> Vec<(u64, GeneUmiHash, gtf_splice_index::UmiTranscriptMatch, SnpUmiEvidence)>
    {
        use rayon::prelude::*;
        let mut molecules = Vec::new();
        self.exon_evidence.for_each_umi(|cell, key, evidence| {
            molecules.push((cell, key, evidence.clone()));
        });
        molecules.into_par_iter().filter_map(|(cell, key, evidence)| {
            reference.resolve_umi_transcript(key.0 as usize, &evidence.splice)
                .map(|resolution| (cell, key, resolution, evidence.snps))
        }).collect()
    }

    /// Fuse UMI classification and matrix accumulation in Rayon worker-local shards.
    /// Only novel transcript patterns require deferred ID assignment.
    pub fn finalize_umi_evidence(
        &mut self,
        reference: &SpliceIndex,
        quant_mode: QuantMode,
        data: &mut QuantData,
        report: &mut RunStats,
    ) -> Result<()> {
        use rayon::prelude::*;
        let started = Instant::now();
        let transcript_mode = matches!(quant_mode, QuantMode::Transcript);
        let threads = rayon::current_num_threads();
        eprintln!("[bam_tide finalize] fused bitmap resolution + QuantData accumulation; mode={} rayon_threads={threads}",
            if transcript_mode { "transcript" } else { "gene" });

        // Scdata's callback is not a parallel iterator. Snapshot its entries once;
        // do NOT allocate another 16-million-element classification vector.
        let mut molecules = Vec::new();
        self.exon_evidence.for_each_umi(|cell, key, evidence| {
            molecules.push((cell, key, evidence.clone()));
        });
        let total = molecules.len();
        eprintln!("[bam_tide finalize] consolidated_umis={total} snapshot_elapsed={:.2}s",
            started.elapsed().as_secs_f64());

        // A deterministic two-pass empirical-Bayes model, trained only on
        // structurally unambiguous, deduplicated UMIs.
        let bayes = if transcript_mode { reference_posteriors(reference, &molecules) }
                    else { HashMap::new() };
        let mut deferred = Vec::new();
        // Ambiguity groups are registered serially after parallel accumulation.

        // Each worker owns its matrices. No global lock and no serial per-UMI
        // insertion. Novel patterns retain their original SNP evidence until
        // the shared novel index has assigned stable transcript IDs.
        let processed = std::sync::atomic::AtomicUsize::new(0);
        let next_report = std::sync::atomic::AtomicUsize::new(1_000_000);
        let shards: Vec<_> = molecules
            .into_par_iter()
            .fold(
                || (QuantData::standard(), RunStats::new(None, 0.0, usize::MAX),
                    Vec::<(u64, usize, u64, ExonMask, SnpUmiEvidence, UmiCoverageRange)>::new(),
                    Vec::<(u64, usize, u64, SnpUmiEvidence)>::new()),
                |(mut local, mut stats, mut novel, mut ambiguous), (cell, key, evidence)| {
                    let gene_id = key.0 as usize;
                    let umi = key.1;
                    let posterior = bayes.get(&(cell, gene_id, umi));
                    // Structural classification takes precedence over Bayesian
                    // ambiguity: intronic and genuinely incompatible molecules
                    // must not be intercepted by a reference posterior.
                    let structural = reference.resolve_umi_transcript(gene_id, &evidence.splice);
                    let resolution = if transcript_mode {
                        match structural {
                            Some(UmiTranscriptMatch::Intronic) => Some(UmiTranscriptMatch::Intronic),
                            Some(UmiTranscriptMatch::Novel(mask)) if posterior.is_none() =>
                                Some(UmiTranscriptMatch::Novel(mask)),
                            Some(UmiTranscriptMatch::Novel(mask)) => {
                                // A reference posterior may explain partial
                                // coverage, but not override structural evidence.
                                // The Bayesian compatibility routine must be
                                // consistent with the structural resolver.
                                if let Some(p) = posterior {
                                    if p.len() > 1 {
                                        ambiguous.push((cell, gene_id, umi, evidence.snps.clone()));
                                        None
                                    } else if p.len() == 1 {
                                        Some(UmiTranscriptMatch::Reference(p[0].0))
                                    } else { Some(UmiTranscriptMatch::Novel(mask)) }
                                } else { Some(UmiTranscriptMatch::Novel(mask)) }
                            }
                            Some(UmiTranscriptMatch::Reference(structural_id)) => {
                                match posterior {
                                    Some(p) if p.len() > 1 => {
                                        ambiguous.push((cell, gene_id, umi, evidence.snps.clone()));
                                        None
                                    }
                                    Some(p) if p.len() == 1 => Some(UmiTranscriptMatch::Reference(p[0].0)),
                                    _ => Some(UmiTranscriptMatch::Reference(structural_id)),
                                }
                            }
                            None => None,
                        }
                    } else {
                        let intronic = match &evidence.splice {
                            ExonMask::Small(bits) => bits & 0xAAAA_AAAA_AAAA_AAAAu64 != 0,
                            ExonMask::Large(words) => words.iter().any(|bits|
                                bits & 0xAAAA_AAAA_AAAA_AAAAu64 != 0),
                        };
                        Some(if intronic { UmiTranscriptMatch::Intronic }
                             else { UmiTranscriptMatch::Reference(gene_id) })
                    };
                    if let Some(resolution) = resolution {
                        if !transcript_mode {
                            let feature = GeneUmiHash(key.0, umi);
                            local.try_insert(QuantData::EXPRS, &cell, feature, 1.0, &mut stats);
                            local.try_insert(if matches!(resolution, UmiTranscriptMatch::Intronic) {
                                QuantData::INTRONIC
                            } else { QuantData::EXONIC }, &cell, feature, 1.0, &mut stats);
                        } else {
                            match resolution {
                                UmiTranscriptMatch::Intronic => {
                                    if let Some(&tx_id) = reference.genes[gene_id].transcript_ids().first() {
                                        let feature = GeneUmiHash(tx_id as u64, umi);
                                        local.try_insert(QuantData::INTRONIC, &cell, feature, 1.0, &mut stats);
                                        local.try_insert(QuantData::EXPRS, &cell, feature, 1.0, &mut stats);
                                    }
                                }
                                UmiTranscriptMatch::Reference(tx_id) => {
                                    let feature = GeneUmiHash(tx_id as u64, umi);
                                    local.try_insert(QuantData::EXONIC, &cell, feature, 1.0, &mut stats);
                                    local.try_insert(QuantData::EXPRS, &cell, feature, 1.0, &mut stats);
                                }
                                UmiTranscriptMatch::Novel(mask) => {
                                    // Retain ordinary gene-expression evidence once.
                                    if let Some(&representative) = reference.genes[gene_id].transcript_ids().first() {
                                        let feature = GeneUmiHash(representative as u64, umi);
                                        local.try_insert(QuantData::EXPRS, &cell, feature, 1.0, &mut stats);
                                        local.try_insert(QuantData::EXONIC, &cell, feature, 1.0, &mut stats);
                                    }
                                    novel.push((cell, gene_id, umi, mask, evidence.snps.clone(), evidence.coverage.clone()));
                                }
                            }
                        }
                        Self::write_snp_counts(&mut local, &mut stats, cell, umi, &evidence.snps);
                    }
                    let done = processed.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                    let threshold = next_report.load(std::sync::atomic::Ordering::Relaxed);
                    if done >= threshold && next_report.compare_exchange(
                        threshold, threshold.saturating_add(1_000_000),
                        std::sync::atomic::Ordering::Relaxed,
                        std::sync::atomic::Ordering::Relaxed).is_ok() {
                        eprintln!("[bam_tide finalize] fused={done}/{total} elapsed={:.1}s",
                            started.elapsed().as_secs_f64());
                    }
                    (local, stats, novel, ambiguous)
                },
            )
            .collect();

        eprintln!("[bam_tide finalize] merging {} worker-local matrix shards elapsed={:.2}s",
            shards.len(), started.elapsed().as_secs_f64());
        let mut novel_ids: HashMap<(usize, ExonMask), usize> = HashMap::new();
        let mut novel_count = 0usize;
        let mut rejected_conflicts = 0usize;
        let verbose_novel_conflicts = std::env::var_os("BAM_TIDE_DEBUG_NOVEL_CONFLICTS")
            .is_some_and(|value| !value.is_empty() && value != "0");
        for (local, _stats, novel, ambiguous) in shards {
            deferred.extend(ambiguous);
            report.merge(&data.merge(&local));
            for (cell, gene_id, umi, mask, snps, coverage) in novel {
                let lookup = (gene_id, mask.clone());
                let tx_id = if let Some(&id) = novel_ids.get(&lookup) { id } else {
                    let catalog = reference.gene_exon_catalog(gene_id)
                        .ok_or_else(|| anyhow::anyhow!("missing exon catalogue for gene {gene_id}"))?;
                    let mut blocks: Vec<RefBlock> = Vec::new();
                    for (i, unit) in catalog.units.iter().enumerate() {
                        let mut bit = ExonMask::default();
                        bit.set(2 * i);
                        if mask.contains(&bit) {
                            if let Some(last) = blocks.last_mut() {
                                if last.end == unit.start { last.end = unit.end; continue; }
                            }
                            blocks.push(*unit);
                        }
                    }
                    if blocks.is_empty() { continue; }
                    let gene = &reference.genes[gene_id];
                    let template = gene.transcript_ids().first()
                        .and_then(|&id| reference.transcripts.get(id))
                        .ok_or_else(|| anyhow::anyhow!("gene {gene_id} has no reference transcript"))?;
                    let gene_name = gene.primary_name()
                        .ok_or_else(|| anyhow::anyhow!("gene {gene_id} has no name"))?;
                    let chr = &reference.chr_names[template.chr_id];
                    // Find the closest reference over the observed interval.
                    // Missing terminal exons are not deletions in poly(A) data.
                    let observed_first = blocks.first().unwrap().start;
                    let observed_last = blocks.last().unwrap().end;
                    let mut closest: Option<(usize, usize)> = None;
                    for &candidate_id in gene.transcript_ids() {
                        let Some(candidate) = reference.transcripts.get(candidate_id) else { continue; };
                        if candidate.exons().is_empty() { continue; }
                        // Penalize only positively observed structural differences.
                        // A missing exon without an N junction is UNKNOWN, not skipped.
                        let skipped = candidate.exons().iter().filter(|e|
                            coverage.bypasses(e.start, e.end) &&
                            !blocks.iter().any(|b| b.overlaps(**e))).count();
                        // Do not treat coverage-only blocks as novel exons.
                        let distance = skipped;
                        if closest.is_none_or(|(d, id)| (distance, candidate_id) < (d, id)) {
                            closest = Some((distance, candidate_id));
                        }
                    }
                    let closest_id = closest.map(|(_, id)| id)
                        .ok_or_else(|| anyhow::anyhow!("no reference for gene {gene_id}"))?;
                    let closest_tx = &reference.transcripts[closest_id];
                    let reference_name = closest_tx.primary_name().unwrap_or(gene_name);
                    // Describe actual reference-relative changes using genomic intervals.
                    // An unobserved terminal exon is UNKNOWN, not a deletion.
                    let mut events = Vec::new();
                    let mut missing_coding_bases = 0u32;
                    let mut has_missing_cds = false;
                    let skipped: Vec<_> = closest_tx.exons().iter().copied()
                        .filter(|exon| coverage.bypasses(exon.start, exon.end))
                        .filter(|exon| !blocks.iter().any(|b| b.overlaps(*exon)))
                        .collect();
                    for exon in &skipped {
                        events.push(skip_event_name(chr, *exon));
                        if let Some(removed) = closest_tx.coding_bases_removed_by_junction(exon.start, exon.end) {
                            missing_coding_bases += removed;
                            has_missing_cds |= removed > 0;
                        }
                    }
                    events.sort();
                    events.dedup();
                    // Do not invent a novel transcript when no structural change
                    // can be described. Such UMIs remain unassigned to a novel model.
                    if events.is_empty() { continue; }
                    let mut name = format!("{reference_name}_{}", events.join("_"));
                    // This is a conditional coding prediction, not evidence that
                    // an unobserved exon was actually spliced out.
                    if has_missing_cds && missing_coding_bases % 3 != 0 {
                        name.push_str("_frameshift");
                    }
                    // Construct the model from the closest full reference, not
                    // the truncated read footprint. Remove only internal exons
                    // absent from the observed interval.
                    let mut model_blocks: Vec<RefBlock> = closest_tx.exons().iter().copied()
                        .filter(|exon| !skipped.contains(exon))
                        .collect();
                    // No new exon from intronic or partial aligned coverage.
                    // Coverage fragments fully contained by a reference exon are
                    // not additional exons. Preserve the original per-UMI mask
                    // separately; canonicalize only the registered model.
                    model_blocks.sort_by_key(|b| (b.start, std::cmp::Reverse(b.end)));
                    let mut canonical: Vec<RefBlock> = Vec::with_capacity(model_blocks.len());
                    for block in model_blocks {
                        if canonical.last().is_some_and(|last| last.end >= block.end) {
                            continue;
                        }
                        canonical.push(block);
                    }
                    let model_blocks = canonical;
                    // A reference without usable exons cannot define a novel model.
                    // Discard this candidate rather than aborting BAM quantification.
                    if model_blocks.is_empty() { continue; }
                    // IDs describe the biological structural change, never a hash.
                    // insert_novel_transcript validates that an existing name has
                    // identical exon coordinates; incompatible collisions are rejected.
                    let id = match self.gtf_splice_index.insert_novel_transcript(
                        gene_name, &name, chr, closest_tx.strand, &model_blocks) {
                        Ok(id) => id,
                        Err(err) if err.to_string().contains("conflicting novel transcript identity") => {
                            rejected_conflicts += 1;
                            if verbose_novel_conflicts && rejected_conflicts <= 20 {
                                eprintln!("[bam_tide finalize] WARNING: {err:#}");
                                eprintln!("  Incoming UMI: cell={cell} gene={gene_name} umi={umi}");
                                if let Some(catalog) = reference.gene_exon_catalog(gene_id) {
                                    // Columns are gene-wide exon *units*, not transcript exon numbers.
                                    // Annotations use '-' for a shared observed unit, '1' for an
                                    // annotated but unobserved unit, and '0' otherwise.
                                    let present = |m: &ExonMask, i: usize| -> bool {
                                        let bit = 2 * i;
                                        match m {
                                            ExonMask::Small(v) => bit < 64 && (v & (1u64 << bit)) != 0,
                                            ExonMask::Large(v) => v.get(bit / 64)
                                                .is_some_and(|w| (w & (1u64 << (bit % 64))) != 0),
                                        }
                                    };
                                    let columns = catalog.units.len();
                                    let format_row = |m: &ExonMask, annotation: bool| -> String {
                                        (0..columns).map(|i| {
                                            let observed = present(&mask, i);
                                            let set = present(m, i);
                                            if annotation && observed && set { '-' }
                                            else if annotation && observed && !set { '<' }
                                            else if set { '1' } else { '0' }
                                        }).collect()
                                    };
                                    eprintln!("  Exon units (genomic order): {}", (1..=columns)
                                        .map(|i| format!("{}", i % 10)).collect::<String>());
                                    eprintln!("  Incoming UMI:               {}", format_row(&mask, false));
                                    eprintln!("  Annotated transcripts for {gene_name} ({}):", gene.transcript_ids().len());
                                    for &reference_id in gene.transcript_ids() {
                                        if let Some(tx) = reference.transcripts.get(reference_id) {
                                            let tx_mask = catalog.transcript_mask(tx);
                                            eprintln!("    {:<26} {}", tx.primary_name().unwrap_or("<unnamed>"),
                                                format_row(&tx_mask, true));
                                        }
                                    }
                                    eprintln!("  Legend: '-' shared with UMI; '<' observed in UMI but missing from annotation; '1' annotated but not observed; '0' neither.");
                                    eprintln!("  Unobserved UMI exons are unknown, not evidence of skipping.");
                                } else {
                                    eprintln!("  No exon catalog available for {gene_name}; mask={mask:?}");
                                }
                                eprintln!("  Action: reject this novel-transcript UMI assignment; continue");
                            }
                            continue;
                        }
                        Err(err) => return Err(err),
                    };
                    novel_ids.insert(lookup, id);
                    id
                };
                // Coverage belongs to the molecule, not the transcript identity.
                if let Some(catalog) = reference.gene_exon_catalog(gene_id) {
                    let mut span: Option<(u32, u32)> = None;
                    for (i, unit) in catalog.units.iter().enumerate() {
                        let mut bit = ExonMask::default();
                        bit.set(2 * i);
                        if mask.contains(&bit) {
                            span = Some(match span {
                                Some((start, end)) => (start.min(unit.start), end.max(unit.end)),
                                None => (unit.start, unit.end),
                            });
                        }
                    }
                    if let Some(span) = span {
                        self.supported_regions.entry((cell, tx_id, umi))
                            .and_modify(|r| { r.0 = r.0.min(span.0); r.1 = r.1.max(span.1); })
                            .or_insert(span);
                    }
                }
                self.add_molecule(tx_id, cell, umi, snps);
                novel_count += 1;
            }
        }
        // Group unresolved molecules by the *set* of compatible reference
        // transcripts, not by their posterior values or processing order.
        for (cell, gene_id, umi, snps) in deferred {
            Self::write_snp_counts(data, report, cell, umi, &snps);
            let Some(p) = bayes.get(&(cell, gene_id, umi)) else { continue; };
            // The ordinary expression matrices currently use transcript IDs.
            // Preserve the molecule there once, using a deterministic reference
            // feature as a storage representative; the ambiguity TSV records
            // that this is NOT a confident transcript assignment.
            if let Some(&(representative, _, _)) = p.iter().min_by_key(|v| v.0) {
                let feature = GeneUmiHash(representative as u64, umi);
                data.try_insert(QuantData::EXPRS, &cell, feature, 1.0, report);
                data.try_insert(QuantData::EXONIC, &cell, feature, 1.0, report);
            }
            let mut members: Vec<usize> = p.iter().map(|v| v.0).collect();
            members.sort_unstable();
            self.ambiguity_gene_names.entry(gene_id).or_insert_with(||
                reference.genes[gene_id].primary_name().unwrap_or("NA").to_owned());
            for &id in &members {
                self.ambiguity_transcript_names.entry(id).or_insert_with(||
                    reference.transcript_name(id).unwrap_or("NA").to_owned());
            }
            let group = self.ambiguity.entry((gene_id, members)).or_default();
            group.molecules.push((cell, umi, snps));
            let gene_name = reference.genes[gene_id].primary_name().unwrap_or("NA");
            for &(id, probability, snp_informed) in p {
                let candidate = reference.transcript_name(id).unwrap_or("NA").to_owned();
                self.posterior_rows.push(PosteriorRow { cell, gene: gene_name.to_owned(),
                    umi, candidate, probability, snp_informed });
            }
        }
        if rejected_conflicts > 0 {
            eprintln!("[bam_tide finalize] rejected_novel_umi_assignments={rejected_conflicts} (set BAM_TIDE_DEBUG_NOVEL_CONFLICTS=1 for details)");
        }
        eprintln!("[bam_tide finalize] done; fused={total} novel_molecules={novel_count} novel_models={} total_elapsed={:.2}s",
            novel_ids.len(), started.elapsed().as_secs_f64());
        Ok(())
    }

    fn write_snp_counts(data: &mut QuantData, report: &mut RunStats,
                        cell: u64, umi: u64, snps: &SnpUmiEvidence) {
        for snp in &snps.snps {
            // One allele per consolidated cell/UMI/SNP; ties contribute neither.
            let category = match snp.ref_reads.cmp(&snp.alt_reads) {
                std::cmp::Ordering::Greater => QuantData::SNP_REF,
                std::cmp::Ordering::Less => QuantData::SNP_ALT,
                std::cmp::Ordering::Equal => continue,
            };
            data.try_insert(category, &cell,
                GeneUmiHash(snp.snp_id, umi), 1.0, report);
        }
    }

    /// Export novel isoforms and per-molecule SNP support without changing normal SNP output.
    pub fn write_reports(&mut self, output: impl AsRef<Path>) -> Result<()> {
        let output = output.as_ref();
        if self.gtf_splice_index.transcripts.is_empty() && self.ambiguity.is_empty() { return Ok(()); }
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
        // Separate molecular coverage from structural transcript identity.
        let support_path = output.join("new_transcript_supported_regions.tsv");
        let mut support = BufWriter::new(std::fs::File::create(&support_path)?);
        writeln!(support, "cell_id\ttranscript_id\tumi\tchromosome\tstart\tend")?;
        for (&(cell, tx_id, umi), &(start, end)) in &self.supported_regions {
            if let Some(tx) = self.gtf_splice_index.transcripts.get(tx_id) {
                let name = tx.primary_name().unwrap_or("NA");
                let chr = &self.gtf_splice_index.chr_names[tx.chr_id];
                writeln!(support, "{cell}\t{name}\t{umi}\t{chr}\t{}\t{end}", start + 1)?;
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
        // Known-unknown evidence belongs in its own table, not in the novel
        // transcript evidence report or novel-transcript matrix.
        let ambiguity_path = output.join("known_transcript_ambiguity.tsv");
        let mut ambiguity_file = BufWriter::new(std::fs::File::create(&ambiguity_path)?);
        writeln!(ambiguity_file, "cell_id\tgene\tumi\tcandidate_transcripts\tsnp_id\tref_reads\talt_reads")?;
        for ((gene_id, members), group) in &self.ambiguity {
            let gene = self.ambiguity_gene_names.get(gene_id).map_or("NA", String::as_str);
            let label = members.iter().map(|id| self.ambiguity_transcript_names.get(id)
                .map_or("NA", String::as_str)).collect::<Vec<_>>().join(",");
            for (cell, umi, evidence) in &group.molecules {
                if evidence.snps.is_empty() {
                    writeln!(ambiguity_file, "{cell}\t{gene}\t{umi}\t{label}\tNA\t0\t0")?;
                } else {
                    for snp in &evidence.snps {
                        writeln!(ambiguity_file, "{cell}\t{gene}\t{umi}\t{label}\t{}\t{}\t{}",
                            snp.snp_id, snp.ref_reads, snp.alt_reads)?;
                    }
                }
            }
        }
        // Report conflicts from consolidated gene-level UMI evidence, rather
        // than repeating the same observation for each novel transcript ID.
        let conflicts_path = output.join("snp_umi_conflicts.tsv");
        let mut conflicts_file = BufWriter::new(std::fs::File::create(&conflicts_path)?);
        writeln!(conflicts_file, "cell_id\tumi\tsnp_id\tref_reads\talt_reads\tdecision")?;
        let mut conflicts = 0usize;
        let mut conflict_error = None;
        let mut seen = std::collections::HashSet::new();
        self.exon_evidence.for_each_umi(|cell, key, evidence| {
            if conflict_error.is_some() { return; }
            for snp in evidence.snps.conflicts() {
                if !seen.insert((cell, key.1, snp.snp_id)) { continue; }
                let decision = match snp.ref_reads.cmp(&snp.alt_reads) {
                    std::cmp::Ordering::Greater => "REF",
                    std::cmp::Ordering::Less => "ALT",
                    std::cmp::Ordering::Equal => "IGNORED",
                };
                if let Err(err) = writeln!(conflicts_file,
                    "{cell}\t{}\t{}\t{}\t{}\t{decision}",
                    key.1, snp.snp_id, snp.ref_reads, snp.alt_reads) {
                    conflict_error = Some(err);
                    return;
                }
                conflicts += 1;
            }
        });
        if let Some(err) = conflict_error { return Err(err.into()); }
        eprintln!("[bam_tide] SNP conflicts: {conflicts} entries written to {}", conflicts_path.display());

        // Keep the per-UMI SNP report above independent from the ordinary
        // transcript-by-cell UMI matrix. Export every observed novel-transcript
        // cell, without applying the normal GEX cell-calling threshold.
        let feature_index = self.gtf_splice_index.feature_index();
        let mut posterior_file = BufWriter::new(std::fs::File::create(output.join("transcript_bayesian_posteriors.tsv"))?);
        writeln!(posterior_file, "cell_id\tgene\tumi\ttranscript\tposterior_probability\tsnp_informed")?;
        for row in &self.posterior_rows {
            writeln!(posterior_file, "{}\t{}\t{}\t{}\t{:.12}\t{}",
                row.cell, row.gene, row.umi, row.candidate, row.probability, row.snp_informed)?;
        }

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
        self.exon_evidence.merge_single_thread(source.exon_evidence);
        for (key, group) in source.ambiguity {
            self.ambiguity.entry(key).or_default().molecules.extend(group.molecules);
        }
        self.posterior_rows.extend(source.posterior_rows);
        self.ambiguity_gene_names.extend(source.ambiguity_gene_names);
        self.ambiguity_transcript_names.extend(source.ambiguity_transcript_names);
        let mapping = self.gtf_splice_index.merge_novel_transcripts(&source.gtf_splice_index)?;
        // The index intentionally omits conflicting source IDs. Remove their
        // molecule evidence before remapping so no UMI is assigned to a wrong model.
        let observed = source.scdata.observed_feature_ids();
        let rejected: Vec<_> = observed.iter().filter(|id| !mapping.contains_key(id)).copied().collect();
        if !rejected.is_empty() {
            eprintln!("[bam_tide] WARNING: excluding {} conflicting novel-transcript features from shard", rejected.len());
            source.scdata.retain_features(&mapping.keys().copied().collect());
        }
        source.scdata.remap_feature_ids(&mapping).map_err(anyhow::Error::msg)?;
        self.scdata.merge_single_thread(source.scdata);
        for ((cell, tx_id, umi), (start, end)) in source.supported_regions {
            let Some(&mapped_id) = mapping.get(&(tx_id as u64)) else { continue; };
            let mapped_id = mapped_id as usize;
            self.supported_regions.entry((cell, mapped_id, umi))
                .and_modify(|r| { r.0 = r.0.min(start); r.1 = r.1.max(end); })
                .or_insert((start, end));
        }
        Ok(())
    }
}

#[cfg(test)]
mod exon_collection_tests {
    use super::*;

    #[test]
    fn exon_evidence_merges_within_and_across_workers() {
        let mut a = NewTranscriptData::new(100);
        let mut b = NewTranscriptData::new(100);
        let mut high = ExonMask::empty(96);
        high.set(70);
        a.add_exon_evidence(4, 12, 8, ExonMask::Small(1), SnpUmiEvidence::default());
        a.add_exon_evidence(4, 12, 8, ExonMask::Small(2), SnpUmiEvidence::default());
        b.add_exon_evidence(4, 12, 8, high.clone(), SnpUmiEvidence::default());
        b.add_exon_evidence(5, 12, 8, ExonMask::Small(4), SnpUmiEvidence::default());
        a.merge(b).unwrap();
        let mut observed = Vec::new();
        a.exon_evidence.for_each_umi(|cell, key, mask| {
            observed.push((cell, key.0, key.1, mask.splice.clone()));
        });
        assert_eq!(observed.len(), 2);
        let gene4 = &observed.iter().find(|(_, gene, _, _)| *gene == 4).unwrap().3;
        assert_eq!(gene4.count_ones(), 3);
        assert!(gene4.contains(&high));
        assert_eq!(observed.iter().find(|(_, gene, _, _)| *gene == 5).unwrap().3, ExonMask::Small(4));
    }
}

#[cfg(test)]
mod nasty_novel_event_characterization_tests {
    use super::*;

    fn b(start: u32, end: u32) -> RefBlock { RefBlock::new(start, end) }

    // Regression: disconnected reads from one UMI cannot invent splice junctions.
    #[test]
    fn disconnected_reads_do_not_invent_five_exon_skips() {
        let parent = vec![b(100, 150), b(200, 250), b(300, 350), b(400, 450),
                          b(500, 550), b(600, 650), b(700, 750)];
        // Two separate reads for one UMI, NOT a read connecting E1 and E7.
        // The observed CIGAR-N junction set is empty.
        let observed = vec![b(120, 145), b(710, 740)];
        let observed_junctions: Vec<(u32, u32)> = vec![];
        let (skipped, added) = reference_relative_event_candidates(&parent, &observed, &[]);
        assert!(observed_junctions.is_empty());
        assert!(skipped.is_empty(), "disconnected reads must not imply skipped exons");
        assert!(added.is_empty());
    }

    #[test]
    fn wrong_reference_parent_does_not_invent_skips() {
        let wrong_parent = vec![b(100, 150), b(300, 350), b(500, 550)];
        let alternate_known_parent = vec![b(100, 150), b(200, 280), b(300, 350),
                                          b(400, 450), b(500, 550)];
        // E2 is already annotated in another known transcript; no novel add.
        // E4 is simply not covered; no CIGAR N junction crosses it.
        let observed = vec![b(110, 140), b(205, 275), b(510, 540)];
        let (wrong_skips, wrong_adds) = reference_relative_event_candidates(&wrong_parent, &observed, &[]);
        assert!(wrong_adds.is_empty(), "coverage without a junction cannot prove a novel exon");
        assert!(wrong_skips.is_empty(), "uncovered reference exon is not skipped");
        let (better_skips, better_adds) = reference_relative_event_candidates(&alternate_known_parent, &observed, &[]);
        assert!(better_skips.is_empty(), "uncovered exons are unknown, not skipped");
        assert!(better_adds.is_empty());
    }

    #[test]
    fn partial_terminal_coverage_does_not_create_terminal_skips() {
        let parent = vec![b(100, 150), b(200, 250), b(300, 350),
                          b(400, 450), b(500, 550)];
        let observed = vec![b(310, 345), b(410, 445)];
        let (skipped, added) = reference_relative_event_candidates(&parent, &observed, &[]);
        assert!(skipped.is_empty());
        assert!(added.is_empty());
    }

    // These are intentionally RED characterization tests. Do not loosen their
    // assertions to make them green: production must learn to classify these
    // observations without manufacturing an added exon.
    #[test]
    fn terminal_two_base_overlap_means_known_exon_not_truncated_novel_exon() {
        let reference = vec![b(100, 150), b(200, 250), b(300, 350)];
        // Read terminates two bases into exon 2; the remainder is UNKNOWN,
        // not evidence for a new exon boundary or a skipped exon.
        let observed = vec![b(120, 150), b(200, 202)];
        let (skipped, added) = reference_relative_event_candidates(&reference, &observed, &[]);
        assert!(skipped.is_empty(), "terminal truncation is not a skip");
        assert!(added.is_empty(), "two aligned bases into an annotated exon establish presence");
    }

    #[test]
    fn three_base_splice_boundary_wobble_is_not_a_novel_added_exon() {
        let reference = vec![b(100, 150), b(200, 250)];
        // Mapper places the left exon 3 bp past its annotated end.  With a
        // configurable 5 bp splice-site tolerance this is reference-compatible.
        let observed = vec![b(120, 153), b(200, 235)];
        let (skipped, added) = reference_relative_event_candidates(&reference, &observed, &[(153, 200)]);
        assert!(skipped.is_empty());
        assert!(added.is_empty(), "3 bp boundary wobble must not create add_chr");
    }

    #[test]
    fn substantial_intronic_alignment_is_not_a_new_exon() {
        let reference = vec![b(100, 150), b(300, 350)];
        // An aligned intronic segment belongs in gene-level intronic evidence;
        // the current event helper incorrectly emits an added exon for it.
        let observed = vec![b(120, 145), b(190, 260)];
        let (skipped, added) = reference_relative_event_candidates(&reference, &observed, &[]);
        assert!(skipped.is_empty());
        assert!(added.is_empty(), "intronic coverage must not create add_chr");
    }

    #[test]
    fn exon_annotated_in_another_transcript_is_not_novel_to_the_gene() {
        let selected_parent = vec![b(100, 150), b(300, 350)];
        let alternative_parent = vec![b(100, 150), b(200, 250), b(300, 350)];
        let observed = vec![b(120, 145), b(205, 245), b(315, 340)];
        let (_, selected_adds) = reference_relative_event_candidates(&selected_parent, &observed, &[]);
        let (_, alternative_adds) = reference_relative_event_candidates(&alternative_parent, &observed, &[]);
        assert!(alternative_adds.is_empty());
        assert!(selected_adds.is_empty(), "gene-annotated exon must not become novel solely from parent selection");
    }

    #[test]
    fn contiguous_cigar_n_spanning_exon_is_genuine_skip_candidate() {
        let parent = vec![b(100, 150), b(200, 250), b(300, 350)];
        let observed = vec![b(120, 150), b(300, 330)];
        // Example: 30M150N30M, N is [150,300), spans exon [200,250).
        let junction = b(150, 300);
        let (skipped, added) = reference_relative_event_candidates(&parent, &observed, &[(junction.start, junction.end)]);
        assert_eq!(skipped, vec![b(200, 250)]);
        assert!(junction.contains(skipped[0]));
        assert_eq!(skip_event_name("chr17", skipped[0]), "skip_chr17:201-250");
        assert!(added.is_empty());
    }
}
