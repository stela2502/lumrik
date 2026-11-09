//! An independent annotation + UMI matrix for newly observed splice isoforms.
use anyhow::{Result, Context};
use std::io::{BufWriter, Write};
use std::path::Path;
use gtf_splice_index::{SpliceIndex, SpliceMatchMode, UmiTranscriptMatch, RefBlock};
use super::cli::QuantMode;
use run_stats::RunStats;
use scdata::QuantData;
use std::collections::HashMap;
use std::time::Instant;
use scdata::{MatrixValueType, Scdata, GeneUmiHash, SnpUmiEvidence, ExonMask, TranscriptUmiEvidence};

pub struct NewTranscriptData {
    pub scdata: Scdata<SnpUmiEvidence>,
    /// Raw gene-wide exon evidence, keyed by (cell, gene, UMI).
    /// Transcript assignment must happen only after collection.
    pub exon_evidence: Scdata<TranscriptUmiEvidence>,
    pub gtf_splice_index: SpliceIndex,
    /// Genomic span supported by each novel-model UMI (0-based, half-open).
    supported_regions: HashMap<(u64, usize, u64), (u32, u32)>,

}

impl NewTranscriptData {
    pub fn new(bin_width: u32) -> Self {
        Self {
            scdata: Scdata::new_with_evidence(1, MatrixValueType::Integer),
            exon_evidence: Scdata::new_with_evidence(1, MatrixValueType::Integer),
            gtf_splice_index: SpliceIndex::new(bin_width).with_match_mode(SpliceMatchMode::Transcript),
            supported_regions: HashMap::new(),
        }
    }

    pub fn add_molecule(&mut self, transcript_id: usize, cell: u64, umi: u64, evidence: SnpUmiEvidence) {
        self.scdata.try_insert_evidence(&cell, GeneUmiHash(transcript_id as u64, umi), evidence);
    }

    pub fn add_exon_evidence(&mut self, gene_id: usize, cell: u64, umi: u64, mask: ExonMask, snps: SnpUmiEvidence) {
        if mask.count_ones() > 0 {
            self.exon_evidence.try_insert_evidence(&cell, GeneUmiHash(gene_id as u64, umi),
                TranscriptUmiEvidence { splice: mask, snps });
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

        // Each worker owns its matrices. No global lock and no serial per-UMI
        // insertion. Novel patterns retain their original SNP evidence until
        // the shared novel index has assigned stable transcript IDs.
        let processed = std::sync::atomic::AtomicUsize::new(0);
        let next_report = std::sync::atomic::AtomicUsize::new(1_000_000);
        let shards: Vec<_> = molecules
            .into_par_iter()
            .fold(
                || (QuantData::standard(), RunStats::new(None, 0.0, usize::MAX),
                    Vec::<(u64, usize, u64, ExonMask, SnpUmiEvidence)>::new()),
                |(mut local, mut stats, mut novel), (cell, key, evidence)| {
                    let gene_id = key.0 as usize;
                    let umi = key.1;
                    let resolution = if transcript_mode {
                        reference.resolve_umi_transcript(gene_id, &evidence.splice)
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
                                    novel.push((cell, gene_id, umi, mask, evidence.snps.clone()));
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
                    (local, stats, novel)
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
        for (local, _stats, novel) in shards {
            report.merge(&data.merge(&local));
            for (cell, gene_id, umi, mask, snps) in novel {
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
                        let reference_mask = catalog.transcript_mask(candidate);
                        let mut distance = 0usize;
                        for i in 0..catalog.units.len() {
                            let unit = catalog.units[i];
                            if unit.end <= observed_first || unit.start >= observed_last { continue; }
                            let mut bit = ExonMask::default();
                            bit.set(2 * i);
                            if mask.contains(&bit) != reference_mask.contains(&bit) { distance += 1; }
                        }
                        if closest.is_none_or(|(d, id)| (distance, candidate_id) < (d, id)) {
                            closest = Some((distance, candidate_id));
                        }
                    }
                    let closest_id = closest.map(|(_, id)| id)
                        .ok_or_else(|| anyhow::anyhow!("no reference for gene {gene_id}"))?;
                    let closest_tx = &reference.transcripts[closest_id];
                    let reference_name = closest_tx.primary_name().unwrap_or(gene_name);
                    // Only *internal* missing reference exons contribute an event.
                    // The current bitmap does not retain splice-junction connectivity;
                    // call these candidates, never confirmed deletions.
                    let mut events = Vec::new();
                    for exon in closest_tx.exons() {
                        if exon.start <= observed_first || exon.end >= observed_last { continue; }
                        if !blocks.iter().any(|b| b.overlaps(*exon)) {
                            events.push(format!("{chr}:{}-{}_internal_missing", exon.start + 1, exon.end));
                        }
                    }
                    // Additional observed regions are coverage, not structural events.
                    // Keep them in per-UMI supported_regions instead of transcript IDs.
                    events.sort();
                    events.dedup();
                    if events.is_empty() {
                        events.push("unresolved".to_owned());
                    }
                    let name = format!("{reference_name}_mod_{}", events.join("_"));
                    // Construct the model from the closest full reference, not
                    // the truncated read footprint. Remove only internal exons
                    // absent from the observed interval.
                    let mut model_blocks: Vec<RefBlock> = closest_tx.exons().iter().copied()
                        .filter(|exon| exon.start <= observed_first || exon.end >= observed_last ||
                            blocks.iter().any(|b| b.overlaps(*exon)))
                        .collect();
                    for block in &blocks {
                        if !model_blocks.iter().any(|exon| exon.contains(*block)) {
                            model_blocks.push(*block);
                        }
                    }
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
                    // Event labels alone are not unique transcript identities:
                    // distinct exon structures can have the same `_observed` event.
                    // Include a deterministic fingerprint of the full structure,
                    // without changing the human-readable genomic event prefix.
                    let mut fingerprint = 0xcbf29ce484222325u64;
                    for value in std::iter::once(gene_id as u64).chain(
                        model_blocks.iter().flat_map(|b| [b.start as u64, b.end as u64])
                    ) {
                        for byte in value.to_le_bytes() {
                            fingerprint ^= byte as u64;
                            fingerprint = fingerprint.wrapping_mul(0x100000001b3);
                        }
                    }
                    let name = format!("{name}_structure_{fingerprint:016x}");
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
