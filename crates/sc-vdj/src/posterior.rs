use crate::align::{
    best_oriented_local_alignment, best_oriented_local_alignment_with_reverse, local_alignment,
    OrientedLocalAlignment,
};
use crate::gex::ExpressionMatrix;
use crate::junction::{measure_junction, JunctionInput, JunctionMeasurement};
use crate::mapper::{CoverageFamilyScore, IdentitySeedScore, VdjMapper};
use crate::reference::VdjReference;
use crate::score::{score_recombination_activity, RecombinationActivityEvidence};
use crate::sterile::{SterileAccumulator, SterileProfile};
use crate::types::{Chain, SegmentKind};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct BamReadEvidence {
    pub cell: String,
    pub umi: String,
    pub read_name: String,
    pub sequence: Vec<u8>,
    /// Raw BAM Phred qualities, one byte per base in `sequence`.
    pub qualities: Vec<u8>,
    pub chr: Option<String>,
    pub ref_start: Option<u32>,
    pub ref_end: Option<u32>,
    /// Exact aligned reference blocks; RefSkip introns are not counted as transcriptional coverage.
    pub ref_blocks: Vec<(u32, u32)>,
    pub mapq: u8,
    pub is_reverse: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
}

#[derive(Debug, Clone)]
pub struct GermlineSegmentSupport {
    pub segment_index: usize,
    pub id: String,
    pub kind: SegmentKind,
    pub local_alignment_score: i32,
    pub supporting_umis: usize,
    pub supporting_reads: usize,
    pub locus_fraction: f64,
    pub distance_to_recombination_center: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecombinationStage {
    None,
    Dj,
    Vj,
    Vdj,
}

#[derive(Debug, Clone)]
pub struct ReadSegmentAlignment {
    pub score: i32,
    /// Coordinates on the read in transcript orientation. When `reverse_complement`
    /// is true these coordinates refer to the reverse-complemented BAM sequence.
    pub query_start: usize,
    pub query_end: usize,
    pub reference_start: usize,
    pub reference_end: usize,
    pub reverse_complement: bool,
}

#[derive(Debug, Clone)]
pub struct RearrangementSupportingRead {
    pub umi: String,
    pub read_name: String,
    /// Sequence exactly as stored in BAM. It is not silently reverse-complemented.
    pub sequence: Vec<u8>,
    /// Raw BAM Phred qualities corresponding to `sequence`.
    pub qualities: Vec<u8>,
    pub bam_is_reverse: bool,
    pub is_supplementary: bool,
    pub supports_v: bool,
    pub supports_j: bool,
    pub supports_d: bool,
    pub supports_c: bool,
    pub v_alignment: Option<ReadSegmentAlignment>,
    pub j_alignment: Option<ReadSegmentAlignment>,
    pub d_alignment: Option<ReadSegmentAlignment>,
    pub c_alignment: Option<ReadSegmentAlignment>,
}

#[derive(Debug, Clone)]
pub struct RearrangementCall {
    pub chain: Chain,
    pub stage: RecombinationStage,
    pub v: Option<GermlineSegmentSupport>,
    pub d: Option<GermlineSegmentSupport>,
    /// For D-bearing chains, true when the selected D was inferred by comparing
    /// all locus-compatible D genes only inside the already-established V/J junction.
    pub d_inferred_from_vj_junction: bool,
    /// Score difference between the best and second-best bounded D hypotheses.
    /// This is deliberately diagnostic rather than an acceptance threshold.
    pub d_hypothesis_margin: Option<i32>,
    pub j: Option<GermlineSegmentSupport>,
    pub c: Option<GermlineSegmentSupport>,
    /// UMIs with coherent support for both the selected V and J segments.
    pub total_supporting_umis: usize,
    /// Original BAM reads from those coherent UMIs that support at least one
    /// segment in the final call. Used only for auditable sequence diagnostics.
    pub supporting_reads: Vec<RearrangementSupportingRead>,
    /// Junction geometry measured on one continuous coherent UMI contig.
    pub junction: Option<JunctionMeasurement>,
    pub notation: String,
}

#[derive(Debug, Clone)]
pub struct CellVdjSummary {
    pub cell: String,
    pub rearrangements: Vec<RearrangementCall>,
    pub sterile: Vec<SterileProfile>,
    pub recombination_activity: RecombinationActivityEvidence,
}

#[derive(Debug, Clone)]
pub struct PosteriorConfig {
    pub sterile_bins: usize,
    pub min_seed_hits: u32,
    pub candidate_segments_per_kind: usize,
    /// Maximum sensitive 9-mer candidates retained per recognized chain/kind
    /// family as a safety fallback around the discriminative identity index.
    pub coverage_candidates_per_family_kind: usize,
    /// Maximum discriminative 13-mer candidates retained per recognized
    /// chain/kind family.
    pub identity_candidates_per_family_kind: usize,
    /// Minimum distinct 13-mer observations before identity evidence nominates
    /// a germline segment.
    pub identity_min_hits: u32,
    /// Dedicated V-3'/J-5' rescue candidates retained per segment kind.
    pub terminal_candidates_per_kind: usize,
    /// Minimum distinct terminal 13-mers required for junction rescue.
    pub terminal_min_hits: u32,
    /// Maximum chain/locus families retained independently for each segment kind.
    pub family_candidates_per_kind: usize,
    /// Minimum relative-position coverage bins for V/J/C recognition.
    pub min_v_coverage_bins: u32,
    pub min_j_coverage_bins: u32,
    pub min_c_coverage_bins: u32,
    pub min_local_score: i32,
    /// Reject C candidates implausibly far from the chain's J recombination center.
    pub max_constant_distance_bp: u64,
    /// Maximum V/J overlap on one oriented read. Large overlap usually means
    /// the two germlines are explaining the same conserved sequence.
    pub max_vj_alignment_overlap: usize,
    /// V evidence must reach this close to the 3' end of the germline V.
    pub max_v_end_distance: usize,
    /// J evidence must begin this close to the 5' end of the germline J.
    pub max_j_start_distance: usize,
    /// Minimum overlap used to assemble raw observations from one UMI into
    /// molecule-level contigs before germline alignment.
    pub min_umi_contig_overlap: usize,
    /// Minimum observed sequence overlap when separate V- and J-anchored reads
    /// from the same UMI are used to establish a continuous V->J path.
    pub min_vj_read_overlap: usize,
}
impl Default for PosteriorConfig {
    fn default() -> Self {
        Self {
            sterile_bins: 64,
            min_seed_hits: 2,
            candidate_segments_per_kind: 8,
            coverage_candidates_per_family_kind: 2,
            identity_candidates_per_family_kind: 3,
            identity_min_hits: 1,
            terminal_candidates_per_kind: 2,
            terminal_min_hits: 2,
            family_candidates_per_kind: 2,
            // Distributed family coverage stays conservative. Dedicated V-3'/J-5'
            // terminal rescue seeds recover short rearrangement junction fragments.
            // V 3' / J 5' bins, so requiring two V bins discards exactly the
            // short junction-spanning molecules needed for low-support chains.
            // C deliberately remains distributed/strict to prevent constant
            // region support inflation.
            min_v_coverage_bins: 2,
            min_j_coverage_bins: 1,
            min_c_coverage_bins: 2,
            min_local_score: 18,
            max_constant_distance_bp: 1_000_000,
            max_vj_alignment_overlap: 12,
            max_v_end_distance: 35,
            max_j_start_distance: 35,
            min_umi_contig_overlap: 20,
            min_vj_read_overlap: 20,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct BestReadMatch {
    segment_index: usize,
    seed_hits: u32,
    local_score: i32,
    alignment: OrientedLocalAlignment,
}

struct UmiContig<'a> {
    umi: &'a str,
    sequence: Vec<u8>,
    members: Vec<&'a BamReadEvidence>,
}

struct SeededRead<'a> {
    umi: &'a str,
    sequence: Vec<u8>,
    members: Vec<&'a BamReadEvidence>,
    multiplicity: usize,
    reverse_sequence: Vec<u8>,
    seed_scores: HashMap<usize, u32>,
    best_by_kind: HashMap<SegmentKind, BestReadMatch>,
    alignments: HashMap<usize, OrientedLocalAlignment>,
    candidate_alignments: usize,
    family_candidates: HashMap<SegmentKind, Vec<Chain>>,
    identity_seed_used: bool,
}

#[derive(Debug, Clone)]
struct SelectedSegment {
    support: GermlineSegmentSupport,
    umis: HashSet<String>,
}

pub struct PosteriorAnalyzer<'a> {
    reference: &'a VdjReference,
    mapper: &'a VdjMapper,
    config: PosteriorConfig,
}
impl<'a> PosteriorAnalyzer<'a> {
    pub fn new(
        reference: &'a VdjReference,
        mapper: &'a VdjMapper,
        config: PosteriorConfig,
    ) -> Self {
        Self {
            reference,
            mapper,
            config,
        }
    }

    pub fn analyze<E: ExpressionMatrix + Sync>(
        &self,
        reads: impl IntoIterator<Item = BamReadEvidence>,
        gex: &E,
    ) -> Vec<CellVdjSummary> {
        self.analyze_for_cells(reads, gex, std::iter::empty::<String>())
    }

    /// Analyze a bounded batch while optionally predeclaring cells that must be
    /// emitted even if receptor routing retained no BAM records for them.
    /// Developmental/recombination GEX scoring remains meaningful for an
    /// otherwise receptor-silent cell.
    pub fn analyze_for_cells<E, I>(
        &self,
        reads: impl IntoIterator<Item = BamReadEvidence>,
        gex: &E,
        required_cells: I,
    ) -> Vec<CellVdjSummary>
    where
        E: ExpressionMatrix + Sync,
        I: IntoIterator<Item = String>,
    {
        let mut ignored = MappingInfo::new(None, 0.0, 0);
        self.analyze_for_cells_with_mapping_info(reads, gex, required_cells, &mut ignored)
    }

    /// Analyze a bounded cell batch and merge detailed per-cell work counters and
    /// timings into the caller's MappingInfo. `cell_work.*` named timings are
    /// summed across Rayon workers and therefore intentionally represent total
    /// parallel work, not wall-clock latency.
    pub fn analyze_for_cells_with_mapping_info<E, I>(
        &self,
        reads: impl IntoIterator<Item = BamReadEvidence>,
        gex: &E,
        required_cells: I,
        report: &mut MappingInfo,
    ) -> Vec<CellVdjSummary>
    where
        E: ExpressionMatrix + Sync,
        I: IntoIterator<Item = String>,
    {
        let mut cells: HashMap<String, Vec<BamReadEvidence>> = HashMap::new();
        for cell in required_cells {
            cells.entry(cell).or_default();
        }
        for r in reads {
            if !r.is_secondary {
                cells.entry(r.cell.clone()).or_default().push(r);
            }
        }
        report.report_n("vdj.cells", cells.len());

        let mut out: Vec<_> = cells
            .into_par_iter()
            .map(|(cell, reads)| {
                let mut local = MappingInfo::new(None, 0.0, 0);
                let summary = self.analyze_cell_profiled(cell, reads, gex, &mut local);
                (summary, local)
            })
            .collect();
        for (_, local) in &out {
            report.merge(local);
        }
        out.sort_by(|a, b| a.0.cell.cmp(&b.0.cell));
        out.into_iter().map(|(summary, _)| summary).collect()
    }

    /// Analyze evidence that was already routed by indexed genomic locus during
    /// the single BAM scan. Each receptor locus is reconstructed independently;
    /// calls from another locus are never allowed to leak out of that bucket.
    pub fn analyze_routed_for_cells_with_mapping_info<E, I>(
        &self,
        routed: crate::bam::RoutedBamEvidence,
        gex: &E,
        required_cells: I,
        report: &mut MappingInfo,
    ) -> Vec<CellVdjSummary>
    where
        E: ExpressionMatrix + Sync,
        I: IntoIterator<Item = String>,
    {
        self.analyze_routed_for_cells_with_progress(
            routed,
            gex,
            required_cells,
            report,
            |_, _| {},
        )
    }

    /// Progress-enabled locus-routed analysis. The callback is invoked from
    /// Rayon worker threads after each cell is complete, so callers can expose
    /// genuine live reconstruction progress without creating artificial shards.
    pub fn analyze_routed_for_cells_with_progress<E, I, G>(
        &self,
        mut routed: crate::bam::RoutedBamEvidence,
        gex: &E,
        required_cells: I,
        report: &mut MappingInfo,
        progress: G,
    ) -> Vec<CellVdjSummary>
    where
        E: ExpressionMatrix + Sync,
        I: IntoIterator<Item = String>,
        G: Fn(&CellVdjSummary, &MappingInfo) + Sync,
    {
        for cell in required_cells {
            routed.entry(cell).or_default();
        }
        report.report_n("vdj.cells", routed.len());

        let mut out: Vec<_> = routed
            .into_par_iter()
            .map(|(cell, loci)| {
                let mut local = MappingInfo::new(None, 0.0, 0);
                let mut rearrangements = Vec::new();
                let mut sterile = Vec::new();

                for (chain, reads) in loci {
                    let summary = self.analyze_cell_profiled(
                        cell.clone(),
                        reads,
                        gex,
                        &mut local,
                    );
                    rearrangements.extend(
                        summary
                            .rearrangements
                            .into_iter()
                            .filter(|call| call.chain == chain),
                    );
                    sterile.extend(
                        summary
                            .sterile
                            .into_iter()
                            .filter(|profile| profile.chain == chain),
                    );
                }

                rearrangements.sort_by_key(|call| call.chain);
                sterile.sort_by_key(|profile| profile.chain);
                let summary = CellVdjSummary {
                    recombination_activity: score_recombination_activity(gex, &cell),
                    cell,
                    rearrangements,
                    sterile,
                };
                progress(&summary, &local);
                (summary, local)
            })
            .collect();
        for (_, local) in &out {
            report.merge(local);
        }
        out.sort_by(|a, b| a.0.cell.cmp(&b.0.cell));
        out.into_iter().map(|(summary, _)| summary).collect()
    }

    fn analyze_cell_profiled<E: ExpressionMatrix>(
        &self,
        cell: String,
        reads: Vec<BamReadEvidence>,
        gex: &E,
        report: &mut MappingInfo,
    ) -> CellVdjSummary {
        report.report_n("vdj.bam_records_in_cells", reads.len());
        report.report_n(
            "vdj.umis",
            reads.iter().map(|read| read.umi.as_str()).collect::<HashSet<_>>().len(),
        );
        report.start_timer("cell_work.sterile_evidence");
        // Sterile/germline transcription remains read-based and is evaluated on
        // the original BAM evidence before UMI assembly. This keeps locus evidence
        // identical while rearrangement inference below becomes molecule-based.
        let mut sterile_acc: HashMap<Chain, SterileAccumulator> = Chain::ALL
            .into_iter()
            .filter_map(|c| {
                SterileAccumulator::new(self.reference, c, self.config.sterile_bins).map(|a| (c, a))
            })
            .collect();
        for read in &reads {
            if self.mapper.map(&read.sequence).is_none() {
                if let Some(chr) = &read.chr {
                    for &(s, e) in &read.ref_blocks {
                        for (chain, acc) in sterile_acc.iter_mut() {
                            if let Some((lchr, ls, le)) = self.reference.locus_bounds(*chain) {
                                if chr == lchr && e > ls && s < le {
                                    acc.observe_n(chr, s, e, &read.umi, 1);
                                }
                            }
                        }
                    }
                }
            }
        }

        report.stop_timer("cell_work.sterile_evidence");

        // Rearrangement inference is performed on a few observed contigs per UMI,
        // not every PCR/read observation. Disconnected/conflicting molecules with
        // the same UMI remain separate contigs rather than being force-merged.
        report.start_timer("cell_work.umi_contig_assembly");
        let umi_contigs = assemble_umi_contigs(&reads, self.config.min_umi_contig_overlap);
        report.report_n("vdj.umi_contigs", umi_contigs.len());
        report.stop_timer("cell_work.umi_contig_assembly");

        report.start_timer("cell_work.seed_and_local_alignment");
        let seeded: Vec<SeededRead<'_>> = umi_contigs
            .into_par_iter()
            .map(|contig| {
                let multiplicity = contig.members.len();
                let reverse_sequence = crate::sequence::reverse_complement(&contig.sequence);
                let ranked = self
                    .mapper
                    .segment_seed_scores_ranked_oriented(&contig.sequence, &reverse_sequence);
                let seed_scores: HashMap<usize, u32> = ranked
                    .iter()
                    .map(|&(idx, hits, _)| (idx, hits))
                    .collect();
                let (
                    best_by_kind,
                    alignments,
                    candidate_alignments,
                    family_candidates,
                    identity_seed_used,
                ) = self.best_read_matches(&contig.sequence, &reverse_sequence, &ranked);
                SeededRead {
                    umi: contig.umi,
                    sequence: contig.sequence,
                    members: contig.members,
                    multiplicity,
                    reverse_sequence,
                    seed_scores,
                    best_by_kind,
                    alignments,
                    candidate_alignments,
                    family_candidates,
                    identity_seed_used,
                }
            })
            .collect();
        report.report_n("vdj.seeded_contigs", seeded.len());
        report.report_n(
            "vdj.initial_local_alignments",
            seeded.iter().map(|evidence| evidence.candidate_alignments).sum(),
        );
        report.report_n(
            "vdj.identity_seed_nominated_contigs",
            seeded.iter().filter(|evidence| evidence.identity_seed_used).count(),
        );
        report.report_n(
            "vdj.coverage_seed_only_contigs",
            seeded.iter().filter(|evidence| !evidence.identity_seed_used).count(),
        );
        report.stop_timer("cell_work.seed_and_local_alignment");

        report.start_timer("cell_work.chain_calls_geometry");
        let mut chain_reads: HashMap<Chain, Vec<&SeededRead<'_>>> = HashMap::new();
        for evidence in &seeded {
            let mut seen_chain = HashSet::new();
            for kind in [SegmentKind::V, SegmentKind::J, SegmentKind::C] {
                if let Some(chains) = evidence.family_candidates.get(&kind) {
                    seen_chain.extend(chains.iter().copied());
                }
            }
            for (&kind, best) in &evidence.best_by_kind {
                let chain = self.reference.segments[best.segment_index].chain;
                let already_supported = evidence
                    .family_candidates
                    .get(&kind)
                    .is_some_and(|chains| chains.contains(&chain));
                if already_supported
                    || self.valid_terminal_extension_rescue(kind, best, evidence.sequence.len())
                {
                    seen_chain.insert(chain);
                }
            }
            // Family coverage is segment-kind specific. A strong IGK-V signal
            // therefore cannot manufacture IGK-C support, and abundant V seeds
            // cannot suppress sparse but coherent J evidence.
            for chain in seen_chain {
                chain_reads.entry(chain).or_default().push(evidence);
            }
        }

        let mut rearrangements: Vec<_> = Chain::ALL
            .par_iter()
            .filter_map(|&chain| {
                chain_reads
                    .get(&chain)
                    .and_then(|rs| self.call_chain(chain, rs))
            })
            .collect();
        rearrangements.sort_by_key(|call| call.chain);
        report.report_n("vdj.rearrangements_called", rearrangements.len());
        report.stop_timer("cell_work.chain_calls_geometry");

        report.start_timer("cell_work.recombination_activity");
        let sterile: Vec<_> = Chain::ALL
            .into_iter()
            .filter_map(|c| sterile_acc.remove(&c).map(|a| a.finish()))
            .collect();
        let recombination_activity = score_recombination_activity(gex, &cell);
        report.stop_timer("cell_work.recombination_activity");
        CellVdjSummary {
            cell,
            rearrangements,
            sterile,
            recombination_activity,
        }
    }

    fn best_read_matches(
        &self,
        sequence: &[u8],
        reverse_sequence: &[u8],
        ranked: &[(usize, u32, u32)],
    ) -> (
        HashMap<SegmentKind, BestReadMatch>,
        HashMap<usize, OrientedLocalAlignment>,
        usize,
        HashMap<SegmentKind, Vec<Chain>>,
        bool,
    ) {
        let coverage = self
            .mapper
            .coverage_family_scores_ranked_oriented(sequence, reverse_sequence);
        let terminal_v = self.mapper.terminal_seed_scores_ranked_oriented(
            sequence,
            reverse_sequence,
            SegmentKind::V,
        );
        let terminal_j = self.mapper.terminal_seed_scores_ranked_oriented(
            sequence,
            reverse_sequence,
            SegmentKind::J,
        );

        let mut family_candidates = HashMap::new();
        for kind in [SegmentKind::V, SegmentKind::J, SegmentKind::C] {
            let mut families = self.family_candidates_for_kind(&coverage, kind);
            if families.is_empty() && kind != SegmentKind::C {
                // Seed-poor fallback remains kind-specific for V/J only. C must
                // satisfy its own distributed coverage criterion; otherwise a
                // couple of accidental V-like kmers can inflate constant support.
                for &(segment_index, seed_hits, _) in ranked {
                    if seed_hits < self.config.min_seed_hits {
                        continue;
                    }
                    let segment = &self.reference.segments[segment_index];
                    if segment.kind != kind || families.contains(&segment.chain) {
                        continue;
                    }
                    families.push(segment.chain);
                    if families.len() >= self.config.family_candidates_per_kind {
                        break;
                    }
                }
            }
            if !families.is_empty() {
                family_candidates.insert(kind, families);
            }
        }

        // Terminal seeds are allowed to broaden candidate SEARCH, but never
        // family/chain EVIDENCE by themselves. This is the critical distinction:
        // conserved V-3'/J-5' kmers are cheap rescue nominators, not proof that
        // an arbitrary contig belongs to that receptor locus.
        let mut search_families = family_candidates.clone();
        for (kind, terminal) in [
            (SegmentKind::V, terminal_v.as_slice()),
            (SegmentKind::J, terminal_j.as_slice()),
        ] {
            let families = search_families.entry(kind).or_default();
            for &(segment_index, hits, _) in terminal {
                if hits < self.config.terminal_min_hits {
                    continue;
                }
                let chain = self.reference.segments[segment_index].chain;
                if !families.contains(&chain) {
                    families.push(chain);
                }
                if families.len() >= self.config.family_candidates_per_kind {
                    break;
                }
            }
        }

        let general_hits: HashMap<usize, u32> = ranked
            .iter()
            .map(|&(segment_index, hits, _)| (segment_index, hits))
            .collect();
        let mut candidate_list: Vec<(usize, u32, Option<(usize, usize)>)> = Vec::new();
        let mut candidate_positions: HashMap<usize, usize> = HashMap::new();
        let mut identity_seed_used = false;

        for kind in [SegmentKind::V, SegmentKind::J] {
            let Some(chains) = search_families.get(&kind) else {
                continue;
            };
            for &chain in chains {
                // Dedicated junction-facing rescue first. It is an additive
                // nominator only: whole-segment identity and 9-mer fallbacks
                // remain in the union and Smith-Waterman still decides.
                let terminal: &[(usize, u32, u32)] = match kind {
                    SegmentKind::V => terminal_v.as_slice(),
                    SegmentKind::J => terminal_j.as_slice(),
                    _ => &[],
                };
                let mut terminal_added = 0usize;
                for &(segment_index, terminal_hits, _) in terminal {
                    if terminal_hits < self.config.terminal_min_hits {
                        continue;
                    }
                    let segment = &self.reference.segments[segment_index];
                    if segment.chain != chain || segment.kind != kind {
                        continue;
                    }
                    let hits = general_hits
                        .get(&segment_index)
                        .copied()
                        .unwrap_or(terminal_hits);
                    if !candidate_positions.contains_key(&segment_index) {
                        candidate_positions.insert(segment_index, candidate_list.len());
                        candidate_list.push((segment_index, hits, None));
                    }
                    terminal_added += 1;
                    if terminal_added >= self.config.terminal_candidates_per_kind {
                        break;
                    }
                }

                // Then nominate by long, discriminative within-family seeds.
                let identity = self.mapper.identity_seed_scores_ranked_oriented(
                    sequence,
                    reverse_sequence,
                    chain,
                    kind,
                );
                let mut identity_added = 0usize;
                for score in identity {
                    if score.distinct_kmers < self.config.identity_min_hits {
                        continue;
                    }
                    let hits = general_hits.get(&score.segment_index).copied().unwrap_or(0);
                    if hits < self.config.min_seed_hits {
                        continue;
                    }
                    identity_seed_used = true;
                    let window = self.identity_alignment_window(kind, score);
                    if let Some(&position) = candidate_positions.get(&score.segment_index) {
                        if candidate_list[position].2.is_none() {
                            candidate_list[position].2 = window;
                        }
                    } else {
                        candidate_positions.insert(score.segment_index, candidate_list.len());
                        candidate_list.push((score.segment_index, hits, window));
                    }
                    identity_added += 1;
                    if identity_added >= self.config.identity_candidates_per_family_kind {
                        break;
                    }
                }

                // Always union a small sensitive 9-mer fallback from the same
                // chain/kind family. The long-seed index is an accelerator and
                // identity discriminator, never an exclusion gate.
                let mut coverage_added = 0usize;
                for &(segment_index, seed_hits, _) in ranked {
                    if seed_hits < self.config.min_seed_hits {
                        continue;
                    }
                    let segment = &self.reference.segments[segment_index];
                    if segment.chain != chain || segment.kind != kind {
                        continue;
                    }
                    if !candidate_positions.contains_key(&segment_index) {
                        candidate_positions.insert(segment_index, candidate_list.len());
                        candidate_list.push((segment_index, seed_hits, None));
                    }
                    coverage_added += 1;
                    if coverage_added >= self.config.coverage_candidates_per_family_kind {
                        break;
                    }
                }
            }
        }

        let candidate_alignment_count = candidate_list.len();
        let aligned: Vec<_> = candidate_list
            .into_iter()
            .map(|(segment_index, seed_hits, window)| {
                let segment = &self.reference.segments[segment_index];
                let alignment = if let Some((start, end)) = window {
                    let mut alignment = best_oriented_local_alignment_with_reverse(
                        sequence,
                        reverse_sequence,
                        &segment.sequence[start..end],
                    );
                    alignment.alignment.reference_start += start;
                    alignment.alignment.reference_end += start;
                    alignment
                } else {
                    best_oriented_local_alignment_with_reverse(
                        sequence,
                        reverse_sequence,
                        &segment.sequence,
                    )
                };
                (segment_index, seed_hits, alignment)
            })
            .collect();

        let mut best = HashMap::new();
        let mut alignments = HashMap::with_capacity(aligned.len());
        for (segment_index, seed_hits, alignment) in aligned {
            alignments.insert(segment_index, alignment);
            let local = alignment.alignment.score;
            if local < self.config.min_local_score {
                continue;
            }
            let segment = &self.reference.segments[segment_index];
            let candidate = BestReadMatch {
                segment_index,
                seed_hits,
                local_score: local,
                alignment,
            };
            let replace = best.get(&segment.kind).map_or(true, |old: &BestReadMatch| {
                candidate.local_score > old.local_score
                    || (candidate.local_score == old.local_score
                        && (candidate.seed_hits > old.seed_hits
                            || (candidate.seed_hits == old.seed_hits
                                && candidate.segment_index < old.segment_index)))
            });
            if replace {
                best.insert(segment.kind, candidate);
            }
        }
        (
            best,
            alignments,
            candidate_alignment_count,
            family_candidates,
            identity_seed_used,
        )
    }

    /// A terminal-only V/J candidate may rescue chain routing only when
    /// local alignment proves the expected recombination-facing germline end
    /// AND the observed molecule continues past that boundary. The alignment is
    /// already in its winning orientation, so this single check also establishes
    /// the direction of the extension without trusting BAM strand flags.
    fn valid_terminal_extension_rescue(
        &self,
        kind: SegmentKind,
        best: &BestReadMatch,
        query_len: usize,
    ) -> bool {
        let segment = &self.reference.segments[best.segment_index];
        match kind {
            SegmentKind::V => {
                self.v_reaches_junction(best, segment.sequence.len())
                    && best.alignment.alignment.query_end < query_len
            }
            SegmentKind::J => {
                self.j_starts_at_junction(best)
                    && best.alignment.alignment.query_start > 0
            }
            SegmentKind::D | SegmentKind::C => false,
        }
    }

    fn minimum_coverage_bins(&self, kind: SegmentKind) -> u32 {
        match kind {
            SegmentKind::V => self.config.min_v_coverage_bins,
            SegmentKind::J => self.config.min_j_coverage_bins,
            SegmentKind::C => self.config.min_c_coverage_bins,
            SegmentKind::D => 1,
        }
    }

    fn family_candidates_for_kind(
        &self,
        scores: &[CoverageFamilyScore],
        kind: SegmentKind,
    ) -> Vec<Chain> {
        let eligible: Vec<_> = scores
            .iter()
            .copied()
            .filter(|score| {
                score.kind == kind
                    && score.distinct_kmers >= self.config.min_seed_hits
                    && score.covered_bins() >= self.minimum_coverage_bins(kind)
            })
            .collect();
        let best_weight = eligible.first().map(|score| score.weighted_score).unwrap_or(0);
        eligible
            .into_iter()
            .filter(|score| {
                best_weight == 0
                    || (score.weighted_score as u64) * 2 >= best_weight as u64
            })
            .take(self.config.family_candidates_per_kind)
            .map(|score| score.chain)
            .collect()
    }

    fn identity_alignment_window(
        &self,
        kind: SegmentKind,
        score: IdentitySeedScore,
    ) -> Option<(usize, usize)> {
        let segment = &self.reference.segments[score.segment_index];
        let bins = self.mapper.coverage_bins() as usize;
        if bins == 0 || segment.sequence.len() < 32 {
            return None;
        }
        let set_bins: Vec<_> = (0..bins)
            .filter(|&bin| score.bin_mask & (1u8 << bin) != 0)
            .collect();
        let (&first, &last) = (set_bins.first()?, set_bins.last()?);
        let bin_start = |bin: usize| segment.sequence.len() * bin / bins;
        let bin_end = |bin: usize| segment.sequence.len() * (bin + 1) / bins;
        match kind {
            // When identity evidence reaches the recombination-facing V end,
            // align only from one bin before the first observed identity bin to
            // the true 3' end. Reference coordinates are restored afterwards.
            SegmentKind::V if last + 1 == bins => {
                Some((bin_start(first.saturating_sub(1)), segment.sequence.len()))
            }
            // Analogously keep the true J 5' end and extend one bin beyond the
            // last observed identity bin.
            SegmentKind::J if first == 0 => {
                Some((0, bin_end((last + 1).min(bins - 1)).min(segment.sequence.len())))
            }
            _ => None,
        }
    }

    fn call_chain(&self, chain: Chain, reads: &[&SeededRead<'_>]) -> Option<RearrangementCall> {
        let (v, j) = rayon::join(
            || self.best_segment(chain, SegmentKind::V, reads, None),
            || self.best_segment(chain, SegmentKind::J, reads, None),
        );
        let (Some(v), Some(j)) = (v, j) else {
            return None;
        };

        // A rearrangement requires molecule-level coherence: the same UMI must
        // support the selected V and J.  This deliberately replaces the old
        // cell/chain-wide UMI count that produced support ~= every UMI in a cell.
        let candidate_vj = umi_intersection(&v.umis, &j.umis);
        let coherent_vj = self.coherent_vj_umis(
            reads,
            &candidate_vj,
            v.support.segment_index,
            j.support.segment_index,
        );
        if coherent_vj.is_empty() {
            return None;
        }

        // D is intrinsically difficult: it is short, frequently clipped and may
        // already carry somatic mutations. V/J coherence establishes the heavy
        // rearrangement. Only then compare every D from the same locus inside the
        // bounded V..J junction, molecule by molecule. D therefore refines the HC
        // identity but never gates the HC call itself.
        let d_hypothesis = if chain.has_d() {
            self.best_d_hypothesis_in_vj_junction(
                chain,
                reads,
                &coherent_vj,
                v.support.segment_index,
                j.support.segment_index,
            )
        } else {
            None
        };
        let d = d_hypothesis.as_ref().map(|x| x.selected.clone());

        // C is supporting context, not evidence that creates a VJ call.  It must
        // be physically plausible and share at least one coherent VJ UMI.
        let c = self
            .best_segment(chain, SegmentKind::C, reads, Some(&coherent_vj))
            .filter(|candidate| {
                candidate.support.distance_to_recombination_center
                    <= self.config.max_constant_distance_bp
                    && candidate.umis.iter().any(|umi| coherent_vj.contains(umi))
            });

        let stage = if chain.has_d() {
            RecombinationStage::Vdj
        } else {
            RecombinationStage::Vj
        };
        let notation = format_call(
            chain,
            Some(&v.support),
            d.as_ref().map(|x| &x.support),
            Some(&j.support),
            c.as_ref().map(|x| &x.support),
        );
        let mut junction = self.measure_call_junction(
            reads,
            &coherent_vj,
            v.support.segment_index,
            j.support.segment_index,
            d.as_ref().map(|x| x.support.segment_index),
        );
        let supporting_reads = self.collect_supporting_reads(
            reads,
            &coherent_vj,
            v.support.segment_index,
            j.support.segment_index,
            d.as_ref().map(|x| x.support.segment_index),
            c.as_ref().map(|x| x.support.segment_index),
        );
        if let Some(measurement) = junction.as_mut() {
            refine_junction_observed_sequences(measurement, &supporting_reads);
        }
        Some(RearrangementCall {
            chain,
            stage,
            v: Some(v.support),
            d: d.map(|x| x.support),
            d_inferred_from_vj_junction: d_hypothesis.is_some(),
            d_hypothesis_margin: d_hypothesis.as_ref().map(|x| x.margin),
            j: Some(j.support),
            c: c.map(|x| x.support),
            total_supporting_umis: coherent_vj.len(),
            supporting_reads,
            junction,
            notation,
        })
    }

    fn measure_call_junction(
        &self,
        reads: &[&SeededRead<'_>],
        coherent_umis: &HashSet<String>,
        v_index: usize,
        j_index: usize,
        d_index: Option<usize>,
    ) -> Option<JunctionMeasurement> {
        let v_ref = &self.reference.segments[v_index].sequence;
        let j_ref = &self.reference.segments[j_index].sequence;
        let d_ref = d_index.map(|idx| self.reference.segments[idx].sequence.as_slice());

        let mut best: Option<(i32, JunctionMeasurement)> = None;
        for evidence in reads {
            if !coherent_umis.contains(evidence.umi) {
                continue;
            }
            let v_alignment = best_oriented_local_alignment_with_reverse(
                &evidence.sequence,
                &evidence.reverse_sequence,
                v_ref,
            );
            let j_alignment = best_oriented_local_alignment_with_reverse(
                &evidence.sequence,
                &evidence.reverse_sequence,
                j_ref,
            );
            if v_alignment.alignment.score < self.config.min_local_score
                || j_alignment.alignment.score < self.config.min_local_score
            {
                continue;
            }

            let d_alignment = d_ref.and_then(|reference| {
                bounded_d_alignment(
                    &evidence.sequence,
                    &evidence.reverse_sequence,
                    v_alignment,
                    j_alignment,
                    reference,
                )
            });
            // If an HC has a selected D hypothesis, junction geometry must be
            // measured with that D inside the V/J interval. There is intentionally
            // no absolute D score cutoff: short/mutated D sequence is expected.
            if d_ref.is_some() && d_alignment.is_none() {
                continue;
            }

            let Some(measurement) = measure_junction(JunctionInput {
                observed: &evidence.sequence,
                v: v_ref,
                v_alignment,
                d: d_ref,
                d_alignment,
                j: j_ref,
                j_alignment,
            }) else {
                continue;
            };
            let score = v_alignment.alignment.score
                + j_alignment.alignment.score
                + d_alignment.map_or(0, |x| x.alignment.score);
            if best
                .as_ref()
                .map_or(true, |(best_score, _)| score > *best_score)
            {
                best = Some((score, measurement));
            }
        }
        best.map(|(_, measurement)| measurement)
    }

    fn best_d_hypothesis_in_vj_junction(
        &self,
        chain: Chain,
        reads: &[&SeededRead<'_>],
        coherent_umis: &HashSet<String>,
        v_index: usize,
        j_index: usize,
    ) -> Option<DHypothesis> {
        let v_ref = &self.reference.segments[v_index].sequence;
        let j_ref = &self.reference.segments[j_index].sequence;
        let d_indices: Vec<usize> = self
            .reference
            .segments_for(chain, SegmentKind::D)
            .map(|(idx, _)| idx)
            .collect();
        if d_indices.is_empty() {
            return None;
        }

        // Each UMI contributes at most its best bounded score for each D gene.
        // This keeps PCR/read multiplicity from turning into false D confidence.
        let mut by_d: HashMap<usize, HashMap<&str, i32>> = HashMap::new();
        let mut reads_by_d: HashMap<usize, usize> = HashMap::new();
        for evidence in reads {
            if !coherent_umis.contains(evidence.umi) {
                continue;
            }
            let v_alignment = best_oriented_local_alignment_with_reverse(
                &evidence.sequence,
                &evidence.reverse_sequence,
                v_ref,
            );
            let j_alignment = best_oriented_local_alignment_with_reverse(
                &evidence.sequence,
                &evidence.reverse_sequence,
                j_ref,
            );
            if v_alignment.alignment.score < self.config.min_local_score
                || j_alignment.alignment.score < self.config.min_local_score
                || v_alignment.reverse_complement != j_alignment.reverse_complement
                || v_alignment.alignment.query_end > j_alignment.alignment.query_start
            {
                continue;
            }

            for &idx in &d_indices {
                let d_ref = &self.reference.segments[idx].sequence;
                let Some(hit) = bounded_d_alignment(
                    &evidence.sequence,
                    &evidence.reverse_sequence,
                    v_alignment,
                    j_alignment,
                    d_ref,
                ) else {
                    continue;
                };
                if hit.alignment.score <= 0 || hit.alignment.reference_end <= hit.alignment.reference_start {
                    continue;
                }
                by_d
                    .entry(idx)
                    .or_default()
                    .entry(evidence.umi)
                    .and_modify(|old| *old = (*old).max(hit.alignment.score))
                    .or_insert(hit.alignment.score);
                *reads_by_d.entry(idx).or_default() += evidence.multiplicity;
            }
        }

        let mut ranked: Vec<(usize, i32, usize)> = by_d
            .into_iter()
            .map(|(idx, by_umi)| {
                let umi_count = by_umi.len();
                let total = by_umi
                    .into_values()
                    .fold(0i32, |acc, score| acc.saturating_add(score));
                (idx, total, umi_count)
            })
            .collect();
        ranked.sort_by(|a, b| {
            b.1.cmp(&a.1)
                .then_with(|| b.2.cmp(&a.2))
                .then_with(|| a.0.cmp(&b.0))
        });
        let &(best_idx, best_score, best_umis) = ranked.first()?;
        let second_score = ranked.get(1).map(|x| x.1).unwrap_or(0);
        let seg = &self.reference.segments[best_idx];
        Some(DHypothesis {
            selected: SelectedSegment {
                support: GermlineSegmentSupport {
                    segment_index: best_idx,
                    id: seg.name.clone(),
                    kind: SegmentKind::D,
                    local_alignment_score: best_score,
                    supporting_umis: best_umis,
                    supporting_reads: reads_by_d.get(&best_idx).copied().unwrap_or(0),
                    locus_fraction: seg.locus_fraction,
                    distance_to_recombination_center: seg.distance_to_recombination_center,
                },
                umis: coherent_umis.clone(),
            },
            margin: best_score.saturating_sub(second_score),
        })
    }

    fn coherent_vj_umis(
        &self,
        reads: &[&SeededRead<'_>],
        candidate_umis: &HashSet<String>,
        v_index: usize,
        j_index: usize,
    ) -> HashSet<String> {
        let mut by_umi: HashMap<&str, Vec<&SeededRead<'_>>> = HashMap::new();
        for evidence in reads {
            if candidate_umis.contains(evidence.umi) {
                by_umi
                    .entry(evidence.umi)
                    .or_default()
                    .push(*evidence);
            }
        }
        by_umi
            .into_iter()
            .filter_map(|(umi, evidence)| {
                self.umi_has_observed_vj_path(&evidence, v_index, j_index)
                    .then(|| umi.to_string())
            })
            .collect()
    }

    fn umi_has_observed_vj_path(
        &self,
        reads: &[&SeededRead<'_>],
        v_index: usize,
        j_index: usize,
    ) -> bool {
        let v_ref = &self.reference.segments[v_index].sequence;
        let j_ref = &self.reference.segments[j_index].sequence;
        let mut v_reads: Vec<(usize, Vec<u8>)> = Vec::new();
        let mut j_reads: Vec<(usize, Vec<u8>)> = Vec::new();

        for (read_idx, evidence) in reads.iter().enumerate() {
            let v_match = evidence
                .best_by_kind
                .get(&SegmentKind::V)
                .filter(|x| x.segment_index == v_index);
            let j_match = evidence
                .best_by_kind
                .get(&SegmentKind::J)
                .filter(|x| x.segment_index == j_index);

            if let (Some(v), Some(j)) = (v_match, j_match) {
                if self.valid_vj_geometry(v, j, v_ref.len()) {
                    return true;
                }
            }
            if let Some(v) = v_match {
                if self.v_reaches_junction(v, v_ref.len()) {
                    v_reads.push((read_idx, oriented_read(&evidence.sequence, v.alignment.reverse_complement)));
                }
            }
            if let Some(j) = j_match {
                if self.j_starts_at_junction(j) {
                    j_reads.push((read_idx, oriented_read(&evidence.sequence, j.alignment.reverse_complement)));
                }
            }
        }

        for (v_idx, v_read) in &v_reads {
            for (j_idx, j_read) in &j_reads {
                if v_idx == j_idx {
                    continue;
                }
                let Some(merged) = overlap_reads(v_read, j_read, self.config.min_vj_read_overlap)
                else {
                    continue;
                };
                let v = best_oriented_local_alignment(&merged, v_ref);
                let j = best_oriented_local_alignment(&merged, j_ref);
                if !v.reverse_complement
                    && !j.reverse_complement
                    && self.valid_oriented_vj_geometry(v, j, v_ref.len())
                {
                    return true;
                }
            }
        }
        false
    }

    fn valid_vj_geometry(&self, v: &BestReadMatch, j: &BestReadMatch, v_len: usize) -> bool {
        if v.alignment.reverse_complement != j.alignment.reverse_complement {
            return false;
        }
        self.valid_oriented_vj_geometry(v.alignment, j.alignment, v_len)
    }

    fn valid_oriented_vj_geometry(
        &self,
        v: OrientedLocalAlignment,
        j: OrientedLocalAlignment,
        v_len: usize,
    ) -> bool {
        if !self.v_alignment_reaches_junction(v, v_len) || !self.j_alignment_starts_at_junction(j) {
            return false;
        }
        let va = v.alignment;
        let ja = j.alignment;
        if va.query_start > ja.query_start {
            return false;
        }
        va.query_end.saturating_sub(ja.query_start) <= self.config.max_vj_alignment_overlap
    }

    fn v_reaches_junction(&self, v: &BestReadMatch, v_len: usize) -> bool {
        self.v_alignment_reaches_junction(v.alignment, v_len)
    }

    fn v_alignment_reaches_junction(&self, v: OrientedLocalAlignment, v_len: usize) -> bool {
        v_len.saturating_sub(v.alignment.reference_end) <= self.config.max_v_end_distance
    }

    fn j_starts_at_junction(&self, j: &BestReadMatch) -> bool {
        self.j_alignment_starts_at_junction(j.alignment)
    }

    fn j_alignment_starts_at_junction(&self, j: OrientedLocalAlignment) -> bool {
        j.alignment.reference_start <= self.config.max_j_start_distance
    }

    fn collect_supporting_reads(
        &self,
        reads: &[&SeededRead<'_>],
        coherent_umis: &HashSet<String>,
        v_index: usize,
        j_index: usize,
        d_index: Option<usize>,
        c_index: Option<usize>,
    ) -> Vec<RearrangementSupportingRead> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for evidence in reads {
            if !coherent_umis.contains(evidence.umi) {
                continue;
            }
            for raw in &evidence.members {
                let key = (raw.read_name.clone(), raw.sequence.clone());
                if !seen.insert(key) {
                    continue;
                }

                // UMI contig coordinates do not necessarily correspond to any one
                // contributing BAM read. Re-align only the final-call audit reads
                // here so FASTA annotations remain coordinates on the raw sequence.
                let reverse = crate::sequence::reverse_complement(&raw.sequence);
                let v_alignment = audit_alignment(
                    &raw.sequence,
                    &reverse,
                    &self.reference.segments[v_index].sequence,
                    self.config.min_local_score,
                );
                let j_alignment = audit_alignment(
                    &raw.sequence,
                    &reverse,
                    &self.reference.segments[j_index].sequence,
                    self.config.min_local_score,
                );
                let d_alignment = d_index.and_then(|idx| {
                    audit_alignment(
                        &raw.sequence,
                        &reverse,
                        &self.reference.segments[idx].sequence,
                        self.config.min_local_score,
                    )
                });
                let c_alignment = c_index.and_then(|idx| {
                    audit_alignment(
                        &raw.sequence,
                        &reverse,
                        &self.reference.segments[idx].sequence,
                        self.config.min_local_score,
                    )
                });
                let supports_v = v_alignment.is_some();
                let supports_j = j_alignment.is_some();
                let supports_d = d_alignment.is_some();
                let supports_c = c_alignment.is_some();
                if !(supports_v || supports_j || supports_d || supports_c) {
                    continue;
                }
                out.push(RearrangementSupportingRead {
                    umi: raw.umi.clone(),
                    read_name: raw.read_name.clone(),
                    sequence: raw.sequence.clone(),
                    qualities: raw.qualities.clone(),
                    bam_is_reverse: raw.is_reverse,
                    is_supplementary: raw.is_supplementary,
                    supports_v,
                    supports_j,
                    supports_d,
                    supports_c,
                    v_alignment,
                    j_alignment,
                    d_alignment,
                    c_alignment,
                });
            }
        }
        out.sort_by(|a, b| {
            a.umi
                .cmp(&b.umi)
                .then_with(|| a.read_name.cmp(&b.read_name))
                .then_with(|| a.sequence.cmp(&b.sequence))
        });
        out
    }

    fn c_bam_locus_compatible(&self, evidence: &SeededRead<'_>, segment_index: usize) -> bool {
        const CONFIDENT_MAPQ: u8 = 20;
        let segment = &self.reference.segments[segment_index];
        let mut has_confident_genomic_mapping = false;

        for member in &evidence.members {
            if member.mapq < CONFIDENT_MAPQ || member.is_supplementary {
                continue;
            }
            let Some(chr) = member.chr.as_deref() else {
                continue;
            };
            has_confident_genomic_mapping = true;
            if chr == segment.chr
                && member
                    .ref_blocks
                    .iter()
                    .any(|&(start, end)| end > segment.start && start < segment.end)
            {
                return true;
            }
        }

        // Preserve sequence rescue for unmapped, low-MAPQ, or supplementary
        // evidence. A confidently mapped read elsewhere in the genome is not
        // allowed to masquerade as a constant-region read merely because it
        // resembles an IG/TR C segment locally.
        !has_confident_genomic_mapping
    }

    fn best_segment(
        &self,
        chain: Chain,
        kind: SegmentKind,
        reads: &[&SeededRead<'_>],
        restrict_umis: Option<&HashSet<String>>,
    ) -> Option<SelectedSegment> {
        // Rank germline candidates by molecule-level evidence. Multiple PCR reads
        // and multiple overlapping contigs from one UMI contribute at most their
        // strongest seed score for that germline segment.
        let mut seed_by_umi: HashMap<usize, HashMap<&str, u32>> = HashMap::new();
        let mut identity_by_umi: HashMap<usize, HashMap<&str, u32>> = HashMap::new();
        for evidence in reads {
            if restrict_umis.is_some_and(|allowed| !allowed.contains(evidence.umi)) {
                continue;
            }
            if kind == SegmentKind::C
                && evidence
                    .family_candidates
                    .get(&SegmentKind::C)
                    .is_some_and(|chains| chains.contains(&chain))
            {
                for score in self.mapper.identity_seed_scores_ranked_oriented(
                    &evidence.sequence,
                    &evidence.reverse_sequence,
                    chain,
                    SegmentKind::C,
                ) {
                    if score.distinct_kmers < self.config.identity_min_hits
                        || !self.c_bam_locus_compatible(evidence, score.segment_index)
                    {
                        continue;
                    }
                    identity_by_umi
                        .entry(score.segment_index)
                        .or_default()
                        .entry(evidence.umi)
                        .and_modify(|old| *old = (*old).max(score.weighted_score))
                        .or_insert(score.weighted_score);
                }
            }
            if matches!(kind, SegmentKind::D | SegmentKind::C) {
                if kind == SegmentKind::C
                    && !evidence
                        .family_candidates
                        .get(&SegmentKind::C)
                        .is_some_and(|chains| chains.contains(&chain))
                {
                    continue;
                }
                for (&idx, &score) in &evidence.seed_scores {
                    let s = &self.reference.segments[idx];
                    if s.chain == chain
                        && s.kind == kind
                        && (kind != SegmentKind::C || self.c_bam_locus_compatible(evidence, idx))
                    {
                        let by_umi = seed_by_umi.entry(idx).or_default();
                        by_umi
                            .entry(evidence.umi)
                            .and_modify(|old| *old = (*old).max(score))
                            .or_insert(score);
                    }
                }
            } else if let Some(best) = evidence.best_by_kind.get(&kind) {
                let s = &self.reference.segments[best.segment_index];
                if s.chain == chain {
                    let by_umi = seed_by_umi.entry(best.segment_index).or_default();
                    by_umi
                        .entry(evidence.umi)
                        .and_modify(|old| *old = (*old).max(best.seed_hits))
                        .or_insert(best.seed_hits);
                }
            }
        }
        let seed_totals: HashMap<usize, u32> = seed_by_umi
            .into_iter()
            .map(|(idx, scores)| {
                let total = scores
                    .into_values()
                    .fold(0u32, |acc, score| acc.saturating_add(score));
                (idx, total)
            })
            .collect();
        let identity_totals: HashMap<usize, u32> = identity_by_umi
            .into_iter()
            .map(|(idx, scores)| {
                let total = scores
                    .into_values()
                    .fold(0u32, |acc, score| acc.saturating_add(score));
                (idx, total)
            })
            .collect();
        let mut cand: Vec<_> = if kind == SegmentKind::C && !identity_totals.is_empty() {
            identity_totals.into_iter().collect()
        } else {
            seed_totals.into_iter().collect()
        };
        cand.sort_by_key(|x| std::cmp::Reverse(x.1));
        cand.truncate(if kind == SegmentKind::C {
            2
        } else {
            self.config.candidate_segments_per_kind
        });
        // D segments are short and can seed poorly; fall back to the complete,
        // usually small D set, but only on coherent V/J UMIs.
        if cand.is_empty() && kind == SegmentKind::D {
            cand = self
                .reference
                .segments_for(chain, kind)
                .map(|(i, _)| (i, 0))
                .collect();
        }

        cand.into_par_iter()
            .filter_map(|(idx, _)| {
                let seg = &self.reference.segments[idx];
                if kind == SegmentKind::C
                    && seg.distance_to_recombination_center > self.config.max_constant_distance_bp
                {
                    return None;
                }

                let mut best_score_by_umi: HashMap<&str, i32> = HashMap::new();
                let mut umis = HashSet::new();
                let mut supporting_reads = 0usize;
                for evidence in reads {
                    if restrict_umis.is_some_and(|allowed| !allowed.contains(evidence.umi)) {
                        continue;
                    }

                    if kind == SegmentKind::C {
                        if !self.c_bam_locus_compatible(evidence, idx) {
                            continue;
                        }
                        let family_supported = evidence
                            .family_candidates
                            .get(&SegmentKind::C)
                            .is_some_and(|chains| chains.contains(&chain));
                        if !family_supported
                            || evidence.seed_scores.get(&idx).copied().unwrap_or(0)
                                < self.config.min_seed_hits
                        {
                            continue;
                        }
                    }

                    // V/J evidence is assigned to the best local-alignment match
                    // for this contig across all receptor chains. Shared k-mers therefore
                    // nominate candidates but cannot make the same read support every
                    // homologous IG/TR segment. D remains a deliberate exception because
                    // its short germline sequence can be seed-poor after junctional editing.
                    let score = if matches!(kind, SegmentKind::D | SegmentKind::C) {
                        evidence
                            .alignments
                            .get(&idx)
                            .copied()
                            .unwrap_or_else(|| {
                                best_oriented_local_alignment_with_reverse(
                                    &evidence.sequence,
                                    &evidence.reverse_sequence,
                                    &seg.sequence,
                                )
                            })
                            .alignment
                            .score
                    } else {
                        let Some(best) = evidence.best_by_kind.get(&kind) else {
                            continue;
                        };
                        if best.segment_index != idx {
                            continue;
                        }
                        best.local_score
                    };
                    if score >= self.config.min_local_score {
                        best_score_by_umi
                            .entry(evidence.umi)
                            .and_modify(|old| *old = (*old).max(score))
                            .or_insert(score);
                        umis.insert(evidence.umi.to_string());
                        supporting_reads += evidence.multiplicity;
                    }
                }
                if supporting_reads == 0 {
                    return None;
                }
                let total = best_score_by_umi
                    .into_values()
                    .fold(0i32, |acc, score| acc.saturating_add(score));
                Some(SelectedSegment {
                    support: GermlineSegmentSupport {
                        segment_index: idx,
                        id: seg.name.clone(),
                        kind,
                        local_alignment_score: total,
                        supporting_umis: umis.len(),
                        supporting_reads,
                        locus_fraction: seg.locus_fraction,
                        distance_to_recombination_center: seg.distance_to_recombination_center,
                    },
                    umis,
                })
            })
            .max_by(|a, b| {
                a.support
                    .local_alignment_score
                    .cmp(&b.support.local_alignment_score)
                    .then_with(|| a.support.supporting_umis.cmp(&b.support.supporting_umis))
                    // Reverse the index comparison so the lower stable reference
                    // index wins exact ties independent of Rayon scheduling.
                    .then_with(|| b.support.segment_index.cmp(&a.support.segment_index))
            })
    }
}

#[derive(Debug, Clone)]
struct DHypothesis {
    selected: SelectedSegment,
    margin: i32,
}

/// Align a D reference only inside an already-coherent V/J junction and project
/// the hit back onto the oriented full molecule. This is the critical difference
/// from a whole-read D search: N/P sequence and flanking V/J cannot attract the D
/// alignment outside the biologically admissible interval.
fn bounded_d_alignment(
    forward: &[u8],
    reverse: &[u8],
    v: OrientedLocalAlignment,
    j: OrientedLocalAlignment,
    d_ref: &[u8],
) -> Option<OrientedLocalAlignment> {
    if v.reverse_complement != j.reverse_complement
        || v.alignment.query_end > j.alignment.query_start
    {
        return None;
    }
    let oriented = if v.reverse_complement { reverse } else { forward };
    let start = v.alignment.query_end;
    let end = j.alignment.query_start;
    if start >= end || end > oriented.len() || d_ref.is_empty() {
        return None;
    }
    let mut hit = local_alignment(&oriented[start..end], d_ref);
    if hit.score <= 0 {
        return None;
    }
    hit.query_start += start;
    hit.query_end += start;
    Some(OrientedLocalAlignment {
        alignment: hit,
        reverse_complement: v.reverse_complement,
    })
}

const CONSENSUS_MIN_BASE_QUAL: u8 = 20;
const CONSENSUS_MIN_OVERLAP: usize = 12;

#[derive(Debug, Clone, Copy)]
struct OverlapProjection {
    reverse_complement: bool,
    target_start: usize,
    read_start: usize,
    overlap: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct UmiBaseCall {
    base: u8,
    qual: u8,
    conflicted: bool,
}

fn refine_junction_observed_sequences(
    measurement: &mut JunctionMeasurement,
    reads: &[RearrangementSupportingRead],
) {
    repair_unresolved_bases(&mut measurement.observed_sequence, reads, |_| true);
    repair_unresolved_bases(&mut measurement.observed_v, reads, |read| read.supports_v);
    if !measurement.observed_d.is_empty() {
        repair_unresolved_bases(&mut measurement.observed_d, reads, |read| read.supports_d);
    }
    repair_unresolved_bases(&mut measurement.observed_j, reads, |read| read.supports_j);
}

/// Resolve only bases that are genuinely unresolved in the selected observed
/// molecule. Germline sequence is deliberately absent from this function: a base
/// can leave `N` only through independent sequencing evidence.
///
/// Reads are first collapsed within each UMI at every target position. One UMI
/// therefore contributes at most one vote, regardless of PCR/read multiplicity.
/// The highest-quality observation wins within a UMI; equal-quality conflicts
/// invalidate that UMI at that position. Across UMIs, a unique winning nucleotide
/// is accepted when its quality-weighted support is at least twice the runner-up.
fn repair_unresolved_bases<F>(
    target: &mut [u8],
    reads: &[RearrangementSupportingRead],
    keep_read: F,
) where
    F: Fn(&RearrangementSupportingRead) -> bool,
{
    let unresolved: Vec<usize> = target
        .iter()
        .enumerate()
        .filter_map(|(idx, base)| (!is_acgt(*base)).then_some(idx))
        .collect();
    if unresolved.is_empty() || target.is_empty() {
        return;
    }

    let mut by_position: HashMap<usize, HashMap<&str, UmiBaseCall>> = HashMap::new();
    for read in reads.iter().filter(|read| keep_read(read)) {
        if read.sequence.len() != read.qualities.len() || read.sequence.is_empty() {
            continue;
        }
        let min_overlap = CONSENSUS_MIN_OVERLAP.min(target.len()).min(read.sequence.len());
        if min_overlap == 0 {
            continue;
        }
        let Some(projection) = best_overlap_projection(target, &read.sequence, min_overlap) else {
            continue;
        };
        let oriented_sequence;
        let oriented_qualities;
        let (sequence, qualities): (&[u8], &[u8]) = if projection.reverse_complement {
            oriented_sequence = crate::sequence::reverse_complement(&read.sequence);
            oriented_qualities = read.qualities.iter().rev().copied().collect::<Vec<_>>();
            (&oriented_sequence, &oriented_qualities)
        } else {
            (&read.sequence, &read.qualities)
        };
        let target_end = projection.target_start + projection.overlap;
        for &target_pos in unresolved
            .iter()
            .filter(|&&pos| pos >= projection.target_start && pos < target_end)
        {
            let read_pos = projection.read_start + (target_pos - projection.target_start);
            let base = sequence[read_pos].to_ascii_uppercase();
            let qual = qualities[read_pos];
            if !is_acgt(base) || qual < CONSENSUS_MIN_BASE_QUAL {
                continue;
            }
            let call = by_position
                .entry(target_pos)
                .or_default()
                .entry(read.umi.as_str())
                .or_default();
            if qual > call.qual {
                call.base = base;
                call.qual = qual;
                call.conflicted = false;
            } else if qual == call.qual && call.base != 0 && call.base != base {
                call.conflicted = true;
            }
        }
    }

    for pos in unresolved {
        let Some(umis) = by_position.get(&pos) else {
            continue;
        };
        let mut support = [0u32; 4];
        for call in umis.values() {
            if call.conflicted || call.qual < CONSENSUS_MIN_BASE_QUAL {
                continue;
            }
            if let Some(idx) = base_index(call.base) {
                // Quality-weighting rewards a well-supported high-quality call,
                // while the UMI collapse above prevents PCR duplication from
                // manufacturing support.
                support[idx] = support[idx].saturating_add(u32::from(call.qual));
            }
        }
        let mut ranked: Vec<(u32, usize)> = support
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, score)| (score, idx))
            .collect();
        ranked.sort_unstable_by(|a, b| b.cmp(a));
        let (best, best_idx) = ranked[0];
        let second = ranked[1].0;
        if best >= u32::from(CONSENSUS_MIN_BASE_QUAL)
            && (second == 0 || best >= second.saturating_mul(2))
        {
            target[pos] = b"ACGT"[best_idx];
        }
    }
}

fn best_overlap_projection(
    target: &[u8],
    read: &[u8],
    min_overlap: usize,
) -> Option<OverlapProjection> {
    let reverse = crate::sequence::reverse_complement(read);
    let forward = best_oriented_overlap_projection(target, read, min_overlap, false);
    let reverse = best_oriented_overlap_projection(target, &reverse, min_overlap, true);
    match (forward, reverse) {
        (Some(a), Some(b)) => Some(if overlap_projection_better(a, b) { a } else { b }),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn best_oriented_overlap_projection(
    target: &[u8],
    read: &[u8],
    min_overlap: usize,
    reverse_complement: bool,
) -> Option<OverlapProjection> {
    if target.is_empty() || read.is_empty() {
        return None;
    }
    let min_overlap = min_overlap.min(target.len()).min(read.len());
    if min_overlap == 0 {
        return None;
    }
    let min_offset = -(read.len() as isize) + min_overlap as isize;
    let max_offset = target.len() as isize - min_overlap as isize;
    let mut best: Option<(usize, usize, OverlapProjection)> = None;
    for offset in min_offset..=max_offset {
        let target_start = offset.max(0) as usize;
        let read_start = (-offset).max(0) as usize;
        let overlap = (target.len() - target_start).min(read.len() - read_start);
        if overlap < min_overlap {
            continue;
        }
        let mut informative = 0usize;
        let mut mismatches = 0usize;
        for i in 0..overlap {
            let a = target[target_start + i].to_ascii_uppercase();
            let b = read[read_start + i].to_ascii_uppercase();
            if !is_acgt(a) || !is_acgt(b) {
                continue;
            }
            informative += 1;
            if a != b {
                mismatches += 1;
            }
        }
        if informative < min_overlap.saturating_sub(2) {
            continue;
        }
        if mismatches > (informative / 20).max(1) {
            continue;
        }
        let projection = OverlapProjection {
            reverse_complement,
            target_start,
            read_start,
            overlap,
        };
        let replace = best.as_ref().map_or(true, |(old_overlap, old_mismatches, old)| {
            overlap > *old_overlap
                || (overlap == *old_overlap
                    && (mismatches < *old_mismatches
                        || (mismatches == *old_mismatches
                            && overlap_projection_better(projection, *old))))
        });
        if replace {
            best = Some((overlap, mismatches, projection));
        }
    }
    best.map(|(_, _, projection)| projection)
}

fn overlap_projection_better(a: OverlapProjection, b: OverlapProjection) -> bool {
    a.overlap > b.overlap
        || (a.overlap == b.overlap
            && (!a.reverse_complement && b.reverse_complement
                || (a.reverse_complement == b.reverse_complement
                    && (a.target_start, a.read_start) < (b.target_start, b.read_start))))
}

fn is_acgt(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn assemble_umi_contigs<'a>(
    reads: &'a [BamReadEvidence],
    min_overlap: usize,
) -> Vec<UmiContig<'a>> {
    let mut by_umi: HashMap<&'a str, Vec<&'a BamReadEvidence>> = HashMap::new();
    for read in reads {
        by_umi.entry(read.umi.as_str()).or_default().push(read);
    }

    let mut umis: Vec<_> = by_umi.into_iter().collect();
    umis.sort_by(|a, b| a.0.cmp(b.0));

    let mut out = Vec::new();
    for (umi, mut members) in umis {
        members.sort_by(|a, b| {
            a.sequence
                .cmp(&b.sequence)
                .then_with(|| a.read_name.cmp(&b.read_name))
        });

        let mut contigs: Vec<UmiContig<'a>> = Vec::new();
        for read in members {
            let mut best: Option<(usize, SequenceMerge)> = None;
            for (idx, contig) in contigs.iter().enumerate() {
                let Some(candidate) = merge_sequences_any_orientation(
                    &contig.sequence,
                    &read.sequence,
                    min_overlap,
                ) else {
                    continue;
                };
                let replace = best.as_ref().map_or(true, |(_, old)| {
                    candidate.better_than(old)
                });
                if replace {
                    best = Some((idx, candidate));
                }
            }

            if let Some((idx, merged)) = best {
                contigs[idx].sequence = merged.sequence;
                contigs[idx].members.push(read);
            } else {
                contigs.push(UmiContig {
                    umi,
                    sequence: read.sequence.clone(),
                    members: vec![read],
                });
            }
        }

        // A read added late can bridge two contigs that were previously disconnected.
        // Repeatedly merge the strongest compatible pair, but never force a UMI into
        // one molecule when no observed sequence overlap exists.
        loop {
            let mut best: Option<(usize, usize, SequenceMerge)> = None;
            for left in 0..contigs.len() {
                for right in left + 1..contigs.len() {
                    let Some(candidate) = merge_sequences_any_orientation(
                        &contigs[left].sequence,
                        &contigs[right].sequence,
                        min_overlap,
                    ) else {
                        continue;
                    };
                    let replace = best.as_ref().map_or(true, |(_, _, old)| {
                        candidate.better_than(old)
                    });
                    if replace {
                        best = Some((left, right, candidate));
                    }
                }
            }
            let Some((left, right, merged)) = best else {
                break;
            };
            contigs[left].sequence = merged.sequence;
            let right_contig = contigs.remove(right);
            contigs[left].members.extend(right_contig.members);
        }

        contigs.sort_by(|a, b| {
            b.members
                .len()
                .cmp(&a.members.len())
                .then_with(|| b.sequence.len().cmp(&a.sequence.len()))
                .then_with(|| a.sequence.cmp(&b.sequence))
        });
        out.extend(contigs);
    }
    out
}

#[derive(Debug, Clone)]
struct SequenceMerge {
    sequence: Vec<u8>,
    overlap: usize,
    mismatches: usize,
    extension: usize,
}

impl SequenceMerge {
    fn better_than(&self, other: &Self) -> bool {
        self.overlap > other.overlap
            || (self.overlap == other.overlap
                && (self.mismatches < other.mismatches
                    || (self.mismatches == other.mismatches
                        && (self.extension > other.extension
                            || (self.extension == other.extension
                                && self.sequence < other.sequence)))))
    }
}

fn merge_sequences_any_orientation(
    left: &[u8],
    right: &[u8],
    min_overlap: usize,
) -> Option<SequenceMerge> {
    let reverse = crate::sequence::reverse_complement(right);
    let forward = merge_oriented_sequences(left, right, min_overlap);
    let reverse = merge_oriented_sequences(left, &reverse, min_overlap);
    match (forward, reverse) {
        (Some(a), Some(b)) => Some(if a.better_than(&b) { a } else { b }),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn merge_oriented_sequences(
    left: &[u8],
    right: &[u8],
    min_overlap: usize,
) -> Option<SequenceMerge> {
    if left.is_empty() || right.is_empty() {
        return None;
    }
    let min_overlap = min_overlap.min(left.len()).min(right.len());
    if min_overlap == 0 {
        return None;
    }

    let min_offset = -(right.len() as isize) + min_overlap as isize;
    let max_offset = left.len() as isize - min_overlap as isize;
    let mut best: Option<SequenceMerge> = None;

    for offset in min_offset..=max_offset {
        let left_start = offset.max(0) as usize;
        let right_start = (-offset).max(0) as usize;
        let overlap = (left.len() - left_start).min(right.len() - right_start);
        if overlap < min_overlap {
            continue;
        }

        let mut mismatches = 0usize;
        for i in 0..overlap {
            let a = left[left_start + i].to_ascii_uppercase();
            let b = right[right_start + i].to_ascii_uppercase();
            if a != b && a != b'N' && b != b'N' {
                mismatches += 1;
            }
        }
        let max_mismatches = (overlap / 20).max(1);
        if mismatches > max_mismatches {
            continue;
        }

        let start = offset.min(0);
        let end = (left.len() as isize).max(offset + right.len() as isize);
        let mut sequence = Vec::with_capacity((end - start) as usize);
        for pos in start..end {
            let a = if pos >= 0 && pos < left.len() as isize {
                Some(left[pos as usize].to_ascii_uppercase())
            } else {
                None
            };
            let right_pos = pos - offset;
            let b = if right_pos >= 0 && right_pos < right.len() as isize {
                Some(right[right_pos as usize].to_ascii_uppercase())
            } else {
                None
            };
            sequence.push(match (a, b) {
                (Some(x), None) | (None, Some(x)) => x,
                (Some(x), Some(y)) if x == y => x,
                (Some(b'N'), Some(y)) => y,
                (Some(x), Some(b'N')) => x,
                (Some(_), Some(_)) => b'N',
                (None, None) => unreachable!(),
            });
        }

        let candidate = SequenceMerge {
            extension: sequence.len().saturating_sub(left.len().max(right.len())),
            sequence,
            overlap,
            mismatches,
        };
        if best.as_ref().map_or(true, |old| candidate.better_than(old)) {
            best = Some(candidate);
        }
    }
    best
}

fn audit_alignment(
    sequence: &[u8],
    reverse: &[u8],
    reference: &[u8],
    min_local_score: i32,
) -> Option<ReadSegmentAlignment> {
    let hit = best_oriented_local_alignment_with_reverse(sequence, reverse, reference);
    (hit.alignment.score >= min_local_score).then(|| read_segment_alignment(hit))
}

fn read_segment_alignment(hit: OrientedLocalAlignment) -> ReadSegmentAlignment {
    ReadSegmentAlignment {
        score: hit.alignment.score,
        query_start: hit.alignment.query_start,
        query_end: hit.alignment.query_end,
        reference_start: hit.alignment.reference_start,
        reference_end: hit.alignment.reference_end,
        reverse_complement: hit.reverse_complement,
    }
}

fn oriented_read(sequence: &[u8], reverse: bool) -> Vec<u8> {
    if reverse {
        crate::sequence::reverse_complement(sequence)
    } else {
        sequence.to_vec()
    }
}

fn overlap_reads(left: &[u8], right: &[u8], min_overlap: usize) -> Option<Vec<u8>> {
    let max_overlap = left.len().min(right.len());
    if max_overlap < min_overlap {
        return None;
    }
    for overlap in (min_overlap..=max_overlap).rev() {
        let a = &left[left.len() - overlap..];
        let b = &right[..overlap];
        let mismatches = a.iter().zip(b).filter(|(x, y)| x != y).count();
        if mismatches > (overlap / 20).max(1) {
            continue;
        }
        let mut merged = Vec::with_capacity(left.len() + right.len() - overlap);
        merged.extend_from_slice(&left[..left.len() - overlap]);
        for (&x, &y) in a.iter().zip(b) {
            merged.push(if x == y { x } else { b'N' });
        }
        merged.extend_from_slice(&right[overlap..]);
        return Some(merged);
    }
    None
}

fn umi_intersection(a: &HashSet<String>, b: &HashSet<String>) -> HashSet<String> {
    if a.len() <= b.len() {
        a.iter().filter(|umi| b.contains(*umi)).cloned().collect()
    } else {
        b.iter().filter(|umi| a.contains(*umi)).cloned().collect()
    }
}

fn support_id(support: Option<&GermlineSegmentSupport>) -> &str {
    support.map_or("?", |support| support.id.as_str())
}

fn format_call(
    chain: Chain,
    v: Option<&GermlineSegmentSupport>,
    d: Option<&GermlineSegmentSupport>,
    j: Option<&GermlineSegmentSupport>,
    c: Option<&GermlineSegmentSupport>,
) -> String {
    if chain.has_d() {
        format!(
            "{}:{}-({})-{}-{}",
            chain,
            support_id(v),
            support_id(d),
            support_id(j),
            support_id(c)
        )
    } else {
        format!(
            "{}:{}-{}-{}",
            chain,
            support_id(v),
            support_id(j),
            support_id(c)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapper::VdjMapperConfig;
    use crate::types::VdjSegment;

    fn segment(name: &str, chain: Chain, sequence: &[u8]) -> VdjSegment {
        VdjSegment {
            name: name.to_string(),
            transcript_id: format!("{name}_tx"),
            gene_id: name.to_string(),
            chain,
            kind: SegmentKind::V,
            chr: "chr1".to_string(),
            start: 0,
            end: sequence.len() as u32,
            strand_minus: false,
            locus_rank: 0,
            locus_fraction: 0.0,
            distance_to_recombination_center: 0,
            sequence: sequence.to_vec(),
        }
    }

    #[test]
    fn umi_intersection_is_not_union() {
        let a = HashSet::from(["u1".to_string(), "u2".to_string()]);
        let b = HashSet::from(["u2".to_string(), "u3".to_string()]);
        assert_eq!(umi_intersection(&a, &b), HashSet::from(["u2".to_string()]));
    }

    #[test]
    fn shared_seed_read_supports_only_best_cross_chain_segment() {
        let reference = VdjReference {
            segments: vec![
                segment("IGHV_best", Chain::Igh, b"AACCGGTTAAAA"),
                segment("TRAV_shared", Chain::Tra, b"AACCGGTTCCCC"),
            ],
        };
        let mapper = VdjMapper::new(
            reference.clone(),
            VdjMapperConfig {
                k: 5,
                min_v_hits: 1,
                min_j_hits: 1,
                min_c_hits: 1,
            },
        );
        let analyzer = PosteriorAnalyzer::new(
            &reference,
            &mapper,
            PosteriorConfig {
                min_seed_hits: 1,
                min_local_score: 4,
                ..PosteriorConfig::default()
            },
        );
        let read = BamReadEvidence {
            cell: "cell".to_string(),
            umi: "umi".to_string(),
            read_name: "read".to_string(),
            sequence: b"AACCGGTTAAAA".to_vec(),
            qualities: vec![30; 12],
            chr: None,
            ref_start: None,
            ref_end: None,
            ref_blocks: Vec::new(),
            mapq: 0,
            is_reverse: false,
            is_secondary: false,
            is_supplementary: false,
        };
        let reverse = crate::sequence::reverse_complement(&read.sequence);
        let ranked = mapper.segment_seed_scores_ranked_oriented(&read.sequence, &reverse);
        let seed_scores: HashMap<usize, u32> = ranked
            .iter()
            .map(|&(idx, hits, _)| (idx, hits))
            .collect();
        assert!(seed_scores.contains_key(&0));
        assert!(seed_scores.contains_key(&1));

        let (best, _, _, _, _) = analyzer.best_read_matches(&read.sequence, &reverse, &ranked);
        assert_eq!(best[&SegmentKind::V].segment_index, 0);
    }
    #[test]
    fn umi_assembly_extends_overlapping_reads_but_keeps_disconnected_molecules() {
        fn read(name: &str, sequence: &[u8]) -> BamReadEvidence {
            BamReadEvidence {
                cell: "cell".to_string(),
                umi: "umi".to_string(),
                read_name: name.to_string(),
                sequence: sequence.to_vec(),
                qualities: vec![30; sequence.len()],
                chr: None,
                ref_start: None,
                ref_end: None,
                ref_blocks: Vec::new(),
                mapq: 0,
                is_reverse: false,
                is_secondary: false,
                is_supplementary: false,
            }
        }

        let reads = vec![
            read("r1", b"AAAACCCCGGGGTTTTAAAA"),
            read("r2", b"GGGGTTTTAAAACCCCGGGG"),
            read("r3", b"TGCATGCATGCATGCATGCA"),
        ];
        let contigs = assemble_umi_contigs(&reads, 12);
        assert_eq!(contigs.len(), 2);
        assert!(contigs.iter().any(|x| x.members.len() == 2 && x.sequence.len() > 20));
        assert!(contigs.iter().any(|x| x.members.len() == 1));
    }

    #[test]
    fn umi_assembly_can_merge_reverse_complemented_observation() {
        let a = b"AAAACCCCGGGGTTTTAAAA";
        let b = crate::sequence::reverse_complement(b"GGGGTTTTAAAACCCCGGGG");
        let merged = merge_sequences_any_orientation(a, &b, 12).expect("compatible overlap");
        assert!(merged.sequence.len() > a.len());
        assert!(merged.overlap >= 12);
    }

    #[test]
    fn unresolved_base_is_repaired_from_quality_weighted_umi_evidence() {
        fn support(umi: &str, name: &str, sequence: &[u8], qualities: &[u8]) -> RearrangementSupportingRead {
            RearrangementSupportingRead {
                umi: umi.to_string(),
                read_name: name.to_string(),
                sequence: sequence.to_vec(),
                qualities: qualities.to_vec(),
                bam_is_reverse: false,
                is_supplementary: false,
                supports_v: true,
                supports_j: true,
                supports_d: false,
                supports_c: false,
                v_alignment: None,
                j_alignment: None,
                d_alignment: None,
                c_alignment: None,
            }
        }

        let mut target = b"AACCGGNTTAACCGG".to_vec();
        let q30 = vec![30; target.len()];
        let q35 = vec![35; target.len()];
        let q40 = vec![40; target.len()];
        let reads = vec![
            // Two PCR observations from one UMI still contribute one molecule vote.
            support("u1", "r1", b"AACCGGATTAACCGG", &q30),
            support("u1", "r2", b"AACCGGATTAACCGG", &q35),
            support("u2", "r3", b"AACCGGATTAACCGG", &q40),
        ];
        repair_unresolved_bases(&mut target, &reads, |_| true);
        assert_eq!(target, b"AACCGGATTAACCGG");
    }

    #[test]
    fn conflicting_umi_evidence_keeps_base_unresolved() {
        fn support(umi: &str, sequence: &[u8]) -> RearrangementSupportingRead {
            RearrangementSupportingRead {
                umi: umi.to_string(),
                read_name: umi.to_string(),
                sequence: sequence.to_vec(),
                qualities: vec![35; sequence.len()],
                bam_is_reverse: false,
                is_supplementary: false,
                supports_v: true,
                supports_j: true,
                supports_d: false,
                supports_c: false,
                v_alignment: None,
                j_alignment: None,
                d_alignment: None,
                c_alignment: None,
            }
        }

        let mut target = b"AACCGGNTTAACCGG".to_vec();
        let reads = vec![
            support("u1", b"AACCGGATTAACCGG"),
            support("u2", b"AACCGGGTTAACCGG"),
        ];
        repair_unresolved_bases(&mut target, &reads, |_| true);
        assert_eq!(target, b"AACCGGNTTAACCGG");
    }

    #[test]
    fn low_quality_base_does_not_repair_unresolved_position() {
        let sequence = b"AACCGGATTAACCGG";
        let mut qualities = vec![35; sequence.len()];
        qualities[6] = CONSENSUS_MIN_BASE_QUAL - 1;
        let read = RearrangementSupportingRead {
            umi: "u1".to_string(),
            read_name: "r1".to_string(),
            sequence: sequence.to_vec(),
            qualities,
            bam_is_reverse: false,
            is_supplementary: false,
            supports_v: true,
            supports_j: true,
            supports_d: false,
            supports_c: false,
            v_alignment: None,
            j_alignment: None,
            d_alignment: None,
            c_alignment: None,
        };
        let mut target = b"AACCGGNTTAACCGG".to_vec();
        repair_unresolved_bases(&mut target, &[read], |_| true);
        assert_eq!(target, b"AACCGGNTTAACCGG");
    }

    #[test]
    fn d_alignment_is_bounded_to_the_established_vj_junction() {
        let observed = b"AAAACCCCACGGGGTTTTAAAA";
        let reverse = crate::sequence::reverse_complement(observed);
        let v_ref = b"AAAACCCC";
        let d_ref = b"GGGG";
        let j_ref = b"TTTTAAAA";
        let v = best_oriented_local_alignment_with_reverse(observed, &reverse, v_ref);
        let j = best_oriented_local_alignment_with_reverse(observed, &reverse, j_ref);
        let d = bounded_d_alignment(observed, &reverse, v, j, d_ref).expect("bounded D hit");
        assert!(!d.reverse_complement);
        assert!(d.alignment.query_start >= v.alignment.query_end);
        assert!(d.alignment.query_end <= j.alignment.query_start);
        assert_eq!(d.alignment.reference_end - d.alignment.reference_start, 4);
        assert_eq!(d.alignment.score, 8);
    }

}
