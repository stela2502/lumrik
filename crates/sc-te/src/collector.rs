use std::collections::HashSet;

use anyhow::{Context, Result};
use int_to_str::int_to_str::IntToStr;
use mapping_info::MappingInfo;
use read_tag_table::ReadTagRecord;
use rust_htslib::bam::HeaderView;
use sc_mapper::process::SamReadCluster;
use scdata::{GeneUmiHash, MatrixValueType, Scdata};

use crate::TeIndex;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct MultiMolecule {
    cell_id: u64,
    umi_id: u64,
    candidates: Vec<u64>,
}

pub struct TeCollector {
    /// Confident TE-overlapping mappings from the original 10x BAM. These are the anchor model.
    pub anchor: Scdata,
    /// Reads that become TE-unambiguous only after permissive remapping. Never used as anchors.
    pub rescued_unique: Scdata,
    anchor_features: HashSet<u64>,
    pending_multi: HashSet<MultiMolecule>,
    report: MappingInfo,
    threads: usize,
}

pub struct TeResult {
    pub anchor: Scdata,
    pub rescued_unique: Scdata,
    pub multi_em: Scdata,
    pub multi_anchored_em: Scdata,
    pub multi_unanchored_em: Scdata,
    pub report: MappingInfo,
}

impl TeResult {
    pub fn combined(mut self) -> Scdata {
        self.anchor.merge_values(self.rescued_unique);
        self.anchor.merge_values(self.multi_em);
        self.anchor
    }
}

impl TeCollector {
    pub fn new(threads: usize) -> Self {
        Self {
            anchor: Scdata::new(threads.max(1), MatrixValueType::Real),
            rescued_unique: Scdata::new(threads.max(1), MatrixValueType::Real),
            anchor_features: HashSet::new(),
            pending_multi: HashSet::new(),
            report: MappingInfo::new(None, 0.0, 0),
            threads: threads.max(1),
        }
    }

    pub fn add_anchor(&mut self, cell_id: u64, umi_id: u64, candidates: &[u64]) -> bool {
        self.report.report("te.original.overlaps_te");
        if candidates.len() != 1 {
            self.report.report("te.original.overlapping_te_features.multiple");
            return false;
        }
        let feature = candidates[0];
        let inserted = self.anchor.try_insert_value(
            &cell_id,
            GeneUmiHash(feature, umi_id),
            1.0,
            &mut self.report,
        );
        if inserted {
            self.anchor_features.insert(feature);
            self.report.report("te.anchor.molecules");
        } else {
            self.report.report("te.anchor.duplicate_molecule");
        }
        inserted
    }

    /// Consume a complete permissive-remap cluster. Candidate identity is based on
    /// distinct spatial TE features, not raw STAR alignment count.
    pub fn push_cluster(
        &mut self,
        cluster: SamReadCluster,
        header: &HeaderView,
        index: &mut TeIndex,
    ) -> Result<bool> {
        let tag = ReadTagRecord::from_qname(&cluster.read_id)
            .with_context(|| format!("failed to decode cell/UMI from {}", cluster.read_id))?;
        let cell_id = IntToStr::new(&tag.cell_seq).into_u64();
        let umi_id = IntToStr::new(&tag.umi_seq).into_u64();
        let mut candidates = HashSet::new();
        for mapped in cluster.records {
            candidates.extend(index.record_overlaps(&mapped.record, header)?);
        }
        if candidates.is_empty() {
            self.report.report("te.remap.no_te_candidate");
            return Ok(false);
        }
        let mut candidates: Vec<_> = candidates.into_iter().collect();
        candidates.sort_unstable();
        if candidates.len() == 1 {
            if self.rescued_unique.try_insert_value(
                &cell_id,
                GeneUmiHash(candidates[0], umi_id),
                1.0,
                &mut self.report,
            ) {
                self.report.report("te.remap.rescued_unique");
            }
        } else {
            self.pending_multi.insert(MultiMolecule {
                cell_id,
                umi_id,
                candidates,
            });
        }
        Ok(true)
    }

    /// Multimapper-only EM. Original anchors are used only to diagnose whether each
    /// multimapper candidate set intersects the observed anchor model; they do not
    /// influence the EM probabilities.
    pub fn finish(mut self, index: &TeIndex, max_iter: usize, epsilon: f64) -> TeResult {
        let molecules: Vec<_> = self.pending_multi.into_iter().collect();
        self.report.report_n("te.multi.molecules", molecules.len());
        let n_features = index.len();
        let mut abundance = vec![0.0f64; n_features];
        for molecule in &molecules {
            let n_anchor = molecule
                .candidates
                .iter()
                .filter(|id| self.anchor_features.contains(id))
                .count();
            if n_anchor == 0 {
                self.report.report("te.multi.unanchored");
            } else {
                self.report.report("te.multi.anchored");
            }
            match n_anchor {
                0 => {}
                1 => self.report.report("te.multi.anchor_candidates.1"),
                _ => self.report.report("te.multi.anchor_candidates.multiple"),
            }
            let w = 1.0 / molecule.candidates.len() as f64;
            for &id in &molecule.candidates {
                abundance[id as usize] += w;
            }
        }

        let mut iterations = 0usize;
        let mut converged = molecules.is_empty();
        for iter in 0..max_iter {
            iterations = iter + 1;
            let mut next = vec![0.0f64; n_features];
            for molecule in &molecules {
                let mass: f64 = molecule
                    .candidates
                    .iter()
                    .map(|&id| abundance[id as usize])
                    .sum();
                if mass > 0.0 {
                    for &id in &molecule.candidates {
                        next[id as usize] += abundance[id as usize] / mass;
                    }
                } else {
                    let w = 1.0 / molecule.candidates.len() as f64;
                    for &id in &molecule.candidates {
                        next[id as usize] += w;
                    }
                }
            }
            let delta = next
                .iter()
                .zip(&abundance)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            abundance = next;
            if delta <= epsilon {
                converged = true;
                break;
            }
        }
        self.report.report_n("te.em.iterations", iterations);
        if converged {
            self.report.report("te.em.converged");
        } else {
            self.report.report("te.em.max_iterations");
        }

        let mut multi_em = Scdata::new(self.threads, MatrixValueType::Real);
        let mut multi_anchored_em = Scdata::new(self.threads, MatrixValueType::Real);
        let mut multi_unanchored_em = Scdata::new(self.threads, MatrixValueType::Real);
        for molecule in molecules {
            let anchored = molecule
                .candidates
                .iter()
                .any(|id| self.anchor_features.contains(id));
            let mass: f64 = molecule
                .candidates
                .iter()
                .map(|&id| abundance[id as usize])
                .sum();
            let uniform = 1.0 / molecule.candidates.len() as f64;
            for &id in &molecule.candidates {
                let weight = if mass > 0.0 {
                    abundance[id as usize] / mass
                } else {
                    uniform
                } as f32;
                let gh = GeneUmiHash(id, molecule.umi_id);
                multi_em.try_insert_value(&molecule.cell_id, gh, weight, &mut self.report);
                if anchored {
                    multi_anchored_em.try_insert_value(
                        &molecule.cell_id,
                        gh,
                        weight,
                        &mut self.report,
                    );
                } else {
                    multi_unanchored_em.try_insert_value(
                        &molecule.cell_id,
                        gh,
                        weight,
                        &mut self.report,
                    );
                }
            }
        }
        TeResult {
            anchor: self.anchor,
            rescued_unique: self.rescued_unique,
            multi_em,
            multi_anchored_em,
            multi_unanchored_em,
            report: self.report,
        }
    }
}
