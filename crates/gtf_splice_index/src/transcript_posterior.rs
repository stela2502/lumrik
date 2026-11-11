//! Bayesian interpretation of incomplete transcript observations.
//!
//! A poly(A)-biased long read is a window onto a transcript, not a complete
//! description of its structure. Within that window, even a few bases of an
//! exon establish its presence, and an absent exon between observed flanks is
//! evidence for skipping. Outside the window, absence means *unknown*.
//!
//! We therefore compare only the first through last observed exon units.
//! Molecules with wider informative windows can supply abundance priors for
//! shorter molecules, but a prior must never rescue a model contradicting
//! observed evidence. This is deliberately separate from transcript discovery:
//! posterior probabilities must not depend on UMI iteration order.
use crate::exon_catalog::{ExonMask, GeneExonCatalog};

/// One candidate model and its informative-molecule support.
#[derive(Debug, Clone)]
pub struct TranscriptPrior<'a> {
    pub id: usize,
    pub exon_mask: &'a ExonMask,
    /// Number of confidently assigned, independently deduplicated molecules.
    pub support: usize,
}

/// Normalized posterior probability for a candidate transcript.
#[derive(Debug, Clone, PartialEq)]
pub struct TranscriptPosterior {
    pub id: usize,
    pub probability: f64,
}

/// Score models given a partially observed molecule and informative-molecule priors.
///
/// Uses a symmetric Dirichlet pseudocount `alpha` for every candidate and a
/// strict likelihood (1 for agreement within the observed span, 0 otherwise).
/// This is a conservative Bayesian baseline, not a sequencing-error model.
/// Incompatible models receive zero probability even when very abundant.
/// The caller must provide priors derived from informative molecules, **not**
/// from the same ambiguous molecules being scored.
///
/// Returns an empty vector for no exon observations, invalid `alpha`, or when
/// no candidate explains the evidence. Results are sorted by posterior then ID.
pub fn posterior_for_partial_observation(
    catalog: &GeneExonCatalog,
    observed: &ExonMask,
    candidates: &[TranscriptPrior<'_>],
    alpha: f64,
) -> Vec<TranscriptPosterior> {
    if !alpha.is_finite() || alpha <= 0.0 { return Vec::new(); }
    let n = catalog.units.len();
    let bit = |mask: &ExonMask, i: usize| {
        let mut probe = ExonMask::default();
        probe.set(2 * i);
        mask.contains(&probe)
    };
    let Some(first) = (0..n).find(|&i| bit(observed, i)) else { return Vec::new(); };
    let last = (first..n).rfind(|&i| bit(observed, i)).unwrap();
    // Intronic evidence must not be mistaken for a mature transcript.
    for i in 0..n.saturating_sub(1) {
        let mut probe = ExonMask::default();
        probe.set(2 * i + 1);
        if observed.contains(&probe) { return Vec::new(); }
    }
    let mut weights: Vec<_> = candidates.iter().filter_map(|candidate| {
        let compatible = (first..=last).all(|i| bit(observed, i) == bit(candidate.exon_mask, i));
        compatible.then_some((candidate.id, candidate.support as f64 + alpha))
    }).collect();
    let total: f64 = weights.iter().map(|(_, weight)| weight).sum();
    if !total.is_finite() || total <= 0.0 { return Vec::new(); }
    let mut result: Vec<_> = weights.drain(..).map(|(id, weight)| TranscriptPosterior {
        id, probability: weight / total,
    }).collect();
    result.sort_by(|a, b| b.probability.total_cmp(&a.probability).then(a.id.cmp(&b.id)));
    result
}

/// Reference/alternative observations at one SNP on one deduplicated molecule.
/// Counts are read-level evidence, not independent molecule counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MoleculeSnpEvidence {
    pub snp_id: usize,
    pub ref_reads: u32,
    pub alt_reads: u32,
}

/// Allele evidence learned for a transcript from *other informative UMIs*.
/// A SNP not listed for a candidate is uninformative, not reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TranscriptSnpSupport {
    pub transcript_id: usize,
    pub snp_id: usize,
    pub ref_molecules: u64,
    pub alt_molecules: u64,
}

/// Combine structural compatibility, molecule-derived transcript priors, and
/// allele evidence without modifying the existing exon-only resolver.
///
/// The structural posterior is computed first: an observed exon contradiction
/// always excludes a transcript. SNP evidence only redistributes probability
/// among structurally compatible models. A SNP outside the molecule's observed
/// region, or one without learned association for *every* surviving candidate,
/// is ignored rather than silently treated as wild type.
///
/// The allele likelihood is a beta-binomial predictive distribution (up to a
/// candidate-independent binomial coefficient), with symmetric `beta` smoothing.
/// Training counts are deduplicated *molecules*, and the queried molecule MUST
/// be excluded from training to prevent self-reinforcement. Multiple reads at
/// the same SNP on one UMI are correlated; we deliberately reduce them to one
/// allele observation (ref, alt or mixed), not independent trials.
///
/// For mixed ref/alt reads, both outcomes are given half weight. This avoids
/// overconfidence from PCR/ONT read depth and does not claim allele phasing
/// when the UMI contains conflicting calls. SNPs are conditionally treated as
/// independent across loci; correlated SNPs require a future haplotype model.
///
/// Invalid smoothing, duplicate training entries, or invalid observations
/// fail closed with an empty result. No SNP evidence yields the exon-only
/// posterior exactly.
pub fn posterior_with_snp_evidence(
    catalog: &GeneExonCatalog,
    observed: &ExonMask,
    candidates: &[TranscriptPrior<'_>],
    alpha: f64,
    snps: &[MoleculeSnpEvidence],
    training: &[TranscriptSnpSupport],
    beta: f64,
) -> Vec<TranscriptPosterior> {
    if !beta.is_finite() || beta <= 0.0 { return Vec::new(); }
    let structural = posterior_for_partial_observation(catalog, observed, candidates, alpha);
    if structural.is_empty() || snps.is_empty() { return structural; }
    let mut seen = std::collections::HashSet::new();
    if training.iter().any(|t| !seen.insert((t.transcript_id, t.snp_id))) {
        return Vec::new();
    }
    let mut observed_ids = std::collections::HashSet::new();
    if snps.iter().any(|s| !observed_ids.insert(s.snp_id)) { return Vec::new(); }
    let mut log_weights: Vec<(usize, f64)> = structural.iter()
        .map(|p| (p.id, p.probability.ln())).collect();
    let mut informative_snp = false;
    for snp in snps {
        if snp.ref_reads == 0 && snp.alt_reads == 0 { continue; }
        // Missing training for any compatible model makes this locus
        // non-discriminating. Never mistake missing evidence for REF.
        let Some(supports) = log_weights.iter().map(|(id, _)|
            training.iter().find(|t| t.transcript_id == *id && t.snp_id == snp.snp_id)
        ).collect::<Option<Vec<_>>>() else { continue; };
        let mut first_likelihood: Option<f64> = None;
        for ((_, weight), t) in log_weights.iter_mut().zip(supports) {
            let denom = t.ref_molecules as f64 + t.alt_molecules as f64 + 2.0 * beta;
            if !denom.is_finite() || denom <= 0.0 { return Vec::new(); }
            let pref = (t.ref_molecules as f64 + beta) / denom;
            let palt = (t.alt_molecules as f64 + beta) / denom;
            let likelihood = match (snp.ref_reads > 0, snp.alt_reads > 0) {
                (true, false) => pref,
                (false, true) => palt,
                (true, true) => 0.5 * (pref + palt),
                (false, false) => 1.0,
            };
            if let Some(first) = first_likelihood {
                if likelihood != first { informative_snp = true; }
            } else {
                first_likelihood = Some(likelihood);
            }
            *weight += likelihood.ln();
        }
    }
    // Non-discriminating SNPs must leave the structural posterior unchanged.
    if !informative_snp { return structural; }
    let max_log = log_weights.iter().map(|(_, w)| *w)
        .fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = log_weights.iter().map(|(_, w)| (*w - max_log).exp()).sum();
    if !total.is_finite() || total <= 0.0 { return Vec::new(); }
    let mut result: Vec<_> = log_weights.into_iter().map(|(id, w)| TranscriptPosterior {
        id, probability: (w - max_log).exp() / total,
    }).collect();
    result.sort_by(|a, b| b.probability.total_cmp(&a.probability).then(a.id.cmp(&b.id)));
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::RefBlock;

    #[test]
    fn unobserved_alternative_first_and_second_exons_remain_ambiguous() {
        // Two mutually exclusive 5' structures: transcript 1 uses exon 1,
        // transcript 2 uses exon 2. No molecule reaches either exon.
        // Even 10,000 shared-region molecules cannot identify their ratio.
        let catalog = GeneExonCatalog { units: (0..100)
            .map(|i| RefBlock::new(100 * i, 100 * i + 30)).collect() };
        let observed = catalog.mask_for_blocks(&catalog.units[40..]);
        let first = catalog.mask_for_blocks(&[&catalog.units[0..1], &catalog.units[40..]].concat());
        let second = catalog.mask_for_blocks(&[&catalog.units[1..2], &catalog.units[40..]].concat());
        let priors = [
            TranscriptPrior { id: 1, exon_mask: &first, support: 0 },
            TranscriptPrior { id: 2, exon_mask: &second, support: 0 },
        ];
        let posterior = posterior_for_partial_observation(&catalog, &observed, &priors, 1.0);
        assert_eq!(posterior.len(), 2, "unobserved 5' exons must not define a deletion");
        assert_eq!(posterior[0].probability, 0.5);
        assert_eq!(posterior[1].probability, 0.5);
        // No SNP associations were learned; an allele observation cannot
        // invent transcript identity without training evidence.
        let with_snp = posterior_with_snp_evidence(&catalog, &observed, &priors,
            1.0, &[MoleculeSnpEvidence { snp_id: 42, ref_reads: 0, alt_reads: 1 }], &[], 1.0);
        assert_eq!(with_snp, posterior);
        // A 10-base overlap with the first exon is positive evidence.
        // Its alternative second exon is now an internal skip only if
        // the molecule also reaches the shared downstream region.
        let informative = catalog.mask_for_blocks(&[
            RefBlock::new(0, 10), catalog.units[40],
        ]);
        let identified = posterior_for_partial_observation(&catalog, &informative, &priors, 1.0);
        assert_eq!(identified.len(), 1);
        assert_eq!(identified[0].id, 1);
    }

    #[test]
    fn ten_base_flank_preserves_internal_skip_and_unknown_terminal_exons() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(100, 130), RefBlock::new(200, 230),
            RefBlock::new(300, 330), RefBlock::new(400, 430),
        ] };
        // Negative-strand gene: the observed window covers only 10 bases of
        // one flank, skips the middle exon, and reaches the other flank.
        let observed = catalog.mask_for_blocks(&[
            RefBlock::new(200, 210), RefBlock::new(400, 430),
        ]);
        let skip = catalog.mask_for_blocks(&[
            catalog.units[0], catalog.units[1], catalog.units[3],
        ]);
        let retained = catalog.mask_for_blocks(&catalog.units);
        let extra_skip = catalog.mask_for_blocks(&[catalog.units[1], catalog.units[3]]);
        let priors = [
            TranscriptPrior { id: 230, exon_mask: &skip, support: 9 },
            TranscriptPrior { id: 231, exon_mask: &extra_skip, support: 90 },
            TranscriptPrior { id: 100, exon_mask: &retained, support: 100 },
        ];
        let result = posterior_for_partial_observation(&catalog, &observed, &priors, 1.0);
        // The unobserved terminal exon cannot distinguish 230 from 231;
        // their informative-molecule priors determine their relative weight.
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].id, 231);
        assert!((result[0].probability - 91.0 / 101.0).abs() < 1e-12);
        assert_eq!(result[1].id, 230);
        // The retained-exon model is contradicted inside the observed window.
        assert!(!result.iter().any(|p| p.id == 100));
    }

    #[test]
    fn no_evidence_has_no_posterior() {
        let catalog = GeneExonCatalog { units: vec![RefBlock::new(10, 20)] };
        let mask = catalog.mask_for_blocks(&catalog.units);
        assert!(posterior_for_partial_observation(&catalog, &ExonMask::default(),
            &[TranscriptPrior { id: 1, exon_mask: &mask, support: 10 }], 1.0).is_empty());
    }

    /// Adversarial negative-strand TP53 analogue: the 5' end is progressively
    /// truncated, but the observed internal skip must remain real. A 10-bp
    /// exon flank is positive evidence, not an absent exon.
    #[test]
    fn nasty_partial_windows_and_very_strong_wrong_prior() {
        let catalog = GeneExonCatalog { units: (0..6)
            .map(|i| RefBlock::new(100 + 100*i, 140 + 100*i)).collect() };
        // 230 retains the unobserved 5' exons; 231 skips one of them.
        let model_230 = catalog.mask_for_blocks(&[
            catalog.units[0], catalog.units[1], catalog.units[2],
            catalog.units[4], catalog.units[5],
        ]);
        let model_231 = catalog.mask_for_blocks(&[
            catalog.units[1], catalog.units[2], catalog.units[4], catalog.units[5],
        ]);
        let full = catalog.mask_for_blocks(&catalog.units);
        let priors = [
            TranscriptPrior { id: 230, exon_mask: &model_230, support: 2 },
            TranscriptPrior { id: 231, exon_mask: &model_231, support: 999_999 },
            TranscriptPrior { id: 100, exon_mask: &full, support: 9_999_999 },
        ];
        // All 40 bases versus just 10 bases of the lower-coordinate flank.
        for flank in [40, 10] {
            let observed = catalog.mask_for_blocks(&[
                RefBlock::new(300, 300 + flank), RefBlock::new(500, 540),
            ]);
            let result = posterior_for_partial_observation(&catalog, &observed, &priors, 1.0);
            assert_eq!(result.len(), 2, "flank={flank}: {result:?}");
            assert_eq!(result[0].id, 231); // prior resolves the unknown terminal exon
            assert!(!result.iter().any(|p| p.id == 100),
                "a huge prior must not rescue a contradicted retained exon");
            assert!((result.iter().map(|p| p.probability).sum::<f64>() - 1.0).abs() < 1e-12);
        }
        // Once exon 0 is actually observed, model 231 is impossible even
        // with a 500,000-fold stronger prior.
        let observed = catalog.mask_for_blocks(&[
            RefBlock::new(100, 110), RefBlock::new(200, 210),
            RefBlock::new(300, 310), RefBlock::new(500, 540),
        ]);
        let result = posterior_for_partial_observation(&catalog, &observed, &priors, 1.0);
        assert_eq!(result, vec![TranscriptPosterior { id: 230, probability: 1.0 }]);
    }

    /// The caller must not turn a one-exon window into an exon-deletion call.
    /// Both distinct structures remain compatible, so their priors decide.
    #[test]
    fn nasty_single_exon_window_is_ambiguous_not_a_deletion() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 50), RefBlock::new(100, 140),
            RefBlock::new(200, 240), RefBlock::new(300, 340),
        ] };
        let retained = catalog.mask_for_blocks(&catalog.units);
        let skipped = catalog.mask_for_blocks(&[catalog.units[0], catalog.units[2], catalog.units[3]]);
        let priors = [
            TranscriptPrior { id: 1, exon_mask: &retained, support: 3 },
            TranscriptPrior { id: 2, exon_mask: &skipped, support: 7 },
        ];
        let observed = catalog.mask_for_blocks(&[RefBlock::new(300, 310)]);
        let result = posterior_for_partial_observation(&catalog, &observed, &priors, 1.0);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].id, 2);
        assert!((result[0].probability - 8.0/12.0).abs() < 1e-12);
        assert!((result[1].probability - 4.0/12.0).abs() < 1e-12);
    }

    /// An observed intronic segment cannot be explained away by a prior.
    #[test]
    fn nasty_intron_and_invalid_hyperparameters_fail_closed() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(30, 40), RefBlock::new(50, 60),
        ] };
        let retained = catalog.mask_for_blocks(&catalog.units);
        let priors = [TranscriptPrior { id: 1, exon_mask: &retained, support: 100 }];
        let intronic = catalog.mask_for_blocks(&[RefBlock::new(20, 21), RefBlock::new(30, 40)]);
        assert!(posterior_for_partial_observation(&catalog, &intronic, &priors, 1.0).is_empty());
        for alpha in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(posterior_for_partial_observation(&catalog, &retained, &priors, alpha).is_empty(),
                "invalid alpha={alpha:?} was accepted");
        }
    }

    /// Cross the 64-bit boundary: large exon masks must not silently truncate
    /// distant exon evidence and turn an incompatible model into a match.
    #[test]
    fn nasty_large_catalog_preserves_high_exon_bits() {
        let catalog = GeneExonCatalog { units: (0..80)
            .map(|i| RefBlock::new(100*i, 100*i + 40)).collect() };
        let mut skip_blocks = catalog.units.clone();
        skip_blocks.remove(65);
        let skipped = catalog.mask_for_blocks(&skip_blocks);
        let retained = catalog.mask_for_blocks(&catalog.units);
        let priors = [
            TranscriptPrior { id: 230, exon_mask: &skipped, support: 1 },
            TranscriptPrior { id: 100, exon_mask: &retained, support: 1_000_000 },
        ];
        let observed = catalog.mask_for_blocks(&[
            RefBlock::new(6400, 6410), RefBlock::new(6600, 6610),
        ]);
        let result = posterior_for_partial_observation(&catalog, &observed, &priors, 1.0);
        assert_eq!(result, vec![TranscriptPosterior { id: 230, probability: 1.0 }]);
    }

    /// Input order must not alter posterior ranking or probability; an
    /// ambiguous molecule must never bootstrap its own informative prior.
    #[test]
    fn nasty_prior_order_independence_and_dirichlet_smoothing() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(30, 40), RefBlock::new(50, 60),
        ] };
        let a = catalog.mask_for_blocks(&catalog.units);
        let b = catalog.mask_for_blocks(&[catalog.units[0], catalog.units[2]]);
        let observed = catalog.mask_for_blocks(&[RefBlock::new(50, 60)]);
        let priors = [
            TranscriptPrior { id: 7, exon_mask: &a, support: 0 },
            TranscriptPrior { id: 2, exon_mask: &b, support: 0 },
        ];
        let forward = posterior_for_partial_observation(&catalog, &observed, &priors, 0.5);
        let reverse = posterior_for_partial_observation(&catalog, &observed,
            &[priors[1].clone(), priors[0].clone()], 0.5);
        assert_eq!(forward, reverse);
        assert_eq!(forward, vec![
            TranscriptPosterior { id: 2, probability: 0.5 },
            TranscriptPosterior { id: 7, probability: 0.5 },
        ]);
    }
    #[test]
    fn snp_rescues_short_window_but_never_overrides_observed_exon() {
        let catalog = GeneExonCatalog { units: (0..6)
            .map(|i| RefBlock::new(100 + 100*i, 140 + 100*i)).collect() };
        let wt = catalog.mask_for_blocks(&catalog.units);
        let skip = catalog.mask_for_blocks(&[
            catalog.units[0], catalog.units[1], catalog.units[2],
            catalog.units[4], catalog.units[5],
        ]);
        let priors = [
            TranscriptPrior { id: 100, exon_mask: &wt, support: 90 },
            TranscriptPrior { id: 230, exon_mask: &skip, support: 10 },
        ];
        let training = [
            TranscriptSnpSupport { transcript_id: 100, snp_id: 7, ref_molecules: 90, alt_molecules: 0 },
            TranscriptSnpSupport { transcript_id: 230, snp_id: 7, ref_molecules: 0, alt_molecules: 10 },
        ];
        // Only the 3'-most exon is observed: structural priors favor WT.
        let short = catalog.mask_for_blocks(&[RefBlock::new(600, 610)]);
        let alt = [MoleculeSnpEvidence { snp_id: 7, ref_reads: 0, alt_reads: 50 }];
        let without = posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &[], &training, 1.0);
        assert_eq!(without[0].id, 100);
        let with = posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &alt, &training, 1.0);
        assert_eq!(with[0].id, 230, "allele-linked evidence should resolve a short read");
        assert!(with[0].probability > 0.5);
        // Even a strong ALT association cannot erase directly observed exon 3.
        let retained = catalog.mask_for_blocks(&[RefBlock::new(400, 410),
            RefBlock::new(500, 510), RefBlock::new(600, 610)]);
        let hard = posterior_with_snp_evidence(&catalog, &retained, &priors, 1.0, &alt, &training, 1.0);
        assert_eq!(hard, vec![TranscriptPosterior { id: 100, probability: 1.0 }]);
    }

    #[test]
    fn snp_missing_unlinked_mixed_and_invalid_inputs() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(100, 140), RefBlock::new(200, 240), RefBlock::new(300, 340),
        ] };
        let a = catalog.mask_for_blocks(&catalog.units);
        let b = catalog.mask_for_blocks(&[catalog.units[0], catalog.units[2]]);
        let priors = [
            TranscriptPrior { id: 1, exon_mask: &a, support: 9 },
            TranscriptPrior { id: 2, exon_mask: &b, support: 1 },
        ];
        let short = catalog.mask_for_blocks(&[RefBlock::new(300, 310)]);
        let base = posterior_for_partial_observation(&catalog, &short, &priors, 1.0);
        let training = [
            TranscriptSnpSupport { transcript_id: 1, snp_id: 7, ref_molecules: 90, alt_molecules: 0 },
            TranscriptSnpSupport { transcript_id: 2, snp_id: 7, ref_molecules: 0, alt_molecules: 10 },
        ];
        for snps in [
            vec![],
            vec![MoleculeSnpEvidence { snp_id: 7, ref_reads: 0, alt_reads: 0 }],
            vec![MoleculeSnpEvidence { snp_id: 999, ref_reads: 0, alt_reads: 2 }],
            vec![MoleculeSnpEvidence { snp_id: 7, ref_reads: 1, alt_reads: 1000 }],
        ] {
            let result = posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &snps, &training, 1.0);
            assert_eq!(result, base, "unavailable or mixed evidence must not change prior");
        }
        let alt = [MoleculeSnpEvidence { snp_id: 7, ref_reads: 0, alt_reads: 1 }];
        assert!(posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &alt, &training, 0.0).is_empty());
        assert!(posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &alt, &[training[0]], 1.0) == base);
        assert!(posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &[alt[0], alt[0]], &training, 1.0).is_empty());
        assert!(posterior_with_snp_evidence(&catalog, &short, &priors, 1.0, &alt, &[training[0], training[0]], 1.0).is_empty());
    }

}
