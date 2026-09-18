//! Receptor-family data model.
//!
//! This module owns biological collection, alignment validation, outlier
//! HC membership validation and LC clone splitting. Reassignment of rejected
//! cells is deliberately owned by the outer analysis/orchestration layer. Plotting consumes the frozen
//! result; it does not decide family membership or recompute mutation depths.

use crate::NovelVRegistries;

#[derive(Debug, Clone)]
pub struct Receptor {
    pub id: String,
    pub chain: String,
    pub v: String,
    pub d: String,
    pub j: String,
    pub c: String,
    pub cdr3_nt: String,
    pub naive: String,
    pub observed: String,
    /// Number of raw sc-vdj reconstructions collapsed into this selected receptor.
    pub alternative_reconstructions: usize,
    /// Read-level support for this receptor reconstruction. Used only to guard
    /// provisional V registration against single-read accidents.
    pub evidence_reads: usize,
}

#[derive(Debug, Clone)]
pub struct CellReceptor {
    pub cell_id: String,
    pub hc: Receptor,
    /// Productive light-chain receptor observations for this cell. Their clone
    /// membership is owned entirely by the finalized HC Family.
    pub lc: Vec<Receptor>,
}

#[derive(Debug, Clone, Copy)]
pub struct FamilyConfig {
    pub max_cdr3_distance: usize,
    /// Permissive first-pass alignment identity. Failure ejects the cell.
    pub min_alignment_identity: f64,
    /// Hard identity required when an ejected cell tries another family/clone.
    pub hard_alignment_identity: f64,
    /// Robust mutation outlier fence: median + max(min_extra, MAD * multiplier).
    pub mutation_mad_multiplier: usize,
    pub mutation_min_extra: usize,
}

impl Default for FamilyConfig {
    fn default() -> Self {
        Self {
            max_cdr3_distance: 3,
            min_alignment_identity: 0.50,
            hard_alignment_identity: 0.75,
            mutation_mad_multiplier: 4,
            mutation_min_extra: 6,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndelEvent {
    pub naive_pos: usize,
    pub inserted: bool,
    pub len: usize,
}

#[derive(Debug, Clone)]
pub struct MutationMeasurement {
    pub substitutions: usize,
    pub indels: Vec<IndelEvent>,
    pub informative_pairs: usize,
    pub matching_pairs: usize,
    pub identity: f64,
}

impl MutationMeasurement {
    pub fn mutation_events(&self) -> usize { self.substitutions + self.indels.len() }
}

#[derive(Debug, Clone)]
pub struct AlignedCell {
    pub cell: CellReceptor,
    pub hc_mutations: MutationMeasurement,
}

#[derive(Debug, Clone)]
pub struct LightMember {
    pub cell: String,
    pub receptor: Receptor,
    pub mutations: MutationMeasurement,
}

#[derive(Debug, Clone)]
pub struct LightClone {
    pub name: String,
    pub root: Receptor,
    pub members: Vec<LightMember>,
}

#[derive(Debug, Clone)]
pub struct Family {
    pub name: String,
    pub root: Receptor,
    candidates: Vec<CellReceptor>,
    pub members: Vec<AlignedCell>,
    pub light_clones: Vec<LightClone>,
    pub max_hc_mutations: Option<usize>,
}

impl Family {
    pub fn new(name: impl Into<String>, seed: CellReceptor) -> Self {
        let root = seed.hc.clone();
        Self { name: name.into(), root, candidates: vec![seed], members: Vec::new(), light_clones: Vec::new(), max_hc_mutations: None }
    }

    /// Initial collection is deliberately structural only. No NW decision is
    /// allowed here. Complete-link CDR3 compatibility is checked against every
    /// cell already collected in this provisional family.
    pub fn add_candidate(&mut self, cell: CellReceptor, cfg: &FamilyConfig) -> Result<(), CellReceptor> {
        if !same_hc_labels(&self.root, &cell.hc)
            || self.candidates.iter().any(|x| edit_distance(&x.hc.cdr3_nt, &cell.hc.cdr3_nt) > cfg.max_cdr3_distance)
        {
            return Err(cell);
        }
        self.candidates.push(cell);
        Ok(())
    }

    /// Align all provisional members to this family's model, report mutation
    /// depths, and eject alignment failures plus robust mutation-count outliers.
    /// The surviving member set is frozen until explicit reassignment.
    pub fn align(&mut self, cfg: &FamilyConfig, novel_v: &NovelVRegistries) -> Vec<CellReceptor> {
        let candidates = std::mem::take(&mut self.candidates);
        let mut measured = Vec::new();
        let mut ejected = Vec::new();
        for cell in candidates {
            match align_fragment(&novel_v.effective_naive(&cell.hc), &cell.hc.observed) {
                Some(m) if m.identity >= cfg.min_alignment_identity => measured.push((cell, m)),
                Some(_) | None => ejected.push(cell),
            }
        }
        let depths: Vec<usize> = measured.iter().map(|(_, m)| m.mutation_events()).collect();
        let cutoff = robust_upper_cutoff(&depths, cfg);
        self.members.clear();
        for (cell, m) in measured {
            if cutoff.is_some_and(|c| m.mutation_events() > c) {
                ejected.push(cell);
            } else {
                self.members.push(AlignedCell { cell, hc_mutations: m });
            }
        }
        self.max_hc_mutations = self.members.iter().map(|x| x.hc_mutations.mutation_events()).max();
        ejected
    }

    /// Second-pass admission for an ejected HC. Structural HC labels must fit;
    /// then the cell is aligned directly into this family's alignment model and
    /// must pass the hard alignment gate. CDR3 is not used as a second rescue
    /// veto: the alignment model is the final arbiter on reassignment.
    pub fn try_integrate(&mut self, cell: CellReceptor, cfg: &FamilyConfig, novel_v: &NovelVRegistries) -> Result<(), CellReceptor> {
        if !same_hc_labels(&self.root, &cell.hc) { return Err(cell); }
        let Some(m) = align_fragment(&novel_v.effective_naive(&cell.hc), &cell.hc.observed) else { return Err(cell) };
        if m.identity < cfg.hard_alignment_identity { return Err(cell); }
        if let Some(cutoff) = robust_upper_cutoff(&self.members.iter().map(|x| x.hc_mutations.mutation_events()).collect::<Vec<_>>(), cfg) {
            if m.mutation_events() > cutoff { return Err(cell); }
        }
        self.members.push(AlignedCell { cell, hc_mutations: m });
        self.max_hc_mutations = self.members.iter().map(|x| x.hc_mutations.mutation_events()).max();
        Ok(())
    }

    /// Split light chains only after HC membership is final.
    ///
    /// LC clone membership is evaluated against the candidate clone's naive
    /// root sequence, not against the receptor's own naive reconstruction.
    /// That makes the alignment informative for clone assignment: a receptor
    /// can fit one LC root well and another poorly even when both share V/J.
    /// After the first structural pass, mutation-count outliers are removed
    /// from their current clone and tried against every other compatible LC
    /// clone in this HC family. A receptor that fits none becomes a new root.
    pub fn split_light_chains(&mut self, cfg: &FamilyConfig, novel_v: &NovelVRegistries) {
        self.light_clones.clear();
        let mut pending: Vec<(String, Receptor)> = self.members.iter()
            .flat_map(|m| m.cell.lc.iter().cloned().map(|lc| (m.cell.cell_id.clone(), lc)))
            .collect();
        pending.sort_by(|a,b| (&a.1.chain,&a.1.v,&a.1.j,&a.1.cdr3_nt,&a.0).cmp(&(&b.1.chain,&b.1.v,&b.1.j,&b.1.cdr3_nt,&b.0)));

        // First collection is observation-driven. A bad/missing genomic V must
        // not prevent mutually compatible receptor observations from meeting.
        for (cell, receptor) in pending {
            if let Some((i, _membership)) = best_light_clone(&self.light_clones, &receptor, cfg, None, false, true, novel_v) {
                let mutations = measure_receptor(&receptor, novel_v);
                self.light_clones[i].members.push(LightMember { cell, receptor, mutations });
            } else {
                self.light_clones.push(new_light_clone(cell, receptor, novel_v));
            }
        }

        self.redistribute_lc_outliers(cfg, novel_v);
    }

    /// Remove mutation-depth outliers from their current LC clone and give
    /// every ejected receptor a genuine second chance against the other LC
    /// roots in the same finalized HC family. The target clone must satisfy
    /// the structural V/J/CDR3 gate, the hard identity gate, and its existing
    /// mutation-depth distribution. If no existing clone explains the receptor,
    /// it seeds a new LC clone rather than being discarded.
    fn redistribute_lc_outliers(&mut self, cfg: &FamilyConfig, novel_v: &NovelVRegistries) {
        let mut ejected: Vec<(String, LightMember)> = Vec::new();
        for clone in &mut self.light_clones {
            let depths: Vec<usize> = clone.members.iter().map(|x| x.mutations.mutation_events()).collect();
            let cutoff = robust_upper_cutoff(&depths, cfg);
            let Some(cutoff) = cutoff else { continue; };

            let source = clone.name.clone();
            let mut keep = Vec::new();
            for member in std::mem::take(&mut clone.members) {
                if member.mutations.mutation_events() > cutoff {
                    ejected.push((source.clone(), member));
                } else {
                    keep.push(member);
                }
            }
            clone.members = keep;
        }
        self.light_clones.retain(|clone| !clone.members.is_empty());

        // Hardest cases first. This prevents an easy low-mutation receptor from
        // broadening a target clone's robust cutoff before a difficult receptor
        // is evaluated against it.
        ejected.sort_by(|a, b| b.1.mutations.mutation_events().cmp(&a.1.mutations.mutation_events())
            .then_with(|| a.1.cell.cmp(&b.1.cell)));

        for (source, member) in ejected {
            if let Some((i, _membership)) = best_light_clone(
                &self.light_clones,
                &member.receptor,
                cfg,
                Some(source.as_str()),
                true,
                false,
                novel_v,
            ) {
                let mutations = measure_receptor(&member.receptor, novel_v);
                self.light_clones[i].members.push(LightMember {
                    cell: member.cell,
                    receptor: member.receptor,
                    mutations,
                });
            } else {
                self.light_clones.push(new_light_clone(member.cell, member.receptor, novel_v));
            }
        }
    }

    pub fn provisional_len(&self) -> usize { self.candidates.len() }

    pub fn worst_lc_alignments(&self, limit: usize, novel_v: &NovelVRegistries) -> Vec<WorstLightAlignment> {
        let mut rows: Vec<WorstLightAlignment> = self.light_clones.iter()
            .flat_map(|clone| clone.members.iter().filter_map(move |member| {
                let (expected, observed, differences) = render_fragment_alignment(&novel_v.effective_naive(&member.receptor), &member.receptor.observed)?;
                Some(WorstLightAlignment {
                    cell: member.cell.clone(),
                    clone_name: clone.name.clone(),
                    receptor_id: member.receptor.id.clone(),
                    mutation_events: member.mutations.mutation_events(),
                    substitutions: member.mutations.substitutions,
                    indel_events: member.mutations.indels.len(),
                    informative_pairs: member.mutations.informative_pairs,
                    matching_pairs: member.mutations.matching_pairs,
                    identity: member.mutations.identity,
                    expected,
                    observed,
                    differences,
                })
            }))
            .collect();
        rows.sort_by(|a, b| b.mutation_events.cmp(&a.mutation_events)
            .then_with(|| a.cell.cmp(&b.cell))
            .then_with(|| a.receptor_id.cmp(&b.receptor_id)));
        rows.truncate(limit);
        rows
    }

    pub fn mutation_report(&self) -> FamilyMutationReport {
        FamilyMutationReport {
            family: self.name.clone(),
            cells: self.members.len(),
            max_hc_mutations: self.max_hc_mutations,
            lc_clones: self.light_clones.len(),
            lc_mutations: mutation_stats(
                self.light_clones
                    .iter()
                    .flat_map(|c| c.members.iter())
                    .map(|x| x.mutations.mutation_events()),
            ),
        }
    }
}

#[derive(Debug, Clone)]
pub struct WorstLightAlignment {
    pub cell: String,
    pub clone_name: String,
    pub receptor_id: String,
    pub mutation_events: usize,
    pub substitutions: usize,
    pub indel_events: usize,
    pub informative_pairs: usize,
    pub matching_pairs: usize,
    pub identity: f64,
    /// Gapped expected (naive) sequence for the selected semi-global placement.
    pub expected: String,
    /// Gapped observed sequence in the exact same alignment columns.
    pub observed: String,
    /// Same-width diagnostic line: spaces are matches, bases are substitutions/
    /// insertions, and '-' marks bases missing from the observation.
    pub differences: String,
}

#[derive(Debug, Clone)]
pub struct FamilyMutationReport {
    pub family: String,
    pub cells: usize,
    pub max_hc_mutations: Option<usize>,
    pub lc_clones: usize,
    pub lc_mutations: MutationStats,
}

#[derive(Debug, Clone, Default)]
pub struct MutationStats {
    pub n: usize,
    pub mean: Option<f64>,
    pub sd: Option<f64>,
    /// The three largest mutation depths, ordered from lower to higher.
    pub max3: Vec<usize>,
}

fn mutation_stats<I>(depths: I) -> MutationStats
where
    I: IntoIterator<Item = usize>,
{
    let mut values: Vec<usize> = depths.into_iter().collect();
    if values.is_empty() { return MutationStats::default(); }

    let n = values.len();
    let mean = values.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
    // Population SD: this describes all LC observations in this finalized HC family.
    let variance = values
        .iter()
        .map(|&x| { let d = x as f64 - mean; d * d })
        .sum::<f64>() / n as f64;
    let sd = variance.sqrt();

    values.sort_unstable();
    let keep_from = values.len().saturating_sub(3);
    let max3 = values[keep_from..].to_vec();

    MutationStats { n, mean: Some(mean), sd: Some(sd), max3 }
}

fn new_light_clone(cell: String, receptor: Receptor, novel_v: &NovelVRegistries) -> LightClone {
    let measurement = measure_receptor(&receptor, novel_v);
    /* fallback retained inside measure_receptor */
    /*
        substitutions: 0,
        indels: Vec::new(),
        informative_pairs: 0,
        matching_pairs: 0,
    */
    let name = format!("LC:{}:{}:{}:{}", receptor.chain, receptor.v, receptor.j, receptor.cdr3_nt);
    LightClone {
        name,
        root: receptor.clone(),
        members: vec![LightMember { cell, receptor, mutations: measurement }],
    }
}

/// Find the best existing LC clone for a receptor. Mutation measurement is
/// clone-relative: candidate observed sequence versus the clone root's naive
/// sequence. This is the key distinction that makes LC redistribution real
/// rather than re-measuring the same receptor against itself for every clone.
fn best_light_clone(
    clones: &[LightClone],
    receptor: &Receptor,
    cfg: &FamilyConfig,
    exclude_name: Option<&str>,
    enforce_target_cutoff: bool,
    observed_seed: bool,
    novel_v: &NovelVRegistries,
) -> Option<(usize, MutationMeasurement)> {
    let mut best: Option<(usize, usize, MutationMeasurement)> = None;
    for (i, clone) in clones.iter().enumerate() {
        if exclude_name.is_some_and(|name| clone.name == name) { continue; }
        if !same_lc_labels(&clone.root, receptor) { continue; }
        let cdr3_distance = edit_distance(&clone.root.cdr3_nt, &receptor.cdr3_nt);
        if cdr3_distance > cfg.max_cdr3_distance { continue; }

        let effective_root_naive;
        let reference = if observed_seed {
            clone.root.observed.as_str()
        } else {
            effective_root_naive = novel_v.effective_naive(&clone.root);
            effective_root_naive.as_str()
        };
        let Some(measurement) = align_fragment(reference, &receptor.observed) else { continue; };
        if measurement.identity < cfg.hard_alignment_identity { continue; }

        if enforce_target_cutoff {
            let depths: Vec<usize> = clone.members.iter().map(|x| x.mutations.mutation_events()).collect();
            if let Some(cutoff) = robust_upper_cutoff(&depths, cfg) {
                if measurement.mutation_events() > cutoff { continue; }
            }
        }

        let replace = best.as_ref().is_none_or(|(_, old_distance, old)| {
            (measurement.mutation_events(), cdr3_distance, std::cmp::Reverse(measurement.matching_pairs))
                < (old.mutation_events(), *old_distance, std::cmp::Reverse(old.matching_pairs))
        });
        if replace { best = Some((i, cdr3_distance, measurement)); }
    }
    best.map(|(i, _, measurement)| (i, measurement))
}

fn measure_receptor(receptor: &Receptor, novel_v: &NovelVRegistries) -> MutationMeasurement {
    align_fragment(&novel_v.effective_naive(receptor), &receptor.observed).unwrap_or(MutationMeasurement {
        substitutions: 0, indels: Vec::new(), informative_pairs: 0, matching_pairs: 0, identity: 0.0,
    })
}

fn same_hc_labels(a: &Receptor, b: &Receptor) -> bool { a.chain == "IGH" && b.chain == "IGH" && a.v == b.v && a.j == b.j }
fn same_lc_labels(a: &Receptor, b: &Receptor) -> bool { a.chain == b.chain && matches!(a.chain.as_str(), "IGK"|"IGL") && a.v == b.v && a.j == b.j }

fn robust_upper_cutoff(depths: &[usize], cfg: &FamilyConfig) -> Option<usize> {
    if depths.len() < 4 { return None; }
    let med = median(depths);
    let deviations: Vec<usize> = depths.iter().map(|x| x.abs_diff(med)).collect();
    let mad = median(&deviations);
    Some(med + cfg.mutation_min_extra.max(mad.saturating_mul(cfg.mutation_mad_multiplier)))
}
fn median(xs: &[usize]) -> usize { let mut v=xs.to_vec(); v.sort_unstable(); v[v.len()/2] }

fn edit_distance(a: &str, b: &str) -> usize {
    let (a,b)=(a.as_bytes(),b.as_bytes());
    let mut prev: Vec<usize>=(0..=b.len()).collect(); let mut cur=vec![0;b.len()+1];
    for (i,&x) in a.iter().enumerate(){ cur[0]=i+1; for (j,&y) in b.iter().enumerate(){ cur[j+1]=(prev[j]+usize::from(!x.eq_ignore_ascii_case(&y))).min(prev[j+1]+1).min(cur[j]+1); } std::mem::swap(&mut prev,&mut cur); }
    prev[b.len()]
}


const MAX_SHORT_INDEL: usize = 3;
const INDEL_RESCUE_ANCHOR: usize = 4;

fn rescue_anchor_matches(a: &[u8], b: &[u8], i: usize, j: usize) -> bool {
    let available = (a.len() - i).min(b.len() - j);
    if available < INDEL_RESCUE_ANCHOR { return false; }
    (0..INDEL_RESCUE_ANCHOR).all(|k| a[i+k].eq_ignore_ascii_case(&b[j+k]))
}

fn walk_in_register(a: &[u8], b: &[u8]) -> Vec<(Option<usize>, Option<usize>)> {
    let (mut i, mut j) = (0usize, 0usize);
    let mut path = Vec::new();
    while i < a.len() && j < b.len() {
        if a[i].eq_ignore_ascii_case(&b[j]) {
            path.push((Some(i), Some(j))); i += 1; j += 1; continue;
        }

        let ins = (1..=MAX_SHORT_INDEL).find(|&d| j+d < b.len() && rescue_anchor_matches(a,b,i,j+d));
        let del = (1..=MAX_SHORT_INDEL).find(|&d| i+d < a.len() && rescue_anchor_matches(a,b,i+d,j));
        match (ins, del) {
            (Some(di), Some(dd)) if di <= dd => {
                for x in 0..di { path.push((None, Some(j+x))); }
                j += di;
            }
            (Some(_), Some(dd)) | (None, Some(dd)) => {
                for x in 0..dd { path.push((Some(i+x), None)); }
                i += dd;
            }
            (Some(di), None) => {
                for x in 0..di { path.push((None, Some(j+x))); }
                j += di;
            }
            (None, None) => {
                path.push((Some(i), Some(j))); i += 1; j += 1;
            }
        }
    }
    path
}

/// Build a mutation-measurement path without allowing a global aligner to
/// scatter terminal sequence through the opposite receptor as isolated gaps.
/// The longest exact common block anchors the two receptor observations. From
/// that block we walk outwards in register. A 1..=3 bp skip is accepted as an
/// indel only when it restores four consecutive matching bases; otherwise the
/// disagreement is a substitution. Sequence left after either side ends is
/// coverage outside the mutually observed span and is ignored.
pub(crate) fn fragment_path(reference: &[u8], observed: &[u8]) -> Option<Vec<(Option<usize>, Option<usize>)>> {
    if reference.is_empty() || observed.is_empty() { return None; }

    // Longest exact common substring: an unambiguous seed inside the receptor.
    let mut prev = vec![0usize; observed.len()+1];
    let mut best_len=0usize; let mut best_a_end=0usize; let mut best_b_end=0usize;
    for i in 1..=reference.len() {
        let mut cur=vec![0usize; observed.len()+1];
        for j in 1..=observed.len() {
            if reference[i-1].eq_ignore_ascii_case(&observed[j-1]) {
                cur[j]=prev[j-1]+1;
                if cur[j] > best_len { best_len=cur[j]; best_a_end=i; best_b_end=j; }
            }
        }
        prev=cur;
    }
    if best_len == 0 { return None; }
    let seed_a=best_a_end-best_len; let seed_b=best_b_end-best_len;

    // Walk the left side in reverse, then map it back to forward coordinates.
    let ar: Vec<u8> = reference[..seed_a].iter().rev().copied().collect();
    let br: Vec<u8> = observed[..seed_b].iter().rev().copied().collect();
    let mut left = walk_in_register(&ar,&br).into_iter().map(|(ai,bj)| {
        (ai.map(|x| seed_a-1-x), bj.map(|x| seed_b-1-x))
    }).collect::<Vec<_>>();
    left.reverse();

    let mut path=left;
    for k in 0..best_len { path.push((Some(seed_a+k),Some(seed_b+k))); }

    let right = walk_in_register(&reference[best_a_end..], &observed[best_b_end..]);
    path.extend(right.into_iter().map(|(ai,bj)| {
        (ai.map(|x| best_a_end+x), bj.map(|x| best_b_end+x))
    }));
    Some(path)
}

/// Render the exact semi-global placement used by mutation measurement.
/// The first line is the expected/naive sequence. The second line intentionally
/// suppresses matching bases: only substitutions, insertions, and deletions are
/// visible, at exactly the same columns as the expected sequence.
fn render_fragment_alignment(reference: &str, observed: &str) -> Option<(String, String, String)> {
    let a = reference.as_bytes();
    let b = observed.as_bytes();
    let path = fragment_path(a, b)?;

    let mut expected = String::with_capacity(path.len());
    let mut observed_line = String::with_capacity(path.len());
    let mut differences = String::with_capacity(path.len());
    for &(ai, bj) in &path {
        match (ai, bj) {
            (Some(ai), Some(bj)) => {
                let x = a[ai].to_ascii_uppercase() as char;
                let y = b[bj].to_ascii_uppercase() as char;
                expected.push(x);
                observed_line.push(y);
                differences.push(if x.eq_ignore_ascii_case(&y) { ' ' } else { y });
            }
            (Some(ai), None) => {
                expected.push(a[ai].to_ascii_uppercase() as char);
                observed_line.push('-');
                differences.push('-');
            }
            (None, Some(bj)) => {
                expected.push('-');
                observed_line.push(b[bj].to_ascii_uppercase() as char);
                differences.push(b[bj].to_ascii_uppercase() as char);
            }
            (None, None) => unreachable!(),
        }
    }
    Some((expected, observed_line, differences))
}

/// Compatibility helper: expected sequence plus differences-only line.
pub fn render_fragment_differences(reference: &str, observed: &str) -> Option<(String, String)> {
    let (expected, _observed, differences) = render_fragment_alignment(reference, observed)?;
    Some((expected, differences))
}

/// Fragment-aware mutation measurement for partial receptor observations.
///
/// A longest exact common block anchors the two sequences. From that anchor we
/// walk towards both ends in register. Short (1..=3 bp) indels are accepted
/// only when they restore four consecutive matches; otherwise disagreement is
/// counted as substitutions. Once either sequence ends, the remaining tail is
/// coverage outside the mutually observed span and is ignored.
///
/// This is mutation measurement only, never the initial structural family gate.
pub fn align_fragment(reference: &str, observed: &str) -> Option<MutationMeasurement> {
    let a=reference.as_bytes(); let b=observed.as_bytes();
    let path=fragment_path(a,b)?;

    let mut informative=0; let mut matches=0; let mut substitutions=0; let mut indels=Vec::new(); let mut k=0;
    while k<path.len(){match path[k]{(Some(ai),Some(bj))=>{let x=a[ai].to_ascii_uppercase();let y=b[bj].to_ascii_uppercase();if x!=b'N'&&y!=b'N'{informative+=1;if x==y{matches+=1}else{substitutions+=1}}k+=1},(Some(ai),None)=>{let start=ai;let mut n=0;while k<path.len()&&matches!(path[k],(Some(_),None)){n+=1;k+=1}indels.push(IndelEvent{naive_pos:start,inserted:false,len:n})},(None,Some(_))=>{let pos=path[..k].iter().rev().find_map(|x|x.0).map_or(0,|x|x+1);let mut n=0;while k<path.len()&&matches!(path[k],(None,Some(_))){n+=1;k+=1}indels.push(IndelEvent{naive_pos:pos,inserted:true,len:n})},(None,None)=>unreachable!()}}
    if informative==0{return None} let identity=matches as f64/informative as f64;
    Some(MutationMeasurement{substitutions,indels,informative_pairs:informative,matching_pairs:matches,identity})
}

#[cfg(test)]
mod tests {
    use super::*;
    fn r(id:&str,chain:&str,v:&str,j:&str,cdr3:&str,naive:&str,obs:&str)->Receptor{Receptor{id:id.into(),chain:chain.into(),v:v.into(),d:String::new(),j:j.into(),c:String::new(),cdr3_nt:cdr3.into(),naive:naive.into(),observed:obs.into(),alternative_reconstructions:0,evidence_reads:2}}
    fn cell(name:&str,cdr3:&str,obs:&str)->CellReceptor{CellReceptor{cell_id:name.into(),hc:r(name,"IGH","V1","J1",cdr3,"AAAACCCCGGGGTTTT",obs),lc:vec![]}}
    #[test] fn initial_collection_is_cdr3_complete_link(){let cfg=FamilyConfig::default();let mut f=Family::new("F",cell("a","AAAA","AAAACCCCGGGGTTTT"));assert!(f.add_candidate(cell("b","AAAT","AAAACCCCGGGGTTTT"),&cfg).is_ok());assert!(f.add_candidate(cell("c","TTTT","AAAACCCCGGGGTTTT"),&cfg).is_err());}
    #[test] fn mutation_stats_reports_mean_sd_and_top_three(){
        let s = mutation_stats([1, 2, 3, 4, 53]);
        assert_eq!(s.n, 5);
        assert!((s.mean.unwrap() - 12.6).abs() < 1e-9);
        assert!((s.sd.unwrap() - 20.22473732833136).abs() < 1e-9);
        assert_eq!(s.max3, vec![3, 4, 53]);
    }

    #[test] fn best_light_clone_uses_target_root_not_receptors_own_naive(){
        let cfg=FamilyConfig::default();
        let root_a=r("a","IGK","VK","JK","AAAA","AAAACCCCGGGGTTTT","AAAACCCCGGGGTTTT");
        let root_b=r("b","IGK","VK","JK","AAAT","TTTTCCCCGGGGAAAA","TTTTCCCCGGGGAAAA");
        let receptor=r("x","IGK","VK","JK","AAAT","AAAACCCCGGGGTTTT","TTTTCCCCGGGGAAAA");
        let clones=vec![
            new_light_clone("a".into(),root_a,&NovelVRegistries::default()),
            new_light_clone("b".into(),root_b,&NovelVRegistries::default()),
        ];
        let (i,m)=best_light_clone(&clones,&receptor,&cfg,None,false,false,&NovelVRegistries::default()).unwrap();
        assert_eq!(i,1);
        assert_eq!(m.mutation_events(),0);
        assert_eq!(m.identity,1.0);
    }

    #[test] fn fragment_alignment_ignores_reference_ends(){let m=align_fragment("TTTTAAAACCCCGGGGAAAA","AAAACCCCGGGG").unwrap();assert_eq!(m.mutation_events(),0);assert_eq!(m.identity,1.0);}
    #[test] fn fragment_alignment_ignores_observed_ends(){let m=align_fragment("AAAACCCCGGGG","TTTTAAAACCCCGGGGAAAA").unwrap();assert_eq!(m.mutation_events(),0);assert_eq!(m.identity,1.0);assert_eq!(m.informative_pairs,12);}
    #[test] fn fragment_alignment_keeps_internal_indels(){let m=align_fragment("AAAACCCCGGGG","AAAACCCCAGGGG").unwrap();assert_eq!(m.indels.len(),1);assert_eq!(m.identity,1.0);}
    #[test] fn fragment_alignment_counts_supported_prefix_disagreement(){
        let m=align_fragment("ATGATGAAAACCCCGGGG", "TATCGATCAAAAACCCCGGGG").unwrap();
        assert_eq!(m.mutation_events(), 6);
        assert_eq!(m.informative_pairs, 18);
    }
    #[test] fn fragment_alignment_keeps_single_internal_substitution(){
        let m=align_fragment("AAAACCCCGGGG","AAAATCCCGGGG").unwrap();
        assert_eq!(m.substitutions,1);
        assert!(m.indels.is_empty());
        assert_eq!(m.informative_pairs,12);
    }
    #[test] fn fragment_alignment_rescues_three_base_internal_insertion(){
        let m=align_fragment("AAAACCCCGGGGTTTT","AAAACCCAAACGGGGTTTT").unwrap();
        assert_eq!(m.substitutions,0);
        assert_eq!(m.indels.len(),1);
        assert_eq!(m.indels[0].len,3);
    }
    #[test] fn fragment_alignment_does_not_fish_matches_from_terminal_extension(){
        let m=align_fragment("AAAACCCCGGGGAAAA","AAAACCCCGGGGTTTTCCCC").unwrap();
        assert_eq!(m.substitutions,4);
        assert!(m.indels.is_empty());
        assert_eq!(m.informative_pairs,16);
    }
    #[test] fn difference_render_hides_matches(){
        let (expected, diff)=render_fragment_differences("AAAACCCCGGGG","AAAATCCCGGGG").unwrap();
        assert_eq!(expected, "AAAACCCCGGGG");
        assert_eq!(diff, "    T       ");
    }
}
