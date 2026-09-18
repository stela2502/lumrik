use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct Call {
    pub cell: String,
    pub id: String,
    pub chain: String,
    pub v: String,
    pub d: String,
    pub j: String,
    pub c: String,
    pub pn_alternative: bool,
    pub v_del_3: u16,
    pub j_del_5: u16,
    pub d_retained_len: u16,
    pub p_total_len: u16,
    pub productivity: String,
    pub support: u64,
    pub rediscovery: u64,
    pub naive: String,
    pub observed: String,
    pub cdr3_nt: String,
    pub cdr3_aa: String,
}
impl Call {
    pub fn productive(&self) -> bool {
        self.productivity == "productive"
    }
    pub fn heavy(&self) -> bool {
        matches!(self.chain.as_str(), "IGH" | "TRB" | "TRD")
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Mutation {
    pub pos: usize,
    pub from: u8,
    pub to: u8,
}

/// Measure substitutions after family membership has already been established.
///
/// Family admission is deliberately strict (V/J plus the configured CDR3 edit
/// gate). Mutation measurement has a different job: establish homologous
/// coordinates between the reconstructed naive receptor and an observed
/// receptor that may start/end at a slightly different position. A compact
/// Needleman-Wunsch edit alignment is therefore appropriate here.
///
/// Gaps participate in the alignment but are NOT counted as somatic mutations.
/// Only aligned, non-N nucleotide substitutions are emitted. This prevents a
/// one-base terminal truncation from shifting the whole receptor and creating
/// hundreds of false substitutions.
#[derive(Debug, Clone)]
pub struct MutationAlignment {
    /// Aligned nucleotide substitutions (SHM SNV events).
    pub mutations: Vec<Mutation>,
    /// Internal contiguous indel runs. Terminal sequence truncation is coverage,
    /// not an indel event, and is deliberately excluded.
    pub indels: Vec<IndelEvent>,
    aligned: Vec<(Option<usize>, Option<usize>)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndelKind {
    Insertion,
    Deletion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndelEvent {
    /// Naive-coordinate boundary immediately before the event.
    pub pos: usize,
    pub kind: IndelKind,
    pub len: usize,
}

impl MutationAlignment {
    /// Biological mutation-event depth: each substitution is one event and each
    /// contiguous internal indel, regardless of length, is one event.
    pub fn event_count(&self) -> usize {
        self.mutations.len() + self.indels.len()
    }

    pub fn inserted_nt(&self) -> usize {
        self.indels.iter().filter(|x| x.kind == IndelKind::Insertion).map(|x| x.len).sum()
    }

    pub fn deleted_nt(&self) -> usize {
        self.indels.iter().filter(|x| x.kind == IndelKind::Deletion).map(|x| x.len).sum()
    }
}

pub fn measure_mutations_nw(naive: &str, observed: &str) -> Option<MutationAlignment> {
    let a = naive.as_bytes();
    let b = observed.as_bytes();
    if a.is_empty() || b.is_empty() {
        return None;
    }

    // MUTATION MEASUREMENT ONLY.
    //
    // `naive` is commonly the complete reconstructed/germline receptor while
    // `observed` may only cover an internal fragment. A global alignment would
    // force the complete germline ends through that fragment and can turn
    // missing coverage into tens/hundreds of apparent mutations. Use a
    // semi-global Needleman-Wunsch alignment instead: the complete observed
    // fragment must align, but unused naive prefix/suffix is free terminal
    // coverage. Family admission is already frozen before this function runs.
    let cols = b.len() + 1;
    let mut score = vec![0u32; (a.len() + 1) * cols];

    // Any naive prefix may be skipped for free; observed sequence may not be
    // skipped. This lets the alignment discover where the fragment starts in
    // the complete germline/reconstructed receptor.
    for i in 0..=a.len() {
        score[i * cols] = 0;
    }
    for j in 0..=b.len() {
        score[j] = j as u32;
    }

    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let subst = score[(i - 1) * cols + (j - 1)]
                + u32::from(!a[i - 1].eq_ignore_ascii_case(&b[j - 1]));
            let delete = score[(i - 1) * cols + j] + 1;
            let insert = score[i * cols + (j - 1)] + 1;
            score[i * cols + j] = subst.min(delete).min(insert);
        }
    }

    // The observed fragment may finish before the complete naive receptor.
    // Pick the best endpoint after consuming all observed bases.
    let mut end_i = 0usize;
    let mut end_score = score[b.len()];
    for i in 1..=a.len() {
        let candidate = score[i * cols + b.len()];
        if candidate < end_score {
            end_score = candidate;
            end_i = i;
        }
    }

    let mut aligned = Vec::with_capacity(a.len().max(b.len()));

    // Naive suffix beyond the observed fragment is terminal missing coverage.
    for ai in (end_i..a.len()).rev() {
        aligned.push((Some(ai), None));
    }

    let (mut i, mut j) = (end_i, b.len());
    while j > 0 {
        if i > 0 {
            let subst_cost = u32::from(!a[i - 1].eq_ignore_ascii_case(&b[j - 1]));
            if score[i * cols + j] == score[(i - 1) * cols + (j - 1)] + subst_cost {
                aligned.push((Some(i - 1), Some(j - 1)));
                i -= 1;
                j -= 1;
                continue;
            }
            if score[i * cols + j] == score[(i - 1) * cols + j] + 1 {
                aligned.push((Some(i - 1), None));
                i -= 1;
                continue;
            }
        }
        aligned.push((None, Some(j - 1)));
        j -= 1;
    }

    // Naive prefix before the fragment is terminal missing coverage.
    while i > 0 {
        aligned.push((Some(i - 1), None));
        i -= 1;
    }
    aligned.reverse();

    // Refuse to manufacture a mutation distance for an unrelated/very poor
    // match. Ns are coverage uncertainty and do not contribute to this check.
    let mut informative = 0usize;
    let mut matches = 0usize;
    for (ai, bj) in aligned.iter().copied() {
        let (Some(ai), Some(bj)) = (ai, bj) else { continue };
        let from = a[ai].to_ascii_uppercase();
        let to = b[bj].to_ascii_uppercase();
        if from == b'N' || to == b'N' {
            continue;
        }
        informative += 1;
        if from == to {
            matches += 1;
        }
    }
    if informative == 0 || matches * 2 < informative {
        return None;
    }

    let mutations = aligned
        .iter()
        .copied()
        .filter_map(|(ai, bj)| {
            let (ai, bj) = (ai?, bj?);
            let from = a[ai].to_ascii_uppercase();
            let to = b[bj].to_ascii_uppercase();
            (from != to && from != b'N' && to != b'N').then_some(Mutation { pos: ai, from, to })
        })
        .collect();

    // Convert internal gap runs into indel events. Leading/trailing gaps are
    // incomplete receptor coverage and must not masquerade as SHM indels.
    let first_paired = aligned.iter().position(|(ai, bj)| ai.is_some() && bj.is_some());
    let last_paired = aligned.iter().rposition(|(ai, bj)| ai.is_some() && bj.is_some());
    let mut indels = Vec::new();
    if let (Some(first), Some(last)) = (first_paired, last_paired) {
        let mut k = first + 1;
        while k < last {
            let kind = match aligned[k] {
                (None, Some(_)) => Some(IndelKind::Insertion),
                (Some(_), None) => Some(IndelKind::Deletion),
                _ => None,
            };
            let Some(kind) = kind else { k += 1; continue };
            let start = k;
            while k < last {
                let same = matches!((kind, aligned[k]),
                    (IndelKind::Insertion, (None, Some(_))) |
                    (IndelKind::Deletion, (Some(_), None)));
                if !same { break; }
                k += 1;
            }
            let pos = aligned[..start]
                .iter()
                .rev()
                .find_map(|(ai, _)| *ai)
                .map(|x| x + 1)
                .unwrap_or(0);
            indels.push(IndelEvent { pos, kind, len: k - start });
        }
    }

    Some(MutationAlignment { mutations, indels, aligned })
}

#[derive(Debug, Clone)]
struct HcCellRejection {
    cell: String,
    reason: &'static str,
    valid_hc_ids: Vec<String>,
}

/// Valkyrn's entry contract for B-cell lineage analysis: exactly one
/// productive IGH reconstruction with a biologically established AIRR CDR3
/// per cell. Cells with zero or multiple such heavy chains are excluded before
/// family construction; their light-chain calls must not leak into an
/// unqualified HC background.
///
/// TODO(later): investigate whether some multiple-HC cells are duplicate
/// reconstructions of one biological HC. Do not collapse/rescue them here;
/// qualification stays deliberately literal until that behavior is proven.
fn qualify_cells_by_hc(calls: Vec<Call>) -> (Vec<Call>, Vec<HcCellRejection>, usize) {
    let mut hc_by_cell: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut all_cells: BTreeSet<String> = BTreeSet::new();
    for c in &calls {
        all_cells.insert(c.cell.clone());
        if c.productive() && c.chain == "IGH" && c.id.starts_with("HC:") && !c.cdr3_nt.is_empty() {
            hc_by_cell
                .entry(c.cell.clone())
                .or_default()
                .push(c.id.clone());
        }
    }

    let mut accepted = BTreeSet::new();
    let mut rejected = Vec::new();
    for cell in &all_cells {
        let ids = hc_by_cell.get(cell).cloned().unwrap_or_default();
        match ids.len() {
            1 => { accepted.insert(cell.clone()); }
            0 => rejected.push(HcCellRejection {
                cell: cell.clone(),
                reason: "no_valid_hc",
                valid_hc_ids: ids,
            }),
            _ => rejected.push(HcCellRejection {
                cell: cell.clone(),
                reason: "multiple_valid_hc",
                valid_hc_ids: ids,
            }),
        }
    }

    let qualified = calls
        .into_iter()
        .filter(|c| accepted.contains(&c.cell))
        .collect();
    (qualified, rejected, all_cells.len())
}

fn write_hc_rejections(out: &Path, rejected: &[HcCellRejection]) -> Result<()> {
    let mut w = writer(out.join("valkyrn_rejected_cells.tsv"))?;
    writeln!(w, "cell\treason\tvalid_hc_count\tvalid_hc_ids")?;
    for r in rejected {
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            r.cell,
            r.reason,
            r.valid_hc_ids.len(),
            r.valid_hc_ids.join(",")
        )?;
    }
    Ok(())
}

fn edit_distance(a: &str, b: &str) -> usize {
    let a = a.as_bytes();
    let b = b.as_bytes();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0usize; b.len() + 1];
    for (i, &x) in a.iter().enumerate() {
        cur[0] = i + 1;
        for (j, &y) in b.iter().enumerate() {
            cur[j + 1] = (prev[j + 1] + 1)
                .min(cur[j] + 1)
                .min(prev[j] + usize::from(x != y));
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

pub fn read_calls(vdj_dir: &Path) -> Result<Vec<Call>> {
    let calls_path = vdj_dir.join("vdj_calls.tsv");
    let airr_path = vdj_dir.join("airr_rearrangements.tsv");
    let airr = read_airr_cdr3(&airr_path)?;
    let file =
        File::open(&calls_path).with_context(|| format!("opening {}", calls_path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines.next().context("vdj_calls.tsv is empty")??;
    let h: Vec<&str> = header.split('\t').collect();
    let ix = |name: &str| {
        h.iter()
            .position(|x| *x == name)
            .with_context(|| format!("vdj_calls.tsv lacks column {name}"))
    };
    let cell = ix("cell")?;
    let id = ix("recombination_id")?;
    let chain = ix("chain")?;
    let v = ix("v")?;
    let d = ix("d")?;
    let j = ix("j")?;
    let c = ix("c")?;
    let pn_alternative = ix("pn_alternative")?;
    let v_del_3 = ix("v_del_3")?;
    let j_del_5 = ix("j_del_5")?;
    let d_retained_len = ix("d_retained_len")?;
    let p_v3_len = ix("p_v3_len")?;
    let p_d5_len = ix("p_d5_len")?;
    let p_d3_len = ix("p_d3_len")?;
    let p_j5_len = ix("p_j5_len")?;
    let productivity = ix("productivity_status")?;
    let support = ix("support_features")?;
    let rediscovery = ix("receptor_rediscovery_reads")?;
    let naive = ix("naive_recombination")?;
    let observed = ix("observed_receptor_sequence")?;
    let mut out = Vec::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let rid = get(id).to_string();
        out.push(Call {
            cell: get(cell).into(),
            id: rid.clone(),
            chain: get(chain).into(),
            v: get(v).into(),
            d: get(d).into(),
            j: get(j).into(),
            c: get(c).into(),
            pn_alternative: get(pn_alternative).eq_ignore_ascii_case("true"),
            v_del_3: get(v_del_3).parse().unwrap_or(0),
            j_del_5: get(j_del_5).parse().unwrap_or(0),
            d_retained_len: get(d_retained_len).parse().unwrap_or(0),
            p_total_len: [p_v3_len, p_d5_len, p_d3_len, p_j5_len]
                .into_iter()
                .filter_map(|col| get(col).parse::<u16>().ok())
                .sum(),
            productivity: get(productivity).into(),
            support: get(support).parse().unwrap_or(0),
            rediscovery: get(rediscovery).parse().unwrap_or(0),
            naive: get(naive).into(),
            observed: get(observed).into(),
            cdr3_nt: airr.get(&rid).map(|x| x.0.clone()).unwrap_or_default(),
            cdr3_aa: airr.get(&rid).map(|x| x.1.clone()).unwrap_or_default(),
        });
    }
    Ok(out)
}

fn read_airr_cdr3(path: &Path) -> Result<HashMap<String, (String, String)>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines.next().context("airr_rearrangements.tsv is empty")??;
    let h: Vec<&str> = header.split('\t').collect();
    let seq = h
        .iter()
        .position(|x| *x == "lumrik_recombination_id")
        .context("AIRR lacks lumrik_recombination_id")?;
    let cdr_nt = h
        .iter()
        .position(|x| *x == "cdr3")
        .context("AIRR lacks cdr3")?;
    let cdr_aa = h
        .iter()
        .position(|x| *x == "cdr3_aa")
        .context("AIRR lacks cdr3_aa")?;
    let mut out = HashMap::new();
    for line in lines {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if let (Some(id), Some(nt), Some(aa)) = (f.get(seq), f.get(cdr_nt), f.get(cdr_aa)) {
            out.insert((*id).into(), ((*nt).into(), (*aa).into()));
        }
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct Family {
    pub name: String,
    pub members: Vec<usize>,
}

/// Hard CDR3 family gate. Structural IDs, P/N decomposition and mutation
/// measurements are provenance only and may never override this distance.
fn cdr3_compatible(a: &str, b: &str, max_cdr3_distance: usize) -> bool {
    !a.is_empty() && !b.is_empty() && edit_distance(a, b) <= max_cdr3_distance
}

/// Deterministic complete-link clustering over a hard compatibility predicate.
///
/// A candidate cluster is merged only when *every* cross-cluster pair is
/// compatible.  This deliberately rejects transitive chaining (A~B, B~C, but
/// A!~C), because that would put a <75%-matching CDR3 pair into the same family.
fn complete_link_groups<F>(indices: &[usize], compatible: F) -> Vec<Vec<usize>>
where
    F: Fn(usize, usize) -> bool,
{
    let mut groups: Vec<Vec<usize>> = indices.iter().copied().map(|i| vec![i]).collect();
    loop {
        let mut merge = None;
        'outer: for left in 0..groups.len() {
            for right in left + 1..groups.len() {
                if groups[left]
                    .iter()
                    .all(|&a| groups[right].iter().all(|&b| compatible(a, b)))
                {
                    merge = Some((left, right));
                    break 'outer;
                }
            }
        }
        let Some((left, right)) = merge else { break };
        let other = groups.remove(right);
        groups[left].extend(other);
    }
    groups
}

/// Fraction of informative aligned nucleotide pairs that agree. Gaps and Ns do
/// not improve this score. `measure_mutations_nw` already applies the permissive
/// >=50% fail-safe; reassignment deliberately uses the harder >=75% gate.
fn mutation_alignment_identity(calls: &[Call], reference_i: usize, call_i: usize) -> Option<f64> {
    let aln = measure_mutations_nw(&calls[reference_i].naive, &calls[call_i].observed)?;
    let a = calls[reference_i].naive.as_bytes();
    let b = calls[call_i].observed.as_bytes();
    let mut informative = 0usize;
    let mut matches = 0usize;
    for (ai, bj) in aln.aligned.iter().copied() {
        let (Some(ai), Some(bj)) = (ai, bj) else { continue };
        let from = a[ai].to_ascii_uppercase();
        let to = b[bj].to_ascii_uppercase();
        if from == b'N' || to == b'N' { continue; }
        informative += 1;
        matches += usize::from(from == to);
    }
    (informative > 0).then_some(matches as f64 / informative as f64)
}

fn family_reference_index(calls: &[Call], members: &[usize]) -> Option<usize> {
    members.iter().copied().max_by_key(|&i| {
        (
            calls[i].observed.len(),
            calls[i].support.saturating_add(calls[i].rediscovery),
        )
    })
}

fn structurally_compatible_with_family(
    calls: &[Call],
    call_i: usize,
    members: &[usize],
    max_cdr3_distance: usize,
) -> bool {
    let c = &calls[call_i];
    members.iter().all(|&j| {
        let x = &calls[j];
        c.chain == x.chain
            && c.v == x.v
            && c.j == x.j
            && cdr3_compatible(&c.cdr3_nt, &x.cdr3_nt, max_cdr3_distance)
    })
}

/// Second-stage HC validation. Initial families are created exclusively by the
/// hard structural gate. A member whose observed receptor cannot be measured
/// against its provisional family root is removed, tried against every other
/// structurally legal HC family, and re-added only when the hard mutation
/// identity gate succeeds. Failure leaves the HC unassigned; it never becomes
/// plot input merely because a permissive NW path exists.
fn refine_heavy_families(
    calls: &[Call],
    mut fams: Vec<Family>,
    max_cdr3_distance: usize,
) -> (Vec<Family>, Vec<usize>) {
    let roots: Vec<Option<usize>> = fams
        .iter()
        .map(|f| family_reference_index(calls, &f.members))
        .collect();
    let mut failed = Vec::new();

    // First pass is deliberately permissive: only a genuinely unmeasurable
    // member is evicted from its provisional family.
    for fi in 0..fams.len() {
        let Some(root) = roots[fi] else { continue };
        let old = std::mem::take(&mut fams[fi].members);
        for i in old {
            if i == root || mutation_alignment_identity(calls, root, i).is_some() {
                fams[fi].members.push(i);
            } else {
                failed.push(i);
            }
        }
    }

    // Retry evicted HCs against every other legal family. Re-admission is hard:
    // complete-link structural compatibility plus >=75% informative identity.
    for i in failed.iter().copied() {
        let mut best: Option<(usize, f64)> = None;
        for fi in 0..fams.len() {
            if fams[fi].members.is_empty()
                || !structurally_compatible_with_family(calls, i, &fams[fi].members, max_cdr3_distance)
            {
                continue;
            }
            let Some(root) = family_reference_index(calls, &fams[fi].members) else { continue };
            let Some(identity) = mutation_alignment_identity(calls, root, i) else { continue };
            if identity < 0.75 { continue; }
            if best.map(|(_, x)| identity > x).unwrap_or(true) {
                best = Some((fi, identity));
            }
        }
        if let Some((fi, _)) = best {
            fams[fi].members.push(i);
        }
    }

    let assigned: BTreeSet<usize> = fams.iter().flat_map(|f| f.members.iter().copied()).collect();
    let unassigned = failed.into_iter().filter(|i| !assigned.contains(i)).collect();
    fams.retain(|f| !f.members.is_empty());
    for f in &mut fams {
        f.members.sort_unstable();
        let c = &calls[f.members[0]];
        f.name = format!("HC:{}:{}:CDR3:{}", c.v, c.j, family_cdr3(calls, &f.members));
    }
    (fams, unassigned)
}

fn family_cdr3(calls: &[Call], members: &[usize]) -> String {
    members
        .iter()
        .map(|&i| calls[i].cdr3_nt.as_str())
        .filter(|x| !x.is_empty())
        .min()
        .unwrap_or("MISSING")
        .to_string()
}

/// Heavy-chain families use only robust receptor labels plus the CDR3 sequence.
///
/// The old HC:<HEX>/P-N family merger is intentionally gone.  HC compact IDs
/// remain provenance in the output, but they no longer decide clonality.  Two
/// productive IGH calls can share a family only when V and J agree and their
/// CDR3 nucleotide sequences have edit distance within the configured maximum.  Complete-link grouping
/// makes the configured CDR3-distance rule true for every pair in the final family.
fn heavy_families(calls: &[Call], max_cdr3_distance: usize) -> Vec<Family> {
    // Bucket first by the hard V/J labels, then collapse identical CDR3s before
    // complete-link clustering.  The previous implementation fed every call to
    // the agglomerator individually.  Large expanded clones therefore contained
    // hundreds/thousands of duplicate CDR3 nodes and turned a tiny biological
    // comparison into a cubic amount of repeated cluster bookkeeping.
    //
    // Collapsing identical CDR3s is semantics-preserving: compatibility inside
    // a V/J bucket depends only on the CDR3 sequence, so all calls carrying the
    // same sequence are interchangeable for the hard configured CDR3-distance invariant.
    let mut buckets: BTreeMap<(String, String), BTreeMap<String, Vec<usize>>> = BTreeMap::new();
    for (i, c) in calls.iter().enumerate() {
        if c.productive() && c.chain == "IGH" && c.id.starts_with("HC:") && !c.cdr3_nt.is_empty() {
            buckets
                .entry((c.v.clone(), c.j.clone()))
                .or_default()
                .entry(c.cdr3_nt.clone())
                .or_default()
                .push(i);
        }
    }

    let mut out = Vec::new();
    for ((v, j), variants) in buckets {
        let variant_members: Vec<Vec<usize>> = variants.into_values().collect();
        let variant_indices: Vec<usize> = (0..variant_members.len()).collect();
        let groups = complete_link_groups(&variant_indices, |a, b| {
            let left = &calls[variant_members[a][0]].cdr3_nt;
            let right = &calls[variant_members[b][0]].cdr3_nt;
            cdr3_compatible(left, right, max_cdr3_distance)
        });

        for group in groups {
            let mut members = Vec::new();
            for variant in group {
                members.extend_from_slice(&variant_members[variant]);
            }
            members.sort_by(|&a, &b| {
                (&calls[a].cdr3_nt, &calls[a].cell, &calls[a].id)
                    .cmp(&(&calls[b].cdr3_nt, &calls[b].cell, &calls[b].id))
            });
            let cdr3 = family_cdr3(calls, &members);
            out.push(Family {
                name: format!("HC:{}:{}:CDR3:{}", v, j, cdr3),
                members,
            });
        }
    }
    out
}

/// Canonicalize productive light chains inside one HC background.
///
/// LC compact recombination IDs are provenance only.  Light-chain families
/// require the same light chain, V and J calls, and pairwise CDR3 edit distance within the configured maximum
/// nucleotide identity.  As for HC, complete-link grouping prevents chaining
/// across the hard CDR3 boundary.
fn canonical_light_ids_with_distance(calls: &[Call], indices: &[usize], max_cdr3_distance: usize) -> HashMap<usize, String> {
    // As for HC, cluster unique CDR3 variants rather than individual receptor
    // calls.  The HC background has already been fixed by the caller; chain/V/J
    // are the remaining hard labels and are used as buckets here.
    let mut buckets: BTreeMap<(String, String, String), BTreeMap<String, Vec<usize>>> =
        BTreeMap::new();
    for &i in indices {
        let c = &calls[i];
        if c.productive()
            && matches!(c.chain.as_str(), "IGK" | "IGL")
            && c.id.starts_with("LC:")
            && !c.cdr3_nt.is_empty()
        {
            buckets
                .entry((c.chain.clone(), c.v.clone(), c.j.clone()))
                .or_default()
                .entry(c.cdr3_nt.clone())
                .or_default()
                .push(i);
        }
    }

    let mut out = HashMap::new();
    for ((chain, v, j), variants) in buckets {
        let variant_members: Vec<Vec<usize>> = variants.into_values().collect();
        let variant_indices: Vec<usize> = (0..variant_members.len()).collect();
        let groups = complete_link_groups(&variant_indices, |a, b| {
            let left = &calls[variant_members[a][0]].cdr3_nt;
            let right = &calls[variant_members[b][0]].cdr3_nt;
            cdr3_compatible(left, right, max_cdr3_distance)
        });

        for group in groups {
            let mut members = Vec::new();
            for variant in group {
                members.extend_from_slice(&variant_members[variant]);
            }
            members.sort_by(|&a, &b| {
                (&calls[a].cdr3_nt, &calls[a].cell, &calls[a].id)
                    .cmp(&(&calls[b].cdr3_nt, &calls[b].cell, &calls[b].id))
            });
            let cdr3 = family_cdr3(calls, &members);
            let name = format!("LC:{}:{}:{}:CDR3:{}", chain, v, j, cdr3);
            for i in members {
                out.insert(i, name.clone());
            }
        }
    }
    out
}

fn canonical_light_ids(calls: &[Call], indices: &[usize]) -> HashMap<usize, String> {
    canonical_light_ids_with_distance(calls, indices, 3)
}

/// Build repertoire families with asymmetric HC/LC trust.
///
/// HC families are defined first.  LC families are then defined only within an
/// HC background; the same LC-like CDR3 on unrelated HC backgrounds is not one
/// clone.  Missing CDR3 calls are not family evidence and therefore stay out of
/// family clustering rather than bypassing the 75% invariant.
pub fn families(calls: &[Call], max_cdr3_distance: usize) -> Vec<Family> {
    let provisional_hc = heavy_families(calls, max_cdr3_distance);
    let (mut out, unassigned_hc) = refine_heavy_families(calls, provisional_hc, max_cdr3_distance);
    if !unassigned_hc.is_empty() {
        eprintln!(
            "Valkyrn: HC mutation validation left {} structurally qualified cells unassigned after hard reassignment",
            unassigned_hc.len()
        );
    }

    // LC is asymmetric by design. It lives only inside a finalized HC family.
    // If an LC does not measure against its provisional LC root, try every other
    // structurally legal LC clone in that HC background with the hard mutation
    // gate. If none fits, the fragment becomes a NEW LC clone root rather than
    // being discarded or contaminating another clone.
    let heavy_families_snapshot: Vec<(String, Vec<usize>)> = out
        .iter()
        .map(|f| (f.name.clone(), f.members.clone()))
        .collect();

    for (hc_name, hc_members) in heavy_families_snapshot {
        let cells: BTreeSet<&str> = hc_members.iter().map(|&i| calls[i].cell.as_str()).collect();
        let indices: Vec<usize> = calls.iter().enumerate().filter(|(_, c)| {
            cells.contains(c.cell.as_str())
                && c.productive()
                && matches!(c.chain.as_str(), "IGK" | "IGL")
                && c.id.starts_with("LC:")
                && !c.cdr3_nt.is_empty()
        }).map(|(i, _)| i).collect();

        let canonical = canonical_light_ids_with_distance(calls, &indices, max_cdr3_distance);
        let mut clones: Vec<Vec<usize>> = {
            let mut grouped: BTreeMap<String, Vec<usize>> = BTreeMap::new();
            for &i in &indices {
                if let Some(name) = canonical.get(&i) {
                    grouped.entry(name.clone()).or_default().push(i);
                }
            }
            grouped.into_values().collect()
        };

        let roots: Vec<Option<usize>> = clones.iter().map(|m| family_reference_index(calls, m)).collect();
        let mut failed = Vec::new();
        for ci in 0..clones.len() {
            let Some(root) = roots[ci] else { continue };
            let old = std::mem::take(&mut clones[ci]);
            for i in old {
                if i == root || mutation_alignment_identity(calls, root, i).is_some() {
                    clones[ci].push(i);
                } else {
                    failed.push(i);
                }
            }
        }

        for i in failed {
            let mut best: Option<(usize, f64)> = None;
            for ci in 0..clones.len() {
                if clones[ci].is_empty()
                    || !structurally_compatible_with_family(calls, i, &clones[ci], max_cdr3_distance)
                {
                    continue;
                }
                let Some(root) = family_reference_index(calls, &clones[ci]) else { continue };
                let Some(identity) = mutation_alignment_identity(calls, root, i) else { continue };
                if identity < 0.75 { continue; }
                if best.map(|(_, x)| identity > x).unwrap_or(true) {
                    best = Some((ci, identity));
                }
            }
            if let Some((ci, _)) = best {
                clones[ci].push(i);
            } else {
                // No existing LC clone is a convincing home: this fragment is
                // the root of a new LC clone in this finalized HC background.
                clones.push(vec![i]);
            }
        }

        for mut members in clones {
            if members.is_empty() { continue; }
            members.sort_unstable();
            let c = &calls[members[0]];
            let lc_name = format!(
                "LC:{}:{}:{}:CDR3:{}",
                c.chain,
                c.v,
                c.j,
                family_cdr3(calls, &members)
            );
            out.push(Family { name: format!("{}+{}", hc_name, lc_name), members });
        }
    }

    // Plot/report thresholds are applied by downstream writers only after this
    // final HC reassignment and LC clone-root pass has frozen family membership.
    out.sort_by_key(|f| std::cmp::Reverse(f.members.len()));
    out
}

pub fn analyze(
    vdj_dir: &Path,
    out_dir: &Path,
    min_structure_family: usize,
    threads: usize,
    min_clonomap_family: usize,
    min_clonomap_paired_family: usize,
    clonomap_k: usize,
    clonomap_radial: bool,
    max_cdr3_distance: usize,
) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("building Valkyrn Rayon pool")?;
    pool.install(|| {
        analyze_inner(
            vdj_dir,
            out_dir,
            min_structure_family,
            min_clonomap_family,
            min_clonomap_paired_family,
            clonomap_k,
            clonomap_radial,
            max_cdr3_distance,
        )
    })
}

fn analyze_inner(
    vdj_dir: &Path,
    out_dir: &Path,
    min_structure_family: usize,
    min_clonomap_family: usize,
    min_clonomap_paired_family: usize,
    clonomap_k: usize,
    clonomap_radial: bool,
    max_cdr3_distance: usize,
) -> Result<()> {
    // Validate the input before touching output, then make every Valkyrn run
    // self-contained. Reusing an output directory otherwise leaves stale plots
    // behind and makes it impossible to tell which renderer produced them.
    eprintln!("Valkyrn: loading VDJ calls from {}", vdj_dir.display());
    let calls = read_calls(vdj_dir)?;
    if calls.is_empty() {
        bail!("no VDJ calls found")
    }
    eprintln!("Valkyrn: loaded {} receptor calls", calls.len());
    let (calls, rejected_cells, cells_examined) = qualify_cells_by_hc(calls);
    let no_hc = rejected_cells.iter().filter(|r| r.reason == "no_valid_hc").count();
    let multiple_hc = rejected_cells.iter().filter(|r| r.reason == "multiple_valid_hc").count();
    let accepted_cells = cells_examined - rejected_cells.len();
    eprintln!("Valkyrn: HC validation");
    eprintln!("  cells examined:                 {cells_examined}");
    eprintln!("  accepted: exactly one valid HC: {accepted_cells}");
    eprintln!("  skipped: no valid HC:           {no_hc}");
    eprintln!("  skipped: multiple valid HCs:    {multiple_hc}");
    if calls.is_empty() {
        bail!("no cells with exactly one valid productive IGH/CDR3 reconstruction")
    }
    if out_dir.exists() {
        let vdj_canonical = fs::canonicalize(vdj_dir)
            .context("canonicalizing Valkyrn VDJ input directory")?;
        let out_canonical = fs::canonicalize(out_dir)
            .context("canonicalizing existing Valkyrn output directory")?;
        if out_canonical == vdj_canonical {
            bail!("refusing to purge Valkyrn output because --out is the VDJ input directory")
        }
        eprintln!("Valkyrn: purging previous output {}", out_dir.display());
        let purge_started = Instant::now();
        fs::remove_dir_all(out_dir).context("purging previous Valkyrn output directory")?;
        eprintln!(
            "Valkyrn: previous output purged in {:.2?}",
            purge_started.elapsed()
        );
    }
    fs::create_dir_all(out_dir)?;
    write_hc_rejections(out_dir, &rejected_cells)?;
    eprintln!("Valkyrn: wrote {} rejected cells to {}", rejected_cells.len(), out_dir.join("valkyrn_rejected_cells.tsv").display());
    eprintln!("Valkyrn: output directory ready; clustering receptor families from {accepted_cells} HC-qualified cells (max CDR3 edit distance = {})", max_cdr3_distance);
    let family_started = Instant::now();
    let fams = families(&calls, max_cdr3_distance);
    eprintln!(
        "Valkyrn: clustered {} receptor families in {:.2?}; computing mutation cache",
        fams.len(),
        family_started.elapsed()
    );
    // Needleman-Wunsch is by far the expensive operation. Compute every
    // rearrangement exactly once, in parallel, then reuse the immutable cache
    // for family summaries, recurrence classification and TSV output.
    let mutation_started = Instant::now();
    let mutation_cache: Vec<Option<Vec<Mutation>>> = calls
        .par_iter()
        .map(|c| {
            measure_mutations_nw(&c.naive, &c.observed)
                .map(|x| x.mutations)
        })
        .collect();
    let rejected_alignments = mutation_cache.iter().filter(|x| x.is_none()).count();
    eprintln!(
        "Valkyrn: mutation cache computed for {} receptor calls in {:.2?} ({} rejected: empty receptor sequence); writing summary tables",
        calls.len(),
        mutation_started.elapsed(),
        rejected_alignments
    );
    write_families(out_dir, &calls, &fams, &mutation_cache)?;
    write_receptors(out_dir, &calls, &fams)?;
    write_mutations(out_dir, &calls, &fams, &mutation_cache)?;
    write_indel_events(out_dir, &calls)?;
    write_cell_qc(out_dir, &calls, &mutation_cache)?;
    write_structure_candidates(out_dir, &calls, &fams, min_structure_family)?;
    write_clonomap(
        out_dir,
        &calls,
        &fams,
        min_clonomap_family,
        min_clonomap_paired_family,
        clonomap_k,
        clonomap_radial,
    )?;
    eprintln!("Valkyrn: writing report");
    write_report(out_dir, &calls, &fams)?;
    eprintln!("Valkyrn: finished -> {}", out_dir.display());
    Ok(())
}

fn pn_reference_index(calls: &[Call], members: &[usize]) -> Option<usize> {
    // Family membership is already fixed by V/J + the hard CDR3 rule.  P/N and
    // retained-D measurements are used only to choose a descriptive ancestral
    // reconstruction for downstream distance plots; they cannot merge families.
    // Ties prefer more explicit P sequence, then stronger receptor evidence.
    members.iter().copied().max_by_key(|&i| {
        (
            calls[i].d_retained_len,
            calls[i].p_total_len,
            calls[i].support.saturating_add(calls[i].rediscovery),
        )
    })
}

fn family_pn_reference<'a>(calls: &'a [Call], f: &Family) -> Option<&'a str> {
    pn_reference_index(calls, &f.members).map(|i| calls[i].naive.as_str())
}

fn pn_distance_for_call(c: &Call, reference: Option<&str>) -> Option<usize> {
    let reference = reference?;
    if reference.is_empty() || c.naive.is_empty() {
        None
    } else {
        Some(edit_distance(reference, &c.naive))
    }
}

fn write_receptors(out: &Path, calls: &[Call], fams: &[Family]) -> Result<()> {
    let mut w = writer(out.join("valkyrn_receptors.tsv"))?;
    writeln!(
        w,
        "family\trecombination_id\tcell\tchain\tcdr3_nt\tcdr3_aa\tpn_alternative\tpn_reconstruction_distance_nt"
    )?;
    for f in fams {
        let reference = family_pn_reference(calls, f);
        for &i in &f.members {
            let c = &calls[i];
            let pn_distance = pn_distance_for_call(c, reference)
                .map(|x| x.to_string())
                .unwrap_or_default();
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                f.name,
                c.id,
                c.cell,
                c.chain,
                c.cdr3_nt,
                c.cdr3_aa,
                c.pn_alternative,
                pn_distance
            )?;
        }
    }
    Ok(())
}

fn write_families(
    out: &Path,
    calls: &[Call],
    fams: &[Family],
    mutation_cache: &[Option<Vec<Mutation>>],
) -> Result<()> {
    let mut w = writer(out.join("valkyrn_families.tsv"))?;
    writeln!(
        w,
        "family\tchain\tcells\tv\td\tj\tstructural_recombination_ids\tcdr3_nt_variants\tcdr3_aa_variants\tpn_alternative_members\tpn_distance_min_nt\tpn_distance_median_nt\tpn_distance_max_nt\tisotypes\tisotype_counts\tlight_partners\tdominant_light_fraction\tshared_mutations\tvariable_mutations"
    )?;
    for f in fams {
        let cells: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        let c = &calls[f.members[0]];
        let cdr_nt: BTreeSet<_> = f
            .members
            .iter()
            .filter_map(|&i| (!calls[i].cdr3_nt.is_empty()).then_some(calls[i].cdr3_nt.as_str()))
            .collect();
        let cdr_aa: BTreeSet<_> = f
            .members
            .iter()
            .filter_map(|&i| (!calls[i].cdr3_aa.is_empty()).then_some(calls[i].cdr3_aa.as_str()))
            .collect();
        let pn_reference = family_pn_reference(calls, f);
        let pn_distances: Vec<usize> = f
            .members
            .iter()
            .filter_map(|&i| pn_distance_for_call(&calls[i], pn_reference))
            .collect();
        let pn_alt_members = f.members.iter().filter(|&&i| calls[i].pn_alternative).count();
        let pn_min = pn_distances.iter().min().map(|x| x.to_string()).unwrap_or_default();
        let pn_med = median_usize(pn_distances.clone())
            .map(|x| format!("{x:.1}"))
            .unwrap_or_default();
        let pn_max = pn_distances.iter().max().map(|x| x.to_string()).unwrap_or_default();
        let structural: BTreeSet<_> = f.members.iter().map(|&i| calls[i].id.as_str()).collect();
        let lights = if c.chain == "IGH" {
            light_partners(calls, &cells)
        } else {
            BTreeMap::new()
        };
        let mut isotypes: BTreeMap<&str, usize> = BTreeMap::new();
        for &i in &f.members {
            if !calls[i].c.is_empty() {
                *isotypes.entry(calls[i].c.as_str()).or_default() += 1;
            }
        }
        let isotype_names = isotypes.keys().copied().collect::<Vec<_>>().join(",");
        let isotype_counts = isotypes
            .iter()
            .map(|(k, v)| format!("{k}:{v}"))
            .collect::<Vec<_>>()
            .join(",");
        let light_total: usize = lights.values().sum();
        let dominant = lights.values().max().copied().unwrap_or(0);
        let (shared, var) = family_mutation_sets(mutation_cache, f);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            f.name,
            c.chain,
            cells.len(),
            c.v,
            c.d,
            c.j,
            structural.len(),
            cdr_nt.into_iter().collect::<Vec<_>>().join(","),
            cdr_aa.into_iter().collect::<Vec<_>>().join(","),
            pn_alt_members,
            pn_min,
            pn_med,
            pn_max,
            isotype_names,
            isotype_counts,
            lights.len(),
            if light_total > 0 {
                dominant as f64 / light_total as f64
            } else {
                0.0
            },
            shared.len(),
            var.len()
        )?;
    }
    Ok(())
}
fn light_partners(calls: &[Call], cells: &BTreeSet<&str>) -> BTreeMap<String, usize> {
    let indices: Vec<usize> = calls.iter().enumerate().filter(|(_, c)| {
        cells.contains(c.cell.as_str())
            && c.productive()
            && matches!(c.chain.as_str(), "IGK" | "IGL")
    }).map(|(i, _)| i).collect();
    let canonical = canonical_light_ids(calls, &indices);
    let mut x = BTreeMap::new();
    for i in indices {
        let id = canonical.get(&i).cloned().unwrap_or_else(|| calls[i].id.clone());
        *x.entry(id).or_default() += 1;
    }
    x
}

fn family_mutation_sets(
    mutation_cache: &[Option<Vec<Mutation>>],
    f: &Family,
) -> (BTreeSet<Mutation>, BTreeSet<Mutation>) {
    let sets: Vec<BTreeSet<Mutation>> = f
        .members
        .iter()
        .map(|&i| mutation_cache[i].as_deref().unwrap_or(&[]).iter().cloned().collect())
        .collect();
    if sets.is_empty() {
        return (Default::default(), Default::default());
    }
    let mut shared = sets[0].clone();
    let mut union = sets[0].clone();
    for s in &sets[1..] {
        shared = shared.intersection(s).cloned().collect();
        union.extend(s.iter().cloned());
    }
    let variable = union.difference(&shared).cloned().collect();
    (shared, variable)
}

fn write_mutations(
    out: &Path,
    calls: &[Call],
    fams: &[Family],
    mutation_cache: &[Option<Vec<Mutation>>],
) -> Result<()> {
    let mut w = writer(out.join("valkyrn_mutations.tsv"))?;
    writeln!(
        w,
        "family\trecombination_id\tcell\tchain\tv\tj\tposition\tnaive_base\tobserved_base\tclass\tindependent_families_same_v"
    )?;

    // Count the same aligned difference across independent receptor families using
    // the same chain and V call. This is deliberately only a recurrence flag:
    // without segment-aware germline coordinates it is not enough to call an
    // unrepresented germline allele.
    let mut recurrence: HashMap<(String, String, Mutation), BTreeSet<String>> = HashMap::new();
    for f in fams {
        let mut seen = BTreeSet::new();
        for &i in &f.members {
            let c = &calls[i];
            for m in mutation_cache[i].as_deref().unwrap_or(&[]).iter().cloned() {
                seen.insert((c.chain.clone(), c.v.clone(), m));
            }
        }
        for (chain, v, m) in seen {
            recurrence
                .entry((chain, v, m))
                .or_default()
                .insert(f.name.clone());
        }
    }

    for f in fams {
        let (shared, _) = family_mutation_sets(mutation_cache, f);
        let cells: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        for &i in &f.members {
            let c = &calls[i];
            for m in mutation_cache[i].as_deref().unwrap_or(&[]).iter().cloned() {
                let n = recurrence
                    .get(&(c.chain.clone(), c.v.clone(), m.clone()))
                    .map_or(1, |x| x.len());
                let class = if cells.len() == 1 {
                    "singleton_observed_difference"
                } else if n > 1 {
                    "recurrent_same_v_candidate"
                } else if shared.contains(&m) {
                    "family_shared_candidate"
                } else {
                    "branch_or_private_candidate"
                };
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    f.name,
                    c.id,
                    c.cell,
                    c.chain,
                    c.v,
                    c.j,
                    m.pos + 1,
                    m.from as char,
                    m.to as char,
                    class,
                    n
                )?;
            }
        }
    }
    Ok(())
}

fn median_usize(mut values: Vec<usize>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    let n = values.len();
    Some(if n % 2 == 1 {
        values[n / 2] as f64
    } else {
        (values[n / 2 - 1] + values[n / 2]) as f64 / 2.0
    })
}

fn write_indel_events(out: &Path, calls: &[Call]) -> Result<()> {
    let mut w = writer(out.join("valkyrn_indels.tsv"))?;
    writeln!(w, "cell\treceptor_id\tchain\tsubstitutions\tindel_events\tinserted_nt\tdeleted_nt\tindels")?;
    for c in calls {
        let Some(aln) = measure_mutations_nw(&c.naive, &c.observed) else { continue };
        if aln.indels.is_empty() {
            continue;
        }
        let events = aln.indels.iter().map(|x| {
            let kind = match x.kind { IndelKind::Insertion => "ins", IndelKind::Deletion => "del" };
            format!("{kind}@{}:{}nt", x.pos, x.len)
        }).collect::<Vec<_>>().join(",");
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            c.cell, c.id, c.chain, aln.mutations.len(), aln.indels.len(),
            aln.inserted_nt(), aln.deleted_nt(), events)?;
    }
    Ok(())
}

fn write_cell_qc(out: &Path, calls: &[Call], mutation_cache: &[Option<Vec<Mutation>>]) -> Result<()> {
    let mut by: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, c) in calls.iter().enumerate() {
        by.entry(&c.cell).or_default().push(i)
    }
    let mut w = writer(out.join("valkyrn_cells.tsv"))?;
    writeln!(
        w,
        "cell\tproductive_heavy\tproductive_light\ttotal_heavy\ttotal_light\tpaired_productive\theavy_ids\tlight_ids\tpaired_ids\thc_depth_median_nt\thc_depth_min_nt\thc_depth_max_nt\tlc_depth_median_nt\tlc_depth_min_nt\tlc_depth_max_nt\tpaired_depth_median_nt"
    )?;
    for (cell, idxs) in by {
        let cs: Vec<&Call> = idxs.iter().map(|&i| &calls[i]).collect();
        let ph = cs.iter().filter(|c| c.productive() && c.heavy()).count();
        let pl = cs.iter().filter(|c| c.productive() && !c.heavy()).count();
        let th = cs.iter().filter(|c| c.heavy()).count();
        let tl = cs.len() - th;
        let heavy_ids: BTreeSet<&str> = cs
            .iter()
            .filter(|c| c.productive() && c.chain == "IGH" && c.id.starts_with("HC:"))
            .map(|c| c.id.as_str())
            .collect();
        let light_ids: BTreeSet<&str> = cs
            .iter()
            .filter(|c| {
                c.productive()
                    && matches!(c.chain.as_str(), "IGK" | "IGL")
                    && c.id.starts_with("LC:")
            })
            .map(|c| c.id.as_str())
            .collect();
        let paired_ids: Vec<String> = heavy_ids
            .iter()
            .flat_map(|h| light_ids.iter().map(move |l| format!("{h}+{l}")))
            .collect();

        let hc_depths: Vec<usize> = idxs
            .iter()
            .copied()
            .filter(|&i| calls[i].productive() && calls[i].chain == "IGH")
            .filter_map(|i| mutation_cache[i].as_ref().map(Vec::len))
            .collect();
        let lc_depths: Vec<usize> = idxs
            .iter()
            .copied()
            .filter(|&i| {
                calls[i].productive() && matches!(calls[i].chain.as_str(), "IGK" | "IGL")
            })
            .filter_map(|i| mutation_cache[i].as_ref().map(Vec::len))
            .collect();
        let hc_med = median_usize(hc_depths.clone());
        let lc_med = median_usize(lc_depths.clone());
        let paired_med = match (hc_med, lc_med) {
            (Some(h), Some(l)) => Some(h + l),
            _ => None,
        };
        let fmt = |x: Option<f64>| x.map(|v| format!("{v:.1}")).unwrap_or_default();
        let minv = |x: &[usize]| x.iter().min().map(|v| v.to_string()).unwrap_or_default();
        let maxv = |x: &[usize]| x.iter().max().map(|v| v.to_string()).unwrap_or_default();

        writeln!(
            w,
            "{cell}\t{ph}\t{pl}\t{th}\t{tl}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            ph > 0 && pl > 0,
            heavy_ids.into_iter().collect::<Vec<_>>().join(","),
            light_ids.into_iter().collect::<Vec<_>>().join(","),
            paired_ids.join(","),
            fmt(hc_med), minv(&hc_depths), maxv(&hc_depths),
            fmt(lc_med), minv(&lc_depths), maxv(&lc_depths), fmt(paired_med)
        )?;
    }
    Ok(())
}

fn write_structure_candidates(
    out: &Path,
    calls: &[Call],
    fams: &[Family],
    min_size: usize,
) -> Result<()> {
    let dir = out.join("structure_candidates");
    fs::create_dir_all(&dir)?;
    let mut m = writer(dir.join("manifest.tsv"))?;
    let mut fa = writer(dir.join("paired_receptors.fasta"))?;
    writeln!(
        m,
        "candidate\tfamily\tcell\theavy_id\tlight_id\theavy_v\theavy_j\tlight_chain\tlight_v\tlight_j\treason"
    )?;
    for f in fams
        .iter()
        .filter(|f| f.members.len() >= min_size && calls[f.members[0]].chain == "IGH")
    {
        let cells: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        let partners = light_partners(calls, &cells);
        if partners.len() < 2 {
            continue;
        }
        for &hi in &f.members {
            let h = &calls[hi];
            for l in calls.iter().filter(|c| {
                c.cell == h.cell && c.productive() && matches!(c.chain.as_str(), "IGK" | "IGL")
            }) {
                let id = format!("{}__{}__{}", f.name, h.id, l.id).replace(':', "_");
                writeln!(
                    m,
                    "{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\texpanded_HC_family_with_diverse_productive_LC",
                    f.name, h.cell, h.id, l.id, h.v, h.j, l.chain, l.v, l.j
                )?;
                writeln!(fa, ">{id}|H\n{}\n>{id}|L\n{}", h.observed, l.observed)?;
            }
        }
    }
    Ok(())
}

fn observed_on_naive_coordinates(naive: &str, observed: &str) -> Option<String> {
    let aln = measure_mutations_nw(naive, observed)?;
    let b = observed.as_bytes();
    let mut out = vec![b'-'; naive.len()];
    for (ai, bj) in aln.aligned {
        if let (Some(ai), Some(bj)) = (ai, bj) {
            out[ai] = b[bj].to_ascii_uppercase();
        }
    }
    String::from_utf8(out).ok()
}

fn xml_escape(x: &str) -> String {
    x.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// Cheap whole-HC overview of productive structural light-chain partners.
///
/// This is deliberately NOT a lineage tree: every LC sector is an independent
/// child of the HC family and no LC is ever connected to another LC. Sector
/// angle is proportional to the number of distinct cells carrying that LC;
/// the outer dot radius also scales with cell count. This remains readable for
/// tens to hundreds of LC partners and requires no PCA, distance matrix or MST.
fn write_lc_constellation(
    path: &Path,
    family: &Family,
    calls: &[Call],
    lights_by_cell: &HashMap<String, Vec<usize>>,
) -> Result<usize> {
    #[derive(Default)]
    struct LcMeta {
        cells: BTreeSet<String>,
        depths: Vec<usize>,
        chain: String,
        v: String,
        j: String,
    }
    let hc_cells: BTreeSet<&str> = family
        .members
        .iter()
        .map(|&i| calls[i].cell.as_str())
        .collect();
    let mut by_lc: BTreeMap<String, LcMeta> = BTreeMap::new();
    for cell in &hc_cells {
        if let Some(ls) = lights_by_cell.get(*cell) {
            for &li in ls {
                let lc = &calls[li];
                let m = by_lc.entry(lc.id.clone()).or_default();
                m.cells.insert((*cell).to_string());
                m.depths.push(measure_mutations_nw(&lc.naive, &lc.observed).map(|x| x.event_count()).unwrap_or(0));
                m.chain = lc.chain.clone();
                m.v = lc.v.clone();
                m.j = lc.j.clone();
            }
        }
    }
    if by_lc.is_empty() {
        return Ok(0);
    }
    let mut items: Vec<_> = by_lc.into_iter().collect();
    items.sort_by(|a, b| {
        b.1.cells
            .len()
            .cmp(&a.1.cells.len())
            .then_with(|| a.0.cmp(&b.0))
    });
    let total: usize = items.iter().map(|(_, m)| m.cells.len()).sum();
    if total == 0 {
        return Ok(0);
    }
    let anchor = &calls[family.members[0]];
    let (w, h) = (1000.0f64, 1000.0f64);
    let (cx, cy) = (500.0, 500.0);
    let inner = 145.0;
    let outer = 365.0;
    let palette = [
        "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
        "#9C755F", "#BAB0AC",
    ];
    let mut svg = writer(path)?;
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img">"#
    )?;
    writeln!(svg, r#"<rect width="100%" height="100%" fill="white"/>"#)?;
    writeln!(
        svg,
        r#"<text x="500" y="38" text-anchor="middle" font-family="sans-serif" font-size="22" font-weight="bold">HC → productive LC repertoire</text>"#
    )?;
    writeln!(
        svg,
        r##"<text x="500" y="66" text-anchor="middle" font-family="sans-serif" font-size="14" fill="#444">{} | {} / {} / {} | {} HC cells | {} structural LCs</text>"##,
        xml_escape(&family.name),
        xml_escape(&anchor.v),
        xml_escape(&anchor.d),
        xml_escape(&anchor.j),
        hc_cells.len(),
        items.len()
    )?;
    let polar = |r: f64, a: f64| (cx + r * a.cos(), cy + r * a.sin());
    let mut angle = -std::f64::consts::FRAC_PI_2;
    for (rank, (id, m)) in items.iter().enumerate() {
        let n = m.cells.len();
        let frac = n as f64 / total as f64;
        let span = frac * std::f64::consts::TAU;
        let gap = (0.012f64).min(span * 0.12);
        let a0 = angle + gap / 2.0;
        let a1 = angle + span - gap / 2.0;
        let amid = angle + span / 2.0;
        angle += span;
        let (x0, y0) = polar(inner, a0);
        let (x1, y1) = polar(outer, a0);
        let (x2, y2) = polar(outer, a1);
        let (x3, y3) = polar(inner, a1);
        let large = if a1 - a0 > std::f64::consts::PI { 1 } else { 0 };
        let color = palette[rank % palette.len()];
        let mut ds = m.depths.clone();
        ds.sort_unstable();
        let med = if ds.is_empty() {
            None
        } else if ds.len() % 2 == 1 {
            Some(ds[ds.len() / 2] as f64)
        } else {
            Some((ds[ds.len() / 2 - 1] + ds[ds.len() / 2]) as f64 / 2.0)
        };
        let tip_r = (4.0 + (n as f64).sqrt() * 1.15).min(28.0);
        let (tx, ty) = polar(outer + 18.0, amid);
        writeln!(
            svg,
            r##"<g><title>{} | {} {} / {} | {} cells ({:.1}%) | median LC depth: {}</title><path d="M {:.2} {:.2} L {:.2} {:.2} A {:.2} {:.2} 0 {} 1 {:.2} {:.2} L {:.2} {:.2} A {:.2} {:.2} 0 {} 0 {:.2} {:.2} Z" fill="{}" fill-opacity="0.78" stroke="white" stroke-width="1"/><circle cx="{:.2}" cy="{:.2}" r="{:.2}" fill="{}" stroke="#222" stroke-width="0.7"/></g>"##,
            xml_escape(id),
            xml_escape(&m.chain),
            xml_escape(&m.v),
            xml_escape(&m.j),
            n,
            100.0 * frac,
            med.map(|x| format!("{x:.1} nt"))
                .unwrap_or_else(|| "NA".into()),
            x0,
            y0,
            x1,
            y1,
            outer,
            outer,
            large,
            x2,
            y2,
            x3,
            y3,
            inner,
            inner,
            large,
            x0,
            y0,
            color,
            tx,
            ty,
            tip_r,
            color
        )?;
        if rank < 12 || frac >= 0.025 {
            let (lx, ly) = polar(outer + 52.0, amid);
            let anchor_txt = if lx < cx { "end" } else { "start" };
            writeln!(
                svg,
                r#"<text x="{:.2}" y="{:.2}" text-anchor="{}" dominant-baseline="middle" font-family="sans-serif" font-size="11">LC{} · {} cells</text>"#,
                lx,
                ly,
                anchor_txt,
                rank + 1,
                n
            )?;
        }
    }
    writeln!(
        svg,
        r##"<circle cx="500" cy="500" r="128" fill="#f7f7f7" stroke="#222" stroke-width="1.2"/>"##
    )?;
    writeln!(
        svg,
        r#"<text x="500" y="477" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold">{}</text>"#,
        xml_escape(&family.name)
    )?;
    writeln!(
        svg,
        r#"<text x="500" y="502" text-anchor="middle" font-family="sans-serif" font-size="14">{} HC cells</text>"#,
        hc_cells.len()
    )?;
    writeln!(
        svg,
        r#"<text x="500" y="525" text-anchor="middle" font-family="sans-serif" font-size="14">{} productive LC rearrangements</text>"#,
        items.len()
    )?;
    writeln!(
        svg,
        r##"<text x="500" y="548" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#555">sectors are independent LC partners — not LC→LC transitions</text>"##
    )?;
    writeln!(
        svg,
        r##"<text x="500" y="965" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#555">sector angle and outer-dot size encode distinct-cell abundance; hover a sector in a browser for LC identity and mutational depth</text>"##
    )?;
    writeln!(svg, "</svg>")?;
    Ok(items.len())
}

fn write_clonomap(
    out: &Path,
    calls: &[Call],
    fams: &[Family],
    min_size: usize,
    min_paired_size: usize,
    k: usize,
    radial_layout: bool,
) -> Result<()> {
    use clonomap::ClonoMap;
    #[derive(Default)]
    struct StateMeta {
        cells: BTreeSet<String>,
        isotypes: BTreeMap<String, usize>,
        lights: BTreeMap<String, usize>,
        hc_depths: Vec<usize>,
        lc_depths: Vec<usize>,
        pn_distances: Vec<usize>,
        original_hc_ids: BTreeSet<String>,
        original_lc_ids: BTreeSet<String>,
    }
    fn fmt_counts(x: &BTreeMap<String, usize>) -> String {
        x.iter()
            .map(|(k, v)| format!("{k}:{v}"))
            .collect::<Vec<_>>()
            .join(",")
    }
    fn dominant(x: &BTreeMap<String, usize>, empty: &str) -> String {
        x.iter()
            .max_by_key(|(_, n)| *n)
            .map(|(k, _)| k.clone())
            .unwrap_or_else(|| empty.to_string())
    }
    fn median(xs: &[usize]) -> Option<f32> {
        if xs.is_empty() {
            return None;
        }
        let mut x = xs.to_vec();
        x.sort_unstable();
        let n = x.len();
        Some(if n % 2 == 1 {
            x[n / 2] as f32
        } else {
            (x[n / 2 - 1] + x[n / 2]) as f32 / 2.0
        })
    }
    fn safe_name(x: &str) -> String {
        x.replace(':', "_").replace('~', "_")
    }

    let mut lights_by_cell: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, c) in calls
        .iter()
        .enumerate()
        .filter(|(_, c)| c.productive() && matches!(c.chain.as_str(), "IGK" | "IGL"))
    {
        lights_by_cell.entry(c.cell.clone()).or_default().push(i);
    }

    fn render_group(
        fdir: &Path,
        label: &str,
        member_indices: &[usize],
        calls: &[Call],
        lights_by_cell: &HashMap<String, Vec<usize>>,
        required_lc: Option<&str>,
        k: usize,
        radial_layout: bool,
    ) -> Result<usize> {
        use ndarray::Array2;

        let Some(&anchor_i) = member_indices.first() else {
            return Ok(0);
        };
        // Plotting must not require every reconstructed HC naive sequence to have
        // exactly the same length. Junction reconstruction can legitimately vary
        // within an already-frozen HC family. Pad mutation feature vectors to the
        // widest member instead of dropping those cells.
        let hc_width = member_indices
            .iter()
            .map(|&i| calls[i].naive.len())
            .max()
            .unwrap_or_else(|| calls[anchor_i].naive.len());
        let member_cells: BTreeSet<&str> = member_indices.iter().map(|&i| calls[i].cell.as_str()).collect();
        let light_indices: Vec<usize> = lights_by_cell.iter()
            .filter(|(cell, _)| member_cells.contains(cell.as_str()))
            .flat_map(|(_, xs)| xs.iter().copied())
            .collect();
        let canonical_lc = canonical_light_ids(calls, &light_indices);
        // Use the compatible reconstruction retaining the most explicit P bases
        // as the P/N reference. Ties prefer stronger receptor support. Distances
        // to this sequence are descriptive; they are never used as a merge cutoff.
        let pn_reference_i = pn_reference_index(calls, member_indices).unwrap_or(anchor_i);
        let pn_reference = calls[pn_reference_i].naive.as_str();

        // Build every HC coordinate against the one family-level HC NAIVE
        // reconstruction.  A family must have one coordinate origin; allowing
        // every member to use its own reconstructed naive sequence hides exactly
        // the P/N reconstruction differences that the PN heat layer reports.
        //
        // LC keeps the original two signals -- mutation state relative to each
        // LC rearrangement's own reconstructed naive plus the weighted structural
        // LC identity -- and ADDS the family-level artificial LC CDR3 coordinate.
        // The latter supplies a shared cross-LC sequence reference without throwing
        // away the local SHM or structural-recombination information.
        #[derive(Clone)]
        struct PairedState {
            label: String,
            hc_seq: String,
            lc_id: String,
            hc_mut: Vec<f32>,
            lc_mut: Vec<f32>,
        }

        let mut lc_ids: BTreeSet<String> = BTreeSet::new();
        let mut lc_local_width = 0usize;
        let mut lc_cdr3s: Vec<&str> = Vec::new();
        for &i in member_indices {
            if let Some(ls) = lights_by_cell.get(&calls[i].cell) {
                for &li in ls {
                    let lc = &calls[li];
                    let lc_id = canonical_lc.get(&li).map(String::as_str).unwrap_or(lc.id.as_str());
                    if required_lc.is_some_and(|want| lc_id != want) {
                        continue;
                    }
                    lc_ids.insert(lc_id.to_string());
                    lc_local_width = lc_local_width.max(lc.naive.len());
                    if !lc.cdr3_nt.is_empty() {
                        lc_cdr3s.push(lc.cdr3_nt.as_str());
                    }
                }
            }
        }
        if lc_ids.is_empty() {
            lc_ids.insert("unpaired".to_string());
        }
        let lc_common_width = lc_cdr3s.iter().map(|x| x.len()).max().unwrap_or(0);
        let mut lc_counts = vec![[0usize; 4]; lc_common_width];
        for seq in &lc_cdr3s {
            for (pos, base) in seq.bytes().enumerate() {
                let slot = match base.to_ascii_uppercase() {
                    b'A' => Some(0),
                    b'C' => Some(1),
                    b'G' => Some(2),
                    b'T' => Some(3),
                    _ => None,
                };
                if let Some(slot) = slot {
                    lc_counts[pos][slot] += 1;
                }
            }
        }
        let lc_consensus: Vec<u8> = lc_counts
            .iter()
            .map(|counts| {
                let max = counts.iter().copied().max().unwrap_or(0);
                if max == 0 || counts.iter().filter(|&&n| n == max).count() != 1 {
                    b'N'
                } else {
                    b"ACGT"[counts.iter().position(|&n| n == max).unwrap()]
                }
            })
            .collect();

        let lc_consensus_vector = |observed: &str| -> Vec<f32> {
            let obs = observed.as_bytes();
            (0..lc_common_width)
                .map(|pos| {
                    let Some(&reference) = lc_consensus.get(pos) else { return 0.0 };
                    if reference == b'N' {
                        0.0
                    } else {
                        match obs.get(pos).map(|b| b.to_ascii_uppercase()) {
                            Some(base) if base == reference => 0.0,
                            // A different base or a shorter CDR3 is a difference
                            // from the artificial family LC reference.
                            _ => 1.0,
                        }
                    }
                })
                .collect()
        };

        // Build one representative CDR3 sequence per canonical LC identity and
        // order identities by sequence proximity.  This ordering is visualization
        // metadata only: it does not alter PCA, MST construction, or LC merging.
        // Start from the most abundant LC, then repeatedly take the nearest
        // unvisited LC by nucleotide edit distance.  Thus adjacent legend symbols
        // (and therefore visually similar split-circle codes) tend to denote close
        // LC CDR3 sequences.
        let mut lc_cdr3_counts: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
        for &li in &light_indices {
            let lc = &calls[li];
            if lc.cdr3_nt.is_empty() {
                continue;
            }
            let lc_id = canonical_lc.get(&li).map(String::as_str).unwrap_or(lc.id.as_str());
            *lc_cdr3_counts
                .entry(lc_id.to_string())
                .or_default()
                .entry(lc.cdr3_nt.clone())
                .or_default() += 1;
        }
        let mut lc_representative: BTreeMap<String, (String, usize)> = BTreeMap::new();
        for (lc_id, seqs) in lc_cdr3_counts {
            let total = seqs.values().sum();
            let representative = seqs
                .into_iter()
                .max_by(|(sa, na), (sb, nb)| na.cmp(nb).then_with(|| sb.cmp(sa)))
                .map(|(seq, _)| seq)
                .unwrap_or_default();
            lc_representative.insert(lc_id, (representative, total));
        }
        let mut lc_visual_order = Vec::with_capacity(lc_representative.len() + 1);
        if let Some((start, _)) = lc_representative
            .iter()
            .max_by(|(ida, (_, na)), (idb, (_, nb))| na.cmp(nb).then_with(|| idb.cmp(ida)))
        {
            let mut current = start.clone();
            lc_visual_order.push(current.clone());
            while lc_visual_order.len() < lc_representative.len() {
                let current_seq = &lc_representative[&current].0;
                let next = lc_representative
                    .iter()
                    .filter(|(id, _)| !lc_visual_order.contains(id))
                    .min_by(|(ida, (sa, na)), (idb, (sb, nb))| {
                        edit_distance(current_seq, sa)
                            .cmp(&edit_distance(current_seq, sb))
                            .then_with(|| nb.cmp(na))
                            .then_with(|| ida.cmp(idb))
                    })
                    .map(|(id, _)| id.clone());
                let Some(next) = next else { break };
                lc_visual_order.push(next.clone());
                current = next;
            }
        }
        if lc_ids.contains("unpaired") {
            lc_visual_order.push("unpaired".to_string());
        }

        let mutation_vector = |naive: &str, observed: &str, width: usize| -> Option<Vec<f32>> {
            let mut v = vec![0.0; width];
            for m in measure_mutations_nw(naive, observed)?.mutations {
                if m.pos < width {
                    v[m.pos] = 1.0;
                }
            }
            Some(v)
        };

        let mut states: BTreeMap<String, PairedState> = BTreeMap::new();
        let mut raw_meta: HashMap<String, StateMeta> = HashMap::new();
        for &i in member_indices {
            let hc = &calls[i];
            if hc.naive.is_empty() || hc.observed.is_empty() {
                continue;
            }
            // Family membership is already frozen. Mutation measurement is local
            // to this receptor: its own reconstructed naive -> observed sequence.
            // Failure here must never change family admission.
            let hc_seq = observed_on_naive_coordinates(&hc.naive, &hc.observed)
                .unwrap_or_else(|| hc.observed.clone());
            let hc_mut = mutation_vector(&hc.naive, &hc.observed, hc_width)
                .unwrap_or_else(|| vec![0.0; hc_width]);
            let mut paired_any = false;

            if let Some(ls) = lights_by_cell.get(&hc.cell) {
                for &li in ls {
                    let lc = &calls[li];
                    let lc_id = canonical_lc.get(&li).map(String::as_str).unwrap_or(lc.id.as_str());
                    if required_lc.is_some_and(|want| lc_id != want) {
                        continue;
                    }
                    if lc.observed.is_empty() {
                        continue;
                    }
                    if lc.cdr3_nt.is_empty() {
                        continue;
                    }
                    paired_any = true;
                    let Some(lc_mut) = mutation_vector(&lc.naive, &lc.observed, lc_local_width) else { continue };
                    let lc_common_mut = lc_consensus_vector(&lc.cdr3_nt);
                    let lc_pattern = lc_mut
                        .iter()
                        .map(|x| if *x > 0.0 { '1' } else { '0' })
                        .collect::<String>();
                    let lc_common_pattern = lc_common_mut
                        .iter()
                        .map(|x| if *x > 0.0 { '1' } else { '0' })
                        .collect::<String>();
                    let label = format!("{}|{}|{}|{}", hc_seq, lc_id, lc_pattern, lc_common_pattern);
                    states.entry(label.clone()).or_insert_with(|| PairedState {
                        label: label.clone(),
                        hc_seq: hc_seq.clone(),
                        lc_id: lc_id.to_string(),
                        hc_mut: hc_mut.clone(),
                        lc_mut: lc_mut.clone(),
                    });
                    let m = raw_meta.entry(label).or_default();
                    m.cells.insert(hc.cell.clone());
                    m.original_hc_ids.insert(hc.id.clone());
                    m.original_lc_ids.insert(lc.id.clone());
                    if let Some(aln) = measure_mutations_nw(&hc.naive, &hc.observed) {
                        m.hc_depths.push(aln.event_count());
                    }
                    m.pn_distances.push(edit_distance(pn_reference, &hc.naive));
                    let iso = if hc.c.is_empty() {
                        "unknown"
                    } else {
                        hc.c.as_str()
                    };
                    *m.isotypes.entry(iso.to_string()).or_default() += 1;
                    *m.lights.entry(lc_id.to_string()).or_default() += 1;
                    if let Some(aln) = measure_mutations_nw(&lc.naive, &lc.observed) {
                        m.lc_depths.push(aln.event_count());
                    }
                }
            }

            // Preserve genuinely unpaired HC cells in whole-family maps. They get
            // no LC mutation features and their own categorical state.
            if !paired_any && required_lc.is_none() {
                let label = format!("{}|unpaired|", hc_seq);
                states.entry(label.clone()).or_insert_with(|| PairedState {
                    label: label.clone(),
                    hc_seq: hc_seq.clone(),
                    lc_id: "unpaired".to_string(),
                    hc_mut: hc_mut.clone(),
                    lc_mut: vec![0.0; lc_local_width],
                });
                let m = raw_meta.entry(label).or_default();
                m.cells.insert(hc.cell.clone());
                m.original_hc_ids.insert(hc.id.clone());
                if let Some(aln) = measure_mutations_nw(&hc.naive, &hc.observed) {
                    m.hc_depths.push(aln.event_count());
                }
                m.pn_distances.push(edit_distance(pn_reference, &hc.naive));
                let iso = if hc.c.is_empty() {
                    "unknown"
                } else {
                    hc.c.as_str()
                };
                *m.isotypes.entry(iso.to_string()).or_default() += 1;
                *m.lights.entry("unpaired".to_string()).or_default() += 1;
            }
        }
        if states.len() < 3 || hc_width < 2 {
            return Ok(0);
        }

        // Add the inferred unmutated heavy-chain state as a real ClonoMap row.
        // It has zero HC mutations, zero LC mutations, and deliberately no LC
        // one-hot identity: the root represents HC ancestry without inventing an
        // ancestral light-chain rearrangement.
        const HC_NAIVE_STATE: &str = "HC NAIVE";
        let mut rows: Vec<PairedState> = states.into_values().collect();
        rows.insert(0, PairedState {
            label: HC_NAIVE_STATE.to_string(),
            hc_seq: pn_reference.to_string(),
            lc_id: String::new(),
            hc_mut: vec![0.0; hc_width],
            lc_mut: vec![0.0; lc_local_width],
        });
        // Topology is inferred independently inside each paired LC lineage.
        // Unrelated LC rearrangements therefore never compete for an MST edge.
        // Within a lineage the metric contains only biologically commensurate
        // mutation coordinates: HC vs the family HC NAIVE plus LC vs that LC's
        // reconstructed naive. The artificial cross-LC CDR3 and categorical LC
        // identity remain state/visualization metadata and do not pull the tree.
        let feature_cols = hc_width + lc_local_width;
        let mut features = Array2::<f32>::zeros((rows.len(), feature_cols));
        let mut topology_groups = Vec::with_capacity(rows.len());
        for (r, state) in rows.iter().enumerate() {
            for (j, x) in state.hc_mut.iter().enumerate() {
                features[[r, j]] = *x;
            }
            for (j, x) in state.lc_mut.iter().enumerate() {
                features[[r, hc_width + j]] = *x;
            }
            topology_groups.push(if r == 0 {
                HC_NAIVE_STATE.to_string()
            } else {
                state.lc_id.clone()
            });
        }
        fs::create_dir_all(fdir)?;
        let n_seqs = rows.len();
        let model = ClonoMap::from_grouped_feature_matrix(
            rows.iter().map(|x| x.label.clone()).collect(),
            features,
            topology_groups,
            0,
            hc_width,
            k,
        )
        .map_err(|e| anyhow::anyhow!("ClonoMap failed for {label}: {e}"))?;
        model
            .pca
            .to_tsv(&model.encoder.sequences, fdir.join("coords.tsv"))
            .map_err(|e| anyhow::anyhow!(e.to_string()))?;
        model.tree.to_tsv(fdir.join("tree.tsv"))?;
        model
            .encoder
            .sequences
            .to_tsv(fdir.join("rows.tsv"))
            .map_err(|e| anyhow::anyhow!("failed to write ClonoMap rows: {e}"))?;

        // HC NAIVE is part of the PCA/MST geometry itself. Because its feature
        // vector is all zeroes, it is the unmutated HC state with an empty LC
        // state rather than an observed LC-bearing receptor chosen after the fact.
        // All observed HC rows are encoded against this same HC reference; all
        // observed LC rows retain their local-naive mutation coordinates and
        // structural identity, with the family artificial LC CDR3 model appended.
        let naive_root = model
            .encoder
            .sequences
            .dna
            .iter()
            .position(|x| x == HC_NAIVE_STATE)
            .ok_or_else(|| anyhow::anyhow!("HC NAIVE state missing from ClonoMap model"))?;

        let mut iso_cat = Vec::new();
        let mut iso_mixed = Vec::new();
        let mut lc_cat = Vec::new();
        let mut lc_mixed = Vec::new();
        let mut abundance = Vec::new();
        let mut hc_depth = Vec::new();
        let mut lc_depth = Vec::new();
        let mut paired_depth = Vec::new();
        let mut pn_distance = Vec::new();
        let mut nodes = writer(fdir.join("nodes.tsv"))?;
        writeln!(
            nodes,
            "node\thc_family\tcell_ids\toriginal_hc_ids\toriginal_lc_ids\tstate\thc_dna\tcell_count\tisotypes\tdominant_isotype\tlight_chains\tdominant_light_chain\thc_depth_nt\tlc_depth_median_nt\tlc_depth_min_nt\tlc_depth_max_nt\tpaired_depth_median_nt\tpn_reconstruction_distance_nt"
        )?;
        for (node, key) in model.encoder.sequences.dna.iter().enumerate() {
            if key == HC_NAIVE_STATE {
                iso_mixed.push(false);
                lc_mixed.push(false);
                iso_cat.push(HC_NAIVE_STATE.to_string());
                lc_cat.push(HC_NAIVE_STATE.to_string());
                abundance.push(1);
                hc_depth.push(None);
                lc_depth.push(None);
                paired_depth.push(None);
                pn_distance.push(Some(0.0));
                writeln!(
                    nodes,
                    "{}\t{}\t\t\t\t{}\t{}\t0\t\tHC NAIVE\t\tLC empty\t\t\t\t\t\t0",
                    node,
                    label,
                    key,
                    pn_reference
                )?;
                continue;
            }
            let m = raw_meta.get(key);
            let empty_iso = BTreeMap::new();
            let empty_lc = BTreeMap::new();
            let isos = m.map(|x| &x.isotypes).unwrap_or(&empty_iso);
            let lcs = m.map(|x| &x.lights).unwrap_or(&empty_lc);
            let ni = dominant(isos, "unknown");
            let nl = dominant(lcs, "unpaired");
            iso_mixed.push(isos.len() > 1);
            lc_mixed.push(lcs.len() > 1);
            iso_cat.push(ni.clone());
            lc_cat.push(nl.clone());
            let n = m.map(|x| x.cells.len()).unwrap_or(1);
            abundance.push(n);
            let hc_depths = m.map(|x| x.hc_depths.as_slice()).unwrap_or(&[]);
            let hd = median(hc_depths);
            hc_depth.push(hd);
            let pnd = median(m.map(|x| x.pn_distances.as_slice()).unwrap_or(&[]));
            pn_distance.push(pnd);
            let depths = m.map(|x| x.lc_depths.as_slice()).unwrap_or(&[]);
            let med = median(depths);
            lc_depth.push(med);
            paired_depth.push(match (hd, med) {
                (Some(h), Some(l)) => Some(h + l),
                (Some(h), None) => Some(h),
                _ => None,
            });
            let min = depths
                .iter()
                .min()
                .map(|x| x.to_string())
                .unwrap_or_default();
            let max = depths
                .iter()
                .max()
                .map(|x| x.to_string())
                .unwrap_or_default();
            let med_s = med.map(|x| format!("{x:.1}")).unwrap_or_default();
            let hd_s = hd.map(|x| format!("{x:.1}")).unwrap_or_default();
            let pair_s = match (hd, med) {
                (Some(h), Some(l)) => format!("{:.1}", h + l),
                (Some(h), None) => format!("{h:.1}"),
                _ => String::new(),
            };
            let pnd_s = pnd.map(|x| format!("{x:.1}")).unwrap_or_default();
            let hc_dna = rows.get(node).map(|x| x.hc_seq.as_str()).unwrap_or("");
            writeln!(
                nodes,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                node,
                label,
                m.map(|x| x.cells.iter().cloned().collect::<Vec<_>>().join(",")).unwrap_or_default(),
                m.map(|x| x.original_hc_ids.iter().cloned().collect::<Vec<_>>().join(",")).unwrap_or_default(),
                m.map(|x| x.original_lc_ids.iter().cloned().collect::<Vec<_>>().join(",")).unwrap_or_default(),
                key,
                hc_dna,
                n,
                fmt_counts(isos),
                ni,
                fmt_counts(lcs),
                nl,
                hd_s,
                med_s,
                min,
                max,
                pair_s,
                pnd_s
            )?;
        }
        let cells: BTreeSet<_> = member_indices
            .iter()
            .map(|&i| calls[i].cell.as_str())
            .collect();
        let title = format!(
            "{} | {} / {} / {} | {} cells | {} paired states",
            label,
            calls[anchor_i].v,
            calls[anchor_i].d,
            calls[anchor_i].j,
            cells.len(),
            model.encoder.sequences.len()
        );
        model
            .tree
            .plot_rooted_annotated_cached(
                model.coords().nrows(),
                naive_root,
                &iso_cat,
                None,
                &iso_mixed,
                &abundance,
                "IGH constant class",
                &title,
                fdir.join("mst_rooted_isotype.svg")
                    .to_string_lossy()
                    .as_ref(),
                radial_layout,
            )
            .map_err(|e| anyhow::anyhow!("failed to plot isotype-rooted MST: {e}"))?;
        model
            .tree
            .plot_rooted_annotated_cached(
                model.coords().nrows(),
                naive_root,
                &lc_cat,
                Some(&lc_visual_order),
                &lc_mixed,
                &abundance,
                "Structural light-chain identity",
                &title,
                fdir.join("mst_rooted_light_chain.svg")
                    .to_string_lossy()
                    .as_ref(),
                radial_layout,
            )
            .map_err(|e| anyhow::anyhow!("failed to plot light-chain-rooted MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                &hc_depth,
                &abundance,
                "HC mutational depth (nt from NAIVE)",
                &title,
                fdir.join("mst_rooted_hc_depth.svg")
                    .to_string_lossy()
                    .as_ref(),
                radial_layout,
            )
            .map_err(|e| anyhow::anyhow!("failed to plot HC-depth MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                &pn_distance,
                &abundance,
                "P/N reconstruction distance (nt)",
                &title,
                fdir.join("mst_rooted_pn_distance.svg")
                    .to_string_lossy()
                    .as_ref(),
                radial_layout,
            )
            .map_err(|e| anyhow::anyhow!("failed to plot P/N-distance MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                &lc_depth,
                &abundance,
                "Linked LC mutational depth (median nt)",
                &title,
                fdir.join("mst_rooted_lc_depth.svg")
                    .to_string_lossy()
                    .as_ref(),
                radial_layout,
            )
            .map_err(|e| anyhow::anyhow!("failed to plot LC-depth MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                &paired_depth,
                &abundance,
                "Paired HC+LC mutational depth (median nt)",
                &title,
                fdir.join("mst_rooted_paired_depth.svg")
                    .to_string_lossy()
                    .as_ref(),
                radial_layout,
            )
            .map_err(|e| anyhow::anyhow!("failed to plot paired-depth MST: {e}"))?;
        Ok(n_seqs)
    }

    let dir = out.join("clonomap");
    fs::create_dir_all(&dir)?;
    let mut readme = writer(dir.join("README.md"))?;
    write!(
        readme,
        r#"# Valkyrn / ClonoMap lineage maps

Each family is drawn on one cached ClonoMap minimum-spanning-tree topology built
from the **paired receptor state**: HC mutation coordinates + LC mutation
coordinates + a weighted one-hot structural LC identity. LC identity therefore
participates in the geometry instead of being painted onto an HC-only tree after
the fact. `HC NAIVE` is a synthetic **but real ClonoMap state** included in the
PCA and MST: its HC mutation coordinates are zero, its LC mutation coordinates
are zero, and its LC identity block is empty. It therefore represents the
reconstructed unmutated HC without assigning any observed LC as ancestral.
`HC NAIVE` is drawn grey, like a missing LC, and is not an observed cell.
Rooted ClonoMaps use a layered
left-to-right layout by default; `--clonomap-radial` restores the radial view. The topology is a paired sequence-state landscape, not a
claim of chronological phylogeny.

## Visual encoding

- **Node size = clone/state abundance**: the number of distinct cells represented
  by that observed HC sequence state. Every SVG includes reference circles drawn
  with the same radius transform, so circle size can be read back as cell count.
- **Edges = the same cached ClonoMap MST** in every view. Annotation never changes
  the topology.
- **Black node border** keeps adjacent states visually distinct. An additional
  outer black ring marks a mixed categorical state: cells with the same HC
  sequence state carry more than one isotype or productive LC identity.
- Categorical views use deterministic **slash-split two-colour circles**. The
  eight Okabe-Ito base colours are deliberately discrete (no near-shade variants);
  ordered pairs provide 56 categorical identities without perceptual shade copies.
- `mst_rooted_isotype.svg`: node fill is IGH constant-region class/isotype.
- `lc_repertoire.svg` (whole-HC families): a cheap radial repertoire overview.
  Every structural LC is an independent sector attached to the HC family; there
  are deliberately no LC-to-LC edges. Sector angle and outer-dot size both encode
  distinct-cell abundance. This is a composition view, not a lineage tree.
- `mst_rooted_hc_depth.svg`: node fill is HC mutational depth.
- `mst_rooted_pn_distance.svg`: node fill is the nucleotide edit distance from
  the family's most P-supported compatible naive junction reconstruction. This
  exposes P/N ambiguity continuously on the **same MST topology**; there is no
  mutation cutoff and the value is not used to decide family membership.
- `mst_rooted_light_chain.svg`: structural LC identity on the paired-receptor
  topology. Different LC identities are separated in the feature space by a
  weighted categorical block rather than inferred from HC sequence alone.
- `mst_rooted_lc_depth.svg`: LC mutational depth on that same paired topology.
- `mst_rooted_paired_depth.svg`: HC + LC mutational depth on that topology.
  Paired HC+LC sub-maps use exactly the same representation after conditioning
  on one structural LC identity.

## Mutational depth

For one receptor call, mutational depth is the number of **aligned nucleotide
substitutions** between Lumrik's reconstructed naive rearrangement and the
error-corrected observed receptor. Valkyrn uses its Needleman-Wunsch alignment
before counting substitutions, so an indel does not shift the remainder of the
sequence into false mismatches. Indels themselves are currently used to establish
the alignment but are **not added to the depth count**.

For an HC sequence state, `hc_depth_nt` is the median HC substitution depth among
cells represented by that state. For LC depth, Valkyrn finds productive IGK/IGL
calls belonging to those cells, computes each LC against **its own reconstructed
naive LC**, and reports the median (`lc_depth_median_nt`) plus minimum and maximum.
In a paired HC+LC sub-map, LC collection is restricted to that structural LC ID.
In a whole-HC map, several productive LCs can contribute to one HC state.

The blue -> yellow -> red scale is normalized to the maximum observed value **in
that plot**. Red therefore means *deepest in this displayed family*, not a fixed
absolute mutation burden across different families. Mutational depth is a measured
molecular divergence; it should not be read directly as chronological cell age.

## Tables

`nodes.tsv` contains the exact abundance, isotype/LC counts and depth values behind
the SVGs. `rows.tsv`, `coords.tsv`, and `tree.tsv` retain the ClonoMap state,
coordinate and MST data for downstream analysis.
"#
    )?;
    let mut summary = writer(dir.join("clonomap_summary.tsv"))?;
    writeln!(
        summary,
        "kind\tfamily\tlight_chain\tcells\tsequences\tstatus\toutput_dir"
    )?;
    let selected: Vec<_> = fams
        .iter()
        .filter(|f| calls[f.members[0]].chain == "IGH")
        .filter(|f| {
            f.members
                .iter()
                .map(|&i| calls[i].cell.as_str())
                .collect::<BTreeSet<_>>()
                .len()
                >= min_size
        })
        .collect();
    eprintln!(
        "Valkyrn: generating {} ClonoMap HC families (>= {} cells; layout={})",
        selected.len(),
        min_size,
        if radial_layout { "radial" } else { "layered" }
    );
    let selected_total = selected.len();
    for (family_no, f) in selected.into_iter().enumerate() {
        let started = std::time::Instant::now();
        let cells: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        eprintln!(
            "  [HC {}/{}] {} — {} cells ...",
            family_no + 1,
            selected_total,
            f.name,
            cells.len()
        );
        let safe = safe_name(&f.name);
        let fdir = dir.join(&safe);
        fs::create_dir_all(&fdir)?;
        let lc_n =
            write_lc_constellation(&fdir.join("lc_repertoire.svg"), f, calls, &lights_by_cell)?;
        match render_group(&fdir, &f.name, &f.members, calls, &lights_by_cell, None, k, radial_layout) {
            Ok(n) => writeln!(
                summary,
                "heavy_family\t{}\t\t{}\t{}\tok;lc_repertoire={}\t{}",
                f.name,
                cells.len(),
                n,
                lc_n,
                safe
            )?,
            Err(e) => {
                writeln!(
                    summary,
                    "heavy_family\t{}\t\t{}\t0\tfailed:{};lc_repertoire={}\t{}",
                    f.name,
                    cells.len(),
                    e.to_string().replace('\t', " "),
                    lc_n,
                    safe
                )?;
            }
        }

        let hc_light_indices: Vec<usize> = cells.iter()
            .filter_map(|cell| lights_by_cell.get(*cell))
            .flat_map(|xs| xs.iter().copied())
            .collect();
        let canonical_lc = canonical_light_ids(calls, &hc_light_indices);
        let mut lc_cells: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for cell in &cells {
            if let Some(ls) = lights_by_cell.get(*cell) {
                for &li in ls {
                    let lc_id = canonical_lc.get(&li).cloned().unwrap_or_else(|| calls[li].id.clone());
                    lc_cells.entry(lc_id).or_default().insert((*cell).to_string());
                }
            }
        }
        for (lc_id, paired_cells) in lc_cells {
            if paired_cells.len() < min_paired_size {
                continue;
            }
            let members: Vec<usize> = f
                .members
                .iter()
                .copied()
                .filter(|&i| paired_cells.contains(&calls[i].cell))
                .collect();
            let psafe = format!("{}__{}", safe, safe_name(&lc_id));
            let pdir = dir.join("paired").join(&psafe);
            let label = format!("{} + {}", f.name, lc_id);
            match render_group(
                &pdir,
                &label,
                &members,
                calls,
                &lights_by_cell,
                Some(&lc_id),
                k,
                radial_layout,
            ) {
                Ok(n) => writeln!(
                    summary,
                    "paired_hc_lc\t{}\t{}\t{}\t{}\tok\tpaired/{}",
                    f.name,
                    lc_id,
                    paired_cells.len(),
                    n,
                    psafe
                )?,
                Err(e) => {
                    writeln!(
                        summary,
                        "paired_hc_lc\t{}\t{}\t{}\t0\tfailed:{}\tpaired/{}",
                        f.name,
                        lc_id,
                        paired_cells.len(),
                        e.to_string().replace('\t', " "),
                        psafe
                    )?;
                }
            }
        }
        eprintln!(
            "  [HC {}/{}] {} done ({:.1}s)",
            family_no + 1,
            selected_total,
            f.name,
            started.elapsed().as_secs_f32()
        );
    }
    Ok(())
}

fn write_report(out: &Path, calls: &[Call], fams: &[Family]) -> Result<()> {
    let cells: BTreeSet<_> = calls.iter().map(|c| c.cell.as_str()).collect();
    let productive = calls.iter().filter(|c| c.productive()).count();
    let hc: Vec<_> = fams
        .iter()
        .filter(|f| calls[f.members[0]].chain == "IGH")
        .collect();
    let expanded = hc.iter().filter(|f| f.members.len() > 1).count();
    let mut diverse = 0;
    for f in &hc {
        let cs: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        if light_partners(calls, &cs).len() > 1 {
            diverse += 1
        }
    }
    let mut w = writer(out.join("README.txt"))?;
    writeln!(
        w,
        "Valkyrn repertoire interpretation\n===============================\n\nCells: {}\nRearrangements: {}\nProductive rearrangements: {}\nProductive receptor families: {}\nIGH families: {}\nExpanded IGH families: {}\nIGH families with >1 productive LC partner: {}\n\nInterpretation notes\n--------------------\nIGH family membership requires the same V and J calls plus a CDR3 nucleotide edit distance no greater than --max-cdr3-distance between every pair of family members. This is a hard boundary: structural HC:<HEX> IDs, P/N decomposition and mutation-distance measurements cannot override it. IGK/IGL families likewise require the same light-chain type, V and J plus at least the configured maximum pairwise CDR3 edit distance, and are only formed inside the same inferred IGH background; a similar light chain across different heavy backgrounds is treated as recurrence, not clonal evidence. Mutation classes compare Lumrik's reconstructed naive rearrangement with its error-corrected observed receptor. A mutation shared by all members of a family is a lineage-shared candidate, not proof of somatic hypermutation: recurrent changes across unrelated families using the same germline segment may indicate an unrepresented germline allele. Structure candidates are deliberately restricted to expanded IGH families with multiple observed productive light-chain partners.\n",
        cells.len(),
        calls.len(),
        productive,
        fams.len(),
        hc.len(),
        expanded,
        diverse
    )?;
    Ok(())
}

fn writer(path: impl AsRef<Path>) -> Result<BufWriter<File>> {
    let p = path.as_ref();
    File::create(p)
        .map(BufWriter::new)
        .with_context(|| format!("creating {}", p.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mutation_diff() {
        let m = measure_mutations_nw("AACCGG", "AATCGA").unwrap().mutations;
        assert_eq!(m.len(), 2);
        assert_eq!(m[0].pos, 2);
        assert_eq!(m[1].pos, 5);
    }
    #[test]
    fn naive_cdr3_is_not_subject_to_family_distance_threshold() {
        // --max-cdr3-distance constrains observed family members, not SHM from
        // naive. Four substitutions in the naive CDR3 must still establish an
        // ungapped coordinate anchor.
        let m = measure_mutations_nw("AAAACCCCGGGG", "AAAATTTTGGGG")
            .unwrap()
            .mutations;
        assert_eq!(m.len(), 4);
    }

    #[test]
    fn left_and_right_are_anchored_independently() {
        // Extra observed prefix must not shift the J-side comparison.
        let m = measure_mutations_nw("AACCGGTT", "XXAACCGGTA")
            .unwrap()
            .mutations;
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].pos, 7);
    }

    #[test]
    fn anchored_comparison_never_shifts_after_a_mismatch() {
        let m = measure_mutations_nw("AACCGGTT", "AATCGGTA")
            .unwrap()
            .mutations;
        assert_eq!(m.len(), 2);
        assert_eq!(m[0].pos, 2);
        assert_eq!(m[1].pos, 7);
    }

    #[test]
    fn internal_three_nt_indel_is_one_mutation_event() {
        let aln = measure_mutations_nw("AAAACCCCGGGG", "AAAATTTCCCCGGGG").unwrap();
        assert_eq!(aln.mutations.len(), 0);
        assert_eq!(aln.indels.len(), 1);
        assert_eq!(aln.indels[0].kind, IndelKind::Insertion);
        assert_eq!(aln.indels[0].len, 3);
        assert_eq!(aln.event_count(), 1);
    }

    #[test]
    fn terminal_truncation_is_not_an_indel_event() {
        let aln = measure_mutations_nw("ATGAACCGGTT", "GAACCGGTT").unwrap();
        assert!(aln.indels.is_empty());
        assert_eq!(aln.event_count(), 0);
    }

    #[test]
    fn complete_germline_is_trimmed_to_observed_internal_fragment() {
        let aln = measure_mutations_nw("TTTTAAAACCCCGGGGAAAA", "AAAACCCCGGGG").unwrap();
        assert!(aln.mutations.is_empty());
        assert!(aln.indels.is_empty());
        assert_eq!(aln.event_count(), 0);
    }

    #[test]
    fn unrelated_fragment_has_no_mutation_measurement() {
        assert!(measure_mutations_nw("AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC").is_none());
    }

    #[test]
    fn cdr3_distance_has_a_hard_configurable_boundary() {
        assert!(cdr3_compatible("AACCGGTT", "AACCGGTT", 3));
        assert!(cdr3_compatible("AACCGGTT", "AATCAGTA", 3));
        assert!(!cdr3_compatible("AACCGGTT", "AATCAATA", 3));
        assert!(cdr3_compatible("AACCGGTT", "AATCAATA", 4));
        assert!(!cdr3_compatible("", "AACCGGTT", 3));
    }

    #[test]
    fn complete_link_prevents_cdr3_chaining() {
        // 0~1 and 1~2 at distance <=2, but 0!~2. Complete-link must not
        // bridge the two endpoints through the middle sequence.
        let seqs = ["AAAAAAAA", "AAAAAACC", "AAAACCCC"];
        let groups = complete_link_groups(&[0, 1, 2], |a, b| {
            cdr3_compatible(seqs[a], seqs[b], 2)
        });
        assert_eq!(groups.len(), 2);
        assert!(groups.iter().all(|g| g.iter().all(|&a| {
            g.iter().all(|&b| cdr3_compatible(seqs[a], seqs[b], 2))
        })));
    }
}
