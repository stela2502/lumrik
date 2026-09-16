use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

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
    pub productivity: String,
    pub support: u64,
    pub rediscovery: u64,
    pub naive: String,
    pub observed: String,
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

/// Align the reconstructed naive receptor to the observed receptor before
/// comparing bases.  The previous implementation zipped the strings directly,
/// so one indel shifted every downstream coordinate and manufactured hundreds
/// of false substitutions.
///
/// This is a compact Needleman-Wunsch edit alignment. Receptor sequences are
/// short enough that the O(n*m) matrix is tiny compared with VDJ reconstruction.
/// Mutation coordinates are always in the naive/ancestral sequence coordinate
/// system. Insertions/deletions establish the alignment but are not emitted as
/// substitutions here; they can become first-class lineage events later.
pub fn mutations(naive: &str, observed: &str) -> Vec<Mutation> {
    let a = naive.as_bytes();
    let b = observed.as_bytes();
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }

    let cols = b.len() + 1;
    let mut score = vec![0u32; (a.len() + 1) * cols];
    for i in 0..=a.len() {
        score[i * cols] = i as u32;
    }
    for j in 0..=b.len() {
        score[j] = j as u32;
    }

    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let subst = score[(i - 1) * cols + (j - 1)] + u32::from(a[i - 1] != b[j - 1]);
            let delete = score[(i - 1) * cols + j] + 1;
            let insert = score[i * cols + (j - 1)] + 1;
            score[i * cols + j] = subst.min(delete).min(insert);
        }
    }

    let mut aligned = Vec::new();
    let (mut i, mut j) = (a.len(), b.len());
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let subst_cost = u32::from(a[i - 1] != b[j - 1]);
            if score[i * cols + j] == score[(i - 1) * cols + (j - 1)] + subst_cost {
                aligned.push((Some(i - 1), Some(j - 1)));
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && score[i * cols + j] == score[(i - 1) * cols + j] + 1 {
            aligned.push((Some(i - 1), None));
            i -= 1;
        } else {
            aligned.push((None, Some(j - 1)));
            j -= 1;
        }
    }
    aligned.reverse();

    aligned
        .into_iter()
        .filter_map(|(ai, bj)| {
            let (ai, bj) = (ai?, bj?);
            let from = a[ai].to_ascii_uppercase();
            let to = b[bj].to_ascii_uppercase();
            (from != to && from != b'N' && to != b'N').then_some(Mutation { pos: ai, from, to })
        })
        .collect()
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
            productivity: get(productivity).into(),
            support: get(support).parse().unwrap_or(0),
            rediscovery: get(rediscovery).parse().unwrap_or(0),
            naive: get(naive).into(),
            observed: get(observed).into(),
            cdr3_aa: airr.get(&rid).cloned().unwrap_or_default(),
        });
    }
    Ok(out)
}

fn read_airr_cdr3(path: &Path) -> Result<HashMap<String, String>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines.next().context("airr_rearrangements.tsv is empty")??;
    let h: Vec<&str> = header.split('\t').collect();
    let seq = h
        .iter()
        .position(|x| *x == "lumrik_recombination_id")
        .context("AIRR lacks lumrik_recombination_id")?;
    let cdr = h
        .iter()
        .position(|x| *x == "cdr3_aa")
        .context("AIRR lacks cdr3_aa")?;
    let mut out = HashMap::new();
    for line in lines {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if let (Some(id), Some(c)) = (f.get(seq), f.get(cdr)) {
            out.insert((*id).into(), (*c).into());
        }
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct Family {
    pub name: String,
    pub members: Vec<usize>,
}

/// Build heavy-chain families from Lumrik's structural recombination ID.
///
/// HC:<HEX> is already the compact, reversible representation of V/D/J plus
/// the measured junction architecture. It is therefore the primary clone key;
/// Valkyrn must not throw that information away and re-cluster on V/J+CDR3 AA.
///
/// `pn_alternative` marks junctions where P/N decomposition is not unique. For
/// those calls only, a conservative fallback can connect otherwise distinct
/// compact IDs when V/D/J, CDR3-AA length and reconstructed naive length agree.
/// The exact compact IDs remain visible on every member in valkyrn_mutations.tsv.
fn heavy_families(calls: &[Call]) -> Vec<Family> {
    let eligible: Vec<usize> = calls
        .iter()
        .enumerate()
        .filter(|(_, c)| c.productive() && c.chain == "IGH" && c.id.starts_with("HC:"))
        .map(|(i, _)| i)
        .collect();

    let mut parent: Vec<usize> = (0..eligible.len()).collect();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }

    let mut by_id: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (local, &global) in eligible.iter().enumerate() {
        by_id
            .entry(calls[global].id.as_str())
            .or_default()
            .push(local);
    }
    for locals in by_id.values() {
        if let Some((&first, rest)) = locals.split_first() {
            for &x in rest {
                let a = root(&mut parent, first);
                let b = root(&mut parent, x);
                if a != b {
                    parent[b] = a;
                }
            }
        }
    }

    // Only ambiguous P/N decompositions get a fallback. Exact structural IDs
    // remain the normal path. Requiring the same V/D/J, CDR3 length and naive
    // rearrangement length prevents the old broad V/J+CDR3-distance clustering.
    let ambiguous: Vec<usize> = eligible
        .iter()
        .enumerate()
        .filter(|(_, g)| calls[**g].pn_alternative)
        .map(|(l, _)| l)
        .collect();
    for a in 0..ambiguous.len() {
        for b in a + 1..ambiguous.len() {
            let la = ambiguous[a];
            let lb = ambiguous[b];
            let ca = &calls[eligible[la]];
            let cb = &calls[eligible[lb]];
            if ca.id != cb.id
                && ca.v == cb.v
                && ca.d == cb.d
                && ca.j == cb.j
                && !ca.cdr3_aa.is_empty()
                && ca.cdr3_aa.len() == cb.cdr3_aa.len()
                && ca.naive.len() == cb.naive.len()
            {
                let ra = root(&mut parent, la);
                let rb = root(&mut parent, lb);
                if ra != rb {
                    parent[rb] = ra;
                }
            }
        }
    }

    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (local, &global) in eligible.iter().enumerate() {
        let r = root(&mut parent, local);
        groups.entry(r).or_default().push(global);
    }
    let mut out = Vec::new();
    for members in groups.into_values() {
        let mut ids: BTreeSet<&str> = members.iter().map(|&i| calls[i].id.as_str()).collect();
        let representative = ids.pop_first().unwrap_or("HC:UNKNOWN");
        let name = if ids.is_empty() {
            representative.to_string()
        } else {
            format!("{}~PNALT", representative)
        };
        out.push(Family { name, members });
    }
    out
}

/// Build repertoire groups with deliberately asymmetric trust.
///
/// Heavy-chain clone identity comes from sc-vdj's compact structural HC:<HEX>
/// recombination ID. Light-chain LC:<HEX> identity is useful, but is only
/// treated as a multi-cell clone in the context of the same HC family.
pub fn families(calls: &[Call], _max_cdr3_distance: usize) -> Vec<Family> {
    let mut out = heavy_families(calls);

    let mut heavy_family_by_cell: HashMap<&str, &str> = HashMap::new();
    for f in &out {
        for &i in &f.members {
            heavy_family_by_cell
                .entry(calls[i].cell.as_str())
                .or_insert(f.name.as_str());
        }
    }

    let mut light: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, c) in calls.iter().enumerate().filter(|(_, c)| {
        c.productive() && matches!(c.chain.as_str(), "IGK" | "IGL") && c.id.starts_with("LC:")
    }) {
        let bg = heavy_family_by_cell
            .get(c.cell.as_str())
            .map(|x| (*x).to_string())
            .unwrap_or_else(|| format!("UNPAIRED:{}", c.cell));
        light.entry((bg, c.id.clone())).or_default().push(i);
    }
    for ((bg, lc), members) in light {
        let name = if bg.starts_with("UNPAIRED:") {
            format!("{}@{}", lc, bg)
        } else {
            format!("{}+{}", bg, lc)
        };
        out.push(Family { name, members });
    }
    out.sort_by_key(|f| std::cmp::Reverse(f.members.len()));
    out
}

pub fn analyze(
    vdj_dir: &Path,
    out_dir: &Path,
    max_cdr3_distance: usize,
    min_structure_family: usize,
    threads: usize,
    min_clonomap_family: usize,
    min_clonomap_paired_family: usize,
    clonomap_k: usize,
) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("building Valkyrn Rayon pool")?;
    pool.install(|| {
        analyze_inner(
            vdj_dir,
            out_dir,
            max_cdr3_distance,
            min_structure_family,
            min_clonomap_family,
            min_clonomap_paired_family,
            clonomap_k,
        )
    })
}

fn analyze_inner(
    vdj_dir: &Path,
    out_dir: &Path,
    max_cdr3_distance: usize,
    min_structure_family: usize,
    min_clonomap_family: usize,
    min_clonomap_paired_family: usize,
    clonomap_k: usize,
) -> Result<()> {
    fs::create_dir_all(out_dir)?;
    let calls = read_calls(vdj_dir)?;
    if calls.is_empty() {
        bail!("no VDJ calls found")
    }
    let fams = families(&calls, max_cdr3_distance);
    // Needleman-Wunsch is by far the expensive operation. Compute every
    // rearrangement exactly once, in parallel, then reuse the immutable cache
    // for family summaries, recurrence classification and TSV output.
    let mutation_cache: Vec<Vec<Mutation>> = calls
        .par_iter()
        .map(|c| mutations(&c.naive, &c.observed))
        .collect();
    write_families(out_dir, &calls, &fams, &mutation_cache)?;
    write_mutations(out_dir, &calls, &fams, &mutation_cache)?;
    write_cell_qc(out_dir, &calls, &mutation_cache)?;
    write_structure_candidates(out_dir, &calls, &fams, min_structure_family)?;
    write_clonomap(
        out_dir,
        &calls,
        &fams,
        min_clonomap_family,
        min_clonomap_paired_family,
        clonomap_k,
    )?;
    write_report(out_dir, &calls, &fams)?;
    Ok(())
}

fn write_families(
    out: &Path,
    calls: &[Call],
    fams: &[Family],
    mutation_cache: &[Vec<Mutation>],
) -> Result<()> {
    let mut w = writer(out.join("valkyrn_families.tsv"))?;
    writeln!(
        w,
        "family\tchain\tcells\tv\td\tj\tstructural_recombination_ids\tcdr3_aa_variants\tisotypes\tisotype_counts\tlight_partners\tdominant_light_fraction\tshared_mutations\tvariable_mutations"
    )?;
    for f in fams {
        let cells: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        let c = &calls[f.members[0]];
        let cdr: BTreeSet<_> = f
            .members
            .iter()
            .map(|&i| calls[i].cdr3_aa.as_str())
            .collect();
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
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            f.name,
            c.chain,
            cells.len(),
            c.v,
            c.d,
            c.j,
            structural.len(),
            cdr.len(),
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
    let mut x = BTreeMap::new();
    for c in calls.iter().filter(|c| {
        cells.contains(c.cell.as_str())
            && c.productive()
            && matches!(c.chain.as_str(), "IGK" | "IGL")
    }) {
        *x.entry(c.id.clone()).or_default() += 1;
    }
    x
}

fn family_mutation_sets(
    mutation_cache: &[Vec<Mutation>],
    f: &Family,
) -> (BTreeSet<Mutation>, BTreeSet<Mutation>) {
    let sets: Vec<BTreeSet<Mutation>> = f
        .members
        .iter()
        .map(|&i| mutation_cache[i].iter().cloned().collect())
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
    mutation_cache: &[Vec<Mutation>],
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
            for m in mutation_cache[i].iter().cloned() {
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
            for m in mutation_cache[i].iter().cloned() {
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

fn write_cell_qc(out: &Path, calls: &[Call], mutation_cache: &[Vec<Mutation>]) -> Result<()> {
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
            .map(|i| mutation_cache[i].len())
            .collect();
        let lc_depths: Vec<usize> = idxs
            .iter()
            .copied()
            .filter(|&i| {
                calls[i].productive() && matches!(calls[i].chain.as_str(), "IGK" | "IGL")
            })
            .map(|i| mutation_cache[i].len())
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

fn observed_on_naive_coordinates(naive: &str, observed: &str) -> String {
    let a = naive.as_bytes();
    let b = observed.as_bytes();
    if a.is_empty() {
        return String::new();
    }
    let cols = b.len() + 1;
    let mut score = vec![0u32; (a.len() + 1) * cols];
    for i in 0..=a.len() {
        score[i * cols] = i as u32;
    }
    for j in 0..=b.len() {
        score[j] = j as u32;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let subst = score[(i - 1) * cols + j - 1] + u32::from(a[i - 1] != b[j - 1]);
            score[i * cols + j] = subst
                .min(score[(i - 1) * cols + j] + 1)
                .min(score[i * cols + j - 1] + 1);
        }
    }
    let mut out = vec![b'-'; a.len()];
    let (mut i, mut j) = (a.len(), b.len());
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = u32::from(a[i - 1] != b[j - 1]);
            if score[i * cols + j] == score[(i - 1) * cols + j - 1] + cost {
                out[i - 1] = b[j - 1].to_ascii_uppercase();
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && score[i * cols + j] == score[(i - 1) * cols + j] + 1 {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    String::from_utf8(out).unwrap_or_default()
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
                m.depths.push(mutations(&lc.naive, &lc.observed).len());
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
) -> Result<()> {
    use clonomap::ClonoMap;
    #[derive(Default)]
    struct StateMeta {
        cells: BTreeSet<String>,
        isotypes: BTreeMap<String, usize>,
        lights: BTreeMap<String, usize>,
        hc_depths: Vec<usize>,
        lc_depths: Vec<usize>,
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
    ) -> Result<usize> {
        use ndarray::Array2;

        const LC_ID_WEIGHT: f32 = 8.0;

        let Some(&anchor_i) = member_indices.first() else {
            return Ok(0);
        };
        let hc_width = calls[anchor_i].naive.len();

        // The old whole-HC map built geometry from HC sequence alone and painted
        // LC identity on afterwards. That necessarily allowed an HC edge to look
        // like an LC-A -> LC-B transition. Build the geometry from the paired
        // receptor state instead: HC mutation state + LC mutation state + a
        // weighted one-hot structural LC identity.
        //
        // LC mutation coordinates are relative to each LC's own reconstructed
        // naive sequence. The categorical block tells the metric that two
        // different structural LC rearrangements are not ordinary mutations of
        // one another.
        #[derive(Clone)]
        struct PairedState {
            label: String,
            hc_seq: String,
            lc_id: String,
            hc_mut: Vec<f32>,
            lc_mut: Vec<f32>,
        }

        let mut lc_ids: BTreeSet<String> = BTreeSet::new();
        let mut lc_width = 0usize;
        for &i in member_indices {
            if let Some(ls) = lights_by_cell.get(&calls[i].cell) {
                for &li in ls {
                    let lc = &calls[li];
                    if required_lc.is_some_and(|want| lc.id != want) {
                        continue;
                    }
                    lc_ids.insert(lc.id.clone());
                    lc_width = lc_width.max(lc.naive.len());
                }
            }
        }
        if lc_ids.is_empty() {
            lc_ids.insert("unpaired".to_string());
        }
        let lc_id_col: BTreeMap<String, usize> = lc_ids
            .iter()
            .cloned()
            .enumerate()
            .map(|(i, id)| (id, i))
            .collect();

        let mutation_vector = |naive: &str, observed: &str, width: usize| -> Vec<f32> {
            let mut v = vec![0.0; width];
            for m in mutations(naive, observed) {
                if m.pos < width {
                    v[m.pos] = 1.0;
                }
            }
            v
        };

        let mut states: BTreeMap<String, PairedState> = BTreeMap::new();
        let mut raw_meta: HashMap<String, StateMeta> = HashMap::new();
        for &i in member_indices {
            let hc = &calls[i];
            if hc.naive.len() != hc_width || hc.observed.is_empty() {
                continue;
            }
            let hc_seq = observed_on_naive_coordinates(&hc.naive, &hc.observed);
            if hc_seq.len() != hc_width {
                continue;
            }
            let hc_mut = mutation_vector(&hc.naive, &hc.observed, hc_width);
            let mut paired_any = false;

            if let Some(ls) = lights_by_cell.get(&hc.cell) {
                for &li in ls {
                    let lc = &calls[li];
                    if required_lc.is_some_and(|want| lc.id != want) {
                        continue;
                    }
                    if lc.observed.is_empty() {
                        continue;
                    }
                    paired_any = true;
                    let lc_mut = mutation_vector(&lc.naive, &lc.observed, lc_width);
                    let lc_pattern = lc_mut
                        .iter()
                        .map(|x| if *x > 0.0 { '1' } else { '0' })
                        .collect::<String>();
                    let label = format!("{}|{}|{}", hc_seq, lc.id, lc_pattern);
                    states.entry(label.clone()).or_insert_with(|| PairedState {
                        label: label.clone(),
                        hc_seq: hc_seq.clone(),
                        lc_id: lc.id.clone(),
                        hc_mut: hc_mut.clone(),
                        lc_mut: lc_mut.clone(),
                    });
                    let m = raw_meta.entry(label).or_default();
                    m.cells.insert(hc.cell.clone());
                    m.hc_depths.push(mutations(&hc.naive, &hc.observed).len());
                    let iso = if hc.c.is_empty() {
                        "unknown"
                    } else {
                        hc.c.as_str()
                    };
                    *m.isotypes.entry(iso.to_string()).or_default() += 1;
                    *m.lights.entry(lc.id.clone()).or_default() += 1;
                    m.lc_depths.push(mutations(&lc.naive, &lc.observed).len());
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
                    lc_mut: vec![0.0; lc_width],
                });
                let m = raw_meta.entry(label).or_default();
                m.cells.insert(hc.cell.clone());
                m.hc_depths.push(mutations(&hc.naive, &hc.observed).len());
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

        let rows: Vec<PairedState> = states.into_values().collect();
        let feature_cols = hc_width + lc_width + lc_id_col.len();
        let mut features = Array2::<f32>::zeros((rows.len(), feature_cols));
        for (r, state) in rows.iter().enumerate() {
            for (j, x) in state.hc_mut.iter().enumerate() {
                features[[r, j]] = *x;
            }
            for (j, x) in state.lc_mut.iter().enumerate() {
                features[[r, hc_width + j]] = *x;
            }
            if let Some(&j) = lc_id_col.get(&state.lc_id) {
                features[[r, hc_width + lc_width + j]] = LC_ID_WEIGHT;
            }
        }
        fs::create_dir_all(fdir)?;
        let n_seqs = rows.len();
        let model = ClonoMap::from_feature_matrix(
            rows.iter().map(|x| x.label.clone()).collect(),
            features,
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

        // Root at the observed paired state closest to its reconstructed naive
        // receptors. This remains a virtual display root; it is not an observed
        // cell and does not create a biological LC-to-LC transition.
        let (naive_root, naive_dist) = model
            .encoder
            .sequences
            .dna
            .iter()
            .enumerate()
            .map(|(i, key)| {
                let m = raw_meta.get(key);
                let h = m
                    .and_then(|x| median(&x.hc_depths))
                    .unwrap_or(f32::INFINITY);
                let l = m.and_then(|x| median(&x.lc_depths)).unwrap_or(0.0);
                (i, h + l)
            })
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, d)| (i, if d.is_finite() { d } else { 0.0 }))
            .unwrap_or((0, 0.0));

        let mut iso_cat = Vec::new();
        let mut iso_mixed = Vec::new();
        let mut lc_cat = Vec::new();
        let mut lc_mixed = Vec::new();
        let mut abundance = Vec::new();
        let mut hc_depth = Vec::new();
        let mut lc_depth = Vec::new();
        let mut paired_depth = Vec::new();
        let mut nodes = writer(fdir.join("nodes.tsv"))?;
        writeln!(
            nodes,
            "node\tstate\thc_dna\tcells\tisotypes\tdominant_isotype\tlight_chains\tdominant_light_chain\thc_depth_nt\tlc_depth_median_nt\tlc_depth_min_nt\tlc_depth_max_nt\tpaired_depth_median_nt"
        )?;
        for (node, key) in model.encoder.sequences.dna.iter().enumerate() {
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
            let hc_dna = rows.get(node).map(|x| x.hc_seq.as_str()).unwrap_or("");
            writeln!(
                nodes,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                node,
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
                pair_s
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
                naive_dist,
                &iso_cat,
                &iso_mixed,
                &abundance,
                "IGH constant class",
                &title,
                fdir.join("mst_rooted_isotype.svg")
                    .to_string_lossy()
                    .as_ref(),
            )
            .map_err(|e| anyhow::anyhow!("failed to plot isotype-rooted MST: {e}"))?;
        model
            .tree
            .plot_rooted_annotated_cached(
                model.coords().nrows(),
                naive_root,
                naive_dist,
                &lc_cat,
                &lc_mixed,
                &abundance,
                "Structural light-chain identity",
                &title,
                fdir.join("mst_rooted_light_chain.svg")
                    .to_string_lossy()
                    .as_ref(),
            )
            .map_err(|e| anyhow::anyhow!("failed to plot light-chain-rooted MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                naive_dist,
                &hc_depth,
                &abundance,
                "HC mutational depth (nt from NAIVE)",
                &title,
                fdir.join("mst_rooted_hc_depth.svg")
                    .to_string_lossy()
                    .as_ref(),
            )
            .map_err(|e| anyhow::anyhow!("failed to plot HC-depth MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                naive_dist,
                &lc_depth,
                &abundance,
                "Linked LC mutational depth (median nt)",
                &title,
                fdir.join("mst_rooted_lc_depth.svg")
                    .to_string_lossy()
                    .as_ref(),
            )
            .map_err(|e| anyhow::anyhow!("failed to plot LC-depth MST: {e}"))?;
        model
            .tree
            .plot_rooted_continuous_cached(
                model.coords().nrows(),
                naive_root,
                naive_dist,
                &paired_depth,
                &abundance,
                "Paired HC+LC mutational depth (median nt)",
                &title,
                fdir.join("mst_rooted_paired_depth.svg")
                    .to_string_lossy()
                    .as_ref(),
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
the fact. The virtual `NAIVE` root is attached to the observed paired state with
the smallest measured HC+LC substitution depth. The topology is a paired
sequence-state landscape, not a claim of chronological phylogeny.

## Visual encoding

- **Node size = clone/state abundance**: the number of distinct cells represented
  by that observed HC sequence state. The same size encoding is used in every
  annotation view.
- **Edges = the same cached ClonoMap MST** in every view. Annotation never changes
  the topology.
- **Thin black ring = mixed state** in categorical views: cells with the same HC
  sequence state carry more than one isotype or productive LC identity.
- `mst_rooted_isotype.svg`: node fill is IGH constant-region class/isotype.
- `lc_repertoire.svg` (whole-HC families): a cheap radial repertoire overview.
  Every structural LC is an independent sector attached to the HC family; there
  are deliberately no LC-to-LC edges. Sector angle and outer-dot size both encode
  distinct-cell abundance. This is a composition view, not a lineage tree.
- `mst_rooted_hc_depth.svg`: node fill is HC mutational depth.
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
    for f in fams.iter().filter(|f| calls[f.members[0]].chain == "IGH") {
        let cells: BTreeSet<_> = f.members.iter().map(|&i| calls[i].cell.as_str()).collect();
        if cells.len() < min_size {
            continue;
        }
        let safe = safe_name(&f.name);
        let fdir = dir.join(&safe);
        fs::create_dir_all(&fdir)?;
        let lc_n =
            write_lc_constellation(&fdir.join("lc_repertoire.svg"), f, calls, &lights_by_cell)?;
        match render_group(&fdir, &f.name, &f.members, calls, &lights_by_cell, None, k) {
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

        let mut lc_cells: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for cell in &cells {
            if let Some(ls) = lights_by_cell.get(*cell) {
                for &li in ls {
                    lc_cells
                        .entry(calls[li].id.clone())
                        .or_default()
                        .insert((*cell).to_string());
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
        "Valkyrn repertoire interpretation\n===============================\n\nCells: {}\nRearrangements: {}\nProductive rearrangements: {}\nProductive receptor families: {}\nIGH families: {}\nExpanded IGH families: {}\nIGH families with >1 productive LC partner: {}\n\nInterpretation notes\n--------------------\nIGH family membership is defined primarily by Lumrik's reversible structural HC:<HEX> recombination identifier. Only calls flagged pn_alternative may use a conservative fallback requiring the same V/D/J, CDR3-AA length and naive rearrangement length. IGK/IGL LC:<HEX> identity is only allowed to form a multi-cell clone inside the same inferred IGH background; a shared light-chain rearrangement across different heavy backgrounds is treated as recurrence, not clonal evidence. Mutation classes compare Lumrik's reconstructed naive rearrangement with its error-corrected observed receptor. A mutation shared by all members of a family is a lineage-shared candidate, not proof of somatic hypermutation: recurrent changes across unrelated families using the same germline segment may indicate an unrepresented germline allele. Structure candidates are deliberately restricted to expanded IGH families with multiple observed productive light-chain partners.\n",
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
        let m = mutations("AACCGG", "AATCGA");
        assert_eq!(m.len(), 2);
        assert_eq!(m[0].pos, 2);
        assert_eq!(m[1].pos, 5);
    }
    #[test]
    fn insertion_does_not_shift_mutations() {
        let m = mutations("AACCGGTT", "AACTCGGTA");
        assert_eq!(m.len(), 1);
        assert_eq!(
            m[0],
            Mutation {
                pos: 7,
                from: b'T',
                to: b'A'
            }
        );
    }
    #[test]
    fn deletion_does_not_shift_mutations() {
        let m = mutations("AACCTGGTT", "AACCGGTA");
        assert_eq!(m.len(), 1);
        assert_eq!(
            m[0],
            Mutation {
                pos: 8,
                from: b'T',
                to: b'A'
            }
        );
    }
}
