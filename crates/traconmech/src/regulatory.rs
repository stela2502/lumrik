//! Deterministic TP53/STRING/MSigDB projection onto experimental expression.
//! STRING associations and MSigDB membership are not proof of direct regulation.
use anyhow::{bail, Context, Result};
use sc_analysis::{SingleCellData, write_group_stats};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub struct Config<'a> {
    pub gene: &'a str,
    pub index: &'a Path,
    pub data: &'a SingleCellData,
    pub out: &'a Path,
    pub string_info: &'a Path,
    pub string_links: &'a Path,
    pub cluster: &'a str,
    pub condition: &'a str,
    pub cases: &'a str,
    pub controls: &'a str,
    pub effect_cluster: bool,
    pub label_angle: Option<f64>,
    pub count_label_layout: &'a str,
    pub co_variant: &'a [String],
    pub min_score: u32,
    pub min_wt_cells: usize,
    pub ai_development: bool,
    pub ai_model: &'a str,
    pub ai_url: &'a str,
    pub ai_primer: Option<&'a str>,
}

fn lines(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().is_some_and(|e| e == "gz") {
        use flate2::read::MultiGzDecoder;
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else { Ok(Box::new(BufReader::new(file))) }
}

fn selected(raw: &str) -> HashSet<String> {
    raw.split(',').map(str::trim).filter(|x| !x.is_empty()).map(str::to_owned).collect()
}

fn stats(data: &SingleCellData, gene: &str, cells: &[usize]) -> (f64, f64, f64) {
    let Some(idx) = data.feature_index(gene) else { return (0.0, 0.0, 0.0); };
    if cells.is_empty() { return (0.0, 0.0, 0.0); }
    let Some(row) = data.matrix.outer_view(idx) else { return (0.0, 0.0, 0.0); };
    let mut sum = 0.0;
    let mut sumsq = 0.0;
    let mut detected = 0;
    for &cell in cells {
        let x = row.get(cell).copied().unwrap_or(0.0) as f64;
        sum += x; sumsq += x*x;
        if x != 0.0 { detected += 1; }
    }
    let n = cells.len() as f64;
    let mean = sum/n;
    (mean, (sumsq/n - mean*mean).max(0.0).sqrt(), detected as f64/n)
}


// Bounded inverse residual variability. Zero-dominated genes receive no weight.
// The expression-dependent floor prevents near-zero SD from dominating.
fn stability_weight(cm: f64, cs: f64, wm: f64, ws: f64, cd: f64, wd: f64) -> f64 {
    if cm.max(wm) < 0.05 || cd.max(wd) < 0.05 { return 0.0; }
    let dispersion = (cs * cs / (cm + 0.1) + ws * ws / (wm + 0.1)) * 0.5;
    1.0 / (1.0 + dispersion)
}

// Cluster genes by the *shape* of their DBA/control effect trajectory.
// Deterministic spherical k-means: normalized profiles, fixed initialization,
// no random seeds. A flat trajectory belongs to module 0.
fn effect_modules(profiles: &BTreeMap<String, Vec<f64>>) -> BTreeMap<String, usize> {
    let mut assignments = BTreeMap::new();
    let mut rows: Vec<(&str, Vec<f64>)> = Vec::new();
    for (gene, values) in profiles {
        if values.len() < 2 {
            assignments.insert(gene.clone(), 0);
            continue;
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let centered: Vec<f64> = values.iter().map(|v| v - mean).collect();
        let norm = centered.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < 1e-9 {
            assignments.insert(gene.clone(), 0);
        } else {
            rows.push((gene, centered.iter().map(|v| v / norm).collect()));
        }
    }
    if rows.is_empty() { return assignments; }
    let k = rows.len().min(8);
    // Farthest-first seeds, with gene-name ordering as a deterministic tie break.
    let mut centers = vec![rows[0].1.clone()];
    while centers.len() < k {
        let mut best = (f64::NEG_INFINITY, 0usize);
        for (i, (_, profile)) in rows.iter().enumerate() {
            let distance = centers.iter().map(|center| {
                profile.iter().zip(center).map(|(a,b)| (a-b)*(a-b)).sum::<f64>()
            }).fold(f64::INFINITY, f64::min);
            if distance > best.0 { best = (distance, i); }
        }
        if best.0 < 1e-12 { break; }
        centers.push(rows[best.1].1.clone());
    }
    let mut labels = vec![0usize; rows.len()];
    for _ in 0..40 {
        let mut changed = false;
        for (i, (_, profile)) in rows.iter().enumerate() {
            let (best, _) = centers.iter().enumerate().map(|(j, center)| {
                (j, profile.iter().zip(center).map(|(a,b)| a*b).sum::<f64>())
            }).max_by(|a,b| a.1.total_cmp(&b.1).then_with(|| b.0.cmp(&a.0))).unwrap();
            if labels[i] != best { labels[i] = best; changed = true; }
        }
        let mut sums = vec![vec![0.0; centers[0].len()]; centers.len()];
        let mut counts = vec![0usize; centers.len()];
        for ((_, profile), &label) in rows.iter().zip(&labels) {
            counts[label] += 1;
            for (sum, value) in sums[label].iter_mut().zip(profile) { *sum += value; }
        }
        for (j, sum) in sums.into_iter().enumerate() {
            if counts[j] == 0 { continue; }
            let norm = sum.iter().map(|v| v*v).sum::<f64>().sqrt();
            if norm > 1e-9 { centers[j] = sum.iter().map(|v| v/norm).collect(); }
        }
        if !changed { break; }
    }
    for ((gene, _), label) in rows.into_iter().zip(labels) {
        assignments.insert(gene.to_owned(), label + 1);
    }
    assignments
}

pub fn run(c: Config<'_>) -> Result<()> {
    fs::create_dir_all(c.out)?;
    // Reuse sc-analysis cluster-vs-rest marker statistics on the ORIGINAL
    // cluster annotations, independent of the experimental grouping.
    // The standard output includes top 20 enriched and depleted markers.
    let marker_dir = c.out.join("cluster_markers");
    if !marker_dir.join("top20_up_down.tsv").is_file()
        || !marker_dir.join("cluster_labels.tsv").is_file()
    {
        eprintln!("[traconmech] sc-analysis markers for {}", c.cluster);
        write_group_stats(&marker_dir, c.data, c.cluster, 10)?;
    }
    let clusters = c.data.annotations.column(c.cluster)
        .with_context(|| format!("missing cluster column {}", c.cluster))?;
    let conditions = c.data.annotations.column(c.condition)
        .with_context(|| format!("missing condition column {}", c.condition))?;
    let cases = selected(c.cases);
    let controls = selected(c.controls);
    if cases.is_empty() || controls.is_empty() || !cases.is_disjoint(&controls) {
        bail!("--case and --control must be nonempty, disjoint lists");
    }
    let mut groups: BTreeMap<String, (Vec<usize>,Vec<usize>)> = BTreeMap::new();
    for (i,(cluster,condition)) in clusters.iter().zip(conditions).enumerate() {
        let entry = groups.entry(cluster.clone()).or_default();
        if cases.contains(condition) { entry.0.push(i); }
        if controls.contains(condition) { entry.1.push(i); }
    }
    // Resolve a WT reference independently for EACH Louvain population.
    // Prefer its own WT cells; only substitute if fewer than min_wt_cells.
    // Rank alternatives by their reference-only profiles, never pooled cases.
    if c.min_wt_cells == 0 { bail!("--min-wt-cells must be positive"); }
    let mut reference_profiles = BTreeMap::new();
    for (cluster, (_, wt_cells)) in &groups {
        if wt_cells.len() >= c.min_wt_cells {
            reference_profiles.insert(cluster.clone(), crate::wt_match::profile(c.data, wt_cells));
        }
    }
    let actual_groups = groups.clone();
    let mut reference_clusters = BTreeMap::new();
    let mut effective_groups = groups.clone();
    let mut matching = BufWriter::new(File::create(c.out.join("wt_reference_clusters.tsv"))?);
    writeln!(matching,"cluster\tcase_n\twithin_cluster_wt_n\twt_reference_cluster\twt_reference_n\tmethod\tpearson_log1p_mean")?;
    let mut unresolved = Vec::new();
    for (cluster, (case_cells, wt_cells)) in &groups {
        if case_cells.is_empty() { continue; }
        let (reference, method, similarity) = if wt_cells.len() >= c.min_wt_cells {
            (cluster.clone(), "same_cluster", Some(1.0))
        } else {
            let source = crate::wt_match::profile(c.data, case_cells);
            let candidate = groups.iter()
                .filter(|(other,(_,wt))| *other != cluster && wt.len() >= c.min_wt_cells)
                .filter_map(|(other,_)| {
                    let score = crate::wt_match::correlation(&source, &reference_profiles[other])?;
                    Some((other,score))
                })
                .max_by(|(a,sa),(b,sb)| sa.total_cmp(sb).then_with(|| b.cmp(a)));
            let Some((other,score)) = candidate else {
                unresolved.push(cluster.clone());
                writeln!(matching,"{cluster}\t{}\t{}\tNA\t0\tunresolved\tNA",
                    case_cells.len(), wt_cells.len())?;
                continue;
            };
            (other.clone(), "nearest_louvain", Some(score))
        };
        reference_clusters.insert(cluster.clone(), reference.clone());
        let reference_cells = &groups[&reference].1;
        effective_groups.get_mut(cluster).expect("cluster present").1 = reference_cells.clone();
        writeln!(matching,"{cluster}\t{}\t{}\t{reference}\t{}\t{method}\t{}",
            case_cells.len(),wt_cells.len(),reference_cells.len(),
            similarity.map(|x|x.to_string()).unwrap_or_else(||"NA".into()))?;
    }
    // All downstream STRING partner statistics, target effects and gene modules
    // use the resolved WT references, never the undersampled original WT cells.
    for cluster in &unresolved { effective_groups.remove(cluster); }
    if groups.values().any(|(cases, _)| !cases.is_empty()) && effective_groups.values().all(|(cases, _)| cases.is_empty()) {
        bail!("no populations have an adequate reference; see wt_reference_clusters.tsv");
    }
    let groups = effective_groups;
    for cluster in c.co_variant {
        if !groups.contains_key(cluster) { bail!("unknown or unresolved co-variant cluster {cluster}"); }
    }
    // STRING info file: protein_external_id, preferred_name, protein_size, annotation.
    let mut proteins = HashMap::<String,String>::new();
    let mut seed_ids = HashSet::new();
    for line in lines(c.string_info)?.lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        let mut parts = line.split_whitespace();
        let (Some(id), Some(name)) = (parts.next(),parts.next()) else { continue };
        proteins.insert(id.into(),name.into());
        if name.eq_ignore_ascii_case(c.gene) { seed_ids.insert(id.to_owned()); }
    }
    if seed_ids.is_empty() { bail!("{} not found in STRING protein info", c.gene); }
    // STRING IDs identify representative proteins. Resolve them through
    // protein.info preferred_name to gene-level associations; do not assign
    // these associations as experimentally established isoform interactions.
    // Both protein.links (3 columns) and protein.links.full (16 columns) are
    // supported. Never interpret `neighborhood` as `combined_score`.
    let mut partners = BTreeMap::<String,u32>::new();
    let mut score_column: Option<usize> = None;
    let mut seed_edges = 0usize;
    let mut passing_edges = 0usize;
    let mut unresolved = 0usize;
    for (line_number, line) in lines(c.string_links)?.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<_> = line.split_whitespace().collect();
        if score_column.is_none() {
            if parts.first() == Some(&"protein1") && parts.get(1) == Some(&"protein2") {
                score_column = Some(parts.iter().position(|v| *v == "combined_score")
                    .with_context(|| format!("STRING links header missing combined_score: {}", c.string_links.display()))?);
                continue;
            }
            // Headerless files are supported only for the standard 3-column layout.
            if parts.len() != 3 {
                bail!("STRING links {}:{} has {} columns without a header; cannot identify combined_score",
                    c.string_links.display(), line_number + 1, parts.len());
            }
            score_column = Some(2);
        }
        let column = score_column.expect("initialized above");
        if parts.len() <= column { bail!("malformed STRING links line {}", line_number + 1); }
        let other = if seed_ids.contains(parts[0]) { parts[1] }
            else if seed_ids.contains(parts[1]) { parts[0] }
            else { continue; };
        seed_edges += 1;
        let score: u32 = parts[column].parse().with_context(||
            format!("invalid combined_score on STRING links line {}", line_number + 1))?;
        if score < c.min_score { continue; }
        passing_edges += 1;
        if let Some(symbol) = proteins.get(other) {
            if !symbol.eq_ignore_ascii_case(c.gene) {
                partners.entry(symbol.clone()).and_modify(|s| *s=(*s).max(score)).or_insert(score);
            }
        } else { unresolved += 1; }
    }
    eprintln!("[traconmech] STRING seed {}: {} IDs; {} incident edges; {} passing score >= {}; {} unmapped partner IDs; {} distinct partner genes",
        c.gene, seed_ids.len(), seed_edges, passing_edges, c.min_score, unresolved, partners.len());
    if seed_edges == 0 {
        bail!("no STRING interactions for {} (STRING IDs: {:?}); check species and links file",
            c.gene, seed_ids);
    }
    if partners.is_empty() {
        bail!("no STRING partner genes for {} survived score threshold {} ({} incident edges)",
            c.gene, c.min_score, seed_edges);
    }
    // MSigDB C3 target collections: use gene-set name as a *candidate* regulator
    // mapping. No generic TF identity is inferred from STRING alone.
    let mut partner_out = BufWriter::new(File::create(c.out.join("partners.tsv"))?);
    writeln!(partner_out,"partner\tstring_score\tcluster\tcase_mean\tcase_sd\tcase_detected\tcontrol_mean\tcontrol_sd\tcontrol_detected\tregulatory_sets")?;
    let mut targets: BTreeMap<String,BTreeSet<String>> = BTreeMap::new();
    let mut sets_out = BufWriter::new(File::create(c.out.join("regulatory_sets.tsv"))?);
    writeln!(sets_out,"regulator\tset\tsource\tgenes")?;
    let regulatory_sets = ommverse::msigdb::all_sets(c.index)?;
    eprintln!("[traconmech] loaded {} MSigDB sets from sidecar", regulatory_sets.len());
    for (partner,score) in &partners {
        let mut matching = Vec::new();
        // The MSigDB reverse lookup is by MEMBER, not by regulator. We instead
        // load all sets via the public collection iterator, matching C3 names.
        for set in &regulatory_sets {
            if !set.source.to_ascii_lowercase().starts_with("c3.") { continue; }
            let upper = set.name.to_ascii_uppercase();
            let symbol = partner.to_ascii_uppercase();
            if upper.split('_').any(|token| token == symbol) {
                matching.push(set);
            }
        }
        for set in &matching {
            writeln!(sets_out,"{}\t{}\t{}\t{}",partner,set.name,set.source,set.genes.len())?;
            for gene in &set.genes { targets.entry(gene.clone()).or_default().insert(set.name.clone()); }
        }
        for (cluster,(case_cells,control_cells)) in &groups {
            if case_cells.is_empty() { continue; }
            let (cm,cs,cd)=stats(c.data,partner,case_cells);
            let (wm,ws,wd)=stats(c.data,partner,control_cells);
            writeln!(partner_out,"{partner}\t{score}\t{cluster}\t{cm}\t{cs}\t{cd}\t{wm}\t{ws}\t{wd}\t{}",matching.len())?;
        }
    }
    // The requested order is a biological trajectory, not BTreeMap key order.
    let ordered: Vec<String> = if c.co_variant.is_empty() {
        groups.iter().filter(|(_, (cases, _))| !cases.is_empty()).map(|(cluster, _)| cluster.clone()).collect()
    } else {
        c.co_variant.to_vec()
    };
    // Native AnnData categorical levels are authoritative for plotting order.
    // For plain string labels, retain the existing deterministic ordering.
    let mut ordered = ordered;
    if c.co_variant.is_empty() {
        if let Some(levels) = c.data.annotations.categorical_levels(c.cluster) {
            let rank: HashMap<&str, usize> = levels.iter().enumerate()
                .map(|(i, level)| (level.as_str(), i)).collect();
            ordered.sort_by_key(|label| rank.get(label.as_str()).copied().unwrap_or(usize::MAX));
        } else if c.cluster == "cell_type" {
            // Numeric developmental_order can also provide ordering for cell_type.
            if let (Some(names), Some(orders)) = (
                c.data.annotations.column(c.cluster),
                c.data.annotations.column("developmental_order"),
            ) {
                let mut rank: HashMap<&str, f64> = HashMap::new();
                for (name, order) in names.iter().zip(orders) {
                    if let Ok(value) = order.parse::<f64>() {
                        if value.is_finite() {
                            rank.entry(name.as_str()).and_modify(|v| *v = v.min(value))
                                .or_insert(value);
                        }
                    }
                }
                ordered.sort_by(|a, b| {
                    let x = rank.get(a.as_str()).copied().unwrap_or(f64::INFINITY);
                    let y = rank.get(b.as_str()).copied().unwrap_or(f64::INFINITY);
                    x.total_cmp(&y).then_with(|| a.cmp(b))
                });
            }
        }
    }
    let mut signatures = BTreeMap::<String,Vec<f64>>::new();
    for gene in targets.keys() {
        let mut effects = Vec::with_capacity(ordered.len());
        for cluster in &ordered {
            let (case_cells,control_cells) = &groups[cluster];
            let (cm,_,_) = stats(c.data,gene,case_cells);
            let (wm,_,_) = stats(c.data,gene,control_cells);
            effects.push(((cm+0.01)/(wm+0.01)).log2());
        }
        signatures.insert(gene.clone(),effects);
    }
    let modules = if c.effect_cluster {
        effect_modules(&signatures)
    } else {
        BTreeMap::new()
    };
    let mut effect_out = BufWriter::new(File::create(c.out.join("effect_clusters.tsv"))?);
    writeln!(effect_out,"gene\teffect_cluster\tmembership_sets\tcluster\tcase_n\tcontrol_n\tcase_mean\tcontrol_mean\tlog2_ratio\tcase_detected\tcontrol_detected\tcase_sd\tcontrol_sd\tstability_weight")?;
    for (gene, member_sets) in &targets {
        let id = modules.get(gene).copied().unwrap_or(0);
        let member = member_sets.iter().cloned().collect::<Vec<_>>().join(";");
        for cluster in &ordered {
            let (case_cells,control_cells) = &groups[cluster];
            let (cm,cs,cd) = stats(c.data,gene,case_cells);
            let (wm,ws,wd) = stats(c.data,gene,control_cells);
            let ratio = ((cm+0.01)/(wm+0.01)).log2();
            let weight = stability_weight(cm, cs, wm, ws, cd, wd);
            writeln!(effect_out,"{gene}\t{id}\t{member}\t{cluster}\t{}\t{}\t{cm}\t{wm}\t{ratio}\t{cd}\t{wd}\t{cs}\t{ws}\t{weight}",case_cells.len(),control_cells.len())?;
        }
    }
    if !c.co_variant.is_empty() {
        let mut out = BufWriter::new(File::create(c.out.join("co_variant_effects.tsv"))?);
        writeln!(out,"gene\t{}",ordered.join("\t"))?;
        for (gene,effects) in &signatures {
            writeln!(out,"{}\t{}",gene,effects.iter().map(ToString::to_string).collect::<Vec<_>>().join("\t"))?;
        }
    }
    if c.effect_cluster {
        // One row per candidate gene; exclusive, deterministic module membership.
        // Use the actual assignments, never infer them from rounded effect summaries.
        let module_ids: Vec<usize> = modules.values().copied()
            .collect::<BTreeSet<_>>().into_iter().collect();
        let mut out = BufWriter::new(File::create(c.out.join("module_membership.tsv"))?);
        write!(out, "gene")?;
        for id in &module_ids { write!(out, "\tmodule_{id}")?; }
        writeln!(out)?;
        for (gene, assigned) in &modules {
            write!(out, "{gene}")?;
            for id in &module_ids { write!(out, "\t{}", u8::from(id == assigned))?; }
            writeln!(out)?;
        }
        out.flush()?;

        let mut out = BufWriter::new(File::create(c.out.join("effect_modules.tsv"))?);
        writeln!(out,"module\tgenes\t{}",ordered.join("\t"))?;
        let mut module_genes: BTreeMap<usize,Vec<&String>> = BTreeMap::new();
        for (gene,module) in &modules { module_genes.entry(*module).or_default().push(gene); }
        for (module,genes) in &module_genes {
            let mut means = vec![0.0;ordered.len()];
            for gene in genes {
                for (sum,value) in means.iter_mut().zip(&signatures[*gene]) { *sum += value; }
            }
            for sum in &mut means { *sum /= genes.len() as f64; }
            writeln!(out,"{}\t{}\t{}",module,genes.len(),means.iter().map(ToString::to_string).collect::<Vec<_>>().join("\t"))?;
        }
        // Per-population module evidence and complete ranked stable-gene lists.
        let mut summary = BufWriter::new(File::create(c.out.join("module_measurements.tsv"))?);
        writeln!(summary, "module\tcluster\tgenes\tcase_n\tcontrol_n\tcase_mean\tcontrol_mean\tmean_effect\tweighted_effect\teffective_genes\tmean_gene_effect_sd")?;
        let mut stable = BufWriter::new(File::create(c.out.join("stable_module_genes.tsv"))?);
        writeln!(stable, "module\tcluster\trank\tgene\tcase_n\tcontrol_n\tcase_mean\tcase_sd\tcontrol_mean\tcontrol_sd\tlog2_ratio\tcase_detected\tcontrol_detected\tstability_weight\tweighted_abs_effect")?;
        for (module, genes) in &module_genes {
            for cluster in &ordered {
                let (case_cells, control_cells) = &groups[cluster];
                let mut rows = Vec::with_capacity(genes.len());
                for gene in genes {
                    let (cm, cs, cd) = stats(c.data, gene, case_cells);
                    let (wm, ws, wd) = stats(c.data, gene, control_cells);
                    let effect = ((cm + 0.01) / (wm + 0.01)).log2();
                    let weight = stability_weight(cm, cs, wm, ws, cd, wd);
                    rows.push((gene.as_str(), cm, cs, cd, wm, ws, wd, effect, weight));
                }
                rows.sort_by(|a,b| (b.7.abs()*b.8).total_cmp(&(a.7.abs()*a.8))
                    .then_with(|| a.0.cmp(b.0)));
                let n = rows.len() as f64;
                let case_mean = rows.iter().map(|r| r.1).sum::<f64>() / n;
                let control_mean = rows.iter().map(|r| r.4).sum::<f64>() / n;
                let effect_mean = rows.iter().map(|r| r.7).sum::<f64>() / n;
                let effect_sd = (rows.iter().map(|r| (r.7-effect_mean).powi(2)).sum::<f64>() / n).sqrt();
                let weight_sum = rows.iter().map(|r| r.8).sum::<f64>();
                let weight_sq = rows.iter().map(|r| r.8*r.8).sum::<f64>();
                let weighted = if weight_sum > 0.0 { rows.iter().map(|r| r.7*r.8).sum::<f64>() / weight_sum } else { 0.0 };
                let effective = if weight_sq > 0.0 { weight_sum*weight_sum/weight_sq } else { 0.0 };
                writeln!(summary, "{module}\t{cluster}\t{}\t{}\t{}\t{case_mean}\t{control_mean}\t{effect_mean}\t{weighted}\t{effective}\t{effect_sd}", rows.len(),case_cells.len(),control_cells.len())?;
                for (i,r) in rows.iter().enumerate() {
                    writeln!(stable, "{module}\t{cluster}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i+1,r.0,case_cells.len(),control_cells.len(),r.1,r.2,r.4,r.5,r.7,r.3,r.6,r.8,r.8*r.7.abs())?;
                }
            }
        }
        summary.flush()?;
        stable.flush()?;
        // Annotation must precede SVG generation so both plot variants have
        // cell-type labels. Reuse saved AI results when rerunning the analysis.
        let ai_order = if c.ai_development {
            eprintln!("[traconmech] developmental AI: loading cached annotation or contacting Ollama ({})", c.ai_model);
            Some(crate::ai::developmental_order_cached(c.out, c.ai_model, c.ai_url, &ordered, c.ai_primer)?)
        } else { None };
        crate::module_plots::write_module_plots(c.out, c.data, &modules, &ordered, &groups, &actual_groups, &reference_clusters, c.cases, c.controls, c.label_angle, c.count_label_layout)?;
        if let Some(ai_order) = ai_order {
            crate::module_plots::write_module_plots_named(c.out, c.data, &modules,
                &ai_order, &groups, &actual_groups, &reference_clusters,
                c.cases, c.controls, "_developmental", c.label_angle, c.count_label_layout)?;
        }

        eprintln!("[traconmech] effect clustering: {} modules (0 = flat trajectory); order: {}",
            module_genes.len(),ordered.join(" -> "));
    }
    // Explain the hypothesis, evidence chain, and interpretation for this run.
    let mut readme = BufWriter::new(File::create(c.out.join("README.md"))?);
    writeln!(readme, "# TraConMech investigation: {}\n", c.gene)?;
    writeln!(readme, "## Hypothesis\n")?;
    writeln!(readme, "A perturbation associated with `{}` may change expression of genes annotated to its STRING-associated partners. Population-specific expression patterns can suggest mechanisms worth testing, but neither STRING association nor gene-set membership establishes causal regulation.\n", c.gene)?;
    writeln!(readme, "## Experimental comparison\n")?;
    writeln!(readme, "Read `experiment.yaml` for the exact reference and experimental labels. WT/control (blue) and EXP/case (red) name comparison roles, not necessarily literal genotypes. Population order: `{}`. Minimum STRING score: `{}`. Effect clustering enabled: `{}`.\n", ordered.join(" -> "), c.min_score, c.effect_cluster)?;
    writeln!(readme, "## Method\n")?;
    writeln!(readme, "1. Find direct STRING partners of `{}` above the configured combined-score threshold. STRING measures association, not physical binding or direction of regulation.\n", c.gene)?;
    writeln!(readme, "2. Match partners to MSigDB regulatory gene sets and collect their genes. One gene may have multiple partner/set paths; these are candidate annotations, not independent causal proofs.\n")?;
    writeln!(readme, "3. For each gene and cell population, measure expression separately in EXP and reference cells and calculate `log2((case_mean + 0.01)/(control_mean + 0.01))`. Positive values indicate higher expression in EXP. Interpret alongside cell counts and detection fractions.\n")?;
    if c.effect_cluster {
        writeln!(readme, "4. Mean-center and L2-normalize each gene's vector of population-specific log2 effects. Deterministic spherical k-means groups similar effect *shapes* into modules 1–8. Genes with fewer than two populations or a flat centered vector enter module 0; it is retained for auditing, not plotted. A flat effect is not necessarily zero expression. Module numbers do not identify regulators.\n")?;
    }
    writeln!(readme, "## Identify the cell populations\n")?;
    writeln!(readme, "Start with `cluster_markers/top20_up_down.tsv`, generated by the existing sc-analysis cluster-versus-rest marker statistics from the original cluster assignments. It contains up to 20 enriched (`up`) and 20 depleted (`down`) genes per population, ranked by adjusted p-value and absolute log2 fold change. `cluster_markers/cluster_labels.tsv` maps sc-analysis identifiers such as `cluster_000` to the original `{}` labels used in the trajectory plots. These are marker candidates, not automatic cell-type annotations. Marker statistics use all cells regardless of case/control label. The `cluster_markers/cluster_XXX.tsv` files retain full ranked results.\n", c.cluster)?;
    writeln!(readme, "## How to explore the results\n")?;
    if c.effect_cluster {
        writeln!(readme, "**Start with `module_N_trajectory.svg` (N=1–8).** Blue and red lines show the median of per-gene jointly normalized log1p population means, with gene-wise 25th–75th percentile bands. Normalization uses both conditions together, preserving their separation. Bands describe variation across genes, not uncertainty across cells. The underlying plot values are in `module_trajectories.tsv`. Cluster labels are categorical: connected lines do not establish a developmental sequence.\n")?;
        writeln!(readme, "**Check signed effects** in `effect_modules.tsv` and `effect_clusters.tsv`: a module groups shapes, so opposing positive and negative gene effects can cancel in a summary.\n")?;
        writeln!(readme, "**Inspect genes** in `stable_module_genes.tsv`. `rank` is within each population, not a global module rank. Examine cell counts, expression means, detection fractions and stability weights.\n")?;
        writeln!(readme, "**Trace candidate mechanisms** from `module_regulatory_sets.tsv` to `gene_regulatory_provenance.tsv` and then the STRING partner of `{}`. Gene-set overlaps are descriptive, not enrichment p-values or FDR; large sets can dominate by size.\n", c.gene)?;
    } else {
        writeln!(readme, "Start with `effect_clusters.tsv` for measured population effects, then inspect `gene_regulatory_provenance.tsv` for candidate partner/set paths. No expression modules were requested.\n")?;
    }
    writeln!(readme, "## Output files\n")?;
    writeln!(readme, "| File | Purpose |\n|---|---|")?;
    for (file, description) in [
        ("cluster_markers/top20_up_down.tsv", "sc-analysis top 20 enriched and depleted genes per existing cluster"),
        ("cluster_markers/cluster_labels.tsv", "Map sc-analysis cluster identifiers to original population labels"),
        ("cluster_markers/marker_summary.tsv", "Compact sc-analysis cluster marker signatures"),
        ("experiment.yaml", "Recorded experiment labels and comparison"),
        ("partners.tsv", "Direct STRING-associated partners"),
        ("regulatory_sets.tsv", "Partner-matched MSigDB gene sets"),
        ("gene_regulatory_provenance.tsv", "Candidate seed-to-partner-to-set-to-gene links"),
        ("effect_clusters.tsv", "Measured gene effects per population"),
        ("module_membership.tsv", "Gene-to-expression-module assignments"),
        ("effect_modules.tsv", "Module-level signed effect summaries"),
        ("module_measurements.tsv", "Module and population measurements"),
        ("stable_module_genes.tsv", "Population-specific ranked genes"),
        ("module_regulatory_sets.tsv", "Module overlap with partner-associated sets"),
        ("module_trajectories.tsv", "Values behind WT/EXP trajectory plots"),
    ] {
        if c.out.join(file).is_file() { writeln!(readme, "| `{file}` | {description} |")?; }
    }
    if c.effect_cluster {
        writeln!(readme, "| `module_N_trajectory.svg` | Blue reference / red experimental trajectory plots for modules 1–8 |")?;
    }
    if !c.co_variant.is_empty() && c.out.join("co_variant_effects.tsv").is_file() {
        writeln!(readme, "| `co_variant_effects.tsv` | Effects in the requested population order |")?;
    }
    for (file, description) in [
        ("wt_reference_clusters.tsv", "Per-population reference selection and fallback"),
        ("closest_wt_clusters.tsv", "Reference candidate ranking"),
        ("dba_vs_closest_wt.tsv", "Case vs matched reference comparisons"),
        ("cluster_counts.tsv", "Cell counts by population and condition"),
        ("cluster_expression.tsv", "Population expression profiles"),
        ("cluster_correlations.tsv", "Population expression correlations"),
    ] {
        if c.out.join(file).is_file() { writeln!(readme, "| `{file}` | {description} |")?; }
    }
    writeln!(readme, "\n## Interpretation and limitations\n")?;
    writeln!(readme, "Expression differences are observations, while STRING/MSigDB paths are hypotheses. Module assignments are based on relative effect trajectories, not shared regulator identity. Module 0 is omitted from plots but preserved in data for auditing. Gene-set overlap is not enrichment, and separately inferred cell clusters are not automatically equivalent across experiments.\n")?;
    readme.flush()?;
    eprintln!("[traconmech] {} STRING partners; {} MSigDB target genes; results in {}",partners.len(),targets.len(),c.out.display());
    Ok(())
}
