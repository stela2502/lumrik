//! Deterministic Seurat AddModuleScore-style expression-matched background scoring.
use anyhow::{ensure, Result};
use sc_analysis::SingleCellData;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};

pub struct ModuleScores {
    pub ids: Vec<usize>,
    pub scores: Vec<Vec<f64>>, // cell x module
    pub module_means: Vec<Vec<f64>>,
    pub coverage: BTreeMap<usize, (usize, usize)>,
}

fn next_u64(state: &mut u64) -> u64 {
    *state ^= *state >> 12;
    *state ^= *state << 25;
    *state ^= *state >> 27;
    state.wrapping_mul(0x2545F4914F6CDD1D)
}

pub fn score(data: &SingleCellData, modules: &BTreeMap<String, usize>) -> Result<ModuleScores> {
    const NBIN: usize = 24;
    const CTRL: usize = 100;
    let n = data.cells.len();
    let ng = data.features.len();
    ensure!(n > 0 && ng > 0, "empty expression matrix");
    let lookup: HashMap<&str, usize> = data.features.iter().enumerate().map(|(i,g)| (g.as_str(),i)).collect();
    let mut ids: Vec<usize> = modules.values().copied().filter(|&id| id != 0).collect();
    ids.sort_unstable(); ids.dedup();
    ensure!(!ids.is_empty(), "no nonzero modules");
    let mut coverage = BTreeMap::new();
    let mut members = Vec::new();
    for &id in &ids {
        let genes: Vec<_> = modules.iter().filter(|(_, mid)| **mid == id).map(|(g,_)| g).collect();
        let indices: Vec<usize> = genes.iter().filter_map(|g| lookup.get(g.as_str()).copied()).collect();
        coverage.insert(id, (genes.len(), indices.len()));
        ensure!(!indices.is_empty(), "module {id} has no genes in dataset");
        members.push(indices);
    }
    // Mean of supplied log1p values, including zeros. Do not transform again.
    let means: Vec<f64> = data.matrix.outer_iterator().collect::<Vec<_>>()
        .into_par_iter()
        .map(|row| row.iter().map(|(_, v)| f64::from(*v)).sum::<f64>() / n as f64)
        .collect();
    // Rank genes into equal-frequency expression bins. Deterministic tie order.
    let mut ranked: Vec<usize> = (0..ng).collect();
    ranked.sort_by(|&a,&b| means[a].total_cmp(&means[b]).then(a.cmp(&b)));
    let mut bin_of = vec![0usize;ng];
    let mut bins = vec![Vec::new(); NBIN];
    for (rank,&g) in ranked.iter().enumerate() {
        let bin = rank * NBIN / ng;
        bin_of[g] = bin;
        bins[bin].push(g);
    }
    let mut selected = Vec::new();
    for (j,genes) in members.iter().enumerate() {
        let mut state = 0x9E3779B97F4A7C15u64 ^ (ids[j] as u64).wrapping_mul(0xD6E8FEB86659FD93);
        let mut controls = HashSet::new();
        for &g in genes {
            let pool = &bins[bin_of[g]];
            let take = CTRL.min(pool.len());
            // Uniform sample without replacement for each gene, then deduplicate
            // across module genes, matching AddModuleScore's control union.
            let mut available = pool.clone();
            for k in 0..take {
                let pick = k + (next_u64(&mut state) as usize % (available.len()-k));
                available.swap(k,pick);
                controls.insert(available[k]);
            }
        }
        ensure!(!controls.is_empty(), "no controls for module {}", ids[j]);
        selected.push(controls);
    }
    let mut module_means = vec![vec![0.0;ids.len()];n];
    let mut background = vec![vec![0.0;ids.len()];n];
    let mut gene_module = vec![Vec::new();ng];
    let mut gene_control = vec![Vec::new();ng];
    for (j,genes) in members.iter().enumerate() {
        for &g in genes { gene_module[g].push(j); }
        for &g in &selected[j] { gene_control[g].push(j); }
    }
    for (g,row) in data.matrix.outer_iterator().enumerate() {
        if gene_module[g].is_empty() && gene_control[g].is_empty() { continue; }
        for (c,v) in row.iter() {
            let x = f64::from(*v);
            for &j in &gene_module[g] { module_means[c][j] += x; }
            for &j in &gene_control[g] { background[c][j] += x; }
        }
    }
    let scores: Vec<Vec<f64>> = module_means.par_iter_mut()
        .zip(background.par_iter())
        .map(|(cell_means, cell_background)| {
            cell_means.iter_mut().enumerate().map(|(j, mean)| {
                *mean /= members[j].len() as f64;
                *mean - cell_background[j] / selected[j].len() as f64
            }).collect()
        }).collect();
    Ok(ModuleScores { ids, scores, module_means, coverage })
}
