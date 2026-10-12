//! Reusable, reference-independent population expression profiles.
//! Uses the same feature-by-cell matrix as the rest of TraConMech.
use anyhow::{bail, Context, Result};
use sc_analysis::SingleCellData;
use crate::experiment::require_annotation;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn mean_and_detection(sums: &[f64], detected: &[usize], n: usize, g: usize) -> (f64,f64) {
    if n == 0 { (0.0,0.0) } else { (sums[g]/n as f64, detected[g] as f64/n as f64) }
}
fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len()!=b.len() || a.len()<2 {return None}
    let n=a.len() as f64;
    let (ma,mb)=(a.iter().sum::<f64>()/n,b.iter().sum::<f64>()/n);
    let (mut aa,mut bb,mut ab)=(0.0,0.0,0.0);
    for (&x,&y) in a.iter().zip(b) {
        let (dx,dy)=(x-ma,y-mb);
        aa+=dx*dx;bb+=dy*dy;ab+=dx*dy;
    }
    if aa<=0.0 || bb<=0.0 {None} else {Some(ab/(aa*bb).sqrt())}
}
fn cell_groups(data:&SingleCellData, cluster:&str, condition:Option<&str>, case:&str, control:&str)
    -> Result<BTreeMap<(String,String),Vec<usize>>> {
    let labels=require_annotation(data, cluster, "population")?;
    let conditions=condition.map(|name|require_annotation(data, name, "condition")).transpose()?;
    let cases:BTreeSet<&str>=case.split(',').filter(|x|!x.is_empty()).collect();
    let controls:BTreeSet<&str>=control.split(',').filter(|x|!x.is_empty()).collect();
    if !cases.is_disjoint(&controls) {bail!("case and control labels overlap");}
    let mut groups=BTreeMap::new();
    for (i,cl) in labels.iter().enumerate() {
        groups.entry((cl.clone(),"all".to_owned())).or_insert_with(Vec::new).push(i);
        if let Some(values)=conditions {
            let state=if cases.contains(values[i].as_str()) {"case"}
                else if controls.contains(values[i].as_str()) {"control"}
                else {"other"};
            groups.entry((cl.clone(),state.to_owned())).or_insert_with(Vec::new).push(i);
        }
    }
    Ok(groups)
}
fn load_modules(path:&Path)->Result<BTreeMap<String,usize>> {
    let reader=BufReader::new(File::open(path).with_context(||format!("opening {}",path.display()))?);
    let mut lines=reader.lines();
    let header=lines.next().context("empty module assignment file")??;
    let columns:Vec<&str>=header.split('\t').collect();
    let gene_col=columns.iter().position(|s|*s=="gene").context("module file missing gene column")?;
    let module_col=columns.iter().position(|s|*s=="effect_cluster").context("module file missing effect_cluster column")?;
    let mut result=BTreeMap::new();
    for line in lines {
        let line=line?;
        let fields:Vec<&str>=line.split('\t').collect();
        if let (Some(&gene),Some(&id))=(fields.get(gene_col),fields.get(module_col)) {
            let id:usize=id.parse().with_context(||format!("bad module ID for {gene}"))?;
            if let Some(previous)=result.insert(gene.to_owned(),id) {
                if previous!=id {bail!("conflicting module assignments for {gene}");}
            }
        }
    }
    Ok(result)
}
pub fn run(data:&SingleCellData,out:&Path,cluster:&str,condition:Option<&str>,
           case:&str,control:&str,reference_modules:Option<&Path>)->Result<()> {
    fs::create_dir_all(out)?;
    let groups=cell_groups(data,cluster,condition,case,control)?;
    let n_genes=data.features.len();
    let mut means=BTreeMap::<(String,String),Vec<f64>>::new();
    let mut detections=BTreeMap::<(String,String),Vec<f64>>::new();
    let mut counts=BufWriter::new(File::create(out.join("cluster_counts.tsv"))?);
    writeln!(counts,"cluster\tcondition\tn_cells")?;
    // Single CSR traversal per group: do not repeatedly scan every cell for every gene.
    for (key,cells) in &groups {
        let mut selected=vec![false;data.cells.len()];
        for &i in cells {selected[i]=true;}
        let mut sums=vec![0.0;n_genes];
        let mut hits=vec![0usize;n_genes];
        for (g,row) in data.matrix.outer_iterator().enumerate() {
            for (cell,value) in row.iter() {
                if selected[cell] {
                    sums[g]+=f64::from(*value);
                    if *value!=0.0 {hits[g]+=1;}
                }
            }
        }
        let mut mean=vec![0.0;n_genes];
        let mut detection=vec![0.0;n_genes];
        for g in 0..n_genes {
            (mean[g],detection[g])=mean_and_detection(&sums,&hits,cells.len(),g);
        }
        writeln!(counts,"{}\t{}\t{}",key.0,key.1,cells.len())?;
        means.insert(key.clone(),mean);
        detections.insert(key.clone(),detection);
    }
    let mut expression=BufWriter::new(File::create(out.join("cluster_expression.tsv"))?);
    writeln!(expression,"gene\tcluster\tcondition\tn_cells\tmean_expression\tdetection_fraction")?;
    for (key,mean) in &means {
        let det=&detections[key];
        let n=groups[key].len();
        for (g,gene) in data.features.iter().enumerate() {
            writeln!(expression,"{gene}\t{}\t{}\t{n}\t{}\t{}",key.0,key.1,mean[g],det[g])?;
        }
    }
    // Correlations are calculated on log1p(mean expression), all measured genes.
    // The same transform and gene universe are used for every within-dataset pair.
    let transformed:BTreeMap<_,_>=means.iter().map(|(k,v)|
        (k.clone(),v.iter().map(|x|x.max(0.0).ln_1p()).collect::<Vec<_>>())).collect();
    let mut corr=BufWriter::new(File::create(out.join("cluster_correlations.tsv"))?);
    writeln!(corr,"cluster_a\tcondition_a\tn_a\tcluster_b\tcondition_b\tn_b\tpearson_log1p_mean\tgenes_compared")?;
    for (a,va) in &transformed {
        for (b,vb) in &transformed {
            let value=pearson(va,vb).map(|x|x.to_string()).unwrap_or_else(||"NA".into());
            writeln!(corr,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",a.0,a.1,groups[a].len(),b.0,b.1,groups[b].len(),value,n_genes)?;
        }
    }
    // Reuse fixed module definitions from an earlier DBA analysis on Johan's data.
    // Never recluster: identical gene-to-module assignments in every dataset.
    if let Some(path)=reference_modules {
        let assignments=load_modules(path)?;
        let lookup:HashMap<&str,usize>=data.features.iter().enumerate().map(|(i,g)|(g.as_str(),i)).collect();
        let mut by_module=BTreeMap::<usize,Vec<usize>>::new();
        let mut missing=BTreeMap::<usize,usize>::new();
        for (gene,module) in &assignments {
            if let Some(&g)=lookup.get(gene.as_str()) {
                by_module.entry(*module).or_default().push(g);
            } else {*missing.entry(*module).or_default()+=1;}
        }
        let mut scores=BufWriter::new(File::create(out.join("projected_modules.tsv"))?);
        writeln!(scores,"module\tcluster\tcondition\tn_cells\tn_reference_genes\tn_measured_genes\tn_missing_genes\tmean_expression\tmean_detection_fraction")?;
        let modules:BTreeSet<usize>=assignments.values().copied().collect();
        for module in modules {
            let genes=by_module.get(&module).map(Vec::as_slice).unwrap_or(&[]);
            for (key,mean) in &means {
                let det=&detections[key];
                let m=if genes.is_empty(){0.0}else{genes.iter().map(|&g|mean[g]).sum::<f64>()/genes.len() as f64};
                let d=if genes.is_empty(){0.0}else{genes.iter().map(|&g|det[g]).sum::<f64>()/genes.len() as f64};
                writeln!(scores,"{module}\t{}\t{}\t{}\t{}\t{}\t{}\t{m}\t{d}",
                    key.0,key.1,groups[key].len(),genes.len()+missing.get(&module).copied().unwrap_or(0),
                    genes.len(),missing.get(&module).copied().unwrap_or(0))?;
            }
        }
    }
    eprintln!("[traconmech] profiled {} populations/condition groups and {} genes in {}",
        groups.len(),n_genes,out.display());
    Ok(())
}
