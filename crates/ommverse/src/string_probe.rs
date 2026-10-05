use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use gtf_splice_index::GeneId;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::Ommverse;

#[derive(Debug, Clone)]
struct Node {
    protein_id: usize,
    gene_id: GeneId,
    features: Vec<String>,
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().is_some_and(|x| x == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn embedded_ensp(string_id: &str) -> &str {
    string_id.split_once('.').map(|(_, id)| id).unwrap_or(string_id)
}

fn feature_vector(omm: &Ommverse, protein_id: usize) -> Vec<String> {
    let mut features = omm.proteins[protein_id]
        .features
        .iter()
        .filter(|feature| feature.source_db.eq_ignore_ascii_case("InterPro"))
        .map(|feature| feature.label.clone())
        .collect::<Vec<_>>();
    features.sort();
    features.dedup();
    features
}

fn protein_gene(omm: &Ommverse, protein_id: usize) -> Option<GeneId> {
    let protein = omm.proteins.get(protein_id)?;
    let mut genes = protein
        .transcript_ids
        .iter()
        .filter_map(|&tx_id| omm.splice.transcripts.get(tx_id).map(|tx| tx.gene_id))
        .collect::<Vec<_>>();
    genes.sort_unstable();
    genes.dedup();
    (genes.len() == 1).then_some(genes[0])
}

pub fn probe_string(index: &Path, info: &Path, links: &Path, out: &Path) -> Result<()> {
    eprintln!("[ommverse] loading Ommverse reference for STRING probe...");
    let omm = Ommverse::load(index)?;
    eprintln!("[ommverse] resolving STRING proteins through Ommverse ENSP aliases...");

    let mut nodes: HashMap<String, Node> = HashMap::new();
    let mut info_records = 0usize;
    let mut resolved = 0usize;
    let mut unresolved = 0usize;
    let mut multi_gene = 0usize;

    for line in open_text(info)?.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        info_records += 1;
        let Some(string_id) = line.split_whitespace().next() else { continue };
        let ensp = embedded_ensp(string_id);
        let Some(protein_id) = omm.protein_id(ensp) else {
            unresolved += 1;
            continue;
        };
        let Some(gene_id) = protein_gene(&omm, protein_id) else {
            multi_gene += 1;
            continue;
        };
        nodes.insert(
            string_id.to_owned(),
            Node { protein_id, gene_id, features: feature_vector(&omm, protein_id) },
        );
        resolved += 1;
    }

    let mut gene_vectors: HashMap<GeneId, BTreeSet<Vec<String>>> = HashMap::new();
    for node in nodes.values() {
        gene_vectors.entry(node.gene_id).or_default().insert(node.features.clone());
    }
    let candidate_genes = gene_vectors
        .iter()
        .filter_map(|(&gene, vectors)| (vectors.len() > 1).then_some(gene))
        .collect::<HashSet<_>>();
    let candidate_nodes = nodes
        .iter()
        .filter_map(|(id, node)| candidate_genes.contains(&node.gene_id).then_some(id.clone()))
        .collect::<HashSet<_>>();

    eprintln!(
        "[ommverse] STRING nodes: {info_records}; resolved: {resolved}; unresolved: {unresolved}; ambiguous-gene: {multi_gene}"
    );
    eprintln!(
        "[ommverse] genes represented: {}; genes with >1 InterPro feature vector: {}",
        gene_vectors.len(), candidate_genes.len()
    );
    eprintln!("[ommverse] streaming STRING links; retaining partner sets only for candidate genes...");

    let mut partners: HashMap<String, BTreeSet<GeneId>> = HashMap::new();
    let mut link_records = 0usize;
    let mut resolved_links = 0usize;
    for line in open_text(links)?.lines() {
        let line = line?;
        if line.starts_with("protein1 ") || line.trim().is_empty() {
            continue;
        }
        link_records += 1;
        let mut fields = line.split_whitespace();
        let Some(a) = fields.next() else { continue };
        let Some(b) = fields.next() else { continue };
        let (Some(na), Some(nb)) = (nodes.get(a), nodes.get(b)) else { continue };
        resolved_links += 1;
        if candidate_nodes.contains(a) {
            partners.entry(a.to_owned()).or_default().insert(nb.gene_id);
        }
        if candidate_nodes.contains(b) {
            partners.entry(b.to_owned()).or_default().insert(na.gene_id);
        }
    }

    let mut by_gene: BTreeMap<GeneId, Vec<(&str, &Node, BTreeSet<GeneId>)>> = BTreeMap::new();
    for (string_id, node) in &nodes {
        if candidate_genes.contains(&node.gene_id) {
            by_gene.entry(node.gene_id).or_default().push((
                string_id.as_str(), node, partners.get(string_id).cloned().unwrap_or_default(),
            ));
        }
    }

    let mut genes_with_multiple_link_sets = 0usize;
    let mut writer = BufWriter::new(File::create(out).with_context(|| format!("creating {}", out.display()))?);
    writeln!(writer, "gene_id\tgene\tn_feature_vectors\tn_link_sets\tstring_id\tprotein_id\ttranscript_ids\tinterpro_features\tn_partner_genes\tpartner_gene_ids")?;
    for (gene_id, rows) in &by_gene {
        let link_sets = rows.iter().map(|(_, _, set)| set.clone()).collect::<BTreeSet<_>>();
        if link_sets.len() > 1 { genes_with_multiple_link_sets += 1; }
        let gene = &omm.splice.genes[*gene_id];
        let gene_name = gene.primary_name().unwrap_or(".");
        let n_feature_vectors = gene_vectors.get(gene_id).map_or(0, BTreeSet::len);
        for (string_id, node, partner_set) in rows {
            let tx = omm.proteins[node.protein_id].transcript_ids.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(",");
            let features = node.features.join(",");
            let partner_ids = partner_set.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(",");
            writeln!(writer, "{gene_id}\t{gene_name}\t{n_feature_vectors}\t{}\t{string_id}\t{}\t{tx}\t{features}\t{}\t{partner_ids}", link_sets.len(), node.protein_id, partner_set.len())?;
        }
    }
    writer.flush()?;

    eprintln!("[ommverse] STRING links streamed: {link_records}; both endpoints resolved: {resolved_links}");
    eprintln!("[ommverse] candidate genes with >1 feature vector: {}", candidate_genes.len());
    eprintln!("[ommverse] candidate genes with >1 partner-gene link set: {genes_with_multiple_link_sets}");
    eprintln!("[ommverse] report: {}", out.display());
    Ok(())
}
