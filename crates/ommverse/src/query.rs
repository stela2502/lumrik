use crate::Ommverse;
use anyhow::{Context, Result, bail};
use rust_htslib::bcf::{self, Read as BcfRead};
use std::collections::HashSet;
use std::fmt;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Entity {
    Gene,
    Variant,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    Text(String),
    Integer(u64),
    InclusiveRange(u64, u64),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Predicate {
    Eq { field: String, value: Value },
    IsNull { field: String, negated: bool },
}

impl Predicate {
    fn field(&self) -> &str {
        match self {
            Self::Eq { field, .. } | Self::IsNull { field, .. } => field,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Query {
    pub entity: Entity,
    pub predicates: Vec<Predicate>,
}

impl Query {
    /// Parse the deliberately small first Ommverse query grammar:
    /// SELECT variant WHERE field = value [AND field = value ...]
    /// Numeric ranges use inclusive biological coordinates: start:end.
    pub fn parse(input: &str) -> Result<Self> {
        let input = input.trim();
        let upper = input.to_ascii_uppercase();
        if !upper.starts_with("SELECT ") {
            bail!("query must start with SELECT");
        }
        let where_at = upper.find(" WHERE ").context("query requires WHERE")?;
        let entity_text = input[7..where_at].trim();
        let entity = match entity_text.to_ascii_lowercase().as_str() {
            "gene" | "genes" => Entity::Gene,
            "variant" | "variants" => Entity::Variant,
            other => bail!("unsupported SELECT entity '{other}' (currently: gene, variant)"),
        };

        let where_text = &input[where_at + 7..];
        let mut predicates = Vec::new();
        for clause in split_and(where_text) {
            predicates.push(parse_predicate(clause)?);
        }
        if predicates.is_empty() {
            bail!("WHERE requires at least one predicate");
        }
        Ok(Self { entity, predicates })
    }
}

fn parse_predicate(clause: &str) -> Result<Predicate> {
    let clause = clause.trim();
    let upper = clause.to_ascii_uppercase();
    for (suffix, negated) in [(" IS NOT NULL", true), (" IS NULL", false)] {
        if upper.ends_with(suffix) {
            let field = clause[..clause.len() - suffix.len()].trim().to_ascii_lowercase();
            if field.is_empty() {
                bail!("empty field in WHERE clause");
            }
            return Ok(Predicate::IsNull { field, negated });
        }
    }
    let (field, raw_value) = clause
        .split_once('=')
        .with_context(|| format!("expected field = value or field IS [NOT] NULL in '{clause}'"))?;
    let field = field.trim().to_ascii_lowercase();
    if field.is_empty() {
        bail!("empty field in WHERE clause");
    }
    Ok(Predicate::Eq { field, value: parse_value(raw_value.trim())? })
}

fn split_and(input: &str) -> Vec<&str> {
    // Values in the first grammar cannot contain AND, so keep the parser tiny
    // while making keyword matching case-insensitive.
    let bytes = input.as_bytes();
    let mut cuts = Vec::new();
    let mut start = 0usize;
    let mut i = 0usize;
    let mut quoted = false;
    while i < bytes.len() {
        if bytes[i] == b'"' {
            quoted = !quoted;
            i += 1;
            continue;
        }
        if !quoted && i + 5 <= bytes.len() && input[i..i + 5].eq_ignore_ascii_case(" AND ") {
            cuts.push(&input[start..i]);
            start = i + 5;
            i = start;
            continue;
        }
        i += 1;
    }
    cuts.push(&input[start..]);
    cuts
}

fn parse_value(raw: &str) -> Result<Value> {
    if raw.starts_with('"') {
        if raw.len() < 2 || !raw.ends_with('"') {
            bail!("unterminated quoted value {raw}");
        }
        return Ok(Value::Text(raw[1..raw.len() - 1].to_owned()));
    }
    if let Some((start, end)) = raw.split_once(':') {
        let start: u64 = start.trim().parse().with_context(|| format!("invalid range start '{start}'"))?;
        let end: u64 = end.trim().parse().with_context(|| format!("invalid range end '{end}'"))?;
        if start == 0 || end == 0 || end < start {
            bail!("invalid inclusive biological range {start}:{end}");
        }
        return Ok(Value::InclusiveRange(start, end));
    }
    let value: u64 = raw.parse().with_context(|| format!("unsupported value '{raw}'"))?;
    Ok(Value::Integer(value))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantRow {
    pub chromosome: String,
    pub position: u64,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
}

impl fmt::Display for VariantRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.chromosome,
            self.position,
            self.id,
            self.reference,
            self.alternates.join(",")
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneRow {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub gene_id: usize,
    pub names: Vec<String>,
    pub strand: String,
}

impl fmt::Display for GeneRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chromosome,
            self.start,
            self.end,
            self.gene_id,
            self.names.join(","),
            self.strand,
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryRow {
    Gene(GeneRow),
    Variant(VariantRow),
}

impl fmt::Display for QueryRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Gene(row) => row.fmt(f),
            Self::Variant(row) => row.fmt(f),
        }
    }
}

pub fn execute(omm: &Ommverse, query: &Query) -> Result<Vec<QueryRow>> {
    match query.entity {
        Entity::Gene => execute_gene_query(omm, query)
            .map(|rows| rows.into_iter().map(QueryRow::Gene).collect()),
        Entity::Variant => execute_variant_query(omm, query)
            .map(|rows| rows.into_iter().map(QueryRow::Variant).collect()),
    }
}

fn execute_gene_query(omm: &Ommverse, query: &Query) -> Result<Vec<GeneRow>> {
    let chromosome = text_predicate(query, "chromosome")
        .or_else(|| text_predicate(query, "chrom"))
        .or_else(|| text_predicate(query, "chr"))
        .context("gene query requires chromosome = \"...\"")?;
    let (start1, end1) = position_predicate_for(query, "gene")?;

    for predicate in &query.predicates {
        match predicate.field() {
            "chromosome" | "chrom" | "chr" | "position" | "pos" => {}
            other => bail!("gene field '{other}' is not searchable yet"),
        }
    }

    let chr_id = omm.splice.chr_id(chromosome)
        .with_context(|| format!("chromosome '{chromosome}' not present in splice index"))?;
    let start0: u32 = (start1 - 1).try_into().context("gene query start exceeds supported genomic coordinate range")?;
    let end0: u32 = end1.try_into().context("gene query end exceeds supported genomic coordinate range")?;

    let mut gene_ids = omm.splice
        .overlapping_genes_with_bins(chr_id, start0, end0)
        .into_iter()
        .map(|(_, gene_id)| gene_id)
        .collect::<Vec<_>>();
    gene_ids.sort_unstable();
    gene_ids.dedup();

    let mut rows = Vec::new();
    for gene_id in gene_ids {
        let gene = &omm.splice.genes[gene_id];
        let mut span_start = u32::MAX;
        let mut span_end = 0u32;
        let mut strand = None;
        for &tx_id in gene.transcript_ids() {
            let tx = &omm.splice.transcripts[tx_id];
            if tx.chr_id != chr_id {
                continue;
            }
            span_start = span_start.min(omm.splice.tx_span_start[tx_id]);
            span_end = span_end.max(omm.splice.tx_span_end[tx_id]);
            strand.get_or_insert(tx.strand);
        }
        if span_start == u32::MAX || span_end <= span_start {
            continue;
        }
        let strand = match strand {
            Some(gtf_splice_index::Strand::Plus) => "+",
            Some(gtf_splice_index::Strand::Minus) => "-",
            _ => ".",
        };
        rows.push(GeneRow {
            chromosome: chromosome.to_owned(),
            start: span_start as u64 + 1,
            end: span_end as u64,
            gene_id,
            names: gene.names.clone(),
            strand: strand.to_owned(),
        });
    }
    rows.sort_by_key(|row| (row.start, row.end, row.gene_id));
    Ok(rows)
}

fn execute_variant_query(omm: &Ommverse, query: &Query) -> Result<Vec<VariantRow>> {
    let chromosome = text_predicate(query, "chromosome")
        .or_else(|| text_predicate(query, "chrom"))
        .or_else(|| text_predicate(query, "chr"))
        .context("variant query requires chromosome = \"...\"")?;
    let (start1, end1) = position_predicate_for(query, "variant")?;

    let clinical_filter = query.predicates.iter().find_map(|predicate| match predicate {
        Predicate::IsNull { field, negated } if field == "clinical_effect" => Some(*negated),
        _ => None,
    });

    for predicate in &query.predicates {
        match predicate.field() {
            "chromosome" | "chrom" | "chr" | "position" | "pos" | "clinical_effect" => {}
            other => bail!("variant field '{other}' is not searchable yet"),
        }
        if matches!(predicate, Predicate::Eq { field, .. } if field == "clinical_effect") {
            bail!("clinical_effect currently supports IS NULL or IS NOT NULL");
        }
    }

    let vcf = variant_source(omm)?;
    let mut reader = bcf::IndexedReader::from_path(&vcf)
        .with_context(|| format!("opening indexed variant source {}", vcf.display()))?;
    let header = reader.header().clone();
    let source_chromosome = resolve_chromosome(&omm.assembly, chromosome);
    let rid = header
        .name2rid(source_chromosome.as_bytes())
        .with_context(|| format!("chromosome '{chromosome}' ({source_chromosome}) not present in {}", vcf.display()))?;

    // User coordinates are 1-based inclusive. htslib fetch uses 0-based,
    // half-open coordinates, so [start1, end1] -> [start1 - 1, end1).
    reader.fetch(rid, start1 - 1, Some(end1))?;

    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record?;
        let pos1 = record.pos() as u64 + 1;
        if pos1 < start1 || pos1 > end1 {
            continue;
        }
        let alleles = record.alleles();
        if alleles.len() < 2 {
            continue;
        }
        let id = String::from_utf8_lossy(&record.id()).into_owned();
        rows.push(VariantRow {
            chromosome: chromosome.to_owned(),
            position: pos1,
            id: if id.is_empty() { ".".to_owned() } else { id },
            reference: String::from_utf8_lossy(alleles[0]).into_owned(),
            alternates: alleles[1..]
                .iter()
                .map(|allele| String::from_utf8_lossy(allele).into_owned())
                .collect(),
        });
    }
    if let Some(require_present) = clinical_filter {
        let clinical = clinical_variant_keys(omm, chromosome, start1, end1)?;
        rows.retain(|row| {
            let present = row.alternates.iter().any(|alt| {
                clinical.contains(&(row.position, row.reference.clone(), alt.clone()))
            });
            present == require_present
        });
    }
    Ok(rows)
}

fn clinical_variant_keys(
    omm: &Ommverse,
    chromosome: &str,
    start1: u64,
    end1: u64,
) -> Result<HashSet<(u64, String, String)>> {
    let vcf = source_for_kind(omm, "clinical_effects")?;
    let mut reader = bcf::IndexedReader::from_path(&vcf)
        .with_context(|| format!("opening indexed clinical-effect source {}", vcf.display()))?;
    let header = reader.header().clone();
    let rid = clinvar_rid(&header, chromosome)
        .with_context(|| format!("chromosome '{chromosome}' not present in {}", vcf.display()))?;
    reader.fetch(rid, start1 - 1, Some(end1))?;

    let mut keys = HashSet::new();
    for record in reader.records() {
        let record = record?;
        let pos1 = record.pos() as u64 + 1;
        if pos1 < start1 || pos1 > end1 {
            continue;
        }
        let alleles = record.alleles();
        if alleles.len() < 2 {
            continue;
        }
        let reference = String::from_utf8_lossy(alleles[0]).into_owned();
        for alt in &alleles[1..] {
            keys.insert((pos1, reference.clone(), String::from_utf8_lossy(alt).into_owned()));
        }
    }
    Ok(keys)
}

fn clinvar_rid(header: &bcf::header::HeaderView, chromosome: &str) -> Option<u32> {
    let bare = chromosome.strip_prefix("chr").unwrap_or(chromosome);
    [bare, chromosome].into_iter().find_map(|name| header.name2rid(name.as_bytes()).ok())
}

fn text_predicate<'a>(query: &'a Query, field: &str) -> Option<&'a str> {
    query.predicates.iter().find_map(|p| match p {
        Predicate::Eq { field: predicate_field, value: Value::Text(value) } if predicate_field == field => {
            Some(value.as_str())
        }
        _ => None,
    })
}

fn position_predicate_for(query: &Query, entity: &str) -> Result<(u64, u64)> {
    let value = query.predicates.iter().find_map(|p| match p {
        Predicate::Eq { field, value } if field == "position" || field == "pos" => Some(value),
        _ => None,
    }).with_context(|| format!("{entity} query requires position = N or position = START:END"))?;
    match value {
        Value::Integer(pos) if *pos > 0 => Ok((*pos, *pos)),
        Value::InclusiveRange(start, end) => Ok((*start, *end)),
        _ => bail!("position must be a positive integer or inclusive START:END range"),
    }
}

fn variant_source(omm: &Ommverse) -> Result<PathBuf> {
    source_for_kind(omm, "variants")
}

fn source_for_kind(omm: &Ommverse, kind: &str) -> Result<PathBuf> {
    let manifest = omm.sources.as_ref().context("Ommverse index has no source manifest")?;
    let resource = manifest.resources.iter()
        .find(|resource| resource.kind == kind)
        .with_context(|| format!("Ommverse source manifest has no {kind} resource"))?;
    let relative = resource.file.as_ref().with_context(|| format!("{kind} resource has no file"))?;
    let path = omm.source_root.join(relative);
    if !path.is_file() {
        bail!("{kind} source is missing: {}", path.display());
    }
    Ok(path)
}

fn resolve_chromosome(assembly: &str, chromosome: &str) -> String {
    if assembly != "hg38" || !chromosome.starts_with("chr") {
        return chromosome.to_owned();
    }
    let bare = &chromosome[3..];
    let number = match bare {
        "X" => 23,
        "Y" => 24,
        _ => match bare.parse::<u32>() {
            Ok(n @ 1..=22) => n,
            _ => return chromosome.to_owned(),
        },
    };
    let version = match number {
        1 => 11, 2 => 12, 3 => 12, 4 => 12, 5 => 10, 6 => 12,
        7 => 14, 8 => 11, 9 => 12, 10 => 11, 11 => 10, 12 => 12,
        13 => 11, 14 => 9, 15 => 10, 16 => 10, 17 => 11, 18 => 10,
        19 => 10, 20 => 11, 21 => 9, 22 => 11, 23 => 11, 24 => 10,
        _ => return chromosome.to_owned(),
    };
    format!("NC_{number:06}.{version}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_gene_range_query() {
        let query = Query::parse(
            r#"SELECT gene WHERE chromosome = "chr17" AND position = 7661779:7687550"#,
        ).unwrap();
        assert_eq!(query.entity, Entity::Gene);
        assert_eq!(query.predicates.len(), 2);
    }

    #[test]
    fn parses_variant_range_query() {
        let query = Query::parse(
            r#"SELECT variant WHERE chromosome = "chr17" AND position = 7674894:7679000"#,
        ).unwrap();
        assert_eq!(query.entity, Entity::Variant);
        assert_eq!(query.predicates.len(), 2);
        assert_eq!(
            query.predicates[1],
            Predicate::Eq { field: "position".to_owned(), value: Value::InclusiveRange(7_674_894, 7_679_000) }
        );
    }

    #[test]
    fn parses_exact_position_query() {
        let query = Query::parse(
            r#"select variant where chromosome = "chr17" and position = 7674894"#,
        ).unwrap();
        assert_eq!(
            query.predicates[1],
            Predicate::Eq { field: "position".to_owned(), value: Value::Integer(7_674_894) }
        );
    }

    #[test]
    fn parses_clinical_effect_presence_query() {
        let query = Query::parse(
            r#"SELECT variant WHERE chromosome = "chr17" AND position = 7674894:7679000 AND clinical_effect IS NOT NULL"#,
        ).unwrap();
        assert_eq!(
            query.predicates[2],
            Predicate::IsNull { field: "clinical_effect".to_owned(), negated: true }
        );
    }

    #[test]
    fn hg38_ucsc_chromosome_maps_to_refseq() {
        assert_eq!(resolve_chromosome("hg38", "chr17"), "NC_000017.11");
    }
}
