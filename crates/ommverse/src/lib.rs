//! Ommverse: Lumrik's integrated genome-to-protein biological reference model.
//!
//! Ommverse deliberately separates biological identity from source formats.
//! UCSC GTF, twoBit and UniProt bigBed files are import formats; callers see
//! genes, transcripts, proteins and protein features.

use anyhow::{bail, Context, Result};
use bigtools::BigBedRead;
use gtf_splice_index::{IdNameKeys, SpliceIndex, Strand, TranscriptId};
use gtf_splice_index::types::RefBlock;
use int_to_dna::{IntToDna, TwoBitReader};
use int_to_prot::IntToProt;
use hmm::{CategoricalEmission, Hmm};
use serde::{Deserialize, Serialize};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 4] = b"OMM1";
pub const OMMVERSE_FORMAT_VERSION: u32 = 5;

pub mod sources;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReviewStatus { SwissProt, Trembl, Other }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ProteinFeatureKind {
    Transmembrane, Domain, SignalPeptide, Cytoplasmic, Extracellular,
    ModifiedResidue, Disulfide, Repeat, Chain, Conflict, Interest,
    Mutagenesis, SpliceVariant, Structure, Other,
    ActiveSite, BindingSite, ConservedSite, ProteinFamily, HomologousSuperfamily, PtmSite,
}

impl ProteinFeatureKind {
    pub const ALL: [Self; 15] = [
        Self::Transmembrane, Self::Domain, Self::SignalPeptide, Self::Cytoplasmic,
        Self::Extracellular, Self::ModifiedResidue, Self::Disulfide, Self::Repeat,
        Self::Chain, Self::Conflict, Self::Interest, Self::Mutagenesis,
        Self::SpliceVariant, Self::Structure, Self::Other,
    ];
}

impl std::str::FromStr for ProteinFeatureKind {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        let normalized = value.to_ascii_lowercase().replace(['-', '_'], "");
        match normalized.as_str() {
            "transmembrane" | "tm" => Ok(Self::Transmembrane),
            "domain" => Ok(Self::Domain),
            "signalpeptide" => Ok(Self::SignalPeptide),
            "cytoplasmic" => Ok(Self::Cytoplasmic),
            "extracellular" => Ok(Self::Extracellular),
            "modifiedresidue" => Ok(Self::ModifiedResidue),
            "disulfide" => Ok(Self::Disulfide),
            "repeat" => Ok(Self::Repeat),
            "chain" => Ok(Self::Chain),
            "conflict" => Ok(Self::Conflict),
            "interest" => Ok(Self::Interest),
            "mutagenesis" => Ok(Self::Mutagenesis),
            "splicevariant" => Ok(Self::SpliceVariant),
            "structure" => Ok(Self::Structure),
            "other" => Ok(Self::Other),
            "activesite" => Ok(Self::ActiveSite),
            "bindingsite" => Ok(Self::BindingSite),
            "conservedsite" => Ok(Self::ConservedSite),
            "family" | "proteinfamily" => Ok(Self::ProteinFamily),
            "homologoussuperfamily" => Ok(Self::HomologousSuperfamily),
            "ptm" | "ptmsite" => Ok(Self::PtmSite),
            _ => bail!("unknown protein feature {value:?}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProteinFeature {
    pub kind: ProteinFeatureKind,
    /// Protein coordinates, 0-based half-open. None when UCSC did not provide
    /// a parseable amino-acid range for this feature.
    pub protein_range: Option<(u32, u32)>,
    pub label: String,
    pub description: String,
    pub chromosome: String,
    pub genomic_start: u32,
    pub genomic_end: u32,
    pub source_db: String,
    pub review_status: ReviewStatus,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Protein {
    pub accession: String,
    pub entry_name: String,
    pub name: String,
    pub gene_symbol: String,
    pub aliases: Vec<String>,
    pub ensembl_gene: Option<String>,
    pub ensembl_protein: Option<String>,
    pub transcript_ids: Vec<TranscriptId>,
    pub review_status: ReviewStatus,
    pub features: Vec<ProteinFeature>,
}


#[derive(Debug, Clone)]
pub struct ProteinTrainingExample {
    pub accession: String,
    pub gene_symbol: String,
    pub sequence: IntToProt,
    /// Residue-level feature truth, parallel to `sequence`.
    pub truth: Vec<bool>,
    pub feature_segments: usize,
}

#[derive(Debug, Clone, Default)]
pub struct DistributionSummary {
    pub count: usize,
    pub mean: f64,
    pub median: f64,
    /// Population standard deviation across the complete accepted corpus.
    pub sd: f64,
    pub q1: f64,
    pub q3: f64,
    pub min: usize,
    pub max: usize,
}

#[derive(Debug, Clone, Default)]
pub struct TrainingCorpusReport {
    pub candidate_proteins: usize,
    pub reconstructed_proteins: usize,
    pub rejected_sequence: usize,
    pub rejected_missing_range: usize,
    pub rejected_out_of_bounds: usize,
    pub total_residues: usize,
    pub feature_residues: usize,
    pub feature_segments: usize,
    /// Lengths of individual curated feature intervals, in residues.
    pub feature_lengths: DistributionSummary,
    /// Background residues between consecutive curated feature intervals on the same protein.
    pub inter_feature_gaps: DistributionSummary,
    /// Number of adjacent feature pairs whose two requested flank windows would overlap.
    pub overlapping_flank_pairs: usize,
    /// Diagnostic flank width used for `overlapping_flank_pairs`.
    pub diagnostic_flank_width: usize,
}

#[derive(Debug, Clone)]
pub struct ProteinTrainingCorpus {
    pub feature: ProteinFeatureKind,
    pub train: Vec<ProteinTrainingExample>,
    pub test: Vec<ProteinTrainingExample>,
    pub report: TrainingCorpusReport,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExactAaFeatureModel {
    pub feature: ProteinFeatureKind,
    pub flank_width: usize,
    /// BACKGROUND, PRE_FEATURE, FEATURE, POST_FEATURE.
    pub initial: [f64; 4],
    /// Row-major four-state transition probabilities.
    pub transition: [f64; 16],
    /// Exact five-bit IntToProt code probabilities for the four states.
    pub emission: [[f64; 32]; 4],
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct FeatureEvaluation {
    pub proteins: usize,
    pub residues: usize,
    pub true_positive: usize,
    pub false_positive: usize,
    pub true_negative: usize,
    pub false_negative: usize,
    pub truth_segments: usize,
    pub predicted_segments: usize,
    pub recovered_segments: usize,
}

impl FeatureEvaluation {
    pub fn precision(&self) -> f64 { ratio(self.true_positive, self.true_positive + self.false_positive) }
    pub fn recall(&self) -> f64 { ratio(self.true_positive, self.true_positive + self.false_negative) }
    pub fn specificity(&self) -> f64 { ratio(self.true_negative, self.true_negative + self.false_positive) }
    pub fn f1(&self) -> f64 {
        let p = self.precision(); let r = self.recall();
        if p + r == 0.0 { 0.0 } else { 2.0 * p * r / (p + r) }
    }
    pub fn segment_recall(&self) -> f64 { ratio(self.recovered_segments, self.truth_segments) }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TopologyState {
    Other = 0,
    Transmembrane = 1,
    CytoShortLoop = 2,
    CytoMediumLoop = 3,
    CytoLongRegion = 4,
    ExtraShortLoop = 5,
    ExtraMediumLoop = 6,
    ExtraLongRegion = 7,
}

impl TopologyState {
    pub const COUNT: usize = 8;
    /// Public biological outputs. Internal cytoplasmic/extracellular substates are
    /// collapsed back to these feature classes for evaluation.
    pub const BIOLOGICAL: [ProteinFeatureKind; 3] = [
        ProteinFeatureKind::Transmembrane,
        ProteinFeatureKind::Cytoplasmic,
        ProteinFeatureKind::Extracellular,
    ];
    pub const ALL: [Self; 8] = [
        Self::Other, Self::Transmembrane, Self::CytoShortLoop, Self::CytoMediumLoop,
        Self::CytoLongRegion, Self::ExtraShortLoop, Self::ExtraMediumLoop,
        Self::ExtraLongRegion,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::Other => "OTHER", Self::Transmembrane => "TM",
            Self::CytoShortLoop => "CYTO_SHORT", Self::CytoMediumLoop => "CYTO_MEDIUM",
            Self::CytoLongRegion => "CYTO_LONG",
            Self::ExtraShortLoop => "EXTRA_SHORT", Self::ExtraMediumLoop => "EXTRA_MEDIUM",
            Self::ExtraLongRegion => "EXTRA_LONG",
        }
    }
    fn index(self) -> usize { self as usize }
    fn feature(self) -> Option<ProteinFeatureKind> {
        match self {
            Self::Other => None,
            Self::Transmembrane => Some(ProteinFeatureKind::Transmembrane),
            Self::CytoShortLoop | Self::CytoMediumLoop | Self::CytoLongRegion => Some(ProteinFeatureKind::Cytoplasmic),
            Self::ExtraShortLoop | Self::ExtraMediumLoop | Self::ExtraLongRegion => Some(ProteinFeatureKind::Extracellular),
        }
    }
    fn region_state(feature: ProteinFeatureKind, start: usize, end: usize, _protein_len: usize) -> Self {
        // Terminality is known geometry, not a sequence property. Do not give the
        // HMM dedicated terminal escape states; terminal regions use the same
        // length-based latent contexts as every other CYTO/EXTRA annotation.
        let len = end.saturating_sub(start);
        match (feature, len) {
            (ProteinFeatureKind::Cytoplasmic, 0..=15) => Self::CytoShortLoop,
            (ProteinFeatureKind::Cytoplasmic, 16..=40) => Self::CytoMediumLoop,
            (ProteinFeatureKind::Cytoplasmic, _) => Self::CytoLongRegion,
            (ProteinFeatureKind::Extracellular, 0..=15) => Self::ExtraShortLoop,
            (ProteinFeatureKind::Extracellular, 16..=40) => Self::ExtraMediumLoop,
            (ProteinFeatureKind::Extracellular, _) => Self::ExtraLongRegion,
            _ => Self::Other,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TopologyTrainingExample {
    pub accession: String,
    pub sequence: IntToProt,
    pub truth: Vec<usize>,
}

#[derive(Debug, Clone, Default)]
pub struct TopologyTrainingReport {
    pub candidate_proteins: usize,
    pub reconstructed_proteins: usize,
    pub rejected_sequence: usize,
    pub rejected_missing_range: usize,
    pub rejected_out_of_bounds: usize,
    pub conflicting_residues: usize,
}

#[derive(Debug, Clone)]
pub struct TopologyTrainingCorpus {
    pub train: Vec<TopologyTrainingExample>,
    pub test: Vec<TopologyTrainingExample>,
    pub report: TopologyTrainingReport,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExactAaTopologyModel {
    pub initial: Vec<f64>,
    pub transition: Vec<f64>,
    pub emission: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemistryTopologyModel {
    pub categories: Vec<u16>,
    pub initial: Vec<f64>,
    pub transition: Vec<f64>,
    pub emission: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TopologyEvaluation {
    /// One binary evaluation per competitive biological state, in
    /// TopologyState::BIOLOGICAL order.
    pub states: Vec<FeatureEvaluation>,
    /// Residue counts for each internal latent state. These make it possible to
    /// see which hidden explanation is consuming sequence before CYTO/EXTRA are
    /// collapsed to their public biological labels.
    pub latent_truth_residues: Vec<usize>,
    pub latent_predicted_residues: Vec<usize>,
    /// Row-major truth x predicted latent-state confusion matrix.
    pub latent_confusion: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemistryFeatureModel {
    pub feature: ProteinFeatureKind,
    pub flank_width: usize,
    pub initial: [f64; 4],
    pub transition: [f64; 16],
    /// Sorted chemistry bitmasks observed in training. The final emission column
    /// is reserved for chemistry combinations absent from training.
    pub categories: Vec<u16>,
    pub emission: Vec<Vec<f64>>,
}

pub const AA_MODEL_VAULT_FORMAT_VERSION: u32 = 5;
const AA_MODEL_VAULT_MAGIC: &[u8; 4] = b"AAV5";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModelVaultMetadata {
    pub assembly: String,
    pub source: String,
    pub ommverse_format_version: u32,
    pub vault_format_version: u32,
    pub flank_width: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ProteinFeatureModel {
    ExactAa { model: ExactAaFeatureModel, evaluation: FeatureEvaluation },
    Chemistry { model: ChemistryFeatureModel, evaluation: FeatureEvaluation },
    /// Opinionated membrane-architecture HMM. Cytoplasmic/extracellular sequence
    /// is represented by latent loop/long/terminal substates and collapsed on output.
    TopologyExactAa { model: ExactAaTopologyModel, evaluation: TopologyEvaluation },
    TopologyChemistry { model: ChemistryTopologyModel, evaluation: TopologyEvaluation },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AaModelVault {
    pub metadata: ModelVaultMetadata,
    pub models: Vec<ProteinFeatureModel>,
}

#[derive(Debug, Clone, Default)]
pub struct ModelVaultTrainingReport {
    pub feature_classes_considered: usize,
    pub feature_classes_trained: usize,
    pub models_trained: usize,
    pub topology_train_proteins: usize,
    pub topology_test_proteins: usize,
    pub topology_conflicting_residues: usize,
    pub skipped: Vec<(ProteinFeatureKind, String)>,
}

#[derive(Debug, Clone, Default)]
pub struct UnannotatedScanReport {
    pub candidate_proteins: usize,
    pub reconstructed_proteins: usize,
    pub rejected_sequence: usize,
    pub predicted_feature_proteins: usize,
    pub predicted_segments: usize,
    pub predicted_residues: usize,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BuildReport {
    pub transcripts: usize,
    pub mapped_proteins: usize,
    pub linked_transcripts: usize,
    pub feature_records: usize,
    pub unlinked_mapping_records: usize,
    pub feature_records_without_protein: usize,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterProEntry {
    pub kind: ProteinFeatureKind,
    pub name: String,
    pub parents: Vec<String>,
    pub children: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV1 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV2 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
    interpro_entries: HashMap<String, InterProEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV3 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
    interpro_entries: HashMap<String, InterProEntry>,
    chromatin: Vec<ChromatinElement>,
}

/// A merged genomic interval where ENCODE has observed at least one
/// DNA-associated protein binding event in any assayed biosample.  Ommverse
/// deliberately stores only this union: factor-, experiment- and biosample-
/// specific evidence remains in the cached ENCODE bigBed for downstream models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProteinBindingRegion {
    pub chromosome_id: u16,
    pub region: RefBlock,
    /// Number of source rPeaks collapsed into this union interval.
    pub source_peak_count: u32,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ProteinBindingUnion {
    pub chromosomes: Vec<String>,
    pub regions: Vec<ProteinBindingRegion>,
    /// Number of ENCODE rPeaks scanned to construct the union.
    pub source_peak_count: usize,
    pub source: String,
    pub source_url: String,
}

// Version-4 compatibility only.  v4 serialized every ENCODE rPeak.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
struct ProteinBindingPeakV4 {
    chromosome_id: u16, region: RefBlock, score: u16, factor_id: u16, peak_id: u32,
    observed_experiments: u16, assayed_experiments: u16, ccre_id: u32,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct ProteinBindingIndexV4 {
    chromosomes: Vec<String>, factors: Vec<String>, ccres: Vec<String>,
    peaks: Vec<ProteinBindingPeakV4>, source: String, source_url: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV4 {
    assembly: String, source_root: PathBuf, genome_twobit: PathBuf, splice: SpliceIndex,
    proteins: Vec<Protein>, report: BuildReport, interpro_entries: HashMap<String, InterProEntry>,
    chromatin: Vec<ChromatinElement>, protein_binding: ProteinBindingIndexV4,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChromatinElement {
    pub chromosome: String,
    pub region: RefBlock,
    /// Source-native stable/display identifier.
    pub name: String,
    /// Source-native score; FANTOM5 enhancer BED uses the BED score column.
    pub score: u32,
    /// Source-native sub-block geometry, stored as absolute genomic intervals.
    pub blocks: Vec<RefBlock>,
    pub source: String,
    pub source_url: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ommverse {
    pub assembly: String,
    pub source_root: PathBuf,
    pub genome_twobit: PathBuf,
    pub splice: SpliceIndex,
    pub proteins: Vec<Protein>,
    pub report: BuildReport,
    /// InterPro vocabulary and parent/child hierarchy retained once per index.
    pub interpro_entries: HashMap<String, InterProEntry>,
    /// Reference regulatory/chromatin observations imported for this assembly.
    pub chromatin: Vec<ChromatinElement>,
    /// Union of loci with ENCODE DNA-associated-protein binding evidence.
    pub protein_binding: ProteinBindingUnion,
    #[serde(skip)]
    protein_by_accession: HashMap<String, usize>,
    #[serde(skip)]
    proteins_by_gene: HashMap<String, Vec<usize>>,
}

fn find_gtf(genes_dir: &Path) -> Result<PathBuf> {
    if !genes_dir.is_dir() {
        bail!("required UCSC Genes directory missing: {}", genes_dir.display());
    }

    let mut gtfs = std::fs::read_dir(genes_dir)
        .with_context(|| format!("reading Genes directory {}", genes_dir.display()))?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
                return false;
            };
            path.is_file() && (name.ends_with(".gtf") || name.ends_with(".gtf.gz"))
        })
        .collect::<Vec<_>>();
    gtfs.sort();

    gtfs.into_iter().next().with_context(|| {
        format!(
            "no GTF annotation (*.gtf or *.gtf.gz) found in {}",
            genes_dir.display()
        )
    })
}

#[derive(Debug, Clone, Default)]
pub struct InterProImportReport {
    pub entries_loaded: usize,
    pub records_streamed: usize,
    pub matched_records: usize,
    pub features_added: usize,
    pub duplicate_features: usize,
    pub malformed_records: usize,
    pub unknown_entries: usize,
}

fn interpro_feature_kind(value: &str) -> ProteinFeatureKind {
    match value.to_ascii_lowercase().replace(['-', '_', ' '], "").as_str() {
        "domain" => ProteinFeatureKind::Domain,
        "family" => ProteinFeatureKind::ProteinFamily,
        "homologoussuperfamily" => ProteinFeatureKind::HomologousSuperfamily,
        "repeat" => ProteinFeatureKind::Repeat,
        "activesite" => ProteinFeatureKind::ActiveSite,
        "bindingsite" => ProteinFeatureKind::BindingSite,
        "conservedsite" => ProteinFeatureKind::ConservedSite,
        "ptm" | "ptmsite" => ProteinFeatureKind::PtmSite,
        _ => ProteinFeatureKind::Other,
    }
}


pub fn ingest_interpro_many<F>(
    ommverses: &mut [Ommverse],
    protein2ipr: &Path,
    entry_list: &Path,
    parent_child_tree: Option<&Path>,
    mut progress: F,
) -> Result<Vec<InterProImportReport>>
where
    F: FnMut(usize, &[InterProImportReport]),
{
    use flate2::read::MultiGzDecoder;
    use std::io::{BufRead, BufReader};

    let mut entries = HashMap::<String, (ProteinFeatureKind, String)>::new();
    let reader = BufReader::new(File::open(entry_list).with_context(|| format!("opening {}", entry_list.display()))?);
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        let mut f = line.split('\t');
        let Some(ipr) = f.next().map(str::trim).filter(|x| !x.is_empty()) else { continue; };
        let kind = interpro_feature_kind(f.next().unwrap_or("").trim());
        let name = f.next().unwrap_or("").trim().to_owned();
        entries.insert(ipr.to_owned(), (kind, name));
    }

    for omm in ommverses.iter_mut() {
        for (ipr, (kind, name)) in &entries {
            omm.interpro_entries.entry(ipr.clone()).or_insert_with(|| InterProEntry {
                kind: *kind, name: name.clone(), parents: Vec::new(), children: Vec::new(),
            });
        }
        if let Some(tree) = parent_child_tree {
            let reader = BufReader::new(File::open(tree).with_context(|| format!("opening {}", tree.display()))?);
            let mut stack: Vec<String> = Vec::new();
            for line in reader.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') { continue; }
                let mut depth = 0usize;
                let bytes = line.as_bytes();
                while bytes.get(depth * 2..depth * 2 + 2) == Some(b"--") { depth += 1; }
                let body = &line[depth * 2..];
                let Some(ipr) = body.split("::").next().map(str::trim).filter(|x| x.starts_with("IPR")) else { continue; };
                stack.truncate(depth);
                if depth > 0 {
                    if let Some(parent) = stack.get(depth - 1).cloned() {
                        if let Some(entry) = omm.interpro_entries.get_mut(ipr) {
                            if !entry.parents.contains(&parent) { entry.parents.push(parent.clone()); }
                        }
                        if let Some(entry) = omm.interpro_entries.get_mut(&parent) {
                            if !entry.children.iter().any(|x| x == ipr) { entry.children.push(ipr.to_owned()); }
                        }
                    }
                }
                stack.push(ipr.to_owned());
            }
        }
    }

    let mut locations = HashMap::<String, Vec<(usize, usize)>>::new();
    for (oi, omm) in ommverses.iter().enumerate() {
        for (pi, protein) in omm.proteins.iter().enumerate() {
            locations.entry(protein.accession.clone()).or_default().push((oi, pi));
        }
    }
    let mut reports = vec![InterProImportReport::default(); ommverses.len()];
    for report in &mut reports { report.entries_loaded = entries.len(); }

    let file = File::open(protein2ipr).with_context(|| format!("opening {}", protein2ipr.display()))?;
    let gz = protein2ipr.extension().is_some_and(|x| x == "gz");
    let input: Box<dyn Read> = if gz { Box::new(MultiGzDecoder::new(file)) } else { Box::new(file) };
    let reader = BufReader::with_capacity(1024 * 1024, input);
    let mut streamed = 0usize;
    for line in reader.lines() {
        let line = line?;
        streamed += 1;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            for r in &mut reports { r.malformed_records += 1; }
            continue;
        }
        let accession = f[0].trim();
        let Some(targets) = locations.get(accession) else {
            if streamed % 1_000_000 == 0 { progress(streamed, &reports); }
            continue;
        };
        let ipr = f[1].trim();
        let signature = f[3].trim();
        let start_1: u32 = match f[4].trim().parse() { Ok(v) if v > 0 => v, _ => { for &(oi, _) in targets { reports[oi].malformed_records += 1; } continue; } };
        let end_1: u32 = match f[5].trim().parse() { Ok(v) if v >= start_1 => v, _ => { for &(oi, _) in targets { reports[oi].malformed_records += 1; } continue; } };
        let known = entries.get(ipr);
        let (kind, entry_name) = known.cloned().unwrap_or_else(|| (ProteinFeatureKind::Other, f[2].trim().to_owned()));
        let description = if signature.is_empty() { entry_name } else if entry_name.is_empty() { format!("member signature {signature}") } else { format!("{entry_name} [{signature}]") };
        for &(oi, pi) in targets {
            let report = &mut reports[oi];
            report.records_streamed = streamed;
            report.matched_records += 1;
            if known.is_none() { report.unknown_entries += 1; }
            let feature = ProteinFeature { kind, protein_range: Some((start_1 - 1, end_1)), label: ipr.to_owned(), description: description.clone(), chromosome: String::new(), genomic_start: 0, genomic_end: 0, source_db: "InterPro".to_owned(), review_status: ReviewStatus::Other };
            if !ommverses[oi].proteins[pi].features.contains(&feature) { ommverses[oi].proteins[pi].features.push(feature); report.features_added += 1; } else { report.duplicate_features += 1; }
        }
        if streamed % 1_000_000 == 0 { progress(streamed, &reports); }
    }
    for (omm, report) in ommverses.iter_mut().zip(reports.iter_mut()) {
        report.records_streamed = streamed;
        omm.report.feature_records += report.features_added;
    }
    progress(streamed, &reports);
    Ok(reports)
}

impl Ommverse {
    /// Fetch all currently supported reference sources for an assembly into a cache,
    /// then build the integrated Ommverse from those source files.
    pub fn fetch_and_build(assembly: &str, cache_root: impl AsRef<Path>) -> Result<Self> {
        Self::fetch_and_build_with_progress(assembly, cache_root, |_, _| {})
    }

    /// Fetch supported sources and build while reporting coarse build stages and
    /// record counters. The callback is intentionally cheap so callers can feed
    /// the shared Lumrik status server without putting it on hot per-record paths.
    pub fn fetch_and_build_with_progress(
        assembly: &str,
        cache_root: impl AsRef<Path>,
        mut progress: impl FnMut(&str, usize),
    ) -> Result<Self> {
        progress("resolving UCSC sources", 0);
        let root = sources::ucsc::fetch(assembly, cache_root.as_ref())?;
        progress("resolving FANTOM5 sources", 0);
        sources::fantom5::fetch(assembly, &root)?;
        progress("resolving ENCODE4 sources", 0);
        sources::encode4::fetch(assembly, &root)?;
        Self::build_ucsc_with_debug_and_progress(root, false, progress)
    }

    fn ingest_fantom5_if_present(&mut self) -> Result<()> {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let Some(url) = sources::fantom5::enhancer_url(&self.assembly) else { return Ok(()); };
        let path = self.source_root.join("Chromatin").join("FANTOM5").join("F5.hg38.enhancers.bed.gz");
        if !path.is_file() { return Ok(()); }
        let file = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
        let reader = BufReader::new(MultiGzDecoder::new(file));
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') { continue; }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 12 { bail!("{}:{}: expected BED12, got {} fields", path.display(), line_no + 1, f.len()); }
            let start: u32 = f[1].parse()?;
            let end: u32 = f[2].parse()?;
            if start >= end { bail!("{}:{}: invalid interval {start}-{end}", path.display(), line_no + 1); }
            let block_count: usize = f[9].parse()?;
            let sizes: Vec<u32> = f[10].trim_end_matches(',').split(',').map(str::parse).collect::<std::result::Result<_, _>>()?;
            let starts: Vec<u32> = f[11].trim_end_matches(',').split(',').map(str::parse).collect::<std::result::Result<_, _>>()?;
            if sizes.len() != block_count || starts.len() != block_count {
                bail!("{}:{}: BED block count mismatch", path.display(), line_no + 1);
            }
            let mut blocks = Vec::with_capacity(block_count);
            for (&size, &offset) in sizes.iter().zip(&starts) {
                let block_start = start.checked_add(offset).context("FANTOM5 block start overflow")?;
                let block_end = block_start.checked_add(size).context("FANTOM5 block end overflow")?;
                if block_start < start || block_end > end || block_start >= block_end {
                    bail!("{}:{}: invalid BED block", path.display(), line_no + 1);
                }
                blocks.push(RefBlock::new(block_start, block_end));
            }
            self.chromatin.push(ChromatinElement {
                chromosome: f[0].to_owned(), region: RefBlock::new(start, end), name: f[3].to_owned(),
                score: f[4].parse().unwrap_or(0), blocks, source: "FANTOM5 enhancer".to_owned(), source_url: url.to_owned(),
            });
        }
        Ok(())
    }

    fn ingest_encode4_binding_if_present(&mut self, progress: &mut dyn FnMut(&str, usize)) -> Result<()> {
        let Some(url) = sources::encode4::tf_rpeaks_url(&self.assembly) else { return Ok(()); };
        let path = self.source_root.join("Chromatin").join("ENCODE4").join("TFrPeakClusters.bb");
        if !path.is_file() { return Ok(()); }

        // The bigBed is coordinate sorted.  Collapse all overlapping/touching
        // rPeaks while streaming, so memory scales with the union rather than
        // with the ~22 million source observations.
        let mut chromosome_ids = HashMap::<String, u16>::new();
        let mut index = ProteinBindingUnion {
            source: "ENCODE4 TF rPeak union".to_owned(),
            source_url: url.to_owned(),
            ..ProteinBindingUnion::default()
        };
        let mut current: Option<ProteinBindingRegion> = None;
        let mut imported = 0usize;

        read_all_bigbed(&path, |chrom, start, end, _rest| {
            let chromosome_id = if let Some(&id) = chromosome_ids.get(chrom) { id } else {
                let id = u16::try_from(index.chromosomes.len()).context("too many ENCODE chromosomes")?;
                index.chromosomes.push(chrom.to_owned());
                chromosome_ids.insert(chrom.to_owned(), id);
                id
            };
            let next = RefBlock::new(start, end);
            match current.as_mut() {
                Some(region) if region.chromosome_id == chromosome_id && next.start <= region.region.end => {
                    region.region.end = region.region.end.max(next.end);
                    region.source_peak_count = region.source_peak_count.saturating_add(1);
                }
                _ => {
                    if let Some(region) = current.take() { index.regions.push(region); }
                    current = Some(ProteinBindingRegion { chromosome_id, region: next, source_peak_count: 1 });
                }
            }
            imported += 1;
            if imported % 1_000_000 == 0 {
                eprintln!("[ommverse] ENCODE4 TF rPeaks: {imported} scanned, {} union regions", index.regions.len());
                progress("building ENCODE4 binding union", imported);
            }
            Ok(())
        })?;
        if let Some(region) = current.take() { index.regions.push(region); }
        index.source_peak_count = imported;
        self.protein_binding = index;
        progress("building ENCODE4 binding union", imported);
        Ok(())
    }

    /// Build an Ommverse directly from the directory layout produced by
    /// scripts/download_ucsc_reference.sh.
    pub fn build_ucsc(root: impl AsRef<Path>) -> Result<Self> {
        Self::build_ucsc_with_debug_and_progress(root, false, |_, _| {})
    }

    /// Build from UCSC and optionally report a bounded sample of feature records
    /// that could not be attached to a mapped protein.
    pub fn build_ucsc_with_debug(root: impl AsRef<Path>, debug_failed_mappings: bool) -> Result<Self> {
        Self::build_ucsc_with_debug_and_progress(root, debug_failed_mappings, |_, _| {})
    }

    pub fn build_ucsc_with_debug_and_progress(
        root: impl AsRef<Path>,
        debug_failed_mappings: bool,
        mut progress: impl FnMut(&str, usize),
    ) -> Result<Self> {
        let root = root.as_ref().canonicalize().with_context(|| format!("reference root {}", root.as_ref().display()))?;
        let assembly = root.file_name().and_then(|x| x.to_str()).context("reference root has no assembly name")?.to_owned();
        let genes_dir = root.join("Genes");
        let gtf = find_gtf(&genes_dir)?;
        let twobit = root.join("Genome").join(format!("{assembly}.2bit"));
        let protein_dir = root.join("Protein");
        for required in [&twobit, &protein_dir] {
            if !required.exists() { bail!("required UCSC reference component missing: {}", required.display()); }
        }

        progress("building transcript index", 0);
        let splice = SpliceIndex::from_path(&gtf, 100_000, IdNameKeys::default())
            .with_context(|| format!("building splice index from {}", gtf.display()))?;
        let mut tx_by_stable = HashMap::<String, TranscriptId>::new();
        for tx in &splice.transcripts {
            for name in &tx.names {
                tx_by_stable.entry(strip_version(name).to_owned()).or_insert(tx.id);
            }
        }

        progress("mapping proteins", 0);
        let (mapping, mapping_schema) = find_swissprot_mapping_bigbed(&protein_dir)?;
        let mut proteins = Vec::<Protein>::new();
        let mut protein_by_accession = HashMap::<String, usize>::new();
        let mut linked = HashSet::new();
        let mut unlinked_mapping_records = 0usize;

        read_all_bigbed(&mapping, |_, _, _, rest| {
            let f: Vec<&str> = rest.split('\t').collect();
            // UCSC mapping schemas vary between assemblies and annotation
            // namespaces.  Do not assume the transcript identifier lives in a
            // fixed column: search from right to left because the most specific
            // cross-reference fields normally live at the end of these records.
            // Only exact GTF-known identifiers (with an optional numeric version
            // stripped) are accepted.
            if f.is_empty() { return Ok(()); }
            // ensGene_*.swissprot.bb is a transcript mapping table, while
            // unipAliSwissprot.bb is a bigPsl protein-to-genome alignment.  In
            // bigPsl the BED name (rest field 0) is the UniProt accession; the
            // later fields are alignment statistics, not UniProt metadata.
            let accession = f[0].trim();
            if accession.is_empty() { return Ok(()); }

            let tx_id = if mapping_schema == MappingSchema::EnsGene {
                let Some(tx_id) = f.iter().rev().find_map(|value| {
                    value
                        .split(|c: char| !(c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '-'))
                        .rev()
                        .filter(|candidate| !candidate.is_empty())
                        .find_map(|candidate| tx_by_stable.get(strip_version(candidate)).copied())
                }) else {
                    unlinked_mapping_records += 1;
                    return Ok(());
                };
                linked.insert(tx_id);
                Some(tx_id)
            } else {
                // The canonical UCSC track already proves that this UniProt
                // protein aligns to this assembly.  Transcript linkage is added
                // below from protMapInfo.tsv, whose job is to bridge the protein
                // and assembly transcript namespaces.
                None
            };

            let idx = *protein_by_accession.entry(accession.to_owned()).or_insert_with(|| {
                let idx = proteins.len();
                proteins.push(Protein {
                    accession: accession.to_owned(),
                    entry_name: if mapping_schema == MappingSchema::EnsGene { f.get(22).copied().unwrap_or("").to_owned() } else { String::new() },
                    review_status: if mapping_schema == MappingSchema::EnsGene { review_status(f.get(23).copied().unwrap_or("")) } else { ReviewStatus::SwissProt },
                    name: if mapping_schema == MappingSchema::EnsGene { f.get(26).copied().unwrap_or("").to_owned() } else { String::new() },
                    gene_symbol: if mapping_schema == MappingSchema::EnsGene { f.get(27).copied().unwrap_or("").to_owned() } else { String::new() },
                    aliases: if mapping_schema == MappingSchema::EnsGene { parse_aliases(f.get(30).copied().unwrap_or(""), f.get(31).copied().unwrap_or("")) } else { Vec::new() },
                    ensembl_gene: if mapping_schema == MappingSchema::EnsGene { nonempty(f.get(f.len().saturating_sub(3)).copied().unwrap_or("")) } else { None },
                    ensembl_protein: if mapping_schema == MappingSchema::EnsGene { nonempty(f.get(f.len().saturating_sub(2)).copied().unwrap_or("")) } else { None },
                    transcript_ids: Vec::new(),
                    features: Vec::new(),
                });
                idx
            });
            if let Some(tx_id) = tx_id {
                if !proteins[idx].transcript_ids.contains(&tx_id) { proteins[idx].transcript_ids.push(tx_id); }
            }
            Ok(())
        })?;

        // UCSC's protMapInfo.tsv is a source-agnostic rescue bridge between
        // UniProt accessions and the transcript namespace used for this
        // assembly.  Some assemblies (notably RefSeq-backed ones) do not have
        // an ensGene Swiss-Prot mapping track, and unipAli*.bb may use isoform
        // accessions as their BED names.  Prefer the mapping BigBed above, but
        // augment it here with canonical UniProt accessions whenever an exact
        // GTF-known transcript identifier is present in protMapInfo.tsv.
        //
        // Example UCSC row:
        // Q9BL78  trembl  ...  NM_058267.6  NM_058267.6  chrI:55339-64021
        //
        // This is intentionally still an exact identifier rescue: no fuzzy
        // gene-name or coordinate inference is performed here.
        let prot_map_info = protein_dir.join("protMapInfo.tsv");
        if prot_map_info.is_file() {
            let text = fs::read_to_string(&prot_map_info)
                .with_context(|| format!("reading {}", prot_map_info.display()))?;
            for line in text.lines() {
                if line.trim().is_empty() || line.starts_with('#') { continue; }
                let f: Vec<&str> = line.split('\t').collect();
                let accession = f.first().copied().unwrap_or("").trim();
                if accession.is_empty() { continue; }

                let Some(tx_id) = f.iter().rev().find_map(|value| {
                    value
                        .split(|c: char| !(c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '-'))
                        .rev()
                        .filter(|candidate| !candidate.is_empty())
                        .find_map(|candidate| tx_by_stable.get(strip_version(candidate)).copied())
                }) else {
                    continue;
                };

                linked.insert(tx_id);
                let idx = *protein_by_accession.entry(accession.to_owned()).or_insert_with(|| {
                    let idx = proteins.len();
                    proteins.push(Protein {
                        accession: accession.to_owned(),
                        entry_name: String::new(),
                        name: String::new(),
                        gene_symbol: String::new(),
                        aliases: Vec::new(),
                        ensembl_gene: None,
                        ensembl_protein: None,
                        transcript_ids: Vec::new(),
                        review_status: review_status(f.get(1).copied().unwrap_or("")),
                        features: Vec::new(),
                    });
                    idx
                });
                if !proteins[idx].transcript_ids.contains(&tx_id) {
                    proteins[idx].transcript_ids.push(tx_id);
                }
            }
        }

        let feature_tracks = [
            ("unipLocTransMemb.bb", ProteinFeatureKind::Transmembrane),
            ("unipDomain.bb", ProteinFeatureKind::Domain),
            ("unipLocSignal.bb", ProteinFeatureKind::SignalPeptide),
            ("unipLocCytopl.bb", ProteinFeatureKind::Cytoplasmic),
            ("unipLocExtra.bb", ProteinFeatureKind::Extracellular),
            ("unipModif.bb", ProteinFeatureKind::ModifiedResidue),
            ("unipDisulfBond.bb", ProteinFeatureKind::Disulfide),
            ("unipRepeat.bb", ProteinFeatureKind::Repeat),
            ("unipChain.bb", ProteinFeatureKind::Chain),
            ("unipConflict.bb", ProteinFeatureKind::Conflict),
            ("unipInterest.bb", ProteinFeatureKind::Interest),
            ("unipMut.bb", ProteinFeatureKind::Mutagenesis),
            ("unipSplice.bb", ProteinFeatureKind::SpliceVariant),
            ("unipStruct.bb", ProteinFeatureKind::Structure),
        ];
        let mut feature_records = 0usize;
        let mut orphan_features = 0usize;
        let mut warnings = Vec::new();
        for (file, kind) in feature_tracks {
            let path = protein_dir.join(file);
            if !path.exists() { warnings.push(format!("optional protein track missing: {file}")); continue; }
            let mut failed_debug_printed = 0usize;
            read_all_bigbed(&path, |chrom, start, end, rest| {
                let f: Vec<&str> = rest.split('\t').collect();
                if f.is_empty() { return Ok(()); }

                // UCSC UniProt feature bigBeds are BED12+14 after BigBed's BED3
                // prefix.  The explicit UniProt accession is therefore f[24].
                // f[19] independently states "... on protein ACCESSION" and is
                // used as a schema/identity cross-check rather than as a rescue.
                let accession = f.get(24).copied().unwrap_or("").trim();
                let aa_text = f.get(19).copied().unwrap_or("");
                let described_accession = protein_accession_from_feature_text(aa_text);
                let identity_consistent = !accession.is_empty()
                    && described_accession.map_or(true, |described| described == accession);

                let Some(&idx) = identity_consistent.then(|| protein_by_accession.get(accession)).flatten() else {
                    orphan_features += 1;
                    if debug_failed_mappings && failed_debug_printed < 10 {
                        eprintln!(
                            "[ommverse] unmapped feature {file} {chrom}:{start}-{end} accession={accession:?} described={described_accession:?} fields={} label={:?}",
                            f.len(),
                            f.get(18).copied().unwrap_or("")
                        );
                        failed_debug_printed += 1;
                    }
                    return Ok(());
                };

                let source_db = f.get(13).copied().unwrap_or("").to_owned();
                let status_text = f.get(17).copied().unwrap_or("");
                let label = f.get(18).copied().unwrap_or("").to_owned();
                let description = aa_text.to_owned();
                // Feature tracks carry the UniProt/gene metadata that the
                // canonical bigPsl alignment deliberately does not.  Use it to
                // enrich canonical protein records without changing identity.
                let protein = &mut proteins[idx];
                if protein.gene_symbol.is_empty() { protein.gene_symbol = f.get(14).copied().unwrap_or("").to_owned(); }
                if protein.name.is_empty() { protein.name = f.get(20).copied().unwrap_or("").to_owned(); }
                if protein.aliases.is_empty() { protein.aliases = parse_aliases(f.get(21).copied().unwrap_or(""), ""); }
                protein.features.push(ProteinFeature {
                    kind, protein_range: parse_amino_acid_range(aa_text), label, description,
                    chromosome: chrom.to_owned(), genomic_start: start, genomic_end: end,
                    source_db, review_status: review_status(status_text),
                });
                feature_records += 1;
                Ok(())
            })?;
        }

        let report = BuildReport {
            transcripts: splice.transcripts.len(), mapped_proteins: proteins.len(), linked_transcripts: linked.len(),
            feature_records, unlinked_mapping_records, feature_records_without_protein: orphan_features, warnings,
        };
        let mut out = Self { assembly, source_root: root, genome_twobit: twobit, splice, proteins, report,
            interpro_entries: HashMap::new(), chromatin: Vec::new(), protein_binding: ProteinBindingUnion::default(), protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new() };
        progress("building lookup indexes", 0);
        out.reindex();
        progress("importing FANTOM5 chromatin", 0);
        out.ingest_fantom5_if_present()?;
        progress("importing ENCODE4 TF rPeaks", 0);
        out.ingest_encode4_binding_if_present(&mut progress)?;
        progress("build complete", out.protein_binding.regions.len());
        Ok(out)
    }

    /// Add assembly-relevant InterPro annotations to an existing Ommverse.
    ///
    /// `protein2ipr.dat.gz` is streamed and only records whose UniProt accession
    /// already exists in this Ommverse are retained. InterPro coordinates are
    /// 1-based inclusive and are converted to Ommverse's 0-based half-open
    /// protein coordinates. The stable IPR accession is stored in `label`, the
    /// human-readable entry name plus contributing member signature in
    /// `description`, and `source_db` is `InterPro`.
    pub fn ingest_interpro(
        &mut self,
        protein2ipr: impl AsRef<Path>,
        entry_list: impl AsRef<Path>,
        parent_child_tree: Option<&Path>,
    ) -> Result<InterProImportReport> {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};

        let entry_list = entry_list.as_ref();
        let mut entries = HashMap::<String, (ProteinFeatureKind, String)>::new();
        let reader = BufReader::new(File::open(entry_list)
            .with_context(|| format!("opening {}", entry_list.display()))?);
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') { continue; }
            let mut f = line.split('\t');
            let Some(ipr) = f.next().map(str::trim).filter(|x| !x.is_empty()) else { continue; };
            let kind_text = f.next().unwrap_or("").trim();
            let name = f.next().unwrap_or("").trim();
            let kind = interpro_feature_kind(kind_text);
            entries.insert(ipr.to_owned(), (kind, name.to_owned()));
            self.interpro_entries.entry(ipr.to_owned()).or_insert_with(|| InterProEntry {
                kind, name: name.to_owned(), parents: Vec::new(), children: Vec::new(),
            });
        }
        if let Some(tree) = parent_child_tree {
            let reader = BufReader::new(File::open(tree).with_context(|| format!("opening {}", tree.display()))?);
            let mut stack: Vec<String> = Vec::new();
            for line in reader.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') { continue; }
                let mut depth = 0usize;
                let bytes = line.as_bytes();
                while bytes.get(depth * 2..depth * 2 + 2) == Some(b"--") { depth += 1; }
                let body = &line[depth * 2..];
                let Some(ipr) = body.split("::").next().map(str::trim).filter(|x| x.starts_with("IPR")) else { continue; };
                stack.truncate(depth);
                if depth > 0 {
                    if let Some(parent) = stack.get(depth - 1).cloned() {
                        if let Some(entry) = self.interpro_entries.get_mut(ipr) {
                            if !entry.parents.contains(&parent) { entry.parents.push(parent.clone()); }
                        }
                        if let Some(entry) = self.interpro_entries.get_mut(&parent) {
                            if !entry.children.iter().any(|x| x == ipr) { entry.children.push(ipr.to_owned()); }
                        }
                    }
                }
                stack.push(ipr.to_owned());
            }
        }

        let mut report = InterProImportReport { entries_loaded: entries.len(), ..Default::default() };
        let path = protein2ipr.as_ref();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let gz = path.extension().is_some_and(|x| x == "gz");
        let input: Box<dyn Read> = if gz { Box::new(MultiGzDecoder::new(file)) } else { Box::new(file) };
        let reader = BufReader::with_capacity(1024 * 1024, input);

        for line in reader.lines() {
            let line = line?;
            report.records_streamed += 1;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 6 { report.malformed_records += 1; continue; }
            let accession = f[0].trim();
            let Some(&idx) = self.protein_by_accession.get(accession) else { continue; };
            report.matched_records += 1;

            let ipr = f[1].trim();
            let signature = f[3].trim();
            let start_1: u32 = match f[4].trim().parse() { Ok(v) if v > 0 => v, _ => { report.malformed_records += 1; continue; } };
            let end_1: u32 = match f[5].trim().parse() { Ok(v) if v >= start_1 => v, _ => { report.malformed_records += 1; continue; } };
            let (kind, entry_name) = entries.get(ipr)
                .cloned()
                .unwrap_or_else(|| (ProteinFeatureKind::Other, f[2].trim().to_owned()));
            if !entries.contains_key(ipr) { report.unknown_entries += 1; }

            let description = if signature.is_empty() {
                entry_name
            } else if entry_name.is_empty() {
                format!("member signature {signature}")
            } else {
                format!("{entry_name} [{signature}]")
            };
            let feature = ProteinFeature {
                kind,
                protein_range: Some((start_1 - 1, end_1)),
                label: ipr.to_owned(),
                description,
                chromosome: String::new(),
                genomic_start: 0,
                genomic_end: 0,
                source_db: "InterPro".to_owned(),
                review_status: ReviewStatus::Other,
            };
            if !self.proteins[idx].features.contains(&feature) {
                self.proteins[idx].features.push(feature);
                report.features_added += 1;
            } else {
                report.duplicate_features += 1;
            }
        }
        self.report.feature_records += report.features_added;
        Ok(report)
    }

    pub fn protein(&self, accession: &str) -> Option<&Protein> {
        self.protein_by_accession.get(accession).map(|&i| &self.proteins[i])
    }

    pub fn proteins_for_gene(&self, symbol: &str) -> impl Iterator<Item=&Protein> {
        self.proteins_by_gene.get(&symbol.to_ascii_lowercase()).into_iter().flatten().map(|&i| &self.proteins[i])
    }

    pub fn proteins_with_feature(&self, kind: ProteinFeatureKind) -> impl Iterator<Item=&Protein> {
        self.proteins.iter().filter(move |p| p.features.iter().any(|f| f.kind == kind))
    }

    /// Reconstruct the protein sequence from the UCSC twoBit genome and the
    /// linked GTF transcript. Protein sequence is deliberately not persisted.
    pub fn protein_sequence(&self, accession: &str) -> Result<IntToProt> {
        let protein = self.protein(accession).with_context(|| format!("unknown protein {accession}"))?;
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        self.protein_sequence_with_reader(protein, &mut genome)
    }

    fn protein_sequence_with_reader(&self, protein: &Protein, genome: &mut TwoBitReader) -> Result<IntToProt> {
        let &tx_id = protein.transcript_ids.first().context("protein has no linked transcript")?;
        let tx = &self.splice.transcripts[tx_id];
        let chr = self.splice.chr_names.get(tx.chr_id).context("transcript chromosome missing")?;
        let mut cdna = Vec::<u8>::with_capacity(tx.transcript_len());
        match tx.strand {
            Strand::Plus | Strand::Unknown => for exon in tx.exons() {
                let dna = genome.sequence(chr, exon.start, exon.end)?;
                cdna.extend_from_slice(dna.to_string(exon.len() as usize).as_bytes());
            },
            Strand::Minus => for exon in tx.exons().iter().rev() {
                let dna = genome.sequence(chr, exon.start, exon.end)?;
                let s = dna.to_string(exon.len() as usize);
                cdna.extend(s.bytes().rev().map(complement));
            },
        }
        let (start, end) = tx.cds_transcript_span().context("transcript has no usable CDS")?;
        if end > cdna.len() { bail!("CDS {}..{} exceeds reconstructed cDNA length {}", start, end, cdna.len()); }
        Ok(IntToDna::try_new(&cdna[start..end]).map_err(anyhow::Error::msg)?.translate())
    }

    /// Build a validated residue-level corpus for one curated protein feature.
    /// The split is deterministic and occurs at protein level, never residue level.
    pub fn training_corpus(&self, feature: ProteinFeatureKind) -> Result<ProteinTrainingCorpus> {
        self.training_corpus_with_flank_diagnostic(feature, 10)
    }

    /// As `training_corpus`, while also describing how often two context windows
    /// of `flank_width` residues would collide between adjacent features.
    pub fn training_corpus_with_flank_diagnostic(&self, feature: ProteinFeatureKind, flank_width: usize) -> Result<ProteinTrainingCorpus> {
        let candidates: Vec<&Protein> = self.proteins_with_feature(feature).collect();
        let mut report = TrainingCorpusReport { candidate_proteins: candidates.len(), diagnostic_flank_width: flank_width, ..Default::default() };
        let mut feature_lengths = Vec::<usize>::new();
        let mut inter_feature_gaps = Vec::<usize>::new();
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let mut examples = Vec::with_capacity(candidates.len());

        for protein in candidates {
            let sequence = match self.protein_sequence_with_reader(protein, &mut genome) {
                Ok(sequence) => sequence,
                Err(_) => { report.rejected_sequence += 1; continue; }
            };
            let mut ranges: Vec<(u32, u32)> = protein.features.iter()
                .filter(|f| f.kind == feature)
                .filter_map(|f| f.protein_range)
                .collect();
            ranges.sort_unstable();
            let selected_features = protein.features.iter().filter(|f| f.kind == feature).count();
            if ranges.len() != selected_features {
                report.rejected_missing_range += 1;
                continue;
            }
            if ranges.iter().any(|&(start, end)| start >= end || end as usize > sequence.len()) {
                report.rejected_out_of_bounds += 1;
                continue;
            }

            let mut truth = vec![false; sequence.len()];
            for &(start, end) in &ranges {
                truth[start as usize..end as usize].fill(true);
            }
            for &(start, end) in &ranges {
                feature_lengths.push((end - start) as usize);
            }
            for pair in ranges.windows(2) {
                let gap = pair[1].0.saturating_sub(pair[0].1) as usize;
                inter_feature_gaps.push(gap);
                if gap < flank_width.saturating_mul(2) { report.overlapping_flank_pairs += 1; }
            }
            report.reconstructed_proteins += 1;
            report.total_residues += sequence.len();
            report.feature_residues += truth.iter().filter(|&&x| x).count();
            report.feature_segments += ranges.len();
            examples.push(ProteinTrainingExample {
                accession: protein.accession.clone(),
                gene_symbol: protein.gene_symbol.clone(),
                sequence,
                truth,
                feature_segments: ranges.len(),
            });
        }

        report.feature_lengths = distribution_summary(&feature_lengths);
        report.inter_feature_gaps = distribution_summary(&inter_feature_gaps);

        examples.sort_by(|a, b| a.accession.cmp(&b.accession));
        let mut train = Vec::with_capacity((examples.len() + 1) / 2);
        let mut test = Vec::with_capacity(examples.len() / 2);
        for (i, example) in examples.into_iter().enumerate() {
            if i % 2 == 0 { train.push(example); } else { test.push(example); }
        }
        Ok(ProteinTrainingCorpus { feature, train, test, report })
    }

    /// Train a supervised four-state exact-amino-acid feature HMM.
    ///
    /// States are BACKGROUND, PRE_FEATURE, FEATURE and POST_FEATURE. PRE/POST
    /// are learned from residues immediately outside curated feature intervals.
    /// The held-out corpus is never consulted here.
    pub fn train_exact_aa_feature_model(&self, corpus: &ProteinTrainingCorpus, flank_width: usize) -> Result<ExactAaFeatureModel> {
        if flank_width == 0 { bail!("flank width must be greater than zero"); }
        let mut initial = [1.0f64; 4];
        let mut transition = [0.0f64; 16];
        let mut emission = [[1.0f64; 32]; 4];

        // Only biologically meaningful transitions receive a pseudocount.
        // POST->PRE permits two nearby features separated by a short loop.
        for (src, dst) in [
            (0,0),(0,1), (1,1),(1,2), (2,2),(2,3),
            (3,3),(3,0),(3,1),
        ] { transition[src * 4 + dst] = 1.0; }

        for example in &corpus.train {
            if example.sequence.is_empty() { continue; }
            let states = feature_context_states(&example.truth, flank_width);
            initial[states[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = states[pos];
                let code = example.sequence.get(pos).context("protein sequence position disappeared")?.code() as usize;
                emission[state][code] += 1.0;
                if pos > 0 { transition[states[pos - 1] * 4 + state] += 1.0; }
            }
        }
        normalize_array(&mut initial);
        for row in 0..4 {
            let sum: f64 = transition[row * 4..row * 4 + 4].iter().sum();
            if sum == 0.0 { bail!("feature HMM state {row} has no outgoing transitions"); }
            for col in 0..4 { transition[row * 4 + col] /= sum; }
        }
        for row in &mut emission { normalize_array(row); }
        Ok(ExactAaFeatureModel { feature: corpus.feature, flank_width, initial, transition, emission })
    }

    pub fn evaluate_exact_aa_feature_model(
        &self,
        model: &ExactAaFeatureModel,
        examples: &[ProteinTrainingExample],
    ) -> Result<FeatureEvaluation> {
        let hmm = model.hmm()?;
        let partials: Result<Vec<FeatureEvaluation>> = examples.par_iter().map(|example| {
            let observations = aa_observations(&example.sequence)?;
            let predicted: Vec<bool> = hmm.infer(&observations)?.viterbi().iter().map(|s| s.0 == 2).collect();
            let mut report = FeatureEvaluation::default();
            accumulate_evaluation(&mut report, &predicted, example);
            Ok(report)
        }).collect();
        Ok(partials?.into_iter().fold(FeatureEvaluation::default(), merge_feature_evaluation))
    }

    /// Apply a frozen model to every reconstructable protein that has no curated
    /// annotation of this feature. This is inference only; these proteins never
    /// contribute to model fitting.
    pub fn scan_unannotated_exact_aa(&self, model: &ExactAaFeatureModel) -> Result<UnannotatedScanReport> {
        let hmm = model.hmm()?;
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let mut report = UnannotatedScanReport::default();
        for protein in &self.proteins {
            if protein.features.iter().any(|f| f.kind == model.feature) { continue; }
            report.candidate_proteins += 1;
            let sequence = match self.protein_sequence_with_reader(protein, &mut genome) {
                Ok(x) if !x.is_empty() => x,
                _ => { report.rejected_sequence += 1; continue; }
            };
            report.reconstructed_proteins += 1;
            let observations = aa_observations(&sequence)?;
            let result = hmm.infer(&observations)?;
            let predicted: Vec<bool> = result.viterbi().iter().map(|s| s.0 == 2).collect();
            let segments = bool_segments(&predicted);
            if !segments.is_empty() {
                report.predicted_feature_proteins += 1;
                report.predicted_segments += segments.len();
                report.predicted_residues += predicted.iter().filter(|&&x| x).count();
            }
        }
        Ok(report)
    }

    pub fn train_chemistry_feature_model(&self, corpus: &ProteinTrainingCorpus, flank_width: usize) -> Result<ChemistryFeatureModel> {
        if flank_width == 0 { bail!("flank width must be greater than zero"); }
        let mut categories = Vec::<u16>::new();
        for example in &corpus.train {
            for pos in 0..example.sequence.len() {
                let bits = example.sequence.get(pos).context("protein sequence position disappeared")?.chemistry();
                categories.push(bits);
            }
        }
        categories.sort_unstable();
        categories.dedup();
        if categories.is_empty() { bail!("chemistry model has no training observations"); }
        let unknown = categories.len();
        let mut initial = [1.0f64; 4];
        let mut transition = [0.0f64; 16];
        let mut emission = vec![vec![1.0f64; categories.len() + 1]; 4];
        for (src, dst) in [(0,0),(0,1),(1,1),(1,2),(2,2),(2,3),(3,3),(3,0),(3,1)] {
            transition[src * 4 + dst] = 1.0;
        }
        for example in &corpus.train {
            if example.sequence.is_empty() { continue; }
            let states = feature_context_states(&example.truth, flank_width);
            initial[states[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = states[pos];
                let bits = example.sequence.get(pos).context("protein sequence position disappeared")?.chemistry();
                let category = categories.binary_search(&bits).unwrap_or(unknown);
                emission[state][category] += 1.0;
                if pos > 0 { transition[states[pos - 1] * 4 + state] += 1.0; }
            }
        }
        normalize_array(&mut initial);
        for row in 0..4 {
            let sum: f64 = transition[row * 4..row * 4 + 4].iter().sum();
            if sum == 0.0 { bail!("feature HMM state {row} has no outgoing transitions"); }
            for col in 0..4 { transition[row * 4 + col] /= sum; }
        }
        for row in &mut emission { normalize_slice(row); }
        Ok(ChemistryFeatureModel { feature: corpus.feature, flank_width, initial, transition, categories, emission })
    }

    pub fn evaluate_chemistry_feature_model(&self, model: &ChemistryFeatureModel, examples: &[ProteinTrainingExample]) -> Result<FeatureEvaluation> {
        let hmm = model.hmm()?;
        let partials: Result<Vec<FeatureEvaluation>> = examples.par_iter().map(|example| {
            let observations = chemistry_observations(&example.sequence, &model.categories)?;
            let predicted: Vec<bool> = hmm.infer(&observations)?.viterbi().iter().map(|s| s.0 == 2).collect();
            let mut report = FeatureEvaluation::default();
            accumulate_evaluation(&mut report, &predicted, example);
            Ok(report)
        }).collect();
        Ok(partials?.into_iter().fold(FeatureEvaluation::default(), merge_feature_evaluation))
    }

    /// Build a supervised membrane-architecture corpus.
    /// Signal peptide is deliberately handled by its independent N-terminal model.
    /// TM is explicit; cytoplasmic/extracellular annotations are expanded into
    /// short-loop, medium-loop, long-region and terminal latent states. These
    /// substates are collapsed back to their biological class during evaluation.
    pub fn topology_training_corpus(&self) -> Result<TopologyTrainingCorpus> {
        // Signal peptide is intentionally excluded here. It is already a strong
        // independent N-terminal detector and should act as upstream evidence rather
        // than compete at every residue of the membrane-architecture HMM.
        let kinds = [
            ProteinFeatureKind::Transmembrane,
            ProteinFeatureKind::Cytoplasmic, ProteinFeatureKind::Extracellular,
        ];
        let mut report = TopologyTrainingReport::default();
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let mut examples = Vec::new();

        for protein in &self.proteins {
            let selected: Vec<&ProteinFeature> = protein.features.iter()
                .filter(|f| kinds.contains(&f.kind)).collect();
            if selected.is_empty() { continue; }
            report.candidate_proteins += 1;
            if selected.iter().any(|f| f.protein_range.is_none()) {
                report.rejected_missing_range += 1;
                continue;
            }
            let sequence = match self.protein_sequence_with_reader(protein, &mut genome) {
                Ok(sequence) if !sequence.is_empty() => sequence,
                _ => { report.rejected_sequence += 1; continue; }
            };
            if selected.iter().any(|f| {
                let (start, end) = f.protein_range.unwrap();
                start >= end || end as usize > sequence.len()
            }) {
                report.rejected_out_of_bounds += 1;
                continue;
            }

            let mut truth = vec![TopologyState::Other.index(); sequence.len()];
            let mut priority = vec![0u8; sequence.len()];
            let mut conflict = false;
            for feature in selected {
                let (start, end) = feature.protein_range.unwrap();
                let state = match feature.kind {
                    ProteinFeatureKind::Transmembrane => TopologyState::Transmembrane,
                    ProteinFeatureKind::Cytoplasmic | ProteinFeatureKind::Extracellular =>
                        TopologyState::region_state(feature.kind, start as usize, end as usize, sequence.len()),
                    _ => unreachable!(),
                };
                // Membrane crossings outrank broad sidedness annotations.
                let p = if state == TopologyState::Transmembrane { 2 } else { 1 };
                for pos in start as usize..end as usize {
                    if p > priority[pos] {
                        truth[pos] = state.index();
                        priority[pos] = p;
                    } else if p == priority[pos] && truth[pos] != state.index() {
                        report.conflicting_residues += 1;
                        conflict = true;
                    }
                }
            }
            if conflict { continue; }
            report.reconstructed_proteins += 1;
            examples.push(TopologyTrainingExample {
                accession: protein.accession.clone(), sequence, truth,
            });
        }
        examples.sort_by(|a, b| a.accession.cmp(&b.accession));
        let mut train = Vec::with_capacity((examples.len() + 1) / 2);
        let mut test = Vec::with_capacity(examples.len() / 2);
        for (i, example) in examples.into_iter().enumerate() {
            if i % 2 == 0 { train.push(example); } else { test.push(example); }
        }
        Ok(TopologyTrainingCorpus { train, test, report })
    }

    pub fn train_exact_aa_topology_model(&self, corpus: &TopologyTrainingCorpus) -> Result<ExactAaTopologyModel> {
        let n = TopologyState::COUNT;
        let mut initial = vec![1.0f64; n];
        // Small pseudocount everywhere: observed biology dominates, but a transition
        // absent from mouse training is not made literally impossible in another species.
        let mut transition = vec![0.1f64; n * n];
        let mut emission = vec![vec![1.0f64; 32]; n];
        for example in &corpus.train {
            if example.sequence.is_empty() { continue; }
            initial[example.truth[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = example.truth[pos];
                let code = example.sequence.get(pos).context("protein sequence position disappeared")?.code() as usize;
                emission[state][code] += 1.0;
                if pos > 0 { transition[example.truth[pos - 1] * n + state] += 1.0; }
            }
        }
        normalize_slice(&mut initial);
        for row in 0..n { normalize_slice(&mut transition[row*n..(row+1)*n]); }
        for row in &mut emission { normalize_slice(row); }
        Ok(ExactAaTopologyModel { initial, transition, emission })
    }

    pub fn train_chemistry_topology_model(&self, corpus: &TopologyTrainingCorpus) -> Result<ChemistryTopologyModel> {
        let n = TopologyState::COUNT;
        let mut categories = Vec::<u16>::new();
        for example in &corpus.train {
            for pos in 0..example.sequence.len() {
                categories.push(example.sequence.get(pos).context("protein sequence position disappeared")?.chemistry());
            }
        }
        categories.sort_unstable(); categories.dedup();
        if categories.is_empty() { bail!("topology chemistry model has no training observations"); }
        let unknown = categories.len();
        let mut initial = vec![1.0f64; n];
        let mut transition = vec![0.1f64; n*n];
        let mut emission = vec![vec![1.0f64; categories.len()+1]; n];
        for example in &corpus.train {
            if example.sequence.is_empty() { continue; }
            initial[example.truth[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = example.truth[pos];
                let bits = example.sequence.get(pos).context("protein sequence position disappeared")?.chemistry();
                let obs = categories.binary_search(&bits).unwrap_or(unknown);
                emission[state][obs] += 1.0;
                if pos > 0 { transition[example.truth[pos-1]*n + state] += 1.0; }
            }
        }
        normalize_slice(&mut initial);
        for row in 0..n { normalize_slice(&mut transition[row*n..(row+1)*n]); }
        for row in &mut emission { normalize_slice(row); }
        Ok(ChemistryTopologyModel { categories, initial, transition, emission })
    }

    pub fn evaluate_exact_aa_topology_model(&self, model: &ExactAaTopologyModel, examples: &[TopologyTrainingExample]) -> Result<TopologyEvaluation> {
        let hmm = model.hmm()?;
        evaluate_topology(examples, |example| {
            let observations = aa_observations(&example.sequence)?;
            Ok(hmm.infer(&observations)?.viterbi().iter().map(|s| s.0).collect())
        })
    }

    pub fn evaluate_chemistry_topology_model(&self, model: &ChemistryTopologyModel, examples: &[TopologyTrainingExample]) -> Result<TopologyEvaluation> {
        let hmm = model.hmm()?;
        evaluate_topology(examples, |example| {
            let observations = chemistry_observations(&example.sequence, &model.categories)?;
            Ok(hmm.infer(&observations)?.viterbi().iter().map(|s| s.0).collect())
        })
    }

    /// Train every currently supported observation model for every feature class
    /// with a non-empty deterministic train/test split. One bad/empty feature does
    /// not abort the vault; it is recorded in the report instead.
    pub fn train_aa_model_vault(&self, flank_width: usize) -> Result<(AaModelVault, ModelVaultTrainingReport)> {
        let mut models = Vec::new();
        let mut report = ModelVaultTrainingReport::default();
        for feature in ProteinFeatureKind::ALL {
            // Chain is effectively whole-protein coverage and Conflict is a
            // curation discrepancy, not a sequence feature. Neither belongs in
            // this biological feature-modelling experiment.
            if matches!(feature, ProteinFeatureKind::Chain | ProteinFeatureKind::Conflict | ProteinFeatureKind::SignalPeptide) { continue; }
            report.feature_classes_considered += 1;
            let corpus = match self.training_corpus_with_flank_diagnostic(feature, flank_width) {
                Ok(corpus) => corpus,
                Err(err) => { report.skipped.push((feature, err.to_string())); continue; }
            };
            if corpus.train.is_empty() || corpus.test.is_empty() || corpus.report.feature_residues == 0 {
                report.skipped.push((feature, format!("insufficient corpus: train {}, test {}, feature residues {}", corpus.train.len(), corpus.test.len(), corpus.report.feature_residues)));
                continue;
            }
            let exact = self.train_exact_aa_feature_model(&corpus, flank_width)?;
            let exact_eval = self.evaluate_exact_aa_feature_model(&exact, &corpus.test)?;
            models.push(ProteinFeatureModel::ExactAa { model: exact, evaluation: exact_eval });
            let chemistry = self.train_chemistry_feature_model(&corpus, flank_width)?;
            let chemistry_eval = self.evaluate_chemistry_feature_model(&chemistry, &corpus.test)?;
            models.push(ProteinFeatureModel::Chemistry { model: chemistry, evaluation: chemistry_eval });
            report.feature_classes_trained += 1;
            report.models_trained += 2;
        }
        // Membrane architecture is deliberately joint. TM competes with multiple
        // latent cytoplasmic/extracellular contexts (short loop, medium loop, long
        // region). Terminality and signal peptide are external/positional evidence,
        // not states in this membrane-architecture model.
        let topology = self.topology_training_corpus()?;
        report.topology_train_proteins = topology.train.len();
        report.topology_test_proteins = topology.test.len();
        report.topology_conflicting_residues = topology.report.conflicting_residues;
        if !topology.train.is_empty() && !topology.test.is_empty() {
            let exact = self.train_exact_aa_topology_model(&topology)?;
            let exact_eval = self.evaluate_exact_aa_topology_model(&exact, &topology.test)?;
            models.push(ProteinFeatureModel::TopologyExactAa { model: exact, evaluation: exact_eval });
            let chemistry = self.train_chemistry_topology_model(&topology)?;
            let chemistry_eval = self.evaluate_chemistry_topology_model(&chemistry, &topology.test)?;
            models.push(ProteinFeatureModel::TopologyChemistry { model: chemistry, evaluation: chemistry_eval });
            report.models_trained += 2;
        }

        let vault = AaModelVault {
            metadata: ModelVaultMetadata {
                assembly: self.assembly.clone(),
                source: "Ommverse/UCSC".to_owned(),
                ommverse_format_version: OMMVERSE_FORMAT_VERSION,
                vault_format_version: AA_MODEL_VAULT_FORMAT_VERSION,
                flank_width,
            },
            models,
        };
        Ok((vault, report))
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut file = File::create(path.as_ref())?;
        file.write_all(MAGIC)?;
        file.write_all(&OMMVERSE_FORMAT_VERSION.to_le_bytes())?;
        bincode::serialize_into(file, self)?;
        Ok(())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let mut magic = [0u8; 4]; file.read_exact(&mut magic)?;
        if &magic != MAGIC { bail!("not an Ommverse index"); }
        let mut version = [0u8; 4]; file.read_exact(&mut version)?;
        let version = u32::from_le_bytes(version);
        let mut out = match version {
            OMMVERSE_FORMAT_VERSION => bincode::deserialize_from(file)?,
            4 => {
                let old: OmmverseV4 = bincode::deserialize_from(file)?;
                // v4 embedded all rPeaks.  Drop that experiment-level payload on
                // load; rebuilding v5 reconstructs the compact union from cache.
                Self { assembly: old.assembly, source_root: old.source_root, genome_twobit: old.genome_twobit,
                    splice: old.splice, proteins: old.proteins, report: old.report, interpro_entries: old.interpro_entries,
                    chromatin: old.chromatin, protein_binding: ProteinBindingUnion::default(), protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new() }
            }
            3 => {
                let old: OmmverseV3 = bincode::deserialize_from(file)?;
                Self { assembly: old.assembly, source_root: old.source_root, genome_twobit: old.genome_twobit,
                    splice: old.splice, proteins: old.proteins, report: old.report, interpro_entries: old.interpro_entries,
                    chromatin: old.chromatin, protein_binding: ProteinBindingUnion::default(), protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new() }
            }
            2 => {
                let old: OmmverseV2 = bincode::deserialize_from(file)?;
                Self { assembly: old.assembly, source_root: old.source_root, genome_twobit: old.genome_twobit,
                    splice: old.splice, proteins: old.proteins, report: old.report, interpro_entries: old.interpro_entries,
                    chromatin: Vec::new(), protein_binding: ProteinBindingUnion::default(), protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new() }
            }
            1 => {
                let old: OmmverseV1 = bincode::deserialize_from(file)?;
                Self {
                    assembly: old.assembly, source_root: old.source_root, genome_twobit: old.genome_twobit,
                    splice: old.splice, proteins: old.proteins, report: old.report, interpro_entries: HashMap::new(), chromatin: Vec::new(), protein_binding: ProteinBindingUnion::default(),
                    protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new(),
                }
            }
            _ => bail!("unsupported Ommverse format version {version}"),
        };
        out.reindex();
        Ok(out)
    }

    fn reindex(&mut self) {
        self.protein_by_accession.clear(); self.proteins_by_gene.clear();
        for (i, p) in self.proteins.iter().enumerate() {
            self.protein_by_accession.insert(p.accession.clone(), i);
            if !p.gene_symbol.is_empty() { self.proteins_by_gene.entry(p.gene_symbol.to_ascii_lowercase()).or_default().push(i); }
        }
    }
}

fn read_all_bigbed(path: &Path, mut visit: impl FnMut(&str, u32, u32, &str) -> Result<()>) -> Result<()> {
    let mut bb = BigBedRead::open_file(path).with_context(|| format!("opening bigBed {}", path.display()))?;
    let chroms: Vec<(String, u32)> = bb.chroms().iter().map(|c| (c.name.clone(), c.length)).collect();
    for (chrom, len) in chroms {
        let entries = bb.get_interval(&chrom, 0, len)?;
        for entry in entries {
            let entry = entry?;
            visit(&chrom, entry.start, entry.end, &entry.rest)?;
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MappingSchema { EnsGene, UnipAliSwissprot }

fn find_swissprot_mapping_bigbed(dir: &Path) -> Result<(PathBuf, MappingSchema)> {
    // unipAliSwissprot.bb is the canonical Swiss-Prot alignment track exposed
    // by UCSC's UniProt trackDb.  Hashed ensGene_*.swissprot.bb files are
    // auxiliary mapping products and multiple variants can coexist in one
    // UniProt release, so they must not override the canonical track.
    let canonical = dir.join("unipAliSwissprot.bb");
    if canonical.is_file() { return Ok((canonical, MappingSchema::UnipAliSwissprot)); }

    // Keep compatibility with older source trees that contain only one of the
    // historical ensGene mapping products.  Never guess when several exist.
    let mut hits = Vec::new();
    for e in fs::read_dir(dir)? {
        let p = e?.path();
        let name = p.file_name().and_then(|x| x.to_str()).unwrap_or("");
        if name.starts_with("ensGene_") && name.ends_with(".swissprot.bb") { hits.push(p); }
    }
    hits.sort();
    match hits.len() {
        1 => Ok((hits.remove(0), MappingSchema::EnsGene)),
        n if n > 1 => bail!("multiple ensGene_*.swissprot.bb mapping bigBeds found in {} and canonical unipAliSwissprot.bb is absent; refusing to guess: {:?}", dir.display(), hits),
        _ => bail!("no Swiss-Prot mapping bigBed found in {} (tried unipAliSwissprot.bb and ensGene_*.swissprot.bb)", dir.display()),
    }
}

fn strip_version(s: &str) -> &str { s.rsplit_once('.').filter(|(_, v)| v.chars().all(|c| c.is_ascii_digit())).map(|(a, _)| a).unwrap_or(s) }
fn nonempty(s: &str) -> Option<String> { let s=s.trim(); (!s.is_empty()).then(|| s.to_owned()) }
fn review_status(s: &str) -> ReviewStatus { if s.contains("Swiss-Prot") || s.eq_ignore_ascii_case("swissprot") { ReviewStatus::SwissProt } else if s.contains("TrEMBL") || s.eq_ignore_ascii_case("trembl") { ReviewStatus::Trembl } else { ReviewStatus::Other } }
fn parse_aliases(a: &str, b: &str) -> Vec<String> { let mut v=Vec::new(); for x in a.split(|c:char| c==',' || c==';' || c.is_whitespace()).chain(b.split(|c:char| c==',' || c==';' || c.is_whitespace())) { let x=x.trim(); if !x.is_empty() && !v.iter().any(|y| y==x) { v.push(x.to_owned()); } } v }
fn complement(b: u8) -> u8 { match b.to_ascii_uppercase() { b'A'=>b'T', b'C'=>b'G', b'G'=>b'C', b'T'=>b'A', x=>x } }

fn protein_accession_from_feature_text(text: &str) -> Option<&str> {
    text.rsplit_once(" on protein ")
        .map(|(_, accession)| accession.trim())
        .filter(|accession| !accession.is_empty())
}

fn parse_amino_acid_range(text: &str) -> Option<(u32,u32)> {
    let tail = text.strip_prefix("amino acids ").or_else(|| text.strip_prefix("amino acid "))?;
    let token = tail.split_whitespace().next()?;
    let (a,b) = token.split_once('-').map(|(a,b)|(a,b)).unwrap_or((token,token));
    let start: u32=a.parse().ok()?; let end:u32=b.parse().ok()?;
    (start > 0 && end >= start).then_some((start-1,end))
}


impl AaModelVault {
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut file = File::create(path)?;
        file.write_all(AA_MODEL_VAULT_MAGIC)?;
        file.write_all(&AA_MODEL_VAULT_FORMAT_VERSION.to_le_bytes())?;
        bincode::serialize_into(file, self)?;
        Ok(())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let mut file = File::open(path)?;
        let mut magic = [0u8; 4]; file.read_exact(&mut magic)?;
        if &magic != AA_MODEL_VAULT_MAGIC { bail!("not an Ommverse AA model vault"); }
        let mut version = [0u8; 4]; file.read_exact(&mut version)?;
        let version = u32::from_le_bytes(version);
        if version != AA_MODEL_VAULT_FORMAT_VERSION { bail!("unsupported AA model vault format version {version}"); }
        Ok(bincode::deserialize_from(file)?)
    }
}

impl ExactAaTopologyModel {
    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self.emission.iter().map(|row| CategoricalEmission::new(row.clone()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(self.initial.clone(), self.transition.clone(), emissions)?)
    }
}

impl ChemistryTopologyModel {
    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self.emission.iter().map(|row| CategoricalEmission::new(row.clone()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(self.initial.clone(), self.transition.clone(), emissions)?)
    }
}

impl ChemistryFeatureModel {
    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self.emission.iter().map(|row| CategoricalEmission::new(row.clone())).collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(self.initial.to_vec(), self.transition.to_vec(), emissions)?)
    }
}

impl ExactAaFeatureModel {
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let file = File::create(path)?;
        bincode::serialize_into(file, self)?;
        Ok(())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        Ok(bincode::deserialize_from(file)?)
    }

    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self.emission.iter()
            .map(|row| CategoricalEmission::new(row.to_vec()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(self.initial.to_vec(), self.transition.to_vec(), emissions)?)
    }
}

fn chemistry_observations(sequence: &IntToProt, categories: &[u16]) -> Result<Vec<usize>> {
    let unknown = categories.len();
    (0..sequence.len()).map(|pos| {
        let bits = sequence.get(pos).context("protein sequence position disappeared")?.chemistry();
        Ok(categories.binary_search(&bits).unwrap_or(unknown))
    }).collect()
}

fn topology_state_from_index(index: usize) -> Option<TopologyState> {
    Some(match index {
        0 => TopologyState::Other,
        1 => TopologyState::Transmembrane,
        2 => TopologyState::CytoShortLoop,
        3 => TopologyState::CytoMediumLoop,
        4 => TopologyState::CytoLongRegion,
        5 => TopologyState::ExtraShortLoop,
        6 => TopologyState::ExtraMediumLoop,
        7 => TopologyState::ExtraLongRegion,
        _ => return None,
    })
}

fn evaluate_topology(
    examples: &[TopologyTrainingExample],
    predict: impl Fn(&TopologyTrainingExample) -> Result<Vec<usize>> + Sync,
) -> Result<TopologyEvaluation> {
    let partials: Result<Vec<TopologyEvaluation>> = examples.par_iter().map(|example| {
        let predicted = predict(example)?;
        if predicted.len() != example.truth.len() { bail!("topology prediction length mismatch"); }
        let mut reports = vec![FeatureEvaluation::default(); TopologyState::BIOLOGICAL.len()];
        let mut truth_residues = vec![0usize; TopologyState::COUNT];
        let mut predicted_residues = vec![0usize; TopologyState::COUNT];
        let mut confusion = vec![0usize; TopologyState::COUNT * TopologyState::COUNT];
        for (&truth, &pred) in example.truth.iter().zip(&predicted) {
            if truth < TopologyState::COUNT && pred < TopologyState::COUNT {
                truth_residues[truth] += 1;
                predicted_residues[pred] += 1;
                confusion[truth * TopologyState::COUNT + pred] += 1;
            }
        }
        for (i, feature) in TopologyState::BIOLOGICAL.iter().enumerate() {
            let truth_mask: Vec<bool> = example.truth.iter().map(|&x| topology_state_from_index(x).and_then(TopologyState::feature) == Some(*feature)).collect();
            let pred_mask: Vec<bool> = predicted.iter().map(|&x| topology_state_from_index(x).and_then(TopologyState::feature) == Some(*feature)).collect();
            let r = &mut reports[i];
            r.proteins += 1; r.residues += truth_mask.len();
            for (&truth, &pred) in truth_mask.iter().zip(&pred_mask) {
                match (truth, pred) {
                    (true,true) => r.true_positive += 1, (false,true) => r.false_positive += 1,
                    (false,false) => r.true_negative += 1, (true,false) => r.false_negative += 1,
                }
            }
            let truth_segments = bool_segments(&truth_mask);
            let pred_segments = bool_segments(&pred_mask);
            r.truth_segments += truth_segments.len();
            r.predicted_segments += pred_segments.len();
            r.recovered_segments += truth_segments.iter().filter(|truth| pred_segments.iter().any(|pred| overlap_fraction(**truth, *pred) >= 0.5)).count();
        }
        Ok(TopologyEvaluation { states: reports, latent_truth_residues: truth_residues, latent_predicted_residues: predicted_residues, latent_confusion: confusion })
    }).collect();
    Ok(partials?.into_iter().fold(empty_topology_evaluation(), merge_topology_evaluation))
}

fn empty_topology_evaluation() -> TopologyEvaluation {
    TopologyEvaluation {
        states: vec![FeatureEvaluation::default(); TopologyState::BIOLOGICAL.len()],
        latent_truth_residues: vec![0; TopologyState::COUNT],
        latent_predicted_residues: vec![0; TopologyState::COUNT],
        latent_confusion: vec![0; TopologyState::COUNT * TopologyState::COUNT],
    }
}

fn merge_topology_evaluation(mut a: TopologyEvaluation, b: TopologyEvaluation) -> TopologyEvaluation {
    for (dst, src) in a.states.iter_mut().zip(b.states) { *dst = merge_feature_evaluation(dst.clone(), src); }
    for (dst, src) in a.latent_truth_residues.iter_mut().zip(b.latent_truth_residues) { *dst += src; }
    for (dst, src) in a.latent_predicted_residues.iter_mut().zip(b.latent_predicted_residues) { *dst += src; }
    for (dst, src) in a.latent_confusion.iter_mut().zip(b.latent_confusion) { *dst += src; }
    a
}

fn merge_feature_evaluation(mut a: FeatureEvaluation, b: FeatureEvaluation) -> FeatureEvaluation {
    a.proteins += b.proteins; a.residues += b.residues;
    a.true_positive += b.true_positive; a.false_positive += b.false_positive;
    a.true_negative += b.true_negative; a.false_negative += b.false_negative;
    a.truth_segments += b.truth_segments; a.predicted_segments += b.predicted_segments;
    a.recovered_segments += b.recovered_segments;
    a
}

fn accumulate_evaluation(report: &mut FeatureEvaluation, predicted: &[bool], example: &ProteinTrainingExample) {
    report.proteins += 1;
    report.residues += predicted.len();
    for (&truth, &pred) in example.truth.iter().zip(predicted) {
        match (truth, pred) {
            (true, true) => report.true_positive += 1,
            (false, true) => report.false_positive += 1,
            (false, false) => report.true_negative += 1,
            (true, false) => report.false_negative += 1,
        }
    }
    let truth_segments = bool_segments(&example.truth);
    let pred_segments = bool_segments(predicted);
    report.truth_segments += truth_segments.len();
    report.predicted_segments += pred_segments.len();
    report.recovered_segments += truth_segments.iter().filter(|truth| pred_segments.iter().any(|pred| overlap_fraction(**truth, *pred) >= 0.5)).count();
}

fn normalize_slice(values: &mut [f64]) {
    let sum: f64 = values.iter().sum();
    if sum > 0.0 { for value in values { *value /= sum; } }
}

fn distribution_summary(values: &[usize]) -> DistributionSummary {
    if values.is_empty() { return DistributionSummary::default(); }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let count = sorted.len();
    let mean = sorted.iter().map(|&x| x as f64).sum::<f64>() / count as f64;
    let variance = sorted.iter().map(|&x| { let d = x as f64 - mean; d * d }).sum::<f64>() / count as f64;
    let percentile = |q: f64| -> f64 {
        if count == 1 { return sorted[0] as f64; }
        let pos = q * (count - 1) as f64;
        let lo = pos.floor() as usize;
        let hi = pos.ceil() as usize;
        let frac = pos - lo as f64;
        sorted[lo] as f64 * (1.0 - frac) + sorted[hi] as f64 * frac
    };
    DistributionSummary {
        count, mean, median: percentile(0.5), sd: variance.sqrt(),
        q1: percentile(0.25), q3: percentile(0.75),
        min: sorted[0], max: sorted[count - 1],
    }
}

fn aa_observations(sequence: &IntToProt) -> Result<Vec<usize>> {
    (0..sequence.len()).map(|pos| {
        sequence.get(pos).map(|aa| aa.code() as usize)
            .context("protein sequence position disappeared")
    }).collect()
}

/// Convert binary feature truth into four supervised states:
/// 0 BACKGROUND, 1 PRE_FEATURE, 2 FEATURE, 3 POST_FEATURE.
///
/// For short gaps between features, each background residue is assigned to the
/// nearest boundary; ties go to PRE_FEATURE. This avoids order-dependent flank
/// overwrites while retaining both sides of nearby features.
fn feature_context_states(truth: &[bool], flank_width: usize) -> Vec<usize> {
    let segments = bool_segments(truth);
    let mut states = vec![0usize; truth.len()];
    for &(start, end) in &segments { states[start..end].fill(2); }
    for pos in 0..truth.len() {
        if states[pos] == 2 { continue; }
        let prev = segments.iter().filter(|(_, end)| *end <= pos).map(|(_, end)| pos + 1 - *end).min();
        let next = segments.iter().filter(|(start, _)| *start > pos).map(|(start, _)| *start - pos).min();
        let prev = prev.filter(|&d| d <= flank_width);
        let next = next.filter(|&d| d <= flank_width);
        states[pos] = match (prev, next) {
            (Some(a), Some(b)) => if b <= a { 1 } else { 3 },
            (None, Some(_)) => 1,
            (Some(_), None) => 3,
            (None, None) => 0,
        };
    }
    states
}

fn bool_segments(values: &[bool]) -> Vec<(usize, usize)> {
    let mut out = Vec::new(); let mut start = None;
    for (i, &value) in values.iter().enumerate() {
        match (start, value) {
            (None, true) => start = Some(i),
            (Some(s), false) => { out.push((s, i)); start = None; }
            _ => {}
        }
    }
    if let Some(s) = start { out.push((s, values.len())); }
    out
}

fn overlap_fraction(a: (usize, usize), b: (usize, usize)) -> f64 {
    let overlap = a.1.min(b.1).saturating_sub(a.0.max(b.0));
    if overlap == 0 { 0.0 } else { overlap as f64 / (a.1 - a.0) as f64 }
}

fn ratio(num: usize, den: usize) -> f64 { if den == 0 { 0.0 } else { num as f64 / den as f64 } }
fn normalize_array<const N: usize>(values: &mut [f64; N]) { let sum: f64 = values.iter().sum(); for x in values { *x /= sum; } }

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn strips_only_numeric_versions() { assert_eq!(strip_version("ENSMUST1.5"), "ENSMUST1"); assert_eq!(strip_version("Q9EST3-1"), "Q9EST3-1"); }
    #[test] fn parses_feature_protein_accession() { assert_eq!(protein_accession_from_feature_text("amino acids 485-507 on protein Q5GH67"), Some("Q5GH67")); assert_eq!(protein_accession_from_feature_text("not a protein range"), None); }
    #[test] fn parses_ucsc_uniprot_coordinates() { assert_eq!(parse_amino_acid_range("amino acids 485-507 on protein Q5GH67"), Some((484,507))); assert_eq!(parse_amino_acid_range("amino acids 42 on protein X"), Some((41,42))); assert_eq!(parse_amino_acid_range("amino acid 197 on protein X"), Some((196,197))); }
    #[test] fn summarizes_geometry_distribution() {
        let s = distribution_summary(&[1, 2, 3, 4, 5]);
        assert_eq!(s.count, 5);
        assert_eq!(s.mean, 3.0);
        assert_eq!(s.median, 3.0);
        assert_eq!(s.q1, 2.0);
        assert_eq!(s.q3, 4.0);
        assert_eq!(s.min, 1);
        assert_eq!(s.max, 5);
        assert!((s.sd - 2.0f64.sqrt()).abs() < 1e-12);
    }
    #[test] fn builds_pre_and_post_feature_context() {
        let truth = [false, false, false, true, true, false, false, false];
        assert_eq!(feature_context_states(&truth, 2), vec![0,1,1,2,2,3,3,0]);
    }
    #[test] fn splits_short_inter_feature_gap_by_nearest_boundary() {
        let truth = [false, true, true, false, false, false, true, true, false];
        assert_eq!(feature_context_states(&truth, 3), vec![1,2,2,3,1,1,2,2,3]);
    }
}
