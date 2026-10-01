
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReviewStatus {
    SwissProt,
    Trembl,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ProteinFeatureKind {
    Transmembrane,
    Domain,
    SignalPeptide,
    Cytoplasmic,
    Extracellular,
    ModifiedResidue,
    Disulfide,
    Repeat,
    Chain,
    Conflict,
    Interest,
    Mutagenesis,
    SpliceVariant,
    Structure,
    Other,
    ActiveSite,
    BindingSite,
    ConservedSite,
    ProteinFamily,
    HomologousSuperfamily,
    PtmSite,
}

impl ProteinFeatureKind {
    pub const ALL: [Self; 15] = [
        Self::Transmembrane,
        Self::Domain,
        Self::SignalPeptide,
        Self::Cytoplasmic,
        Self::Extracellular,
        Self::ModifiedResidue,
        Self::Disulfide,
        Self::Repeat,
        Self::Chain,
        Self::Conflict,
        Self::Interest,
        Self::Mutagenesis,
        Self::SpliceVariant,
        Self::Structure,
        Self::Other,
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
    /// Protein identity aliases (UniProt, ENSP, entry names and source-provided protein IDs).
    pub identifiers: Vec<String>,
    /// Gene-level aliases carried by protein source metadata.
    pub aliases: Vec<String>,
    pub ensembl_gene: Option<String>,
    pub ensembl_protein: Option<String>,
    pub transcript_ids: Vec<TranscriptId>,
    pub review_status: ReviewStatus,
    pub features: Vec<ProteinFeature>,
}

impl Identifiable for Protein {
    fn aliases(&self) -> &[String] { &self.identifiers }
}

impl Plottable for ProteinFeature {
    fn axis(&self) -> Axis { Axis::Proteomic }

    fn blocks(&self) -> Vec<RefBlock> {
        self.protein_range
            .filter(|(start, end)| start < end)
            .map(|(start, end)| vec![RefBlock::new(start, end)])
            .unwrap_or_default()
    }
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
    pub fn precision(&self) -> f64 {
        ratio(self.true_positive, self.true_positive + self.false_positive)
    }
    pub fn recall(&self) -> f64 {
        ratio(self.true_positive, self.true_positive + self.false_negative)
    }
    pub fn specificity(&self) -> f64 {
        ratio(self.true_negative, self.true_negative + self.false_positive)
    }
    pub fn f1(&self) -> f64 {
        let p = self.precision();
        let r = self.recall();
        if p + r == 0.0 {
            0.0
        } else {
            2.0 * p * r / (p + r)
        }
    }
    pub fn segment_recall(&self) -> f64 {
        ratio(self.recovered_segments, self.truth_segments)
    }
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
        Self::Other,
        Self::Transmembrane,
        Self::CytoShortLoop,
        Self::CytoMediumLoop,
        Self::CytoLongRegion,
        Self::ExtraShortLoop,
        Self::ExtraMediumLoop,
        Self::ExtraLongRegion,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::Other => "OTHER",
            Self::Transmembrane => "TM",
            Self::CytoShortLoop => "CYTO_SHORT",
            Self::CytoMediumLoop => "CYTO_MEDIUM",
            Self::CytoLongRegion => "CYTO_LONG",
            Self::ExtraShortLoop => "EXTRA_SHORT",
            Self::ExtraMediumLoop => "EXTRA_MEDIUM",
            Self::ExtraLongRegion => "EXTRA_LONG",
        }
    }
    fn index(self) -> usize {
        self as usize
    }
    fn feature(self) -> Option<ProteinFeatureKind> {
        match self {
            Self::Other => None,
            Self::Transmembrane => Some(ProteinFeatureKind::Transmembrane),
            Self::CytoShortLoop | Self::CytoMediumLoop | Self::CytoLongRegion => {
                Some(ProteinFeatureKind::Cytoplasmic)
            }
            Self::ExtraShortLoop | Self::ExtraMediumLoop | Self::ExtraLongRegion => {
                Some(ProteinFeatureKind::Extracellular)
            }
        }
    }
    fn region_state(
        feature: ProteinFeatureKind,
        start: usize,
        end: usize,
        _protein_len: usize,
    ) -> Self {
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
    ExactAa {
        model: ExactAaFeatureModel,
        evaluation: FeatureEvaluation,
    },
    Chemistry {
        model: ChemistryFeatureModel,
        evaluation: FeatureEvaluation,
    },
    /// Opinionated membrane-architecture HMM. Cytoplasmic/extracellular sequence
    /// is represented by latent loop/long/terminal substates and collapsed on output.
    TopologyExactAa {
        model: ExactAaTopologyModel,
        evaluation: TopologyEvaluation,
    },
    TopologyChemistry {
        model: ChemistryTopologyModel,
        evaluation: TopologyEvaluation,
    },
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

impl Identifiable for InterProEntry {
    fn aliases(&self) -> &[String] {
        // The owning InterPro HashMap retains the IPR accession; the entry's
        // human-readable name is the alias stored on the object itself.
        std::slice::from_ref(&self.name)
    }
}

impl Connected for InterProEntry {
    fn connections(&self) -> Vec<Connection<'_>> {
        let mut out = Vec::with_capacity(self.parents.len() + self.children.len());
        out.extend(self.parents.iter().map(|target| Connection { target, relation: "parent" }));
        out.extend(self.children.iter().map(|target| Connection { target, relation: "child" }));
        out
    }
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

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV7 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
    interpro_entries: HashMap<String, InterProEntry>,
    chromatin: Vec<ChromatinElement>,
    protein_binding: ProteinBindingUnion,
    ctcf: CtcfArchitecture,
    experimental_loops: ExperimentalLoopArchitecture,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV6 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
    interpro_entries: HashMap<String, InterProEntry>,
    chromatin: Vec<ChromatinElement>,
    protein_binding: ProteinBindingUnion,
    ctcf: CtcfArchitecture,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV5 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
    interpro_entries: HashMap<String, InterProEntry>,
    chromatin: Vec<ChromatinElement>,
    protein_binding: ProteinBindingUnion,
}

/// A merged ENCODE CTCF rPeak anchor. Multiple overlapping CTCF rPeaks are
/// collapsed into one candidate anchor; experiment/biosample detail stays in
/// the cached ENCODE bigBed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct CtcfAnchor {
    pub chromosome_id: u16,
    pub region: RefBlock,
    pub source_peak_count: u32,
}

impl Plottable for CtcfAnchor {
    fn axis(&self) -> Axis { Axis::Genomic }
    fn blocks(&self) -> Vec<RefBlock> { vec![self.region] }
}

/// The interval between two adjacent candidate CTCF anchors. This is a
/// reference architectural unit, not a claim that a physical chromatin loop
/// has been observed in a particular cell type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct CtcfDomain {
    pub chromosome_id: u16,
    pub region: RefBlock,
    pub left_anchor: u32,
    pub right_anchor: u32,
}

impl Plottable for CtcfDomain {
    fn axis(&self) -> Axis { Axis::Genomic }
    fn blocks(&self) -> Vec<RefBlock> { vec![self.region] }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CtcfArchitecture {
    pub chromosomes: Vec<String>,
    pub anchors: Vec<CtcfAnchor>,
    pub domains: Vec<CtcfDomain>,
    pub source_peak_count: usize,
    pub source: String,
    pub source_url: String,
}

/// A 25-kb genomic endpoint bin supported by one or more independent loop
/// experiments. Multiple callers/resolutions from the same experiment count once.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExperimentalLoopAnchor {
    pub chromosome_id: u16,
    pub region: RefBlock,
    pub experiment_count: u32,
}

impl Plottable for ExperimentalLoopAnchor {
    fn axis(&self) -> Axis { Axis::Genomic }
    fn blocks(&self) -> Vec<RefBlock> { vec![self.region] }
}

/// A loop observed in one or more independent experiments. Endpoints are the
/// compact anchor bins above; Ommverse stores recurrence, not per-experiment evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExperimentalLoop {
    pub left_anchor: u32,
    pub right_anchor: u32,
    pub experiment_count: u32,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ExperimentalLoopArchitecture {
    pub chromosomes: Vec<String>,
    pub anchors: Vec<ExperimentalLoopAnchor>,
    pub loops: Vec<ExperimentalLoop>,
    pub experiment_count: usize,
    pub source_loop_count: usize,
    pub anchor_bin_size: u32,
    pub source: String,
    pub source_url: String,
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

impl Plottable for ProteinBindingRegion {
    fn axis(&self) -> Axis { Axis::Genomic }
    fn blocks(&self) -> Vec<RefBlock> { vec![self.region] }
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
    chromosome_id: u16,
    region: RefBlock,
    score: u16,
    factor_id: u16,
    peak_id: u32,
    observed_experiments: u16,
    assayed_experiments: u16,
    ccre_id: u32,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct ProteinBindingIndexV4 {
    chromosomes: Vec<String>,
    factors: Vec<String>,
    ccres: Vec<String>,
    peaks: Vec<ProteinBindingPeakV4>,
    source: String,
    source_url: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
struct OmmverseV4 {
    assembly: String,
    source_root: PathBuf,
    genome_twobit: PathBuf,
    splice: SpliceIndex,
    proteins: Vec<Protein>,
    report: BuildReport,
    interpro_entries: HashMap<String, InterProEntry>,
    chromatin: Vec<ChromatinElement>,
    protein_binding: ProteinBindingIndexV4,
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

impl Identifiable for ChromatinElement {
    fn aliases(&self) -> &[String] { std::slice::from_ref(&self.name) }
}

impl Plottable for ChromatinElement {
    fn axis(&self) -> Axis { Axis::Genomic }
    fn blocks(&self) -> Vec<RefBlock> {
        if self.blocks.is_empty() { vec![self.region] } else { self.blocks.clone() }
    }
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
    /// Candidate CTCF-bounded architectural units derived from ENCODE4 rPeaks.
    pub ctcf: CtcfArchitecture,
    /// Compact recurrence map from experimentally observed chromatin loops.
    pub experimental_loops: ExperimentalLoopArchitecture,
    /// Resolved biological sources and retrieval recipes used for this build.
    pub sources: Option<sources::SourceManifest>,
    #[serde(skip)]
    protein_by_accession: HashMap<String, usize>,
    #[serde(skip)]
    proteins_by_gene: HashMap<String, Vec<usize>>,
    /// Internal TranscriptId -> ProteinId adjacency, rebuilt on load.
    #[serde(skip)]
    proteins_by_transcript: HashMap<TranscriptId, Vec<usize>>,
    /// Case-insensitive gene symbols, aliases and stable accessions -> internal GeneId.
    #[serde(skip)]
    gene_by_name: HashMap<String, GeneId>,
}

