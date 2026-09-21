//! Ommverse: Lumrik's integrated genome-to-protein biological reference model.
//!
//! Ommverse deliberately separates biological identity from source formats.
//! UCSC GTF, twoBit and UniProt bigBed files are import formats; callers see
//! genes, transcripts, proteins and protein features.

use anyhow::{bail, Context, Result};
use bigtools::BigBedRead;
use gtf_splice_index::{IdNameKeys, SpliceIndex, Strand, TranscriptId};
use int_to_dna::{IntToDna, TwoBitReader};
use int_to_prot::IntToProt;
use hmm::{CategoricalEmission, Hmm};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 4] = b"OMM1";
pub const OMMVERSE_FORMAT_VERSION: u32 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReviewStatus { SwissProt, Trembl, Other }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ProteinFeatureKind {
    Transmembrane, Domain, SignalPeptide, Cytoplasmic, Extracellular,
    ModifiedResidue, Disulfide, Repeat, Chain, Conflict, Interest,
    Mutagenesis, SpliceVariant, Structure, Other,
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

pub const AA_MODEL_VAULT_FORMAT_VERSION: u32 = 1;
const AA_MODEL_VAULT_MAGIC: &[u8; 4] = b"AAV1";

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
pub struct Ommverse {
    pub assembly: String,
    pub source_root: PathBuf,
    pub genome_twobit: PathBuf,
    pub splice: SpliceIndex,
    pub proteins: Vec<Protein>,
    pub report: BuildReport,
    #[serde(skip)]
    protein_by_accession: HashMap<String, usize>,
    #[serde(skip)]
    proteins_by_gene: HashMap<String, Vec<usize>>,
}

impl Ommverse {
    /// Build an Ommverse directly from the directory layout produced by
    /// scripts/download_ucsc_reference.sh.
    pub fn build_ucsc(root: impl AsRef<Path>) -> Result<Self> {
        Self::build_ucsc_with_debug(root, false)
    }

    /// Build from UCSC and optionally report a bounded sample of feature records
    /// that could not be attached to a mapped protein.
    pub fn build_ucsc_with_debug(root: impl AsRef<Path>, debug_failed_mappings: bool) -> Result<Self> {
        let root = root.as_ref().canonicalize().with_context(|| format!("reference root {}", root.as_ref().display()))?;
        let assembly = root.file_name().and_then(|x| x.to_str()).context("reference root has no assembly name")?.to_owned();
        let gtf = root.join("Genes").join(format!("{assembly}.knownGene.gtf.gz"));
        let twobit = root.join("Genome").join(format!("{assembly}.2bit"));
        let protein_dir = root.join("Protein");
        for required in [&gtf, &twobit, &protein_dir] {
            if !required.exists() { bail!("required UCSC reference component missing: {}", required.display()); }
        }

        let splice = SpliceIndex::from_path(&gtf, 100_000, IdNameKeys::default())
            .with_context(|| format!("building splice index from {}", gtf.display()))?;
        let mut tx_by_stable = HashMap::<String, TranscriptId>::new();
        for tx in &splice.transcripts {
            for name in &tx.names {
                tx_by_stable.entry(strip_version(name).to_owned()).or_insert(tx.id);
            }
        }

        let mapping = find_mapping_bigbed(&protein_dir, "ensGene_", ".swissprot.bb")?;
        let mut proteins = Vec::<Protein>::new();
        let mut protein_by_accession = HashMap::<String, usize>::new();
        let mut linked = HashSet::new();
        let mut unlinked_mapping_records = 0usize;

        read_all_bigbed(&mapping, |_, _, _, rest| {
            let f: Vec<&str> = rest.split('\t').collect();
            // UCSC ensGene UniProt mapping is bigGenePred plus UniProt fields.
            // The stable identifiers occupy the final three columns.
            if f.len() < 41 { return Ok(()); }
            let accession = f[0].trim();
            let ensembl_tx = f.last().copied().unwrap_or("").trim();
            let Some(&tx_id) = tx_by_stable.get(strip_version(ensembl_tx)) else {
                unlinked_mapping_records += 1;
                return Ok(());
            };
            linked.insert(tx_id);
            let idx = *protein_by_accession.entry(accession.to_owned()).or_insert_with(|| {
                let idx = proteins.len();
                proteins.push(Protein {
                    accession: accession.to_owned(),
                    entry_name: f.get(22).copied().unwrap_or("").to_owned(),
                    review_status: review_status(f.get(23).copied().unwrap_or("")),
                    name: f.get(26).copied().unwrap_or("").to_owned(),
                    gene_symbol: f.get(27).copied().unwrap_or("").to_owned(),
                    aliases: parse_aliases(f.get(30).copied().unwrap_or(""), f.get(31).copied().unwrap_or("")),
                    ensembl_gene: nonempty(f.get(f.len().saturating_sub(3)).copied().unwrap_or("")),
                    ensembl_protein: nonempty(f.get(f.len().saturating_sub(2)).copied().unwrap_or("")),
                    transcript_ids: Vec::new(),
                    features: Vec::new(),
                });
                idx
            });
            if !proteins[idx].transcript_ids.contains(&tx_id) { proteins[idx].transcript_ids.push(tx_id); }
            Ok(())
        })?;

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
                proteins[idx].features.push(ProteinFeature {
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
            protein_by_accession: HashMap::new(), proteins_by_gene: HashMap::new() };
        out.reindex();
        Ok(out)
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
        let mut report = FeatureEvaluation::default();
        for example in examples {
            let observations = aa_observations(&example.sequence)?;
            let result = hmm.infer(&observations)?;
            let predicted: Vec<bool> = result.viterbi().iter().map(|s| s.0 == 2).collect();
            report.proteins += 1;
            report.residues += predicted.len();
            for (&truth, &pred) in example.truth.iter().zip(&predicted) {
                match (truth, pred) {
                    (true, true) => report.true_positive += 1,
                    (false, true) => report.false_positive += 1,
                    (false, false) => report.true_negative += 1,
                    (true, false) => report.false_negative += 1,
                }
            }
            let truth_segments = bool_segments(&example.truth);
            let pred_segments = bool_segments(&predicted);
            report.truth_segments += truth_segments.len();
            report.predicted_segments += pred_segments.len();
            report.recovered_segments += truth_segments.iter().filter(|truth| {
                pred_segments.iter().any(|pred| overlap_fraction(**truth, *pred) >= 0.5)
            }).count();
        }
        Ok(report)
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
        let mut report = FeatureEvaluation::default();
        for example in examples {
            let observations = chemistry_observations(&example.sequence, &model.categories)?;
            accumulate_evaluation(&mut report, &hmm.infer(&observations)?.viterbi().iter().map(|s| s.0 == 2).collect::<Vec<_>>(), example);
        }
        Ok(report)
    }

    /// Train every currently supported observation model for every feature class
    /// with a non-empty deterministic train/test split. One bad/empty feature does
    /// not abort the vault; it is recorded in the report instead.
    pub fn train_aa_model_vault(&self, flank_width: usize) -> Result<(AaModelVault, ModelVaultTrainingReport)> {
        let mut models = Vec::new();
        let mut report = ModelVaultTrainingReport::default();
        for feature in ProteinFeatureKind::ALL {
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
        if version != OMMVERSE_FORMAT_VERSION { bail!("unsupported Ommverse format version {version}"); }
        let mut out: Self = bincode::deserialize_from(file)?;
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

fn find_mapping_bigbed(dir: &Path, prefix: &str, suffix: &str) -> Result<PathBuf> {
    let mut hits = Vec::new();
    for e in fs::read_dir(dir)? {
        let p = e?.path();
        let name = p.file_name().and_then(|x| x.to_str()).unwrap_or("");
        if name.starts_with(prefix) && name.ends_with(suffix) { hits.push(p); }
    }
    hits.sort();
    match hits.len() {
        0 => bail!("no {prefix}*{suffix} mapping bigBed found in {}", dir.display()),
        1 => Ok(hits.remove(0)),
        _ => bail!("multiple {prefix}*{suffix} mapping bigBeds found in {}; refusing to guess: {:?}", dir.display(), hits),
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
fn normalize_pair(values: &mut [f64; 2]) { let sum = values[0] + values[1]; values[0] /= sum; values[1] /= sum; }
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
