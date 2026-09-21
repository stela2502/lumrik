use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use ommverse::{AaModelVault, ExactAaFeatureModel, Ommverse, ProteinFeatureKind, ProteinFeatureModel};
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Parser)]
#[command(name="ommverse", about="Build and explore Lumrik's integrated biological reference model")]
struct Cli { #[command(subcommand)] command: Command }
#[derive(Subcommand)]
enum Command {
    Build { #[arg(long)] reference: PathBuf, #[arg(long)] out: PathBuf, #[arg(long)] debug_failed_mappings: bool },
    Protein { #[arg(long)] index: PathBuf, accession: String, #[arg(long)] sequence: bool },
    Gene { #[arg(long)] index: PathBuf, symbol: String },
    /// Build and validate a supervised protein-feature dataset from an Ommverse index.
    ///
    /// Complete protein sequences are reconstructed from the reference genome and
    /// transcript annotation. Curated feature intervals are projected onto those
    /// sequences and proteins with inconsistent annotations are rejected. The
    /// surviving proteins are split deterministically 50:50 at protein level.
    /// Residues outside the selected feature remain in the examples as background.
    /// This command validates/describes the corpus; it does not train a model.
    TrainingCorpus {
        #[arg(long)] index: PathBuf,
        #[arg(long)] feature: String,
        /// Context width used only to report collisions between adjacent feature flanks.
        #[arg(long, default_value_t = 10)] flank_width: usize,
    },

    /// Train and evaluate a feature HMM from exact amino-acid identity.
    ///
    /// The selected curated feature supplies supervised residue labels. Training
    /// uses only the deterministic training half; the other half is held out for
    /// evaluation. Observations are exact IntToProt residue identities: no chemistry,
    /// substitution matrix, or protein-language representation is supplied.
    ///
    /// The HMM has four biological states: BACKGROUND -> PRE_FEATURE -> FEATURE ->
    /// POST_FEATURE -> BACKGROUND. PRE/POST describe the immediately adjacent
    /// residues outside each curated feature. --flank-width controls that context
    /// window (default 10 aa). With --scan-unannotated the frozen model is then
    /// applied to proteins lacking the selected annotation; they never train it.
    TrainExactAa {
        #[arg(long)] index: PathBuf,
        #[arg(long)] feature: String,
        #[arg(long)] out: PathBuf,
        /// Number of residues immediately outside each curated feature assigned to PRE/POST states.
        #[arg(long, default_value_t = 10)] flank_width: usize,
        #[arg(long)] scan_unannotated: bool,
    },

    /// Train a portable collection of all currently supported amino-acid feature HMMs.
    ///
    /// Every curated feature class with a usable deterministic train/test split is
    /// trained independently. At present the vault contains both exact-AA and
    /// intrinsic-chemistry HMMs. Empty/unsupported feature classes are reported and
    /// skipped rather than aborting the run.
    TrainModelVault {
        #[arg(long)] index: PathBuf,
        #[arg(long)] out: PathBuf,
        #[arg(long, default_value_t = 10)] flank_width: usize,
    },

    /// Show provenance and held-out metrics stored with a model vault.
    ModelVaultInfo { #[arg(long)] vault: PathBuf },

    /// Apply every frozen model in a vault to the matching curated feature in another Ommverse index.
    ///
    /// Nothing is fitted or updated: the vault is loaded exactly as trained and each
    /// model is evaluated against the target index's curated annotations. This is the
    /// cross-reference / cross-species transfer test.
    EvaluateModelVault {
        #[arg(long)] index: PathBuf,
        #[arg(long)] vault: PathBuf,
    },

    /// Apply a previously trained exact-AA feature model to unannotated proteins.
    ScanExactAa { #[arg(long)] index: PathBuf, #[arg(long)] model: PathBuf },
}

fn percent(n: usize, d: usize) -> f64 {
    if d == 0 { 0.0 } else { 100.0 * n as f64 / d as f64 }
}

fn print_distribution(label: &str, s: &ommverse::DistributionSummary) {
    println!("  {label}:");
    println!("    n:                       {}", s.count);
    println!("    mean:                    {:.2} aa", s.mean);
    println!("    median:                  {:.2} aa", s.median);
    println!("    SD:                      {:.2} aa", s.sd);
    println!("    Q1 / Q3:                 {:.2} / {:.2} aa", s.q1, s.q3);
    println!("    min / max:               {} / {} aa", s.min, s.max);
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Build { reference, out, debug_failed_mappings } => {
            let model=Ommverse::build_ucsc_with_debug(&reference, debug_failed_mappings)?; model.save(&out)?;
            println!("Ommverse {}", model.assembly);
            println!("  transcripts:             {}", model.report.transcripts);
            println!("  proteins:                {}", model.report.mapped_proteins);
            println!("  linked transcripts:      {}", model.report.linked_transcripts);
            println!("  protein features:        {}", model.report.feature_records);
            println!("  unlinked mappings:       {}", model.report.unlinked_mapping_records);
            println!("  orphan features:         {}", model.report.feature_records_without_protein);
            println!("  transmembrane proteins:  {}", model.proteins_with_feature(ProteinFeatureKind::Transmembrane).count());
            for w in &model.report.warnings { eprintln!("warning: {w}"); }
            println!("  saved: {}", out.display());
        }
        Command::Protein { index, accession, sequence } => {
            let omm=Ommverse::load(index)?; let p=omm.protein(&accession).with_context(|| format!("protein {accession} not found"))?;
            println!("{}  {}", p.accession, p.name); println!("  gene: {}", p.gene_symbol); println!("  entry: {}", p.entry_name);
            println!("  transcripts: {:?}", p.transcript_ids); println!("  features: {}", p.features.len());
            for f in &p.features { println!("    {:?}\t{:?}\t{}", f.kind, f.protein_range, f.description); }
            if sequence { println!("  sequence: {}", omm.protein_sequence(&accession)?); }
        }
        Command::Gene { index, symbol } => {
            let omm=Ommverse::load(index)?; let proteins:Vec<_>=omm.proteins_for_gene(&symbol).collect();
            if proteins.is_empty() { anyhow::bail!("gene {symbol} not found"); }
            println!("{}", symbol); for p in proteins { println!("  {}\t{}\t{} features", p.accession,p.name,p.features.len()); }
        }
        Command::TrainingCorpus { index, feature, flank_width } => {
            let omm = Ommverse::load(index)?;
            let feature = ProteinFeatureKind::from_str(&feature)?;
            let corpus = omm.training_corpus_with_flank_diagnostic(feature, flank_width)?;
            let r = &corpus.report;
            println!("Ommverse supervised protein-feature corpus");
            println!("  feature:                   {:?}", feature);
            println!("  sequence:                  complete reconstructed proteins");
            println!("  labels:                    curated Ommverse feature intervals");
            println!("  background:                all residues outside the feature");
            println!("  split:                     deterministic protein-level 50:50");
            println!("  model training:            none (corpus validation only)");
            println!("  candidate proteins:        {}", r.candidate_proteins);
            println!("  reconstructed:             {}", r.reconstructed_proteins);
            println!("  rejected sequence:         {}", r.rejected_sequence);
            println!("  rejected missing range:    {}", r.rejected_missing_range);
            println!("  rejected out of bounds:    {}", r.rejected_out_of_bounds);
            println!("  training proteins:         {}", corpus.train.len());
            println!("  test proteins:             {}", corpus.test.len());
            println!("  total residues:            {}", r.total_residues);
            println!("  feature residues:          {}", r.feature_residues);
            println!("  feature segments:          {}", r.feature_segments);
            println!();
            println!("Feature geometry");
            print_distribution("feature length", &r.feature_lengths);
            println!();
            println!("Adjacent feature geometry");
            print_distribution("inter-feature gap", &r.inter_feature_gaps);
            println!("  flank diagnostic width:    {} aa per side", r.diagnostic_flank_width);
            println!("  overlapping flank pairs:   {} / {} ({:.2}%)",
                r.overlapping_flank_pairs,
                r.inter_feature_gaps.count,
                percent(r.overlapping_flank_pairs, r.inter_feature_gaps.count));
        }
        Command::TrainExactAa { index, feature, out, flank_width, scan_unannotated } => {
            let omm = Ommverse::load(index)?;
            let feature = ProteinFeatureKind::from_str(&feature)?;
            println!("Ommverse protein-feature training");
            println!("  feature:          {:?}", feature);
            println!("  observation:      exact amino-acid identity");
            println!("  states:           BACKGROUND / PRE_FEATURE / FEATURE / POST_FEATURE");
            println!("  flank width:      {} aa", flank_width);
            println!("  supervision:      curated Ommverse protein-feature coordinates");
            println!("  split:            deterministic protein-level 50:50");
            println!("  held-out labels:  never used for fitting");
            println!("  unannotated scan: {}", if scan_unannotated { "enabled; inference only" } else { "disabled" });
            println!();
            println!("Preparing supervised corpus...");
            let corpus = omm.training_corpus(feature)?;
            println!("  accepted proteins: {} (train {}, test {})", corpus.train.len() + corpus.test.len(), corpus.train.len(), corpus.test.len());
            println!("Training four-state exact-AA HMM...");
            let model = omm.train_exact_aa_feature_model(&corpus, flank_width)?;
            let evaluation = omm.evaluate_exact_aa_feature_model(&model, &corpus.test)?;
            model.save(&out)?;
            println!("Evaluating held-out proteins...");
            println!("Exact-AA {:?} model", feature);
            println!("  training proteins:       {}", corpus.train.len());
            println!("  held-out proteins:       {}", corpus.test.len());
            println!("  held-out residues:       {}", evaluation.residues);
            println!("  precision:               {:.4}", evaluation.precision());
            println!("  recall:                  {:.4}", evaluation.recall());
            println!("  specificity:             {:.4}", evaluation.specificity());
            println!("  F1:                      {:.4}", evaluation.f1());
            println!("  truth segments:          {}", evaluation.truth_segments);
            println!("  predicted segments:      {}", evaluation.predicted_segments);
            println!("  segments recovered >=50%:{} ({:.4})", evaluation.recovered_segments, evaluation.segment_recall());
            println!("  mean truth segment len:  {:.2} aa", ratio_usize(evaluation.true_positive + evaluation.false_negative, evaluation.truth_segments));
            println!("  mean predicted seg len:  {:.2} aa", ratio_usize(evaluation.true_positive + evaluation.false_positive, evaluation.predicted_segments));
            println!("  initial:                 {:?}", model.initial);
            println!("  transition rows:         BACKGROUND / PRE_FEATURE / FEATURE / POST_FEATURE");
            for row in 0..4 { println!("    {:?}", &model.transition[row * 4..row * 4 + 4]); }
            println!("  learned AA emissions (BACKGROUND / PRE_FEATURE / FEATURE / POST_FEATURE):");
            for residue in b"ACDEFGHIKLMNPQRSTVWY" {
                // Canonical IntToProt codes are not alphabetical; ask IntToProt for the code.
                let aa = int_to_prot::IntToProt::new([*residue]).get(0).unwrap();
                let c = aa.code() as usize;
                println!("    {}  {:.6}  {:.6}  {:.6}  {:.6}", *residue as char, model.emission[0][c], model.emission[1][c], model.emission[2][c], model.emission[3][c]);
            }
            println!("  saved: {}", out.display());
            if scan_unannotated {
                print_unannotated(&omm.scan_unannotated_exact_aa(&model)?);
            }
        }
        Command::TrainModelVault { index, out, flank_width } => {
            let omm = Ommverse::load(&index)?;
            println!("Ommverse AA model-vault training");
            println!("  assembly:          {}", omm.assembly);
            println!("  observations:      exact-AA + intrinsic chemistry");
            println!("  feature classes:   all curated Ommverse feature kinds");
            println!("  split:             deterministic protein-level 50:50 per feature");
            println!("  flank width:       {} aa", flank_width);
            println!("  held-out labels:   never used for fitting");
            println!();
            let (vault, report) = omm.train_aa_model_vault(flank_width)?;
            vault.save(&out)?;
            println!("Model vault complete");
            println!("  feature classes considered: {}", report.feature_classes_considered);
            println!("  feature classes trained:    {}", report.feature_classes_trained);
            println!("  models trained:             {}", report.models_trained);
            println!("  models stored:              {}", vault.models.len());
            for (feature, reason) in &report.skipped { println!("  skipped {:?}: {}", feature, reason); }
            println!("  saved: {}", out.display());
        }
        Command::ModelVaultInfo { vault } => {
            let vault = AaModelVault::load(vault)?;
            println!("Ommverse AA model vault");
            println!("  assembly:                  {}", vault.metadata.assembly);
            println!("  source:                    {}", vault.metadata.source);
            println!("  Ommverse format version:   {}", vault.metadata.ommverse_format_version);
            println!("  vault format version:      {}", vault.metadata.vault_format_version);
            println!("  flank width:               {} aa", vault.metadata.flank_width);
            println!("  models:                    {}", vault.models.len());
            for entry in &vault.models {
                match entry {
                    ProteinFeatureModel::ExactAa { model, evaluation } => println!("    {:?} / exact-AA: precision {:.4}, recall {:.4}, F1 {:.4}", model.feature, evaluation.precision(), evaluation.recall(), evaluation.f1()),
                    ProteinFeatureModel::Chemistry { model, evaluation } => println!("    {:?} / chemistry: precision {:.4}, recall {:.4}, F1 {:.4}", model.feature, evaluation.precision(), evaluation.recall(), evaluation.f1()),
                }
            }
        }
        Command::EvaluateModelVault { index, vault } => {
            let omm = Ommverse::load(&index)?;
            let vault = AaModelVault::load(&vault)?;
            println!("Ommverse model-vault transfer evaluation");
            println!("  trained assembly:          {}", vault.metadata.assembly);
            println!("  target assembly:           {}", omm.assembly);
            println!("  models:                    {}", vault.models.len());
            println!("  fitting on target:         none");
            println!("  target labels:             evaluation only");
            println!();

            for entry in &vault.models {
                match entry {
                    ProteinFeatureModel::ExactAa { model, evaluation: source } => {
                        let corpus = omm.training_corpus(model.feature)?;
                        let target_examples: Vec<_> = corpus.train.iter().chain(corpus.test.iter()).cloned().collect();
                        let target = omm.evaluate_exact_aa_feature_model(model, &target_examples)?;
                        print_transfer_evaluation(model.feature, "exact-AA", source, &target);
                    }
                    ProteinFeatureModel::Chemistry { model, evaluation: source } => {
                        let corpus = omm.training_corpus(model.feature)?;
                        let target_examples: Vec<_> = corpus.train.iter().chain(corpus.test.iter()).cloned().collect();
                        let target = omm.evaluate_chemistry_feature_model(model, &target_examples)?;
                        print_transfer_evaluation(model.feature, "chemistry", source, &target);
                    }
                }
            }
        }
        Command::ScanExactAa { index, model } => {
            let omm = Ommverse::load(index)?;
            let model = ExactAaFeatureModel::load(model)?;
            print_unannotated(&omm.scan_unannotated_exact_aa(&model)?);
        }
    }
    Ok(())
}

fn print_unannotated(r: &ommverse::UnannotatedScanReport) {
    println!("Unannotated-protein scan");
    println!("  candidate proteins:      {}", r.candidate_proteins);
    println!("  reconstructed:           {}", r.reconstructed_proteins);
    println!("  rejected sequence:       {}", r.rejected_sequence);
    println!("  predicted proteins:      {}", r.predicted_feature_proteins);
    println!("  predicted segments:      {}", r.predicted_segments);
    println!("  predicted residues:      {}", r.predicted_residues);
}

fn ratio_usize(num: usize, den: usize) -> f64 { if den == 0 { 0.0 } else { num as f64 / den as f64 } }

fn print_transfer_evaluation(
    feature: ProteinFeatureKind,
    observation: &str,
    source: &ommverse::FeatureEvaluation,
    target: &ommverse::FeatureEvaluation,
) {
    println!("{:?} / {}", feature, observation);
    println!("  source held-out: precision {:.4}, recall {:.4}, F1 {:.4}", source.precision(), source.recall(), source.f1());
    println!("  target corpus:   {} proteins / {} residues / {} truth segments", target.proteins, target.residues, target.truth_segments);
    println!("  target:          precision {:.4}, recall {:.4}, specificity {:.4}, F1 {:.4}", target.precision(), target.recall(), target.specificity(), target.f1());
    println!("  target counts:   TP {} / FP {} / TN {} / FN {}", target.true_positive, target.false_positive, target.true_negative, target.false_negative);
    println!("  segments:        predicted {}, recovered {} / {} ({:.4})", target.predicted_segments, target.recovered_segments, target.truth_segments, target.segment_recall());
    println!();
}
