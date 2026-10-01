impl Ommverse {

    /// Build a validated residue-level corpus for one curated protein feature.
    /// The split is deterministic and occurs at protein level, never residue level.
    pub fn training_corpus(&self, feature: ProteinFeatureKind) -> Result<ProteinTrainingCorpus> {
        self.training_corpus_with_flank_diagnostic(feature, 10)
    }

    /// As `training_corpus`, while also describing how often two context windows
    /// of `flank_width` residues would collide between adjacent features.
    pub fn training_corpus_with_flank_diagnostic(
        &self,
        feature: ProteinFeatureKind,
        flank_width: usize,
    ) -> Result<ProteinTrainingCorpus> {
        let candidates: Vec<&Protein> = self.proteins_with_feature(feature).collect();
        let mut report = TrainingCorpusReport {
            candidate_proteins: candidates.len(),
            diagnostic_flank_width: flank_width,
            ..Default::default()
        };
        let mut feature_lengths = Vec::<usize>::new();
        let mut inter_feature_gaps = Vec::<usize>::new();
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let mut examples = Vec::with_capacity(candidates.len());

        for protein in candidates {
            let sequence = match self.protein_sequence_with_reader(protein, &mut genome) {
                Ok(sequence) => sequence,
                Err(_) => {
                    report.rejected_sequence += 1;
                    continue;
                }
            };
            let mut ranges: Vec<(u32, u32)> = protein
                .features
                .iter()
                .filter(|f| f.kind == feature)
                .filter_map(|f| f.protein_range)
                .collect();
            ranges.sort_unstable();
            let selected_features = protein
                .features
                .iter()
                .filter(|f| f.kind == feature)
                .count();
            if ranges.len() != selected_features {
                report.rejected_missing_range += 1;
                continue;
            }
            if ranges
                .iter()
                .any(|&(start, end)| start >= end || end as usize > sequence.len())
            {
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
                if gap < flank_width.saturating_mul(2) {
                    report.overlapping_flank_pairs += 1;
                }
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
            if i % 2 == 0 {
                train.push(example);
            } else {
                test.push(example);
            }
        }
        Ok(ProteinTrainingCorpus {
            feature,
            train,
            test,
            report,
        })
    }

    /// Train a supervised four-state exact-amino-acid feature HMM.
    ///
    /// States are BACKGROUND, PRE_FEATURE, FEATURE and POST_FEATURE. PRE/POST
    /// are learned from residues immediately outside curated feature intervals.
    /// The held-out corpus is never consulted here.
    pub fn train_exact_aa_feature_model(
        &self,
        corpus: &ProteinTrainingCorpus,
        flank_width: usize,
    ) -> Result<ExactAaFeatureModel> {
        if flank_width == 0 {
            bail!("flank width must be greater than zero");
        }
        let mut initial = [1.0f64; 4];
        let mut transition = [0.0f64; 16];
        let mut emission = [[1.0f64; 32]; 4];

        // Only biologically meaningful transitions receive a pseudocount.
        // POST->PRE permits two nearby features separated by a short loop.
        for (src, dst) in [
            (0, 0),
            (0, 1),
            (1, 1),
            (1, 2),
            (2, 2),
            (2, 3),
            (3, 3),
            (3, 0),
            (3, 1),
        ] {
            transition[src * 4 + dst] = 1.0;
        }

        for example in &corpus.train {
            if example.sequence.is_empty() {
                continue;
            }
            let states = feature_context_states(&example.truth, flank_width);
            initial[states[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = states[pos];
                let code = example
                    .sequence
                    .get(pos)
                    .context("protein sequence position disappeared")?
                    .code() as usize;
                emission[state][code] += 1.0;
                if pos > 0 {
                    transition[states[pos - 1] * 4 + state] += 1.0;
                }
            }
        }
        normalize_array(&mut initial);
        for row in 0..4 {
            let sum: f64 = transition[row * 4..row * 4 + 4].iter().sum();
            if sum == 0.0 {
                bail!("feature HMM state {row} has no outgoing transitions");
            }
            for col in 0..4 {
                transition[row * 4 + col] /= sum;
            }
        }
        for row in &mut emission {
            normalize_array(row);
        }
        Ok(ExactAaFeatureModel {
            feature: corpus.feature,
            flank_width,
            initial,
            transition,
            emission,
        })
    }

    pub fn evaluate_exact_aa_feature_model(
        &self,
        model: &ExactAaFeatureModel,
        examples: &[ProteinTrainingExample],
    ) -> Result<FeatureEvaluation> {
        let hmm = model.hmm()?;
        let partials: Result<Vec<FeatureEvaluation>> = examples
            .par_iter()
            .map(|example| {
                let observations = aa_observations(&example.sequence)?;
                let predicted: Vec<bool> = hmm
                    .infer(&observations)?
                    .viterbi()
                    .iter()
                    .map(|s| s.0 == 2)
                    .collect();
                let mut report = FeatureEvaluation::default();
                accumulate_evaluation(&mut report, &predicted, example);
                Ok(report)
            })
            .collect();
        Ok(partials?
            .into_iter()
            .fold(FeatureEvaluation::default(), merge_feature_evaluation))
    }

    /// Apply a frozen model to every reconstructable protein that has no curated
    /// annotation of this feature. This is inference only; these proteins never
    /// contribute to model fitting.
    pub fn scan_unannotated_exact_aa(
        &self,
        model: &ExactAaFeatureModel,
    ) -> Result<UnannotatedScanReport> {
        let hmm = model.hmm()?;
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let mut report = UnannotatedScanReport::default();
        for protein in &self.proteins {
            if protein.features.iter().any(|f| f.kind == model.feature) {
                continue;
            }
            report.candidate_proteins += 1;
            let sequence = match self.protein_sequence_with_reader(protein, &mut genome) {
                Ok(x) if !x.is_empty() => x,
                _ => {
                    report.rejected_sequence += 1;
                    continue;
                }
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

    pub fn train_chemistry_feature_model(
        &self,
        corpus: &ProteinTrainingCorpus,
        flank_width: usize,
    ) -> Result<ChemistryFeatureModel> {
        if flank_width == 0 {
            bail!("flank width must be greater than zero");
        }
        let mut categories = Vec::<u16>::new();
        for example in &corpus.train {
            for pos in 0..example.sequence.len() {
                let bits = example
                    .sequence
                    .get(pos)
                    .context("protein sequence position disappeared")?
                    .chemistry();
                categories.push(bits);
            }
        }
        categories.sort_unstable();
        categories.dedup();
        if categories.is_empty() {
            bail!("chemistry model has no training observations");
        }
        let unknown = categories.len();
        let mut initial = [1.0f64; 4];
        let mut transition = [0.0f64; 16];
        let mut emission = vec![vec![1.0f64; categories.len() + 1]; 4];
        for (src, dst) in [
            (0, 0),
            (0, 1),
            (1, 1),
            (1, 2),
            (2, 2),
            (2, 3),
            (3, 3),
            (3, 0),
            (3, 1),
        ] {
            transition[src * 4 + dst] = 1.0;
        }
        for example in &corpus.train {
            if example.sequence.is_empty() {
                continue;
            }
            let states = feature_context_states(&example.truth, flank_width);
            initial[states[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = states[pos];
                let bits = example
                    .sequence
                    .get(pos)
                    .context("protein sequence position disappeared")?
                    .chemistry();
                let category = categories.binary_search(&bits).unwrap_or(unknown);
                emission[state][category] += 1.0;
                if pos > 0 {
                    transition[states[pos - 1] * 4 + state] += 1.0;
                }
            }
        }
        normalize_array(&mut initial);
        for row in 0..4 {
            let sum: f64 = transition[row * 4..row * 4 + 4].iter().sum();
            if sum == 0.0 {
                bail!("feature HMM state {row} has no outgoing transitions");
            }
            for col in 0..4 {
                transition[row * 4 + col] /= sum;
            }
        }
        for row in &mut emission {
            normalize_slice(row);
        }
        Ok(ChemistryFeatureModel {
            feature: corpus.feature,
            flank_width,
            initial,
            transition,
            categories,
            emission,
        })
    }

    pub fn evaluate_chemistry_feature_model(
        &self,
        model: &ChemistryFeatureModel,
        examples: &[ProteinTrainingExample],
    ) -> Result<FeatureEvaluation> {
        let hmm = model.hmm()?;
        let partials: Result<Vec<FeatureEvaluation>> = examples
            .par_iter()
            .map(|example| {
                let observations = chemistry_observations(&example.sequence, &model.categories)?;
                let predicted: Vec<bool> = hmm
                    .infer(&observations)?
                    .viterbi()
                    .iter()
                    .map(|s| s.0 == 2)
                    .collect();
                let mut report = FeatureEvaluation::default();
                accumulate_evaluation(&mut report, &predicted, example);
                Ok(report)
            })
            .collect();
        Ok(partials?
            .into_iter()
            .fold(FeatureEvaluation::default(), merge_feature_evaluation))
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
            ProteinFeatureKind::Cytoplasmic,
            ProteinFeatureKind::Extracellular,
        ];
        let mut report = TopologyTrainingReport::default();
        let mut genome = TwoBitReader::open(&self.genome_twobit)?;
        let mut examples = Vec::new();

        for protein in &self.proteins {
            let selected: Vec<&ProteinFeature> = protein
                .features
                .iter()
                .filter(|f| kinds.contains(&f.kind))
                .collect();
            if selected.is_empty() {
                continue;
            }
            report.candidate_proteins += 1;
            if selected.iter().any(|f| f.protein_range.is_none()) {
                report.rejected_missing_range += 1;
                continue;
            }
            let sequence = match self.protein_sequence_with_reader(protein, &mut genome) {
                Ok(sequence) if !sequence.is_empty() => sequence,
                _ => {
                    report.rejected_sequence += 1;
                    continue;
                }
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
                    ProteinFeatureKind::Cytoplasmic | ProteinFeatureKind::Extracellular => {
                        TopologyState::region_state(
                            feature.kind,
                            start as usize,
                            end as usize,
                            sequence.len(),
                        )
                    }
                    _ => unreachable!(),
                };
                // Membrane crossings outrank broad sidedness annotations.
                let p = if state == TopologyState::Transmembrane {
                    2
                } else {
                    1
                };
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
            if conflict {
                continue;
            }
            report.reconstructed_proteins += 1;
            examples.push(TopologyTrainingExample {
                accession: protein.accession.clone(),
                sequence,
                truth,
            });
        }
        examples.sort_by(|a, b| a.accession.cmp(&b.accession));
        let mut train = Vec::with_capacity((examples.len() + 1) / 2);
        let mut test = Vec::with_capacity(examples.len() / 2);
        for (i, example) in examples.into_iter().enumerate() {
            if i % 2 == 0 {
                train.push(example);
            } else {
                test.push(example);
            }
        }
        Ok(TopologyTrainingCorpus {
            train,
            test,
            report,
        })
    }

    pub fn train_exact_aa_topology_model(
        &self,
        corpus: &TopologyTrainingCorpus,
    ) -> Result<ExactAaTopologyModel> {
        let n = TopologyState::COUNT;
        let mut initial = vec![1.0f64; n];
        // Small pseudocount everywhere: observed biology dominates, but a transition
        // absent from mouse training is not made literally impossible in another species.
        let mut transition = vec![0.1f64; n * n];
        let mut emission = vec![vec![1.0f64; 32]; n];
        for example in &corpus.train {
            if example.sequence.is_empty() {
                continue;
            }
            initial[example.truth[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = example.truth[pos];
                let code = example
                    .sequence
                    .get(pos)
                    .context("protein sequence position disappeared")?
                    .code() as usize;
                emission[state][code] += 1.0;
                if pos > 0 {
                    transition[example.truth[pos - 1] * n + state] += 1.0;
                }
            }
        }
        normalize_slice(&mut initial);
        for row in 0..n {
            normalize_slice(&mut transition[row * n..(row + 1) * n]);
        }
        for row in &mut emission {
            normalize_slice(row);
        }
        Ok(ExactAaTopologyModel {
            initial,
            transition,
            emission,
        })
    }

    pub fn train_chemistry_topology_model(
        &self,
        corpus: &TopologyTrainingCorpus,
    ) -> Result<ChemistryTopologyModel> {
        let n = TopologyState::COUNT;
        let mut categories = Vec::<u16>::new();
        for example in &corpus.train {
            for pos in 0..example.sequence.len() {
                categories.push(
                    example
                        .sequence
                        .get(pos)
                        .context("protein sequence position disappeared")?
                        .chemistry(),
                );
            }
        }
        categories.sort_unstable();
        categories.dedup();
        if categories.is_empty() {
            bail!("topology chemistry model has no training observations");
        }
        let unknown = categories.len();
        let mut initial = vec![1.0f64; n];
        let mut transition = vec![0.1f64; n * n];
        let mut emission = vec![vec![1.0f64; categories.len() + 1]; n];
        for example in &corpus.train {
            if example.sequence.is_empty() {
                continue;
            }
            initial[example.truth[0]] += 1.0;
            for pos in 0..example.sequence.len() {
                let state = example.truth[pos];
                let bits = example
                    .sequence
                    .get(pos)
                    .context("protein sequence position disappeared")?
                    .chemistry();
                let obs = categories.binary_search(&bits).unwrap_or(unknown);
                emission[state][obs] += 1.0;
                if pos > 0 {
                    transition[example.truth[pos - 1] * n + state] += 1.0;
                }
            }
        }
        normalize_slice(&mut initial);
        for row in 0..n {
            normalize_slice(&mut transition[row * n..(row + 1) * n]);
        }
        for row in &mut emission {
            normalize_slice(row);
        }
        Ok(ChemistryTopologyModel {
            categories,
            initial,
            transition,
            emission,
        })
    }

    pub fn evaluate_exact_aa_topology_model(
        &self,
        model: &ExactAaTopologyModel,
        examples: &[TopologyTrainingExample],
    ) -> Result<TopologyEvaluation> {
        let hmm = model.hmm()?;
        evaluate_topology(examples, |example| {
            let observations = aa_observations(&example.sequence)?;
            Ok(hmm
                .infer(&observations)?
                .viterbi()
                .iter()
                .map(|s| s.0)
                .collect())
        })
    }

    pub fn evaluate_chemistry_topology_model(
        &self,
        model: &ChemistryTopologyModel,
        examples: &[TopologyTrainingExample],
    ) -> Result<TopologyEvaluation> {
        let hmm = model.hmm()?;
        evaluate_topology(examples, |example| {
            let observations = chemistry_observations(&example.sequence, &model.categories)?;
            Ok(hmm
                .infer(&observations)?
                .viterbi()
                .iter()
                .map(|s| s.0)
                .collect())
        })
    }

    /// Train every currently supported observation model for every feature class
    /// with a non-empty deterministic train/test split. One bad/empty feature does
    /// not abort the vault; it is recorded in the report instead.
    pub fn train_aa_model_vault(
        &self,
        flank_width: usize,
    ) -> Result<(AaModelVault, ModelVaultTrainingReport)> {
        let mut models = Vec::new();
        let mut report = ModelVaultTrainingReport::default();
        for feature in ProteinFeatureKind::ALL {
            // Chain is effectively whole-protein coverage and Conflict is a
            // curation discrepancy, not a sequence feature. Neither belongs in
            // this biological feature-modelling experiment.
            if matches!(
                feature,
                ProteinFeatureKind::Chain
                    | ProteinFeatureKind::Conflict
                    | ProteinFeatureKind::SignalPeptide
            ) {
                continue;
            }
            report.feature_classes_considered += 1;
            let corpus = match self.training_corpus_with_flank_diagnostic(feature, flank_width) {
                Ok(corpus) => corpus,
                Err(err) => {
                    report.skipped.push((feature, err.to_string()));
                    continue;
                }
            };
            if corpus.train.is_empty()
                || corpus.test.is_empty()
                || corpus.report.feature_residues == 0
            {
                report.skipped.push((
                    feature,
                    format!(
                        "insufficient corpus: train {}, test {}, feature residues {}",
                        corpus.train.len(),
                        corpus.test.len(),
                        corpus.report.feature_residues
                    ),
                ));
                continue;
            }
            let exact = self.train_exact_aa_feature_model(&corpus, flank_width)?;
            let exact_eval = self.evaluate_exact_aa_feature_model(&exact, &corpus.test)?;
            models.push(ProteinFeatureModel::ExactAa {
                model: exact,
                evaluation: exact_eval,
            });
            let chemistry = self.train_chemistry_feature_model(&corpus, flank_width)?;
            let chemistry_eval = self.evaluate_chemistry_feature_model(&chemistry, &corpus.test)?;
            models.push(ProteinFeatureModel::Chemistry {
                model: chemistry,
                evaluation: chemistry_eval,
            });
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
            models.push(ProteinFeatureModel::TopologyExactAa {
                model: exact,
                evaluation: exact_eval,
            });
            let chemistry = self.train_chemistry_topology_model(&topology)?;
            let chemistry_eval =
                self.evaluate_chemistry_topology_model(&chemistry, &topology.test)?;
            models.push(ProteinFeatureModel::TopologyChemistry {
                model: chemistry,
                evaluation: chemistry_eval,
            });
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

}
