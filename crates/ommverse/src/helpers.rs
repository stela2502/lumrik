fn add_identifier(identifiers: &mut Vec<String>, value: &str) {
    let value = value.trim();
    if value.is_empty() { return; }
    if !identifiers.iter().any(|known| known == value) { identifiers.push(value.to_owned()); }
    let base = strip_version(value);
    if base != value && !identifiers.iter().any(|known| known == base) { identifiers.push(base.to_owned()); }
}

fn read_all_bigbed(
    path: &Path,
    mut visit: impl FnMut(&str, u32, u32, &str) -> Result<()>,
) -> Result<()> {
    let mut bb = BigBedRead::open_file(path)
        .with_context(|| format!("opening bigBed {}", path.display()))?;
    let chroms: Vec<(String, u32)> = bb
        .chroms()
        .iter()
        .map(|c| (c.name.clone(), c.length))
        .collect();
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
enum MappingSchema {
    EnsGene,
    UnipAliSwissprot,
}

fn find_swissprot_mapping_bigbed(dir: &Path) -> Result<(PathBuf, MappingSchema)> {
    // unipAliSwissprot.bb is the canonical Swiss-Prot alignment track exposed
    // by UCSC's UniProt trackDb.  Hashed ensGene_*.swissprot.bb files are
    // auxiliary mapping products and multiple variants can coexist in one
    // UniProt release, so they must not override the canonical track.
    let canonical = dir.join("unipAliSwissprot.bb");
    if canonical.is_file() {
        return Ok((canonical, MappingSchema::UnipAliSwissprot));
    }

    // Keep compatibility with older source trees that contain only one of the
    // historical ensGene mapping products.  Never guess when several exist.
    let mut hits = Vec::new();
    for e in fs::read_dir(dir)? {
        let p = e?.path();
        let name = p.file_name().and_then(|x| x.to_str()).unwrap_or("");
        if name.starts_with("ensGene_") && name.ends_with(".swissprot.bb") {
            hits.push(p);
        }
    }
    hits.sort();
    match hits.len() {
        1 => Ok((hits.remove(0), MappingSchema::EnsGene)),
        n if n > 1 => bail!(
            "multiple ensGene_*.swissprot.bb mapping bigBeds found in {} and canonical unipAliSwissprot.bb is absent; refusing to guess: {:?}",
            dir.display(),
            hits
        ),
        _ => bail!(
            "no Swiss-Prot mapping bigBed found in {} (tried unipAliSwissprot.bb and ensGene_*.swissprot.bb)",
            dir.display()
        ),
    }
}

fn strip_version(s: &str) -> &str {
    s.rsplit_once('.')
        .filter(|(_, v)| v.chars().all(|c| c.is_ascii_digit()))
        .map(|(a, _)| a)
        .unwrap_or(s)
}
fn nonempty(s: &str) -> Option<String> {
    let s = s.trim();
    (!s.is_empty()).then(|| s.to_owned())
}
fn review_status(s: &str) -> ReviewStatus {
    if s.contains("Swiss-Prot") || s.eq_ignore_ascii_case("swissprot") {
        ReviewStatus::SwissProt
    } else if s.contains("TrEMBL") || s.eq_ignore_ascii_case("trembl") {
        ReviewStatus::Trembl
    } else {
        ReviewStatus::Other
    }
}
fn parse_aliases(a: &str, b: &str) -> Vec<String> {
    let mut v = Vec::new();
    for x in a
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .chain(b.split(|c: char| c == ',' || c == ';' || c.is_whitespace()))
    {
        let x = x.trim();
        if !x.is_empty() && !v.iter().any(|y| y == x) {
            v.push(x.to_owned());
        }
    }
    v
}
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

fn protein_accession_from_feature_text(text: &str) -> Option<&str> {
    text.rsplit_once(" on protein ")
        .map(|(_, accession)| accession.trim())
        .filter(|accession| !accession.is_empty())
}

fn parse_amino_acid_range(text: &str) -> Option<(u32, u32)> {
    let tail = text
        .strip_prefix("amino acids ")
        .or_else(|| text.strip_prefix("amino acid "))?;
    let token = tail.split_whitespace().next()?;
    let (a, b) = token
        .split_once('-')
        .map(|(a, b)| (a, b))
        .unwrap_or((token, token));
    let start: u32 = a.parse().ok()?;
    let end: u32 = b.parse().ok()?;
    (start > 0 && end >= start).then_some((start - 1, end))
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
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != AA_MODEL_VAULT_MAGIC {
            bail!("not an Ommverse AA model vault");
        }
        let mut version = [0u8; 4];
        file.read_exact(&mut version)?;
        let version = u32::from_le_bytes(version);
        if version != AA_MODEL_VAULT_FORMAT_VERSION {
            bail!("unsupported AA model vault format version {version}");
        }
        Ok(bincode::deserialize_from(file)?)
    }
}

impl ExactAaTopologyModel {
    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self
            .emission
            .iter()
            .map(|row| CategoricalEmission::new(row.clone()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(
            self.initial.clone(),
            self.transition.clone(),
            emissions,
        )?)
    }
}

impl ChemistryTopologyModel {
    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self
            .emission
            .iter()
            .map(|row| CategoricalEmission::new(row.clone()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(
            self.initial.clone(),
            self.transition.clone(),
            emissions,
        )?)
    }
}

impl ChemistryFeatureModel {
    fn hmm(&self) -> Result<Hmm<CategoricalEmission>> {
        let emissions = self
            .emission
            .iter()
            .map(|row| CategoricalEmission::new(row.clone()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(
            self.initial.to_vec(),
            self.transition.to_vec(),
            emissions,
        )?)
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
        let emissions = self
            .emission
            .iter()
            .map(|row| CategoricalEmission::new(row.to_vec()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Hmm::new(
            self.initial.to_vec(),
            self.transition.to_vec(),
            emissions,
        )?)
    }
}

fn chemistry_observations(sequence: &IntToProt, categories: &[u16]) -> Result<Vec<usize>> {
    let unknown = categories.len();
    (0..sequence.len())
        .map(|pos| {
            let bits = sequence
                .get(pos)
                .context("protein sequence position disappeared")?
                .chemistry();
            Ok(categories.binary_search(&bits).unwrap_or(unknown))
        })
        .collect()
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
    let partials: Result<Vec<TopologyEvaluation>> = examples
        .par_iter()
        .map(|example| {
            let predicted = predict(example)?;
            if predicted.len() != example.truth.len() {
                bail!("topology prediction length mismatch");
            }
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
                let truth_mask: Vec<bool> = example
                    .truth
                    .iter()
                    .map(|&x| {
                        topology_state_from_index(x).and_then(TopologyState::feature)
                            == Some(*feature)
                    })
                    .collect();
                let pred_mask: Vec<bool> = predicted
                    .iter()
                    .map(|&x| {
                        topology_state_from_index(x).and_then(TopologyState::feature)
                            == Some(*feature)
                    })
                    .collect();
                let r = &mut reports[i];
                r.proteins += 1;
                r.residues += truth_mask.len();
                for (&truth, &pred) in truth_mask.iter().zip(&pred_mask) {
                    match (truth, pred) {
                        (true, true) => r.true_positive += 1,
                        (false, true) => r.false_positive += 1,
                        (false, false) => r.true_negative += 1,
                        (true, false) => r.false_negative += 1,
                    }
                }
                let truth_segments = bool_segments(&truth_mask);
                let pred_segments = bool_segments(&pred_mask);
                r.truth_segments += truth_segments.len();
                r.predicted_segments += pred_segments.len();
                r.recovered_segments += truth_segments
                    .iter()
                    .filter(|truth| {
                        pred_segments
                            .iter()
                            .any(|pred| overlap_fraction(**truth, *pred) >= 0.5)
                    })
                    .count();
            }
            Ok(TopologyEvaluation {
                states: reports,
                latent_truth_residues: truth_residues,
                latent_predicted_residues: predicted_residues,
                latent_confusion: confusion,
            })
        })
        .collect();
    Ok(partials?
        .into_iter()
        .fold(empty_topology_evaluation(), merge_topology_evaluation))
}

fn empty_topology_evaluation() -> TopologyEvaluation {
    TopologyEvaluation {
        states: vec![FeatureEvaluation::default(); TopologyState::BIOLOGICAL.len()],
        latent_truth_residues: vec![0; TopologyState::COUNT],
        latent_predicted_residues: vec![0; TopologyState::COUNT],
        latent_confusion: vec![0; TopologyState::COUNT * TopologyState::COUNT],
    }
}

fn merge_topology_evaluation(
    mut a: TopologyEvaluation,
    b: TopologyEvaluation,
) -> TopologyEvaluation {
    for (dst, src) in a.states.iter_mut().zip(b.states) {
        *dst = merge_feature_evaluation(dst.clone(), src);
    }
    for (dst, src) in a
        .latent_truth_residues
        .iter_mut()
        .zip(b.latent_truth_residues)
    {
        *dst += src;
    }
    for (dst, src) in a
        .latent_predicted_residues
        .iter_mut()
        .zip(b.latent_predicted_residues)
    {
        *dst += src;
    }
    for (dst, src) in a.latent_confusion.iter_mut().zip(b.latent_confusion) {
        *dst += src;
    }
    a
}

fn merge_feature_evaluation(mut a: FeatureEvaluation, b: FeatureEvaluation) -> FeatureEvaluation {
    a.proteins += b.proteins;
    a.residues += b.residues;
    a.true_positive += b.true_positive;
    a.false_positive += b.false_positive;
    a.true_negative += b.true_negative;
    a.false_negative += b.false_negative;
    a.truth_segments += b.truth_segments;
    a.predicted_segments += b.predicted_segments;
    a.recovered_segments += b.recovered_segments;
    a
}

fn accumulate_evaluation(
    report: &mut FeatureEvaluation,
    predicted: &[bool],
    example: &ProteinTrainingExample,
) {
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
    report.recovered_segments += truth_segments
        .iter()
        .filter(|truth| {
            pred_segments
                .iter()
                .any(|pred| overlap_fraction(**truth, *pred) >= 0.5)
        })
        .count();
}

fn normalize_slice(values: &mut [f64]) {
    let sum: f64 = values.iter().sum();
    if sum > 0.0 {
        for value in values {
            *value /= sum;
        }
    }
}

fn distribution_summary(values: &[usize]) -> DistributionSummary {
    if values.is_empty() {
        return DistributionSummary::default();
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let count = sorted.len();
    let mean = sorted.iter().map(|&x| x as f64).sum::<f64>() / count as f64;
    let variance = sorted
        .iter()
        .map(|&x| {
            let d = x as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / count as f64;
    let percentile = |q: f64| -> f64 {
        if count == 1 {
            return sorted[0] as f64;
        }
        let pos = q * (count - 1) as f64;
        let lo = pos.floor() as usize;
        let hi = pos.ceil() as usize;
        let frac = pos - lo as f64;
        sorted[lo] as f64 * (1.0 - frac) + sorted[hi] as f64 * frac
    };
    DistributionSummary {
        count,
        mean,
        median: percentile(0.5),
        sd: variance.sqrt(),
        q1: percentile(0.25),
        q3: percentile(0.75),
        min: sorted[0],
        max: sorted[count - 1],
    }
}

fn aa_observations(sequence: &IntToProt) -> Result<Vec<usize>> {
    (0..sequence.len())
        .map(|pos| {
            sequence
                .get(pos)
                .map(|aa| aa.code() as usize)
                .context("protein sequence position disappeared")
        })
        .collect()
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
    for &(start, end) in &segments {
        states[start..end].fill(2);
    }
    for pos in 0..truth.len() {
        if states[pos] == 2 {
            continue;
        }
        let prev = segments
            .iter()
            .filter(|(_, end)| *end <= pos)
            .map(|(_, end)| pos + 1 - *end)
            .min();
        let next = segments
            .iter()
            .filter(|(start, _)| *start > pos)
            .map(|(start, _)| *start - pos)
            .min();
        let prev = prev.filter(|&d| d <= flank_width);
        let next = next.filter(|&d| d <= flank_width);
        states[pos] = match (prev, next) {
            (Some(a), Some(b)) => {
                if b <= a {
                    1
                } else {
                    3
                }
            }
            (None, Some(_)) => 1,
            (Some(_), None) => 3,
            (None, None) => 0,
        };
    }
    states
}

fn bool_segments(values: &[bool]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut start = None;
    for (i, &value) in values.iter().enumerate() {
        match (start, value) {
            (None, true) => start = Some(i),
            (Some(s), false) => {
                out.push((s, i));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(s) = start {
        out.push((s, values.len()));
    }
    out
}

fn overlap_fraction(a: (usize, usize), b: (usize, usize)) -> f64 {
    let overlap = a.1.min(b.1).saturating_sub(a.0.max(b.0));
    if overlap == 0 {
        0.0
    } else {
        overlap as f64 / (a.1 - a.0) as f64
    }
}

fn ratio(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}
fn normalize_array<const N: usize>(values: &mut [f64; N]) {
    let sum: f64 = values.iter().sum();
    for x in values {
        *x /= sum;
    }
}

