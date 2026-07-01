use anyhow::{bail, Context, Result};
use int_to_str::IntToStr;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TagEntry {
    pub id: usize,
    pub name: String,
    pub sequence: Vec<u8>,
}

impl TagEntry {
    pub fn new(id: usize, name: impl Into<String>, sequence: impl AsRef<[u8]>) -> Self {
        Self {
            id,
            name: name.into(),
            sequence: sequence.as_ref().to_ascii_uppercase(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagCall {
    Match {
        tag_id: usize,
        tag_name: String,
        matches: usize,
        scanned_kmers: usize,
    },

    NoMatch {
        best_matches: usize,
        scanned_kmers: usize,
    },

    MultiMatch {
        tag_ids: Vec<usize>,
        tag_names: Vec<String>,
        matches: usize,
        scanned_kmers: usize,
    },
}

impl TagCall {
    pub fn is_match(&self) -> bool {
        matches!(self, Self::Match { .. })
    }

    pub fn status(&self) -> &'static str {
        match self {
            Self::Match { .. } => "match",
            Self::NoMatch { .. } => "no_match",
            Self::MultiMatch { .. } => "multi_match",
        }
    }

    pub fn tag_id(&self) -> Option<usize> {
        match self {
            Self::Match { tag_id, .. } => Some(*tag_id),
            _ => None,
        }
    }

    pub fn tag_name(&self) -> Option<&str> {
        match self {
            Self::Match { tag_name, .. } => Some(tag_name.as_str()),
            _ => None,
        }
    }

    pub fn matches(&self) -> usize {
        match self {
            Self::Match { matches, .. } => *matches,
            Self::NoMatch { best_matches, .. } => *best_matches,
            Self::MultiMatch { matches, .. } => *matches,
        }
    }
}

/// Generic fast k-mer tag mapper.
///
/// This is the reusable extraction of the old Rustody sample-id idea:
///
/// - build a `kmer -> tag_id` map from known tag sequences
/// - remove kmers that occur in more than one tag
/// - scan a read and count votes per tag
/// - return the unique best tag, no match, or multimatch
///
/// K-mer encoding uses `int_to_str::IntToStr`.
#[derive(Debug, Clone)]
pub struct FastTagMapper {
    name: String,
    kmer_size: usize,
    tags: Vec<TagEntry>,
    kmers: HashMap<u64, usize>,
    bad_kmers: HashSet<u64>,
    min_matches: usize,
}

impl FastTagMapper {
    pub fn new(name: impl Into<String>, kmer_size: usize) -> Self {
        Self {
            name: name.into(),
            kmer_size,
            tags: Vec::new(),
            kmers: HashMap::new(),
            bad_kmers: HashSet::new(),
            min_matches: 3,
        }
    }

    pub fn with_min_matches(mut self, min_matches: usize) -> Self {
        self.min_matches = min_matches;
        self
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }

    pub fn min_matches(&self) -> usize {
        self.min_matches
    }

    pub fn tags(&self) -> &[TagEntry] {
        &self.tags
    }

    pub fn len(&self) -> usize {
        self.tags.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tags.is_empty()
    }

    pub fn add_tag(&mut self, name: impl Into<String>, sequence: impl AsRef<[u8]>) -> Result<usize> {
        let id = self.tags.len();
        let entry = TagEntry::new(id, name, sequence);

        if self.kmer_size == 0 {
            bail!("kmer_size must be > 0");
        }

        if entry.sequence.len() < self.kmer_size {
            bail!(
                "tag '{}' is shorter than kmer_size {}",
                entry.name,
                self.kmer_size
            );
        }

        self.add_tag_kmers(id, &entry.sequence);
        self.tags.push(entry);

        Ok(id)
    }

    pub fn from_entries(
        name: impl Into<String>,
        kmer_size: usize,
        entries: impl IntoIterator<Item = (String, Vec<u8>)>,
    ) -> Result<Self> {
        let mut mapper = Self::new(name, kmer_size);

        for (tag_name, seq) in entries {
            mapper.add_tag(tag_name, seq)?;
        }

        Ok(mapper)
    }

    /// Read a simple two-column TSV:
    ///
    /// name<TAB>sequence
    ///
    /// A header is allowed if the first row contains `name` and `sequence`.
    pub fn from_tsv(path: &Path, name: impl Into<String>, kmer_size: usize) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut mapper = Self::new(name, kmer_size);
        let mut first_data_line = true;

        for line in reader.lines() {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 2 {
                bail!("tag TSV line needs at least two columns: {line}");
            }

            if first_data_line {
                first_data_line = false;

                let c0 = cols[0].trim().to_ascii_lowercase();
                let c1 = cols[1].trim().to_ascii_lowercase();

                if c0 == "name" && (c1 == "sequence" || c1 == "seq") {
                    continue;
                }
            }

            mapper.add_tag(cols[0].trim(), cols[1].trim().as_bytes())?;
        }

        Ok(mapper)
    }

    /// Read FASTA records as tag entries.
    ///
    /// Header up to first whitespace becomes the tag name.
    pub fn from_fasta(path: &Path, name: impl Into<String>, kmer_size: usize) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut mapper = Self::new(name, kmer_size);
        let mut current_name: Option<String> = None;
        let mut current_seq = Vec::<u8>::new();

        for line in reader.lines() {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            let line = line.trim();

            if line.is_empty() {
                continue;
            }

            if let Some(rest) = line.strip_prefix('>') {
                if let Some(tag_name) = current_name.take() {
                    mapper.add_tag(tag_name, &current_seq)?;
                    current_seq.clear();
                }

                let tag_name = rest
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();

                if tag_name.is_empty() {
                    bail!("empty FASTA tag name in {}", path.display());
                }

                current_name = Some(tag_name);
            } else {
                current_seq.extend_from_slice(line.as_bytes());
            }
        }

        if let Some(tag_name) = current_name.take() {
            mapper.add_tag(tag_name, &current_seq)?;
        }

        Ok(mapper)
    }

    pub fn bd_human(kmer_size: usize) -> Result<Self> {
        let mut mapper = Self::new("bd-human", kmer_size);

        for (i, seq) in BD_HUMAN_SAMPLE_TAGS.iter().enumerate() {
            mapper.add_tag(format!("BD_Human_SampleTag{:02}", i + 1), *seq)?;
        }

        Ok(mapper)
    }

    pub fn bd_mouse(kmer_size: usize) -> Result<Self> {
        let mut mapper = Self::new("bd-mouse", kmer_size);

        for (i, seq) in BD_MOUSE_SAMPLE_TAGS.iter().enumerate() {
            mapper.add_tag(format!("BD_Mouse_SampleTag{:02}", i + 1), *seq)?;
        }

        Ok(mapper)
    }

    /// Match a read using every kmer after `start`, stepping by `jump`.
    ///
    /// This mirrors old Rustody `SampleIds::get(seq, jump, start)`,
    /// but returns a structured `TagCall`.
    pub fn call(&self, seq: &[u8], jump: usize, start: usize) -> TagCall {
        if self.tags.is_empty() || self.kmer_size == 0 || seq.len() < self.kmer_size {
            return TagCall::NoMatch {
                best_matches: 0,
                scanned_kmers: 0,
            };
        }

        let jump = jump.max(1);
        let mut counts = vec![0usize; self.tags.len()];
        let mut scanned_kmers = 0usize;

        let mut pos = start;
        while pos + self.kmer_size <= seq.len() {
            let kmer = &seq[pos..pos + self.kmer_size];

            if let Some(encoded) = encode_kmer(kmer) {
                scanned_kmers += 1;

                if let Some(&tag_id) = self.kmers.get(&encoded) {
                    counts[tag_id] += 1;
                }
            }

            pos += jump;
        }

        let best_matches = counts.iter().copied().max().unwrap_or(0);

        if best_matches < self.min_matches {
            return TagCall::NoMatch {
                best_matches,
                scanned_kmers,
            };
        }

        let best_ids: Vec<usize> = counts
            .iter()
            .enumerate()
            .filter_map(|(tag_id, &count)| (count == best_matches).then_some(tag_id))
            .collect();

        if best_ids.len() == 1 {
            let tag = &self.tags[best_ids[0]];

            TagCall::Match {
                tag_id: tag.id,
                tag_name: tag.name.clone(),
                matches: best_matches,
                scanned_kmers,
            }
        } else {
            let tag_names = best_ids
                .iter()
                .map(|&id| self.tags[id].name.clone())
                .collect();

            TagCall::MultiMatch {
                tag_ids: best_ids,
                tag_names,
                matches: best_matches,
                scanned_kmers,
            }
        }
    }

    fn add_tag_kmers(&mut self, tag_id: usize, seq: &[u8]) {
        let mut seen_in_this_tag = HashSet::<u64>::new();

        for pos in 0..=seq.len() - self.kmer_size {
            let kmer = &seq[pos..pos + self.kmer_size];
            let Some(encoded) = encode_kmer(kmer) else {
                continue;
            };

            if !seen_in_this_tag.insert(encoded) {
                continue;
            }

            if self.bad_kmers.contains(&encoded) {
                continue;
            }

            match self.kmers.get(&encoded).copied() {
                None => {
                    self.kmers.insert(encoded, tag_id);
                }
                Some(existing_tag_id) if existing_tag_id == tag_id => {}
                Some(_) => {
                    self.kmers.remove(&encoded);
                    self.bad_kmers.insert(encoded);
                }
            }
        }
    }
}

impl fmt::Display for FastTagMapper {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "FastTagMapper: {}", self.name)?;
        writeln!(f, "  tags       : {}", self.tags.len())?;
        writeln!(f, "  kmer_size  : {}", self.kmer_size)?;
        writeln!(f, "  kmers      : {}", self.kmers.len())?;
        writeln!(f, "  bad_kmers  : {}", self.bad_kmers.len())?;
        writeln!(f, "  min_matches: {}", self.min_matches)
    }
}

fn encode_kmer(kmer: &[u8]) -> Option<u64> {
    if kmer.iter().any(|b| {
        let b = b.to_ascii_uppercase();
        b == b'N' || !(b == b'A' || b == b'C' || b == b'G' || b == b'T')
    }) {
        return None;
    }

    Some(IntToStr::new(kmer.to_ascii_uppercase()).into_u64())
}

pub const BD_HUMAN_SAMPLE_TAGS: [&[u8]; 12] = [
    b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG",
    b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG",
    b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT",
    b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT",
    b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG",
    b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC",
    b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT",
    b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT",
    b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG",
    b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC",
    b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC",
    b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG",
];

pub const BD_MOUSE_SAMPLE_TAGS: [&[u8]; 12] = [
    b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG",
    b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC",
    b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC",
    b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
    b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT",
    b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
    b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC",
    b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
    b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC",
    b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
    b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT",
    b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC",
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bd_human_sample_tag_exact_match() {
        let mapper = FastTagMapper::bd_human(9).unwrap();
        let call = mapper.call(BD_HUMAN_SAMPLE_TAGS[0], 1, 0);

        match call {
            TagCall::Match { tag_id, tag_name, matches, .. } => {
                assert_eq!(tag_id, 0);
                assert_eq!(tag_name, "BD_Human_SampleTag01");
                assert!(matches >= 3);
            }
            other => panic!("unexpected call: {other:?}"),
        }
    }

    #[test]
    fn no_match_for_junk() {
        let mapper = FastTagMapper::bd_mouse(9).unwrap();
        let call = mapper.call(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 1, 0);

        assert!(matches!(call, TagCall::NoMatch { .. }));
    }

    #[test]
    fn custom_mapper() {
        let mut mapper = FastTagMapper::new("custom", 4).with_min_matches(2);
        mapper.add_tag("tag_a", b"AAAACCCC").unwrap();
        mapper.add_tag("tag_b", b"GGGGTTTT").unwrap();

        let call = mapper.call(b"NNNAAAACCCCNNN", 1, 0);

        match call {
            TagCall::Match { tag_name, .. } => assert_eq!(tag_name, "tag_a"),
            other => panic!("unexpected call: {other:?}"),
        }
    }
}

