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
        mismatches: usize,
        matched_bases: usize,
        expected_bases: usize,
        offset: usize,
        scanned_offsets: usize,
    },
    NoMatch {
        best_mismatches: usize,
        best_matched_bases: usize,
        expected_bases: usize,
        scanned_offsets: usize,
    },
    MultiMatch {
        tag_ids: Vec<usize>,
        tag_names: Vec<String>,
        mismatches: usize,
        matched_bases: usize,
        expected_bases: usize,
        scanned_offsets: usize,
    },
}

impl TagCall {
    pub fn status(&self) -> &'static str {
        match self {
            Self::Match { .. } => "match",
            Self::NoMatch { .. } => "no_match",
            Self::MultiMatch { .. } => "multi_match",
        }
    }

    pub fn success(&self) -> bool{
        match self {
          Self::Match{..} => true,
          _ => false,
        }
    }

    pub fn tag_id_string(&self) -> String {
        match self {
            Self::Match { tag_id, .. } => tag_id.to_string(),
            Self::MultiMatch { tag_ids, .. } => tag_ids.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            Self::NoMatch { .. } => String::new(),
        }
    }

    pub fn tag_name_string(&self) -> String {
        match self {
            Self::Match { tag_name, .. } => tag_name.clone(),
            Self::MultiMatch { tag_names, .. } => tag_names.join(","),
            Self::NoMatch { .. } => String::new(),
        }
    }

    pub fn best_tag_id(&self) -> Option<u64> {
        match self {
            Self::Match { tag_id, .. } => Some(*tag_id as u64),
            Self::NoMatch { .. } => None,
            Self::MultiMatch { .. } => None,
        }
    }
    pub fn diagnostic(&self) -> String {
        match self {
            Self::Match {
                tag_name,
                mismatches,
                matched_bases,
                expected_bases,
                offset,
                scanned_offsets,
                ..
            } => format!(
                "match: {tag_name}; offset={offset}; mismatches={mismatches}; overlap={matched_bases}/{expected_bases}; scanned_offsets={scanned_offsets}"
            ),

            Self::NoMatch {
                best_mismatches,
                best_matched_bases,
                expected_bases,
                scanned_offsets,
            } => format!(
                "no_match: best_mismatches={best_mismatches}; best_overlap={best_matched_bases}/{expected_bases}; scanned_offsets={scanned_offsets}"
            ),

            Self::MultiMatch {
                tag_names,
                mismatches,
                matched_bases,
                expected_bases,
                scanned_offsets,
                ..
            } => format!(
                "multi_match: {}; mismatches={mismatches}; overlap={matched_bases}/{expected_bases}; scanned_offsets={scanned_offsets}",
                tag_names.join(",")
            ),
        }
    }
}

/// Very strict tag mapper for BD sample tags / HTO-like reads.
///
/// This does NOT do whole-read permissive kmer-voting.
/// It first finds candidate tags with a small kmer seed, then validates by
/// base-level comparison over the expected tag-length window.
///
/// Acceptance:
/// - the tag must start close to the beginning of the read (`max_leading_bases`)
/// - the compared tag/read overlap must cover at least `min_overlap_fraction`
/// - at most `max_mismatches` over the overlap
/// - trailing junk/Ns after the tag are allowed
/// - large leading junk before the tag is rejected
#[derive(Debug, Clone)]
pub struct FastTagMapper {
    name: String,
    kmer_size: usize,
    tags: Vec<TagEntry>,
    seed_to_tags: HashMap<u64, Vec<usize>>,
    name_to_id: HashMap<String, u64>,
    max_mismatches: usize,
    min_overlap_fraction: f32,
    max_leading_bases: usize,
}

impl FastTagMapper {
    pub fn new(name: impl Into<String>, kmer_size: usize) -> Self {
        Self {
            name: name.into(),
            kmer_size,
            
            tags: Vec::new(),
            name_to_id: HashMap::new(),

            seed_to_tags: HashMap::new(),
            max_mismatches: 2,
            min_overlap_fraction: 0.75,
            max_leading_bases: 4,
        }
    }
    pub fn tags(&self) -> &[TagEntry] {
        &self.tags
    }

    pub fn tag(&self, id: usize) -> Option<&TagEntry> {
        self.tags.get(id)
    }

    pub fn tag_count(&self) -> usize {
        self.tags.len()
    }

    pub fn with_max_mismatches(mut self, max_mismatches: usize) -> Self {
        self.max_mismatches = max_mismatches;
        self
    }

    pub fn with_min_overlap_fraction(mut self, min_overlap_fraction: f32) -> Self {
        self.min_overlap_fraction = min_overlap_fraction;
        self
    }

    pub fn with_max_leading_bases(mut self, max_leading_bases: usize) -> Self {
        self.max_leading_bases = max_leading_bases;
        self
    }

    pub fn name(&self) -> &str {
        &self.name
    }


    pub fn add_tag(&mut self, name: impl Into<String>, sequence: impl AsRef<[u8]>) -> Result<usize> {
        if self.kmer_size == 0 {
            bail!("kmer_size must be > 0");
        }

        let id = self.tags.len();
        let entry = TagEntry::new(id, name, sequence);

        if self
            .name_to_id
            .insert(entry.name.clone(), id as u64)
            .is_some()
        {
            bail!("duplicate tag name '{}'", entry.name);
        }

        if entry.sequence.len() < self.kmer_size {
            bail!("tag '{}' is shorter than kmer_size {}", entry.name, self.kmer_size);
        }

        let mut seen = HashSet::new();
        for pos in 0..=entry.sequence.len() - self.kmer_size {
            if let Some(seed) = encode_kmer(&entry.sequence[pos..pos + self.kmer_size]) {
                if seen.insert(seed) {
                    self.seed_to_tags.entry(seed).or_default().push(id);
                }
            }
        }

        self.tags.push(entry);
        Ok(id)
    }

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

                let tag_name = rest.split_whitespace().next().unwrap_or("").to_string();
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

    /// Strict call for sample/tag reads.
    ///
    /// This tries candidate offsets from 0..=max_leading_bases.
    /// The tag is expected to begin close to the read start. Trailing sequence is ignored.
    pub fn call(&self, read: &[u8]) -> TagCall {
        if self.tags.is_empty() || read.is_empty() {
            return TagCall::NoMatch {
                best_mismatches: usize::MAX,
                best_matched_bases: 0,
                expected_bases: 0,
                scanned_offsets: 0,
            };
        }

        let read = read.to_ascii_uppercase();
        let max_offset = self.max_leading_bases.min(read.len().saturating_sub(1));
        let mut candidates = HashSet::<(usize, usize)>::new(); // (tag_id, offset)
        let mut scanned_offsets = 0usize;

        for offset in 0..=max_offset {
            scanned_offsets += 1;
            if offset + self.kmer_size > read.len() {
                continue;
            }

            let scan_end = (offset + 16).min(read.len().saturating_sub(self.kmer_size) + 1);
            for pos in offset..scan_end {
                if pos + self.kmer_size > read.len() {
                    break;
                }

                let Some(seed) = encode_kmer(&read[pos..pos + self.kmer_size]) else {
                    continue;
                };

                if let Some(tag_ids) = self.seed_to_tags.get(&seed) {
                    for &tag_id in tag_ids {
                        candidates.insert((tag_id, offset));
                    }
                }
            }
        }

        let mut best: Vec<(usize, usize, usize, usize)> = Vec::new();
        // tag_id, offset, mismatches, matched_bases
        let mut best_mismatches = usize::MAX;
        let mut best_matched_bases = 0usize;
        let mut best_expected = 0usize;

        for (tag_id, offset) in candidates {
            let tag = &self.tags[tag_id];
            let available = read.len().saturating_sub(offset);
            let overlap = available.min(tag.sequence.len());

            if overlap == 0 {
                continue;
            }

            let min_required = ((tag.sequence.len() as f32) * self.min_overlap_fraction).ceil() as usize;
            if overlap < min_required {
                if overlap > best_matched_bases {
                    best_matched_bases = overlap;
                    best_expected = tag.sequence.len();
                }
                continue;
            }

            let mut mismatches = 0usize;
            for i in 0..overlap {
                let a = read[offset + i];
                let b = tag.sequence[i];

                if a == b'N' {
                    // N in read = unknown, but not an exact match. Count as mismatch
                    // because sample tags should be strict.
                    mismatches += 1;
                } else if a != b {
                    mismatches += 1;
                }

                if mismatches > self.max_mismatches {
                    break;
                }
            }

            if mismatches < best_mismatches
                || (mismatches == best_mismatches && overlap > best_matched_bases)
            {
                best.clear();
                best.push((tag_id, offset, mismatches, overlap));
                best_mismatches = mismatches;
                best_matched_bases = overlap;
                best_expected = tag.sequence.len();
            } else if mismatches == best_mismatches && overlap == best_matched_bases {
                best.push((tag_id, offset, mismatches, overlap));
            }
        }

        if best.is_empty()
            || best_mismatches > self.max_mismatches
            || best_matched_bases < ((best_expected as f32) * self.min_overlap_fraction).ceil() as usize
        {
            return TagCall::NoMatch {
                best_mismatches,
                best_matched_bases,
                expected_bases: best_expected,
                scanned_offsets,
            };
        }

        let unique_tag_ids: HashSet<usize> = best.iter().map(|x| x.0).collect();
        if unique_tag_ids.len() == 1 {
            let (tag_id, offset, mismatches, matched_bases) = best[0];
            let tag = &self.tags[tag_id];

            TagCall::Match {
                tag_id,
                tag_name: tag.name.clone(),
                mismatches,
                matched_bases,
                expected_bases: tag.sequence.len(),
                offset,
                scanned_offsets,
            }
        } else {
            let mut tag_ids: Vec<usize> = unique_tag_ids.into_iter().collect();
            tag_ids.sort_unstable();
            let tag_names = tag_ids.iter().map(|&id| self.tags[id].name.clone()).collect();

            TagCall::MultiMatch {
                tag_ids,
                tag_names,
                mismatches: best_mismatches,
                matched_bases: best_matched_bases,
                expected_bases: best_expected,
                scanned_offsets,
            }
        }
    }
}

impl fmt::Display for FastTagMapper {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "FastTagMapper: {}", self.name)?;
        writeln!(f, "  tags                : {}", self.tags.len())?;
        writeln!(f, "  kmer_size           : {}", self.kmer_size)?;
        writeln!(f, "  seeds               : {}", self.seed_to_tags.len())?;
        writeln!(f, "  max_mismatches      : {}", self.max_mismatches)?;
        writeln!(f, "  min_overlap_fraction: {:.3}", self.min_overlap_fraction)?;
        writeln!(f, "  max_leading_bases   : {}", self.max_leading_bases)
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

    fn concat(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    fn mutate_many(seq: &[u8], positions: &[usize]) -> Vec<u8> {
        let mut out = seq.to_vec();
        for &pos in positions {
            out[pos] = match out[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                x => x,
            };
        }
        out
    }

    fn expect_match(call: TagCall, id: usize) {
        match call {
            TagCall::Match {
                tag_id,
                mismatches,
                matched_bases,
                expected_bases,
                offset,
                ..
            } => {
                assert_eq!(tag_id, id);
                assert!(mismatches <= 2);
                assert!(matched_bases * 4 >= expected_bases * 3);
                assert!(offset <= 4);
            }
            other => panic!("expected match to {id}, got {other:?}"),
        }
    }

    #[test]
    fn all_bd_mouse_exact_tags_match_with_trailing_ns() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        for (i, tag) in BD_MOUSE_SAMPLE_TAGS.iter().enumerate() {
            let read = concat(&[*tag, b"NNNNNNNNNNNNNNNNNNN"]);
            expect_match(mapper.call(&read), i);
        }
    }

    #[test]
    fn all_bd_human_exact_tags_match_with_trailing_ns() {
        let mapper = FastTagMapper::bd_human(4).unwrap();
        for (i, tag) in BD_HUMAN_SAMPLE_TAGS.iter().enumerate() {
            let read = concat(&[*tag, b"NNNNNNNNNNNNNNNNNNN"]);
            expect_match(mapper.call(&read), i);
        }
    }

    #[test]
    fn leading_junk_before_full_tag_is_rejected() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        let tag = BD_MOUSE_SAMPLE_TAGS[4];

        let read = concat(&[b"NNNNNNNNNNNNNNNNNNN", tag]);
        match mapper.call(&read) {
            TagCall::NoMatch { .. } => {}
            other => panic!("leading junk must reject, got {other:?}"),
        }
    }

    #[test]
    fn truncated_after_most_of_tag_is_accepted() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        let tag = BD_MOUSE_SAMPLE_TAGS[4];

        let read = concat(&[&tag[..40], b"NNNNNNNNNNNNNNNNNNN"]);
        match mapper.call(&read) {
            TagCall::NoMatch { .. } => {}
            other => panic!("too short insert is accepted"),
        }
    }

    #[test]
    fn small_leading_damage_is_rejected() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        let tag = BD_MOUSE_SAMPLE_TAGS[4];

        let read = concat(&[b"NNNN", &tag[4..], b"NNNNNNNNNNNNNNNNNNN"]);

        match mapper.call(&read) {
            TagCall::NoMatch { .. } => {}
            other => panic!("leading unknown bases before partial tag must reject, got {other:?}"),
        }
    }

    #[test]
    fn too_short_partial_tag_is_rejected() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        let tag = BD_MOUSE_SAMPLE_TAGS[4];

        let read = concat(&[&tag[..24], b"NNNNNNNNNNNNNNNNNNN"]);
        match mapper.call(&read) {
            TagCall::NoMatch { .. } => {}
            other => panic!("too short partial tag must reject, got {other:?}"),
        }
    }

    #[test]
    fn one_mismatches_are_accepted() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        let tag = BD_MOUSE_SAMPLE_TAGS[4];

        let read = concat(&[&mutate_many(tag, &[10]), b"NNNNNNNNNNNN"]);
        expect_match(mapper.call(&read), 4);
    }

    #[test]
    fn three_mismatches_are_rejected() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        let tag = BD_MOUSE_SAMPLE_TAGS[4];

        let read = concat(&[&mutate_many(tag, &[10, 20, 32]), b"NNNNNNNNNNNN"]);
        match mapper.call(&read) {
            TagCall::NoMatch { .. } => {}
            other => panic!("three mismatches must reject, got {other:?}"),
        }
    }

    #[test]
    fn exact_human_tags_do_not_match_mouse_mapper() {
        let mapper = FastTagMapper::bd_mouse(4).unwrap();
        for tag in BD_HUMAN_SAMPLE_TAGS {
            match mapper.call(tag) {
                TagCall::NoMatch { .. } => {}
                other => panic!("human tag matched mouse mapper: {other:?}"),
            }
        }
    }

    #[test]
    fn exact_mouse_tags_do_not_match_human_mapper() {
        let mapper = FastTagMapper::bd_human(4).unwrap();
        for tag in BD_MOUSE_SAMPLE_TAGS {
            match mapper.call(tag) {
                TagCall::NoMatch { .. } => {}
                other => panic!("mouse tag matched human mapper: {other:?}"),
            }
        }
    }
}
