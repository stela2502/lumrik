//clone_data.rs

use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

pub struct RegionBounds {
    pub cdr1_start: usize,
    pub cdr1_end: usize,
    pub cdr2_start: usize,
    pub cdr2_end: usize,
    pub cdr3_start: usize,
    pub cdr3_end: usize,
}

pub struct CloneData {
    pub dna: Vec<String>,  // unique DNA sequences
    pub aa: Vec<String>,   // AA translated from DNA
    pub count: Vec<usize>, // counts per DNA sequence

    // CDR1-3 positions
    pub bounds: Option<RegionBounds>,
}

impl Default for CloneData {
    fn default() -> Self {
        Self {
            dna: Vec::new(),
            aa: Vec::new(),
            count: Vec::new(),
            bounds: None,
        }
    }
}

impl CloneData {
    pub fn from_raw(raw_dna: &[String]) -> Self {
        use std::collections::HashMap;

        let mut map: HashMap<String, usize> = HashMap::new();

        // Collapse duplicate DNA sequences
        for seq in raw_dna {
            *map.entry(seq.clone()).or_insert(0) += 1;
        }

        // Build dna + count vectors
        let mut dna = Vec::new();
        let mut count = Vec::new();

        for (seq, cnt) in map.into_iter() {
            dna.push(seq);
            count.push(cnt);
        }

        // Translate DNA → AA
        let aa: Vec<String> = dna.iter().map(|s| translate_imgt_alignment(s)).collect();

        CloneData {
            dna,
            aa,
            count,
            bounds: None,
        }
    }

    /// Construct row metadata for caller-defined, already unique feature states.
    /// No deduplication is performed; row order is preserved exactly.
    pub fn from_unique_labels(labels: Vec<String>) -> Self {
        let count = vec![1; labels.len()];
        let aa = vec![String::new(); labels.len()];
        Self {
            dna: labels,
            aa,
            count,
            bounds: None,
        }
    }

    /// Write dna, aa, count as a TSV (3-column table)
    pub fn to_tsv<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn Error>> {
        self.to_delimited('\t', path)?;
        Ok(())
    }

    /// General 3-column delimited writer (dna, aa, count)
    pub fn to_delimited<P: AsRef<Path>>(&self, sep: char, path: P) -> Result<(), Box<dyn Error>> {
        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        let rows = self.dna.len();

        // -----------------------------
        // Validate row-aligned fields
        // -----------------------------
        if self.aa.len() != rows || self.count.len() != rows {
            return Err(format!(
                "Row mismatch: dna={}, aa={}, count={}",
                rows,
                self.aa.len(),
                self.count.len(),
            )
            .into());
        }

        // -----------------------------
        // Header
        // -----------------------------
        writeln!(w, "dna{}aa{}count", sep, sep)?;

        // -----------------------------
        // Rows
        // -----------------------------
        for i in 0..rows {
            writeln!(
                w,
                "{}{}{}{}{}",
                self.dna[i], sep, self.aa[i], sep, self.count[i]
            )?;
        }

        Ok(())
    }

    pub fn len(&self) -> usize {
        self.dna.len()
    }

    /// Collapse DNA sequences by their translated AA sequence.
    /// After this call, CloneData contains one DNA sequence per unique AA sequence,
    /// with counts summed across all DNA sequences that translated to that AA.
    pub fn make_unique_aa(&mut self) {
        use std::collections::HashMap;

        // Map: AA -> (total_count, representative_DNA)
        let mut map: HashMap<String, (usize, String)> = HashMap::new();

        for i in 0..self.aa.len() {
            let aa_seq = &self.aa[i];
            let dna_seq = &self.dna[i];
            let cnt = self.count[i];

            map.entry(aa_seq.clone())
                .and_modify(|(total, _rep_dna)| {
                    *total += cnt; // accumulate counts
                })
                .or_insert((cnt, dna_seq.clone())); // keep first DNA as representative
        }

        // Rebuild CloneData vectors in AA-unique form
        let mut new_dna = Vec::new();
        let mut new_aa = Vec::new();
        let mut new_count = Vec::new();

        for (aa_seq, (total_count, rep_dna)) in map.into_iter() {
            new_dna.push(rep_dna);
            new_aa.push(aa_seq);
            new_count.push(total_count);
        }

        self.dna = new_dna;
        self.aa = new_aa;
        self.count = new_count;
        self.bounds = None;
    }

    /// Identify Cys104 and J118 anchor using IMGT-like consensus across whole clone set.
    ///
    /// Returns (cdr3_start, cdr3_end) indices for *each* sequence.
    /// If detection fails, returns None for that sequence.
    pub fn detect_cdr3_bounds_multi(&self) -> (usize, usize) {
        let aa_seqs = &self.aa;
        let n = aa_seqs.len();
        let max_len = aa_seqs.iter().map(|s| s.len()).max().unwrap();

        // -------------------------------------------------------------
        // 1. Find stable Cys (C104) = position with highest C-frequency
        // -------------------------------------------------------------
        let mut c_freq = vec![0usize; max_len];

        for seq in aa_seqs {
            for (i, ch) in seq.chars().enumerate() {
                if ch == 'C' {
                    c_freq[i] += 1;
                }
            }
        }

        let (cdr3_start, c_score) = c_freq
            .iter()
            .enumerate()
            .max_by_key(|&(_, &count)| count)
            .expect("no positions?");

        // Require this to be a stable cysteine
        if *c_score < n / 4 {
            panic!("No stable Cys104-like position found across sequences");
        }

        // -------------------------------------------------------------
        // 2. Find J118 anchor: look from the RIGHT side for W/F/G cluster
        // -------------------------------------------------------------
        let mut anchor_freq = vec![0usize; max_len];

        for seq in aa_seqs {
            let bytes = seq.as_bytes();
            for i in 0..bytes.len() {
                let aa = bytes[bytes.len() - 1 - i] as char;
                if aa == 'W' || aa == 'F' || aa == 'G' {
                    anchor_freq[bytes.len() - 1 - i] += 1;
                }
            }
        }

        let (cdr3_end, a_score) = anchor_freq
            .iter()
            .enumerate()
            .max_by_key(|&(_, &count)| count)
            .expect("no anchor residues found");

        if *a_score < n / 5 {
            panic!("No stable J-anchor residue (W/F/G) found near end");
        }

        if cdr3_end <= cdr3_start {
            panic!("CDR3 end < start: alignment broken?");
        }

        (cdr3_start, cdr3_end)
    }
    // the keys for the Newton tree
    pub fn aa_with_count_labels(&self) -> Vec<String> {
        self.aa
            .iter()
            .zip(self.count.iter())
            .map(|(aa, count)| format!("{}.{}", aa, count))
            .collect()
    }
}

fn translate_maybe_gap_codon(codon: &[char]) -> char {
    // If one or more positions was a gap / missing nucleotide → X
    if codon.contains(&'N') {
        return 'X';
    }

    let codon_string: String = codon.iter().collect();
    translate_codon(&codon_string)
}

/// Translate an IMGT-gapped nucleotide alignment into an aligned AA sequence.
///
/// Rules:
/// - Codons consist ONLY of nucleotides {A,C,G,T,N}.
/// - Any gap '.' or '-' is *ignored* for codon construction.
/// - Every 3 real nucleotides produce one amino acid.
/// - If the alignment has gap characters at the amino-acid level,
///   they must be inserted explicitly by caller (not here).
/// - If we encounter less than 3 nucleotides at the end → emit 'X'.
/// - If a region contains only gaps → emit 'X'.
pub fn translate_imgt_alignment(seq: &str) -> String {
    let mut aa = String::new();
    let mut codon: Vec<char> = Vec::new();

    for c in seq.chars() {
        match c {
            'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't' => {
                codon.push(c.to_ascii_uppercase());
            }
            '.' | '-' => {
                // gap = missing nucleotide, use placeholder N
                codon.push('N');
            }
            _ => {
                // unknown symbol → treat as N
                codon.push('N');
            }
        }

        if codon.len() == 3 {
            aa.push(translate_maybe_gap_codon(&codon));
            codon.clear();
        }
    }

    // leftover nucleotides (<3) → missing codon → X
    if !codon.is_empty() {
        aa.push('X');
    }

    aa
}

/// Translate a codon to an amino acid
fn translate_codon(c: &str) -> char {
    match &c.to_ascii_uppercase()[..] {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        "TAA" | "TAG" | "TGA" => '*', // stop
        _ => 'X',                     // unknown / incomplete codon
    }
}

#[test]
fn test_imgt_long_gap_yields_long_x_prefix() {
    // The IMGT-gapped alignment
    let nt = "\
................................................................................................................................................................................................................................................................AAGCAGTGGTATCAACGCAGAGTTCAGTGGGGGAGGACACAGCCCTTTATTACTGTGCAAGACGGGATTACTACGGTAGCCACTACTTTGACTACTGGGGCCAAGACACCACTCTCACAGTCTCCTCAG";

    let AA_TAIL: &str = "AVVSTQSSVGEDTALYYCARRDYYGSHYFDYWGQDTTLTVSS";
    // Count leading gaps
    let n_gaps = nt.chars().take_while(|c| *c == '.').count();

    // 3 gaps = 1 codon = 1 missing AA = X
    let expected_x_codons = ((n_gaps as f32) / 3.0).ceil() as usize;

    // Build expected aligned AA sequence
    let mut expected = "X".repeat(expected_x_codons);
    expected.push_str(AA_TAIL);

    // Perform translation
    let aa = translate_imgt_alignment(nt);

    // 1. Leading X's are correct
    assert!(
        aa.starts_with(&"X".repeat(expected_x_codons)),
        "AA does not start with expected X prefix.\nExpected prefix: {}\nGot: {}",
        "X".repeat(expected_x_codons),
        &aa[..expected_x_codons.min(aa.len())]
    );

    // 2. Tail AA matches the known biological translation
    let tail_start = expected_x_codons;
    let aa_tail = &aa[tail_start..tail_start + AA_TAIL.len()];

    assert_eq!(
        aa_tail, AA_TAIL,
        "Translated AA tail does not match expected.\nAA_TAIL: {}\nExpect: {}",
        aa_tail, AA_TAIL
    );
}
