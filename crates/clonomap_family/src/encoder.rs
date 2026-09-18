use crate::CloneData;
use ndarray::Array2;
use std::error::Error;

/// OneHotEncoder:
///  - Keeps existing functionality (same API!).
///  - Now also computes:
///       * aligned DNA
///       * ungapped coding DNA
///       * best reading frame
///       * amino-acid sequence
///       * mutation maps
///       * DNA & AA mutation lists
///
/// No behavior visible to your existing code is broken.
pub struct OneHotEncoder {
    pub seq_len: usize,

    // Consensus base for each aligned position
    pub consensus: Vec<u8>,

    // Original aligned sequences
    pub sequences: CloneData,

    // Relative 0/1/-1 PCA-friendly encoding
    pub encoded_states: Vec<Vec<f32>>,

    // NEW: Ungapped coding sequences (one per input sequence)
    pub coding_sequences: Vec<String>,

    // NEW: Best frame for each sequence (0,1,2)
    pub best_frames: Vec<usize>,
}

impl OneHotEncoder {
    pub fn new() -> Self {
        Self {
            seq_len: 0,
            consensus: Vec::new(),
            sequences: CloneData::default(),
            encoded_states: Vec::new(),
            coding_sequences: Vec::new(),
            best_frames: Vec::new(),
        }
    }

    // ========================================================================
    //  THE IMPORTANT ONE — called by PCA pipeline
    // ========================================================================
    pub fn encode_relative(
        &mut self,
        sequences: &[String],
        aa_based: bool,
    ) -> Result<Array2<f32>, Box<dyn Error>> {
        if sequences.is_empty() {
            return Err("No sequences provided".into());
        }

        // --------------------------------------------------------------------
        // Validate and store aligned sequences
        // --------------------------------------------------------------------
        let alignment_width = sequences[0].len();
        self.seq_len = alignment_width;

        for (i, s) in sequences.iter().enumerate() {
            if s.len() != alignment_width {
                return Err(format!(
                    "Sequence length mismatch at index {}: expected {}, got {}",
                    i,
                    alignment_width,
                    s.len()
                )
                .into());
            }
        }
        let mut cdata = CloneData::from_raw(sequences);

        // collapse to unique AA set
        if aa_based {
            cdata.make_unique_aa();
        }
        self.sequences = cdata;

        self.consensus.clear();
        self.encoded_states.clear();
        self.coding_sequences.clear();
        self.best_frames.clear();

        let n = self.sequences.len();
        let mut x = Array2::<f32>::zeros((n, alignment_width));

        // --------------------------------------------------------------------
        // Build DNA consensus per position
        // --------------------------------------------------------------------
        for col in 0..alignment_width {
            let mut counts = [0u32; 4];

            for s in &self.sequences.dna {
                match s.as_bytes()[col] {
                    b'A' => counts[0] += 1,
                    b'C' => counts[1] += 1,
                    b'G' => counts[2] += 1,
                    b'T' => counts[3] += 1,
                    _ => {} // ignore gaps
                }
            }

            let (idx, _) = counts.iter().enumerate().max_by_key(|(_, c)| *c).unwrap();
            self.consensus.push(b"ACGT"[idx]);
        }

        // --------------------------------------------------------------------
        // Build PCA-friendly relative encoding  (mutation = 1, consensus = 0)
        // --------------------------------------------------------------------
        for (i, seq) in self.sequences.dna.iter().enumerate() {
            let mut row_state = Vec::with_capacity(alignment_width);

            for (j, base) in seq.as_bytes().iter().enumerate() {
                let val = match *base {
                    b'.' | b'-' => -1.0, // gap
                    b => {
                        if b == self.consensus[j] {
                            0.0
                        } else {
                            1.0
                        }
                    }
                };
                x[[i, j]] = val;
                row_state.push(val);
            }

            self.encoded_states.push(row_state);
        }

        // --------------------------------------------------------------------
        // Produce ungapped coding DNA for frame detection
        // --------------------------------------------------------------------
        for seq in &self.sequences.dna {
            let coding = seq
                .chars()
                .filter(|c| *c != '-' && *c != '.')
                .collect::<String>();
            self.coding_sequences.push(coding);
        }

        Ok(x)
    }

    /// Install already encoded states while keeping the legacy row metadata
    /// used by ClonoMap plotting/export code.
    ///
    /// Unlike `encode_relative`, this does not derive features from nucleotide
    /// consensus. The caller owns the feature semantics. This is intentionally
    /// useful for composite receptor states such as HC mutations + LC mutations
    /// + a weighted one-hot structural LC identity.
    pub fn set_external_states(
        &mut self,
        row_labels: Vec<String>,
        features: &Array2<f32>,
    ) -> Result<(), Box<dyn Error>> {
        if features.nrows() != row_labels.len() {
            return Err(format!(
                "External state mismatch: {} labels for {} feature rows",
                row_labels.len(),
                features.nrows()
            )
            .into());
        }
        if row_labels.iter().any(|x| x.is_empty()) {
            return Err("External state labels must not be empty".into());
        }

        self.seq_len = 0;
        self.consensus.clear();
        self.coding_sequences.clear();
        self.best_frames.clear();
        self.encoded_states = features.outer_iter().map(|row| row.to_vec()).collect();

        // CloneData is used as row metadata by the existing PCA/tree exporters.
        // Do not collapse these labels: each feature row is already a unique
        // caller-defined state.
        self.sequences = CloneData::from_unique_labels(row_labels);
        Ok(())
    }

    // ========================================================================
    //  Helper functions for PCA annotation (optional use by caller)
    // ========================================================================

    /// Return human-readable mutation pattern ".X..X..."
    pub fn mutation_pattern(&self, row: usize) -> String {
        self.encoded_states[row]
            .iter()
            .map(|v| {
                if *v == 0.0 {
                    '.'
                } else if *v == 1.0 {
                    'X'
                } else {
                    '-'
                }
            })
            .collect()
    }

    /// Return number of mutated DNA positions
    pub fn mutation_count(&self, row: usize) -> usize {
        self.encoded_states[row]
            .iter()
            .filter(|v| **v == 1.0)
            .count()
    }
}
