//! Ommverse: Lumrik's integrated genome-to-protein biological reference model.
//!
//! Ommverse deliberately separates biological identity from source formats.
//! UCSC GTF, twoBit and UniProt bigBed files are import formats; callers see
//! genes, transcripts, proteins and protein features.

pub mod protein_index;
pub use protein_index::{ProteinFeatureId, ProteinFeatureIndex, ProjectedProteinFeature};

use anyhow::{Context, Result, bail};
use bigtools::BigBedRead;
use gtf_splice_index::types::RefBlock;
use gtf_splice_index::{Axis, Connected, Connection, Gene, GeneId, IdNameKeys, Identifiable, Plottable, SpliceIndex, Strand, Transcript, TranscriptId};
use hmm::{CategoricalEmission, Hmm};
use int_to_dna::{IntToDna, TwoBitReader};
use int_to_prot::IntToProt;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 4] = b"OMM1";
pub const OMMVERSE_FORMAT_VERSION: u32 = 10;

pub mod sources;
pub mod query;
pub mod schema;

// Source layout is split by responsibility while deliberately sharing the crate-root
// namespace. `include!` keeps this reorganisation API- and logic-neutral.
include!("model.rs");
include!("build.rs");
include!("lookup.rs");
include!("sequence.rs");
include!("training.rs");
include!("storage.rs");
include!("helpers.rs");
include!("tests.rs");
