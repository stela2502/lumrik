pub mod int_to_dna;

#[cfg(feature = "alignment")]
pub mod alignment;

#[allow(unused_imports)]
pub use int_to_dna::IntToDna;

#[cfg(feature = "alignment")]
pub use alignment::{AlignmentConfig, AlignmentResult, ExactOverlap};
