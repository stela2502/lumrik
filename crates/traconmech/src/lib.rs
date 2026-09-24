//! TraConMech: trace biological context to candidate mechanisms.
//!
//! The crate deliberately starts small. Ommverse supplies the stable biological
//! reference; experimental observations are consumed through thin adapters.
//! Missing modalities must reduce what can be concluded, not make an
//! investigation impossible.

/// Human-readable name used by CLIs and reports.
pub const NAME: &str = "TraConMech";

/// Short description of the project's core job.
pub const TAGLINE: &str = "Trace Context to Mechanism";
