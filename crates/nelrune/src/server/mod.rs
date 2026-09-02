//! Backward-compatible Nelrune names for the workspace-wide Lumrik status server.
//!
//! New code should import `lumrik_status` directly.  Keeping these aliases
//! avoids breaking older callers while the server itself is no longer owned by
//! the Nelrune application.

pub use lumrik_status::{StatusServer as HealthServer, spawn_status_server as spawn_health_server};
