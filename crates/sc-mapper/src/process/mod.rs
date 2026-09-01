// crates/sc-mapper/src/process/mod.rs

mod mapper_process;
mod mapper_record;
mod options;
mod sam_cluster_buffer;

mod bwa;
mod minimap2;
mod star;

pub use mapper_process::{MapperProcess, check_binary};

pub use sam_cluster_buffer::{
    SamClusterBuffer, SamClusterReceiver, SamClusterSender, sam_cluster_channel,
};

pub use mapper_record::{MapperRecord, SamReadCluster};

pub use bwa::Bwa;
pub use minimap2::Minimap2;
pub use star::Star;

pub(crate) use options::{has_option, remove_option};
