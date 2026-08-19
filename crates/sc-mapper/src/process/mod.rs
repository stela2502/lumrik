// crates/sc-mapper/src/process/mod.rs

mod mapper_process;
mod sam_cluster_buffer;
mod mapper_record;

mod bwa;
mod minimap2;
mod star;

pub use mapper_process::{check_binary, MapperProcess};

pub use sam_cluster_buffer::{
    sam_cluster_channel,
    SamClusterBuffer,
    SamClusterReceiver,
    SamClusterSender,
};

pub use mapper_record::{
    MapperRecord,
    SamReadCluster,
};

pub use bwa::Bwa;
pub use minimap2::Minimap2;
pub use star::Star;
