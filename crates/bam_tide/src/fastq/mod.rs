pub mod reader;
pub mod record;
pub mod writer;

pub use reader::{FastqPairReader, FastqRead, SimpleFastqReader};
pub use record::FastqRecord;
pub use writer::FastqWriter;
