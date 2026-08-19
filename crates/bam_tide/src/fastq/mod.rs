pub mod record;
pub mod writer;
pub mod reader;

pub use record::FastqRecord;
pub use writer::FastqWriter;
pub use reader::{FastqPairReader, SimpleFastqReader, FastqRead};