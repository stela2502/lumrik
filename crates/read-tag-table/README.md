# read-tag-table

A lightweight Rust crate for storing, serializing, and exchanging single-cell read tag information.

`read-tag-table` provides a common representation for read-level cell barcode and UMI assignments and was extracted from the bam_tide ecosystem to serve as a reusable interchange format between tools.

The crate supports:

* Binary serialization for fast loading and storage
* TSV and TSV.GZ import
* TSV export
* Read-level lookup by read identifier
* Cell/UMI pair statistics
* Integration with single-cell sequencing workflows

---

## Installation

### From crates.io

```bash
cargo add read-tag-table
```

or manually:

```toml
[dependencies]
read-tag-table = "0.1"
```

### From GitHub

```toml
[dependencies]
read-tag-table = { git = "https://github.com/stela2502/read-tag-table.git" }
```

---

## Concepts

### ReadTagRecord

A `ReadTagRecord` stores the barcode and UMI assignment for a single sequencing read.

Each record contains:

* read identifier
* optional original read identifier
* cell barcode sequence
* cell barcode quality values
* UMI sequence
* UMI quality values

Example:

```rust
use read_tag_table::ReadTagRecord;

let record = ReadTagRecord::new(
    "read1".to_string(),
    None,
    b"CELL1",
    vec![40,40,40,40,40],
    b"UMI1",
    vec![40,40,40,40],
);
```

---

### ReadTagTable

A `ReadTagTable` stores many `ReadTagRecord` objects and provides lookup and summary functionality.

```rust
use read_tag_table::{ReadTagRecord, ReadTagTable};

let mut table = ReadTagTable::new();

table.insert(
    ReadTagRecord::new(
        "read1".to_string(),
        None,
        b"CELL1",
        Vec::new(),
        b"UMI1",
        Vec::new(),
    )
);
```

## Command Line Integration

The crate includes a ready-to-use Clap integration.

```rust
use clap::Parser;
use read_tag_table::ReadTagTableCli;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    read_tag_table: ReadTagTableCli,
}
```

Load a table directly:

```rust
let table = cli.read_tag_table.load()?;
```

For workflows using multiple read-tag tables:

```rust
let table = cli.read_tag_table.load_for_id(0)?;
```

Supported inputs:

- Binary `.bin` tables
- TSV files
- TSV.GZ files

The appropriate reader is selected automatically.

---

## Binary Storage

The preferred storage format is the native binary format.

### Save

```rust
table.save_binary("read_tags.bin")?;
```

or:

```rust
table.save("read_tags.bin")?;
```

### Load

```rust
let table = ReadTagTable::load_binary("read_tags.bin")?;
```

The binary format is considerably faster to load than TSV and is recommended for large datasets.

---

## Loading TSV Files

Read-tag tables can be loaded from TSV files.

Required columns:

| Column  |
| ------- |
| read_id |
| raw_cb  |
| raw_umi |

Optional columns:

| Column           |
| ---------------- |
| original_read_id |
| quality_cb       |
| quality_umi      |

Example:

```text
read_id raw_cb  raw_umi
read1   CELL1   UMI1
read2   CELL2   UMI2
```

Load using a configuration:

```rust
use read_tag_table::{ReadTagTable, ReadTagTableConfig};

let config = ReadTagTableConfig {
    path: "read_tags.tsv".into(),
    read_id_column: "read_id".to_string(),
    original_read_id_column: "original_read_id".to_string(),
    cell_column: "raw_cb".to_string(),
    cell_qual_column: "quality_cb".to_string(),
    umi_column: "raw_umi".to_string(),
    umi_qual_column: "quality_umi".to_string(),
};

let table = ReadTagTable::from_config(&config)?;
```

Both plain TSV and TSV.GZ files are supported.

---

## Lookup by Read ID

Retrieve the barcode and UMI associated with a read:

```rust
if let Some((cell, umi)) = table.cell_umi_for_read("read1") {
    println!("{cell} {umi}");
}
```

Or access the full record:

```rust
if let Some(record) = table.get("read1") {
    println!("{}", record.cell_string());
}
```

---

## Merging Tables

Multiple tables can be combined:

```rust
table_a.merge(table_b);
```

Records from the second table overwrite records with identical read identifiers.

---

## Exporting TSV

Write a table as TSV:

```rust
let mut file = std::fs::File::create("read_tags.tsv")?;
table.write_tsv(&mut file)?;
```

The output format is compatible with the TSV reader.

---

## Cell / UMI Statistics

The crate can summarize cell-barcode / UMI observations.

```rust
let stats = table.summarize_pairs(
    2, // minimum observations per cell/UMI pair
    2, // minimum UMIs per cell
);
```

Reported statistics include:

* number of valid cells
* number of unique cell/UMI combinations
* total observations
* UMIs per cell
* observations per cell/UMI pair

This functionality is useful for quality control and library complexity analysis.

---

## Typical Workflow

```text
FASTQ/BAM processing
        ↓
ReadTagRecord generation
        ↓
ReadTagTable
        ↓
save_binary()
        ↓
read_tags.bin
        ↓
load_binary()
        ↓
downstream analysis
```

---

## Intended Use

`read-tag-table` is designed as a lightweight exchange format between tools that generate or consume single-cell barcode assignments.

Examples include:

* FASTQ normalizers
* Primer stripping workflows
* Primer restoration workflows
* Single-cell barcode translation
* Cell/UMI quality control
* VDJ and targeted sequencing pipelines

---

## License

See the repository LICENSE file.
