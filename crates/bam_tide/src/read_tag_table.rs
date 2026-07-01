use anyhow::{Context, Result};
use clap::Args;
use csv::WriterBuilder;
use flate2::read::MultiGzDecoder;
use mapping_info::MappingInfo;
use serde::{Deserialize, Serialize};

use std::{
    collections::{HashMap, HashSet},
    fmt,
    fs::File,
    io::{BufReader, Read, Write},
    path::{Path, PathBuf},
};

#[derive(Debug, Clone, Args)]
pub struct ReadTagTableCli {
    /// Optional external read-tag table.
    ///
    /// Accepted formats:
    /// - binary `.bin` files written with `ReadTagTable::save_binary`
    /// - TSV files
    /// - TSV.GZ files
    #[arg(long = "read-tag-table", num_args = 1..)]
    pub read_tag_table: Vec<PathBuf>,

    #[arg(long = "rt-read-id-column", default_value = "read_id")]
    pub rt_read_id_column: String,

    #[arg(long = "rt-cell-column", default_value = "raw_cb")]
    pub rt_cell_column: String,

    #[arg(long = "rt-cell-qual-column", default_value = "quality_cb")]
    pub rt_cell_qual_column: String,

    #[arg(long = "rt-umi-column", default_value = "raw_umi")]
    pub rt_umi_column: String,

    #[arg(long = "rt-umi-qual-column", default_value = "quality_umi")]
    pub rt_umi_qual_column: String,

    #[arg(long = "rt-original-read-id-column", default_value = "original_read_id")]
    pub rt_original_read_id_column: String,
}

impl ReadTagTableCli {
    pub fn to_config_for_id(&self, id: usize) -> Result<ReadTagTableConfig> {
        if self.read_tag_table.is_empty() {
            anyhow::bail!("No --read-tag-table files supplied");
        }

        let path = self
            .read_tag_table
            .get(id)
            .with_context(|| format!("No --read-tag-table for id {id}"))?
            .clone();

        Ok(self.config_for_path(path))
    }

    pub fn to_config(&self) -> Result<ReadTagTableConfig> {
        if self.read_tag_table.is_empty() {
            anyhow::bail!("No --read-tag-table files supplied");
        }

        if self.read_tag_table.len() > 1 {
            anyhow::bail!(
                "Number of --read-tag-table files (seen {}) > 1 is not supported by this function",
                self.read_tag_table.len()
            );
        }

        Ok(self.config_for_path(self.read_tag_table[0].clone()))
    }

    pub fn load(&self) -> Result<ReadTagTable> {
        let config = self.to_config()?;
        ReadTagTable::from_path_or_config(&config)
    }

    pub fn load_for_id(&self, id: usize) -> Result<ReadTagTable> {
        let config = self.to_config_for_id(id)?;
        ReadTagTable::from_path_or_config(&config)
    }

    fn config_for_path(&self, path: PathBuf) -> ReadTagTableConfig {
        ReadTagTableConfig {
            path,
            read_id_column: self.rt_read_id_column.clone(),
            original_read_id_column: self.rt_original_read_id_column.clone(),
            cell_column: self.rt_cell_column.clone(),
            cell_qual_column: self.rt_cell_qual_column.clone(),
            umi_column: self.rt_umi_column.clone(),
            umi_qual_column: self.rt_umi_qual_column.clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ReadTagTableConfig {
    pub path: PathBuf,
    pub read_id_column: String,
    pub original_read_id_column: String,
    pub cell_column: String,
    pub cell_qual_column: String,
    pub umi_column: String,
    pub umi_qual_column: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadTagRecord {
    pub read_id: String,
    pub original_read_id: Option<String>,
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
}

impl ReadTagRecord {
    pub fn new(
        read_id: String,
        original_read_id: Option<String>,
        cell_seq: impl AsRef<[u8]>,
        cell_qual: impl AsRef<[u8]>,
        umi_seq: impl AsRef<[u8]>,
        umi_qual: impl AsRef<[u8]>,
    ) -> Self {
        Self {
            read_id,
            original_read_id,
            cell_seq: cell_seq.as_ref().to_vec(),
            cell_qual: cell_qual.as_ref().to_vec(),
            umi_seq: umi_seq.as_ref().to_vec(),
            umi_qual: umi_qual.as_ref().to_vec(),
        }
    }

    pub fn from_slices(
        read_id: impl Into<String>,
        original_read_id: Option<String>,
        cell_seq: &[u8],
        cell_qual: &[u8],
        umi_seq: &[u8],
        umi_qual: &[u8],
    ) -> Self {
        Self::new(
            read_id.into(),
            original_read_id,
            cell_seq.to_vec(),
            cell_qual.to_vec(),
            umi_seq.to_vec(),
            umi_qual.to_vec(),
        )
    }

    pub fn from_tsv_fields(
        read_id: impl Into<String>,
        original_read_id: Option<String>,
        cell: &str,
        cell_qual: Option<&str>,
        umi: &str,
        umi_qual: Option<&str>,
    ) -> Self {
        Self::new(
            read_id.into(),
            original_read_id,
            cell.as_bytes().to_vec(),
            cell_qual.map(phred_from_ascii).unwrap_or_default(),
            umi.as_bytes().to_vec(),
            umi_qual.map(phred_from_ascii).unwrap_or_default(),
        )
    }

    pub fn cell_string(&self) -> String {
        seq_to_string(&self.cell_seq)
    }

    pub fn umi_string(&self) -> String {
        seq_to_string(&self.umi_seq)
    }

    pub fn cell_qual_string(&self) -> String {
        phred_to_ascii(&self.cell_qual)
    }

    pub fn umi_qual_string(&self) -> String {
        phred_to_ascii(&self.umi_qual)
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReadTagTable {
    records: HashMap<String, ReadTagRecord>,
}

impl ReadTagTable {
    pub fn new() -> Self {
        Self {
            records: HashMap::new(),
        }
    }

    pub fn from_path_or_config(config: &ReadTagTableConfig) -> Result<Self> {
        if is_binary_read_tag_table(&config.path) {
            Self::load_binary(&config.path)
        } else {
            Self::from_config(config)
        }
    }

    pub fn from_config(config: &ReadTagTableConfig) -> Result<Self> {
        let mut mapping_info = MappingInfo::new(None, 0.0, 0);
        let reader = open_maybe_gz(&config.path)?;

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_reader(reader);

        let headers = rdr.headers()?.clone();

        let read_id_ix = required_column_ix(&headers, &config.read_id_column)?;
        let cell_ix = required_column_ix(&headers, &config.cell_column)?;
        let umi_ix = required_column_ix(&headers, &config.umi_column)?;

        let original_read_id_ix = optional_column_ix(&headers, &config.original_read_id_column);
        let cell_qual_ix = optional_column_ix(&headers, &config.cell_qual_column);
        let umi_qual_ix = optional_column_ix(&headers, &config.umi_qual_column);

        let mut records = HashMap::new();

        for rec in rdr.records() {
            let rec = rec?;

            let read_id = rec.get(read_id_ix).unwrap_or("").trim();
            let cell = rec.get(cell_ix).unwrap_or("").trim();
            let umi = rec.get(umi_ix).unwrap_or("").trim();

            if read_id.is_empty() || cell.is_empty() || umi.is_empty() {
                continue;
            }

            let record = ReadTagRecord::from_tsv_fields(
                read_id.to_string(),
                get_optional(&rec, original_read_id_ix),
                cell,
                get_optional_ref(&rec, cell_qual_ix),
                umi,
                get_optional_ref(&rec, umi_qual_ix),
            );

            records.insert(read_id.to_string(), record);
        }

        mapping_info.stop_file_io_time();

        let (h, m, s, ms) = MappingInfo::split_duration(mapping_info.file_io_time);
        println!(
            "Read-tag table loaded from TSV: {} entries in {}:{:02}:{:02}.{:03}",
            records.len(),
            h,
            m,
            s,
            ms,
        );

        Ok(Self { records })
    }

    pub fn save_binary<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;

        bincode::serialize_into(file, self)
            .with_context(|| format!("writing binary read-tag table {}", path.display()))?;

        Ok(())
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        self.save_binary(path)
    }

    pub fn load_binary<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;

        let table: Self = bincode::deserialize_from(file)
            .with_context(|| format!("reading binary read-tag table {}", path.display()))?;

        Ok(table)
    }

    pub fn get(&self, read_id: &str) -> Option<&ReadTagRecord> {
        self.records.get(read_id)
    }

    pub fn cell_umi_for_read(&self, read_id: &str) -> Option<(String, String)> {
        let rec = self.records.get(read_id)?;
        Some((rec.cell_string(), rec.umi_string()))
    }

    pub fn cell_umi_bytes_for_read(&self, read_id: &str) -> Option<(&[u8], &[u8])> {
        let rec = self.records.get(read_id)?;
        Some((rec.cell_seq.as_slice(), rec.umi_seq.as_slice()))
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn pair_counts(&self) -> HashMap<(Vec<u8>, Vec<u8>), u64> {
        let mut counts = HashMap::new();

        for rec in self.records.values() {
            if rec.cell_seq.is_empty() || rec.umi_seq.is_empty() {
                continue;
            }

            *counts
                .entry((rec.cell_seq.clone(), rec.umi_seq.clone()))
                .or_insert(0) += 1;
        }

        counts
    }

    pub fn insert(&mut self, record: ReadTagRecord) -> Option<ReadTagRecord> {
        self.records.insert(record.read_id.clone(), record)
    }

    pub fn merge(&mut self, other: Self) {
        self.records.extend(other.records);
    }

    pub fn summarize_pairs(&self, min_pair_count: u64, min_cell_umis: u64) -> PairStats {
        let counts = self.pair_counts();
        PairStats::from_counts(&counts, min_pair_count, min_cell_umis)
    }

    pub fn write_tsv<W: Write>(&self, inner: W) -> Result<()> {
        let mut writer = ReadTagTableWriter::new(inner)?;

        for record in self.records.values() {
            writer.write_record(record)?;
        }

        writer.flush()
    }
}

fn is_binary_read_tag_table(path: &Path) -> bool {
    path.extension()
        .and_then(|x| x.to_str())
        .is_some_and(|x| x.eq_ignore_ascii_case("bin"))
}

#[derive(Debug, Clone)]
pub struct PairStats {
    pub cell_entries: usize,
    pub unique_cell_umi_combos: usize,
    pub total_pair_observations: u64,
    pub umis_per_cell: Summary,
    pub detections_per_cell_umi: Summary,
}

impl PairStats {
    pub fn from_counts(
        cb_umi_counts: &HashMap<(Vec<u8>, Vec<u8>), u64>,
        min_pair_count: u64,
        min_cell_umis: u64,
    ) -> Self {
        let mut cell_to_umis: HashMap<&[u8], HashSet<&[u8]>> = HashMap::new();

        for ((cell, umi), count) in cb_umi_counts {
            if *count < min_pair_count {
                continue;
            }

            cell_to_umis
                .entry(cell.as_slice())
                .or_default()
                .insert(umi.as_slice());
        }

        let valid_cells: HashSet<&[u8]> = cell_to_umis
            .iter()
            .filter_map(|(cell, umis)| {
                if umis.len() as u64 >= min_cell_umis {
                    Some(*cell)
                } else {
                    None
                }
            })
            .collect();

        let mut final_cell_to_umi_count: HashMap<&[u8], u64> = HashMap::new();
        let mut final_pair_counts = Vec::new();
        let mut total_pair_observations = 0_u64;

        for ((cell, _umi), count) in cb_umi_counts {
            if *count < min_pair_count {
                continue;
            }

            if !valid_cells.contains(cell.as_slice()) {
                continue;
            }

            *final_cell_to_umi_count.entry(cell.as_slice()).or_insert(0) += 1;
            final_pair_counts.push(*count);
            total_pair_observations += *count;
        }

        let umis_per_cell: Vec<u64> = final_cell_to_umi_count.values().copied().collect();

        Self {
            cell_entries: final_cell_to_umi_count.len(),
            unique_cell_umi_combos: final_pair_counts.len(),
            total_pair_observations,
            umis_per_cell: summarize(umis_per_cell),
            detections_per_cell_umi: summarize(final_pair_counts),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Summary {
    pub n: usize,
    pub mean: f64,
    pub median: f64,
    pub min: u64,
    pub max: u64,
}

fn summarize(mut values: Vec<u64>) -> Summary {
    if values.is_empty() {
        return Summary {
            n: 0,
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
        };
    }

    values.sort_unstable();

    let n = values.len();
    let min = values[0];
    let max = values[n - 1];

    let sum: u128 = values.iter().map(|&x| x as u128).sum();
    let mean = sum as f64 / n as f64;

    let median = if n % 2 == 0 {
        let a = values[n / 2 - 1];
        let b = values[n / 2];
        (a as f64 + b as f64) / 2.0
    } else {
        values[n / 2] as f64
    };

    Summary {
        n,
        mean,
        median,
        min,
        max,
    }
}

fn required_column_ix(headers: &csv::StringRecord, name: &str) -> Result<usize> {
    headers
        .iter()
        .position(|h| h == name)
        .with_context(|| format!("Could not find required column '{name}'"))
}

fn optional_column_ix(headers: &csv::StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|h| h == name)
}

fn get_optional(rec: &csv::StringRecord, ix: Option<usize>) -> Option<String> {
    let value = rec.get(ix?)?.trim();

    if value.is_empty() {
        None
    } else {
        Some(value.to_string())
    }
}

fn get_optional_ref<'a>(rec: &'a csv::StringRecord, ix: Option<usize>) -> Option<&'a str> {
    let value = rec.get(ix?)?.trim();

    if value.is_empty() {
        None
    } else {
        Some(value)
    }
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(reader)))
    } else {
        Ok(Box::new(reader))
    }
}

fn seq_to_string(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).to_string()
}

fn phred_from_ascii(text: &str) -> Vec<u8> {
    text.as_bytes()
        .iter()
        .map(|q| q.saturating_sub(33))
        .collect()
}

fn phred_to_ascii(qual: &[u8]) -> String {
    qual.iter()
        .map(|q| q.saturating_add(33) as char)
        .collect()
}

fn display_qual(qual: &[u8]) -> String {
    if qual.is_empty() {
        "-".to_string()
    } else {
        phred_to_ascii(qual)
    }
}

pub const READ_TAG_TABLE_COLUMNS: [&str; 8] = [
    "read_id",
    "original_read_id",
    "orientation",
    "raw_cb",
    "quality_cb",
    "raw_umi",
    "quality_umi",
    "status",
];

#[derive(Debug, Clone)]
pub struct ReadTagWriteRecord<'a> {
    pub read_id: &'a str,
    pub original_read_id: Option<&'a str>,
    pub orientation: Option<&'a str>,
    pub raw_cb: String,
    pub quality_cb: String,
    pub raw_umi: String,
    pub quality_umi: String,
    pub status: &'a str,
}

pub struct ReadTagTableWriter<W: Write> {
    writer: csv::Writer<W>,
}

impl<W: Write> ReadTagTableWriter<W> {
    pub fn new(inner: W) -> Result<Self> {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(inner);

        writer
            .write_record(READ_TAG_TABLE_COLUMNS)
            .context("writing read-tag table header")?;

        Ok(Self { writer })
    }

    pub fn write_record(&mut self, rec: &ReadTagRecord) -> Result<()> {
        self.writer
            .write_record([
                rec.read_id.as_str(),
                rec.original_read_id.as_deref().unwrap_or(""),
                "",
                rec.cell_string().as_str(),
                rec.cell_qual_string().as_str(),
                rec.umi_string().as_str(),
                rec.umi_qual_string().as_str(),
                "ok",
            ])
            .with_context(|| format!("writing read-tag table row for read_id '{}'", rec.read_id))?;

        Ok(())
    }

    /// Compatibility shim for older TSV-streaming code.
    /// Normalizers should prefer `ReadTagTable::insert()` + `ReadTagTable::save()`.
    pub fn write_tsv_record(&mut self, rec: &ReadTagWriteRecord<'_>) -> Result<()> {
        self.writer
            .write_record([
                rec.read_id,
                rec.original_read_id.unwrap_or(""),
                rec.orientation.unwrap_or(""),
                &rec.raw_cb,
                &rec.quality_cb,
                &rec.raw_umi,
                &rec.quality_umi,
                rec.status,
            ])
            .with_context(|| format!("writing read-tag table row for read_id '{}'", rec.read_id))?;

        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("flushing read-tag table writer")?;
        Ok(())
    }
}

impl fmt::Display for ReadTagRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "read_id={}, original_read_id={}, cell={}, cell_qual={}, umi={}, umi_qual={}",
            self.read_id,
            self.original_read_id.as_deref().unwrap_or("-"),
            self.cell_string(),
            display_qual(&self.cell_qual),
            self.umi_string(),
            display_qual(&self.umi_qual),
        )
    }
}

impl fmt::Display for ReadTagTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "ReadTagTable: {} records", self.records.len())?;

        for (i, record) in self.records.values().take(3).enumerate() {
            writeln!(f, "  [{}] {}", i, record)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tmp_path(name: &str) -> PathBuf {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();

        std::env::temp_dir().join(format!("read_tag_table_test_{stamp}_{name}"))
    }

    fn write_text(path: &Path, text: &str) {
        fs::write(path, text).unwrap();
    }

    fn default_config(path: PathBuf) -> ReadTagTableConfig {
        ReadTagTableConfig {
            path,
            read_id_column: "read_id".to_string(),
            original_read_id_column: "original_read_id".to_string(),
            cell_column: "raw_cb".to_string(),
            cell_qual_column: "quality_cb".to_string(),
            umi_column: "raw_umi".to_string(),
            umi_qual_column: "quality_umi".to_string(),
        }
    }

    #[test]
    fn reads_plain_tsv_table() {
        let path = tmp_path("table.tsv");

        write_text(
            &path,
            concat!(
                "read_id\toriginal_read_id\traw_cb\tquality_cb\traw_umi\tquality_umi\n",
                "read1\torig1\tCELL1\tIIII\tUMI1\tJJJJ\n",
                "read2\torig2\tCELL1\tIIII\tUMI2\tJJJJ\n",
                "read3\t\tCELL2\t\tUMI3\t\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        assert_eq!(table.len(), 3);
        assert!(!table.is_empty());

        let rec = table.get("read1").unwrap();
        assert_eq!(rec.read_id, "read1");
        assert_eq!(rec.original_read_id.as_deref(), Some("orig1"));
        assert_eq!(rec.cell_seq, b"CELL1");
        assert_eq!(rec.cell_qual, vec![40, 40, 40, 40]);
        assert_eq!(rec.umi_seq, b"UMI1");
        assert_eq!(rec.umi_qual, vec![41, 41, 41, 41]);

        let rec = table.get("read3").unwrap();
        assert_eq!(rec.original_read_id, None);
        assert!(rec.cell_qual.is_empty());
        assert!(rec.umi_qual.is_empty());

        fs::remove_file(path).ok();
    }

    #[test]
    fn skips_rows_with_missing_required_values() {
        let path = tmp_path("missing.tsv");

        write_text(
            &path,
            concat!(
                "read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "\tCELL2\tUMI2\n",
                "read3\t\tUMI3\n",
                "read4\tCELL4\t\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        assert_eq!(table.len(), 1);
        assert!(table.get("read1").is_some());
        assert!(table.get("read3").is_none());

        fs::remove_file(path).ok();
    }

    #[test]
    fn missing_required_column_returns_error() {
        let path = tmp_path("bad_columns.tsv");

        write_text(&path, concat!("read_id\traw_cb\n", "read1\tCELL1\n",));

        let err = ReadTagTable::from_config(&default_config(path.clone()))
            .unwrap_err()
            .to_string();

        assert!(
            err.contains("Could not find required column 'raw_umi'"),
            "unexpected error: {err}"
        );

        fs::remove_file(path).ok();
    }

    #[test]
    fn reads_gzipped_tsv_table() {
        let path = tmp_path("table.tsv.gz");

        let file = File::create(&path).unwrap();
        let mut gz = GzEncoder::new(file, Compression::default());

        gz.write_all(
            concat!(
                "read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "read2\tCELL2\tUMI2\n",
            )
            .as_bytes(),
        )
        .unwrap();

        gz.finish().unwrap();

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        assert_eq!(table.len(), 2);
        assert_eq!(
            table.cell_umi_for_read("read1"),
            Some(("CELL1".to_string(), "UMI1".to_string()))
        );
        assert_eq!(table.cell_umi_for_read("missing"), None);

        fs::remove_file(path).ok();
    }

    #[test]
    fn loads_binary_when_extension_is_bin() {
        let path = tmp_path("read_tags.bin");

        let mut table = ReadTagTable::new();
        table.insert(ReadTagRecord::new(
            "read1".to_string(),
            Some("orig1".to_string()),
            b"CELL1".to_vec(),
            vec![40, 40, 40, 40],
            b"UMI1".to_vec(),
            vec![41, 41, 41, 41],
        ));

        table.save_binary(&path).unwrap();

        let loaded = ReadTagTable::from_path_or_config(&default_config(path.clone())).unwrap();

        assert_eq!(loaded.len(), 1);
        assert_eq!(
            loaded.cell_umi_for_read("read1"),
            Some(("CELL1".to_string(), "UMI1".to_string()))
        );

        fs::remove_file(path).ok();
    }

    #[test]
    fn cli_load_uses_binary_for_bin_extension() {
        let path = tmp_path("cli_read_tags.bin");

        let mut table = ReadTagTable::new();
        table.insert(ReadTagRecord::new(
            "read1".to_string(),
            None,
            b"CELL1".to_vec(),
            Vec::new(),
            b"UMI1".to_vec(),
            Vec::new(),
        ));

        table.save_binary(&path).unwrap();

        let cli = ReadTagTableCli {
            read_tag_table: vec![path.clone()],
            rt_read_id_column: "read_id".to_string(),
            rt_cell_column: "raw_cb".to_string(),
            rt_cell_qual_column: "quality_cb".to_string(),
            rt_umi_column: "raw_umi".to_string(),
            rt_umi_qual_column: "quality_umi".to_string(),
            rt_original_read_id_column: "original_read_id".to_string(),
        };

        let loaded = cli.load().unwrap();

        assert_eq!(loaded.len(), 1);
        assert_eq!(
            loaded.cell_umi_for_read("read1"),
            Some(("CELL1".to_string(), "UMI1".to_string()))
        );

        fs::remove_file(path).ok();
    }

    #[test]
    fn pair_counts_counts_cell_umi_observations() {
        let path = tmp_path("pairs.tsv");

        write_text(
            &path,
            concat!(
                "read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "read2\tCELL1\tUMI1\n",
                "read3\tCELL1\tUMI2\n",
                "read4\tCELL2\tUMI3\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();
        let counts = table.pair_counts();

        assert_eq!(
            counts.get(&(b"CELL1".to_vec(), b"UMI1".to_vec())),
            Some(&2)
        );
        assert_eq!(
            counts.get(&(b"CELL1".to_vec(), b"UMI2".to_vec())),
            Some(&1)
        );
        assert_eq!(
            counts.get(&(b"CELL2".to_vec(), b"UMI3".to_vec())),
            Some(&1)
        );

        fs::remove_file(path).ok();
    }

    #[test]
    fn summarize_pairs_applies_pair_and_cell_thresholds() {
        let path = tmp_path("summary.tsv");

        write_text(
            &path,
            concat!(
                "read_id\traw_cb\traw_umi\n",
                "r1\tCELL1\tUMI1\n",
                "r2\tCELL1\tUMI1\n",
                "r3\tCELL1\tUMI2\n",
                "r4\tCELL1\tUMI2\n",
                "r5\tCELL2\tUMI3\n",
                "r6\tCELL2\tUMI3\n",
                "r7\tCELL3\tUMI4\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        let stats = table.summarize_pairs(2, 2);

        assert_eq!(stats.cell_entries, 1);
        assert_eq!(stats.unique_cell_umi_combos, 2);
        assert_eq!(stats.total_pair_observations, 4);

        assert_eq!(stats.umis_per_cell.n, 1);
        assert_eq!(stats.umis_per_cell.min, 2);
        assert_eq!(stats.umis_per_cell.max, 2);

        assert_eq!(stats.detections_per_cell_umi.n, 2);
        assert_eq!(stats.detections_per_cell_umi.min, 2);
        assert_eq!(stats.detections_per_cell_umi.max, 2);

        fs::remove_file(path).ok();
    }

    #[test]
    fn display_for_record_is_informative() {
        let rec = ReadTagRecord::new(
            "read1".to_string(),
            Some("orig1".to_string()),
            b"CELL1".to_vec(),
            Vec::new(),
            b"UMI1".to_vec(),
            vec![41, 41, 41, 41],
        );

        let text = rec.to_string();

        assert!(text.contains("read_id=read1"));
        assert!(text.contains("original_read_id=orig1"));
        assert!(text.contains("cell=CELL1"));
        assert!(text.contains("cell_qual=-"));
        assert!(text.contains("umi=UMI1"));
        assert!(text.contains("umi_qual=JJJJ"));
    }

    #[test]
    fn display_for_table_shows_count_and_at_most_three_records() {
        let path = tmp_path("display.tsv");

        write_text(
            &path,
            concat!(
                "read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "read2\tCELL2\tUMI2\n",
                "read3\tCELL3\tUMI3\n",
                "read4\tCELL4\tUMI4\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();
        let text = table.to_string();

        assert!(text.contains("ReadTagTable: 4 records"));

        let shown_records = text
            .lines()
            .filter(|line| line.trim_start().starts_with('['))
            .count();

        assert_eq!(shown_records, 3);

        fs::remove_file(path).ok();
    }

    #[test]
    fn writer_writes_header_and_rows() {
        let mut out = Vec::new();

        {
            let mut writer = ReadTagTableWriter::new(&mut out).unwrap();

            writer
                .write_record(&ReadTagRecord::new(
                    "read1".to_string(),
                    Some("orig1".to_string()),
                    b"CELL1".to_vec(),
                    vec![40, 40, 40, 40],
                    b"UMI1".to_vec(),
                    vec![41, 41, 41, 41],
                ))
                .unwrap();

            writer.flush().unwrap();
        }

        let text = String::from_utf8(out).unwrap();

        assert!(text.starts_with(&READ_TAG_TABLE_COLUMNS.join("\t")));
        assert!(text.contains(
            "read1\torig1\t\tCELL1\tIIII\tUMI1\tJJJJ\tok"
        ));
    }
}