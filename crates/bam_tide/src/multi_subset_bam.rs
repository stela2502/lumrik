use anyhow::{Context, Result};
use rust_htslib::bam::record::{Aux, Record};

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

pub struct Subsetter {
    tags: BTreeMap<String, usize>,
    pub ofile_names: Vec<PathBuf>,
    group_names: Vec<String>,
}

impl Default for Subsetter {
    fn default() -> Self {
        Self::new()
    }
}

impl Subsetter {
    pub fn new() -> Self {
        Self {
            tags: BTreeMap::new(),
            ofile_names: Vec::new(),
            group_names: Vec::new(),
        }
    }

    /// Read one value-per-line tag list.
    ///
    /// Each list file creates one output BAM.
    pub fn read_simple_list(&mut self, bc_file: &Path, prefix: &str) -> Result<()> {
        let writer_id = self.ofile_names.len();

        let file = File::open(bc_file)
            .with_context(|| format!("opening tag list {}", bc_file.display()))?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let tag_value =
                line.with_context(|| format!("reading line from {}", bc_file.display()))?;
            let tag_value = tag_value.trim();

            if tag_value.is_empty() || tag_value.starts_with('#') {
                continue;
            }

            self.tags.insert(tag_value.to_string(), writer_id);
        }

        let stem = bc_file
            .file_stem()
            .and_then(|x| x.to_str())
            .with_context(|| format!("could not derive output name from {}", bc_file.display()))?;

        self.group_names.push(stem.to_string());

        self.ofile_names
            .push(PathBuf::from(format!("{prefix}{stem}.bam")));

        Ok(())
    }

    pub fn group_name(&self, id: usize) -> &str {
        &self.group_names[id]
    }

    pub fn process_record_with_name<'a>(
        &'a self,
        record: &Record,
        tag: &[u8; 2],
    ) -> Option<(usize, &'a str)> {
        let id = self.process_record(record, tag)?;
        Some((id, &self.group_names[id]))
    }

    pub fn process_record(&self, record: &Record, tag: &[u8; 2]) -> Option<usize> {
        let tag_value = get_tag_value(record, tag)?;
        self.tags.get(&tag_value).copied()
    }
}

fn get_tag_value(record: &Record, tag: &[u8; 2]) -> Option<String> {
    match record.aux(tag).ok()? {
        Aux::String(s) => Some(s.to_string()),
        Aux::Char(c) => Some((c as char).to_string()),
        Aux::I8(v) => Some(v.to_string()),
        Aux::U8(v) => Some(v.to_string()),
        Aux::I16(v) => Some(v.to_string()),
        Aux::U16(v) => Some(v.to_string()),
        Aux::I32(v) => Some(v.to_string()),
        Aux::U32(v) => Some(v.to_string()),
        Aux::Float(v) => Some(v.to_string()),
        _ => None,
    }
}
