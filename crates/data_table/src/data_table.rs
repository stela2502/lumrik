use crate::{Factor, factor::FactorJson};
use anyhow::Result;
use csv::WriterBuilder;
use ndarray::{Array2, Axis, concatenate, s};
use rand::seq::SliceRandom;
use rand::rng;
use serde_json;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct DataTable {
    pub headers: Vec<String>,
    pub numeric_data: Array2<f64>,
    pub factors: HashMap<String, Factor>,
    pub exclude: HashSet<String>,
    header_lookup: HashMap<String, usize>,
    next_order: f64,
    //pub max_levels: usize,
}

impl Default for DataTable {
    fn default() -> Self {
        DataTable {
            headers: Vec::new(),
            numeric_data: Array2::zeros((0, 0)),
            factors: HashMap::new(),
            exclude: HashSet::new(),
            header_lookup: HashMap::new(),
            next_order: 0.0,
        }
    }
}

impl fmt::Display for DataTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "DataTable Summary:")?;
        writeln!(
            f,
            "Rows: {}, Columns: {}",
            self.numeric_data.nrows(),
            self.numeric_data.ncols()
        )?;
        writeln!(f, "Headers:")?;
        for (idx, header) in self.headers.iter().enumerate() {
            writeln!(f, "  {}: {}", idx, header)?;
        }
        writeln!(f, "Factors:")?;
        for factor in self.factors.values() {
            let _ = writeln!(f, "{}", factor);
        }
        writeln!(f, "Excluded columns: {:?}", self.exclude)
    }
}

impl DataTable {
    /// ------------------------------------------------------------------------
    /// Core loader: reads TSV/CSV file, applies factors if available.
    /// If the `.factors.json` exists, it will be used.
    /// If not, it loads data heuristically and continues silently.
    /// ------------------------------------------------------------------------
    pub fn from_tsv<P: AsRef<Path> + std::fmt::Debug, FF: AsRef<Path> + std::fmt::Debug>(
        file_path: P,
        delimiter: u8,
        categorical_cols: HashSet<String>,
        factors_file: FF,
    ) -> Result<Self> {
        let mut ret = DataTable::default();

        // --- 1. Load existing factors if the file exists ---
        let factors_path = factors_file.as_ref();
        if factors_path.exists() {
            println!("🧬 Using existing factors file: {:?}", factors_path);
            ret.load_factors(factors_path)?;
        } else {
            println!("ℹ️ No factors file found, proceeding without it (new dataset?)");
        }

        // --- 2. Open file and prepare CSV reader ---
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .from_path(&file_path)?;

        // --- 3. Expand headers based on existing factors (for one-hot) ---
        let headers: Vec<String> = rdr
            .headers()?
            .iter()
            .flat_map(|s| {
                // 💡 handle R-style empty header (rownames column)
                let _ = if s.trim().is_empty() {
                    println!("⚠️ Detected unlabeled first(?) column -> treating as 'barcode'");
                    "barcode".to_string()
                } else {
                    s.to_string()
                };
                if let Some(fact) = ret.factors.get(s) {
                    if fact.one_hot {
                        let mut c = vec![s.to_string()];
                        c.extend(fact.all_column_names());
                        c
                    } else {
                        vec![s.to_string()]
                    }
                } else {
                    vec![s.to_string()]
                }
            })
            .collect();

        ret.headers = headers.clone();
        ret.header_lookup = headers
            .iter()
            .enumerate()
            .map(|(id, name)| (name.clone(), id))
            .collect();

        // --- 4. Mark forced categorical columns as factors ---
        for header in &headers {
            if categorical_cols.contains(header) {
                println!("Forcing header {header} to be a factor");
                ret.factors
                    .entry(header.clone())
                    .or_insert_with(|| Factor::new(header, false));
            }
        }

        // --- 5. Parse rows (one-hot aware logic) ---
        let n_cols = headers.len();
        let mut raw_rows: Vec<Vec<f64>> = Vec::new();

        for result in rdr.records() {
            let record = result?;
            let mut row: Vec<f64> = Vec::with_capacity(n_cols);
            let mut expanded = 0;

            for (i, value) in record.iter().enumerate() {
                let trimmed = value.trim().trim_matches('"');
                match trimmed.replace(',', ".").parse::<f64>() {
                    Ok(num) => {
                        if let Some(factor) = ret.factors.get_mut(&headers[i + expanded]) {
                            let (idx, col_to_add, all_cols) = factor.push(&num.to_string());
                            let alt = if idx.is_nan() { f64::NAN } else { 0.0 };
                            match all_cols {
                                Some(cols) => {
                                    expanded += cols.len();
                                    row.push(factor.get_f64(trimmed));
                                    for cname in cols.iter().cloned() {
                                        row.push(if cname == col_to_add { idx } else { alt });
                                    }
                                }
                                None => row.push(idx),
                            }
                        } else {
                            row.push(num);
                        }
                    }
                    Err(_) => {
                        if (trimmed.is_empty() || trimmed.eq_ignore_ascii_case("NA"))
                            && !ret.factors.contains_key(&headers[i + expanded]) {
                                row.push(f64::NAN);
                                continue;
                            }
                        let factor = ret
                            .factors
                            .entry(headers[i + expanded].clone())
                            .or_insert_with(|| Factor::new(&headers[i + expanded], false));
                        let (idx, col_to_add, all_cols) = factor.push(trimmed);
                        let alt = if idx.is_nan() { f64::NAN } else { 0.0 };

                        match all_cols {
                            Some(cols) => {
                                expanded += cols.len();
                                row.push(factor.get_f64(trimmed));
                                for cname in cols.iter() {
                                    row.push(if *cname == col_to_add { idx } else { alt });
                                }
                            }
                            None => row.push(if trimmed.is_empty() { f64::NAN } else { idx }),
                        }
                    }
                }
            }

            raw_rows.push(row);
        }

        // --- 6. Build Array2 ---
        let n_rows = raw_rows.len();
        let mut numeric_data = Array2::<f64>::zeros((n_rows, n_cols));
        for (i, row) in raw_rows.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                numeric_data[[i, j]] = val;
            }
        }

        ret.numeric_data = numeric_data;
        Ok(ret)
    }

    /// High-level loader: handles factor file and performs cell-name consistency checks.
    /// Reuses `from_tsv()` internally.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug, FF: AsRef<Path> + std::fmt::Debug>(
        file_path: P,
        delimiter: u8,
        categorical_cols: HashSet<String>,
        factors_file: FF,
    ) -> Result<Self> {
        let file_path_ref = file_path.as_ref();
        let factors_file_ref = factors_file.as_ref();

        // --- 1. Load the numeric/factor data from TSV ---
        let mut ret = Self::from_tsv(file_path_ref, delimiter, categorical_cols, &factors_file)?;

        // --- 2. Load or create factors file ---
        if factors_file_ref.exists() {
            println!("Factors are loaded from file");
            ret.load_factors(factors_file_ref)?;
        } else {
            println!(
                "Saved a new factors file to fine tune the factors: '{:?}'",
                &factors_file_ref
            );
            ret.save_factors(factors_file_ref)?;
            panic!("
Please review and update the factors file so that it accurately reflects the logic in the data.

The factors file is a JSON-formatted file, for example:

[
  {{
    'column': 'status2',
    'levels': [
      '1',
      '0'
    ],
    'numeric': [
      0.0,
      1.0
    ],
    'matching': null,
    'one_hot': false
  }}
]
In this example, there is an error: the numeric values do not match the actual data. They should be [1.0, 0.0]. Once corrected, the factor will work as expected.


The one_hot option allows the factor to be expanded into multiple 0.0/1.0 columns—two in this case. This is particularly useful when the factor levels have no inherent numeric order or relationship:

[
  {{
    'column': 'tp53_mutation_type',
    'levels': [
      'Missense',
      'Nonsense',
      'Frameshift',
      'Splice_site',
      'Silent',
    ],
    'numeric': [
      0.0,
      1.0,
      2.0,
      3.0,
      4.0,
    ],
    'matching': null,
    'one_hot': true
  }}
]

This factor is categorical with no inherent order, so it’s a good candidate for one-hot encoding in a model. 
Each level will be represented as a separate binary column (0/1) if one-hot encoding is used.

                ");
        }

        Ok(ret)
    }

    /// Add a new data row (group) initialized with NaN values across all columns (cells).
    ///
    /// Each row represents a dataset or feature.
    /// If this is the first row added, you must specify the expected number of columns via `size`.
    /// If `is_factor` is true, a corresponding `Factor` is also created and registered.
    ///
    /// # Panics
    /// Panics if the data table is empty and no `size` is provided.
    ///
    /// # Example
    /// ```
    /// use data_table::DataTable;
    /// let mut data = DataTable::default();
    /// data.add_dataset("age", false, Some(10));     // first row → needs size
    /// data.add_dataset("treatment", true, None);    // second row → inferred columns
    /// ```
    pub fn add_dataset(&mut self, dataset_name: &str, is_factor: bool, size: Option<usize>) {
        // Don't add the same dataset twice
        if self.headers.contains(&dataset_name.to_string()) {
            #[cfg(debug_assertions)]
            println!("Row '{}' already exists, skipping.", dataset_name);
            return;
        }

        let n_rows = self.numeric_data.nrows();
        let n_cols = self.numeric_data.ncols();

        // Case 1: Empty table — initialize with one row and known number of columns
        if n_rows == 0 {
            let init_cols = size.unwrap_or_else(|| {
                panic!("🛑 Must provide 'size' (number of cells) when adding the first dataset!");
            });
            self.numeric_data = Array2::<f64>::from_elem((1, init_cols), f64::NAN);
        } else {
            // Case 2: Append a new row of NaN across existing columns
            let new_row = Array2::<f64>::from_elem((1, n_cols), f64::NAN);
            self.numeric_data = concatenate(Axis(0), &[self.numeric_data.view(), new_row.view()])
                .expect("❌ Failed to concatenate new row");
        }

        // Register dataset name (row name)
        self.headers.push(dataset_name.to_string());

        // Register mapping for fast lookup
        let row_idx = self.numeric_data.nrows() - 1;

        self.header_lookup.insert(dataset_name.to_string(), row_idx);

        // Optionally register a Factor
        if is_factor {
            self.factors
                .insert(dataset_name.to_string(), Factor::new(dataset_name, false));
        }

        #[cfg(debug_assertions)]
        println!(
            "✅ Added dataset row '{}' (factor={}) → total rows: {}, columns: {}",
            dataset_name,
            is_factor,
            self.numeric_data.nrows(),
            self.numeric_data.ncols()
        );
    }

    /// Update a single cell by dataset (row) name and sample (column) index.
    /// Returns true if successful, false if the dataset or index is invalid.
    pub fn update_value(&mut self, dataset_name: &str, col: usize, value: f64) -> bool {
        // Look up the row index for this dataset
        let row_idx = match self.header_lookup.get(dataset_name) {
            Some(idx) => idx,
            None => {
                eprintln!(
                    "⚠️ update_value(): dataset '{}' not found. Headers: {:?}",
                    dataset_name, self.headers
                );
                return false;
            }
        };

        let (n_rows, n_cols) = self.numeric_data.dim();
        if *row_idx >= n_rows {
            eprintln!(
                "🛑 Row index {} out of bounds (rows={}) for dataset '{}'",
                row_idx, n_rows, dataset_name
            );
            return false;
        }
        if col >= n_cols {
            eprintln!(
                "🛑 Column index {} out of bounds (cols={}) for dataset '{}'",
                col, n_cols, dataset_name
            );
            return false;
        }

        if self.numeric_data[[*row_idx, col]] == value {
            false
        } else {
            self.numeric_data[[*row_idx, col]] = value;
            true
        }
    }

    /// Update a single cell by dataset (row) name and sample (column) index.
    /// The value will represent the order in which this function is called.
    /// reset the order with self.reset_order().
    pub fn update_order(&mut self, dataset_name: &str, col: usize) -> bool {
        // Look up the row index for this dataset
        let row_idx = match self.header_lookup.get(dataset_name) {
            Some(idx) => idx,
            None => {
                eprintln!(
                    "⚠️ update_value(): dataset '{}' not found. Headers: {:?}",
                    dataset_name, self.headers
                );
                return false;
            }
        };
        let (n_rows, n_cols) = self.numeric_data.dim();
        if *row_idx >= n_rows {
            eprintln!(
                "🛑 Row index {} out of bounds (rows={}) for dataset '{}'",
                row_idx, n_rows, dataset_name
            );
            return false;
        }
        if col >= n_cols {
            eprintln!(
                "🛑 Column index {} out of bounds (cols={}) for dataset '{}'",
                col, n_cols, dataset_name
            );
            return false;
        }

        self.numeric_data[[*row_idx, col]] = self.next_order;
        self.next_order += 1.0;
        true
    }
    pub fn reset_order(&mut self) {
        self.next_order = 0.0;
    }

    /// Update a single cell by dataset name and row index.
    /// Returns true if successful, false if the column or index is invalid.
    pub fn update_value_str(&mut self, dataset_name: &str, col: usize, value: &str) -> bool {
        if let Some(fact) = self.factors.get_mut(dataset_name) {
            let val = fact.level_to_index(value);
            self.update_value(dataset_name, col, val)
        } else {
            // Numeric column → try to parse directly
            match value.trim().parse::<f64>() {
                Ok(v) => self.update_value(dataset_name, col, v),
                Err(_) => {
                    eprintln!(
                        "Value '{}' not recognized for column '{}'",
                        value, dataset_name
                    );
                    false
                }
            }
        }
    }

    /// Split into train and test by a fraction (e.g. 0.7 = 70% train, 30% test).
    pub fn train_test_split(&self, train_fraction: f64) -> (DataTable, DataTable) {
        assert!(train_fraction > 0.0 && train_fraction < 1.0);

        let n_rows = self.numeric_data.nrows();
        let mut indices: Vec<usize> = (0..n_rows).collect();
        let mut rng = rng();
        indices.shuffle(&mut rng);

        let train_size = (n_rows as f64 * train_fraction).round() as usize;

        let train_idx = &indices[..train_size];
        let test_idx = &indices[train_size..];

        let train_data = self.numeric_data.select(Axis(0), train_idx);
        let test_data = self.numeric_data.select(Axis(0), test_idx);

        let make_subset = |data: Array2<f64>| DataTable {
            headers: self.headers.clone(),
            numeric_data: data,
            factors: self.factors.clone(),
            exclude: self.exclude.clone(),
            header_lookup: self.header_lookup.clone(),
            ..Default::default()
        };

        (make_subset(train_data), make_subset(test_data))
    }


    pub fn data_summary(&self) {
        println!(
            "Shape: {} rows x {} columns",
            self.numeric_data.nrows(),
            self.numeric_data.ncols()
        );

        // 2️⃣ Check first few rows
        println!("First 5 rows:");
        for i in 0..self.numeric_data.nrows().min(5) {
            println!("{:?}", self.numeric_data.row(i).to_vec());
        }

        // 3️⃣ Column-wise summaries
        for j in 0..self.numeric_data.ncols().min(10) {
            let col = self.numeric_data.column(j);
            let n_total = col.len();
            let n_na = col.iter().filter(|v| v.is_nan()).count();
            let mean = col.iter().filter(|v| !v.is_nan()).sum::<f64>() / (n_total - n_na) as f64;
            println!(
                "Col {}: NA fraction = {:.2}, mean of non-NA = {:.2}",
                j,
                n_na as f64 / n_total as f64,
                mean
            );
        }
    }



    /// Remove all rows that contain any NaN in numeric_data
    pub fn filter_all_na_rows(&mut self, usable: &[String] ) {
        let n_rows = self.numeric_data.nrows();
        let n_cols = self.numeric_data.ncols();
        println!(
            "filter_all_na_rows got {} rows and {} columns and checks {} of these columns for na's",
            n_rows,
            n_cols,
            usable.len()
        );

        // Build a lookup of usable column indices
        let usable_indices: Vec<usize> = usable
            .iter()
            .filter_map(|col| self.headers.iter().position(|h| h == col))
            .collect();
        if usable_indices.len() == n_rows {
            println!("No rows containing NA values found");
            return;
        }
        // Determine which rows to keep
        let keep_rows: Vec<usize> = (0..n_rows)
            .filter(|&i| {
                !usable_indices
                    .iter()
                    .any(|&j| self.numeric_data[[i, j]].is_nan())
            })
            .collect();
        if keep_rows.len() == n_rows {
            println!("No na's found in the {} columns", usable.len());
            return;
        }

        // Rebuild numeric_data with only the kept rows
        let mut filtered = Array2::<f64>::zeros((keep_rows.len(), n_cols));
        for (new_i, &old_i) in keep_rows.iter().enumerate() {
            filtered
                .row_mut(new_i)
                .assign(&self.numeric_data.slice(s![old_i, ..]));
        }
        self.numeric_data = filtered;

        // Filter factors: keep only entries for kept rows
        let mut new_factors = HashMap::<String, Factor>::new();

        for factor in self.factors.values() {
            if let Some(col_id) = self.headers.iter().position(|h| h == &factor.column_name) {
                new_factors.insert(
                    factor.column_name.clone(),
                    factor.subset(&self.numeric_data, col_id),
                );
            } else {
                panic!(
                    "a previousely known column has vanished?! {}",
                    factor.column_name
                );
            }
        }
        self.factors = new_factors;

        println!(
            "Filtered out {} rows containing NaNs. Remaining rows: {}",
            n_rows - keep_rows.len(),
            keep_rows.len()
        );
    }
    #[allow(clippy::doc_overindented_list_items)]
    /// Impute missing values in `numeric_data` using K-nearest neighbours.
    ///
    /// * `k` – number of neighbours to use (e.g. 3 for “mean of 3 closest”).
    /// * `min_common` – minimum number of shared non-NA features required
    ///                  between two rows to consider them neighbours.
    /// * `weighted` – if `true`, use distance-weighted mean;
    ///                if `false`, simple mean of the k neighbours.
    ///
    /// Missing values are represented as `f64::NAN`.
    pub fn impute_knn(&mut self, k: usize, min_common: usize, weighted: bool) {
        let eps = 1e-8_f64;
        let (n_rows, n_cols) = (self.numeric_data.nrows(), self.numeric_data.ncols());

        // ---------- z-score scale copy for distance calculations ----------
        let scaled = self.numeric_data.clone();
        let mut col_means = vec![0.0; n_cols];
        let mut col_stds = vec![1.0; n_cols];

        for j in 0..n_cols {
            let col = self.numeric_data.column(j);
            let vals: Vec<f64> = col.iter().copied().filter(|v| !v.is_nan()).collect();
            if vals.is_empty() {
                continue;
            }
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / vals.len() as f64;
            let std = var.sqrt().max(eps);
            col_means[j] = mean;
            col_stds[j] = std;
            #[allow(clippy::needless_range_loop)]
            for j in 0..n_cols {
                if self.factors.contains_key(&self.headers[j]) {
                    // factor column → use mode of non-NA values
                    let mut counts = std::collections::HashMap::new();
                    for &v in &vals {
                        *counts.entry(v as usize).or_insert(0) += 1;
                    }
                    let &mode = counts.iter().max_by_key(|(_, c)| *c).unwrap().0;
                    col_means[j] = mode as f64; // for later distance scaling, can still normalize if needed
                } else {
                    // numeric → mean
                    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
                    col_means[j] = mean;
                }
            }
        }

        // helper to compute Euclidean distance using only shared non-NA cols
        let distance = |a: usize, b: usize| -> Option<f64> {
            let mut sum = 0.0;
            let mut cnt = 0;
            for j in 0..n_cols {
                let va = scaled[[a, j]];
                let vb = scaled[[b, j]];
                if va.is_nan() || vb.is_nan() {
                    continue;
                }
                let d = va - vb;
                sum += d * d;
                cnt += 1;
            }
            if cnt >= min_common {
                Some((sum / cnt as f64).sqrt())
            } else {
                None
            }
        };

        // ---------- main loop over rows with NAs ----------
        for i in 0..n_rows {
            // which columns are missing for this row?
            let missing: Vec<_> = (0..n_cols)
                .filter(|&j| self.numeric_data[[i, j]].is_nan())
                .collect();
            if missing.is_empty() {
                continue;
            }

            // collect neighbours and their distances
            let mut neigh = Vec::new();
            for j in 0..n_rows {
                if j == i {
                    continue;
                }
                if let Some(d) = distance(i, j) {
                    neigh.push((j, d));
                }
            }
            neigh.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            if neigh.is_empty() {
                continue;
            }

            // impute each missing column
            for &col in &missing {
                let mut vals = Vec::new();
                for &(row_idx, dist) in &neigh {
                    let v = self.numeric_data[[row_idx, col]];
                    if !v.is_nan() {
                        vals.push((v, dist));
                        if vals.len() >= k {
                            break;
                        }
                    }
                }
                if vals.is_empty() {
                    continue;
                }

                let imputed = if self.factors.contains_key(&self.headers[col]) {
                    // Factor column → mode of neighbors
                    let mut counts = std::collections::HashMap::new();
                    for &(v, _) in &vals {
                        let level = v.round() as usize; // round to nearest valid factor index
                        *counts.entry(level).or_insert(0) += 1;
                    }
                    *counts.iter().max_by_key(|(_, c)| *c).unwrap().0 as f64
                } else {
                    // Numeric → weighted or regular mean
                    if weighted {
                        let mut num = 0.0;
                        let mut den = 0.0;
                        for (v, d) in vals {
                            let w = 1.0 / (d + eps);
                            num += w * v;
                            den += w;
                        }
                        num / den
                    } else {
                        vals.iter().map(|(v, _)| *v).sum::<f64>() / vals.len() as f64
                    }
                };
                self.numeric_data[[i, col]] = imputed;
            }
        }
    }

    /// Remove any rows where `column` has NaN in the numeric data
    pub fn filter_na(&mut self, column: &str) {
        // Find the column index
        let col_idx = self
            .headers
            .iter()
            .position(|h| h == column)
            .expect("Column not found");

        // Determine which rows to keep
        let keep_rows: Vec<usize> = self
            .numeric_data
            .column(col_idx)
            .indexed_iter()
            .filter_map(|(i, &val)| if !val.is_nan() { Some(i) } else { None })
            .collect();

        // Rebuild numeric_data with only the kept rows
        let mut filtered = Array2::<f64>::zeros((keep_rows.len(), self.numeric_data.ncols()));
        for (new_i, &old_i) in keep_rows.iter().enumerate() {
            filtered
                .row_mut(new_i)
                .assign(&self.numeric_data.slice(s![old_i, ..]));
        }
        self.numeric_data = filtered;
    }
    /*
    /// Remove columns with variance below `threshold`
    pub fn filter_low_var(&mut self, threshold: f64) -> usize {
        return 0;
        let mut keep_cols = Vec::new();
        let mut filtered = 0;
        for j in 0..self.numeric_data.ncols() {
            let col = self.numeric_data.column(j);
            let vals: Vec<f64> = col.iter().copied().filter(|v| !v.is_nan()).collect();
            if vals.is_empty() {
                filtered += 1;
                println!("Dropping empty column: {}", self.headers[j]);
                continue;
            }

            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / vals.len() as f64;

            if var > threshold {
                keep_cols.push(j);
            } else {
                filtered += 1;
                println!(
                    "Dropping low-variance column (mean {mean:2e}, var {var:2e}): {}",
                    self.headers[j]
                );
            }
        }

        // Rebuild numeric_data and headers
        let mut filtered_data = Array2::<f64>::zeros((self.numeric_data.nrows(), keep_cols.len()));
        let mut new_headers = Vec::with_capacity(keep_cols.len());

        for (new_j, &old_j) in keep_cols.iter().enumerate() {
            filtered_data
                .column_mut(new_j)
                .assign(&self.numeric_data.column(old_j));
            new_headers.push(self.headers[old_j].clone());
        }

        self.numeric_data = filtered_data;
        self.headers = new_headers;

        // Remove factors whose columns were dropped
        self.factors.retain(|name, _| self.headers.contains(name));
        filtered
    }
    */

    /// Compute fraction of NAs per feature column
    pub fn feature_na_fraction(&self) -> HashMap<String, f64> {
        let mut na_frac: HashMap<String, f64> = HashMap::new();
        for (idx, name) in self.headers.clone().into_iter().enumerate() {
            let col = self.numeric_data.column(idx);
            let n_total = col.len();
            let n_na = col.iter().filter(|v| v.is_nan()).count();
            na_frac.insert(name.clone(), n_na as f64 / n_total as f64);
        }
        na_frac
    }

    /// Return features filtered by max allowed NA fraction
    pub fn filter_features_by_na(&self, max_na_frac: f64) -> Vec<String> {
        self.feature_na_fraction()
            .into_iter()
            .filter(|(_, frac)| *frac <= max_na_frac)
            .map(|(name, _)| name)
            .collect()
    }
    /*
    /// Return filtered numeric data and corresponding indices for allowed features
    pub fn as_array2_filtered(&self, max_na_frac: f64) -> (Array2<f64>, Vec<usize>) {
        let allowed_features = self.filter_features_by_na( max_na_frac);
        let mut indices = Vec::new();
        for name in &allowed_features {
            if let Some(idx) = self.headers.iter().position(|h| h == name) {
                indices.push(idx);
            }
        }

        let n_rows = self.numeric_data.nrows();
        let n_cols = indices.len();
        let mut arr = Array2::<f64>::zeros((n_rows, n_cols));

        for (j, &col_idx) in indices.iter().enumerate() {
            arr.column_mut(j).assign(&self.numeric_data.column(col_idx));
        }

        (arr, indices)
    }*/

    /// Return numeric data as ndarray, optionally selecting columns
    #[allow(dead_code)]
    pub fn as_array2(&self, columns: Option<&[String]>) -> Array2<f64> {
        match columns {
            Some(cols) => {
                let indices: Vec<usize> = cols
                    .iter()
                    .map(|c| {
                        self.headers
                            .iter()
                            .position(|h| h == c)
                            .expect("Column not found")
                    })
                    .collect();
                let mut arr = Array2::<f64>::zeros((self.numeric_data.nrows(), indices.len()));
                for (j, &col_idx) in indices.iter().enumerate() {
                    arr.column_mut(j).assign(&self.numeric_data.column(col_idx));
                }
                arr
            }
            None => self.numeric_data.clone(),
        }
    }

    /// Return a single column as Vec<f64>
    pub fn as_vec_f64(&self, column: &str) -> Vec<f64> {
        let idx = self
            .headers
            .iter()
            .position(|h| h == column)
            .unwrap_or_else(|| {
                panic!(
                    "Column '{}' not found in the dataset; all columns: \n{}",
                    column,
                    self.headers(20).join("\n")
                )
            });
        self.numeric_data.column(idx).to_vec()
    }

    /// Return a single column as Vec<u8>
    pub fn as_vec_u8(&self, column: &str) -> Vec<u8> {
        let idx = self.headers.iter().position(|h| h == column)
            .unwrap_or_else(|| 
                panic!("Column '{}' not found in the dataset; all columns: \n{}\nColumn '{}' not found in the dataset", 
                    column, self.headers(20).join("\n") ,column ));
        self.numeric_data
            .column(idx)
            .iter()
            .map(|&v| v as u8)
            .collect()
    }

    /// Return a single column as Option<Vec<String>> - if it is a factor
    pub fn as_vec_string(&self, column: &str) -> Option<Vec<String>> {
        let idx = self
            .headers
            .iter()
            .position(|h| h == column)
            .expect("Column not found");
        if self.factors.contains_key(column) {
            let fact = self.factors.get(column).unwrap();
            //println!("I will use the factor {fact:?} to translate the columns ids like {:?}",  self.numeric_data.column(idx).iter().take(10).collect::<Vec<_>>() );
            Some(
                self.numeric_data
                    .column(idx)
                    .iter()
                    .map(|&v| fact.get_string(v))
                    .collect(),
            )
        } else {
            //println!("Column '{column}' is no Factor here\n{:?}\nFactors I have: \n{:?}\n", self.headers, self.factors.keys() );
            None
        }
    }

    pub fn headers(&self, max: usize) -> Vec<&str> {
        self.headers
            .iter() // borrow each String
            .take(self.headers.len().min(max)) // take at most `max`
            .map(|s| s.as_str()) // convert &String → &str
            .collect() // collect into Vec<&str>
    }

    /// Write data (numeric + factor) to CSV
    pub fn to_file<P: AsRef<std::path::Path>>(&self, file_path: P, delimiter: u8) -> Result<()> {
        let mut wtr = WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(file_path)?;

        // Write header
        wtr.write_record(&self.headers)?;

        for i in 0..self.numeric_data.nrows() {
            let mut record: Vec<String> = Vec::with_capacity(self.headers.len());

            for (j, _col_name) in self.headers.iter().enumerate() {
                let val = self.numeric_data[[i, j]];
                if val.is_nan() {
                    record.push("NA".to_string());
                } else {
                    record.push(val.to_string());
                }
            }

            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        Ok(())
    }

    /// Save all factors to a JSON file
    pub fn save_factors<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let defs: Vec<FactorJson> = self
            .factors
            .iter()
            .map(|(col, factor)| factor.as_json(col))
            .collect();
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &defs)?;
        Ok(())
    }

    /// Load factors from a JSON file
    pub fn load_factors<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let defs: Vec<FactorJson> = serde_json::from_reader(reader)?;
        for def in defs {
            self.factors
                .insert(def.column.clone(), Factor::from_def(&def));
        }
        Ok(())
    }
}



#[cfg(test)]
mod tests{
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;

    use crate::{DataTable, Factor};

    #[test]
    fn test_survivaldata_one_hot_factors() -> Result<(), Box<dyn std::error::Error>> {
        // --- Prepare temporary directory ---
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let csv_path = dir.as_path().join("test_factors.csv");
        let factors_path = dir.as_path().join("test_factors_factors.json");

        // --- Write CSV file ---
        // Two categorical columns: Color, Number
        // Values: Color = Red Blue Blue Green Blue Red Na Na
        //         Number = 1 2 1 3 2 1 Na Na
        let csv_content = "\
Color,Number
Red,1
Blue,2
Blue,1
Green,3
Blue,2
Red,1
Na,Na
Na,Na
";
        let mut file = File::create(&csv_path)?;
        file.write_all(csv_content.as_bytes())?;

        let mut color = Factor::new("Color", false);
        let _ = color.push("Red");
        let _ = color.push("Blue");
        let _ = color.push("Green");
        color.one_hot = true;

        let mut number = Factor::new("Number", false);
        let _ = number.push("1");
        let _ = number.push("2");
        let _ = number.push("3");
        number.one_hot = true;
        let mut temp = DataTable::default();
        temp.factors.insert("Color".to_string(), color);
        temp.factors.insert("Number".to_string(), number);

        let _ = temp.save_factors(&factors_path); //default location
        let categorical_cols = HashSet::<String>::new();
        // --- Load DataTable ---
        let data = DataTable::from_file(&csv_path, b',', categorical_cols, &factors_path)?;

        // --- Assertions ---

        // 1️⃣ Check original factor columns contain indices / NaN
        println!("The data object: {data}");
        let color_idx_col = data.headers.iter().position(|h| h == "Color").unwrap();
        let number_idx_col = data.headers.iter().position(|h| h == "Number").unwrap();
        assert_eq!(
            color_idx_col, 0,
            "expected color_idx_col 0 - have these cols: {:?}",
            &data.headers
        );
        assert_eq!(
            number_idx_col, 4,
            "expected number_idx_col 4 - have these cols: {:?}",
            &data.headers
        );

        assert_eq!(data.numeric_data[[0, color_idx_col]], 0.0); // Red first index 0
        assert_eq!(data.numeric_data[[1, color_idx_col]], 1.0); // Blue index 1
        assert!(data.numeric_data[[6, color_idx_col]].is_nan()); // Na

        assert_eq!(data.numeric_data[[0, number_idx_col]], 0.0); // 1 -> index 0
        assert_eq!(data.numeric_data[[1, number_idx_col]], 1.0); // 2 -> index 1
        assert!(data.numeric_data[[6, number_idx_col]].is_nan()); // Na

        assert_eq!(data.headers[1], "Color_Red");

        // Find the indices of the one-hot columns
        let color_red_idx = data.headers.iter().position(|h| h == "Color_Red").unwrap();
        let color_blue_idx = data.headers.iter().position(|h| h == "Color_Blue").unwrap();
        let color_green_idx = data
            .headers
            .iter()
            .position(|h| h == "Color_Green")
            .unwrap();

        let number_1_idx = data.headers.iter().position(|h| h == "Number_1").unwrap();
        let number_2_idx = data.headers.iter().position(|h| h == "Number_2").unwrap();
        let number_3_idx = data.headers.iter().position(|h| h == "Number_3").unwrap();

        // --- Expected one-hot values per row ---
        // Rows: Color = Red, Blue, Blue, Green, Blue, Red, NA, NA
        //       Number = 1,2,1,3,2,1,NA,NA
        let expected_color = [
            (1.0, 0.0, 0.0),                // Red
            (0.0, 1.0, 0.0),                // Blue
            (0.0, 1.0, 0.0),                // Blue
            (0.0, 0.0, 1.0),                // Green
            (0.0, 1.0, 0.0),                // Blue
            (1.0, 0.0, 0.0),                // Red
            (f64::NAN, f64::NAN, f64::NAN), // NA
            (f64::NAN, f64::NAN, f64::NAN), // NA
        ];

        let expected_number = [
            (1.0, 0.0, 0.0),                // 1
            (0.0, 1.0, 0.0),                // 2
            (1.0, 0.0, 0.0),                // 1
            (0.0, 0.0, 1.0),                // 3
            (0.0, 1.0, 0.0),                // 2
            (1.0, 0.0, 0.0),                // 1
            (f64::NAN, f64::NAN, f64::NAN), // NA
            (f64::NAN, f64::NAN, f64::NAN), // NA
        ];

        for row in 0..data.numeric_data.nrows() {
            // Color
            let val_red = data.numeric_data[[row, color_red_idx]];
            let val_blue = data.numeric_data[[row, color_blue_idx]];
            let val_green = data.numeric_data[[row, color_green_idx]];

            let (exp_red, exp_blue, exp_green) = expected_color[row];
            if exp_red.is_nan() {
                assert!(val_red.is_nan());
                assert!(val_blue.is_nan());
                assert!(val_green.is_nan());
            } else {
                assert_eq!(val_red, exp_red);
                assert_eq!(val_blue, exp_blue);
                assert_eq!(val_green, exp_green);
            }

            // Number
            let val_1 = data.numeric_data[[row, number_1_idx]];
            let val_2 = data.numeric_data[[row, number_2_idx]];
            let val_3 = data.numeric_data[[row, number_3_idx]];

            let (exp_1, exp_2, exp_3) = expected_number[row];
            if exp_1.is_nan() {
                assert!(
                    val_1.is_nan(),
                    "expected NA value for Number_1 at row {row}; all data {}",
                    data.numeric_data
                );
                assert!(
                    val_2.is_nan(),
                    "expected NA value for Number_2 at row {row}; all data {}",
                    data.numeric_data
                );
                assert!(
                    val_3.is_nan(),
                    "expected NA value for Number_3 at row {row}; all data {}",
                    data.numeric_data
                );
            } else {
                assert_eq!(
                    val_1, exp_1,
                    "Number_1 mismatch at row {}/{}: {} != {} all data {}",
                    number_1_idx, row, val_1, exp_1, data.numeric_data
                );
                assert_eq!(
                    val_2, exp_2,
                    "Number_2 mismatch at row {}/{}: {} != {} all data {}",
                    number_2_idx, row, val_2, exp_2, data.numeric_data
                );
                assert_eq!(
                    val_3, exp_3,
                    "Number_3 mismatch at row {}/{}: {} != {} all data {}",
                    number_3_idx, row, val_3, exp_3, data.numeric_data
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_survivaldata_train_test_split_unique() -> Result<(), Box<dyn std::error::Error>> {
        use std::fs::File;
        use std::io::Write;
        use tempfile::tempdir;

        // --- Prepare temporary directory and files ---
        let dir = tempdir()?;
        let csv_path = dir.path().join("unique_data.csv");
        let factors_path = dir.path().join("unique_factors.json");

        // --- Generate CSV data with 10 rows, 10 columns, unique values 1..100 ---
        let nrows = 10;
        let ncols = 10;
        let mut csv_content = String::new();
        // header
        csv_content.push_str(
            &(0..ncols)
                .map(|i| format!("Col{}", i))
                .collect::<Vec<_>>()
                .join(","),
        );
        csv_content.push('\n');
        // data
        let mut val = 1;
        for _ in 0..nrows {
            let row = (0..ncols)
                .map(|_| {
                    let s = val.to_string();
                    val += 1;
                    s
                })
                .collect::<Vec<_>>()
                .join(",");
            csv_content.push_str(&row);
            csv_content.push('\n');
        }

        // write CSV
        let mut f = File::create(&csv_path)?;
        f.write_all(csv_content.as_bytes())?;

        // --- Create empty factors JSON ---
        let mut f_factors = File::create(&factors_path)?;
        f_factors.write_all(b"[]")?;

        // --- Load data ---
        let categorical_cols = HashSet::<String>::new();
        let data = DataTable::from_file(&csv_path, b',', categorical_cols, &factors_path)?;

        // --- Split ---
        let (train, test) = data.train_test_split(0.7);

        // --- Assertions ---
        assert_eq!(
            train.numeric_data.nrows() + test.numeric_data.nrows(),
            nrows,
            "Row counts must add up"
        );
        assert_eq!(train.headers, data.headers);
        assert_eq!(test.headers, data.headers);

        // Check that all rows are unique across train and test
        let train_rows: Vec<Vec<f64>> = train
            .numeric_data
            .rows()
            .into_iter()
            .map(|r| r.to_vec())
            .collect();
        let test_rows: Vec<Vec<f64>> = test
            .numeric_data
            .rows()
            .into_iter()
            .map(|r| r.to_vec())
            .collect();

        for row in &train_rows {
            assert!(
                !test_rows.contains(row),
                "Row appeared in both train and test!"
            );
        }

        Ok(())
    }

    #[test]
    fn test_update_value_and_update_value_str() {
        // --- Setup ---
        let mut data = DataTable::default();

        // Add numeric and factor columns
        data.add_dataset("age", false, Some(10));
        data.add_dataset("treatment", true, None);

        // --- 1️⃣ Normal numeric update ---
        assert!(data.update_value("age", 0, 42.0));
        assert_eq!(data.numeric_data[[0, 0]], 42.0);

        // --- 2️⃣ Out-of-bounds numeric update ---
        assert!(!data.update_value("age", 10, 46.0)); // Should print error, return false

        // --- 3️⃣ Factor update (new level insertion) ---
        assert!(data.update_value_str("treatment", 0, "DrugA"));
        assert!(data.update_value_str("treatment", 1, "DrugB"));

        // Factor should contain both levels now
        let factor = data.factors.get("treatment").unwrap();
        let levels = factor.get_levels();
        assert_eq!(levels, vec!["DrugA", "DrugB"]);

        // Check numeric encoding values (0.0, 1.0)
        let row_idx = data.headers.iter().position(|h| h == "treatment").unwrap();
        assert_eq!(data.numeric_data[[row_idx, 0]], 0.0);
        assert_eq!(data.numeric_data[[row_idx, 1]], 1.0);

        // --- 4️⃣ Numeric update via update_value_str ---
        assert!(data.update_value_str("age", 1, "39.5"));
        assert_eq!(data.numeric_data[[0, 1]], 39.5);

        // --- 5️⃣ Non-factor string value to numeric column (should fail gracefully) ---
        assert!(!data.update_value_str("age", 2, "NaNish"));

        // --- 6️⃣ Out-of-bounds update via update_value_str ---
        assert!(!data.update_value_str("age", 10, "46"));

        // --- 7️⃣ Factor update with existing value (should not create new level) ---
        assert!(
            data.update_value_str("treatment", 2, "DrugA"),
            "update_value_str() failed for dataset='treatment', column=2, value='DrugA'. \
             Factors present: {:?}, headers: {:?}, numeric shape: {:?}",
            data.factors.keys().collect::<Vec<_>>(),
            data.headers,
            data.numeric_data.dim()
        );
        let factor = data.factors.get("treatment").unwrap();
        let levels = factor.get_levels();
        assert_eq!(levels, vec!["DrugA", "DrugB"]); // still only two levels

        // --- 8️⃣ Unknown column name ---
        assert!(!data.update_value_str("unknown", 0, "123.0"));

        println!("{}", data);
    }

    #[test]
    fn test_order() {
        let mut sd = DataTable::default();

        // 1️⃣ create a new dataset column with 10 entries
        sd.add_dataset("group_000", true, Some(10));
        assert!(sd.factors.contains_key("group_000"));
        let _ = sd.factors.get("group_000").unwrap();
        assert!(
            sd.numeric_data[[0, 3]].is_nan(),
            "initially we have NaN here"
        );

        // 2️⃣ update one entry (name, col, value)
        let changed = sd.update_value("group_000", 3, 1.0);
        assert!(changed);

        // 3️⃣ updating again with the same value should return false
        let changed_again = sd.update_value("group_000", 3, 1.0);
        assert!(!changed_again);

        // 4️⃣ add the order column and update order
        sd.add_dataset("group_000_order", false, None);
        assert!(sd.headers.contains(&"group_000_order".to_string()));

        sd.update_order("group_000_order", 3);
        assert_eq!(sd.numeric_data[[1, 3]], 0.0, "First order gets 0.0");
        sd.update_order("group_000_order", 1);
        assert_eq!(sd.numeric_data[[1, 1]], 1.0, "Second order gets 1.0");

        // 5️⃣ verify reset_order() resets all to "0"
        sd.reset_order();
        assert_eq!(sd.next_order, 0.0, "reset resets to zero");
    }
}
