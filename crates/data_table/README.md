# DataTable Rust Library

## Overview

The `DataTable` crate provides a flexible framework for reading, processing, and managing survival or biomedical data tables that contain both **numeric** and **categorical** (factor) variables.  
It is particularly suited for bioinformatics, epidemiology, and statistical modeling pipelines where tabular datasets must be cleaned, imputed, split, and exported with precise factor-level tracking.

The library handles automatic detection and management of categorical variables through *factor definitions*, supports one-hot encoding, and includes utilities for data filtering, imputation, and persistence to CSV and JSON formats.

---

## Features

- **CSV Input Parsing**
  - Reads numeric and categorical columns from delimited files.
  - Automatically infers factor definitions from provided JSON or generates them if missing.
  - One-hot expansion for categorical variables with flexible level mapping.

- **Factor Management**
  - `Factor` objects track categorical levels, numeric encodings, and one-hot states.
  - Load and save factors to JSON for reproducibility.

- **Data Cleaning**
  - Remove rows or columns with missing values (`NaN`).
  - Impute missing values via K-Nearest Neighbours (`impute_knn`).
  - Filter low-variance columns or features with excessive missingness.

- **Data Export**
  - Write numeric and factor-expanded data to CSV.
  - Save and reload factor definitions to/from JSON.

- **Utilities**
  - Split datasets into training and test sets with random shuffling.
  - Access individual columns as `Vec<f64>`, `Vec<u8>`, or categorical strings.
  - Print concise dataset summaries for debugging and inspection.

---

## Example Usage

```rust
use std::collections::HashSet;
use survival_data::DataTable;

fn main() -> anyhow::Result<()> {
    let file_path = "data/example.csv";
    let factors_path = "data/factors.json";
    let categorical: HashSet<String> = ["status".into(), "sex".into()].into();

    // Load data, automatically building factors if needed
    let data = DataTable::from_file(file_path, b',', categorical, factors_path)?;

    // Clean and impute data
    let usable_features = data.filter_features_by_na(0.1);
    let mut filtered_data = data.clone();
    filtered_data.filter_all_na_rows(&usable_features);
    filtered_data.impute_knn(3, 5, true);

    // Split for training/testing
    let (train, test) = filtered_data.train_test_split(0.7);

    // Save cleaned data
    train.to_file("data/train.csv", b',')?;
    test.to_file("data/test.csv", b',');
    Ok(())
}
```

---

## Example Factor Definition (JSON)

```json
[
  {
    "column": "status",
    "levels": ["alive", "dead"],
    "numeric": [0.0, 1.0],
    "matching": null,
    "one_hot": false
  },
  {
    "column": "treatment",
    "levels": ["control", "drugA", "drugB"],
    "numeric": [0.0, 1.0, 2.0],
    "matching": null,
    "one_hot": true
  }
]
```

---

## Create Win and Linux binaries on Linux:

On Ubuntu you can simply:

```
rustup target add x86_64-unknown-linux-musl
rustup target add x86_64-pc-windows-gnu
```

And then install the missing linkers with
```
sudo apt update
sudo apt install musl-tools mingw-w64
```

Later on you can compile for both systems using

```
cargo build --release --target x86_64-unknown-linux-musl
cargo build --release --target x86_64-pc-windows-gnu
```


If you need to explain this crate quickly in a chat session, say:

> “This Rust library (`DataTable`) helps parse and preprocess survival or biomedical datasets. It automatically handles categorical factors (with optional one-hot encoding), performs data cleaning and imputation, and supports consistent factor definitions via JSON. It’s ideal for preparing structured data for downstream modeling or statistical analysis.”

---

## License

MIT License © 2025
