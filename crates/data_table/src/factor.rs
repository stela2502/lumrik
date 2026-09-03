use anyhow::Result;
use ndarray::Array2;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

/// Only save labels in JSON
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FactorJson {
    pub column: String,
    pub levels: Vec<String>,
    pub numeric: Option<Vec<f64>>,
    pub matching: Option<Vec<String>>,
    pub one_hot: bool, // NEW
}

#[derive(Debug, Clone)]
pub struct Factor {
    pub column_name: String,
    levels: Vec<String>,                                // level labels
    level_to_index: HashMap<String, f64>,               // fast lookup
    index_to_level: HashMap<OrderedFloat<f64>, String>, // <-- instead of a indices vector
    matching: Option<Vec<String>>, // this could match to multiple column names. Like SNP or something
    pub one_hot: bool,             // NEW
}

impl fmt::Display for Factor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Factor '{}':", self.column_name)?;
        writeln!(f, "  One-hot: {}", self.one_hot)?;
        writeln!(f, "  Levels: {:?}", self.levels)?;
        writeln!(f, "  Matching: {:?}", self.matching)?;
        writeln!(f, "  All columns: {:?}", self.all_column_names())
    }
}

impl Factor {
    /// Create a new Factor with a column name and one_hot flag
    pub fn new(column_name: &str, one_hot: bool) -> Self {
        Factor {
            column_name: column_name.to_string(),
            levels: Vec::new(),
            level_to_index: HashMap::new(),
            index_to_level: HashMap::new(),
            matching: None,
            one_hot,
        }
    }

    /*
    /// Subset a factor to only the levels actually observed in `row_data` (numeric values).
    /// Returns a new Factor with:
    /// - levels: filtered in order of first appearance
    /// - level_to_index: numeric values preserved from the original factor
    pub fn subset(&self, data: &Array2<f64>, col_id: usize ) -> Self {
        let mut new_factor = Self::new(&self.column_name, self.one_hot);

        new_factor.matching = self.matching.clone();

        // Reconstruct levels vector, level_to_index, and numeric for new factor
        for &val in data.column(col_id).iter() {
            if val.is_nan() {
                continue;
            }

            let key = val.to_string();

            if !new_factor.level_to_index.contains_key(&key) {
                if let Some(num_val) = self.index_to_level.get(&OrderedFloat(val)) {
                    new_factor.levels.push(key.clone());
                    new_factor.level_to_index.insert(key.clone(), val);
                    new_factor.index_to_level.insert(OrderedFloat(val), key.clone());
                }else {
                    panic!("Error: Trying to insert a previousely not existing key '{}'\nfrom\n{}\nTo\n{}",
                        key, self, new_factor );
                }
            }
        }
        new_factor
    }
    */

    /// Compare level names between two Factors.
    /// Returns Ok(()) if identical, or Err(msg) describing the mismatch.
    pub fn compare_levels(&self, other: &Self) -> Result<(), String> {
        if self.levels.len() != other.levels.len() {
            return Err(format!(
                "Level count mismatch: expected {} vs found {}",
                self.levels.len(),
                other.levels.len()
            ));
        }

        for (i, (a, b)) in self.levels.iter().zip(other.levels.iter()).enumerate() {
            if a != b {
                return Err(format!(
                    "Mismatch at index {}: expected '{}' but found '{}'",
                    i, a, b
                ));
            }
        }

        Ok(())
    }

    /// Subset factor with **extreme verbosity**.
    pub fn subset(&self, data: &Array2<f64>, col_id: usize) -> Self {
        //println!("➡️ Starting subset for factor '{}', column id {}", self.column_name, col_id);

        let mut new_factor = Self::new(&self.column_name, self.one_hot);
        new_factor.matching = self.matching.clone();

        /*println!("    Original factor levels: {:?}", self.levels);
        println!("    Original level_to_index: {:?}", self.level_to_index);
        println!("    Original index_to_level: {:?}", self.index_to_level);*/

        for (i, &val) in data.column(col_id).iter().enumerate() {
            //println!("📍 Row {}: numeric value = {}", i, val);

            if val.is_nan() {
                //println!("    ❌ Skipping NaN value");
                continue;
            }

            let lvl = match self.index_to_level.get(&OrderedFloat(val)) {
                Some(s) => {
                    //println!("    ✅ Found original level '{}'", s);
                    s.clone()
                }
                None => panic!(
                    "❌ Unknown numeric value {} in column '{}' at row {}",
                    val, self.column_name, i
                ),
            };

            if !new_factor.index_to_level.contains_key(&OrderedFloat(val)) {
                //println!("    ➕ Adding new level '{}' with numeric {}", lvl, val);
                new_factor.levels.push(lvl.clone());
                new_factor.level_to_index.insert(lvl.clone(), val);
                new_factor
                    .index_to_level
                    .insert(OrderedFloat(val), lvl.clone());
            }
        }

        /*println!("🏁 Subset complete!");
        println!("    New factor levels: {:?}", new_factor.levels);
        println!("    New level_to_index: {:?}", new_factor.level_to_index);
        println!("    New index_to_level: {:?}", new_factor.index_to_level);
        */
        new_factor
    }

    pub fn extra_columns(&self) -> usize {
        if self.one_hot { self.levels.len() } else { 0 }
    }

    /// returns all one_hot column names or the original colum, name only
    pub fn all_column_names(&self) -> Vec<String> {
        if self.one_hot {
            self.levels
                .iter()
                .map(|lvl| self.build_one_hot_column(lvl))
                .collect()
        } else {
            //vec![]
            vec![self.column_name.to_string()]
        }
    }

    fn build_one_hot_column(&self, value: &str) -> String {
        if self.one_hot {
            format!("{}_{}", self.column_name, value)
        } else {
            self.column_name.clone()
        }
    }

    // Return (and possibly insert) the numeric index for a given level.
    ///
    /// If the level does not exist yet, it is added to `levels`,
    /// and the internal maps (`level_to_index` and `index_to_level`) are updated.
    ///
    /// Special cases:
    /// - Empty strings or "NA" → return NaN.
    pub fn level_to_index(&mut self, value: &str) -> f64 {
        let trimmed = value.trim();

        // Handle missing values
        if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("NA") {
            return f64::NAN;
        }

        // For one-hot encoding, presence means 1.0 — but we still must register the level
        if self.one_hot {
            if !self.level_to_index.contains_key(trimmed) {
                let new_idx = self.levels.len() as f64;
                self.levels.push(trimmed.to_string());
                self.level_to_index.insert(trimmed.to_string(), new_idx);
                self.index_to_level
                    .insert(OrderedFloat(new_idx), trimmed.to_string());
            }
            return 1.0;
        }

        // For dense encoding: assign incremental indices
        if let Some(&idx) = self.level_to_index.get(trimmed) {
            idx
        } else {
            let new_idx = self.levels.len() as f64;
            self.levels.push(trimmed.to_string());
            self.level_to_index.insert(trimmed.to_string(), new_idx);
            self.index_to_level
                .insert(OrderedFloat(new_idx), trimmed.to_string());
            new_idx
        }
    }

    /// Push a value for this factor.
    ///
    /// Returns:
    /// - f64: numeric index (or binary for one-hot)
    /// - String: name of the column to update
    /// - Option<Vec<String>>: all one-hot columns (if applicable)
    pub fn push(&mut self, value: &str) -> (f64, String, Option<Vec<String>>) {
        let trimmed = value.trim();

        if self.one_hot {
            // Ensure level is registered and retrieve index behavior
            let idx = self.level_to_index(trimmed);

            // all other levels become zero columns
            let zero_cols = self.all_column_names();
            let zero_cols_option = if zero_cols.len() == 1 {
                None
            } else {
                Some(zero_cols)
            };

            (idx, self.build_one_hot_column(trimmed), zero_cols_option)
        } else {
            let idx = self.level_to_index(trimmed);
            (idx, self.column_name.to_string(), None)
        }
    }

    pub fn get_f64(&self, trimmed: &str) -> f64 {
        *self.level_to_index.get(trimmed).unwrap_or(&f64::NAN)
    }

    pub fn get_levels(&self) -> &[String] {
        &self.levels
    }
    pub fn get_string(&self, value: f64) -> String {
        match self.index_to_level.get(&OrderedFloat(value)) {
            Some(string) => string.to_string(),
            None => "NA".to_string(),
        }
    }

    /// Create Factor from FactorDef
    pub fn from_def(def: &FactorJson) -> Self {
        let mut level_to_index: std::collections::HashMap<String, f64> = HashMap::new();
        let mut index_to_level: std::collections::HashMap<OrderedFloat<f64>, String> =
            HashMap::new();

        match &def.numeric {
            Some(numbers) => {
                for (lvl, &num) in def.levels.iter().zip(numbers.iter()) {
                    level_to_index.insert(lvl.clone(), num);
                    index_to_level.insert(OrderedFloat(num), lvl.clone());
                }
            }
            None => {
                for (i, lvl) in def.levels.iter().enumerate() {
                    let num = i as f64;
                    level_to_index.insert(lvl.clone(), num);
                    index_to_level.insert(OrderedFloat(num), lvl.clone());
                }
            }
        }

        Factor {
            column_name: def.column.to_string(),
            levels: def.levels.clone(),
            level_to_index,
            index_to_level,
            matching: def.matching.clone(),
            one_hot: def.one_hot,
        }
    }

    /// Convert Factor into JSON representation
    pub fn as_json(&self, column: &str) -> FactorJson {
        let numeric: Vec<f64> = self
            .levels
            .iter()
            .map(|name| {
                self.level_to_index.get(name).copied().unwrap_or_else(|| {
                    panic!(
                        "Factor::as_json error: level '{}' not found in mapping for column '{}'",
                        name, column
                    )
                })
            })
            .collect();

        FactorJson {
            column: self.column_name.to_string(),
            levels: self.levels.clone(),
            numeric: Some(numeric),
            matching: self.matching.clone(),
            one_hot: self.one_hot,
        }
    }
    /// Modify numeric values for levels in this factor.
    ///
    /// # Arguments
    /// * `numeric_values` - New numeric values to assign.
    /// * `levels_to_change` - Optional subset of levels to modify. If `None`, all levels are overwritten in order.
    pub fn modify_levels(
        &mut self,
        numeric_values: &[f64],
        levels_to_change: Option<&[String]>,
    ) -> Result<(), String> {
        match levels_to_change {
            Some(levels) => {
                // Modify only the specified levels
                for (lvl, &val) in levels.iter().zip(numeric_values.iter()) {
                    if self.levels.contains(lvl) {
                        self.level_to_index.insert(lvl.clone(), val);
                        self.index_to_level.insert(OrderedFloat(val), lvl.clone());
                    } else {
                        return Err(format!(
                            "Level '{}' does not exist in factor '{}'",
                            lvl, self.column_name
                        ));
                    }
                }
            }
            None => {
                // Overwrite all levels in order
                if numeric_values.len() != self.levels.len() {
                    return Err(format!(
                        "Length mismatch: {} numeric values, but factor '{}' has {} levels",
                        numeric_values.len(),
                        self.column_name,
                        self.levels.len()
                    ));
                }
                for (lvl, &val) in self.levels.iter().zip(numeric_values.iter()) {
                    self.level_to_index.insert(lvl.clone(), val);
                    self.index_to_level.insert(OrderedFloat(val), lvl.clone());
                }
            }
        }
        Ok(())
    }
}


#[cfg(test)]
mod tests {

    use super::*;
    use ndarray::array;
    use ordered_float::OrderedFloat;
    use std::collections::HashMap;

    use std::path::PathBuf;

    use crate::{DataTable, Factor};
    

    #[test]
    fn test_save_and_load_factors() -> Result<(), Box<dyn std::error::Error>> {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let factors_path = dir.as_path().join("factors.json");

        // --- Create Factors ---
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

        // --- Insert into DataTable ---
        let mut data = DataTable::default();
        data.factors.insert("Color".to_string(), color);
        data.factors.insert("Number".to_string(), number);

        // --- Save to JSON ---
        data.save_factors(&factors_path)?;

        // --- Load from JSON ---
        data.factors = HashMap::new();
        data.load_factors(&factors_path)?;
        let loaded_factors = &data.factors;

        // Check factors exist
        assert!(loaded_factors.contains_key("Color"));
        assert!(loaded_factors.contains_key("Number"));

        // Check levels
        let color = loaded_factors.get("Color").unwrap();
        assert_eq!(color.get_levels(), vec!["Red", "Blue", "Green"]);
        assert!(color.one_hot);

        let number = loaded_factors.get("Number").unwrap();
        assert_eq!(number.get_levels(), vec!["1", "2", "3"]);
        assert!(number.one_hot);

        Ok(())
    }


    #[test]
    fn test_modify_levels() {
        // --- Setup ---
        let mut factor = Factor::new("Color", false);
        let _ = factor.push("Red");
        let _ = factor.push("Blue");
        let _ = factor.push("Green");

        // --- Modify subset of levels ---
        let levels_to_change = vec!["Red".to_string(), "Green".to_string()];
        let new_values = vec![10.0, 30.0];
        let _ = factor.modify_levels(&new_values, Some(&levels_to_change));

        assert_eq!(*factor.level_to_index.get("Red").unwrap(), 10.0);
        assert_eq!(*factor.level_to_index.get("Green").unwrap(), 30.0);
        assert_eq!(*factor.level_to_index.get("Blue").unwrap(), 1.0); // unchanged

        // --- Overwrite all levels ---
        let all_new_values = vec![100.0, 200.0, 300.0];
        let _ = factor.modify_levels(&all_new_values, None);

        assert_eq!(*factor.level_to_index.get("Red").unwrap(), 100.0);
        assert_eq!(*factor.level_to_index.get("Blue").unwrap(), 200.0);
        assert_eq!(*factor.level_to_index.get("Green").unwrap(), 300.0);

        // --- Test panic for non-existing level ---

        assert!(
            factor
                .modify_levels(&[1.0], Some(&["Yellow".to_string()]))
                .is_err(),
            "Should panic when modifying non-existing level"
        );

        assert!(
            factor.modify_levels(&[1.0, 2.0], None).is_err(),
            "Should panic for length mismatch"
        );
    }

    #[test]
    fn test_push_one_hot() {
        let mut f = Factor::new("Color", false);
        assert_eq!(
            f.push("Red"),
            (0.0, "Color".to_string(), None),
            "new factor Red gets 0"
        );
        assert_eq!(
            f.push("Blue"),
            (1.0, "Color".to_string(), None),
            "new factor Blue gets 1"
        );
        let (idx, col, opt) = f.push("NA");
        assert!(idx.is_nan(), "new factor NA should be NaN");
        assert_eq!(col, "Color");
        assert_eq!(opt, None);
        f.one_hot = true;

        // Push first value "Red"
        let (idx, col_to_add, all_cols) = f.push("Red");
        assert_eq!(idx, 1.0);
        assert_eq!(col_to_add, "Color_Red");
        assert_eq!(
            all_cols.unwrap(),
            vec!["Color_Red".to_string(), "Color_Blue".to_string()]
        );

        // Push another value "Blue"
        let (idx, col_to_add, all_cols) = f.push("Blue");
        assert_eq!(idx, 1.0);
        assert_eq!(col_to_add, "Color_Blue");
        assert_eq!(
            all_cols.unwrap(),
            vec!["Color_Red".to_string(), "Color_Blue".to_string()]
        );

        // Now push a missing value
        let (idx, col_to_add, all_cols) = f.push("NA");
        assert!(idx.is_nan());
        assert_eq!(col_to_add, "Color_NA");
        assert_eq!(
            all_cols.unwrap(),
            vec!["Color_Red".to_string(), "Color_Blue".to_string()]
        );
    }

    #[test]
    fn test_push_categorical_indexed() {
        let mut f = Factor::new("Color", false);
        assert_eq!(
            f.push("Red"),
            (0.0, "Color".to_string(), None),
            "new factor Red gets 0"
        );
        assert_eq!(
            f.push("Blue"),
            (1.0, "Color".to_string(), None),
            "new factor Blue gets 1"
        );
        let (idx, col, opt) = f.push("NA");
        assert!(idx.is_nan(), "new factor NA should be NaN");
        assert_eq!(col, "Color");
        assert_eq!(opt, None);

        assert_eq!(
            f.push("Red"),
            (0.0, "Color".to_string(), None),
            "2# new factor Red gets 0"
        );
        assert_eq!(
            f.push("Blue"),
            (1.0, "Color".to_string(), None),
            "2# new factor Blue gets 1"
        );
        let (idx, col, opt) = f.push("NA");
        assert!(idx.is_nan(), "2# new factor NA should be NaN");
        assert_eq!(col, "Color");
        assert_eq!(opt, None);
    }

    #[test]
    fn test_factor_subset_verbose() {
        // --- Step 1: Create original factor ---
        let mut factor = Factor::new("Color", false);
        factor.levels = vec!["Red".into(), "Blue".into(), "Green".into()];
        factor.level_to_index = vec![
            ("Red".into(), 0.0),
            ("Blue".into(), 1.0),
            ("Green".into(), 2.0),
        ]
        .into_iter()
        .collect();
        factor.index_to_level = vec![
            (OrderedFloat(0.0), "Red".into()),
            (OrderedFloat(1.0), "Blue".into()),
            (OrderedFloat(2.0), "Green".into()),
        ]
        .into_iter()
        .collect();

        // --- Step 2: Create row data (Array2) ---
        let data = array![
            [0.0],      // Red
            [2.0],      // Green
            [2.0],      // Green
            [f64::NAN], // Missing
            [0.0],      // Red
        ];

        // --- Step 3: Subset factor using column 0 ---
        let mut new_factor = factor.subset(&data, 0);

        // --- Step 4: Assertions ---
        let expected_levels = vec!["Red".to_string(), "Green".to_string()];
        assert_eq!(
            new_factor.levels, expected_levels,
            "Levels should preserve first appearance order"
        );

        let mut expected_level_to_index: HashMap<String, f64> = HashMap::new();
        expected_level_to_index.insert("Red".into(), 0.0);
        expected_level_to_index.insert("Green".into(), 2.0);
        //expected_level_to_index.insert("Green".into(), 2.0);
        assert_eq!(
            new_factor.level_to_index, expected_level_to_index,
            "level_to_index mapping should be preserved"
        );

        let mut expected_index_to_level: HashMap<OrderedFloat<f64>, String> = HashMap::new();
        expected_index_to_level.insert(OrderedFloat(0.0), "Red".into());
        expected_index_to_level.insert(OrderedFloat(2.0), "Green".into());
        //expected_index_to_level.insert(OrderedFloat(2.0), "Green".into());
        assert_eq!(
            new_factor.index_to_level, expected_index_to_level,
            "index_to_level mapping should be preserved"
        );

        assert_eq!(
            new_factor.push("Green"),
            factor.push("Green"),
            "Get the same value for 'Blue' from both factors"
        );

        // --- Step 5: Print for braggy verification ---
        println!("Original factor: {:?}", factor);
        println!("Subset factor: {:?}", new_factor);
    }
}
