// matrix_value.rs
use core::fmt;

/// Supported MatrixMarket numeric/value classes.
///
/// This is used to decide how matrix entries are interpreted during
/// import/export, and which MatrixMarket header should be written.
#[derive(Debug, PartialEq, Clone)]
pub enum MatrixValueType {
    /// Integer-valued matrix entries.
    Integer,

    /// Real-valued matrix entries.
    Real,

    /// Complex-valued matrix entries.
    Complex,

    /// Pattern matrix without explicit stored values.
    Pattern,

    /// Any not yet recognized MatrixMarket type string.
    Unknown(String),
}

impl MatrixValueType {
    /// Return the canonical MatrixMarket type keyword.
    pub fn as_matrix_market_str(&self) -> &str {
        match self {
            MatrixValueType::Integer => "integer",
            MatrixValueType::Real => "real",
            MatrixValueType::Complex => "complex",
            MatrixValueType::Pattern => "pattern",
            MatrixValueType::Unknown(_) => "unknown",
        }
    }

    /// Parse the MatrixMarket type keyword from a header token.
    pub fn from_matrix_market_str(value: &str) -> Self {
        match value.trim().to_ascii_lowercase().as_str() {
            "integer" => MatrixValueType::Integer,
            "real" => MatrixValueType::Real,
            "complex" => MatrixValueType::Complex,
            "pattern" => MatrixValueType::Pattern,
            other => MatrixValueType::Unknown(other.to_string()),
        }
    }

    /// Return true if this type is currently supported for sparse export.
    pub fn supports_sparse_write(&self) -> bool {
        matches!(self, MatrixValueType::Integer | MatrixValueType::Real)
    }

    /// Return true if this type is currently supported for MatrixMarket import.
    pub fn supports_matrix_market_read(&self) -> bool {
        matches!(self, MatrixValueType::Integer | MatrixValueType::Real)
    }

    /// format a internal f32 the right way.
    pub fn format_matrix_market_value(&self, value: f32) -> Result<String, String> {
        match self {
            MatrixValueType::Integer => {
                if !value.is_finite() {
                    return Err("non-finite value in integer matrix".to_string());
                }
                if value.fract() != 0.0 {
                    return Err(format!(
                        "non-integer value {value} encountered in MatrixMarket integer matrix"
                    ));
                }
                Ok((value as i64).to_string())
            }
            MatrixValueType::Real => {
                if !value.is_finite() {
                    return Err("non-finite value in real matrix".to_string());
                }
                Ok(value.to_string())
            }
            MatrixValueType::Complex => {
                Err("complex MatrixMarket writing is not implemented".to_string())
            }
            MatrixValueType::Pattern => Ok(String::new()),
            MatrixValueType::Unknown(s) => {
                Err(format!("cannot write unknown MatrixMarket type: {s}"))
            }
        }
    }
}

impl fmt::Display for MatrixValueType {
    /// Print the MatrixMarket keyword form of the value type.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_matrix_market_str())
    }
}

#[cfg(test)]
mod matrix_value_type_tests {
    use super::MatrixValueType;

    #[test]
    fn format_integer_from_f32() {
        let ty = MatrixValueType::Integer;
        assert_eq!(ty.format_matrix_market_value(20.0).unwrap(), "20");
        assert_eq!(ty.format_matrix_market_value(0.0).unwrap(), "0");
    }

    #[test]
    fn reject_fractional_value_for_integer_matrix() {
        let ty = MatrixValueType::Integer;
        let err = ty.format_matrix_market_value(2.25).unwrap_err();
        assert!(err.contains("non-integer value"));
    }

    #[test]
    fn reject_non_finite_for_integer_matrix() {
        let ty = MatrixValueType::Integer;
        assert!(ty.format_matrix_market_value(f32::NAN).is_err());
        assert!(ty.format_matrix_market_value(f32::INFINITY).is_err());
    }

    #[test]
    fn format_real_from_f32() {
        let ty = MatrixValueType::Real;
        assert_eq!(ty.format_matrix_market_value(2.5).unwrap(), "2.5");
    }

    #[test]
    fn reject_non_finite_for_real_matrix() {
        let ty = MatrixValueType::Real;
        assert!(ty.format_matrix_market_value(f32::NAN).is_err());
        assert!(ty.format_matrix_market_value(f32::NEG_INFINITY).is_err());
    }

    #[test]
    fn pattern_returns_empty_string() {
        let ty = MatrixValueType::Pattern;
        assert_eq!(ty.format_matrix_market_value(1.0).unwrap(), "");
    }

    #[test]
    fn complex_is_not_implemented() {
        let ty = MatrixValueType::Complex;
        assert!(ty.format_matrix_market_value(1.0).is_err());
    }

    #[test]
    fn unknown_is_rejected() {
        let ty = MatrixValueType::Unknown("weird".to_string());
        assert!(ty.format_matrix_market_value(1.0).is_err());
    }
}
