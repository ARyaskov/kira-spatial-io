//! Canonical JSON serialisation for metadata.

use std::io::Write;

use serde_json::{Map, Value};

use crate::error::SpatialIoError;

/// Returns recursively canonicalized JSON with sorted object keys.
pub fn canonicalize_json(v: &Value) -> Value {
    match v {
        Value::Null => Value::Null,
        Value::Bool(b) => Value::Bool(*b),
        Value::Number(n) => Value::Number(n.clone()),
        Value::String(s) => Value::String(s.clone()),
        Value::Array(arr) => Value::Array(arr.iter().map(canonicalize_json).collect()),
        Value::Object(obj) => {
            let mut keys: Vec<&str> = obj.keys().map(|k| k.as_str()).collect();
            keys.sort_unstable();

            let mut canonical = Map::with_capacity(obj.len());
            for key in keys {
                canonical.insert(key.to_string(), canonicalize_json(&obj[key]));
            }
            Value::Object(canonical)
        }
    }
}

/// Writes canonical JSON, rejecting NaN/Inf in numeric values.
pub fn write_canonical_json<W: Write>(w: &mut W, v: &Value) -> Result<(), SpatialIoError> {
    validate_json_numbers(v)?;
    serde_json::to_writer(w, v).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!("failed to serialize canonical json: {e}"))
    })
}

fn validate_json_numbers(v: &Value) -> Result<(), SpatialIoError> {
    match v {
        Value::Null | Value::Bool(_) | Value::String(_) => Ok(()),
        Value::Number(n) => {
            if n.is_i64() || n.is_u64() {
                return Ok(());
            }
            if let Some(f) = n.as_f64()
                && f.is_finite()
            {
                return Ok(());
            }
            Err(SpatialIoError::InvalidFloat(
                "metadata contains non-finite number (NaN or Infinity)".to_string(),
            ))
        }
        Value::Array(arr) => {
            for item in arr {
                validate_json_numbers(item)?;
            }
            Ok(())
        }
        Value::Object(obj) => {
            for value in obj.values() {
                validate_json_numbers(value)?;
            }
            Ok(())
        }
    }
}
