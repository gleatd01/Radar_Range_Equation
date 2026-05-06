use crate::{RadarError, Result};

pub fn ensure_positive(value: f64, name: &str) -> Result<()> {
    if value <= 0.0 {
        return Err(RadarError::InvalidParameter(format!(
            "{} must be positive",
            name
        )));
    }
    Ok(())
}

pub fn ensure_unit_interval(value: f64, name: &str) -> Result<()> {
    if value <= 0.0 || value > 1.0 {
        return Err(RadarError::InvalidParameter(format!(
            "{} must be between 0 and 1",
            name
        )));
    }
    Ok(())
}
