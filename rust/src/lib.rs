//! Radar Range Equation Calculator
//!
//! This library provides tools for calculating radar range equations
//! and related radar system parameters.

use std::f64::consts::PI;

/// Error type for radar range equation calculations
#[derive(Debug, Clone, PartialEq)]
pub enum RadarError {
    /// One or more parameters are non-positive
    InvalidParameter(String),
}

impl std::fmt::Display for RadarError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            RadarError::InvalidParameter(msg) => write!(f, "Invalid parameter: {}", msg),
        }
    }
}

impl std::error::Error for RadarError {}

/// Calculate the maximum detection range of a radar system.
///
/// Uses the radar range equation to determine the maximum range at which
/// a target can be detected.
///
/// # Arguments
///
/// * `transmit_power` - Transmitted power in Watts
/// * `antenna_gain` - Antenna gain (dimensionless, linear scale)
/// * `wavelength` - Radar wavelength in meters
/// * `radar_cross_section` - Target radar cross section in square meters
/// * `min_detectable_signal` - Minimum detectable signal power in Watts
///
/// # Returns
///
/// Maximum range in meters
///
/// # Errors
///
/// Returns `RadarError::InvalidParameter` if any input parameters are non-positive
///
/// # Examples
///
/// ```
/// use radar_range_equation::calculate_max_range;
///
/// let max_range = calculate_max_range(
///     1000.0,    // 1 kW
///     1000.0,    // 30 dBi
///     0.03,      // 10 GHz
///     1.0,       // 1 m²
///     1e-13,     // -100 dBm
/// ).unwrap();
///
/// assert!(max_range > 0.0);
/// ```
pub fn calculate_max_range(
    transmit_power: f64,
    antenna_gain: f64,
    wavelength: f64,
    radar_cross_section: f64,
    min_detectable_signal: f64,
) -> Result<f64, RadarError> {
    if transmit_power <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "transmit_power must be positive".to_string(),
        ));
    }
    if antenna_gain <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "antenna_gain must be positive".to_string(),
        ));
    }
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "wavelength must be positive".to_string(),
        ));
    }
    if radar_cross_section <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "radar_cross_section must be positive".to_string(),
        ));
    }
    if min_detectable_signal <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "min_detectable_signal must be positive".to_string(),
        ));
    }

    let numerator = transmit_power
        * antenna_gain.powi(2)
        * wavelength.powi(2)
        * radar_cross_section;
    let denominator = (4.0 * PI).powi(3) * min_detectable_signal;

    let max_range = (numerator / denominator).powf(0.25);

    Ok(max_range)
}

/// Calculate the received power from a target at a given range.
///
/// # Arguments
///
/// * `transmit_power` - Transmitted power in Watts
/// * `antenna_gain` - Antenna gain (dimensionless, linear scale)
/// * `wavelength` - Radar wavelength in meters
/// * `radar_cross_section` - Target radar cross section in square meters
/// * `range_m` - Range to target in meters
///
/// # Returns
///
/// Received power in Watts
///
/// # Errors
///
/// Returns `RadarError::InvalidParameter` if any input parameters are non-positive
///
/// # Examples
///
/// ```
/// use radar_range_equation::calculate_received_power;
///
/// let power = calculate_received_power(
///     1000.0,    // 1 kW
///     1000.0,    // 30 dBi
///     0.03,      // 10 GHz
///     1.0,       // 1 m²
///     10000.0,   // 10 km
/// ).unwrap();
///
/// assert!(power > 0.0);
/// ```
pub fn calculate_received_power(
    transmit_power: f64,
    antenna_gain: f64,
    wavelength: f64,
    radar_cross_section: f64,
    range_m: f64,
) -> Result<f64, RadarError> {
    if transmit_power <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "transmit_power must be positive".to_string(),
        ));
    }
    if antenna_gain <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "antenna_gain must be positive".to_string(),
        ));
    }
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "wavelength must be positive".to_string(),
        ));
    }
    if radar_cross_section <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "radar_cross_section must be positive".to_string(),
        ));
    }
    if range_m <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "range_m must be positive".to_string(),
        ));
    }

    let numerator = transmit_power
        * antenna_gain.powi(2)
        * wavelength.powi(2)
        * radar_cross_section;
    let denominator = (4.0 * PI).powi(3) * range_m.powi(4);

    let received_power = numerator / denominator;

    Ok(received_power)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_calculate_max_range_basic() {
        let max_range = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
        assert!(max_range > 0.0);
    }

    #[test]
    fn test_calculate_max_range_negative_power() {
        let result = calculate_max_range(-1000.0, 1000.0, 0.03, 1.0, 1e-13);
        assert!(result.is_err());
    }

    #[test]
    fn test_calculate_max_range_zero_wavelength() {
        let result = calculate_max_range(1000.0, 1000.0, 0.0, 1.0, 1e-13);
        assert!(result.is_err());
    }

    #[test]
    fn test_calculate_max_range_larger_power_increases_range() {
        let range1 = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
        let range2 = calculate_max_range(2000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
        assert!(range2 > range1);
    }

    #[test]
    fn test_calculate_received_power_basic() {
        let power = calculate_received_power(1000.0, 1000.0, 0.03, 1.0, 10000.0).unwrap();
        assert!(power > 0.0);
    }

    #[test]
    fn test_calculate_received_power_negative_range() {
        let result = calculate_received_power(1000.0, 1000.0, 0.03, 1.0, -1000.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_inverse_fourth_power_law() {
        let power1 = calculate_received_power(1000.0, 1000.0, 0.03, 1.0, 10000.0).unwrap();
        let power2 = calculate_received_power(1000.0, 1000.0, 0.03, 1.0, 20000.0).unwrap();

        // Power at 2x range should be 1/16 (2^4) of power at 1x range
        assert_relative_eq!(power2 * 16.0, power1, epsilon = power1 * 0.01);
    }

    #[test]
    fn test_consistency_with_max_range() {
        let transmit_power = 1000.0;
        let antenna_gain = 1000.0;
        let wavelength = 0.03;
        let radar_cross_section = 1.0;
        let min_detectable_signal = 1e-13;

        let max_range = calculate_max_range(
            transmit_power,
            antenna_gain,
            wavelength,
            radar_cross_section,
            min_detectable_signal,
        )
        .unwrap();

        let received_power = calculate_received_power(
            transmit_power,
            antenna_gain,
            wavelength,
            radar_cross_section,
            max_range,
        )
        .unwrap();

        // Should be very close to min_detectable_signal
        assert_relative_eq!(
            received_power,
            min_detectable_signal,
            epsilon = min_detectable_signal * 0.01
        );
    }
}
