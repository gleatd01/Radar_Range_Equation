//! # Radar Range Equation
//!
//! A Rust library for calculating radar range equations.
//!
//! This library provides functions for computing:
//! - Maximum detection range
//! - Received power
//! - Wavelength and frequency conversions
//! - Antenna gain calculations
//! - Various unit conversions
//!
//! ## Example
//!
//! ```
//! use radar_range_equation::*;
//!
//! // Calculate wavelength from frequency
//! let frequency = 10e9; // 10 GHz
//! let wavelength = calculate_wavelength(frequency).unwrap();
//!
//! // Calculate maximum range
//! let max_range = calculate_max_range(
//!     1000.0,  // transmit power (W)
//!     1000.0,  // antenna gain (linear)
//!     wavelength,
//!     1.0,     // radar cross-section (m²)
//!     1e-13,   // min detectable signal (W)
//! ).unwrap();
//! println!("Maximum range: {:.2} meters", max_range);
//! ```

/// Physical constants used in radar calculations
pub mod constants {
    /// Speed of light in m/s
    pub const SPEED_OF_LIGHT: f64 = 299792458.0;
    
    /// Boltzmann constant in J/K
    pub const BOLTZMANN: f64 = 1.380649e-23;
    
    /// Pi
    pub const PI: f64 = std::f64::consts::PI;
    
    /// Beamwidth coefficient for Gaussian approximation
    pub const BEAMWIDTH_COEFFICIENT: f64 = 65.0;
}

/// Errors that can occur during radar calculations
#[derive(Debug, Clone, PartialEq)]
pub enum RadarError {
    /// Invalid input parameter (e.g., negative value where positive expected)
    InvalidParameter(String),
    /// Division by zero or other mathematical error
    MathError(String),
}

impl std::fmt::Display for RadarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RadarError::InvalidParameter(msg) => write!(f, "Invalid parameter: {}", msg),
            RadarError::MathError(msg) => write!(f, "Math error: {}", msg),
        }
    }
}

impl std::error::Error for RadarError {}

/// Result type for radar calculations
pub type Result<T> = std::result::Result<T, RadarError>;

/// Calculate wavelength from frequency
///
/// # Arguments
///
/// * `frequency` - Frequency in Hz
///
/// # Returns
///
/// Wavelength in meters
///
/// # Example
///
/// ```
/// use radar_range_equation::calculate_wavelength;
/// let wavelength = calculate_wavelength(10e9).unwrap(); // 10 GHz
/// assert!((wavelength - 0.0299792458).abs() < 1e-6);
/// ```
pub fn calculate_wavelength(frequency: f64) -> Result<f64> {
    if frequency <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Frequency must be positive".to_string(),
        ));
    }
    Ok(constants::SPEED_OF_LIGHT / frequency)
}

/// Calculate frequency from wavelength
///
/// # Arguments
///
/// * `wavelength` - Wavelength in meters
///
/// # Returns
///
/// Frequency in Hz
pub fn calculate_frequency(wavelength: f64) -> Result<f64> {
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    Ok(constants::SPEED_OF_LIGHT / wavelength)
}

/// Calculate transmit antenna gain
///
/// # Arguments
///
/// * `effective_aperture` - Effective aperture in m²
/// * `wavelength` - Wavelength in meters
///
/// # Returns
///
/// Antenna gain (linear, not dB)
pub fn calculate_antenna_gain(effective_aperture: f64, wavelength: f64) -> Result<f64> {
    if effective_aperture <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Effective aperture must be positive".to_string(),
        ));
    }
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    Ok((4.0 * constants::PI * effective_aperture) / (wavelength * wavelength))
}

/// Calculate effective aperture for rectangular antenna
///
/// # Arguments
///
/// * `efficiency` - Antenna efficiency (0-1)
/// * `horizontal_dimension` - Horizontal dimension in meters
/// * `vertical_dimension` - Vertical dimension in meters
///
/// # Returns
///
/// Effective aperture in m²
pub fn calculate_effective_aperture_rect(
    efficiency: f64,
    horizontal_dimension: f64,
    vertical_dimension: f64,
) -> Result<f64> {
    if efficiency <= 0.0 || efficiency > 1.0 {
        return Err(RadarError::InvalidParameter(
            "Efficiency must be between 0 and 1".to_string(),
        ));
    }
    if horizontal_dimension <= 0.0 || vertical_dimension <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Dimensions must be positive".to_string(),
        ));
    }
    Ok(efficiency * horizontal_dimension * vertical_dimension)
}

/// Calculate effective aperture for circular antenna
///
/// # Arguments
///
/// * `efficiency` - Antenna efficiency (0-1)
/// * `diameter` - Antenna diameter in meters
///
/// # Returns
///
/// Effective aperture in m²
pub fn calculate_effective_aperture_circ(efficiency: f64, diameter: f64) -> Result<f64> {
    if efficiency <= 0.0 || efficiency > 1.0 {
        return Err(RadarError::InvalidParameter(
            "Efficiency must be between 0 and 1".to_string(),
        ));
    }
    if diameter <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Diameter must be positive".to_string(),
        ));
    }
    Ok(efficiency * constants::PI * diameter * diameter * 0.25)
}

/// Calculate maximum detection range using the radar range equation
///
/// # Arguments
///
/// * `transmit_power` - Transmit power in watts
/// * `antenna_gain` - Antenna gain (linear, not dB)
/// * `wavelength` - Wavelength in meters
/// * `radar_cross_section` - Target radar cross-section in m²
/// * `min_detectable_signal` - Minimum detectable signal in watts
///
/// # Returns
///
/// Maximum range in meters
///
/// # Example
///
/// ```
/// use radar_range_equation::calculate_max_range;
/// let max_range = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
/// assert!(max_range > 0.0);
/// ```
pub fn calculate_max_range(
    transmit_power: f64,
    antenna_gain: f64,
    wavelength: f64,
    radar_cross_section: f64,
    min_detectable_signal: f64,
) -> Result<f64> {
    if transmit_power <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Transmit power must be positive".to_string(),
        ));
    }
    if antenna_gain <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Antenna gain must be positive".to_string(),
        ));
    }
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    if radar_cross_section <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Radar cross-section must be positive".to_string(),
        ));
    }
    if min_detectable_signal <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Min detectable signal must be positive".to_string(),
        ));
    }

    let numerator = transmit_power
        * antenna_gain.powi(2)
        * wavelength.powi(2)
        * radar_cross_section;
    let denominator = (4.0 * constants::PI).powi(3) * min_detectable_signal;
    Ok((numerator / denominator).powf(0.25))
}

/// Calculate received power at a given range
///
/// # Arguments
///
/// * `transmit_power` - Transmit power in watts
/// * `transmit_gain` - Transmit antenna gain (linear)
/// * `receive_gain` - Receive antenna gain (linear)
/// * `wavelength` - Wavelength in meters
/// * `radar_cross_section` - Target radar cross-section in m²
/// * `range` - Range to target in meters
///
/// # Returns
///
/// Received power in watts
pub fn calculate_received_power(
    transmit_power: f64,
    transmit_gain: f64,
    receive_gain: f64,
    wavelength: f64,
    radar_cross_section: f64,
    range: f64,
) -> Result<f64> {
    if transmit_power <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Transmit power must be positive".to_string(),
        ));
    }
    if transmit_gain <= 0.0 || receive_gain <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Antenna gains must be positive".to_string(),
        ));
    }
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    if radar_cross_section <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Radar cross-section must be positive".to_string(),
        ));
    }
    if range <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Range must be positive".to_string(),
        ));
    }

    let numerator = transmit_power
        * transmit_gain
        * receive_gain
        * wavelength.powi(2)
        * radar_cross_section;
    let denominator = (4.0 * constants::PI).powi(3) * range.powi(4);
    Ok(numerator / denominator)
}

/// Calculate beamwidth (Gaussian approximation)
///
/// # Arguments
///
/// * `wavelength` - Wavelength in meters
/// * `horizontal_dimension` - Horizontal antenna dimension in meters
///
/// # Returns
///
/// Beamwidth in radians
pub fn calculate_beamwidth(wavelength: f64, horizontal_dimension: f64) -> Result<f64> {
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    if horizontal_dimension <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Horizontal dimension must be positive".to_string(),
        ));
    }
    Ok((constants::BEAMWIDTH_COEFFICIENT * constants::PI / 180.0) * (wavelength / horizontal_dimension))
}

/// Calculate Doppler frequency shift
///
/// # Arguments
///
/// * `velocity` - Target velocity in m/s (negative for closing)
/// * `wavelength` - Wavelength in meters
///
/// # Returns
///
/// Doppler frequency shift in Hz
pub fn calculate_doppler_frequency(velocity: f64, wavelength: f64) -> Result<f64> {
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    Ok(-2.0 * velocity / wavelength)
}

/// Calculate velocity from Doppler shift
///
/// # Arguments
///
/// * `doppler_frequency` - Doppler frequency shift in Hz
/// * `wavelength` - Wavelength in meters
///
/// # Returns
///
/// Velocity in m/s (negative for closing)
pub fn calculate_velocity_from_doppler(doppler_frequency: f64, wavelength: f64) -> Result<f64> {
    if wavelength <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Wavelength must be positive".to_string(),
        ));
    }
    Ok(-wavelength * doppler_frequency / 2.0)
}

/// Calculate radar cross section of a sphere
///
/// # Arguments
///
/// * `radius` - Sphere radius in meters
///
/// # Returns
///
/// Radar cross-section in m²
pub fn calculate_sphere_rcs(radius: f64) -> Result<f64> {
    if radius <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Radius must be positive".to_string(),
        ));
    }
    Ok(constants::PI * radius * radius)
}

/// Calculate sphere radius from radar cross section
///
/// # Arguments
///
/// * `radar_cross_section` - Radar cross-section in m²
///
/// # Returns
///
/// Sphere radius in meters
pub fn calculate_sphere_radius(radar_cross_section: f64) -> Result<f64> {
    if radar_cross_section <= 0.0 {
        return Err(RadarError::InvalidParameter(
            "Radar cross-section must be positive".to_string(),
        ));
    }
    Ok((radar_cross_section / constants::PI).sqrt())
}

/// Unit conversion utilities
pub mod convert {
    use super::*;

    /// Convert linear value to dB
    pub fn linear_to_db(linear: f64) -> f64 {
        if linear <= 0.0 {
            f64::NEG_INFINITY
        } else {
            10.0 * linear.log10()
        }
    }

    /// Convert dB to linear value
    pub fn db_to_linear(db: f64) -> f64 {
        10.0_f64.powf(db / 10.0)
    }

    /// Convert radians to degrees
    pub fn rad_to_deg(radians: f64) -> f64 {
        radians * 180.0 / constants::PI
    }

    /// Convert degrees to radians
    pub fn deg_to_rad(degrees: f64) -> f64 {
        degrees * constants::PI / 180.0
    }

    /// Convert meters to nautical miles
    pub fn meters_to_nautical_miles(meters: f64) -> f64 {
        meters / 1852.0
    }

    /// Convert nautical miles to meters
    pub fn nautical_miles_to_meters(nautical_miles: f64) -> f64 {
        nautical_miles * 1852.0
    }

    /// Convert meters to miles
    pub fn meters_to_miles(meters: f64) -> f64 {
        meters / 1609.34
    }

    /// Convert miles to meters
    pub fn miles_to_meters(miles: f64) -> f64 {
        miles * 1609.34
    }

    /// Convert meters to kilometers
    pub fn meters_to_kilometers(meters: f64) -> f64 {
        meters / 1000.0
    }

    /// Convert kilometers to meters
    pub fn kilometers_to_meters(kilometers: f64) -> f64 {
        kilometers * 1000.0
    }

    /// Convert feet to meters
    pub fn feet_to_meters(feet: f64) -> f64 {
        feet * 0.3048
    }

    /// Convert meters to feet
    pub fn meters_to_feet(meters: f64) -> f64 {
        meters / 0.3048
    }

    /// Convert Hz to MHz
    pub fn hz_to_mhz(hz: f64) -> f64 {
        hz / 1e6
    }

    /// Convert MHz to Hz
    pub fn mhz_to_hz(mhz: f64) -> f64 {
        mhz * 1e6
    }

    /// Convert Hz to GHz
    pub fn hz_to_ghz(hz: f64) -> f64 {
        hz / 1e9
    }

    /// Convert GHz to Hz
    pub fn ghz_to_hz(ghz: f64) -> f64 {
        ghz * 1e9
    }

    /// Convert watts to kilowatts
    pub fn watts_to_kilowatts(watts: f64) -> f64 {
        watts / 1000.0
    }

    /// Convert kilowatts to watts
    pub fn kilowatts_to_watts(kilowatts: f64) -> f64 {
        kilowatts * 1000.0
    }

    /// Convert watts to milliwatts
    pub fn watts_to_milliwatts(watts: f64) -> f64 {
        watts * 1000.0
    }

    /// Convert milliwatts to watts
    pub fn milliwatts_to_watts(milliwatts: f64) -> f64 {
        milliwatts / 1000.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_wavelength() {
        let wavelength = calculate_wavelength(10e9).unwrap();
        assert!((wavelength - 0.0299792458).abs() < 1e-6);
    }

    #[test]
    fn test_calculate_frequency() {
        let frequency = calculate_frequency(0.03).unwrap();
        assert!((frequency - 9.99308193e9).abs() < 1e3);
    }

    #[test]
    fn test_calculate_max_range() {
        let max_range = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
        assert!(max_range > 0.0);
        // Verify the value is in a reasonable range
        assert!(max_range > 8000.0 && max_range < 9000.0);
    }

    #[test]
    fn test_calculate_received_power() {
        let power = calculate_received_power(1000.0, 1000.0, 1000.0, 0.03, 1.0, 10000.0).unwrap();
        assert!(power > 0.0);
    }

    #[test]
    fn test_calculate_doppler_frequency() {
        let doppler = calculate_doppler_frequency(-100.0, 0.03).unwrap();
        assert!((doppler - 6666.67).abs() < 0.01);
    }

    #[test]
    fn test_calculate_velocity_from_doppler() {
        let velocity = calculate_velocity_from_doppler(6666.67, 0.03).unwrap();
        assert!((velocity - (-100.0)).abs() < 0.01);
    }

    #[test]
    fn test_sphere_rcs() {
        let rcs = calculate_sphere_rcs(1.0).unwrap();
        assert!((rcs - constants::PI).abs() < 1e-6);
    }

    #[test]
    fn test_sphere_radius() {
        let radius = calculate_sphere_radius(constants::PI).unwrap();
        assert!((radius - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_linear_to_db() {
        assert!((convert::linear_to_db(10.0) - 10.0).abs() < 1e-6);
        assert!((convert::linear_to_db(100.0) - 20.0).abs() < 1e-6);
        assert!((convert::linear_to_db(1.0) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_db_to_linear() {
        assert!((convert::db_to_linear(10.0) - 10.0).abs() < 1e-6);
        assert!((convert::db_to_linear(20.0) - 100.0).abs() < 1e-6);
        assert!((convert::db_to_linear(0.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_angle_conversions() {
        assert!((convert::rad_to_deg(constants::PI) - 180.0).abs() < 1e-6);
        assert!((convert::deg_to_rad(180.0) - constants::PI).abs() < 1e-6);
    }

    #[test]
    fn test_distance_conversions() {
        assert!((convert::meters_to_nautical_miles(1852.0) - 1.0).abs() < 1e-6);
        assert!((convert::feet_to_meters(1.0) - 0.3048).abs() < 1e-6);
    }

    #[test]
    fn test_frequency_conversions() {
        assert!((convert::hz_to_mhz(1e6) - 1.0).abs() < 1e-6);
        assert!((convert::hz_to_ghz(1e9) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_invalid_parameters() {
        assert!(calculate_wavelength(-1.0).is_err());
        assert!(calculate_max_range(-1.0, 1000.0, 0.03, 1.0, 1e-13).is_err());
        assert!(calculate_sphere_rcs(-1.0).is_err());
    }
}
