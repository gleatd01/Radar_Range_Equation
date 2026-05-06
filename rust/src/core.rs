use crate::constants;
use crate::validation::{ensure_positive, ensure_unit_interval};
use crate::Result;

/// Calculate wavelength from frequency.
pub fn calculate_wavelength(frequency: f64) -> Result<f64> {
    ensure_positive(frequency, "Frequency")?;
    Ok(constants::SPEED_OF_LIGHT / frequency)
}

/// Calculate frequency from wavelength.
pub fn calculate_frequency(wavelength: f64) -> Result<f64> {
    ensure_positive(wavelength, "Wavelength")?;
    Ok(constants::SPEED_OF_LIGHT / wavelength)
}

/// Calculate transmit antenna gain.
pub fn calculate_antenna_gain(effective_aperture: f64, wavelength: f64) -> Result<f64> {
    ensure_positive(effective_aperture, "Effective aperture")?;
    ensure_positive(wavelength, "Wavelength")?;
    Ok((4.0 * constants::PI * effective_aperture) / (wavelength * wavelength))
}

/// Calculate effective aperture for rectangular antenna.
pub fn calculate_effective_aperture_rect(
    efficiency: f64,
    horizontal_dimension: f64,
    vertical_dimension: f64,
) -> Result<f64> {
    ensure_unit_interval(efficiency, "Efficiency")?;
    ensure_positive(horizontal_dimension, "Horizontal dimension")?;
    ensure_positive(vertical_dimension, "Vertical dimension")?;
    Ok(efficiency * horizontal_dimension * vertical_dimension)
}

/// Calculate effective aperture for circular antenna.
pub fn calculate_effective_aperture_circ(efficiency: f64, diameter: f64) -> Result<f64> {
    ensure_unit_interval(efficiency, "Efficiency")?;
    ensure_positive(diameter, "Diameter")?;
    Ok(efficiency * constants::PI * diameter * diameter * 0.25)
}

/// Calculate maximum detection range using the radar range equation.
pub fn calculate_max_range(
    transmit_power: f64,
    antenna_gain: f64,
    wavelength: f64,
    radar_cross_section: f64,
    min_detectable_signal: f64,
) -> Result<f64> {
    ensure_positive(transmit_power, "Transmit power")?;
    ensure_positive(antenna_gain, "Antenna gain")?;
    ensure_positive(wavelength, "Wavelength")?;
    ensure_positive(radar_cross_section, "Radar cross-section")?;
    ensure_positive(min_detectable_signal, "Min detectable signal")?;

    let numerator = transmit_power
        * antenna_gain.powi(2)
        * wavelength.powi(2)
        * radar_cross_section;
    let denominator = (4.0 * constants::PI).powi(3) * min_detectable_signal;
    Ok((numerator / denominator).powf(0.25))
}

/// Calculate received power at a given range.
pub fn calculate_received_power(
    transmit_power: f64,
    transmit_gain: f64,
    receive_gain: f64,
    wavelength: f64,
    radar_cross_section: f64,
    range: f64,
) -> Result<f64> {
    ensure_positive(transmit_power, "Transmit power")?;
    ensure_positive(transmit_gain, "Transmit gain")?;
    ensure_positive(receive_gain, "Receive gain")?;
    ensure_positive(wavelength, "Wavelength")?;
    ensure_positive(radar_cross_section, "Radar cross-section")?;
    ensure_positive(range, "Range")?;

    let numerator = transmit_power
        * transmit_gain
        * receive_gain
        * wavelength.powi(2)
        * radar_cross_section;
    let denominator = (4.0 * constants::PI).powi(3) * range.powi(4);
    Ok(numerator / denominator)
}

/// Calculate beamwidth (Gaussian approximation).
pub fn calculate_beamwidth(wavelength: f64, horizontal_dimension: f64) -> Result<f64> {
    ensure_positive(wavelength, "Wavelength")?;
    ensure_positive(horizontal_dimension, "Horizontal dimension")?;
    Ok((constants::BEAMWIDTH_COEFFICIENT * constants::PI / 180.0)
        * (wavelength / horizontal_dimension))
}

/// Calculate Doppler frequency shift.
pub fn calculate_doppler_frequency(velocity: f64, wavelength: f64) -> Result<f64> {
    ensure_positive(wavelength, "Wavelength")?;
    Ok(-2.0 * velocity / wavelength)
}

/// Calculate velocity from Doppler shift.
pub fn calculate_velocity_from_doppler(doppler_frequency: f64, wavelength: f64) -> Result<f64> {
    ensure_positive(wavelength, "Wavelength")?;
    Ok(-wavelength * doppler_frequency / 2.0)
}

/// Calculate radar cross section of a sphere.
pub fn calculate_sphere_rcs(radius: f64) -> Result<f64> {
    ensure_positive(radius, "Radius")?;
    Ok(constants::PI * radius * radius)
}

/// Calculate sphere radius from radar cross section.
pub fn calculate_sphere_radius(radar_cross_section: f64) -> Result<f64> {
    ensure_positive(radar_cross_section, "Radar cross-section")?;
    Ok((radar_cross_section / constants::PI).sqrt())
}
