use crate::constants;

/// Convert linear value to dB.
pub fn linear_to_db(linear: f64) -> f64 {
    if linear <= 0.0 {
        f64::NEG_INFINITY
    } else {
        10.0 * linear.log10()
    }
}

/// Convert dB to linear value.
pub fn db_to_linear(db: f64) -> f64 {
    10.0_f64.powf(db / 10.0)
}

/// Convert radians to degrees.
pub fn rad_to_deg(radians: f64) -> f64 {
    radians * 180.0 / constants::PI
}

/// Convert degrees to radians.
pub fn deg_to_rad(degrees: f64) -> f64 {
    degrees * constants::PI / 180.0
}

/// Convert meters to nautical miles.
pub fn meters_to_nautical_miles(meters: f64) -> f64 {
    meters / 1852.0
}

/// Convert nautical miles to meters.
pub fn nautical_miles_to_meters(nautical_miles: f64) -> f64 {
    nautical_miles * 1852.0
}

/// Convert meters to miles.
pub fn meters_to_miles(meters: f64) -> f64 {
    meters / 1609.34
}

/// Convert miles to meters.
pub fn miles_to_meters(miles: f64) -> f64 {
    miles * 1609.34
}

/// Convert meters to kilometers.
pub fn meters_to_kilometers(meters: f64) -> f64 {
    meters / 1000.0
}

/// Convert kilometers to meters.
pub fn kilometers_to_meters(kilometers: f64) -> f64 {
    kilometers * 1000.0
}

/// Convert feet to meters.
pub fn feet_to_meters(feet: f64) -> f64 {
    feet * 0.3048
}

/// Convert meters to feet.
pub fn meters_to_feet(meters: f64) -> f64 {
    meters / 0.3048
}

/// Convert Hz to MHz.
pub fn hz_to_mhz(hz: f64) -> f64 {
    hz / 1e6
}

/// Convert MHz to Hz.
pub fn mhz_to_hz(mhz: f64) -> f64 {
    mhz * 1e6
}

/// Convert Hz to GHz.
pub fn hz_to_ghz(hz: f64) -> f64 {
    hz / 1e9
}

/// Convert GHz to Hz.
pub fn ghz_to_hz(ghz: f64) -> f64 {
    ghz * 1e9
}

/// Convert watts to kilowatts.
pub fn watts_to_kilowatts(watts: f64) -> f64 {
    watts / 1000.0
}

/// Convert kilowatts to watts.
pub fn kilowatts_to_watts(kilowatts: f64) -> f64 {
    kilowatts * 1000.0
}

/// Convert watts to milliwatts.
pub fn watts_to_milliwatts(watts: f64) -> f64 {
    watts * 1000.0
}

/// Convert milliwatts to watts.
pub fn milliwatts_to_watts(milliwatts: f64) -> f64 {
    milliwatts / 1000.0
}
