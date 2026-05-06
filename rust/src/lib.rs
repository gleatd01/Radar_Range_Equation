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

pub mod constants;
mod core;
pub mod convert;
mod errors;
mod validation;

pub use core::{
    calculate_antenna_gain, calculate_beamwidth, calculate_doppler_frequency,
    calculate_effective_aperture_circ, calculate_effective_aperture_rect,
    calculate_frequency, calculate_max_range, calculate_received_power,
    calculate_sphere_radius, calculate_sphere_rcs, calculate_velocity_from_doppler,
    calculate_wavelength,
};

pub use errors::{RadarError, Result};
