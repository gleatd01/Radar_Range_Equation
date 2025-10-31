# Radar Range Equation - Rust Crate

A Rust library for calculating radar range equations. This crate provides functions for computing maximum detection range, received power, Doppler frequency shifts, and other radar-related calculations.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
radar_range_equation = { path = "../rust" }
```

Or, if published to crates.io:

```toml
[dependencies]
radar_range_equation = "0.0.1"
```

## Features

- **Radar Range Calculations**: Calculate maximum detection range using the radar range equation
- **Power Calculations**: Compute received power at a given range
- **Frequency/Wavelength**: Convert between frequency and wavelength
- **Antenna Calculations**: Calculate antenna gain and effective aperture
- **Doppler Calculations**: Compute Doppler frequency shifts and velocities
- **Unit Conversions**: Comprehensive unit conversion utilities
- **Cross-section Calculations**: Calculate radar cross-section for spheres
- **Error Handling**: Robust error handling with custom error types
- **No Dependencies**: Pure Rust implementation with no external dependencies

## Quick Start

```rust
use radar_range_equation::*;

fn main() {
    // Calculate wavelength from frequency
    let frequency = 10e9; // 10 GHz
    let wavelength = calculate_wavelength(frequency).unwrap();
    println!("Wavelength: {:.4} m", wavelength);

    // Calculate maximum radar range
    let max_range = calculate_max_range(
        1000.0,  // transmit power (W)
        1000.0,  // antenna gain (linear)
        wavelength,
        1.0,     // radar cross-section (m²)
        1e-13,   // min detectable signal (W)
    ).unwrap();
    println!("Maximum range: {:.2} meters", max_range);

    // Calculate Doppler frequency
    let velocity = -100.0; // 100 m/s closing
    let doppler_freq = calculate_doppler_frequency(velocity, wavelength).unwrap();
    println!("Doppler shift: {:.3} MHz", convert::hz_to_mhz(doppler_freq));
}
```

## API Reference

### Core Functions

#### `calculate_wavelength(frequency: f64) -> Result<f64>`
Calculate wavelength from frequency.
- **Parameters**: `frequency` - Frequency in Hz
- **Returns**: `Result<f64>` - Wavelength in meters

#### `calculate_max_range(...) -> Result<f64>`
Calculate maximum detection range using the radar range equation.
- **Parameters**:
  - `transmit_power: f64` - Transmit power in watts
  - `antenna_gain: f64` - Antenna gain (linear, not dB)
  - `wavelength: f64` - Wavelength in meters
  - `radar_cross_section: f64` - Target RCS in m²
  - `min_detectable_signal: f64` - Minimum detectable signal in watts
- **Returns**: `Result<f64>` - Maximum range in meters

#### `calculate_received_power(...) -> Result<f64>`
Calculate received power at a given range.
- **Parameters**:
  - `transmit_power: f64` - Transmit power in watts
  - `transmit_gain: f64` - Transmit antenna gain (linear)
  - `receive_gain: f64` - Receive antenna gain (linear)
  - `wavelength: f64` - Wavelength in meters
  - `radar_cross_section: f64` - Target RCS in m²
  - `range: f64` - Range to target in meters
- **Returns**: `Result<f64>` - Received power in watts

#### `calculate_doppler_frequency(velocity: f64, wavelength: f64) -> Result<f64>`
Calculate Doppler frequency shift.
- **Parameters**:
  - `velocity: f64` - Target velocity in m/s (negative for closing)
  - `wavelength: f64` - Wavelength in meters
- **Returns**: `Result<f64>` - Doppler frequency shift in Hz

#### `calculate_antenna_gain(effective_aperture: f64, wavelength: f64) -> Result<f64>`
Calculate transmit antenna gain.
- **Parameters**:
  - `effective_aperture: f64` - Effective aperture in m²
  - `wavelength: f64` - Wavelength in meters
- **Returns**: `Result<f64>` - Antenna gain (linear, not dB)

### Unit Conversions

The `convert` module provides various unit conversion utilities:

```rust
use radar_range_equation::convert;

// Power conversions
convert::linear_to_db(100.0);        // Convert to dB
convert::db_to_linear(20.0);         // Convert from dB

// Distance conversions
convert::meters_to_kilometers(1000.0);
convert::meters_to_nautical_miles(1852.0);
convert::feet_to_meters(100.0);

// Frequency conversions
convert::hz_to_mhz(1e6);
convert::hz_to_ghz(1e9);

// Angle conversions
convert::rad_to_deg(std::f64::consts::PI);
convert::deg_to_rad(180.0);
```

## Examples

Run the basic example with:

```bash
cargo run --example basic
```

See the [examples](examples/) directory for more comprehensive examples.

## Testing

Run the tests with:

```bash
cargo test
```

Run tests with output:

```bash
cargo test -- --nocapture
```

## Radar Range Equation

The radar range equation relates the range of a radar to the characteristics of the transmitter, receiver, antenna, target, and environment:

```
         ┌─────────────────────────────────────────────┐
         │    P_t × G² × λ² × σ                       │
R_max = ⁴│  ──────────────────────────────────         │
         │    (4π)³ × P_min                           │
         └─────────────────────────────────────────────┘
```

Where:
- `R_max` = Maximum range
- `P_t` = Transmit power
- `G` = Antenna gain
- `λ` = Wavelength
- `σ` = Radar cross-section
- `P_min` = Minimum detectable signal

## Constants

The crate provides common physical constants through the `constants` module:

- `constants::SPEED_OF_LIGHT` - Speed of light (299,792,458 m/s)
- `constants::BOLTZMANN` - Boltzmann constant (1.380649×10⁻²³ J/K)
- `constants::PI` - Pi (π)

## Error Handling

The crate uses a custom `RadarError` type for error handling:

```rust
use radar_range_equation::*;

match calculate_wavelength(-1.0) {
    Ok(wavelength) => println!("Wavelength: {}", wavelength),
    Err(RadarError::InvalidParameter(msg)) => println!("Error: {}", msg),
    Err(e) => println!("Error: {}", e),
}
```

Error types:
- `RadarError::InvalidParameter` - Invalid input parameter (e.g., negative value)
- `RadarError::MathError` - Mathematical error (e.g., division by zero)

## Performance

This crate is optimized for performance with:
- Zero-cost abstractions
- Inline functions where appropriate
- No heap allocations in core calculations
- No external dependencies

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## References

- Skolnik, M. I. (2008). Radar Handbook, Third Edition. McGraw-Hill.
- Richards, M. A. (2014). Fundamentals of Radar Signal Processing, Second Edition. McGraw-Hill.
