# Changelog

All notable changes to the Radar Range Equation project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- CONTRIBUTING.md guide for contributors
- CHANGELOG.md to track version history
- Development requirements files for easier setup

## [2025.10.27] - 2025-10-27

### Added
- Support for theta_B calculations
- Calculations for spheres, RCS and radius from RCS
- Plotting capabilities for visualizing radar signals including:
  - Pulsed radar signals
  - CW Doppler signals
  - CWFM (Continuous Wave Frequency Modulated) signals
  - Pulse compression signals
  - Range profiles
  - Doppler spectra
- Comprehensive plotting API via `plot` module
- Example files demonstrating plotting capabilities

### Python Package
- Added `plot.py` module with visualization functions
- Added `example_plotting.py` with comprehensive examples
- Added `test_plotting.py` for testing plot functionality
- Added `PLOTTING_SUMMARY.md` documentation

## Previous Versions

### Core Features (Initial Release)
- Multi-language support: Python, Flutter/Dart, and Rust implementations
- Radar range equation calculations:
  - Maximum detection range
  - Received power at a given range
  - Wavelength from frequency
  - Antenna gain calculations
  - Doppler frequency shifts
  - Radar cross-section calculations
- Comprehensive equation library (72 equations)
- Extensive solver functions (76 solvers) covering:
  - Base/Common radar equations
  - Doppler CW Radar
  - CWFM Radar
  - Pulsed Radar (range ambiguity, pulse repetition, integration)
  - Direction Finding (angle estimation and accuracy)
  - Pulse Compression (range resolution, compression ratio)
  - Electronic Warfare (chaff, jamming, false targets, countermeasures)
- Unit conversion utilities for:
  - Power (linear to dB, dB to linear)
  - Distance (meters, kilometers, nautical miles, feet)
  - Frequency (Hz, MHz, GHz)
  - Angles (radians, degrees)
- Full test suites for all implementations
- Comprehensive documentation and examples

### Python Package
- SymPy-based symbolic equation system
- NumPy/SciPy numeric solvers
- Type hints for better IDE support
- Variable namespace system (`vars`)
- Analysis helpers for pulse parsing and integration gains
- Example scripts for various use cases

### Flutter/Dart Package
- Pure Dart implementation
- Comprehensive unit conversion utilities
- Full API documentation
- Example applications
- Test suite

### Rust Crate
- Pure Rust implementation with no external dependencies
- Custom error handling with `RadarError` type
- Inline optimizations for performance
- Comprehensive documentation
- Example programs
- Test suite

---

## Version History Notes

- Versions are dated following the format YYYY.MM.DD for clarity
- Python package available on PyPI as `radar-range-equation`
- Flutter/Dart package available as local package
- Rust crate available as local crate (not yet published to crates.io)

[Unreleased]: https://github.com/gleatd01/Radar_Range_Equation/compare/v2025.10.27...HEAD
[2025.10.27]: https://github.com/gleatd01/Radar_Range_Equation/releases/tag/v2025.10.27
