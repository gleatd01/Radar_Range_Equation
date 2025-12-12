# Radar Range Equation

![Python Package Tests](https://img.shields.io/badge/tests-100%25%20passing-brightgreen)
![Equations](https://img.shields.io/badge/equations-72-blue)
![Solvers](https://img.shields.io/badge/solvers-76-blue)

A multi-language library providing radar range equation calculations for Python, Flutter/Dart, and Rust.

## Overview

This repository contains implementations of radar range equation calculators in multiple programming languages. The library calculates:

- **Maximum Detection Range**: Determines how far a radar can detect a target
- **Received Power**: Calculates the power received from a target at a given range

All implementations provide the same core functionality with language-specific APIs and conventions.

## Package Structure

```
Radar_Range_Equation/
├── python/          # Python package
├── flutter/         # Flutter/Dart package
├── rust/            # Rust crate
├── README.md        # This file
└── LICENSE          # MIT License
```

## Language-Specific Documentation

Each language implementation has its own README with detailed installation and usage instructions:

- [Python Package](python/README.md) - [PyPI](https://pypi.org/project/Radar-Range-Equation/)
- [Flutter/Dart Package](flutter/README.md)
- [Rust Crate](rust/README.md)

## Quick Start

### Python

```bash
cd python
pip install -e .
```

```python
from radar_range_equation import calculate_max_range

max_range = calculate_max_range(
    transmit_power=1000,
    antenna_gain=1000,
    wavelength=0.03,
    radar_cross_section=1.0,
    min_detectable_signal=1e-13
)
print(f"Maximum range: {max_range:.2f} meters")
```

### Flutter/Dart

```bash
cd flutter
dart pub get
```

```dart
import 'package:radar_range_equation/radar_range_equation.dart';

final maxRange = calculateMaxRange(
  transmitPower: 1000,
  antennaGain: 1000,
  wavelength: 0.03,
  radarCrossSection: 1.0,
  minDetectableSignal: 1e-13,
);
print('Maximum range: ${maxRange.toStringAsFixed(2)} meters');
```

### Rust

```bash
cd rust
cargo build
```

```rust
use radar_range_equation::calculate_max_range;

let max_range = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13)?;
println!("Maximum range: {:.2} meters", max_range);
```

## Features

- **Multi-language Support**: Python, Flutter/Dart, and Rust implementations
- **Consistent API**: Similar function signatures across all languages
- **Well-tested**: Comprehensive test suites for each implementation
- **Documentation**: Full API documentation and examples
- **Type Safety**: Strong typing where applicable (Python type hints, Dart, Rust)

## Radar Range Equation

The radar range equation relates the range of a radar to the characteristics of the transmitter, receiver, antenna, target, and environment. This library implements the basic radar equation:

### Maximum Range Equation

```math
R_{max} = \sqrt[4]{\frac{P_t \cdot G^2 \cdot \lambda^2 \cdot \sigma}{(4\pi)^3 \cdot S_{min}}}
```

Where:
- $R_{max}$ = Maximum detection range (meters)
- $P_t$ = Transmit power (watts)
- $G$ = Antenna gain (dimensionless)
- $\lambda$ = Wavelength (meters)
- $\sigma$ = Radar cross-section (m²)
- $S_{min}$ = Minimum detectable signal (watts)

### Wavelength Equation

```math
\lambda = \frac{c}{f}
```

Where:
- $\lambda$ = Wavelength (meters)
- $c$ = Speed of light (3×10⁸ m/s)
- $f$ = Frequency (Hz)

### Available Equations and Solvers

This package includes **72 equations** and **76 solver functions** covering:

- **Base/Common**: Fundamental radar equations (wavelength, gain, range)
- **Doppler CW Radar**: Doppler frequency shift and velocity calculations
- **CWFM Radar**: Continuous Wave FM radar equations
- **Pulsed Radar**: Range ambiguity, pulse repetition, and integration
- **Direction Finding**: Angle estimation and accuracy
- **Pulse Compression**: Range resolution and pulse compression ratio
- **Electronic Warfare**: Chaff, jamming, false targets, and countermeasures

For detailed documentation of all equations with LaTeX rendering, see:
- [Complete Equations Reference](EQUATIONS.md) - Full list of all 72 equations
- [Python Package README](python/README.md) - Python-specific usage and examples

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. When adding support for a new language:

1. Create a new directory for the language
2. Follow the existing structure (src, tests, examples)
3. Include comprehensive tests
4. Add language-specific README
5. Update this main README

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Future Expansion

The project structure is designed to easily accommodate additional language implementations. To add a new language:

1. Create a new directory at the root level
2. Implement the core functions following your language's conventions
3. Add tests and examples
4. Document installation and usage in a language-specific README
5. Update the main README with your language

## References

- Skolnik, M. I. (2008). Radar Handbook, Third Edition. McGraw-Hill.
- Richards, M. A. (2014). Fundamentals of Radar Signal Processing, Second Edition. McGraw-Hill.
