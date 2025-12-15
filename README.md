# Radar Range Equation

[![Python Package Tests](https://img.shields.io/badge/tests-100%25%20passing-brightgreen)](https://github.com/gleatd01/Radar_Range_Equation/actions)
[![Equations](https://img.shields.io/badge/equations-72-blue)](EQUATIONS.md)
[![Solvers](https://img.shields.io/badge/solvers-76-blue)](SOLVERS.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![PyPI](https://img.shields.io/badge/PyPI-radar--range--equation-blue)](https://pypi.org/project/radar-range-equation/)

A multi-language library providing radar range equation calculations for Python, Flutter/Dart, and Rust.

**📖 [Quickstart Guide](QUICKSTART.md)** | **📚 [Documentation](python/README.md)** | **🤝 [Contributing](CONTRIBUTING.md)** | **📝 [Changelog](CHANGELOG.md)** | **📋 [Code of Conduct](CODE_OF_CONDUCT.md)**

## Overview

This repository contains implementations of radar range equation calculators in multiple programming languages. The library calculates:

- **Maximum Detection Range**: Determines how far a radar can detect a target
- **Received Power**: Calculates the power received from a target at a given range

All implementations provide the same core functionality with language-specific APIs and conventions.

## 🚀 Installation

Choose your preferred language:

### Python
```bash
pip install radar-range-equation
```

### Flutter/Dart
```yaml
# Add to pubspec.yaml
dependencies:
  radar_range_equation:
    git:
      url: https://github.com/gleatd01/Radar_Range_Equation.git
      path: flutter
```

### Rust
```toml
# Add to Cargo.toml
[dependencies]
radar_range_equation = { git = "https://github.com/gleatd01/Radar_Range_Equation.git", path = "rust" }
```

For detailed installation instructions, see the [Quickstart Guide](QUICKSTART.md).

## 📦 Package Structure

```
Radar_Range_Equation/
├── python/              # Python package (PyPI: radar-range-equation)
├── flutter/             # Flutter/Dart package
├── rust/                # Rust crate
├── .github/             # CI/CD workflows and issue templates
├── README.md            # This file
├── QUICKSTART.md        # Quick start guide
├── CONTRIBUTING.md      # Contribution guidelines
├── CHANGELOG.md         # Version history
├── EQUATIONS.md         # Complete equations reference
├── SOLVERS.md           # Complete solver reference
└── LICENSE              # MIT License
```

## 📚 Documentation

Each language implementation has its own comprehensive documentation:

- **Python**: [README](python/README.md) • [PyPI Package](https://pypi.org/project/radar-range-equation/)
- **Flutter/Dart**: [README](flutter/README.md)
- **Rust**: [README](rust/README.md)

## ⚡ Quick Start

### Python

```bash
pip install radar-range-equation
```

```python
import radar_range_equation as RRE

# Set variables
RRE.vars.P_t = 1000         # Transmit power (W)
RRE.vars.G_t = 1000         # Antenna gain
RRE.vars.G_r = 1000
RRE.vars.wavelength = 0.03  # 30 cm wavelength
RRE.vars.sigma = 1.0        # Radar cross-section (m²)
RRE.vars.S_min = 1e-13      # Minimum detectable signal (W)

# Calculate maximum range
max_range = RRE.solve.R_max()
print(f"Maximum range: {max_range/1000:.2f} km")

# Visualize radar signal
RRE.plot.pulsed_radar_signal()
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

For detailed documentation, see:
- [Complete Equations Reference](EQUATIONS.md) - Full list of all 72 equations with LaTeX
- [Complete Solver Reference](SOLVERS.md) - Full list of all 76 solver functions
- [Python Package README](python/README.md) - Python-specific usage and examples

## 🤝 Contributing

We welcome contributions! Whether you're fixing bugs, adding features, improving documentation, or implementing support for a new language, your help is appreciated.

**Getting Started:**
1. Read the [Contributing Guide](CONTRIBUTING.md) and [Code of Conduct](CODE_OF_CONDUCT.md)
2. Run the development setup script: `./setup-dev.sh`
3. Check [open issues](https://github.com/gleatd01/Radar_Range_Equation/issues)
4. Fork the repository and create a feature branch
5. Submit a Pull Request

**Adding a New Language:**
1. Create a new directory for the language
2. Follow the existing structure (src, tests, examples)
3. Include comprehensive tests
4. Add language-specific README
5. Update this main README

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## 📋 Version History

See [CHANGELOG.md](CHANGELOG.md) for detailed version history and release notes.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🎯 Roadmap

Future plans for the project:

- [ ] Publish Rust crate to crates.io
- [ ] Publish Flutter package to pub.dev
- [ ] Add more radar equation variants
- [ ] Expand plotting capabilities
- [ ] Add interactive web calculator
- [ ] Additional language implementations (JavaScript, Java, C++)

## 🆘 Getting Help

- **Documentation**: Start with [QUICKSTART.md](QUICKSTART.md)
- **Examples**: Check language-specific example files
- **Issues**: Report bugs or ask questions on [GitHub Issues](https://github.com/gleatd01/Radar_Range_Equation/issues)
- **Discussions**: Use GitHub Discussions for general questions

## 📚 References

- Skolnik, M. I. (2008). *Radar Handbook, Third Edition*. McGraw-Hill.
- Richards, M. A. (2014). *Fundamentals of Radar Signal Processing, Second Edition*. McGraw-Hill.

## ⭐ Star History

If you find this project useful, please consider giving it a star! It helps others discover the project.

---

Made with ❤️ for the radar engineering community
