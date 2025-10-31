# Radar Range Equation

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

```
         ┌─────────────────────────────────────────────┐
         │    P_t × G² × λ² × σ                       │
R_max = ⁴│  ──────────────────────────────             │
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
