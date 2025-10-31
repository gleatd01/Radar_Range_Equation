# Radar Range Equation - Flutter/Dart Package

A Dart package for calculating radar range equations. This package provides functions for computing maximum detection range, received power, Doppler frequency shifts, and other radar-related calculations.

## Installation

Add this to your package's `pubspec.yaml` file:

```yaml
dependencies:
  radar_range_equation:
    path: ../flutter
```

Then run:

```bash
dart pub get
```

## Features

- **Radar Range Calculations**: Calculate maximum detection range using the radar range equation
- **Power Calculations**: Compute received power at a given range
- **Frequency/Wavelength**: Convert between frequency and wavelength
- **Antenna Calculations**: Calculate antenna gain and effective aperture
- **Doppler Calculations**: Compute Doppler frequency shifts and velocities
- **Unit Conversions**: Comprehensive unit conversion utilities
- **Cross-section Calculations**: Calculate radar cross-section for spheres

## Quick Start

```dart
import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  // Calculate wavelength from frequency
  final frequency = 10e9; // 10 GHz
  final wavelength = calculateWavelength(frequency);
  print('Wavelength: ${wavelength * 100} cm');

  // Calculate maximum radar range
  final maxRange = calculateMaxRange(
    transmitPower: 1000.0,      // 1 kW
    antennaGain: 1000.0,        // Linear gain
    wavelength: wavelength,
    radarCrossSection: 1.0,     // 1 m²
    minDetectableSignal: 1e-13, // 0.1 pW
  );
  print('Maximum range: ${maxRange.toStringAsFixed(2)} meters');

  // Calculate Doppler frequency
  final velocity = -100.0; // 100 m/s closing
  final dopplerFreq = calculateDopplerFrequency(velocity, wavelength);
  print('Doppler shift: ${RadarConvert.hzToMhz(dopplerFreq).toStringAsFixed(3)} MHz');
}
```

## API Reference

### Core Functions

#### `calculateWavelength(double frequency)`
Calculate wavelength from frequency.
- **Parameters**: `frequency` - Frequency in Hz
- **Returns**: Wavelength in meters

#### `calculateMaxRange(...)`
Calculate maximum detection range using the radar range equation.
- **Parameters**:
  - `transmitPower` - Transmit power in watts
  - `antennaGain` - Antenna gain (linear, not dB)
  - `wavelength` - Wavelength in meters
  - `radarCrossSection` - Target RCS in m²
  - `minDetectableSignal` - Minimum detectable signal in watts
- **Returns**: Maximum range in meters

#### `calculateReceivedPower(...)`
Calculate received power at a given range.
- **Parameters**:
  - `transmitPower` - Transmit power in watts
  - `transmitGain` - Transmit antenna gain (linear)
  - `receiveGain` - Receive antenna gain (linear)
  - `wavelength` - Wavelength in meters
  - `radarCrossSection` - Target RCS in m²
  - `range` - Range to target in meters
- **Returns**: Received power in watts

#### `calculateDopplerFrequency(double velocity, double wavelength)`
Calculate Doppler frequency shift.
- **Parameters**:
  - `velocity` - Target velocity in m/s (negative for closing)
  - `wavelength` - Wavelength in meters
- **Returns**: Doppler frequency shift in Hz

#### `calculateAntennaGain(double effectiveAperture, double wavelength)`
Calculate transmit antenna gain.
- **Parameters**:
  - `effectiveAperture` - Effective aperture in m²
  - `wavelength` - Wavelength in meters
- **Returns**: Antenna gain (linear, not dB)

### Unit Conversions

The `RadarConvert` class provides various unit conversion utilities:

```dart
// Power conversions
RadarConvert.linearToDb(100)        // Convert to dB
RadarConvert.dbToLinear(20)          // Convert from dB

// Distance conversions
RadarConvert.metersToKilometers(1000)
RadarConvert.metersToNauticalMiles(1852)
RadarConvert.feetToMeters(100)

// Frequency conversions
RadarConvert.hzToMhz(1e6)
RadarConvert.hzToGhz(1e9)

// Angle conversions
RadarConvert.radToDeg(3.14159)
RadarConvert.degToRad(180)
```

## Examples

See the [example](example/example.dart) directory for more comprehensive examples.

## Testing

Run the tests with:

```bash
dart test
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

The package provides common physical constants through the `RadarConstants` class:

- `RadarConstants.speedOfLight` - Speed of light (299,792,458 m/s)
- `RadarConstants.boltzmann` - Boltzmann constant (1.380649×10⁻²³ J/K)
- `RadarConstants.pi` - Pi (π)

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## References

- Skolnik, M. I. (2008). Radar Handbook, Third Edition. McGraw-Hill.
- Richards, M. A. (2014). Fundamentals of Radar Signal Processing, Second Edition. McGraw-Hill.
