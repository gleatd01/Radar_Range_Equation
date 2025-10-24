# Flutter/Dart Package

This directory contains the Flutter/Dart implementation of the Radar Range Equation tools.

## Installation

Add this to your package's `pubspec.yaml` file:

```yaml
dependencies:
  radar_range_equation:
    path: ../flutter  # Adjust path as needed
```

Or for development:

```yaml
dependencies:
  radar_range_equation:
    git:
      url: https://github.com/gleatd01/Radar_Range_Equation.git
      path: flutter
```

Then run:

```bash
flutter pub get
```

## Usage

```dart
import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  // Calculate maximum detection range
  final maxRange = calculateMaxRange(
    transmitPower: 1000,           // 1 kW
    antennaGain: 1000,              // 30 dBi (linear scale)
    wavelength: 0.03,               // 10 GHz
    radarCrossSection: 1.0,         // 1 m²
    minDetectableSignal: 1e-13,     // -100 dBm
  );

  print('Maximum range: ${maxRange.toStringAsFixed(2)} meters');

  // Calculate received power at a specific range
  final receivedPower = calculateReceivedPower(
    transmitPower: 1000,
    antennaGain: 1000,
    wavelength: 0.03,
    radarCrossSection: 1.0,
    rangeM: 10000,  // 10 km
  );

  print('Received power: ${receivedPower.toStringAsExponential(2)} Watts');
}
```

## Testing

```bash
cd flutter
dart test
```

## Linting

```bash
cd flutter
dart analyze
```

## Project Structure

```
flutter/
├── lib/
│   ├── radar_range_equation.dart
│   └── src/
│       └── radar_range.dart
├── test/
│   └── radar_range_test.dart
├── example/
│   └── basic_example.dart
├── pubspec.yaml
└── README.md
```
