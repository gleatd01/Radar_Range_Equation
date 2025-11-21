# Radar Range Equation - Node.js/TypeScript

A Node.js/TypeScript library providing radar range equation calculations for determining maximum detection range and received power for radar systems.

## Installation

```bash
npm install radar-range-equation
```

## Features

- **Maximum Range Calculation**: Determine how far a radar can detect a target
- **Received Power Calculation**: Calculate the power received from a target at a given range
- **Wavelength and Frequency Conversions**: Convert between frequency and wavelength
- **Antenna Calculations**: Calculate antenna gain and effective aperture
- **Unit Conversions**: Comprehensive utilities for converting between units
- **TypeScript Support**: Full TypeScript definitions included
- **Zero Dependencies**: No runtime dependencies
- **Well Tested**: Comprehensive test suite with >95% coverage

## Quick Start

```typescript
import { calculateMaxRange, convert } from 'radar-range-equation';

// Calculate maximum detection range
const maxRange = calculateMaxRange(
  1000,    // transmitPower (W)
  1000,    // antennaGain
  0.03,    // wavelength (m)
  1.0,     // radarCrossSection (m²)
  1e-13    // minDetectableSignal (W)
);

console.log(`Maximum range: ${maxRange.toFixed(2)} meters`);
console.log(`              ${convert.metersToKm(maxRange).toFixed(2)} km`);
```

## API Reference

### Core Functions

#### `calculateMaxRange(transmitPower, antennaGain, wavelength, radarCrossSection, minDetectableSignal)`

Calculates the maximum detection range of a radar system using the radar range equation.

**Parameters:**
- `transmitPower` (number): Transmit power in watts
- `antennaGain` (number): Antenna gain (dimensionless)
- `wavelength` (number): Wavelength in meters
- `radarCrossSection` (number): Radar cross-section in square meters
- `minDetectableSignal` (number): Minimum detectable signal in watts

**Returns:** Maximum detection range in meters

**Example:**
```typescript
const maxRange = calculateMaxRange(1000, 1000, 0.03, 1.0, 1e-13);
// Returns: ~46148.04 meters
```

#### `calculateReceivedPower(transmitPower, antennaGain, wavelength, radarCrossSection, range)`

Calculates the power received by a radar from a target at a specified range.

**Parameters:**
- `transmitPower` (number): Transmit power in watts
- `antennaGain` (number): Antenna gain (dimensionless)
- `wavelength` (number): Wavelength in meters
- `radarCrossSection` (number): Radar cross-section in square meters
- `range` (number): Range to target in meters

**Returns:** Received power in watts

**Example:**
```typescript
const receivedPower = calculateReceivedPower(1000, 1000, 0.03, 1.0, 10000);
```

#### `calculateWavelength(frequency, speedOfLight?)`

Calculates wavelength from frequency.

**Parameters:**
- `frequency` (number): Frequency in Hz
- `speedOfLight` (number, optional): Speed of light in m/s (default: 299792458)

**Returns:** Wavelength in meters

**Example:**
```typescript
const wavelength = calculateWavelength(10e9); // 10 GHz
// Returns: ~0.02998 meters
```

#### `calculateAntennaGain(effectiveAperture, wavelength)`

Calculates antenna gain from effective aperture.

**Parameters:**
- `effectiveAperture` (number): Effective aperture in square meters
- `wavelength` (number): Wavelength in meters

**Returns:** Antenna gain (dimensionless)

**Example:**
```typescript
const gain = calculateAntennaGain(1.0, 0.03);
// Returns: ~13962.63
```

#### `calculateEffectiveApertureRect(efficiency, horizontalDimension, verticalDimension)`

Calculates effective aperture for a rectangular antenna.

**Parameters:**
- `efficiency` (number): Antenna efficiency (0 to 1)
- `horizontalDimension` (number): Horizontal dimension in meters
- `verticalDimension` (number): Vertical dimension in meters

**Returns:** Effective aperture in square meters

**Example:**
```typescript
const aperture = calculateEffectiveApertureRect(0.6, 2.0, 1.5);
// Returns: 1.8 m²
```

#### `calculateEffectiveApertureCirc(efficiency, diameter)`

Calculates effective aperture for a circular antenna.

**Parameters:**
- `efficiency` (number): Antenna efficiency (0 to 1)
- `diameter` (number): Antenna diameter in meters

**Returns:** Effective aperture in square meters

**Example:**
```typescript
const aperture = calculateEffectiveApertureCirc(0.6, 2.0);
// Returns: ~1.885 m²
```

### Conversion Utilities

The `convert` object provides various unit conversion functions:

#### Power/Signal Conversions
- `linearToDb(value)` - Convert linear value to decibels
- `dbToLinear(valueDb)` - Convert decibels to linear value

#### Angle Conversions
- `radToDeg(radians)` - Convert radians to degrees
- `degToRad(degrees)` - Convert degrees to radians

#### Distance Conversions
- `metersToKm(meters)` - Convert meters to kilometers
- `kmToMeters(km)` - Convert kilometers to meters
- `metersToNauticalMiles(meters)` - Convert meters to nautical miles
- `nauticalMilesToMeters(nauticalMiles)` - Convert nautical miles to meters
- `metersToMiles(meters)` - Convert meters to miles
- `milesToMeters(miles)` - Convert miles to meters

#### Frequency Conversions
- `hzToGhz(hz)` - Convert hertz to gigahertz
- `ghzToHz(ghz)` - Convert gigahertz to hertz

#### Power Conversions
- `wattsToKw(watts)` - Convert watts to kilowatts
- `kwToWatts(kw)` - Convert kilowatts to watts

**Example:**
```typescript
import { convert } from 'radar-range-equation';

const rangeKm = convert.metersToKm(46148);
console.log(`${rangeKm} km`); // 45.403 km

const powerDb = convert.linearToDb(1000);
console.log(`${powerDb} dB`); // 30 dB

const angleRad = convert.degToRad(90);
console.log(`${angleRad} rad`); // 1.5708 rad
```

### Physical Constants

The `constants` object provides commonly used physical constants:

- `SPEED_OF_LIGHT` - Speed of light in m/s (299792458)
- `BOLTZMANN` - Boltzmann constant in J/K (1.380649e-23)
- `PI` - Mathematical constant π

**Example:**
```typescript
import { constants } from 'radar-range-equation';

console.log(`c = ${constants.SPEED_OF_LIGHT} m/s`);
```

## Usage Examples

### Example 1: Basic Range Calculation

```typescript
import { calculateMaxRange, calculateWavelength, convert } from 'radar-range-equation';

// System parameters
const frequency = 10e9; // 10 GHz
const wavelength = calculateWavelength(frequency);
const transmitPower = 1000; // 1 kW
const antennaGain = 1000;
const radarCrossSection = 1.0; // 1 m²
const minDetectableSignal = 1e-13; // W

// Calculate maximum range
const maxRange = calculateMaxRange(
  transmitPower,
  antennaGain,
  wavelength,
  radarCrossSection,
  minDetectableSignal
);

console.log(`Maximum Detection Range: ${maxRange.toFixed(2)} m`);
console.log(`                         ${convert.metersToKm(maxRange).toFixed(2)} km`);
```

### Example 2: Received Power Analysis

```typescript
import { calculateReceivedPower, convert } from 'radar-range-equation';

const ranges = [5000, 10000, 20000, 40000]; // meters

console.log('Range (km) | Received Power (dBW)');
console.log('-----------|--------------------');

ranges.forEach(range => {
  const power = calculateReceivedPower(1000, 1000, 0.03, 1.0, range);
  const powerDb = convert.linearToDb(power);
  console.log(`${convert.metersToKm(range).toFixed(1).padEnd(10)} | ${powerDb.toFixed(2)}`);
});
```

### Example 3: Antenna Design

```typescript
import {
  calculateEffectiveApertureCirc,
  calculateAntennaGain,
  calculateWavelength,
  convert
} from 'radar-range-equation';

// Design parameters
const frequency = 35e9; // 35 GHz (Ka-band)
const wavelength = calculateWavelength(frequency);
const efficiency = 0.65; // 65% efficiency
const diameter = 0.5; // 0.5 meter dish

// Calculate antenna characteristics
const effectiveAperture = calculateEffectiveApertureCirc(efficiency, diameter);
const gain = calculateAntennaGain(effectiveAperture, wavelength);
const gainDb = convert.linearToDb(gain);

console.log(`Antenna Design for ${convert.hzToGhz(frequency)} GHz:`);
console.log(`  Diameter: ${diameter} m`);
console.log(`  Efficiency: ${(efficiency * 100).toFixed(0)}%`);
console.log(`  Effective Aperture: ${effectiveAperture.toFixed(4)} m²`);
console.log(`  Gain: ${gain.toFixed(2)} (${gainDb.toFixed(2)} dBi)`);
```

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

## Development

### Building

```bash
npm install
npm run build
```

### Testing

```bash
npm test
```

### Running Examples

```bash
npx ts-node examples/basic-usage.ts
```

## Related Packages

This is part of a multi-language radar range equation library:

- [Python Package](https://pypi.org/project/Radar-Range-Equation/) - Python implementation
- [Rust Crate](../rust/) - Rust implementation
- [Flutter/Dart Package](../flutter/) - Flutter/Dart implementation

## License

MIT License - see the [LICENSE](../LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## References

- Skolnik, M. I. (2008). *Radar Handbook, Third Edition*. McGraw-Hill.
- Richards, M. A. (2014). *Fundamentals of Radar Signal Processing, Second Edition*. McGraw-Hill.
