import 'dart:math';

/// Calculate the maximum detection range of a radar system.
///
/// Uses the radar range equation to determine the maximum range at which
/// a target can be detected.
///
/// Parameters:
/// - [transmitPower]: Transmitted power in Watts
/// - [antennaGain]: Antenna gain (dimensionless, linear scale)
/// - [wavelength]: Radar wavelength in meters
/// - [radarCrossSection]: Target radar cross section in square meters
/// - [minDetectableSignal]: Minimum detectable signal power in Watts
///
/// Returns:
/// Maximum range in meters
///
/// Throws [ArgumentError] if any input parameters are non-positive
///
/// Example:
/// ```dart
/// final maxRange = calculateMaxRange(
///   transmitPower: 1000,           // 1 kW
///   antennaGain: 1000,              // 30 dB
///   wavelength: 0.03,               // 10 GHz
///   radarCrossSection: 1.0,         // 1 m²
///   minDetectableSignal: 1e-13,     // -100 dBm
/// );
/// ```
double calculateMaxRange({
  required double transmitPower,
  required double antennaGain,
  required double wavelength,
  required double radarCrossSection,
  required double minDetectableSignal,
}) {
  if (transmitPower <= 0 ||
      antennaGain <= 0 ||
      wavelength <= 0 ||
      radarCrossSection <= 0 ||
      minDetectableSignal <= 0) {
    throw ArgumentError('All parameters must be positive');
  }

  final numerator = transmitPower *
      pow(antennaGain, 2) *
      pow(wavelength, 2) *
      radarCrossSection;
  final denominator = pow(4 * pi, 3) * minDetectableSignal;

  final maxRange = pow(numerator / denominator, 0.25);

  return maxRange as double;
}

/// Calculate the received power from a target at a given range.
///
/// Parameters:
/// - [transmitPower]: Transmitted power in Watts
/// - [antennaGain]: Antenna gain (dimensionless, linear scale)
/// - [wavelength]: Radar wavelength in meters
/// - [radarCrossSection]: Target radar cross section in square meters
/// - [rangeM]: Range to target in meters
///
/// Returns:
/// Received power in Watts
///
/// Throws [ArgumentError] if any input parameters are non-positive
///
/// Example:
/// ```dart
/// final power = calculateReceivedPower(
///   transmitPower: 1000,
///   antennaGain: 1000,
///   wavelength: 0.03,
///   radarCrossSection: 1.0,
///   rangeM: 10000,
/// );
/// ```
double calculateReceivedPower({
  required double transmitPower,
  required double antennaGain,
  required double wavelength,
  required double radarCrossSection,
  required double rangeM,
}) {
  if (transmitPower <= 0 ||
      antennaGain <= 0 ||
      wavelength <= 0 ||
      radarCrossSection <= 0 ||
      rangeM <= 0) {
    throw ArgumentError('All parameters must be positive');
  }

  final numerator = transmitPower *
      pow(antennaGain, 2) *
      pow(wavelength, 2) *
      radarCrossSection;
  final denominator = pow(4 * pi, 3) * pow(rangeM, 4);

  final receivedPower = numerator / denominator;

  return receivedPower as double;
}
