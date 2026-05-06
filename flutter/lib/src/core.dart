import 'dart:math' as math;

import 'constants.dart';

/// Calculate wavelength from frequency.
///
/// Parameters:
/// - [frequency]: Frequency in Hz.
///
/// Returns wavelength in meters.
double calculateWavelength(double frequency) {
  return RadarConstants.speedOfLight / frequency;
}

/// Calculate frequency from wavelength.
///
/// Parameters:
/// - [wavelength]: Wavelength in meters.
///
/// Returns frequency in Hz.
double calculateFrequency(double wavelength) {
  return RadarConstants.speedOfLight / wavelength;
}

/// Calculate transmit antenna gain.
///
/// Parameters:
/// - [effectiveAperture]: Effective aperture in m².
/// - [wavelength]: Wavelength in meters.
///
/// Returns antenna gain (linear, not dB).
double calculateAntennaGain(double effectiveAperture, double wavelength) {
  return (4 * RadarConstants.pi * effectiveAperture) / (wavelength * wavelength);
}

/// Calculate effective aperture for rectangular antenna.
///
/// Parameters:
/// - [efficiency]: Antenna efficiency (0-1).
/// - [horizontalDimension]: Horizontal dimension in meters.
/// - [verticalDimension]: Vertical dimension in meters.
///
/// Returns effective aperture in m².
double calculateEffectiveApertureRect(
    double efficiency, double horizontalDimension, double verticalDimension) {
  return efficiency * horizontalDimension * verticalDimension;
}

/// Calculate effective aperture for circular antenna.
///
/// Parameters:
/// - [efficiency]: Antenna efficiency (0-1).
/// - [diameter]: Antenna diameter in meters.
///
/// Returns effective aperture in m².
double calculateEffectiveApertureCirc(double efficiency, double diameter) {
  return efficiency * RadarConstants.pi * diameter * diameter * 0.25;
}

/// Calculate maximum detection range using the radar range equation.
///
/// Parameters:
/// - [transmitPower]: Transmit power in watts.
/// - [antennaGain]: Antenna gain (linear, not dB).
/// - [wavelength]: Wavelength in meters.
/// - [radarCrossSection]: Target radar cross-section in m².
/// - [minDetectableSignal]: Minimum detectable signal in watts.
///
/// Returns maximum range in meters.
double calculateMaxRange(
  double transmitPower,
  double antennaGain,
  double wavelength,
  double radarCrossSection,
  double minDetectableSignal,
) {
  final antennaGainSquared = math.pow(antennaGain, 2).toDouble();
  final wavelengthSquared = math.pow(wavelength, 2).toDouble();
  final denominatorFactor = math.pow(4 * RadarConstants.pi, 3).toDouble();
  final numerator =
      transmitPower * antennaGainSquared * wavelengthSquared * radarCrossSection;
  final denominator = denominatorFactor * minDetectableSignal;
  return math.pow(numerator / denominator, 0.25).toDouble();
}

/// Calculate received power at a given range.
///
/// Parameters:
/// - [transmitPower]: Transmit power in watts.
/// - [transmitGain]: Transmit antenna gain (linear).
/// - [receiveGain]: Receive antenna gain (linear).
/// - [wavelength]: Wavelength in meters.
/// - [radarCrossSection]: Target radar cross-section in m².
/// - [range]: Range to target in meters.
///
/// Returns received power in watts.
double calculateReceivedPower(
  double transmitPower,
  double transmitGain,
  double receiveGain,
  double wavelength,
  double radarCrossSection,
  double range,
) {
  final wavelengthSquared = math.pow(wavelength, 2).toDouble();
  final rangeFourth = math.pow(range, 4).toDouble();
  final numerator =
      transmitPower * transmitGain * receiveGain * wavelengthSquared * radarCrossSection;
  final denominator =
      math.pow(4 * RadarConstants.pi, 3).toDouble() * rangeFourth;
  return numerator / denominator;
}

/// Calculate beamwidth (Gaussian approximation).
///
/// Parameters:
/// - [wavelength]: Wavelength in meters.
/// - [horizontalDimension]: Horizontal antenna dimension in meters.
///
/// Returns beamwidth in radians.
double calculateBeamwidth(double wavelength, double horizontalDimension) {
  return (RadarConstants.beamwidthCoefficient * RadarConstants.pi / 180) *
      (wavelength / horizontalDimension);
}

/// Calculate Doppler frequency shift.
///
/// Parameters:
/// - [velocity]: Target velocity in m/s (negative for closing).
/// - [wavelength]: Wavelength in meters.
///
/// Returns Doppler frequency shift in Hz.
double calculateDopplerFrequency(double velocity, double wavelength) {
  return -2 * velocity / wavelength;
}

/// Calculate velocity from Doppler shift.
///
/// Parameters:
/// - [dopplerFrequency]: Doppler frequency shift in Hz.
/// - [wavelength]: Wavelength in meters.
///
/// Returns velocity in m/s (negative for closing).
double calculateVelocityFromDoppler(double dopplerFrequency, double wavelength) {
  return -wavelength * dopplerFrequency / 2;
}

/// Calculate radar cross section of a sphere.
///
/// Parameters:
/// - [radius]: Sphere radius in meters.
///
/// Returns radar cross-section in m².
double calculateSphereRCS(double radius) {
  return RadarConstants.pi * radius * radius;
}

/// Calculate sphere radius from radar cross section.
///
/// Parameters:
/// - [radarCrossSection]: Radar cross-section in m².
///
/// Returns sphere radius in meters.
double calculateSphereRadius(double radarCrossSection) {
  return math.sqrt(radarCrossSection / RadarConstants.pi);
}
