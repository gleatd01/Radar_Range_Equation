/// A Dart library for calculating radar range equations.
///
/// This library provides functions for computing:
/// - Maximum detection range
/// - Received power
/// - Wavelength and frequency conversions
/// - Antenna gain calculations
/// - Various unit conversions
library radar_range_equation;

import 'dart:math' as math;

/// Physical constants
class RadarConstants {
  /// Speed of light in m/s
  static const double speedOfLight = 299792458.0;
  
  /// Boltzmann constant in J/K
  static const double boltzmann = 1.380649e-23;
  
  /// Pi
  static const double pi = math.pi;
}

/// Calculate wavelength from frequency
///
/// Parameters:
/// - [frequency]: Frequency in Hz
///
/// Returns wavelength in meters
double calculateWavelength(double frequency) {
  return RadarConstants.speedOfLight / frequency;
}

/// Calculate frequency from wavelength
///
/// Parameters:
/// - [wavelength]: Wavelength in meters
///
/// Returns frequency in Hz
double calculateFrequency(double wavelength) {
  return RadarConstants.speedOfLight / wavelength;
}

/// Calculate transmit antenna gain
///
/// Parameters:
/// - [effectiveAperture]: Effective aperture in m²
/// - [wavelength]: Wavelength in meters
///
/// Returns antenna gain (linear, not dB)
double calculateAntennaGain(double effectiveAperture, double wavelength) {
  return (4 * RadarConstants.pi * effectiveAperture) / (wavelength * wavelength);
}

/// Calculate effective aperture for rectangular antenna
///
/// Parameters:
/// - [efficiency]: Antenna efficiency (0-1)
/// - [horizontalDimension]: Horizontal dimension in meters
/// - [verticalDimension]: Vertical dimension in meters
///
/// Returns effective aperture in m²
double calculateEffectiveApertureRect(
    double efficiency, double horizontalDimension, double verticalDimension) {
  return efficiency * horizontalDimension * verticalDimension;
}

/// Calculate effective aperture for circular antenna
///
/// Parameters:
/// - [efficiency]: Antenna efficiency (0-1)
/// - [diameter]: Antenna diameter in meters
///
/// Returns effective aperture in m²
double calculateEffectiveApertureCirc(double efficiency, double diameter) {
  return efficiency * RadarConstants.pi * math.pow(diameter / 2, 2);
}

/// Calculate maximum detection range using the radar range equation
///
/// Parameters:
/// - [transmitPower]: Transmit power in watts
/// - [antennaGain]: Antenna gain (linear, not dB)
/// - [wavelength]: Wavelength in meters
/// - [radarCrossSection]: Target radar cross-section in m²
/// - [minDetectableSignal]: Minimum detectable signal in watts
///
/// Returns maximum range in meters
double calculateMaxRange(
  double transmitPower,
  double antennaGain,
  double wavelength,
  double radarCrossSection,
  double minDetectableSignal,
) {
  final numerator = transmitPower *
      math.pow(antennaGain, 2) *
      math.pow(wavelength, 2) *
      radarCrossSection;
  final denominator =
      math.pow(4 * RadarConstants.pi, 3) * minDetectableSignal;
  return math.pow(numerator / denominator, 0.25);
}

/// Calculate received power at a given range
///
/// Parameters:
/// - [transmitPower]: Transmit power in watts
/// - [transmitGain]: Transmit antenna gain (linear)
/// - [receiveGain]: Receive antenna gain (linear)
/// - [wavelength]: Wavelength in meters
/// - [radarCrossSection]: Target radar cross-section in m²
/// - [range]: Range to target in meters
///
/// Returns received power in watts
double calculateReceivedPower(
  double transmitPower,
  double transmitGain,
  double receiveGain,
  double wavelength,
  double radarCrossSection,
  double range,
) {
  final numerator = transmitPower *
      transmitGain *
      receiveGain *
      math.pow(wavelength, 2) *
      radarCrossSection;
  final denominator = math.pow(4 * RadarConstants.pi, 3) * math.pow(range, 4);
  return numerator / denominator;
}

/// Calculate beamwidth (Gaussian approximation)
///
/// Parameters:
/// - [wavelength]: Wavelength in meters
/// - [horizontalDimension]: Horizontal antenna dimension in meters
///
/// Returns beamwidth in radians
double calculateBeamwidth(double wavelength, double horizontalDimension) {
  return (65 * RadarConstants.pi / 180) * (wavelength / horizontalDimension);
}

/// Calculate Doppler frequency shift
///
/// Parameters:
/// - [velocity]: Target velocity in m/s (negative for closing)
/// - [wavelength]: Wavelength in meters
///
/// Returns Doppler frequency shift in Hz
double calculateDopplerFrequency(double velocity, double wavelength) {
  return -2 * velocity / wavelength;
}

/// Calculate velocity from Doppler shift
///
/// Parameters:
/// - [dopplerFrequency]: Doppler frequency shift in Hz
/// - [wavelength]: Wavelength in meters
///
/// Returns velocity in m/s (negative for closing)
double calculateVelocityFromDoppler(double dopplerFrequency, double wavelength) {
  return -wavelength * dopplerFrequency / 2;
}

/// Calculate radar cross section of a sphere
///
/// Parameters:
/// - [radius]: Sphere radius in meters
///
/// Returns radar cross-section in m²
double calculateSphereRCS(double radius) {
  return RadarConstants.pi * radius * radius;
}

/// Calculate sphere radius from radar cross section
///
/// Parameters:
/// - [radarCrossSection]: Radar cross-section in m²
///
/// Returns sphere radius in meters
double calculateSphereRadius(double radarCrossSection) {
  return math.sqrt(radarCrossSection / RadarConstants.pi);
}

/// Unit conversion utilities
class RadarConvert {
  /// Convert linear value to dB
  static double linearToDb(double linear) {
    if (linear <= 0) {
      return double.negativeInfinity;
    }
    return 10 * math.log(linear) / math.ln10;
  }

  /// Convert dB to linear value
  static double dbToLinear(double db) {
    return math.pow(10, db / 10).toDouble();
  }

  /// Convert radians to degrees
  static double radToDeg(double radians) {
    return radians * 180 / RadarConstants.pi;
  }

  /// Convert degrees to radians
  static double degToRad(double degrees) {
    return degrees * RadarConstants.pi / 180;
  }

  /// Convert meters to nautical miles
  static double metersToNauticalMiles(double meters) {
    return meters / 1852.0;
  }

  /// Convert nautical miles to meters
  static double nauticalMilesToMeters(double nauticalMiles) {
    return nauticalMiles * 1852.0;
  }

  /// Convert meters to miles
  static double metersToMiles(double meters) {
    return meters / 1609.34;
  }

  /// Convert miles to meters
  static double milesToMeters(double miles) {
    return miles * 1609.34;
  }

  /// Convert meters to kilometers
  static double metersToKilometers(double meters) {
    return meters / 1000.0;
  }

  /// Convert kilometers to meters
  static double kilometersToMeters(double kilometers) {
    return kilometers * 1000.0;
  }

  /// Convert feet to meters
  static double feetToMeters(double feet) {
    return feet * 0.3048;
  }

  /// Convert meters to feet
  static double metersToFeet(double meters) {
    return meters / 0.3048;
  }

  /// Convert Hz to MHz
  static double hzToMhz(double hz) {
    return hz / 1e6;
  }

  /// Convert MHz to Hz
  static double mhzToHz(double mhz) {
    return mhz * 1e6;
  }

  /// Convert Hz to GHz
  static double hzToGhz(double hz) {
    return hz / 1e9;
  }

  /// Convert GHz to Hz
  static double ghzToHz(double ghz) {
    return ghz * 1e9;
  }

  /// Convert watts to kilowatts
  static double wattsToKilowatts(double watts) {
    return watts / 1000.0;
  }

  /// Convert kilowatts to watts
  static double kilowattsToWatts(double kilowatts) {
    return kilowatts * 1000.0;
  }

  /// Convert watts to milliwatts
  static double wattsToMilliwatts(double watts) {
    return watts * 1000.0;
  }

  /// Convert milliwatts to watts
  static double milliwattsToWatts(double milliwatts) {
    return milliwatts / 1000.0;
  }
}
