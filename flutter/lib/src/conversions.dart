import 'constants.dart';
import 'dart:math' as math;

/// Unit conversion utilities.
class RadarConvert {
  /// Convert linear value to dB.
  static double linearToDb(double linear) {
    if (linear <= 0) {
      return double.negativeInfinity;
    }
    return 10 * math.log(linear) / math.ln10;
  }

  /// Convert dB to linear value.
  static double dbToLinear(double db) {
    return math.pow(10, db / 10).toDouble();
  }

  /// Convert radians to degrees.
  static double radToDeg(double radians) {
    return radians * 180 / RadarConstants.pi;
  }

  /// Convert degrees to radians.
  static double degToRad(double degrees) {
    return degrees * RadarConstants.pi / 180;
  }

  /// Convert meters to nautical miles.
  static double metersToNauticalMiles(double meters) {
    return meters / 1852.0;
  }

  /// Convert nautical miles to meters.
  static double nauticalMilesToMeters(double nauticalMiles) {
    return nauticalMiles * 1852.0;
  }

  /// Convert meters to miles.
  static double metersToMiles(double meters) {
    return meters / 1609.34;
  }

  /// Convert miles to meters.
  static double milesToMeters(double miles) {
    return miles * 1609.34;
  }

  /// Convert meters to kilometers.
  static double metersToKilometers(double meters) {
    return meters / 1000.0;
  }

  /// Convert kilometers to meters.
  static double kilometersToMeters(double kilometers) {
    return kilometers * 1000.0;
  }

  /// Convert feet to meters.
  static double feetToMeters(double feet) {
    return feet * 0.3048;
  }

  /// Convert meters to feet.
  static double metersToFeet(double meters) {
    return meters / 0.3048;
  }

  /// Convert Hz to MHz.
  static double hzToMhz(double hz) {
    return hz / 1e6;
  }

  /// Convert MHz to Hz.
  static double mhzToHz(double mhz) {
    return mhz * 1e6;
  }

  /// Convert Hz to GHz.
  static double hzToGhz(double hz) {
    return hz / 1e9;
  }

  /// Convert GHz to Hz.
  static double ghzToHz(double ghz) {
    return ghz * 1e9;
  }

  /// Convert watts to kilowatts.
  static double wattsToKilowatts(double watts) {
    return watts / 1000.0;
  }

  /// Convert kilowatts to watts.
  static double kilowattsToWatts(double kilowatts) {
    return kilowatts * 1000.0;
  }

  /// Convert watts to milliwatts.
  static double wattsToMilliwatts(double watts) {
    return watts * 1000.0;
  }

  /// Convert milliwatts to watts.
  static double milliwattsToWatts(double milliwatts) {
    return milliwatts / 1000.0;
  }
}
