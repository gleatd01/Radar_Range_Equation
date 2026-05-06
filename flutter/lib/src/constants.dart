import 'dart:math' as math;

/// Physical constants used in radar calculations.
class RadarConstants {
  /// Speed of light in m/s.
  static const double speedOfLight = 299792458.0;

  /// Boltzmann constant in J/K.
  static const double boltzmann = 1.380649e-23;

  /// Pi.
  static const double pi = math.pi;

  /// Beamwidth coefficient for Gaussian approximation.
  static const double beamwidthCoefficient = 65.0;
}
