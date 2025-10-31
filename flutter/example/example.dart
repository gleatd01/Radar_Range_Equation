import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  print('Radar Range Equation - Flutter/Dart Example\n');

  // Example 1: Calculate wavelength from frequency
  print('Example 1: Calculate wavelength');
  final frequency = 10e9; // 10 GHz
  final wavelength = calculateWavelength(frequency);
  print('Frequency: ${RadarConvert.hzToGhz(frequency)} GHz');
  print('Wavelength: ${wavelength * 100} cm\n');

  // Example 2: Calculate maximum range
  print('Example 2: Calculate maximum radar range');
  final transmitPower = 1000.0; // 1 kW
  final antennaGain = 1000.0; // Linear gain
  final radarCrossSection = 1.0; // 1 m² RCS
  final minDetectableSignal = 1e-13; // 0.1 pW

  final maxRange = calculateMaxRange(
    transmitPower,
    antennaGain,
    wavelength,
    radarCrossSection,
    minDetectableSignal,
  );

  print('Transmit Power: ${RadarConvert.wattsToKilowatts(transmitPower)} kW');
  print('Antenna Gain: $antennaGain (linear)');
  print('Wavelength: ${wavelength * 100} cm');
  print('Radar Cross Section: $radarCrossSection m²');
  print('Min Detectable Signal: $minDetectableSignal W');
  print('Maximum Range: ${maxRange.toStringAsFixed(2)} meters');
  print(
      'Maximum Range: ${RadarConvert.metersToKilometers(maxRange).toStringAsFixed(2)} km\n');

  // Example 3: Calculate Doppler frequency
  print('Example 3: Calculate Doppler frequency shift');
  final velocity = -100.0; // 100 m/s closing velocity
  final dopplerFreq = calculateDopplerFrequency(velocity, wavelength);
  print('Target velocity: ${-velocity} m/s (closing)');
  print('Wavelength: ${wavelength * 100} cm');
  print(
      'Doppler frequency shift: ${RadarConvert.hzToMhz(dopplerFreq).toStringAsFixed(3)} MHz\n');

  // Example 4: Calculate antenna properties
  print('Example 4: Calculate antenna effective aperture');
  final efficiency = 0.6;
  final diameter = 10.0; // 10 meters
  final effectiveAperture = calculateEffectiveApertureCirc(efficiency, diameter);
  print('Antenna diameter: $diameter m');
  print('Antenna efficiency: $efficiency');
  print(
      'Effective aperture: ${effectiveAperture.toStringAsFixed(2)} m²\n');

  // Example 5: Unit conversions
  print('Example 5: Unit conversions');
  final rangeMeters = 50000.0;
  print('Range: $rangeMeters meters');
  print(
      '  = ${RadarConvert.metersToKilometers(rangeMeters).toStringAsFixed(2)} km');
  print(
      '  = ${RadarConvert.metersToNauticalMiles(rangeMeters).toStringAsFixed(2)} nautical miles');
  print(
      '  = ${RadarConvert.metersToMiles(rangeMeters).toStringAsFixed(2)} miles\n');

  final powerDb = 30.0; // 30 dB
  final powerLinear = RadarConvert.dbToLinear(powerDb);
  print('Power: $powerDb dB = ${powerLinear.toStringAsFixed(2)} (linear)');
}
