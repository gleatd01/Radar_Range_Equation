import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  print('Radar Range Equation Example');
  print('=' * 40);

  // System parameters
  const transmitPower = 1000.0; // 1 kW
  const antennaGain = 1000.0; // 30 dBi (linear scale)
  const wavelength = 0.03; // 10 GHz (c/f = 3e8/10e9)
  const radarCrossSection = 1.0; // 1 m²
  const minDetectableSignal = 1e-13; // -100 dBm

  print('Transmit Power: $transmitPower W');
  print('Antenna Gain: $antennaGain (linear)');
  print('Wavelength: $wavelength m');
  print('Target RCS: $radarCrossSection m²');
  print('Min Detectable Signal: $minDetectableSignal W');
  print('');

  // Calculate maximum range
  final maxRange = calculateMaxRange(
    transmitPower: transmitPower,
    antennaGain: antennaGain,
    wavelength: wavelength,
    radarCrossSection: radarCrossSection,
    minDetectableSignal: minDetectableSignal,
  );

  print(
      'Maximum Detection Range: ${maxRange.toStringAsFixed(2)} m (${(maxRange / 1000).toStringAsFixed(2)} km)');
  print('');

  // Calculate received power at various ranges
  print('Received Power at Different Ranges:');
  print('-' * 40);
  for (final rangeKm in [5, 10, 15, 20]) {
    final rangeM = rangeKm * 1000.0;
    final power = calculateReceivedPower(
      transmitPower: transmitPower,
      antennaGain: antennaGain,
      wavelength: wavelength,
      radarCrossSection: radarCrossSection,
      rangeM: rangeM,
    );
    print('  $rangeKm km: ${power.toStringAsExponential(2)} W');
  }
}
