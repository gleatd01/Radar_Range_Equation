import 'package:test/test.dart';
import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  group('Radar Range Equation Tests', () {
    test('Calculate wavelength from frequency', () {
      final wavelength = calculateWavelength(10e9); // 10 GHz
      expect(wavelength, closeTo(0.0299792458, 1e-6));
    });

    test('Calculate frequency from wavelength', () {
      final frequency = calculateFrequency(0.03);
      expect(frequency, closeTo(9.99308193e9, 1e3));
    });

    test('Calculate antenna gain', () {
      final effectiveAperture = 10.0; // 10 m²
      final wavelength = 0.03; // 0.03 m
      final gain = calculateAntennaGain(effectiveAperture, wavelength);
      expect(gain, greaterThan(0));
      expect(gain, closeTo(139626.34, 1));
    });

    test('Calculate effective aperture for rectangular antenna', () {
      final efficiency = 0.6;
      final horizontal = 5.0;
      final vertical = 4.0;
      final aperture =
          calculateEffectiveApertureRect(efficiency, horizontal, vertical);
      expect(aperture, closeTo(12.0, 1e-6));
    });

    test('Calculate effective aperture for circular antenna', () {
      final efficiency = 0.6;
      final diameter = 10.0;
      final aperture = calculateEffectiveApertureCirc(efficiency, diameter);
      expect(aperture, closeTo(47.123889803, 1e-6));
    });

    test('Calculate maximum range', () {
      final transmitPower = 1000.0; // 1 kW
      final antennaGain = 1000.0;
      final wavelength = 0.03; // 3 cm
      final radarCrossSection = 1.0; // 1 m²
      final minDetectableSignal = 1e-13; // 0.1 pW

      final maxRange = calculateMaxRange(
        transmitPower,
        antennaGain,
        wavelength,
        radarCrossSection,
        minDetectableSignal,
      );

      expect(maxRange, greaterThan(0));
      // The exact value depends on the formula
      expect(maxRange, closeTo(40278.97, 1));
    });

    test('Calculate received power', () {
      final transmitPower = 1000.0;
      final transmitGain = 1000.0;
      final receiveGain = 1000.0;
      final wavelength = 0.03;
      final radarCrossSection = 1.0;
      final range = 10000.0; // 10 km

      final receivedPower = calculateReceivedPower(
        transmitPower,
        transmitGain,
        receiveGain,
        wavelength,
        radarCrossSection,
        range,
      );

      expect(receivedPower, greaterThan(0));
    });

    test('Calculate beamwidth', () {
      final wavelength = 0.03;
      final horizontalDimension = 5.0;
      final beamwidth = calculateBeamwidth(wavelength, horizontalDimension);
      expect(beamwidth, greaterThan(0));
      expect(beamwidth, closeTo(0.00679607, 1e-6));
    });

    test('Calculate Doppler frequency', () {
      final velocity = -100.0; // 100 m/s closing
      final wavelength = 0.03;
      final dopplerFreq = calculateDopplerFrequency(velocity, wavelength);
      expect(dopplerFreq, closeTo(6666.67, 0.01));
    });

    test('Calculate velocity from Doppler', () {
      final dopplerFreq = 6666.67;
      final wavelength = 0.03;
      final velocity = calculateVelocityFromDoppler(dopplerFreq, wavelength);
      expect(velocity, closeTo(-100.0, 0.01));
    });

    test('Calculate sphere RCS', () {
      final radius = 1.0;
      final rcs = calculateSphereRCS(radius);
      expect(rcs, closeTo(RadarConstants.pi, 1e-6));
    });

    test('Calculate sphere radius from RCS', () {
      final rcs = RadarConstants.pi;
      final radius = calculateSphereRadius(rcs);
      expect(radius, closeTo(1.0, 1e-6));
    });
  });

  group('Unit Conversion Tests', () {
    test('Linear to dB conversion', () {
      expect(RadarConvert.linearToDb(10), closeTo(10, 1e-6));
      expect(RadarConvert.linearToDb(100), closeTo(20, 1e-6));
      expect(RadarConvert.linearToDb(1), closeTo(0, 1e-6));
    });

    test('dB to Linear conversion', () {
      expect(RadarConvert.dbToLinear(10), closeTo(10, 1e-6));
      expect(RadarConvert.dbToLinear(20), closeTo(100, 1e-6));
      expect(RadarConvert.dbToLinear(0), closeTo(1, 1e-6));
    });

    test('Radian to Degree conversion', () {
      expect(RadarConvert.radToDeg(RadarConstants.pi), closeTo(180, 1e-6));
      expect(RadarConvert.radToDeg(RadarConstants.pi / 2), closeTo(90, 1e-6));
    });

    test('Degree to Radian conversion', () {
      expect(RadarConvert.degToRad(180), closeTo(RadarConstants.pi, 1e-6));
      expect(RadarConvert.degToRad(90), closeTo(RadarConstants.pi / 2, 1e-6));
    });

    test('Meters to Nautical Miles conversion', () {
      expect(RadarConvert.metersToNauticalMiles(1852), closeTo(1, 1e-6));
    });

    test('Feet to Meters conversion', () {
      expect(RadarConvert.feetToMeters(1), closeTo(0.3048, 1e-6));
    });

    test('Hz to MHz conversion', () {
      expect(RadarConvert.hzToMhz(1e6), closeTo(1, 1e-6));
    });

    test('Hz to GHz conversion', () {
      expect(RadarConvert.hzToGhz(1e9), closeTo(1, 1e-6));
    });
  });
}
