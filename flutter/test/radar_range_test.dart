import 'package:test/test.dart';
import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  group('calculateMaxRange', () {
    test('basic calculation with typical values', () {
      final maxRange = calculateMaxRange(
        transmitPower: 1000,
        antennaGain: 1000,
        wavelength: 0.03,
        radarCrossSection: 1.0,
        minDetectableSignal: 1e-13,
      );

      expect(maxRange, greaterThan(0));
      expect(maxRange, isA<double>());
    });

    test('negative transmit power throws ArgumentError', () {
      expect(
        () => calculateMaxRange(
          transmitPower: -1000,
          antennaGain: 1000,
          wavelength: 0.03,
          radarCrossSection: 1.0,
          minDetectableSignal: 1e-13,
        ),
        throwsArgumentError,
      );
    });

    test('zero wavelength throws ArgumentError', () {
      expect(
        () => calculateMaxRange(
          transmitPower: 1000,
          antennaGain: 1000,
          wavelength: 0,
          radarCrossSection: 1.0,
          minDetectableSignal: 1e-13,
        ),
        throwsArgumentError,
      );
    });

    test('larger transmit power increases range', () {
      final range1 = calculateMaxRange(
        transmitPower: 1000,
        antennaGain: 1000,
        wavelength: 0.03,
        radarCrossSection: 1.0,
        minDetectableSignal: 1e-13,
      );

      final range2 = calculateMaxRange(
        transmitPower: 2000,
        antennaGain: 1000,
        wavelength: 0.03,
        radarCrossSection: 1.0,
        minDetectableSignal: 1e-13,
      );

      expect(range2, greaterThan(range1));
    });
  });

  group('calculateReceivedPower', () {
    test('basic calculation', () {
      final power = calculateReceivedPower(
        transmitPower: 1000,
        antennaGain: 1000,
        wavelength: 0.03,
        radarCrossSection: 1.0,
        rangeM: 10000,
      );

      expect(power, greaterThan(0));
      expect(power, isA<double>());
    });

    test('negative range throws ArgumentError', () {
      expect(
        () => calculateReceivedPower(
          transmitPower: 1000,
          antennaGain: 1000,
          wavelength: 0.03,
          radarCrossSection: 1.0,
          rangeM: -1000,
        ),
        throwsArgumentError,
      );
    });

    test('follows inverse fourth power law', () {
      final power1 = calculateReceivedPower(
        transmitPower: 1000,
        antennaGain: 1000,
        wavelength: 0.03,
        radarCrossSection: 1.0,
        rangeM: 10000,
      );

      final power2 = calculateReceivedPower(
        transmitPower: 1000,
        antennaGain: 1000,
        wavelength: 0.03,
        radarCrossSection: 1.0,
        rangeM: 20000,
      );

      // Power at 2x range should be 1/16 (2^4) of power at 1x range
      expect((power2 * 16 - power1).abs(), lessThan(power1 * 0.01));
    });

    test('consistency with max range', () {
      const transmitPower = 1000.0;
      const antennaGain = 1000.0;
      const wavelength = 0.03;
      const radarCrossSection = 1.0;
      const minDetectableSignal = 1e-13;

      final maxRange = calculateMaxRange(
        transmitPower: transmitPower,
        antennaGain: antennaGain,
        wavelength: wavelength,
        radarCrossSection: radarCrossSection,
        minDetectableSignal: minDetectableSignal,
      );

      final receivedPower = calculateReceivedPower(
        transmitPower: transmitPower,
        antennaGain: antennaGain,
        wavelength: wavelength,
        radarCrossSection: radarCrossSection,
        rangeM: maxRange,
      );

      // Should be very close to minDetectableSignal
      final relativeError =
          (receivedPower - minDetectableSignal).abs() / minDetectableSignal;
      expect(relativeError, lessThan(0.01));
    });
  });
}
