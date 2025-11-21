import {
  calculateMaxRange,
  calculateReceivedPower,
  calculateWavelength,
  calculateAntennaGain,
  calculateEffectiveApertureRect,
  calculateEffectiveApertureCirc,
  convert,
  constants
} from '../src/index';

describe('Radar Range Equation', () => {
  describe('calculateMaxRange', () => {
    it('should calculate maximum range correctly', () => {
      const transmitPower = 1000; // W
      const antennaGain = 1000;
      const wavelength = 0.03; // m
      const radarCrossSection = 1.0; // m²
      const minDetectableSignal = 1e-13; // W

      const maxRange = calculateMaxRange(
        transmitPower,
        antennaGain,
        wavelength,
        radarCrossSection,
        minDetectableSignal
      );

      expect(maxRange).toBeGreaterThan(0);
      expect(maxRange).toBeCloseTo(46148.04, 0); // Expected value for monostatic radar with G³
    });

    it('should throw error for invalid transmit power', () => {
      expect(() => calculateMaxRange(-1, 1000, 0.03, 1.0, 1e-13)).toThrow('Transmit power must be positive');
      expect(() => calculateMaxRange(0, 1000, 0.03, 1.0, 1e-13)).toThrow('Transmit power must be positive');
    });

    it('should throw error for invalid antenna gain', () => {
      expect(() => calculateMaxRange(1000, -1, 0.03, 1.0, 1e-13)).toThrow('Antenna gain must be positive');
      expect(() => calculateMaxRange(1000, 0, 0.03, 1.0, 1e-13)).toThrow('Antenna gain must be positive');
    });

    it('should throw error for invalid wavelength', () => {
      expect(() => calculateMaxRange(1000, 1000, -0.03, 1.0, 1e-13)).toThrow('Wavelength must be positive');
      expect(() => calculateMaxRange(1000, 1000, 0, 1.0, 1e-13)).toThrow('Wavelength must be positive');
    });

    it('should throw error for invalid radar cross-section', () => {
      expect(() => calculateMaxRange(1000, 1000, 0.03, -1.0, 1e-13)).toThrow('Radar cross-section must be positive');
      expect(() => calculateMaxRange(1000, 1000, 0.03, 0, 1e-13)).toThrow('Radar cross-section must be positive');
    });

    it('should throw error for invalid minimum detectable signal', () => {
      expect(() => calculateMaxRange(1000, 1000, 0.03, 1.0, -1e-13)).toThrow('Minimum detectable signal must be positive');
      expect(() => calculateMaxRange(1000, 1000, 0.03, 1.0, 0)).toThrow('Minimum detectable signal must be positive');
    });
  });

  describe('calculateReceivedPower', () => {
    it('should calculate received power correctly', () => {
      const transmitPower = 1000; // W
      const antennaGain = 1000;
      const wavelength = 0.03; // m
      const radarCrossSection = 1.0; // m²
      const range = 10000; // m

      const receivedPower = calculateReceivedPower(
        transmitPower,
        antennaGain,
        wavelength,
        radarCrossSection,
        range
      );

      expect(receivedPower).toBeGreaterThan(0);
      expect(receivedPower).toBeLessThan(transmitPower);
    });

    it('should throw error for invalid range', () => {
      expect(() => calculateReceivedPower(1000, 1000, 0.03, 1.0, -10000)).toThrow('Range must be positive');
      expect(() => calculateReceivedPower(1000, 1000, 0.03, 1.0, 0)).toThrow('Range must be positive');
    });

    it('should show inverse fourth power relationship', () => {
      const transmitPower = 1000;
      const antennaGain = 1000;
      const wavelength = 0.03;
      const radarCrossSection = 1.0;

      const power1 = calculateReceivedPower(transmitPower, antennaGain, wavelength, radarCrossSection, 1000);
      const power2 = calculateReceivedPower(transmitPower, antennaGain, wavelength, radarCrossSection, 2000);

      // Power should decrease by factor of 16 when range doubles (R^4)
      expect(power1 / power2).toBeCloseTo(16, 1);
    });
  });

  describe('calculateWavelength', () => {
    it('should calculate wavelength from frequency', () => {
      const frequency = 10e9; // 10 GHz
      const wavelength = calculateWavelength(frequency);

      expect(wavelength).toBeCloseTo(0.0299792458, 5);
    });

    it('should use custom speed of light if provided', () => {
      const frequency = 10e9;
      const customSpeed = 3e8;
      const wavelength = calculateWavelength(frequency, customSpeed);

      expect(wavelength).toBeCloseTo(0.03, 5);
    });

    it('should throw error for invalid frequency', () => {
      expect(() => calculateWavelength(-10e9)).toThrow('Frequency must be positive');
      expect(() => calculateWavelength(0)).toThrow('Frequency must be positive');
    });
  });

  describe('calculateAntennaGain', () => {
    it('should calculate antenna gain from effective aperture', () => {
      const effectiveAperture = 1.0; // m²
      const wavelength = 0.03; // m

      const gain = calculateAntennaGain(effectiveAperture, wavelength);

      expect(gain).toBeGreaterThan(0);
      expect(gain).toBeCloseTo(13962.634, 2);
    });

    it('should throw error for invalid effective aperture', () => {
      expect(() => calculateAntennaGain(-1.0, 0.03)).toThrow('Effective aperture must be positive');
      expect(() => calculateAntennaGain(0, 0.03)).toThrow('Effective aperture must be positive');
    });

    it('should throw error for invalid wavelength', () => {
      expect(() => calculateAntennaGain(1.0, -0.03)).toThrow('Wavelength must be positive');
      expect(() => calculateAntennaGain(1.0, 0)).toThrow('Wavelength must be positive');
    });
  });

  describe('calculateEffectiveApertureRect', () => {
    it('should calculate effective aperture for rectangular antenna', () => {
      const efficiency = 0.6;
      const horizontal = 2.0; // m
      const vertical = 1.5; // m

      const aperture = calculateEffectiveApertureRect(efficiency, horizontal, vertical);

      expect(aperture).toBeCloseTo(1.8, 5);
    });

    it('should throw error for invalid efficiency', () => {
      expect(() => calculateEffectiveApertureRect(-0.1, 2.0, 1.5)).toThrow('Efficiency must be between 0 and 1');
      expect(() => calculateEffectiveApertureRect(0, 2.0, 1.5)).toThrow('Efficiency must be between 0 and 1');
      expect(() => calculateEffectiveApertureRect(1.5, 2.0, 1.5)).toThrow('Efficiency must be between 0 and 1');
    });

    it('should throw error for invalid dimensions', () => {
      expect(() => calculateEffectiveApertureRect(0.6, -2.0, 1.5)).toThrow('Horizontal dimension must be positive');
      expect(() => calculateEffectiveApertureRect(0.6, 2.0, -1.5)).toThrow('Vertical dimension must be positive');
    });
  });

  describe('calculateEffectiveApertureCirc', () => {
    it('should calculate effective aperture for circular antenna', () => {
      const efficiency = 0.6;
      const diameter = 2.0; // m

      const aperture = calculateEffectiveApertureCirc(efficiency, diameter);

      expect(aperture).toBeCloseTo(1.88495559, 5);
    });

    it('should throw error for invalid efficiency', () => {
      expect(() => calculateEffectiveApertureCirc(-0.1, 2.0)).toThrow('Efficiency must be between 0 and 1');
      expect(() => calculateEffectiveApertureCirc(0, 2.0)).toThrow('Efficiency must be between 0 and 1');
      expect(() => calculateEffectiveApertureCirc(1.5, 2.0)).toThrow('Efficiency must be between 0 and 1');
    });

    it('should throw error for invalid diameter', () => {
      expect(() => calculateEffectiveApertureCirc(0.6, -2.0)).toThrow('Diameter must be positive');
      expect(() => calculateEffectiveApertureCirc(0.6, 0)).toThrow('Diameter must be positive');
    });
  });

  describe('convert utilities', () => {
    describe('linearToDb and dbToLinear', () => {
      it('should convert between linear and dB', () => {
        expect(convert.linearToDb(10)).toBeCloseTo(10, 5);
        expect(convert.linearToDb(100)).toBeCloseTo(20, 5);
        expect(convert.linearToDb(1000)).toBeCloseTo(30, 5);

        expect(convert.dbToLinear(10)).toBeCloseTo(10, 5);
        expect(convert.dbToLinear(20)).toBeCloseTo(100, 5);
        expect(convert.dbToLinear(30)).toBeCloseTo(1000, 5);
      });

      it('should return -Infinity for zero or negative linear values', () => {
        expect(convert.linearToDb(0)).toBe(-Infinity);
        expect(convert.linearToDb(-10)).toBe(-Infinity);
      });
    });

    describe('angle conversions', () => {
      it('should convert between radians and degrees', () => {
        expect(convert.radToDeg(Math.PI)).toBeCloseTo(180, 5);
        expect(convert.radToDeg(Math.PI / 2)).toBeCloseTo(90, 5);
        expect(convert.degToRad(180)).toBeCloseTo(Math.PI, 5);
        expect(convert.degToRad(90)).toBeCloseTo(Math.PI / 2, 5);
      });
    });

    describe('distance conversions', () => {
      it('should convert between meters and kilometers', () => {
        expect(convert.metersToKm(1000)).toBe(1);
        expect(convert.kmToMeters(1)).toBe(1000);
      });

      it('should convert between meters and nautical miles', () => {
        expect(convert.metersToNauticalMiles(1852)).toBe(1);
        expect(convert.nauticalMilesToMeters(1)).toBe(1852);
      });

      it('should convert between meters and miles', () => {
        expect(convert.metersToMiles(1609.34)).toBeCloseTo(1, 5);
        expect(convert.milesToMeters(1)).toBeCloseTo(1609.34, 2);
      });
    });

    describe('frequency conversions', () => {
      it('should convert between Hz and GHz', () => {
        expect(convert.hzToGhz(10e9)).toBe(10);
        expect(convert.ghzToHz(10)).toBe(10e9);
      });
    });

    describe('power conversions', () => {
      it('should convert between watts and kilowatts', () => {
        expect(convert.wattsToKw(1000)).toBe(1);
        expect(convert.kwToWatts(1)).toBe(1000);
      });
    });
  });

  describe('constants', () => {
    it('should provide correct physical constants', () => {
      expect(constants.SPEED_OF_LIGHT).toBe(299792458);
      expect(constants.BOLTZMANN).toBeCloseTo(1.380649e-23, 30);
      expect(constants.PI).toBe(Math.PI);
    });
  });
});
