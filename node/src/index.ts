/**
 * Radar Range Equation Package
 * 
 * Provides radar range equation calculations for determining maximum detection range
 * and received power for radar systems.
 */

/**
 * Calculate the maximum detection range of a radar system.
 * 
 * The radar range equation relates the range of a radar to the characteristics
 * of the transmitter, receiver, antenna, target, and environment.
 * 
 * Formula: R_max = ((P_t × G² × λ² × σ) / ((4π)³ × P_min))^(1/4)
 * 
 * @param transmitPower - Transmit power in watts (P_t)
 * @param antennaGain - Antenna gain (dimensionless) (G)
 * @param wavelength - Wavelength in meters (λ)
 * @param radarCrossSection - Radar cross-section in square meters (σ)
 * @param minDetectableSignal - Minimum detectable signal in watts (P_min)
 * @returns Maximum detection range in meters
 * 
 * @example
 * ```typescript
 * const maxRange = calculateMaxRange(1000, 1000, 0.03, 1.0, 1e-13);
 * console.log(`Maximum range: ${maxRange.toFixed(2)} meters`);
 * ```
 */
export function calculateMaxRange(
  transmitPower: number,
  antennaGain: number,
  wavelength: number,
  radarCrossSection: number,
  minDetectableSignal: number
): number {
  if (transmitPower <= 0) {
    throw new Error('Transmit power must be positive');
  }
  if (antennaGain <= 0) {
    throw new Error('Antenna gain must be positive');
  }
  if (wavelength <= 0) {
    throw new Error('Wavelength must be positive');
  }
  if (radarCrossSection <= 0) {
    throw new Error('Radar cross-section must be positive');
  }
  if (minDetectableSignal <= 0) {
    throw new Error('Minimum detectable signal must be positive');
  }

  const fourPi = 4 * Math.PI;
  const numerator = transmitPower * Math.pow(antennaGain, 2) * Math.pow(wavelength, 2) * radarCrossSection;
  const denominator = Math.pow(fourPi, 3) * minDetectableSignal;
  
  return Math.pow(numerator / denominator, 0.25);
}

/**
 * Calculate the received power at a given range.
 * 
 * This function calculates the power received by a radar from a target
 * at a specified range.
 * 
 * Formula: P_r = (P_t × G² × λ² × σ) / ((4π)³ × R⁴)
 * 
 * @param transmitPower - Transmit power in watts (P_t)
 * @param antennaGain - Antenna gain (dimensionless) (G)
 * @param wavelength - Wavelength in meters (λ)
 * @param radarCrossSection - Radar cross-section in square meters (σ)
 * @param range - Range to target in meters (R)
 * @returns Received power in watts
 * 
 * @example
 * ```typescript
 * const receivedPower = calculateReceivedPower(1000, 1000, 0.03, 1.0, 10000);
 * console.log(`Received power: ${receivedPower} watts`);
 * ```
 */
export function calculateReceivedPower(
  transmitPower: number,
  antennaGain: number,
  wavelength: number,
  radarCrossSection: number,
  range: number
): number {
  if (transmitPower <= 0) {
    throw new Error('Transmit power must be positive');
  }
  if (antennaGain <= 0) {
    throw new Error('Antenna gain must be positive');
  }
  if (wavelength <= 0) {
    throw new Error('Wavelength must be positive');
  }
  if (radarCrossSection <= 0) {
    throw new Error('Radar cross-section must be positive');
  }
  if (range <= 0) {
    throw new Error('Range must be positive');
  }

  const fourPi = 4 * Math.PI;
  const numerator = transmitPower * Math.pow(antennaGain, 2) * Math.pow(wavelength, 2) * radarCrossSection;
  const denominator = Math.pow(fourPi, 3) * Math.pow(range, 4);
  
  return numerator / denominator;
}

/**
 * Calculate wavelength from frequency.
 * 
 * @param frequency - Frequency in Hz
 * @param speedOfLight - Speed of light in m/s (default: 299792458)
 * @returns Wavelength in meters
 * 
 * @example
 * ```typescript
 * const wavelength = calculateWavelength(10e9); // 10 GHz
 * console.log(`Wavelength: ${wavelength} meters`);
 * ```
 */
export function calculateWavelength(frequency: number, speedOfLight: number = 299792458): number {
  if (frequency <= 0) {
    throw new Error('Frequency must be positive');
  }
  return speedOfLight / frequency;
}

/**
 * Calculate antenna gain from effective aperture.
 * 
 * @param effectiveAperture - Effective aperture in square meters
 * @param wavelength - Wavelength in meters
 * @returns Antenna gain (dimensionless)
 * 
 * @example
 * ```typescript
 * const gain = calculateAntennaGain(1.0, 0.03);
 * console.log(`Antenna gain: ${gain}`);
 * ```
 */
export function calculateAntennaGain(effectiveAperture: number, wavelength: number): number {
  if (effectiveAperture <= 0) {
    throw new Error('Effective aperture must be positive');
  }
  if (wavelength <= 0) {
    throw new Error('Wavelength must be positive');
  }
  return (4 * Math.PI * effectiveAperture) / Math.pow(wavelength, 2);
}

/**
 * Calculate effective aperture for a rectangular antenna.
 * 
 * @param efficiency - Antenna efficiency (0 to 1)
 * @param horizontalDimension - Horizontal dimension in meters
 * @param verticalDimension - Vertical dimension in meters
 * @returns Effective aperture in square meters
 * 
 * @example
 * ```typescript
 * const aperture = calculateEffectiveApertureRect(0.6, 2.0, 1.5);
 * console.log(`Effective aperture: ${aperture} m²`);
 * ```
 */
export function calculateEffectiveApertureRect(
  efficiency: number,
  horizontalDimension: number,
  verticalDimension: number
): number {
  if (efficiency <= 0 || efficiency > 1) {
    throw new Error('Efficiency must be between 0 and 1');
  }
  if (horizontalDimension <= 0) {
    throw new Error('Horizontal dimension must be positive');
  }
  if (verticalDimension <= 0) {
    throw new Error('Vertical dimension must be positive');
  }
  return efficiency * horizontalDimension * verticalDimension;
}

/**
 * Calculate effective aperture for a circular antenna.
 * 
 * @param efficiency - Antenna efficiency (0 to 1)
 * @param diameter - Antenna diameter in meters
 * @returns Effective aperture in square meters
 * 
 * @example
 * ```typescript
 * const aperture = calculateEffectiveApertureCirc(0.6, 2.0);
 * console.log(`Effective aperture: ${aperture} m²`);
 * ```
 */
export function calculateEffectiveApertureCirc(efficiency: number, diameter: number): number {
  if (efficiency <= 0 || efficiency > 1) {
    throw new Error('Efficiency must be between 0 and 1');
  }
  if (diameter <= 0) {
    throw new Error('Diameter must be positive');
  }
  return efficiency * Math.PI * Math.pow(diameter / 2, 2);
}

/**
 * Conversion utilities for radar calculations
 */
export const convert = {
  /**
   * Convert linear value to decibels (dB)
   */
  linearToDb(value: number): number {
    if (value <= 0) {
      return -Infinity;
    }
    return 10 * Math.log10(value);
  },

  /**
   * Convert decibels (dB) to linear value
   */
  dbToLinear(valueDb: number): number {
    return Math.pow(10, valueDb / 10);
  },

  /**
   * Convert radians to degrees
   */
  radToDeg(radians: number): number {
    return radians * (180 / Math.PI);
  },

  /**
   * Convert degrees to radians
   */
  degToRad(degrees: number): number {
    return degrees * (Math.PI / 180);
  },

  /**
   * Convert meters to kilometers
   */
  metersToKm(meters: number): number {
    return meters / 1000;
  },

  /**
   * Convert kilometers to meters
   */
  kmToMeters(km: number): number {
    return km * 1000;
  },

  /**
   * Convert meters to nautical miles
   */
  metersToNauticalMiles(meters: number): number {
    return meters / 1852;
  },

  /**
   * Convert nautical miles to meters
   */
  nauticalMilesToMeters(nauticalMiles: number): number {
    return nauticalMiles * 1852;
  },

  /**
   * Convert meters to miles
   */
  metersToMiles(meters: number): number {
    return meters / 1609.34;
  },

  /**
   * Convert miles to meters
   */
  milesToMeters(miles: number): number {
    return miles * 1609.34;
  },

  /**
   * Convert hertz to gigahertz
   */
  hzToGhz(hz: number): number {
    return hz / 1e9;
  },

  /**
   * Convert gigahertz to hertz
   */
  ghzToHz(ghz: number): number {
    return ghz * 1e9;
  },

  /**
   * Convert watts to kilowatts
   */
  wattsToKw(watts: number): number {
    return watts / 1000;
  },

  /**
   * Convert kilowatts to watts
   */
  kwToWatts(kw: number): number {
    return kw * 1000;
  }
};

/**
 * Physical constants used in radar calculations
 */
export const constants = {
  /** Speed of light in m/s */
  SPEED_OF_LIGHT: 299792458,
  
  /** Boltzmann constant in J/K */
  BOLTZMANN: 1.380649e-23,
  
  /** Pi */
  PI: Math.PI
};
