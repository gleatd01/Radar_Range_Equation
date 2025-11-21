/**
 * Basic usage example for the radar-range-equation package
 */

import {
  calculateMaxRange,
  calculateReceivedPower,
  calculateWavelength,
  calculateAntennaGain,
  calculateEffectiveApertureCirc,
  convert,
  constants
} from '../src/index';

console.log('=== Radar Range Equation Examples ===\n');

// Example 1: Calculate maximum detection range
console.log('Example 1: Calculate Maximum Detection Range');
console.log('-------------------------------------------');

const transmitPower = 1000; // W (1 kW)
const antennaGain = 1000; // dimensionless
const wavelength = 0.03; // m (10 GHz)
const radarCrossSection = 1.0; // m²
const minDetectableSignal = 1e-13; // W

const maxRange = calculateMaxRange(
  transmitPower,
  antennaGain,
  wavelength,
  radarCrossSection,
  minDetectableSignal
);

console.log(`Transmit Power: ${transmitPower} W`);
console.log(`Antenna Gain: ${antennaGain}`);
console.log(`Wavelength: ${wavelength} m`);
console.log(`Radar Cross Section: ${radarCrossSection} m²`);
console.log(`Minimum Detectable Signal: ${minDetectableSignal} W`);
console.log(`\nMaximum Detection Range: ${maxRange.toFixed(2)} meters`);
console.log(`                         ${convert.metersToKm(maxRange).toFixed(2)} km`);
console.log(`                         ${convert.metersToNauticalMiles(maxRange).toFixed(2)} nautical miles\n`);

// Example 2: Calculate received power at a specific range
console.log('Example 2: Calculate Received Power at Range');
console.log('-------------------------------------------');

const range = 10000; // m (10 km)
const receivedPower = calculateReceivedPower(
  transmitPower,
  antennaGain,
  wavelength,
  radarCrossSection,
  range
);

console.log(`Range to Target: ${range} meters (${convert.metersToKm(range)} km)`);
console.log(`Received Power: ${receivedPower.toExponential(3)} W`);
console.log(`Received Power: ${convert.linearToDb(receivedPower).toFixed(2)} dBW\n`);

// Example 3: Calculate wavelength from frequency
console.log('Example 3: Calculate Wavelength from Frequency');
console.log('-------------------------------------------');

const frequency = 10e9; // 10 GHz
const calculatedWavelength = calculateWavelength(frequency);

console.log(`Frequency: ${convert.hzToGhz(frequency)} GHz`);
console.log(`Wavelength: ${calculatedWavelength.toFixed(4)} meters`);
console.log(`            ${(calculatedWavelength * 100).toFixed(2)} cm\n`);

// Example 4: Calculate antenna gain from physical parameters
console.log('Example 4: Calculate Antenna Gain');
console.log('-------------------------------------------');

const efficiency = 0.6; // 60% efficiency
const diameter = 2.0; // 2 meter diameter
const effectiveAperture = calculateEffectiveApertureCirc(efficiency, diameter);
const gain = calculateAntennaGain(effectiveAperture, calculatedWavelength);

console.log(`Antenna Diameter: ${diameter} m`);
console.log(`Antenna Efficiency: ${(efficiency * 100).toFixed(0)}%`);
console.log(`Effective Aperture: ${effectiveAperture.toFixed(3)} m²`);
console.log(`Antenna Gain: ${gain.toFixed(2)}`);
console.log(`Antenna Gain: ${convert.linearToDb(gain).toFixed(2)} dBi\n`);

// Example 5: Compare ranges at different frequencies
console.log('Example 5: Range vs Frequency Comparison');
console.log('-------------------------------------------');

const frequencies = [1e9, 3e9, 10e9, 35e9]; // 1, 3, 10, 35 GHz
console.log('Frequency (GHz) | Wavelength (cm) | Max Range (km)');
console.log('----------------|-----------------|---------------');

frequencies.forEach(freq => {
  const wl = calculateWavelength(freq);
  const rng = calculateMaxRange(transmitPower, antennaGain, wl, radarCrossSection, minDetectableSignal);
  console.log(
    `${convert.hzToGhz(freq).toString().padEnd(15)} | ` +
    `${(wl * 100).toFixed(2).padEnd(15)} | ` +
    `${convert.metersToKm(rng).toFixed(2)}`
  );
});

console.log('\n=== All Examples Complete ===');
