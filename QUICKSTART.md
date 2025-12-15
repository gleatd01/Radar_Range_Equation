# Quickstart Guide

Get started with Radar Range Equation in 5 minutes! This guide will walk you through installation and your first calculations.

## Choose Your Language

<details open>
<summary><b>Python</b></summary>

### Installation

```bash
pip install radar-range-equation
```

### Your First Calculation

```python
import radar_range_equation as RRE

# Set physical constants
RRE.vars.c = 3.0e8        # Speed of light (m/s)
RRE.vars.f = 10e9         # Frequency: 10 GHz
RRE.vars.P_t = 1000       # Transmit power: 1 kW
RRE.vars.sigma = 1.0      # Radar cross-section: 1 m²
RRE.vars.S_min = 1e-13    # Minimum detectable signal: 0.1 pW

# Calculate wavelength
wavelength = RRE.solve.wavelength()
print(f"Wavelength: {wavelength*100:.2f} cm")

# Calculate antenna effective aperture (60 ft diameter circular antenna)
RRE.vars.D = RRE.convert.ft_to_m(60)  
RRE.vars.eta = 0.6  # 60% efficiency
RRE.vars.A_e = RRE.solve.A_e_circ()
print(f"Effective aperture: {RRE.vars.A_e:.2f} m²")

# Calculate antenna gain
RRE.vars.G_t = RRE.solve.G_t()
RRE.vars.G_r = RRE.vars.G_t  # Same antenna for transmit and receive
print(f"Antenna gain: {RRE.vars.G_t:.2f} (linear)")
print(f"Antenna gain: {RRE.convert.lin_to_db(RRE.vars.G_t):.2f} dB")

# Calculate maximum detection range
R_max = RRE.solve.R_max()
print(f"Maximum range: {R_max/1000:.2f} km")

# Visualize a pulsed radar signal
RRE.plot.pulsed_radar_signal(
    amplitude=20,
    frequency=0.5e6,    # 0.5 MHz
    pulse_width=15e-6,  # 15 microseconds
    pri=50e-6,          # 50 microseconds
    num_pulses=3
)
```

### What's Next?

- Explore [example scripts](python/README.md#testing) in the `python/` directory
- Read the [full Python documentation](python/README.md)
- Check out [EQUATIONS.md](EQUATIONS.md) for all available equations
- See [SOLVERS.md](SOLVERS.md) for all solver functions

</details>

<details>
<summary><b>Flutter/Dart</b></summary>

### Installation

Add to your `pubspec.yaml`:

```yaml
dependencies:
  radar_range_equation:
    git:
      url: https://github.com/gleatd01/Radar_Range_Equation.git
      path: flutter
```

Then run:
```bash
dart pub get
```

### Your First Calculation

```dart
import 'package:radar_range_equation/radar_range_equation.dart';

void main() {
  // Calculate wavelength from frequency
  final frequency = 10e9; // 10 GHz
  final wavelength = calculateWavelength(frequency);
  print('Wavelength: ${(wavelength * 100).toStringAsFixed(2)} cm');

  // Calculate maximum radar range
  final maxRange = calculateMaxRange(
    transmitPower: 1000.0,      // 1 kW
    antennaGain: 1000.0,        // Linear gain
    wavelength: wavelength,
    radarCrossSection: 1.0,     // 1 m²
    minDetectableSignal: 1e-13, // 0.1 pW
  );
  print('Maximum range: ${(maxRange / 1000).toStringAsFixed(2)} km');

  // Calculate Doppler frequency shift
  final velocity = -100.0; // 100 m/s closing velocity
  final dopplerFreq = calculateDopplerFrequency(velocity, wavelength);
  print('Doppler shift: ${RadarConvert.hzToMhz(dopplerFreq).toStringAsFixed(3)} MHz');

  // Unit conversions
  print('10 dB = ${RadarConvert.dbToLinear(10).toStringAsFixed(2)} (linear)');
  print('1000 m = ${RadarConvert.metersToKilometers(1000)} km');
}
```

### What's Next?

- Run the [example program](flutter/example/example.dart)
- Read the [full Flutter documentation](flutter/README.md)
- Run tests with `dart test`

</details>

<details>
<summary><b>Rust</b></summary>

### Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
radar_range_equation = { git = "https://github.com/gleatd01/Radar_Range_Equation.git", path = "rust" }
```

### Your First Calculation

```rust
use radar_range_equation::*;

fn main() -> Result<(), RadarError> {
    // Calculate wavelength from frequency
    let frequency = 10e9; // 10 GHz
    let wavelength = calculate_wavelength(frequency)?;
    println!("Wavelength: {:.2} cm", wavelength * 100.0);

    // Calculate maximum radar range
    let max_range = calculate_max_range(
        1000.0,     // transmit power (W)
        1000.0,     // antenna gain (linear)
        wavelength, // wavelength (m)
        1.0,        // radar cross-section (m²)
        1e-13,      // min detectable signal (W)
    )?;
    println!("Maximum range: {:.2} km", max_range / 1000.0);

    // Calculate Doppler frequency shift
    let velocity = -100.0; // 100 m/s closing velocity
    let doppler_freq = calculate_doppler_frequency(velocity, wavelength)?;
    println!("Doppler shift: {:.3} MHz", convert::hz_to_mhz(doppler_freq));

    // Unit conversions
    println!("10 dB = {:.2} (linear)", convert::db_to_linear(10.0));
    println!("1000 m = {} km", convert::meters_to_kilometers(1000.0));

    Ok(())
}
```

### What's Next?

- Run the [example program](rust/examples/basic.rs) with `cargo run --example basic`
- Read the [full Rust documentation](rust/README.md)
- Run tests with `cargo test`

</details>

## Common Use Cases

### Calculate Detection Range

Determine how far your radar can detect a target:

**Python:**
```python
import radar_range_equation as RRE
RRE.vars.P_t = 1000  # 1 kW transmit power
RRE.vars.G_t = 1000  # Antenna gain
RRE.vars.G_r = 1000
RRE.vars.wavelength = 0.03  # 30 cm (1 GHz)
RRE.vars.sigma = 1.0  # 1 m² target
RRE.vars.S_min = 1e-13
range_km = RRE.solve.R_max() / 1000
print(f"Max range: {range_km:.2f} km")
```

### Doppler Shift Calculation

Calculate the frequency shift due to target motion:

**Python:**
```python
import radar_range_equation as RRE
wavelength = 0.03  # 30 cm
velocity = -50  # 50 m/s closing (negative = towards radar)
RRE.vars.wavelength = wavelength
RRE.vars.v = velocity
doppler_hz = RRE.solve.f_d()
print(f"Doppler shift: {doppler_hz:.2f} Hz")
```

### Unit Conversions

Convert between different units:

**Python:**
```python
import radar_range_equation as RRE

# Power conversions
linear = RRE.convert.db_to_lin(30)  # 30 dB to linear
db = RRE.convert.lin_to_db(1000)    # 1000 to dB

# Distance conversions
km = RRE.convert.m_to_km(5000)      # 5000 m to km
nm = RRE.convert.m_to_nmi(10000)    # 10 km to nautical miles

# Frequency conversions
mhz = RRE.convert.hz_to(1e9, 'MHz')  # 1 GHz to MHz
```

## Troubleshooting

### Python: "lambda is a reserved keyword"

The wavelength variable uses Python's reserved word `lambda`. Use `getattr`/`setattr`:

```python
import radar_range_equation as RRE

# Setting wavelength
setattr(RRE.vars, 'lambda', 0.03)

# Getting wavelength
wavelength = getattr(RRE.vars, 'lambda')
```

Or use the `wavelength` variable instead:
```python
RRE.vars.wavelength = 0.03
```

### Installation Issues

If you encounter installation issues:

1. **Check Python version**: Python 3.9+ required
   ```bash
   python --version
   ```

2. **Upgrade pip**:
   ```bash
   pip install --upgrade pip
   ```

3. **Install from source**:
   ```bash
   git clone https://github.com/gleatd01/Radar_Range_Equation.git
   cd Radar_Range_Equation
   pip install -e .
   ```

## Getting Help

- **Documentation**: Check [README.md](README.md) and language-specific docs
- **Examples**: Browse example files in each language directory
- **Issues**: Report bugs or ask questions on [GitHub Issues](https://github.com/gleatd01/Radar_Range_Equation/issues)
- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md)

## Next Steps

1. **Explore Examples**: Each language has example files demonstrating various features
2. **Read Full Docs**: Language-specific READMEs have comprehensive API documentation
3. **Review Equations**: [EQUATIONS.md](EQUATIONS.md) lists all 72 available equations
4. **Try Solvers**: [SOLVERS.md](SOLVERS.md) documents all 76 solver functions
5. **Contribute**: Help improve the library by contributing code or documentation

Happy calculating! 🎯📡
