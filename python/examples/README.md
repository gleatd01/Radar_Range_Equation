# Python Examples

This directory contains example scripts demonstrating various features of the Radar Range Equation Python package.

## Available Examples

### Basic Examples

Located in the `python/` directory (one level up):

- **`example.py`** - Basic usage example showing fundamental calculations
- **`example_plotting.py`** - Comprehensive plotting examples (pulsed radar, CW Doppler, CWFM, pulse compression, etc.)
- **`example_analysis.py`** - Analysis features including pulse timing and integration gains
- **`example_angle_estimation.py`** - Direction finding and angle estimation examples
- **`example_interactive_solver.py`** - Interactive solver demonstrations
- **`example_tactical_scenario.py`** - Tactical radar scenario examples
- **`example_problem_statement_tactical.py`** - Problem statement and solution examples

## Running Examples

### Prerequisites

Make sure you have installed the package:

```bash
pip install radar-range-equation
```

Or install in development mode from the repository root:

```bash
pip install -e .
```

### Run an Example

```bash
# From the repository root
python python/example.py
python python/example_plotting.py
python python/example_analysis.py
```

### Example Output Location

Some examples generate plots that are saved in the `python/examples/plots/` directory.

## Example Categories

### 1. Basic Calculations (`example.py`)
- Setting variables
- Basic radar equation calculations
- Wavelength calculations
- Maximum range computation

### 2. Visualization (`example_plotting.py`)
- Pulsed radar signals
- CW Doppler signals
- CWFM radar signals
- Pulse compression visualization
- Range profiles
- Doppler spectra

### 3. Analysis (`example_analysis.py`)
- Pulse timing analysis
- Integration gain calculations
- Jammer analysis
- Detection probability

### 4. Direction Finding (`example_angle_estimation.py`)
- Monopulse angle estimation
- Angle accuracy calculations
- Beamwidth effects

### 5. Interactive Solvers (`example_interactive_solver.py`)
- Solving for different variables
- Equation manipulation
- Multi-variable solutions

### 6. Tactical Scenarios (`example_tactical_scenario.py`)
- Real-world scenarios
- Multiple target tracking
- Countermeasure effects

## Creating Your Own Examples

To create your own example:

1. Import the package:
   ```python
   import radar_range_equation as RRE
   ```

2. Set your variables using `RRE.vars`:
   ```python
   RRE.vars.P_t = 1000  # Transmit power in watts
   RRE.vars.f = 10e9    # Frequency in Hz
   ```

3. Use solver functions:
   ```python
   max_range = RRE.solve.R_max()
   ```

4. Visualize results:
   ```python
   RRE.plot.pulsed_radar_signal()
   ```

## Tips

- **Variable Names**: See `python/README.md` for complete list of available variables
- **Solver Functions**: Check `SOLVERS.md` for all 76 solver functions
- **Equations**: Refer to `EQUATIONS.md` for the 72 available equations
- **Lambda Variable**: Remember that `lambda` is a Python keyword, use `getattr(RRE.vars, 'lambda')` or use `RRE.vars.wavelength` instead

## Need Help?

- Check the main [README.md](../../README.md)
- See the [QUICKSTART.md](../../QUICKSTART.md) guide
- Read [CONTRIBUTING.md](../../CONTRIBUTING.md)
- Open an issue on [GitHub](https://github.com/gleatd01/Radar_Range_Equation/issues)

## Contributing Examples

Have a useful example? We'd love to include it! 

1. Create your example script
2. Add comments explaining what it does
3. Test it thoroughly
4. Submit a pull request
5. Update this README

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.
