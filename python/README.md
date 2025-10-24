# Python Package

This directory contains the Python implementation of the Radar Range Equation tools.

## Installation

### From Source

```bash
cd python
pip install -e .
```

### For Development

```bash
cd python
pip install -e ".[dev]"
```

## Usage

```python
from radar_range_equation import calculate_max_range, calculate_received_power

# Calculate maximum detection range
max_range = calculate_max_range(
    transmit_power=1000,           # 1 kW
    antenna_gain=1000,              # 30 dBi (linear scale)
    wavelength=0.03,                # 10 GHz
    radar_cross_section=1.0,        # 1 m²
    min_detectable_signal=1e-13     # -100 dBm
)

print(f"Maximum range: {max_range:.2f} meters")

# Calculate received power at a specific range
received_power = calculate_received_power(
    transmit_power=1000,
    antenna_gain=1000,
    wavelength=0.03,
    radar_cross_section=1.0,
    range_m=10000  # 10 km
)

print(f"Received power: {received_power:.2e} Watts")
```

## Testing

```bash
cd python
pytest
```

## Linting and Formatting

```bash
# Format code
black src tests

# Lint code
ruff check src tests

# Type checking
mypy src
```

## Project Structure

```
python/
├── src/
│   └── radar_range_equation/
│       ├── __init__.py
│       ├── radar_range.py
│       └── py.typed
├── tests/
│   └── test_radar_range.py
├── examples/
│   └── basic_example.py
└── pyproject.toml
```
