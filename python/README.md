# Radar Range Equation
A basic toolbox for solving radar range equations

## Testing

After building and installing the package, you can run the test script to verify functionality:

```bash
python python/test_package.py
```

This test script verifies that:
- The package can be imported successfully
- Variables can be set dynamically (e.g., `c`, `f`, `lambda`)
- The `redefine_variable` function works correctly
- Calculations work as expected (e.g., `lambda = c/f`)

### Example Usage

```python
import radar_range_equation as RRE

# Set the speed of light (m/s)
RRE.vars.c = 3.00 * 10**8

# Set the frequency (Hz)
RRE.vars.f = 10

# Calculate and set wavelength (m)
# Note: 'lambda' is a reserved keyword in Python, so use setattr
setattr(RRE.vars, 'lambda', RRE.vars.c / RRE.vars.f)

# Print the wavelength
print(getattr(RRE.vars, 'lambda'))  # Output: 30000000.0
```
