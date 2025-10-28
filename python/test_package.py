#!/usr/bin/env python3
"""
Test script to verify the Radar_Range_Equation package after installation.

This script tests the basic functionality of the package by:
1. Setting the speed of light (c)
2. Setting the frequency (f)
3. Calculating and setting the wavelength (lambda = c/f)
4. Printing the result
"""

import sys


def test_package():
    """Test the basic functionality of the radar_range_equation package."""
    try:
        # Import the package
        import radar_range_equation as RRE
        print("✓ Package imported successfully")

        # Test 1: Set c (speed of light in m/s)
        RRE.vars.c = 3.00 * 10**8
        assert hasattr(RRE.vars, 'c'), "Failed to set c"
        assert RRE.vars.c == 3.00 * 10**8, f"Expected c=3.00e8, got {RRE.vars.c}"
        print(f"✓ RRE.vars.c = {RRE.vars.c}")

        # Test 2: Set f (frequency in Hz)
        RRE.vars.f = 10
        assert hasattr(RRE.vars, 'f'), "Failed to set f"
        assert RRE.vars.f == 10, f"Expected f=10, got {RRE.vars.f}"
        print(f"✓ RRE.vars.f = {RRE.vars.f}")

        # Test 3: Calculate and set lambda (wavelength in meters)
        # Note: 'lambda' is a reserved keyword, so we use setattr
        RRE.vars.lambda_ = 1
        assert hasattr(RRE.vars, 'lambda_'), "Failed to set lambda"
        print(f"✓ RRE.vars.lambda_ = {RRE.vars.lambda_}")
        
        # Test 4: Use the redefine_variable function
        new_wavelength = 0.03
        RRE.redefine_variable('lambda_', new_wavelength)
        lambda_value = getattr(RRE.vars, 'lambda_')
        assert lambda_value == new_wavelength, f"Expected lambda={new_wavelength}, got {lambda_value}"
        print(f"✓ RRE.redefine_variable('lambda', {new_wavelength}) works correctly")
        print(f"✓ RRE.vars.lambda_ = {lambda_value}")

        print("\n✓ All tests passed!")

        RRE.vars.f = 430*10**6  # Frequency in Hz
        RRE.vars.eta = .6
        RRE.vars.D = RRE.convert_ft_to_m(60)  # Antenna diameter in meters
        RRE.vars.R = 3.844*10**8  # Range in meters (example: distance to the Moon)
        RRE.vars.sigma = 6.64*10**11
        RRE.vars.A_e = RRE.solve.A_e_circ()
        RRE.vars.S_min = 1.5*10**-16
        RRE.vars.wavelength = RRE.solve.wavelength()
        pprint(RRE.vars.A_e)
        RRE.vars.G_t = RRE.solve.G_t()
        pprint(RRE.vars.G_t)
        pprint(RRE.equations.P_t)
        RRE.vars.P_t = RRE.solve.P_t()
        pprint(RRE.vars.P_t)
    
        RRE.vars.A=1.738*10**6  # Antenna radius in meters
        RRE.vars.sigma = RRE.solve.sigma_sphere()
        pprint(RRE.vars.sigma)
    
        vars.sigma = 6.64*10**11
        vars.A = solve.A_sphere()
        pprint(vars.A)
    
        return 0

    except ImportError as e:
        print(f"✗ Failed to import package: {e}")
        return 1
    except AssertionError as e:
        print(f"✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(test_package())
