#!/usr/bin/env python3
"""Test script for angle estimation functions.

This script validates the new angle estimation functions added to the
radar_range_equation package, including amplitude, phase, and time comparison methods.
"""

import sys
import math
from math import isclose


def test_angle_estimation():
    """Run tests for the angle estimation functions.

    Returns 0 on success, 1 on failure.
    """
    try:
        import radar_range_equation as RRE

        print("✓ Imported radar_range_equation")

        # Test that new functions are exported
        for name in ("v_phi", "calculate_Theta", "v_phi_full", "estimate_phi_hat",
                     "sigma_phi_amplitude", "sigma_phi_phase", "sigma_phi_time", "db_to_linear"):
            assert hasattr(RRE, name), f"Missing exported function: {name}"
        print("✓ All new functions are exported")

        # Test v_phi (Gaussian beam approximation)
        result = RRE.v_phi(5.0, 5.0, 0.0277)
        assert isclose(result, 1.0, rel_tol=1e-9), f"v_phi at center should be ~1.0, got {result}"
        print(f"✓ v_phi(5.0, 5.0, 0.0277) = {result}")

        # Test calculate_Theta
        theta_B = 10.0
        Theta = RRE.calculate_Theta(theta_B)
        expected_Theta = (4 * math.log(2)) / (theta_B**2)
        assert isclose(Theta, expected_Theta, rel_tol=1e-9), f"calculate_Theta({theta_B}) = {Theta}, expected {expected_Theta}"
        print(f"✓ calculate_Theta({theta_B}) = {Theta}")

        # Test v_phi_full
        result = RRE.v_phi_full(5.0, 5.0, 10.0)
        assert isclose(result, 1.0, rel_tol=1e-9), f"v_phi_full at center should be ~1.0, got {result}"
        print(f"✓ v_phi_full(5.0, 5.0, 10.0) = {result}")

        # Test estimate_phi_hat
        phi_hat = RRE.estimate_phi_hat(1.0, 1.0, 10.0, 5.0)
        expected = (1.0 / 1.0) * (10.0**2 / (8 * math.log(2) * 5.0))
        assert isclose(phi_hat, expected, rel_tol=1e-9), f"estimate_phi_hat = {phi_hat}, expected {expected}"
        print(f"✓ estimate_phi_hat(1.0, 1.0, 10.0, 5.0) = {phi_hat}")

        # Test db_to_linear
        S_N_dB = 10
        S_N_linear = RRE.db_to_linear(S_N_dB)
        expected_linear = 10.0
        assert isclose(S_N_linear, expected_linear, rel_tol=1e-9), f"db_to_linear({S_N_dB}) = {S_N_linear}, expected {expected_linear}"
        print(f"✓ db_to_linear({S_N_dB}) = {S_N_linear}")

        # Test sigma_phi_amplitude (from problem statement page 2)
        theta_B_1 = 10.0
        phi_s_1 = 5.0
        S_N_1 = 10.0
        sigma_phi_amp = RRE.sigma_phi_amplitude(theta_B_1, S_N_1, phi_s_1)
        # Expected: approx 0.8 degrees according to problem statement
        expected_amp = (theta_B_1**2 * math.sqrt(1 / S_N_1)) / (8 * math.sqrt(2) * phi_s_1 * math.log(2))
        assert isclose(sigma_phi_amp, expected_amp, rel_tol=1e-9), f"sigma_phi_amplitude = {sigma_phi_amp}, expected {expected_amp}"
        print(f"✓ sigma_phi_amplitude({theta_B_1}, {S_N_1}, {phi_s_1}) = {sigma_phi_amp:.4f} degrees")

        # Test calculation for target sigma_phi (from problem statement)
        target_sigma_phi = 0.5
        phi_s_calc = (theta_B_1**2 * math.sqrt(1 / S_N_1)) / (8 * math.sqrt(2) * target_sigma_phi * math.log(2))
        print(f"✓ Calculated phi_s for target sigma_phi of {target_sigma_phi} = {phi_s_calc:.4f} degrees")

        # Test sigma_phi_phase (from problem statement page 3)
        c = 3e8
        f_o = 10e9
        lambda_ = c / f_o  # 0.03 m
        S_N_dB_2 = 8
        S_N_linear_2 = RRE.db_to_linear(S_N_dB_2)  # approx 6.3
        d_1 = 2.0  # 2 meters
        sigma_phi_phase_val = RRE.sigma_phi_phase(lambda_, d_1, S_N_linear_2)
        expected_phase = (lambda_ / (2 * math.pi * d_1)) * math.sqrt(1 / S_N_linear_2)
        assert isclose(sigma_phi_phase_val, expected_phase, rel_tol=1e-9), f"sigma_phi_phase = {sigma_phi_phase_val}, expected {expected_phase}"
        print(f"✓ sigma_phi_phase({lambda_}, {d_1}, {S_N_linear_2:.2f}) = {sigma_phi_phase_val:.6f} radians")
        sigma_phi_phase_deg = sigma_phi_phase_val * (180 / math.pi)
        print(f"  = {sigma_phi_phase_deg:.4f} degrees")

        # Test calculation of d for target sigma_phi (from problem statement)
        target_sigma_phi_deg_2 = 0.5
        target_sigma_phi_rad_2 = target_sigma_phi_deg_2 * (math.pi / 180)
        d_calc_2 = (lambda_ / (2 * math.pi * target_sigma_phi_rad_2)) * math.sqrt(1 / S_N_linear_2)
        print(f"✓ Calculated d for target sigma_phi of {target_sigma_phi_deg_2} deg = {d_calc_2:.4f} meters")

        # Test sigma_phi_time (from problem statement page 3)
        B_1 = 200e6  # 200 MHz
        sigma_phi_time_val = RRE.sigma_phi_time(c, d_1, B_1)
        expected_time = c / (d_1 * B_1)
        assert isclose(sigma_phi_time_val, expected_time, rel_tol=1e-9), f"sigma_phi_time = {sigma_phi_time_val}, expected {expected_time}"
        print(f"✓ sigma_phi_time({c:.2e}, {d_1}, {B_1:.2e}) = {sigma_phi_time_val:.4f} radians")
        sigma_phi_time_deg = sigma_phi_time_val * (180 / math.pi)
        print(f"  = {sigma_phi_time_deg:.2f} degrees")

        # Test calculation of B for target sigma_phi (from problem statement)
        target_sigma_phi_deg_3 = 0.5
        target_sigma_phi_rad_3 = target_sigma_phi_deg_3 * (math.pi / 180)
        B_calc_2 = c / (d_1 * target_sigma_phi_rad_3)
        print(f"✓ Calculated B for target sigma_phi of {target_sigma_phi_deg_3} deg = {B_calc_2:.2e} Hz")
        B_calc_2_GHz = B_calc_2 / 1e9
        print(f"  = {B_calc_2_GHz:.2f} GHz")

        print("\n✓ All angle estimation tests passed")
        return 0

    except ImportError as e:
        print(f"✗ ImportError: {e}")
        return 1
    except AssertionError as e:
        print(f"✗ Assertion failed: {e}")
        return 1
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(test_angle_estimation())
