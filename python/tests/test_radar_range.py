"""Tests for radar range equation calculations."""

import pytest
import math
from radar_range_equation import calculate_max_range, calculate_received_power


class TestCalculateMaxRange:
    """Tests for the calculate_max_range function."""
    
    def test_basic_calculation(self):
        """Test basic range calculation with typical values."""
        max_range = calculate_max_range(
            transmit_power=1000,
            antenna_gain=1000,
            wavelength=0.03,
            radar_cross_section=1.0,
            min_detectable_signal=1e-13
        )
        assert max_range > 0
        assert isinstance(max_range, float)
    
    def test_negative_power_raises_error(self):
        """Test that negative transmit power raises ValueError."""
        with pytest.raises(ValueError):
            calculate_max_range(
                transmit_power=-1000,
                antenna_gain=1000,
                wavelength=0.03,
                radar_cross_section=1.0,
                min_detectable_signal=1e-13
            )
    
    def test_zero_wavelength_raises_error(self):
        """Test that zero wavelength raises ValueError."""
        with pytest.raises(ValueError):
            calculate_max_range(
                transmit_power=1000,
                antenna_gain=1000,
                wavelength=0,
                radar_cross_section=1.0,
                min_detectable_signal=1e-13
            )
    
    def test_larger_power_increases_range(self):
        """Test that larger transmit power increases range."""
        range_1 = calculate_max_range(
            transmit_power=1000,
            antenna_gain=1000,
            wavelength=0.03,
            radar_cross_section=1.0,
            min_detectable_signal=1e-13
        )
        range_2 = calculate_max_range(
            transmit_power=2000,
            antenna_gain=1000,
            wavelength=0.03,
            radar_cross_section=1.0,
            min_detectable_signal=1e-13
        )
        assert range_2 > range_1


class TestCalculateReceivedPower:
    """Tests for the calculate_received_power function."""
    
    def test_basic_calculation(self):
        """Test basic received power calculation."""
        power = calculate_received_power(
            transmit_power=1000,
            antenna_gain=1000,
            wavelength=0.03,
            radar_cross_section=1.0,
            range_m=10000
        )
        assert power > 0
        assert isinstance(power, float)
    
    def test_negative_range_raises_error(self):
        """Test that negative range raises ValueError."""
        with pytest.raises(ValueError):
            calculate_received_power(
                transmit_power=1000,
                antenna_gain=1000,
                wavelength=0.03,
                radar_cross_section=1.0,
                range_m=-1000
            )
    
    def test_inverse_fourth_power_law(self):
        """Test that received power follows inverse fourth power law."""
        power_1 = calculate_received_power(
            transmit_power=1000,
            antenna_gain=1000,
            wavelength=0.03,
            radar_cross_section=1.0,
            range_m=10000
        )
        power_2 = calculate_received_power(
            transmit_power=1000,
            antenna_gain=1000,
            wavelength=0.03,
            radar_cross_section=1.0,
            range_m=20000
        )
        # Power at 2x range should be 1/16 (2^4) of power at 1x range
        assert abs(power_2 * 16 - power_1) < power_1 * 0.01  # Within 1%
    
    def test_consistency_with_max_range(self):
        """Test that received power at max range equals min detectable signal."""
        transmit_power = 1000
        antenna_gain = 1000
        wavelength = 0.03
        radar_cross_section = 1.0
        min_detectable_signal = 1e-13
        
        max_range = calculate_max_range(
            transmit_power=transmit_power,
            antenna_gain=antenna_gain,
            wavelength=wavelength,
            radar_cross_section=radar_cross_section,
            min_detectable_signal=min_detectable_signal
        )
        
        received_power = calculate_received_power(
            transmit_power=transmit_power,
            antenna_gain=antenna_gain,
            wavelength=wavelength,
            radar_cross_section=radar_cross_section,
            range_m=max_range
        )
        
        # Should be very close to min_detectable_signal
        assert abs(received_power - min_detectable_signal) / min_detectable_signal < 0.01
