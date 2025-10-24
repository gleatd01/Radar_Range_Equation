"""Core radar range equation calculations."""

import math
from typing import Optional


def calculate_max_range(
    transmit_power: float,
    antenna_gain: float,
    wavelength: float,
    radar_cross_section: float,
    min_detectable_signal: float,
) -> float:
    """Calculate the maximum detection range of a radar system.

    Uses the radar range equation to determine the maximum range at which
    a target can be detected.

    Args:
        transmit_power: Transmitted power in Watts
        antenna_gain: Antenna gain (dimensionless, linear scale)
        wavelength: Radar wavelength in meters
        radar_cross_section: Target radar cross section in square meters
        min_detectable_signal: Minimum detectable signal power in Watts

    Returns:
        Maximum range in meters

    Raises:
        ValueError: If any input parameters are non-positive

    Example:
        >>> range_m = calculate_max_range(
        ...     transmit_power=1000,      # 1 kW
        ...     antenna_gain=1000,         # 30 dBi
        ...     wavelength=0.03,           # 10 GHz
        ...     radar_cross_section=1.0,   # 1 m²
        ...     min_detectable_signal=1e-13
        ... )
    """
    if any(x <= 0 for x in [transmit_power, antenna_gain, wavelength, 
                             radar_cross_section, min_detectable_signal]):
        raise ValueError("All parameters must be positive")

    numerator = (transmit_power * antenna_gain**2 * wavelength**2 * 
                 radar_cross_section)
    denominator = (4 * math.pi)**3 * min_detectable_signal
    
    max_range = (numerator / denominator) ** 0.25
    
    return max_range


def calculate_received_power(
    transmit_power: float,
    antenna_gain: float,
    wavelength: float,
    radar_cross_section: float,
    range_m: float,
) -> float:
    """Calculate the received power from a target at a given range.

    Args:
        transmit_power: Transmitted power in Watts
        antenna_gain: Antenna gain (dimensionless, linear scale)
        wavelength: Radar wavelength in meters
        radar_cross_section: Target radar cross section in square meters
        range_m: Range to target in meters

    Returns:
        Received power in Watts

    Raises:
        ValueError: If any input parameters are non-positive

    Example:
        >>> power = calculate_received_power(
        ...     transmit_power=1000,
        ...     antenna_gain=1000,
        ...     wavelength=0.03,
        ...     radar_cross_section=1.0,
        ...     range_m=10000
        ... )
    """
    if any(x <= 0 for x in [transmit_power, antenna_gain, wavelength, 
                             radar_cross_section, range_m]):
        raise ValueError("All parameters must be positive")

    numerator = (transmit_power * antenna_gain**2 * wavelength**2 * 
                 radar_cross_section)
    denominator = (4 * math.pi)**3 * range_m**4
    
    received_power = numerator / denominator
    
    return received_power
