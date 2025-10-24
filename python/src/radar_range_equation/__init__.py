"""Radar Range Equation Calculator.

This package provides tools for calculating radar range equations
and related radar system parameters.
"""

__version__ = "0.1.0"

from .radar_range import calculate_max_range, calculate_received_power

__all__ = ["calculate_max_range", "calculate_received_power"]
