#!/usr/bin/env python3
"""
Example demonstrating the LFM spectrogram generation capability.

This example shows how to generate a spectrogram for Linear Frequency Modulated
(LFM) radar pulses, showing the frequency sweep from 7850 MHz to 8150 MHz.

This matches the exact usage requested in the problem statement.
"""

import sys
sys.path.insert(0, 'python/src')  # For running from repo without installation
import radar_range_equation as RRE

# Example Usage - This generates the graph with the parameters from the problem statement
print("="*70)
print("LFM Spectrogram Example")
print("="*70)
print("\nGenerating spectrogram for LFM radar pulses...")
print("Parameters:")
print("  - Start frequency: 7850 MHz")
print("  - End frequency: 8150 MHz")
print("  - Pulse duration: 40 µs")
print("  - Pulse start times: [0, 120] µs")
print("  - Total duration: 180 µs")
print("  - Sampling frequency: 500 MHz")
print("  - Noise factor: 0.05")

# Generate the spectrogram
f, t, Sxx, fig = RRE.plot.generate_lfm_spectrogram(
    f_start_MHz=7850,
    f_end_MHz=8150,
    pulse_duration_us=40,
    pulse_start_times_us=[0, 120],
    total_duration_us=180
)

print(f"\n✓ Spectrogram generated successfully!")
print(f"  - Frequency bins: {len(f)}")
print(f"  - Time bins: {len(t)}")
print(f"  - Spectrogram shape: {Sxx.shape}")
print("\nThe spectrogram shows two LFM chirp pulses:")
print("  - First pulse: 0-40 µs, sweeping from 7850-8150 MHz")
print("  - Second pulse: 120-160 µs, sweeping from 7850-8150 MHz")
print("\n" + "="*70)

# Note: The plot will display automatically when show=True (default)
