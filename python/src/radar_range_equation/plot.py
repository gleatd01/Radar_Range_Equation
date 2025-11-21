"""Plotting utilities for visualizing radar signals and equations.

This module provides functions to visualize various radar signal types including:
    - Pulsed radar transmit signals
    - CW (Continuous Wave) signals
    - CWFM (Continuous Wave Frequency Modulated) signals
    - Range profiles
    - Doppler spectra
    - Tactical scenarios (radar, target, and jammer positions)
    - LFM (Linear Frequency Modulated) spectrograms

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.plot.pulsed_radar_signal(
    ...     amplitude=20,
    ...     frequency=0.5e6,
    ...     pulse_width=15e-6,
    ...     pri=50e-6,
    ...     num_pulses=3,
    ...     time_span=150e-6
    ... )
    >>> 
    >>> # Generate LFM spectrogram
    >>> RRE.plot.generate_lfm_spectrogram(
    ...     f_start_MHz=7850,
    ...     f_end_MHz=8150,
    ...     pulse_duration_us=40
    ... )
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List
from scipy.signal import spectrogram


def pulsed_radar_signal(
    amplitude: float = 20,
    frequency: float = 0.5e6,
    pulse_width: float = 15e-6,
    pri: float = 50e-6,
    num_pulses: int = 3,
    time_span: float = 150e-6,
    num_points: int = 10000,
    figsize: Tuple[float, float] = (10, 5),
    show: bool = True
) -> Tuple[np.ndarray, np.ndarray, plt.Figure]:
    """Plot a pulsed radar transmit signal.
    
    Creates a visualization of a pulsed radar signal showing the carrier wave
    modulated by rectangular pulses. The signal consists of periodic bursts
    of a sinusoidal carrier wave.
    
    Args:
        amplitude: Peak amplitude of the carrier signal (default: 20)
        frequency: Carrier frequency in Hz (default: 0.5e6 = 0.5 MHz)
        pulse_width: Duration of each pulse in seconds (default: 15e-6 = 15 µs)
        pri: Pulse Repetition Interval in seconds (default: 50e-6 = 50 µs)
        num_pulses: Number of pulses to display (default: 3)
        time_span: Total time span to display in seconds (default: 150e-6 = 150 µs)
        num_points: Number of sample points for smooth plotting (default: 10000)
        figsize: Figure size as (width, height) in inches (default: (10, 5))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        Tuple containing:
            - t (np.ndarray): Time vector in microseconds
            - signal (np.ndarray): The modulated signal values
            - fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> # Create a pulsed radar signal with default parameters
        >>> t, signal, fig = RRE.plot.pulsed_radar_signal()
        >>> 
        >>> # Customize the signal parameters
        >>> t, signal, fig = RRE.plot.pulsed_radar_signal(
        ...     amplitude=30,
        ...     frequency=1e6,
        ...     pulse_width=10e-6,
        ...     pri=40e-6,
        ...     num_pulses=5,
        ...     show=False
        ... )
        >>> # Save to file
        >>> fig.savefig('pulsed_radar.png')
    
    Notes:
        - The carrier signal is a cosine wave with the specified frequency
        - Pulses are rectangular (on/off) modulation of the carrier
        - Time axis is displayed in microseconds for readability
        - Angular frequency ω = 2π * frequency
    """
    # Create high-resolution time vector (in seconds)
    t_sec = np.linspace(0, time_span, num_points)
    
    # Convert to microseconds for display
    t = t_sec * 1e6  # Convert to microseconds
    
    # Create the carrier signal (cosine wave)
    # Angular frequency: ω = 2π * f
    carrier = amplitude * np.cos(2 * np.pi * frequency * t_sec)
    
    # Create the pulse mask (on/off pattern)
    mask = np.zeros_like(t_sec, dtype=bool)
    
    # Add each pulse to the mask
    for i in range(num_pulses):
        pulse_start = i * pri
        pulse_end = pulse_start + pulse_width
        pulse_mask = (t_sec >= pulse_start) & (t_sec < pulse_end)
        mask = mask | pulse_mask
    
    # Apply the pulse mask to create the final signal
    # NumPy treats True as 1 and False as 0 in calculations
    signal = carrier * mask
    
    # Create the plot
    fig = plt.figure(figsize=figsize)
    plt.plot(t, signal, linewidth=1)
    
    # Customize the plot
    plt.xlabel('Time ($\\mu$s)', fontsize=12)
    plt.ylabel('Transmit signal', fontsize=12)
    plt.title('Pulsed Radar Transmit Signal', fontsize=14)
    
    # Set axis limits and ticks
    plt.xlim(0, time_span * 1e6)  # Convert to microseconds
    # Create evenly spaced x-axis ticks
    num_x_ticks = int(time_span * 1e6 / 10) + 1
    plt.xticks(np.linspace(0, time_span * 1e6, num_x_ticks))
    
    # Set y-axis limits with padding
    y_max = amplitude * 1.1
    plt.ylim(-y_max, y_max)
    # Create y-axis ticks
    y_tick_step = amplitude / 4  # Create 5 ticks per side
    plt.yticks(np.arange(-amplitude, amplitude + y_tick_step, y_tick_step))
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Display the plot if requested
    if show:
        plt.show()
    
    return t, signal, fig


def cw_doppler_signal(
    amplitude: float = 1.0,
    carrier_freq: float = 10e9,
    doppler_shift: float = 6667,
    time_span: float = 1e-3,
    num_points: int = 10000,
    figsize: Tuple[float, float] = (10, 6),
    show: bool = True
) -> Tuple[np.ndarray, np.ndarray, plt.Figure]:
    """Plot a Continuous Wave (CW) Doppler radar signal.
    
    Visualizes the transmitted and received signals for a CW Doppler radar,
    showing the frequency shift caused by target motion.
    
    Args:
        amplitude: Signal amplitude (default: 1.0)
        carrier_freq: Carrier frequency in Hz (default: 10e9 = 10 GHz)
        doppler_shift: Doppler frequency shift in Hz (default: 6667 Hz)
        time_span: Time duration to display in seconds (default: 1e-3 = 1 ms)
        num_points: Number of sample points (default: 10000)
        figsize: Figure size as (width, height) in inches (default: (10, 6))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        Tuple containing:
            - t (np.ndarray): Time vector in milliseconds
            - signals (dict): Dictionary with 'tx' and 'rx' signal arrays
            - fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> t, signals, fig = RRE.plot.cw_doppler_signal(
        ...     carrier_freq=10e9,
        ...     doppler_shift=5000
        ... )
    """
    # Create time vector (in seconds)
    t_sec = np.linspace(0, time_span, num_points)
    
    # Convert to milliseconds for display
    t = t_sec * 1e3
    
    # Create transmitted signal (carrier only)
    # For visualization, we'll show a lower frequency representation
    # since the actual carrier is too fast to visualize
    vis_freq = 1e4  # Visualization frequency (10 kHz)
    tx_signal = amplitude * np.cos(2 * np.pi * vis_freq * t_sec)
    
    # Create received signal (carrier + Doppler shift)
    rx_signal = amplitude * np.cos(2 * np.pi * (vis_freq + doppler_shift) * t_sec)
    
    # Create the plot with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    
    # Plot transmitted signal
    ax1.plot(t, tx_signal, 'b-', linewidth=1)
    ax1.set_ylabel('Transmitted Signal', fontsize=11)
    ax1.set_title('CW Doppler Radar Signals', fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_ylim(-amplitude * 1.2, amplitude * 1.2)
    
    # Plot received signal
    ax2.plot(t, rx_signal, 'r-', linewidth=1)
    ax2.set_xlabel('Time (ms)', fontsize=12)
    ax2.set_ylabel('Received Signal', fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_ylim(-amplitude * 1.2, amplitude * 1.2)
    
    plt.tight_layout()
    
    if show:
        plt.show()
    
    signals = {'tx': tx_signal, 'rx': rx_signal}
    return t, signals, fig


def cwfm_signal(
    amplitude: float = 1.0,
    carrier_freq: float = 35e9,
    modulation_freq: float = 100,
    frequency_deviation: float = 30e6,
    time_span: float = 0.02,
    num_points: int = 10000,
    figsize: Tuple[float, float] = (10, 6),
    show: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, plt.Figure]:
    """Plot a Continuous Wave Frequency Modulated (CWFM) radar signal.
    
    Visualizes the transmitted CWFM signal and its frequency modulation pattern.
    
    Args:
        amplitude: Signal amplitude (default: 1.0)
        carrier_freq: Carrier frequency in Hz (default: 35e9 = 35 GHz)
        modulation_freq: Modulation frequency in Hz (default: 100 Hz)
        frequency_deviation: Peak frequency deviation in Hz (default: 30e6 = 30 MHz)
        time_span: Time duration in seconds (default: 0.02 = 20 ms)
        num_points: Number of sample points (default: 10000)
        figsize: Figure size as (width, height) in inches (default: (10, 6))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        Tuple containing:
            - t (np.ndarray): Time vector in milliseconds
            - signal (np.ndarray): The modulated signal
            - freq_modulation (np.ndarray): Frequency modulation pattern
            - fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> t, signal, freq_mod, fig = RRE.plot.cwfm_signal(
        ...     modulation_freq=200,
        ...     frequency_deviation=50e6
        ... )
    """
    # Create time vector (in seconds)
    t_sec = np.linspace(0, time_span, num_points)
    
    # Convert to milliseconds for display
    t = t_sec * 1e3
    
    # Create frequency modulation (triangular wave)
    freq_modulation = frequency_deviation * np.sin(2 * np.pi * modulation_freq * t_sec)
    
    # Create phase modulation (integral of frequency modulation)
    phase = 2 * np.pi * np.cumsum(freq_modulation) * (time_span / num_points)
    
    # Create the CWFM signal
    # For visualization, use a lower carrier frequency
    vis_freq = 1e3  # 1 kHz for visualization
    signal = amplitude * np.cos(2 * np.pi * vis_freq * t_sec + phase / 1e4)
    
    # Create the plot with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    
    # Plot frequency modulation
    ax1.plot(t, freq_modulation / 1e6, 'g-', linewidth=2)
    ax1.set_ylabel('Frequency Deviation (MHz)', fontsize=11)
    ax1.set_title('CWFM Radar Signal', fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    
    # Plot modulated signal
    ax2.plot(t, signal, 'b-', linewidth=1)
    ax2.set_xlabel('Time (ms)', fontsize=12)
    ax2.set_ylabel('Transmit Signal', fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_ylim(-amplitude * 1.2, amplitude * 1.2)
    
    plt.tight_layout()
    
    if show:
        plt.show()
    
    return t, signal, freq_modulation, fig


def pulse_compression_signal(
    pulse_width: float = 20e-6,
    bandwidth: float = 200e6,
    amplitude: float = 1.0,
    num_points: int = 10000,
    figsize: Tuple[float, float] = (10, 8),
    show: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, plt.Figure]:
    """Plot pulse compression radar signals (before and after compression).
    
    Visualizes the transmitted chirp pulse and the compressed pulse after
    matched filtering.
    
    Args:
        pulse_width: Pulse duration in seconds (default: 20e-6 = 20 µs)
        bandwidth: Chirp bandwidth in Hz (default: 200e6 = 200 MHz)
        amplitude: Signal amplitude (default: 1.0)
        num_points: Number of sample points (default: 10000)
        figsize: Figure size as (width, height) in inches (default: (10, 8))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        Tuple containing:
            - t (np.ndarray): Time vector in microseconds
            - chirp (np.ndarray): Transmitted chirp pulse
            - compressed (np.ndarray): Compressed pulse
            - fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> t, chirp, compressed, fig = RRE.plot.pulse_compression_signal(
        ...     pulse_width=10e-6,
        ...     bandwidth=100e6
        ... )
    """
    # Create time vector (in seconds)
    t_sec = np.linspace(-pulse_width, pulse_width, num_points)
    
    # Convert to microseconds for display
    t = t_sec * 1e6
    
    # Calculate chirp rate
    chirp_rate = bandwidth / pulse_width
    
    # Create transmitted chirp pulse
    chirp = np.zeros_like(t_sec)
    pulse_mask = (t_sec >= -pulse_width/2) & (t_sec <= pulse_width/2)
    chirp[pulse_mask] = amplitude * np.cos(2 * np.pi * chirp_rate * 0.5 * t_sec[pulse_mask]**2)
    
    # Create compressed pulse (sinc function approximation)
    # Time constant for compressed pulse
    t_c = 1 / bandwidth
    # Scale factor: divide by 2 * num_points for normalization from integration approximation
    compressed = amplitude * bandwidth * pulse_width * np.sinc(t_sec / t_c) / (2 * num_points)
    
    # Normalize compressed pulse to match input amplitude
    compressed = compressed / np.max(np.abs(compressed)) * amplitude
    
    # Create the plot with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize)
    
    # Plot transmitted chirp
    ax1.plot(t, chirp, 'b-', linewidth=1)
    ax1.set_ylabel('Amplitude', fontsize=11)
    ax1.set_title('Pulse Compression Radar Signal', fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_ylim(-amplitude * 1.2, amplitude * 1.2)
    ax1.text(0.02, 0.95, f'Transmitted Chirp\n(τ = {pulse_width*1e6:.1f} µs, B = {bandwidth/1e6:.0f} MHz)',
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Plot instantaneous frequency
    inst_freq = np.zeros_like(t_sec)
    inst_freq[pulse_mask] = chirp_rate * t_sec[pulse_mask]
    ax2.plot(t, inst_freq / 1e6, 'g-', linewidth=2)
    ax2.set_ylabel('Frequency (MHz)', fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax2.text(0.02, 0.95, 'Instantaneous Frequency',
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    # Plot compressed pulse
    ax3.plot(t, compressed, 'r-', linewidth=1.5)
    ax3.set_xlabel('Time (µs)', fontsize=12)
    ax3.set_ylabel('Amplitude', fontsize=11)
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.set_ylim(-amplitude * 1.2, amplitude * 1.2)
    pcr = pulse_width * bandwidth
    ax3.text(0.02, 0.95, f'Compressed Pulse\n(PCR = {pcr:.0f})',
             transform=ax3.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    
    plt.tight_layout()
    
    if show:
        plt.show()
    
    return t, chirp, compressed, fig


def range_profile(
    ranges: List[float],
    amplitudes: List[float],
    max_range: float = 300,
    figsize: Tuple[float, float] = (10, 5),
    show: bool = True
) -> plt.Figure:
    """Plot a radar range profile showing target detections.
    
    Creates a bar chart or stem plot showing detected targets at various ranges
    with their corresponding amplitudes.
    
    Args:
        ranges: List of target ranges in meters
        amplitudes: List of corresponding signal amplitudes
        max_range: Maximum range to display in meters (default: 300)
        figsize: Figure size as (width, height) in inches (default: (10, 5))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> ranges = [50, 100, 150, 200]
        >>> amplitudes = [0.8, 1.0, 0.6, 0.4]
        >>> fig = RRE.plot.range_profile(ranges, amplitudes)
    """
    fig = plt.figure(figsize=figsize)
    
    # Create stem plot
    markerline, stemlines, baseline = plt.stem(ranges, amplitudes, basefmt=' ')
    plt.setp(markerline, marker='o', markersize=8, color='red')
    plt.setp(stemlines, linewidth=2, color='blue')
    
    plt.xlabel('Range (m)', fontsize=12)
    plt.ylabel('Normalized Amplitude', fontsize=12)
    plt.title('Radar Range Profile', fontsize=14)
    plt.xlim(0, max_range)
    plt.ylim(0, max(amplitudes) * 1.2 if amplitudes else 1.0)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    if show:
        plt.show()
    
    return fig


def doppler_spectrum(
    velocities: List[float],
    amplitudes: List[float],
    max_velocity: float = 100,
    figsize: Tuple[float, float] = (10, 5),
    show: bool = True
) -> plt.Figure:
    """Plot a Doppler spectrum showing velocity distribution of targets.
    
    Creates a frequency spectrum plot showing targets at different velocities
    with their corresponding signal strengths.
    
    Args:
        velocities: List of target velocities in m/s (negative = approaching)
        amplitudes: List of corresponding signal amplitudes
        max_velocity: Maximum velocity magnitude to display (default: 100 m/s)
        figsize: Figure size as (width, height) in inches (default: (10, 5))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> velocities = [-50, -20, 10, 30]
        >>> amplitudes = [0.9, 0.7, 0.5, 0.8]
        >>> fig = RRE.plot.doppler_spectrum(velocities, amplitudes)
    """
    fig = plt.figure(figsize=figsize)
    
    # Create stem plot
    markerline, stemlines, baseline = plt.stem(velocities, amplitudes, basefmt=' ')
    plt.setp(markerline, marker='o', markersize=8, color='red')
    plt.setp(stemlines, linewidth=2, color='green')
    
    plt.xlabel('Velocity (m/s)', fontsize=12)
    plt.ylabel('Normalized Amplitude', fontsize=12)
    plt.title('Doppler Spectrum', fontsize=14)
    plt.xlim(-max_velocity, max_velocity)
    plt.ylim(0, max(amplitudes) * 1.2 if amplitudes else 1.0)
    plt.axvline(x=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add legend for velocity direction
    plt.text(0.02, 0.98, 'Negative: Approaching\nPositive: Receding',
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3),
             fontsize=9)
    
    if show:
        plt.show()
    
    return fig


def tactical_scenario(
    radar_pos: Tuple[float, float] = (10, 20),
    target_pos: Tuple[float, float] = (50, -20),
    jammer_pos: Tuple[float, float] = (70, -10),
    xlim: Tuple[float, float] = (0, 75),
    ylim: Tuple[float, float] = (-25, 25),
    figsize: Tuple[float, float] = (8, 6),
    show: bool = True
) -> plt.Figure:
    """Plot a tactical scenario showing radar, target, and support jammer positions.
    
    Creates a 2D plot showing the positions of radar, target, and support jammer
    on an X-Y coordinate system with customizable markers and colors.
    
    Args:
        radar_pos: (x, y) position of the radar in km (default: (10, 20))
        target_pos: (x, y) position of the target in km (default: (50, -20))
        jammer_pos: (x, y) position of the support jammer in km (default: (70, -10))
        xlim: X-axis limits as (min, max) in km (default: (0, 75))
        ylim: Y-axis limits as (min, max) in km (default: (-25, 25))
        figsize: Figure size as (width, height) in inches (default: (8, 6))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> # Use default positions
        >>> fig = RRE.plot.tactical_scenario()
        >>> 
        >>> # Custom positions
        >>> fig = RRE.plot.tactical_scenario(
        ...     radar_pos=(20, 15),
        ...     target_pos=(60, -15),
        ...     jammer_pos=(80, -5)
        ... )
    
    Notes:
        - Radar is displayed as a blue 'x' marker
        - Target is displayed as a green solid square marker
        - Support jammer is displayed as a red hollow circle marker
        - Grid is displayed with dashed lines for easy position reading
    """
    # Create the figure and axes
    fig = plt.figure(figsize=figsize)
    
    # Plot the individual data points
    # Radar (Blue 'x')
    plt.scatter(radar_pos[0], radar_pos[1], marker='x', color='blue', s=80, label='Radar')
    
    # Target (Green Square)
    plt.scatter(target_pos[0], target_pos[1], marker='s', color='green', s=60, label='Target')
    
    # Support jammer (Red hollow circle)
    plt.scatter(jammer_pos[0], jammer_pos[1], marker='o', facecolors='none', 
                edgecolors='red', s=80, label='Support jammer')
    
    # Customize the plot
    plt.xlabel('X (km)', fontsize=12)
    plt.ylabel('Y (km)', fontsize=12)
    
    # Set axis limits
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    
    # Set axis ticks
    x_ticks = np.arange(xlim[0], xlim[1] + 1, 10)
    y_ticks = np.arange(ylim[0], ylim[1] + 1, 5)
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)
    
    # Add the grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add the legend
    plt.legend()
    
    # Ensure the layout is clean
    plt.tight_layout()
    
    # Display the plot
    if show:
        plt.show()
    
    return fig


def generate_lfm_spectrogram(
    f_start_MHz: float = 7850, 
    f_end_MHz: float = 8150, 
    pulse_duration_us: float = 40, 
    pulse_start_times_us: List[float] = None, 
    total_duration_us: float = 180, 
    Fs_MHz: float = 500, 
    noise_factor: float = 0.05,
    figsize: Tuple[float, float] = (10, 5),
    show: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, plt.Figure]:
    """
    Generates and plots a spectrogram for Linear Frequency Modulated (LFM) radar pulses.
    
    Creates a time-frequency representation of LFM chirp pulses, showing how the
    frequency changes linearly over the duration of each pulse. This is useful for
    visualizing radar pulse compression waveforms.
    
    Parameters:
        f_start_MHz: Start frequency of the chirp in MHz (default: 7850)
        f_end_MHz: End frequency of the chirp in MHz (default: 8150)
        pulse_duration_us: Duration of a single pulse in microseconds (default: 40)
        pulse_start_times_us: List of start times for each pulse in microseconds.
            If None, defaults to [0, 120] (default: None)
        total_duration_us: Total duration of the simulation in microseconds (default: 180)
        Fs_MHz: Sampling frequency in MHz. Must be > (f_end - f_start) (default: 500)
        noise_factor: Magnitude of the complex Gaussian noise added to the signal (default: 0.05)
        figsize: Figure size as (width, height) in inches (default: (10, 5))
        show: Whether to display the plot immediately (default: True)
    
    Returns:
        Tuple containing:
            - f_plot_RF (np.ndarray): Frequency vector in MHz
            - time_sec (np.ndarray): Time vector in seconds
            - Sxx_db (np.ndarray): Spectrogram power spectral density in dB
            - fig (plt.Figure): The matplotlib Figure object
    
    Example:
        >>> import radar_range_equation as RRE
        >>> # Generate spectrogram with default parameters
        >>> f, t, Sxx, fig = RRE.plot.generate_lfm_spectrogram()
        >>> 
        >>> # Custom parameters
        >>> f, t, Sxx, fig = RRE.plot.generate_lfm_spectrogram(
        ...     f_start_MHz=7850,
        ...     f_end_MHz=8150,
        ...     pulse_duration_us=40,
        ...     pulse_start_times_us=[0, 120],
        ...     total_duration_us=180,
        ...     show=False
        ... )
    
    Notes:
        - The signal is generated at baseband (centered at 0 Hz) to avoid aliasing artifacts
        - Complex exponential representation is used to create clean chirp signals
        - The spectrogram uses overlapping windows for better time-frequency resolution
        - Dynamic range of the plot is set to 40 dB for better visualization
    """
    # Default pulse start times
    if pulse_start_times_us is None:
        pulse_start_times_us = [0, 120]
    
    # --- 1. Parameter Setup ---
    # Convert units to SI (Hz, seconds)
    Fs = Fs_MHz * 1e6
    T_total = total_duration_us * 1e-6
    T_pulse = pulse_duration_us * 1e-6
    
    f_start_RF = f_start_MHz * 1e6
    f_end_RF = f_end_MHz * 1e6
    
    # Calculate Bandwidth and Center Frequency
    bandwidth = f_end_RF - f_start_RF
    f_center_RF = (f_start_RF + f_end_RF) / 2
    
    # Check Nyquist
    if Fs <= bandwidth:
        print(f"Warning: Fs ({Fs_MHz} MHz) is too low for bandwidth ({bandwidth/1e6} MHz). Aliasing will occur.")

    # Time Vector for the total duration
    N_total = int(Fs * T_total)
    t = np.linspace(0, T_total, N_total, endpoint=False)
    
    # --- 2. Pulse Generation (Baseband) ---
    # We generate the pulse at Baseband (centered at 0 Hz) using complex exponential
    # to avoid "V-shape" artifacts at DC.
    f_start_BB = -bandwidth / 2
    f_end_BB = bandwidth / 2
    
    # Time vector for a single pulse
    N_pulse = int(Fs * T_pulse)
    t_pulse_vec = np.linspace(0, T_pulse, N_pulse, endpoint=False)
    
    # Generate Phase for Linear Chirp: phi(t) = 2*pi * (f0*t + k/2 * t^2)
    k_chirp = (f_end_BB - f_start_BB) / T_pulse
    phase = 2 * np.pi * (f_start_BB * t_pulse_vec + 0.5 * k_chirp * t_pulse_vec**2)
    signal_pulse_BB = np.exp(1j * phase)
    
    # --- 3. Signal Assembly ---
    signal_total = np.zeros(N_total, dtype=complex)
    
    for start_time_us in pulse_start_times_us:
        # Convert start time to index
        t_start = start_time_us * 1e-6
        idx_start = int(t_start * Fs)
        idx_end = idx_start + N_pulse
        
        # Bounds check to ensure we don't write past the array
        if idx_end <= N_total:
            signal_total[idx_start:idx_end] = signal_pulse_BB
        else:
            print(f"Warning: Pulse starting at {start_time_us}us exceeds total duration.")

    # Add Noise
    noise = noise_factor * (np.random.randn(N_total) + 1j * np.random.randn(N_total))
    signal_with_noise = signal_total + noise
    
    # --- 4. Spectrogram Calculation ---
    nperseg = int(0.5e-6 * Fs) # Window size (approx 0.5 us)
    noverlap = int(nperseg * 0.9)
    
    # Compute spectrogram (complex mode, full spectrum)
    f, time_sec, Sxx = spectrogram(
        signal_with_noise,
        Fs,
        nperseg=nperseg,
        noverlap=noverlap,
        scaling='spectrum',
        mode='complex',
        return_onesided=False
    )
    
    # Shift zero frequency to the center
    f = np.fft.fftshift(f)
    Sxx = np.fft.fftshift(Sxx, axes=0)
    
    # Convert to dB
    Sxx_db = 10 * np.log10(np.abs(Sxx) + 1e-12)
    
    # --- 5. Plotting ---
    fig = plt.figure(figsize=figsize)
    
    # Shift frequency axis back to RF for display
    f_plot_RF = (f + f_center_RF) / 1e6
    
    # Define view limits slightly wider than the chirp
    view_margin = 0 # MHz margin
    f_view_min = f_start_MHz - view_margin
    f_view_max = f_end_MHz + view_margin
    
    # Slice data for better contrast in plotting
    idx_min = np.searchsorted(f_plot_RF, f_view_min)
    idx_max = np.searchsorted(f_plot_RF, f_view_max)
    
    # Handle edge case where indices might be out of bounds
    idx_min = max(0, idx_min)
    idx_max = min(len(f_plot_RF), idx_max)
    
    if idx_max <= idx_min:
        # Fallback if search fails or range is invalid
        idx_min, idx_max = 0, len(f_plot_RF)

    f_plot = f_plot_RF[idx_min:idx_max]
    Sxx_plot = Sxx_db[idx_min:idx_max, :]
    
    # Calculate edges for pcolormesh
    dt = time_sec[1] - time_sec[0] if len(time_sec) > 1 else 1e-6
    time_plot_edges = np.append(time_sec, time_sec[-1] + dt) * 1e6
    
    if len(f_plot) > 1:
        df = f_plot[1] - f_plot[0]
        f_plot_edges = np.append(f_plot, f_plot[-1] + df)
    else:
        # Fallback for very small slices
        f_plot_edges = np.array([f_plot[0], f_plot[0] + 1.0])

    # Plot
    plt.pcolormesh(
        time_plot_edges, 
        f_plot_edges, 
        Sxx_plot,
        cmap='jet',
        vmin=np.max(Sxx_db) - 40, # Dynamic range of 40 dB
        vmax=np.max(Sxx_db)
    )
    
    plt.xlim(0, total_duration_us)
    plt.ylim(f_view_min, f_view_max)
    
    plt.xlabel(r'Time ($\mu$s)')
    plt.ylabel('Frequency (MHz)')
    plt.title(f'Spectrogram of LFM Radar Pulses ({f_start_MHz}-{f_end_MHz} MHz)')
    plt.colorbar(label='Power Spectral Density (dB)')
    plt.grid(False)
    
    if show:
        plt.show()
    
    return f_plot_RF, time_sec, Sxx_db, fig
