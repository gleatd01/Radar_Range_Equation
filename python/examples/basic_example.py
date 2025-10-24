"""Example usage of the radar range equation package."""

from radar_range_equation import calculate_max_range, calculate_received_power


def main():
    """Demonstrate basic radar calculations."""
    
    # System parameters
    transmit_power = 1000  # 1 kW
    antenna_gain = 1000    # 30 dBi (linear scale)
    wavelength = 0.03      # 10 GHz (c/f = 3e8/10e9)
    radar_cross_section = 1.0  # 1 m²
    min_detectable_signal = 1e-13  # -100 dBm
    
    print("Radar Range Equation Example")
    print("=" * 40)
    print(f"Transmit Power: {transmit_power} W")
    print(f"Antenna Gain: {antenna_gain} (linear)")
    print(f"Wavelength: {wavelength} m")
    print(f"Target RCS: {radar_cross_section} m²")
    print(f"Min Detectable Signal: {min_detectable_signal} W")
    print()
    
    # Calculate maximum range
    max_range = calculate_max_range(
        transmit_power=transmit_power,
        antenna_gain=antenna_gain,
        wavelength=wavelength,
        radar_cross_section=radar_cross_section,
        min_detectable_signal=min_detectable_signal
    )
    
    print(f"Maximum Detection Range: {max_range:.2f} m ({max_range/1000:.2f} km)")
    print()
    
    # Calculate received power at various ranges
    print("Received Power at Different Ranges:")
    print("-" * 40)
    for range_km in [5, 10, 15, 20]:
        range_m = range_km * 1000
        power = calculate_received_power(
            transmit_power=transmit_power,
            antenna_gain=antenna_gain,
            wavelength=wavelength,
            radar_cross_section=radar_cross_section,
            range_m=range_m
        )
        print(f"  {range_km} km: {power:.2e} W")


if __name__ == "__main__":
    main()
