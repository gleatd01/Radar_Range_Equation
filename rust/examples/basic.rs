//! Example demonstrating the usage of the radar_range_equation library

use radar_range_equation::*;

fn main() {
    println!("Radar Range Equation - Rust Example\n");

    // Example 1: Calculate wavelength from frequency
    println!("Example 1: Calculate wavelength");
    let frequency = 10e9; // 10 GHz
    let wavelength = calculate_wavelength(frequency).unwrap();
    println!("Frequency: {:.2} GHz", convert::hz_to_ghz(frequency));
    println!("Wavelength: {:.4} m ({:.2} cm)\n", wavelength, wavelength * 100.0);

    // Example 2: Calculate maximum range
    println!("Example 2: Calculate maximum radar range");
    let transmit_power = 1000.0; // 1 kW
    let antenna_gain = 1000.0; // Linear gain
    let radar_cross_section = 1.0; // 1 m² RCS
    let min_detectable_signal = 1e-13; // 0.1 pW

    let max_range = calculate_max_range(
        transmit_power,
        antenna_gain,
        wavelength,
        radar_cross_section,
        min_detectable_signal,
    )
    .unwrap();

    println!(
        "Transmit Power: {:.2} kW",
        convert::watts_to_kilowatts(transmit_power)
    );
    println!("Antenna Gain: {} (linear)", antenna_gain);
    println!("Wavelength: {:.4} m", wavelength);
    println!("Radar Cross Section: {} m²", radar_cross_section);
    println!("Min Detectable Signal: {} W", min_detectable_signal);
    println!("Maximum Range: {:.2} meters", max_range);
    println!(
        "Maximum Range: {:.2} km\n",
        convert::meters_to_kilometers(max_range)
    );

    // Example 3: Calculate Doppler frequency
    println!("Example 3: Calculate Doppler frequency shift");
    let velocity = -100.0; // 100 m/s closing velocity
    let doppler_freq = calculate_doppler_frequency(velocity, wavelength).unwrap();
    println!("Target velocity: {} m/s (closing)", -velocity);
    println!("Wavelength: {:.4} m", wavelength);
    println!(
        "Doppler frequency shift: {:.3} MHz\n",
        convert::hz_to_mhz(doppler_freq)
    );

    // Example 4: Calculate antenna properties
    println!("Example 4: Calculate antenna effective aperture");
    let efficiency = 0.6;
    let diameter = 10.0; // 10 meters
    let effective_aperture = calculate_effective_aperture_circ(efficiency, diameter).unwrap();
    println!("Antenna diameter: {} m", diameter);
    println!("Antenna efficiency: {}", efficiency);
    println!("Effective aperture: {:.2} m²\n", effective_aperture);

    // Example 5: Unit conversions
    println!("Example 5: Unit conversions");
    let range_meters = 50000.0;
    println!("Range: {} meters", range_meters);
    println!(
        "  = {:.2} km",
        convert::meters_to_kilometers(range_meters)
    );
    println!(
        "  = {:.2} nautical miles",
        convert::meters_to_nautical_miles(range_meters)
    );
    println!("  = {:.2} miles\n", convert::meters_to_miles(range_meters));

    let power_db = 30.0; // 30 dB
    let power_linear = convert::db_to_linear(power_db);
    println!(
        "Power: {} dB = {:.2} (linear)",
        power_db, power_linear
    );

    // Example 6: Sphere radar cross-section
    println!("\nExample 6: Sphere radar cross-section");
    let sphere_radius = 0.5; // 0.5 meter radius
    let sphere_rcs = calculate_sphere_rcs(sphere_radius).unwrap();
    println!("Sphere radius: {} m", sphere_radius);
    println!("Sphere RCS: {:.4} m²", sphere_rcs);
    
    // Verify by calculating radius back
    let calculated_radius = calculate_sphere_radius(sphere_rcs).unwrap();
    println!("Verified radius: {:.4} m", calculated_radius);
}
