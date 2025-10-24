use radar_range_equation::{calculate_max_range, calculate_received_power};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Radar Range Equation Example");
    println!("{}", "=".repeat(40));

    // System parameters
    let transmit_power = 1000.0; // 1 kW
    let antenna_gain = 1000.0; // 30 dBi (linear scale)
    let wavelength = 0.03; // 10 GHz (c/f = 3e8/10e9)
    let radar_cross_section = 1.0; // 1 m²
    let min_detectable_signal = 1e-13; // -100 dBm

    println!("Transmit Power: {} W", transmit_power);
    println!("Antenna Gain: {} (linear)", antenna_gain);
    println!("Wavelength: {} m", wavelength);
    println!("Target RCS: {} m²", radar_cross_section);
    println!("Min Detectable Signal: {} W", min_detectable_signal);
    println!();

    // Calculate maximum range
    let max_range = calculate_max_range(
        transmit_power,
        antenna_gain,
        wavelength,
        radar_cross_section,
        min_detectable_signal,
    )?;

    println!(
        "Maximum Detection Range: {:.2} m ({:.2} km)",
        max_range,
        max_range / 1000.0
    );
    println!();

    // Calculate received power at various ranges
    println!("Received Power at Different Ranges:");
    println!("{}", "-".repeat(40));
    for range_km in [5, 10, 15, 20] {
        let range_m = (range_km * 1000) as f64;
        let power = calculate_received_power(
            transmit_power,
            antenna_gain,
            wavelength,
            radar_cross_section,
            range_m,
        )?;
        println!("  {} km: {:.2e} W", range_km, power);
    }

    Ok(())
}
