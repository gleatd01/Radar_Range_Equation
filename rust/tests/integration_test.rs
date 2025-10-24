use radar_range_equation::{calculate_max_range, calculate_received_power};

#[test]
fn test_integration_max_range() {
    let max_range = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
    assert!(max_range > 0.0);
    // Reasonable range for these parameters (should be around 8 km)
    assert!(max_range > 5000.0);
    assert!(max_range < 15000.0);
}

#[test]
fn test_integration_received_power() {
    let power = calculate_received_power(1000.0, 1000.0, 0.03, 1.0, 10000.0).unwrap();
    assert!(power > 0.0);
    // Power should be very small at this range
    assert!(power < 1e-10);
}

#[test]
fn test_integration_roundtrip() {
    let transmit_power = 1000.0;
    let antenna_gain = 1000.0;
    let wavelength = 0.03;
    let radar_cross_section = 1.0;
    let min_detectable_signal = 1e-13;

    // Calculate max range
    let max_range = calculate_max_range(
        transmit_power,
        antenna_gain,
        wavelength,
        radar_cross_section,
        min_detectable_signal,
    )
    .unwrap();

    // Calculate received power at max range
    let received_power = calculate_received_power(
        transmit_power,
        antenna_gain,
        wavelength,
        radar_cross_section,
        max_range,
    )
    .unwrap();

    // Should match the minimum detectable signal (within 1%)
    let relative_error = (received_power - min_detectable_signal).abs() / min_detectable_signal;
    assert!(relative_error < 0.01);
}
