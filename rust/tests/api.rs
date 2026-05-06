use radar_range_equation::*;

#[test]
fn test_calculate_wavelength() {
    let wavelength = calculate_wavelength(10e9).unwrap();
    assert!((wavelength - 0.0299792458).abs() < 1e-6);
}

#[test]
fn test_calculate_frequency() {
    let frequency = calculate_frequency(0.03).unwrap();
    assert!((frequency - 9.99308193e9).abs() < 1e3);
}

#[test]
fn test_calculate_max_range() {
    let max_range = calculate_max_range(1000.0, 1000.0, 0.03, 1.0, 1e-13).unwrap();
    assert!(max_range > 0.0);
    assert!(max_range > 8000.0 && max_range < 9000.0);
}

#[test]
fn test_calculate_received_power() {
    let power = calculate_received_power(1000.0, 1000.0, 1000.0, 0.03, 1.0, 10000.0).unwrap();
    assert!(power > 0.0);
}

#[test]
fn test_calculate_doppler_frequency() {
    let doppler = calculate_doppler_frequency(-100.0, 0.03).unwrap();
    assert!((doppler - 6666.67).abs() < 0.01);
}

#[test]
fn test_calculate_velocity_from_doppler() {
    let velocity = calculate_velocity_from_doppler(6666.67, 0.03).unwrap();
    assert!((velocity - (-100.0)).abs() < 0.01);
}

#[test]
fn test_sphere_rcs() {
    let rcs = calculate_sphere_rcs(1.0).unwrap();
    assert!((rcs - constants::PI).abs() < 1e-6);
}

#[test]
fn test_sphere_radius() {
    let radius = calculate_sphere_radius(constants::PI).unwrap();
    assert!((radius - 1.0).abs() < 1e-6);
}

#[test]
fn test_linear_to_db() {
    assert!((convert::linear_to_db(10.0) - 10.0).abs() < 1e-6);
    assert!((convert::linear_to_db(100.0) - 20.0).abs() < 1e-6);
    assert!((convert::linear_to_db(1.0) - 0.0).abs() < 1e-6);
}

#[test]
fn test_db_to_linear() {
    assert!((convert::db_to_linear(10.0) - 10.0).abs() < 1e-6);
    assert!((convert::db_to_linear(20.0) - 100.0).abs() < 1e-6);
    assert!((convert::db_to_linear(0.0) - 1.0).abs() < 1e-6);
}

#[test]
fn test_angle_conversions() {
    assert!((convert::rad_to_deg(constants::PI) - 180.0).abs() < 1e-6);
    assert!((convert::deg_to_rad(180.0) - constants::PI).abs() < 1e-6);
}

#[test]
fn test_distance_conversions() {
    assert!((convert::meters_to_nautical_miles(1852.0) - 1.0).abs() < 1e-6);
    assert!((convert::feet_to_meters(1.0) - 0.3048).abs() < 1e-6);
}

#[test]
fn test_frequency_conversions() {
    assert!((convert::hz_to_mhz(1e6) - 1.0).abs() < 1e-6);
    assert!((convert::hz_to_ghz(1e9) - 1.0).abs() < 1e-6);
}

#[test]
fn test_invalid_parameters() {
    assert!(calculate_wavelength(-1.0).is_err());
    assert!(calculate_max_range(-1.0, 1000.0, 0.03, 1.0, 1e-13).is_err());
    assert!(calculate_sphere_rcs(-1.0).is_err());
}
