# Rust Package

This directory contains the Rust implementation of the Radar Range Equation tools.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
radar-range-equation = { git = "https://github.com/gleatd01/Radar_Range_Equation.git", subdir = "rust" }
```

Or for local development:

```toml
[dependencies]
radar-range-equation = { path = "../rust" }
```

## Usage

```rust
use radar_range_equation::{calculate_max_range, calculate_received_power};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Calculate maximum detection range
    let max_range = calculate_max_range(
        1000.0,    // 1 kW
        1000.0,    // 30 dBi (linear scale)
        0.03,      // 10 GHz
        1.0,       // 1 m²
        1e-13,     // -100 dBm
    )?;

    println!("Maximum range: {:.2} meters", max_range);

    // Calculate received power at a specific range
    let received_power = calculate_received_power(
        1000.0,
        1000.0,
        0.03,
        1.0,
        10000.0,  // 10 km
    )?;

    println!("Received power: {:.2e} Watts", received_power);

    Ok(())
}
```

## Building

```bash
cd rust
cargo build
```

## Testing

```bash
cd rust
cargo test
```

## Running Examples

```bash
cd rust
cargo run --example basic_example
```

## Documentation

Generate and view the documentation:

```bash
cd rust
cargo doc --open
```

## Linting

```bash
cd rust
cargo clippy
```

## Formatting

```bash
cd rust
cargo fmt
```

## Project Structure

```
rust/
├── src/
│   └── lib.rs
├── tests/
│   └── integration_test.rs
├── examples/
│   └── basic_example.rs
├── Cargo.toml
└── README.md
```
