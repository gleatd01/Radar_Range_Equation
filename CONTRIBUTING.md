# Contributing to Radar Range Equation

Thank you for your interest in contributing to the Radar Range Equation project! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [How to Contribute](#how-to-contribute)
- [Code Style](#code-style)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Adding a New Language Implementation](#adding-a-new-language-implementation)

## Getting Started

Before you begin:
1. Read the [README.md](README.md) to understand the project
2. Check existing [issues](https://github.com/gleatd01/Radar_Range_Equation/issues) and [pull requests](https://github.com/gleatd01/Radar_Range_Equation/pulls)
3. For major changes, open an issue first to discuss your ideas

## Development Setup

### Python Package

1. Clone the repository:
   ```bash
   git clone https://github.com/gleatd01/Radar_Range_Equation.git
   cd Radar_Range_Equation
   ```

2. Install Python dependencies (Python 3.9+ required):
   ```bash
   pip install -e .
   ```

3. Install development dependencies:
   ```bash
   pip install -r python/requirements-dev.txt
   ```

4. Run tests to verify setup:
   ```bash
   python python/test_package.py
   python python/test_equations_solvers.py
   ```

### Flutter/Dart Package

1. Navigate to the Flutter directory:
   ```bash
   cd flutter
   ```

2. Get dependencies:
   ```bash
   dart pub get
   ```

3. Run tests:
   ```bash
   dart test
   ```

### Rust Crate

1. Navigate to the Rust directory:
   ```bash
   cd rust
   ```

2. Build the project:
   ```bash
   cargo build
   ```

3. Run tests:
   ```bash
   cargo test
   ```

## How to Contribute

### Reporting Bugs

When reporting bugs, please include:
- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior vs actual behavior
- Code samples or test cases
- Your environment (OS, language version, package version)

### Suggesting Enhancements

For feature requests:
- Use a clear, descriptive title
- Provide a detailed description of the proposed feature
- Explain why this feature would be useful
- Include code examples if applicable

### Code Contributions

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add or update tests as needed
5. Update documentation as needed
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to your branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

## Code Style

### Python
- Follow [PEP 8](https://pep8.org/) style guide
- Use type hints where applicable
- Write docstrings for functions and classes
- Keep functions focused and concise
- Use meaningful variable names

### Dart/Flutter
- Follow [Effective Dart](https://dart.dev/guides/language/effective-dart) guidelines
- Use camelCase for variable and function names
- Document public APIs

### Rust
- Follow [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Run `cargo fmt` before committing
- Run `cargo clippy` and address warnings
- Document public APIs with doc comments

## Testing

All contributions should include appropriate tests:

### Python
- Add unit tests in the `python/` directory
- Test files should start with `test_`
- Run all tests before submitting:
  ```bash
  python python/test_package.py
  python python/test_equations_solvers.py
  python python/test_plotting.py
  python python/test_analysis.py
  ```

### Flutter/Dart
- Add tests in the `flutter/test/` directory
- Run tests with: `dart test`

### Rust
- Add tests in the `rust/src/` or `rust/tests/` directory
- Run tests with: `cargo test`

## Submitting Changes

### Pull Request Guidelines

1. **Title**: Use a clear, descriptive title
2. **Description**: Include:
   - What changes were made
   - Why the changes were made
   - Any breaking changes
   - Related issue numbers (if applicable)
3. **Tests**: Ensure all tests pass
4. **Documentation**: Update relevant documentation
5. **Commits**: Keep commits focused and write clear commit messages

### Pull Request Process

1. Ensure your code follows the style guidelines
2. Update documentation as needed
3. Add tests for new functionality
4. Ensure all tests pass
5. Update the README.md if needed
6. Submit the pull request
7. Respond to review feedback

## Adding a New Language Implementation

We welcome implementations in additional programming languages! To add a new language:

1. **Create a new directory** at the root level (e.g., `javascript/`, `java/`, etc.)

2. **Implement core functionality**:
   - Calculate maximum detection range
   - Calculate received power
   - Calculate wavelength from frequency
   - Calculate Doppler frequency
   - Unit conversion utilities

3. **Follow the existing structure**:
   ```
   language_name/
   ├── src/          # Source code
   ├── test/         # Tests
   ├── examples/     # Usage examples
   └── README.md     # Language-specific documentation
   ```

4. **Include comprehensive tests** matching the functionality in other implementations

5. **Add documentation**:
   - Create a detailed README.md in your language directory
   - Include installation instructions
   - Provide usage examples
   - Document the API

6. **Update the main README**:
   - Add your language to the package structure
   - Add a quick start example
   - Link to your language-specific README

7. **Maintain consistency**:
   - Use similar function/method names (adjusted for language conventions)
   - Implement the same core calculations
   - Return the same units

## Questions?

If you have questions or need help, please:
- Open an issue on GitHub
- Check existing documentation
- Review examples in the repository

Thank you for contributing to Radar Range Equation!
