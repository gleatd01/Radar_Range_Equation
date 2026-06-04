#!/bin/bash
# Development environment setup script for Radar Range Equation

set -e  # Exit on error

echo "========================================"
echo "Radar Range Equation - Development Setup"
echo "========================================"
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_status() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_error() {
    echo -e "${RED}[✗]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    print_error "Please run this script from the root of the repository (the directory containing pyproject.toml)"
    exit 1
fi

# Function to setup Python environment
setup_python() {
    echo ""
    echo "Setting up Python environment..."
    echo "--------------------------------"
    
    # Check Python version
    if ! command -v python3 &> /dev/null; then
        print_error "Python 3 is not installed. Please install Python 3.9 or higher."
        return 1
    fi
    
    python_version=$(python3 --version | cut -d' ' -f2 | cut -d'.' -f1,2)
    print_status "Found Python $python_version"
    
    # Install package in editable mode
    echo "Installing Python package in editable mode..."
    pip install -e . || {
        print_error "Failed to install Python package"
        return 1
    }
    print_status "Python package installed"
    
    # Install development dependencies
    if [ -f "python/requirements-dev.txt" ]; then
        echo "Installing development dependencies..."
        pip install -r python/requirements-dev.txt || {
            print_warning "Some development dependencies failed to install"
        }
        print_status "Development dependencies installed"
    fi
    
    # Run tests
    echo "Running Python tests..."
    python python/test_package.py && print_status "Package tests passed" || print_warning "Some tests failed"
    
    print_status "Python setup complete!"
    return 0
}

# Function to setup Flutter/Dart environment
setup_flutter() {
    echo ""
    echo "Setting up Flutter/Dart environment..."
    echo "--------------------------------------"
    
    if ! command -v dart &> /dev/null; then
        print_warning "Dart is not installed. Skipping Flutter setup."
        echo "To install Dart, visit: https://dart.dev/get-dart"
        return 1
    fi
    
    dart_version=$(dart --version 2>&1 | head -n 1)
    print_status "Found $dart_version"
    
    cd flutter
    echo "Getting Flutter/Dart dependencies..."
    dart pub get || {
        print_error "Failed to get Dart dependencies"
        cd ..
        return 1
    }
    print_status "Flutter/Dart dependencies installed"
    
    echo "Running Dart tests..."
    dart test && print_status "Dart tests passed" || print_warning "Some tests failed"
    cd ..
    
    print_status "Flutter/Dart setup complete!"
    return 0
}

# Function to setup Rust environment
setup_rust() {
    echo ""
    echo "Setting up Rust environment..."
    echo "------------------------------"
    
    if ! command -v cargo &> /dev/null; then
        print_warning "Rust is not installed. Skipping Rust setup."
        echo "To install Rust, visit: https://rustup.rs/"
        return 1
    fi
    
    rust_version=$(rustc --version)
    print_status "Found $rust_version"
    
    cd rust
    echo "Building Rust crate..."
    cargo build || {
        print_error "Failed to build Rust crate"
        cd ..
        return 1
    }
    print_status "Rust crate built"
    
    echo "Running Rust tests..."
    cargo test && print_status "Rust tests passed" || print_warning "Some tests failed"
    cd ..
    
    print_status "Rust setup complete!"
    return 0
}

# Main setup logic
main() {
    echo "This script will set up your development environment for:"
    echo "  - Python"
    echo "  - Flutter/Dart (if installed)"
    echo "  - Rust (if installed)"
    echo ""
    
    read -p "Continue? (y/N) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Setup cancelled."
        exit 0
    fi
    
    # Setup Python (required)
    setup_python || {
        print_error "Python setup failed. This is required."
        exit 1
    }
    
    # Setup Flutter/Dart (optional)
    setup_flutter || print_warning "Flutter/Dart setup skipped or failed"
    
    # Setup Rust (optional)
    setup_rust || print_warning "Rust setup skipped or failed"
    
    echo ""
    echo "========================================"
    echo "Development environment setup complete!"
    echo "========================================"
    echo ""
    echo "Next steps:"
    echo "  1. Read CONTRIBUTING.md for contribution guidelines"
    echo "  2. Check out QUICKSTART.md for usage examples"
    echo "  3. Explore the example files in python/, flutter/, and rust/"
    echo "  4. Run tests before making changes:"
    echo "     - Python: python python/test_package.py"
    echo "     - Dart: cd flutter && dart test"
    echo "     - Rust: cd rust && cargo test"
    echo ""
    print_status "Happy coding!"
}

# Run main function
main
