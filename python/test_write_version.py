#!/usr/bin/env python3
"""
Test script for write_version.py functionality.
"""
import sys
import tempfile
from pathlib import Path
from datetime import datetime

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from write_version import get_caldev_version, write_version_file


def test_get_caldev_version():
    """Test that calendar version format is correct."""
    # Test with a known datetime
    test_dt = datetime(2025, 10, 27, 16, 30, 45)
    version = get_caldev_version(test_dt)
    
    expected = "2025.10.27.16.dev"
    assert version == expected, f"Expected {expected}, got {version}"
    print(f"✓ get_caldev_version() format test passed: {version}")


def test_write_version_file():
    """Test that version file is written correctly."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_path = Path(tmpdir) / "test_package" / "__version__.py"
        test_version = "2025.10.27.16.dev"
        
        write_version_file(test_version, test_path)
        
        # Verify file exists
        assert test_path.exists(), f"File was not created at {test_path}"
        
        # Verify content
        content = test_path.read_text(encoding="utf-8")
        assert test_version in content, f"Version {test_version} not found in content"
        assert '__version__ = "2025.10.27.16.dev"' in content, "Incorrect version format in file"
        assert '__all__ = ["__version__"]' in content, "__all__ declaration missing"
        
        print(f"✓ write_version_file() test passed")
        print(f"  Created file at: {test_path}")
        print(f"  Content preview: {content[:100]}...")


def test_version_format_components():
    """Test that version format has all required components."""
    test_dt = datetime(2025, 1, 5, 9, 0, 0)  # Test with single-digit month/day/hour
    version = get_caldev_version(test_dt)
    
    parts = version.replace(".dev", "").split(".")
    assert len(parts) == 4, f"Expected 4 parts in version, got {len(parts)}: {parts}"
    
    year, month, day, hour = parts
    assert year == "2025", f"Year should be 2025, got {year}"
    assert month == "01", f"Month should be zero-padded 01, got {month}"
    assert day == "05", f"Day should be zero-padded 05, got {day}"
    assert hour == "09", f"Hour should be zero-padded 09, got {hour}"
    assert version.endswith(".dev"), f"Version should end with .dev, got {version}"
    
    print(f"✓ version format components test passed: {version}")


def main():
    """Run all tests."""
    try:
        print("Testing write_version.py functionality...\n")
        
        test_get_caldev_version()
        test_write_version_file()
        test_version_format_components()
        
        print("\n✓ All write_version.py tests passed!")
        return 0
        
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
