#!/usr/bin/env python3
"""
Write calendar-dev version into python/src/radar_range_equation/__version__.py.

Default timezone: America/New_York (handles EDT/EST via zoneinfo).
Format: YYYY.MM.DD.HH.dev
Example: 2025.10.27.16.dev
"""
from datetime import datetime
from pathlib import Path
import argparse
import sys

try:
    from zoneinfo import ZoneInfo
except Exception:  # pragma: no cover - Python < 3.9
    ZoneInfo = None

# Path relative to repository root
VERSION_FILE = Path("python/src/radar_range_equation/__version__.py")


def get_caldev_version(now: datetime) -> str:
    # YYYY.MM.DD.HH.dev using provided datetime (tz-aware or naive)
    return now.strftime("%Y.%m.%d.%H") + ".dev"


def write_version_file(version: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    content = f'''# This file is managed by scripts/write_version.py.
# Do not edit manually if you use this script to set versions.
__all__ = ["__version__"]

__version__ = "{version}"
'''
    path.write_text(content, encoding="utf-8")
    print(f"Wrote version {version!r} to {path}")


def parse_args():
    p = argparse.ArgumentParser(description="Write calendar-dev version to package __version__.py")
    p.add_argument("--tz", default="America/New_York", help="Timezone name for calendar version (IANA), default: America/New_York")
    p.add_argument("--utc", action="store_true", help="Force UTC time instead of local timezone")
    p.add_argument("--out", default=str(VERSION_FILE), help="Path to write version file")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.utc:
        now = datetime.utcnow()
    else:
        if ZoneInfo is None:
            print("zoneinfo not available on this Python; falling back to UTC. (Use Python >= 3.9 for zoneinfo)", file=sys.stderr)
            now = datetime.utcnow()
        else:
            try:
                tz = ZoneInfo(args.tz)
            except Exception as exc:
                print(f"Failed to load timezone '{args.tz}': {exc}. Falling back to UTC.", file=sys.stderr)
                now = datetime.utcnow()
            else:
                now = datetime.now(tz)

    version = get_caldev_version(now)
    write_version_file(version, Path(args.out))


if __name__ == "__main__":
    main()
