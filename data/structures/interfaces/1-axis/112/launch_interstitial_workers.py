#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    target = Path(__file__).with_name("launch_tetra_workers.py")
    print("Legacy wrapper: use launch_tetra_workers.py for canonical tetra interstitial jobs.")
    argv = [sys.executable, str(target), *sys.argv[1:]]
    raise SystemExit(__import__("subprocess").call(argv))


if __name__ == "__main__":
    raise SystemExit(main())
