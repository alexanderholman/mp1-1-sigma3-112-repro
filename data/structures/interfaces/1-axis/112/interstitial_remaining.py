#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Report tetra n=1 interstitial completion/remaining counts.")
    p.add_argument("--root", default="__variations__/defects_tetra/interstitial/1-site")
    return p.parse_args()


def is_done(case: Path) -> bool:
    return (
        (case / "CONTCAR_MACE").exists()
        and (case / "initial-energies.txt").exists()
        and (case / "relaxed-energies.txt").exists()
    )


def main() -> int:
    root = Path(parse_args().root)
    total = 0
    done = 0
    for sp in sorted([p for p in root.iterdir() if p.is_dir()]):
        st = 0
        sd = 0
        for c in sorted([d for d in sp.iterdir() if d.is_dir() and d.name.startswith("v")]):
            total += 1
            st += 1
            if is_done(c):
                done += 1
                sd += 1
        print(f"{sp.name} {sd}/{st}")

    print(f"total {done}/{total} remaining={total-done}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
