#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build complete-only render tasks for substitutional/interstitial species. "
            "Tasks are consumed by run_render_task_worker.py"
        )
    )
    p.add_argument("--defects-root", default="__variations__/defects")
    p.add_argument("--tetra-defects-root", default="__variations__/defects_tetra")
    p.add_argument("--queue-root", default="plots/render_queue")
    p.add_argument("--modes", nargs="+", default=["substitutional", "interstitial", "vacancy"])  # one or more
    p.add_argument("--refresh-substitutional-maps", action="store_true")
    p.add_argument("--refresh-vacancy-maps", action="store_true")
    p.add_argument("--overwrite-tasks", action="store_true")
    return p.parse_args()


def is_case_done(case: Path) -> bool:
    return (
        (case / "CONTCAR_MACE").exists()
        and (case / "initial-energies.txt").exists()
        and (case / "relaxed-energies.txt").exists()
    )


def species_completion(defects_root: Path, tetra_defects_root: Path, mode: str) -> list[tuple[str, int, int]]:
    if mode == "vacancy":
        root = defects_root / "vacancy" / "1-site"
        if not root.exists():
            return []
        cases = sorted([c for c in root.iterdir() if c.is_dir() and c.name.startswith("v")])
        total = len(cases)
        done = sum(1 for c in cases if is_case_done(c))
        return [("Si", done, total)]

    root_base = tetra_defects_root if mode == "interstitial" else defects_root
    root = root_base / mode / "1-site"
    rows: list[tuple[str, int, int]] = []
    if not root.exists():
        return rows

    for sp in sorted([p for p in root.iterdir() if p.is_dir()]):
        cases = sorted([c for c in sp.iterdir() if c.is_dir() and c.name.startswith("v")])
        total = len(cases)
        done = sum(1 for c in cases if is_case_done(c))
        rows.append((sp.name, done, total))
    return rows


def build_task(mode: str, species: str) -> dict:
    if mode == "substitutional":
        cmd = [
            sys.executable,
            "render_site_maps_ase.py",
            "--root",
            "plots/substitutional_site_maps",
            "--species",
            species,
            "--output-name",
            "concat.png",
            "--backup-mode",
            "move",
        ]
        output_path = f"plots/substitutional_site_maps/{species}/concat.png"
    elif mode == "vacancy":
        cmd = [
            sys.executable,
            "render_site_maps_ase.py",
            "--root",
            "plots/vacancy_site_maps",
            "--species",
            species,
            "--output-name",
            "concat.png",
            "--backup-mode",
            "move",
        ]
        output_path = f"plots/vacancy_site_maps/{species}/concat.png"
    else:
        cmd = [
            sys.executable,
            "render_interstitial_maps_ase.py",
            "--outdir",
            "plots/interstitial_tetra_site_maps",
            "--species",
            species,
            "--backup-mode",
            "move",
            "--min-done",
            "1",
        ]
        output_path = f"plots/interstitial_tetra_site_maps/{species}/concat.png"

    return {
        "mode": mode,
        "species": species,
        "command": cmd,
        "output_path": output_path,
    }


def maybe_refresh_substitutional_maps(enabled: bool) -> None:
    if not enabled:
        return
    subprocess.run([sys.executable, "build_substitutional_site_maps.py"], check=True)


def maybe_refresh_vacancy_maps(enabled: bool) -> None:
    if not enabled:
        return
    subprocess.run([sys.executable, "build_vacancy_site_map.py"], check=True)


def main() -> int:
    args = parse_args()
    defects_root = Path(args.defects_root).resolve()
    tetra_defects_root = Path(args.tetra_defects_root).resolve()
    queue_root = Path(args.queue_root).resolve()
    tasks_dir = queue_root / "tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)

    maybe_refresh_substitutional_maps(args.refresh_substitutional_maps)
    maybe_refresh_vacancy_maps(args.refresh_vacancy_maps)

    created = 0
    skipped = 0
    summary = []

    for mode in args.modes:
        rows = species_completion(defects_root, tetra_defects_root, mode)
        for species, done, total in rows:
            summary.append({"mode": mode, "species": species, "done": done, "total": total})
            if total == 0 or done != total:
                continue

            task = build_task(mode, species)
            task_name = f"{mode}__{species}.json"
            task_path = tasks_dir / task_name
            done_path = queue_root / "done" / task_name
            if done_path.exists() and not args.overwrite_tasks:
                skipped += 1
                continue
            if task_path.exists() and not args.overwrite_tasks:
                skipped += 1
                continue
            task_path.write_text(json.dumps(task, indent=2), encoding="ascii")
            created += 1

    (queue_root / "summary.json").write_text(json.dumps(summary, indent=2), encoding="ascii")
    print(f"tasks_created={created}")
    print(f"tasks_skipped={skipped}")
    print(f"tasks_dir={tasks_dir}")
    print(f"summary={queue_root / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
