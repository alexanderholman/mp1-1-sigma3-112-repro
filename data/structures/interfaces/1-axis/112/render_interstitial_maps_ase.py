#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.atoms import Atoms as AtomsType
from ase.data import atomic_numbers, covalent_radii
from ase.io import read
from ase.visualize.plot import plot_atoms
from PIL import Image, ImageDraw, ImageFont


VIEW_SPECS = {
    "top": "0x,0y,0z",
    "front": "90x,0y,0z",
    "left": "0x,90y,0z",
    "perspective": "-25x,35y,10z",
}


RAINBOW_ANCHORS = [
    (0.00, (215, 48, 39)),   # red   (negative)
    (0.16, (244, 109, 67)),  # orange
    (0.32, (254, 224, 139)), # yellow
    (0.50, (46, 181, 76)),   # green (zero)
    (0.66, (102, 194, 165)), # cyan-teal
    (0.80, (50, 136, 189)),  # blue
    (0.90, (94, 79, 162)),   # indigo
    (1.00, (123, 50, 148)),  # violet (positive)
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render interstitial maps with host atoms in gray and interstitial points in color. "
            "Color meaning: red=negative, green=zero, violet=positive."
        )
    )
    parser.add_argument("--manifest", default="__variations__/defects_tetra/manifest_interstitial_tetra.json")
    parser.add_argument("--base-poscar", default="CALC/VASP/unrelaxed/relax/POSCAR")
    parser.add_argument("--site", default="1-site")
    parser.add_argument("--outdir", default="plots/interstitial_tetra_site_maps")
    parser.add_argument("--species", nargs="*", default=["all"])
    parser.add_argument("--output-name", default="concat.png")
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--legend-width-ratio", type=float, default=0.38)
    parser.add_argument("--host-radius", type=float, default=0.55)
    parser.add_argument("--host-alpha", type=float, default=0.25)
    parser.add_argument("--dopant-radius-scale", type=float, default=0.85)
    parser.add_argument("--dopant-alpha", type=float, default=0.95)
    parser.add_argument("--backup-mode", choices=["move", "copy", "none"], default="move")
    parser.add_argument("--min-done", type=int, default=1)
    return parser.parse_args()


def ensure_atoms(obj) -> AtomsType:
    if isinstance(obj, AtomsType):
        return obj
    if isinstance(obj, list) and obj and isinstance(obj[0], AtomsType):
        return obj[0]
    raise ValueError("Could not load Atoms")


def rainbow_color_u8(t: float) -> tuple[int, int, int]:
    t = max(0.0, min(1.0, t))
    for i in range(len(RAINBOW_ANCHORS) - 1):
        t0, c0 = RAINBOW_ANCHORS[i]
        t1, c1 = RAINBOW_ANCHORS[i + 1]
        if t0 <= t <= t1:
            u = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
            return (
                int(c0[0] * (1.0 - u) + c1[0] * u),
                int(c0[1] * (1.0 - u) + c1[1] * u),
                int(c0[2] * (1.0 - u) + c1[2] * u),
            )
    return RAINBOW_ANCHORS[-1][1]


def value_to_rgba(v: float, vabs: float, alpha: float) -> tuple[float, float, float, float]:
    if vabs <= 1e-12:
        c = rainbow_color_u8(0.5)
    else:
        t = (v + vabs) / (2.0 * vabs)
        c = rainbow_color_u8(t)
    return (c[0] / 255.0, c[1] / 255.0, c[2] / 255.0, alpha)


def make_backup(species_dir: Path, targets: list[Path], mode: str) -> None:
    if mode == "none":
        return
    existing = [p for p in targets if p.exists()]
    if not existing:
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    bdir = species_dir / f"_backup_{stamp}"
    bdir.mkdir(parents=True, exist_ok=True)
    for src in existing:
        dst = bdir / src.name
        if mode == "move":
            shutil.move(str(src), str(dst))
        else:
            shutil.copy2(str(src), str(dst))


def load_manifest(path: Path) -> list[dict]:
    return json.loads(path.read_text(encoding="ascii"))


def read_delta(case_dir: Path) -> float | None:
    e0 = case_dir / "initial-energies.txt"
    e1 = case_dir / "relaxed-energies.txt"
    if not (e0.exists() and e1.exists()):
        return None
    try:
        return float(e1.read_text(encoding="utf-8").strip()) - float(e0.read_text(encoding="utf-8").strip())
    except Exception:
        return None


def as_float(value: object, default: float = 0.0) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except Exception:
            return default
    return default


def as_frac3(value: object) -> tuple[float, float, float]:
    if isinstance(value, (list, tuple)) and len(value) >= 3:
        return (as_float(value[0]), as_float(value[1]), as_float(value[2]))
    return (0.0, 0.0, 0.0)


def collect_interstitial_data(manifest_rows: list[dict], site: str) -> dict[str, dict[int, dict[str, object]]]:
    data: dict[str, dict[int, dict[str, object]]] = defaultdict(dict)
    target_n = int(site.split("-site", 1)[0])
    for row in manifest_rows:
        if row.get("defect_type") != "interstitial":
            continue
        if int(row.get("n_sites", 0)) != target_n:
            continue
        idxs = row.get("selected_indices", [])
        fracs = row.get("selected_frac", [])
        if len(idxs) != 1 or len(fracs) != 1:
            continue
        species = row.get("element")
        idx = int(idxs[0])
        frac = tuple(float(x) for x in fracs[0])
        case_dir = Path(row["output_dir"]).resolve()
        delta = read_delta(case_dir)
        done = delta is not None
        data[str(species)][idx] = {
            "frac": frac,
            "delta": float(delta) if done else 0.0,
            "done": done,
            "case_dir": str(case_dir),
        }
    return data


def render_view_png(atoms: AtomsType, colors: list[tuple[float, float, float, float]], radii: list[float], rotation: str, out_path: Path, dpi: int) -> None:
    fig = plt.figure(figsize=(4, 4), dpi=dpi)
    ax = fig.add_subplot(111)
    ax.set_facecolor((0.95, 0.95, 0.95))
    plot_atoms(
        atoms,
        ax,
        rotation=rotation,
        radii=radii,
        colors=colors,
        show_unit_cell=2,
    )
    ax.set_axis_off()
    fig.tight_layout(pad=0)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def draw_concat_with_legend(species: str, image_paths: dict[str, Path], out_path: Path, vmin: float, vmax: float, legend_width_ratio: float) -> None:
    imgs = {k: Image.open(p).convert("RGB") for k, p in image_paths.items()}
    tw, th = imgs["top"].size
    for k in imgs:
        if imgs[k].size != (tw, th):
            imgs[k] = imgs[k].resize((tw, th), resample=Image.Resampling.BICUBIC)

    grid_w = 2 * tw
    grid_h = 2 * th
    legend_w = int(max(240, legend_width_ratio * tw))
    canvas = Image.new("RGB", (grid_w + legend_w, grid_h), (245, 245, 245))

    canvas.paste(imgs["front"], (0, 0))
    canvas.paste(imgs["top"], (tw, 0))
    canvas.paste(imgs["left"], (0, th))
    canvas.paste(imgs["perspective"], (tw, th))

    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    draw.text((8, 8), "FRONT", fill=(10, 10, 10), font=font)
    draw.text((tw + 8, 8), "TOP", fill=(10, 10, 10), font=font)
    draw.text((8, th + 8), "LEFT", fill=(10, 10, 10), font=font)
    draw.text((tw + 8, th + 8), "PERSPECTIVE", fill=(10, 10, 10), font=font)

    x0 = grid_w
    y0 = 0
    w = legend_w
    h = grid_h
    bar_x0 = x0 + int(0.40 * w)
    bar_x1 = x0 + int(0.64 * w)
    bar_y0 = y0 + int(0.08 * h)
    bar_y1 = y0 + int(0.92 * h)

    for y in range(bar_y0, bar_y1 + 1):
        t = 1.0 - (y - bar_y0) / max(1, (bar_y1 - bar_y0))
        c = rainbow_color_u8(t)
        draw.line((bar_x0, y, bar_x1, y), fill=c)

    draw.rectangle((bar_x0, bar_y0, bar_x1, bar_y1), outline=(25, 25, 25), width=2)
    txt_x = bar_x1 + 10
    draw.text((txt_x, bar_y0 - 8), f"{vmax:+.3f} eV", fill=(20, 20, 20), font=font)
    draw.text((txt_x, (bar_y0 + bar_y1) // 2 - 8), f"{0.0:+.3f} eV", fill=(20, 20, 20), font=font)
    draw.text((txt_x, bar_y1 - 8), f"{vmin:+.3f} eV", fill=(20, 20, 20), font=font)

    draw.text((x0 + 8, y0 + 8), f"{species} key", fill=(20, 20, 20), font=font)
    draw.text((x0 + 8, y0 + 30), "violet=positive", fill=(20, 20, 20), font=font)
    draw.text((x0 + 8, y0 + 50), "green=zero", fill=(20, 20, 20), font=font)
    draw.text((x0 + 8, y0 + 70), "red=negative", fill=(20, 20, 20), font=font)
    draw.text((x0 + 8, y0 + 90), "gray=host/missing", fill=(20, 20, 20), font=font)

    canvas.save(out_path)


def render_species(species_dir: Path, species: str, host_atoms: AtomsType, point_data: dict[int, dict[str, object]], args: argparse.Namespace) -> tuple[int, int]:
    # Build candidate list sorted by index for stable visuals.
    items = sorted(point_data.items(), key=lambda x: x[0])
    if not items:
        return 0, 0

    done_vals = [as_float(v.get("delta", 0.0)) for _, v in items if bool(v.get("done", False))]
    done_count = len(done_vals)
    total_count = len(items)
    if done_count < args.min_done:
        return done_count, total_count

    vabs = max(abs(min(done_vals)), abs(max(done_vals))) if done_vals else 1.0
    if vabs <= 1e-12:
        vabs = 1.0
    vmin, vmax = -vabs, vabs

    interstitial_fracs = [as_frac3(v.get("frac", (0.0, 0.0, 0.0))) for _, v in items]
    interstitial_atoms = Atoms(
        symbols=[species] * len(interstitial_fracs),
        scaled_positions=interstitial_fracs,
        cell=host_atoms.cell,
        pbc=True,
    )

    combo = host_atoms.copy() + interstitial_atoms

    host_n = len(host_atoms)
    point_n = len(interstitial_atoms)
    host_alpha = max(0.0, min(1.0, float(args.host_alpha)))
    dopant_alpha = max(0.0, min(1.0, float(args.dopant_alpha)))
    host_color = (0.62, 0.62, 0.62, host_alpha)
    missing_color = (0.45, 0.45, 0.45, dopant_alpha)

    colors: list[tuple[float, float, float, float]] = [host_color] * host_n
    for _, info in items:
        if bool(info.get("done", False)):
            colors.append(value_to_rgba(as_float(info.get("delta", 0.0)), vabs, dopant_alpha))
        else:
            colors.append(missing_color)

    z = atomic_numbers.get(species, 14)
    dopant_radius = float(covalent_radii[z]) * float(args.dopant_radius_scale)
    radii = [float(args.host_radius)] * host_n + [dopant_radius] * point_n

    targets = [species_dir / n for n in ["top.png", "front.png", "left.png", "perspective.png", args.output_name]]
    make_backup(species_dir, targets, args.backup_mode)

    for view_name, rot in VIEW_SPECS.items():
        render_view_png(combo, colors, radii, rot, species_dir / f"{view_name}.png", args.dpi)

    draw_concat_with_legend(
        species=species,
        image_paths={k: species_dir / f"{k}.png" for k in ["top", "front", "left", "perspective"]},
        out_path=species_dir / args.output_name,
        vmin=vmin,
        vmax=vmax,
        legend_width_ratio=args.legend_width_ratio,
    )
    return done_count, total_count


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.manifest).resolve()
    base_poscar = Path(args.base_poscar).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not manifest_path.exists():
        print(f"ERROR: missing manifest: {manifest_path}")
        return 1
    if not base_poscar.exists():
        print(f"ERROR: missing base POSCAR: {base_poscar}")
        return 1

    host_atoms = ensure_atoms(read(str(base_poscar)))
    manifest = load_manifest(manifest_path)
    by_species = collect_interstitial_data(manifest, args.site)

    if args.species == ["all"]:
        species_list = sorted(by_species)
    else:
        species_list = [s for s in args.species if s in by_species]

    summary_rows = []
    created = 0
    for species in species_list:
        species_dir = outdir / species
        species_dir.mkdir(parents=True, exist_ok=True)
        done_count, total_count = render_species(species_dir, species, host_atoms, by_species[species], args)
        if total_count == 0 or done_count < args.min_done:
            summary_rows.append(
                {
                    "species": species,
                    "done": done_count,
                    "total": total_count,
                    "status": "skipped",
                    "concat": "",
                }
            )
            continue

        created += 1
        summary_rows.append(
            {
                "species": species,
                "done": done_count,
                "total": total_count,
                "status": "created",
                "concat": str((species_dir / args.output_name).resolve()),
            }
        )

    with (outdir / "summary.csv").open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=["species", "done", "total", "status", "concat"])
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"created={created}")
    print(f"summary={outdir / 'summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
