# Defect Generator

`generate_defects.py` creates defect variants from a base POSCAR and prepares each case for VASP.

## Defaults

- Base structure: `CALC/VASP/unrelaxed/relax/POSCAR`
- Output root: `__variations__/defects`
- Defect types: `vacancy`, `substitutional`
- Canonical interstitial workflow: tetra-only via `generate_interstitial_tetra.py`
- Default multiplicity: `n=1`
- Default combination mode: exhaustive (`--max-variants-per-case -1`)
- Per-case preparation (always):
  1. `potcar.py`
  2. `nbands.py`
  3. `incar.py`
  4. `kpoints.py`

## Output Layout

- Vacancy: `__variations__/defects/vacancy/<n>-site/vXXXX/`
- Substitutional: `__variations__/defects/substitutional/<n>-site/<Element>/vXXXX/`

Tetra interstitial outputs are generated separately:

- Interstitial (tetra): `__variations__/defects_tetra/interstitial/<n>-site/<Element>/vXXXX/`

Each variant directory contains at least:

- `POSCAR`
- `POTCAR`
- `nbands.txt`
- `nelect.txt`
- `INCAR`
- `KPOINTS`

Run-level files are written to `__variations__/defects/`:

- `manifest.json`
- `manifest.csv`
- `summary.json`

Tetra interstitial run-level files are written to `__variations__/defects_tetra/`:

- `manifest_interstitial_tetra.json`
- `manifest_interstitial_tetra.csv`
- `summary_interstitial_tetra.json`

## Typical Runs

Exhaustive `n=1` (default):

```bash
python3 generate_defects.py
```

Exhaustive `n=1,2`:

```bash
python3 generate_defects.py --n-defects 1 2 --max-variants-per-case -1
```

Dry run only:

```bash
python3 generate_defects.py --dry-run
```

## Notes

- The script filters elements against direct POTCAR availability in
  `~/research/projects/poly_axis_interfaces/src/tools/vasp_utilities/common/POTCARS/<Element>/POTCAR`.
- If a required output is missing for a variant, that variant is marked as failed in the manifest.
- Legacy grid interstitial outputs under `__variations__/defects/interstitial/` are archived; use the tetra workflow instead.
