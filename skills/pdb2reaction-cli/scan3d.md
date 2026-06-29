# `pdb2reaction scan3d`

## Purpose

3D distance scan with harmonic restraints. Drives three bonds toward
target distances and produces a 3D grid of relaxed geometries. Rare in
practice; usually a sequence of 1D / 2D scans (or `all-scan-list.md`)
captures the chemistry with less compute. Provided for inherently
3D-coupled mechanisms (e.g. proton-coupled electron transfer with two
protons + one redox-donor distance).

## Synopsis

```bash
pdb2reaction scan3d -i input.pdb \
    -s '[(a1,b1,low_1,high_1), (a2,b2,low_2,high_2), (a3,b3,low_3,high_3)]' \
    [--csv surface.csv] \
    [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] [-o ./result_scan3d/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required (unless `--csv`) | Reactant `.pdb` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required (unless `--csv`) | Python literal with **three** 4-tuples `(i, j, low, high)` |
| `--csv` | path | none | Skip the scan; load a precomputed `surface.csv` for downstream plotting |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan3d/` | Output directory |
| `--ref-pdb` / `--config` / `--help-advanced` | — | — | Standard |
| `--out-json / --no-out-json` | flag | `--no-out-json` | Emit `result.json` summary |

## Examples

```bash
pdb2reaction scan3d -i 1.R.pdb -l 'SAM:1' \
    -s '[("OH TYR 100","HC TYR 100",1.50,2.40), ("HC TYR 100","O ASP 50",1.20,2.20), ("FE 200","O ASP 50",2.10,3.10)]' \
    -b uma -o result_scan3d
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/grid/point_i<d1Å>_j<d2Å>_k<d3Å>.xyz` | always | final relaxed grid point |
| `<out_dir>/grid/preopt_i<d1Å>_j<d2Å>_k<d3Å>.{xyz,pdb,gjf}` | always | starting-structure snapshot (preoptimized only when `--preopt` is set) |
| `<out_dir>/grid/inner_path_d1_NNN_d2_MMM_trj.xyz` | `--dump` | inner-loop trajectory |
| `<out_dir>/scan3d_density.html` | always (with input csv when `--csv` is passed) | interactive 3D iso-surface |
| `<out_dir>/surface.csv` | unless `--csv` (post-mortem reuses the input csv) | 3D energy surface (i, j, k, d1_A, d2_A, d3_A, energy_hartree, bias_converged, energy_kcal, plus axis labels) |

`result.json` stores grid metadata and energy values; `surface.csv`
holds the per-grid-point tabulation ready for slicing / contour plotting
in pandas or matplotlib.

## Caveats

- Output volume scales as `n_steps1 × n_steps2 × n_steps3`. Keep grid
  sizes modest (a 5×5×5 grid is already 125 points).
- For inherently 1D or 2D mechanisms, `scan.md` / `scan2d.md` are
  drastically cheaper.
- `--csv` is the post-mortem entry point: rerun energy diagrams or
  contour plots without redoing the scan.

## See also

- `scan.md`, `scan2d.md` — lower-dim analogs.
- `all-scan-list.md` — staged sequential scans inside the full
  pipeline (avoids the 3D grid blow-up when stages are decoupled).
- Defaults: `import pdb2reaction.core.defaults as d; print(d.OUT_DIR_SCAN3D)`