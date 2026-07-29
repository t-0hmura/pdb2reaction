# `pdb2reaction scan3d`

## Purpose

3D distance scan with harmonic restraints. Drives three bonds toward
target distances and produces a 3D grid of restrained optimized geometries.
Use it when a three-coordinate grid is required; lower-dimensional scans need
fewer attempted optimizations but answer a different question. A distance grid does not by itself model
electronic-state changes or establish a mechanism.

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
| `-i, --input` | path | required (unless `--csv`) | Reactant `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required (unless `--csv`) | Python literal with **three** 4-tuples `(i, j, low, high)` |
| `--csv` | path | none | Skip the scan; load a precomputed `surface.csv` for downstream plotting |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan3d/` | Output directory |
| `--ref-pdb` / `--config` / `--dry-run` / `--help-advanced` | — | — | Standard |
| `--out-json / --no-out-json` | flag | `--no-out-json` | Emit `result.json` summary |

## Examples

```bash
pdb2reaction scan3d -i 1.R.pdb -l 'SAM:1' \
    -s '[("OH TYR 100","HC TYR 100",1.50,2.40), ("HC TYR 100","O ASP 50",1.20,2.20), ("FE HEM 200","O ASP 50",2.10,3.10)]' \
    -b uma -o result_scan3d
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/grid/point_iDDD_jDDD_kDDD.{xyz,pdb,cif,gjf}` (DDD = round(d×100), Å) | each attempted grid point reaches output | final attempted grid geometry; format companions depend on the input topology and `--convert-files`; `surface.csv` records `bias_converged` |
| `<out_dir>/grid/preopt_iDDD_jDDD_kDDD.{xyz,pdb,cif,gjf}` (DDD = round(d×100), Å) | a new scan is run | starting-structure snapshot (actually preoptimized only with `--preopt`); format companions depend on the input topology and `--convert-files` |
| `<out_dir>/grid/inner_path_d1_NNN_d2_MMM_trj.xyz` | `--dump` | inner-loop trajectory |
| `<out_dir>/scan3d_density.html` | scan/CSV data pass validation and interpolation | interactive 3D iso-surface |
| `<out_dir>/surface.csv` | unless `--csv` (post-mortem reuses the input csv) | 3D energy surface (i, j, k, d1_A, d2_A, d3_A, energy_hartree, bias_converged, energy_kcal, plus axis labels) |

`result.json` stores grid metadata and file paths; `surface.csv` holds the
per-grid-point energies and `bias_converged` values. In `--csv` post-mortem
mode, the input CSV is read but not copied to `<out_dir>`, and the result JSON
does not advertise a new `surface_csv` file.
When rounded distance tags collide, later point filenames append
`_grid_III_JJJ_KKK` with zero-based grid indices.

## Caveats

- Output volume scales as `n_steps1 × n_steps2 × n_steps3`; a 5×5×5 grid
  attempts 125 restrained geometry optimizations.
- Use `scan.md` / `scan2d.md` when fewer coordinates answer the question;
  their number of attempted grid optimizations scales in fewer dimensions.
- `--csv` is the post-mortem entry point for regenerating the 3D density
  visualization without redoing the scan.

## See also

- `scan.md`, `scan2d.md` — lower-dim analogs.
- `all-scan-list.md` — staged sequential scans inside the full
  pipeline (avoids the 3D grid blow-up when stages are decoupled).
- Defaults: `import pdb2reaction.core.defaults as d; print(d.OUT_DIR_SCAN3D)`
