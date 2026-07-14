# `pdb2reaction scan2d`

## Purpose

2D distance scan with harmonic restraints. Drives two bonds toward
target distances simultaneously and produces a 2D grid of relaxed
geometries. It can help inspect coupling between two chosen coordinates, but
the grid alone does not prove a concerted or stepwise mechanism.

## Synopsis

```bash
pdb2reaction scan2d -i input.pdb \
    -s '[(a1, b1, low_A, high_A), (a2, b2, low_B, high_B)]' \
    [-l 'RES:Q,...'] \
    [-b uma|orb|mace|aimnet2] [-o ./result_scan2d/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Reactant `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required | Inline Python literal containing **two** 4-tuples `(i, j, low, high)` (one per axis), or a YAML/JSON spec file. |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan2d/` | Output directory |
| `--ref-pdb` | path | none | Residue context for XYZ/GJF |
| `--config` / `--dry-run` / `--help-advanced` | — | — | Standard |

Each 4-tuple defines one scan axis: atoms `i,j` (1-based indices or
PDB-style atom specs) and the distance range `[low, high]` in Å. Both
bonds are driven across their ranges simultaneously, generating the
grid.

## Examples

```bash
pdb2reaction scan2d -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60,3.10), ("H11 GPP 321","OE2 GLU 186",0.90,2.40)]' \
    -b uma -o result_scan2d
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/grid/point_iDDD_jDDD.{xyz,pdb,cif,gjf}` (DDD = round(d×100), Å) | each attempted grid point reaches output | final attempted grid geometry; format companions depend on the input topology and `--convert-files`; `surface.csv` records `bias_converged` |
| `<out_dir>/grid/preopt_iDDD_jDDD.{xyz,pdb,cif,gjf}` (DDD = round(d×100), Å) | `--preopt` | pre-relaxation snapshot; format companions depend on the input topology and `--convert-files` |
| `<out_dir>/grid/inner_path_d1_NNN_trj.xyz` | `--dump` | inner-loop trajectory |
| `<out_dir>/scan2d_map.png` | successful grid + interpolation/PNG export | 2D energy surface heatmap |
| `<out_dir>/scan2d_landscape.html` | successful grid + interpolation | interactive 3D landscape |
| `<out_dir>/surface.csv` | at least one grid record | 2D energy surface (i, j, d1_A, d2_A, energy_hartree, bias_converged, energy_kcal, plus axis labels) |

`result.json` stores grid metadata and file paths; per-point energies and
convergence flags are in `surface.csv`. `result.json["status"] == "completed"`
means the grid and plots were written, not that every `bias_converged` value is
true.

## Caveats

- `-s` literal must contain **exactly two** tuples for `scan2d`.
- Output volume scales as `n_steps1 × n_steps2`; a 10×10 grid attempts 100
  restrained geometry optimizations.
- Atom-spec syntax is identical to `scan.md`.

## See also

- `scan.md`, `scan3d.md` — 1D / 3D analogs.
- `all-scan-list.md` — sequential scans (different from coupled grid).
- Defaults: `import pdb2reaction.core.defaults as d; print(d.BIAS_KW, d.OUT_DIR_SCAN2D)`
