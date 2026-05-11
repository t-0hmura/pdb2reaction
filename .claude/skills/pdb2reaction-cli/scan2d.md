# `pdb2reaction scan2d`

## Purpose

2D distance scan with harmonic restraints. Drives two bonds toward
target distances simultaneously and produces a 2D grid of relaxed
geometries. Useful for mapping concerted-vs-stepwise reaction
surfaces (e.g. SN2 attack + leaving-group departure).

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
| `-i, --input` | path | required | Reactant `.pdb` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required | Inline Python literal containing **two** 4-tuples `(i, j, low, high)` (one per axis), or a YAML/JSON spec file. |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan2d/` | Output directory |
| `--ref-pdb` | path | none | Residue context for XYZ/GJF |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

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
| `<out_dir>/grid/point_i<d1Å>_j<d2Å>.xyz` | always | final relaxed grid point |
| `<out_dir>/grid/preopt_i<d1Å>_j<d2Å>.{xyz,pdb,gjf}` | `--preopt` | pre-relaxation snapshot |
| `<out_dir>/grid/inner_path_d1_NNN_trj.xyz` | `--dump` | inner-loop trajectory |
| `<out_dir>/scan2d_map.png` | always | 2D energy surface heatmap |
| `<out_dir>/scan2d_landscape.html` | always | interactive 3D landscape |
| `<out_dir>/surface.csv` | always | 2D energy surface (d1_A, d2_A, energy_hartree, energy_kcal, plus axis labels) |

`result.json` stores grid metadata and energy values; `surface.csv` is
ready for downstream contour plotting; `scan2d_map.png` is the static
2D heatmap; `scan2d_landscape.html` is an interactive 3D rendering.

## Caveats

- `-s` literal must contain **exactly two** tuples for `scan2d`.
- Output volume scales as `n_steps1 × n_steps2`. Keep grids modest;
  10×10 = 100 single-point optimizations.
- Atom-spec syntax is identical to `scan.md`.

## See also

- `scan.md`, `scan3d.md` — 1D / 3D analogs.
- `all-scan-list.md` — sequential scans (different from coupled grid).
- Defaults: `import pdb2reaction.defaults as d; print(d.BIAS_KW, d.OUT_DIR_SCAN2D)`