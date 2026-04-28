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
| `--ref-pdb` / `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

```bash
pdb2reaction scan3d -i 1.R.pdb -l 'SAM:1' \
    -s '[("OH TYR 100","HC TYR 100",1.50,2.40), ("HC TYR 100","O ASP 50",1.20,2.20), ("FE 200","O ASP 50",2.10,3.10)]' \
    -b uma -o result_scan3d
```

## Output

```
result_scan3d/
├── result.json
├── grid_NN_MM_LL/              # per grid-point relaxed geometries
│   └── final.xyz
├── surface.csv                 # 3D energy surface (axis_1, axis_2, axis_3, energy)
└── scan3d.log
```

`result.json` stores grid metadata and energy values; `surface.csv`
holds the four-column tabulation ready for slicing / contour plotting
in pandas or matplotlib.

## Caveats

- Output volume scales as `n_steps1 × n_steps2 × n_steps3`. Keep grid
  sizes modest (5×5×5 = 125 points already).
- For inherently 1D or 2D mechanisms, `scan.md` / `scan2d.md` are
  drastically cheaper.
- `--csv` is the post-mortem entry point: rerun energy diagrams or
  contour plots without redoing the scan.

## See also

- `scan.md`, `scan2d.md` — lower-dim analogs.
- `all-scan-list.md` — staged sequential scans inside the full
  pipeline (avoids the 3D grid blow-up when stages are decoupled).
- Defaults: `import pdb2reaction.defaults as d; print(d.OUT_DIR_SCAN3D)`