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
    -s '[(a1,b1,T1), (a2,b2,T2), (a3,b3,T3)]' \
    [--csv surface.csv] \
    [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] [-o ./result_scan3d/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required (unless `--csv`) | Reactant `.pdb` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required (unless `--csv`) | Python literal with **three** tuples |
| `--csv` | path | none | Skip the scan; load a precomputed `surface.csv` for downstream plotting |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan3d/` | Output directory |
| `--ref-pdb` / `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

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