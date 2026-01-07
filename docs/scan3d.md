# `scan3d` subcommand

## Overview
`scan3d` performs a three-distance grid scan with harmonic restraints using the
UMA calculator. You provide exactly one `--scan-list` literal containing three
quadruples `(i, j, lowÅ, highÅ)`. The tool builds linear grids for each distance
with `--max-step-size`, reorders the values so that those nearest to the
(pre‑optimized) starting structure are visited first, and then nests the loops
(outer d₁, middle d₂, inner d₃). Each grid point is relaxed with the
appropriate restraints active; unbiased energies are recorded so you can compare
points directly. A precomputed `surface.csv` can also be visualized without
rerunning the scan.

## Usage
```bash
pdb2reaction scan3d -i INPUT.{pdb|xyz|trj|...} -q CHARGE [-m MULT] \
                    --scan-list '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]' [options] \
                    [--convert-files {True|False}]
```

### Examples
```bash
# Minimal three-distance scan
pdb2reaction scan3d -i input.pdb -q 0 \
    --scan-list '[("TYR,285,CA",1.30,3.10),("MMT,309,C10",1.20,3.20),("TYR,285,CB",1.10,3.00)]'

# LBFGS relaxations, dumped inner trajectories, and HTML isosurface plot
pdb2reaction scan3d -i input.pdb -q 0 \
    --scan-list '[("TYR,285,CA",1.30,3.10),("MMT,309,C10",1.20,3.20),("TYR,285,CB",1.10,3.00)]' \
    --max-step-size 0.20 --dump True --out-dir ./result_scan3d/ --opt-mode light \
    --preopt True --baseline min

# Plot only from an existing surface.csv (skip new energy evaluation)
pdb2reaction scan3d -i input.pdb -q 0 \
    --scan-list '[("TYR,285,CA",1.30,3.10),("MMT,309,C10",1.20,3.20),("TYR,285,CB",1.10,3.00)]' \
    --csv ./result_scan3d/surface.csv --out-dir ./result_scan3d/
```

## Workflow
1. Load the structure through `geom_loader`, resolve charge/spin from CLI or
   embedded Gaussian templates, and optionally run an unbiased preoptimization
   when `--preopt True`. If `-q` is omitted but `--ligand-charge` is provided, the
   structure is treated as an enzyme–substrate complex and `extract.py`’s charge
   summary derives the total charge before scanning.
2. Parse the single `--scan-list` literal (default 1-based indices unless
   `--one-based False` is passed) into three quadruples. For PDB inputs, each
   atom entry can be an integer index or a selector string like `'TYR,285,CA'`;
   delimiters may be spaces, commas, slashes, backticks, or backslashes, and
   token order is flexible (fallback assumes resname, resseq, atom). Build each linear grid using
   `h = --max-step-size` and reorder the values so the ones closest to the
   starting distances are visited first.
3. Outer loop over `d1[i]`: relax with only the d₁ restraint active, starting
   from the previously scanned geometry whose d₁ value is closest. Snapshot that
   structure.
4. Middle loop over `d2[j]`: relax with d₁ and d₂ restraints, starting from the
   closest (d₁, d₂) geometry. Snapshot that result.
5. Inner loop over `d3[k]`: relax with all three restraints, measure the
   unbiased energy (bias removed for evaluation), and write the constrained
   geometry and convergence flag.
6. After the scan completes, assemble `surface.csv`, apply the kcal/mol
   baseline shift (`--baseline {min|first}`), and generate a 3D RBF-interpolated
   isosurface plot (`scan3d_density.html`) honoring `--zmin/--zmax`. When
   `--csv` is provided, only this plotting step runs.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge (CLI > template > 0). Overrides `--ligand-charge` when both are set. | Required when not in template |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1. | `1` |
| `--scan-list TEXT` | **Single** Python literal with three quadruples `(i,j,lowÅ,highÅ)`. `i`/`j` can be integer indices or PDB atom selectors like `'TYR,285,CA'`. | Required |
| `--one-based {True|False}` | Interpret `(i, j)` indices as 1- or 0-based. | `True` |
| `--max-step-size FLOAT` | Maximum change allowed per distance increment (Å). Controls grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Maximum optimizer cycles during each biased relaxation. Overrides `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS, `heavy` → RFOptimizer. | `light` |
| `--freeze-links BOOL` | When the input is PDB, freeze parents of link hydrogens. | `True` |
| `--dump BOOL` | Write `inner_path_d1_###_d2_###.trj` for each (d₁, d₂). | `False` |
| `--convert-files {True|False}` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--out-dir TEXT` | Output directory root for grids and plots. | `./result_scan3d/` |
| `--csv PATH` | Load an existing `surface.csv` and only plot it (no new scan). | _None_ |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--args-yaml FILE` | YAML overrides for `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`. | _None_ |
| `--preopt BOOL` | Run an unbiased optimization before scanning. | `True` |
| `--baseline {min,first}` | Shift kcal/mol energies so the global min or `(i,j,k)=(0,0,0)` is zero. | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | Manual limits for the isosurface color bands (kcal/mol). | Autoscaled |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical knobs to those documented for
  [`opt`](opt.md#yaml-configuration-args-yaml). `opt.dump` is forced to `False`
  so trajectory control stays on the CLI.

More YAML options about `opt` are available in [docs/opt.md](opt.md#yaml-
configuration-args-yaml).

## YAML configuration (`--args-yaml`)
A minimal example (extend using the keys documented for [`opt`](opt.md#yaml-
configuration-args-yaml)):

```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  device: auto               # UMA device selection
opt:
  thresh: gau                # convergence preset (Gaussian/Baker-style)
  max_cycles: 10000          # optimizer cycle cap
  dump: false                # trajectory dumping disabled (CLI controls dumping)
  out_dir: ./result_scan3d/  # output directory
lbfgs:
  max_step: 0.3              # maximum step length
  out_dir: ./result_scan3d/  # LBFGS-specific output directory
rfo:
  trust_radius: 0.3          # trust-region radius
  out_dir: ./result_scan3d/  # RFO-specific output directory
bias:
  k: 100.0                  # harmonic bias strength (eV·Å⁻²)
```

More YAML options about `opt` are available in [docs/opt.md](opt.md).

### Section `bias`
- `k` (`100`): Harmonic strength in eV·Å⁻². Overridden by `--bias-k`.

## Outputs
```
out_dir/ (default: ./result_scan3d/)
├─ surface.csv                     # Grid metadata; may include a reference row (i=j=k=-1)
├─ scan3d_density.html             # 3D energy isosurface visualization
├─ grid/point_i###_j###_k###.xyz   # Relaxed geometry for each grid point (Å×100 tags)
├─ grid/point_i###_j###_k###.pdb   # PDB companions when conversion is enabled and templates exist
├─ grid/point_i###_j###_k###.gjf   # Gaussian companions when templates exist and conversion is enabled
├─ grid/preopt_i###_j###_k###.xyz  # Starting structure saved before scanning (preoptimized when --preopt is True)
└─ grid/inner_path_d1_###_d2_###.trj # Present only when --dump is True (mirrored to .pdb/.gjf with conversion)
```

## Notes
- UMA via `uma_pysis` is the only calculator backend and reuses the same
  `HarmonicBiasCalculator` as the 1D/2D scans.
- Ångström limits are converted to Bohr internally to cap LBFGS steps and RFO
  trust radii; optimizer scratch files live under temporary directories.
- `--baseline` defaults to the global minimum; `--baseline first` anchors the
  `(i,j,k)=(0,0,0)` grid point when present.
- 3D visualization uses RBF interpolation on a 50×50×50 grid with
  semi-transparent step-colored isosurfaces (no cross-sectional planes).
- `--freeze-links` merges user `freeze_atoms` with detected link-H parents for
  PDB inputs, keeping extracted pockets rigid.
