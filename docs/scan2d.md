# `scan2d` subcommand

## Overview
`scan2d` performs a two-distance (d₁, d₂) grid scan using harmonic restraints
and UMA-based relaxations. You supply one `--scan-list` literal with two
quadruples `(i, j, lowÅ, highÅ)`; the tool constructs linear grids for both
ranges using `--max-step-size`, then **reorders each axis so that points closest
to the (pre)optimized structure are visited first**. Each grid point is relaxed
and written alongside a ready-to-plot CSV/figure bundle. Energies reported in
`surface.csv` are always evaluated **without bias** so you can compare grid
points directly. Optimizations use LBFGS when `--opt-mode light` (default)
or RFOptimizer when `--opt-mode heavy`.

## Usage
```bash
pdb2reaction scan2d -i INPUT.{pdb|xyz|trj|...} -q CHARGE [-m MULT] \
                    --scan-list "[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]" [options]
                    [--convert-files {True|False}]
```

### Examples
```bash
# Minimal two-distance scan
pdb2reaction scan2d -i input.pdb -q 0 \
    --scan-list "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

# LBFGS, dumped inner trajectories, and Plotly outputs
pdb2reaction scan2d -i input.pdb -q 0 \
    --scan-list "[(12,45,1.30,3.10),(10,55,1.20,3.20)]" \
    --max-step-size 0.20 --dump True --out-dir ./result_scan2d/ --opt-mode light \
    --preopt True --baseline min
```

## Workflow
1. Load the input geometry via `geom_loader`, resolve charge/spin, and optionally
   run an unbiased preoptimization when `--preopt True`. If `-q` is omitted but
   `--ligand-charge` is provided, the structure is treated as an enzyme–substrate
   complex and `extract.py`’s charge summary derives the total charge before the
   scan. The preoptimized structure is saved under `grid/preopt_i###_j###.*` and
   its unbiased energy is stored in `surface.csv` with indices `i = j = -1`.
2. Parse the single `--scan-list` literal into two quadruples, normalize indices
   (1-based by default), and construct linear grids: `N = ceil(|high − low| / h)`
   with `h = --max-step-size`. Zero-length spans collapse to a single point.
   Each axis is then reordered so that the distance closest to the preoptimized
   geometry is indexed as `i = 0` / `j = 0`.
3. Iterate over every `d1[i]` (nearest-first ordering). For each value, relax
   the system with **only the d₁ restraint** active, snapshot that geometry, then
   run the inner loop over `d2[j]` with **both restraints** applied starting from
   the nearest previously converged structure.
4. At each `(i, j)` pair, store the biased-optimization result under
   `<out-dir>/grid/point_i###_j###.xyz`, record whether the bias converged, and
   evaluate the UMA energy without bias. Optional per-outer-step inner
   trajectories are saved as `inner_path_d1_###.trj` when `--dump True`.
5. After all points are visited, write `<out-dir>/surface.csv` with columns
   `i,j,d1_A,d2_A,energy_hartree,energy_kcal,bias_converged`, shifting the kcal
   reference via `--baseline {min|first}`. With `--baseline first`, the reference
   is the first grid entry (`i = j = 0` after reordering), not necessarily
   `(low₁, low₂)`. Generate `scan2d_map.png` (2D contour) and
   `scan2d_landscape.html` (3D surface) in `<out-dir>/`. Use `--zmin/--zmax` to
   clamp the color scale.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge (CLI > template > 0). Overrides `--ligand-charge` when both are set. | Required when not in template |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1 (CLI > template > 1). | `.gjf` template value or `1` |
| `--scan-list TEXT` | **Single** Python literal with two quadruples `(i,j,lowÅ,highÅ)`. | Required |
| `--one-based {True|False}` | Interpret `(i, j)` indices as 1- or 0-based. | `True` |
| `--max-step-size FLOAT` | Maximum change allowed for either distance per increment (Å). Determines the grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Maximum optimizer cycles during each biased relaxation. Overrides `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS, `heavy` → RFOptimizer. | `light` |
| `--freeze-links BOOL` | When the input is PDB, freeze parents of link hydrogens. | `True` |
| `--dump BOOL` | Write `inner_path_d1_###.trj` for each outer step. | `False` |
| `--convert-files {True|False}` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--out-dir TEXT` | Output directory root for grids and plots. | `./result_scan2d/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--args-yaml FILE` | YAML overrides for `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`. | _None_ |
| `--preopt BOOL` | Run an unbiased optimization before scanning. | `True` |
| `--baseline {min,first}` | Shift kcal/mol energies so the global min or first grid point is zero. | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | Manual limits for the contour/surface color scale (kcal/mol). | Autoscaled |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical knobs to those documented for
  [`opt`](opt.md#yaml-configuration-args-yaml). `opt.dump` is forced to `False`
  so trajectory control stays on the CLI.

### Section `bias`
- `k` (`100`): Harmonic strength in eV·Å⁻². Overridden by `--bias-k`.

## Outputs
```
out_dir/ (default: ./result_scan2d/)
├─ surface.csv                # Structured grid table
├─ scan2d_map.png             # 2D contour (requires Kaleido; the run stops if PNG export fails)
├─ scan2d_landscape.html      # 3D surface visualization
├─ grid/point_i###_j###.xyz   # Relaxed geometries for every (i, j) pair
├─ grid/point_i###_j###.pdb   # PDB companions when conversion is enabled and templates exist
├─ grid/point_i###_j###.gjf   # Gaussian companions when templates exist and conversion is enabled
└─ grid/inner_path_d1_###.trj # Present only when --dump is True (mirrored to .pdb for PDB inputs with conversion)
```

## Notes
- UMA via `uma_pysis` is the only calculator backend and reuses the same
  `HarmonicBiasCalculator` as the 1D scan.
- Ångström limits are converted to Bohr internally to cap LBFGS steps and RFO
  trust radii; Optimizer scratch files live under temporary directories.
- The bias is always removed before final energies are recorded so you can reuse
  `surface.csv` in downstream fitting or visualization scripts.
- `--freeze-links` merges user `freeze_atoms` with detected link-H parents for
  PDB inputs, keeping extracted pockets rigid.
- Charge inherits Gaussian template metadata when available. If `-q` is omitted
  but `--ligand-charge` is provided, the full structure is treated as an
  enzyme–substrate complex and `extract.py`’s charge summary computes the total
  charge; explicit `-q` still overrides. Otherwise charge defaults to `0` and
  spin to `1`.

## YAML configuration (`--args-yaml`)
A minimal example (extend with the same keys documented in [`opt`](opt.md#yaml-
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
  out_dir: ./result_scan2d/  # output directory
lbfgs:
  max_step: 0.3              # maximum step length
  out_dir: ./result_scan2d/  # LBFGS-specific output directory
rfo:
  trust_radius: 0.3          # trust-region radius
  out_dir: ./result_scan2d/  # RFO-specific output directory
bias:
  k: 100.0                  # harmonic bias strength (eV·Å⁻²)
```

More YAML options about `opt` are available in [docs/opt.md](opt.md).
