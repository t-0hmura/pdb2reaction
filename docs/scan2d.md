# `2d-scan` subcommand

## Purpose
Performs a two-distance, two-dimensional scan with harmonic restraints on each
pair while using UMA-based LBFGS/RFO relaxations under the hood. The inner loop
freezes the first distance and explores the second one so you obtain a PES grid
whose energies are always evaluated without bias.

## Usage
```bash
pdb2reaction 2d-scan -i INPUT -q CHARGE \
                   --scan-list "[(i1,j1,low1,high1),(i2,j2,low2,high2)]" \
                   [--one-based/--zero-based] [--max-step-size ΔÅ] [--bias-k k] \
                   [--relax-max-cycles N] [--opt-mode light|lbfgs|heavy|rfo] \
                   [--freeze-links BOOL] [--dump BOOL] [--out-dir DIR] \
                   [--thresh PRESET] [--args-yaml FILE] [--preopt BOOL] \
                   [--baseline min|first] [--zmin FLOAT] [--zmax FLOAT]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | `.gjf` template value or `0` |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--scan-list TEXT` | **Single** Python-like literal with two `(i, j, lowÅ, highÅ)` quadruples, one per distance. | Required |
| `--one-based / --zero-based` | Interpret `(i, j)` indices in `--scan-list` as 1-based or 0-based. | `--one-based` |
| `--max-step-size FLOAT` | Largest change allowed for either distance per grid increment. Determines the number of points. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` (eV·Å⁻²). Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Maximum optimizer cycles for every biased relaxation. | `10000` |
| `--opt-mode TEXT` | Relaxation backend (`light|lbfgs` for LBFGS, `heavy|rfo` for RFOptimizer). | `light` |
| `--freeze-links BOOL` | When the input is PDB, also freeze the parent atoms of link hydrogens. | `True` |
| `--dump BOOL` | If `True`, write `inner_path_d1_###.trj` files with the inner (d2) scan snapshots. | `False` |
| `--out-dir TEXT` | Output directory for grids/plots. | `./result_scan2d/` |
| `--thresh TEXT` | Override the UMA convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--args-yaml FILE` | YAML overrides for `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias` sections. | _None_ |
| `--preopt BOOL` | Perform an unbiased optimization before scanning. | `True` |
| `--baseline {min,first}` | Reference for `energy_kcal` in `surface.csv`. `min` zeroes the global minimum, `first` zeroes (i=0,j=0). | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | Manually clamp the color scale used for the contour/surface plots (kcal/mol). | Autoscaled |

## Algorithm overview
1. Load the structure via `geom_loader` and resolve charge/spin (CLI overrides >
   embedded template > defaults).
2. Optionally run an unbiased LBFGS/RFO pre-optimization with the UMA
   calculator configured by `--args-yaml`.
3. Parse the two quadruples from `--scan-list`, normalize them to 0-based indices
   (unless `--zero-based`), and build two linear grids according to
   `--max-step-size`.
4. Loop over each value of the first distance (`d1`), relaxing the geometry with
   only that bias applied. For every relaxed `d1`, snapshot the structure and
   scan the second distance (`d2`) while keeping `d1` fixed.
5. At every `(d1_i, d2_j)` pair, run a biased optimization, record the unbiased
   UMA single-point energy, and store the relaxed coordinates as
   `grid/point_i###_j###.xyz`. Optional inner trajectories are dumped as
   `grid/inner_path_d1_###.trj` when `--dump True`.
6. Write all grid records into `<out-dir>/surface.csv` with columns
   `i,j,d1_A,d2_A,energy_hartree,energy_kcal,bias_converged`. Energies are
   converted to kcal/mol relative to the baseline requested via `--baseline`.
7. Generate square Plotly figures inside `<out-dir>/plots/`:
   - `scan2d_contour.png` (or `.html` fallback) for the 2D contour heatmap.
   - `scan2d_surface.html` for the 3D surface plus projected contour.

## YAML configuration (`--args-yaml`)
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical keys to the ones documented
  for [`opt`](opt.md#yaml-configuration-args-yaml). `opt.dump` is internally
  forced to `False`; use `--dump` to control trajectory output.
- `bias.k`: Harmonic strength in eV·Å⁻². Overwritten by `--bias-k` when present.

A minimal example:
```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
  model: uma-s-1p1
  device: auto
opt:
  thresh: gau
  max_cycles: 10000
  out_dir: ./result_scan2d/
lbfgs:
  max_step: 0.3
rfo:
  trust_radius: 0.3
bias:
  k: 100.0
```
