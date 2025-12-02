# `opt` subcommand

## Overview
`pdb2reaction opt` performs a single-structure geometry optimization with the pysisyphus LBFGS ("light") or RFOptimizer ("heavy") engines while UMA provides energies, gradients, and Hessians. Input structures can be `.pdb`, `.xyz`, `.trj`, or any format supported by `geom_loader`. Settings are applied in the order **built-in defaults → CLI overrides → `--args-yaml` overrides** (YAML has the highest precedence), making it easy to keep lightweight defaults while selectively overriding options. The optimizer preset now defaults to the LBFGS-based **`light`** mode.

When the starting structure is a PDB or Gaussian template, format-aware conversion mirrors outputs into `.pdb` or multi-geometry `.gjf` companions, controlled by `--convert-files/--no-convert-files` (enabled by default). PDB-specific conveniences include:
- With `--freeze-links` (default `True`), parent atoms of link hydrogens are detected and merged into `geom.freeze_atoms` (0-based indices).
- Output conversion produces `final_geometry.pdb` (and `optimization.pdb` when dumping trajectories) using the input PDB as the topology reference.

A Gaussian `.gjf` template seeds the charge/spin defaults and enables automatic export of the optimized structure and trajectories as `.gjf` when conversion is enabled.

## Usage
```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} -q CHARGE -m MULT \
                 [--opt-mode light|heavy] [--freeze-links BOOL] \
                 [--dist-freeze "[(i,j,target_A), ...]"] [--one-based|--zero-based] \
                 [--bias-k K_eV_per_A2] [--dump BOOL] [--out-dir DIR] \
                 [--max-cycles N] [--thresh PRESET] [--args-yaml FILE] \
                 [--convert-files/--no-convert-files]
```

## Workflow
- **Optimizers**: `--opt-mode light` (default) → L-BFGS; `--opt-mode heavy` → rational-function optimizer with trust-region control.
- **Restraints**: `--dist-freeze` consumes Python-literal tuples `(i, j, target_A)`; omitting the third element restrains the starting distance. `--bias-k` sets a global harmonic strength (eV·Å⁻²). Indices default to 1-based but can be flipped to 0-based with `--zero-based`.
- **Charge/spin resolution**: CLI `-q/-m` override `.gjf` template metadata, which in turn override the `calc` defaults. If no template exists, the fallback is `0/1`. Always pass the physically correct values explicitly.
- **Freeze atoms**: CLI freeze-link logic is merged with YAML `geom.freeze_atoms`, then propagated to the UMA calculator (`calc.freeze_atoms`).
- **Dumping & conversion**: `--dump True` mirrors `opt.dump=True` and writes `optimization.trj`; when conversion is enabled, trajectories are mirrored to `.pdb` for PDB inputs. `opt.dump_restart` can emit restart YAML snapshots.
- **Exit codes**: `0` success, `2` zero step (step norm < `min_step_norm`), `3` optimizer failure, `130` keyboard interrupt, `1` unexpected error.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Required unless the input is a `.gjf` template that already encodes charge. | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Falls back to `.gjf` template or `1`. | Template/`1` |
| `--dist-freeze TEXT` | Repeatable string parsed as Python literal describing `(i,j,target_A)` tuples for harmonic restraints. | _None_ |
| `--one-based / --zero-based` | Interpret `--dist-freeze` indices as 1-based (default) or 0-based. | `--one-based` |
| `--bias-k FLOAT` | Harmonic bias strength applied to every `--dist-freeze` tuple (eV·Å⁻²). | `10.0` |
| `--freeze-links BOOL` | Toggle link-hydrogen parent freezing (PDB inputs only). | `True` |
| `--max-cycles INT` | Hard limit on optimization iterations (`opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Choose optimizer: `light` (LBFGS) or `heavy` (RFO). | `light` |
| `--dump BOOL` | Emit trajectory dumps (`optimization.trj`). | `False` |
| `--convert-files/--no-convert-files` | Enable or disable XYZ/TRJ → PDB/GJF companions for inputs with PDB/Gaussian templates. | `--convert-files` |
| `--out-dir TEXT` | Output directory for all files. | `./result_opt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--args-yaml FILE` | Supply YAML overrides (sections `geom`, `calc`, `opt`, `lbfgs`, `rfo`). | _None_ |

## Outputs
```
out_dir/
├─ final_geometry.xyz          # Always written
├─ final_geometry.pdb          # Only when the input was a PDB and conversion is enabled
├─ final_geometry.gjf          # When a Gaussian template was detected and conversion is enabled
├─ optimization.trj            # Only if dumping is enabled
├─ optimization.pdb            # PDB conversion of the trajectory (PDB inputs, conversion enabled)
├─ optimization.gjf            # Multi-geometry Gaussian dump (Gaussian templates, conversion enabled)
└─ restart*.yml                # Optional restarts when opt.dump_restart is set
```
The console prints the resolved `geom`, `calc`, `opt`, `lbfgs`/`rfo` blocks plus cycle-by-cycle progress and total runtime.

## YAML sections (`--args-yaml`)
YAML values override CLI, which override the defaults below.

### `geom`
- `coord_type` (`"cart"`): Cartesian vs. `"dlc"` delocalized internal coordinates.
- `freeze_atoms` (`[]`): Base 0-based frozen indices; automatically merged with CLI link detection.

### `calc`
- UMA configuration (`model`, `task_name`, device selection, neighbor radii, Hessian format, etc.).
- `charge`/`spin` mirror the CLI options; defaults come from `.gjf` when present.

### `opt`
Shared optimizer controls used by both LBFGS and RFO:
- `thresh` presets (Gaussian-like or Baker rule). Presets translate to the force/step thresholds documented in `pdb2reaction/opt.py`.
- `max_cycles`, `print_every` (`100`), `min_step_norm` (`1e-8`), `assert_min_step`, convergence toggles (`rms_force`, etc.), RMSD-based `converge_to_geom_rms_thresh`, `overachieve_factor`, `check_eigval_structure`, `line_search`.
- Dumping/bookkeeping fields (`dump`, `dump_restart`, `prefix`, `out_dir`).

### `lbfgs`
Extends `opt` with L-BFGS specifics: `keep_last`, `beta`, `gamma_mult`, `max_step`, `control_step`, `double_damp`, and the optional regularisation parameters `mu_reg`/`max_mu_reg_adaptions`.

### `rfo`
Extends `opt` with RFOptimizer fields: trust-region sizing (`trust_radius`, `trust_min`, `trust_max`, `trust_update`), `max_energy_incr`, Hessian management (`hessian_update`, `hessian_init`, `hessian_recalc`, `hessian_recalc_adapt`, `small_eigval_thresh`), micro-iteration controls (`alpha0`, `max_micro_cycles`, `rfo_overlaps`), DIIS helpers (`gdiis`, `gediis`, thresholds, `gdiis_test_direction`), and `adapt_step_func`.

### Example YAML
```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
  model: uma-s-1p1
  task_name: omol
  device: auto
  max_neigh: null
  radius: null
  r_edges: false
  out_hess_torch: true
  freeze_atoms: null
  hessian_calc_mode: Analytical
  return_partial_hessian: true
opt:
  thresh: baker
  max_cycles: 10000
  print_every: 100
  min_step_norm: 1.0e-08
  assert_min_step: true
  rms_force: null
  rms_force_only: false
  max_force_only: false
  force_only: false
  converge_to_geom_rms_thresh: 0.05
  overachieve_factor: 0.0
  check_eigval_structure: false
  line_search: true
  dump: false
  dump_restart: false
  prefix: ""
  out_dir: ./result_opt/
lbfgs:
  thresh: baker
  max_cycles: 10000
  print_every: 100
  min_step_norm: 1.0e-08
  assert_min_step: true
  rms_force: null
  rms_force_only: false
  max_force_only: false
  force_only: false
  converge_to_geom_rms_thresh: 0.05
  overachieve_factor: 0.0
  check_eigval_structure: false
  line_search: true
  dump: false
  dump_restart: false
  prefix: ""
  out_dir: ./result_opt/
  keep_last: 7
  beta: 1.0
  gamma_mult: false
  max_step: 0.3
  control_step: true
  double_damp: true
  mu_reg: null
  max_mu_reg_adaptions: 10
rfo:
  thresh: baker
  max_cycles: 10000
  print_every: 10
  min_step_norm: 1.0e-08
  assert_min_step: true
  rms_force: null
  rms_force_only: false
  max_force_only: false
  force_only: false
  converge_to_geom_rms_thresh: 0.05
  overachieve_factor: 0.0
  check_eigval_structure: false
  line_search: true
  dump: false
  dump_restart: false
  prefix: ""
  out_dir: ./result_opt/
  trust_radius: 0.1
  trust_update: true
  trust_min: 0.0
  trust_max: 0.1
  max_energy_incr: null
  hessian_update: bfgs
  hessian_init: calc
  hessian_recalc: 200
  hessian_recalc_adapt: null
  small_eigval_thresh: 1.0e-08
  alpha0: 1.0
  max_micro_cycles: 50
  rfo_overlaps: false
  gediis: false
  gdiis: true
  gdiis_thresh: 0.0025
  gediis_thresh: 0.01
  gdiis_test_direction: true
  adapt_step_func: true
```
