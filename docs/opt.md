# `opt` subcommand

## Overview
`pdb2reaction opt` performs a single-structure geometry optimization with the pysisyphus LBFGS ('light') or RFOptimizer ('heavy') engines while UMA provides energies, gradients, and Hessians. Input structures can be `.pdb`, `.xyz`, `.trj`, or any format supported by `geom_loader`. Settings are applied in the order **built-in defaults → CLI overrides → `--args-yaml` overrides** (YAML has the highest precedence), making it easy to keep lightweight defaults while selectively overriding options. The optimizer preset now defaults to the LBFGS-based **`light`** mode.

When the starting structure is a PDB or Gaussian template, format-aware conversion mirrors the optimized structure into `.pdb` (PDB inputs) and `.gjf` (Gaussian templates) companions, controlled by `--convert-files {True|False}` (enabled by default). PDB-specific conveniences include:
- With `--freeze-links` (default `True`), parent atoms of link hydrogens are detected and merged into `geom.freeze_atoms` (0-based indices).
- Output conversion produces `final_geometry.pdb` (and `optimization.pdb` when dumping trajectories) using the input PDB as the topology reference.

A Gaussian `.gjf` template seeds the charge/spin defaults and enables automatic export of the optimized structure as `.gjf` when conversion is enabled.

## Usage
```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} -q CHARGE -m MULT \
                 [--opt-mode light|heavy] [--freeze-links BOOL] \
                 [--dist-freeze '[(i,j,target_A), ...]'] [--one-based {True|False}] \
                 [--bias-k K_eV_per_A2] [--dump BOOL] [--out-dir DIR] \
                 [--max-cycles N] [--thresh PRESET] [--args-yaml FILE] \
                 [--convert-files {True|False}]
```

## Workflow
- **Optimizers**: `--opt-mode light` (default) → L-BFGS; `--opt-mode heavy` → rational-function optimizer with trust-region control.
- **Restraints**: `--dist-freeze` consumes Python-literal tuples `(i, j, target_A)`; omitting the third element restrains the starting distance. `--bias-k` sets a global harmonic strength (eV·Å⁻²). Indices default to 1-based but can be flipped to 0-based with `--one-based False`.
- **Charge/spin resolution**: CLI `-q/-m` override `.gjf` template metadata, which in turn override the `calc` defaults. If `-q` is omitted but `--ligand-charge` is provided, the input is treated as an enzyme–substrate complex and `extract.py`’s charge summary derives the total charge; explicit `-q` still overrides. If no template or ligand-charge hint exists, the fallback is `0/1`. Always pass the physically correct values explicitly.
- **Freeze atoms**: CLI freeze-link logic is merged with YAML `geom.freeze_atoms`, then propagated to the UMA calculator (`calc.freeze_atoms`).
- **Dumping & conversion**: `--dump True` mirrors `opt.dump=True` and writes `optimization.trj`; when conversion is enabled, trajectories are mirrored to `.pdb` for PDB inputs. `opt.dump_restart` can emit restart YAML snapshots.
- **Exit codes**: `0` success, `2` zero step (step norm < `min_step_norm`), `3` optimizer failure, `130` keyboard interrupt, `1` unexpected error.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Required unless the input is a `.gjf` template that already encodes charge. Overrides `--ligand-charge` when both are set. | Required when not in template |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Falls back to `.gjf` template or `1`. | Template/`1` |
| `--dist-freeze TEXT` | Repeatable string parsed as Python literal describing `(i,j,target_A)` tuples for harmonic restraints. | _None_ |
| `--one-based {True|False}` | Interpret `--dist-freeze` indices as 1-based (default) or 0-based. | `True` |
| `--bias-k FLOAT` | Harmonic bias strength applied to every `--dist-freeze` tuple (eV·Å⁻²). | `10.0` |
| `--freeze-links BOOL` | Toggle link-hydrogen parent freezing (PDB inputs only). | `True` |
| `--max-cycles INT` | Hard limit on optimization iterations (`opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Choose optimizer: `light` (LBFGS) or `heavy` (RFO). | `light` |
| `--dump BOOL` | Emit trajectory dumps (`optimization.trj`). | `False` |
| `--convert-files {True|False}` | Enable or disable XYZ/TRJ → PDB companions for PDB inputs and XYZ → GJF companions for Gaussian templates. | `True` |
| `--out-dir TEXT` | Output directory for all files. | `./result_opt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--args-yaml FILE` | Supply YAML overrides (sections `geom`, `calc`, `opt`, `lbfgs`, `rfo`). | _None_ |

## Outputs
```
out_dir/
├─ final_geometry.xyz          # Always written
├─ final_geometry.pdb          # Only when the input was a PDB and conversion is enabled
├─ final_geometry.gjf          # When a Gaussian template was detected and conversion is enabled
├─ optimization.trj            # Only if dumping is enabled
├─ optimization.pdb            # PDB conversion of the trajectory (PDB inputs, conversion enabled)
└─ restart*.yml                # Optional restarts when opt.dump_restart is set
```
The console prints the resolved `geom`, `calc`, `opt`, `lbfgs`/`rfo` blocks plus cycle-by-cycle progress and total runtime.

## YAML sections (`--args-yaml`)
YAML values override CLI, which override the defaults below.

### `geom`
- `coord_type` (`'cart'`): Cartesian vs. `'dlc'` delocalized internal coordinates.
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
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  task_name: omol            # UMA task name
  device: auto               # UMA device selection
  max_neigh: null            # maximum neighbors for graph construction
  radius: null               # cutoff radius for neighbor search
  r_edges: false             # store radial edges
  out_hess_torch: true       # request torch-form Hessian
  freeze_atoms: null         # calculator-level frozen atoms
  hessian_calc_mode: FiniteDifference   # Hessian mode selection
  return_partial_hessian: false         # full Hessian (avoids shape mismatches)
opt:
  thresh: gau                # convergence preset (Gaussian/Baker-style)
  max_cycles: 10000          # optimizer cycle cap
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum norm for step acceptance
  assert_min_step: true      # stop if steps fall below threshold
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # geom RMS threshold when converging to ref
  overachieve_factor: 0.0    # factor to tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ''                 # filename prefix
  out_dir: ./result_opt/     # output directory
lbfgs:
  thresh: gau                # LBFGS convergence preset
  max_cycles: 10000          # iteration limit
  print_every: 100           # logging stride
  min_step_norm: 1.0e-08     # minimum accepted step norm
  assert_min_step: true      # assert when steps stagnate
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
  overachieve_factor: 0.0    # tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ''                 # filename prefix
  out_dir: ./result_opt/     # output directory
  keep_last: 7               # history size for LBFGS buffers
  beta: 1.0                  # initial damping beta
  gamma_mult: false          # multiplicative gamma update toggle
  max_step: 0.3              # maximum step length
  control_step: true         # control step length adaptively
  double_damp: true          # double damping safeguard
  mu_reg: null               # regularization strength
  max_mu_reg_adaptions: 10   # cap on mu adaptations
rfo:
  thresh: gau                # RFOptimizer convergence preset
  max_cycles: 10000          # iteration cap
  print_every: 100           # logging stride (matches shared opt defaults)
  min_step_norm: 1.0e-08     # minimum accepted step norm
  assert_min_step: true      # assert when steps stagnate
  rms_force: null            # explicit RMS force target
  rms_force_only: false      # rely only on RMS force convergence
  max_force_only: false      # rely only on max force convergence
  force_only: false          # skip displacement checks
  converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
  overachieve_factor: 0.0    # tighten thresholds
  check_eigval_structure: false   # validate Hessian eigenstructure
  line_search: true          # enable line search
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  prefix: ''                 # filename prefix
  out_dir: ./result_opt/     # output directory
  trust_radius: 0.1          # trust-region radius
  trust_update: true         # enable trust-region updates
  trust_min: 0.0             # minimum trust radius
  trust_max: 0.1             # maximum trust radius
  max_energy_incr: null      # allowed energy increase per step
  hessian_update: bfgs       # Hessian update scheme
  hessian_init: calc         # Hessian initialization source
  hessian_recalc: 200        # rebuild Hessian every N steps
  hessian_recalc_adapt: null # adaptive Hessian rebuild limit
  small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
  alpha0: 1.0                # initial micro step
  max_micro_cycles: 50       # micro-iteration limit
  rfo_overlaps: false        # enable RFO overlaps
  gediis: false              # enable GEDIIS
  gdiis: true                # enable GDIIS
  gdiis_thresh: 0.0025       # GDIIS acceptance threshold
  gediis_thresh: 0.01        # GEDIIS acceptance threshold
  gdiis_test_direction: true # test descent direction before DIIS
  adapt_step_func: true      # adaptive step scaling
```
