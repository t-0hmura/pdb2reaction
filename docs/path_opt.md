# `path-opt` subcommand

## Purpose
Optimizes a minimum-energy path between two endpoints using the pysisyphus Growing String method with UMA providing energies, gradients, and Hessians.

## Usage
```bash
pdb2reaction path-opt -i REACTANT PRODUCT -q CHARGE [--spin 2S+1]
                      [--freeze-links BOOL]
                      [--max-nodes N] [--max-cycles N] [--climb BOOL]
                      [--dump BOOL] [--out-dir DIR] [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Two endpoint structures (reactant, product). | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-links BOOL` | Explicit `True`/`False`. For PDB inputs, freeze link-hydrogen parents. | `True` |
| `--max-nodes INT` | Internal nodes in the string (total images = `max_nodes + 2`). | `30` |
| `--max-cycles INT` | Maximum optimizer cycles. | `1000` |
| `--climb BOOL` | Explicit `True`/`False`. Enable climbing-image refinement. | `True` |
| `--dump BOOL` | Explicit `True`/`False`. Dump optimizer trajectories and restarts. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

### Shared sections
- `geom`, `calc`: same keys as [`opt`](opt.md#yaml-configuration-args-yaml). `--freeze-links` augments `geom.freeze_atoms` for PDB inputs.

### Section `gs`
Growing String controls (defaults shown in parentheses).

- `max_nodes` (`30`): Internal nodes (overridden by `--max-nodes`).
- `perp_thresh` (`5e-3`): Growth convergence threshold on perpendicular forces.
- `reparam_check` (`"rms"`), `reparam_every` (`1`), `reparam_every_full` (`1`): Reparametrisation cadence and metric.
- `param` (`"equi"`): Node spacing (`"equi"` or `"energy"`).
- `max_micro_cycles` (`10`): Micro-iterations per macro step.
- `reset_dlc` (`True`): Reset DLC coordinates when appropriate.
- `climb` (`True`), `climb_rms` (`5e-4`), `climb_lanczos` (`True`), `climb_lanczos_rms` (`5e-4`), `climb_fixed` (`False`): Climbing image controls (overridden by `--climb`).
- `scheduler` (`None`): Execution scheduler (normally left `None` for shared calculator use).

### Section `opt`
StringOptimizer controls (defaults in parentheses).

- `type` (`"string"`): Label for bookkeeping.
- `stop_in_when_full` (`1000`): Extra cycles allowed after full growth (overridden by `--max-cycles`).
- `align` (`False`): Internal alignment disabled (external Kabsch alignment is used).
- `scale_step` (`"global"`): Step scaling policy.
- `max_cycles` (`1000`): Macro-iteration cap (overridden by `--max-cycles`).
- `dump` (`False`), `dump_restart` (`False`): Trajectory/restart dumping (dump toggled by CLI).
- `reparam_thresh` (`1e-3`), `coord_diff_thresh` (`0.0`): Reparametrisation and pruning thresholds.
- `out_dir` (`"./result_path_opt/"`), `print_every` (`1`): Output location and logging cadence.

## Outputs
- `<out-dir>/final_geometries.trj` (+ `.pdb` when a PDB reference is available).
- `<out-dir>/gsm_hei.xyz` (+ `.pdb` for PDB inputs) for the highest-energy image.
- `<out-dir>/align_refine/` alignment diagnostics.
- Optional optimizer dumps/restarts when `--dump` or YAML toggles them.

## Notes
- Inputs are rigidly aligned via Kabsch before optimisation; if freeze atoms are present, only those atoms guide the fit.
- UMA calculator instances are shared across all images for efficiency.
- When the input is PDB, link-hydrogen parents are merged into `geom.freeze_atoms` before loading.
- Exit codes: `0` success, `3` optimization failure, `4`/`5` write failures, `130` interrupt, `1` unexpected error.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. CLI values override YAML. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

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
gs:
  max_nodes: 30
  perp_thresh: 0.005
  reparam_check: rms
  reparam_every: 1
  reparam_every_full: 1
  param: equi
  max_micro_cycles: 10
  reset_dlc: true
  climb: true
  climb_rms: 0.0005
  climb_lanczos: true
  climb_lanczos_rms: 0.0005
  climb_fixed: false
  scheduler: null
opt:
  type: string
  stop_in_when_full: 1000
  align: false
  scale_step: global
  max_cycles: 1000
  dump: false
  dump_restart: false
  reparam_thresh: 0.001
  coord_diff_thresh: 0.0
  out_dir: ./result_path_opt/
  print_every: 1
```