# `opt` subcommand

## Purpose
Performs single-structure geometry optimizations with pysisyphus using either the LBFGS ("light") or RFOptimizer ("heavy") algorithms backed by the UMA machine-learning calculator.

## Usage
```bash
pdb2reaction opt -i INPUT -q CHARGE [--spin 2S+1] [--opt-mode light|lbfgs|heavy|rfo]
                 [--freeze-links/--no-freeze-links] [--dump/--no-dump]
                 [--out-dir DIR] [--max-cycles N] [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader` (`.pdb`, `.xyz`, `.trj`, …). | Required |
| `-q, --charge INT` | Total charge passed to UMA. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-links / --no-freeze-links` | When the input is PDB, detect link hydrogens and freeze their parent atoms (merged with `geom.freeze_atoms`). | `--freeze-links` |
| `--max-cycles INT` | Maximum optimization cycles (`opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Select optimizer: `light`/`lbfgs` → LBFGS, `heavy`/`rfo` → RFO. | `light` |
| `--dump / --no-dump` | Emit `optimization.trj`. | `--no-dump` |
| `--out-dir TEXT` | Output directory. | `./result_opt/` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## YAML configuration (`--args-yaml`)
Pass a YAML mapping; CLI values override YAML, which override the defaults below.

### Section `geom`
Geometry loader options (see also `pdb2reaction.uma_pysis.GEOM_KW_DEFAULT`).

- `coord_type` (`"cart"`): Coordinate representation for `geom_loader` (`"cart"` or `"dlc"`).
- `freeze_atoms` (`[]`): 0-based indices to freeze; merged with link parents when `--freeze-links`.

### Section `calc`
UMA calculator options (see also `pdb2reaction/uma_pysis.py`).

- `charge` / `spin`: Overridden by `-q/-s`.
- `model` (`"uma-s-1p1"`) and `task_name` (`"omol"`): UMA checkpoint and dataset tag.
- `device` (`"auto"`), `max_neigh`, `radius`, `r_edges`: Graph/device controls.
- `out_hess_torch` (`True`): Request a Torch Hessian (GPU-friendly).
- `freeze_atoms`: Auto-filled from `geom.freeze_atoms`.
- `hessian_calc_mode` (`"Analytical"`) and `return_partial_hessian` (`True`): UMA Hessian controls.

### Section `opt`
Shared optimizer controls.

- `thresh` (`"gau"`): Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`).
- `max_cycles` (`10000`), `print_every` (`1`), `dump` (`False`), `dump_restart` (`False`), `prefix` (`""`), `out_dir` (`"./result_opt/"`).
- Safeguards: `min_step_norm` (`1e-8`), `assert_min_step` (`True`).
- Convergence flags: `rms_force`, `rms_force_only`, `max_force_only`, `force_only`.
- Extras: `converge_to_geom_rms_thresh`, `overachieve_factor`, `check_eigval_structure`, `line_search`.

### Section `lbfgs`
Specific to the LBFGS optimizer (`--opt-mode light|lbfgs`).

- `keep_last` (`7`): Memory depth.
- `beta` (`1.0`), `gamma_mult` (`False`): Initial scaling.
- Step control: `max_step` (`0.30`), `control_step` (`True`).
- Safeguards: `double_damp` (`True`).
- Regularization: `mu_reg` (`None`), `max_mu_reg_adaptions` (`10`).

### Section `rfo`
Specific to RFOptimizer (`--opt-mode heavy|rfo`).

- Trust region: `trust_radius` (`0.30`), `trust_min` (`0.01`), `trust_max` (`0.30`), `trust_update` (`True`).
- Energy guard: `max_energy_incr` (`None`).
- Hessian management: `hessian_update` (`"bfgs"`), `hessian_init` (`"calc"`), `hessian_recalc` (`100`), `hessian_recalc_adapt` (`2.0`).
- Numerics: `small_eigval_thresh` (`1e-8`), `line_search` (`True`), `alpha0` (`1.0`), `max_micro_cycles` (`25`), `rfo_overlaps` (`False`).
- DIIS controls: `gdiis` (`True`), `gediis` (`False`), `gdiis_thresh` (`2.5e-3`), `gediis_thresh` (`1.0e-2`), `gdiis_test_direction` (`True`).
- `adapt_step_func` (`False`).

## Outputs
- `final_geometry.xyz` (+ `.pdb` when the input was PDB).
- Optional `optimization.trj` (+ `.pdb` when dumping and PDB input).
- Optional restart YAML files when `opt.dump_restart` is set.
- Console blocks summarising resolved `geom`, `calc`, `opt`, and LBFGS/RFO settings.

## Notes
- Always provide the chemically correct charge and multiplicity.
- `--freeze-links` is PDB-only and merges link parents into `geom.freeze_atoms`.
- CLI precedence: CLI > YAML > built-in defaults.
- Exit codes: `0` success, `2` zero-step error, `3` optimizer failure, `130` interrupt, `1` unexpected error.
