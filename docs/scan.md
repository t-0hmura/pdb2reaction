# `scan` subcommand

## Purpose
Performs staged bond-length scans with harmonic restraints and single-structure relaxations after each step, using UMA plus pysisyphus optimizers.

## Usage
```bash
pdb2reaction scan -i INPUT -q CHARGE --scan-lists "[(i,j,target), ...]" [...]
                  [--one-based/--zero-based] [--max-step-size ΔÅ] [--bias-k k]
                  [--relax-max-cycles N] [--opt-mode light|lbfgs|heavy|rfo]
                  [--freeze-links BOOL] [--dump BOOL]
                  [--out-dir DIR] [--preopt BOOL] [--endopt BOOL]
                  [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--scan-lists TEXT` | Repeatable Python-like list literal of `(i, j, targetÅ)` triples defining each scan stage. | Required |
| `--one-based / --zero-based` | Interpret `(i, j)` indices as 1-based (default) or 0-based. | `--one-based` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per integration step (Å). | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` (eV·Å⁻²). Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Maximum optimizer cycles per step (overrides `opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Relaxation optimizer (`light|lbfgs` or `heavy|rfo`). | `light` |
| `--freeze-links BOOL` | Explicit `True`/`False`. Freeze link-hydrogen parents for PDB inputs. | `True` |
| `--dump BOOL` | Explicit `True`/`False`. Dump stage trajectories. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_scan/` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |
| `--preopt BOOL` | Explicit `True`/`False`. Pre-optimize the initial structure before scanning. | `True` |
| `--endopt BOOL` | Explicit `True`/`False`. Unbiased relaxation after each stage. | `True` |

### Shared sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: see [`opt`](opt.md#yaml-configuration-args-yaml) for all keys and defaults. `opt.dump` is internally forced to `False`; dumping is controlled by `--dump`.

### Section `bias`
Controls the harmonic restraint applied during each scan stage.

- `k` (`100`): Harmonic strength in eV·Å⁻².

### Section `bond`
Parameters for UMA-based bond-change detection (mirrors `path_search`).

- `device` (`"cuda"`): UMA device for bond analysis.
- `bond_factor` (`1.20`): Covalent-radius scale factor for bond cutoff.
- `margin_fraction` (`0.05`): Fractional tolerance when comparing bond lengths.
- `delta_fraction` (`0.05`): Minimum fractional change to classify a bond as forming/breaking.

## Outputs
- `<out-dir>/stage_XX/scan.trj` (and `.pdb` when the input was PDB and dumping is enabled).
- Final unbiased geometries after `--endopt` stored per stage.
- Console summaries of resolved `geom`, `calc`, `opt`, `bias`, `bond`, and optimizer blocks.
- Final human-readable stage summary printed to the terminal (no separate YAML file).

## Notes
- `--scan-lists` accepts multiple literals; each defines one stage, and tuple indices are normalized to 0-based internally.
- When `--one-based` is used (default), indices follow PDB conventions; invalid indices or non-positive target distances raise `click.BadParameter`.
- Pre- and end-of-stage optimizations share UMA calculator instances for efficiency.
- Stage dumping writes one trajectory per stage; `--dump` also triggers `.pdb` conversion for PDB inputs.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. CLI parameters override YAML. Shared sections reuse the definitions documented for [`opt`](opt.md#yaml-configuration-args-yaml).

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
  thresh: gau
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
  out_dir: ./result_scan/
lbfgs:
  thresh: gau
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
  out_dir: ./result_scan/
  keep_last: 7
  beta: 1.0
  gamma_mult: false
  max_step: 0.3
  control_step: true
  double_damp: true
  mu_reg: null
  max_mu_reg_adaptions: 10
rfo:
  thresh: gau
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
  out_dir: ./result_scan/
  trust_radius: 0.3
  trust_update: true
  trust_min: 0.01
  trust_max: 0.3
  max_energy_incr: null
  hessian_update: bfgs
  hessian_init: calc
  hessian_recalc: 100
  hessian_recalc_adapt: 2.0
  small_eigval_thresh: 1.0e-08
  alpha0: 1.0
  max_micro_cycles: 25
  rfo_overlaps: false
  gediis: false
  gdiis: true
  gdiis_thresh: 0.0025
  gediis_thresh: 0.01
  gdiis_test_direction: true
  adapt_step_func: false
bias:
  k: 100
bond:
  device: cuda
  bond_factor: 1.2
  margin_fraction: 0.05
  delta_fraction: 0.05
```
