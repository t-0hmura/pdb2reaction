# `scan` subcommand

## Overview
`scan` performs a staged, bond-length–driven scan using the UMA calculator and
harmonic restraints. Each tuple `(i, j, targetÅ)` defines a distance target. At
every integration step the temporary targets are updated, the restraint wells
are applied, and the entire structure is relaxed with RFOptimizer (`--opt-mode` heavy, default)
or LBFGS (`--opt-mode` light). After the biased walk, you can optionally
run unbiased pre-/post-optimizations to clean up the geometries that get written
to disk.

## Usage
```bash
pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} -q CHARGE [-m MULT] \
                  --scan-lists "[(i,j,targetÅ), ...]" [options]
                  [--convert-files/--no-convert-files]
```

### Examples
```bash
# Single-stage, minimal inputs (PDB)
pdb2reaction scan -i input.pdb -q 0 --scan-lists "[(12,45,1.35)]"

# Two stages, LBFGS relaxations, and trajectory dumping
pdb2reaction scan -i input.pdb -q 0 \
    --scan-lists "[(12,45,1.35)]" \
    --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
    --max-step-size 0.20 --dump True --out-dir ./result_scan/ --opt-mode light \
    --preopt True --endopt True
```

## Workflow
1. Load the structure through `geom_loader`, resolving charge/spin from the CLI
   overrides, the embedded Gaussian template (if present), or defaults.
2. Optionally run an unbiased preoptimization (`--preopt True`) before any
   biasing so the starting point is relaxed.
3. For each stage literal supplied via `--scan-lists`, parse and normalize the
   `(i, j)` indices (1-based by default). Compute the per-bond displacement
   `Δ = target − current` and split it into `N = ceil(max(|Δ|) / h)` steps using
   `h = --max-step-size`. Every bond receives its own `δ = Δ / N` increment.
4. March through all steps, updating the temporary targets, applying the
   harmonic wells `E = Σ ½ k (|ri − rj| − target)²`, and minimizing with UMA.
   Optimizer cycles are capped by `--relax-max-cycles` (overriding YAML).
5. After the last step of each stage, optionally run an unbiased relaxation
   (`--endopt True`) before reporting covalent bond changes and writing the
   `result.*` files.
6. Repeat for every stage; optional trajectories are dumped only when `--dump`
   is `True`.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge (CLI > template > 0). | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1 (CLI > template > 1). | `.gjf` template value or `1` |
| `--scan-lists TEXT` | Repeatable Python literal with `(i,j,targetÅ)` tuples. Each literal is one stage. | Required |
| `--one-based / --zero-based` | Interpret atom indices as 1- or 0-based. | `--one-based` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per step (Å). Controls the number of integration steps. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Cap on optimizer cycles during each biased step. Overrides `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS, `heavy` → RFOptimizer. | `heavy` |
| `--freeze-links BOOL` | When the input is PDB, freeze the parents of link hydrogens. | `True` |
| `--dump BOOL` | Dump concatenated biased trajectories (`scan.trj`/`scan.pdb`). | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `--convert-files` |
| `--out-dir TEXT` | Output directory root. | `./result_scan/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | Inherit YAML |
| `--args-yaml FILE` | YAML overrides for `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond`. | _None_ |
| `--preopt BOOL` | Run an unbiased optimization before scanning. | `True` |
| `--endopt BOOL` | Run an unbiased optimization after each stage. | `True` |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical keys to those documented in
  [`opt`](opt.md#yaml-configuration-args-yaml). `opt.dump` is internally forced
  to `False`; use `--dump` to control stage trajectories.

### Section `bias`
- `k` (`100`): Harmonic strength in eV·Å⁻².

### Section `bond`
UMA-based bond-change detection mirrored from `path_search`:
- `device` (`"cuda"`): UMA device for graph analysis.
- `bond_factor` (`1.20`): Covalent-radius scaling for cutoff.
- `margin_fraction` (`0.05`): Fractional tolerance for comparisons.
- `delta_fraction` (`0.05`): Minimum relative change to flag formation/breaking.

## Outputs
- `<out-dir>/preopt/` when `--preopt True`:
  - `result.xyz`, `result.gjf` (if the input provided a Gaussian template and conversion is enabled), and
    `result.pdb` (for PDB inputs when conversion is enabled).
- `<out-dir>/stage_XX/` for each stage:
  - `result.xyz` and optional `result.gjf`/`result.pdb` mirrors of the final
    structure (after `--endopt`, conversion enabled).
  - `scan.trj` when `--dump True` plus `scan.pdb`/`scan.gjf` companions for
    PDB/Gaussian inputs when conversion is enabled.
- Console summaries of the resolved `geom`, `calc`, `opt`, `bias`, `bond`, and
  optimizer blocks plus per-stage bond-change reports.

## Notes
- `--scan-lists` may be repeated; each literal becomes one stage. Tuples must
  have positive targets. Atom indices are normalized to 0-based internally.
- `--freeze-links` augments user `freeze_atoms` by adding parents of link-H
  atoms in PDB files so pockets stay rigid.
- UMA is the only supported calculator; energies are not re-queried for every
  biased frame to avoid redundant evaluations.
- Charge and spin inherit Gaussian template metadata when available; otherwise
  `-q/--charge` is required and spin defaults to `1`.
- Trajectories are written only when `--dump` is `True`; this also triggers PDB
  or Gaussian conversion when `--convert-files` is enabled.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. CLI parameters override YAML. Shared sections
reuse the definitions documented for [`opt`](opt.md#yaml-configuration-args-yaml).

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
