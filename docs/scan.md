# `scan` subcommand

## Overview
`scan` performs a staged, bond-length–driven scan using the UMA calculator and
harmonic restraints. Each tuple `(i, j, targetÅ)` defines a distance target. At
every integration step the temporary targets are updated, the restraint wells
are applied, and the entire structure is relaxed with LBFGS (`--opt-mode` light, default)
or RFOptimizer (`--opt-mode` heavy). After the biased walk, you can optionally
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
   overrides, the embedded Gaussian template (if present), or defaults. If `-q`
   is omitted but `--ligand-charge` is provided, the input is treated as an
   enzyme–substrate complex and `extract.py`’s charge summary derives the total
   charge before any scans.
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
| `-q, --charge INT` | Total charge (CLI > template > 0). When omitted, charge can be inferred from `--ligand-charge`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge` supplies it |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1 (CLI > template > 1). | `.gjf` template value or `1` |
| `--scan-lists TEXT` | Repeatable Python literal with `(i,j,targetÅ)` tuples. Each literal is one stage. | Required |
| `--one-based / --zero-based` | Interpret atom indices as 1- or 0-based. | `--one-based` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per step (Å). Controls the number of integration steps. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². Overrides `bias.k`. | `100` |
| `--relax-max-cycles INT` | Cap on optimizer cycles during preopt, each biased step, and end-of-stage cleanups. Overrides `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `light` → LBFGS, `heavy` → RFOptimizer. | `light` |
| `--freeze-links BOOL` | When the input is PDB, freeze the parents of link hydrogens. | `True` |
| `--dump BOOL` | Dump concatenated biased trajectories (`scan.trj`/`scan.pdb`). | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `--convert-files` |
| `--out-dir TEXT` | Output directory root. | `./result_scan/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
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
```
out_dir/ (default: ./result_scan/)
├─ preopt/                   # Present when --preopt is True
│  ├─ result.xyz
│  ├─ result.pdb             # PDB companion for PDB inputs when conversion is enabled
│  └─ result.gjf             # When a Gaussian template exists and conversion is enabled
└─ stage_XX/                 # One folder per stage
    ├─ result.xyz
    ├─ result.pdb             # PDB mirror of the final structure (conversion enabled)
    ├─ result.gjf             # Gaussian mirror when templates exist and conversion is enabled
    ├─ scan.trj               # Written when --dump is True
    ├─ scan.pdb               # Trajectory companion for PDB inputs when conversion is enabled
    └─ scan.gjf               # Trajectory companion when a Gaussian template exists and conversion is enabled
```
- Console summaries of the resolved `geom`, `calc`, `opt`, `bias`, `bond`, and optimizer blocks plus per-stage bond-change reports.

## Notes
- `--scan-lists` may be repeated; each literal becomes one stage. Tuples must
  have positive targets. Atom indices are normalized to 0-based internally.
- `--freeze-links` augments user `freeze_atoms` by adding parents of link-H
  atoms in PDB files so pockets stay rigid.
- UMA is the only supported calculator; energies are not re-queried for every
  biased frame to avoid redundant evaluations.
- Charge and spin inherit Gaussian template metadata when available. If `-q` is
  omitted but `--ligand-charge` is provided, the full structure is treated as an
  enzyme–substrate complex and `extract.py`’s charge summary computes the total
  charge; explicit `-q` still overrides. Otherwise charge defaults to `0` and
  spin to `1`.
- Stage results (`result.xyz` plus optional PDB/GJF companions) are written
  regardless of `--dump`; trajectories are written only when `--dump` is `True`
  and converted to `scan.pdb`/`scan.gjf` when conversion is enabled.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. YAML parameters override CLI. Shared sections
reuse the definitions documented for [`opt`](opt.md#yaml-configuration-args-yaml).

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
  prefix: ""                 # filename prefix
  out_dir: ./result_scan/    # output directory
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
  prefix: ""                 # filename prefix
  out_dir: ./result_scan/    # output directory
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
  prefix: ""                 # filename prefix
  out_dir: ./result_scan/    # output directory
  trust_radius: 0.3          # trust-region radius
  trust_update: true         # enable trust-region updates
  trust_min: 0.01            # minimum trust radius
  trust_max: 0.3             # maximum trust radius
  max_energy_incr: null      # allowed energy increase per step
  hessian_update: bfgs       # Hessian update scheme
  hessian_init: calc         # Hessian initialization source
  hessian_recalc: 100        # rebuild Hessian every N steps
  hessian_recalc_adapt: 2.0  # adaptive Hessian rebuild factor
  small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
  alpha0: 1.0                # initial micro step
  max_micro_cycles: 25       # micro-iteration limit
  rfo_overlaps: false        # enable RFO overlaps
  gediis: false              # enable GEDIIS
  gdiis: true                # enable GDIIS
  gdiis_thresh: 0.0025       # GDIIS acceptance threshold
  gediis_thresh: 0.01        # GEDIIS acceptance threshold
  gdiis_test_direction: true # test descent direction before DIIS
  adapt_step_func: false     # adaptive step scaling toggle
bias:
  k: 100                    # harmonic bias strength (eV·Å⁻²)
bond:
  device: cuda               # UMA device for bond analysis
  bond_factor: 1.2           # covalent-radius scaling
  margin_fraction: 0.05      # tolerance margin for comparisons
  delta_fraction: 0.05       # minimum relative change to flag bonds
```
