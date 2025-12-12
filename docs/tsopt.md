# `tsopt` subcommand

## Overview
`pdb2reaction tsopt` optimizes transition states using two complementary workflows:

- **light** mode: Hessian Dimer search with periodic exact-Hessian refreshes, a
  memory-conscious flatten loop to remove surplus imaginary modes, and PHVA-aware
  Hessian updates for the active degrees of freedom.
- **heavy** mode: RS-I-RFO optimizer with configurable trust-region safeguards.

Both modes use the UMA calculator for energies/gradients/Hessians, inherit `geom`/`calc`/`opt`
settings from YAML, and always write the final imaginary mode in `.trj`. When
`--convert-files` is enabled (default), PDB inputs mirror trajectories into `.pdb`
companions and Gaussian templates receive multi-geometry `.gjf` exports.

## Usage
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-m 2S+1] \
                    [--opt-mode light|heavy] \
                    [--freeze-links {True|False}] [--max-cycles N] [--thresh PRESET] \
                    [--dump {True|False}] [--out-dir DIR] [--args-yaml FILE] \
                    [--hessian-calc-mode Analytical|FiniteDifference] \
                    [--convert-files/--no-convert-files]
```

### Examples
```bash
# Recommended baseline: specify charge/multiplicity and pick the light workflow
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode light --out-dir ./result_tsopt/

# Light mode with YAML overrides, finite-difference Hessian, and freeze-links handling
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links True \
    --opt-mode light --max-cycles 10000 --dump False \
    --out-dir ./result_tsopt/ --args-yaml ./args.yaml \
    --hessian-calc-mode FiniteDifference

# Heavy mode (RS-I-RFO) driven entirely by YAML
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode heavy \
    --args-yaml ./args.yaml --out-dir ./result_tsopt/
```

## Workflow
- **Charge/spin resolution**: when the input is `.gjf`, charge and multiplicity inherit the
  template values; otherwise `-q/--charge` is required and multiplicity defaults to `1`.
  Override them explicitly to ensure UMA runs on the intended state.
- **Geometry loading & freeze-links**: structures are read via
  `pysisyphus.helpers.geom_loader`. On PDB inputs, `--freeze-links True` finds link hydrogens
  and freezes their parent atoms. The merged set is echoed, stored in `geom.freeze_atoms`, and
  forwarded to UMA's `calc.freeze_atoms`.
- **UMA Hessians**: `--hessian-calc-mode` toggles between analytical and finite-difference
  evaluations; both honor active (PHVA) subspaces. UMA may return only the active block when
  frozen atoms are present.
- **Light mode details**:
  - The Hessian Dimer stage periodically refreshes the dimer direction by evaluating an exact
    Hessian (active subspace, TR-projected) and prefers `torch.lobpcg` for the lowest
    eigenpair when `root == 0` (falling back to `torch.linalg.eigh`).
  - A flatten loop updates the stored active Hessian via Bofill (SR1/MS ↔ PSB blend) using
    displacements Δx and gradient differences Δg. Each loop estimates imaginary modes, flattens
    once, refreshes the dimer direction, runs a dimer+LBFGS micro-segment, and performs a final
    Bofill update. Once only one imaginary mode remains, a final exact Hessian is computed for
    frequency analysis.
  - If `root != 0`, that root seeds only the initial dimer direction; subsequent refreshes
    follow the most negative mode (`root = 0`).
- **Heavy mode (RS-I-RFO)**: runs the RS-I-RFO optimizer with optional Hessian reference files,
  R+S splitting safeguards, and micro-cycle controls defined in the `rsirfo` YAML section.
- **Mode export & conversion**: the converged imaginary mode is always written to `vib/final_imag_mode_*.trj`
  and mirrored to `.pdb`/`.gjf` when conversion is enabled. When the input was PDB, the optimization
  trajectory and final geometry are also converted to PDB via the input template; Gaussian templates
  receive multi-geometry `.gjf` companions.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Required unless the input is a `.gjf` template with charge metadata. | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | PDB-only. Freeze parents of link hydrogens (merged into `geom.freeze_atoms`). | `True` |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | Light/Heavy aliases listed above. | `light` |
| `--dump BOOL` | Explicit `True`/`False`. Dump trajectories. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ (uses YAML/default of `FiniteDifference`) |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB or Gaussian inputs. | `--convert-files` |
| `--args-yaml FILE` | YAML overrides (`geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo`). | _None_ |

## Outputs (& directory layout)
```
out_dir/ (default: ./result_tsopt/)
├─ final_geometry.xyz            # Always written
├─ final_geometry.pdb            # When the input was PDB (conversion enabled)
├─ final_geometry.gjf            # When the input was Gaussian (conversion enabled)
├─ optimization_all.trj          # Light-mode dump when --dump is True
├─ optimization_all.pdb          # Light-mode companion for PDB inputs (conversion enabled)
├─ optimization.trj              # Heavy-mode trajectory
├─ optimization.pdb              # Heavy-mode PDB companion when conversion is enabled
├─ vib/
│  ├─ final_imag_mode_±XXXX.Xcm-1.trj
│  └─ final_imag_mode_±XXXX.Xcm-1.pdb/.gjf
└─ .dimer_mode.dat               # Light-mode orientation seed
```

## Notes
- `--opt-mode` aliases map exactly to the workflows described above; pick one for the intended
  algorithm rather than adjusting YAML keys manually (default: `light`).
- Imaginary-mode detection defaults to ~5 cm⁻¹ (configurable via
  `hessian_dimer.neg_freq_thresh_cm`). The selected `root` determines which imaginary mode is
  exported when multiple remain.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after YAML merging, mirroring the
  behavior of other subcommands.
- PHVA translation/rotation projection mirrors the implementation in `freq`, reducing GPU
  memory consumption while preserving correct eigenvectors in the active space.

## YAML configuration (`--args-yaml`)
Provide a mapping; YAML values override CLI. Shared sections reuse
[`opt`](opt.md#yaml-configuration-args-yaml). Keep the full block below intact if it already
matches your workflow—adjust only the values you need to change.

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
  out_dir: ./result_tsopt/   # output directory
hessian_dimer:
  thresh_loose: gau_loose    # loose convergence preset
  thresh: gau                # main convergence preset
  update_interval_hessian: 1000   # Hessian rebuild cadence
  neg_freq_thresh_cm: 5.0    # negative frequency threshold (cm^-1)
  flatten_amp_ang: 0.1       # flattening amplitude (Å)
  flatten_max_iter: 0        # flattening iteration cap (default: Disabled)
  flatten_sep_cutoff: 2.0    # minimum distance between representative atoms (Å)
  flatten_k: 10              # representative atoms sampled per mode
  mem: 100000                # memory limit for solver
  use_lobpcg: true           # enable LOBPCG eigen solver
  device: auto               # device selection for eigensolver
  root: 0                    # targeted TS root index
  dimer:
    length: 0.0189           # dimer separation (Bohr)
    rotation_max_cycles: 15  # max rotation iterations
    rotation_method: fourier # rotation optimizer method
    rotation_thresh: 0.0001  # rotation convergence threshold
    rotation_tol: 1          # rotation tolerance factor
    rotation_max_element: 0.001   # max rotation matrix element
    rotation_interpolate: true    # interpolate rotation steps
    rotation_disable: false   # disable rotations entirely
    rotation_disable_pos_curv: true   # disable when positive curvature detected
    rotation_remove_trans: true   # remove translational components
    trans_force_f_perp: true  # project forces perpendicular to translation
    bonds: null               # bond list for constraints
    N_hessian: null           # Hessian size override
    bias_rotation: false      # bias rotational search
    bias_translation: false   # bias translational search
    bias_gaussian_dot: 0.1    # Gaussian bias dot product
    seed: null                # RNG seed for rotations
    write_orientations: true  # write rotation orientations
    forward_hessian: true     # propagate Hessian forward
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
    out_dir: ./result_opt/     # output directory
    keep_last: 7               # history size for LBFGS buffers
    beta: 1.0                  # initial damping beta
    gamma_mult: false          # multiplicative gamma update toggle
    max_step: 0.3              # maximum step length
    control_step: true         # control step length adaptively
    double_damp: true          # double damping safeguard
    mu_reg: null               # regularization strength
    max_mu_reg_adaptions: 10   # cap on mu adaptations
rsirfo:
  thresh: gau                # RS-IRFO convergence preset
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
  out_dir: ./result_opt/     # output directory
  roots: [0]                 # target root indices
  hessian_ref: null          # reference Hessian
  rx_modes: null             # reaction-mode definitions for projection
  prim_coord: null           # primary coordinates to monitor
  rx_coords: null            # reaction coordinates to monitor
  hessian_update: bofill     # Hessian update scheme override
  hessian_recalc_reset: true # reset recalc counter after exact Hessian
  max_micro_cycles: 50       # micro-iterations per macro cycle
  augment_bonds: false       # augment reaction path based on bond analysis
  min_line_search: true      # enforce minimum line-search step
  max_line_search: true      # enforce maximum line-search step
  assert_neg_eigval: false   # require a negative eigenvalue at convergence
```

