# `tsopt` subcommand

## Overview
`pdb2reaction tsopt` optimizes transition states using two complementary workflows:

- **light** mode: Hessian Dimer search with periodic exact-Hessian refreshes, a
  memory-conscious flatten loop to remove surplus imaginary modes, and PHVA-aware
  Hessian updates for the active degrees of freedom.
- **heavy** mode: RS-I-RFO optimizer with configurable trust-region safeguards.

Both modes use the UMA calculator for energies/gradients/Hessians, inherit `geom`/`calc`/`opt`
settings from YAML, and always write the final imaginary mode in `.trj` and `.pdb` formats.

## Usage
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-m 2S+1] \
                    [--opt-mode light|heavy] \
                    [--freeze-links {True|False}] [--max-cycles N] [--thresh PRESET] \
                    [--dump {True|False}] [--out-dir DIR] [--args-yaml FILE] \
                    [--hessian-calc-mode Analytical|FiniteDifference]
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
- **Mode export**: the converged imaginary mode is always written to `vib/final_imag_mode_*.trj`
  and `.pdb`. When the input was PDB, the optimization trajectory and final geometry are also
  converted to PDB via the input template.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Required unless the input is a `.gjf` template with charge metadata. | Required when not in template |
| `-m, --mult INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | PDB-only. Freeze parents of link hydrogens (merged into `geom.freeze_atoms`). | `True` |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | Light/Heavy aliases listed above. | `light` |
| `--dump BOOL` | Explicit `True`/`False`. Dump trajectories. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (use YAML/default) |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ (use YAML/default) |
| `--args-yaml FILE` | YAML overrides (`geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo`). | _None_ |

## Outputs (& directory layout)
```
out-dir/ (default: ./result_tsopt/)
  ├─ final_geometry.xyz                # Always written
  ├─ final_geometry.pdb                # When the input was PDB
  ├─ optimization_all.trj/.pdb         # Light-mode dump (requires --dump True)
  ├─ optimization.trj/.pdb             # Heavy-mode trajectory
  ├─ vib/
  │   ├─ final_imag_mode_±XXXX.Xcm-1.trj
  │   └─ final_imag_mode_±XXXX.Xcm-1.pdb
  └─ .dimer_mode.dat                   # Light-mode orientation seed
```

## Notes
- `--opt-mode` aliases map exactly to the workflows described above; pick one for the intended
  algorithm rather than adjusting YAML keys manually.
- Imaginary-mode detection defaults to ~5 cm⁻¹ (configurable via
  `hessian_dimer.neg_freq_thresh_cm`). The selected `root` determines which imaginary mode is
  exported when multiple remain.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after YAML merging, mirroring the
  behavior of other subcommands.
- PHVA translation/rotation projection mirrors the implementation in `freq`, reducing GPU
  memory consumption while preserving correct eigenvectors in the active space.

## YAML configuration (`--args-yaml`)
Provide a mapping; CLI values override YAML. Shared sections reuse
[`opt`](opt.md#yaml-configuration-args-yaml). Keep the full block below intact if it already
matches your workflow—adjust only the values you need to change.

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
  out_dir: ./result_tsopt/
hessian_dimer:
  thresh_loose: gau_loose
  thresh: gau
  update_interval_hessian: 1000
  neg_freq_thresh_cm: 5.0
  flatten_amp_ang: 0.1
  flatten_max_iter: 20
  mem: 100000
  use_lobpcg: true
  device: auto
  root: 0
  dimer:
    length: 0.0189
    rotation_max_cycles: 15
    rotation_method: fourier
    rotation_thresh: 0.0001
    rotation_tol: 1
    rotation_max_element: 0.001
    rotation_interpolate: true
    rotation_disable: false
    rotation_disable_pos_curv: true
    rotation_remove_trans: true
    trans_force_f_perp: true
    bonds: null
    N_hessian: null
    bias_rotation: false
    bias_translation: false
    bias_gaussian_dot: 0.1
    seed: null
    write_orientations: true
    forward_hessian: true
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
    out_dir: ./result_opt/
    keep_last: 7
    beta: 1.0
    gamma_mult: false
    max_step: 0.3
    control_step: true
    double_damp: true
    mu_reg: null
    max_mu_reg_adaptions: 10
rsirfo:
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
  out_dir: ./result_opt/
  roots: [0]
  hessian_ref: null
  rx_modes: null
  prim_coord: null
  rx_coords: null
  hessian_update: bofill
  hessian_recalc_reset: true
  max_micro_cycles: 50
  augment_bonds: false
  min_line_search: true
  max_line_search: true
  assert_neg_eigval: false
```
