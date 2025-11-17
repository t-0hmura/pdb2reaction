# `tsopt` subcommand

## Purpose
Optimizes transition states using either the Hessian Dimer method ("light") or RS-I-RFO ("heavy"), leveraging UMA for energies, gradients, and Hessians, and exporting the final imaginary mode.

## Usage
```bash
pdb2reaction tsopt -i INPUT -q CHARGE [--spin 2S+1]
                    [--freeze-links BOOL] [--thresh PRESET]
                    [--max-cycles N]
                    [--opt-mode light|lbfgs|dimer|simple|simpledimer|hessian_dimer|
                               heavy|rfo|rsirfo|rs-i-rfo]
                    [--dump BOOL] [--outdir DIR]
                    [--hessian-calc-mode Analytical|FiniteDifference]
                    [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | `.gjf` template value or `0` |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | Explicit `True`/`False`. For PDB inputs, freeze link-hydrogen parents (propagated to UMA). | `True` |
| `--max-cycles INT` | Maximum macro cycles (forwarded to `opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Hessian Dimer aliases: `light`/`lbfgs`/`dimer`/`simple`/`simpledimer`/`hessian_dimer`. RS-I-RFO aliases: `heavy`/`rfo`/`rsirfo`/`rs-i-rfo`. | `light` |
| `--dump BOOL` | Explicit `True`/`False`. Dump optimization trajectories. | `False` |
| `--outdir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Override the convergence preset for both workflows (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (use YAML/default) |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ (use YAML/default) |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

### Shared sections
- `geom`, `calc`, `opt`: same keys as [`opt`](opt.md#yaml-configuration-args-yaml). `--freeze-links` augments `geom.freeze_atoms` and pushes the list into `calc.freeze_atoms`.

### Section `hessian_dimer`
Controls the light-mode TS workflow. Top-level keys:

- `thresh_loose` (`"gau_loose"`): First-pass convergence preset.
- `thresh` (`"gau"`): Main convergence preset.
- `update_interval_hessian` (`1000`): LBFGS cycles between exact Hessian refreshes.
- `neg_freq_thresh_cm` (`5.0`): Threshold (cm⁻¹) below which a frequency counts as imaginary.
- `flatten_amp_ang` (`0.10` Å), `flatten_max_iter` (`20`): Flattening displacement and iteration cap.
- `mem` (`100000`): Scratch memory forwarded to UMA during Hessian evaluations.
- `use_lobpcg` (`True`): Attempt LOBPCG for the most negative eigenpair.
- `device` (`"auto"`): Torch device for Hessian/mode algebra.
- `root` (`0`): Imaginary mode index to follow.
- `dimer`: Nested dictionary forwarded to the pysisyphus `Dimer` calculator. Key highlights (defaults in parentheses):
  - `length` (`0.0189` Bohr), rotation controls (`rotation_max_cycles: 15`, `rotation_method: "fourier"`, `rotation_thresh: 1e-4`, etc.).
  - Bias flags: `bias_rotation` (`False`), `bias_translation` (`False`), `bias_gaussian_dot` (`0.1`).
  - Optional initial guesses: `bonds`, `N_hessian`.
- `lbfgs`: Nested dictionary for the inner LBFGS pass; defaults align with [`opt`](opt.md#section-lbfgs) but tuned for TS search (e.g., `keep_last`, `max_step`, `double_damp`).

### Section `rsirfo`
Controls the heavy-mode RS-I-RFO optimizer. Defaults derive from [`opt`](opt.md#section-rfo) with TS-specific additions:

- `roots` (`[0]`): Mode indices to follow uphill.
- `hessian_ref` (`null`): Optional reference Hessian path.
- `rx_modes`, `prim_coord`, `rx_coords`: Optional monitoring definitions.
- `hessian_update` (`"bofill"`), `hessian_recalc_reset` (`True`), `max_micro_cycles` (`50`).
- Step safeguards: `augment_bonds` (`False`), `min_line_search` (`True`), `max_line_search` (`True`).
- `assert_neg_eigval` (`False`): Require a negative eigenvalue at convergence.
- Other trust-region parameters inherit from [`opt`](opt.md#section-rfo).

## Outputs
- `<outdir>/final_geometry.xyz` (+ `.pdb` when the input was PDB).
- `<outdir>/optimization.trj` or `optimization_all.trj` (mode-dependent) and PDB conversions when applicable.
- `<outdir>/vib/final_imag_mode_±XXXX.Xcm-1.(trj|pdb)` for the converged imaginary mode.
- Optional `.dimer_mode.dat` (light mode) and dump files when `--dump` is enabled.

## Notes
- Charge/spin inherit `.gjf` template metadata when present; otherwise the CLI defaults to `0`/`1`. Override them explicitly to
  ensure the UMA calculations use the intended state.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after YAML merging.
- In light mode, the flattened active-subspace Hessian is maintained and updated between LBFGS passes.
- Imaginary-mode analysis uses PHVA and translation/rotation projection consistent with the frequency module.
- The `--opt-mode` aliases shown above mirror the CLI behavior: light-side names select the Hessian Dimer workflow, while heavy-side names run RS-I-RFO.

## YAML configuration (`--args-yaml`)
Provide a mapping; CLI overrides YAML. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

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
  outdir: ./result_tsopt/
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
    outdir: ./result_opt/
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
  outdir: ./result_opt/
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