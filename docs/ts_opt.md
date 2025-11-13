# `ts_opt` subcommand

## Purpose
Optimizes transition states using either the Hessian Dimer method ("light") or RS-I-RFO ("heavy"), leveraging UMA for energies, gradients, and Hessians, and exporting the final imaginary mode.

## Usage
```bash
pdb2reaction ts_opt -i INPUT -q CHARGE [--spin 2S+1]
                    [--freeze-links/--no-freeze-links]
                    [--max-cycles N] [--opt-mode light|heavy]
                    [--dump/--no-dump] [--out-dir DIR]
                    [--hessian-calc-mode Analytical|FiniteDifference]
                    [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-links / --no-freeze-links` | For PDB inputs, freeze link-hydrogen parents (propagated to UMA). | `--freeze-links` |
| `--max-cycles INT` | Maximum macro cycles (forwarded to `opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | `light` → Hessian Dimer, `heavy` → RS-I-RFO. | `light` |
| `--dump / --no-dump` | Dump optimization trajectories. | `--no-dump` |
| `--out-dir TEXT` | Output directory. | `./result_ts_opt/` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ (use YAML/default) |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## YAML configuration (`--args-yaml`)
Provide a mapping; CLI overrides YAML. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

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
- `<out-dir>/final_geometry.xyz` (+ `.pdb` when the input was PDB).
- `<out-dir>/optimization.trj` or `optimization_all.trj` (mode-dependent) and PDB conversions when applicable.
- `<out-dir>/vib/final_imag_mode_±XXXX.Xcm-1.(trj|pdb)` for the converged imaginary mode.
- Optional `.dimer_mode.dat` (light mode) and dump files when `--dump` is enabled.

## Notes
- Always provide accurate charge and multiplicity; they are propagated to UMA.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after YAML merging.
- In light mode, the flattened active-subspace Hessian is maintained and updated between LBFGS passes.
- Imaginary-mode analysis uses PHVA and translation/rotation projection consistent with the frequency module.
