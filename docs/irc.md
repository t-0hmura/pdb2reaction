# `irc` subcommand

## Purpose
Runs Intrinsic Reaction Coordinate (IRC) integrations with the EulerPC predictor–corrector method using UMA, writing forward/backward trajectories and optional HDF5 dumps.

## Usage
```bash
pdb2reaction irc -i INPUT -q CHARGE [--spin 2S+1]
                 [--max-cycles N] [--step-size Δs] [--root k]
                 [--forward BOOL] [--backward BOOL]
                 [--freeze-links/--no-freeze-links]
                 [--out-dir DIR]
                 [--hessian-calc-mode Analytical|FiniteDifference]
                 [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--max-cycles INT` | Maximum IRC steps (overrides `irc.max_cycles`). | _None_ (use YAML/default `125`) |
| `--step-size FLOAT` | Step length in mass-weighted coordinates (overrides `irc.step_length`). | _None_ (default `0.10`) |
| `--root INT` | Imaginary-mode index for the initial displacement (`irc.root`). | _None_ (default `0`) |
| `--forward BOOL` | Run forward branch (`irc.forward`). Explicit `True` or `False`. | _None_ (default `True`) |
| `--backward BOOL` | Run backward branch (`irc.backward`). Explicit `True` or `False`. | _None_ (default `True`) |
| `--freeze-links / --no-freeze-links` | For PDB inputs, freeze link-hydrogen parents (merged with `geom.freeze_atoms`). | `--freeze-links` |
| `--out-dir TEXT` | Output directory (`irc.out_dir`). | `./result_irc/` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode. | _None_ |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## YAML configuration (`--args-yaml`)
Provide a mapping; CLI overrides YAML. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml) for geometry/calculator keys.

### Shared sections
- `geom`: same keys as [`opt`](opt.md#section-geom). `--freeze-links` augments `geom.freeze_atoms`.
- `calc`: same keys as [`opt`](opt.md#section-calc). `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after merging.

### Section `irc`
EulerPC / IRC controls (defaults in parentheses).

- `step_length` (`0.10`): Step length in mass-weighted coordinates (overridden by `--step-size`).
- `max_cycles` (`125`): Maximum IRC steps (overridden by `--max-cycles`).
- `downhill` (`False`): Follow downhill potential energy.
- `forward` (`True`), `backward` (`True`): Integrate forward/backward branches.
- `root` (`0`): Imaginary-mode index for initial displacement.
- `hessian_init` (`"calc"`): Source for initial Hessian.
- `displ` (`"energy"`), `displ_energy` (`1.0e-3`), `displ_length` (`0.10`): Displacement control.
- Convergence thresholds: `rms_grad_thresh` (`1.0e-3`), `hard_rms_grad_thresh` (`null`), `energy_thresh` (`1.0e-6`), `imag_below` (`0.0`).
- Termination flags: `force_inflection` (`True`), `check_bonds` (`False`).
- Output: `out_dir` (`"./result_irc/"`), `prefix` (`""`), `dump_fn` (`"irc_data.h5"`), `dump_every` (`5`).
- EulerPC specifics: `hessian_update` (`"bofill"`), `hessian_recalc` (`null`), `max_pred_steps` (`500`), `loose_cycles` (`3`), `corr_func` (`"mbs"`).

## Outputs
- `<out-dir>/irc_data.h5` (written every `dump_every` steps).
- `<out-dir>/<prefix>finished_irc.trj` plus optional forward/backward `.trj` and `.pdb` files when the input was PDB.
- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus timing.

## Notes
- CLI boolean options `--forward` and `--backward` require explicit `True` or `False` (e.g., `--forward True`).
- When the input is PDB, trajectory files are converted to PDB using the input topology.
- IRC calculations reuse the UMA calculator across steps; large `step_length` values may destabilise the integration.
