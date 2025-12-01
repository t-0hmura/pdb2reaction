# `irc` subcommand

## Overview
Run EulerPC-based Intrinsic Reaction Coordinate (IRC) integrations with UMA. The CLI is intentionally narrow: anything not listed below must be provided through YAML so that geometry handling, calculator settings, and low-level EulerPC knobs remain explicit and reproducible.

## Usage
```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-m 2S+1]
                 [--max-cycles N] [--step-size Δs] [--root k]
                 [--forward True|False] [--backward True|False]
                 [--freeze-links True|False]
                 [--out-dir DIR]
                 [--convert-files/--no-convert-files]
                 [--hessian-calc-mode Analytical|FiniteDifference]
                 [--args-yaml FILE]
```

### Examples
```bash
# Forward-only branch, finite-difference Hessian, larger step size
pdb2reaction irc -i ts.xyz -q -1 -m 2 --forward True --backward False \
                --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB input so finished and directional trajectories are also exported as PDB
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## Workflow
1. **Input preparation** – Any format supported by `geom_loader` is accepted. If the source is `.pdb`, EulerPC trajectories are automatically converted to PDB using the original topology, and `--freeze-links` augments `geom.freeze_atoms` by freezing parents of link hydrogens.
2. **Configuration merge** – Defaults → CLI → YAML (`geom`, `calc`, `irc`). Charge/multiplicity inherit `.gjf` template metadata when available; otherwise `-q/--charge` is required and multiplicity defaults to `1`. Always set them explicitly to remain on the intended PES.
3. **IRC integration** – EulerPC integrates forward/backward branches according to `irc.forward/backward`, `irc.step_length`, `irc.root`, and the Hessian workflow configured through UMA (`calc.*`, `--hessian-calc-mode`). Hessians are updated with the configured scheme (`bofill` by default) and can be recalculated periodically.
4. **Outputs** – Trajectories (`finished`, `forward`, `backward`) are written as `.trj` and, for PDB inputs, mirrored to `.pdb`. Optional HDF5 dumps capture per-step frames when `dump_every` > 0.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge; overrides `calc.charge`. Required unless the input is a `.gjf` template with charge metadata. | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); overrides `calc.spin`. | `.gjf` template value or `1` |
| `--max-cycles INT` | Maximum IRC steps; overrides `irc.max_cycles`. | _None_ (use YAML/default `125`) |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; overrides `irc.step_length`. | _None_ (default `0.10`) |
| `--root INT` | Imaginary-mode index for the initial displacement; overrides `irc.root`. | _None_ (default `0`) |
| `--forward BOOL` | Run forward branch (`irc.forward`). Requires explicit `True`/`False`. | _None_ (default `True`) |
| `--backward BOOL` | Run backward branch (`irc.backward`). Requires explicit `True`/`False`. | _None_ (default `True`) |
| `--freeze-links BOOL` | For PDB inputs, freeze link-H parents (merged with `geom.freeze_atoms`). | `True` |
| `--out-dir TEXT` | Output directory (`irc.out_dir`). | `./result_irc/` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `--convert-files` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`calc.hessian_calc_mode`). | _None_ |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## Outputs
- `<out-dir>/irc_data.h5` (written every `dump_every` steps).
- `<out-dir>/<prefix>finished_irc.trj` plus optional forward/backward `.trj` files.
- When the input is PDB or Gaussian, matching `.pdb`/`.gjf` trajectories are emitted using the original templates when `--convert-files` is enabled.
- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus wall-clock timing.

## Notes
- CLI booleans (`--forward`, `--backward`) must be spelled out (`True`/`False`) to be merged into YAML when desired.
- UMA is reused throughout the IRC; aggressive `step_length` values can destabilise EulerPC.
- Charge/spin inherit `.gjf` metadata when possible; override them explicitly for non-standard states.
- `--freeze-links` only applies to PDB inputs, keeping parent atoms of link hydrogens frozen during Hessian construction.

## YAML configuration (`--args-yaml`)
Provide a mapping; YAML overrides CLI. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml) for geometry/calculator keys: `--freeze-links` augments `geom.freeze_atoms` for PDB inputs, and `--hessian-calc-mode` plus CLI charge/spin values supplement the merged `calc` block.

`irc` keys (defaults in parentheses):
- `step_length` (`0.10`), `max_cycles` (`125`): primary integration controls surfaced via `--step-size`/`--max-cycles`.
- `downhill` (`False`), `forward` (`True`), `backward` (`True`), `root` (`0`): branch selection and initial displacement index (`--forward`, `--backward`, `--root`).
- `hessian_init` (`"calc"`), `hessian_update` (`"bofill"`), `hessian_recalc` (`null`): Hessian initialization/update cadence.
- `displ`, `displ_energy`, `displ_length`: displacement construction; keep defaults unless debugging.
- Convergence thresholds: `rms_grad_thresh` (`1.0e-3`), `hard_rms_grad_thresh` (`null`), `energy_thresh` (`1.0e-6`), `imag_below` (`0.0`).
- Output & diagnostics: `force_inflection` (`True`), `check_bonds` (`False`), `out_dir` (`"./result_irc/"`), `prefix` (`""`), `dump_fn` (`"irc_data.h5"`), `dump_every` (`5`), `max_pred_steps` (`500`), `loose_cycles` (`3`), `corr_func` (`"mbs"`).

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
irc:
  step_length: 0.1
  max_cycles: 125
  downhill: false
  forward: true
  backward: true
  root: 0
  hessian_init: calc
  displ: energy
  displ_energy: 0.001
  displ_length: 0.1
  rms_grad_thresh: 0.001
  hard_rms_grad_thresh: null
  energy_thresh: 0.000001
  imag_below: 0.0
  force_inflection: true
  check_bonds: false
  out_dir: ./result_irc/
  prefix: ""
  dump_fn: irc_data.h5
  dump_every: 5
  hessian_update: bofill
  hessian_recalc: null
  max_pred_steps: 500
  loose_cycles: 3
  corr_func: mbs
```
