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
                 [--convert-files {True|False}]
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
2. **Configuration merge** – Defaults → CLI → YAML (`geom`, `calc`, `irc`). Charge/multiplicity inherit `.gjf` template metadata when available. If `-q` is omitted but `--ligand-charge` is provided, the structure is treated as an enzyme–substrate complex and `extract.py`’s charge summary derives the total charge; explicit `-q` still overrides. Otherwise charge defaults to `0` and multiplicity to `1`. Always set them explicitly to remain on the intended PES.
3. **IRC integration** – EulerPC integrates forward/backward branches according to `irc.forward/backward`, `irc.step_length`, `irc.root`, and the Hessian workflow configured through UMA (`calc.*`, `--hessian-calc-mode`). Hessians are updated with the configured scheme (`bofill` by default) and can be recalculated periodically.
4. **Outputs** – Trajectories (`finished`, `forward`, `backward`) are written as `.trj` and, for PDB inputs, mirrored to `.pdb`. Optional HDF5 dumps capture per-step frames when `dump_every` > 0.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge; overrides `calc.charge`. Required unless the input is a `.gjf` template with charge metadata. Overrides `--ligand-charge` when both are set. | Required when not in template |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); overrides `calc.spin`. | `.gjf` template value or `1` |
| `--max-cycles INT` | Maximum IRC steps; overrides `irc.max_cycles`. | _None_ (use YAML/default `125`) |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; overrides `irc.step_length`. | _None_ (default `0.10`) |
| `--root INT` | Imaginary-mode index for the initial displacement; overrides `irc.root`. | _None_ (default `0`) |
| `--forward BOOL` | Run forward branch (`irc.forward`). Requires explicit `True`/`False`. | _None_ (default `True`) |
| `--backward BOOL` | Run backward branch (`irc.backward`). Requires explicit `True`/`False`. | _None_ (default `True`) |
| `--freeze-links BOOL` | For PDB inputs, freeze link-H parents (merged with `geom.freeze_atoms`). | `True` |
| `--out-dir TEXT` | Output directory (`irc.out_dir`). | `./result_irc/` |
| `--convert-files {True|False}` | Toggle XYZ/TRJ → PDB companions for PDB inputs. | `True` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`calc.hessian_calc_mode`). | `FiniteDifference` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## Outputs
```
out_dir/ (default: ./result_irc/)
├─ irc_data.h5                # Written every `dump_every` steps when enabled
├─ <prefix>finished_irc.trj   # Complete IRC trajectory
├─ <prefix>forward_irc.trj    # Present when the forward branch runs
├─ <prefix>backward_irc.trj   # Present when the backward branch runs
└─ *.pdb                      # Trajectory companions for PDB inputs (when conversion is enabled)
```
- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus wall-clock timing.

## Notes
- CLI booleans (`--forward`, `--backward`) must be spelled out (`True`/`False`) to be merged into YAML when desired.
- UMA is reused throughout the IRC; aggressive `step_length` values can destabilise EulerPC.
- When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.
- Charge/spin inherit `.gjf` metadata when possible. If `-q` is omitted but `--ligand-charge` is provided, the input is treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge; explicit `-q` still overrides. Otherwise charge defaults to `0` and multiplicity to `1`. Override them explicitly for non-standard states.
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
irc:
  step_length: 0.1           # integration step length
  max_cycles: 125            # maximum steps along IRC
  downhill: false            # follow downhill direction only
  forward: true              # propagate in forward direction
  backward: true             # propagate in backward direction
  root: 0                    # normal-mode root index
  hessian_init: calc         # Hessian initialization source
  displ: energy              # displacement construction method
  displ_energy: 0.001        # energy-based displacement scaling
  displ_length: 0.1          # length-based displacement fallback
  rms_grad_thresh: 0.001     # RMS gradient convergence threshold
  hard_rms_grad_thresh: null # hard RMS gradient stop
  energy_thresh: 0.000001    # energy change threshold
  imag_below: 0.0            # imaginary frequency cutoff
  force_inflection: true     # enforce inflection detection
  check_bonds: false         # check bonds during propagation
  out_dir: ./result_irc/     # output directory
  prefix: ""                 # filename prefix
  dump_fn: irc_data.h5       # IRC data filename
  dump_every: 5              # dump stride
  hessian_update: bofill     # Hessian update scheme
  hessian_recalc: null       # Hessian rebuild cadence
  max_pred_steps: 500        # predictor-corrector max steps
  loose_cycles: 3            # loose cycles before tightening
  corr_func: mbs             # correlation function choice
```
