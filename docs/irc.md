# `irc`

## Overview

> **Summary:** Trace the intrinsic reaction coordinate (IRC) from a transition state toward reactant and product. By default it runs both forward and backward directions. When VRAM permits, `--hessian-calc-mode Analytical` is recommended.

### At a glance
- **Input:** A TS structure (ideally already optimized and validated).
- **Branches:** Runs both directions by default (`--forward True`, `--backward True`).
- **Key knobs:** `--step-size` (mass-weighted step length) and `--max-cycles` (number of steps).
- **Hard overrides:** IRC forces `geom.coord_type = cart` and `calc.return_partial_hessian = false` after merge (even if YAML sets them).
- **Outputs:** `finished_irc.trj` plus `forward_irc.trj`/`backward_irc.trj` (and `.pdb` companions when a PDB reference exists and conversion is enabled).

`pdb2reaction irc` runs EulerPC-based IRC integrations with UMA. The CLI is intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible.

For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB output conversion. A typical workflow is `tsopt` → `freq` (confirm **one** imaginary mode) → `irc`.

## Usage
```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] \
                 [--workers N] [--workers-per-node N] [-m 2S+1]
                 [--max-cycles N] [--step-size Δs] [--root k]
                 [--forward True|False] [--backward True|False]
                 [--freeze-links True|False]
                 [--out-dir DIR]
                 [--convert-files {True\|False}] [--ref-pdb FILE]
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
1. **Input preparation** – Any format supported by `geom_loader` is accepted. When a reference PDB is available (input is `.pdb` or `--ref-pdb` is supplied), EulerPC trajectories are converted to PDB using that topology, and `--freeze-links` augments `geom.freeze_atoms` by freezing parents of link hydrogens for PDB inputs.
2. **Configuration merge** – Defaults → CLI → YAML (`geom`, `calc`, `irc`). Charge/multiplicity inherit `.gjf` template metadata when available. For non-`.gjf` inputs, `-q/--charge` is required unless `--ligand-charge` is provided (supported for PDB inputs or XYZ/GJF with `--ref-pdb`); explicit `-q` still overrides. Multiplicity defaults to `1` when omitted. Always set them explicitly to remain on the intended PES. **IRC always forces Cartesian coordinates (`geom.coord_type = cart`) and `calc.return_partial_hessian = false`, regardless of YAML.**
3. **IRC integration** – EulerPC integrates forward/backward branches according to `irc.forward/backward`, `irc.step_length`, `irc.root`, and the Hessian workflow configured through UMA (`calc.*`, `--hessian-calc-mode`). Hessians are updated with the configured scheme (`bofill` by default) and can be recalculated periodically.
4. **Outputs** – Trajectories (`finished`, `forward`, `backward`) are written as `.trj` and, when a reference PDB is available, mirrored to `.pdb`.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge; used unless YAML sets `calc.charge`. Required unless a `.gjf` template or `--ligand-charge` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Explicit `-q` still overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); used unless YAML sets `calc.spin`. | `.gjf` template value or `1` |
| `--max-cycles INT` | Maximum IRC steps; used unless YAML sets `irc.max_cycles`. | `125` |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; used unless YAML sets `irc.step_length`. | `0.10` |
| `--root INT` | Imaginary-mode index for the initial displacement; used unless YAML sets `irc.root`. | `0` |
| `--forward {True\|False}` | Run forward branch (`irc.forward`), used unless YAML sets `irc.forward`. | `True` |
| `--backward {True\|False}` | Run backward branch (`irc.backward`), used unless YAML sets `irc.backward`. | `True` |
| `--freeze-links {True\|False}` | For PDB inputs, freeze link-H parents (merged with `geom.freeze_atoms`). | `True` |
| `--out-dir TEXT` | Output directory (`irc.out_dir`), used unless YAML sets `irc.out_dir`. | `./result_irc/` |
| `--convert-files {True\|False}` | Toggle XYZ/TRJ → PDB companions when a reference PDB is available. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`calc.hessian_calc_mode`), used unless YAML sets `calc.hessian_calc_mode`. | `FiniteDifference` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## Outputs
```
out_dir/ (default: ./result_irc/)
├─ <prefix>finished_irc.trj   # Complete IRC trajectory
├─ <prefix>forward_irc.trj    # Present when the forward branch runs
├─ <prefix>backward_irc.trj   # Present when the backward branch runs
└─ *.pdb                      # Trajectory companions when a reference PDB is available (conversion enabled)
```
- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus wall-clock timing.

## Notes
- CLI booleans (`--forward`, `--backward`) must be spelled out (`True`/`False`) to be merged into YAML when desired.
- UMA is reused throughout the IRC; aggressive `step_length` values can destabilize EulerPC.
- When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.
- Charge/spin inherit `.gjf` metadata when possible. For non-`.gjf` inputs, `-q/--charge` is required unless `--ligand-charge` is provided (supported for PDB inputs or XYZ/GJF with `--ref-pdb`); explicit `-q` still overrides. Multiplicity defaults to `1` when omitted. Override them explicitly for non-standard states.
- `--freeze-links` only applies to PDB inputs, keeping parent atoms of link hydrogens frozen during Hessian construction.

## YAML configuration (`--args-yaml`)
Provide a mapping; YAML overrides CLI. Shared sections reuse [YAML Reference](yaml-reference.md) for geometry/calculator keys: `--freeze-links` augments `geom.freeze_atoms` for PDB inputs, and `--hessian-calc-mode` plus CLI charge/spin values supplement the merged `calc` block. For `irc`, `geom.coord_type` is forced to `cart` and `calc.return_partial_hessian` is forced to `false` after YAML/CLI merging.

`irc` keys (defaults in parentheses):
- `step_length` (`0.10`), `max_cycles` (`125`): primary integration controls surfaced via `--step-size`/`--max-cycles`.
- `downhill` (`False`), `forward` (`True`), `backward` (`True`), `root` (`0`): branch selection and initial displacement index (`--forward`, `--backward`, `--root`).
- `hessian_init` (`"calc"`), `hessian_update` (`"bofill"`), `hessian_recalc` (`null`): Hessian initialization/update cadence.
- `displ`, `displ_energy`, `displ_length`: displacement construction; keep defaults unless debugging.
- Convergence thresholds: `rms_grad_thresh` (`1.0e-3`), `hard_rms_grad_thresh` (`null`), `energy_thresh` (`1.0e-6`), `imag_below` (`0.0`).
- Output & diagnostics: `force_inflection` (`True`), `check_bonds` (`False`), `out_dir` (`"./result_irc/"`), `prefix` (`""`), `max_pred_steps` (`500`), `loose_cycles` (`3`), `corr_func` (`"mbs"`).

```yaml
geom:
  coord_type: cart           # forced to cart for irc (YAML value ignored)
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
  return_partial_hessian: false         # forced false for irc (full Hessian)
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
  hessian_update: bofill     # Hessian update scheme
  hessian_recalc: null       # Hessian rebuild cadence
  max_pred_steps: 500        # predictor-corrector max steps
  loose_cycles: 3            # loose cycles before tightening
  corr_func: mbs             # correlation function choice
```

---

## See Also

- [tsopt](tsopt.md) — Optimize the TS before running IRC
- [freq](freq.md) — Verify the TS candidate has one imaginary frequency; analyze IRC endpoints
- [opt](opt.md) — Optimize IRC endpoints to true minima
- [all](all.md) — End-to-end workflow that runs IRC after tsopt
- [YAML Reference](yaml-reference.md) — Full `irc` configuration options
- [Glossary](glossary.md) — Definition of IRC (Intrinsic Reaction Coordinate)
