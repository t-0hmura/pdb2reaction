# `irc`

## Overview

> **Summary:** Runs EulerPC (Euler Predictor-Corrector)-based intrinsic reaction coordinate (IRC) integration from a transition state toward reactants and products. By default both forward and backward branches are computed. Setting `--hessian-calc-mode Analytical` is strongly recommended when VRAM permits.

### At a glance
- **Input:** A TS structure (ideally already optimized and validated).
- **Key knobs:** `--step-size` (mass-weighted step length) and `--max-cycles` (number of steps).
- **Hard overrides:** IRC forces `geom.coord_type = cart` after merge (even if YAML sets it). `calc.return_partial_hessian` is forced to `true` (partial Hessian with active-DOF processing in pysisyphus).

`pdb2reaction irc` runs EulerPC-based IRC integrations with an MLIP backend (UMA by default). The CLI is intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible.

For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB output conversion. A typical workflow is `tsopt` (which includes an imaginary-frequency check; confirm **one** imaginary frequency) → `irc`.

## Minimal example

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## Output checklist

- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`

## Common examples

1. Run only the forward branch.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward \
 --out-dir ./result_irc_forward
```

2. Increase step size and use analytical Hessians.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.20 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

3. Keep both branches and raise the step limit.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 150 \
 --out-dir ./result_irc_long
```

## Usage
```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] \
 [--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] [-m 2S+1]
 [--max-cycles N] [--step-size Δs] [--root k]
 [--freeze-links/--no-freeze-links]
 [--out-dir DIR]
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
 [--hessian-calc-mode Analytical|FiniteDifference]
 [--show-config] [--dry-run]
```

### Examples
```bash
# Forward-only branch, finite-difference Hessian, larger step size
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward \
 --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB input so finished and directional trajectories are also exported as PDB
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## Workflow
1. **Input preparation** – Any format supported by `geom_loader` is accepted. When a reference PDB is available (input is `.pdb` or `--ref-pdb` is supplied), EulerPC trajectories are converted to PDB using that topology, and `--freeze-links` augments `geom.freeze_atoms` by freezing parents of link hydrogens for PDB inputs. Note: `geom.coord_type` is forced to `cart` (Cartesian) regardless of YAML/CLI settings, and `calc.return_partial_hessian` is forced to `true` (partial Hessian with active-DOF processing).
2. **EulerPC integration** – The EulerPC predictor-corrector integrator traces the IRC path from the transition state. Forward and/or backward branches are run according to `--forward`/`--backward` flags. Each step uses a mass-weighted steepest-descent predictor followed by a corrector step.
3. **Trajectory output** – Finished, forward, and backward IRC trajectories are written as XYZ files. When a reference PDB is available, PDB companions are also generated (`--convert-files`).

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge; used unless YAML sets `calc.charge`. Required unless a `.gjf` template or `--ligand-charge` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Explicit `-q` still overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `--ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); used unless YAML sets `calc.spin`. | `.gjf` template value or `1` |
| `--max-cycles INT` | Maximum IRC steps; used unless YAML sets `irc.max_cycles`. | `125` |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; used unless YAML sets `irc.step_length`. | `0.10` |
| `--root INT` | Index of the imaginary vibrational mode for the initial displacement; used unless YAML sets `irc.root`. | `0` |
| `--forward/--no-forward` | Run forward branch (`irc.forward`), used unless YAML sets `irc.forward`. | `True` |
| `--backward/--no-backward` | Run backward branch (`irc.backward`), used unless YAML sets `irc.backward`. | `True` |
| `--freeze-links/--no-freeze-links` | For PDB inputs, freeze link-H parents (merged with `geom.freeze_atoms`). See [extract](extract.md) for link-hydrogen details. | `True` |
| `--out-dir TEXT` | Output directory (`irc.out_dir`), used unless YAML sets `irc.out_dir`. | `./result_irc/` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB companions when a reference PDB is available. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`calc.hessian_calc_mode`), used unless YAML sets `calc.hessian_calc_mode`. | `FiniteDifference` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. | `False` |

## Outputs
```
out_dir/ (default:./result_irc/)
├─ <prefix>finished_irc_trj.xyz # Complete IRC trajectory
├─ <prefix>forward_irc_trj.xyz # Present when the forward branch runs
└─ *.pdb # Trajectory companions when a reference PDB is available (conversion enabled)
```
- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus wall-clock timing.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The MLIP backend is reused throughout the IRC; aggressive `step_length` values can destabilize EulerPC.
- For Hessian evaluation modes, see [MLIP Calculator](uma_pysis.md#hessian-evaluation).
- When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see [Concepts: Link hydrogen](concepts.md#link-hydrogen-and-frozen-atoms)).

See [CLI Conventions: Configuration precedence](cli_conventions.md#configuration-precedence) for the full resolution order.
Shared sections reuse [YAML Reference](yaml_reference.md) for geometry/calculator keys: `--freeze-links` augments `geom.freeze_atoms` for PDB inputs, and `--hessian-calc-mode` plus CLI charge/spin values supplement the merged `calc` block. For `irc`, `geom.coord_type` is forced to `cart` and `calc.return_partial_hessian` is forced to `true` after YAML/CLI merging.

`irc` keys (defaults in parentheses):
- `step_length` (`0.10`), `max_cycles` (`125`): primary integration controls surfaced via `--step-size`/`--max-cycles`.
- `hessian_init` (`"calc"`), `hessian_update` (`"bofill"`), `hessian_recalc` (`null`): Hessian initialization/update cadence.
- `displ`, `displ_energy`, `displ_length`: displacement construction; keep defaults unless debugging.
- Convergence thresholds: `rms_grad_thresh` (`1.0e-3`), `hard_rms_grad_thresh` (`null`), `energy_thresh` (`1.0e-6`), `imag_below` (`0.0`).
- Output & diagnostics: `force_inflection` (`True`), `check_bonds` (`False`), `out_dir` (`"./result_irc/"`), `prefix` (`""`), `max_pred_steps` (`500`), `loose_cycles` (`3`), `corr_func` (`"mbs"`).

```yaml
geom:
 coord_type: cart # forced to cart for irc (YAML value ignored)
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # UMA model tag
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true # forced true for irc (partial Hessian with active-DOF processing)
irc:
 step_length: 0.1 # integration step length
 max_cycles: 125 # maximum steps along IRC
 downhill: false # follow downhill direction only
 forward: true # propagate in forward direction
 backward: true # propagate in backward direction
 root: 0 # normal-mode root index
 hessian_init: calc # Hessian initialization source
 displ: energy # displacement construction method
 displ_energy: 0.001 # energy-based displacement scaling
 displ_length: 0.1 # length-based displacement fallback
 rms_grad_thresh: 0.001 # RMS gradient convergence threshold
 hard_rms_grad_thresh: null # hard RMS gradient stop
 energy_thresh: 0.000001 # energy change threshold
 imag_below: 0.0 # imaginary frequency cutoff
 force_inflection: true # enforce inflection detection
 check_bonds: false # check bonds during propagation
 out_dir: ./result_irc/ # output directory
 prefix: "" # filename prefix
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: null # Hessian rebuild cadence
 max_pred_steps: 500 # predictor-corrector max steps
 loose_cycles: 3 # loose cycles before tightening
 corr_func: mbs # correlation function choice
```

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [tsopt](tsopt.md) — Optimize the TS before running IRC
- [freq](freq.md) — Full vibrational analysis and thermochemistry (imaginary-frequency check is already included in `tsopt`)
- [opt](opt.md) — Optimize IRC endpoints to true minima
- [all](all.md) — End-to-end workflow that runs IRC after tsopt
- [YAML Reference](yaml_reference.md) — Full `irc` configuration options
- [Glossary](glossary.md) — Definition of IRC (Intrinsic Reaction Coordinate)
