# `tsopt`

## Overview

> **Summary:** Optimize a transition-state *candidate* using Dimer (`--opt-mode light`), RS‑I‑RFO (`--opt-mode heavy`, default), or hybrid (`--opt-mode hybrid`: Dimer then RS‑I‑RFO flatten stage). A validated TS should show **exactly one** imaginary frequency; always confirm the mode/connectivity with freq/IRC.

### At a glance
- **Input:** A TS guess (HEI from `path-opt`/`path-search`, or your own structure) in any `geom_loader`-supported format.
- **Modes:** `heavy` = RS‑I‑RFO (default, generally more robust). `light` = Hessian Guided Dimer (often cheaper per step). `hybrid` = Dimer to convergence, then RS‑I‑RFO flatten loop only.
- **Quality control:** The optimized structure is still a *candidate* until [freq](freq.md) and [irc](irc.md) confirm the expected mode and connectivity.
- **Optional cleanup:** `--flatten` (default enabled) controls surplus-imaginary-mode cleanup.
- **Output conversion:** With `--convert-files` (default), PDB inputs can be mirrored to `.pdb` (when `--dump`), and Gaussian templates write a `.gjf` for the final geometry.

### Choosing `--opt-mode`
- Use **`--opt-mode heavy` (RS‑I‑RFO)** when you want the default, conservative optimizer and you can afford Hessian work.
- Use **`--opt-mode light` (Dimer)** when you want a lighter-weight search, or when you plan to iterate quickly from several TS guesses.
- Use **`--opt-mode hybrid`** when you want Dimer search behavior but heavy-style RS‑I‑RFO flattening.

For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB/GJF output conversion. If you need a TS guess first, run [path-opt](path_opt.md) (two structures) or [path-search](path_search.md) (two or more structures) and then validate/optimize the HEI with `tsopt` → `freq` → `irc`.

## Minimal example

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## Output checklist

- `result_tsopt/final_geometry.pdb` (or `final_geometry.xyz`)
- `result_tsopt/vib/final_imag_mode_*_trj.xyz`
- `result_tsopt/vib/final_imag_mode_*.pdb` (for PDB inputs)

## Common examples

1. Use light mode with analytical Hessian when VRAM is sufficient.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt_light
```

2. Keep the optimization trajectory for inspection.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. Run heavy mode with YAML overrides.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode heavy --config tsopt.yaml --out-dir ./result_tsopt_heavy
```

4. Run hybrid mode (Dimer + RS-I-RFO flatten stage).

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hybrid --flatten --out-dir ./result_tsopt_hybrid
```

## Usage
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [--opt-mode light|heavy|hybrid] [--flatten/--no-flatten] \
 [--freeze-links/--no-freeze-links] [--max-cycles N] [--thresh PRESET] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### Examples
```bash
# Recommended baseline: specify charge/multiplicity and pick the light workflow
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode light --out-dir ./result_tsopt/

# Light mode with YAML overrides, finite-difference Hessian, and freeze-links handling
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links \
 --opt-mode light --max-cycles 10000 --no-dump \
 --out-dir ./result_tsopt/ --config ./args.yaml \
 --hessian-calc-mode FiniteDifference

# Heavy mode (RS-I-RFO) driven entirely by YAML
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode heavy \
 --config ./args.yaml --out-dir ./result_tsopt/
```

## Workflow
- **Charge/spin resolution**: when the input is `.gjf`, charge and multiplicity inherit the
 template values. If `-q` is omitted but `--ligand-charge` is provided, the structure is
 treated as an enzyme–substrate complex and `extract.py`’s charge summary derives the total
 charge for PDB inputs (or XYZ/GJF with `--ref-pdb`); explicit `-q` still overrides. Otherwise
 `-q/--charge` is required and multiplicity defaults to `1`. Override them explicitly to ensure
 UMA runs on the intended state.
- **Geometry loading & freeze-links**: structures are read via
 `pysisyphus.helpers.geom_loader`. On PDB inputs, `--freeze-links` finds link hydrogens
 and freezes their parent atoms. The merged set is echoed, stored in `geom.freeze_atoms`, and
 forwarded to UMA's `calc.freeze_atoms`.
- **UMA Hessians**: `--hessian-calc-mode` toggles between analytical and finite-difference
 evaluations; both honor active (PHVA) subspaces. UMA may return only the active block when
 frozen atoms are present.
 When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.
- **Light mode details**:
 - The Hessian Guided Dimer stage periodically refreshes the dimer direction by evaluating an exact
 Hessian (active subspace, TR-projected) and prefers `torch.lobpcg` for the lowest
 eigenpair when `root == 0` (falling back to `torch.linalg.eigh`).
 - When enabled (`--flatten`), the flatten loop updates the stored active Hessian via
 Bofill (SR1/MS ↔ PSB blend; toggle via `hessian_dimer.flatten_loop_bofill`) using
 displacements Δx and gradient differences Δg. Each loop estimates imaginary modes, flattens
 once, refreshes the dimer direction, runs a dimer+LBFGS micro-segment, and (optionally)
 performs a Bofill update. Once only one imaginary mode remains, a final exact Hessian is
 computed for frequency analysis.
 - If `root != 0`, that root seeds only the initial dimer direction; subsequent refreshes
 follow the most negative mode (`root = 0`).
- **Heavy mode (RS-I-RFO)**: runs the RS-I-RFO optimizer with optional Hessian reference files,
 R+S splitting safeguards, and micro-cycle controls defined in the `rsirfo` YAML section.
 When `--flatten` is enabled and more than one imaginary mode remains after
 convergence, the workflow flattens extra modes and reruns RS-I-RFO until only one
 imaginary mode remains or the flatten iteration cap is reached.
- **Hybrid mode**: runs Dimer first (like light mode) and then runs only the heavy-style RS-I-RFO flatten stage.
- **Mode export & conversion**: the converged imaginary mode is always written to `vib/final_imag_mode_*_trj.xyz`
 and mirrored to `.pdb` when the input was PDB and conversion is enabled. The optimization
 trajectory and final geometry are also converted to PDB via the input template when `--dump`;
 Gaussian templates receive a `.gjf` companion for the final geometry only.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Required unless a `.gjf` template or `--ligand-charge` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `--ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | PDB-only. Freeze parents of link hydrogens (merged into `geom.freeze_atoms`). | `True` |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | Choose optimizer: `light` (Dimer), `heavy` (RSIRFO), or `hybrid` (Dimer then RSIRFO flatten stage). | `heavy` |
| `--dump/--no-dump` | Dump trajectories. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--flatten/--no-flatten` | Enable the extra-imaginary-mode flattening loop (`False` forces `flatten_max_iter=0`). Applies to light (dimer loop), heavy (post-RSIRFO), and hybrid (post-dimer RSIRFO stage). | `True` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB or Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--dry-run/--no-dry-run` | Validate inputs/config and print the execution plan without running TS optimization. | `False` |

## Outputs
```
out_dir/ (default:./result_tsopt/)
├─ final_geometry.xyz # Always written
├─ final_geometry.pdb # When the input was PDB (conversion enabled)
├─ final_geometry.gjf # When the input was Gaussian (conversion enabled)
├─ optimization_all_trj.xyz # Light-mode dump when --dump is True
├─ optimization_all.pdb # Light-mode companion for PDB inputs (conversion enabled, --dump)
├─ optimization_trj.xyz # Heavy-mode trajectory when --dump is True
├─ optimization.pdb # Heavy-mode PDB companion when conversion is enabled and --dump is True
├─ vib/
│ ├─ final_imag_mode_±XXXX.Xcm-1_trj.xyz
│ └─ final_imag_mode_±XXXX.Xcm-1.pdb
└─.dimer_mode.dat # Light-mode orientation seed
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.
- Imaginary-mode detection defaults to ~5 cm⁻¹ (configurable via
 `hessian_dimer.neg_freq_thresh_cm`). The selected `root` determines which imaginary mode is
 exported when multiple remain.
- Use `--opt-mode` to choose the algorithm workflow directly (`heavy` by default), instead of
 manually editing YAML mode mappings.
- PHVA translation/rotation projection follows the same implementation as `freq`, while reducing
 memory usage and preserving correct active-space eigenvectors.
- Config merge precedence is `defaults < config < explicit CLI < override`.

Shared sections reuse
[YAML Reference](yaml_reference.md). Keep the full block below intact if it already
matches your workflow—adjust only the values you need to change.

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
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
 return_partial_hessian: false # full Hessian (avoids shape mismatches)
opt:
 thresh: baker # convergence preset (Gaussian/Baker-style)
 max_cycles: 10000 # optimizer cycle cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum norm for step acceptance
 assert_min_step: true # stop if steps fall below threshold
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # geom RMS threshold when converging to ref
 overachieve_factor: 0.0 # factor to tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_tsopt/ # output directory
hessian_dimer:
 thresh_loose: gau_loose # loose convergence preset
 thresh: baker # main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # negative frequency threshold (cm^-1)
 flatten_amp_ang: 0.1 # flattening amplitude (Å)
 flatten_max_iter: 50 # flattening iteration cap (disabled when --no-flatten)
 flatten_sep_cutoff: 0.0 # minimum distance between representative atoms (Å)
 flatten_k: 10 # representative atoms sampled per mode
 flatten_loop_bofill: false # Bofill update for flatten displacements
 mem: 100000 # memory limit for solver
 device: auto # device selection for eigensolver
 root: 0 # targeted TS root index
 dimer:
 length: 0.0189 # dimer separation (Bohr)
 rotation_max_cycles: 15 # max rotation iterations
 rotation_method: fourier # rotation optimizer method
 rotation_thresh: 0.0001 # rotation convergence threshold
 rotation_tol: 1 # rotation tolerance factor
 rotation_max_element: 0.001 # max rotation matrix element
 rotation_interpolate: true # interpolate rotation steps
 rotation_disable: false # disable rotations entirely
 rotation_disable_pos_curv: true # disable when positive curvature detected
 rotation_remove_trans: true # remove translational components
 trans_force_f_perp: true # project forces perpendicular to translation
 bonds: null # bond list for constraints
 N_hessian: null # Hessian size override
 bias_rotation: false # bias rotational search
 bias_translation: false # bias translational search
 bias_gaussian_dot: 0.1 # Gaussian bias dot product
 seed: null # RNG seed for rotations
 write_orientations: true # write rotation orientations
 forward_hessian: true # propagate Hessian forward
 lbfgs:
 thresh: baker # LBFGS convergence preset
 max_cycles: 10000 # iteration limit
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_tsopt/ # output directory
 keep_last: 7 # history size for LBFGS buffers
 beta: 1.0 # initial damping beta
 gamma_mult: false # multiplicative gamma update toggle
 max_step: 0.3 # maximum step length
 control_step: true # control step length adaptively
 double_damp: true # double damping safeguard
 mu_reg: null # regularization strength
 max_mu_reg_adaptions: 10 # cap on mu adaptations
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 10000 # iteration cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir:./result_tsopt/ # output directory
 roots: [0] # target root indices
 hessian_ref: null # reference Hessian
 rx_modes: null # reaction-mode definitions for projection
 prim_coord: null # primary coordinates to monitor
 rx_coords: null # reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme override
 hessian_recalc_reset: true # reset recalc counter after exact Hessian
 max_micro_cycles: 50 # micro-iterations per macro cycle
 augment_bonds: false # augment reaction path based on bond analysis
 min_line_search: true # enforce minimum line-search step
 max_line_search: true # enforce maximum line-search step
 assert_neg_eigval: false # require a negative eigenvalue at convergence
```

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [path-search](path_search.md) — MEP search that identifies TS candidates (HEI)
- [irc](irc.md) — Trace the reaction path from an optimized TS
- [freq](freq.md) — Confirm a single imaginary frequency (expected for a validated TS)
- [all](all.md) — End-to-end workflow that chains extraction → MEP → tsopt → IRC → freq
- [YAML Reference](yaml_reference.md) — Full `hessian_dimer` (Hessian Guided Dimer) and `rsirfo` configuration options
- [Glossary](glossary.md) — Definitions of TS, Dimer, RS-I-RFO, Hessian
