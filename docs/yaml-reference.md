# YAML Reference

This page lists, section by section, the keys and default values you can set in a YAML configuration file (`--config`). The list of sections, the order of precedence, and the mapping from CLI flags to YAML keys come first.

## Overview

| Section | Description | Used by |
|---------|-------------|---------|
| [`geom`](#geom) | Geometry and coordinate settings | opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search, dft, sp |
| [`calc`](#calc) | Machine-learning interatomic potential (MLIP) backend configuration | opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search, sp, dft (`charge` and `spin` only) |
| [`opt`](#opt) | Shared optimizer settings | opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGS optimizer settings | opt, scan, scan2d, scan3d, path-search, path-opt |
| [`rfo`](#rfo) | RFO optimizer settings | opt, scan, scan2d, scan3d, path-search, path-opt |
| [`gs`](#gs) | Growing String Method (GSM) string: nodes, climbing image, reparameterization | path-opt, path-search |
| [`dmf`](#dmf) | Direct Max Flux settings | path-opt, path-search |
| [`stopt`](#stopt) | StringOptimizer, which moves the GSM string; its `thresh` and `max_cycles` decide when GSM stops | path-opt, path-search |
| [`irc`](#irc-section) | IRC integration settings | irc |
| [`freq`](#freq-section) | Vibrational analysis settings | freq (`zero_cutoff_cm` also opt, tsopt) |
| [`thermo`](#thermo) | Thermochemistry settings | freq |
| [`dft`](#dft-section) | DFT calculation settings of the `dft` command | dft |
| [`bias`](#bias) | Harmonic bias settings | scan, scan2d, scan3d |
| [`bond`](#bond) | Bond-change detection settings | scan, path-search |
| [`search`](#search) | Recursive path search settings | path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian Guided Dimer TS optimization | tsopt |
| [`rsirfo`](#rsirfo) | RS-P-RFO / RS-I-RFO TS optimization | tsopt |
| `sp` | Single-point options (`hess`, default `false`; `hessian_calc_mode`; same as `--hess` and `--hessian-calc-mode`); see [sp](sp.md) | sp |

(yaml-configuration-precedence)=
## Configuration precedence

Settings are applied in the following order (later sources override earlier ones):

```
built-in defaults  <  --config (YAML)  <  CLI flags
```

1. **Built-in defaults** — the values shown by `pdb2reaction <subcmd> --help-advanced` and as `[default: …]` in the [Command Reference](reference/commands/index.md).
2. **`--config`** — a YAML file that overrides defaults (e.g., `--config my_settings.yaml`).
3. **CLI flags** — explicit command-line options (e.g., `-q -1`, `--thresh gau_loose`). Options left at their CLI default do not mask YAML values.

For example, if the YAML sets `charge: 0` but the CLI passes `-q -1`, the charge will be `-1`.

This precedence applies uniformly to `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `dft`, and `sp`.

To check the values a run uses, run it with {ref}`-v 3 <verbosity-levels>`, which prints each section as its name, a dashed underline, and the values in effect:

```text
opt
---
thresh: gau
max_cycles: 100000
…
```

`--show-config` (not on `scan`, `scan2d`, or `scan3d`) prints the loaded YAML file and its top-level keys.

A misspelled section name prints `[config] WARNING: YAML section(s) … are not recognized and were ignored.` A misspelled key inside a section is handled by section:

- `calc` prints `[backend] WARNING: … ignored calc setting(s) …` and goes on.
- `freq`, `thermo`, `bias`, `bond`, `search`, `sp`, the keys directly under `hessian_dimer`, and `geom` in the scan and path commands ignore it without a message. With `-v 3` it appears as an extra line in that section's block.
- Every other section stops with an error that names the key.

(common-cli-to-yaml-mapping)=
## Common CLI-to-YAML mapping

| CLI flag | YAML key | Section |
|----------|----------|---------|
| `-q` / `--charge` | `charge` | `calc` |
| `-m` / `--multiplicity` | `spin` | `calc` |
| `-b` / `--backend` | `backend` | `calc` |
| `--backend-model` | `model` | `calc` |
| `--solvent` | `solvent` | `calc` |
| _(YAML only)_ | `device` | `calc` |
| `--thresh` | `thresh` | `opt` |
| `--max-cycles` | `max_cycles` | Command-specific: `opt` for `opt`/`tsopt` and `irc` for `irc` |
| `--max-cycles-gsm` | `max_cycles` | `stopt` (also sets `stopt.stop_in_when_full`) |
| `--dmf-max-iterations` | `max_cycles` | `dmf` |
| `--gsm-param` | `param` | `gs` |
| `--dump` | `dump` | `opt` (opt, tsopt, scan), `stopt` (path-opt, path-search), `thermo` (freq) |
| `--step-size` (irc) | `step_length` | `irc` |
| `--opt-mode` | _(CLI only)_ | — |
| `--freeze-atoms` | `freeze_atoms` | `geom` |
| `--coord-type` | `coord_type` | `geom` |
| `--temperature` (freq, `all --freq-temperature`) | `temperature` | `thermo` |
| `--pressure` (freq, `all --freq-pressure`) | `pressure_atm` | `thermo` |
| `--dft-engine` | `engine` | `dft` |

```{note}
**Name mismatch — `--pressure` vs `pressure_atm`.** Both take atm (converted to Pa internally); only the YAML key names the unit.
```

### Default `--thresh` per subcommand

| Subcommand | Default `--thresh` |
|------------|-------------------|
| `opt` | `gau` |
| `tsopt` (Hessian Dimer) | `baker` |
| `tsopt` (RS-P-RFO / RS-I-RFO) | `baker` |
| `scan` | `gau` |
| `scan2d`, `scan3d` | `baker` |
| `path-opt`, `path-search` (single-structure optimizations) | `gau` |
| `path-opt`, `path-search` (GSM string: `--thresh-gsm`, `stopt.thresh`) | `gau_loose` |
| `all` (pre-opt, post-opt min) | `gau` |
| `all` (post-opt TS stage) | `baker` |

Accepted values: `gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`. Override per run with `--thresh <preset>` or under `opt.thresh` in YAML.

```{note}
**Subcommands without `--thresh`.** `irc`, `freq`, and `dft` do **not** expose `--thresh`:

- `irc` — convergence is governed by `irc.rms_grad_thresh`, `irc.energy_thresh`, and `irc.max_cycles`. The optimizer preset family does not apply because IRC follows a predictor–corrector integrator, not a force-based minimizer.
- `freq` — there is no optimization step, so no `--thresh`. Numerical accuracy is governed by `--hessian-calc-mode` and the underlying MLIP precision.
- `dft` — SCF convergence uses `dft.conv_tol` and `dft.max_cycle`, not the `gau`/`baker` preset family.
```

## Shared Sections

### `geom`

Geometry loading and coordinate handling.

```yaml
geom:
 coord_type: cart # "cart" (Cartesian), "redund" (redundant internals), "dlc" (delocalized internals), or "tric" (translation-rotation internals) for opt, tsopt, scan, scan2d, and scan3d; all, path-opt, and path-search accept cart and dlc only
 freeze_atoms: [] # 1-based atom indices to freeze; if `--freeze-links` is on (PDB/mmCIF input, or XYZ/GJF with `--ref-pdb`), the auto-detected cap-H parent indices are merged in
```

**Notes:**
- Frozen atoms have zeroed forces. With the default `return_partial_hessian: true`, the Hessian covers only the movable atoms; setting it false returns a full matrix with frozen rows and columns zeroed
- In Cartesian PHVA (partial Hessian vibrational analysis), only the rigid motions of the whole system that keep the frozen atoms in place are removed; see [freq](freq.md#rigid-modes-with-frozen-boundaries)
- For `irc`, `geom.coord_type` is always `cart`, whatever YAML or the CLI sets

---

### `calc`

Energy/force calculator configuration.

```yaml
calc:
 backend: uma           # uma, orb, mace, aimnet2, dft, or auto
 precision: auto # auto (uma/aimnet2 fp32, orb/mace fp64) | fp32 | fp64; aimnet2 accepts auto/fp32 and rejects fp64
 charge: 0 # Total charge; used only when written here (no built-in default); -q and -l override it
 spin: 1 # Spin multiplicity 2S+1 (overridden by CLI -m)
 model: uma-s-1p2 # UMA: uma-s-1p2 | uma-m-1p1. Without a model set, backend orb / mace / aimnet2 uses orb_v3_conservative_omol / MACE-OMOL-0 / aimnet2
 task_name: omol # Task tag recorded in UMA batches
 device: auto # Device: "cuda", "cpu", or "auto"
 max_neigh: null # Maximum neighbors for graph construction
 radius: null # Cutoff radius for neighbor search
 r_edges: false # Store radial edges
 workers: 1 # UMA inference workers
 workers_per_node: 1 # Workers per node for parallel predictor
 out_hess_torch: true # Return Hessian as torch.Tensor
 hessian_double: true # Assemble/return Hessian in float64
 # freeze_atoms: null # Inherited from geom.freeze_atoms; do not set directly
 hessian_calc_mode: FiniteDifference # Hessian mode: "Analytical" or "FiniteDifference"
 return_partial_hessian: true  # Return the Hessian of the movable atoms only
 print_timing: true # Print Hessian timing breakdown
 print_vram: true # Print CUDA VRAM usage during Hessian (UMA backend only)
 # xTB solvent correction (computationally expensive)
 solvent: none           # none, water, methanol, acetonitrile, dmso, thf, or toluene
 solvent_model: alpb     # xTB solvent model: "alpb" or "cpcmx"
 xtb_cmd: xtb            # xTB command plus optional arguments, e.g. "xtb --etemp 1000"
 xtb_acc: 0.2            # xTB accuracy parameter
 # Used only when backend: dft
 dft:
  func_basis: wb97m-v/def2-svp
  engine: gpu             # gpu (GPU4PySCF) | cpu (PySCF)
  lowmem: true             # direct JK without a persistent DF tensor
  density_fit: false       # default: the opposite of lowmem
  nprocs: auto             # PySCF/OpenMP threads from scheduler/affinity
  memory: auto             # host RAM limit, e.g. 64GB (not GPU VRAM)
  solvent: none
  solvent_model: smd      # pcm | smd
  save_scf_checkpoint: false
  checkpoint_path: null   # default when enabled (commands other than all): <out-dir>/_work/dft_scf/state.chk
  pyscf:
   mol: {}
   mf: {}
   grids: {}
   density_fit: {}
   with_df: {}
   with_solvent: {}
```

`backend: dft` computes energies, forces, and Hessians with DFT (PySCF or GPU4PySCF) through the `calc.dft` block above; see [Refine an MLIP TS with DFT](dft-backend.md). The top-level [`dft` section](#dft-section) configures the separate `dft` command (also run by `all --dft`). Both take the same SCF keys (`conv_tol`, `max_cycle`, `grid_level`, …); `save_scf_checkpoint` and `checkpoint_path` exist only in `calc.dft`. `hessian_calc_mode: Analytical` avoids the finite-displacement error, but its time and memory depend on the backend and system, so try it on your system first. For how `charge` and `spin` combine with the CLI and `.gjf` templates, see {ref}`Charge specification <charge-specification>`.

---

### `opt`

Shared single-structure optimizer controls used by both L-BFGS and RFO.

```yaml
opt:
 thresh: gau # Convergence preset: gau_loose, gau, gau_tight, gau_vtight, baker, never
 max_cycles: 100000 # Maximum optimizer iterations
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum step norm for acceptance
 assert_min_step: true # Stop if steps fall below threshold
 rms_force: null # Explicit RMS force target
 rms_force_only: false # Rely only on RMS force convergence
 max_force_only: false # Rely only on max force convergence
 force_only: false # Skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when converging to reference geometry
 overachieve_factor: 0.0 # 0.0 = off; >0: converge when forces < thresh/factor, ignoring step (not used by baker)
 check_eigval_structure: false # Validate Hessian eigenstructure
 line_search: true # Enable line search
 energy_plateau: false # Opt-in: stop as stalled when the energy range flattens (see note below)
 energy_plateau_thresh: 1.0e-04 # au (~0.06 kcal/mol); stalled-state threshold for the plateau check
 energy_plateau_window: 50 # Number of most recent steps inspected for the plateau check
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 prefix: "" # Filename prefix
 out_dir: ./result_opt/ # Output directory
```

**Energy plateau stop (off by default):**
`--stop-plateau` on `opt` / `tsopt` / `all` turns `energy_plateau` on, and
`--stop-plateau-thresh` / `--stop-plateau-window` set the two values below.
The optimizer then stops as `stalled`, not converged, when the energy range
(max − min) over the last `energy_plateau_window` steps falls below
`energy_plateau_thresh`. Use it when MLIP force noise keeps the force above the
`baker` threshold after the energy has flattened. When `tsopt` stops on a
plateau, it still computes the Hessian and reports n_imag; a TS optimization
that reaches `max_cycles` without converging skips the Hessian. GSM and DMF
path optimizations do not use this stop.

**Convergence Presets** (forces in Hartree/Bohr and steps in Bohr in Cartesian coordinates; Hartree/rad and rad for angles):

| Preset | Max Force | RMS Force | Max Step | RMS Step |
|--------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` requires all four columns plus `|delta E| < 1e-6` hartree against the
previous cycle, which is stricter than the Baker criterion as stated by Bakken and
Helgaker (*J. Chem. Phys.* **117**, 9160 (2002)): `max(|force|) <= 3e-4`
**and** (`|delta E| < 1e-6` **or** `max(|step|) <= 3e-4`). We use the stricter
form because the published one can accept geometries with a remaining RMS force,
which on machine-learned surfaces end on higher-order saddle points.

---

### `lbfgs`

L-BFGS optimizer settings (extends `opt`).

```yaml
lbfgs:
  # Inherits all opt settings, plus:
 keep_last: 7 # History size for L-BFGS buffers
 beta: 1.0 # Initial damping beta
 gamma_mult: false # Multiplicative gamma update toggle
 max_step: 0.3 # Maximum step length
 control_step: true # Control step length adaptively
 double_damp: true # Double damping safeguard
 mu_reg: null # Regularization strength
 max_mu_reg_adaptions: 10 # Cap on mu adaptations
 reject_uphill: false # Opt in to rejecting energy rises above the tolerance
 uphill_tolerance: 0.0001 # Energy-rise tolerance (Hartree)
 rejection_step_floor: 1.0e-07 # Smallest retry step
 max_rejections_at_floor: 3 # Stop after repeated rejection at the floor
```

---

### `rfo`

Rational Function Optimizer settings (extends `opt`).

```yaml
rfo:
  # Inherits all opt settings, plus:
 trust_radius: 0.10 # Trust-region radius
 trust_update: true # Enable trust-region updates
 trust_min: 0.0001 # Minimum trust radius
 trust_max: 0.10 # Maximum trust radius (bohr)
 max_energy_incr: null # Allowed energy increase per step
 reject_uphill: false # Opt in to rejecting energy rises above the tolerance
 uphill_tolerance: 0.0001 # Energy-rise tolerance (Hartree)
 rejection_trust_floor: 1.0e-07 # Smallest retry trust radius
 max_rejections_at_floor: 3 # Stop after repeated rejection at the floor
 hessian_update: ts_bfgs # Hessian update scheme: ts_bfgs, bfgs, bofill, etc.
 hessian_init: calc # Hessian initialization: calc, unit, etc.
 hessian_recalc: 500 # Rebuild Hessian every N steps
 hessian_recalc_adapt: null # Adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # Eigenvalue threshold for stability
 alpha0: 1.0 # Initial micro step
 max_micro_cycles: 50 # RS iteration limit per step
 rfo_overlaps: false # Enable RFO overlaps
 gediis: false # Enable GEDIIS
 gdiis: true # Enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # Test descent direction before DIIS
 adapt_step_func: true # Adaptive step scaling
```

## Path Optimization Sections

### `gs`

Growing String Method settings.

```yaml
gs:
 fix_first: true # Keep first endpoint fixed
 fix_last: true # Keep last endpoint fixed
 max_nodes: 20 # Maximum string nodes (internal images); for GSM the total path has +2 endpoints
 perp_thresh: 0.005 # Perpendicular displacement threshold
 reparam_check: rms # Reparameterization check metric
 reparam_every: 1 # Reparameterization stride
 reparam_every_full: 1 # Full reparameterization stride
 param: equi # Parameterization scheme
 max_micro_cycles: 10 # RS iteration limit per step
 reset_dlc: true # Rebuild delocalized coordinates each step
 climb: true # Enable climbing image
 climb_rms: 0.0005 # Climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # Keep climbing image fixed
 scheduler: null # Optional scheduler backend
```

```{note}
`gs.max_nodes` / `--max-nodes` is the number of movable internal images for both **GSM** and **DMF**. Both engines retain two endpoints, so the complete path contains `max_nodes + 2` images. See [`path-opt`](path-opt.md).

`gs.param` accepts `equi` or `energy`. Energy weighting is applied only after the GSM string is fully grown and shifts node density toward high-energy regions.
```

---

### `dmf`

Direct Max Flux settings for MEP optimization. DMF builds its initial path with FB-ENM (flat-bottom elastic network model) or, with `correlated: true`, CFB-ENM (correlated FB-ENM).

```yaml
dmf:
 backend: gpu # gpu (dmf.torch / CUDA, default) | cpu (dmf / NumPy)
 max_cycles: 3000 # Maximum DMF/IPOPT iterations (overridden by --dmf-max-iterations)
 tol: tight # IPOPT dual_inf_tol: tight (0.04) | middle (0.10) | loose (0.20) or a positive float (overridden by --dmf-tol)
 correlated: true # Build the initial path with CFB-ENM instead of FB-ENM
 sequential: true # Sequential DMF execution
 fbenm_only_endpoints: false # Run FB-ENM beyond endpoints
 fbenm_options:
   delta_scale: 0.2 # FB-ENM displacement scaling
   bond_scale: 1.25 # Bond cutoff scaling
   fix_planes: true # Enforce planar constraints
 cfbenm_options:
   bond_scale: 1.25 # CFB-ENM bond cutoff scaling
   corr0_scale: 1.1 # Correlation scale for corr0
   corr1_scale: 1.5 # Correlation scale for corr1
   corr2_scale: 1.6 # Correlation scale for corr2
   eps: 0.05 # Correlation epsilon
   pivotal: true # Pivotal residue handling
   single: true # Single-atom pivots
   remove_fourmembered: true # Prune four-membered rings
 dmf_options:
   remove_rotation_and_translation: false # Keep rigid-body motions
   mass_weighted: false # Toggle mass weighting
   parallel: false # Enable parallel DMF
   eps_vel: 0.01 # Velocity tolerance
   eps_rot: 0.01 # Rotational tolerance
   beta: 10.0 # Beta parameter for DMF
   update_teval: false # Update transition evaluation
 ipopt_options: {} # Raw IPOPT options, e.g. {dual_inf_tol: 0.04}
 k_fix: 300.0 # Harmonic constant for restraints (top-level dmf key, NOT under dmf_options)
```

`dmf.tol` is the tolerance the DMF solve applies last, so it takes precedence over an `ipopt_options.dual_inf_tol` set in the same file. Set only `ipopt_options.dual_inf_tol` (and leave `dmf.tol` unset) to pin the raw IPOPT option instead. Gaussian presets such as `gau_tight` are rejected here; they belong to `--thresh` and `--thresh-gsm`.

---

### `search`

Recursive path search settings (path-search only).

```yaml
search:
 max_depth: 10 # Recursive subdivision levels allowed (0 = no subdivision)
 stitch_rmsd_thresh: 0.0001 # RMSD threshold (Bohr) for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold (Bohr) for bridging nodes
 max_nodes_segment: 20 # Max nodes per segment
 max_nodes_bridge: 5 # Max nodes per bridge
 kink_max_nodes: 3 # Max nodes for kink optimizations
 max_seq_kink: 2 # Max sequential kinks
 refine_mode: null # Refinement strategy: peak, minima, or null (auto)
```

---

### `stopt`

StringOptimizer settings for chain-of-states path optimization. `stopt.lbfgs` and `stopt.rfo` set the single-structure optimizers, in the same way as `opt.lbfgs` and `opt.rfo`.

```yaml
stopt:
 type: string # Optimizer type label
 thresh: gau_loose # StringOptimizer convergence preset
 stop_in_when_full: 300 # Cycles allowed after the string is fully grown; then the run stops unconverged
 align: false # Always false in path-opt/path-search; the images are superposed separately by a Kabsch fit
 scale_step: global # Step scaling mode
 max_cycles: 300 # Maximum StringOptimizer iterations
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 reparam_thresh: 0.0 # Reparameterization threshold
 coord_diff_thresh: 0.0 # Coordinate-difference threshold
 out_dir: ./result_path_opt/ # Output directory
 print_every: 10 # Logging stride
```

## TS Optimization Sections

TS optimization uses **two mutually exclusive** algorithm sections, selected by `--opt-mode`:
- `--opt-mode dimer` (or `grad`) → uses `hessian_dimer` section
- `--opt-mode rsprfo` (or `hess`, default), `rsirfo`, or `trim` → uses `rsirfo` section

A key set in only one of `opt` and the active section applies to both;
different values in the two stop the run. With `thresh` in neither, TS optimization uses `baker`,
not the `gau` shown under `opt`.

### `hessian_dimer`

Hessian Guided Dimer TS optimization settings.

```yaml
hessian_dimer:
 thresh_loose: gau_loose # Loose convergence preset
 thresh: baker # Main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # n_imag cutoff (cm⁻¹); see the freq section
 flatten_amp_ang: 0.1 # Flattening amplitude (Å)
 flatten_max_iter: 0 # Flatten rounds (0 = off); --flatten uses 50 when this is 0
 flatten_sep_cutoff: 0.0 # Minimum distance between representative atoms
 flatten_k: 10 # Representative atoms sampled per mode
 flatten_loop_bofill: false # Bofill update for flatten displacements
 mem: 100000 # Memory limit for solver
 device: auto # Device selection for eigensolver
 root: 0 # Targeted TS root index
 dimer:
   length: 0.0189 # Dimer separation (Bohr)
   rotation_max_cycles: 15 # Max rotation iterations
   rotation_method: fourier # Rotation optimizer method
   rotation_thresh: 0.0001 # Rotation convergence threshold
   rotation_tol: 1 # Rotation tolerance factor
   rotation_max_element: 0.001 # Max rotation matrix element
   rotation_interpolate: true # Interpolate rotation steps
   rotation_disable: false # Disable rotations entirely
   rotation_disable_pos_curv: true # Disable when positive curvature detected
   rotation_remove_trans: true # Remove the selected rigid-null components
   trans_force_f_perp: true # Project forces perpendicular to translation
   bonds: null # Bond list for constraints
   N_hessian: null # Hessian size override
   bias_rotation: false # Bias rotational search
   bias_translation: false # Bias translational search
   bias_gaussian_dot: 0.1 # Gaussian bias dot product
   seed: null # RNG seed for rotations
   write_orientations: false # Write rotation orientations (explicit true is allowed)
   forward_hessian: true # Propagate Hessian forward
 lbfgs:                    # sibling of `dimer` under `hessian_dimer`, not nested inside it
   # Same keys as the top-level lbfgs section
   thresh: baker
   line_search: false # Required: Dimer effective force is not energy-conjugate
```

Inner L-BFGS settings live under `hessian_dimer.lbfgs`, not the top-level
`lbfgs` section. Shared `print_every` and `energy_plateau*` values follow the
conflict rule above. `line_search` is fixed to `false`; setting it to `true` is
rejected because Dimer's projected/inverted effective force is not the gradient
of the reported physical energy. `max_cycles` is not set here: each L-BFGS run
between Hessian updates takes at most the cycles left in `opt.max_cycles`.

```{note}
**`flatten_max_iter`.** `--flatten` removes extra imaginary modes in up to
`flatten_max_iter` rounds (50 when the value is 0), and `--no-flatten` turns it
off. With neither flag, a positive value set here turns it on. See
{ref}`When --flatten is on <flatten-precedence-caveat>`.
```

---

### `rsirfo`

RS-I-RFO / RS-P-RFO TS optimization settings.

```yaml
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 100000 # Shared with opt.max_cycles; conflicting explicit values are rejected
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum accepted step norm
 assert_min_step: true # Assert when steps stagnate
 roots: [0] # Exactly one target root index (first-order TS only)
 hessian_ref: null # Reference Hessian
 rx_modes: null # Reaction-mode definitions
 prim_coord: null # Primary coordinates to monitor
 rx_coords: null # Reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: 500 # Rebuild exact Hessian every N macro steps (inherited from rfo)
 hessian_recalc_reset: true # Reset recalc counter after exact Hessian
 max_micro_cycles: 50 # RS iteration limit per step
 augment_bonds: false # Augment reaction path based on bond analysis
 min_line_search: false # Always false: RS-P-RFO does not use line searches
 max_line_search: false # Always false: RS-P-RFO does not use line searches
 assert_neg_eigval: false # Require negative eigenvalue at convergence
 track_mode_by_overlap: false # Track the selected TS mode by overlap with the previous Hessian
 reject_mode_loss: false # Once a TS mode is found, reject steps that lose it and retry with a smaller trust radius
 mode_loss_trust_floor: 1.0e-05 # Smallest trust radius for those retries
 max_mode_loss_rejections: 5 # Rejections allowed at that floor before stopping
 verify_saddle: true # At convergence, count n_imag with an exact Hessian; n_imag = 0 is not accepted
 saddle_imaginary_threshold_cm: 5.0 # n_imag cutoff (cm⁻¹); see the freq section
 saddle_recovery_step: 0.01 # Uphill step used to leave a minimum (n_imag = 0)
 saddle_recovery_check_interval: 50 # Steps between exact-Hessian checks during that recovery
 saddle_recovery_max_cycles: 0 # Maximum recovery steps; 0 turns recovery off
 out_dir: ./result_tsopt/ # Output directory
 # Also inherits rfo-like settings: trust_radius, trust_update, etc.
```

In RS-P-RFO, an explicit `true` for `min_line_search` or `max_line_search`
prints a warning and falls back to `false`. RS-I-RFO and TRIM ignore both keys,
and Dimer uses `hessian_dimer.lbfgs.line_search`.

```{note}
**`--flatten` precedence.** The flatten loop of the RS-P-RFO, RS-I-RFO, and TRIM
paths also reads `hessian_dimer.flatten_max_iter`, with the same rules as in the
note under `hessian_dimer`. See {ref}`When --flatten is on <flatten-precedence-caveat>`.
```

## IRC Section

(irc-section)=
### `irc` (section)

IRC integration settings.

```yaml
irc:
 step_length: 0.1 # Integration step length (Bohr, unweighted Cartesian; --step-size)
 never_stop: false # Ignore physical endpoint criteria and trace to max_cycles
 max_cycles: 125 # Maximum steps along IRC
 forward: true # Propagate in forward direction
 backward: true # Propagate in backward direction
 root: 0 # Normal-mode root index
 hessian_init: calc # Hessian initialization source
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: null # Hessian rebuild cadence
 energy_increase_thresh: 0.0   # Stop on any one-step rise in ordinary mode
 dump_every: null # Disabled; positive cadence writes a coordinate/energy/gradient checkpoint without a Hessian
 dump_fn: irc_data.h5 # Checkpoint filename used only when dump_every is set
 displ: energy # Displacement construction method
 displ_energy: 0.001 # Energy-based displacement scaling
 displ_length: 0.1 # Length-based displacement fallback
 rms_grad_thresh: 0.001 # RMS gradient convergence threshold
 hard_rms_grad_thresh: null # Hard RMS gradient stop
 energy_thresh: 0.000001 # Energy change threshold
 imag_below: 0.0 # IRC starts only if the root mode has ν ≤ this value (cm⁻¹)
 force_inflection: true # Enforce inflection detection
 check_bonds: false # Check bonds during propagation
 out_dir: ./result_irc/ # Output directory
 prefix: "" # Filename prefix
 max_pred_steps: 500 # Predictor-corrector max steps
 loose_cycles: 3 # Loose cycles before tightening
 corr_func: mbs # EulerPC corrector function
```

`corr_func` selects the corrector step of the predictor–corrector IRC integrator (EulerPC). Only `"mbs"` (modified Bulirsch–Stoer) is registered; other values raise a construction error.

## Vibrational Analysis Sections

(freq-section)=
### `freq` (section)

Vibrational frequency analysis settings.

```yaml
freq:
 zero_cutoff_cm: 5.0 # Imaginary modes satisfy ν < -zero_cutoff_cm
 amplitude_ang: 0.8 # Displacement amplitude for modes (Å)
 n_frames: 20 # Number of frames per mode trajectory
 max_write: 10 # Maximum number of modes to write
 sort: value # Sort order: "value" or "abs"
 out_dir: ./result_freq/ # Output directory
```

The default imaginary-mode criterion is ν < −5.00 cm⁻¹. `freq`, `opt`
(flattening), and `tsopt` count n_imag with this `freq.zero_cutoff_cm`;
`irc` does not read it. `tsopt` also takes this cutoff from
`hessian_dimer.neg_freq_thresh_cm` or `rsirfo.saddle_imaginary_threshold_cm`;
if two of the three keys are set explicitly to different values, the run stops
with an error. `n_negative_modes` also counts the negative
frequencies inside the cutoff. Neither n_imag nor `n_negative_modes` decides
whether an optimization has converged. Whatever the cutoff, the output keeps
every signed frequency, and thermochemistry uses every positive mode.

---

### `thermo`

Thermochemistry settings.

```yaml
thermo:
 temperature: 298.15 # Thermochemistry temperature (K)
 pressure_atm: 1.0 # Thermochemistry pressure (atm)
 symmetry_number: null # Auto-detect; a positive integer is an advanced override
 dump: false # Write thermoanalysis.yaml
```

## DFT Section

(dft-section)=
### `dft` (section)

DFT calculation settings.

```yaml
dft:
 func: wb97m-v # Exchange-correlation functional
 basis: def2-svp # Basis set name
 func_basis: null # Combined "FUNC/BASIS" string (overrides func/basis)
 conv_tol: 1.0e-09 # SCF convergence tolerance (hartree)
 max_cycle: 100 # Maximum SCF iterations
 grid_level: 3 # PySCF grid level
 engine: gpu # SCF backend: "gpu" (GPU4PySCF) or "cpu" (PySCF)
 solvent: none # Native PySCF solvent name
 solvent_model: smd # pcm | smd
 pyscf: {} # PySCF object-name attribute forwarding
 lowmem: true # Low-memory direct JK; false enables density fitting
 nprocs: auto # PySCF/OpenMP threads from scheduler/affinity
 memory: auto # Host RAM limit, e.g. 64GB (not GPU VRAM)
 verbose: 0 # PySCF verbosity (0-9); applies at -v 0/1; at the default -v 2 and at -v 3 it is raised to at least 4
 out_dir: ./result_dft/ # Output directory root
```

## Scan Sections

Scan coordinates go in `-s/--scan-lists`, **not** in the `--config` YAML.
See {ref}`Scan-list spec <scan-list-spec>` for the syntax.

(bias-section)=
### `bias`

Harmonic bias settings for `scan`, `scan2d`, and `scan3d`.

```yaml
bias:
 k: 300 # Harmonic bias strength (eV·Å⁻²)
```

**Shared spring constant across subcommands.** The same physical harmonic penalty (`k`, in eV·Å⁻²) appears in the following places with the same default of `300`:

| YAML key | Used by | CLI flag |
|----------|---------|----------|
| `bias.k` | `scan`, `scan2d`, `scan3d` | `--restraint-k` |
| `dmf.k_fix` | `path-opt` / `path-search` with `--mep-mode dmf` | — (YAML only) |

`opt` also applies `--restraint-k` (same default) to its `--distance-restraint` pairs, but reads it only from the CLI, not from the `bias:` section.

Override any of these to tune how stiff the harmonic restraint is. A smaller value (e.g. `20.0`) is appropriate when the geometry should relax against a soft guidance term; the default enforces near-rigid pinning.

---

### `bond`

Bond-change detection from element covalent radii.

```yaml
bond:
 device: auto # Device for the distance calculation: "cuda", "cpu", or "auto"
 bond_factor: 1.2 # Covalent-radius scaling for cutoff
 margin_fraction: 0.05 # Fractional tolerance for comparisons
 delta_fraction: 0.05 # Minimum relative change to flag bond formation/breaking
```

## Example: Complete Configuration File

```yaml
# pdb2reaction configuration example

geom:
 coord_type: cart
 freeze_atoms: []

calc:
 backend: uma
 model: uma-s-1p2 # Model name for the selected backend (UMA: uma-s-1p2 | uma-m-1p1)
 device: auto
 hessian_calc_mode: FiniteDifference # Portable default; benchmark Analytical before opting in

gs:
 max_nodes: 12
 climb: true
 climb_lanczos: true

stopt:
 thresh: gau_loose
 max_cycles: 300
 dump: false

lbfgs:
 max_cycles: 100000

rfo:
 max_cycles: 100000

bond:
 bond_factor: 1.2
 delta_fraction: 0.05

search:
 max_depth: 10
 max_nodes_segment: 20

freq:
 max_write: 10
 amplitude_ang: 0.8

thermo:
 temperature: 298.15
 pressure_atm: 1.0
 symmetry_number: null

dft:
 func: wb97m-v
 basis: def2-svp
 grid_level: 3
```

## Notes

- `workers` and `workers_per_node` take effect only with the UMA backend.
- With `workers > 1`, UMA cannot compute analytical Hessians: an explicit `hessian_calc_mode: Analytical` stops the run with an error. Use `workers: 1` or `FiniteDifference`; see {ref}`Workers and analytical Hessians <workers-analytical-error>`.
- `freq` and `irc` always use the partial Hessian, whatever `calc.return_partial_hessian` is set to.
- `all` passes the same file to each stage it runs, and each stage reads the sections listed for that command in the **Used by** column of the Overview table. For example, the TS stage reads the `tsopt` sections, including `opt`, `hessian_dimer`, and `rsirfo`, and the `dft` section takes effect with `all --dft`.
- `opt.lbfgs` and `opt.rfo` are other names for `lbfgs` and `rfo`, and `freq.thermo` is another name for `thermo`. Two different values for the same setting stop the run with an error, for example `lbfgs.max_cycles` and `opt.lbfgs.max_cycles`, or `opt.max_cycles` and `lbfgs.max_cycles` when L-BFGS is the selected optimizer. `-o/--out-dir` overrides the `out_dir` keys, and `all` sets the output directory of each stage itself.

## See Also

- [all](all.md) - End-to-end workflow
- [opt](opt.md) - Single-structure optimization
- [tsopt](tsopt.md) - Transition state optimization
- [path-search](path-search.md) - Recursive MEP search
- [freq](freq.md) - Vibrational analysis
- [dft](dft.md) - DFT calculations
- [Backends](backends.md) - MLIP backend details
- [Troubleshooting](troubleshooting.md) - Common errors and fixes
