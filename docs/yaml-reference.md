# YAML Reference

```{tip}
**Looking for the `all` end-to-end command?** The `all` command consumes every section described on this page and forwards them to the appropriate subcommand. See [`all`](all.md) (in particular its "Subcommand → YAML Sections" mapping) for which sections are read by which stage.
```

(yaml-configuration-precedence)=
## Configuration precedence

Settings are resolved in the following order (later sources override earlier ones):

```
built-in defaults  <  --config (YAML)  <  CLI flags
```

1. **Built-in defaults** — hard-coded values in `pdb2reaction/core/defaults.py`.
2. **`--config`** — a YAML file that overrides defaults (e.g., `--config my_settings.yaml`).
3. **CLI flags** — explicit command-line options (e.g., `-q -1`, `--thresh gau_loose`). Only *explicitly supplied* flags override YAML; options left at their CLI default do not mask YAML values.

For example, if the YAML sets `charge: 0` but the CLI passes `-q -1`, the charge will be `-1`.

This precedence applies uniformly to `all`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `dft`. See also {ref}`CLI Conventions: Configuration precedence <configuration-precedence>`.

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
| `--max-cycles` | `max_cycles` | `opt` |
| `--dump` | `dump` | `opt` |
| `--opt-mode` | `opt_mode` | `opt`, `tsopt`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `all` |
| `--freeze-atoms` | `freeze_atoms` | `geom` |
| `--coord-type` | `coord_type` | `geom` |
| `--temperature` (freq, `all --freq-temperature`) | `temperature` | `thermo` |
| `--pressure` (freq, `all --freq-pressure`) | `pressure_atm` | `thermo` |
| `--engine` (`dft` subcommand) / `--dft-engine` (`all` wrapper) | `engine` | `dft` |

```{note}
**Name mismatch — `--pressure` vs `pressure_atm`.** On the CLI the flag is `--pressure` (units implicit: atm); the matching YAML key under `thermo:` is `pressure_atm` with an explicit unit suffix. Both carry atm values and get converted to Pa internally.
```

```{note}
**Name mismatch — `--engine` vs `--dft-engine`.** The standalone `dft` subcommand exposes the backend selector as `--engine` (gpu / cpu). In `pdb2reaction all`, to avoid colliding with other engines, the same flag is renamed `--dft-engine` — see {ref}`the --engine vs --dft-engine note in CLI Conventions <engine-vs-dft-engine>`.
```

### Default `--thresh` per subcommand

`--thresh` defaults differ per subcommand because TS optimizers use a tighter "baker" preset while minimizers use the standard "gau" preset.

| Subcommand | Default `--thresh` | Backing defaults block |
|------------|-------------------|------------------------|
| `opt` | `gau` | `OPT_BASE_KW` (→ `lbfgs` / `rfo`) |
| `tsopt` (Hessian Dimer) | `baker` | `HESSIAN_DIMER_KW`, inner `LBFGS_TS_KW` |
| `tsopt` (RS-P-RFO / RS-I-RFO) | `baker` | `RSIRFO_KW` |
| `scan` | `gau` | `OPT_BASE_KW` |
| `scan2d`, `scan3d` | `baker` | `scan_common.py` (`thresh_default="baker"`) |
| `path-search` (per-step opt) | `gau` | `OPT_BASE_KW` |
| `path-opt` / StringOptimizer | `gau_loose` | `STOPT_KW` |
| `all` (pre-opt, post-opt min) | `gau` | `OPT_BASE_KW` |
| `all` (post-opt TS stage) | `baker` | `HESSIAN_DIMER_KW` / `RSIRFO_KW` |

Accepted values: `gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`. Override per run with `--thresh <preset>` or under `opt.thresh` in YAML.

```{note}
**Subcommands without `--thresh`.** `irc`, `freq`, and `dft` do **not** expose `--thresh`:

- `irc` — convergence is governed by `irc.rms_grad_thresh`, `irc.energy_thresh`, and `irc.max_cycles` (see [`irc` section](#irc-section)). The optimizer preset family does not apply because IRC follows a predictor–corrector integrator, not a force-based minimizer.
- `freq` — there is no optimization step, so no `--thresh`. Numerical accuracy is governed by `--hessian-calc-mode` and the underlying MLIP precision.
- `dft` — SCF convergence uses `dft.conv_tol` (default `1e-9` hartree) and `dft.max_cycle`, not the `gau`/`baker` preset family. See the [`dft` section](#dft-section).
```

## Overview


| Section | Description | Used by |
|---------|-------------|---------|
| [`geom`](#geom) | Geometry and coordinate settings | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search, dft |
| [`calc`](#calc) | MLIP backend configuration | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`opt`](#opt) | Shared optimizer settings | opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGS optimizer settings | opt, scan, scan2d, scan3d, path-search, path-opt |
| [`rfo`](#rfo) | RFO optimizer settings | opt, scan, scan2d, scan3d, path-search, path-opt |
| [`gs`](#gs) | Growing String Method settings | path-opt, path-search |
| [`dmf`](#dmf) | Direct Max Flux settings | path-opt, path-search |
| [`stopt`](#stopt) | StringOptimizer settings | path-opt, path-search |
| [`irc`](#irc-section) | IRC integration settings | irc |
| [`freq`](#freq-section) | Vibrational analysis settings | freq |
| [`thermo`](#thermo) | Thermochemistry settings | freq |
| [`dft`](#dft-section) | DFT calculation settings | dft |
| [`bias`](#bias) | Harmonic bias settings | scan, scan2d, scan3d |
| [`bond`](#bond) | Bond-change detection settings | scan, path-search |
| [`search`](#search) | Recursive path search settings | path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian Guided Dimer TS optimization | tsopt |
| [`rsirfo`](#rsirfo) | RS-P-RFO / RS-I-RFO TS optimization | tsopt, all |

## Shared Sections

### `geom`

Geometry loading and coordinate handling.

```yaml
geom:
 coord_type: cart # Coordinate type: "cart" (Cartesian) or "dlc" (delocalized internals)
 freeze_atoms: [] # 1-based atom indices to freeze; if `--freeze-links` is on (PDB input, or XYZ/GJF with `--ref-pdb`), the auto-detected cap-H parent indices are merged in
```

**Notes:**
- `freeze_atoms` from YAML is merged with atoms detected via `--freeze-links` for PDB inputs
- Frozen atoms have zeroed forces; their Hessian columns are also zeroed
- For `irc`, `geom.coord_type` is forced to `cart` after YAML/CLI merging

---

### `calc`

MLIP backend configuration. Supports multiple backends (UMA, ORB, MACE, AIMNet2) and optional xTB solvent corrections.

```yaml
calc:
 backend: uma           # MLIP backend: "uma", "orb", "mace", or "aimnet2"
 precision: auto # auto (per backend: uma fp32, orb/mace fp64) | fp32 | fp64; routed to the backend-native kwarg (ORB precision, MACE default_dtype)
 charge: 0 # Total system charge (overridden by CLI -q)
 spin: 1 # Spin multiplicity 2S+1 (overridden by CLI -m)
 model: uma-s-1p2 # uma-s-1p2 | uma-m-1p1
 task_name: omol # Task tag recorded in UMA batches
 device: auto # Device: "cuda", "cpu", or "auto"
 max_neigh: null # Maximum neighbors for graph construction
 radius: null # Cutoff radius for neighbor search
 r_edges: false # Store radial edges
 workers: 1 # UMA inference workers (workers>1 + explicit Analytical is an error)
 workers_per_node: 1 # Workers per node for parallel predictor
 out_hess_torch: true # Return Hessian as torch.Tensor
 hessian_double: true # Assemble/return Hessian in float64
 # freeze_atoms: null # Inherited from geom.freeze_atoms; do not set directly
 hessian_calc_mode: FiniteDifference # Hessian mode: "Analytical" or "FiniteDifference"
 return_partial_hessian: true  # Return active-DOF block Hessian
 print_timing: true # Print Hessian timing breakdown
 print_vram: true # Print CUDA VRAM usage during Hessian (UMA backend only)
 # Solvent correction (xTB)
 solvent: none           # Implicit solvent name (e.g. "water", "methanol") or "none" to disable
 solvent_model: alpb     # xTB solvent model: "alpb" or "cpcmx"
 xtb_cmd: xtb            # Path to xTB executable
 xtb_acc: 0.2            # xTB accuracy parameter
```

**Notes:**
- `backend` selects the MLIP engine. All backends (UMA, ORB, MACE, AIMNet2) support both analytical (autograd) and finite-difference Hessians; multi-worker inference is UMA-only.
- `workers` / `workers_per_node` are effective with the UMA backend only.
- `solvent` enables xTB-based implicit solvent corrections (delta correction approach). Requires `xtb` to be installed.
- `hessian_calc_mode: Analytical` is recommended when sufficient VRAM is available
- `workers > 1` disables analytical Hessians for the UMA parallel predictor. An explicit `hessian_calc_mode: Analytical` request raises `RuntimeError`; use `workers = 1` or select `FiniteDifference`. See {ref}`the MLIP Calculator hessian-evaluation note <hessian-evaluation>` for details.
- Charge/spin inherit `.gjf` template metadata when available
- `freq` sets `calc.return_partial_hessian = true` by default (PHVA); YAML can override.
- IRC forces `geom.coord_type = cart` and `calc.return_partial_hessian = true` regardless of YAML (partial Hessian with active-DOF processing).

---

### `opt`

Shared single-structure optimizer controls used by both L-BFGS and RFO.

```yaml
opt:
 thresh: gau # Convergence preset: gau_loose, gau, gau_tight, gau_vtight, baker, never
 max_cycles: 10000 # Maximum optimizer iterations
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum step norm for acceptance
 assert_min_step: true # Stop if steps fall below threshold
 rms_force: null # Explicit RMS force target
 rms_force_only: false # Rely only on RMS force convergence
 max_force_only: false # Rely only on max force convergence
 force_only: false # Skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when converging to reference geometry
 overachieve_factor: 0.0 # Factor to tighten thresholds
 check_eigval_structure: false # Validate Hessian eigenstructure
 line_search: true # Enable line search
 energy_plateau: true # Declare convergence when the energy range flattens (see note below)
 energy_plateau_thresh: 1.0e-04 # au (~0.06 kcal/mol); fallback convergence threshold for the plateau check
 energy_plateau_window: 50 # Number of most recent steps inspected for the plateau check
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 prefix: "" # Filename prefix
 out_dir: ./result_opt/ # Output directory
```

**Energy plateau fallback:**
When `energy_plateau: true`, the optimizer declares convergence if the energy range
(max − min) over the last `energy_plateau_window` steps falls below
`energy_plateau_thresh` (default `1e-4` au ≈ 0.06 kcal/mol over 50 steps). This prevents
wasted cycles when MLIP force noise can exceed the `baker` threshold
(`max_force = 3×10⁻⁴ au`), so the force criterion may never be satisfied
even after the energy has flattened.
The fallback is **skipped** for chain-of-states (COS) optimizers such as `stopt`,
`gs`, and DMF, because those store per-image energy arrays rather than a single scalar
trace.

**Convergence Presets:**

| Preset | Max Force | RMS Force | Max Step | RMS Step |
|--------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

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
 hessian_update: bfgs # Hessian update scheme: bfgs, bofill, etc.
 hessian_init: calc # Hessian initialization: calc, unit, etc.
 hessian_recalc: 500 # Rebuild Hessian every N steps
 hessian_recalc_adapt: null # Adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # Eigenvalue threshold for stability
 alpha0: 1.0 # Initial micro step
 max_micro_cycles: 50 # Micro-iteration limit
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
 max_micro_cycles: 10 # Micro-iteration limit
 reset_dlc: true # Rebuild delocalized coordinates each step
 climb: true # Enable climbing image
 climb_rms: 0.0005 # Climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # Keep climbing image fixed
 scheduler: null # Optional scheduler backend
```

```{note}
**`max_nodes` semantics differ between GSM and DMF.** For **GSM** (`mep-mode gsm`), `gs.max_nodes` is the number of **internal images** — the total path length is `max_nodes + 2` because the two endpoints are fixed. For **DMF** (`mep-mode dmf`), the CLI flag `--max-nodes` instead counts the **number of movable images** with no implicit endpoint expansion. See [`path-opt`](path-opt.md) for per-method details.
```

---

### `dmf`

Direct Max Flux settings for MEP optimization.

```{note}
**`--max-nodes` for DMF means "number of movable images"** — unlike GSM, DMF does not add +2 for fixed endpoints. See the `gs.max_nodes` note above for the contrast.
```

```yaml
dmf:
 max_cycles: 300 # Maximum DMF/IPOPT iterations (overridden by --max-cycles)
 correlated: true # Correlated DMF propagation
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
 k_fix: 300.0 # Harmonic constant for restraints (top-level dmf key, NOT under dmf_options)
```

---

### `search`

Recursive path search settings (path-search only).

```yaml
search:
 max_depth: 10 # Recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 20 # Max nodes per segment
 max_nodes_bridge: 5 # Max nodes per bridge
 kink_max_nodes: 3 # Max nodes for kink optimizations
 max_seq_kink: 2 # Max sequential kinks
 refine_mode: null # Refinement strategy: peak, minima, or null (auto)
```

---

### `stopt`

StringOptimizer settings for chain-of-states path optimization (`path-opt`, `path-search`).

```yaml
stopt:
 type: string # Optimizer type label
 thresh: gau_loose # StringOptimizer convergence preset
 stop_in_when_full: 300 # Early stop threshold when the string is full
 align: false # Alignment toggle (forced to False in path-opt/path-search; external Kabsch alignment is used instead)
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

Shared optimizer settings (`thresh`, `max_cycles`, `dump`) are read from the `opt` section above.

### `hessian_dimer`

Hessian Guided Dimer TS optimization settings (tsopt --opt-mode grad).

```yaml
hessian_dimer:
 thresh_loose: gau_loose # Loose convergence preset
 thresh: baker # Main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # Imaginary-frequency detection threshold (cm⁻¹)
 flatten_amp_ang: 0.1 # Flattening amplitude (Å)
 flatten_max_iter: 50 # Flattening iteration cap (see note below)
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
   rotation_remove_trans: true # Remove translational components
   trans_force_f_perp: true # Project forces perpendicular to translation
   bonds: null # Bond list for constraints
   N_hessian: null # Hessian size override
   bias_rotation: false # Bias rotational search
   bias_translation: false # Bias translational search
   bias_gaussian_dot: 0.1 # Gaussian bias dot product
   seed: null # RNG seed for rotations
   write_orientations: true # Write rotation orientations
   forward_hessian: true # Propagate Hessian forward
   lbfgs:
     # Same keys as lbfgs section
     thresh: baker
     max_cycles: 10000
```

```{note}
**`flatten_max_iter` CLI precedence exception.** The YAML value `hessian_dimer.flatten_max_iter: 50` (and the analogous `rsirfo.flatten_max_iter` override) is **overridden to `0` by the CLI unless `--flatten` is explicitly passed**, regardless of the normal `defaults < YAML < CLI` ordering. See {ref}`flatten-precedence-caveat` for the full behavior table.
```

---

### `rsirfo`

RS-I-RFO / RS-P-RFO TS optimization settings (used by tsopt `--opt-mode rsirfo`, `rsprfo` (the `hess` default), and `trim`).

```yaml
rsirfo:
 thresh: baker # RS-IRFO convergence preset
 max_cycles: 10000 # Iteration cap
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum accepted step norm
 assert_min_step: true # Assert when steps stagnate
 roots: [0] # Target root indices
 hessian_ref: null # Reference Hessian
 rx_modes: null # Reaction-mode definitions
 prim_coord: null # Primary coordinates to monitor
 rx_coords: null # Reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: 500 # Rebuild exact Hessian every N macro steps (inherited from rfo)
 hessian_recalc_reset: true # Reset recalc counter after exact Hessian
 max_micro_cycles: 50 # Micro-iterations per macro cycle
 augment_bonds: false # Augment reaction path based on bond analysis
 min_line_search: true # Enforce minimum line-search step
 max_line_search: true # Enforce maximum line-search step
 assert_neg_eigval: false # Require negative eigenvalue at convergence
 track_mode_by_overlap: false # Track the selected TS mode by overlap with the previous Hessian
 saddle_recovery_check_interval: 50 # Exact PHVA cadence during n_imag=0 recovery
 saddle_recovery_max_cycles: 200 # Bounded n_imag=0 recovery cap
 # Also inherits rfo-like settings: trust_radius, trust_update, etc.
```

```{note}
**`--flatten` precedence.** The flatten loop for both Hessian-Dimer and RS-I-RFO paths is configured under the `hessian_dimer:` YAML section (key `flatten_max_iter`, default 50); `rsirfo:` does not define its own flatten counter. The CLI overrides `flatten_max_iter` to `0` unless `--flatten` is explicitly passed on the command line. See {ref}`flatten-precedence-caveat`.
```

## IRC Section

(irc-section)=
### `irc` (section)

IRC integration settings.

```yaml
irc:
 step_length: 0.1 # Integration step length
 never_stop: false # Ignore transient energy rise/plateau stops; convergence/max_cycles still apply
 max_cycles: 125 # Maximum steps along IRC
 downhill: false # Follow downhill direction only
 forward: true # Propagate in forward direction
 backward: true # Propagate in backward direction
 root: 0 # Normal-mode root index
 hessian_init: calc # Hessian initialization source
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: null # Hessian rebuild cadence
 displ: energy # Displacement construction method
 displ_energy: 0.001 # Energy-based displacement scaling
 displ_length: 0.1 # Length-based displacement fallback
 rms_grad_thresh: 0.001 # RMS gradient convergence threshold
 hard_rms_grad_thresh: null # Hard RMS gradient stop
 energy_thresh: 0.000001 # Energy change threshold
 imag_below: 0.0 # Imaginary frequency cutoff
 force_inflection: true # Enforce inflection detection
 check_bonds: false # Check bonds during propagation
 out_dir: ./result_irc/ # Output directory
 prefix: "" # Filename prefix
 max_pred_steps: 500 # Predictor-corrector max steps
 loose_cycles: 3 # Loose cycles before tightening
 corr_func: mbs # EulerPC corrector function (only "mbs" is currently registered)
```

The `corr_func` key selects the corrector step used by the predictor–corrector IRC integrator (EulerPC). Only `"mbs"` (the pysisyphus-native modified Bulirsch–Stoer implementation, default) is currently registered; other values raise a construction error.

## Vibrational Analysis Sections

(freq-section)=
### `freq` (section)

Vibrational frequency analysis settings.

```yaml
freq:
 amplitude_ang: 0.8 # Displacement amplitude for modes (Å)
 n_frames: 20 # Number of frames per mode animation
 max_write: 10 # Maximum number of modes to write
 sort: value # Sort order: "value" or "abs"
 out_dir: ./result_freq/ # Output directory
```

---

### `thermo`

Thermochemistry settings.

```yaml
thermo:
 temperature: 298.15 # Thermochemistry temperature (K)
 pressure_atm: 1.0 # Thermochemistry pressure (atm)
 dump: false # Write thermoanalysis.yaml
```

## DFT Section

(dft-section)=
### `dft` (section)

DFT calculation settings.

```yaml
dft:
 func: wb97m-v # Exchange-correlation functional
 basis: def2-tzvpd # Basis set name
 func_basis: null # Combined "FUNC/BASIS" string (overrides func/basis)
 conv_tol: 1.0e-09 # SCF convergence tolerance (hartree)
 max_cycle: 100 # Maximum SCF iterations
 grid_level: 3 # PySCF grid level
 engine: gpu # SCF backend: "gpu" (GPU4PySCF) or "cpu" (PySCF)
 lowmem: true # Use gpu4pyscf rks_lowmem.RKS for closed-shell GPU runs
 verbose: 0 # PySCF verbosity (0-9); CLI -v 2/3 raises runtime PySCF verbosity to >=4
 out_dir: ./result_dft/ # Output directory root
```

## Scan Sections

Scan coordinates are specified via `--scan-lists/-s` (inline or YAML file), **not** in the main YAML config.
See [Quickstart: scan workflow](quickstart-scan.md) for scan coordinate syntax (PDB selectors, multi-stage).

(bias-section)=
### `bias`

Harmonic bias settings for scans and restraint-based optimizations.

```yaml
bias:
 k: 300 # Harmonic bias strength (eV·Å⁻²)
```

**Shared spring constant across subcommands.** The same physical harmonic penalty (`k`, in eV·Å⁻²) appears in the following places with the same default of `300`:

| YAML key | Used by | CLI flag |
|----------|---------|----------|
| `bias.k` | `scan`, `scan2d`, `scan3d` | `--bias-k` |
| `dmf.k_fix` | `path-opt` / `path-search` when `mep_mode: dmf` | — (YAML only) |

`opt` also accepts `--bias-k` (applied to `--dist-freeze` pairs) but reads it only from the CLI flag, which defaults to the same `300.0` constant; it does not honor the `bias:` YAML section.

Override any of these to tune how stiff the harmonic restraint is. A smaller value (e.g. `20.0`) is appropriate when the geometry should relax against a soft guidance term; the default `300.0` enforces near-rigid pinning.

---

### `bond`

MLIP-based bond-change detection.

```yaml
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # Covalent-radius scaling for cutoff
 margin_fraction: 0.05 # Fractional tolerance for comparisons
 delta_fraction: 0.05 # Minimum relative change to flag bond formation/breaking
```

## Example: Complete Configuration File

Below is a complete example combining multiple sections:

```yaml
# pdb2reaction configuration example

geom:
 coord_type: cart
 freeze_atoms: []

calc:
 backend: uma
 charge: 0
 spin: 1
 model: uma-s-1p2 # uma-s-1p2 | uma-m-1p1
 device: auto
 hessian_calc_mode: Analytical # Recommended when VRAM permits
 solvent: none                 # Set to e.g. "water" for implicit solvent

gs:
 max_nodes: 12
 climb: true
 climb_lanczos: true

stopt:
 thresh: gau_loose
 max_cycles: 300
 dump: false
 out_dir: ./result_all/

opt:
 thresh: gau
 lbfgs:                     # also accepted at top-level as `lbfgs:`
   max_cycles: 10000
 rfo:                       # also accepted at top-level as `rfo:`
   max_cycles: 10000

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

dft:
 func: wb97m-v
 basis: def2-tzvpd
 grid_level: 3
```

## See Also

- [all](all.md) - End-to-end workflow
- [opt](opt.md) - Single-structure optimization
- [tsopt](tsopt.md) - Transition state optimization
- [path-search](path-search.md) - Recursive MEP search
- [freq](freq.md) - Vibrational analysis
- [dft](dft.md) - DFT calculations
- [uma-pysis](uma-pysis.md) - MLIP backend details
