# `path-search`

## Overview

> **Summary:** Build a continuous MEP from **two or more** structures with GSM (default) or DMF (`--mep-mode dmf`). Automatically refines only regions with bond changes and exports the highest-energy image (HEI) as a TS candidate (validate with tsopt + IRC).

### At a glance
- **Use when:** You have R → … → P structures (2+ inputs) and want a single stitched MEP with automatic refinement.
- **Method:** Chains GSM/DMF segments and recursively refines only sub-intervals that still contain covalent changes.
- **Outputs:** `mep_trj.xyz` (main trajectory), `summary.json` (segment-by-segment results), and optional plots/merged PDBs when enabled.
- **Defaults:** `--mep-mode gsm`, `--opt-mode grad` (LBFGS), `--no-preopt`, `--align`, `--thresh gau`, `--thresh-stopt gau_loose`.
- **Next step:** HEI output alone does **not** validate a TS. Follow with [tsopt](tsopt.md) (includes imaginary-frequency check) and [irc](irc.md).

`pdb2reaction path-search` builds a continuous minimum-energy path (MEP) across two or more structures using GSM (default) or DMF (`--mep-mode dmf`). It selectively refines only those regions where covalent bond changes are detected, then stitches the resolved subpaths into a single trajectory.


When `--convert-files` is enabled (default), the command mirrors trajectories into `.pdb` companions when PDB references exist, and writes `.gjf` companions for HEI snapshots when Gaussian templates exist. For XYZ/GJF inputs, `--ref-pdb` supplies an active site model (binding pocket)-level PDB topology while keeping XYZ coordinates, and `--ref-full-pdb` enables full-template merges (XYZ/GJF inputs still do not produce PDB companions).

The recursive decomposition automatically detects multistep reactions and builds a detailed MEP for each elementary step.  However, complex multistep mechanisms may require manual trial-and-error—adjusting input intermediates, scan specifications, or convergence thresholds—to obtain a satisfactory pathway.

If you only have **two** endpoints and do not need recursive refinement, [path-opt](path-opt.md) is the simpler option.

## Minimal example

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_search
```

## Output checklist

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.json`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png` (when plotting succeeds)

## Common examples

1. Provide explicit intermediates for a multistep path.

```bash
pdb2reaction path-search -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 -m 1 \
 --out-dir ./result_path_search_multi
```

2. Enable merged full-system outputs with template references.

```bash
pdb2reaction path-search -i R.pdb IM1.pdb P.pdb -q 0 -m 1 \
 --ref-full-pdb holo_template.pdb --out-dir ./result_path_search_merge
```

3. Use DMF mode with minima refinement.

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --refine-mode minima --out-dir ./result_path_search_dmf
```

## Usage
```bash
pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1]
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx]
 [--workers N] [--workers-per-node N]
 [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--thresh PRESET] [--thresh-stopt PRESET]
 [--refine-mode {peak|minima}]
 [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
 [--opt-mode grad|hess] [--dump/--no-dump]
 [--out-dir DIR] [--preopt/--no-preopt]
 [--align/--no-align] [--ref-full-pdb FILE...] [--ref-pdb FILE...]
 [--convert-files/--no-convert-files]
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

### Examples
- **Active site model-only** MEP between two endpoints:
  ```bash
  pdb2reaction path-search -i reactant.pdb product.pdb -q 0
  ```
- **Multistep** search with YAML overrides and merged full-system output:
  ```bash
  pdb2reaction path-search \
  -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 \
  --ref-full-pdb holo_template.pdb --out-dir ./run_ps
  ```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more structures in reaction order (reactant → product). Pass all files after a single `-i`/`--input`. | Required |
| `-q, --charge INT` | Net charge. Required for non-`.gjf` inputs unless `--ligand-charge/-l` derivation succeeds (PDB inputs). Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only; `workers_per_node` forwarded to the parallel predictor). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | When loading PDB active site models, freeze the parent atoms of link hydrogens. See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-nodes INT` | Internal nodes per MEP segment (GSM string images or DMF images). | `20` |
| `--max-cycles INT` | Maximum MEP optimization cycles (GSM/DMF). | `300` |
| `--climb/--no-climb` | Enable climbing image for GSM segments (bridge segments always run without climbing). | `True` |
| `--opt-mode TEXT` | Single-structure optimizer for HEI±1/kink nodes. `grad` maps to LBFGS; `hess` maps to RFO. See {ref}`opt-mode-semantics` for how the same token maps across subcommands (tsopt uses Dimer/RS-I-RFO, not L-BFGS/RFO). | `grad` |
| `--mep-mode {gsm\|dmf}` | Segment generator: GSM (string-based) or DMF (direct flux). | `gsm` |
| `--refine-mode {peak\|minima}` | Seeds for refinement: `peak` optimizes HEI±1; `minima` searches outward from the HEI toward the nearest local minima on each side. Defaults to `peak` for GSM and `minima` for DMF when omitted. | _Auto_ |
| `--dump/--no-dump` | Dump MEP (GSM/DMF) and single-structure trajectories. Restart YAML is written only when enabled in YAML. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB or Gaussian inputs. | `True` |
| `-o, --out-dir TEXT` | Output directory. | `./result_path_search/` |
| `--thresh TEXT` | Override convergence preset for single-structure optimizations only (`opt.lbfgs/rfo.thresh`). | `gau` |
| `--thresh-stopt TEXT` | Override convergence preset for the string optimizer (`stopt.thresh`). | `gau_loose` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layer metadata) and continue. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running path search. | `False` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint before MEP search. **Scope-dependent default:** `False` under standalone `path-search`; **flipped to `True` when invoked via `pdb2reaction all`** (see [`all` → MEP Search Options](all.md#mep-search-options)). | `False` |
| `--align/--no-align` | Align all inputs to the first structure before searching. | `True` |
| `--ref-full-pdb PATH...` | Full-size template PDBs (one per input, unless `--align` lets you reuse the first). | _None_ |
| `--ref-pdb PATH...` | Active site model reference PDBs used for the final full-system merge when inputs are XYZ/GJF (one per input, matching input order). | _None_ |

## Workflow
1. **Initial segment per pair (GSM/DMF)** – run `GrowingString` or DMF between each adjacent input (A→B) to obtain a coarse MEP and identify the highest-energy image (HEI).
2. **Local relaxation around HEI** – refine either HEI ± 1 (`refine-mode=peak`) or the nearest local minima on each side of the HEI (`refine-mode=minima`) with the chosen single-structure optimizer (`opt-mode`) to recover nearby minima (`End1`, `End2`).
    > **Default:** When `--refine-mode` is omitted, it defaults to `peak` for GSM and `minima` for DMF.
3. **Decide between kink vs. refinement**:
 - If no covalent bond change is detected between `End1` and `End2`, treat the region as a *kink* -- a conformational rearrangement with no bond breaking or formation (see [Glossary](glossary.md)): insert `search.kink_max_nodes` linear nodes and optimize each individually.
 - Otherwise, the region is a *reactive segment* -- a segment in which covalent bond changes are detected between the endpoints (see [Glossary](glossary.md)). Launch a **refinement segment (GSM/DMF)** between `End1` and `End2` to sharpen the barrier.
4. **Selective recursion** – compare bond changes for `(A→End1)` and `(End2→B)` using the `bond` thresholds. Recurse only on sub-intervals that still contain covalent updates. Recursion depth is capped by `search.max_depth`.
5. **Stitching & bridging** – concatenate resolved subpaths, dropping duplicate endpoints when RMSD ≤ `search.stitch_rmsd_thresh`. If the RMSD gap between two stitched pieces exceeds `search.bridge_rmsd_thresh`, insert a *bridge segment* -- a connecting segment between two non-adjacent intermediates (see [Glossary](glossary.md)) -- using GSM/DMF. When the interface itself shows a bond change, a brand-new recursive segment replaces the bridge.
6. **Alignment & merging (optional)** – with `--align` (default), pre-optimized structures are rigidly aligned to the first input and `freeze_atoms` are reconciled. Provide `--ref-full-pdb` to merge active site model trajectories back into full-size PDB templates (one template per input unless alignment allows reuse of the first file).

Bond-change detection relies on `bond_changes.compare_structures` with thresholds surfaced under the `bond` YAML section. MLIP backends are constructed once and shared across all structures for efficiency.

## Outputs
```
out_dir/ (default:./result_path_search/)
├─ mep_trj.xyz # Primary MEP trajectory
├─ mep.pdb # PDB companion when inputs were PDB templates and conversion is enabled
├─ mep_w_ref.pdb # Merged full-system MEP (requires ref PDB/template)
├─ mep_w_ref_seg_XX.pdb # Merged per-segment paths when covalent changes exist (requires ref PDB)
├─ summary.json # Barrier and classification summary for every recursive segment
├─ summary.log # Text summary
├─ mep_plot.png # ΔE profile generated via `trj2fig` (kcal/mol, reactant reference)
├─ energy_diagram_MEP.png # Static export of the MEP state-energy diagram (relative to reactant)
└─ seg_000_*/ # GSM/DMF dumps, HEI snapshots, kink/refinement diagnostics per segment
```
- Console reports covering resolved configuration blocks (`geom`, `calc`, `gs`, `stopt`, `opt.*`, `bond`, `search`).

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Provide at least two inputs; otherwise the command exits with an "invalid value" error for `-i/--input`.
- Repeat `--ref-full-pdb` once per file when providing multiple templates; with `--align`, only the first template is reused for merges.
- All MLIP backends are shared across structures for efficiency.
- When `--dump` is set, MEP (GSM/DMF) and single-structure optimizations emit trajectories. Restart YAML is written only when `dump_restart` is enabled in YAML.

See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.
The YAML root must be a mapping. Shared sections reuse [YAML Reference](yaml-reference.md): `geom`/`calc` mirror single-structure options (with `--freeze-links` augmenting `geom.freeze_atoms` for PDBs), and `stopt` inherits the StringOptimizer knobs documented for `path-opt` (see [path-opt.md](path-opt.md)).

```{note}
**Reference duplication.** The YAML keys for `geom`, `calc`, `gs`, `dmf`, `stopt`, `opt.lbfgs`, and `opt.rfo` listed below mirror the canonical definitions in [YAML Reference](yaml-reference.md). When the two pages disagree, the canonical [YAML Reference](yaml-reference.md) entries (and `pdb2reaction/defaults.py`) take precedence; the appendix on this page is reproduced inline only for `path-search`-specific defaults (e.g. `out_dir: ./result_path_search/`) and convenience lookup.
```

`gs` (Growing String) inherits defaults from `pdb2reaction.path_opt.GS_KW` with overrides for `max_nodes` (internal nodes per segment), climb behavior (`climb`, `climb_rms`, `climb_fixed`), and reparameterization cadence (`reparam_every_full`, `reparam_check`).

`opt` houses the single-structure optimizers used for HEI±1 and kink nodes, split into `lbfgs` and `rfo` subsections. Each subsection mirrors [YAML Reference](yaml-reference.md) but defaults to `out_dir: ./result_path_search/` and `dump: False`.

`bond` carries the MLIP-based bond-change detection parameters shared with {ref}`scan <section-bond>`: `device`, `bond_factor`, `margin_fraction`, and `delta_fraction`.


`dmf` bundles Direct Max Flux + (C)FB-ENM controls applied whenever `--mep-mode dmf` is selected. The defaults mirror the shared `DMF_KW` dictionary and can be overridden per run:

### Example YAML (default value)
```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # net charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 task_name: omol # UMA task name
 device: auto # MLIP device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true  # partial Hessian over active DOFs
gs:
 fix_first: true # keep the first endpoint fixed during optimization
 fix_last: true # keep the last endpoint fixed during optimization
 max_nodes: 20 # maximum string nodes
 perp_thresh: 0.005 # perpendicular displacement threshold
 reparam_check: rms # reparametrization check metric
 reparam_every: 1 # reparametrization stride
 reparam_every_full: 1 # full reparametrization stride
 param: equi # parametrization scheme
 max_micro_cycles: 10 # micro-iteration limit
 reset_dlc: true # rebuild delocalized coordinates each step
 climb: true # enable climbing image
 climb_rms: 0.0005 # climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # keep climbing image fixed
 scheduler: null # optional scheduler backend
stopt:
 type: string # optimizer type label
 thresh: gau_loose # StringOptimizer convergence preset
 stop_in_when_full: 300 # early stop threshold when string is full
 align: false # alignment toggle (kept off)
 scale_step: global # step scaling mode
 max_cycles: 300 # maximum optimization cycles
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 reparam_thresh: 0.0 # reparametrization threshold
 coord_diff_thresh: 0.0 # coordinate difference threshold
 out_dir: ./result_path_search/ # output directory
 print_every: 10 # logging stride
dmf:
 max_cycles: 300 # maximum DMF/IPOPT iterations
 correlated: true # correlated DMF propagation
 sequential: true # sequential DMF execution
 fbenm_only_endpoints: false # run FB-ENM beyond endpoints
 fbenm_options:
   delta_scale: 0.2 # FB-ENM displacement scaling
   bond_scale: 1.25 # bond cutoff scaling
   fix_planes: true # enforce planar constraints
 cfbenm_options:
   bond_scale: 1.25 # CFB-ENM bond cutoff scaling
   corr0_scale: 1.1 # correlation scale for corr0
   corr1_scale: 1.5 # correlation scale for corr1
   corr2_scale: 1.6 # correlation scale for corr2
   eps: 0.05 # correlation epsilon
   pivotal: true # pivotal residue handling
   single: true # single-atom pivots
   remove_fourmembered: true # prune four-membered rings
 dmf_options:
   remove_rotation_and_translation: false # keep rigid-body motions
   mass_weighted: false # toggle mass weighting
   parallel: false # enable parallel DMF
   eps_vel: 0.01 # velocity tolerance
   eps_rot: 0.01 # rotational tolerance
   beta: 10.0 # beta parameter for DMF
   update_teval: false # update transition evaluation
   k_fix: 300.0 # harmonic constant for restraints
opt:
 lbfgs:
   thresh: gau # LBFGS convergence preset
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
   out_dir: ./result_path_search/ # output directory
   keep_last: 7 # history size for LBFGS buffers
   beta: 1.0 # initial damping beta
   gamma_mult: false # multiplicative gamma update toggle
   max_step: 0.3 # maximum step length
   control_step: true # control step length adaptively
   double_damp: true # double damping safeguard
   mu_reg: null # regularization strength
   max_mu_reg_adaptions: 10 # cap on mu adaptations
 rfo:
   thresh: gau # RFOptimizer convergence preset
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
   out_dir: ./result_path_search/ # output directory
   trust_radius: 0.10 # trust-region radius
   trust_update: true # enable trust-region updates
   trust_min: 0.0001 # minimum trust radius
   trust_max: 0.10 # maximum trust radius
   max_energy_incr: null # allowed energy increase per step
   hessian_update: bfgs # Hessian update scheme
   hessian_init: calc # Hessian initialization source
   hessian_recalc: 500 # rebuild Hessian every N steps
   hessian_recalc_adapt: null # adaptive Hessian rebuild factor
   small_eigval_thresh: 1.0e-08 # eigenvalue threshold for stability
   alpha0: 1.0 # initial micro step
   max_micro_cycles: 50 # micro-iteration limit
   rfo_overlaps: false # enable RFO overlaps
   gediis: false # enable GEDIIS
   gdiis: true # enable GDIIS
   gdiis_thresh: 0.0025 # GDIIS acceptance threshold
   gediis_thresh: 0.01 # GEDIIS acceptance threshold
   gdiis_test_direction: true # test descent direction before DIIS
   adapt_step_func: true # adaptive step scaling toggle
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
search:
 max_depth: 10 # recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 10 # max nodes per segment
 max_nodes_bridge: 5 # max nodes per bridge
 kink_max_nodes: 3 # max nodes for kink optimizations
 max_seq_kink: 2 # max sequential kinks
 refine_mode: null # optional refinement strategy (auto-chooses when null)
```

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [path-opt](path-opt.md) — Single-pass MEP optimization (no recursive refinement)
- [tsopt](tsopt.md) — Optimize the HEI as a transition state
- [extract](extract.md) — Generate active site model PDBs for path-search inputs
- [all](all.md) — End-to-end workflow (uses recursive path-search by default; `--refine-path False` for single-pass path-opt)
- [YAML Reference](yaml-reference.md) — Full `gs`, `dmf`, `bond`, `search` configuration options
- [Glossary](glossary.md) — Definitions of MEP, GSM, DMF, HEI
