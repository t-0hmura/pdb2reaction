# `path-search` subcommand

## Overview
Construct a continuous minimum-energy path (MEP) across **two or more** structures ordered along a reaction coordinate. `path-search` chains together Growing String Method (GSM) segments, selectively refines only those regions with covalent changes, and (optionally) merges PDB pockets back into full-size templates. With `--mep-mode dmf`, the same recursive workflow runs using DMF-generated segments instead of GSM, and **GSM is now the default**. Format-aware conversions mirror trajectories and HEI snapshots into `.pdb` or multi-geometry `.gjf` companions when `--convert-files` is enabled (default) and matching templates exist.

## Usage
```bash
pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb -q CHARGE [--multiplicity 2S+1]
                         [--mep-mode {gsm|dmf}] [--freeze-links BOOL] [--thresh PRESET]
                         [--refine-mode {peak|minima}]
                         [--max-nodes N] [--max-cycles N] [--climb BOOL]
                         [--opt-mode light|heavy] [--dump BOOL]
                         [--out-dir DIR] [--preopt BOOL]
                         [--align/--no-align] [--ref-pdb FILE ...]
                         [--convert-files/--no-convert-files]
                         [--args-yaml FILE]
```

### Examples
- **Pocket-only** MEP between two endpoints:
  ```bash
  pdb2reaction path-search -i reactant.pdb product.pdb -q 0
  ```
- **Multistep** search with YAML overrides and merged full-system output:
  ```bash
  pdb2reaction path-search \
      -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 \
      --args-yaml params.yaml --ref-pdb holo_template.pdb --out-dir ./run_ps
  ```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more structures in reaction order (reactant → product). Repeat `-i` or pass multiple paths after one flag. | Required |
| `-q, --charge INT` | Total charge. Required unless the first input is a `.gjf` template that already stores charge. Overrides `--ligand-charge` when both are set. | Required when not in template |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex even when pockets are skipped. | `None` |
| `--workers`, `--workers-per-nodes` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_nodes` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | Explicit `True`/`False`. When loading PDB pockets, freeze the parent atoms of link hydrogens. | `True` |
| `--max-nodes INT` | Internal nodes for GSM segments (`String` has `max_nodes + 2` images). | `10` |
| `--max-cycles INT` | Maximum GSM optimization cycles. | `300` |
| `--climb BOOL` | Explicit `True`/`False`. Enable climbing image for the first segment in each pair. | `True` |
| `--opt-mode TEXT` | Single-structure optimizer for HEI±1/kink nodes. `light` maps to LBFGS; `heavy` maps to RFO. | `light` |
| `--mep-mode {gsm\|dmf}` | Segment generator: GSM (string-based) or DMF (direct flux). | `gsm` |
| `--refine-mode {peak\|minima}` | Seeds for refinement: `peak` optimizes HEI±1; `minima` searches outward from the HEI toward the nearest local minima on each side. Defaults to `peak` for GSM and `minima` for DMF when omitted. | _Auto_ |
| `--dump BOOL` | Explicit `True`/`False`. Dump GSM and single-structure trajectories/restarts. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB or Gaussian inputs. | `--convert-files` |
| `--out-dir TEXT` | Output directory. | `./result_path_search/` |
| `--thresh TEXT` | Override convergence preset for GSM and per-image optimizations (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |
| `--preopt BOOL` | Explicit `True`/`False`. Pre-optimise each endpoint before GSM (recommended). | `True` |
| `--align / --no-align` | Flag toggle. Align all inputs to the first structure before searching. | `--align` |
| `--ref-pdb PATH...` | Full-size template PDBs (one per input, unless `--align` lets you reuse the first). | _None_ |

## Workflow
1. **Initial GSM per pair** – run `GrowingString` between each adjacent input (A→B) to obtain a coarse MEP and identify the highest-energy image (HEI).
2. **Local relaxation around HEI** – refine either HEI ± 1 (`refine-mode=peak`) or the nearest local minima on each side of the HEI (`refine-mode=minima`) with the chosen single-structure optimizer (`opt-mode`) to recover nearby minima (`End1`, `End2`).
3. **Decide between kink vs. refinement**:
   - If no covalent bond change is detected between `End1` and `End2`, treat the region as a *kink*: insert `search.kink_max_nodes` linear nodes and optimize each individually.
   - Otherwise, launch a **refinement GSM** between `End1` and `End2` to sharpen the barrier.
4. **Selective recursion** – compare bond changes for `(A→End1)` and `(End2→B)` using the `bond` thresholds. Recurse only on sub-intervals that still contain covalent updates. Recursion depth is capped by `search.max_depth`.
5. **Stitching & bridging** – concatenate resolved subpaths, dropping duplicate endpoints when RMSD ≤ `search.stitch_rmsd_thresh`. If the RMSD gap between two stitched pieces exceeds `search.bridge_rmsd_thresh`, insert a bridge GSM. When the interface itself shows a bond change, a brand-new recursive segment replaces the bridge.
6. **Alignment & merging (optional)** – with `--align` (default), pre-optimized structures are rigidly aligned to the first input and `freeze_atoms` are reconciled. Provide `--ref-pdb` to merge pocket trajectories back into full-size PDB templates (one template per input unless alignment allows reuse of the first file).

Bond-change detection relies on `bond_changes.compare_structures` with thresholds surfaced under the `bond` YAML section. UMA calculators are constructed once and shared across all structures for efficiency.

## Outputs
```
out_dir/ (default: ./result_path_search/)
├─ mep.trj                  # Primary MEP trajectory
├─ mep.pdb                  # PDB companion when inputs were PDB templates and conversion is enabled
├─ mep_w_ref.pdb            # Merged full-system MEP (requires ref PDB/template)
├─ mep_w_ref_seg_XX.pdb     # Merged per-segment paths when covalent changes exist (requires ref PDB)
├─ summary.yaml             # Barrier and classification summary for every recursive segment
├─ mep_plot.png             # ΔE profile generated via `trj2fig` (kcal/mol, reactant reference)
├─ energy_diagram.html      # Plotly state-energy diagram (relative to reactant)
├─ energy_diagram.png       # Static export of the energy diagram
└─ segments/seg_000_*/      # GSM dumps, HEI snapshots, kink/refinement diagnostics per segment
```
- Console reports covering resolved configuration blocks (`geom`, `calc`, `gs`, `opt`, `sopt.*`, `bond`, `search`).

## Notes
- Provide at least two inputs; `click.BadParameter` is raised otherwise.
- `--ref-pdb` can be given once followed by multiple filenames; with `--align`, only the first template is reused for merges.
- All UMA calculators are shared across structures for efficiency.
- When `--dump` is set, GSM and single-structure optimizations emit trajectories and restart YAML files.
- Charge/spin inherit `.gjf` template metadata when available. If `-q` is omitted but `--ligand-charge` is provided, the inputs are treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge; explicit `-q` still overrides. Otherwise charge defaults to 0 and multiplicity to `1`.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. YAML parameters override the CLI values. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml): `geom`/`calc` mirror single-structure options (with `--freeze-links` augmenting `geom.freeze_atoms` for PDBs), and `opt` inherits the StringOptimizer knobs documented for `path_opt`.

`gs` (Growing String) inherits defaults from `pdb2reaction.path_opt.GS_KW` with overrides for `max_nodes` (internal nodes per segment), climb behavior (`climb`, `climb_rms`, `climb_fixed`), and reparameterization cadence (`reparam_every_full`, `reparam_check`).

`sopt` houses the single-structure optimizers used for HEI±1 and kink nodes, split into `lbfgs` and `rfo` subsections. Each subsection mirrors [`opt`](opt.md#yaml-configuration-args-yaml) but defaults to `out_dir: ./result_path_search/` and `dump: False`.

`bond` carries the UMA-based bond-change detection parameters shared with [`scan`](scan.md#section-bond): `device`, `bond_factor`, `margin_fraction`, and `delta_fraction`.

`search` governs the recursion logic: `max_depth`, `stitch_rmsd_thresh`, `bridge_rmsd_thresh`, `max_nodes_segment`, `max_nodes_bridge`, `kink_max_nodes`, `max_seq_kink`, and `refine_mode` (auto-selects `peak` for GSM or `minima` for DMF when left `null`). The legacy `rmsd_align` flag is ignored but kept for compatibility.

`dmf` bundles Direct Max Flux + (C)FB-ENM controls applied whenever `--mep-mode dmf` is selected. The defaults mirror the shared `DMF_KW` dictionary and can be overridden per run:

### Example YAML (default value)
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
gs:
  fix_first: true            # keep the first endpoint fixed during optimization
  fix_last: true             # keep the last endpoint fixed during optimization
  max_nodes: 10              # maximum string nodes
  perp_thresh: 0.005         # perpendicular displacement threshold
  reparam_check: rms         # reparametrization check metric
  reparam_every: 1           # reparametrization stride
  reparam_every_full: 1      # full reparametrization stride
  param: equi                # parametrization scheme
  max_micro_cycles: 10       # micro-iteration limit
  reset_dlc: true            # rebuild delocalized coordinates each step
  climb: true                # enable climbing image
  climb_rms: 0.0005          # climbing RMS threshold
  climb_lanczos: true        # Lanczos refinement for climbing
  climb_lanczos_rms: 0.0005  # Lanczos RMS threshold
  climb_fixed: false         # keep climbing image fixed
  scheduler: null            # optional scheduler backend
opt:
  type: string               # optimizer type label
  stop_in_when_full: 300     # early stop threshold when string is full
  align: false               # alignment toggle (kept off)
  scale_step: global         # step scaling mode
  max_cycles: 300            # maximum optimization cycles
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  reparam_thresh: 0.0        # reparametrization threshold
  coord_diff_thresh: 0.0     # coordinate difference threshold
  out_dir: ./result_path_search/   # output directory
  print_every: 10            # logging stride
dmf:
  correlated: true           # correlated DMF propagation
  sequential: true           # sequential DMF execution
  fbenm_only_endpoints: false   # run FB-ENM beyond endpoints
  fbenm_options:
    delta_scale: 0.2         # FB-ENM displacement scaling
    bond_scale: 1.25         # bond cutoff scaling
    fix_planes: true         # enforce planar constraints
    two_hop_mode: sparse     # neighbor traversal strategy
  cfbenm_options:
    bond_scale: 1.25         # CFB-ENM bond cutoff scaling
    corr0_scale: 1.1         # correlation scale for corr0
    corr1_scale: 1.5         # correlation scale for corr1
    corr2_scale: 1.6         # correlation scale for corr2
    eps: 0.05                # correlation epsilon
    pivotal: true            # pivotal residue handling
    single: true             # single-atom pivots
    remove_fourmembered: true   # prune four-membered rings
    two_hop_mode: sparse     # neighbor traversal strategy
  dmf_options:
    remove_rotation_and_translation: false  # keep rigid-body motions
    mass_weighted: false     # toggle mass weighting
    parallel: false          # enable parallel DMF
    eps_vel: 0.01            # velocity tolerance
    eps_rot: 0.01            # rotational tolerance
    beta: 10.0               # beta parameter for DMF
    update_teval: false      # update transition evaluation
  k_fix: 100.0               # harmonic constant for restraints
sopt:
  lbfgs:
    thresh: gau                # LBFGS convergence preset
    max_cycles: 10000          # iteration limit
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_path_search/   # output directory
    keep_last: 7               # history size for LBFGS buffers
    beta: 1.0                  # initial damping beta
    gamma_mult: false          # multiplicative gamma update toggle
    max_step: 0.3              # maximum step length
    control_step: true         # control step length adaptively
    double_damp: true          # double damping safeguard
    mu_reg: null               # regularization strength
    max_mu_reg_adaptions: 10   # cap on mu adaptations
  rfo:
    thresh: gau                # RFOptimizer convergence preset
    max_cycles: 10000          # iteration cap
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_path_search/   # output directory
    trust_radius: 0.3          # trust-region radius
    trust_update: true         # enable trust-region updates
    trust_min: 0.01            # minimum trust radius
    trust_max: 0.3             # maximum trust radius
    max_energy_incr: null      # allowed energy increase per step
    hessian_update: bfgs       # Hessian update scheme
    hessian_init: calc         # Hessian initialization source
    hessian_recalc: 100        # rebuild Hessian every N steps
    hessian_recalc_adapt: 2.0  # adaptive Hessian rebuild factor
    small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
    alpha0: 1.0                # initial micro step
    max_micro_cycles: 25       # micro-iteration limit
    rfo_overlaps: false        # enable RFO overlaps
    gediis: false              # enable GEDIIS
    gdiis: true                # enable GDIIS
    gdiis_thresh: 0.0025       # GDIIS acceptance threshold
    gediis_thresh: 0.01        # GEDIIS acceptance threshold
    gdiis_test_direction: true # test descent direction before DIIS
    adapt_step_func: false     # adaptive step scaling toggle
bond:
  device: cuda                # UMA device for bond analysis
  bond_factor: 1.2            # covalent-radius scaling
  margin_fraction: 0.05       # tolerance margin for comparisons
  delta_fraction: 0.05        # minimum relative change to flag bonds
search:
  max_depth: 10               # recursion depth limit
  stitch_rmsd_thresh: 0.0001  # RMSD threshold for stitching segments
  bridge_rmsd_thresh: 0.0001  # RMSD threshold for bridging nodes
  rmsd_align: true            # legacy alignment flag (ignored)
  max_nodes_segment: 10       # max nodes per segment
  max_nodes_bridge: 5         # max nodes per bridge
  kink_max_nodes: 3           # max nodes for kink optimizations
  max_seq_kink: 2             # max sequential kinks
  refine_mode: null           # optional refinement strategy (auto-chooses when null)
```
