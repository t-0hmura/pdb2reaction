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
| `-q, --charge INT` | Total charge. Required unless the first input is a `.gjf` template that already stores charge. | Required when not in template |
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
- `<out-dir>/mep.trj` (and `.pdb` companions when the inputs were PDB templates and conversion is enabled).
- `<out-dir>/mep_w_ref.pdb` merged full-system MEP (requires `--ref-pdb` or auto-provided templates; obeys conversion flag for generated companions).
- `<out-dir>/mep_w_ref_seg_XX.pdb` merged per-segment paths for segments with covalent changes (requires `--ref-pdb`).
- `<out-dir>/summary.yaml` summarising barriers and classification for every recursive segment.
- `<out-dir>/mep_plot.png` ΔE profile generated via `trj2fig` (`kcal/mol`, reference = reactant).
- `<out-dir>/energy_diagram.html` and `.png` Plotly diagrams of state energies (relative to reactant, kcal/mol).
- Per-segment folders (`segments/seg_000_*`) containing GSM dumps, HEI snapshots, merged HEI files, linear kink optimisations, and diagnostic energy plots.
- Console reports covering resolved configuration blocks (`geom`, `calc`, `gs`, `opt`, `sopt.*`, `bond`, `search`).

## Notes
- Provide at least two inputs; `click.BadParameter` is raised otherwise.
- `--ref-pdb` can be given once followed by multiple filenames; with `--align`, only the first template is reused for merges.
- All UMA calculators are shared across structures for efficiency.
- When `--dump` is set, GSM and single-structure optimizations emit trajectories and restart YAML files.
- Charge/spin inherit `.gjf` template metadata when available; otherwise `-q/--charge` is required and multiplicity defaults to
  `1`. Override them explicitly when you need a different electronic state.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. YAML parameters override the CLI values. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml): `geom`/`calc` mirror single-structure options (with `--freeze-links` augmenting `geom.freeze_atoms` for PDBs), and `opt` inherits the StringOptimizer knobs documented for `path_opt`.

`gs` (Growing String) inherits defaults from `pdb2reaction.path_opt.GS_KW` with overrides for `max_nodes` (internal nodes per segment), climb behavior (`climb`, `climb_rms`, `climb_fixed`), and reparameterization cadence (`reparam_every_full`, `reparam_check`).

`sopt` houses the single-structure optimizers used for HEI±1 and kink nodes, split into `lbfgs` and `rfo` subsections. Each subsection mirrors [`opt`](opt.md#yaml-configuration-args-yaml) but defaults to `out_dir: ./result_path_search/` and `dump: False`.

`bond` carries the UMA-based bond-change detection parameters shared with [`scan`](scan.md#section-bond): `device`, `bond_factor`, `margin_fraction`, and `delta_fraction`.

`search` governs the recursion logic: `max_depth`, `stitch_rmsd_thresh`, `bridge_rmsd_thresh`, `max_nodes_segment`, `max_nodes_bridge`, `kink_max_nodes`, and `refine_mode` (auto-selects `peak` for GSM or `minima` for DMF when left `null`). The legacy `rmsd_align` flag is ignored but kept for compatibility.

`dmf` bundles Direct Max Flux + (C)FB-ENM controls applied whenever `--mep-mode dmf` is selected. The defaults mirror the shared `DMF_KW` dictionary and can be overridden per run:

### Example YAML (default value)
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
gs:
  max_nodes: 10
  perp_thresh: 0.005
  reparam_check: rms
  reparam_every: 1
  reparam_every_full: 1
  param: equi
  max_micro_cycles: 10
  reset_dlc: true
  climb: true
  climb_rms: 0.0005
  climb_lanczos: true
  climb_lanczos_rms: 0.0005
  climb_fixed: false
  scheduler: null
opt:
  type: string
  stop_in_when_full: 100
  align: false
  scale_step: global
  max_cycles: 100
  dump: false
  dump_restart: false
  reparam_thresh: 0.001
  coord_diff_thresh: 0.0
  out_dir: ./result_path_search/
  print_every: 10
dmf:
  correlated: true
  sequential: true
  fbenm_only_endpoints: false
  fbenm_options:
    delta_scale: 0.2
    bond_scale: 1.25
    fix_planes: true
    two_hop_mode: sparse
  cfbenm_options:
    bond_scale: 1.25
    corr0_scale: 1.1
    corr1_scale: 1.5
    corr2_scale: 1.6
    eps: 0.05
    pivotal: true
    single: true
    remove_fourmembered: true
    two_hop_mode: dense
  dmf_options:
    remove_rotation_and_translation: false
    mass_weighted: false
    parallel: false
    eps_vel: 0.01
    eps_rot: 0.01
    beta: 10.0
    update_teval: false
  k_fix: 100.0
sopt:
  lbfgs:
    thresh: gau
    max_cycles: 10000
    print_every: 100
    min_step_norm: 1.0e-08
    assert_min_step: true
    rms_force: null
    rms_force_only: false
    max_force_only: false
    force_only: false
    converge_to_geom_rms_thresh: 0.05
    overachieve_factor: 0.0
    check_eigval_structure: false
    line_search: true
    dump: false
    dump_restart: false
    prefix: ""
    out_dir: ./result_path_search/
    keep_last: 7
    beta: 1.0
    gamma_mult: false
    max_step: 0.3
    control_step: true
    double_damp: true
    mu_reg: null
    max_mu_reg_adaptions: 10
  rfo:
    thresh: gau
    max_cycles: 10000
    print_every: 100
    min_step_norm: 1.0e-08
    assert_min_step: true
    rms_force: null
    rms_force_only: false
    max_force_only: false
    force_only: false
    converge_to_geom_rms_thresh: 0.05
    overachieve_factor: 0.0
    check_eigval_structure: false
    line_search: true
    dump: false
    dump_restart: false
    prefix: ""
    out_dir: ./result_path_search/
    trust_radius: 0.3
    trust_update: true
    trust_min: 0.01
    trust_max: 0.3
    max_energy_incr: null
    hessian_update: bfgs
    hessian_init: calc
    hessian_recalc: 100
    hessian_recalc_adapt: 2.0
    small_eigval_thresh: 1.0e-08
    alpha0: 1.0
    max_micro_cycles: 25
    rfo_overlaps: false
    gediis: false
    gdiis: true
    gdiis_thresh: 0.0025
    gediis_thresh: 0.01
    gdiis_test_direction: true
    adapt_step_func: false
bond:
  device: cuda
  bond_factor: 1.2
  margin_fraction: 0.05
  delta_fraction: 0.05
search:
  max_depth: 10
  stitch_rmsd_thresh: 0.0001
  bridge_rmsd_thresh: 0.0001
  rmsd_align: true
  max_nodes_segment: 10
  max_nodes_bridge: 5
  kink_max_nodes: 3
```
