# `path-search` subcommand

## Purpose
Runs a recursive Growing String (GSM) search across multiple structures (reactant → intermediates → product), stitches segment paths, and optionally merges pocket trajectories back into full PDB templates.

## Usage
```bash
pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb -q CHARGE [--spin 2S+1]
                         [--freeze-links BOOL] [--thresh PRESET]
                         [--max-nodes N] [--max-cycles N] [--climb BOOL]
                         [--sopt-mode lbfgs|rfo|light|heavy] [--dump BOOL]
                         [--out-dir DIR] [--pre-opt BOOL]
                         [--align/--no-align] [--ref-pdb FILE ...]
                         [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more structures in reaction order (reactant → product). A single `-i` may be followed by multiple paths. | Required |
| `-q, --charge INT` | Total charge. | `.gjf` template value or `0` |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | Explicit `True`/`False`. When loading PDB pockets, freeze the parent atoms of link hydrogens. | `True` |
| `--max-nodes INT` | Internal nodes for GSM segments (`String` has `max_nodes + 2` images). | `10` |
| `--max-cycles INT` | Maximum GSM optimization cycles. | `100` |
| `--climb BOOL` | Explicit `True`/`False`. Enable climbing image for the first segment in each pair. | `True` |
| `--sopt-mode TEXT` | Single-structure optimizer for HEI±1/kink nodes (`light|lbfgs` or `heavy|rfo`). | `lbfgs` |
| `--dump BOOL` | Explicit `True`/`False`. Dump GSM and single-structure trajectories. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_path_search/` |
| `--thresh TEXT` | Override convergence preset for GSM and per-image optimizations (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (use YAML/default) |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |
| `--pre-opt BOOL` | Explicit `True`/`False`. Pre-optimise each endpoint before the GSM search. | `True` |
| `--align / --no-align` | Flag toggle. Align all inputs to the first structure before searching. | `--align` |
| `--ref-pdb PATH...` | Full-size template PDBs (one per input, unless `--align` lets you reuse the first). | _None_ |

### Shared sections
- `geom`, `calc`: same keys as [`opt`](opt.md#yaml-configuration-args-yaml). `--freeze-links` augments `geom.freeze_atoms` when inputs are PDB.

### Section `gs`
Growing String controls for main segments. Defaults derive from `pdb2reaction.path_opt.GS_KW` with overrides:

- `max_nodes` (`10`): Internal nodes per GSM segment (CLI override).
- `reparam_every_full` (`1`), `climb_rms` (`5e-4`), `climb_fixed` (`False`): Growth/climb behaviour.
- Additional keys match [`path_opt`](path_opt.md#section-gs).

### Section `opt`
StringOptimizer controls for GSM (defaults in parentheses).

- `stop_in_when_full` (`100`) and `max_cycles` (`100`): Cycle caps (CLI overrides `--max-cycles`).
- `dump` (`False`), `dump_restart` (`False`), `out_dir` (`"./result_path_search/"`), `print_every` (`10`), `align` (`False`).
- Other keys mirror [`path_opt`](path_opt.md#section-opt).

### Section `sopt`
Single-structure optimization defaults for HEI±1 and kink nodes. Split into `lbfgs` and `rfo` subsections.

- `sopt.lbfgs`: Same keys as [`opt`](opt.md#section-lbfgs) but with defaults adapted for path search (`out_dir: ./result_path_search/`, `dump: False`).
- `sopt.rfo`: Same keys as [`opt`](opt.md#section-rfo) with analogous defaults.

### Section `bond`
Bond-change detection parameters (identical to [`scan`](scan.md#section-bond)).

- `device` (`"cuda"`), `bond_factor` (`1.20`), `margin_fraction` (`0.05`), `delta_fraction` (`0.05`).

### Section `search`
Recursive path-building controls.

- `max_depth` (`10`): Recursion depth for segment refinement.
- `stitch_rmsd_thresh` (`1.0e-4`): Maximum RMSD to consider endpoints duplicates during stitching.
- `bridge_rmsd_thresh` (`1.0e-4`): RMSD threshold to trigger insertion of bridge GSMs.
- `rmsd_align` (`True`): Retained for compatibility (ignored internally).
- `max_nodes_segment` (`10`): Max nodes per recursive segment (defaults to `--max-nodes` if unspecified).
- `max_nodes_bridge` (`5`): Nodes for bridge GSMs.
- `kink_max_nodes` (`3`): Linear interpolation nodes for skipped GSM at kinks.

## Outputs
- `<out-dir>/mep.trj` (and `.pdb` when the inputs were PDB pockets).
- `<out-dir>/mep_w_ref.pdb` merged full-system MEP (requires `--ref-pdb` or auto-provided templates).
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
- Charge/spin inherit `.gjf` template metadata when available; otherwise the CLI defaults to `0`/`1`. Override them explicitly
  when you need a different electronic state.

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. CLI parameters override YAML values. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

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