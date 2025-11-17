# `all` subcommand

## Purpose
Run an end-to-end enzymatic reaction workflow on pocket models: extract pockets, optionally perform a staged scan for a single structure, run the recursive GSM minimum-energy-path search, merge the pocket path back into the original full systems, and (optionally) execute TS optimisation, thermochemistry, and DFT post-processing per reactive segment. When exactly one structure is supplied without `--scan-lists`, enabling `--tsopt True` triggers a TSOPT-only pocket workflow (no path search).

## Usage
```bash
# Multi-structure ensemble (reaction order)
pdb2reaction all -i R.pdb [I.pdb ...] P.pdb -c SUBSTRATE_SPEC \
                 [--ligand-charge MAP_OR_NUMBER] [--spin 2S+1] \
                 [--freeze-links True|False] [--max-nodes N] [--max-cycles N] \
                 [--climb True|False] [--sopt-mode lbfgs|rfo|light|heavy] \
                 [--dump True|False] [--preopt True|False] \
                 [--args-yaml FILE] [--out-dir DIR] \
                 [--tsopt True|False] [--thermo True|False] [--dft True|False]

# Single-structure + staged scan (pocket scan results become intermediates)
pdb2reaction all -i SINGLE.pdb -c SUBSTRATE_SPEC \
                 --scan-lists "[(i,j,target_Å), ...]" [--scan-lists "..."] \
                 [other options as above]

# Single-structure TSOPT-only mode (no path_search)
pdb2reaction all -i SINGLE.pdb -c SUBSTRATE_SPEC --tsopt True [other toggles]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full PDBs in reaction order, or a single PDB when used with `--scan-lists` or `--tsopt True`. A single `-i` may be followed by multiple files. | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs like `123,124` / `A:123,B:456`, or residue names like `GPP,MMT`). | Required |
| `--out-dir PATH` | Top-level output directory. | `./result_all/` |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). | `0.0` |
| `--include-H2O BOOLEAN` | Include waters (set `False` to drop HOH/WAT/TIP3/SOL). | `True` |
| `--exclude-backbone BOOLEAN` | Remove backbone atoms on non-substrate amino acids. | `True` |
| `--add-linkH BOOLEAN` | Add link hydrogens for severed bonds (carbon-only). | `True` |
| `--selected_resn TEXT` | Residues to force include (comma/space separated; chain/insertion codes allowed). | `""` |
| `--ligand-charge TEXT` | Total charge or mapping for unknown residues (recommended). | `None` |
| `--verbose BOOLEAN` | Enable INFO-level logging in the extractor. | `True` |
| `-s, --spin INT` | Spin multiplicity forwarded to GSM, scan, and post-processing. | `1` |
| `--freeze-links BOOLEAN` | Freeze link parents in pocket PDBs during GSM/scan (also reused by scan/tsopt/freq). | `True` |
| `--max-nodes INT` | GSM internal nodes per segment. | `10` |
| `--max-cycles INT` | GSM maximum optimisation cycles. | `100` |
| `--climb BOOLEAN` | Enable transition-state climbing for the first segment in each pair. | `True` |
| `--sopt-mode [lbfgs|rfo|light|heavy]` | Single-structure optimiser for HEI±1/kink nodes. | `lbfgs` |
| `--opt-mode [light|lbfgs|heavy|rfo]` | UMA optimiser preset for scan/tsopt. Use `light`/`lbfgs` (LBFGS) or `heavy`/`rfo` (RFO). When omitted, scan inherits LBFGS or RFO from `--sopt-mode`, and tsopt falls back to `light`. | `None` |
| `--dump BOOLEAN` | Dump GSM and single-structure trajectories, propagating the same setting to scan/tsopt/freq. | `False` |
| `--args-yaml FILE` | YAML forwarded unchanged to `path_search`, `scan`, `tsopt`, `freq`, and `dft`. | _None_ |
| `--preopt BOOLEAN` | preoptimise pocket endpoints before GSM (also used as the default for scan preoptimisation). | `True` |
| `--hessian-calc-mode [Analytical|FiniteDifference]` | Shared UMA Hessian engine forwarded to tsopt and freq. | _None_ |
| `--tsopt BOOLEAN` | Run TS optimisation + pseudo-IRC per reactive segment, or enable TSOPT-only mode (single input, no scan). | `False` |
| `--thermo BOOLEAN` | Run vibrational analysis (freq) on R/TS/P and build UMA Gibbs diagram. | `False` |
| `--dft BOOLEAN` | Run single-point DFT on R/TS/P and build DFT energy diagram (adds DFT//UMA when `--thermo True`). | `False` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles` for every TS refinement. | _None_ |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory (relative paths are resolved against the default). | _None_ |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Override `freq --max-write`. | _None_ |
| `--freq-amplitude-ang FLOAT` | Override `freq --amplitude-ang` (Å). | _None_ |
| `--freq-n-frames INT` | Override `freq --n-frames`. | _None_ |
| `--freq-sort [value|abs]` | Override freq mode sorting behaviour. | _None_ |
| `--freq-temperature FLOAT` | Override freq thermochemistry temperature (K). | _None_ |
| `--freq-pressure FLOAT` | Override freq thermochemistry pressure (atm). | _None_ |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Override `dft --func-basis`. | _None_ |
| `--dft-max-cycle INT` | Override `dft --max-cycle`. | _None_ |
| `--dft-conv-tol FLOAT` | Override `dft --conv-tol`. | _None_ |
| `--dft-grid-level INT` | Override `dft --grid-level`. | _None_ |
| `--scan-lists TEXT...` | One or more Python-like lists describing staged scans on the extracted pocket (single-input runs only). Each list element is `(i,j,target\_Å)` (values are parsed from a Python-like literal). Atom indices come from the original `all` input PDB (1-based) and are remapped internally to the pocket ordering. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory (`<out-dir>/scan` by default). | _None_ |
| `--scan-one-based BOOLEAN` | Force 1-based (`True`) or 0-based (`False`) scan indexing; `None` keeps the scan default (1-based). | _None_ |
| `--scan-max-step-size FLOAT` | Override scan `--max-step-size` (Å). | _None_ |
| `--scan-bias-k FLOAT` | Override the harmonic bias strength `k` (eV/Å²). | _None_ |
| `--scan-relax-max-cycles INT` | Override scan relaxation max cycles per step. | _None_ |
| `--scan-preopt BOOLEAN` | Override the scan preoptimisation toggle (otherwise `--preopt` is reused). | _None_ |
| `--scan-endopt BOOLEAN` | Override the scan end-of-stage optimisation toggle. | _None_ |

## Outputs
- `<out-dir>/pockets/`: Pocket PDBs for each input.
- `<out-dir>/scan/`: Present when `--scan-lists` is used; contains staged pocket scan results (`stage_XX/result.pdb`).
- `<out-dir>/path_search/`: GSM results (trajectory, merged full-system PDBs, energy diagrams, `summary.yaml`, per-segment folders).
- `<out-dir>/path_search/tsopt_seg_XX/`: Present when post-processing is enabled; includes TS optimisation, pseudo-IRC, frequency, and DFT outputs plus diagrams.
- `<out-dir>/tsopt_single/`: Present only in TSOPT-only mode (single structure, `--tsopt True`, no scan); contains TS optimisation outputs and diagrams.
- Console logs covering pocket charge summary, resolved configuration blocks, scan stages, and per-stage timing.

## Notes
- The total pocket charge from the extractor (first model) is rounded to the nearest integer and propagated to scan/GSM/TSOPT (a console note is emitted if rounding occurs).
- Reference PDB templates for merging are taken automatically from the original inputs; the explicit `--ref-pdb` option of `path_search` is intentionally hidden in this wrapper.
- When both `--thermo` and `--dft` are enabled, the post-processing stage also produces a DFT//UMA Gibbs diagram (DFT energy + UMA thermal correction).
- Always provide `--ligand-charge` when formal charges are not inferable to ensure the correct total charge propagates through the pipeline.
- Single-input runs require either `--scan-lists` (staged scan feeding GSM) or `--tsopt True` (TSOPT-only mode). All boolean toggles expect an explicit `True`/`False` value.

## YAML configuration (`--args-yaml`)
The same YAML file is forwarded unchanged to **every** invoked subcommand. Each tool reads the sections described in its own documentation:

- [`path_search`](path_search.md#yaml-configuration-args-yaml): `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search`.
- [`scan`](scan.md#yaml-configuration-args-yaml): `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond`.
- [`tsopt`](tsopt.md#yaml-configuration-args-yaml): `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo`.
- [`freq`](freq.md#yaml-configuration-args-yaml): `geom`, `calc`, `freq`, `thermo`.
- [`dft`](dft.md#yaml-configuration-args-yaml): `dft`.

Include whichever sections you need at the YAML root; overlapping names such as `geom`, `calc`, or `opt` will be shared across modules, while module-specific blocks (for example `freq` or `dft`) apply only where supported. CLI values always take precedence over the YAML contents.

Example snippet combining shared and module-specific sections:

```yaml
geom:
  coord_type: cart
calc:
  model: uma-s-1p1
  hessian_calc_mode: Analytical
gs:
  max_nodes: 12
freq:
  max_write: 8
dft:
  grid_level: 6
```

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