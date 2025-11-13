# `all` subcommand

## Purpose
Runs the full pipeline: pocket extraction across multiple full PDBs, minimum-energy-path search on the pockets, automatic merge back into the originals, and optional per-segment post-processing (TS optimisation, thermochemistry, DFT).

## Usage
```bash
pdb2reaction all -i R.pdb [I.pdb ...] P.pdb -c SUBSTRATE_SPEC
                 [--ligand-charge MAP_OR_NUMBER]
                 [--spin 2S+1]
                 [--freeze-links/--no-freeze-links]
                 [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
                 [--sopt-mode lbfgs|rfo|light|heavy]
                 [--dump/--no-dump] [--out-dir DIR]
                 [--pre-opt/--no-pre-opt]
                 [--args-yaml FILE]
                 [--tsopt/--no-tsopt] [--thermo/--no-thermo] [--dft/--no-dft]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full PDBs in reaction order. A single `-i` may be followed by multiple files. | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names) passed to the extractor. | Required |
| `--out-dir PATH` | Top-level output directory. | `./result_all/` |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). | `0.0` |
| `--include-H2O / --no-include-H2O` | Include waters. | `--include-H2O` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate residues. | `--exclude-backbone` |
| `--add-linkH / --no-add-linkH` | Add carbon-only link hydrogens. | `--add-linkH` |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--ligand-charge TEXT` | Total charge or mapping for unknown residues (recommended). | `None` |
| `--verbose / --no-verbose` | Extractor logging. | `--verbose` |
| `-s, --spin INT` | Spin multiplicity forwarded to path search and post-processing. | `1` |
| `--freeze-links / --no-freeze-links` | Freeze link parents in pocket PDBs. | `--freeze-links` |
| `--max-nodes INT` | GSM internal nodes per segment. | `10` |
| `--max-cycles INT` | GSM max cycles. | `100` |
| `--climb / --no-climb` | Enable climbing image for the first segment per pair. | `--climb` |
| `--sopt-mode TEXT` | Single-structure optimiser for HEI±1/kink nodes. | `lbfgs` |
| `--dump / --no-dump` | Dump GSM and single-structure trajectories. | `--no-dump` |
| `--args-yaml FILE` | YAML forwarded to `path_search` (see below). | _None_ |
| `--pre-opt / --no-pre-opt` | Pre-optimise pocket endpoints before GSM. | `--pre-opt` |
| `--tsopt / --no-tsopt` | Run TS optimisation and pseudo-IRC per reactive segment. | `--no-tsopt` |
| `--thermo / --no-thermo` | Run vibrational analysis (freq) on R/TS/P and build Gibbs diagram. | `--no-thermo` |
| `--dft / --no-dft` | Run single-point DFT on R/TS/P and build DFT energy diagram. | `--no-dft` |

## Outputs
- `<out-dir>/pockets/`: Pocket PDBs for each input.
- `<out-dir>/path_search/`: GSM results (trajectory, merged full-system PDBs, energy diagrams, `summary.yaml`, segment folders).
- Optional `<out-dir>/path_search/tsopt_seg_XX/` subtrees with TS optimisation, pseudo-IRC, frequency, and DFT results depending on toggles.
- Console logs covering pocket charge summary, resolved configuration blocks, and per-stage timing.

## Notes
- The total pocket charge from the extractor (first model) is rounded to the nearest integer and used as the GSM charge.
- Reference PDB templates for merging are taken automatically from the original inputs; the explicit `--ref-pdb` option of `path_search` is intentionally hidden in this wrapper.
- When both `--thermo` and `--dft` are enabled, the post-processing stage also produces a DFT//UMA Gibbs diagram (DFT energy + UMA thermal correction).
- Always provide `--ligand-charge` when formal charges are not inferable to ensure the correct total charge propagates through the pipeline.

## YAML configuration (`--args-yaml`)
The YAML file is forwarded unchanged to `path_search`. See [`path_search`](path_search.md#yaml-configuration-args-yaml) for accepted sections (`geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search`).

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
  stop_in_when_full: 1000
  align: false
  scale_step: global
  max_cycles: 1000
  dump: false
  dump_restart: false
  reparam_thresh: 0.001
  coord_diff_thresh: 0.0
  out_dir: ./result_path_search/
  print_every: 1
sopt:
  base:
    thresh: gau
    max_cycles: 10000
    print_every: 1
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
  lbfgs:
    thresh: gau
    max_cycles: 10000
    print_every: 1
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
    print_every: 1
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