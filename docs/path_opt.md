# `path-opt` subcommand

## Overview
`pdb2reaction path-opt` searches for a minimum-energy path (MEP) between two endpoint structures using pysisyphus' Growing String method (GSM). UMA supplies energies/gradients/Hessians for every image, while an external rigid-body alignment routine keeps the string well behaved before the optimizer begins. Configuration follows the precedence **CLI > `--args-yaml` > defaults** across the `geom`, `calc`, `gs`, and `opt` sections.

## Usage
```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} -q CHARGE -m MULT \
                      [--freeze-links BOOL] [--max-nodes N] [--max-cycles N] \
                      [--climb BOOL] [--dump BOOL] [--thresh PRESET] \
                      [--out-dir DIR] [--args-yaml FILE]
```

## Workflow
1. **Pre-alignment & freeze resolution**
   - All endpoints after the first are Kabsch-aligned to the first structure. If either endpoint defines `freeze_atoms`, only those atoms participate in the RMSD fit and the resulting transform is applied to every atom.
   - For PDB inputs with `--freeze-links=True` (default), parent atoms of link hydrogens are detected and merged into `freeze_atoms`.
   - `StringOptimizer.align` remains disabled because the explicit alignment/refinement already handles the superposition.
2. **String growth and HEI export**
   - After the path is grown and refined, the tool searches for the highest-energy internal local maximum (preferred). If none exists, it falls back to the maximum among internal nodes; if no internal nodes are present, the global maximum is exported.
   - The highest-energy image (HEI) is written both as `.xyz` and `.pdb` when a PDB reference exists, and as `.gjf` when a Gaussian template is available.

### Key behaviours
- **Endpoints**: Exactly two structures are required. Formats follow `geom_loader`. PDB inputs also enable trajectory/HEI PDB exports.
- **Charge/spin**: CLI overrides `.gjf` template metadata; otherwise defaults to `0/1`. Always set them explicitly for correct states.
- **Growing string**: `--max-nodes` controls the number of *internal* nodes (total images = `max_nodes + 2`). GSM growth and the optional climbing-image refinement share a convergence threshold preset supplied via `--thresh` or YAML (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`).
- **Climbing image**: `--climb` toggles both the standard climbing step and the Lanczos-based tangent refinement.
- **Dumping**: `--dump True` mirrors `opt.dump=True` for the StringOptimizer, producing trajectory/restart dumps inside `out_dir`.
- **Exit codes**: `0` success, `3` optimizer failure, `4` trajectory write error, `5` HEI export error, `130` interrupt, `1` unexpected error.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product structures (`.pdb`/`.xyz`). | Required |
| `-q, --charge INT` | Total charge (`calc.charge`). Required unless the input is a `.gjf` template that already encodes charge. | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity (`calc.spin`). | Template/`1` |
| `--freeze-links BOOL` | PDB-only: freeze link-H parents (merged with YAML). | `True` |
| `--max-nodes INT` | Number of internal nodes (string images = `max_nodes + 2`). | `10` |
| `--max-cycles INT` | Optimizer macro-iteration cap (`opt.max_cycles`). | `300` |
| `--climb BOOL` | Enable climbing-image refinement (and Lanczos tangent). | `True` |
| `--dump BOOL` | Dump GSM trajectories/restarts. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--thresh TEXT` | Override convergence preset for GSM/string optimizer. | YAML/default |
| `--args-yaml FILE` | YAML overrides (sections `geom`, `calc`, `gs`, `opt`). | _None_ |
| `--preopt BOOL` | Pre-optimise each endpoint with the selected single-structure optimizer before alignment/GSM. | `False` |
| `--preopt-max-cycles INT` | Cap for endpoint preoptimization cycles. | `10000` |
| `--fix-ends BOOL` | Keep the endpoint geometries fixed during GSM growth/refinement. | `False` |

## Outputs
```
out_dir/
├─ final_geometries.trj        # XYZ path; comment line holds energies when provided
├─ final_geometries.pdb        # Only when the first endpoint was a PDB
├─ gsm_hei.xyz                 # Highest-energy image with its energy on the comment line
├─ gsm_hei.pdb                 # HEI converted to PDB when the first endpoint was a PDB
├─ gsm_hei.gjf                 # HEI written using a detected Gaussian template
├─ align_refine/               # Intermediate files from the rigid alignment/refinement stage (created only when scan output exists)
└─ <optimizer dumps/restarts>  # Present when dumping is enabled
```
Console output echoes the resolved YAML blocks and prints cycle-by-cycle GSM progress with timing information.

## YAML configuration (`--args-yaml`)
CLI inputs override YAML, which override the defaults listed below.

### `geom`
- Same keys as [`opt`](opt.md) (`coord_type`, `freeze_atoms`, etc.); `--freeze-links` augments `freeze_atoms` for PDBs.

### `calc`
- UMA calculator setup identical to the single-structure optimization (`model`, `device`, neighbour radii, Hessian options, etc.).

### `gs`
- Controls the Growing String representation: `max_nodes`, `perp_thresh`, reparametrisation cadence (`reparam_check`, `reparam_every`, `reparam_every_full`, `param`), `max_micro_cycles`, DLC resets, climb toggles/thresholds, and optional scheduler hooks.

### `opt`
- StringOptimizer settings: type labels, `stop_in_when_full`, `align=False` (kept off), `scale_step`, `max_cycles`, dumping flags, `reparam_thresh`, `coord_diff_thresh`, `out_dir`, and `print_every`.

### Example YAML
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
  max_nodes: 30
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
  out_dir: ./result_path_opt/
  print_every: 5
```
