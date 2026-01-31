# `path-opt`

## Overview

`pdb2reaction path-opt` finds the minimum-energy path (MEP) between exactly two structures using GSM (default) or DMF (`--mep-mode dmf`). It outputs the path trajectory and identifies the highest-energy image (HEI). For multi-structure workflows with automatic refinement, use `path-search` instead. UMA supplies energies/gradients/Hessians for every image, while an external rigid-body alignment routine keeps the string well-behaved before the optimizer begins. Configuration follows the precedence **defaults → CLI → `--args-yaml`** across the `geom`, `calc`, `gs`, and `opt` sections. When `--convert-files` is enabled (default), trajectories are mirrored to `.pdb` companions when PDB references exist, and XYZ snapshots (for example the HEI) are mirrored to `.gjf` companions when Gaussian templates exist. GSM is the default path generator. The default `--opt-mode` is **light** (LBFGS); use `--opt-mode heavy` for RFO.

## Usage
```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
                      [--workers N] [--workers-per-node N] \
                      [--mep-mode {gsm|dmf}] [--freeze-links {True\|False}] [--max-nodes N] [--max-cycles N] \
                      [--climb {True\|False}] [--dump {True\|False}] [--thresh PRESET] \
                      [--preopt {True\|False}] [--preopt-max-cycles N] [--opt-mode light|heavy] [--fix-ends {True\|False}] \
                      [--out-dir DIR] [--args-yaml FILE] \
                      [--convert-files {True\|False}] [--ref-pdb FILE]
```

## Workflow
1. **Pre-alignment & freeze resolution**
   - All endpoints after the first are Kabsch-aligned to the first structure. If either endpoint defines `freeze_atoms`, only those atoms participate in the RMSD fit and the resulting transform is applied to every atom.
   - For PDB inputs with `--freeze-links=True` (default), parent atoms of link hydrogens are detected and merged into `freeze_atoms`.
   - `StringOptimizer.align` remains disabled because the explicit alignment/refinement already handles the superposition.
2. **String growth and HEI export**
   - After the path is grown and refined, the tool searches for the highest-energy internal local maximum (preferred). If none exists, it falls back to the maximum among internal nodes; if no internal nodes are present, the global maximum is exported.
   - The highest-energy image (HEI) is written both as `.xyz` and `.pdb` when a PDB reference exists, and as `.gjf` when a Gaussian template is available; these conversions honor `--convert-files`.

### Key behaviors
- **Endpoints**: Exactly two structures are required. Formats follow `geom_loader`. PDB inputs (or XYZ/GJF with `--ref-pdb`) enable trajectory/HEI PDB exports.
- **Charge/spin**: CLI overrides `.gjf` template metadata. If `-q` is omitted but `--ligand-charge` is provided, the endpoints are treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge when the inputs are PDBs; explicit `-q` still overrides. When no template/derivation applies, the charge defaults to `0` (spin defaults to `1`). Always set them explicitly for correct states.
- **MEP segments**: `--max-nodes` controls the number of *internal* nodes/images for the GSM string or DMF path (total images = `max_nodes + 2` for GSM). GSM growth and the optional climbing-image refinement share a convergence threshold preset supplied via `--thresh` or YAML (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`).
- **Climbing image**: `--climb` toggles both the standard climbing step and the Lanczos-based tangent refinement.
- **Dumping**: `--dump True` mirrors `opt.dump=True` for the StringOptimizer, producing trajectory/restart dumps inside `out_dir`.
- **Exit codes**: `0` success, `3` optimizer failure, `4` trajectory write error, `5` HEI export error, `130` interrupt, `1` unexpected error.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product structures (`.pdb`/`.xyz`). | Required |
| `-q, --charge INT` | Total charge (`calc.charge`). Optional; when omitted, `.gjf` templates or `--ligand-charge` (PDB inputs) may supply it, otherwise it falls back to `0`. Overrides `--ligand-charge` when both are set. | Template/`0` |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex for PDB inputs; otherwise the charge falls back to `0`. | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (`calc.spin`). | Template/`1` |
| `--freeze-links {True\|False}` | PDB-only: freeze link-H parents (merged with YAML). | `True` |
| `--max-nodes INT` | Number of internal nodes (string images = `max_nodes + 2`). | `10` |
| `--mep-mode {gsm\|dmf}` | Select GSM (string-based) or DMF (direct flux) path generator. | `gsm` |
| `--max-cycles INT` | Optimizer macro-iteration cap (`opt.max_cycles`). | `300` |
| `--climb {True\|False}` | Enable climbing-image refinement (and Lanczos tangent). | `True` |
| `--dump {True\|False}` | Dump MEP trajectories/restarts (GSM/DMF). | `False` |
| `--opt-mode TEXT` | Single-structure optimizer for endpoint preoptimization (`light` = LBFGS, `heavy` = RFO). | `light` |
| `--convert-files {True\|False}` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology for XYZ/GJF inputs (keeps XYZ coordinates) to enable PDB conversions. | _None_ |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--thresh TEXT` | Override convergence preset for GSM/string optimizer. | `gau` |
| `--args-yaml FILE` | YAML overrides (sections `geom`, `calc`, `gs`, `opt`). | _None_ |
| `--preopt {True\|False}` | Pre-optimize each endpoint with the selected single-structure optimizer before alignment/MEP search (GSM/DMF). | `False` |
| `--preopt-max-cycles INT` | Cap for endpoint preoptimization cycles. | `10000` |
| `--fix-ends {True\|False}` | Keep the endpoint geometries fixed during GSM growth/refinement. | `False` |

## Outputs
```
out_dir/
├─ final_geometries.trj        # XYZ path; comment line holds energies when provided
├─ final_geometries.pdb        # When a PDB reference is available (input PDB or --ref-pdb) and conversion enabled
├─ hei.xyz                     # Highest-energy image with its energy on the comment line
├─ hei.pdb                     # HEI converted to PDB when a PDB reference is available (conversion enabled)
├─ hei.gjf                     # HEI written using a detected Gaussian template (conversion enabled)
├─ align_refine/               # Intermediate files from the rigid alignment/refinement stage (created when alignment runs)
└─ <optimizer dumps/restarts>  # Present when dumping is enabled
```
Console output echoes the resolved YAML blocks and prints cycle-by-cycle MEP progress (GSM/DMF) with timing information.

## YAML configuration (`--args-yaml`)
YAML inputs override CLI, which override the defaults listed below.

### `geom`
- Same keys as [`opt`](opt.md) (`coord_type`, `freeze_atoms`, etc.); `--freeze-links` augments `freeze_atoms` for PDBs.

### `calc`
- UMA calculator setup identical to the single-structure optimization (`model`, `device`, neighbor radii, Hessian options, etc.).

### `dmf`
- Direct Max Flux + (C)FB-ENM interpolation controls. Keys mirror the CLI-accessible `dmf` block:

### `gs`
- Controls the Growing String representation: `max_nodes`, `perp_thresh`, reparameterization cadence (`reparam_check`, `reparam_every`, `reparam_every_full`, `param`), `max_micro_cycles`, DLC resets, climb toggles/thresholds, and optional scheduler hooks.

### `opt`
- StringOptimizer settings: type labels, `stop_in_when_full`, `align=False` (kept off), `scale_step`, `max_cycles`, dumping flags, `reparam_thresh`, `coord_diff_thresh`, `out_dir`, and `print_every`.

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
  out_dir: ./result_path_opt/   # output directory
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
  k_fix: 300.0               # harmonic constant for restraints
```

---

## See Also

- [path-search](path_search.md) — Recursive MEP search with automatic refinement (for 2+ structures)
- [tsopt](tsopt.md) — Optimize the HEI as a transition state
- [extract](extract.md) — Generate pocket PDBs for path-opt inputs
- [all](all.md) — End-to-end workflow (uses path-search by default)
- [YAML Reference](yaml-reference.md) — Full `gs`, `dmf`, `opt` configuration options
- [Glossary](glossary.md) — Definitions of MEP, GSM, DMF, HEI