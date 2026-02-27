# `path-opt`

## Overview

> **Summary:** Find an MEP between **exactly two** structures with GSM (default) or DMF (`--mep-mode dmf`). Writes the path trajectory and exports the highest-energy image (HEI) as a TS candidate.

### At a glance
- **Use when:** You have reactant and product endpoints (R → P) and want a first-pass MEP.
- **Method:** GSM by default; switch to DMF with `--mep-mode dmf`.
- **Outputs:** `final_geometries_trj.xyz` (path) and `hei.xyz` (HEI), plus optional `.pdb`/`.gjf` companions when conversion is enabled.
- **Defaults:** `--opt-mode grad` (LBFGS), `--climb`, `--max-nodes 10`, `--thresh gau`.
- **Next step:** Validate the HEI with `tsopt` → `freq` (expect **one** imaginary mode) → `irc`.

`pdb2reaction path-opt` searches for a minimum-energy path (MEP) between two endpoints and reports the highest-energy image (HEI). Treat the HEI as a *candidate* transition state until it is validated with [freq](freq.md) and [irc](irc.md). For workflows that start from **two or more** structures and automatically refine only the reactive region, use [path-search](path_search.md).

UMA provides energies, gradients, and Hessians for every image. Before optimization starts, a rigid-body alignment step keeps the string well-behaved; if you define `freeze_atoms`, only those atoms are used for the RMSD fit (the transform is still applied to all atoms).


## Minimal example

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_opt
```

## Output checklist

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb` (when PDB conversion is available)

## Common examples

1. Pre-optimize endpoints before MEP search.

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

2. Use DMF mode instead of GSM.

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

3. Freeze link parents and disable climb for a quick pass.

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --freeze-links --no-climb --out-dir ./result_path_opt_fast
```

## Usage
```bash
pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [--workers N] [--workers-per-node N] \
 [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--max-nodes N] [--max-cycles N] \
 [--climb/--no-climb] [--dump/--no-dump] [--thresh PRESET] \
 [--preopt/--no-preopt] [--preopt-max-cycles N] [--opt-mode grad|hess] [--fix-ends/--no-fix-ends] \
 [--show-config/--no-show-config] [--dry-run/--no-dry-run] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

## Workflow
1. **Pre-alignment & freeze resolution**
 - All endpoints after the first are Kabsch-aligned to the first structure. If either endpoint defines `freeze_atoms`, only those atoms participate in the RMSD fit and the resulting transform is applied to every atom.
 - For PDB inputs with `--freeze-links=True` (default), parent atoms of link hydrogens are detected and merged into `freeze_atoms`.
2. **String growth and HEI export**
 - After the path is grown and refined, the tool searches for the highest-energy internal local maximum (preferred). If none exists, it falls back to the maximum among internal nodes; if no internal nodes are present, the global maximum is exported.
 - The highest-energy image (HEI) is written both as `.xyz` and `.pdb` when a PDB reference exists, and as `.gjf` when a Gaussian template is available; these conversions honor `--convert-files`.

### Key behaviors
- **Endpoints**: Exactly two structures are required. Formats follow `geom_loader`. PDB inputs (or XYZ/GJF with `--ref-pdb`) enable trajectory/HEI PDB exports.
- **Charge/spin**: CLI overrides `.gjf` template metadata. If `-q` is omitted but `--ligand-charge` is provided, the endpoints are treated as an enzyme–substrate complex and `extract.py`’s charge summary computes the total charge for PDB inputs (or XYZ/GJF when `--ref-pdb` is supplied); explicit `-q` still overrides. For non-`.gjf` inputs, omitting `-q` aborts unless derivation succeeds. If `.gjf` inputs lack charge metadata and `-q` is not provided, the command aborts; multiplicity defaults to `1` when omitted. Always set them explicitly for correct states.
- **MEP segments**: `--max-nodes` controls the number of *internal* nodes/images for the GSM string or DMF path (total images = `max_nodes + 2` for GSM). GSM growth and the optional climbing-image refinement share a convergence threshold preset supplied via `--thresh` or YAML (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`).
- **Climbing image**: `--climb` toggles both the standard climbing step and the Lanczos-based tangent refinement.
- **Dumping**: `--dump` mirrors `opt.dump=True` for the StringOptimizer, producing trajectory dumps inside `out_dir`. Restart YAML is written only when enabled in YAML.
- **Exit codes**: `0` success, `3` optimizer failure, `4` trajectory write error, `5` HEI export error, `130` interrupt, `1` unexpected error.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product structures (`.pdb`/`.xyz`). | Required |
| `-q, --charge INT` | Total charge (`calc.charge`). Required for non-`.gjf` inputs unless `--ligand-charge` derivation succeeds (PDB inputs or XYZ/GJF with `--ref-pdb`). `.gjf` templates can supply it; if `.gjf` inputs lack charge metadata, the run aborts unless `-q` is provided. Overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex for PDB inputs (or XYZ/GJF when `--ref-pdb` is supplied). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (`calc.spin`). | Template/`1` |
| `--freeze-links/--no-freeze-links` | PDB-only: freeze link-H parents (merged with YAML). | `True` |
| `--max-nodes INT` | Number of internal nodes (string images = `max_nodes + 2`). | `10` |
| `--mep-mode {gsm\|dmf}` | Select GSM (string-based) or DMF (direct flux) path generator. | `gsm` |
| `--max-cycles INT` | Optimizer macro-iteration cap (`opt.max_cycles`). | `300` |
| `--climb/--no-climb` | Enable climbing-image refinement (and Lanczos tangent). | `True` |
| `--dump/--no-dump` | Dump MEP trajectories (GSM/DMF). Restart YAML is written only when enabled in YAML. | `False` |
| `--opt-mode TEXT` | Single-structure optimizer for endpoint preoptimization (`grad` = LBFGS, `hess` = RFO). | `grad` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology for XYZ/GJF inputs (keeps XYZ coordinates) to enable PDB conversions. | _None_ |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--thresh TEXT` | Override convergence preset for GSM/string optimizer. | `gau` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layers) and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running optimization. | `False` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint with the selected single-structure optimizer before alignment/MEP search (GSM/DMF). | `False` |
| `--preopt-max-cycles INT` | Cap for endpoint preoptimization cycles. | `10000` |
| `--fix-ends/--no-fix-ends` | Keep the endpoint geometries fixed during GSM growth/refinement. | `False` |

## Outputs
```
out_dir/
├─ final_geometries_trj.xyz # XYZ path; comment line holds energies when provided
├─ final_geometries.pdb # When a PDB reference is available (input PDB or --ref-pdb) and conversion enabled
├─ hei.xyz # Highest-energy image with its energy on the comment line
├─ hei.pdb # HEI converted to PDB when a PDB reference is available (conversion enabled)
├─ hei.gjf # HEI written using a detected Gaussian template (conversion enabled)
├─ align_refine/ # Intermediate files from the rigid alignment/refinement stage (created when alignment runs)
└─ <optimizer dumps> # Trajectory dumps when --dump (restart YAML only via YAML dump_restart)
```
Console output echoes the resolved YAML blocks and prints cycle-by-cycle MEP progress (GSM/DMF) with timing information.

Merge order is **defaults < config < explicit CLI < override**.
### `geom`
- Same keys as [`opt`](opt.md) (`coord_type`, `freeze_atoms`, etc.); `--freeze-links` augments `freeze_atoms` for PDBs.

### `calc`
- UMA calculator setup identical to the single-structure optimization (`model`, `device`, neighbor radii, Hessian options, etc.).

### `dmf`
- Direct Max Flux + (C)FB-ENM interpolation controls. Keys mirror the CLI-accessible `dmf` block:

### `gs`
- Controls the Growing String representation: `max_nodes`, `perp_thresh`, reparameterization cadence (`reparam_check`, `reparam_every`, `reparam_every_full`, `param`), `max_micro_cycles`, DLC resets, climb toggles/thresholds, and optional scheduler hooks.

### `opt`
- StringOptimizer settings: type labels, `stop_in_when_full`, `scale_step`, `max_cycles`, dumping flags, `reparam_thresh`, `coord_diff_thresh`, `out_dir`, and `print_every`.

### `sopt.lbfgs` / `sopt.rfo`
- Single-structure preoptimization settings for endpoints. Keys mirror the `lbfgs`/`rfo` sections in [YAML Reference](yaml_reference.md). YAML overrides CLI `--preopt-max-cycles`.

### Example YAML (default value)
```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # UMA model tag
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: false # full Hessian (avoids shape mismatches)
gs:
 fix_first: true # keep the first endpoint fixed during optimization
 fix_last: true # keep the last endpoint fixed during optimization
 max_nodes: 10 # maximum string nodes
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
opt:
 type: string # optimizer type label
 stop_in_when_full: 300 # early stop threshold when string is full
 scale_step: global # step scaling mode
 max_cycles: 300 # maximum optimization cycles
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 reparam_thresh: 0.0 # reparametrization threshold
 coord_diff_thresh: 0.0 # coordinate difference threshold
 out_dir:./result_path_opt/ # output directory
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
 two_hop_mode: sparse # neighbor traversal strategy
 cfbenm_options:
 bond_scale: 1.25 # CFB-ENM bond cutoff scaling
 corr0_scale: 1.1 # correlation scale for corr0
 corr1_scale: 1.5 # correlation scale for corr1
 corr2_scale: 1.6 # correlation scale for corr2
 eps: 0.05 # correlation epsilon
 pivotal: true # pivotal residue handling
 single: true # single-atom pivots
 remove_fourmembered: true # prune four-membered rings
 two_hop_mode: sparse # neighbor traversal strategy
 dmf_options:
 remove_rotation_and_translation: false # keep rigid-body motions
 mass_weighted: false # toggle mass weighting
 parallel: false # enable parallel DMF
 eps_vel: 0.01 # velocity tolerance
 eps_rot: 0.01 # rotational tolerance
 beta: 10.0 # beta parameter for DMF
 update_teval: false # update transition evaluation
 k_fix: 300.0 # harmonic constant for restraints
```

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [path-search](path_search.md) — Recursive MEP search with automatic refinement (for 2+ structures)
- [tsopt](tsopt.md) — Optimize the HEI as a TS candidate (validate with freq/IRC)
- [extract](extract.md) — Generate pocket PDBs for path-opt inputs
- [all](all.md) — End-to-end workflow (uses path-search by default)
- [YAML Reference](yaml_reference.md) — Full `gs`, `dmf`, `opt` configuration options
- [Glossary](glossary.md) — Definitions of MEP, GSM, DMF, HEI
