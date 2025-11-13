# `path_search` subcommand

## Purpose
Runs a recursive Growing String (GSM) search across multiple structures (reactant → intermediates → product), stitches segment paths, and optionally merges pocket trajectories back into full PDB templates.

## Usage
```bash
pdb2reaction path_search -i R.pdb [I.pdb ...] P.pdb -q CHARGE [--spin 2S+1]
                         [--freeze-links/--no-freeze-links]
                         [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
                         [--sopt-mode lbfgs|rfo|light|heavy] [--dump/--no-dump]
                         [--out-dir DIR] [--pre-opt/--no-pre-opt]
                         [--align/--no-align] [--ref-pdb FILE ...]
                         [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more structures in reaction order (reactant → product). A single `-i` may be followed by multiple paths. | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-links / --no-freeze-links` | For PDB inputs, freeze link-hydrogen parents when building pockets. | `--freeze-links` |
| `--max-nodes INT` | Internal nodes for GSM segments (`String` has `max_nodes + 2` images). | `10` |
| `--max-cycles INT` | Maximum GSM optimization cycles. | `100` |
| `--climb / --no-climb` | Enable climbing image for the first segment in each pair. | `--climb` |
| `--sopt-mode TEXT` | Single-structure optimizer for HEI±1/kink nodes (`light|lbfgs` or `heavy|rfo`). | `lbfgs` |
| `--dump / --no-dump` | Dump GSM and single-structure trajectories. | `--no-dump` |
| `--out-dir TEXT` | Output directory. | `./result_path_search/` |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |
| `--pre-opt / --no-pre-opt` | Pre-optimise each endpoint before the GSM search. | `--pre-opt` |
| `--align / --no-align` | Align all inputs to the first structure before searching. | `--align` |
| `--ref-pdb PATH...` | Full-size template PDBs (one per input, unless `--align` lets you reuse the first). | _None_ |

## YAML configuration (`--args-yaml`)
The YAML root must be a mapping. CLI parameters override YAML values. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

### Shared sections
- `geom`, `calc`: same keys as [`opt`](opt.md#yaml-configuration-args-yaml). `--freeze-links` augments `geom.freeze_atoms` when inputs are PDB.

### Section `gs`
Growing String controls for main segments. Defaults derive from `pdb2reaction.path_opt.GS_KW` with overrides:

- `max_nodes` (`10`): Internal nodes per GSM segment (CLI override).
- `reparam_every_full` (`1`), `climb_rms` (`5e-4`), `climb_fixed` (`False`): Growth/climb behaviour.
- Additional keys match [`path_opt`](path_opt.md#section-gs).

### Section `opt`
StringOptimizer controls for GSM (defaults in parentheses).

- `stop_in_when_full` (`1000`) and `max_cycles` (`1000`): Cycle caps (CLI overrides `--max-cycles`).
- `dump` (`False`), `dump_restart` (`False`), `out_dir` (`"./result_path_search/"`), `print_every` (`1`), `align` (`False`).
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
- `<out-dir>/mep.trj` (+ `.pdb` when pockets were PDB inputs).
- `<out-dir>/mep_w_ref.pdb` merged full-system trajectory (requires `--ref-pdb` or auto-supplied templates).
- `<out-dir>/summary.yaml` summarising segment barriers and classification.
- Per-segment folders containing GSM dumps, merged structures, and diagnostic energy plots.
- Console reports covering resolved configuration blocks (`geom`, `calc`, `gs`, `opt`, `sopt.*`, `bond`, `search`).

## Notes
- Provide at least two inputs; `click.BadParameter` is raised otherwise.
- `--ref-pdb` can be given once followed by multiple filenames; with `--align`, only the first template is reused for merges.
- All UMA calculators are shared across structures for efficiency.
- When `--dump` is set, GSM and single-structure optimizations emit trajectories and restart YAML files.
