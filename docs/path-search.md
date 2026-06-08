# `path-search`

Build a continuous minimum-energy path (MEP) from **two or more** structures (R → … → P) via GSM (default, `--mep-mode gsm`, string-based) or DMF (`--mep-mode dmf`, direct flux). Use it when you need a single stitched MEP with automatic refinement. It refines only those regions where covalent bond changes are detected (`--refine-mode peak` optimizes HEI±1, `--refine-mode minima` searches outward toward the nearest local minima, defaulting to `peak` for GSM and `minima` for DMF), then stitches the resolved subpaths into a single trajectory and exports the highest-energy image (HEI) of each segment as a TS candidate (validate with tsopt + IRC). The recursive decomposition automatically detects multistep reactions and builds a detailed MEP for each elementary step. Complex multistep mechanisms may require manual trial-and-error—adjusting input intermediates, scan specifications, or convergence thresholds—to obtain a satisfactory pathway. If you only have **two** endpoints and do not need recursive refinement, [path-opt](path-opt.md) is the simpler option.

## Examples

Command form:

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

Two endpoints (reactant → product):

```bash
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_search
```

Provide explicit intermediates for a multistep path:

```bash
# Provide explicit intermediates for a multistep path
pdb2reaction path-search -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 -m 1 \
 --out-dir ./result_path_search_multi
```

Enable merged full-system outputs with template references:

```bash
# Enable merged full-system outputs with template references
pdb2reaction path-search -i R.pdb IM1.pdb P.pdb -q 0 -m 1 \
 --ref-full-pdb holo_template.pdb --out-dir ./result_path_search_merge
```

Use DMF mode with minima refinement:

```bash
# Use DMF mode with minima refinement
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --refine-mode minima --out-dir ./result_path_search_dmf
```

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

Bond-change detection relies on `bond_changes.compare_structures` with thresholds surfaced under the `bond` YAML section. All MLIP backends are constructed once and shared across structures for efficiency.

## Outputs

```text
out_dir/ (default:./result_path_search/)
├─ mep_trj.xyz # Primary MEP trajectory
├─ mep.pdb # PDB companion when inputs were PDB templates and conversion is enabled
├─ mep.gjf # Gaussian companion when a Gaussian template is detected
├─ mep_w_ref.pdb # Merged full-system MEP (requires ref PDB/template)
├─ mep_seg_XX_trj.xyz # Per-segment MEP trajectory (XYZ)
├─ mep_seg_XX.pdb # Per-segment PDB companion (when conversion is enabled)
├─ mep_seg_XX.gjf # Per-segment Gaussian companion (when a template is detected)
├─ mep_w_ref_seg_XX.pdb # Merged per-segment paths when covalent changes exist (requires ref PDB)
├─ hei_seg_XX.xyz # Per-segment highest-energy image
├─ hei_seg_XX.pdb # HEI PDB companion (when conversion is enabled)
├─ hei_seg_XX.gjf # HEI Gaussian companion (when a template is detected)
├─ hei_w_ref_seg_XX.pdb # Merged HEI in full-system context (requires ref PDB)
├─ summary.json # Barrier and classification summary for every recursive segment
├─ summary.log # Text summary
├─ mep_plot.png # ΔE profile generated via `trj2fig` (kcal/mol, reactant reference)
├─ energy_diagram_MEP.png # Static export of the MEP state-energy diagram (relative to reactant)
└─ seg_NNN_*/ # GSM/DMF dumps, HEI snapshots, kink/refinement diagnostics per segment
```

- Console reports covering resolved configuration blocks (`geom`, `calc`, `gs`, `stopt`, `opt.*`, `bond`, `search`); see {ref}`verbosity-levels`.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation — do not hand-duplicate it here.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more structures in reaction order (reactant → product). Pass all files after a single `-i`/`--input`. | Required |
| `-q, --charge INT` | Net charge. Required for non-`.gjf` inputs unless `--ligand-charge/-l` derivation succeeds (PDB inputs). Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted (PDB inputs only — XYZ/GJF must supply `-q` explicitly). | _None_ |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only; `workers_per_node` forwarded to the parallel predictor). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | When loading PDB active site models, freeze the parent atoms of link hydrogens. See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-nodes INT` | Internal nodes per MEP segment (GSM string images or DMF images). | `20` |
| `--max-cycles INT` | Maximum MEP optimization cycles (GSM/DMF). | `300` |
| `--climb/--no-climb` | Enable climbing image for GSM segments (bridge segments always run without climbing). | `True` |
| `--opt-mode TEXT` | Single-structure optimizer for HEI±1/kink nodes. `grad` maps to L-BFGS; `hess` maps to RFO. See {ref}`opt-mode-semantics` for how the same token maps across subcommands (tsopt uses Dimer/RS-I-RFO, not L-BFGS/RFO). | `grad` |
| `--mep-mode {gsm\|dmf}` | Segment generator: GSM (string-based) or DMF (direct flux). | `gsm` |
| `--refine-mode {peak\|minima}` | Seeds for refinement: `peak` optimizes HEI±1; `minima` searches outward from the HEI toward the nearest local minima on each side. Defaults to `peak` for GSM and `minima` for DMF when omitted. | _Auto_ |
| `--dump/--no-dump` | Dump MEP (GSM/DMF) and single-structure trajectories. Restart YAML is written only when enabled in YAML. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB or Gaussian inputs. XYZ/GJF inputs do not produce a PDB companion of their own primary trajectory. | `True` |
| `-o, --out-dir TEXT` | Output directory. | `./result_path_search/` |
| `--thresh TEXT` | Override convergence preset for single-structure optimizations only (`opt.lbfgs/rfo.thresh`). | `gau` |
| `--thresh-stopt TEXT` | Override convergence preset for the string optimizer (`stopt.thresh`). | `gau_loose` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layer metadata) and continue. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running path search. | `False` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint with the selected single-structure optimizer (L-BFGS/RFO) before MEP search. | `True` |
| `--align/--no-align` | Align all inputs to the first structure before searching. | `True` |
| `--ref-full-pdb PATH...` | Full-size template PDBs (one per input, unless `--align` lets you reuse the first). | _None_ |
| `--ref-pdb PATH...` | Active site model reference PDBs used for the final full-system merge when inputs are XYZ/GJF (one per input, matching input order). | _None_ |

See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

## YAML configuration

The YAML root must be a mapping. Shared sections reuse [YAML Reference](yaml-reference.md): `geom`/`calc` mirror single-structure options (with `--freeze-links` augmenting `geom.freeze_atoms` for PDBs), and `stopt` inherits the StringOptimizer knobs documented for `path-opt` (see [path-opt.md](path-opt.md)).

`bond` and `search` are central to the recursion logic and shown below; `gs`, `dmf`, `stopt`, `opt.lbfgs`, and `opt.rfo` are reproduced only for the `path-search`-specific `out_dir` overrides.

`bond` carries the MLIP-based bond-change detection parameters shared with {ref}`scan <section-bond>`: `device`, `bond_factor`, `margin_fraction`, and `delta_fraction`.

### `path-search`-specific overrides

```yaml
stopt:
 out_dir: ./result_path_search/ # path-search override (canonical default: ./result_path_opt/)
opt:
 lbfgs:
   out_dir: ./result_path_search/ # path-search override (canonical default: ./result_opt/)
 rfo:
   out_dir: ./result_path_search/ # path-search override (canonical default: ./result_opt/)
```

`bond` and `search` are kept here because they are central to the `path-search` recursion logic:

```yaml
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
search:
 max_depth: 10 # recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 20 # max nodes per segment
 max_nodes_bridge: 5 # max nodes per bridge
 kink_max_nodes: 3 # max nodes for kink optimizations
 max_seq_kink: 2 # max sequential kinks
 refine_mode: null # optional refinement strategy (auto-chooses when null)
```

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [path-opt](path-opt.md) — Single-pass MEP optimization (no recursive refinement)
- [tsopt](tsopt.md) — Optimize the HEI as a transition state
- [extract](extract.md) — Generate active site model PDBs for path-search inputs
- [all](all.md) — End-to-end workflow (uses recursive path-search by default; `--refine-path False` for single-pass path-opt)
- [YAML Reference](yaml-reference.md) — Full `gs`, `dmf`, `bond`, `search` configuration options
- [Glossary](glossary.md) — Definitions of MEP, GSM, DMF, HEI
