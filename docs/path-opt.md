# `path-opt`

`pdb2reaction path-opt` searches for a minimum-energy path (MEP) between **exactly two** structures with GSM (default) or DMF (`--mep-mode dmf`). It writes the path trajectory and exports the highest-energy image (HEI) as a TS candidate. Treat the HEI as a *candidate* transition state until it is validated with [tsopt](tsopt.md) (which includes an imaginary-frequency check) and [irc](irc.md). For workflows that start from **two or more** structures and automatically refine only the reactive region, use [path-search](path-search.md).

Use it when you have exactly two endpoint structures (R → P) and need a first-pass MEP without recursive refinement. Choose GSM (default) for a string-based path generator, or switch to DMF with `--mep-mode dmf` for the Direct Max Flux generator.

An MLIP backend (UMA by default; switch with `-b/--backend` to ORB, MACE, or AIMNet2) provides energies and forces for MEP images. Hessians are used only by selected single-structure optimization steps, not for GSM/DMF path scoring. Before optimization starts, a rigid-body alignment step keeps the string stable.

```{note}
**Frozen atoms in DMF mode** use `HarmonicFixAtoms` (harmonic restraints with k=300 eV/Å²) instead of pysisyphus's hard coordinate freeze used by GSM. This means frozen atoms in DMF can move slightly from their reference positions, which differs from the rigid freeze in GSM mode.
```

## Examples

Command form:

```bash
pdb2reaction path-opt -i REACTANT.{pdb|cif|mmcif|xyz|gjf} PRODUCT.{pdb|cif|mmcif|xyz|gjf} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] \
 [--workers N] [--workers-per-node N] \
 [--mep-mode {gsm|dmf}] [--freeze-links/--no-freeze-links] [--max-nodes N] [--max-cycles-gsm N] [--max-cycles-dmf N] \
 [--climb/--no-climb] [--dump/--no-dump] [--thresh PRESET] [--thresh-gsm PRESET] [--thresh-dmf TOL] \
 [--preopt/--no-preopt] [--preopt-max-cycles N] [--opt-mode grad|hess] [--fix-ends/--no-fix-ends] \
 [--show-config/--no-show-config] [--dry-run/--no-dry-run] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

MEP search between two endpoints:

```bash
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --out-dir ./result_path_opt
```

Pre-optimize endpoints before MEP search:

```bash
# Pre-optimize endpoints before MEP search
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

Use DMF mode instead of GSM:

```bash
# Use DMF mode instead of GSM
pdb2reaction path-opt -i reactant.pdb product.pdb -q 0 -m 1 \
 --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

```{note}
DMF mode additionally requires `cyipopt` (install from conda-forge before running with `--mep-mode dmf`). `pydmf` ships with `pdb2reaction` as a dependency. The default `--dmf-backend gpu` uses the PyTorch/CUDA `dmf.torch` backend; pass `--dmf-backend cpu` (`dmf`/NumPy) on a GPU out-of-memory error.

```

A quick pass that freezes cap parents and disables climb: add `--freeze-links --no-climb`.

## Workflow

1. **Pre-alignment & freeze resolution**
 - All endpoints after the first are Kabsch-aligned to the first structure. If either endpoint defines `freeze_atoms`, only those atoms participate in the RMSD fit and the resulting transform is applied to every atom.
 - When `--freeze-links` is active, cap-hydrogen parent atoms are automatically frozen (see {ref}`Cap hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
2. **String growth and HEI export**
 - After the path is grown and refined, the global maximum-energy image is exported as the HEI. An endpoint HEI is not an internal transition-state candidate and requires a separate TS optimization before TS-specific analysis.
 - The highest-energy image (HEI) is written both as `.xyz` and `.pdb` when a PDB reference exists, and as `.gjf` when a Gaussian template is available; these conversions honor `--convert-files`.

## Outputs

```text
out_dir/
├─ final_geometries_trj.xyz # XYZ path; comment line holds energies when provided
├─ final_geometries.pdb # PDB of every image when PDB/mmCIF topology is available and conversion enabled
├─ final_geometries.cif # Bridged mmCIF/oversized-PDB input: original IDs restored
├─ final_geometries.gjf # GSM only: Gaussian companion when a template is detected (conversion enabled)
├─ hei.xyz # Highest-energy image with its energy on the comment line
├─ hei.pdb # HEI converted to PDB when a PDB reference is available (conversion enabled)
├─ hei.gjf # HEI written using a detected Gaussian template (conversion enabled)
├─ align_refine/ # Intermediate files from the rigid alignment/refinement stage (created when alignment runs)
└─ <optimizer dumps> # Trajectory dumps when --dump (restart YAML only via YAML dump_restart)
```

Console output echoes the resolved YAML blocks and prints cycle-by-cycle MEP progress (GSM/DMF) with timing information.

See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product structures (`.pdb`/`.cif`/`.mmcif`/`.xyz`/`.gjf`). | Required |
| `-q, --charge INT` | Total charge (`calc.charge`). Required for non-`.gjf` inputs unless `--ligand-charge/-l` derivation succeeds (PDB/mmCIF inputs or XYZ/GJF with `--ref-pdb`). `.gjf` templates can supply it; if `.gjf` inputs lack charge metadata, the run aborts unless `-q` is provided. Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex for PDB/mmCIF inputs (or XYZ/GJF when `--ref-pdb` is supplied). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism; `workers_per_node` is forwarded to the parallel predictor. `workers > 1` cannot be combined with an explicit analytical Hessian request. See {ref}`workers-analytical-error`. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (`calc.spin`). | Template/`1` |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF input (or XYZ/GJF with `--ref-pdb`): freeze cap-H parents (merged with YAML). See [extract](extract.md) for cap-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-nodes INT` | Number of movable internal images for GSM or DMF. Both engines retain two endpoints, so total images = `max_nodes + 2`. | `20` |
| `--mep-mode {gsm\|dmf}` | Select GSM (string-based) or DMF (Direct Max Flux) path generator. | `gsm` |
| `--dmf-backend {cpu\|gpu}` | DMF compute backend (`--mep-mode dmf` only): `gpu` (`dmf.torch`/CUDA) or `cpu` (`dmf`/NumPy). On a GPU out-of-memory error, retry with `cpu`. | `gpu` |
| `--max-cycles-gsm INT` | GSM string-optimizer cycle cap (sets `stopt.max_cycles` and `stopt.stop_in_when_full`). | `300` |
| `--max-cycles-dmf INT` | DMF IPOPT iteration cap (sets `dmf.max_cycles`). | `300` |
| `--climb/--no-climb` | Enable GSM climbing-image refinement (and Lanczos tangent). Accepted but unused with DMF. | `True` |
| `--dump/--no-dump` | Dump GSM/single-optimizer trajectories. Accepted but unused by the DMF path solver. Restart YAML is written only when enabled in YAML. | `False` |
| `--opt-mode TEXT` | Single-structure optimizer for endpoint preoptimization (`grad` = L-BFGS, `hess` = RFO). | `grad` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/CIF/GJF companions according to the input topology/template. | `True` |
| `--ref-pdb FILE` | Reference PDB topology for XYZ/GJF inputs (keeps XYZ coordinates) to enable PDB conversions. | _None_ |
| `-o, --out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--thresh TEXT` | Override convergence preset for endpoint preoptimization only (`opt.lbfgs/rfo.thresh`). | `gau` |
| `--thresh-gsm TEXT` | Override convergence preset for the GSM string optimizer (`stopt.thresh`). | `gau_loose` |
| `--thresh-dmf TEXT` | Override the IPOPT dual-infeasibility tolerance of the DMF optimizer (`dmf.tol`): `tight` (0.04), `middle` (0.10), `loose` (0.20), or a positive float. Gaussian presets are rejected. | `tight` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layers) and continue. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running optimization. | `False` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint with the selected single-structure optimizer before alignment/MEP search (GSM/DMF). | `True` |
| `--preopt-max-cycles INT` | Cap for endpoint preoptimization cycles. | `100000` |
| `--fix-ends/--no-fix-ends` | Keep endpoint geometries fixed during GSM growth/refinement. Accepted but unused with DMF. | `True` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |

## YAML configuration

### YAML sections used by `path-opt`

See [YAML Reference](yaml-reference.md) for full key listings:

- [`geom`](yaml-reference.md#geom) — `--freeze-links` augments `freeze_atoms` for PDB/mmCIF topology inputs.
- [`calc`](yaml-reference.md#calc) — MLIP backend setup.
- [`gs`](yaml-reference.md#gs) — Growing String representation (GSM mode).
- [`dmf`](yaml-reference.md#dmf) — Direct Max Flux + (C)FB-ENM interpolation (DMF mode).
- [`stopt`](yaml-reference.md#stopt) — StringOptimizer settings.
- [`opt.lbfgs`](yaml-reference.md#lbfgs) / [`opt.rfo`](yaml-reference.md#rfo) — Endpoint single-structure preoptimization. An explicitly supplied CLI `--preopt-max-cycles` overrides the YAML value; when omitted, YAML remains effective.

### `path-opt`-specific defaults

The following keys differ from the canonical defaults when invoked via `path-opt`:

```yaml
stopt:
 out_dir: ./result_path_opt/ # output directory (path-opt default)
opt:
 lbfgs:
   out_dir: ./result_path_opt/ # output directory (path-opt default)
 rfo:
   out_dir: ./result_path_opt/ # output directory (path-opt default)
```

## Exit codes

See {ref}`exit-codes` in CLI Conventions.

## See Also

- [path-search](path-search.md) — Recursive MEP search with automatic refinement (for 2+ structures)
- [tsopt](tsopt.md) — Optimize the HEI as a TS candidate (includes imaginary-frequency check; follow with IRC)
- [extract](extract.md) — Generate active site model (binding pocket) PDBs for path-opt inputs
- [all](all.md) — End-to-end workflow (defaults to single-pass path-opt; add `--refine-path` for recursive path-search. The `--refine-path` flag lives on `pdb2reaction all` only — see [all.md → MEP search](all.md#mep-search) for its definition.)
- [YAML Reference](yaml-reference.md) — Full `gs`, `dmf`, `stopt`, `opt` configuration options
- [Glossary](glossary.md) — Definitions of MEP, GSM, DMF, HEI
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
