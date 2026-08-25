# Architecture: pdb2reaction

## 1. Overview

`pdb2reaction` is a Python CLI for enzymatic reaction-path analysis on an
active-site cluster model. It uses built-in MLIPs or a custom ASE calculator
for geometry/path stages and can add an optional PySCF/GPU4PySCF DFT
single-point stage (extract → MEP → tsopt → IRC → freq → dft).


Two bundled forks (`pysisyphus/`, `thermoanalysis/`) live at the repo top as repo-internal modules. They are deliberately **not** the upstream PyPI distributions; reinstalling them from PyPI alongside this package silently breaks the local extensions. See §6.

---

## 2. Layered structure (6 core layers)

> The package also ships a 7th, optional-integration subpackage `pdb2reaction/mcp/` (MCP server: `server.py`, `_tools.py`, `_runner.py`; entry points `p2r-mcp` / `pdb2reaction-mcp`). It sits outside the 6-layer core graph below and is documented in [mcp_server.md](mcp_server.md).

### 2.1 Layer table

| layer | dir | responsibility | may depend on |
|---|---|---|---|
| **L1 Interface** | `pdb2reaction/cli/` | Click root group, shared option-decorator factories (`common_options.py`), `--help-advanced`, bool flag normalization, subcommand resolver | `workflows/`, `core/` |
| **L2 Application** | `pdb2reaction/workflows/` | per-subcommand orchestration; one file per stage runner (`all.py`, `path_search.py`, `tsopt.py`, `extract.py`, `irc.py`, `freq.py`, `dft.py`, …) | `domain/`, `backends/`, `io/`, `core/` |
| **L3 Domain** | `pdb2reaction/domain/` | chemistry-aware helper logic (bond change detection, bond summary, element-info propagation) | `core/` |
| **L4a Infra (MLIP)** | `pdb2reaction/backends/` | MLIP backend dispatcher + per-backend adapter (UMA / Orb / MACE / AIMNet2) | `core/` |
| **L4b Infra (I/O)** | `pdb2reaction/io/` | output layout, summary, trajectory, PDB fix, energy diagram, Hessian cache | `core/` |
| **L5 Foundation** | `pdb2reaction/core/` | defaults (primary source for shared defaults), utilities, logging, output, and result publication | (none, by design intent) |
| (bundle, not a layer) | `<repo>/pysisyphus/`, `<repo>/thermoanalysis/` | repo-internal forks (optimizer / thermochemistry) | (sibling, layer-external) |

**Dependency direction (design goal)**: `L1 → L2 → {L3, L4} → L5`. The full top-to-bottom direction is not mechanically enforced. The enforced subset keeps the product module graph **acyclic** (no strongly connected component among `pdb2reaction/*` modules) and prevents any `core`/`domain` module from importing a `workflows/*` module. Two CI gates cover different invariants: `.github/scripts/check_engineering_markers.py` verifies `# CHEMISTRY-RULE:{4,5,7}` marker coverage, the `# DOMAIN_PURE` markers on `workflows/dft.py`, `workflows/tsopt.py`, and `workflows/sp.py`, and that the MLIP SDKs (`fairchem` / `orb_models` / `mace` / `aimnet`) are imported only under `backends/`; `.github/scripts/check_import_graph.py` (a static AST import graph) asserts there is no product cycle, no `pysisyphus/**` import of `pdb2reaction`, and no `core → workflows` / `domain → workflows` back-edge. Remaining allowed downward edges / one-way facades are: `workflows/* → cli.common_options` / `cli.decorators` / `cli.help_pages`; `core/utils.py → domain.add_elem_info`, `io.structure_formats`, and `io.charge`; `io/charge.py → domain.residue_data` and `io.structure_formats`; `io/structure_formats.py → domain.add_elem_info`; and `io/trj2fig.py → backends`. The canonical residue tables (`domain/residue_data.py`), the charge engine (`io/charge.py`), and the console-gated charge-summary logger (`core/utils.py`) are re-exported from `workflows/extract.py` so existing import paths keep working. Bundled forks sit outside the layer graph and may be imported from any layer via their absolute package path (`from pysisyphus.X import Y`).

### 2.2 ASCII map of the package tree

```
pdb2reaction/ [GH: t-0hmura/pdb2reaction]
├── pyproject.toml packages.find = ["pdb2reaction*", ...] (glob, frozen)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── architecture.md ← this file
│ └──... (Sphinx documentation site)
├── pdb2reaction/ ← package body, 6-layer physical dir
│ ├── __init__.py PEP 562 lazy: _LAZY_SYMBOLS / _LAZY_MODULES + __getattr__
│ ├── __main__.py `from .cli import cli`
│ ├── _version.py / py.typed
│ │
│ ├── cli/ # === L1 Interface ===
│ │ ├── app.py Click group + _LAZY_SUBCOMMANDS registry (absolute paths)
│ │ ├── common_options.py @add_print_every_option / @add_irc_pos_def_option / @add_precision_option / @add_coord_type_option / @add_ml_charge_spin_options
│ │ ├── decorators.py resolve_yaml_sources / load_merged_yaml_cfg / _write_error_json
│ │ ├── help_pages.py --help-advanced pager
│ │ ├── bool_compat.py --flag / --no-flag normalization
│ │ └── default_group.py subcommand resolver, lazy module import
│ │
│ ├── workflows/ # === L2 Application ===
│ │ ├── all.py full pipeline orchestrator (extract → … → DFT)
│ │ ├── path_search.py / path_opt.py MEP search / COS wrapper
│ │ ├── tsopt.py / freq.py / irc.py / dft.py per-stage runners
│ │ ├── opt.py / sp.py / scan.py / scan2d.py /
│ │ │ scan3d.py / scan_common.py geometry opt / single point / scans
│ │ ├── extract.py active-site extraction CLI
│ │ ├── restraints.py restraint helpers
│ │ └── align_freeze.py Kabsch + frozen-subset rmsd
│ │
│ ├── domain/ # === L3 Domain ===
│ │ ├── bond_changes.py R↔P bond detection
│ │ ├── bond_summary.py post-IRC diagnostic
│ │ ├── add_elem_info.py PDB element column normalizer
│ │ └── residue_data.py residue/ion charge tables
│ │
│ ├── backends/ # === L4a Infra (MLIP) ===
│ │ ├── __init__.py backend dispatch + registry
│ │ ├── base.py MLIPCalculator protocol
│ │ ├── custom.py custom ASE-calculator adapter
│ │ ├── _determinism.py deterministic reduction shim
│ │ └── uma.py / orb.py / mace.py / aimnet2.py per-backend adapters
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.log writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / SVG / PDF / HTML / CSV
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── altloc.py altloc identity/selection helpers
│ │ ├── charge.py residue-aware charge resolution
│ │ ├── structure_formats.py PDB/mmCIF bridge + identifier restoration
│ │ └── hessian_cache.py in-memory Hessian cache
│ │
│ └── core/ # === L5 Foundation ===
│   ├── defaults.py primary source for shared defaults
│   ├── utils.py PDB / XYZ / plot helpers
│   ├── output.py / result_commit.py output/result ownership
│   └── pes_composition.py energy-component composition
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (CI, release, engineering, and documentation checks)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ repo-internal fork (slimmed to the numerical features used here)
 thermoanalysis/ repo-internal fork
```

### 2.3 Per-layer responsibility detail

**L1 `cli/`** owns the root group, lazy registry, argv normalization, progressive help, and shared option-decorator factories. Workflow and utility modules own their actual `@click.command()` objects and may define command-local options inline; the root resolves and invokes those objects. Every `_LAZY_SUBCOMMANDS` entry uses an **absolute module path** (`pdb2reaction.workflows.all`, `pdb2reaction.io.trj2fig`, …).

**L2 `workflows/`** contains compute-command modules plus shared stage helpers such as `scan_common.py`, `restraints.py`, and `align_freeze.py`. A command module normally exposes one `@click.command()` named `cli`; utility commands also live in `io/` and `domain/`, so neither “one file per subcommand” nor “every workflow file is a command” is an invariant.

**L3 `domain/`**. Chemistry-aware helper logic that may import `torch` / `numpy` / `pysisyphus.constants` (numeric backends), but **may not import** MLIP runtimes (`fairchem`, `orb_models`, `mace`, `aimnet`). `.github/scripts/check_engineering_markers.py` enforces this deny list via an external-library import-scope check across non-`backends/` files. (The `# DOMAIN_PURE` docstring marker itself lives on selected workflow modules — `workflows/dft.py`, `tsopt.py`, `sp.py` — not on `domain/`.) Domain helpers are reusable by any L2 stage runner.

**L4a `backends/`**. MLIP backend dispatcher (`__init__.py` + `base.py`) plus one adapter per supported MLIP (`uma.py`, `orb.py`, `mace.py`, `aimnet2.py`). The dispatcher also exposes the custom ASE-calculator boundary.

**L4b `io/`**. Output-side I/O concerns: per-stage summary writer, energy diagram, trajectory rendering, PDB altloc fix, PDB/mmCIF bridge/template restoration, and in-memory Hessian cache. Output format is owned here and consumed by stage runners; current non-foundation imports are listed in the back-edge inventory above.

**L5 `core/`** is the lowest layer. `defaults.py` is the **primary source** for shared numerical and CLI defaults; grep it first, then inspect justified command-local defaults (for example, path-engine choices). `utils.py` contains shared configuration, structure, coordinate, and plotting helpers.

### 2.4 Lazy-import mechanism (conceptual diagram)

```text
External consumer                        Package root             Layer dir
---------------------------------------  ----------------------   ---------

from pdb2reaction.core.utils import x ──► (direct dotted import) ──► pdb2reaction/core/utils.py
import pdb2reaction.io.trj2fig        ──► (direct dotted import) ──► pdb2reaction/io/trj2fig.py

from pdb2reaction import <Symbol>     ──► pdb2reaction/__init__.py
                                          __getattr__("<Symbol>")
                                          └─► _LAZY_SYMBOLS["<Symbol>"]
                                              = "pdb2reaction.<layer>.<module>"
                                              └─► importlib.import_module(...)

from pdb2reaction import <module>     ──► pdb2reaction/__init__.py
        (= module attr)                   __getattr__("<module>")
                                          └─► _LAZY_MODULES["<module>"]
                                              = "pdb2reaction.<layer>.<module>"
                                              └─► importlib.import_module(...) returns module

pdb2reaction myaction                 ──► pdb2reaction/cli/app.py
                                          _LAZY_SUBCOMMANDS["myaction"]
                                          = ("pdb2reaction.workflows.myaction", "cli", "...")
                                          └─► importlib.import_module(absolute path)
                                              └─► getattr(module, "cli") → Click command
```

Two layers of lazy-import compatibility plus CLI dispatch:

1. **Root symbol attribute** (`from pdb2reaction import <Symbol>`) — handled by `pdb2reaction/__init__.py:_LAZY_SYMBOLS` + PEP 562 `__getattr__`. Symbols are loaded on first access from the layer-dir path; import cost stays zero at `pdb2reaction` import time.
2. **Root module attribute** (`from pdb2reaction import <module>`) — handled by `_LAZY_MODULES`. `__getattr__` returns the module object itself via `importlib.import_module`. The registry is currently empty, so no root module-attribute paths are exported.

The CLI subcommand resolver (`cli/app.py:_LAZY_SUBCOMMANDS`) uses **absolute** module paths (e.g. `"pdb2reaction.workflows.all"`) so that moving `default_group.py` into `cli/` does not silently break subcommand discovery (the registry no longer depends on `__package__`).

---

## 3. Fresh-eyes 5-step navigation (≈ 40 min total)

For a contributor opening the repo for the first time, follow this path top-to-bottom; each step closes one concern.

| step | minutes | open | what you learn |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/README.md) | one-paragraph elevator pitch + single-command usage |
| 2 | 5 | this file (`docs/architecture.md`) §2 + §4 | 6-layer dir tree, dependency direction, where each concern lives |
| 3 | 5 | [`pdb2reaction/cli/app.py`](../pdb2reaction/cli/app.py) | Click root group, `_LAZY_SUBCOMMANDS` registry, absolute-path resolution |
| 4 | 20 | [`pdb2reaction/workflows/all.py`](../pdb2reaction/workflows/all.py) (skim) | one full subcommand top-to-bottom; trace `extract → MEP → tsopt → IRC → freq → dft` |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/CONTRIBUTING.md) §3 + §4 | 5 add-a-X recipes + the "do not touch" hidden constraints |

After step 5 you can read any other file by following the file index in §4. The package is intentionally **flat-within-each-layer** — there is no nested package below `pdb2reaction/<layer>/`, so you never need to navigate more than two directories deep.

---

## 4. File index — "where does this concern live?"

### 4.1 CLI / entry (L1 `cli/`)

| concern | file |
|---|---|
| Click root group + subcommand dispatch | `pdb2reaction/cli/app.py` |
| Subcommand resolver (lazy import) | `pdb2reaction/cli/default_group.py` |
| `python -m pdb2reaction` entry | `pdb2reaction/__main__.py` |
| YAML source resolution + standardized exception handling | `pdb2reaction/cli/decorators.py` |
| `--help-advanced` pager | `pdb2reaction/cli/help_pages.py` |
| Bool flag compat (`--flag` / `--no-flag` + value style) | `pdb2reaction/cli/bool_compat.py` |
| Shared option-decorator factories (`--print-every`, `--irc-pos-def`, `--precision`, `--coord-type`, `--charge / --ligand-charge / --multiplicity`) | `pdb2reaction/cli/common_options.py` |

### 4.2 Workflow stage runners (L2 `workflows/`)

| concern | file |
|---|---|
| Full pipeline orchestrator | `pdb2reaction/workflows/all.py` |
| Geometry optimization (L-BFGS / RFO) | `pdb2reaction/workflows/opt.py` |
| Scan and 2D/3D energy-landscape grids + shared | `pdb2reaction/workflows/scan{,2d,3d,_common}.py` |
| MEP search (GSM) | `pdb2reaction/workflows/path_search.py` |
| MEP optimizer core (pysisyphus COS) | `pdb2reaction/workflows/path_opt.py` |
| TS optimization (RS-P-RFO + Bofill + macro/micro) | `pdb2reaction/workflows/tsopt.py` |
| Vibrational analysis (PHVA + backend-agnostic active-DOF Hessian handling) | `pdb2reaction/workflows/freq.py` |
| IRC integration (macro / micro) | `pdb2reaction/workflows/irc.py` |
| Single-point DFT (gpu4pyscf subprocess) | `pdb2reaction/workflows/dft.py` |
| Active-site extraction (cluster cap) | `pdb2reaction/workflows/extract.py` |
| Restraint helpers | `pdb2reaction/workflows/restraints.py` |
| Kabsch / frozen-subset alignment | `pdb2reaction/workflows/align_freeze.py` |

### 4.3 Chemistry helpers (L3 `domain/`)

| concern | file |
|---|---|
| R↔P bond change detection | `pdb2reaction/domain/bond_changes.py` |
| Post-IRC bond summary | `pdb2reaction/domain/bond_summary.py` |
| PDB element column normalizer | `pdb2reaction/domain/add_elem_info.py` |
| Residue and ion charge tables | `pdb2reaction/domain/residue_data.py` |

### 4.4 MLIP backends (L4a `backends/`)

| concern | file |
|---|---|
| Backend dispatch + registry | `pdb2reaction/backends/__init__.py` |
| `MLIPCalculator` protocol + base | `pdb2reaction/backends/base.py` |
| Custom ASE-calculator adapter | `pdb2reaction/backends/custom.py` |
| Deterministic reduction shim | `pdb2reaction/backends/_determinism.py` |
| Per-backend adapters | `pdb2reaction/backends/{uma, orb, mace, aimnet2}.py` |

See [Backends](backends.md) for the add-a-backend recipe.

### 4.5 I/O (L4b `io/`)

| concern | file |
|---|---|
| `summary.json` / `summary.log` writer | `pdb2reaction/io/summary.py` |
| Plotly energy diagram | `pdb2reaction/io/energy_diagram.py` |
| Trajectory → PNG / SVG / PDF / HTML / CSV | `pdb2reaction/io/trj2fig.py` |
| PDB altloc resolution | `pdb2reaction/io/pdb_fix.py` |
| Altloc identity and coherent selection | `pdb2reaction/io/altloc.py` |
| Residue-aware charge resolution | `pdb2reaction/io/charge.py` |
| PDB/mmCIF bridge + identifier restoration | `pdb2reaction/io/structure_formats.py` |
| In-memory Hessian cache (per-run TTL) | `pdb2reaction/io/hessian_cache.py` |
| Harmonic restraint setup | `pdb2reaction/workflows/restraints.py` (L2 stage helper) |

### 4.6 Foundation (L5 `core/`)

| concern | file |
|---|---|
| **Shared/default numerical settings (primary source; verify command-local exceptions)** | `pdb2reaction/core/defaults.py` |
| PDB / XYZ / plot helpers | `pdb2reaction/core/utils.py` |
| `-v` / `--verbose LEVEL` logging wiring | `pdb2reaction/core/logging.py` |
| Output/result ownership helpers | `pdb2reaction/core/output.py`, `pdb2reaction/core/result_commit.py` |
| Energy-component composition | `pdb2reaction/core/pes_composition.py` |

### 4.7 Repo-internal bundled forks

| dir | role | selected divergent files (see the directory README for the complete list) |
|---|---|---|
| `pysisyphus/` | optimizer / TS / IRC engine | `irc/IRC.py` (opt-in `require_pos_def_hessian` PSD convergence guard), `optimizers/hessian_updates.py` (GPU-resident in-place rank-two Bofill update with explicit `PYSIS_BOFILL_CPU_OFFLOAD=1` fallback), `tsoptimizers/{RSIRFOptimizer,RSPRFOptimizer,TRIM,TSHessianOptimizer}.py`, `calculators/{Calculator,Dimer}.py`, `_array.py` (torch/numpy backend shim) |
| `thermoanalysis/` | thermochemistry (ΔG, ZPE, partition functions) | `QCData.py` (branding diff vs upstream) |

See each dir's `README.md` for the touch-restriction boundary.

---

## 5. Hidden constraints (read this before any patch)

### 5.1 Chemistry rules (grep recipe)

Three correctness-critical rules are implemented in `workflows/dft.py` and `workflows/tsopt.py`. They are **not** detected by smoke tests — silent drift here breaks reaction-path accuracy. Inline `# CHEMISTRY-RULE:N` markers and `# DOMAIN_PURE` module-docstring markers identify the rules; `.github/scripts/check_engineering_markers.py` enforces marker completeness in CI.

To find every chemistry rule before editing:

```bash
# List all rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' pdb2reaction/

# List every # DOMAIN_PURE marker (= chemistry-rule host modules)
grep -rn '# DOMAIN_PURE' pdb2reaction/
```

The three rules (marker IDs are non-contiguous) are:

| marker | rule | host file |
|---|---|---|
| 4 | gpu4pyscf `rks_lowmem` triple-guard | `pdb2reaction/workflows/dft.py` |
| 5 | def2 family auto-ECP injection | `pdb2reaction/workflows/dft.py` |
| 7 | `bofill_update` advanced-indexing scatter | `pdb2reaction/workflows/tsopt.py` |

Editing any of these requires the `[CHEMISTRY-RULE:N]` commit prefix, lab
sign-off, the relevant focused regression tests, and a HEAVY benchmark; see
`CONTRIBUTING.md` §4.1 and §5.

### 5.2 VRAM-management invariant (do not refactor `del` chains)

The IRC / TSopt / Freq stages explicitly `del` GPU-resident objects (`calc`, `geom`, `hess`) between stages to free CUDA memory; stage boundaries also use `gc.collect()` and, where CUDA allocations are present, `torch.cuda.empty_cache()`. **Do not refactor these release operations out** — long-running `all` jobs with large active-site models OOM without them.

### 5.3 Bundled forks: do NOT install upstream alongside

The bundled `pysisyphus/` and `thermoanalysis/` packages are **forks**. Reinstalling `pip install pysisyphus` or `pip install thermoanalysis` next to this package silently breaks:

- `pysisyphus/irc/IRC.py` — initial-displacement memory hygiene + opt-in `require_pos_def_hessian` kwarg
- `pysisyphus/optimizers/hessian_updates.py` — GPU-resident in-place rank-two Bofill update; opt-in `PYSIS_BOFILL_CPU_OFFLOAD=1` fallback
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — RSIRFO kwargs (host package import path divergent between forks)
- `pysisyphus/calculators/{Calculator,Dimer}.py` — GPU-aware backend hooks (the 30+ QM backends have been removed; only the abstract base + Dimer TS calculator remain)
- `pysisyphus/_array.py` — `get_xp` / `_outer` / `_dot` / `_eigh` shim used by `optimizers/hessian_updates.py` (and progressively by other hot-path files)
- `thermoanalysis/QCData.py` — branding / I/O diff vs upstream

### 5.4 Packaging changes require an isolated-install gate

The `include` glob (`pdb2reaction*`) already auto-discovers new layer
subpackages, so ordinary internal files need no package-discovery edit. If a
feature really changes package discovery or runtime dependencies, build both
sdist and wheel, inspect wheel contents, install the wheel in a clean
environment, and run the CLI smoke suite. Treat this as a release-contract
change, not as an incidental refactor.

### 5.5 `_LAZY_SUBCOMMANDS` registry must use absolute paths

`pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` resolves every subcommand through an **absolute** module path. Switching any entry back to a relative dotted import (`".all"` etc.) silently breaks subcommand discovery whenever `default_group.py` moves, because the resolver's `__package__` then drifts away from the package root.

---

## 6. Bundled forks (repo-internal)

`pdb2reaction` ships **two** repo-internal modules at the repo top:

| dir | upstream PyPI? | purpose | scope of edits allowed |
|---|---|---|---|
| `pysisyphus/` | NO — fork, do not `pip install pysisyphus` alongside | optimizer, TS, IRC, COS, calculators | maintained locally; behavior changes need a focused regression test and matching numerical benchmark |
| `thermoanalysis/` | NO — fork (branding/I/O diff) | ΔG, ZPE, partition functions, `QCData` | maintained locally; numerical/I/O changes need thermochemistry golden tests |

Each directory carries a `README.md` listing its divergent files and required
validation. From the layer model these forks live **outside** the L1..L5
graph: any layer may import them via the absolute package path
(`from pysisyphus.X import Y`) without breaking the layer direction.

---

## 7. Recommended deeper reading order

After the Fresh-eyes tour (§3), follow this depth-first reading order:

1. `pdb2reaction/core/defaults.py` — internalise the default-value table; everything downstream reads from here.
2. `pdb2reaction/cli/app.py` — Click root + `_LAZY_SUBCOMMANDS` registry.
3. `pdb2reaction/workflows/all.py` — one full pipeline top-to-bottom.
4. `pdb2reaction/workflows/extract.py` — active-site cluster cap.
5. `pdb2reaction/backends/__init__.py` + `base.py` — MLIP dispatcher and per-backend adapter contract.
6. `pdb2reaction/workflows/tsopt.py` — RS-P-RFO + Bofill scatter (CHEMISTRY-RULE:7).
7. `pdb2reaction/workflows/freq.py` — vibrational analysis on the cluster model.
8. `pdb2reaction/workflows/irc.py` — VRAM hygiene + IRC integration.
9. `pdb2reaction/workflows/dft.py` — single-point DFT with gpu4pyscf (CHEMISTRY-RULE:4 + :5).
10. `pdb2reaction/core/utils.py` — shared PDB / XYZ / plot helpers.
