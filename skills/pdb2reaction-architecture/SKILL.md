---
name: pdb2reaction-architecture
description: Where the source code lives in `pdb2reaction`. 6 physical layer directories (`cli` / `workflows` / `domain` / `backends` / `io` / `core`) + 2 repo-internal forks (`pysisyphus` / `thermoanalysis`). Tells an agent which directory to grep for a given concern (Click option, stage runner, MLIP backend, output writer, chemistry default, cap-atom math) before touching code. TRIGGER on questions like "where is X implemented", "which file defines flag Y", "how is the repo organized", "what's safe to refactor". SKIP for usage questions — those belong to `pdb2reaction-cli` / `-overview`.
---

# pdb2reaction architecture (one-screen map)

## 6 layers + bundled forks

```
pdb2reaction/                          ← the package body, one folder per layer
├── cli/        # L1 — Click root group, --help-advanced, bool normalization,
│               #      shared option-decorator factories, subcommand resolver.
├── workflows/  # L2 — the 13 stage subcommands, one module each (`all.py`,
│               #      `scan.py`, `opt.py`, `path_opt.py`, `path_search.py`,
│               #      `tsopt.py`, `freq.py`, `irc.py`, `dft.py`, `sp.py`,
│               #      `scan2d.py`, `scan3d.py`, `extract.py`), plus shared
│               #      stage helpers (`scan_common.py`, `restraints.py`,
│               #      `align_freeze.py`, ...). The remaining 5 subcommands
│               #      live in `io/` and `domain/` — see the table below.
├── domain/     # L3 — chemistry-aware helpers (bond changes, bond summary,
│               #      element-info repair). Hosts the `bond-summary` and
│               #      `add-elem-info` subcommands. May use torch/numpy; no
│               #      MLIP SDK (fairchem/orb_models/mace/aimnet) dependency.
├── backends/   # L4a — MLIP backend dispatcher + per-backend adapters
│               #       (`uma.py`, `orb.py`, `mace.py`, `aimnet2.py`).
├── io/         # L4b — summary writer, energy diagram, trajectory plot,
│               #       Hessian cache, PDB altloc fix. Hosts the `trj2fig`,
│               #       `energy-diagram` and `fix-altloc` subcommands.
└── core/       # L5 — `defaults.py` (primary source of truth for most
                #       CLI defaults), `utils.py` (PDB / XYZ / plot helpers),
                #       `logging.py`.

pysisyphus/        ← bundled fork of the optimizer / TS / IRC engine.
                     Slimmed to the subset pdb2reaction actually uses; its own
                     README's "Divergent files" table is the authoritative list
                     of locally maintained numerical code. Read that table
                     before editing: numerical changes need focused unit tests
                     and the relevant optimizer/IRC benchmark, while upstream
                     replacement is never safe.

thermoanalysis/    ← bundled fork for ΔG / ZPE / partition functions.
                     `workflows/freq.py` and `pysisyphus/Geometry.py` consume
                     it. Numerical or I/O changes need thermochemistry golden
                     tests; do not replace it with the upstream package.
```

The intended dependency direction is `L1 → L2 → {L3, L4} → L5`. The full
top-to-bottom direction is still a refactoring goal rather than a mechanically
enforced rule, but a subset **is** now enforced: `.github/scripts/check_import_graph.py`
(a static AST import graph) proves the product graph is acyclic and that no
`core`/`domain` module imports a `workflows/*` module (and no `pysisyphus/**`
file imports `pdb2reaction`). Remaining allowed downward edges include
`workflows/* → cli.common_options/cli.decorators/cli.help_pages`,
`core/utils.py → domain.add_elem_info/io.structure_formats/io.charge`, and
`io/charge.py → domain.residue_data`. The residue tables (`domain/residue_data.py`)
and charge engine (`io/charge.py`) are re-exported from `workflows/extract.py`,
so inspect current imports before moving code. The bundled forks sit outside the
layer graph and may be imported by any layer.

## Where to look first

| concern | open this |
|---|---|
| Default for any CLI flag | `pdb2reaction/core/defaults.py` (primary source of truth for most defaults — grep here first; a few workflow-local defaults live inline, e.g. path-opt `--mep-mode`) |
| Subcommand body / orchestration | Canonical registry = `pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` — grep the subcommand name there first; it maps every subcommand to its module. Stage subcommands live in `pdb2reaction/workflows/<subcmd>.py` (hyphens become underscores: `path-opt` → `path_opt.py`); the 5 utility subcommands live elsewhere — `trj2fig` → `io/trj2fig.py`, `fix-altloc` → `io/pdb_fix.py`, `energy-diagram` → `io/energy_diagram.py`, `add-elem-info` → `domain/add_elem_info.py`, `bond-summary` → `domain/bond_summary.py` |
| New MLIP backend | `pdb2reaction/backends/<backend>.py` + register in `backends/__init__.py:BACKEND_REGISTRY` |
| `--help` / option decorator | `pdb2reaction/cli/common_options.py` (shared) or the subcommand file (inline) |
| Aggregate `summary.json` schema | `pdb2reaction/workflows/all.py`, `pdb2reaction/workflows/_all_helpers.py`, `pdb2reaction/workflows/path_search.py`; common leaf/error envelope writer in `pdb2reaction/core/utils.py:write_result_json` |
| Human summary log / trajectory conversion / energy diagram | `pdb2reaction/io/` (the machine JSON schema is not owned solely by this layer) |
| Chemistry rule (#4 gpu4pyscf `rks_lowmem`, #5 def2 auto-ECP, #7 Bofill advanced-indexing) | grep `# CHEMISTRY-RULE:` — the only three markers in the package: `workflows/dft.py` (#4, #5) and `workflows/tsopt.py` (#7); maintainer approval and a scheduled numerical benchmark are required |
| Cap-H geometry | `pdb2reaction/workflows/extract.py` (`compute_linkH_atoms`) — chemistry-sensitive, and reached by opening the file (it carries no `CHEMISTRY-RULE` marker) |
| TS / IRC / optimizer internals | `pysisyphus/` — read `pysisyphus/README.md`, add a focused regression test, and run the relevant numerical benchmark for behavior changes |
| MCP server / agent integration | `pdb2reaction/mcp/` — see [`pdb2reaction-mcp`](../pdb2reaction-mcp/SKILL.md) |

## Hidden constraints to remember

1. **`pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS`** entries MUST use absolute module paths (`"pdb2reaction.workflows.all"`, never `".all"`). Relative dotted paths silently break the resolver if `default_group.py` moves.
2. **VRAM hygiene**: the explicit `del calc; gc.collect(); torch.cuda.empty_cache()` sequence between stages releases retained calculators before the next full-protein stage.
3. **Packaging changes are release-sensitive**: edits to `pyproject.toml` package discovery or runtime dependencies require build, wheel-content, isolated-install, and CLI smoke tests. Do not add a dependency merely to simplify an internal helper.
4. **Bundled-fork edits are supported but high-risk**: preserve upstream attribution, add a regression test for the exact failure mode, and run the matching optimizer/IRC/thermochemistry benchmark. A `[CHEMISTRY-RULE:N]` prefix is required only when an actual marked chemistry rule changes.

## See also
- Full architecture (~320 lines): [`docs/architecture.md`](../../docs/architecture.md)
- Contributor recipe + per-step gate cycle: [`CONTRIBUTING.md`](../../CONTRIBUTING.md)
- Engineering-marker coverage check: [`.github/scripts/check_engineering_markers.py`](../../.github/scripts/check_engineering_markers.py)
- Import-graph gate (no product cycle / no `core`·`domain` → `workflows` / no `pysisyphus` → product): [`.github/scripts/check_import_graph.py`](../../.github/scripts/check_import_graph.py) (tests in `tests/test_import_graph.py`)
