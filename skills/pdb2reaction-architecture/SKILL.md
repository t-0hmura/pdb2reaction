---
name: pdb2reaction-architecture
description: Where the source code lives in `pdb2reaction`. 6 physical layer directories (`cli` / `workflows` / `domain` / `backends` / `io` / `core`) + 2 repo-internal forks (`pysisyphus` / `thermoanalysis`). Tells an agent which directory to grep for a given concern (Click option, stage runner, MLIP backend, output writer, chemistry default, cap-atom math) before touching code. TRIGGER on questions like "where is X implemented", "which file defines flag Y", "how is the repo organised", "what's safe to refactor". SKIP for usage questions — those belong to `pdb2reaction-cli` / `-overview`.
---

# pdb2reaction architecture (one-screen map)

## 6 layers + bundled forks

```
pdb2reaction/                          ← the package body, one folder per layer
├── cli/        # L1 — Click root group, --help-advanced, bool normalisation,
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
│               #       (`uma.py`, `orb.py`, `mace.py`, `aimnet2.py`) +
│               #       xTB ALPB solvent correction.
├── io/         # L4b — summary writer, energy diagram, trajectory plot,
│               #       Hessian cache, PDB altloc fix. Hosts the `trj2fig`,
│               #       `energy-diagram` and `fix-altloc` subcommands.
└── core/       # L5 — `defaults.py` (primary source of truth for most
                #       CLI defaults), `utils.py` (PDB / XYZ / plot helpers),
                #       `logging.py`.

pysisyphus/        ← bundled fork of the optimiser / TS / IRC engine.
                     Slimmed to the subset pdb2reaction actually uses; its own
                     README's "Divergent files" table is the authoritative list
                     of files to keep as-is (chemistry-rule load-bearing), and
                     CONTRIBUTING.md 4.3 names the ones whose logic is frozen
                     for this release line. Only annotation edits are expected
                     in the normal workflow.

thermoanalysis/    ← bundled fork for ΔG / ZPE / partition functions.
                     QCData.py is the only consumer; same touch restriction
                     as pysisyphus.
```

Dependency direction is one-way: `L1 → L2 → {L3, L4} → L5`. The bundled forks sit outside the layer graph and may be imported from any layer.

## Where to look first

| concern | open this |
|---|---|
| Default for any CLI flag | `pdb2reaction/core/defaults.py` (primary source of truth for most defaults — grep here first; a few workflow-local defaults live inline, e.g. path-opt `--mep-mode`) |
| Subcommand body / orchestration | Canonical registry = `pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` — grep the subcommand name there first; it maps every subcommand to its module. Stage subcommands live in `pdb2reaction/workflows/<subcmd>.py` (hyphens become underscores: `path-opt` → `path_opt.py`); the 5 utility subcommands live elsewhere — `trj2fig` → `io/trj2fig.py`, `fix-altloc` → `io/pdb_fix.py`, `energy-diagram` → `io/energy_diagram.py`, `add-elem-info` → `domain/add_elem_info.py`, `bond-summary` → `domain/bond_summary.py` |
| New MLIP backend | `pdb2reaction/backends/<backend>.py` + register in `backends/__init__.py:BACKEND_REGISTRY` |
| `--help` / option decorator | `pdb2reaction/cli/common_options.py` (shared) or the subcommand file (inline) |
| Output schema (summary.json, trajectory, energy diagram) | `pdb2reaction/io/` |
| Chemistry rule (#4 gpu4pyscf `rks_lowmem`, #5 def2 auto-ECP, #7 Bofill advanced-indexing) | grep `# CHEMISTRY-RULE:` — the only three markers in the package: `workflows/dft.py` (#4, #5) and `workflows/tsopt.py` (#7); lab sign-off required to edit |
| Cap-atom (link-H) geometry | `pdb2reaction/workflows/extract.py` (`compute_linkH_atoms`) — chemistry-sensitive, and reached by opening the file (it carries no `CHEMISTRY-RULE` marker) |
| TS / IRC / optimiser internals | `pysisyphus/` (annotation-only — chemistry-rule risk) |
| MCP server / agent integration | `pdb2reaction/mcp/` — see [`pdb2reaction-mcp`](../pdb2reaction-mcp/SKILL.md) |

## Hidden constraints to remember

1. **`pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS`** entries MUST use absolute module paths (`"pdb2reaction.workflows.all"`, never `".all"`). Relative dotted paths silently break the resolver if `default_group.py` moves.
2. **VRAM hygiene**: `# DO NOT INLINE` markers around `del calc; gc.collect(); torch.cuda.empty_cache()` between stages are load-bearing — removing them OOMs the next stage on full-protein systems.
3. **`pyproject.toml [tool.setuptools.packages.find].include`** and `dependencies` arrays are treated as 0-diff for this release line. Adding a vendor / internal dir or pinning a new runtime dep breaks behaviour-level guarantees and is out of scope.
4. **Bundled-fork edits** to `pysisyphus/` / `thermoanalysis/` outside the 5 divergent files require `[CHEMISTRY-RULE:N]` commit prefix and a HEAVY benchmark.

## See also
- Full architecture (~320 lines): [`docs/architecture.md`](../../docs/architecture.md)
- Contributor recipe + per-step gate cycle: [`CONTRIBUTING.md`](../../CONTRIBUTING.md)
- Engineering-marker coverage check: [`.github/scripts/check_engineering_markers.py`](../../.github/scripts/check_engineering_markers.py)
