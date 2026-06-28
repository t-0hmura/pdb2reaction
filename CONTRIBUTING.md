# Contributing to pdb2reaction

Thank you for your interest in contributing to **pdb2reaction**.

This document is for **contributors and maintainers**. For end-user usage, see [`README.md`](README.md) and [`docs/getting-started.md`](docs/getting-started.md).

---

## 1. Before you start

`pdb2reaction` follows a **per-step gate cycle** that every change must pass before merge. Internalising this cycle prevents accidentally breaking the behaviour-level guarantees this release line carries.

### 1.1 Gate cycle

| stage | what runs | how to invoke locally | failure means |
|---|---|---|---|
| 1. Unit tests | `pytest tests/ -q` | `pytest tests/ -q` | logic regression; **never delete or skip the failing test** — root-cause it |
| 2. Engineering markers | `# CHEMISTRY-RULE:N` coverage, `# DOMAIN_PURE` coverage, external-library import scope | `python .github/scripts/check_engineering_markers.py` | a required marker is missing, or an MLIP SDK is imported outside `backends/` |
| 3. Help registry drift | CLI `--help` and `--help-advanced` compliance with registry | `python .github/scripts/check_help_registry.py` | CLI option mismatch — re-run after CLI changes |
| 4. Smoke | `tests/smoke/run.sh` exercises the canonical CLI surface (`extract` → `path-search` → `tsopt` → `irc` → `freq` → `all`) on a representative cluster system | qsub `tests/smoke/run.sh` on HPC | functional regression — root-cause before merge |

### 1.2 Before any patch

1. Read [`docs/architecture.md`](docs/architecture.md) §5 "Hidden constraints" once per session — VRAM `del` invariant, chemistry rules, repo-internal fork policy, `pyproject.toml` 0-diff arrays, `_LAZY_SUBCOMMANDS` absolute-path rule.
2. Grep [`pdb2reaction/core/defaults.py`](pdb2reaction/core/defaults.py) for any default value you are about to touch — that file is the single source of truth.
3. If you are about to edit any file under `pysisyphus/` or `thermoanalysis/` — re-read the per-dir `README.md` in that dir to confirm logic edits are forbidden in this release line (annotation-only is allowed).
4. Identify which layer your change belongs to (`cli/`, `workflows/`, `domain/`, `backends/`, `io/`, `core/`). Stay inside one layer per commit when possible; the dependency direction is `L1 → L2 → {L3, L4} → L5` and must not be inverted.

### 1.3 `--dump` / `-v` usage examples

```bash
# Default run — INFO-level logging, no diagnostic dump
pdb2reaction all -i R.pdb P.pdb -q 0 -b uma --out-dir ./result_all

# -v LEVEL (0-3) is a per-subcommand option; -v 3 is full DEBUG
pdb2reaction all -i R.pdb P.pdb -q 0 -b uma -v 3 --out-dir ./result_all

# --dump on freq: write thermoanalysis.yaml alongside the standard outputs
pdb2reaction freq -i opt.pdb -q 0 --dump
```

Use `--dump` when reproducing a HEAVY-tier regression or attaching artefacts to a bug report. Use `-v 3` when diagnosing an import-time or stage-bridge issue; the additional log volume is acceptable for short runs.

### 1.4 Downstream parser freeze rule

The contents of `summary.log` and `summary.json` are **frozen byte-for-byte**. Even cosmetic edits ("Step 1" → "step 1", added trailing periods, reflowed JSON keys) break downstream regex parsers. If a log line genuinely needs to change, treat it as an intentional behaviour change and name the reason in the commit body so downstream maintainers can update their parsers.

---

## 2. Project layout

See [`docs/architecture.md`](docs/architecture.md) for the full 6-layer dir tree, file index, dependency direction, and recommended reading order. Short version:

- `pdb2reaction/cli/` — L1 Interface (Click group, decorator factories, help, bool compat, subcommand resolver).
- `pdb2reaction/workflows/` — L2 Application (one file per subcommand stage runner).
- `pdb2reaction/domain/` — L3 Domain (chemistry-aware helpers: bond changes, bond summary, element info).
- `pdb2reaction/backends/` — L4a Infra (MLIP dispatcher + per-backend adapter).
- `pdb2reaction/io/` — L4b Infra (summary, trajectory, diagram, PDB fix, Hessian cache).
- `pdb2reaction/core/` — L5 Foundation (defaults, utils, future errors / types / logging).
- `pysisyphus/`, `thermoanalysis/` — bundled forks at the repo top; **not** upstream PyPI.
- `tests/` — golden gates that the per-step cycle checks.
- `tests/smoke/` — short representative job covering the canonical CLI surface.

---

## 3. Recipes

Five "add-a-X" recipes cover ~90 % of contributor changes. Each names the exact files you touch (with the correct layer dir) and the gate that catches mistakes.

### 3.1 Add a subcommand

**Goal**: expose a new CLI subcommand `pdb2reaction myaction --opt1 X --opt2 Y`.

| step | action | file |
|---|---|---|
| 1 | Add a Python module `pdb2reaction/workflows/myaction.py` with a top-level `@click.command(...)` named `cli` | new file in L2 |
| 2 | Register defaults (every default value) in `pdb2reaction/core/defaults.py` under a new `MYACTION_DEFAULTS` dict | `pdb2reaction/core/defaults.py` (L5) |
| 3 | Wire the command into the lazy registry — add `"myaction": ("pdb2reaction.workflows.myaction", "cli", "<short description>")` to `_LAZY_SUBCOMMANDS` | `pdb2reaction/cli/app.py` (L1) |
| 4 | If your subcommand has value-style bool flags (`--flag True`), add it to `_COMMAND_BOOL_VALUE_OPTIONS` in the same file | `pdb2reaction/cli/app.py` |
| 5 | Add a docs page `docs/myaction.md` (and `docs/ja/myaction.md` if you maintain the JP set); add a unit test in `tests/test_myaction.py` | new files |

**Gate that catches mistakes**: the unit-test suite (step 3) will fail if the subcommand cannot be discovered or instantiated; the engineering-marker check (step 4) will fail if a new backend SDK leaks outside `backends/`.

**Note on absolute paths**: `_LAZY_SUBCOMMANDS` entries MUST use **absolute** module paths (`pdb2reaction.workflows.myaction`). Relative dotted strings (`".myaction"`) silently break subcommand discovery if `default_group.py` ever moves; see `docs/architecture.md` §5.5.

### 3.2 Add an MLIP backend

**Goal**: introduce a new MLIP backend `XYZModel` consumable as `--backend xyz`.

| step | action | file |
|---|---|---|
| 1 | Create `pdb2reaction/backends/xyz.py` with `XYZCalculator(MLIPCalculator)` (pysisyphus path) and `XYZASECalculator(...)` (ASE path) | new file in L4a |
| 2 | Subclass `MLIPCalculator` (`backends/base.py:120`) and implement `_compute_energy_forces_ev(elem, coord_ang)`; the base class supplies the finite-difference Hessian assembly and ASE adapter (see `backends/uma.py:314` for the reference implementation) | `pdb2reaction/backends/base.py` |
| 3 | Register in `BACKEND_REGISTRY` dict with `module / pysis_cls / ase_cls` keys, and add the accepted-kwargs set to `_BACKEND_ACCEPTED_KEYS` and `_ASE_ACCEPTED_KEYS` | `pdb2reaction/backends/__init__.py` |
| 4 | Add `xyz` to `resolve_backend` fallback order if it should participate in `--backend auto` | `pdb2reaction/backends/__init__.py` |
| 5 | Document model identifiers, install command, accepted kwargs in `docs/backends.md`; add a smoke entry in `tests/smoke/run.sh` | `docs/backends.md`, `tests/smoke/run.sh` |

**Gate that catches mistakes**: the smoke run (step 5) will exercise the new backend end-to-end and the engineering-marker check (step 4) will confirm the new backend's external SDK import stays inside `backends/`.

### 3.3 Add an output format

**Goal**: emit a new artefact `summary_v2.csv` alongside the existing `summary.json` / `summary.log`.

| step | action | file |
|---|---|---|
| 1 | Add a writer function in `pdb2reaction/io/summary.py` that consumes the same in-memory summary dict | `pdb2reaction/io/summary.py` (L4b) |
| 2 | Default emit path / on-or-off flag lives in `pdb2reaction/core/defaults.py` | `pdb2reaction/core/defaults.py` (L5) |
| 3 | Wire into `@add_common_dump_options` factory if the user can opt out (the factory is **future**, expected to land in a later release; until then, attach a per-subcommand `@click.option("--dump-<artefact>",...)` directly to the L2 stage runner) | `pdb2reaction/cli/decorators.py` (future) |
| 4 | Update the `files` map written into `summary.json` to advertise the new artefact (per `docs/json-output.md`) | `pdb2reaction/io/summary.py` |
| 5 | Add docs in `docs/json-output.md` + a unit test for round-trip serialisation | `docs/json-output.md`, new test |

**Gate that catches mistakes**: the smoke run (step 5) will exercise the new artefact end-to-end; any change to a downstream-parser-visible log line is governed by §1.4.

### 3.4 Add a workflow stage

**Goal**: insert a new stage (e.g. an intermediate `validate` step between TSOpt and Freq) into the `all` workflow.

| step | action | file |
|---|---|---|
| 1 | Implement the stage as a standalone subcommand first (Recipe 3.1) | `pdb2reaction/workflows/validate.py` |
| 2 | Add an internal entry to the `all` pipeline orchestrator, preserving the VRAM `del` pattern between stages | `pdb2reaction/workflows/all.py` |
| 3 | Add a `_StageContext` field if the stage needs persistent context (future, `core/types.py`) | `pdb2reaction/core/types.py` (future) |
| 4 | Update `pdb2reaction/io/summary.py` to record the new stage's entry in `summary.json` | `pdb2reaction/io/summary.py` |
| 5 | Update `tests/smoke/run.sh` to include the new stage in the representative run | `tests/smoke/run.sh` |

**Gate that catches mistakes**: the smoke run (step 5) — a new stage in `all` is a behaviour change and should be exercised end-to-end before merge.

### 3.5 Add a test

**Goal**: add a unit test for new behaviour or a regression test for a fixed bug.

| step | action | file |
|---|---|---|
| 1 | Pick the right tier: pure-Python logic → `tests/test_<feature>.py`; multi-stage smoke → `tests/smoke/`; chemistry-rule regression → `tests/domain_golden/` (future) | as appropriate |
| 2 | Use `pytest` style: one assertion per logical thing; name the test for the symptom (`test_irc_initial_displacement_does_not_oom`) | new test |
| 3 | If the test consumes a fixture, prefer the existing `tests/data/` directory; do **not** add large binary fixtures (> 100 KB) — use generators | `tests/data/`, `tests/conftest.py` |
| 4 | Run `pytest tests/test_<feature>.py -q -x` until green, then `pytest tests/ -q` to confirm no cross-test breakage | local |
| 5 | If the test depends on a new public Click command or symbol, land Recipe 3.1 / 3.3 first so the golden gate stays green | sequencing |

**Gate that catches mistakes**: `pytest` itself (step 3 of the gate cycle); CI will block merge.

---

## 4. Do not touch list

These are **hard constraints** enforced by the release process. Violating them either breaks correctness (chemistry rules), behaviour-level guarantees (`pyproject.toml` arrays, downstream-parser log lines), or upstream-fork compatibility.

### 4.1 Chemistry rules

The reaction-path correctness rules listed in [`docs/architecture.md`](docs/architecture.md) §5.1 must not be reordered, simplified, or factored out. They are marked with `# CHEMISTRY-RULE:N` inline comments and `# DOMAIN_PURE` module-docstring markers. The CI gate `.github/scripts/check_engineering_markers.py` enforces marker completeness and confines MLIP-only SDK imports (`fairchem`, `orb_models`, `mace`, `aimnet`) to the `backends/` layer. The three rules are #4 (gpu4pyscf `rks_lowmem`) and #5 (def2 auto-ECP) in `workflows/dft.py`, and #7 (`bofill_update` advanced-indexing) in `workflows/tsopt.py`.

Use the grep recipe before any patch:

```bash
grep -rnE '# CHEMISTRY-RULE:[0-9]+' pdb2reaction/
grep -rn '# DOMAIN_PURE' pdb2reaction/
```

### 4.2 VRAM-management invariant (`del` chains)

The IRC / TSopt / Freq stages explicitly `del calc`, `del geom`, `del hess` between stages and the `all` workflow runs `gc.collect()` at stage boundaries. **Do not refactor those `del` / `gc.collect()` statements out** — long-running jobs OOM without them.

### 4.3 Divergent files in bundled forks

Logic edits to these files are forbidden in this release line (annotation-only is allowed: docstring + type hints):

- `pysisyphus/irc/IRC.py`
- `pysisyphus/optimizers/hessian_updates.py`
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py`
- `thermoanalysis/QCData.py`

Do not `pip install pysisyphus` or `pip install thermoanalysis` from PyPI alongside this package — silent runtime breakage.

### 4.4 `pyproject.toml` arrays are 0-diff

`pyproject.toml [tool.setuptools.packages.find].include` and `dependencies` arrays must not change in this release line. Adding a `vendor/` or `internal/` container directory, or pinning a new runtime dependency, breaks the install contract and is forbidden by the release scope. Reflow / comment edits to `pyproject.toml` are fine; **array contents** are frozen. (The `pdb2reaction*` glob already auto-discovers any new layer subpackage — no `include`-array edit is needed when adding a file inside an existing layer dir.)

### 4.5 `_LAZY_SUBCOMMANDS` absolute-path rule

Entries in `pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` MUST use absolute module paths (`"pdb2reaction.workflows.all"`, never `".all"`). Relative dotted strings silently break the resolver when `default_group.py` moves; see [`docs/architecture.md`](docs/architecture.md) §5.5.

### 4.6 Chemistry default choices

Default basis set (def2-TZVPD), default functional (ωB97M-V), default convergence thresholds, default ECP handling, default solvent models — **none** of these are open for change without a `[CHEMISTRY-RULE]` commit and explicit lab decision. Grep `pdb2reaction/workflows/dft.py` (`DFT_DEFAULT_FUNC` / `DFT_DEFAULT_BASIS`) for the DFT functional/basis and `pdb2reaction/core/defaults.py` for the other defaults; if you think a change is justified, open an issue first.

### 4.7 Downstream-parser-visible log lines

Any `summary.log` or `summary.json` line that downstream parsers consume is **frozen byte-for-byte**. See §1.4 above.

---

## 5. Commit prefix conventions

The prefix tells the reviewer what to expect and which gate cycle stage will be exercised.

| prefix | meaning | typical pattern |
|---|---|---|
| `[CHEMISTRY FREEZE]` | Explicit "no chemistry change" marker on a polish-only edit; reviewer must verify | `[CHEMISTRY FREEZE] docstring polish on IRC.py — no logic change` |
| `[CHEMISTRY-RULE]` | Modifies an actual chemistry-correctness rule — requires lab sign-off + HEAVY benchmark | `[CHEMISTRY-RULE:4] dft.py adjust rks_lowmem cutoff after gpu4pyscf update` |
| `[DOMAIN_PURE]` | Adjusts the `# DOMAIN_PURE` marker or the import-deny gate | `[DOMAIN_PURE] add tsopt.py to deny-gate scope` |

---

## 6. Where to ask

| forum | best for |
|---|---|
| GitHub issue | reproducible bug, feature request, design question with a concrete proposal |
| GitHub discussion | open-ended design / chemistry-method discussion, "what is the right way to..." |

---

## License

By contributing, you agree that your contributions will be licensed under the
[GPL-3.0 License](LICENSE).
