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
| 2. Engineering markers | `# CHEMISTRY-RULE:N` coverage, `# DOMAIN_PURE` coverage, external-library import scope, import-layer direction | `python .github/scripts/check_engineering_markers.py` and `python .github/scripts/check_import_graph.py` | a required marker is missing, an MLIP SDK is imported outside `backends/`, or an import crosses a forbidden layer |
| 3. Documentation quality | generated command references, EN/JA structure, links, and authored command contracts | `python .github/scripts/run_docs_quality.py` | public guidance or generated CLI documentation has drifted |
| 4. Help registry drift | CLI `--help` and `--help-advanced` compliance with registry | `python .github/scripts/check_help_registry.py` | CLI option mismatch — re-run after CLI changes |
| 5. Smoke | `tests/smoke/run.sh` exercises the canonical CLI surface (`extract` → `path-search` → `tsopt` → `irc` → `freq` → `all`) on a representative cluster system | copy `tests/smoke/` to scratch, then invoke `bash run.sh` from a site-specific scheduler wrapper | functional regression — root-cause before merge |

### 1.2 Before any patch

1. Read [`docs/architecture.md`](docs/architecture.md) §5 "Hidden constraints" once per session — VRAM `del` invariant, chemistry rules, repo-internal fork policy, packaging-change gates, and the `_LAZY_SUBCOMMANDS` absolute-path rule.
2. Grep [`pdb2reaction/core/defaults.py`](pdb2reaction/core/defaults.py) first for any default you are about to touch; it is the primary source for shared defaults, with a small number of justified command-local exceptions.
3. Before editing `pysisyphus/` or `thermoanalysis/`, read that directory's `README.md`, identify the local divergence you are changing, and select a focused regression test plus the matching numerical benchmark. Never replace a divergent file wholesale with upstream.
4. Identify which layer your change belongs to (`cli/`, `workflows/`, `domain/`, `backends/`, `io/`, `core/`). Stay inside one layer per commit when possible; the dependency direction is `L1 → L2 → {L3, L4} → L5` and must not be inverted.

### 1.3 `--dump` / `-v` usage examples

```bash
# Default run — INFO-level logging, no diagnostic dump
pdb2reaction all -i R.pdb P.pdb -q 0 -b uma --out-dir ./result_all

# -v LEVEL (0-3) is a per-subcommand option; -v 3 is full DEBUG
pdb2reaction all -i R.pdb P.pdb -q 0 -b uma -v 3 --out-dir ./result_all

# the freq --dump flag: write thermoanalysis.yaml alongside the standard outputs
pdb2reaction freq -i opt.pdb -q 0 --dump
```

Use `--dump` when reproducing a resource-intensive regression or attaching artefacts to a bug report. Use `-v 3` when diagnosing an import-time or stage-bridge issue; the additional log volume is acceptable for short runs.

### 1.4 Downstream output compatibility

Treat `summary.json` as a versioned public schema and `summary.log` as
human-readable output. Do not rename/remove machine-readable keys or paths
accidentally: document an intentional schema change, update both English and
Japanese references, add a compatibility/migration test when needed, and bump
the schema version for a structural break. Human-readable labels may be fixed,
but tests must verify that backend/model provenance remains accurate.

---

## 2. Project layout

See [`docs/architecture.md`](docs/architecture.md) for the full 6-layer dir tree, file index, dependency direction, and recommended reading order. Short version:

- `pdb2reaction/cli/` — L1 Interface (root group, registry, shared decorators, help, bool compatibility, resolver).
- `pdb2reaction/workflows/` — L2 Application (compute-command modules plus shared stage helpers; utility commands also live in `io/` / `domain/`).
- `pdb2reaction/domain/` — L3 Domain (chemistry-aware helpers: bond changes, bond summary, element info).
- `pdb2reaction/backends/` — L4a Infra (MLIP dispatcher + per-backend adapter).
- `pdb2reaction/io/` — L4b Infra (summary, trajectory, diagram, PDB fix, Hessian cache).
- `pdb2reaction/core/` — L5 Foundation (defaults, utilities, logging, output, and result publication).
- `pysisyphus/`, `thermoanalysis/` — bundled forks at the repo top; **not** upstream PyPI.
- `tests/` — golden gates that the per-step cycle checks.
- `tests/smoke/` — short representative job covering the canonical CLI surface.

---

## 3. Recipes

The following "add-a-X" recipes name the files to touch and the checks that
cover common contributor changes.

### 3.1 Add a subcommand

**Goal**: expose a new CLI subcommand `pdb2reaction myaction --opt1 X --opt2 Y`.

| step | action | file |
|---|---|---|
| 1 | Add a Python module `pdb2reaction/workflows/myaction.py` with a top-level `@click.command(...)` named `cli` | new file in L2 |
| 2 | Register shared/numerical defaults in `pdb2reaction/core/defaults.py`; keep a command-local default only when its scope is genuinely local and document the exception | `pdb2reaction/core/defaults.py` or the new module |
| 3 | Wire the command into the lazy registry — add `"myaction": ("pdb2reaction.workflows.myaction", "cli", "<short description>")` to `_LAZY_SUBCOMMANDS` | `pdb2reaction/cli/app.py` (L1) |
| 4 | Register every boolean option in the matching `_COMMAND_BOOL_*_OPTIONS` table so `--flag` / `--no-flag` and legacy input normalization stay covered | `pdb2reaction/cli/app.py` |
| 5 | Add a docs page `docs/myaction.md` (and `docs/ja/myaction.md` if you maintain the JP set); add a unit test in `tests/test_myaction.py` | new files |

**Gates that catch mistakes**: gate stage 1 exercises the unit/regression test;
stage 3 catches registry/help drift (including boolean registration); stage 4
must cover the new command when it belongs to the canonical smoke surface.

**Note on absolute paths**: `_LAZY_SUBCOMMANDS` entries MUST use **absolute** module paths (`pdb2reaction.workflows.myaction`). Relative dotted strings (`".myaction"`) silently break subcommand discovery if `default_group.py` ever moves; see `docs/architecture.md` §5.5.

### 3.2 Add an MLIP backend

**Goal**: introduce a new MLIP backend `XYZModel` consumable as `--backend xyz`.

| step | action | file |
|---|---|---|
| 1 | Create `pdb2reaction/backends/xyz.py` with `XYZCalculator(MLIPCalculator)` (pysisyphus path) and `XYZASECalculator(...)` (ASE path) | new file in L4a |
| 2 | Subclass `MLIPCalculator` (`backends/base.py:120`) and implement `_compute_energy_forces_ev(elem, coord_ang)`; the base supplies the pysis calculator contract and finite-difference Hessian assembly. Implement the separate ASE adapter in the backend module (see `backends/uma.py`) | `pdb2reaction/backends/base.py`, `pdb2reaction/backends/xyz.py` |
| 3 | Register in `BACKEND_REGISTRY` dict with `module / pysis_cls / ase_cls` keys, and add the accepted-kwargs set to `_BACKEND_ACCEPTED_KEYS` and `_ASE_ACCEPTED_KEYS` | `pdb2reaction/backends/__init__.py` |
| 4 | For `--backend auto`, add the `xyz` import probe to `_BACKEND_AVAILABILITY_MODULES` and add `xyz` to the `resolve_backend` fallback tuple | `pdb2reaction/backends/__init__.py` |
| 5 | Document model identifiers, install command, accepted kwargs in `docs/backends.md`; add a smoke entry in `tests/smoke/run.sh` | `docs/backends.md`, `tests/smoke/run.sh` |

**Gates that catch mistakes**: gate stage 2 confirms the external SDK import
stays inside `backends/`; gate stage 4 exercises the backend end to end after
the recipe's smoke entry (step 5) is added.

### 3.3 Add an output format

**Goal**: emit a new artefact `summary_v2.csv` alongside the existing `summary.json` / `summary.log`.

| step | action | file |
|---|---|---|
| 1 | Add a writer function in `pdb2reaction/io/summary.py` that consumes the same in-memory summary dict | `pdb2reaction/io/summary.py` (L4b) |
| 2 | Default emit path / on-or-off flag lives in `pdb2reaction/core/defaults.py` | `pdb2reaction/core/defaults.py` (L5) |
| 3 | Attach a per-subcommand `@click.option("--dump-<artefact>",...)` to the L2 stage runner when users can opt out | the stage runner |
| 4 | Update the aggregate `files` / `key_output_files` enrichment to advertise the new artefact (per `docs/json-output.md`); update the human-readable tree separately if needed | `pdb2reaction/workflows/all.py` (`_enrich_summary`), and `pdb2reaction/io/summary.py` for `summary.log` parity |
| 5 | Add docs in `docs/json-output.md` + a unit test for round-trip serialisation | `docs/json-output.md`, new test |

**Gates that catch mistakes**: the unit test in step 5 checks the writer; add
gate-stage-4 smoke coverage when the artefact belongs to the canonical smoke
surface. Any downstream-parser-visible change is governed by §1.4.

### 3.4 Add a workflow stage

**Goal**: insert a new stage (e.g. an intermediate `validate` step between TSOpt and Freq) into the `all` workflow.

| step | action | file |
|---|---|---|
| 1 | Implement the stage as a standalone subcommand first (Recipe 3.1) | `pdb2reaction/workflows/validate.py` |
| 2 | Add an internal entry to the `all` pipeline orchestrator, preserving the VRAM `del` pattern between stages | `pdb2reaction/workflows/all.py` |
| 3 | Pass the stage result explicitly to the aggregate producer and any later consumer | `pdb2reaction/workflows/all.py` |
| 4 | Record the new stage in the aggregate producer; update the human-readable writer separately if the stage belongs in `summary.log` | `pdb2reaction/workflows/all.py` (`_enrich_summary`), `pdb2reaction/workflows/path_search.py` when relevant, and `pdb2reaction/io/summary.py` for log parity |
| 5 | Update `tests/smoke/run.sh` to include the new stage in the representative run | `tests/smoke/run.sh` |

**Gate that catches mistakes**: the smoke run (step 5) — a new stage in `all` is a behaviour change and should be exercised end-to-end before merge.

### 3.5 Add a test

**Goal**: add a unit test for new behaviour or a regression test for a fixed bug.

| step | action | file |
|---|---|---|
| 1 | Pick the right tier: pure-Python or chemistry-rule logic → `tests/test_<feature>.py`; multi-stage smoke → `tests/smoke/` | as appropriate |
| 2 | Use `pytest` style: one assertion per logical thing; name the test for the symptom (`test_irc_initial_displacement_does_not_oom`) | new test |
| 3 | If the test consumes a fixture, prefer the existing `tests/data/` directory; do **not** add large binary fixtures (> 100 KB) — use generators | `tests/data/`, `tests/conftest.py` |
| 4 | Run `pytest tests/test_<feature>.py -q -x` until green, then `pytest tests/ -q` to confirm no cross-test breakage | local |
| 5 | If the test depends on a new public Click command or symbol, land Recipe 3.1 / 3.3 first so the golden gate stays green | sequencing |

**Gate that catches mistakes**: `pytest` is gate stage 1; CI blocks a failing merge.

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

IRC, TSopt, and Freq release their own large calculator, geometry, Hessian, or temporary objects at stage-specific boundaries; `all` also invokes garbage collection between child stages. Preserve those existing release/GC boundaries when changing an owner—there is no single identical `del` sequence shared by all three workflows.

### 4.3 Divergent files in bundled forks

The files listed in `pysisyphus/README.md` and
`thermoanalysis/README.md` are maintained parts of pdb2reaction. Logic edits
are allowed when they fix a demonstrated defect, but require all of the
following:

1. a regression test reproducing the defect;
2. numerical validation appropriate to the change (optimizer/IRC GPU matrix
   or thermochemistry golden data);
3. preservation of upstream attribution and the public symbols used by the
   workflows; and
4. an updated divergence table when the local behavior changes.

Do not install upstream `pysisyphus` or `thermoanalysis` alongside this
package and never overwrite a divergent file wholesale.

### 4.4 Packaging-change gate

Changes to `[tool.setuptools.packages.find].include` or runtime dependencies
are allowed only when required by the feature. Verify an sdist and wheel,
inspect wheel contents, install the wheel into a clean environment, and run
the CLI smoke suite. The existing `pdb2reaction*` glob already discovers new
subpackages below the package, so most internal files need no discovery edit.

### 4.5 `_LAZY_SUBCOMMANDS` absolute-path rule

Entries in `pdb2reaction/cli/app.py:_LAZY_SUBCOMMANDS` MUST use absolute module paths (`"pdb2reaction.workflows.all"`, never `".all"`). Relative dotted strings silently break the resolver when `default_group.py` moves; see [`docs/architecture.md`](docs/architecture.md) §5.5.

### 4.6 Chemistry default choices

Default basis set (def2-TZVPD), default functional (ωB97M-V), default convergence thresholds, default ECP handling, default solvent models — **none** of these are open for change without a `[CHEMISTRY-RULE]` commit, a documented numerical comparison, and maintainer approval. Grep `pdb2reaction/workflows/dft.py` (`DFT_DEFAULT_FUNC` / `DFT_DEFAULT_BASIS`) for the DFT functional/basis and `pdb2reaction/core/defaults.py` for the other defaults; if you think a change is justified, open an issue first.

### 4.7 Downstream-parser-visible output

Apply the versioned-schema and provenance rules in §1.4. Do not make a
machine-readable compatibility break as a side effect of wording cleanup.

---

## 5. Commit prefix conventions

The prefix identifies the expected scope and the gate cycle stage to exercise.

| prefix | meaning | typical pattern |
|---|---|---|
| `[CHEMISTRY FREEZE]` | Explicit "no chemistry change" marker on a polish-only edit; maintainers verify the scope | `[CHEMISTRY FREEZE] docstring polish on IRC.py — no logic change` |
| `[CHEMISTRY-RULE]` | Modifies a chemistry-correctness rule; requires maintainer approval and a documented scheduled numerical benchmark | `[CHEMISTRY-RULE:4] dft.py adjust rks_lowmem cutoff after gpu4pyscf update` |
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
