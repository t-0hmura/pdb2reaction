# Changelog

All notable changes to **pdb2reaction** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## Unreleased

## [0.4.4] — 2026-07-06

### Changed
- **Behavior change (default): the default TS optimizer is once again RS-I-RFO** (Restricted-Step
  Image RFO), reverting the RS-P-RFO default introduced in 0.4.2. The 0.4.2 choice was based on a
  benchmark that predates the 0.4.1 charge/spin fix; on the corrected (charge-honoring) benchmark,
  RS-I-RFO produces clean first-order saddles at least as reliably as RS-P-RFO across the reliable
  backends, and it keeps the default consistent with the mlmm_toolkit microiteration path.
  `hess`/`heavy` resolve to `rsirfo` again; RS-P-RFO stays available via `--opt-mode rsprfo` and
  TRIM via `--opt-mode trim`.
- **Behavior change (default): the default UMA model is now `uma-s-1p2`** (was `uma-s-1p1`). At the
  same small-model cost it is more robust on the benchmark (fewer optimization/frequency errors and
  a few more clean saddles). Other models (`uma-s-1p1`, `uma-m-1p1`, MACE-OMOL, Orb-v3-omol) remain
  selectable via `-b` / `--backend-model` / config.
- **Centralized the default UMA model** in a single constant `DEFAULT_UMA_MODEL`
  (`pdb2reaction/core/defaults.py`); the `uma-s-1p1` defaults previously hardcoded across backends
  and workflows now all reference it. Docstrings, the `--opt-mode` CLI help, the summary-log legend,
  docs, and skills were updated for consistency with both default changes.

## [0.4.3] — 2026-07-06
### Fixed
- **ORB analytical Hessian failed with a donated-buffer error on some environments.**
  `OrbCalculator._compute_analytical_hessian_ev` builds the Hessian via a double backward
  (`torch.autograd.functional.hessian`). Where Orb's conservative model / the torch aot_autograd
  stack enables the donated-buffer optimization, the second backward raised *"This backward
  function was compiled with non-empty donated buffers which requires create_graph=False and
  retain_graph=False."* The precision guard (`float32-high`/`float64`) was **not** sufficient — the
  optimization fires even at float64 on affected torch/orb builds. Now explicitly sets
  `torch._functorch.config.donated_buffer=False` for the duration of the Hessian and restores it
  afterward. It is a no-op where the optimization is not enabled, and fixes
  `--hessian-calc-mode Analytical` with `-b orb`. The failure is environment-dependent (it appears
  only with certain torch/orb builds; it does not reproduce with torch 2.8.0 / orb 0.5.5).

## [0.4.2] — 2026-07-05

### Changed
- **Behavior change (default):** the default TS optimizer is now **RS-P-RFO** (Restricted-Step
  Partitioned RFO, Banerjee), changed from RS-I-RFO. This affects `tsopt --opt-mode hess` (the
  default) and the `all` TSOPT / post-IRC stage. At equal wall-time RS-P-RFO converges to a clean
  first-order saddle (single imaginary mode) more robustly on backends where RS-I-RFO tends to land
  on high-order saddles. RS-I-RFO remains available via `--opt-mode rsirfo`. The `hess`/`heavy`
  aliases now resolve to `rsprfo`; `rsirfo` is a distinct explicit alias. All three RFO-family
  optimizers still share the `rsirfo` YAML block. Docs updated throughout.

## [0.4.1] — 2026-07-05

### Fixed
- **Charge/spin were silently dropped on the pysisyphus MLIP path.** `AtomicData.from_ase`
  was called without `r_data_keys`, so fairchem-core ≥2.x ran every UMA calculation at
  charge=0/spin=0 regardless of `-q`/`-m` (a regression against older fairchem, which read
  `atoms.info` unconditionally). Now passes `r_data_keys=["spin","charge"]`: charge is honored
  and `spin` is the spin multiplicity (2S+1). The ASE/DMF path was unaffected.
- `--opt-mode trim` crashed with frozen atoms (partial Hessian); TRIM now reduces the gradient
  to and expands the step from the active subspace, like RS-I-RFO / RS-P-RFO.
- The DMF path optimizer was pinned to fp32 while the HEI was ranked at the requested precision;
  `create_ase_calculator` now honors `--precision`. `--mep-mode dmf` now errors under `--solvent`
  (the ASE path has no implicit-solvent wrapper) and points to `--mep-mode gsm`.

### Changed
- Documentation/CLI corrections to match the code (tsopt `--opt-mode` aliases include
  `trim`/`rsprfo`; MACE install note; `create_ase_calculator` kwargs; `--radius-het2het` help).

## [0.4.0] — 2026-06-28

### Changed
- **Behavior change (default):** `all --refine-path` now defaults to `False`.
  The `all` pipeline's MEP stage runs a single-pass `path-opt` by default; pass
  `--refine-path True` to run the recursive `path-search` (automatic multi-step
  bond-change segmentation), which was the previous default. The default MEP
  work directory is now `_work/path_opt/` (was `_work/path_search/`). The
  standalone `path-search` subcommand is unchanged. Docs/skills updated
  throughout to reflect the new default.
- `--dft-func-basis` is now surfaced in the primary `pdb2reaction <subcmd>
  --help` (previously only under `--help-advanced`), so the DFT//MLIP
  functional/basis is discoverable without the advanced listing.
- `--precision fp64` now also forces the Hessian to fp64 (`hessian_double`)
  so the optimiser / eigen linear algebra cannot silently run in a lower
  precision than the model; a config that set `hessian_double=False` under
  fp64 is overridden with a warning.
- AIMNet2 now rejects both `--precision fp64` (its model inputs are cast to
  float32 upstream) and `--deterministic` (its forces come from a custom
  CUDA kernel outside torch's deterministic-algorithms control), with clear
  errors instead of running misleadingly.
- The `all`-pipeline determinism check in the smoke suite is now an
  informational monitor (reports drift, does not gate); bit-exactness is an
  opt-in via `--deterministic`, not a default guarantee.
- `--help` of `pdb2reaction` now groups subcommands under semantic
  sections ("Pipelines" / "Pipeline stages" / "Inputs & topology" /
  "Analysis") in a configurable, deterministic order; subcommands not
  listed in any section fall through to a trailing "Other" bucket so
  we never hide an entry silently.
- CLI exception renderer appends `Try 'pdb2reaction <subcmd> -h for
  help.` to every user-input-style error so first-time users see a
  recovery path, and routes the full traceback through
  `logging.getLogger(...).exception` so log scrapers / `-v` users get
  the structured record alongside the human-readable terminal echo.
- Repo-wide ruff `F401` sweep: removed 29 unused imports including the
  unused `warnings` shadow in `backends/__init__.py`.
- `docs/cli-conventions.md` now spells out the four permanent boolean
  forms (`--flag` / `--no-flag` / `--flag True/yes/1/on` /
  `--flag False/no/0/off`) and adds a "Contributing a new bool flag"
  section pointing at the `add_*_option` factory + `_COMMAND_BOOL_*`
  registries.
- `docs/architecture.md` lists the actual `pdb2reaction/workflows/`
  files (the previous text referenced workflow modules that do not
  exist in this repo) and refreshes the per-file LOC numbers cited
  in §2.3.

### Fixed
- Bond-change detection (`domain/bond_changes.compare_structures`) is now
  row-chunked instead of building dense N×N distance matrices, removing a CUDA
  out-of-memory failure on large solvated clusters (~20k+ atoms) during
  `path-search` / `scan` kink detection on 16–24 GB GPUs.
- OPC / TIP4P 4-point water with a virtual site (Amber `EPW`, element `EP`) is
  now tolerated by the structure readers instead of crashing on the massless
  extra point.
- MCP `run_single_point_dft` emitted `--functional` / `--basis` (two flags
  the `dft` CLI does not accept) so every call failed; it now passes the
  single combined `--func-basis FUNC/BASIS` argument.
- `MLIPCalculator.__init__` no longer swallows malformed `freeze_atoms`
  input with `except Exception: freeze_iter=[]`. A typo in the CLI
  string used to silently un-freeze every atom and let TS opt walk
  through what the user thought was frozen. Narrowed to (TypeError,
  ValueError) and raises with the offending value.
- `MLIPCalculator._compute_full_hessian_au` referenced `torch` but the
  module never imported it; any analytical-Hessian backend that
  returned a torch tensor would crash with `NameError: torch is not
  defined`. Added the missing `import torch`.
- `path_search._run.freeze_atoms_for_log` pre-initialised before its
  inner try block; the comprehension is also assigned later in the
  nested function, so Python treats the name as local and the previous
  raise hit `UnboundLocalError` instead of falling back to the empty
  list from the cli() scope.
- `workflows/all.py` GPU-release block before the DFT subprocess fork
  rebinds names to `None` explicitly instead of calling
  `del locals()[name]` — the latter is a CPython no-op (locals()
  returns a copy of the frame namespace) and the torch.nn.Module
  references stayed pinned through the subsequent gc.collect() +
  empty_cache(), so the subprocess started with the GPU still occupied.

### Added
- `--calc-file PATH` (with `--calc-factory NAME`): load an arbitrary ASE
  Calculator from a user Python file as a `custom` backend — usable on every
  subcommand and forwarded through the `all` pipeline. Couple GFN-xTB, DFTB+,
  ORCA, or any ASE-compatible engine without modifying the package; energy /
  forces follow the ASE eV / eV·Å⁻¹ contract and Hessians use the
  finite-difference path. See `docs/backends.md`.
- `--backend-model NAME` flag on every backend-using subcommand (`opt`,
  `tsopt`, `freq`, `irc`, `scan` / `scan2d` / `scan3d`, `path-opt`,
  `path-search`, `sp`, `all`) to override the model variant for the selected
  `--backend` (e.g. `--backend uma --backend-model uma-s-1p2`), routed to the
  backend's `model` kwarg. Previously settable only through `--config` YAML.
- `--deterministic` flag on every compute subcommand (`opt`, `tsopt`,
  `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`,
  `all`, `sp`) for bit-reproducible GPU runs (deterministic algorithms +
  an `index_reduce_` shim). Process-global, slower, and fails loud if the
  build cannot honour it; `PDB2REACTION_STRICT_DETERMINISTIC=1` is the env-var
  equivalent. Verified bit-identical energy and forces on uma / orb / mace.
- `docs/reproducibility.md` documenting the determinism / precision model
  and the per-backend reproducibility guarantees.
- `tests/test_help_grouping.py` locks the four-bucket `--help` section
  rendering + order.
- `result.json` / `summary.json` envelope now carries
  `schema_version: "1.0"` and `write_result_json` mirrors every per-stage
  `result.json` payload to a sibling `summary.json` so agents can
  converge on a single filename across every subcommand. `result.json`
  is preserved for back-compat. `RESULT_JSON_STATUS_VALUES` enumerates
  the allowed `status` strings.
- Structured error envelope when a subcommand fails: the JSON envelope
  now carries `error_class_chain` (MRO names), `error_module`, and
  `error_label` alongside the legacy `error` / `error_type` / `status`
  keys, so MCP clients can pattern-match the exception hierarchy
  without parsing text.
- `docs/output-layout.md` (new): single-page reference for the filename
  conventions per subcommand + agent recipe for reading `summary.json`
  with class-chain pattern matching. `docs/json-output.md` and
  `docs/mcp_server.md` updated with the schema_version, summary.json
  mirror, and error envelope semantics.
- `tests/test_write_result_json.py`, `tests/test_error_envelope.py`
  (~6 tier-1 assertions) lock in the new envelope contracts.
- `pdb2reaction.workflows._all_helpers` (new module) provides the
  landing zone for `pdb2reaction.workflows.all.cli()` decomposition.
  Currently exposes:
    * `build_energy_level_dict` (factors the 4-way R/TS/P energy-level
      dict duplication out of cli() — same shape for UMA / Gibbs /
      DFT / Gibbs-DFT entries, now in one place).
    * `build_pipeline_summary_payload` (factors the inner body of
      `_write_pipeline_summary_log` for unit-testability).
- `tests/test_all_helpers.py` (5 cases) pins the extracted helper
  contracts (energy-level kcal projection + no-input-mutation, summary
  payload shape, DFT-disabled path drops dft_func_basis, AllContext
  signature drift guard).
- `AllContext` frozen dataclass added to `_all_helpers`: bundles the
  65 `pdb2reaction all` CLI parameters in declaration order; foundation
  for future cli() decomposition steps.
- `pdb2reaction.workflows._path_yaml_helpers` (new module):
  `apply_single_opt_yaml_layer` extracted from the 44-LOC nested
  closure that previously lived (verbatim) in both `path_opt.cli()`
  and `path_search.cli()`. Both sites now delegate to the shared
  helper.
- `pdb2reaction.mcp._runner` exposes `SubcmdResultDict` (TypedDict),
  `MCP_SUBCMD_RESULT_SCHEMA_VERSION = "1.0"`, and
  `MCP_SUBCMD_RESULT_STATUSES` enum; `SubcmdResult.to_dict` now emits
  `schema_version` so MCP clients can pin the contract.
- Smoke `tests/smoke/run.sh` expanded with per-stage `--coord-type
  {dlc,redund,tric}` + `--precision fp64` test coverage
  (test53a/d/g/j/k/m); test53 itself capped at `--max-cycles 5 --no-
  tsopt/thermo/dft` so the DLC code path is exercised without the
  multi-hour convergence the uncapped run requires.
- Smoke model-based complex scans now pass explicit boundary
  `--freeze-atoms` for the generated `p_complex_model.pdb`, and add
  complex-model DLC regressions (test69/70) so frozen-DOF internal-coordinate
  paths are not covered only by small-molecule inputs.
- `--precision fp32|fp64` accepted on every calculator-constructing
  subcommand. The flag was previously available only on `tsopt / freq /
  irc / sp`; it now also covers `opt / all / path-opt / path-search /
  scan / scan2d / scan3d`. For `all`, the value propagates to every
  child stage through the shared args YAML, so a single top-level
  switch covers the full pipeline.
- `--irc-pos-def` (IRC convergence guard requiring PSD mass-weighted
  Hessian) is opt-in on `irc`; blocks the IRC "shoulder" false
  convergence where the rms-only criterion calls success before
  reaching the local minimum.

### Removed
- **BREAKING:** Flat-top compatibility shim layer removed. The package now
  lives under 6 layer directories (`cli/`, `workflows/`, `domain/`,
  `backends/`, `io/`, `core/`); the shims at `pdb2reaction/<file>.py` that
  re-exported the new locations have been deleted in this release. External
  code must migrate dotted imports to the layered paths:

  | Old (removed)              | New                                |
  |----------------------------|------------------------------------|
  | `pdb2reaction.{all,opt,tsopt,freq,irc,scan,scan2d,scan3d,path_opt,path_search,extract,dft}` | `pdb2reaction.workflows.<same>` |
  | `pdb2reaction.align_freeze_atoms` | `pdb2reaction.workflows.align_freeze` |
  | `pdb2reaction.scan_common` | `pdb2reaction.workflows.scan_common` |
  | `pdb2reaction.harmonic_constraints` | `pdb2reaction.workflows.restraints` |
  | `pdb2reaction.{defaults,utils}` | `pdb2reaction.core.<same>`     |
  | `pdb2reaction.uma_pysis`   | `pdb2reaction.backends.uma`        |
  | `pdb2reaction.{bond_changes,bond_summary,add_elem_info}` | `pdb2reaction.domain.<same>` |
  | `pdb2reaction.{energy_diagram,trj2fig,hessian_cache}` | `pdb2reaction.io.<same>` |
  | `pdb2reaction.fix_altloc`  | `pdb2reaction.io.pdb_fix`          |
  | `pdb2reaction.summary_log` | `pdb2reaction.io.summary`          |
  | `pdb2reaction.cli_utils`   | `pdb2reaction.cli.decorators`      |
  | `pdb2reaction.{bool_compat,default_group}` | `pdb2reaction.cli.<same>` |
  | `pdb2reaction.advanced_help` | `pdb2reaction.cli.help_pages`    |

  The `pdb2reaction` console-script CLI is unaffected — only Python imports change.
- `--trust-band` / `--hessian-window` / `--weighted-trust` CLI flags
  (and their `add_*_option` factories). The trust-radius / multistep
  Hessian-update knobs they exposed showed no benefit on small TS
  benchmarks and actively slowed convergence (rho-band trust update
  −33 %, hessian_window > 1 −47 % cycles on a 20-atom TS), with no
  evidence of speed-up on production-scale systems. The vendored
  pysisyphus `HessianOptimizer` kwargs are left dormant; no
  behaviour change since defaults were always legacy.

### Documentation
- Documented that the `[orb]` extra's `torch_scatter` has no PyPI binary wheel
  (sdist only, fails under PEP517 build isolation): install from PyG's
  prebuilt-wheel index, e.g.
  `pip install "pdb2reaction[orb]" -f https://data.pyg.org/whl/torch-2.8.0+cu129.html`.
- Switched the README overview image to an absolute URL so it renders on PyPI.
- Cleanup pass on in-source comments and docs: removed personal-name
  attribution, internal-channel references, private memo filenames, and
  internal review-process markers (phase IDs, dated user-correspondence
  tags) that had leaked into shipped sources. Technical rationale is
  preserved verbatim — no runtime behaviour change.

## [0.3.10] — 2026-05-17

### Fixed
- `scan` / `scan2d` config precedence now matches the other subcommands: `defaults < --config (YAML) < CLI`. Previously `build_scan_configs` applied the YAML configuration *after* the CLI-derived values, so a `--config` file silently overrode explicit CLI options (e.g. `--thresh`, `--bias-k`, `--workers`) for scans. Runs that pass options only on the CLI (or only via YAML) are unaffected. Added `tests/test_scan_precedence.py`.

### Documentation
- Corrected the `--workers` help string across all subcommands (and the generated CLI reference): `>1` does not make Hessian computation unsupported — the analytical (autograd) Hessian is unavailable and pdb2reaction silently uses the FiniteDifference Hessian instead.
- `energy-diagram` docs (EN/JP): the renderer is Plotly, not Matplotlib.
- `quickstart-all` (JP): removed the nonexistent `--irc` toggle; IRC validation runs automatically as part of `--tsopt` (matches the EN page).

### Changed
- AI-agent skill bundle moved from `.claude/skills/` to top-level `skills/` so non-Claude agents (Codex, Cursor, aider, …) can read the same instructions. Copy the directory into your project (e.g.\ as `.claude/skills/` for Claude Code) to activate. README / docs / drift-check scripts updated.

## [0.3.9] — 2026-05-11

Default-value alignment plus removal of the `--resume` flag.

### Removed
- **BREAKING:** `pdb2reaction all --resume / --no-resume`. The resumed-run path silently dropped per-segment TSOPT energies, UMA reference energies, freq thermal corrections, and DFT results from aggregate diagrams and `summary.json`; sentinel checks were existence-only (no integrity verification, no parameter-identity comparison) and the TS-only branch consumed prior outputs without checking the input PDB. To pick up a walltime-truncated run, invoke the standalone subcommands (`pdb2reaction tsopt / irc / freq / dft`) against the segment outputs `all` already produced.

### Changed
- `SEARCH_KW.max_nodes_segment` 10 → 20.
- `path-search`/`path-opt` `--preopt` default `False` → `True` (matches `all`).
- `path-opt` `--fix-ends` default `False` → `True` (matches `GS_KW`).
- `MLIPCalculator` (and all backend subclasses) `return_partial_hessian` / `out_hess_torch` class-kwarg defaults `False` → `True` (matches `CALC_KW_DEFAULT`).
- Revised documents (EN+JA).

### Fixed
- `path-search` writes all stage artefacts (stopt/lbfgs/rfo) to `-o` instead of leaking preopt to `./result_opt/`.
- `all` no longer overrides YAML `stopt.max_cycles` / `gs.climb` when user did not pass the corresponding flag.
- `OrbCalculator.__init__` exposes `out_hess_torch` explicitly (parity with MACE/AIMNet2).

## [0.3.8] — 2026-05-01

### Added
- `pdb2reaction dft --lowmem/--no-lowmem` (default `True`): closed-shell GPU SCF now uses `gpu4pyscf.dft.rks_lowmem.RKS` (direct-JK, no density fitting). Open-shell / CPU / older `gpu4pyscf` paths fall back to standard `RKS`/`UKS`. YAML key `dft.lowmem`.

### Changed
- Closed-shell GPU DFT defaults switch from DF + standard `RKS` to direct-JK `rks_lowmem.RKS`. Absolute energies shift by sub-mHa; pass `--no-lowmem` to reproduce v0.3.6/v0.3.7 numbers.
- Documentation pruned (EN+JA): per-command pages drop redundant Summary↔At a glance↔intro repetitions, version-stamp annotations are moved to this CHANGELOG, and decorative `---` separators between H2 headings are removed (~130 across docs/).

### Fixed
- `OrbASECalculator` default precision changed from `'float32'` (silent slow path, blocks autograd Hessian) to `'float32-high'`, matching `OrbCalculator`.
- `dft.py` exception handler: `out_dir_path` is now pre-bound before YAML override resolution, so `apply_yaml_overrides` failures no longer mask the original exception with `NameError`.
- `freq.py` analytical-Hessian path now clones the result before mass-weighted projection, preventing in-place mutation of the cached Geometry Hessian.
- `freq.py` `prepared_input` cleanup span widened so the temp file is removed even when an exception fires before the main try block.
- `solvent.py`: when the wrapped backend returns a torch GPU tensor for the Hessian, it is detached to CPU/numpy before adding the xTB delta correction (previously raised on UMA + analytical + GPU + solvent).
- `path_opt.py` finally block: `shared_calc = gs = geoms = None` indent corrected so it runs once outside the cleanup loop.
- `summary_log.py`: `nu_imag` resolution no longer treats `0.0` as missing (was `or`-chain falsy fallback); `seg_NN` tag fallback width unified to `:03d`.
- Documentation realigned with source: `uma-pysis.md` `return_partial_hessian` default is `False` (not `True`); `scan2d` output filenames are distance-tag (`point_iDDD_jDDD`, `DDD = round(d × 100)` Å), not step indices; `all.md` output tree expands `models/model_<input>.pdb`, `freq/{R,TS,P}/`, and `tsopt_single/` subdirectories; `path-opt.md` lists `final_geometries.pdb` (no `_trj` suffix) and adds `.gjf` companion; `path-search.md` adds `mep_seg_NN`, `hei_seg_NN`, `hei_w_ref_seg_NN.pdb`, `mep.gjf` to the output tree; JA `troubleshooting.md` charge-mapping example now includes the missing `extract` subcommand; JA `getting-started.md` adds the missing Agent Skills section.
- Internal cleanup: removed unused `CliArgBuilder` class and dead YAML-merge trampolines in `all.py`; removed legacy argparse `main()` in `add_elem_info.py`; made the `defaults.py` module docstring self-contained; `extract.py` `if api==True:` → `if api:`.
- Tests: `test_summary_log.py` no longer skips silently when `fairchem` is absent; unused imports removed from 4 test files.

## [0.3.7] — 2026-04-28

### Changed
- Docs restructured for newcomer onboarding (EN+JA): goal-based 4-card "Start here" map on `docs/index.md`; unified 5-bullet "At a glance" block on all 11 calculation command pages; unified Quickstart template; `getting-started.md` workflow-modes table + flag matrix trimmed; `installation.md` split into Required vs Optional with the MACE+UMA conflict promoted to a `{warning}`; `dft.md` size/OOM caveats consolidated into a single "Practical limits" subsection.
- Reference duplication eliminated — each canonical fact has exactly one home and is cross-linked from every other page: `--selected-resn` ID-vs-name (`cli-conventions.md` `(selected-resn-takes-ids)`), `--opt-mode` polysemy (`(opt-mode-semantics)`), `--engine` / `--dft-engine` alias (`(engine-vs-dft-engine)`), exit-codes table (`(exit-codes)`), 5 cm⁻¹ vs 100 cm⁻¹ TS thresholds (`glossary.md` `(imaginary-mode-thresholds)`), scan-list spec (`(scan-list-spec)`).
- Command pages (`tsopt`, `path-search`, `path-opt`, `opt`, `freq`, `irc`) no longer reproduce the YAML schema; canonical schema lives in `yaml-reference.md`. ~1500 lines of duplication removed; only command-specific overrides remain inline.
- `recipes-common-errors.md` is now a symptom→page router; detailed fixes live in `troubleshooting.md`.

### Added
- `glossary.md` (EN+JA): `DFT//MLIP` entry.
- `yaml-reference.md` overview: `rsirfo` "Used by" lists `tsopt, all` (was `tsopt`).
- `getting-started.md` "Important CLI options": `--opt-mode-post grad|hess` (`all`-only).

### Fixed
- `docs/ja/extract.md`: typo `:electric` → `:charge`.
- `docs/ja/cli-conventions.md` `--selected-resn`: previously claimed silent no-match; corrected to match the actual `ValueError` raised when residue-name tokens are passed.
- `docs/yaml-reference.md` IRC: removed redundant `calc.return_partial_hessian` bullet duplicating the preceding line (EN+JA).
- `docs/ja/getting-started.md`: corrected relative path to `reference/commands/index.md`.

## [0.3.6] — 2026-04-21

### Added
- GPU-resident analytical Hessian for all four backends (Orb, MACE, AIMNet2 in addition to UMA); previously only UMA provided a native analytical Hessian and the other backends silently fell back to finite differences when `--hessian-calc-mode Analytical` was requested. The backend-level silent fallback is removed — Orb / MACE / AIMNet2 now either produce an analytical Hessian or raise `BackendError`. The long-standing worker-level downgrade (UMA multi-worker path uses finite differences regardless of `hessian_calc_mode`) is unchanged and remains documented in `docs/uma-pysis.md` and `docs/yaml-reference.md`.

### Changed
- Orb backend default precision: `float32` → `float32-high` (higher-precision matmul on Ampere+).
- Documentation: EN/JA synchronized (analytical-Hessian backend coverage, Orb description, MyST cross-references) across `README.md`, `docs/*.md`, and `docs/ja/*.md`. Sphinx HTML now builds with zero warnings.
- Output tree: per-segment `tsopt/` → `ts/`; `structures/` subdirectory added under each `post_seg_NN/`.

### Removed
- `examples/benchmark/` and `scripts/validate_benchmark.py` / `scripts/validate_summary.py`. The 6-enzyme / 23-step cluster-model benchmark now ships as a separate Zenodo data bundle, not as part of the software repository.

### Fixed
- Zenodo DOI typo in `README.md` and `CITATION.cff`: `10.5281/zenodo.19197878` (unrelated record by another author) → `10.5281/zenodo.19197865` (pdb2reaction concept DOI).
- TS optimization: reverted the TR (translation/rotation) projection that destabilized convergence on link-hydrogen-capped clusters.
- Orb backend description in `README.md` and `troubleshooting.md` (EN + JA): the old "higher failure rate / SVD failures" wording did not describe the current post-analytical-Hessian behavior; reworded to "correctly identifies the reaction coordinate but TS typically carries extra small imaginary modes".
- Sphinx cross-reference warnings: several `file.md#anchor` call sites converted to `{ref}...<label>` form.
- `tests/smoke/test.md`: test count 35 → 41 (run.sh has `test1` .. `test41`); rows 36–41 and the dry-run block realigned with actual indices.

### Upgrade notes
- Users who relied on `examples/benchmark/` or `scripts/validate_*` should pull the benchmark set from the separate Zenodo data bundle, or keep a copy of the 0.3.5 tarball.
- Runs that implicitly depended on `--hessian-calc-mode Analytical` silently falling back to finite differences on Orb / MACE / AIMNet2 will now compute true analytical Hessians. Set `--hessian-calc-mode FiniteDifference` explicitly to restore the old behavior.

## [0.3.5] — 2026-04-13

### Added
- Energy plateau convergence fallback for optimizers stuck on flat PES (range-based criterion, threshold 1e-4 au).
- `_all_mw_freqs_cm` helper for TS imaginary-mode tracking (commented debug prints).

### Changed
- `--refine-path` default reverted to `True` (recursive `path-search` is again the primary MEP mode under `pdb2reaction all`).
- `trust_max` lowered from 0.20 to 0.10 for RFO / RS-I-RFO optimizers for MLIP stability.
- Energy-plateau convergence criterion switched from mean to range with threshold 2e-4 au, and skipped for chain-of-states optimizers.
- DFT stage now runs as a subprocess to avoid libcusolver conflicts with PyTorch; releases calculator/result refs and frees GPU memory before the subprocess.
- DFT output shows GPU device name and CPU thread count.
- Doc version bumped to v0.3.5; `def2-svp` recommended as OOM workaround for large systems.

### Fixed
- PDB trajectory conversion: `MODEL`/`ENDMDL` missing on first frame.
- Preopt output directory overwrite across segments (include segment name in tag).
- Error `result.json` now written on every CLI subcommand failure.
- Graceful handling of DFT failure (skip diagrams, expose status in summary).
- EulerPC corrector integration loop: safety guards against NaN / zero-norm.
- Tangent normalization and SVD align guarded against NaN / zero.
- Removed internal pysisyphus reference from `--input` help text.
- Docs: quickstart paths, troubleshooting, exit codes, `track_mode_by_overlap`, `dft/` in output tree, `scan`/`path-opt`/`irc` in See Also, `path_opt → path_search`, `--model → calc.model`, `--print-parsed` scope, JA toctree captions, `yaml-reference --config` flag.

## [0.3.4] — 2026-04-05

### Added
- Global pre-alignment stage and expanded smoke tests covering it.
- Auto-ECP selection for def2 basis sets; removed `--engine auto`.

### Fixed
- Exception-safe `AMINO_ACIDS` restore via `try/finally` in `extract()`.
- `bond-summary` PDB loading: `geom_from_pdb_str → geom_from_pdb`.
- Improved xTB not-found error message with install instructions.
- Documented GPU4PySCF Blackwell OOM limitation in `dft.md` (EN/JA).

## [0.3.3] — 2026-04-05

### Added
- JSON Output Reference page (EN/JA) covering `result.json` across every subcommand.
- `--out-json` to every subcommand; `summary.yaml` migrated to `summary.json`; `result.json` enriched with status, backend, and config fingerprint.
- `--modified-residue` option to `extract` and `all` commands plus troubleshooting entry.
- Python API Reference page; bond-summary / uma-pysis integrated with API references; standalone API page retired.
- Tsutsumi et al. 2022 citation for the bezA example system.

### Changed
- `--refine-path` promoted to the default `--help` display (previously `--help-advanced`); temporarily defaulted to `False` with `path-opt` as the primary MEP path (reverted in 0.3.5).
- Renamed `pocket → model` throughout filenames, directories, identifiers, CLI help, and user messages.
- Rewrote quickstart-scan as an `all --scan-lists` workflow guide.
- `reference/yaml.md` renamed to `reference/api-reference.md`; API page rewritten.
- `-s/--scan-lists` accepts multiple values (`-s a -s b`) instead of requiring re-specification.
- README examples migrated to bezA; examples/ cross-referenced from docs.
- EN and JA documentation pages resynced (structure and content).

### Fixed
- `_resolve_device`: handle `'auto' → cuda/cpu`.
- bezA description corrected (bornyl diphosphate synthase → methyltransferase; now sourced from Tsutsumi et al. 2022).
- Resume guard, first-input handling, zero-mass, `dir()` branches, scan guard, kink HEI, `_to_json` numpy/torch support, `BOND_KW` device, bare assert, dead `dft` check, missing `tabulate` dependency.
- MACE install docs: clarify that `fairchem-core` must be uninstalled first (e3nn conflict).
- Shell quoting in `pip install` extras examples (use double quotes).
- Documentation polish (EN + JA): JA UMA hard-code, benchmark GPU spec, pipeline diagram, path-search See Also, Stage 4 guard, `UMA → MLIP` label.

## [0.3.2] — 2026-03-24

### Added
- CITATION.cff with author + software co-author metadata for v0.3.2.
- `--verbose` flag to smoke test commands.
- Leading blank line before config blocks for readability; defaults filtering on bond/DMF blocks.
- Auto-shortening of absolute paths in all CLI output via `click.echo` patch.

### Changed
- RFO / RS-I-RFO trust-radius defaults lowered for MLIP stability; docs updated accordingly.
- Removed `--verbose` flag from command middle positions (now group-level / end-of-command).
- Shortened imaginary-mode filenames: `final_imag_mode → imag`.
- IRC initial-displacement clamp increased from 0.5 to 3.0 au.

### Fixed
- `click.echo` double blank lines: `_patched_echo` suppresses consecutive blanks; `pretty_block` spacing adjusted; section banners de-paded.
- IRC bisection: eliminated in-place mutation of initial displacement.
- JA doc defaults synced; `--scan-preopt` / `--scan-endopt` defaults; `--preopt` default (False → True).
- Theme toggle: document-level capture; 3-state → 2-state (light ↔ dark).
- `_patch_click_echo` definition order (NameError on import).
- `--verbose` flag semantics reconciled; always show full config dump.
- `reference/commands` and `reference/yaml` regenerated from code to match shipped CLI.

## [0.3.1] — 2026-03-18

### Added
- `logging.getLogger(__name__)` in all CLI modules for structured logging.
- `CONTRIBUTING.md` and `CHANGELOG.md`.
- `show_default=True` on all `--backend`, `--solvent`, and `--solvent-model` CLI options.
- Coverage reporting (`pytest-cov`) in CI.
- Additional unit tests for `summary_log`, `align`, `cli_utils`, and `defaults`.
- Bidirectional scan (4-tuple) documentation for `scan` command.
- `-b/--backend` and `-o/--out-dir` options added to all subcommands.
- `--dry-run` moved to `--help-advanced` across all subcommands.
- `--solvent` promoted to primary `--help` display.
- `-s/--scan-lists` unified with `--spec`: auto-detects inline literals vs YAML/JSON file paths.

### Fixed
- `_mep_skipped_by_resume` variable used before definition in `all.py`.
- Temporary directory leaks in `scan2d` and `scan3d` (added `finally: shutil.rmtree`).
- All `WARNING` messages in `all.py` now write to stderr (`err=True`).
- `ja/all.md` documented wrong default for `--opt-mode` (`hess` → `grad`).
- Duplicated `_is_param_explicit` helper replaced with `cli_param_overridden` from `utils`.
- Version banner no longer printed during `--help` tab completion (`ctx.resilient_parsing` guard).
- Documentation: `--preopt` / `--endopt` default values corrected from `True` to `False` in scan/scan2d/scan3d docs (EN/JA).

### Changed
- Centralized Click parameter-source checking via `cli_param_overridden(ctx, name)`.
- Bool options: `--flag/--no-flag` style promoted in docs and help (legacy `--flag True/False` still supported).
