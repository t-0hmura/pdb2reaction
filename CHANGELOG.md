# Changelog

All notable changes to **pdb2reaction** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

## [0.3.8] — 2026-04-29

### Added
- `pdb2reaction dft --lowmem/--no-lowmem` (default `True`): closed-shell GPU runs now use `gpu4pyscf.dft.rks_lowmem.RKS`, which performs SCF with a memory-efficient direct-JK pipeline (no density fitting). Open-shell, CPU, and pre-`rks_lowmem` `gpu4pyscf` installs auto-fall back to standard `RKS`/`UKS`. Selectable via the new `--lowmem/--no-lowmem` CLI flag and `dft.lowmem` YAML key. Population analysis on lowmem builds a CPU surrogate from the converged MOs because `rks_lowmem.RKS.to_cpu()` is not implemented upstream.
- Agent skills: cluster-boundary frozen-atoms reference (`pdb2reaction-cli/freeze-atoms.md`); 1-line input→output cheatsheet on `pdb2reaction-cli/SKILL.md`; `summary.json` schema split out to `pdb2reaction-workflows-output/summary-json.md`.
- `scripts/check_skill_drift.py`: warning-only prose-vs-source consistency check for the agent skills, wired into Docs Quality.

### Changed
- Closed-shell GPU DFT defaults switch from DF + standard `RKS` to direct-JK `rks_lowmem.RKS`. Absolute energies shift by sub-mHa relative to v0.3.6/v0.3.7; pass `--no-lowmem` (or set `dft.lowmem: false`) to reproduce earlier numbers.
- HPC PBS / SLURM preambles in the agent skills now also load `gcc` and an optional `<OPENMPI_MODULE>` (the latter only when running multi-node Ray); inline notes describe when each module is actually required.
- Sphinx config: drop the PyTorch `intersphinx_mapping` entry that was causing build hangs and remove the now-unused `intersphinx_timeout`.
- Documentation tone (EN+JA): `docs/index.md` "Start here" and per-page "At a glance" hooks rephrased from first-person scenarios to declarative goal phrases.

### Fixed
- TS-only mode under `pdb2reaction all`: docs and agent skills now correctly state that a single `-i` requires `--tsopt` (in addition to no `--scan-lists`); without `--tsopt`, `BadParameter` is raised.
- `pdb2reaction-install-backends/env-cuda.md`: DFT does not auto-fall back to CPU when GPU4PySCF is missing; `--engine cpu` must be passed explicitly.
- `pdb2reaction-install-backends/xtb.md`: `--solvent` is accepted by `scan`, `scan2d`, and `scan3d`; only `dft` and `extract` reject it.
- `pdb2reaction-cli/freq.md` and `pdb2reaction-structure-io/pdb.md`: PDB B-factors are not a freeze flag for pdb2reaction; the freeze set is assembled from `--freeze-atoms`, `--freeze-links`, and YAML `geom.freeze_atoms`.
- Output-tree, JSON schema, and scan-list grammar references in the skills realigned with the source.

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
