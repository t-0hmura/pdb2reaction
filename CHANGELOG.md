# Changelog

All notable changes to **pdb2reaction** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

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
