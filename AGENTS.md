# Working on pdb2reaction

This file is for coding agents that change the code, docs, or skills in this repository.
To run pdb2reaction for a user, load the skills in `skills/` instead (see `skills/README.md`).
The full rules and their reasons are in `CONTRIBUTING.md`; § numbers below refer to it.

## Before you finish

1. `pytest tests/ -q`. Fix a failing test at its cause; never delete or skip it.
2. `python .github/scripts/check_engineering_markers.py` and `python .github/scripts/check_import_graph.py`.
3. `python .github/scripts/run_docs_quality.py`. It stops at the first failing step; to see every failure, run the scripts it lists one by one.
4. `python .github/scripts/check_help_registry.py` after CLI changes.
5. Smoke: copy `tests/smoke/` outside the repository and run `bash run.sh` through your scheduler.

After changing CLI options, regenerate `docs/reference/commands/` with `python .github/scripts/generate_reference.py`.
Write boolean options as `--flag/--no-flag`.

## Do not touch

Read §4 before changing any of these:

- Chemistry rules marked `# CHEMISTRY-RULE:N`, and `# DOMAIN_PURE` modules (§4.1).
- The `del` chains that release GPU memory (§4.2).
- Divergent files in the bundled forks (§4.3).
- Packaging settings (§4.4).
- Absolute paths in `_LAZY_SUBCOMMANDS` (§4.5).
- Chemistry default choices (§4.6).
- Output that downstream parsers read, including `summary.json` (§4.7, §1.4).

## Editing skills

Each skill is `skills/<name>/SKILL.md`.
After editing skills, run `python .github/scripts/run_docs_quality.py` (it runs the skill checkers) and `pytest tests/test_skill_command_smoke.py tests/test_cli_completion.py -q`.
The checkers require:

- Frontmatter with exactly `name` and `description`; `name` is hyphen-case and equals the folder name.
- A `description` of at most 1024 characters, without `<` or `>`, that states TRIGGER and SKIP.
- A `SKILL.md` of at most 500 lines.
- Real subcommands and flags in shell code blocks (no language, `bash`, `console`, `sh`, or `shell`) and in inline `pdb2reaction …` (`check_skill_commands.py`). Write placeholders as `<...>`.
- Every subcommand is assigned to a page in `skills/pdb2reaction-cli/` through `CLI_COMMAND_PAGES` in `check_skill_quality.py`, and every page there is assigned; add a new page to that map. Flags in the first column of a page's tables belong to its assigned commands.
- One `-s/--scan-lists` per `all` or `scan` command.
- Real flags in backticks in prose in the folders listed in `PDB2REACTION_CLI_DIRS` (`check_skill_drift.py`). Add a new skill folder that documents the CLI there.
- The sentences that `check_skill_quality.py` pins, and none of the phrases it bans. When you move or reword a pinned sentence, update the checker in the same change.

Do not add flag tables to skills; point to `--help-advanced` and `docs/reference/commands/`. Describe current behaviour only.
When you add, remove, or rename a skill, update `skills/README.md` and the `[Unreleased]` section of `CHANGELOG.md` in the same change.

## More

[CONTRIBUTING.md](CONTRIBUTING.md) · [Architecture](docs/architecture.md) · [Skills index](skills/README.md)
