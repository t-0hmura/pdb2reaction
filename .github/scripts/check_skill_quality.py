#!/usr/bin/env python3
"""Validate agent-skill structure and high-risk pdb2reaction semantics.

This complements the live-CLI flag checker. It validates skill frontmatter,
reference coverage, and a small set of behavioral statements
whose ambiguity has caused real operator/agent mistakes.
"""

from __future__ import annotations

import re
import shlex
import sys
from pathlib import Path

import click
import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_DIR = REPO_ROOT / "skills"
FRONTMATTER_RE = re.compile(r"\A---\r?\n(.*?)\r?\n---(?:\r?\n|\Z)", re.DOTALL)
NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
LONG_FLAG_RE = re.compile(r"--[a-z][a-z0-9-]*")


def _issue(errors: list[str], path: Path, message: str, line: int | None = None) -> None:
    rel = path.relative_to(REPO_ROOT)
    where = f"{rel}:{line}" if line is not None else str(rel)
    errors.append(f"{where}: {message}")


def _line_of(text: str, token: str) -> int:
    pos = text.find(token)
    return 1 if pos < 0 else text.count("\n", 0, pos) + 1


def _validate_root_skill(path: Path, errors: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    if len(text.splitlines()) > 500:
        _issue(errors, path, "SKILL.md exceeds the 500-line progressive-disclosure limit")

    match = FRONTMATTER_RE.match(text)
    if match is None:
        _issue(errors, path, "missing or malformed YAML frontmatter")
        return
    try:
        frontmatter = yaml.safe_load(match.group(1))
    except yaml.YAMLError as exc:
        _issue(errors, path, f"invalid YAML frontmatter: {exc}")
        return
    if not isinstance(frontmatter, dict):
        _issue(errors, path, "frontmatter must be a mapping")
        return

    expected_keys = {"name", "description"}
    if set(frontmatter) != expected_keys:
        _issue(
            errors,
            path,
            f"frontmatter keys must be exactly {sorted(expected_keys)}, got {sorted(frontmatter)}",
        )

    name = frontmatter.get("name")
    description = frontmatter.get("description")
    if not isinstance(name, str) or not NAME_RE.fullmatch(name):
        _issue(errors, path, f"invalid hyphen-case skill name: {name!r}")
    elif name != path.parent.name:
        _issue(errors, path, f"skill name {name!r} does not match directory {path.parent.name!r}")

    if not isinstance(description, str) or not description.strip():
        _issue(errors, path, "description must be a non-empty string")
        return
    if len(description) > 1024:
        _issue(errors, path, "description exceeds 1024 characters")
    if "<" in description or ">" in description:
        _issue(errors, path, "description must not contain angle brackets")
    if "TRIGGER" not in description and "use only when" not in description.lower():
        _issue(errors, path, "description must state when the skill should trigger")
    if "SKIP" not in description and "use only when" not in description.lower():
        _issue(errors, path, "description must state a skip boundary")


def _live_subcommands() -> set[str]:
    sys.path.insert(0, str(REPO_ROOT))
    from pdb2reaction.cli import cli as root_cli

    ctx = click.Context(root_cli)
    return set(root_cli.list_commands(ctx))


def _live_subcommand_flags() -> dict[str, set[str]]:
    """Return the registered long options for each live Click subcommand."""
    sys.path.insert(0, str(REPO_ROOT))
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.cli.app import _COMMAND_BOOL_VALUE_OPTIONS

    flags_by_command: dict[str, set[str]] = {}
    ctx = click.Context(root_cli)
    for name in root_cli.list_commands(ctx):
        command = root_cli.get_command(ctx, name)
        if command is None:
            continue
        flags = {"--help", "--help-advanced"}
        for parameter in command.params:
            flags.update(
                option
                for option in (
                    *(getattr(parameter, "opts", ()) or ()),
                    *(getattr(parameter, "secondary_opts", ()) or ()),
                )
                if option.startswith("--")
            )
        # ``all`` keeps legacy Click BOOL parameters for compatibility; the
        # argv normalizer also exposes their canonical synthetic --no-* form.
        for option in _COMMAND_BOOL_VALUE_OPTIONS.get(name, ()):
            flags.add(option)
            flags.add(f"--no-{option[2:]}")
        flags_by_command[name] = flags
    return flags_by_command


def _validate_cli_page_coverage(errors: list[str]) -> None:
    commands = _live_subcommands()
    cli_dir = SKILLS_DIR / "pdb2reaction-cli"
    pages = {path.stem for path in cli_dir.glob("*.md") if path.name != "SKILL.md"}
    missing = sorted(commands - pages)
    if missing:
        _issue(errors, cli_dir / "SKILL.md", f"missing per-command pages: {missing}")


def _validate_cli_option_table_ownership(errors: list[str]) -> None:
    """Ensure a per-command skill table does not borrow another command's flag.

    The broad drift check catches flags that exist nowhere. This check catches
    the subtler case where a valid flag is documented on the wrong command page.
    Only a table's first cell is inspected so cross-command prose remains legal.
    """
    cli_dir = SKILLS_DIR / "pdb2reaction-cli"
    for command_name, valid_flags in _live_subcommand_flags().items():
        path = cli_dir / f"{command_name}.md"
        if not path.exists():
            continue
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.startswith("|"):
                continue
            first_cell = line.split("|", 2)[1]
            for flag in LONG_FLAG_RE.findall(first_cell):
                if flag not in valid_flags:
                    _issue(
                        errors,
                        path,
                        f"option-table flag {flag} is not registered on {command_name!r}",
                        lineno,
                    )


def _iter_logical_commands(path: Path):
    """Yield shell command text from fenced blocks, joining continuations."""
    in_fence = False
    pending: list[str] = []
    start = 0
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = raw.strip()
        if stripped.startswith("```"):
            if in_fence and pending:
                yield start, " ".join(pending)
                pending = []
            in_fence = not in_fence
            continue
        if not in_fence:
            continue
        if stripped.endswith("\\"):
            if not pending:
                start = lineno
            pending.append(stripped[:-1].strip())
            continue
        if pending:
            pending.append(stripped)
            yield start, " ".join(pending)
            pending = []
        elif stripped.startswith("pdb2reaction "):
            yield lineno, stripped


def _validate_scan_flag_occurrences(path: Path, errors: list[str]) -> None:
    for lineno, command in _iter_logical_commands(path):
        try:
            tokens = shlex.split(command)
        except ValueError:
            continue
        if len(tokens) < 2 or tokens[:2] not in (["pdb2reaction", "all"], ["pdb2reaction", "scan"]):
            continue
        count = sum(token in {"-s", "--scan-lists"} for token in tokens)
        if count > 1:
            _issue(
                errors,
                path,
                "repeat of -s/--scan-lists in one command; use one flag followed by all stage values",
                lineno,
            )


def _require(path: Path, fragments: tuple[str, ...], errors: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for fragment in fragments:
        if fragment not in text:
            _issue(errors, path, f"required high-risk guidance missing: {fragment!r}")


def _validate_high_risk_semantics(errors: list[str]) -> None:
    all_page = SKILLS_DIR / "pdb2reaction-cli" / "all.md"
    scan_page = SKILLS_DIR / "pdb2reaction-cli" / "all-scan-list.md"
    tsopt_page = SKILLS_DIR / "pdb2reaction-cli" / "tsopt.md"
    irc_page = SKILLS_DIR / "pdb2reaction-cli" / "irc.md"
    path_opt_page = SKILLS_DIR / "pdb2reaction-cli" / "path-opt.md"
    freq_page = SKILLS_DIR / "pdb2reaction-cli" / "freq.md"
    hpc_page = SKILLS_DIR / "pdb2reaction-hpc" / "SKILL.md"
    uma_page = SKILLS_DIR / "pdb2reaction-install-backends" / "uma.md"
    cuda_page = SKILLS_DIR / "pdb2reaction-install-backends" / "env-cuda.md"
    xtb_page = SKILLS_DIR / "pdb2reaction-install-backends" / "xtb.md"
    structure_page = SKILLS_DIR / "pdb2reaction-structure-io" / "SKILL.md"
    cif_page = SKILLS_DIR / "pdb2reaction-structure-io" / "cif.md"
    charge_page = SKILLS_DIR / "pdb2reaction-structure-io" / "charge-multiplicity.md"
    pdb_page = SKILLS_DIR / "pdb2reaction-structure-io" / "pdb.md"
    extract_page = SKILLS_DIR / "pdb2reaction-cli" / "extract.md"
    output_page = SKILLS_DIR / "pdb2reaction-workflows-output" / "SKILL.md"
    summary_page = SKILLS_DIR / "pdb2reaction-workflows-output" / "summary-json.md"
    ts_strategy = SKILLS_DIR / "pdb2reaction-ts-strategy" / "SKILL.md"

    _require(all_page, ("temporary directory", "one `-s` occurrence"), errors)
    _require(
        scan_page,
        (
            "Use exactly one `--scan-lists` flag",
            "Repeating the flag is rejected",
            "`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`",
        ),
        errors,
    )
    _require(tsopt_page, ("`--ref-mode`", "ordinary standalone `tsopt` runs should omit it"), errors)
    _require(irc_page, ("reduce `--step-size` first", "`--never-stop`"), errors)
    _require(
        irc_page,
        (
            "`completed` is not an IRC convergence verdict",
            "`*_converged`",
            "`never_stop_energy_bypasses`",
            "inserts one underscore",
        ),
        errors,
    )
    _require(
        freq_page,
        ("every strictly negative frequency", "`tsopt` deliberately applies"),
        errors,
    )
    _require(path_opt_page, ("| `--max-nodes` | int | 20 |",), errors)
    _require(
        hpc_page,
        ("BackendError", "FiniteDifference", "resource syntax is not interchangeable"),
        errors,
    )
    _require(uma_page, ("BackendError", "rather than changing the explicitly requested method"), errors)
    _require(
        cuda_page,
        (
            "prebuilt PyTorch wheel contains its CUDA runtime dependencies",
            "Do not use\n`PYTORCH_NO_CUDA_PRELOAD`",
        ),
        errors,
    )
    _require(
        xtb_page,
        (
            "E_xTB(solvent) - E_xTB(vacuum)",
            "`trj2fig` also accepts it, but only uses it when",
            "It is not accepted by `dft`, `extract`,",
        ),
        errors,
    )
    _require(
        structure_page,
        ("standard amino acids and recognized ions use internal tables",),
        errors,
    )
    _require(
        cif_page,
        (
            "619,938 residues",
            "`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`",
            "synthetic identifiers are implementation details",
            "With `--convert-files` enabled",
        ),
        errors,
    )
    _require(
        extract_page,
        ("`'B:SAM'`", "`'B:SAM:321'`", "same stem as `.cif`"),
        errors,
    )
    _require(
        output_page,
        ("mmCIF/oversized-PDB inputs also receive `.cif` companions",),
        errors,
    )
    _require(
        charge_page,
        (
            "`-q` is an assertion",
            "`-q` does not silently replace extraction's charge result",
            "only chemically correct",
        ),
        errors,
    )
    _require(
        pdb_page,
        ("Within one R/IM/P reaction path", "WT/mutant or other cross-variant models"),
        errors,
    )
    _require(
        summary_page,
        (
            "top-level `mlip_backend` / `mlip_model`",
            "`mlip_precision`",
            "`mlip`, `gibbs_mlip`, and `gibbs_dft_mlip` are the only emitted identifiers",
            "no compatibility aliases are written",
            "filenames use `MLIP`",
        ),
        errors,
    )
    _require(
        ts_strategy,
        (
            "--refine-path",
            "deliberately off by default",
            "`--ref-mode` is not a normal standalone",
            "or guarantee identical\n  output across PyTorch/backend versions",
        ),
        errors,
    )

    banned = {
        "auto-downgrades to finite differences": "an explicit UMA analytical+multi-worker request is an error",
        "only works for UMA; other backends fall back": "all four built-in backends support explicit Analytical Hessians",
        "each segment crosses exactly one TS": "path segmentation proposes candidates; TS/frequency/IRC validation is required",
        "torch_scatter": "current orb-models does not require torch_scatter; diagnose actual package metadata",
        "closed-shell systems; use `-m 1` only": "AIMNet2 has model-dependent open-shell support",
        "Both are accepted by most": "Torque and PBSPro resource syntax is site-specific and not interchangeable",
        "E_MLIP_or_DFT": "xTB correction wraps the base MLIP/custom calculator, not standalone dft",
        "derived total always": "charge derivation is mechanically consistent but still depends on chemically correct residue states",
        "Never re-extract compared states independently": "distinguish one-path identity requirements from cross-variant comparisons",
        "chain IDs are not part of the spec": "chain-qualified atom selectors are supported and required for ambiguous repeated residues",
    }
    for path in sorted(SKILLS_DIR.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        for token, replacement in banned.items():
            if token in text:
                _issue(errors, path, f"stale/ambiguous phrase {token!r}; {replacement}", _line_of(text, token))


def main() -> int:
    errors: list[str] = []
    root_skills = sorted(SKILLS_DIR.glob("*/SKILL.md"))
    markdown_files = sorted(SKILLS_DIR.rglob("*.md"))
    if not root_skills:
        print("[skill-quality] no skills found")
        return 1

    for path in root_skills:
        _validate_root_skill(path, errors)
    for path in markdown_files:
        _validate_scan_flag_occurrences(path, errors)
    _validate_cli_page_coverage(errors)
    _validate_cli_option_table_ownership(errors)
    _validate_high_risk_semantics(errors)

    if errors:
        print(f"[skill-quality] FAILED: {len(errors)} issue(s)")
        for error in errors:
            print(f"- {error}")
        return 1
    print(
        f"[skill-quality] OK: {len(root_skills)} root skills, "
        f"{len(markdown_files)} markdown files"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
