#!/usr/bin/env python3
"""Validate agent-skill structure and high-risk pdb2reaction semantics.

This complements the live-CLI flag checker. It validates skill frontmatter,
reference coverage, and a small set of high-risk behavioral statements.
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

# Command -> page in skills/pdb2reaction-cli/. The all-* pages also belong to `all`.
CLI_COMMAND_PAGES = {
    "all": "all.md",
    "extract": "extract.md",
    "opt": "opt.md",
    "path-opt": "path.md",
    "path-search": "path.md",
    "scan": "scan.md",
    "scan2d": "scan.md",
    "scan3d": "scan.md",
    "tsopt": "tsopt.md",
    "irc": "irc.md",
    "freq": "freq.md",
    "dft": "dft.md",
    "sp": "utilities.md",
    "fix-altloc": "utilities.md",
    "add-elem-info": "utilities.md",
    "bond-summary": "utilities.md",
    "trj2fig": "utilities.md",
    "energy-diagram": "utilities.md",
}
ALL_MODE_PAGES = ("all-endpoint-mep.md", "all-scan-list.md", "all-ts-only.md")


def _issue(errors: list[str], path: Path, message: str, line: int | None = None) -> None:
    rel = path.relative_to(REPO_ROOT)
    where = f"{rel}:{line}" if line is not None else str(rel)
    entry = f"{where}: {message}"
    if entry not in errors:
        errors.append(entry)


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


def _cli_page_commands() -> dict[str, set[str]]:
    """Map each cli page name to the commands it covers."""
    owners: dict[str, set[str]] = {}
    for command, page in CLI_COMMAND_PAGES.items():
        owners.setdefault(page, set()).add(command)
    for page in ALL_MODE_PAGES:
        owners.setdefault(page, set()).add("all")
    return owners


def _validate_cli_page_coverage(errors: list[str]) -> None:
    commands = _live_subcommands()
    cli_dir = SKILLS_DIR / "pdb2reaction-cli"
    pages = {path.name for path in cli_dir.glob("*.md") if path.name != "SKILL.md"}
    unmapped = sorted(commands - set(CLI_COMMAND_PAGES))
    if unmapped:
        _issue(errors, cli_dir / "SKILL.md", f"commands missing from CLI_COMMAND_PAGES: {unmapped}")
    owners = _cli_page_commands()
    absent = sorted(set(owners) - pages)
    if absent:
        _issue(errors, cli_dir / "SKILL.md", f"cli pages listed in the checker but missing: {absent}")
    unknown = sorted(pages - set(owners))
    if unknown:
        _issue(errors, cli_dir / "SKILL.md", f"cli pages not assigned to a command: {unknown}")


def _validate_cli_option_table_ownership(errors: list[str]) -> None:
    """Ensure a per-command skill table does not borrow another command's flag.

    The broad drift check catches flags that exist nowhere. This check catches
    the subtler case where a valid flag is documented on the wrong command page.
    Only a table's first cell is inspected so cross-command prose remains legal.
    """
    cli_dir = SKILLS_DIR / "pdb2reaction-cli"
    flags_by_command = _live_subcommand_flags()
    for page, commands in sorted(_cli_page_commands().items()):
        path = cli_dir / page
        if not path.exists():
            _issue(errors, path, "cli page listed in the checker is missing")
            continue
        valid_flags = set().union(*(flags_by_command.get(c, set()) for c in commands))
        command_name = ", ".join(sorted(commands))
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


def _normalize(text: str) -> str:
    return " ".join(text.split())


def _require(path: Path, fragments: tuple[str, ...], errors: list[str]) -> None:
    if not path.exists():
        _issue(errors, path, "missing file with required high-risk guidance")
        return
    text = _normalize(path.read_text(encoding="utf-8"))
    for fragment in fragments:
        if _normalize(fragment) not in text:
            _issue(errors, path, f"required high-risk guidance missing: {fragment!r}")


def _validate_high_risk_semantics(errors: list[str]) -> None:
    cli = SKILLS_DIR / "pdb2reaction-cli"
    install = SKILLS_DIR / "pdb2reaction-install"
    structure = SKILLS_DIR / "pdb2reaction-model-setup"
    overview = SKILLS_DIR / "pdb2reaction-overview"
    all_page = cli / "all.md"
    scan_page = cli / "all-scan-list.md"
    tsopt_page = cli / "tsopt.md"
    opt_page = cli / "opt.md"
    irc_page = cli / "irc.md"
    freq_page = cli / "freq.md"
    hpc_page = SKILLS_DIR / "pdb2reaction-hpc" / "SKILL.md"
    backends_page = install / "backends.md"
    structure_page = structure / "SKILL.md"
    formats_page = structure / "formats.md"
    model_setup_page = SKILLS_DIR / "pdb2reaction-model-setup" / "SKILL.md"
    extract_page = cli / "extract.md"
    output_page = overview / "SKILL.md"
    summary_page = overview / "outputs.md"
    ts_strategy = overview / "ts-strategy.md"

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
    _require(
        tsopt_page,
        (
            "`--ref-mode`",
            "Ordinary standalone `tsopt` runs should omit it",
            "`tsopt` always forces `reject_uphill=False`",
        ),
        errors,
    )
    _require(
        opt_page,
        (
            "final convergence check on the retained geometry",
            "convergence requires ALL of `max(|force|) <= 3e-4`",
            "deliberately tightened variant of the published",
        ),
        errors,
    )
    _require(irc_page, ("reduce `--step-size` first", "`--never-stop`"), errors)
    _require(
        irc_page,
        (
            "`completed` is not an IRC convergence verdict",
            "`*_integration_converged`",
            "IRC has no independent scientific success verdict",
            "`never_stop_energy_bypasses`",
            "inserts one underscore",
        ),
        errors,
    )
    _require(
        freq_page,
        (
            "`freq` retains every signed physical mode",
            "ν < −5.00 cm⁻¹",
            "Raw negative counts are diagnostic",
            "freq.zero_cutoff_cm",
            "E + G_corr = G",
        ),
        errors,
    )
    for page in (irc_page, output_page, summary_page):
        if not page.exists():
            _issue(errors, page, "missing file checked for removed IRC verdicts")
            continue
        text = page.read_text(encoding="utf-8")
        for obsolete in ('d["forward_status"]', 'd["backward_status"]',
                         '`*_status == "stopped"`'):
            if obsolete in text:
                _issue(errors, page, f"removed IRC verdict used: {obsolete}")
    _require(
        summary_page,
        (
            "`references`",
            "{method, citation, doi}",
            "immediately before elapsed time",
        ),
        errors,
    )
    _require(
        hpc_page,
        ("BackendError", "FiniteDifference", "resource syntax is not interchangeable"),
        errors,
    )
    _require(backends_page, ("BackendError", "rather than changing the explicitly requested method"), errors)
    _require(
        backends_page,
        (
            "prebuilt PyTorch wheel contains its CUDA runtime dependencies",
            "Do not use `PYTORCH_NO_CUDA_PRELOAD`",
            "torch==2.13.0",
            "`cu126`, `cu130`, `cu132`, and `cpu`",
        ),
        errors,
    )
    _require(
        SKILLS_DIR / "pdb2reaction-install" / "SKILL.md",
        ("torch==2.13.0",),
        errors,
    )
    _require(
        backends_page,
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
        formats_page,
        (
            "619,938 residues",
            "`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`",
            "the `.cif` keeps the original chain IDs and residue numbers",
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
        structure_page,
        (
            "Explicit `-q` sets the total",
            "a mismatch produces a warning",
            "only chemically correct",
        ),
        errors,
    )
    _require(
        model_setup_page,
        ("Within one R/IM/P reaction path", "WT/mutant or other cross-variant models"),
        errors,
    )
    _require(
        summary_page,
        (
            "top-level `mlip_backend` / `mlip_model`",
            "`mlip_precision`",
            "`mlip`, `gibbs_mlip`, and `gibbs_dft_mlip` are the only emitted identifiers",
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
            "or guarantee identical output across PyTorch/backend versions",
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
