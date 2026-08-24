#!/usr/bin/env python3
"""Detect drift between ``skills/`` markdown and pdb2reaction
source for content that ``check_skill_commands.py`` does not cover —
prose tables, JSON snippets, and output trees rather than bash blocks.

Catches three classes of staleness:

1. **Unknown CLI flags** — any ``--flag-name`` token that exists nowhere
   in the Click subcommand graph. Either a typo or a removed flag.
2. **Unknown status literals** — result status is command-specific:
   optimizers use ``converged`` / ``not_converged``; completed stage
   runners use ``completed`` or ``ok``; aggregate workflows use
   ``success`` / ``partial`` / ``failed``; failures use ``error``.
3. **Renamed strings** — file names and JSON keys that were renamed in
   the source (``opt_trj.xyz`` → ``optimization_trj.xyz``,
   ``"n_cycles"`` → ``"n_opt_cycles"``, etc.). Maintained as a small
   allowlist of (old, new) pairs.

Exits non-zero on drift so incorrect agent instructions cannot merge.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import click

REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_DIR = REPO_ROOT / "skills"

sys.path.insert(0, str(REPO_ROOT))
from pdb2reaction.cli import cli as root_cli  # noqa: E402
from pdb2reaction.core.utils import RESULT_JSON_STATUS_VALUES  # noqa: E402


# Renamed file/key/string literals — (old, new, note).
# NOTE: pdb2reaction keeps ``final_energy_hartree`` as a canonical
# per-stage scan key (see pdb2reaction/workflows/scan.py and
# docs/json-output.md), so it is intentionally NOT listed as a rename.
RENAMED_STRINGS: list[tuple[str, str, str]] = [
    ("opt_trj.xyz", "optimization_trj.xyz", "trajectory file renamed"),
    ("tsopt_trj.xyz", "optimization_trj.xyz", "trajectory file renamed"),
    ("\"n_cycles\"", "\"n_opt_cycles\"", "result.json key renamed"),
    ("\"gradient_max\"", "\"final_max_force\"", "result.json key renamed"),
    ("\"structure_path\"", "files.final_geometry_xyz",
     "structure_path is no longer a top-level result.json key"),
    ("energy_diagram_UMA", "energy_diagram_MLIP", "backend-neutral diagram filename"),
    ("energy_diagram_G_UMA", "energy_diagram_G_MLIP", "backend-neutral diagram filename"),
    ("energy_diagram_G_DFT_plus_UMA", "energy_diagram_G_DFT_plus_MLIP", "backend-neutral diagram filename"),
    ("DFT//UMA", "DFT//MLIP", "backend-neutral human-readable method label"),
    ("uma_model", "mlip_model", "backend-neutral provenance field"),
    ("opt.log", "(none — pysisyphus loggers are silenced; use --dump)",
     "skill referenced a log file that is not produced"),
    ("tsopt.log", "(none — pysisyphus loggers are silenced; use --dump)",
     "skill referenced a log file that is not produced"),
    ("freq.log", "(none — pysisyphus loggers are silenced; use --dump)",
     "skill referenced a log file that is not produced"),
    ("irc.log", "(none — pysisyphus loggers are silenced; use --dump)",
     "skill referenced a log file that is not produced"),
]

# Status is command-specific; this is the union of public values documented in
# docs/json-output.md and emitted by the workflows.
CANONICAL_STATUS: set[str] = set(RESULT_JSON_STATUS_VALUES)

# MCP runner-level status vocabulary (pdb2reaction/mcp/_runner.py): a
# distinct enum from the CLI summary.json status above, valid only inside
# the pdb2reaction-mcp skill.
MCP_STATUS: set[str] = {"ok", "summary_missing", "summary_parse_error"}
MCP_SKILL_DIRS: set[str] = {"pdb2reaction-mcp"}

# Match a backtick-quoted flag in prose (`--foo`) — the typical
# "documented CLI option" callout in skill markdown tables. Bash
# blocks that begin with ``pdb2reaction`` are already covered by
# .github/scripts/check_skill_commands.py; this script targets prose where
# bare ``--xxx`` could be a non-pdb2reaction system tool.
FLAG_TOKEN_RE = re.compile(r"`(--[a-z][a-z0-9-]*)`")
STATUS_RE = re.compile(r'"status"\s*:\s*"([a-zA-Z_]+)"')

# Skill subdirs that document the pdb2reaction CLI (not external tools
# like pip / conda / nvidia-smi / sbatch). Unknown-flag check is
# restricted to these so legitimate ``--cpus-per-task`` etc. references
# in HPC / env-detect / install-backends docs are not flagged.
PDB2REACTION_CLI_DIRS = {
    "pdb2reaction-cli",
    "pdb2reaction-overview",
    "pdb2reaction-workflows-output",
    "pdb2reaction-structure-io",
}


def _collect_flag_union() -> set[str]:
    """Return the union of all ``--flag`` tokens registered on any
    Click subcommand. Used to catch typos / removed flags."""
    union: set[str] = set()
    ctx = click.Context(root_cli)
    for name in root_cli.list_commands(ctx):
        cmd = root_cli.get_command(ctx, name)
        if cmd is None:
            continue
        for p in cmd.params:
            for opt in getattr(p, "opts", []) or []:
                if opt.startswith("--"):
                    union.add(opt)
            for opt in getattr(p, "secondary_opts", []) or []:
                if opt.startswith("--"):
                    union.add(opt)
    # universal flags every subcommand inherits
    union.update({"--help", "--help-advanced"})
    return union


def _scan_file(path: Path, flag_union: set[str]) -> list[str]:
    rel = path.relative_to(REPO_ROOT)
    warnings: list[str] = []
    try:
        text = path.read_text()
    except UnicodeDecodeError:
        return warnings

    in_pdb2reaction_cli_dir = (
        len(rel.parts) >= 2 and rel.parts[1] in PDB2REACTION_CLI_DIRS
    )
    in_mcp_dir = len(rel.parts) >= 2 and rel.parts[1] in MCP_SKILL_DIRS

    for lineno, line in enumerate(text.splitlines(), start=1):
        # 1. Unknown flag tokens (only in pdb2reaction-cli skill dirs)
        if in_pdb2reaction_cli_dir:
            for m in FLAG_TOKEN_RE.finditer(line):
                flag = m.group(1)
                # accept --no-<x> if --<x> is registered (Click bool toggle)
                if flag.startswith("--no-"):
                    base = "--" + flag[5:]
                    if base in flag_union:
                        continue
                if flag in flag_union:
                    continue
                warnings.append(
                    f"{rel}:{lineno}: unknown flag {flag} "
                    f"(not registered on any subcommand)"
                )

        # 2. Status enum literals
        allowed_status = CANONICAL_STATUS | (MCP_STATUS if in_mcp_dir else set())
        for m in STATUS_RE.finditer(line):
            v = m.group(1)
            if v not in allowed_status:
                warnings.append(
                    f"{rel}:{lineno}: status {v!r} not in canonical "
                    f"{{{', '.join(sorted(allowed_status))}}}"
                )

        # 3. Renamed strings
        for old, new, note in RENAMED_STRINGS:
            if old in line:
                warnings.append(
                    f"{rel}:{lineno}: stale {old!r} ({note}); use {new!r}"
                )

    return warnings


def main() -> int:
    flag_union = _collect_flag_union()
    all_warnings: list[str] = []
    n_files = 0
    for path in sorted(SKILLS_DIR.rglob("*.md")):
        all_warnings.extend(_scan_file(path, flag_union))
        n_files += 1

    if all_warnings:
        print(
            f"[skill-drift] WARN: {len(all_warnings)} warning(s) "
            f"across {n_files} file(s):\n"
        )
        for w in all_warnings:
            print(w)
        print(
            "\nAddress the items above (or update RENAMED_STRINGS / "
            "CANONICAL_STATUS in .github/scripts/check_skill_drift.py if a "
            "documented schema change is intentional)."
        )
    else:
        print(f"[skill-drift] OK: no warnings across {n_files} file(s)")
    return 1 if all_warnings else 0


if __name__ == "__main__":
    raise SystemExit(main())
