#!/usr/bin/env python3
"""Lightweight smoke tests for commands embedded in docs markdown files."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL_NAME = "pdb2reaction"
CLI_MODULE = "pdb2reaction"
DOCS_SMOKE_COMMAND_TIMEOUT_SEC = float(os.environ.get("DOCS_SMOKE_COMMAND_TIMEOUT_SEC", "120"))

sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pdb2reaction.cli import cli as root_cli  # noqa: E402

import docs_command_contract as contract  # noqa: E402
from docs_command_contract import (  # noqa: E402
    AuthoredCommand,
    subcommand_from_tokens as _subcommand_from_tokens,
)


_ALL_ONLY_PATH_EXTS = {".pdb", ".xyz", ".gjf", ".yaml", ".yml", ".json"}


def _extract_docs_authored_commands() -> list[AuthoredCommand]:
    return contract.extract_markdown_commands(contract.docs_markdown_paths())


def _extract_docs_commands() -> list[str]:
    return [cmd.text for cmd in _extract_docs_authored_commands()]


def _prepare_fixture_files(tmp: Path) -> dict[str, Path]:
    pdb_text = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  LIG A   1       1.400   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    r_pdb = tmp / "R.pdb"
    p_pdb = tmp / "P.pdb"
    xyz = tmp / "input.xyz"
    gjf = tmp / "input.gjf"
    cfg = tmp / "config.yaml"
    calc_file = tmp / "my_calc.py"
    out_dir = tmp / "result_all"

    r_pdb.write_text(pdb_text, encoding="utf-8")
    p_pdb.write_text(pdb_text, encoding="utf-8")
    xyz.write_text("1\n\nC 0.0 0.0 0.0\n", encoding="utf-8")
    gjf.write_text("%chk=test\n#p hf/3-21g\n\nTitle\n\n0 1\nC 0.0 0.0 0.0\n\n", encoding="utf-8")
    cfg.write_text("extract:\n  radius: 2.6\n", encoding="utf-8")
    calc_file.write_text(
        "from ase.calculators.emt import EMT\n\n"
        "def get_calculator(**kwargs):\n"
        "    return EMT()\n",
        encoding="utf-8",
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    return {
        "r_pdb": r_pdb,
        "p_pdb": p_pdb,
        "xyz": xyz,
        "gjf": gjf,
        "config": cfg,
        "calc_file": calc_file,
        "out_dir": out_dir,
    }


def _sanitize_all_args(args: list[str], fixture: dict[str, Path]) -> list[str]:
    out: list[str] = []
    staged_scan = any(tok in {"-s", "--scan-lists"} for tok in args)
    fixture_inputs = [str(fixture["r_pdb"])]
    if not staged_scan:
        fixture_inputs.append(str(fixture["p_pdb"]))
    saw_input = False
    saw_dry_run = False
    saw_center = False
    i = 0
    while i < len(args):
        tok = args[i]
        if tok in {"--version", "-h", "--help"}:
            i += 1
            continue
        if tok in {"-c", "--center"}:
            # Doc examples center on real cofactor residues (SAM, GPP, ...);
            # the synthetic fixture only has a LIG residue. Normalize to LIG so
            # the dry-run extract pre-check exercises command structure, not the
            # specific residues, which the smoke fixture cannot carry.
            saw_center = True
            out.extend([tok, "LIG"])
            i += 2
            continue
        if tok in {"-l", "--ligand-charge", "-q", "--charge"}:
            # Drop system-specific charge args: they reference the doc example's
            # real residues/total charge, which would not match the LIG fixture.
            # Omitting them lets the extractor derive a consistent charge.
            i += 2
            continue
        if tok in {"-s", "--scan-lists"}:
            out.append(tok)
            i += 1
            n_values = 0
            while i < len(args) and not args[i].startswith("-"):
                out.append(f"[(1,2,{1.5 + 0.1 * n_values:.1f})]")
                n_values += 1
                i += 1
            continue
        if tok in {"-i", "--input"}:
            saw_input = True
            out.extend([tok, *fixture_inputs])
            i += 1
            while i < len(args) and not args[i].startswith("-"):
                i += 1
            continue
        if tok == "--config":
            out.extend([tok, str(fixture["config"])])
            i += 2
            continue
        if tok == "--calc-file":
            out.extend([tok, str(fixture["calc_file"])])
            i += 2
            continue
        if tok == "--out-dir":
            out.extend([tok, str(fixture["out_dir"])])
            i += 2
            continue
        if tok == "--dry-run":
            saw_dry_run = True
            out.append(tok)
            i += 1
            continue
        if tok == "--no-dry-run":
            saw_dry_run = True
            out.append("--dry-run")
            i += 1
            continue
        if (not tok.startswith("-")) and Path(tok).suffix.lower() in _ALL_ONLY_PATH_EXTS:
            ext = Path(tok).suffix.lower()
            if ext == ".pdb":
                out.append(str(fixture["r_pdb"]))
            elif ext == ".xyz":
                out.append(str(fixture["xyz"]))
            elif ext == ".gjf":
                out.append(str(fixture["gjf"]))
            else:
                out.append(str(fixture["config"]))
            i += 1
            continue

        out.append(tok)
        i += 1

    if not saw_input:
        out.extend(["-i", *fixture_inputs])
    if not saw_center:
        out.extend(["-c", "LIG"])
    if not saw_dry_run:
        out.append("--dry-run")
    if "--out-dir" not in out:
        out.extend(["--out-dir", str(fixture["out_dir"])])
    return out


def _run_help_smoke(commands: list[str]) -> None:
    runner = CliRunner()
    subcommands = sorted({_subcommand_from_tokens(shlex.split(cmd)) for cmd in commands})
    for subcmd in subcommands:
        result = runner.invoke(root_cli, [subcmd, "--help"], catch_exceptions=False)
        if result.exit_code != 0:
            raise RuntimeError(
                f"[help-smoke] failed for '{TOOL_NAME} {subcmd} --help':\n{result.output}"
            )
    print(f"[help-smoke] validated {len(subcommands)} subcommands from docs.")


def _validate_authored_commands(commands: list[AuthoredCommand]) -> None:
    """Validate every authored option token via the shared contract module.

    Retains examples containing ``<placeholder>`` / ``[...]`` notation for
    static validation and reports each unknown option with ``file:line``.
    """
    errors = contract.validate_option_names(commands)
    if errors:
        raise RuntimeError("[option-smoke] docs option validation failed:\n" + "\n".join(errors))
    print(f"[option-smoke] validated option names in {len(commands)} docs examples.")


def _validate_option_names(commands: list[str]) -> None:
    """Backwards-compatible string entrypoint used by the regression tests."""
    _validate_authored_commands([contract._make(REPO_ROOT, 0, cmd) for cmd in commands])


def _select_executable_all_commands(commands: list[str]) -> list[str]:
    """Select concrete ``all`` examples, retaining synopsis for static checks only."""

    all_cmds: set[str] = set()
    for cmd in commands:
        if not contract.is_executable(cmd):
            continue
        tokens = shlex.split(cmd)
        if not tokens or tokens[0] != TOOL_NAME:
            continue
        if len(tokens) >= 2 and tokens[1] == "all":
            all_cmds.add(cmd)
            continue
        if len(tokens) >= 2 and tokens[1].startswith("-"):
            if any(tok in {"-i", "--input"} for tok in tokens[1:]):
                all_cmds.add(cmd)
    return sorted(all_cmds)


def _run_all_dry_run_smoke(commands: list[str]) -> None:
    try:
        probe = subprocess.run(
            [sys.executable, "-m", CLI_MODULE, "all", "--help"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=DOCS_SMOKE_COMMAND_TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            "[dry-run-smoke] timeout while probing availability for "
            f"'{TOOL_NAME} all --help' ({DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s)."
        ) from exc
    probe_output = f"{probe.stdout}\n{probe.stderr}"
    if "Command 'all' is unavailable" in probe_output or "Missing dependency:" in probe_output:
        raise RuntimeError(
            "[dry-run-smoke] required 'all' command is unavailable in this environment."
        )

    all_cmds = _select_executable_all_commands(commands)
    if not all_cmds:
        raise RuntimeError("No 'all' command examples found in docs.")

    with tempfile.TemporaryDirectory(prefix=f"{TOOL_NAME}_docs_smoke_") as tmpdir:
        fixture = _prepare_fixture_files(Path(tmpdir))
        for raw in all_cmds:
            tokens = shlex.split(raw)
            if not tokens or tokens[0] != TOOL_NAME:
                continue
            args = tokens[1:]
            if not args or args[0].startswith("-"):
                args = ["all", *args]
            if args[0] != "all":
                continue
            dry_args = _sanitize_all_args(args, fixture)
            try:
                completed = subprocess.run(
                    [sys.executable, "-m", CLI_MODULE, *dry_args],
                    cwd=REPO_ROOT,
                    text=True,
                    capture_output=True,
                    timeout=DOCS_SMOKE_COMMAND_TIMEOUT_SEC,
                )
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(
                    f"[dry-run-smoke] timeout for docs command ({DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s):\n"
                    f"  {raw}\n"
                    f"sanitized args:\n"
                    f"  {TOOL_NAME} {' '.join(dry_args)}"
                ) from exc
            if completed.returncode != 0:
                raise RuntimeError(
                    f"[dry-run-smoke] failed for docs command:\n  {raw}\n"
                    f"sanitized args:\n  {TOOL_NAME} {' '.join(dry_args)}\n\n"
                    f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
                )
    print(
        f"[dry-run-smoke] validated {len(all_cmds)} docs examples "
        f"(timeout={DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s)."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    authored = _extract_docs_authored_commands()
    if not authored:
        raise RuntimeError("No commands were extracted from docs markdown code fences.")

    commands = [cmd.text for cmd in authored]
    _run_help_smoke(commands)
    _validate_authored_commands(authored)
    _run_all_dry_run_smoke(commands)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
