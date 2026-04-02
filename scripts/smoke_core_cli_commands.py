#!/usr/bin/env python3
"""Smoke test core CLI subcommands with lightweight fixtures."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_ROOT = REPO_ROOT / "tests" / "smoke"


def _fixture_paths() -> dict[str, Path]:
    paths = {
        "r_pdb": FIXTURE_ROOT / "r.pdb",
        "p_pdb": FIXTURE_ROOT / "p.pdb",
        "r_complex_pdb": FIXTURE_ROOT / "r_complex.pdb",
        "ts_gjf": FIXTURE_ROOT / "ts.gjf",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise RuntimeError(
            "[core-cli-smoke] missing required fixtures:\n  - " + "\n  - ".join(missing)
        )
    return {k: v.resolve() for k, v in paths.items()}


def _invoke_or_raise(runner: CliRunner, args: list[str], label: str) -> None:
    from pdb2reaction.cli import cli as root_cli

    result = runner.invoke(root_cli, args, catch_exceptions=False)
    if result.exit_code != 0:
        raise RuntimeError(
            f"[core-cli-smoke] failed: {label}\n"
            f"command: pdb2reaction {' '.join(args)}\n"
            f"exit_code: {result.exit_code}\n"
            f"output:\n{result.output}"
        )


def _run_help_smoke(runner: CliRunner) -> int:
    subcommands = ["extract", "path-opt", "path-search", "tsopt", "freq", "irc", "dft"]
    for subcmd in subcommands:
        _invoke_or_raise(runner, [subcmd, "--help"], f"{subcmd} --help")
        _invoke_or_raise(runner, [subcmd, "--help-advanced"], f"{subcmd} --help-advanced")
    return len(subcommands)


def _run_dry_run_smoke(runner: CliRunner, fixtures: dict[str, Path]) -> int:
    commands: list[list[str]] = [
        [
            "path-opt",
            "-i",
            str(fixtures["r_pdb"]),
            str(fixtures["p_pdb"]),
            "-q",
            "0",
            "--dry-run",
        ],
        [
            "path-search",
            "-i",
            str(fixtures["r_pdb"]),
            "-i",
            str(fixtures["p_pdb"]),
            "-q",
            "0",
            "--dry-run",
        ],
        ["tsopt", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["freq", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["irc", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["dft", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
    ]
    for cmd in commands:
        _invoke_or_raise(runner, cmd, " ".join(cmd[:1]) + " --dry-run")
    return len(commands)


def _run_extract_smoke(runner: CliRunner, fixtures: dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory(prefix="pdb2reaction_extract_smoke_") as tmpdir:
        output_pdb = Path(tmpdir) / "model.pdb"
        _invoke_or_raise(
            runner,
            [
                "extract",
                "-i",
                str(fixtures["r_complex_pdb"]),
                "-c",
                "PRE",
                "-o",
                str(output_pdb),
                "--no-verbose",
            ],
            "extract execution smoke",
        )
        if not output_pdb.exists():
            raise RuntimeError(
                "[core-cli-smoke] extract succeeded but output model was not created: "
                f"{output_pdb}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    fixtures = _fixture_paths()
    runner = CliRunner()

    n_help = _run_help_smoke(runner)
    n_dry = _run_dry_run_smoke(runner, fixtures)
    _run_extract_smoke(runner, fixtures)

    print(f"[core-cli-smoke] help validated for {n_help} subcommands.")
    print(f"[core-cli-smoke] dry-run validated for {n_dry} subcommands.")
    print("[core-cli-smoke] extract execution smoke passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
