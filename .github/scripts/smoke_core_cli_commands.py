#!/usr/bin/env python3
"""Smoke test core CLI subcommands with lightweight fixtures."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = REPO_ROOT / "tests" / "smoke"
sys.path.insert(0, str(REPO_ROOT))


def _fixture_paths() -> dict[str, Path]:
    paths = {
        "r_pdb": FIXTURE_ROOT / "r.pdb",
        "p_pdb": FIXTURE_ROOT / "p.pdb",
        "r_complex_pdb": FIXTURE_ROOT / "r_complex.pdb",
        "r_complex_cif": FIXTURE_ROOT / "r_complex.cif",
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

    # Some grouped multi-value options intentionally inspect raw argv to
    # preserve stage boundaries. CliRunner does not update sys.argv itself, so
    # mirror a real process invocation while the command runs.
    original_argv = sys.argv[:]
    sys.argv = ["pdb2reaction", *args]
    try:
        result = runner.invoke(root_cli, args, catch_exceptions=False)
    finally:
        sys.argv = original_argv
    if result.exit_code != 0:
        raise RuntimeError(
            f"[core-cli-smoke] failed: {label}\n"
            f"command: pdb2reaction {' '.join(args)}\n"
            f"exit_code: {result.exit_code}\n"
            f"output:\n{result.output}"
        )


def _run_help_smoke(runner: CliRunner) -> int:
    subcommands = [
        "add-elem-info", "all", "bond-summary", "dft", "energy-diagram",
        "extract", "fix-altloc", "freq", "irc", "opt", "path-opt",
        "path-search", "scan", "scan2d", "scan3d", "sp", "trj2fig", "tsopt",
    ]
    advanced = {
        "all", "dft", "extract", "freq", "irc", "opt", "path-opt",
        "path-search", "scan", "scan2d", "scan3d", "sp", "tsopt",
    }
    for subcmd in subcommands:
        _invoke_or_raise(runner, [subcmd, "--help"], f"{subcmd} --help")
        if subcmd in advanced:
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
            "-1",
            "--dry-run",
        ],
        [
            "path-search",
            "-i",
            str(fixtures["r_pdb"]),
            "-i",
            str(fixtures["p_pdb"]),
            "-q",
            "-1",
            "--dry-run",
        ],
        ["tsopt", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["freq", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["irc", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["dft", "-i", str(fixtures["ts_gjf"]), "--dry-run"],
        ["opt", "-i", str(fixtures["r_pdb"]), "-q", "-1", "--dry-run"],
        ["sp", "-i", str(fixtures["r_pdb"]), "-q", "-1", "--dry-run"],
        [
            "scan", "-i", str(fixtures["r_pdb"]), "-q", "-1",
            "--scan-lists", "[(1,5,1.8)]", "--dry-run",
        ],
        [
            "scan2d", "-i", str(fixtures["r_pdb"]), "-q", "-1",
            "--scan-lists", "[(1,5,1.8,2.2),(1,6,1.7,2.1)]", "--dry-run",
        ],
        [
            "scan3d", "-i", str(fixtures["r_pdb"]), "-q", "-1",
            "--scan-lists", "[(1,5,1.8,2.2),(1,6,1.7,2.1),(2,3,1.0,1.2)]",
            "--dry-run",
        ],
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
                "-v",
                "0",
            ],
            "extract execution smoke",
        )
        if not output_pdb.exists():
            raise RuntimeError(
                "[core-cli-smoke] extract succeeded but output model was not created: "
                f"{output_pdb}"
            )

        cif_output_pdb = Path(tmpdir) / "model_from_cif.pdb"
        _invoke_or_raise(
            runner,
            [
                "extract",
                "-i",
                str(fixtures["r_complex_cif"]),
                "-c",
                "LONG_CHAIN:PRE:10001",
                "-o",
                str(cif_output_pdb),
                "-v",
                "0",
            ],
            "mmCIF extract execution smoke",
        )
        cif_output = cif_output_pdb.with_suffix(".cif")
        if not cif_output_pdb.exists() or not cif_output.exists():
            raise RuntimeError(
                "[core-cli-smoke] mmCIF extract did not create both internal PDB "
                f"and public CIF outputs: {cif_output_pdb}, {cif_output}"
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
    print("[core-cli-smoke] PDB and mmCIF extract execution smokes passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
