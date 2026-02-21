#!/usr/bin/env python3
"""Smoke-test trajectory dump behavior for opt/tsopt commands."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOL_NAME = "pdb2reaction"


@dataclass(frozen=True)
class Case:
    name: str
    args: list[str]
    expect_present: tuple[str, ...]
    expect_absent: tuple[str, ...] = ()


def _write_tiny_xyz(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "3",
                "water guess",
                "O 0.000000 0.000000 0.000000",
                "H 0.957200 0.000000 0.000000",
                "H -0.239987 0.927297 0.000000",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _frame_count(xyz_path: Path) -> int:
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    i = 0
    n = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        nat = int(lines[i].strip())
        i += nat + 2
        n += 1
    return n


def _run_cli(args: list[str]) -> None:
    cmd = [TOOL_NAME, *args]
    proc = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        tail = (proc.stdout + "\n" + proc.stderr)[-3000:]
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{tail}")


def _check_runnable() -> tuple[bool, str]:
    exe = shutil.which(TOOL_NAME)
    if exe is None:
        return False, f"'{TOOL_NAME}' is not in PATH."

    probe = subprocess.run(
        [TOOL_NAME, "opt", "--help"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )
    probe_out = f"{probe.stdout}\n{probe.stderr}"
    if "Command 'opt' is unavailable" in probe_out or "Missing dependency:" in probe_out:
        return False, "opt command is unavailable in this environment."
    if probe.returncode != 0:
        return False, f"probe command failed with exit code {probe.returncode}."
    return True, "ok"


def _validate_case(case: Case) -> None:
    _run_cli(case.args)
    out_dir = Path(case.args[case.args.index("--out-dir") + 1]).resolve()
    for rel in case.expect_present:
        p = out_dir / rel
        if not p.exists():
            raise RuntimeError(f"[{case.name}] expected '{rel}' to exist.")
        if p.suffix in {".trj", ".xyz"}:
            n_frames = _frame_count(p)
            if n_frames <= 0:
                raise RuntimeError(f"[{case.name}] '{rel}' has no frames.")
    for rel in case.expect_absent:
        p = out_dir / rel
        if p.exists():
            raise RuntimeError(f"[{case.name}] expected '{rel}' to be absent.")


def main() -> int:
    runnable, reason = _check_runnable()
    if not runnable:
        print(f"[dump-smoke] skipped: {reason}")
        return 0

    with tempfile.TemporaryDirectory(prefix="pdb2_dump_smoke_") as td:
        base = Path(td)
        tiny_xyz = base / "tiny.xyz"
        _write_tiny_xyz(tiny_xyz)

        cases = [
            Case(
                name="opt_light_dump",
                args=[
                    "opt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    "--opt-mode",
                    "light",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "opt_light_dump"),
                ],
                expect_present=("optimization.trj",),
            ),
            Case(
                name="opt_heavy_dump",
                args=[
                    "opt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    "--opt-mode",
                    "heavy",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "opt_heavy_dump"),
                ],
                expect_present=("optimization.trj",),
            ),
            Case(
                name="tsopt_light_dump",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    "--opt-mode",
                    "light",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "tsopt_light_dump"),
                ],
                expect_present=("optimization.trj", "optimization_all.trj"),
            ),
            Case(
                name="tsopt_heavy_dump",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    "--opt-mode",
                    "heavy",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "tsopt_heavy_dump"),
                ],
                expect_present=("optimization.trj",),
                expect_absent=("optimization_all.trj",),
            ),
            Case(
                name="tsopt_heavy_nodump",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    "--opt-mode",
                    "heavy",
                    "--max-cycles",
                    "1",
                    "--no-dump",
                    "--out-dir",
                    str(base / "tsopt_heavy_nodump"),
                ],
                expect_present=(),
                expect_absent=("optimization.trj", "optimization_all.trj"),
            ),
            Case(
                name="tsopt_heavy_dump_legacy_bool",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    "--opt-mode",
                    "heavy",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "True",
                    "--out-dir",
                    str(base / "tsopt_heavy_dump_legacy_bool"),
                ],
                expect_present=("optimization.trj",),
            ),
        ]

        for case in cases:
            _validate_case(case)

    print(f"[dump-smoke] validated {len(cases)} cases.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
