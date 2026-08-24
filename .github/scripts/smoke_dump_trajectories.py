#!/usr/bin/env python3
"""Smoke-test trajectory dump behavior for opt/tsopt commands."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CLI_MODULE = "pdb2reaction"
CALC_FILE = Path(__file__).with_name("docs_harmonic_calc.py")
TIMEOUT_ENV = "PDB2REACTION_DUMP_CASE_TIMEOUT_SEC"
GENERIC_TIMEOUT_ENV = "DOCS_DUMP_CASE_TIMEOUT_SEC"
DEFAULT_CASE_TIMEOUT_SEC = 300.0


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


def _run_cli(args: list[str], timeout_sec: float | None = None) -> None:
    cmd = [sys.executable, "-m", CLI_MODULE, *args]
    try:
        proc = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=timeout_sec,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Command timed out after {timeout_sec} sec: {' '.join(cmd)}"
        ) from exc
    if proc.returncode != 0:
        tail = (proc.stdout + "\n" + proc.stderr)[-3000:]
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{tail}")


def _check_runnable() -> tuple[bool, str]:
    probe = subprocess.run(
        [sys.executable, "-m", CLI_MODULE, "opt", "--help"],
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


def _validate_case(case: Case, timeout_sec: float | None = None) -> None:
    _run_cli(case.args, timeout_sec=timeout_sec)
    out_dir = Path(case.args[case.args.index("--out-dir") + 1]).resolve()
    for rel in case.expect_present:
        p = out_dir / rel
        if not p.exists():
            raise RuntimeError(f"[{case.name}] expected '{rel}' to exist.")
        if p.suffix == ".xyz":
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
        raise RuntimeError(f"[dump-smoke] required trajectory smoke cannot run: {reason}")

    timeout_raw = ""
    timeout_env = None
    for env_name in (GENERIC_TIMEOUT_ENV, TIMEOUT_ENV):
        raw = os.environ.get(env_name, "").strip()
        if raw:
            timeout_raw = raw
            timeout_env = env_name
            break
    timeout_sec = float(timeout_raw) if timeout_raw else DEFAULT_CASE_TIMEOUT_SEC
    if timeout_sec <= 0:
        timeout_sec = None

    if timeout_sec is None:
        print("[dump-smoke] per-case timeout: disabled")
    else:
        timeout_source = f"env: {timeout_env}" if timeout_env else "default"
        print(f"[dump-smoke] per-case timeout: {timeout_sec:.1f}s ({timeout_source})")

    with tempfile.TemporaryDirectory(prefix="pdb2_dump_smoke_") as td:
        base = Path(td)
        tiny_xyz = base / "tiny.xyz"
        _write_tiny_xyz(tiny_xyz)
        calc_args = ["--calc-file", str(CALC_FILE)]

        cases = [
            Case(
                name="opt_grad_dump",
                args=[
                    "opt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    *calc_args,
                    "--opt-mode",
                    "grad",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "opt_grad_dump"),
                ],
                expect_present=("optimization_trj.xyz",),
            ),
            Case(
                name="opt_hess_dump",
                args=[
                    "opt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    *calc_args,
                    "--opt-mode",
                    "hess",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "opt_hess_dump"),
                ],
                expect_present=("optimization_trj.xyz",),
            ),
            Case(
                name="tsopt_grad_dump",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    *calc_args,
                    "--opt-mode",
                    "grad",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "tsopt_grad_dump"),
                ],
                expect_present=("optimization_trj.xyz", "optimization_all_trj.xyz"),
            ),
            Case(
                name="tsopt_hess_dump",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    *calc_args,
                    "--opt-mode",
                    "hess",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "--out-dir",
                    str(base / "tsopt_hess_dump"),
                ],
                expect_present=("optimization_trj.xyz",),
                expect_absent=("optimization_all_trj.xyz",),
            ),
            Case(
                name="tsopt_hess_nodump",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    *calc_args,
                    "--opt-mode",
                    "hess",
                    "--max-cycles",
                    "1",
                    "--no-dump",
                    "--out-dir",
                    str(base / "tsopt_hess_nodump"),
                ],
                expect_present=(),
                expect_absent=("optimization_trj.xyz", "optimization_all_trj.xyz"),
            ),
            Case(
                name="tsopt_hess_dump_legacy_bool",
                args=[
                    "tsopt",
                    "-i",
                    str(tiny_xyz),
                    "-q",
                    "0",
                    "-m",
                    "1",
                    *calc_args,
                    "--opt-mode",
                    "hess",
                    "--max-cycles",
                    "1",
                    "--dump",
                    "True",
                    "--out-dir",
                    str(base / "tsopt_hess_dump_legacy_bool"),
                ],
                expect_present=("optimization_trj.xyz",),
            ),
        ]

        for idx, case in enumerate(cases, start=1):
            print(f"[dump-smoke] case {idx}/{len(cases)}: {case.name}")
            started = time.perf_counter()
            _validate_case(case, timeout_sec=timeout_sec)
            elapsed = time.perf_counter() - started
            print(f"[dump-smoke] case ok: {case.name} ({elapsed:.1f}s)")

    print(f"[dump-smoke] validated {len(cases)} cases.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
