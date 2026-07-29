#!/usr/bin/env python3
"""Strict semantic checks for release-smoke outputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def require_finite(value: object, path: str = "root") -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise SystemExit(f"non-finite value at {path}")
        return
    if isinstance(value, dict):
        for key, child in value.items():
            require_finite(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            require_finite(child, f"{path}[{index}]")


def count_xyz_frames(path: Path) -> int:
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    frames = 0
    while index < len(lines):
        try:
            atoms = int(lines[index].strip())
        except ValueError as exc:
            raise SystemExit(f"invalid XYZ frame in {path}") from exc
        index += atoms + 2
        frames += 1
    if index != len(lines):
        raise SystemExit(f"truncated XYZ trajectory: {path}")
    return frames


def check_all(root: Path, require_thermo: bool, require_dft: bool) -> None:
    summary_path = root / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    require_finite(summary)

    # The pipeline must have run to the end.  This is the part a smoke can
    # guarantee, and the part a broken pipeline fails.
    execution = summary.get("execution_status")
    if execution != "completed":
        raise SystemExit(f"all pipeline did not complete: execution_status={execution!r}")

    # The scientific outcome is throttled here on purpose (`--max-cycles`,
    # `--thresh gau_loose`), so demanding `success` would demand that the MLIP
    # find publication-grade minima on a deliberately cheap run.  What IS
    # required is that the reported outcome be TRUE: either a real success, or a
    # `partial` that says what is missing.  A silent degradation, or a `failed`,
    # still fails the lane -- and so does a `partial` with no stated reason,
    # which is how a missing outcome record would look.
    status = summary.get("status")
    scientific = summary.get("scientific_status")
    if status not in ("success", "partial") or scientific not in ("success", "partial"):
        raise SystemExit(
            f"all summary is neither success nor an explained partial: "
            f"status={status!r} scientific_status={scientific!r}"
        )
    reasons = summary.get("scientific_status_reasons") or []
    if scientific == "partial" and not reasons:
        raise SystemExit("all summary is partial but states no reason")
    if scientific == "success" and reasons:
        raise SystemExit(f"all summary claims success yet states reasons: {reasons}")

    segments = summary.get("post_segments") or []
    if not segments:
        raise SystemExit("all summary has no post-TS segments")
    for segment in segments:
        if int((segment.get("ts_imag") or {}).get("n_imag", -1)) != 1:
            raise SystemExit("post-TS segment is not a first-order saddle")
        tag = str(segment["tag"])
        # Use the directory the producer published for this segment. ``tag`` is
        # the path-search segment id (``seg_%03d``) and is NOT a directory name;
        # the per-segment deliverable dir is indexed separately (``seg_%02d``).
        # Rebuilding the path from ``tag`` asserts an invariant the product has
        # never promised, and fails on the very first segment.
        published = segment.get("post_dir")
        if not published:
            raise SystemExit(f"segment {tag!r} published no post_dir")
        segment_root = Path(published)
        if not segment_root.is_dir():
            raise SystemExit(f"published post_dir is not a directory: {segment_root}")
        for direction in ("forward", "backward"):
            trajectory = segment_root / "irc" / f"{direction}_irc_trj.xyz"
            if not trajectory.is_file() or count_xyz_frames(trajectory) < 2:
                raise SystemExit(f"missing/nontrivial IRC branch: {trajectory}")
        if require_thermo:
            # Thermochemistry needs R and P to be clean minima, which a throttled
            # run cannot promise.  So: present, or explicitly accounted for in
            # the reasons. Absent and unexplained is an outcome-record failure.
            missing = [
                state
                for state in ("R", "TS", "P")
                if not (segment_root / "freq" / state / "thermoanalysis.yaml").is_file()
                or (segment_root / "freq" / state / "thermoanalysis.yaml").stat().st_size == 0
            ]
            if missing:
                explained = any("thermochemistry" in str(r).lower() for r in reasons)
                if not explained:
                    raise SystemExit(
                        f"thermochemistry missing for {missing} and no reason says so; "
                        f"reasons={reasons}"
                    )
        if require_dft:
            dft_files = [path for path in segment_root.rglob("*") if path.is_file() and "dft" in path.as_posix().lower()]
            if not dft_files:
                raise SystemExit(f"missing DFT artifacts below {segment_root}")


def check_tsopt(root: Path) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("status") != "converged":
        raise SystemExit(f"TS optimization did not converge: {payload.get('status')!r}")
    if int(payload.get("n_imaginary_modes", -1)) != 1:
        raise SystemExit("TS optimization did not produce exactly one imaginary mode")


def check_opt_config(
    root: Path,
    *,
    expected_charge: int,
    expected_model: str,
    expected_max_cycles: int,
    expected_precision: str,
) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("command") != "opt":
        raise SystemExit(f"unexpected result command: {payload.get('command')!r}")
    if int(payload.get("charge", 10**9)) != expected_charge:
        raise SystemExit(f"runtime charge did not come from YAML: {payload.get('charge')!r}")
    if payload.get("model") != expected_model:
        raise SystemExit(f"runtime model did not come from YAML: {payload.get('model')!r}")
    if int(payload.get("max_cycles", -1)) != expected_max_cycles:
        raise SystemExit(f"CLI max-cycles did not override YAML: {payload.get('max_cycles')!r}")
    if payload.get("mlip_precision") != expected_precision:
        raise SystemExit(
            f"runtime precision did not come from YAML: {payload.get('mlip_precision')!r}"
        )


def check_sp_hessian(root: Path) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("status") != "ok" or payload.get("backend") != "custom":
        raise SystemExit("custom calculator SP did not complete successfully")
    if payload.get("custom_calculator") != "harmonic_calc.py:get_calculator":
        raise SystemExit(
            f"custom calculator provenance is wrong: {payload.get('custom_calculator')!r}"
        )
    if payload.get("mlip_model") != payload.get("custom_calculator"):
        raise SystemExit("custom calculator was not recorded as mlip_model")
    if payload.get("mlip_precision") is not None:
        raise SystemExit("custom calculator incorrectly reports a built-in MLIP precision")
    hessian = np.load(root / "hessian.npy", allow_pickle=False)
    forces = np.load(root / "forces.npy", allow_pickle=False)
    if forces.ndim != 2 or forces.shape[1] != 3 or not np.isfinite(forces).all():
        raise SystemExit(f"custom calculator forces are invalid: {forces.shape}")
    if Path(str(payload.get("forces_path", ""))).resolve() != (root / "forces.npy").resolve():
        raise SystemExit("custom calculator result does not identify forces.npy")
    if Path(str(payload.get("hessian_path", ""))).resolve() != (root / "hessian.npy").resolve():
        raise SystemExit("custom calculator result does not identify hessian.npy")
    expected = 3 * int(payload.get("n_atoms", 0) or 0)
    if expected == 0:
        # The SP schema predates n_atoms; infer it from forces while retaining
        # an exact square-matrix contract.
        expected = 3 * int(np.atleast_2d(forces).shape[0])
    if hessian.shape != (expected, expected):
        raise SystemExit(f"unexpected Hessian shape: {hessian.shape}, expected {(expected, expected)}")
    if not np.isfinite(hessian).all():
        raise SystemExit("custom finite-difference Hessian contains non-finite values")
    if not np.allclose(hessian, hessian.T, atol=1.0e-8, rtol=0.0):
        raise SystemExit("custom finite-difference Hessian is not symmetric")
    if not np.any(np.abs(hessian) > 1.0e-8):
        raise SystemExit("custom finite-difference Hessian is identically zero")
    from ase.units import Bohr, Hartree

    atoms = expected // 3
    projector = np.eye(atoms) - np.ones((atoms, atoms)) / atoms
    exact = np.kron(projector, np.eye(3)) * (Bohr**2 / Hartree)
    if not np.allclose(hessian, exact, atol=1.0e-6, rtol=1.0e-5):
        error = float(np.max(np.abs(hessian - exact)))
        raise SystemExit(f"custom finite-difference Hessian has wrong scale/sign (max error={error:.3e})")


def check_irc_never_stop(root: Path) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("status") != "completed" or payload.get("never_stop") is not True:
        raise SystemExit("IRC never-stop run did not complete with the requested mode")
    if int(payload.get("never_stop_energy_bypasses", 0)) < 1:
        raise SystemExit("IRC never-stop did not bypass an actual energy stop")
    if int(payload.get("n_frames_forward", 0)) < 2 or int(payload.get("n_frames_backward", 0)) < 2:
        raise SystemExit("IRC never-stop did not produce both nontrivial branches")


def check_provenance(
    root: Path,
    *,
    expected_backend: str,
    expected_model: str,
    expected_precision: str,
) -> None:
    summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
    expected = {
        "mlip_backend": expected_backend,
        "mlip_model": expected_model,
        "mlip_precision": expected_precision,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise SystemExit(
                f"incorrect {key} provenance: {summary.get(key)!r}, expected {value!r}"
            )
    text = (root / "summary.log").read_text(encoding="utf-8")
    if f"MLIP backend       : {expected_backend}" not in text:
        raise SystemExit("summary.log omits the resolved MLIP backend")
    if f"MLIP model         : {expected_model}" not in text:
        raise SystemExit("summary.log omits the resolved MLIP model")
    for stale in ("UMA model", "DFT//UMA", "energy_diagram_UMA"):
        if stale in text:
            raise SystemExit(f"summary.log contains stale backend label: {stale}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "kind",
        choices=("all", "tsopt", "opt-config", "sp-hessian", "irc-never-stop", "provenance"),
    )
    parser.add_argument("root", type=Path)
    parser.add_argument("--require-thermo", action="store_true")
    parser.add_argument("--require-dft", action="store_true")
    parser.add_argument("--expected-charge", type=int)
    parser.add_argument("--expected-model")
    parser.add_argument("--expected-max-cycles", type=int)
    parser.add_argument("--expected-precision")
    parser.add_argument("--expected-backend")
    args = parser.parse_args()
    if args.kind == "all":
        check_all(args.root, args.require_thermo, args.require_dft)
    elif args.kind == "tsopt":
        check_tsopt(args.root)
    elif args.kind == "opt-config":
        if None in (
            args.expected_charge,
            args.expected_model,
            args.expected_max_cycles,
            args.expected_precision,
        ):
            parser.error(
                "opt-config requires --expected-charge, --expected-model, "
                "--expected-max-cycles, and --expected-precision"
            )
        check_opt_config(
            args.root,
            expected_charge=args.expected_charge,
            expected_model=args.expected_model,
            expected_max_cycles=args.expected_max_cycles,
            expected_precision=args.expected_precision,
        )
    elif args.kind == "sp-hessian":
        check_sp_hessian(args.root)
    elif args.kind == "irc-never-stop":
        check_irc_never_stop(args.root)
    else:
        if None in (args.expected_backend, args.expected_model, args.expected_precision):
            parser.error(
                "provenance requires --expected-backend, --expected-model, "
                "and --expected-precision"
            )
        check_provenance(
            args.root,
            expected_backend=args.expected_backend,
            expected_model=args.expected_model,
            expected_precision=args.expected_precision,
        )


if __name__ == "__main__":
    main()
