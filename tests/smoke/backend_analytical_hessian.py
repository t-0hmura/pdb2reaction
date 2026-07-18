#!/usr/bin/env python3
"""Strict real-model analytical-Hessian smoke for one MLIP backend."""

from __future__ import annotations

import argparse
import json

import numpy as np
import torch

from pysisyphus.constants import ANG2BOHR


def _as_numpy(value) -> np.ndarray:
    if isinstance(value, torch.Tensor):
        value = value.detach().cpu().numpy()
    return np.asarray(value, dtype=np.float64)


def _metrics(analytical: np.ndarray, finite_difference: np.ndarray) -> dict[str, float]:
    difference = analytical - finite_difference
    fd_norm = max(float(np.linalg.norm(finite_difference)), 1.0e-12)
    fd_max = max(float(np.max(np.abs(finite_difference))), 1.0e-12)
    return {
        "relative_frobenius_error": float(np.linalg.norm(difference)) / fd_norm,
        "relative_max_error": float(np.max(np.abs(difference))) / fd_max,
        "analytical_asymmetry": float(np.max(np.abs(analytical - analytical.T))),
        "fd_asymmetry": float(np.max(np.abs(finite_difference - finite_difference.T))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("backend", choices=("uma", "orb", "mace", "aimnet2"))
    parser.add_argument("--max-relative-frobenius", type=float, default=0.20)
    parser.add_argument("--max-relative-element", type=float, default=0.35)
    args = parser.parse_args()

    if not torch.cuda.is_available():
        raise SystemExit("CUDA is required for the real-model analytical-Hessian lane.")

    common = dict(
        charge=0,
        spin=1,
        device="cuda",
        hessian_calc_mode="Analytical",
        return_partial_hessian=False,
        hessian_double=True,
        out_hess_torch=False,
        print_timing=False,
    )
    if args.backend == "uma":
        from pdb2reaction.backends.uma import UMACalculator

        calculator = UMACalculator(
            model="uma-s-1p2",
            task_name="omol",
            precision="fp32",
            workers=1,
            **common,
        )
    elif args.backend == "orb":
        from pdb2reaction.backends.orb import OrbCalculator

        calculator = OrbCalculator(
            model="orb_v3_conservative_omol",
            precision="float64",
            compile_model=False,
            **common,
        )
    elif args.backend == "mace":
        from pdb2reaction.backends.mace import MACECalculator

        calculator = MACECalculator(
            model="MACE-OMOL-0",
            default_dtype="float64",
            **common,
        )
    else:
        from pdb2reaction.backends.aimnet2 import AIMNet2Calculator

        calculator = AIMNet2Calculator(model="aimnet2", **common)

    elements = ["H", "H"]
    coordinates_bohr = np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]], dtype=np.float64
    ).reshape(-1) * ANG2BOHR

    analytical_result = calculator.get_hessian(elements, coordinates_bohr)
    analytical = _as_numpy(analytical_result["hessian"]).reshape(6, 6)
    calculator.hessian_calc_mode = "FiniteDifference"
    fd_result = calculator.get_hessian(elements, coordinates_bohr)
    finite_difference = _as_numpy(fd_result["hessian"]).reshape(6, 6)

    if not np.isfinite(analytical).all() or not np.isfinite(finite_difference).all():
        raise SystemExit("Hessian contains non-finite values.")
    metrics = _metrics(analytical, finite_difference)
    if metrics["analytical_asymmetry"] > 1.0e-10:
        raise SystemExit(f"Analytical Hessian is not symmetric: {metrics}")
    if metrics["relative_frobenius_error"] > args.max_relative_frobenius:
        raise SystemExit(f"Analytical/FD Frobenius mismatch: {metrics}")
    if metrics["relative_max_error"] > args.max_relative_element:
        raise SystemExit(f"Analytical/FD element mismatch: {metrics}")

    print(json.dumps({"backend": args.backend, "shape": [6, 6], **metrics}, sort_keys=True))


if __name__ == "__main__":
    main()
