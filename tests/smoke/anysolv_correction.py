#!/usr/bin/env python3
"""Required real-model AnySolv lane; run in a scratch copy with the anysolv extra."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from ase import Atoms
from ase.optimize import BFGS
from pysisyphus.constants import AU2EV, BOHR2ANG

from pdb2reaction.backends import create_calculator, create_ase_calculator
from pdb2reaction.backends.anysolv import AnySolvCorrection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["uma", "mace"], default="uma")
    parser.add_argument("--moe-checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("anysolv_gate.json"))
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for the AnySolv real-model gate")
    if not args.moe_checkpoint.is_file():
        raise FileNotFoundError(args.moe_checkpoint)
    symbols = ["O", "H", "H"]
    x = np.array([[0.,0.,0.],[.96,.02,0.],[-.26,.93,.03]])
    records = []
    for checkpoint in ("model_compact", str(args.moe_checkpoint.resolve())):
        correction = AnySolvCorrection(solvent="water", checkpoint=checkpoint,
                                        device="cuda", dtype="float64")
        e,f,_ = correction.compute(symbols,x,need_forces=True)
        fd = np.zeros_like(x)
        step = 1e-4
        for i in range(x.size):
            plus,minus = x.copy(),x.copy()
            plus.flat[i] += step
            minus.flat[i] -= step
            fd.flat[i] = -(correction.compute(symbols,plus)[0]-correction.compute(symbols,minus)[0])/(2*step)
        error = float(np.max(np.abs(f-fd)))
        if error > 1e-4:
            raise AssertionError(f"AnySolv energy/force mismatch: {error}")
        record = {"checkpoint":checkpoint,"energy_ev":e,"force_error_ev_angstrom":error}
        if checkpoint == "model_compact":
            _,_,h = correction.compute(symbols,x,need_hessian=True)
            correction.hessian_step *= .5
            _,_,h_half = correction.compute(symbols,x,need_hessian=True)
            error_h = float(np.max(np.abs(h-h_half)))
            if error_h > 1e-3:
                raise AssertionError(f"AnySolv Hessian step sensitivity: {error_h}")
            record["hessian_step_error_ev_angstrom2"] = error_h
        records.append(record)

    common = dict(device="cuda",charge=0,spin=1,compile_model=False,
                  hessian_calc_mode="FiniteDifference",return_partial_hessian=False,
                  out_hess_torch=False,print_timing=False)
    if args.backend == "mace":
        common.update(model="mace-polar-1-s",default_dtype="float64")
    else:
        common.update(precision="fp64")
    base = create_calculator(args.backend,**common)
    corrected = create_calculator(args.backend,**common,solvent="water",solvent_model="anysolv")
    coords = (x/BOHR2ANG).reshape(-1)
    gas = base.get_hessian(symbols,coords)
    solution = corrected.get_hessian(symbols,coords)
    delta = AnySolvCorrection(solvent="water",device="cuda",dtype="float64")
    de,df,dh = delta.compute(symbols,x,need_hessian=True)
    np.testing.assert_allclose(solution["energy"]-gas["energy"],de/AU2EV,atol=1e-10)
    np.testing.assert_allclose(solution["forces"]-gas["forces"],df.reshape(-1)*BOHR2ANG/AU2EV,atol=1e-9)
    np.testing.assert_allclose(solution["hessian"]-gas["hessian"],dh*BOHR2ANG**2/AU2EV,atol=1e-8)
    atoms = Atoms(symbols,positions=x)
    atoms.calc = create_ase_calculator(args.backend,**common,solvent="water",solvent_model="anysolv")
    np.testing.assert_allclose(atoms.get_potential_energy(),solution["energy"]*AU2EV,atol=1e-7)
    np.testing.assert_allclose(atoms.get_forces(),solution["forces"].reshape(-1,3)*AU2EV/BOHR2ANG,atol=1e-6)
    optimizer = BFGS(atoms,logfile="anysolv_opt.log",trajectory="anysolv_opt.traj")
    optimizer.run(steps=2)
    if not np.isfinite(atoms.get_forces()).all():
        raise AssertionError("Non-finite force after AnySolv optimization steps")
    # The total Hessian is the one consumed by frequency workflows.
    frequencies = np.linalg.eigvalsh(np.asarray(solution["hessian"]))
    if not np.isfinite(frequencies).all():
        raise AssertionError("Non-finite total Hessian eigenvalues")
    args.output.write_text(json.dumps({"status":"passed","backend":args.backend,
                                      "models":records,"ase_optimizer_steps":optimizer.nsteps},indent=2)+"\n")
    print("PASS AnySolv compact/MoE derivatives, total PES, Hessian, ASE and optimizer")


if __name__ == "__main__":
    main()
