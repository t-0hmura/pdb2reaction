# pdb2reaction/ts_opt.py

"""
pdb2reaction ts_opt — Transition-state optimization CLI
====================================================================

Usage (CLI)
-----
    pdb2reaction ts_opt -i INPUT.(pdb|xyz|trj) -q CHARGE -s SPIN
                        [--opt-mode light|heavy]
                        [--freeze-links True|False]
                        [--max-cycles N]
                        [--dump True|False]
                        [--out-dir DIR]
                        [--args-yaml args.yaml]
                        [--hessian-calc-mode Analytical|FiniteDifference]

Examples::
    # Minimal (recommended to always specify charge and spin)
    pdb2reaction ts_opt -i ts_cand.pdb -q 0 -s 1 --opt-mode light --out-dir ./result_ts_opt/

    # Light mode (HessianDimer) with YAML overrides and finite-difference Hessian
    pdb2reaction ts_opt -i ts_cand.pdb -q 0 -s 1 \
      --freeze-links True --opt-mode light --max-cycles 10000 --dump False \
      --out-dir ./result_ts_opt/ --args-yaml ./args.yaml \
      --hessian-calc-mode FiniteDifference

    # Heavy mode (RS-I-RFO) using YAML
    pdb2reaction ts_opt -i ts_cand.pdb -q 0 -s 1 --opt-mode heavy \
      --args-yaml ./args.yaml --out-dir ./result_ts_opt/


Description
-----
Transition-state optimization with two modes:
- light: HessianDimer TS search with periodic Hessian updates and a memory-efficient
  flatten loop to remove excess imaginary modes.
- heavy: RS-I-RFO Hessian-based TS optimizer.

Configuration is driven by YAML overrides for sections: geom, calc, opt, hessian_dimer,
and rsirfo. The hessian_dimer section accepts nested dimer and lbfgs dictionaries that
are forwarded to the respective pysisyphus components.

Structures are loaded via pysisyphus.helpers.geom_loader (PDB/XYZ/TRJ/etc.). The UMA
calculator (pdb2reaction.uma_pysis) provides energies, gradients, and Hessians. For PDB
inputs, optimization trajectories and final_geometry.xyz are also converted to PDB. The
final imaginary mode is always written as both .trj and .pdb.

Key behaviors and algorithmic notes:
- Direction selection: choose which imaginary mode to follow using root
  (0 = most negative, 1 = second-most, ...), consistent with RS-I-RFO. For root == 0,
  the implementation prefers torch.lobpcg for the lowest eigenpair and falls back to
  torch.linalg.eigh as needed.
- PHVA and TR projection: an active-degree-of-freedom (PHVA) subspace with translation/
  rotation (TR) projection reduces GPU memory use. This mirrors freq.py and respects
  freeze_atoms (particularly in light mode).
- Flatten loop (light mode): a single exact active-subspace Hessian is kept in memory.
  After the initial Dimer stage, one exact Hessian is evaluated and the PHVA block
  extracted. Between geometry updates a Bofill update (SR1/MS and PSB blend) is applied
  to the active-space Cartesian Hessian using displacements delta = Dx and gradient
  differences y = Dg. Each iteration: estimate imaginary modes from the updated active
  Hessian (mass-weighted, TR-projected), perform a flatten step, update the Hessian via
  Bofill, refresh the Dimer direction, run a Dimer LBFGS segment, and apply a final
  Bofill update. Continue until only one imaginary mode remains, then compute a final
  exact Hessian for frequency analysis.
- UMA integration: freeze_atoms are passed explicitly, and UMA defaults to
  return_partial_hessian=True so PHVA returns an active-subspace block when applicable.
  --hessian-calc-mode selects Analytical or FiniteDifference Hessians.


Outputs (& Directory Layout)
-----
DIR (default: ./result_ts_opt/)
  ├── final_geometry.xyz
  ├── final_geometry.pdb                 # written if input was PDB
  ├── optimization_all.trj               # light mode, when --dump True
  ├── optimization_all.pdb               # PDB conversion of the above (PDB input)
  ├── optimization.trj                   # heavy mode (RS-I-RFO)
  ├── optimization.pdb                   # PDB conversion of the above (PDB input)
  ├── vib/
  │     ├── final_imag_mode_+/-XXXX.Xcm-1.trj
  │     └── final_imag_mode_+/-XXXX.Xcm-1.pdb
  └── .dimer_mode.dat                    # internal dimer orientation seed (light mode)


Notes:
-----
- Always set charge (-q) and spin (-s) to avoid unphysical conditions.
- --opt-mode light runs the HessianDimer with periodic Hessian-based direction refresh;
  --opt-mode heavy runs RS-I-RFO.
- --freeze-links is PDB-only and freezes parents of link hydrogens; these indices are
  merged into geom.freeze_atoms and propagated to calc.freeze_atoms for UMA.
- --hessian-calc-mode overrides calc.hessian_calc_mode from YAML. Finite-difference
  Hessians will honor the active subspace when freeze_atoms are present.
- Convergence and stepping behavior are configurable via YAML in hessian_dimer.lbfgs,
  hessian_dimer.dimer, opt, and rsirfo sections (e.g., thresholds, trust radii, memory).
- Imaginary-mode detection uses a default threshold of ~5 cm^-1; the primary mode written
  at the end is chosen via root.
"""

from __future__ import annotations

import sys
import math
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List

import click
import numpy as np
import torch
from ase import Atoms
from ase.io import write
from ase.data import atomic_masses
import ase.units as units
import yaml
import time

# ---------------- pysisyphus / pdb2reaction imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AMU2AU, AU2EV
from pysisyphus.calculators.Dimer import Dimer  # Dimer calculator (orientation-projected forces)

# RS-I-RFO optimizer (transition state)
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer  # type: ignore

# local helpers from pdb2reaction
from .uma_pysis import uma_pysis
from .utils import (
    convert_xyz_to_pdb,
    freeze_links as detect_freeze_links,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    merge_freeze_atom_indices,
)


def _norm_opt_mode(mode: str) -> str:
    m = (mode or "").strip().lower()
    if m in ("light", "lbfgs", "dimer", "simple", "simpledimer", "hessian_dimer"):
        return "light"
    if m in ("heavy", "rfo", "rsirfo", "rs-i-rfo"):
        return "heavy"
    raise click.BadParameter(f"Unknown --opt-mode '{mode}'. Use: light|heavy")


# ===================================================================
#               Mass-weighted projection & vib analysis
# ===================================================================

def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)


def _build_tr_basis(coords_bohr_t: torch.Tensor,
                    masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Mass‐weighted translation/rotation basis (Tx,Ty,Tz,Rx,Ry,Rz), shape (3N, r<=6).
    """
    device, dtype = coords_bohr_t.device, coords_bohr_t.dtype
    N = coords_bohr_t.shape[0]
    m_au = masses_au_t.to(dtype=dtype, device=device)
    m_sqrt = torch.sqrt(m_au).reshape(-1, 1)

    # COM
    com = (m_au.reshape(-1, 1) * coords_bohr_t).sum(0) / m_au.sum()
    x = coords_bohr_t - com

    eye3 = torch.eye(3, dtype=dtype, device=device)
    cols = []
    # Translations
    for i in range(3):
        cols.append((eye3[i].repeat(N, 1) * m_sqrt).reshape(-1, 1))
    # Rotations
    for i in range(3):
        rot = torch.cross(x, eye3[i].expand_as(x), dim=1) * m_sqrt
        cols.append(rot.reshape(-1, 1))
    return torch.cat(cols, dim=1)


def _tr_orthonormal_basis(coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor,
                          rtol: float = 1e-12) -> Tuple[torch.Tensor, int]:
    """
    Orthonormalize TR basis in mass-weighted space by SVD. Returns (Q, rank).
    """
    B = _build_tr_basis(coords_bohr_t, masses_au_t)
    U, S, Vh = torch.linalg.svd(B, full_matrices=False)
    r = int((S > rtol * S.max()).sum().item())
    Q = U[:, :r]
    del B, S, Vh, U
    return Q, r


# ---- in-place mass weighting ----
def _mass_weighted_hessian(H_t: torch.Tensor,
                           masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Return Hmw = M^{-1/2} H M^{-1/2} (in-place on H_t).
    """
    dtype, device = H_t.dtype, H_t.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)
        H_t.mul_(inv_sqrt_m_row)
        H_t.mul_(inv_sqrt_m_col)
        del masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row
        return H_t


def _mw_projected_hessian_inplace(H_t: torch.Tensor,
                                  coords_bohr_t: torch.Tensor,
                                  masses_au_t: torch.Tensor,
                                  freeze_idx: Optional[List[int]] = None) -> torch.Tensor:
    """
    Mass-weight H in-place, optionally restrict to active DOF subspace (PHVA) and
    project out TR motions (in that subspace), also in-place. No explicit symmetrization.
    Returns the (possibly reduced) Hessian to be diagonalized.
    """
    dtype, device = H_t.dtype, H_t.device
    with torch.no_grad():
        N = coords_bohr_t.shape[0]
        if freeze_idx is not None and len(freeze_idx) > 0:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen]
            if len(active_idx) == 0:
                raise RuntimeError("All atoms are frozen; no active DOF left for TR projection.")
            # mass-weight first
            H_t = _mass_weighted_hessian(H_t, masses_au_t)
            # take active DOF submatrix
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            H_t = H_t[mask_dof][:, mask_dof]
            # TR basis and projection in active subspace (in-place)
            coords_act = coords_bohr_t[active_idx, :]
            masses_act = masses_au_t[active_idx]
            Q, _ = _tr_orthonormal_basis(coords_act, masses_act)  # (3N_act, r)
            Qt = Q.T
            QtH = Qt @ H_t
            H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
            H_t.addmm_((QtH.T), Qt, beta=1.0, alpha=-1.0)
            QtHQ = QtH @ Q
            H_t.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
            del Q, Qt, QtH, QtHQ, mask_dof, coords_act, masses_act, active_idx, frozen
        else:
            # Full DOF: mass-weight + TR projection in-place
            H_t = _mass_weighted_hessian(H_t, masses_au_t)
            Q, _ = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)  # (3N, r)
            Qt = Q.T
            QtH = Qt @ H_t
            H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
            H_t.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)
            QtHQ = QtH @ Q
            H_t.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
            del Q, Qt, QtH, QtHQ
        if torch.cuda.is_available() and device.type == "cuda":
            torch.cuda.empty_cache()
        return H_t


def _self_check_tr_projection(H_t: torch.Tensor,
                              coords_bohr_t: torch.Tensor,
                              masses_au_t: torch.Tensor,
                              freeze_idx: Optional[List[int]] = None) -> Tuple[float, float, int]:
    """
    (Low-VRAM) Skipped heavy clone-based check. Keep signature for compatibility.
    """
    # To conserve VRAM, do not allocate a full clone for checks.
    # Return zeros and rank estimate based on current TR basis only.
    with torch.no_grad():
        N = coords_bohr_t.shape[0]
        if freeze_idx and len(freeze_idx) > 0:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen]
            if len(active_idx) == 0:
                return 0.0, 0.0, 0
            coords_act = coords_bohr_t[active_idx, :]
            masses_act = masses_au_t[active_idx]
            _, r = _tr_orthonormal_basis(coords_act, masses_act)
            return 0.0, 0.0, r
        else:
            _, r = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)
            return 0.0, 0.0, r


def _mode_direction_by_root(H_t: torch.Tensor,
                            coords_bohr_t: torch.Tensor,
                            masses_au_t: torch.Tensor,
                            root: int = 0,
                            freeze_idx: Optional[List[int]] = None) -> np.ndarray:
    """
    Get the eigenvector (Cartesian space) corresponding to the `root`-th most negative
    eigenvalue (root=0: most negative) of the mass-weighted, TR-projected Hessian.
    PHVA (active-subspace) is applied if freeze_idx is provided: frozen DOFs are zero.
    root==0 prefers torch.lobpcg; fallback to eigh (UPLO='U').
    """
    with torch.no_grad():
        # In-place: mass weight + (active-subspace) TR projection
        Hmw_proj = _mw_projected_hessian_inplace(H_t, coords_bohr_t, masses_au_t, freeze_idx=freeze_idx)

        # Solve eigenproblem in the (possibly reduced) space
        if int(root) == 0:
            try:
                w, v_mw_sub = torch.lobpcg(Hmw_proj, k=1, largest=False)
                u_mw_sub = v_mw_sub[:, 0]
            except Exception:
                evals_f, evecs_f = torch.linalg.eigh(Hmw_proj, UPLO="U")
                u_mw_sub = evecs_f[:, torch.argmin(evals_f)]
                del evals_f, evecs_f
        else:
            evals, evecs_mw = torch.linalg.eigh(Hmw_proj, UPLO="U")  # ascending
            neg = (evals < 0)
            neg_inds = torch.nonzero(neg, as_tuple=False).view(-1)
            if neg_inds.numel() == 0:
                pick = int(torch.argmin(evals).item())
            else:
                k = max(0, min(int(root), neg_inds.numel() - 1))
                pick = int(neg_inds[k].item())
            u_mw_sub = evecs_mw[:, pick]
            del evals, evecs_mw

        # Embed back to full 3N (frozen DOF as zeros) if we solved in subspace
        N = coords_bohr_t.shape[0]
        if freeze_idx and len(freeze_idx) > 0:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Hmw_proj.device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            u_mw_full = torch.zeros(3 * N, dtype=Hmw_proj.dtype, device=Hmw_proj.device)
            u_mw_full[mask_dof] = u_mw_sub
            u_mw = u_mw_full
            del mask_dof, frozen
        else:
            u_mw = u_mw_sub

        # Convert mass-weighted → Cartesian & normalize
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=Hmw_proj.dtype, device=Hmw_proj.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        v = inv_sqrt_m * u_mw
        v = v / torch.linalg.norm(v)
        mode = v.reshape(-1, 3).detach().cpu().numpy()

        del masses_amu_t, m3, inv_sqrt_m, v, u_mw, u_mw_sub
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return mode


def _calc_full_hessian_torch(geom, uma_kwargs: dict, device: torch.device) -> torch.Tensor:
    """
    UMA calculator producing (possibly partial) Hessian as torch.Tensor in Hartree/Bohr^2 (3N or 3N_act).
    """
    # respect caller's out_hess_torch preference (force True here)
    kw = dict(uma_kwargs or {})
    kw["out_hess_torch"] = True
    calc = uma_pysis(**kw)
    H_t = calc.get_hessian(geom.atoms, geom.coords)["hessian"].to(device=device)
    return H_t


def _calc_energy(geom, uma_kwargs: dict) -> float:
    calc = uma_pysis(out_hess_torch=False, **uma_kwargs)
    geom.set_calculator(calc)
    E = geom.energy
    geom.set_calculator(None)
    return E


def _calc_gradient(geom, uma_kwargs: dict) -> np.ndarray:
    """
    Return true Cartesian gradient (shape 3N,) in Hartree/Bohr.
    """
    calc = uma_pysis(out_hess_torch=False, **uma_kwargs)
    geom.set_calculator(calc)
    g = np.array(geom.gradient, dtype=float).reshape(-1)
    geom.set_calculator(None)
    return g


def _frequencies_cm_and_modes(H_t: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              tol: float = 1e-6,
                              freeze_idx: Optional[List[int]] = None) -> Tuple[np.ndarray, torch.Tensor]:
    """
    In-place PHVA/TR projection (active-subspace if freeze_idx) and diagonalization.
    Returns:
      freqs_cm : (nmode,) numpy (negatives are imaginary)
      modes    : (nmode, 3N) torch (mass-weighted eigenvectors embedded to full 3N)
    """
    with torch.no_grad():
        Z = np.array(atomic_numbers, dtype=int)
        N = int(len(Z))
        masses_amu = np.array([atomic_masses[z] for z in Z])  # amu
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        # in-place mass-weight + (active-subspace) TR projection
        Hmw = _mw_projected_hessian_inplace(H_t, coords_bohr_t, masses_au_t, freeze_idx=freeze_idx)

        # eigensolve (upper triangle)
        omega2, Vsub = torch.linalg.eigh(Hmw, UPLO="U")

        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]
        Vsub = Vsub[:, sel]  # (3N_act or 3N, nsel)

        # embed modes to full 3N
        if freeze_idx and len(freeze_idx) > 0:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Hmw.device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Hmw.dtype, device=Hmw.device)
            modes[:, mask_dof] = Vsub.T
            del mask_dof, frozen
        else:
            modes = Vsub.T  # (nsel, 3N)

        # convert to cm^-1
        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()

        del omega2, hnu, Vsub
        if torch.cuda.is_available() and H_t.is_cuda:
            torch.cuda.empty_cache()
        return freqs_cm, modes


def _frequencies_cm_only(H_t: torch.Tensor,
                         atomic_numbers: List[int],
                         coords_bohr: np.ndarray,
                         device: torch.device,
                         tol: float = 1e-6,
                         freeze_idx: Optional[List[int]] = None) -> np.ndarray:
    """
    Frequencies only (PHVA/TR in-place; no eigenvectors) for quick checks.
    """
    with torch.no_grad():
        Z = np.array(atomic_numbers, dtype=int)
        masses_amu = np.array([atomic_masses[z] for z in Z])  # amu
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        Hmw = _mw_projected_hessian_inplace(H_t, coords_bohr_t, masses_au_t, freeze_idx=freeze_idx)
        omega2 = torch.linalg.eigvalsh(Hmw, UPLO="U")

        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]

        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()

        del omega2, hnu, sel
        if torch.cuda.is_available() and H_t.is_cuda:
            torch.cuda.empty_cache()
        return freqs_cm


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            out_pdb: Path,
                            amplitude_ang: float = 0.25,
                            n_frames: int = 20,
                            comment: str = "imag mode") -> None:
    """
    Write a single imaginary mode animation both as .trj (XYZ-like) and .pdb.
    """
    ref_ang = geom.coords.reshape(-1, 3) * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # .trj (XYZ-like concatenation)
    try:
        from pysisyphus.xyzloader import make_trj_str  # type: ignore
        amp_ang = amplitude_ang
        steps = np.sin(2.0 * np.pi * np.arange(n_frames) / n_frames)[:, None, None] * (amp_ang * mode[None, :, :])
        traj_ang = ref_ang[None, :, :] + steps  # (T,N,3) in Å
        traj_bohr = traj_ang.reshape(n_frames, -1, 3) * ANG2BOHR
        comments = [f"{comment}  frame={i+1}/{n_frames}" for i in range(n_frames)]
        trj_str = make_trj_str(geom.atoms, traj_bohr, comments=comments)
        out_trj.write_text(trj_str, encoding="utf-8")
    except Exception:
        with out_trj.open("w", encoding="utf-8") as f:
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                coords = ref_ang + phase * amplitude_ang * mode
                f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
                for sym, (x, y, z) in zip(geom.atoms, coords):
                    f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")

    # .pdb (MODEL/ENDMDL)
    atoms0 = Atoms(geom.atoms, positions=ref_ang, pbc=False)
    for i in range(n_frames):
        phase = np.sin(2.0 * np.pi * i / n_frames)
        ai = atoms0.copy()
        ai.set_positions(ref_ang + phase * amplitude_ang * mode)
        write(out_pdb, ai, append=(i != 0))


# ===================================================================
#            Active-subspace helpers & Bofill update
# ===================================================================

def _active_indices(N: int, freeze_idx: Optional[List[int]]) -> List[int]:
    if not freeze_idx:
        return list(range(N))
    fz = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
    return [i for i in range(N) if i not in fz]


def _active_mask_dof(N: int, freeze_idx: Optional[List[int]]) -> np.ndarray:
    mask = np.ones(3 * N, dtype=bool)
    if freeze_idx:
        for i in freeze_idx:
            if 0 <= int(i) < N:
                mask[3 * int(i):3 * int(i) + 3] = False
    return mask


def _extract_active_block(H_full: torch.Tensor, mask_dof: np.ndarray) -> torch.Tensor:
    """
    Return the active-DOF block as a torch.Tensor sharing device/dtype.
    """
    device = H_full.device
    m = torch.as_tensor(mask_dof, device=device, dtype=torch.bool)
    return H_full[m][:, m].clone()


def _embed_active_vector(vec_act: torch.Tensor,
                         mask_dof: np.ndarray,
                         total_3N: int) -> torch.Tensor:
    """
    Embed a (3N_act,) vector back to full (3N,) with zeros on frozen DOFs.
    """
    device = vec_act.device
    dtype = vec_act.dtype
    full = torch.zeros(total_3N, device=device, dtype=dtype)
    m = torch.as_tensor(mask_dof, device=device, dtype=torch.bool)
    full[m] = vec_act
    return full


def _mw_tr_project_active_inplace(H_act: torch.Tensor,
                                  coords_act_t: torch.Tensor,
                                  masses_act_au_t: torch.Tensor) -> torch.Tensor:
    """
    Mass-weight & project TR in the *active* subspace (in-place; no explicit symmetrization).
    """
    with torch.no_grad():
        # mass-weight
        masses_amu_t = (masses_act_au_t / AMU2AU).to(dtype=H_act.dtype, device=H_act.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)
        H_act.mul_(inv_sqrt_m_row)
        H_act.mul_(inv_sqrt_m_col)
        # TR basis & projection
        Q, _ = _tr_orthonormal_basis(coords_act_t, masses_act_au_t)  # (3N_act, r)
        Qt = Q.T
        QtH = Qt @ H_act
        H_act.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
        H_act.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)
        QtHQ = QtH @ Q
        H_act.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
        del masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row, Q, Qt, QtH, QtHQ
        return H_act


def _frequencies_from_Hact(H_act: torch.Tensor,
                           atomic_numbers: List[int],
                           coords_bohr: np.ndarray,
                           active_idx: List[int],
                           device: torch.device,
                           tol: float = 1e-6) -> np.ndarray:
    """
    Frequencies (cm^-1) computed from active-block Hessian with active-space TR projection.
    """
    with torch.no_grad():
        coords_act = torch.as_tensor(coords_bohr.reshape(-1, 3)[active_idx, :], dtype=H_act.dtype, device=device)
        masses_act_au = torch.as_tensor([atomic_masses[int(z)] * AMU2AU
                                         for z in np.array(atomic_numbers, int)[active_idx]],
                                        dtype=H_act.dtype, device=device)
        Hmw = H_act.clone()
        _mw_tr_project_active_inplace(Hmw, coords_act, masses_act_au)
        omega2 = torch.linalg.eigvalsh(Hmw, UPLO="U")
        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]
        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()
        del coords_act, masses_act_au, Hmw, omega2, hnu, sel
        if torch.cuda.is_available() and H_act.is_cuda:
            torch.cuda.empty_cache()
        return freqs_cm


def _modes_from_Hact_embedded(H_act: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              active_idx: List[int],
                              device: torch.device,
                              tol: float = 1e-6) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize active-block Hessian with mass-weight/TR in active space and return:
      freqs_cm : (nmode,)
      modes    : (nmode, 3N) mass-weighted eigenvectors embedded to full 3N (torch)
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_act = torch.as_tensor(coords_bohr.reshape(-1, 3)[active_idx, :], dtype=H_act.dtype, device=device)
        masses_act_au = torch.as_tensor([atomic_masses[int(z)] * AMU2AU
                                         for z in np.array(atomic_numbers, int)[active_idx]],
                                        dtype=H_act.dtype, device=device)
        Hmw = H_act.clone()
        _mw_tr_project_active_inplace(Hmw, coords_act, masses_act_au)
        omega2, Vsub = torch.linalg.eigh(Hmw, UPLO="U")
        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]
        Vsub = Vsub[:, sel]  # (3N_act, nsel)

        # Embed to full 3N (mass-weighted eigenvectors)
        modes_full = torch.zeros((Vsub.shape[1], 3 * N), dtype=Hmw.dtype, device=device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))  # give frozen list
        mask_t = torch.as_tensor(mask_dof, dtype=torch.bool, device=device)
        modes_full[:, mask_t] = Vsub.T
        # frequencies
        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()

        del coords_act, masses_act_au, Hmw, omega2, Vsub, mask_t
        if torch.cuda.is_available() and H_act.is_cuda:
            torch.cuda.empty_cache()
        return freqs_cm, modes_full


def _mode_direction_by_root_from_Hact(H_act: torch.Tensor,
                                      coords_bohr: np.ndarray,
                                      atomic_numbers: List[int],
                                      masses_au_t: torch.Tensor,
                                      active_idx: List[int],
                                      device: torch.device,
                                      root: int = 0) -> np.ndarray:
    """
    TS direction from the *active* Hessian block. Mass-weighting/TR are done in the
    active space. Result is embedded back to full 3N in Cartesian space.
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_act = torch.as_tensor(coords_bohr.reshape(-1, 3)[active_idx, :], dtype=H_act.dtype, device=device)
        masses_act_au = masses_au_t[active_idx].to(device=device, dtype=H_act.dtype)
        # mass-weight + TR in active space
        Hmw = H_act.clone()
        _mw_tr_project_active_inplace(Hmw, coords_act, masses_act_au)

        # eigenvector for requested root
        if int(root) == 0:
            try:
                w, V = torch.lobpcg(Hmw, k=1, largest=False)
                u_mw = V[:, 0]
            except Exception:
                vals, vecs = torch.linalg.eigh(Hmw, UPLO="U")
                u_mw = vecs[:, torch.argmin(vals)]
                del vals, vecs
        else:
            vals, vecs = torch.linalg.eigh(Hmw, UPLO="U")
            neg = (vals < 0)
            neg_inds = torch.nonzero(neg, as_tuple=False).view(-1)
            if neg_inds.numel() == 0:
                pick = int(torch.argmin(vals).item())
            else:
                k = max(0, min(int(root), neg_inds.numel() - 1))
                pick = int(neg_inds[k].item())
            u_mw = vecs[:, pick]
            del vals, vecs

        # Mass un-weight to Cartesian in the active space, then embed to full 3N
        masses_act_amu = (masses_act_au / AMU2AU).to(dtype=H_act.dtype, device=device)
        m3 = torch.repeat_interleave(masses_act_amu, 3)
        v_cart_act = u_mw / torch.sqrt(m3)
        v_cart_act = v_cart_act / torch.linalg.norm(v_cart_act)

        full = torch.zeros(3 * N, dtype=H_act.dtype, device=device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))
        mask_t = torch.as_tensor(mask_dof, dtype=torch.bool, device=device)
        full[mask_t] = v_cart_act
        mode = full.reshape(-1, 3).detach().cpu().numpy()

        del coords_act, masses_act_au, masses_act_amu, m3, v_cart_act, full, mask_t, Hmw, u_mw
        if torch.cuda.is_available() and H_act.is_cuda:
            torch.cuda.empty_cache()
        return mode


def _bofill_update_active(H_act: torch.Tensor,
                          delta_act: np.ndarray,
                          g_new_act: np.ndarray,
                          g_old_act: np.ndarray,
                          eps: float = 1e-12) -> torch.Tensor:
    """
    Memory-efficient Bofill update on the *active* Cartesian Hessian block.
    Apply symmetric rank-1/2 updates directly **in place** using only the **upper triangle**
    index set (and mirror to the lower) to avoid allocating large NxN temporaries.
    No explicit (H+H^T)/2 symmetrization step is performed.
    """
    device = H_act.device
    dtype = H_act.dtype

    # as torch vectors
    d = torch.as_tensor(delta_act, dtype=dtype, device=device).reshape(-1)
    g0 = torch.as_tensor(g_old_act, dtype=dtype, device=device).reshape(-1)
    g1 = torch.as_tensor(g_new_act, dtype=dtype, device=device).reshape(-1)
    y = g1 - g0

    # Use current symmetric H_act for matvec (no extra allocation)
    Hd = H_act @ d
    xi = y - Hd

    d_dot_xi = torch.dot(d, xi)
    d_norm2 = torch.dot(d, d)
    xi_norm2 = torch.dot(xi, xi)

    # guards
    denom_ms = d_dot_xi if torch.abs(d_dot_xi) > eps else torch.sign(d_dot_xi + 0.0) * eps
    denom_psb_d4 = d_norm2 * d_norm2 if d_norm2 > eps else eps
    denom_psb_d2 = d_norm2 if d_norm2 > eps else eps
    denom_phi = d_norm2 * xi_norm2 if (d_norm2 > eps and xi_norm2 > eps) else (1.0)

    phi = 1.0 - (d_dot_xi * d_dot_xi) / denom_phi
    phi = torch.clamp(phi, 0.0, 1.0)

    # coefficients for rank updates
    alpha = (1.0 - phi) / denom_ms                      # for xi xi^T
    beta  = - phi * (d_dot_xi / denom_psb_d4)           # for d d^T
    gamma = phi / denom_psb_d2                          # for d xi^T + xi d^T

    n = H_act.shape[0]
    iu0, iu1 = torch.triu_indices(n, n, device=device)
    is_diag = (iu0 == iu1)
    off = ~is_diag

    # Diagonal contributions (i == j): alpha*xi_i^2 + beta*d_i^2 + 2*gamma*d_i*xi_i
    if is_diag.any():
        idx = iu0[is_diag]
        H_act[idx, idx].add_(alpha * xi[idx] * xi[idx]
                             + beta * d[idx] * d[idx]
                             + 2.0 * gamma * d[idx] * xi[idx])

    # Off-diagonal (i < j): symmetric update
    if off.any():
        i = iu0[off]; j = iu1[off]
        inc = (alpha * xi[i] * xi[j]
               + beta * d[i] * d[j]
               + gamma * (d[i] * xi[j] + xi[i] * d[j]))
        H_act[i, j].add_(inc)
        H_act[j, i].add_(inc)

    return H_act


# ===================================================================
#                        HessianDimer (extended)
# ===================================================================

class HessianDimer:
    """
    Dimer-based TS search with periodic Hessian updates.

    Extensions in this implementation:
      - `root` parameter: choose which imaginary mode to follow (0 = most negative).
      - Pass-through kwargs: `dimer_kwargs` and `lbfgs_kwargs` to tune internals.
      - Hard cap on total LBFGS steps across segments: `max_total_cycles`.
      - PHVA (active DOF subspace) + TR projection for mode picking,
        respecting ``freeze_atoms``, with in-place operations. When ``root == 0`` the
        implementation prefers LOBPCG.
      - The flatten loop uses a *Bofill*-updated active Hessian block, so the
        expensive exact Hessian is evaluated only once before the flatten loop and
        once at the end for the final frequency analysis.
      - UMA calculator kwargs accept ``freeze_atoms`` and ``hessian_calc_mode`` and
        default to ``return_partial_hessian=True`` (active-block Hessian when frozen).
    """

    def __init__(self,
                 fn: str,
                 out_dir: str = "./result_dimer",
                 thresh_loose: str = "gau_loose",
                 thresh: str = "baker",
                 update_interval_hessian: int = 50,
                 neg_freq_thresh_cm: float = 5.0,
                 flatten_amp_ang: float = 0.20,
                 flatten_max_iter: int = 20,
                 mem: int = 100000,
                 use_lobpcg: bool = True,  # kept for backward compat (not used when root!=0)
                 uma_kwargs: Optional[dict] = None,
                 device: str = "auto",
                 dump: bool = False,
                 #
                 # New:
                 root: int = 0,
                 dimer_kwargs: Optional[Dict[str, Any]] = None,
                 lbfgs_kwargs: Optional[Dict[str, Any]] = None,
                 max_total_cycles: int = 10000,
                 #
                # Pass geom kwargs so freeze-links and YAML geometry overrides apply on the light path (fix #1)
                 geom_kwargs: Optional[Dict[str, Any]] = None,
                 ) -> None:

        self.fn = fn
        self.out_dir = Path(out_dir); self.out_dir.mkdir(parents=True, exist_ok=True)
        self.vib_dir = self.out_dir / "vib"; self.vib_dir.mkdir(parents=True, exist_ok=True)

        self.thresh_loose = thresh_loose
        self.thresh = thresh
        self.update_interval_hessian = int(update_interval_hessian)
        self.neg_freq_thresh_cm = float(neg_freq_thresh_cm)
        self.flatten_amp_ang = float(flatten_amp_ang)
        self.flatten_max_iter = int(flatten_max_iter)
        self.mem = int(mem)
        self.use_lobpcg = bool(use_lobpcg)  # used only when root==0 shortcut
        self.root = int(root)
        self.dimer_kwargs = dict(dimer_kwargs or {})
        self.lbfgs_kwargs = dict(lbfgs_kwargs or {})
        self.max_total_cycles = int(max_total_cycles)

        # Track total cycles globally across ALL loops/segments (fix #2)
        self._cycles_spent = 0

        # UMA settings
        self.uma_kwargs = dict(charge=0, spin=1, model="uma-s-1p1",
                               task_name="omol", device="auto") if uma_kwargs is None else dict(uma_kwargs)

        # Geometry & masses (use provided geom kwargs so freeze_atoms etc. apply)
        gkw = dict(geom_kwargs or {})
        coord_type = gkw.pop("coord_type", "cart")
        self.geom = geom_loader(fn, coord_type=coord_type, **gkw)
        self.masses_amu = np.array([atomic_masses[z] for z in self.geom.atomic_numbers])
        self.masses_au_t = torch.as_tensor(self.masses_amu * AMU2AU, dtype=torch.float32)

        # --- Preserve freeze list (for PHVA) ---
        self.freeze_atoms: List[int] = list(gkw.get("freeze_atoms", [])) if "freeze_atoms" in gkw else []

        # Device
        self.device = _torch_device(device)
        self.masses_au_t = self.masses_au_t.to(self.device)

        # temp file for Dimer orientation (N_raw)
        self.mode_path = self.out_dir / ".dimer_mode.dat"

        self.dump = bool(dump)
        self.optim_all_path = self.out_dir / "optimization_all.trj"

    # ----- One dimer segment for up to n_steps; returns (steps_done, converged) -----
    def _dimer_segment(self, threshold: str, n_steps: int) -> Tuple[int, bool]:
        # Dimer calculator using current mode as initial N
        calc_sp = uma_pysis(out_hess_torch=False, **self.uma_kwargs)

        # Merge user dimer kwargs (but enforce N_raw & write_orientations)
        dimer_kwargs = dict(self.dimer_kwargs)
        dimer_kwargs.update({
            "calculator": calc_sp,
            "N_raw": str(self.mode_path),
            "write_orientations": False,  # runner override to reduce IO
            "seed": 0,                    # runner override for determinism
            "mem": self.mem,              # accepted by Calculator base through **kwargs
        })
        dimer = Dimer(**dimer_kwargs)

        self.geom.set_calculator(dimer)

        # LBFGS kwargs: enforce thresh/max_cycles/out_dir/dump; allow others
        lbfgs_kwargs = dict(self.lbfgs_kwargs)
        lbfgs_kwargs.update({
            "max_cycles": n_steps,
            "thresh": threshold,
            "out_dir": str(self.out_dir),
            "dump": self.dump,
        })
        opt = LBFGS(self.geom, **lbfgs_kwargs)
        opt.run()
        steps = opt.cur_cycle
        converged = opt.is_converged
        self.geom.set_calculator(None)

        # Append to concatenated trajectory if dump enabled
        if self.dump:
            part_path = self.out_dir / "optimization.trj"
            if part_path.exists():
                with self.optim_all_path.open("a", encoding="utf-8") as f_all, \
                     part_path.open("r", encoding="utf-8") as f_part:
                    f_all.write(f_part.read())
        return steps, converged

    # ----- Loop dimer segments, updating mode from Hessian every interval -----
    def _dimer_loop(self, threshold: str) -> int:
        """
        Run multiple LBFGS segments separated by periodic Hessian-based mode updates.
        Consumes from a *global* cycle budget self.max_total_cycles (fix #2).
        """
        steps_in_this_call = 0
        while True:
            remaining_global = max(0, self.max_total_cycles - self._cycles_spent)
            if remaining_global == 0:
                break
            steps_this = min(self.update_interval_hessian, remaining_global)
            steps, ok = self._dimer_segment(threshold, steps_this)
            self._cycles_spent += steps
            steps_in_this_call += steps
            if ok:
                break
            # If budget exhausted after this segment, stop before doing a Hessian update
            if (self.max_total_cycles - self._cycles_spent) <= 0:
                break
            # Update mode from Hessian (respect freeze atoms via PHVA)
            H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
            N = len(self.geom.atomic_numbers)
            coords_bohr_t = torch.as_tensor(self.geom.coords.reshape(-1, 3),
                                            dtype=H_t.dtype, device=H_t.device)
            # full vs active-block Hessian
            if H_t.size(0) == 3 * N:
                mode_xyz = _mode_direction_by_root(
                    H_t, coords_bohr_t, self.masses_au_t,
                    root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
                )
            else:
                # partial (active) Hessian returned by UMA
                active_idx = _active_indices(N, self.freeze_atoms if len(self.freeze_atoms) > 0 else [])
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_t, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root
                )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del H_t, coords_bohr_t, mode_xyz
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        return steps_in_this_call

    # ----- Flatten (resolve multiple imaginary modes) -----
    def _flatten_once(self) -> bool:
        """
        Legacy: exact-Hessian-based flattening (kept for reference / fallback).
        """
        H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        freqs_cm, modes = _frequencies_cm_and_modes(
            H_t, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3), self.device,
            freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
        )
        del H_t
        neg_idx_all = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx_all) <= 1:
            del modes
            return False

        # Identify the "primary" imaginary by root among negative modes
        order = np.argsort(freqs_cm[neg_idx_all])  # ascending (more negative first)
        root_clamped = max(0, min(self.root, len(order) - 1))
        primary_idx = neg_idx_all[order[root_clamped]]

        targets = [i for i in neg_idx_all if i != primary_idx]
        if not targets:
            del modes
            return False

        # Reference structure and energy
        ref = self.geom.coords.reshape(-1, 3).copy()
        _ = _calc_energy(self.geom, self.uma_kwargs)  # E_ref (unused, but keeps semantics)

        # mass scaling so that carbon ~ amplitude
        mass_scale = np.sqrt(12.011 / self.masses_amu)[:, None]
        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        disp_total = np.zeros_like(ref)
        for idx in targets:
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)  # mass-weighted eigenvector embedded to 3N
            # Convert to Cartesian step direction already done downstream in writer,
            # but for flattening we only need a normalized direction in Cartesian:
            # use masses to unweight:
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)
            disp0 = amp_bohr * mass_scale * v_cart

            self.geom.coords = (ref +  disp0).reshape(-1)
            E_plus  = _calc_energy(self.geom, self.uma_kwargs)
            self.geom.coords = (ref -  disp0).reshape(-1)
            E_minus = _calc_energy(self.geom, self.uma_kwargs)
            self.geom.coords = ref.reshape(-1)

            disp_total += (disp0 if E_plus <= E_minus else -disp0)

        del modes
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        self.geom.coords = (ref + disp_total).reshape(-1)
        return True

    def _flatten_once_with_modes(self, freqs_cm: np.ndarray, modes: torch.Tensor) -> bool:
        """
        Flatten using precomputed (approximate) modes (mass-weighted, embedded).
        """
        neg_idx_all = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx_all) <= 1:
            return False

        order = np.argsort(freqs_cm[neg_idx_all])
        root_clamped = max(0, min(self.root, len(order) - 1))
        primary_idx = neg_idx_all[order[root_clamped]]
        targets = [i for i in neg_idx_all if i != primary_idx]
        if not targets:
            return False

        ref = self.geom.coords.reshape(-1, 3).copy()
        _ = _calc_energy(self.geom, self.uma_kwargs)

        mass_scale = np.sqrt(12.011 / self.masses_amu)[:, None]
        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        disp_total = np.zeros_like(ref)
        for idx in targets:
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)
            disp0 = amp_bohr * mass_scale * v_cart

            self.geom.coords = (ref + disp0).reshape(-1)
            E_plus = _calc_energy(self.geom, self.uma_kwargs)
            self.geom.coords = (ref - disp0).reshape(-1)
            E_minus = _calc_energy(self.geom, self.uma_kwargs)
            self.geom.coords = ref.reshape(-1)

            disp_total += (disp0 if E_plus <= E_minus else -disp0)

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        self.geom.coords = (ref + disp_total).reshape(-1)
        return True

    # ----- Run full procedure -----
    def run(self) -> None:
        if self.dump and self.optim_all_path.exists():
            self.optim_all_path.unlink()

        N = len(self.geom.atomic_numbers)
        active_idx = _active_indices(N, self.freeze_atoms if len(self.freeze_atoms) > 0 else [])
        mask_dof = _active_mask_dof(N, self.freeze_atoms if len(self.freeze_atoms) > 0 else [])

        # (1) Initial Hessian → pick direction by `root`
        H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        coords_bohr_t = torch.as_tensor(self.geom.coords.reshape(-1, 3),
                                        dtype=H_t.dtype, device=H_t.device)

        if H_t.size(0) == 3 * N:
            # Skip heavy TR-projection residual check to conserve VRAM.
            print("[CHECK] TR-projection residual check skipped to conserve VRAM.")
            mode_xyz = _mode_direction_by_root(
                H_t, coords_bohr_t, self.masses_au_t,
                root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
            )
        else:
            print("[CHECK] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
            mode_xyz = _mode_direction_by_root_from_Hact(
                H_t, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
            )
        np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
        del mode_xyz, coords_bohr_t, H_t
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # (2) Loose loop
        if self.root!=0:
            print("[WARNING] root != 0. Use this 'root' in first dimer loop")
            print(f"Dimer Loop with initial direction from mode {self.root}...")
            self.root==0
            self.thresh_loose = self.thresh
        else:
            print("Loose Dimer Loop...")

        self._dimer_loop(self.thresh_loose)

        # (3) Update mode & normal loop
        H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        coords_bohr_t = torch.as_tensor(self.geom.coords.reshape(-1, 3),
                                        dtype=H_t.dtype, device=H_t.device)
        if H_t.size(0) == 3 * N:
            print("[CHECK] TR-projection residual check skipped to conserve VRAM.")
            mode_xyz = _mode_direction_by_root(
                H_t, coords_bohr_t, self.masses_au_t,
                root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
            )
        else:
            print("[CHECK] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
            mode_xyz = _mode_direction_by_root_from_Hact(
                H_t, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
            )
        np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
        del mode_xyz, coords_bohr_t, H_t
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        print("Normal Dimer Loop...")
        self._dimer_loop(self.thresh)

        # (4) Flatten Loop — *reduced* exact Hessian calls via Bofill updates (active DOF only)
        print("Flatten Loop with Bofill-updated active Hessian...")

        # (4.1) Evaluate one exact Hessian at the loop start and prepare the active block
        H_any = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        if H_any.size(0) == 3 * N:
            # full → extract active
            H_act = _extract_active_block(H_any, mask_dof)  # torch (3N_act,3N_act)
        else:
            # UMA already returned active-block Hessian
            H_act = H_any
        del H_any
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # Gradient & coordinates snapshot for quasi-Newton updates
        x_prev = self.geom.coords.copy().reshape(-1)             # (3N,)
        g_prev = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)  # (3N,)

        # Flatten iterations with *approximate* Hessian updates
        for it in range(self.flatten_max_iter):
            if (self.max_total_cycles - self._cycles_spent) <= 0:
                break

            # (a) Estimate current imaginary modes using the *active* Hessian
            freqs_est = _frequencies_from_Hact(H_act, self.geom.atomic_numbers,
                                               self.geom.coords.reshape(-1, 3), active_idx, self.device)
            n_imag = int(np.sum(freqs_est < -abs(self.neg_freq_thresh_cm)))
            print(f"[IMAG~] n≈{n_imag}  (approx imag: {[float(x) for x in freqs_est if x < -abs(self.neg_freq_thresh_cm)]})")
            if n_imag <= 1:
                break

            # (b) Get approximate modes for flattening (embedded, mass-weighted)
            freqs_cm_approx, modes_embedded = _modes_from_Hact_embedded(
                H_act, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3), active_idx, self.device
            )

            # (c) Do flatten step using the approximate modes
            x_before_flat = self.geom.coords.copy().reshape(-1)
            did_flatten = self._flatten_once_with_modes(freqs_cm_approx, modes_embedded)
            if not did_flatten:
                break
            x_after_flat = self.geom.coords.copy().reshape(-1)

            # (d) Bofill update using UMA gradients across the flatten displacement
            g_after_flat = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)
            delta_flat_full = x_after_flat - x_before_flat
            delta_flat_act = delta_flat_full[mask_dof]
            g_old_act = g_prev[mask_dof]
            g_new_act = g_after_flat[mask_dof]
            H_act = _bofill_update_active(H_act, delta_flat_act, g_new_act, g_old_act)

            # (e) Refresh dimer direction from updated active Hessian
            mode_xyz = _mode_direction_by_root_from_Hact(
                H_act, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
            )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")

            # (f) Re-optimize with Dimer (consumes global cycle budget)
            self._dimer_loop(self.thresh)

            # (g) Bofill update again across the optimization displacement
            x_after_opt = self.geom.coords.copy().reshape(-1)
            g_after_opt = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)
            delta_opt_full = x_after_opt - x_after_flat
            delta_opt_act = delta_opt_full[mask_dof]
            g_old_act2 = g_after_flat[mask_dof]
            g_new_act2 = g_after_opt[mask_dof]
            H_act = _bofill_update_active(H_act, delta_opt_act, g_new_act2, g_old_act2)

            # (h) Prepare for next iteration
            x_prev = x_after_opt
            g_prev = g_after_opt

        # (5) Final outputs
        final_xyz = self.out_dir / "final_geometry.xyz"
        atoms_final = Atoms(self.geom.atoms, positions=(self.geom.coords.reshape(-1, 3) * BOHR2ANG), pbc=False)
        write(final_xyz, atoms_final)

        # Final Hessian → imaginary mode animation
        H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        if H_t.size(0) == 3 * N:
            freqs_cm, modes = _frequencies_cm_and_modes(
                H_t, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3), self.device,
                freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
            )
        else:
            freqs_cm, modes = _modes_from_Hact_embedded(
                H_t, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3),
                active_idx, self.device
            )

        del H_t
        neg_idx = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx) == 0:
            print("[INFO] No imaginary mode found at the end (ν_min = %.2f cm^-1)." % (freqs_cm.min(),))
            del modes
        else:
            # primary (by root)
            order = np.argsort(freqs_cm[neg_idx])
            primary_idx = neg_idx[order[max(0, min(self.root, len(order)-1))]]
            # convert mass-weighted mode to Cartesian here for writing
            mode_mw = modes[primary_idx]  # (3N,)
            masses_amu_t = torch.as_tensor(self.masses_amu, dtype=mode_mw.dtype, device=mode_mw.device)
            m3 = torch.repeat_interleave(masses_amu_t, 3)
            v_cart = (mode_mw / torch.sqrt(m3)).detach().cpu().numpy()
            v_cart = v_cart / np.linalg.norm(v_cart)
            del modes, masses_amu_t, m3
            out_trj = self.vib_dir / f"final_imag_mode_{freqs_cm[primary_idx]:+.2f}cm-1.trj"
            out_pdb = self.vib_dir / f"final_imag_mode_{freqs_cm[primary_idx]:+.2f}cm-1.pdb"
            _write_mode_trj_and_pdb(self.geom, v_cart, out_trj, out_pdb,
                                    amplitude_ang=0.8, n_frames=20,
                                    comment=f"imag {freqs_cm[primary_idx]:+.2f} cm^-1")

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        print(f"[DONE] Saved final geometry → {final_xyz}")
        print(f"[DONE] Mode files → {self.vib_dir}")


# ===================================================================
#                         Defaults for CLI
# ===================================================================

# Geometry defaults
GEOM_KW = {
    "coord_type": "cart",     # str, coordinate representation for geom_loader ("cart" recommended here)
    "freeze_atoms": [],       # List[int], 0-based indices to freeze
}

# UMA calculator defaults
CALC_KW = {
    "charge": 0,              # int, total charge
    "spin": 1,                # int, multiplicity (2S+1)
    "model": "uma-s-1p1",     # str, UMA pretrained model ID
    "task_name": "omol",      # str, dataset/task tag carried into UMA's AtomicData
    "device": "auto",         # "cuda" | "cpu" | "auto"
    "max_neigh": None,        # Optional[int], override model's neighbor cap
    "radius": None,           # Optional[float], cutoff radius (Å)
    "r_edges": False,         # bool, store edge vectors in graph (UMA option)
    "out_hess_torch": True,   # bool, RSIRFO sets this to True for torch Hessian; HessianDimer sets per-call
    # --- Hessian interfaces to UMA ---
    "hessian_calc_mode": "Analytical",        # "Analytical" | "FiniteDifference"
    "return_partial_hessian": True,           # Default: receive the active-DOF Hessian block from UMA
    # freeze_atoms is injected below from the geometry CLI/YAML configuration
}

# Optimizer base (common) — used by both RSIRFO and the inner LBFGS of HessianDimer
OPT_BASE_KW = {
    "thresh": "gau",          # str, convergence preset: "gau_loose"|"gau"|"gau_tight"|"gau_vtight"|"baker"|"never"
    "max_cycles": 10000,      # int, hard cap
    "print_every": 1,         # int, logging cadence
    "min_step_norm": 1e-8,    # float, if ||step|| < this and assert_min_step, raise
    "assert_min_step": True,  # bool, enforce min_step_norm check

    # Convergence criteria toggles
    "rms_force": None,        # Optional[float], override RMS(F) threshold (None → preset)
    "rms_force_only": False,  # bool
    "max_force_only": False,  # bool
    "force_only": False,      # bool

    # Extra mechanisms (not used by inner LBFGS unless supported)
    "converge_to_geom_rms_thresh": 0.05,
    "overachieve_factor": 0.0,
    "check_eigval_structure": False,

    # Dumping / restart
    "dump": False,
    "dump_restart": False,
    "prefix": "",
    "out_dir": "./result_ts_opt/",
}

DIMER_KW = {
    # --- Geometry / length ---
    "length": 0.0189,                 # float, dimer half-length in Bohr (≈ 0.01 Å)

    # --- Rotation loop settings ---
    "rotation_max_cycles": 15,        # int, max rotation cycles per update
    "rotation_method": "fourier",     # "fourier" (robust) | "direct" (steepest-like)
    "rotation_thresh": 1e-4,          # float, threshold on ||rot_force||
    "rotation_tol": 1,                # float (deg), skip rotation if |angle| < tol
    "rotation_max_element": 0.001,    # float, max element for direct rotation step (Bohr)
    "rotation_interpolate": True,     # bool, interpolate f1 during rotation (Fourier method)
    "rotation_disable": False,        # bool, disable rotation (use given N as-is)
    "rotation_disable_pos_curv": True,# bool, if curvature positive after rotation, restore previous N
    "rotation_remove_trans": True,    # bool, remove net translation from N

    # --- Translational part of projected force ---
    "trans_force_f_perp": True,       # bool, include perpendicular force into translational component when C<0

    # --- Initial guess helpers (mutually exclusive) ---
    "bonds": None,                    # Optional[List[Tuple[int,int]]], use weighted-bond mode as N_raw
    "N_hessian": None,                # Optional[str], path to HDF5 Hessian; use 1st imag. mode as N_raw

    # --- Optional bias potentials for robustness ---
    "bias_rotation": False,           # bool, quadratic bias along initial N during rotation
    "bias_translation": False,        # bool, add Gaussians along the path to stabilize translation
    "bias_gaussian_dot": 0.1,         # float, target dot product when tuning Gaussian height

    # --- Stochastic control / IO ---
    "seed": None,                     # Optional[int], RNG seed; Runner sets to 0 for determinism
    "write_orientations": True,       # bool, write N.trj each call; Runner overrides to False to reduce IO

    # --- Hessian forwarding ---
    "forward_hessian": True,          # bool, allow forwarding get_hessian to wrapped calculator
}

# ---------- Reference: internal L-BFGS defaults for TS optimization highlighting deviations from OPT_BASE_KW ----------
LBFGS_TS_KW = {
    **OPT_BASE_KW,

    # History / memory
    "keep_last": 7,                   # int, number of (s, y) pairs

    # Preconditioner / initial scaling
    "beta": 1.0,                      # float, β in -(H + βI)^{-1} g
    "gamma_mult": False,              # bool, estimate β from previous cycle (Nocedal 7.20)

    # Step-size control
    "max_step": 0.30,                 # float, max component magnitude of step (Bohr, working coords)
    "control_step": True,             # bool, scale step to satisfy max_component ≤ max_step

    # Safeguards & line search
    "double_damp": True,              # bool, ensure s·y > 0
    "line_search": False,             # bool, polynomial line search on last step (off by default)

    # Regularized L-BFGS (μ_reg)
    "mu_reg": None,                   # Optional[float], initial regularization; if set, disables double_damp etc.
    "max_mu_reg_adaptions": 10,       # int, max trial steps for μ adaptation
}

# HessianDimer defaults (CLI-level)
hessian_dimer_KW = {
    "thresh_loose": "gau_loose",      # str, loose threshold preset for first pass
    "thresh": "gau",                  # str, main threshold preset for TS search
    "update_interval_hessian": 1000,  # int, LBFGS cycles per Hessian refresh for direction
    "neg_freq_thresh_cm": 5.0,        # float, treat ν < -this as imaginary (cm^-1)
    "flatten_amp_ang": 0.10,          # float, mass-scaled displacement amplitude for flattening (Å)
    "flatten_max_iter": 20,           # int, max flattening iterations
    "mem": 100000,                    # int, scratch/IO memory passed through Calculator (**kwargs)
    "use_lobpcg": True,               # kept for compat (root==0 uses LOBPCG regardless; fallback if fails)
    "device": "auto",                 # str, "cuda"|"cpu"|"auto" for torch-side ops
    "root": 0,                        # int, 0: follow the most negative mode
    "dimer": {**DIMER_KW},            # Dict, default kwargs forwarded to Dimer (Runner may override some)
    "lbfgs": {**LBFGS_TS_KW},         # Dict, default kwargs forwarded to inner LBFGS
}

# RSIRFO (TS Hessian optimizer) defaults (subset; additional keys may be provided)
RSIRFO_KW = {
    "roots": [0],                # which mode(s) to follow uphill
    "hessian_ref": None,
    "rx_modes": None,
    "prim_coord": None,
    "rx_coords": None,
    "hessian_init": "calc",
    "hessian_update": "bofill",
    "hessian_recalc": 100,       # Optional[int], recalc exact Hessian every N cycles
    "hessian_recalc_adapt": 2.0, # If norm(force) become 1/hessian_recalc_adapt, recalc hessian
    "hessian_recalc_reset": True,
    "max_micro_cycles": 50,
    "trust_radius": 0.30,
    "trust_max": 0.30,
    "augment_bonds": False,
    "min_line_search": False,
    "max_line_search": False,
    "assert_neg_eigval": False,
}


# ===================================================================
#                            CLI
# ===================================================================

@click.command(
    help="Transition state optimization (Dimer or RS-I-RFO).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure (.pdb, .xyz, .trj, ...)",
)
@click.option("-q", "--charge", type=int, required=True, help="Total charge")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1)")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze parent atoms of link hydrogens (PDB only).")
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Max cycles / steps cap")
@click.option("--opt-mode", type=str, default="light", show_default=True,
              help="light (=Dimer) or heavy (=RSIRFO)")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump optimization trajectory")
@click.option("--out-dir", type=str, default="./result_ts_opt/", show_default=True, help="Output directory")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, opt, dimer, rsirfo).",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None,
    help="Choose UMA Hessian evaluation mode (overrides YAML/calc.hessian_calc_mode).",
)
def cli(
    input_path: Path,
    charge: int,
    spin: int,
    freeze_links: bool,
    max_cycles: int,
    opt_mode: str,
    dump: bool,
    out_dir: str,
    args_yaml: Optional[Path],
    hessian_calc_mode: Optional[str],
) -> None:
    time_start = time.perf_counter()

    # --------------------------
    # 1) Assemble configuration
    # --------------------------
    yaml_cfg = load_yaml_dict(args_yaml)
    geom_cfg = dict(GEOM_KW)
    calc_cfg = dict(CALC_KW)
    opt_cfg  = dict(OPT_BASE_KW)
    simple_cfg = dict(hessian_dimer_KW)
    rsirfo_cfg = dict(RSIRFO_KW)

    # YAML → defaults
    apply_yaml_overrides(
        yaml_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (opt_cfg, (("opt",),)),
            (simple_cfg, (("hessian_dimer",),)),
            (rsirfo_cfg, (("rsirfo",),)),
        ],
    )

    # CLI overrides
    calc_cfg["charge"] = int(charge)
    calc_cfg["spin"]   = int(spin)
    opt_cfg["max_cycles"] = int(max_cycles)
    opt_cfg["dump"]       = bool(dump)
    opt_cfg["out_dir"]    = out_dir

    # Hessian mode override from CLI
    if hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

    # Freeze links (PDB only): merge with existing list
    if freeze_links and input_path.suffix.lower() == ".pdb":
        try:
            detected = detect_freeze_links(input_path)
        except Exception as e:
            click.echo(f"[freeze-links] WARNING: Could not detect link parents: {e}", err=True)
            detected = []
        merged = merge_freeze_atom_indices(geom_cfg, detected)
        if merged:
            click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, merged))}")

    # ---- Paess freeze_atoms from geom → calc (so UMA knows active DOF for FD Hessian) ----
    if "freeze_atoms" in geom_cfg:
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))

    kind = _norm_opt_mode(opt_mode)
    out_dir_path = Path(opt_cfg["out_dir"]).resolve()

    # Pretty-print config summary
    click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
    click.echo(pretty_block("calc", calc_cfg))
    click.echo(pretty_block("opt",  {**opt_cfg, "out_dir": str(out_dir_path)}))
    if kind == "light":
        # Split out pass-through dicts for readability
        sd_cfg_for_echo = dict(simple_cfg)
        sd_cfg_for_echo["dimer"] = dict(simple_cfg.get("dimer", {}))
        sd_cfg_for_echo["lbfgs"] = dict(simple_cfg.get("lbfgs", {}))
        click.echo(pretty_block("hessian_dimer", sd_cfg_for_echo))
    else:
        click.echo(pretty_block("rsirfo", rsirfo_cfg))

    # --------------------------
    # 2) Prepare geometry dir
    # --------------------------
    out_dir_path.mkdir(parents=True, exist_ok=True)

    # --------------------------
    # 3) Run
    # --------------------------
    try:
        if kind == "light":
            # HessianDimer runner construction
            uma_kwargs_for_sd = dict(calc_cfg)
            # device comes from calc_cfg; out_hess_torch/return_partial_hessian handled in wrapper
            runner = HessianDimer(
                fn=str(input_path),
                out_dir=str(out_dir_path),
                thresh_loose=simple_cfg.get("thresh_loose", "gau_loose"),
                thresh=simple_cfg.get("thresh", "gau"),
                update_interval_hessian=int(simple_cfg.get("update_interval_hessian", 200)),
                neg_freq_thresh_cm=float(simple_cfg.get("neg_freq_thresh_cm", 5.0)),
                flatten_amp_ang=float(simple_cfg.get("flatten_amp_ang", 0.10)),
                flatten_max_iter=int(simple_cfg.get("flatten_max_iter", 20)),
                mem=int(simple_cfg.get("mem", 100000)),
                use_lobpcg=bool(simple_cfg.get("use_lobpcg", True)),
                uma_kwargs=uma_kwargs_for_sd,
                device=str(simple_cfg.get("device", calc_cfg.get("device", "auto"))),
                dump=bool(opt_cfg["dump"]),
                root=int(simple_cfg.get("root", 0)),
                dimer_kwargs=dict(simple_cfg.get("dimer", {})),
                lbfgs_kwargs=dict(simple_cfg.get("lbfgs", {})),
                max_total_cycles=int(opt_cfg["max_cycles"]),
                # Propagate geometry settings (freeze_atoms, coord_type, ...) to the HessianDimer runner (fix #1)
                geom_kwargs=dict(geom_cfg),
            )

            click.echo("\n=== TS optimization (HessianDimer) started ===\n")
            runner.run()
            click.echo("\n=== TS optimization (HessianDimer) finished ===\n")

            # Convert outputs if input is PDB
            if input_path.suffix.lower() == ".pdb":
                ref_pdb = input_path.resolve()
                # final_geometry.xyz → final_geometry.pdb
                final_xyz = out_dir_path / "final_geometry.xyz"
                final_pdb = out_dir_path / "final_geometry.pdb"
                try:
                    convert_xyz_to_pdb(final_xyz, ref_pdb, final_pdb)
                    click.echo(f"[convert] Wrote '{final_pdb}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert final geometry to PDB: {e}", err=True)
                # optimization_all.trj → optimization_all.pdb
                all_trj = out_dir_path / "optimization_all.trj"
                if all_trj.exists():
                    try:
                        opt_pdb = out_dir_path / "optimization_all.pdb"
                        convert_xyz_to_pdb(all_trj, ref_pdb, opt_pdb)
                        click.echo(f"[convert] Wrote '{opt_pdb}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)

        else:
            # RS-I-RFO (heavy)
            #  - Build geometry now and attach UMA calculator
            coord_type = geom_cfg.get("coord_type", "cart")
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)
            geometry = geom_loader(input_path, coord_type=coord_type, **coord_kwargs)

            calc_builder_or_instance = uma_pysis(**calc_cfg)  # includes freeze_atoms / hessian_calc_mode / return_partial_hessian
            try:
                geometry.set_calculator(calc_builder_or_instance())
            except TypeError:
                geometry.set_calculator(calc_builder_or_instance)

            # Construct RSIRFO optimizer
            rs_args = dict(rsirfo_cfg)
            opt_base = dict(opt_cfg)
            opt_base["out_dir"] = str(out_dir_path)

            # Normalize roots/root
            roots = rs_args.get("roots", None)
            root_single = rs_args.get("root", None)
            if roots is None and root_single is not None:
                roots = [int(root_single)]
            if roots is None:
                roots = [0]
            rs_args["roots"] = [int(x) for x in roots]
            if "root" in rs_args:
                rs_args.pop("root")

            optimizer = RSIRFOptimizer(geometry, **rs_args, **opt_base)

            click.echo("\n=== TS optimization (RS-I-RFO) started ===\n")
            optimizer.run()
            click.echo("\n=== TS optimization (RS-I-RFO) finished ===\n")

            # Convert outputs if input is PDB
            if input_path.suffix.lower() == ".pdb":
                ref_pdb = input_path.resolve()
                # final_geometry.xyz (from optimizer.final_fn) → final_geometry.pdb
                final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
                final_pdb = out_dir_path / "final_geometry.pdb"
                try:
                    convert_xyz_to_pdb(final_xyz_path, ref_pdb, final_pdb)
                    click.echo(f"[convert] Wrote '{final_pdb}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert final geometry to PDB: {e}", err=True)
                # optimization.trj → optimization.pdb
                trj_path = out_dir_path / "optimization.trj"
                if trj_path.exists():
                    try:
                        opt_pdb = out_dir_path / "optimization.pdb"
                        convert_xyz_to_pdb(trj_path, ref_pdb, opt_pdb)
                        click.echo(f"[convert] Wrote '{opt_pdb}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)

            # --- RSIRFO: write final imaginary mode like HessianDimer (PHVA/in-place or active) ---
            geometry.set_calculator(None)
            uma_kwargs_for_heavy = dict(calc_cfg)
            uma_kwargs_for_heavy["out_hess_torch"] = True
            device = _torch_device(calc_cfg.get("device", "auto"))

            H_t = _calc_full_hessian_torch(geometry, uma_kwargs_for_heavy, device)
            N = len(geometry.atomic_numbers)
            if H_t.size(0) == 3 * N:
                freqs_cm, modes = _frequencies_cm_and_modes(
                    H_t, geometry.atomic_numbers, geometry.coords.reshape(-1, 3), device,
                    freeze_idx=list(geom_cfg.get("freeze_atoms", [])) if len(geom_cfg.get("freeze_atoms", [])) > 0 else None
                )
            else:
                act_idx = _active_indices(N, list(geom_cfg.get("freeze_atoms", [])))
                freqs_cm, modes = _modes_from_Hact_embedded(
                    H_t, geometry.atomic_numbers, geometry.coords.reshape(-1, 3), act_idx, device
                )

            neg_idx = np.where(freqs_cm < -5.0)[0]  # use same threshold default as HessianDimer
            if len(neg_idx) == 0:
                click.echo("[INFO] No imaginary mode found at the end for RSIRFO.", err=True)
            else:
                roots = rs_args.get("roots", [0])
                main_root = int(roots[0]) if roots else 0
                order = np.argsort(freqs_cm[neg_idx])
                pick_idx = neg_idx[order[max(0, min(main_root, len(order)-1))]]
                mode_mw = modes[pick_idx]
                masses_amu_t = torch.as_tensor([atomic_masses[z] for z in geometry.atomic_numbers],
                                               dtype=mode_mw.dtype, device=mode_mw.device)
                m3 = torch.repeat_interleave(masses_amu_t, 3)
                v_cart = (mode_mw / torch.sqrt(m3)).detach().cpu().numpy()
                v_cart = v_cart / np.linalg.norm(v_cart)

                vib_dir = out_dir_path / "vib"
                vib_dir.mkdir(parents=True, exist_ok=True)
                out_trj = vib_dir / f"final_imag_mode_{freqs_cm[pick_idx]:+.2f}cm-1.trj"
                out_pdb = vib_dir / f"final_imag_mode_{freqs_cm[pick_idx]:+.2f}cm-1.pdb"

                @dataclass
                class _GProxy:
                    atoms: List[str]
                    coords: np.ndarray
                gproxy = _GProxy(atoms=geometry.atoms, coords=geometry.coords.copy())
                _write_mode_trj_and_pdb(gproxy, v_cart, out_trj, out_pdb,
                                        amplitude_ang=0.8, n_frames=20,
                                        comment=f"imag {freqs_cm[pick_idx]:+.2f} cm-1")

            if torch.cuda.is_available():
                torch.cuda.empty_cache()

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for TS Opt: {hh:02d}:{mm:02d}:{ss:06.3f}")

    except ZeroStepLength:
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        import traceback
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


# Allow `python -m pdb2reaction.ts_opt` direct execution
if __name__ == "__main__":
    cli()
