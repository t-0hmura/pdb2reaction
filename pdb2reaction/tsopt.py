# pdb2reaction/tsopt.py

"""
tsopt — Transition-state optimization CLI
====================================================================

Usage (CLI)
-----------
    pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q <charge>] [-m <multiplicity>] \
        [--opt-mode {light|heavy}] \
        [--freeze-links {True|False}] [--max-cycles <int>] [--dump {True|False}] \
        [--out-dir <dir>] [--args-yaml <file>] \
        [--hessian-calc-mode {Analytical|FiniteDifference}]

Examples
--------
    # Minimal (recommended: always specify charge and spin)
    pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode light \
        --out-dir ./result_tsopt/

    # Light mode (Hessian Dimer) with YAML overrides and finite-difference Hessian
    pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
        --freeze-links True --opt-mode light --max-cycles 10000 --dump False \
        --out-dir ./result_tsopt/ --args-yaml ./args.yaml \
        --hessian-calc-mode FiniteDifference

    # Heavy mode (RS-I-RFO) using YAML
    pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode heavy \
        --args-yaml ./args.yaml --out-dir ./result_tsopt/


Description
-----------
Transition-state optimization with two modes:

- **light**: Hessian Dimer TS search with periodic Hessian-based direction refresh and a
  memory-efficient flatten loop to eliminate excess imaginary modes. After the initial dimer
  stage, one exact Hessian is evaluated and its active (PHVA) block is kept and updated by
  Bofill (SR1/MS and PSB blend) between geometry updates in the flatten loop. Each
  flatten iteration:
    * estimate approximate imaginary modes using the current active Hessian (mass-weighted,
      TR-projected),
    * select only **spatially separated** extra imaginary modes (using representative atoms
      and a distance cutoff) and perform a mass-scaled flatten step in those modes,
    * apply a Bofill update **only for the flatten displacement** (no update for the dimer step),
    * refresh the dimer direction from the updated active Hessian,
    * run a dimer–LBFGS segment (consuming the global cycle budget),
    * recompute an exact Hessian at the end of the dimer segment for the next iteration.
  Once only one imaginary mode remains, a final exact Hessian is computed for frequency analysis.
  *If `root != 0`, the specified root is used only to seed the initial dimer direction; the
  algorithm then follows the most negative mode (`root = 0`) on subsequent updates.*

- **heavy**: RS-I-RFO Hessian-based TS optimizer.

The CLI `--opt-mode` accepts two modes:
`light` maps to the Hessian Dimer workflow, and `heavy` selects RS-I-RFO.
The default is `light`.

Configuration is driven by YAML overrides for sections: `geom`, `calc`, `opt`, `hessian_dimer`,
and `rsirfo`. The `hessian_dimer` section accepts nested `dimer` and `lbfgs` dictionaries
forwarded to the respective pysisyphus components. The optional key `use_lobpcg` in
`hessian_dimer` is **deprecated and ignored**; the implementation always tries LOBPCG for the
lowest eigenpair when `root == 0`, falling back to `torch.linalg.eigh` as needed.

Structures are loaded via `pysisyphus.helpers.geom_loader` (PDB/XYZ/TRJ/etc.). The UMA
calculator (`pdb2reaction.uma_pysis`) provides energies, gradients, and Hessians. UMA may return
an active-subspace (partial) Hessian block when `freeze_atoms` are present; finite-difference
Hessians honor the active subspace. `--hessian-calc-mode` selects Analytical or FiniteDifference.

For PDB inputs, optimization trajectories and the final geometry are also converted to PDB.
The final imaginary mode is written as `.pdb` only when the input is PDB (the `.trj` is always written).

Key behaviors and algorithmic notes
-----------------------------------
- **Direction selection**: choose which imaginary mode to follow using `root`
  (0 = most negative). For `root == 0`, the code prefers `torch.lobpcg` for the lowest eigenpair
  and falls back to `torch.linalg.eigh` when necessary. In *light* mode with `root != 0`,
  the initial dimer direction uses that root once and subsequent updates default to `root = 0`.

- **PHVA and TR projection**: an active-degree-of-freedom (PHVA) subspace with translation/
  rotation (TR) projection reduces GPU memory use. This respects `freeze_atoms`. A heavy
  clone-based TR self-check is disabled to conserve VRAM.

- **Flatten loop (light mode)**: one exact active-subspace Hessian is evaluated at the start of
  the flatten loop. Each iteration:
    * estimate approximate imaginary modes using the current active Hessian,
    * select only **spatially separated** extra imaginary modes (using representative atoms
      and a distance cutoff) and perform a mass-scaled flatten step in those modes,
    * apply a Bofill update **only for the flatten displacement** (no update for the dimer step),
    * refresh the dimer direction from the updated active Hessian,
    * run a dimer–LBFGS segment,
    * recompute an exact Hessian at the end of the dimer segment for the next iteration.
  Continue until only one imaginary mode remains, then compute a final exact Hessian for
  frequency analysis.

- **UMA integration**: `freeze_atoms` are propagated to UMA. Finite-difference Hessians honor
  the active subspace. UMA defaults to returning a partial (active) Hessian when applicable.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_tsopt/)
  ├─ final_geometry.xyz              # Optimized TS in XYZ format
  ├─ final_geometry.pdb              # Written when the input was PDB
  ├─ optimization_all.trj            # Light-mode trajectory when --dump True
  ├─ optimization_all.pdb            # PDB conversion of the above (PDB input)
  ├─ optimization.trj                # Heavy-mode (RS-I-RFO) trajectory
  ├─ optimization.pdb                # PDB conversion of the heavy-mode trajectory (PDB input)
  ├─ vib/
  │   ├─ final_imag_mode_±XXXX.Xcm-1.trj  # Final imaginary mode animation (.trj)
  │   └─ final_imag_mode_±XXXX.Xcm-1.pdb  # PDB animation (only when the input is PDB)
  └─ .dimer_mode.dat                 # Internal dimer orientation seed (light mode)

Notes
-----
- **Charge/spin**: `-q/--charge` and `-m/--mult` inherit `.gjf` template values when the input
  is `.gjf`; otherwise `-q/--charge` is required and the CLI aborts if omitted (multiplicity
  still defaults to 1). Override explicitly to avoid unphysical conditions.

- `--opt-mode light` runs Hessian Dimer with periodic Hessian-based direction refresh;
  `--opt-mode heavy` runs RS-I-RFO.

- `--freeze-links` is PDB-only and freezes parents of link hydrogens; these indices are merged
  into `geom.freeze_atoms` and also propagated to `calc.freeze_atoms` for UMA.

- Convergence and stepping behavior are configurable via YAML in
  `hessian_dimer.lbfgs`, `hessian_dimer.dimer`, `opt`, and `rsirfo` sections
  (e.g., thresholds, trust radii, memory).

- Imaginary-mode detection uses a default threshold of ~5 cm⁻¹; the primary mode written
  at the end is chosen via `root`.
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
from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .opt import (
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
)
from .utils import (
    convert_xyz_to_pdb,
    detect_freeze_links,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    normalize_choice,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
)
from .freq import (
    _torch_device,
    _build_tr_basis,
    _tr_orthonormal_basis,
    _mass_weighted_hessian,
    _calc_full_hessian_torch,
    _calc_energy,
    _write_mode_trj_and_pdb,
)


# Normalization helper
_OPT_MODE_ALIASES = (
    (("light",), "light"),
    (("heavy",), "heavy"),
)


# ===================================================================
#               Mass-weighted projection & vib analysis
# ===================================================================

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
    Lightweight TR-projection diagnostic. Returns zeros for residuals (clone-based
    checks disabled to conserve VRAM) and the rank estimate of the TR basis used.
    """
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
    For root==0 the routine prefers LOBPCG and falls back to eigh as needed.
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


def _calc_gradient(geom, uma_kwargs: dict) -> np.ndarray:
    """
    Return true Cartesian gradient (shape 3N,) in Hartree/Bohr.
    """
    calc = uma_pysis(**uma_kwargs)
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
        modes_full = torch.zeros((Vsub.shape[1], 3 * N), dtype=Hmw.dtype, device=Hmw.device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))
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
    Apply symmetric rank-1/2 updates directly in place using only the upper triangle
    index set (and mirror to the lower) to avoid allocating large temporaries.
    No explicit (H+H^T)/2 symmetrization step is performed.
    """
    device = H_act.device
    dtype = H_act.dtype

    # as torch vectors
    d = torch.as_tensor(delta_act, dtype=dtype, device=device).reshape(-1)
    g0 = torch.as_tensor(g_old_act, dtype=dtype, device=device).reshape(-1)
    g1 = torch.as_tensor(g_new_act, dtype=dtype, device=device).reshape(-1)
    y = g1 - g0

    # Use current symmetric H_act for matvec
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
    beta = -phi * (d_dot_xi / denom_psb_d4)             # for d d^T
    gamma = phi / denom_psb_d2                          # for d xi^T + xi d^T

    n = H_act.shape[0]
    iu0, iu1 = torch.triu_indices(n, n, device=device)
    is_diag = (iu0 == iu1)
    off = ~is_diag

    # Diagonal contributions (i == j): alpha*xi_i^2 + beta*d_i^2 + 2*gamma*d_i*xi_i
    if is_diag.any():
        idx = iu0[is_diag]
        diag_inc = (alpha * xi[idx] * xi[idx]
                    + beta * d[idx] * d[idx]
                    + 2.0 * gamma * d[idx] * xi[idx])
        # NOTE: use assignment so that advanced indexing actually updates H_act
        H_act[idx, idx] = H_act[idx, idx] + diag_inc

    # Off-diagonal (i < j): symmetric update
    if off.any():
        i = iu0[off]
        j = iu1[off]
        inc = (alpha * xi[i] * xi[j]
               + beta * d[i] * d[j]
               + gamma * (d[i] * xi[j] + xi[i] * d[j]))
        H_act[i, j] = H_act[i, j] + inc
        H_act[j, i] = H_act[j, i] + inc

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
        respecting `freeze_atoms`. For `root == 0` the implementation prefers LOBPCG.
      - The flatten loop uses a Bofill-updated active Hessian block in the active DOF
        subspace, but Bofill is applied *only* for the flatten displacements; after
        each dimer segment in the flatten loop, a fresh exact Hessian is recomputed.
      - Only **spatially separated** extra imaginary modes (based on representative
        atoms and a distance cutoff) are used for flattening to avoid overly large
        displacements when clustered imaginary modes are present.
      - UMA calculator kwargs accept `freeze_atoms` and `hessian_calc_mode` and
        default to returning a partial (active) Hessian when applicable.

    Note: the `use_lobpcg` argument is deprecated and ignored.
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
                 use_lobpcg: bool = True,  # deprecated; ignored
                 uma_kwargs: Optional[dict] = None,
                 device: str = "auto",
                 dump: bool = False,
                 #
                 root: int = 0,
                 dimer_kwargs: Optional[Dict[str, Any]] = None,
                 lbfgs_kwargs: Optional[Dict[str, Any]] = None,
                 max_total_cycles: int = 10000,
                 #
                 # Multi-mode flatten control
                 flatten_sep_cutoff: float = 2.0,
                 flatten_k: int = 10,
                 #
                 # Propagate geometry kwargs so freeze-links and YAML geometry overrides
                 # also apply in light mode.
                 geom_kwargs: Optional[Dict[str, Any]] = None,
                 ) -> None:

        self.fn = fn
        self.out_dir = Path(out_dir); self.out_dir.mkdir(parents=True, exist_ok=True)
        self.vib_dir = self.out_dir / "vib"; self.vib_dir.mkdir(parents=True, exist_ok=True)
        self.ref_pdb: Optional[Path] = Path(fn) if Path(fn).suffix.lower() == ".pdb" else None

        self.thresh_loose = thresh_loose
        self.thresh = thresh
        self.update_interval_hessian = int(update_interval_hessian)
        self.neg_freq_thresh_cm = float(neg_freq_thresh_cm)
        self.flatten_amp_ang = float(flatten_amp_ang)
        self.flatten_max_iter = int(flatten_max_iter)
        self.mem = int(mem)
        self.use_lobpcg = bool(use_lobpcg)  # deprecated; ignored
        self.root = int(root)
        self.dimer_kwargs = dict(dimer_kwargs or {})
        self.lbfgs_kwargs = dict(lbfgs_kwargs or {})
        self.max_total_cycles = int(max_total_cycles)

        # multi-mode flatten controls
        self.flatten_sep_cutoff = float(flatten_sep_cutoff)
        self.flatten_k = int(flatten_k)

        # Track total cycles globally across all loops/segments
        self._cycles_spent = 0

        # UMA settings
        self.uma_kwargs = dict(charge=0, spin=1, model="uma-s-1p1",
                               task_name="omol", device="auto") if uma_kwargs is None else dict(uma_kwargs)

        # Geometry & masses (use provided geom kwargs so freeze_atoms etc. apply)
        gkw = dict(geom_kwargs or {})
        coord_type = gkw.pop("coord_type", GEOM_KW_DEFAULT["coord_type"])
        self.geom = geom_loader(fn, coord_type=coord_type, **gkw)
        self.masses_amu = np.array([atomic_masses[z] for z in self.geom.atomic_numbers])
        self.masses_au_t = torch.as_tensor(self.masses_amu * AMU2AU, dtype=torch.float32)

        # Preserve freeze list (for PHVA)
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
        calc_sp = uma_pysis(**self.uma_kwargs)

        # Merge user dimer kwargs (but enforce N_raw & write_orientations)
        dimer_kwargs = dict(self.dimer_kwargs)
        dimer_kwargs.update({
            "calculator": calc_sp,
            "N_raw": str(self.mode_path),
            "write_orientations": False,
            "seed": 0,
            "mem": self.mem,
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
        Consumes from the global cycle budget self.max_total_cycles.
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
            coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
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
                    H_t, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root
                )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del H_t, coords_bohr_t, mode_xyz
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        return steps_in_this_call

    # ----- helpers for flatten selection -----
    def _representative_atoms_for_mode(self, mode_vec: torch.Tensor) -> np.ndarray:
        """
        Return indices of atoms with the largest displacements for a given mode.
        """
        vec = mode_vec.view(-1, 3)
        norms = torch.linalg.norm(vec, dim=1)
        k = min(self.flatten_k, vec.shape[0])
        if k <= 0:
            return np.zeros(0, dtype=int)
        topk = torch.topk(norms, k=k, largest=True)
        return topk.indices.detach().cpu().numpy()

    def _select_flatten_targets(self,
                                freqs_cm: np.ndarray,
                                modes: torch.Tensor) -> List[int]:
        """
        Select a subset of imaginary modes to flatten:
          * exclude the primary (TS) mode selected by `root`,
          * from the remaining imaginary modes, greedily pick modes whose
            representative atoms are separated by at least `flatten_sep_cutoff`
            (Å) from previously selected modes.
        """
        neg_idx_all = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx_all) <= 1:
            return []

        # sort imaginary modes (most negative first)
        order = np.argsort(freqs_cm[neg_idx_all])
        sorted_neg = neg_idx_all[order]

        # primary mode index in freqs_cm (TS mode to keep)
        root_clamped = max(0, min(self.root, len(order) - 1))
        primary_idx = sorted_neg[root_clamped]

        # candidates = all other imaginary modes
        candidates = [int(i) for i in sorted_neg if int(i) != int(primary_idx)]
        if not candidates:
            return []

        # atomic coordinates in Å
        coords_ang = torch.as_tensor(
            self.geom.cart_coords.reshape(-1, 3) * BOHR2ANG,
            dtype=modes.dtype,
            device=modes.device,
        )

        targets: List[int] = []
        reps_list: List[np.ndarray] = []

        for idx in candidates:
            rep = self._representative_atoms_for_mode(modes[idx])
            if rep.size == 0:
                continue
            rep_coords = coords_ang[rep]  # (k, 3)
            if not reps_list:
                targets.append(idx)
                reps_list.append(rep)
                continue

            # check distance to all previously selected representative sets
            accept = True
            for prev_rep in reps_list:
                prev_coords = coords_ang[prev_rep]  # (k, 3)
                dmat = torch.cdist(rep_coords, prev_coords)
                min_dist = float(torch.min(dmat).item())
                if min_dist < self.flatten_sep_cutoff:
                    accept = False
                    break
            if accept:
                targets.append(idx)
                reps_list.append(rep)

        return targets

    def _flatten_once_with_modes(self, freqs_cm: np.ndarray, modes: torch.Tensor) -> bool:
        """
        Flatten using precomputed (approximate) modes (mass-weighted, embedded).

        Only spatially separated extra imaginary modes are used, following a
        greedy selection similar to the PartialHessianDimer implementation.
        Modes are applied sequentially from the current geometry, so that each
        mode's displacement is decided from a 1D energy scan along that mode.
        """
        neg_idx_all = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx_all) <= 1:
            return False

        # choose targets based on spatial separation
        targets = self._select_flatten_targets(freqs_cm, modes)
        if not targets:
            return False

        # mass scaling (C moves exactly flatten_amp_ang Å)
        mass_scale = np.sqrt(12.011 / self.masses_amu)[:, None]
        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        # ensure energy reference is set up
        _ = _calc_energy(self.geom, self.uma_kwargs)

        # work in Bohr coordinates
        for idx in targets:
            # mode is currently mass-weighted & embedded to 3N
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
            # convert to Cartesian (Å-scale amplitude, but coords are Bohr)
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)

            disp = amp_bohr * mass_scale * v_cart  # Bohr
            ref = self.geom.cart_coords.reshape(-1, 3)

            plus = ref + disp
            minus = ref - disp

            self.geom.coords = plus.reshape(-1)
            E_plus = _calc_energy(self.geom, self.uma_kwargs)

            self.geom.coords = minus.reshape(-1)
            E_minus = _calc_energy(self.geom, self.uma_kwargs)

            # keep lower-energy side and continue from there
            self.geom.coords = (plus if E_plus <= E_minus else minus).reshape(-1)

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
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
        coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
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
                H_t, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
            )
        np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
        del mode_xyz, coords_bohr_t, H_t
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # (2) Loose loop (or initial pass)
        if self.root != 0:
            print("[WARNING] root != 0. Using this root only for the initial dimer loop.")
            print(f"Dimer loop with initial direction from mode {self.root}...")
            self.root = 0
            self.thresh_loose = self.thresh
        else:
            print("Loose dimer loop...")

        self._dimer_loop(self.thresh_loose)

        # (3) Update mode & normal loop
        H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
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
                H_t, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
            )
        np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
        del mode_xyz, coords_bohr_t, H_t
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        print("Normal dimer loop...")
        self._dimer_loop(self.thresh)

        # (4) Flatten loop — exact Hessian each iteration & Bofill only for flatten
        print("Flatten loop with Bofill-updated active Hessian (flatten displacements only)...")

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

        # Flatten iterations with approximate Hessian updates
        for it in range(self.flatten_max_iter):
            if (self.max_total_cycles - self._cycles_spent) <= 0:
                break

            # (a) Approximate frequencies & modes from the active Hessian
            freqs_cm_approx, modes_embedded = _modes_from_Hact_embedded(
                H_act, self.geom.atomic_numbers,
                self.geom.cart_coords.reshape(-1, 3), active_idx, self.device
            )
            n_imag = int(np.sum(freqs_cm_approx < -abs(self.neg_freq_thresh_cm)))
            approx_ims = [float(x) for x in freqs_cm_approx if x < -abs(self.neg_freq_thresh_cm)]
            print(f"[IMAG~] n≈{n_imag}  (approx imag: {approx_ims})")
            if n_imag <= 1:
                break

            # (b) Flatten using the approximate modes
            x_before_flat = self.geom.cart_coords.copy().reshape(-1)
            g_before_flat = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)

            did_flatten = self._flatten_once_with_modes(freqs_cm_approx, modes_embedded)
            if not did_flatten:
                break

            x_after_flat = self.geom.cart_coords.copy().reshape(-1)
            g_after_flat = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)

            # (c) Bofill update using UMA gradients across the flatten displacement
            delta_flat_full = x_after_flat - x_before_flat
            delta_flat_act = delta_flat_full[mask_dof]
            g_old_act = g_before_flat[mask_dof]
            g_new_act = g_after_flat[mask_dof]
            H_act = _bofill_update_active(H_act, delta_flat_act, g_new_act, g_old_act)

            # (d) Refresh dimer direction from updated active Hessian
            mode_xyz = _mode_direction_by_root_from_Hact(
                H_act, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
            )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")

            # (e) Re-optimize with Dimer (consumes global cycle budget)
            self._dimer_loop(self.thresh)

            # (f) After dimer optimization, recompute an exact Hessian for the next iteration
            H_any = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
            if H_any.size(0) == 3 * N:
                H_act = _extract_active_block(H_any, mask_dof)
            else:
                H_act = H_any
            del H_any
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

        # (5) Final outputs
        final_xyz = self.out_dir / "final_geometry.xyz"
        atoms_final = Atoms(self.geom.atoms, positions=(self.geom.cart_coords.reshape(-1, 3) * BOHR2ANG), pbc=False)
        write(final_xyz, atoms_final)

        # Final Hessian → imaginary mode animation
        H_t = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        if H_t.size(0) == 3 * N:
            freqs_cm, modes = _frequencies_cm_and_modes(
                H_t, self.geom.atomic_numbers, self.geom.cart_coords.reshape(-1, 3), self.device,
                freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
            )
        else:
            freqs_cm, modes = _modes_from_Hact_embedded(
                H_t, self.geom.atomic_numbers, self.geom.cart_coords.reshape(-1, 3),
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
            mode_mw = modes[primary_idx]
            masses_amu_t = torch.as_tensor(self.masses_amu, dtype=mode_mw.dtype, device=mode_mw.device)
            m3 = torch.repeat_interleave(masses_amu_t, 3)
            v_cart = (mode_mw / torch.sqrt(m3)).detach().cpu().numpy()
            v_cart = v_cart / np.linalg.norm(v_cart)
            del modes, masses_amu_t, m3
            out_trj = self.vib_dir / f"final_imag_mode_{freqs_cm[primary_idx]:+.2f}cm-1.trj"
            out_pdb = self.vib_dir / f"final_imag_mode_{freqs_cm[primary_idx]:+.2f}cm-1.pdb"
            _write_mode_trj_and_pdb(
                self.geom,
                v_cart,
                out_trj,
                amplitude_ang=0.8,
                n_frames=20,
                comment=f"imag {freqs_cm[primary_idx]:+.2f} cm-1",
                ref_pdb=self.ref_pdb,
                write_pdb=self.ref_pdb is not None,
                out_pdb=out_pdb,
            )

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        print(f"[DONE] Saved final geometry → {final_xyz}")
        print(f"[DONE] Mode files → {self.vib_dir}")


# ===================================================================
#                         Defaults for CLI
# ===================================================================

# Geometry defaults
GEOM_KW = dict(GEOM_KW_DEFAULT)

# UMA calculator defaults
CALC_KW = dict(_UMA_CALC_KW)

# Optimizer base (common) — used by both RSIRFO and the inner LBFGS of HessianDimer
OPT_BASE_KW = dict(_OPT_BASE_KW)
OPT_BASE_KW.update({
    "out_dir": "./result_tsopt/",
})

DIMER_KW = {
    # --- Geometry / length ---
    "length": 0.0189,                 # dimer half-length in Bohr (≈ 0.01 Å)

    # --- Rotation loop settings ---
    "rotation_max_cycles": 15,        # max rotation cycles per update
    "rotation_method": "fourier",     # "fourier" (robust) | "direct" (steepest-like)
    "rotation_thresh": 1e-4,          # threshold on ||rot_force||
    "rotation_tol": 1,                # degrees; skip rotation if |angle| < tol
    "rotation_max_element": 0.001,    # max element for direct rotation step (Bohr)
    "rotation_interpolate": True,     # interpolate f1 during rotation (Fourier method)
    "rotation_disable": False,        # disable rotation (use given N as-is)
    "rotation_disable_pos_curv": True,# if curvature positive after rotation, restore previous N
    "rotation_remove_trans": True,    # remove net translation from N

    # --- Translational part of projected force ---
    "trans_force_f_perp": True,       # include perpendicular force into translational component when C<0

    # --- Initial guess helpers (mutually exclusive) ---
    "bonds": None,                    # Optional[List[Tuple[int,int]]], use weighted-bond mode as N_raw
    "N_hessian": None,                # Optional[str], path to HDF5 Hessian; use 1st imag. mode as N_raw

    # --- Optional bias potentials for robustness ---
    "bias_rotation": False,           # quadratic bias along initial N during rotation
    "bias_translation": False,        # add Gaussians along the path to stabilize translation
    "bias_gaussian_dot": 0.1,         # target dot product when tuning Gaussian height

    # --- Stochastic control / IO ---
    "seed": None,                     # RNG seed; Runner sets to 0 for determinism
    "write_orientations": True,       # write N.trj each call; Runner overrides to False to reduce IO

    # --- Hessian forwarding ---
    "forward_hessian": True,          # allow forwarding get_hessian to wrapped calculator
}

# Reference: internal L-BFGS defaults for TS optimization highlighting deviations from OPT_BASE_KW
LBFGS_TS_KW: Dict[str, Any] = dict(_LBFGS_KW)

# HessianDimer defaults (CLI-level)
hessian_dimer_KW = {
    "thresh_loose": "gau_loose",      # loose threshold preset for first pass
    "thresh": "gau",                  # main threshold preset for TS search
    "update_interval_hessian": 1000,  # LBFGS cycles per Hessian refresh for direction
    "neg_freq_thresh_cm": 5.0,        # treat ν < -this as imaginary (cm^-1)
    "flatten_amp_ang": 0.10,          # mass-scaled displacement amplitude for flattening (Å)
    "flatten_max_iter": 20,           # max flattening iterations
    "flatten_sep_cutoff": 2.0,        # minimum distance between representative atoms (Å)
    "flatten_k": 10,                  # number of representative atoms per mode
    "mem": 100000,                    # scratch/IO memory passed through Calculator (**kwargs)
    "use_lobpcg": True,               # deprecated (ignored)
    "device": "auto",                 # "cuda"|"cpu"|"auto" for torch-side ops
    "root": 0,                        # 0: follow the most negative mode
    "dimer": {**DIMER_KW},            # default kwargs forwarded to Dimer (Runner may override some)
    "lbfgs": {**LBFGS_TS_KW},         # default kwargs forwarded to inner LBFGS
}

# RSIRFO (TS Hessian optimizer) defaults (subset; additional keys may be provided)
RSIRFO_KW: Dict[str, Any] = dict(_RFO_KW)
RSIRFO_KW.update({
    "roots": [0],               # mode indices to follow uphill
    "hessian_ref": None,        # reference Hessian file (HDF5)
    "rx_modes": None,           # reaction-mode definitions for projection
    "prim_coord": None,         # primary coordinates for monitoring
    "rx_coords": None,          # reaction coordinates for monitoring
    "hessian_update": "bofill", # override base "bfgs"
    "hessian_recalc_reset": True,# reset recalc counter after exact Hessian
    "max_micro_cycles": 50,     # micro-iterations per macro cycle
    "augment_bonds": False,     # augment reaction path based on bond analysis
    "min_line_search": True,    # enforce minimum line-search step
    "max_line_search": True,    # enforce maximum line-search step
    "assert_neg_eigval": False, # ensure a negative eigenvalue at convergence
})


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
@click.option("-q", "--charge", type=int, required=False, help="Charge of the ML region.")
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1) for the ML region.")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze parent atoms of link hydrogens (PDB only).")
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Max cycles / steps cap")
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="light (=Dimer) or heavy (=RSIRFO)",
)
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump optimization trajectory")
@click.option("--out-dir", type=str, default="./result_tsopt/", show_default=True, help="Output directory")
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help="Convergence preset for the active optimizer (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, opt, dimer, rsirfo).",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="Choose UMA Hessian evaluation mode (overrides YAML/calc.hessian_calc_mode). Defaults to 'FiniteDifference'.",
)
def cli(
    input_path: Path,
    charge: Optional[int],
    spin: Optional[int],
    freeze_links: bool,
    convert_files: bool,
    max_cycles: int,
    opt_mode: str,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    hessian_calc_mode: Optional[str],
) -> None:
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)
    time_start = time.perf_counter()

    # --------------------------
    # 1) Assemble configuration (defaults ← CLI ← YAML)
    # --------------------------
    yaml_cfg = load_yaml_dict(args_yaml)
    geom_cfg = dict(GEOM_KW)
    calc_cfg = dict(CALC_KW)
    opt_cfg  = dict(OPT_BASE_KW)
    simple_cfg = dict(hessian_dimer_KW)
    rsirfo_cfg = dict(RSIRFO_KW)

    # CLI overrides
    calc_cfg["charge"] = int(charge)
    calc_cfg["spin"]   = int(spin)
    opt_cfg["max_cycles"] = int(max_cycles)
    opt_cfg["dump"]       = bool(dump)
    opt_cfg["out_dir"]    = out_dir
    if thresh is not None:
        opt_cfg["thresh"] = str(thresh)
        simple_cfg["thresh"] = str(thresh)
        rsirfo_cfg["thresh"] = str(thresh)

    # Hessian mode override from CLI
    if hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

    # YAML overrides (highest precedence)
    apply_yaml_overrides(
        yaml_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (opt_cfg,  (("opt",),)),
            (simple_cfg, (("hessian_dimer",),)),
            (rsirfo_cfg, (("rsirfo",),)),
        ],
    )

    # Freeze links (PDB only): merge with existing list
    if freeze_links and source_path.suffix.lower() == ".pdb":
        try:
            detected = detect_freeze_links(source_path)
        except Exception as e:
            click.echo(f"[freeze-links] WARNING: Could not detect link parents: {e}", err=True)
            detected = []
        merged = merge_freeze_atom_indices(geom_cfg, detected)
        if merged:
            click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, merged))}")

    # Pass freeze_atoms from geom → calc (so UMA knows active DOF for FD Hessian)
    if "freeze_atoms" in geom_cfg:
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))

    kind = normalize_choice(
        opt_mode,
        param="--opt-mode",
        alias_groups=_OPT_MODE_ALIASES,
        allowed_hint="light|heavy",
    )
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
            runner = HessianDimer(
                fn=str(geom_input_path),
                out_dir=str(out_dir_path),
                thresh_loose=simple_cfg.get("thresh_loose", "gau_loose"),
                thresh=simple_cfg.get("thresh", "gau"),
                update_interval_hessian=int(simple_cfg.get("update_interval_hessian", 200)),
                neg_freq_thresh_cm=float(simple_cfg.get("neg_freq_thresh_cm", 5.0)),
                flatten_amp_ang=float(simple_cfg.get("flatten_amp_ang", 0.10)),
                flatten_max_iter=int(simple_cfg.get("flatten_max_iter", 20)),
                mem=int(simple_cfg.get("mem", 100000)),
                use_lobpcg=bool(simple_cfg.get("use_lobpcg", True)),  # deprecated; ignored
                uma_kwargs=uma_kwargs_for_sd,
                device=str(simple_cfg.get("device", calc_cfg.get("device", "auto"))),
                dump=bool(opt_cfg["dump"]),
                root=int(simple_cfg.get("root", 0)),
                dimer_kwargs=dict(simple_cfg.get("dimer", {})),
                lbfgs_kwargs=dict(simple_cfg.get("lbfgs", {})),
                max_total_cycles=int(opt_cfg["max_cycles"]),
                flatten_sep_cutoff=float(simple_cfg.get("flatten_sep_cutoff", 2.0)),
                flatten_k=int(simple_cfg.get("flatten_k", 10)),
                # Propagate geometry settings (freeze_atoms, coord_type, ...) to the HessianDimer runner
                geom_kwargs=dict(geom_cfg),
            )

            click.echo("\n=== TS optimization (Hessian Dimer) started ===\n")
            runner.run()
            click.echo("\n=== TS optimization (Hessian Dimer) finished ===\n")

            needs_pdb = source_path.suffix.lower() == ".pdb"
            needs_gjf = prepared_input.is_gjf
            ref_pdb = source_path.resolve() if needs_pdb else None
            final_xyz = out_dir_path / "final_geometry.xyz"
            try:
                convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=out_dir_path / "final_geometry.pdb" if needs_pdb else None,
                    out_gjf_path=out_dir_path / "final_geometry.gjf" if needs_gjf else None,
                )
                click.echo("[convert] Wrote 'final_geometry' outputs.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final geometry: {e}", err=True)

            if bool(opt_cfg.get("dump", False)) and needs_pdb:
                all_trj = out_dir_path / "optimization_all.trj"
                if all_trj.exists():
                    try:
                        convert_xyz_like_outputs(
                            all_trj,
                            prepared_input,
                            ref_pdb_path=ref_pdb,
                            out_pdb_path=out_dir_path / "optimization_all.pdb" if needs_pdb else None,
                        )
                        click.echo("[convert] Wrote 'optimization_all' outputs.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert optimization trajectory: {e}", err=True)
                else:
                    click.echo("[convert] WARNING: 'optimization_all.trj' not found; skipping conversion.", err=True)

        else:
            # RS-I-RFO (heavy)
            #  - Build geometry now and attach UMA calculator
            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)
            geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

            calc_builder_or_instance = uma_pysis(**calc_cfg)  # includes freeze_atoms / hessian_calc_mode / partial Hessian
            try:
                geometry.set_calculator(calc_builder_or_instance())
            except TypeError:
                geometry.set_calculator(calc_builder_or_instance)

            # Construct RSIRFO optimizer
            rs_args = dict(rsirfo_cfg)
            opt_base = dict(opt_cfg)
            opt_base["out_dir"] = str(out_dir_path)
            rs_args["out_dir"] = str(out_dir_path)

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

            rsirfo_kwargs = {**opt_base, **rs_args}
    
            for diis_kw in ("gediis", "gdiis", "gdiis_thresh", "gediis_thresh", "gdiis_test_direction", "adapt_step_func"):
                rsirfo_kwargs.pop(diis_kw, None)

            optimizer = RSIRFOptimizer(geometry, **rsirfo_kwargs)

            click.echo("\n=== TS optimization (RS-I-RFO) started ===\n")
            optimizer.run()
            click.echo("\n=== TS optimization (RS-I-RFO) finished ===\n")

            needs_pdb = source_path.suffix.lower() == ".pdb"
            needs_gjf = prepared_input.is_gjf
            ref_pdb = source_path.resolve() if needs_pdb else None
            final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
            try:
                convert_xyz_like_outputs(
                    final_xyz_path,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=out_dir_path / "final_geometry.pdb" if needs_pdb else None,
                    out_gjf_path=out_dir_path / "final_geometry.gjf" if needs_gjf else None,
                )
                click.echo("[convert] Wrote 'final_geometry' outputs.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final geometry: {e}", err=True)

            if bool(opt_cfg.get("dump", False)) and needs_pdb:
                trj_path = out_dir_path / "optimization.trj"
                if trj_path.exists():
                    try:
                        convert_xyz_like_outputs(
                            trj_path,
                            prepared_input,
                            ref_pdb_path=ref_pdb,
                            out_pdb_path=out_dir_path / "optimization.pdb" if needs_pdb else None,
                        )
                        click.echo("[convert] Wrote 'optimization' outputs.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert optimization trajectory: {e}", err=True)
                else:
                    click.echo("[convert] WARNING: 'optimization.trj' not found; skipping conversion.", err=True)

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

            # Use configurable neg_freq_thresh_cm (same default as light mode)
            neg_freq_thresh_cm = float(simple_cfg.get("neg_freq_thresh_cm", 5.0))
            neg_idx = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
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
                ref_pdb = source_path if source_path.suffix.lower() == ".pdb" else None
                _write_mode_trj_and_pdb(
                    gproxy,
                    v_cart,
                    out_trj,
                    amplitude_ang=0.8,
                    n_frames=20,
                    comment=f"imag {freqs_cm[pick_idx]:+.2f} cm-1",
                    ref_pdb=ref_pdb,
                    write_pdb=ref_pdb is not None,
                    out_pdb=out_pdb,
                )

            if torch.cuda.is_available():
                torch.cuda.empty_cache()

        click.echo(format_elapsed("[time] Elapsed Time for TS Opt", time_start))

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
    finally:
        prepared_input.cleanup()
