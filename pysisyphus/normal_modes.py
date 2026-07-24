# pysisyphus/normal_modes.py

"""Pure normal-mode kernel: mass weighting, rigid projection, diagonalization.

This is a *lower* bundled-engine module: a sibling of ``pysisyphus.tr_projection``
that sits BELOW the product workflows. It must import neither ``pdb2reaction`` nor
``mlmm`` so that bundled-engine code (e.g.
``pysisyphus.tsoptimizers.TSHessianOptimizer``) can consume the mass/mode kernel
without an upward product dependency.

The bounded-peak Hessian symmetrizer lives here as well so the module is
self-contained; ``pdb2reaction.core.utils`` re-exports it for backward
compatibility.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
import torch
from ase.data import atomic_masses
import ase.units as units

from pysisyphus.constants import BOHR2ANG, AMU2AU, AU2EV
from pysisyphus._array import active_square
from pysisyphus.tr_projection import active_tr_basis, project_hessian_inplace


def symmetrize_inplace(H, chunk: int = 512):
    """Symmetrize a square Hessian-like tensor in place with bounded peak VRAM.

    Replaces the 2x-peak idiom ``_t = H.T.clone(); H.add_(_t).mul_(0.5); del _t``
    with a chunked average that writes BOTH triangles symmetrically (no
    upper-triangle-only tricks). Peak extra allocation is bounded by
    ``chunk * chunk`` elements (vs ``N * N`` for the naive form).
    """

    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(
            f"symmetrize_inplace expects a square 2-D tensor, got shape {tuple(H.shape)}"
        )
    N = H.shape[0]
    if N == 0:
        return H

    if N <= chunk:
        tmp = H.T.contiguous()
        H.add_(tmp).mul_(0.5)
        del tmp
        return H

    for i in range(0, N, chunk):
        ie = min(i + chunk, N)
        diag = H[i:ie, i:ie]
        diag_tmp = diag.T.contiguous()
        # Assign out-of-place: an in-place add_ on this strided diagonal view is
        # self-overlapping when the block is 1x1 (N % chunk == 1), which torch
        # rejects with a RuntimeError. Use the same write pattern as the
        # off-diagonal blocks below.
        H[i:ie, i:ie] = diag.add(diag_tmp).mul(0.5)
        del diag_tmp
        for j in range(ie, N, chunk):
            je = min(j + chunk, N)
            upper = H[i:ie, j:je]
            lower_T = H[j:je, i:ie].T
            avg = upper.add(lower_T).mul_(0.5)
            upper.copy_(avg)
            H[j:je, i:ie].copy_(avg.T)
            del avg
    return H


def _safe_masses_amu(atomic_numbers: Sequence[int]) -> np.ndarray:
    """Look up atomic masses with a clear error for unknown atomic numbers."""
    max_z = len(atomic_masses) - 1
    bad = [z for z in atomic_numbers if z < 0 or z > max_z or atomic_masses[z] == 0.0]
    if bad:
        raise ValueError(
            f"Unknown or unsupported atomic number(s): {sorted(set(bad))}. "
            "Check that all elements in the input structure are valid."
        )
    return np.array([atomic_masses[z] for z in atomic_numbers])


def _mw_projected_hessian(H: torch.Tensor,
                          coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor,
                          tr_projection: str = "constrained",
                          projection_info: Optional[dict] = None) -> torch.Tensor:
    """
    Project out translations/rotations in mass-weighted space:
    Hmw = M^{-1/2} H M^{-1/2};  P = I - QQ^T;  Hmw_proj = P Hmw P

    To save memory, update **H in-place** (no clone) and return it.
    The output is explicitly symmetrized after TR projection.
    """
    if H.dtype != torch.float64:
        H = H.to(dtype=torch.float64)
    dtype, device = H.dtype, H.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        inv_sqrt_m_col = inv_sqrt_m.view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m.view(-1, 1)

        H.mul_(inv_sqrt_m_row)
        H.mul_(inv_sqrt_m_col)

        Q, info = active_tr_basis(
            coords_bohr_t,
            masses_au_t,
            list(range(int(coords_bohr_t.shape[0]))),
            mode=tr_projection,
        )
        project_hessian_inplace(H, Q)
        if projection_info is not None:
            projection_info.clear()
            projection_info.update(info.as_dict())

        # Bounded-peak symmetrization (helper writes both triangles; peak temp
        # <= chunk^2 instead of full N×N clone).
        symmetrize_inplace(H)

        del masses_amu_t, m3, inv_sqrt_m, inv_sqrt_m_col, inv_sqrt_m_row
        del Q

        if torch.cuda.is_available() and device.type == "cuda":
            torch.cuda.empty_cache()
        return H


# PHVA helper: mass-weighted Hessian without TR projection (for active subspace)
def _mass_weighted_hessian(H: torch.Tensor,
                           masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Return Hmw = M^{-1/2} H M^{-1/2} (no symmetrization/TR projection; in-place).
    """
    dtype, device = H.dtype, H.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        inv_sqrt_m_col = inv_sqrt_m.view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m.view(-1, 1)
        H.mul_(inv_sqrt_m_row)
        H.mul_(inv_sqrt_m_col)
        del masses_amu_t, m3, inv_sqrt_m, inv_sqrt_m_col, inv_sqrt_m_row
        return H


def _frequencies_cm_and_modes(H: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              # tol is a mass-weighted eigenvalue (ω², au) floor for numerical-noise
                              # removal: 1e-6 au(ω²) ≈ 5.14 cm⁻¹, so near-zero modes (FD/analytical
                              # Hessian noise, translation/rotation residue) are dropped before the
                              # imaginary count and thermochemistry. NOT a physical soft-mode cutoff;
                              # tsopt counts imaginary at neg_freq_thresh_cm=5.0 (≈ the same floor).
                              tol: float = 1e-6,
                              freeze_idx: Optional[List[int]] = None,
                              tr_projection: str = "constrained",
                              projection_info: Optional[dict] = None) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize a (possibly PHVA/active-subspace) TR-projected mass-weighted Hessian
    to obtain frequencies (cm^-1) and mass-weighted eigenvectors (modes).

    If `freeze_idx` is provided (list of 0-based atom indices), perform
    Partial Hessian Vibrational Analysis (PHVA). Supports two cases:

      A) Full Hessian given (3N×3N):
         1) mass-weight the full Hessian
         2) take the active subspace by removing DOF of frozen atoms
         3) remove only constrained-system rigid null modes represented in the active subspace
         4) diagonalize and embed eigenvectors back to 3N by zero-filling frozen DOF

      B) Already-reduced (active-block) Hessian given (3N_act×3N_act):
         1) mass-weight using only active masses
         2) apply the same constrained-system rigid-null treatment in active space
         3) diagonalize and embed back to 3N by zero-filling frozen DOF

    Returns:
      freqs_cm : (nmode,) numpy, negatives are imaginary
      modes    : (nmode, 3N) torch (mass-weighted eigenvectors)
    """
    with torch.no_grad():
        if H.dtype != torch.float64:
            H = H.to(dtype=torch.float64)
        Z = np.array(atomic_numbers, dtype=int)
        N = int(len(Z))
        masses_amu = _safe_masses_amu(Z)
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H.dtype, device=device)

        if freeze_idx is not None and len(freeze_idx) > 0:
            frozen_set = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen_set]
            n_active = len(active_idx)
            if n_active == 0:
                raise ValueError("PHVA requires at least one active atom")

            expected_act_dim = 3 * n_active
            is_partial = (H.shape[0] == expected_act_dim and H.shape[1] == expected_act_dim)
            is_full = (H.shape[0] == 3 * N and H.shape[1] == 3 * N)
            if not (is_partial or is_full):
                raise ValueError(
                    "Hessian shape is inconsistent with the full and active PHVA spaces: "
                    f"got {tuple(H.shape)}, expected {(3 * N, 3 * N)} or "
                    f"{(expected_act_dim, expected_act_dim)}"
                )
            Q, info = active_tr_basis(
                coords_bohr_t,
                masses_au_t,
                active_idx,
                mode=tr_projection,
            )
            if projection_info is not None:
                projection_info.clear()
                projection_info.update(info.as_dict())

            if is_partial:
                masses_act = masses_au_t[active_idx]
                Hmw_act = _mass_weighted_hessian(H, masses_act)
                project_hessian_inplace(Hmw_act, Q)
                # Bounded-peak symmetrization (helper writes both triangles).
                symmetrize_inplace(Hmw_act)

                omega2, Vsub = torch.linalg.eigh(Hmw_act)

                del Hmw_act
                del H
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                sel = torch.abs(omega2) > tol
                omega2 = omega2[sel]
                Vsub = Vsub[:, sel]

                modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Vsub.dtype, device=Vsub.device)
                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Vsub.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False
                modes[:, mask_dof] = Vsub.T
                del Q, mask_dof

            else:
                H = _mass_weighted_hessian(H, masses_au_t)

                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=H.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False

                active_dof = torch.nonzero(mask_dof, as_tuple=False).flatten()
                H = active_square(H, active_dof)
                del active_dof
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                project_hessian_inplace(H, Q)
                # Bounded-peak symmetrization (helper writes both triangles).
                symmetrize_inplace(H)

                omega2, Vsub = torch.linalg.eigh(H)

                del H
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                sel = torch.abs(omega2) > tol
                omega2 = omega2[sel]
                Vsub = Vsub[:, sel]

                modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Vsub.dtype, device=Vsub.device)
                modes[:, mask_dof] = Vsub.T
                del Vsub, mask_dof, Q

        else:
            H = _mw_projected_hessian(
                H,
                coords_bohr_t,
                masses_au_t,
                tr_projection=tr_projection,
                projection_info=projection_info,
            )
            omega2, V = torch.linalg.eigh(H)

            del H
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            sel = torch.abs(omega2) > tol
            omega2 = omega2[sel]
            modes = V[:, sel].T
            del V

        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()

        del omega2, hnu, sel
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return freqs_cm, modes


def _mw_mode_to_cart(mode_mw_3N_t: torch.Tensor,
                     masses_au_t: torch.Tensor) -> np.ndarray:
    """
    Convert one mass-weighted eigenvector (3N,) to Cartesian (3N,) and L2-normalize.
    """
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=mode_mw_3N_t.dtype, device=mode_mw_3N_t.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        v_cart = torch.sqrt(1.0 / m3) * mode_mw_3N_t
        v_cart.div_(torch.linalg.norm(v_cart))
        arr = v_cart.detach().cpu().numpy()
        del masses_amu_t, m3, v_cart
        return arr
