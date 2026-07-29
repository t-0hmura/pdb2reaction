"""Transition-state optimization with Dimer and rational-function methods."""
# DOMAIN_PURE

from __future__ import annotations

import gc
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import click
from pdb2reaction.core.output import emit
import numpy as np
import torch
from ase import Atoms
import ase.units as units
from ase.io import write

# ---------------- pysisyphus / pdb2reaction imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.intcoords.exceptions import RebuiltInternalsException
from pysisyphus.constants import BOHR2ANG, AMU2AU, AU2EV
from pysisyphus._array import active_square
from pysisyphus.calculators.Dimer import Dimer  # Dimer calculator (orientation-projected forces)
from pysisyphus.tr_projection import (
    active_tr_basis,
    allows_saddle_certification,
    compact_project_hessian,
    full_cartesian_tr_basis,
    normalize_tr_projection_mode,
    project_hessian_inplace,
)

# Additional TS optimizer classes for --opt-mode trim, rsprfo
from pysisyphus.tsoptimizers.TRIM import TRIM
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer

# RS-I-RFO optimizer (transition state)
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer  # type: ignore

# local helpers from pdb2reaction
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    THRESH_CHOICES,
    UMA_CALC_KW,
    OPT_BASE_KW,
    TSOPT_MODE_ALIASES,
    TS_IMAG_SOFT_WARN_CM,
    OUT_DIR_TSOPT,
    HESSIAN_DIMER_CLI_KW,
    RSIRFO_KW,
    DEFAULT_UMA_MODEL,
    apply_backend_defaults,
)
from pdb2reaction.core.utils import (
    resolve_freeze_atoms,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    strip_inherited_keys,
    cli_param_overridden,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
    PreparedInputStructure,
    emit_optimizer_terminal_status,
    optimizer_cycle_count,
    optimizer_terminal_status,
)
from pdb2reaction.cli.common_options import (
    add_ml_charge_spin_options,
    add_precision_option, add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_coord_type_option,
    add_print_every_option,
)
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, _write_error_json, render_cli_exception
from pdb2reaction.workflows.freq import (
    _torch_device,
    _safe_masses_amu,
    _mass_weighted_hessian,
    _calc_full_hessian_torch,
    _calc_energy,
    _write_mode_trj_and_pdb,
    _frequencies_cm_and_modes,
)

logger = logging.getLogger(__name__)


def _set_cartesian_flatten_coords(geom, cart_coords: np.ndarray) -> None:
    """Install a Cartesian TS trial and accept a completed internals rebuild."""

    try:
        geom.cart_coords = np.asarray(cart_coords, dtype=float).reshape(-1)
    except RebuiltInternalsException:
        # The Cartesian setter installs the requested geometry before it uses
        # this exception to report that its primitive internals were rebuilt.
        geom.clear()


# TS optimizer class map. All three classes inherit from TSHessianOptimizer
# and share the kwargs surface built by `_build_rsirfo_kwargs`, so they are
# interchangeable at construction time.
TSOPT_CLASS_MAP = {"rsirfo": RSIRFOptimizer, "trim": TRIM, "rsprfo": RSPRFOptimizer}
TSOPT_DISPLAY_NAME = {"rsirfo": "RS-I-RFO", "rsprfo": "RS-P-RFO", "trim": "TRIM"}
# Total Cartesian displacement norms.  These remain small compared with a
# typical neighboring-image separation, but the second shell is large enough
# to escape a shallow minimum reached from a coarse HEI.
PATH_MODE_RESTART_AMPLITUDES_ANG = (-0.10, 0.10, -0.20, 0.20)
# A flatten displacement perturbs both the surplus mode and the saddle mode it
# must retain.  One exact retry cycle can stop before the latter is
# re-established and was observed to soften both modes together.  Reuse the
# normal three-check persistence window before advancing to another
# displacement.
FLATTEN_RETRY_HIGHER_ORDER_CHECKS = 3


def _mirrored_flatten_start(
    saddle_coords: np.ndarray,
    primary_start: np.ndarray,
) -> np.ndarray:
    """Return the opposite signed displacement about a saddle geometry."""
    saddle = np.asarray(saddle_coords, dtype=float)
    primary = np.asarray(primary_start, dtype=float)
    if saddle.shape != primary.shape:
        raise ValueError("saddle and flatten-start coordinates must have equal shapes")
    return 2.0 * saddle - primary


def _flatten_branch_needs_alternate(result: Dict[str, Any]) -> bool:
    """Whether a signed flatten trial failed the requested first-order test."""
    optimizer = result["optimizer"]
    target_negative = (
        getattr(optimizer, "_last_exact_target_mode_is_negative", None)
        is True
    )
    return int(result["n_imag"]) != 1 or not target_negative


def _effective_flatten_iterations(
    configured: int,
    *,
    has_reference_mode: bool,
    n_imag: int,
    target_mode_is_negative: Optional[bool],
) -> Tuple[int, bool]:
    """Return the explicit flatten budget and whether path safety vetoed it."""
    iterations = max(int(configured), 0)
    vetoed = (
        iterations > 0
        and has_reference_mode
        and int(n_imag) > 1
        and target_mode_is_negative is not True
    )
    return (0 if vetoed else iterations), vetoed


def _transported_path_mode_full(
    optimizer,
    geometry,
    fallback: Optional[np.ndarray],
) -> Optional[np.ndarray]:
    """Return the latest path mode as a normalized full Cartesian vector.

    Post-Hessian flattening creates a fresh optimizer. Carrying the previous
    optimizer's PHVA/overlap-transported mode into that restart prevents the
    new instance from jumping back to a different root selected by the older
    HEI tangent. Non-Cartesian modes cannot be expanded this way and retain the
    public Cartesian fallback.
    """
    mode = None
    if geometry.coord_type in ("cart", "cartesian"):
        mode = getattr(optimizer, "_last_exact_physical_mode", None)
        if mode is None:
            mode = getattr(optimizer, "_physical_ts_mode", None)
        if mode is None:
            mode = getattr(optimizer, "_saddle_recovery_mode", None)
        if mode is None and len(getattr(optimizer, "ts_modes", ())):
            mode = optimizer.ts_modes[0]
        if isinstance(mode, torch.Tensor):
            mode = mode.detach().cpu().numpy()
        if mode is not None:
            mode = np.asarray(mode, dtype=float).reshape(-1)
            if mode.size != geometry.cart_coords.size:
                try:
                    mode = optimizer.full_from_active(mode)
                except (AttributeError, IndexError, ValueError):
                    mode = None

    if mode is None and fallback is not None:
        mode = np.asarray(fallback, dtype=float).reshape(-1)
    if mode is None or mode.size != geometry.cart_coords.size:
        return None
    norm = float(np.linalg.norm(mode))
    if not np.isfinite(norm) or norm <= 0.0:
        return None
    return mode / norm


def _initial_path_root_mode_full(optimizer, geometry) -> Optional[np.ndarray]:
    """Return the HEI's initially selected soft root in full Cartesian space."""
    if geometry.coord_type not in ("cart", "cartesian"):
        return None
    mode = getattr(optimizer, "_initial_reference_root_mode", None)
    if isinstance(mode, torch.Tensor):
        mode = mode.detach().cpu().numpy()
    if mode is None:
        return None
    mode = np.asarray(mode, dtype=float).reshape(-1)
    if mode.size != geometry.cart_coords.size:
        try:
            mode = optimizer.full_from_active(mode)
        except (AttributeError, IndexError, ValueError):
            return None
    if mode.size != geometry.cart_coords.size:
        return None
    norm = float(np.linalg.norm(mode))
    if not np.all(np.isfinite(mode)) or norm <= 0.0:
        return None
    return mode / norm


def _path_restart_mode_candidates(
    optimizer,
    geometry,
    reference_mode: np.ndarray,
) -> List[Tuple[str, np.ndarray]]:
    """Build bounded restart directions, preferring the explicit MEP tangent.

    A coarse or kinked discrete path can distribute its tangent over a soft
    reaction root and a much stiffer vibration.  The MEP tangent remains the
    first choice.  If it differs materially from the soft root selected by the
    initial overlap-band/curvature rule, retain that root as a second bounded
    shell rather than increasing the displacement without limit.
    """
    reference = np.asarray(reference_mode, dtype=float).reshape(-1)
    reference_norm = float(np.linalg.norm(reference))
    if (
        reference.size != geometry.cart_coords.size
        or not np.all(np.isfinite(reference))
        or reference_norm <= 0.0
    ):
        return []
    reference = reference / reference_norm
    candidates = [("mep-tangent", reference)]

    soft_root = _initial_path_root_mode_full(optimizer, geometry)
    if soft_root is not None:
        overlap = abs(float(np.dot(reference, soft_root)))
        if overlap < 0.95:
            candidates.append(("initial-soft-root", soft_root))
    return candidates


def _mw_projected_hessian_inplace(H: torch.Tensor,
                                  coords_bohr_t: torch.Tensor,
                                  masses_au_t: torch.Tensor,
                                  freeze_idx: Optional[List[int]] = None,
                                  tr_projection: str = "constrained",
                                  projection_info: Optional[dict] = None,
                                  compact: bool = False):
    """
    Mass-weight H in-place, optionally restrict to active DOF subspace (PHVA) and
    project out TR motions (in that subspace), also in-place. Symmetrizes after TR projection.
    Returns the active Hessian to be diagonalized.  With ``compact=True``,
    diagonalization is performed in the true orthogonal complement of the
    rigid basis and ``(H_compact, lift)`` is returned.  This avoids selecting
    artificial zero eigenvectors introduced by a same-size ``P H P`` matrix.
    """
    if H.dtype != torch.float64:
        H = H.to(dtype=torch.float64)
    dtype, device = H.dtype, H.device
    with torch.no_grad():
        N = coords_bohr_t.shape[0]
        if freeze_idx is not None and len(freeze_idx) > 0:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen]
            if len(active_idx) == 0:
                raise RuntimeError("All atoms are frozen; no active DOF left for TR projection.")
            # mass-weight first
            H = _mass_weighted_hessian(H, masses_au_t.to(dtype=dtype, device=device))
            # take active DOF submatrix
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            active_dof = torch.nonzero(mask_dof, as_tuple=False).flatten()
            H = active_square(H, active_dof)
            del active_dof
            Q, info = active_tr_basis(
                coords_bohr_t.to(dtype=dtype, device=device),
                masses_au_t.to(dtype=dtype, device=device),
                active_idx,
                mode=tr_projection,
            )
            if compact:
                H, lift = compact_project_hessian(H, Q)
            else:
                project_hessian_inplace(H, Q)
                lift = None
            if projection_info is not None:
                projection_info.clear()
                projection_info.update(info.as_dict())
            # Bounded-peak symmetrization (helper writes both triangles).
            from pdb2reaction.core.utils import symmetrize_inplace
            symmetrize_inplace(H)
            del Q, mask_dof, active_idx, frozen
        else:
            # Full DOF: mass-weight + TR projection in-place
            H = _mass_weighted_hessian(H, masses_au_t.to(dtype=dtype, device=device))
            coords_bohr_t = coords_bohr_t.to(dtype=dtype, device=device)
            masses_au_t = masses_au_t.to(dtype=dtype, device=device)
            Q, info = active_tr_basis(
                coords_bohr_t,
                masses_au_t,
                list(range(int(N))),
                mode=tr_projection,
            )
            if compact:
                H, lift = compact_project_hessian(H, Q)
            else:
                project_hessian_inplace(H, Q)
                lift = None
            if projection_info is not None:
                projection_info.clear()
                projection_info.update(info.as_dict())
            # Bounded-peak symmetrization (helper writes both triangles).
            from pdb2reaction.core.utils import symmetrize_inplace
            symmetrize_inplace(H)
            del Q
        if torch.cuda.is_available() and device.type == "cuda":
            torch.cuda.empty_cache()
        return (H, lift) if compact else H


def _omega2_to_freq_cm(omega2: torch.Tensor) -> float:
    """
    Convert a (possibly signed) omega^2 eigenvalue to cm^-1, matching freq.py.
    """
    s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
    hnu = s_new * torch.sqrt(torch.abs(omega2))
    hnu = torch.where(omega2 < 0, -hnu, hnu)
    return float((hnu / units.invcm).item())

def _mode_direction_by_root(H: torch.Tensor,
                            coords_bohr_t: torch.Tensor,
                            masses_au_t: torch.Tensor,
                            root: int = 0,
                            freeze_idx: Optional[List[int]] = None,
                            tr_projection: str = "constrained",
                            projection_info: Optional[dict] = None) -> Tuple[np.ndarray, float]:
    """
    Get the eigenvector (Cartesian space) corresponding to the `root`-th most negative
    eigenvalue (root=0: most negative) of the mass-weighted, TR-projected Hessian.
    PHVA (active-subspace) is applied if freeze_idx is provided: frozen DOFs are zero.
    For root==0 the routine prefers LOBPCG and falls back to eigh as needed.
    """
    with torch.no_grad():
        # In-place: mass weight + (active-subspace) TR projection
        Hmw_proj, lift = _mw_projected_hessian_inplace(
            H,
            coords_bohr_t,
            masses_au_t,
            freeze_idx=freeze_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )
        if Hmw_proj.shape[0] == 0:
            raise RuntimeError(
                "No Dimer orientation remains after rigid-null projection."
            )

        # Solve eigenproblem in the (possibly reduced) space
        if int(root) == 0:
            try:
                w, v_mw_sub = torch.lobpcg(Hmw_proj, k=1, largest=False)
                u_mw_reduced = v_mw_sub[:, 0]
                omega2 = w[0]
            except Exception:
                evals_f, evecs_f = torch.linalg.eigh(Hmw_proj)
                pick = int(torch.argmin(evals_f).item())
                u_mw_reduced = evecs_f[:, pick]
                omega2 = evals_f[pick]
                del evals_f, evecs_f
        else:
            evals, evecs_mw = torch.linalg.eigh(Hmw_proj)  # ascending
            neg = (evals < 0)
            neg_inds = torch.nonzero(neg, as_tuple=False).view(-1)
            if neg_inds.numel() == 0:
                pick = int(torch.argmin(evals).item())
            else:
                k = max(0, min(int(root), neg_inds.numel() - 1))
                pick = int(neg_inds[k].item())
            u_mw_reduced = evecs_mw[:, pick]
            omega2 = evals[pick]
            del evals, evecs_mw

        u_mw_sub = (
            lift.T @ u_mw_reduced if lift is not None else u_mw_reduced
        )

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
        m3 = torch.repeat_interleave(masses_amu_t, 3).clamp(min=1e-10)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        v = inv_sqrt_m * u_mw
        v = v / torch.linalg.norm(v)
        mode = v.reshape(-1, 3).detach().cpu().numpy()
        freq_cm = _omega2_to_freq_cm(omega2)

        del masses_amu_t, m3, inv_sqrt_m, v, u_mw, u_mw_sub, u_mw_reduced, omega2
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return mode, freq_cm


def _calc_gradient(geom, calc_kwargs: dict) -> np.ndarray:
    """
    Return true Cartesian gradient (shape 3N,) in Hartree/Bohr.
    """
    calc = create_calculator(**calc_kwargs)
    geom.set_calculator(calc)
    echo_resolved_device()
    g = np.array(geom.gradient, dtype=float).reshape(-1)
    geom.set_calculator(None)
    return g



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
    idx = np.flatnonzero(np.asarray(mask_dof, dtype=bool))
    return active_square(H_full, idx)


def _representative_atoms_for_mode(mode_vec: torch.Tensor, flatten_k: int) -> np.ndarray:
    """
    Return indices of atoms with the largest displacements for a given mode.
    """
    vec = mode_vec.view(-1, 3)
    norms = torch.linalg.norm(vec, dim=1)
    k = min(int(flatten_k), vec.shape[0])
    if k <= 0:
        return np.zeros(0, dtype=int)
    topk = torch.topk(norms, k=k, largest=True)
    return topk.indices.detach().cpu().numpy()


def _select_flatten_targets_for_geom(
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    coords_bohr: np.ndarray,
    neg_freq_thresh_cm: float,
    root: int,
    flatten_sep_cutoff: float,
    flatten_k: int,
    primary_idx: Optional[int] = None,
) -> List[int]:
    """
    Select a subset of imaginary modes to flatten for a geometry.
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) <= 1:
        return []

    order = np.argsort(freqs_cm[neg_idx_all])
    sorted_neg = neg_idx_all[order]
    if primary_idx is None or int(primary_idx) not in set(map(int, sorted_neg)):
        root_clamped = max(0, min(int(root), len(order) - 1))
        primary_idx = int(sorted_neg[root_clamped])
    else:
        primary_idx = int(primary_idx)
    candidates = [int(i) for i in sorted_neg if int(i) != int(primary_idx)]
    if not candidates:
        return []

    coords_ang = torch.as_tensor(
        coords_bohr.reshape(-1, 3) * BOHR2ANG,
        dtype=modes.dtype,
        device=modes.device,
    )

    targets: List[int] = []
    reps_list: List[np.ndarray] = []

    for idx in candidates:
        rep = _representative_atoms_for_mode(modes[idx], flatten_k)
        if rep.size == 0:
            continue
        rep_coords = coords_ang[rep]
        if not reps_list:
            targets.append(idx)
            reps_list.append(rep)
            continue

        accept = True
        for prev_rep in reps_list:
            prev_coords = coords_ang[prev_rep]
            dmat = torch.cdist(rep_coords, prev_coords)
            min_dist = float(torch.min(dmat).item())
            if min_dist < float(flatten_sep_cutoff):
                accept = False
                break
        if accept:
            targets.append(idx)
            reps_list.append(rep)

    return targets


def _flatten_once_with_modes_for_geom(
    geom,
    masses_amu: np.ndarray,
    calc_kwargs: dict,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    flatten_amp_ang: float,
    flatten_sep_cutoff: float,
    flatten_k: int,
    root: int,
    reference_mode: Optional[np.ndarray] = None,
) -> bool:
    """
    Flatten extra imaginary modes for a geometry (single pass).
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) <= 1:
        return False

    m3 = np.repeat(masses_amu, 3).reshape(-1, 3)
    primary_idx = None
    if reference_mode is not None:
        reference = np.asarray(reference_mode, dtype=float).reshape(-1)
        reference /= np.linalg.norm(reference)
        overlaps = []
        for idx in neg_idx_all:
            v_mw = modes[int(idx)].detach().cpu().numpy().reshape(-1, 3)
            v_cart = (v_mw / np.sqrt(m3)).reshape(-1)
            v_cart /= np.linalg.norm(v_cart)
            overlaps.append(abs(float(np.dot(reference, v_cart))))
        primary_idx = int(neg_idx_all[int(np.argmax(overlaps))])

    targets = _select_flatten_targets_for_geom(
        freqs_cm,
        modes,
        geom.cart_coords,
        neg_freq_thresh_cm,
        root,
        flatten_sep_cutoff,
        flatten_k,
        primary_idx=primary_idx,
    )
    if not targets:
        return False

    amp_bohr = float(flatten_amp_ang) / BOHR2ANG
    E_ref = _calc_energy(geom, calc_kwargs)

    for idx in targets:
        v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
        v_cart = v_mw / np.sqrt(m3)
        v_cart /= np.linalg.norm(v_cart)

        # ``modes`` are eigenvectors of the mass-weighted Hessian. Dividing by
        # sqrt(mass) above already converts the selected eigenvector to its
        # Cartesian direction. Applying a second atom-wise mass factor here
        # rotates that direction (most strongly toward H atoms), so the trial
        # is no longer guaranteed to follow the negative-curvature mode.
        disp = amp_bohr * v_cart
        ref = geom.cart_coords.reshape(-1, 3)

        plus = ref + disp
        minus = ref - disp

        _set_cartesian_flatten_coords(geom, plus)
        E_plus = _calc_energy(geom, calc_kwargs)

        _set_cartesian_flatten_coords(geom, minus)
        E_minus = _calc_energy(geom, calc_kwargs)

        use_plus = E_plus <= E_minus
        _set_cartesian_flatten_coords(geom, plus if use_plus else minus)
        E_keep = E_plus if use_plus else E_minus
        delta_e = E_keep - E_ref
        click.echo(
            f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
            f"E_disp={E_keep:.8f} Ha \u0394E={delta_e:+.8f} Ha"
        )

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True


def _mw_tr_project_active_inplace(H: torch.Tensor,
                                  coords_all_t: torch.Tensor,
                                  masses_all_au_t: torch.Tensor,
                                  active_idx: Sequence[int],
                                  tr_projection: str = "constrained",
                                  projection_info: Optional[dict] = None,
                                  compact: bool = False):
    """
    Mass-weight & project TR in the *active* subspace (in-place; no explicit symmetrization).
    """
    with torch.no_grad():
        # mass-weight
        masses_act_au_t = masses_all_au_t[list(active_idx)]
        masses_amu_t = (masses_act_au_t / AMU2AU).to(dtype=H.dtype, device=H.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3).clamp(min=1e-10)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)
        H.mul_(inv_sqrt_m_row)
        H.mul_(inv_sqrt_m_col)
        Q, info = active_tr_basis(
            coords_all_t.to(dtype=H.dtype, device=H.device),
            masses_all_au_t.to(dtype=H.dtype, device=H.device),
            active_idx,
            mode=tr_projection,
        )
        if compact:
            H, lift = compact_project_hessian(H, Q)
        else:
            project_hessian_inplace(H, Q)
            lift = None
        if projection_info is not None:
            projection_info.clear()
            projection_info.update(info.as_dict())
        del masses_act_au_t, masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row, Q
        return (H, lift) if compact else H




def _mode_direction_by_root_from_Hact(H: torch.Tensor,
                                      coords_bohr: np.ndarray,
                                      atomic_numbers: List[int],
                                      masses_au_t: torch.Tensor,
                                      active_idx: List[int],
                                      device: torch.device,
                                      root: int = 0,
                                      tr_projection: str = "constrained",
                                      projection_info: Optional[dict] = None) -> Tuple[np.ndarray, float]:
    """
    TS direction from the *active* Hessian block. Mass-weighting/TR are done in the
    active space. Result is embedded back to full 3N in Cartesian space.
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_all = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H.dtype, device=device)
        masses_act_au = masses_au_t[active_idx].to(device=device, dtype=H.dtype)
        # mass-weight + TR in active space
        Hmw, lift = _mw_tr_project_active_inplace(
            H.clone(),
            coords_all,
            masses_au_t,
            active_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )
        if Hmw.shape[0] == 0:
            raise RuntimeError(
                "No Dimer orientation remains after rigid-null projection."
            )
        # Bounded-peak symmetrization (helper writes both triangles).
        from pdb2reaction.core.utils import symmetrize_inplace
        symmetrize_inplace(Hmw)

        # eigenvector for requested root
        if int(root) == 0:
            try:
                w, V = torch.lobpcg(Hmw, k=1, largest=False)
                u_mw_reduced = V[:, 0]
                omega2 = w[0]
            except Exception:
                vals, vecs = torch.linalg.eigh(Hmw)
                pick = int(torch.argmin(vals).item())
                u_mw_reduced = vecs[:, pick]
                omega2 = vals[pick]
                del vals, vecs
        else:
            vals, vecs = torch.linalg.eigh(Hmw)
            neg = (vals < 0)
            neg_inds = torch.nonzero(neg, as_tuple=False).view(-1)
            if neg_inds.numel() == 0:
                pick = int(torch.argmin(vals).item())
            else:
                k = max(0, min(int(root), neg_inds.numel() - 1))
                pick = int(neg_inds[k].item())
            u_mw_reduced = vecs[:, pick]
            omega2 = vals[pick]
            del vals, vecs

        u_mw = lift.T @ u_mw_reduced if lift is not None else u_mw_reduced

        # Mass un-weight to Cartesian in the active space, then embed to full 3N
        masses_act_amu = (masses_act_au / AMU2AU).to(dtype=H.dtype, device=device)
        m3 = torch.repeat_interleave(masses_act_amu, 3)
        v_cart_act = u_mw / torch.sqrt(m3)
        v_cart_act = v_cart_act / torch.linalg.norm(v_cart_act)

        full = torch.zeros(3 * N, dtype=H.dtype, device=device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))
        mask_t = torch.as_tensor(mask_dof, dtype=torch.bool, device=device)
        full[mask_t] = v_cart_act
        mode = full.reshape(-1, 3).detach().cpu().numpy()
        freq_cm = _omega2_to_freq_cm(omega2)

        del coords_all, masses_act_au, masses_act_amu, m3, v_cart_act
        del full, mask_t, Hmw, u_mw, u_mw_reduced, omega2, lift
        if torch.cuda.is_available() and H.is_cuda:
            torch.cuda.empty_cache()
        return mode, freq_cm


def _bofill_update_active(H: torch.Tensor,
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
    device = H.device
    dtype = H.dtype

    # as torch vectors
    d = torch.as_tensor(delta_act, dtype=dtype, device=device).reshape(-1)
    g0 = torch.as_tensor(g_old_act, dtype=dtype, device=device).reshape(-1)
    g1 = torch.as_tensor(g_new_act, dtype=dtype, device=device).reshape(-1)
    y = g1 - g0

    # Use current symmetric H for matvec
    Hd = H @ d
    xi = y - Hd

    d_dot_xi = torch.dot(d, xi)
    d_norm2 = torch.dot(d, d)
    xi_norm2 = torch.dot(xi, xi)

    # guards
    if torch.abs(d_dot_xi) > eps:
        denom_ms = d_dot_xi
    else:
        sign = torch.sign(d_dot_xi)
        denom_ms = (sign if sign != 0 else torch.tensor(1.0, device=device)) * eps
    denom_psb_d4 = d_norm2 * d_norm2 if d_norm2 > eps else eps
    denom_psb_d2 = d_norm2 if d_norm2 > eps else eps
    denom_phi = d_norm2 * xi_norm2 if (d_norm2 > eps and xi_norm2 > eps) else (1.0)

    phi = 1.0 - (d_dot_xi * d_dot_xi) / denom_phi
    phi = torch.clamp(phi, 0.0, 1.0)

    # coefficients for rank updates
    alpha = (1.0 - phi) / denom_ms                      # for xi xi^T
    beta = -phi * (d_dot_xi / denom_psb_d4)             # for d d^T
    gamma = phi / denom_psb_d2                          # for d xi^T + xi d^T

    n = H.shape[0]
    iu0, iu1 = torch.triu_indices(n, n, device=device)
    is_diag = (iu0 == iu1)
    off = ~is_diag

    # Diagonal contributions (i == j): alpha*xi_i^2 + beta*d_i^2 + 2*gamma*d_i*xi_i
    if is_diag.any():
        idx = iu0[is_diag]
        diag_inc = (alpha * xi[idx] * xi[idx]
                    + beta * d[idx] * d[idx]
                    + 2.0 * gamma * d[idx] * xi[idx])
        # Explicit read-modify-write preserves duplicate-index accumulation;
        # indexed ``+=`` does not provide those semantics in PyTorch.
        # CHEMISTRY-RULE:7 Bofill advanced-indexing (assignment, NOT in-place +=).
        # NOTE: use assignment so that advanced indexing actually updates H
        H[idx, idx] = H[idx, idx] + diag_inc

    # Off-diagonal (i < j): symmetric update
    if off.any():
        i = iu0[off]
        j = iu1[off]
        inc = (alpha * xi[i] * xi[j]
               + beta * d[i] * d[j]
               + gamma * (d[i] * xi[j] + xi[i] * d[j]))
        H[i, j] = H[i, j] + inc
        H[j, i] = H[j, i] + inc

    return H


def _imaginary_mode_indices_and_values(
    freqs_cm: np.ndarray,
    neg_freq_thresh_cm: float,
) -> Tuple[np.ndarray, List[float]]:
    """
    Return indices and values of frequencies regarded as imaginary under threshold.
    """
    thresh = abs(float(neg_freq_thresh_cm))
    neg_idx = np.where(freqs_cm < -thresh)[0]
    ims = [float(freqs_cm[i]) for i in neg_idx]
    return neg_idx, ims


def _warn_if_leading_imaginary_mode_is_soft(ims: List[float]) -> None:
    """
    Warn when the imaginary mode that certifies the saddle is very soft.

    Certification counts imaginary modes; it does not assess their character.
    The warning does not change the terminal status or counting rule.
    """
    if not ims:
        return
    leading = min(float(x) for x in ims)
    if abs(leading) >= TS_IMAG_SOFT_WARN_CM:
        return
    emit(
        f"[tsopt] WARNING: the leading imaginary mode is {leading:.2f} cm^-1, "
        f"below {TS_IMAG_SOFT_WARN_CM:.0f} cm^-1. Visualize the mode and "
        f"confirm IRC connectivity before treating this as a transition state.",
        narrative=True,
    )


def _echo_imaginary_modes(
    freqs_cm: np.ndarray,
    neg_freq_thresh_cm: float,
) -> np.ndarray:
    """
    Print and return imaginary-mode indices.
    """
    neg_idx, ims = _imaginary_mode_indices_and_values(freqs_cm, neg_freq_thresh_cm)
    emit(f"[Imaginary modes] n={len(neg_idx)}  ({ims})", narrative=True)
    _warn_if_leading_imaginary_mode_is_soft(ims)
    return neg_idx


def _write_all_imaginary_modes(
    geom,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    masses_amu: Sequence[float],
    vib_dir: Path,
    prepared_input: Optional["PreparedInputStructure"] = None,
    ref_pdb: Optional[Path] = None,
    amplitude_ang: float = 0.8,
    n_frames: int = 20,
) -> int:
    """
    Write all imaginary modes as trajectory (and optional PDB) files.
    """
    neg_idx, _ = _imaginary_mode_indices_and_values(freqs_cm, neg_freq_thresh_cm)
    if len(neg_idx) == 0:
        return 0

    vib_dir.mkdir(parents=True, exist_ok=True)
    order = neg_idx[np.argsort(freqs_cm[neg_idx])]

    masses_amu_t = torch.as_tensor(masses_amu, dtype=modes.dtype, device=modes.device)
    m3 = torch.repeat_interleave(masses_amu_t, 3)
    freq_name_counts: Dict[str, int] = {}
    n_written = 0

    for mode_idx in order:
        i = int(mode_idx)
        freq = float(freqs_cm[i])
        mode_mw = modes[i]
        v_cart = (mode_mw / torch.sqrt(m3)).detach().cpu().numpy()
        norm = np.linalg.norm(v_cart)
        if norm <= 0.0:
            click.echo(
                f"[Mode write] WARNING: Skip mode {i} ({freq:+.2f} cm^-1) due to near-zero norm.",
                err=True,
            )
            continue
        v_cart = v_cart / norm

        freq_key = f"{freq:+.2f}"
        freq_name_counts[freq_key] = freq_name_counts.get(freq_key, 0) + 1
        serial = freq_name_counts[freq_key]
        suffix = "" if serial == 1 else f"_{serial:02d}"
        stem = f"imag_{freq:+.2f}cm-1{suffix}"

        out_trj = vib_dir / f"{stem}_trj.xyz"
        out_pdb = vib_dir / f"{stem}.pdb"
        _write_mode_trj_and_pdb(
            geom,
            v_cart,
            out_trj,
            amplitude_ang=amplitude_ang,
            n_frames=n_frames,
            comment=f"imag {freq:+.2f} cm-1",
            ref_pdb=ref_pdb,
            write_pdb=ref_pdb is not None,
            prepared_input=prepared_input,
            out_pdb=out_pdb,
        )
        n_written += 1

    del masses_amu_t, m3
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return n_written


#                   HessianDimer (Hessian Guided Dimer)

def _tsopt_terminal_status(
    optimizer: Any,
    *,
    saddle_verified: bool,
    projection_certifiable: bool = True,
) -> str:
    """Compose a TS optimizer's public status.

    Precedence: a stall wins (energy-plateau outcome is
    never a converged saddle); otherwise a genuine first-order saddle
    (``is_converged`` and ``saddle_verified``) is ``converged``; anything else
    is ``not_converged``.  n_imag / stop_reason are recorded separately so a
    stall does not hide saddle-order evidence.
    """
    if getattr(optimizer, "is_stalled", False):
        return "stalled"
    if (
        getattr(optimizer, "is_converged", False)
        and saddle_verified
        and projection_certifiable
    ):
        return "converged"
    return "not_converged"


def _finalize_dimer_saddle_status(
    runner: Any,
    freqs_cm: np.ndarray,
    neg_freq_thresh_cm: float,
) -> np.ndarray:
    """Record the final exact-Hessian verdict on a dimer runner."""

    neg_idx, imaginary = _imaginary_mode_indices_and_values(
        freqs_cm, neg_freq_thresh_cm
    )
    certifiable = allows_saddle_certification(
        getattr(runner, "tr_projection", "constrained"),
        getattr(runner, "freeze_atoms", ()),
    )
    runner.n_imaginary_modes = len(neg_idx)
    runner.imaginary_frequencies_cm = imaginary
    runner.saddle_order_verified = bool(certifiable and len(neg_idx) == 1)
    if not runner.saddle_order_verified:
        runner.is_converged = False
    if not certifiable:
        runner.stop_reason = (
            "legacy-active projection is comparison-only for frozen systems"
        )
        click.echo(f"[tsopt] WARNING: {runner.stop_reason}.", err=True)
    return neg_idx


class HessianDimer:
    """
    Hessian Guided Dimer TS search with periodic Hessian updates.

    Extensions in this implementation:
      - `root` parameter: choose which imaginary mode to follow (0 = most negative).
      - Pass-through kwargs: `dimer_kwargs` and `lbfgs_kwargs` to tune internals.
      - Hard cap on total LBFGS steps across segments: `max_total_cycles`.
      - PHVA (active DOF subspace) + TR projection for mode picking,
        respecting `freeze_atoms`. For `root == 0` the implementation prefers LOBPCG.
      - The flatten loop can optionally use a Bofill-updated active Hessian block in the
        active DOF subspace (controlled by `flatten_loop_bofill`), but Bofill is applied
        *only* for the flatten displacements; after each dimer segment in the flatten
        loop, a fresh exact Hessian is recomputed.
      - Only **spatially separated** extra imaginary modes (based on representative
        atoms and a distance cutoff) are used for flattening to avoid overly large
        displacements when clustered imaginary modes are present.
      - Calculator kwargs accept `freeze_atoms` and `hessian_calc_mode` and
        default to returning a partial (active) Hessian when applicable.

    """

    def __init__(self,
                 fn: str,
                 out_dir: str = OUT_DIR_TSOPT,
                 thresh_loose: str = "gau_loose",
                 thresh: str = "baker",
                 update_interval_hessian: int = 500,
                 # neg_freq_thresh_cm: a mode counts as imaginary iff freqs_cm < -5.0. It is
                 # SUBORDINATE to the ~5.14 cm⁻¹ tol floor in _frequencies_cm_and_modes
                 # (|ω²| > 1e-6 au drops near-zero modes first), so survivors already exceed
                 # ~5.14 cm⁻¹ and 5.0 never fires independently. Do NOT raise it above ~5.14
                 # expecting a looser cutoff — that would discard genuine soft modes.
                 neg_freq_thresh_cm: float = 5.0,
                 flatten_amp_ang: float = 0.10,
                 flatten_max_iter: int = 0,
                 mem: int = 100000,
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
                 flatten_sep_cutoff: float = 0.0,
                 flatten_k: int = 10,
                 flatten_loop_bofill: bool = False,
                 #
                 # Propagate geometry kwargs so freeze-links and YAML geometry overrides
                 # also apply in grad mode.
                 geom_kwargs: Optional[Dict[str, Any]] = None,
                 prepared_input: Optional["PreparedInputStructure"] = None,
                 ) -> None:

        self.fn = fn
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.vib_dir = self.out_dir / "vib"
        self.vib_dir.mkdir(parents=True, exist_ok=True)
        # Use prepared_input.source_path if available (may be overridden by --ref-pdb),
        # otherwise fall back to fn directly.
        _src = (prepared_input.source_path if prepared_input is not None else Path(fn))
        self.ref_pdb: Optional[Path] = _src if _src.suffix.lower() == ".pdb" else None
        self.prepared_input = prepared_input

        self.thresh_loose = thresh_loose
        self.thresh = thresh
        self.update_interval_hessian = int(update_interval_hessian)
        self.neg_freq_thresh_cm = float(neg_freq_thresh_cm)
        self.flatten_amp_ang = float(flatten_amp_ang)
        self.flatten_max_iter = int(flatten_max_iter)
        self.mem = int(mem)
        self.root = int(root)
        self.dimer_kwargs = dict(dimer_kwargs or {})
        self.lbfgs_kwargs = dict(lbfgs_kwargs or {})
        self.max_total_cycles = int(max_total_cycles)

        # multi-mode flatten controls
        self.flatten_sep_cutoff = float(flatten_sep_cutoff)
        self.flatten_k = int(flatten_k)
        self.flatten_loop_bofill = bool(flatten_loop_bofill)

        # Track total cycles globally across all loops/segments
        self._cycles_spent = 0

        # Real convergence state of the last optimization loop that ran.
        # True iff that loop terminated because the optimizer reported
        # convergence (not because the global cycle budget was exhausted).
        # Mirrors the RS-I-RFO branch's `last_optimizer.is_converged` so the
        # result.json status reflects whether the TS-opt actually converged.
        self.is_converged = False

        # Propagate the stall state from a child LBFGS whose
        # energy plateaued (energy stopped decreasing while its force/step
        # criteria stayed unmet).  A stall stops all later segments/loops and
        # is never reported as a converged TS.
        self.is_stalled = False
        self.stop_reason = ""
        self.saddle_order_verified = False
        self.n_imaginary_modes: Optional[int] = None
        self.imaginary_frequencies_cm: List[float] = []

        # UMA settings
        self.uma_kwargs = dict(charge=0, spin=1, model=DEFAULT_UMA_MODEL,
                               task_name="omol", device="auto") if uma_kwargs is None else dict(uma_kwargs)

        # Geometry & masses (use provided geom kwargs so freeze_atoms etc. apply)
        gkw = dict(geom_kwargs or {})
        coord_type = gkw.pop("coord_type", GEOM_KW_DEFAULT["coord_type"])
        self.geom = geom_loader(fn, coord_type=coord_type, **gkw)
        self.tr_projection = self.geom.tr_projection
        self.rigid_projection_info: Dict[str, Any] = {}
        self.masses_amu = _safe_masses_amu(self.geom.atomic_numbers)
        self.masses_au_t = torch.as_tensor(self.masses_amu * AMU2AU, dtype=torch.float32)

        # Preserve freeze list (for PHVA)
        self.freeze_atoms: List[int] = [int(i) for i in self.geom.freeze_atoms]

        # Device
        self.device = _torch_device(device)
        self.masses_au_t = self.masses_au_t.to(self.device)

        # temp file for Dimer orientation (N_raw)
        self.mode_path = self.out_dir / ".dimer_mode.dat"

        self.dump = bool(dump)
        self.optim_all_path = self.out_dir / "optimization_all_trj.xyz"

        self._raw_hessian_cache_cpu: Optional[torch.Tensor] = None
        self._raw_hessian_coords_cpu: Optional[np.ndarray] = None

    @property
    def termination_status(self) -> str:
        """Public terminal outcome: ``stalled`` > ``converged`` > ``not_converged``."""
        if self.is_stalled:
            return "stalled"
        if self.is_converged:
            return "converged"
        return "not_converged"

    def _cache_raw_hessian_cpu(self, H: torch.Tensor) -> None:
        """
        Cache the *raw* Hessian on CPU for the current geometry.
        IMPORTANT:
          - Must be called before any in-place mass-weight / TR projection happens.
          - This cache is CPU-only by design (no persistent GPU cache).
        """
        self._raw_hessian_cache_cpu = H.detach().cpu().clone()
        self._raw_hessian_coords_cpu = self.geom.cart_coords.copy()

    def _sync_geom_hessian_cache(self, H: torch.Tensor) -> None:
        """
        Share a full-space raw Hessian with Geometry-level cache so RS-I-RFO/RFO
        can reuse the same Hessian on unchanged coordinates.

        Uses the CPU cache if available to avoid an extra GPU clone.
        """
        try:
            n_dof = int(self.geom.cart_coords.size)
            if tuple(H.shape) != (n_dof, n_dof):
                return
            # Prefer the CPU cache to avoid duplicating GPU memory
            if self._raw_hessian_cache_cpu is not None:
                H_cache = self._raw_hessian_cache_cpu
            else:
                H_cache = H.detach().cpu().clone()
            if hasattr(self.geom, "within_partial_hessian"):
                self.geom.within_partial_hessian = None
            self.geom.cart_hessian = H_cache
        except Exception:
            # Cache hand-off is best-effort; optimization must continue even if it fails.
            pass

    def _reuse_cached_hessian(self) -> Optional[torch.Tensor]:
        """
        If the cached geometry matches the current geometry, return the cached *raw*
        Hessian transferred to self.device. No persistent GPU caching is performed.

        Note:
          - On CPU, `.to('cpu')` can share storage; we clone to protect cache against
            later in-place ops (e.g., mass-weight/TR projection).
        """
        if self._raw_hessian_cache_cpu is None or self._raw_hessian_coords_cpu is None:
            return None
        if not np.array_equal(self.geom.cart_coords, self._raw_hessian_coords_cpu):
            return None
        H_dev = self._raw_hessian_cache_cpu.to(self.device)
        if self.device.type == "cpu":
            H_dev = H_dev.clone()
        return H_dev

    def _calc_full_hessian_cached(self, allow_reuse: bool) -> torch.Tensor:
        """
        Compute an exact Hessian, caching a CPU copy of the *raw* Hessian.
        If allow_reuse is True and the current geometry matches the cached one,
        reuse the cached Hessian (CPU→device transfer) instead of recalculating.
        """
        if allow_reuse:
            cached = self._reuse_cached_hessian()
            if cached is not None:
                click.echo("[Hessian] Reusing cached raw Hessian (0-step convergence).")
                self._sync_geom_hessian_cache(cached)
                return cached
        H = _calc_full_hessian_torch(self.geom, self.uma_kwargs, self.device)
        self._cache_raw_hessian_cpu(H)
        self._sync_geom_hessian_cache(H)
        return H

    # ----- One dimer segment for up to n_steps; returns (steps_done, converged) -----
    def _dimer_segment(self, threshold: str, n_steps: int) -> Tuple[int, bool]:
        # Dimer calculator using current mode as initial N
        calc_sp = create_calculator(**self.uma_kwargs)

        n_atoms = len(self.geom.atomic_numbers)
        active_idx = _active_indices(n_atoms, self.freeze_atoms)
        rigid_masses = self.masses_au_t.detach().to(device="cpu", dtype=torch.float64)
        projection_mode = self.tr_projection
        rigid_basis, projection = full_cartesian_tr_basis(
            torch.as_tensor(self.geom.cart_coords.reshape(-1, 3), dtype=torch.float64),
            rigid_masses,
            active_idx,
            mode=projection_mode,
        )

        def rigid_basis_at(coords_flat: np.ndarray) -> np.ndarray:
            basis, _ = full_cartesian_tr_basis(
                torch.as_tensor(coords_flat.reshape(-1, 3), dtype=torch.float64),
                rigid_masses,
                active_idx,
                mode=projection_mode,
            )
            return basis.detach().cpu().numpy()

        self.rigid_projection_info.clear()
        self.rigid_projection_info.update(projection.as_dict())

        # Merge user dimer kwargs, while enforcing the geometry's exact Cartesian
        # constraints and its selected rigid-null treatment.
        dimer_kwargs = dict(self.dimer_kwargs)
        dimer_kwargs.update({
            "calculator": calc_sp,
            "N_raw": str(self.mode_path),
            "frozen_atoms": self.freeze_atoms,
            "rigid_basis": rigid_basis.detach().cpu().numpy(),
            "rigid_basis_getter": rigid_basis_at,
            "write_orientations": False,
            "seed": 0,
            "mem": self.mem,
        })
        dimer = Dimer(**dimer_kwargs)

        self.geom.set_calculator(dimer)
        echo_resolved_device()

        # LBFGS kwargs: enforce thresh/max_cycles/out_dir/dump; allow others
        lbfgs_kwargs = _force_ts_reject_uphill_off(self.lbfgs_kwargs)
        lbfgs_kwargs.update({
            "max_cycles": n_steps,
            "thresh": threshold,
            "out_dir": str(self.out_dir),
            "dump": self.dump,
        })
        opt = LBFGS(self.geom, **lbfgs_kwargs)
        opt.run()
        steps = opt.cur_cycle + 1
        converged = opt.is_converged
        # Propagate an energy-plateau stall from the child LBFGS. A
        # stalled child is not converged; the caller stops all later segments.
        if getattr(opt, "is_stalled", False):
            self.is_stalled = True
            self.stop_reason = getattr(opt, "stop_reason", "") or self.stop_reason
        self.geom.set_calculator(None)

        # Append to concatenated trajectory if dump enabled
        if self.dump:
            part_path = self.out_dir / "optimization_trj.xyz"
            if part_path.exists():
                with self.optim_all_path.open("a", encoding="utf-8") as f_all, \
                     part_path.open("r", encoding="utf-8") as f_part:
                    f_all.write(f_part.read())
        return steps, converged

    # ----- Loop dimer segments, updating mode from Hessian every interval -----
    def _dimer_loop(self, threshold: str) -> Tuple[int, bool, bool]:
        """
        Run multiple LBFGS segments separated by periodic Hessian-based mode updates.
        Consumes from the global cycle budget self.max_total_cycles.

        Returns:
          (steps_in_this_call, zero_step_converged, loop_converged)
        where `zero_step_converged` is True iff the loop terminated by convergence
        without changing the geometry (i.e., 0-step convergence; used as a
        Hessian-reuse hint), and `loop_converged` is True iff the loop terminated
        because the underlying optimizer reported convergence (as opposed to the
        global cycle budget being exhausted). `loop_converged` is the honest
        convergence signal threaded out to the result.json status.
        """
        steps_in_this_call = 0
        zero_step_converged = False
        loop_converged = False
        while True:
            remaining_global = max(0, self.max_total_cycles - self._cycles_spent)
            if remaining_global == 0:
                break
            steps_this = min(self.update_interval_hessian, remaining_global)
            coords_before = self.geom.cart_coords.copy()
            steps, ok = self._dimer_segment(threshold, steps_this)
            self._cycles_spent += steps
            steps_in_this_call += steps
            # A stalled child stops all further segments in this loop.
            if self.is_stalled:
                break
            if ok:
                loop_converged = True
                if np.array_equal(self.geom.cart_coords, coords_before):
                    zero_step_converged = True
                break
            # If budget exhausted after this segment, stop before doing a Hessian update
            if (self.max_total_cycles - self._cycles_spent) <= 0:
                break
            # Update mode from Hessian (respect freeze atoms via PHVA)
            H = self._calc_full_hessian_cached(allow_reuse=False)
            N = len(self.geom.atomic_numbers)
            coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
                                            dtype=H.dtype, device=H.device)
            # full vs active-block Hessian
            if H.size(0) == 3 * N:
                mode_xyz, mode_freq_cm = _mode_direction_by_root(
                    H, coords_bohr_t, self.masses_au_t,
                    root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            else:
                # partial (active) Hessian returned by UMA
                active_idx = _active_indices(N, self.freeze_atoms if len(self.freeze_atoms) > 0 else [])
                mode_xyz, mode_freq_cm = _mode_direction_by_root_from_Hact(
                    H, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            emit(f"[Dimer mode] root={self.root} freq={mode_freq_cm:+.2f} cm^-1", narrative=True)
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del H, coords_bohr_t, mode_xyz, mode_freq_cm
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        return steps_in_this_call, zero_step_converged, loop_converged

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
        Flatten using precomputed modes (mass-weighted, embedded).

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

        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        # ensure energy reference is set up
        E_ref = _calc_energy(self.geom, self.uma_kwargs)

        # work in Bohr coordinates
        for idx in targets:
            # mode is currently mass-weighted & embedded to 3N
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
            # convert to Cartesian (Å-scale amplitude, but coords are Bohr)
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)

            # v_cart is already the un-mass-weighted Cartesian normal-mode
            # direction. Preserve it; an additional per-atom mass factor
            # would rotate the flatten displacement away from the mode.
            disp = amp_bohr * v_cart  # Bohr
            ref = self.geom.cart_coords.reshape(-1, 3)

            plus = ref + disp
            minus = ref - disp

            _set_cartesian_flatten_coords(self.geom, plus)
            E_plus = _calc_energy(self.geom, self.uma_kwargs)

            _set_cartesian_flatten_coords(self.geom, minus)
            E_minus = _calc_energy(self.geom, self.uma_kwargs)

            # keep lower-energy side and continue from there
            use_plus = E_plus <= E_minus
            _set_cartesian_flatten_coords(
                self.geom, plus if use_plus else minus
            )
            E_keep = E_plus if use_plus else E_minus
            delta_e = E_keep - E_ref
            click.echo(
                f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
                f"E_disp={E_keep:.8f} Ha \u0394E={delta_e:+.8f} Ha"
            )

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
        H = self._calc_full_hessian_cached(allow_reuse=False)
        coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
                                        dtype=H.dtype, device=H.device)

        if H.size(0) == 3 * N:
            mode_xyz, mode_freq_cm = _mode_direction_by_root(
                H, coords_bohr_t, self.masses_au_t,
                root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                tr_projection=self.tr_projection,
                projection_info=self.rigid_projection_info,
            )
        else:
            click.echo("[CHECK] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
            mode_xyz, mode_freq_cm = _mode_direction_by_root_from_Hact(
                H, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root,
                tr_projection=self.tr_projection,
                projection_info=self.rigid_projection_info,
            )
        emit(f"[Dimer mode] root={self.root} freq={mode_freq_cm:+.2f} cm^-1", narrative=True)
        np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
        del mode_xyz, coords_bohr_t, H, mode_freq_cm
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # (2) Loose loop (or initial pass)
        if self.root != 0:
            click.echo("[WARNING] root != 0. Using this root only for the initial dimer loop.", err=True)
            click.echo(f"Dimer loop with initial direction from mode {self.root}...")
            self.root = 0
            self.thresh_loose = self.thresh
        else:
            click.echo("Loose dimer loop...")

        _steps_loose, zero_step_loose, conv_loose = self._dimer_loop(self.thresh_loose)
        self.is_converged = conv_loose

        # (3) Update mode & normal loop. A stalled loose loop stops all further
        # optimization work: skip the Hessian/mode update and the
        # normal + flatten loops so a stalled TS search is never retried.
        zero_step_normal = zero_step_loose
        if not self.is_stalled:
            H = self._calc_full_hessian_cached(allow_reuse=zero_step_loose)
            coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
                                            dtype=H.dtype, device=H.device)
            if H.size(0) == 3 * N:
                mode_xyz, mode_freq_cm = _mode_direction_by_root(
                    H, coords_bohr_t, self.masses_au_t,
                    root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            else:
                click.echo("[CHECK] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
                mode_xyz, mode_freq_cm = _mode_direction_by_root_from_Hact(
                    H, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            emit(f"[Dimer mode] root={self.root} freq={mode_freq_cm:+.2f} cm^-1", narrative=True)
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del mode_xyz, coords_bohr_t, H, mode_freq_cm
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            click.echo("Normal dimer loop...")
            _steps_normal, zero_step_normal, conv_normal = self._dimer_loop(self.thresh)
            self.is_converged = conv_normal

        # (4) Flatten loop — exact Hessian each iteration & optional Bofill update.
        # A stalled optimization never enters the flatten/retry loop.
        if self.flatten_max_iter > 0 and not self.is_stalled:
            if self.flatten_loop_bofill:
                click.echo("Flatten loop with Bofill-updated active Hessian (flatten displacements only)...")
            else:
                click.echo("Flatten loop without Bofill updates (flatten displacements only)...")

            # (4.1) Evaluate one exact Hessian at the loop start and prepare the active block
            H = self._calc_full_hessian_cached(allow_reuse=zero_step_normal)
            if H.size(0) == 3 * N:
                # full → extract active
                H = _extract_active_block(H, mask_dof)  # torch (3N_act,3N_act)
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            # Flatten iterations with Hessian updates
            for _it in range(self.flatten_max_iter):
                if (self.max_total_cycles - self._cycles_spent) <= 0:
                    break

                # (a) Frequencies & modes from the active Hessian
                freqs_cm, modes_embedded = _frequencies_cm_and_modes(
                    H.clone(),
                    self.geom.atomic_numbers,
                    self.geom.cart_coords.reshape(-1, 3),
                    self.device,
                    freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
                n_imag = int(np.sum(freqs_cm < -abs(self.neg_freq_thresh_cm)))
                ims = [float(x) for x in freqs_cm if x < -abs(self.neg_freq_thresh_cm)]
                emit(f"[Imaginary modes] n={n_imag}  ({ims})", narrative=True)
                _warn_if_leading_imaginary_mode_is_soft(ims)
                if n_imag <= 1:
                    break

                # (b) Flatten other imaginary modes
                if self.flatten_loop_bofill:
                    x_before_flat = self.geom.cart_coords.copy().reshape(-1)
                    g_before_flat = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)

                did_flatten = self._flatten_once_with_modes(freqs_cm, modes_embedded)
                if not did_flatten:
                    break

                if self.flatten_loop_bofill:
                    x_after_flat = self.geom.cart_coords.copy().reshape(-1)
                    g_after_flat = _calc_gradient(self.geom, self.uma_kwargs).reshape(-1)

                # (c) Bofill update using UMA gradients across the flatten displacement
                if self.flatten_loop_bofill:
                    delta_flat_full = x_after_flat - x_before_flat
                    delta_flat_act = delta_flat_full[mask_dof]
                    g_old_act = g_before_flat[mask_dof]
                    g_new_act = g_after_flat[mask_dof]
                    H = _bofill_update_active(H, delta_flat_act, g_new_act, g_old_act)

                # (d) Refresh dimer direction
                mode_xyz, mode_freq_cm = _mode_direction_by_root_from_Hact(
                    H, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
                emit(f"[Dimer mode] root={self.root} freq={mode_freq_cm:+.2f} cm^-1", narrative=True)
                np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")

                # (e) Re-optimize with Dimer (consumes global cycle budget)
                _steps_flat, zero_step_flat, conv_flat = self._dimer_loop(self.thresh)
                self.is_converged = conv_flat

                # A stall inside the flatten loop stops the remaining iterations
                # Do not keep retrying a stalled optimization.
                if self.is_stalled:
                    break

                # (f) After dimer optimization, recompute an exact Hessian for the next iteration
                H = self._calc_full_hessian_cached(allow_reuse=zero_step_flat)
                if H.size(0) == 3 * N:
                    H = _extract_active_block(H, mask_dof)
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

        # Surface non-convergence: if the last optimization loop exhausted the
        # global cycle budget without the optimizer reporting convergence, the
        # final geometry is NOT a converged TS. Emit a visible WARNING (the
        # result.json status is set from self.termination_status at the call
        # site). A stall is a distinct energy-plateau outcome.
        if self.is_stalled:
            click.echo(
                f"[tsopt] WARNING: TS optimization stalled (energy plateau): "
                f"{self.stop_reason}",
                err=True,
            )
        elif not self.is_converged:
            click.echo(
                f"[tsopt] WARNING: max cycles reached without convergence "
                f"(max_total_cycles={self.max_total_cycles}, "
                f"cycles_spent={self._cycles_spent}).",
                err=True,
            )

        # (5) Final outputs
        final_xyz = self.out_dir / "final_geometry.xyz"
        atoms_final = Atoms(self.geom.atoms, positions=(self.geom.cart_coords.reshape(-1, 3) * BOHR2ANG), pbc=False)
        write(final_xyz, atoms_final)

        # Final Hessian → imaginary mode animation
        H = self._calc_full_hessian_cached(allow_reuse=True)
        freqs_cm, modes = _frequencies_cm_and_modes(
            H,
            self.geom.atomic_numbers,
            self.geom.cart_coords.reshape(-1, 3),
            self.device,
            freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
            tr_projection=self.tr_projection,
            projection_info=self.rigid_projection_info,
        )

        del H
        neg_idx = _finalize_dimer_saddle_status(
            self, freqs_cm, self.neg_freq_thresh_cm
        )
        _echo_imaginary_modes(freqs_cm, self.neg_freq_thresh_cm)
        if len(neg_idx) == 0:
            click.echo("[INFO] No imaginary mode found at the end (ν_min = %.2f cm^-1)." % (freqs_cm.min(),))
            del modes
        else:
            _write_all_imaginary_modes(
                self.geom,
                freqs_cm,
                modes,
                self.neg_freq_thresh_cm,
                self.masses_amu,
                self.vib_dir,
                prepared_input=self.prepared_input,
                ref_pdb=self.ref_pdb,
                amplitude_ang=0.8,
                n_frames=20,
            )
            del modes

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        emit(f"[DONE] Saved final geometry → {final_xyz}", detail=True)
        click.echo(f"[DONE] Mode files → {self.vib_dir}")



# Geometry defaults
GEOM_KW = dict(GEOM_KW_DEFAULT)

CALC_KW = dict(UMA_CALC_KW)

# Optimizer base (common) — used by both RSIRFO and the inner LBFGS of Hessian Guided Dimer
OPT_BASE_KW_LOCAL = {
    **OPT_BASE_KW,
    "out_dir": OUT_DIR_TSOPT,
    "thresh": "baker",
    "check_eigval_structure": True,
}


def _force_ts_reject_uphill_off(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Return TS optimizer kwargs with physical-energy rejection disabled."""
    effective = dict(kwargs)
    effective["reject_uphill"] = False
    return effective


def _build_rsirfo_kwargs(
    opt_cfg: Dict[str, Any],
    rsirfo_cfg: Dict[str, Any],
    out_dir_path: Path,
) -> Dict[str, Any]:
    """Return the exact kwargs used for RSIRFO after normalization."""
    rs_args = dict(rsirfo_cfg)
    opt_base = dict(opt_cfg)
    opt_base["out_dir"] = str(out_dir_path)
    rs_args["out_dir"] = str(out_dir_path)

    roots = rs_args.get("roots", None)
    root_single = rs_args.get("root", None)
    if roots is None and root_single is not None:
        roots = [int(root_single)]
    if roots is None:
        roots = [0]
    rs_args["roots"] = [int(x) for x in roots]
    rs_args.pop("root", None)

    # Keep top-level opt knobs (max_cycles/dump/thresh/...) authoritative unless
    # rsirfo.* explicitly overrides them to a non-default value.
    for k in list(rs_args.keys()):
        if k in opt_base and k in RSIRFO_KW and rs_args[k] == RSIRFO_KW[k]:
            rs_args.pop(k, None)

    rsirfo_kwargs = {**opt_base, **rs_args}

    # RSIRFO ignores these DIIS-related knobs; drop them for clarity.
    for diis_kw in ("gediis", "gdiis", "gdiis_thresh", "gediis_thresh", "gdiis_test_direction", "adapt_step_func"):
        rsirfo_kwargs.pop(diis_kw, None)

    return _force_ts_reject_uphill_off(rsirfo_kwargs)


def _validate_reference_mode_optimizer(
    mode: str,
    reference_mode_path: Optional[Path],
) -> None:
    """Reject path-mode guidance for optimizers that do not consume it."""
    if reference_mode_path is not None and mode == "dimer":
        raise click.BadParameter(
            "--ref-mode requires a Hessian TS optimizer; use "
            "--opt-mode hess, rsirfo, rsprfo, or trim.",
            param_hint="--ref-mode",
        )



@click.command(
    help="Transition state optimization (Dimer or RS-P-RFO).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .cif, .mmcif, .xyz, .gjf, _trj.xyz, ...).",
)
@click.option(
    "--ref-mode",
    "reference_mode_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "Advanced/internal path-mode hint for Hessian TS recovery "
        "(.npy or whitespace Cartesian 3N text). The all workflow supplies "
        "this automatically from its MEP; ordinary standalone tsopt runs "
        "normally omit it."
    ),
)
@click.option(
    "--workers",
    type=int,
    default=CALC_KW["workers"],
    show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel MLIP predictor (workers>1).",
)
@click.option(
    "--freeze-links/--no-freeze-links",
    "freeze_links",
    default=True,
    show_default=True,
    help="Freeze parent atoms of cap hydrogens (PDB/mmCIF input or XYZ/GJF with --ref-pdb).",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--tr-projection",
    type=click.Choice(["constrained", "legacy-active"], case_sensitive=False),
    default=GEOM_KW["tr_projection"],
    show_default=True,
    help=(
        "Rigid translation/rotation treatment for Cartesian PHVA. "
        "'constrained' removes only full-system rigid motions compatible with "
        "the frozen atoms; 'legacy-active' is deprecated comparison-only "
        "behavior and must not be used for pass/HOSP transition-state certification."
    ),
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/CIF/GJF companions based on the input format.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB/mmCIF topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Maximum number of optimization cycles.")
@click.option(
    "--flatten/--no-flatten",
    "flatten",
    default=False,
    show_default=True,
    help="Enable the extra-imaginary-mode flattening loop (grad: dimer loop, hess: post-RS-P-RFO).",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess", "dimer", "rsirfo", "trim", "rsprfo"], case_sensitive=False),
    default="hess",
    show_default=True,
    help=(
        "TS optimizer: 'grad'/'dimer' → Hessian Guided Dimer; "
        "'hess'/'rsprfo' → RS-P-RFO (Banerjee, default); 'trim' → TRIM (Helgaker); 'rsirfo' → RS-I-RFO."
    ),
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write the per-cycle optimization trajectory ('optimization_trj.xyz' for rsirfo/hess, 'optimization_all_trj.xyz' for grad/dimer) in the output directory.",
)
@click.option("-o", "--out-dir", type=str, default=OUT_DIR_TSOPT, show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default=False,
    help=(
        "Convergence preset for the active optimizer "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'baker' when not provided."
    ),
)
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Base YAML configuration file applied before explicit CLI options.",
)
@click.option(
    "--show-config/--no-show-config",
    "show_config",
    default=False,
    show_default=True,
    help="Print resolved configuration and continue execution.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running TS optimization.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="Choose MLIP Hessian evaluation mode (used unless YAML sets calc.hessian_calc_mode). Defaults to 'FiniteDifference'.",
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              show_default=True, help="xTB solvent model.")
@add_ml_charge_spin_options()
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@add_coord_type_option()
@add_print_every_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    reference_mode_path: Optional[Path],
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links: bool,
    freeze_atoms_text: Optional[str],
    tr_projection: str,
    convert_files: bool,
    ref_pdb: Optional[Path],
    max_cycles: int,
    flatten: bool,
    opt_mode: str,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    hessian_calc_mode: Optional[str],
    backend: str,
    solvent: str,
    solvent_model: str,
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
    cli_coord_type: Optional[str],
    print_every: Optional[int],
) -> None:
    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )
    from pdb2reaction.core.utils import resolve_configured_charge_spin
    charge, spin = resolve_configured_charge_spin(
        merged_yaml_cfg, charge=charge, spin=spin, ligand_charge=ligand_charge,
    )

    set_convert_file_enabled(convert_files)
    with prepared_cli_input(
        input_path,
        ref_pdb=ref_pdb,
        charge=charge,
        spin=spin,
        ligand_charge=ligand_charge,
        prefix="[tsopt]",
    ) as (prepared_input, resolved_charge, resolved_spin):
        validate_charge_spin_for_prepared(prepared_input, resolved_charge, resolved_spin)
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path
        time_start = time.perf_counter()

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg = dict(OPT_BASE_KW_LOCAL)
        import copy
        simple_cfg = copy.deepcopy(HESSIAN_DIMER_CLI_KW)
        # tsopt default: keep flatten loop off unless enabled via config or explicit --flatten.
        simple_cfg["flatten_max_iter"] = 0
        rsirfo_cfg = dict(RSIRFO_KW)

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (opt_cfg, (("opt",),)),
                (simple_cfg, (("hessian_dimer",),)),
                (rsirfo_cfg, (("rsirfo",),)),
            ],
        )

        # resolved_charge / resolved_spin already incorporate explicit CLI
        # -q/-m, gjf metadata fill, and --ligand-charge derivation. Assign
        # directly; an earlier .get("charge", resolved) idiom silently
        # returned the UMA_CALC_KW default 0 when -q was not passed,
        # dropping the derived value.
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        if cli_param_overridden(ctx, "workers"):
            calc_cfg["workers"] = int(workers)
        if cli_param_overridden(ctx, "workers_per_node"):
            calc_cfg["workers_per_node"] = int(workers_per_node)
        if cli_param_overridden(ctx, "backend"):
            calc_cfg["backend"] = backend
        if cli_param_overridden(ctx, "solvent"):
            calc_cfg["solvent"] = solvent
        if cli_param_overridden(ctx, "solvent_model"):
            calc_cfg["solvent_model"] = solvent_model
        from pdb2reaction.backends import apply_effective_precision
        apply_effective_precision(calc_cfg, precision)
        from pdb2reaction.backends import apply_backend_model_to_calc_cfg
        # Unconditional: also pops a raw backend_model token from a --config YAML
        # (the helper no-ops when neither the CLI arg nor the YAML names one).
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        from pdb2reaction.backends import apply_calc_file_to_calc_cfg
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        apply_backend_defaults(calc_cfg)
        if cli_param_overridden(ctx, "max_cycles"):
            opt_cfg["max_cycles"] = int(max_cycles)
        if cli_param_overridden(ctx, "dump"):
            opt_cfg["dump"] = bool(dump)
        if cli_param_overridden(ctx, "out_dir"):
            opt_cfg["out_dir"] = out_dir
        if cli_param_overridden(ctx, "thresh") and thresh is not None:
            opt_cfg["thresh"] = str(thresh)
            simple_cfg["thresh"] = str(thresh)
            rsirfo_cfg["thresh"] = str(thresh)
        if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
            geom_cfg["coord_type"] = str(cli_coord_type).lower()
        if cli_param_overridden(ctx, "tr_projection"):
            geom_cfg["tr_projection"] = str(tr_projection).lower()
        if cli_param_overridden(ctx, "hessian_calc_mode") and hessian_calc_mode is not None:
            calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
        # CLI knobs → cfg. Downstream plumbing into rsirfo_cfg / simple_cfg
        # below already understands these keys.
        if cli_param_overridden(ctx, "print_every") and print_every is not None:
            opt_cfg["print_every"] = int(print_every)
        if cli_param_overridden(ctx, "flatten"):
            if flatten:
                # --flatten explicitly enables flattening even when defaults/config disable it.
                default_flatten_iter = int(HESSIAN_DIMER_CLI_KW.get("flatten_max_iter", 0))
                if int(simple_cfg.get("flatten_max_iter", 0)) <= 0 and default_flatten_iter > 0:
                    simple_cfg["flatten_max_iter"] = default_flatten_iter
            else:
                simple_cfg["flatten_max_iter"] = 0

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (opt_cfg, (("opt",),)),
                (simple_cfg, (("hessian_dimer",),)),
                (rsirfo_cfg, (("rsirfo",),)),
            ],
        )

        # A TS search follows a saddle-search direction, so physical energy is
        # not required to decrease. Keep this invariant after every YAML merge.
        opt_cfg["reject_uphill"] = False
        rsirfo_cfg["reject_uphill"] = False
        simple_cfg["lbfgs"] = _force_ts_reject_uphill_off(
            simple_cfg.get("lbfgs", {})
        )

        try:
            geom_cfg["tr_projection"] = normalize_tr_projection_mode(
                geom_cfg.get("tr_projection")
            )
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc

        if "print_every" in opt_cfg:
            try:
                opt_print_every = int(opt_cfg["print_every"])
                if opt_print_every >= 1:
                    rsirfo_cfg.setdefault("print_every", opt_print_every)
                    simple_cfg.setdefault("lbfgs", {})
                    if isinstance(simple_cfg["lbfgs"], dict):
                        simple_cfg["lbfgs"].setdefault("print_every", opt_print_every)
            except (TypeError, ValueError):
                pass

        # Optimizer-side exact validation and the final reported PHVA must use
        # the same definition of an imaginary mode.
        rsirfo_cfg["saddle_imaginary_threshold_cm"] = float(
            simple_cfg.get("neg_freq_thresh_cm", 5.0)
        )

        # Convert 1-based YAML freeze_atoms to 0-based internal
        if geom_cfg.get("freeze_atoms"):
            geom_cfg["freeze_atoms"] = yaml_freeze_to_internal(geom_cfg["freeze_atoms"])
        # Merge CLI --freeze-atoms (already 0-based)
        try:
            freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)
        if freeze_atoms_cli:
            merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)
        # Normalize freeze_atoms and optionally add link-parent indices for PDB inputs
        resolve_freeze_atoms(geom_cfg, source_path, freeze_links)

        # Pass freeze_atoms from geom → calc (so UMA knows active DOF for FD Hessian)
        if "freeze_atoms" in geom_cfg:
            calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
        calc_cfg["return_partial_hessian"] = True

        kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=TSOPT_MODE_ALIASES,
            allowed_hint="grad|hess|dimer|rsirfo|trim|rsprfo",
        )
        _validate_reference_mode_optimizer(kind, reference_mode_path)
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()

        # Default-verbosity entry summary (skipped in child mode; under -v
        # the full pretty_block dump below restores everything).
        from pdb2reaction.core.utils import echo_run_summary
        _backend_disp = calc_cfg.get("backend", "?")
        _model = calc_cfg.get("model")
        _precision = calc_cfg.get("precision", "fp32")
        _max_cycles = opt_cfg.get("max_cycles", "?")
        echo_run_summary({
            "input": str(input_path),
            "backend": f"{_backend_disp} ({_model}, {_precision})" if _model else _backend_disp,
            "opt": f"{kind}, max_cycles={_max_cycles}",
            "out": str(out_dir_path),
        })

        # Pretty-print config summary
        click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
        click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
        echo_opt = {**opt_cfg, "out_dir": str(out_dir_path)}
        click.echo(pretty_block("opt", echo_opt))
        if kind == "dimer":
            sd_cfg_for_echo = dict(simple_cfg)
            sd_cfg_for_echo["dimer"] = dict(simple_cfg.get("dimer", {}))
            sd_cfg_for_echo["lbfgs"] = strip_inherited_keys(
                dict(simple_cfg.get("lbfgs", {})),
                echo_opt,
                mode="same",
            )
            click.echo(pretty_block("hessian_dimer", sd_cfg_for_echo))
        else:
            rsirfo_kwargs_for_echo = _build_rsirfo_kwargs(opt_cfg, rsirfo_cfg, out_dir_path)
            rsirfo_kwargs_for_echo = strip_inherited_keys(
                rsirfo_kwargs_for_echo,
                echo_opt,
                mode="same",
            )
            click.echo(pretty_block("rsirfo", rsirfo_kwargs_for_echo))

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override_yaml": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                force=True)
            )
        if dry_run:
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_geometry": str(geom_input_path),
                        "output_dir": str(out_dir_path),
                        "optimizer_mode": str(kind),
                        "convert_files": bool(convert_files),
                        "freeze_links": bool(freeze_links),
                        "tr_projection": geom_cfg["tr_projection"],
                        "flatten": int(simple_cfg.get("flatten_max_iter", 0)) > 0,
                        "will_run_tsopt": True,
                        "will_write_summary": True,
                    },
                )
            )
            click.echo("[dry-run] Validation complete. TS optimization execution was skipped.")
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)
        _tsopt_result_data: Optional[Dict[str, Any]] = None

        try:
            if kind == "dimer":
                # Hessian Guided Dimer runner construction
                calc_kwargs_for_dimer = dict(calc_cfg)
                runner = HessianDimer(
                    fn=str(geom_input_path),
                    out_dir=str(out_dir_path),
                    thresh_loose=simple_cfg.get("thresh_loose", "gau_loose"),
                    thresh=simple_cfg.get("thresh", "baker"),
                    update_interval_hessian=int(simple_cfg.get("update_interval_hessian", 500)),
                    neg_freq_thresh_cm=float(simple_cfg.get("neg_freq_thresh_cm", 5.0)),
                    flatten_amp_ang=float(simple_cfg.get("flatten_amp_ang", 0.10)),
                    flatten_max_iter=int(simple_cfg.get("flatten_max_iter", 0)),
                    mem=int(simple_cfg.get("mem", 100000)),
                    uma_kwargs=calc_kwargs_for_dimer,
                    device=str(simple_cfg.get("device", calc_cfg.get("device", "auto"))),
                    dump=bool(opt_cfg["dump"]),
                    root=int(simple_cfg.get("root", 0)),
                    dimer_kwargs=dict(simple_cfg.get("dimer", {})),
                    lbfgs_kwargs=dict(simple_cfg.get("lbfgs", {})),
                    max_total_cycles=int(opt_cfg["max_cycles"]),
                    flatten_sep_cutoff=float(simple_cfg.get("flatten_sep_cutoff", 0.0)),
                    flatten_k=int(simple_cfg.get("flatten_k", 10)),
                    flatten_loop_bofill=bool(simple_cfg.get("flatten_loop_bofill", False)),
                    # Propagate geometry settings (freeze_atoms, coord_type, ...) to the HessianDimer runner
                    geom_kwargs=dict(geom_cfg),
                    prepared_input=prepared_input,
                )

                emit("\n====== TS optimization (Hessian Guided Dimer) ======\n", narrative=True)
                runner.run()
                emit_optimizer_terminal_status(
                    "tsopt",
                    converged=getattr(runner, "is_converged", None),
                    cycles=getattr(runner, "_cycles_spent", None),
                    max_cycles=int(opt_cfg.get("max_cycles", 0)) or None,
                    stalled=getattr(runner, "is_stalled", False),
                    stop_reason=getattr(runner, "stop_reason", None) or None,
                )

                needs_pdb = source_path.suffix.lower() == ".pdb"
                needs_gjf = prepared_input.is_gjf
                ref_pdb = source_path.resolve() if needs_pdb else None

                final_xyz = out_dir_path / "final_geometry.xyz"
                if convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=out_dir_path / "final_geometry.pdb" if needs_pdb else None,
                    out_gjf_path=out_dir_path / "final_geometry.gjf" if needs_gjf else None,
                    context="final geometry",
                ):
                    click.echo("[convert] Wrote 'final_geometry' outputs.")

                if bool(opt_cfg.get("dump", False)) and needs_pdb:
                    all_trj = out_dir_path / "optimization_all_trj.xyz"
                    if all_trj.exists():
                        if convert_xyz_like_outputs(
                            all_trj,
                            prepared_input,
                            ref_pdb_path=ref_pdb,
                            out_pdb_path=out_dir_path / "optimization_all.pdb" if needs_pdb else None,
                            context="optimization trajectory",
                        ):
                            click.echo("[convert] Wrote 'optimization_all' outputs.")
                    else:
                        click.echo("[convert] WARNING: 'optimization_all_trj.xyz' not found; skipping conversion.", err=True)

                # Collect result data for --out-json (dimer mode)
                if out_json:
                    _dimer_energy = _calc_energy(runner.geom, dict(calc_cfg))
                    _dimer_imag = list(runner.imaginary_frequencies_cm)
                    _tsopt_result_data = {
                        "status": _tsopt_terminal_status(
                            runner,
                            saddle_verified=runner.saddle_order_verified,
                        ),
                        "energy_hartree": _dimer_energy,
                        "n_imaginary_modes": runner.n_imaginary_modes,
                        "imaginary_frequencies_cm": _dimer_imag,
                        "opt_mode": "dimer",
                        "n_atoms": len(runner.geom.atomic_numbers),
                        "n_opt_cycles": runner._cycles_spent,
                        "backend": calc_cfg.get("backend", backend),
                        "charge": calc_cfg["charge"],
                        "spin": calc_cfg["spin"],
                        "model": calc_cfg.get("model"),
                        "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                        "solvent": calc_cfg.get("solvent", "none"),
                        "thresh": opt_cfg.get("thresh", simple_cfg.get("thresh")),
                        "max_cycles": opt_cfg.get("max_cycles", simple_cfg.get("max_cycles")),
                        "input_file": str(input_path),
                        "files": {"final_geometry_xyz": "final_geometry.xyz"},
                        "rigid_projection": dict(runner.rigid_projection_info),
                    }
                    # Additive stop_reason, present only for a non-converged
                    # stop so a converged TS run's JSON stays byte-compatible.
                    _dimer_stop_reason = getattr(runner, "stop_reason", "") or ""
                    if _dimer_stop_reason:
                        _tsopt_result_data["stop_reason"] = _dimer_stop_reason
                    for ext in (".pdb", ".cif", ".gjf"):
                        f = out_dir_path / f"final_geometry{ext}"
                        if f.exists():
                            _tsopt_result_data["files"][f"final_geometry_{ext[1:]}"] = f.name
                    # Add trajectory files if they exist
                    for _trj_name in ("optimization_all_trj.xyz", "optimization_all.pdb", "optimization_all.cif"):
                        _tf = out_dir_path / _trj_name
                        if _tf.exists():
                            _key = _trj_name.replace(".", "_").replace("-", "_")
                            _tsopt_result_data["files"][_key] = _trj_name
                    # List imaginary mode vib files
                    _vib_dir = out_dir_path / "vib"
                    if _vib_dir.exists():
                        _tsopt_result_data["files"]["imaginary_mode_files"] = sorted([
                            f"vib/{f.name}" for f in _vib_dir.iterdir() if f.suffix in ('.xyz', '.pdb', '.cif')
                        ])

            else:
                # hess-family optimizer (RS-P-RFO default / RS-I-RFO / TRIM)
                # Build the geometry and attach the configured MLIP calculator.
                coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
                coord_kwargs = dict(geom_cfg)
                coord_kwargs.pop("coord_type", None)
                geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)
                initial_ts_cart_coords = geometry.cart_coords.copy()

                calc = create_calculator(**calc_cfg)  # includes freeze_atoms / hessian_calc_mode / partial Hessian
                geometry.set_calculator(calc)
                echo_resolved_device()

                # Construct RSIRFO optimizer
                rsirfo_kwargs = _build_rsirfo_kwargs(opt_cfg, rsirfo_cfg, out_dir_path)
                reference_mode = None
                if reference_mode_path is not None:
                    try:
                        if reference_mode_path.suffix.lower() == ".npy":
                            reference_mode = np.load(reference_mode_path)
                        else:
                            reference_mode = np.loadtxt(reference_mode_path)
                    except (OSError, ValueError) as exc:
                        raise click.ClickException(
                            f"Failed to read --ref-mode: {exc}"
                        ) from exc
                    reference_mode = np.asarray(reference_mode, dtype=float).reshape(-1)
                    expected = int(geometry.cart_coords.size)
                    if reference_mode.size != expected:
                        raise click.ClickException(
                            "--ref-mode must contain one Cartesian value per "
                            f"degree of freedom ({reference_mode.size} read, {expected} expected)."
                        )
                    rsirfo_kwargs["reference_mode"] = reference_mode

                # TRIM / RS-P-RFO are pure-numpy, single-pass TS optimizers and
                # do not support the torch-backed line_search used by RS-I-RFO.
                if kind != "rsirfo":
                    rsirfo_kwargs["min_line_search"] = False
                    rsirfo_kwargs["max_line_search"] = False

                optimizer = TSOPT_CLASS_MAP[kind](geometry, **rsirfo_kwargs)

                emit("\n====== TS optimization (" + TSOPT_DISPLAY_NAME.get(kind, kind.upper()) + ") ======\n", narrative=True)
                optimizer.run()
                emit_optimizer_terminal_status(
                    "tsopt",
                    converged=getattr(optimizer, "is_converged", None),
                    cycles=optimizer_cycle_count(optimizer),
                    max_cycles=int(opt_cfg.get("max_cycles", 0)) or None,
                    stalled=getattr(optimizer, "is_stalled", False),
                    stop_reason=getattr(optimizer, "stop_reason", None) or None,
                )
                last_optimizer = optimizer

                # --- RSIRFO: count imaginary modes and optional flatten loop ---
                geometry.set_calculator(None)
                calc_kwargs_for_rsirfo = dict(calc_cfg)
                calc_kwargs_for_rsirfo["out_hess_torch"] = True
                device = _torch_device(calc_cfg.get("device", "auto"))

                def _attach_rsirfo_calc() -> None:
                    geometry.set_calculator(calc)

                rigid_projection_info: Dict[str, Any] = {}

                def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                    H = _calc_full_hessian_torch(geometry, calc_kwargs_for_rsirfo, device)
                    from pdb2reaction.io.hessian_cache import (
                        store as _hess_store,
                        identity_from_context as _hess_identity,
                    )
                    freeze_atoms = list(geom_cfg.get("freeze_atoms", []))
                    if freeze_atoms and H.shape[0] < 3 * len(geometry.atomic_numbers):
                        all_dofs = set(range(3 * len(geometry.atomic_numbers)))
                        frozen_dofs = set()
                        for _fi in freeze_atoms:
                            frozen_dofs.update([3 * _fi, 3 * _fi + 1, 3 * _fi + 2])
                        _active_dofs = sorted(all_dofs - frozen_dofs)
                    else:
                        _active_dofs = None
                    _ts_calc_cfg = dict(calc_cfg)
                    _ts_calc_cfg.setdefault("freeze_atoms", freeze_atoms)
                    _hess_store(
                        "ts",
                        H,
                        active_dofs=_active_dofs,
                        meta={
                            "cart_coords": geometry.cart_coords,
                            "source": "tsopt_exact",
                        },
                        identity=_hess_identity(geometry, _ts_calc_cfg, role="ts"),
                    )
                    freqs_local, modes_local = _frequencies_cm_and_modes(
                        H,
                        geometry.atomic_numbers,
                        geometry.cart_coords.reshape(-1, 3),
                        device,
                        freeze_idx=list(geom_cfg.get("freeze_atoms", [])) if len(geom_cfg.get("freeze_atoms", [])) > 0 else None,
                        tr_projection=geom_cfg["tr_projection"],
                        projection_info=rigid_projection_info,
                    )
                    rigid_projection_info.update({
                        "hessian_space": (
                            "full" if H.shape[0] == 3 * len(geometry.atomic_numbers)
                            else "active"
                        ),
                        "raw_hessian_shape": list(H.shape),
                        "source": "tsopt_exact",
                    })
                    del H
                    return freqs_local, modes_local

                freqs_cm, modes = _calc_freqs_and_modes()
                click.echo(pretty_block("rigid_projection", rigid_projection_info))

                neg_freq_thresh_cm = float(simple_cfg.get("neg_freq_thresh_cm", 5.0))
                neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm if x < -abs(neg_freq_thresh_cm)]
                emit(f"[Imaginary modes] n={n_imag}  ({ims})", narrative=True)
                _warn_if_leading_imaginary_mode_is_soft(ims)

                # A local minimum plus its Hessian does not identify the
                # intended neighboring saddle.  When the path supplied that
                # missing direction, retry a failed Hessian TS run from small
                # bounded displacements on both sides of the original HEI.
                # Keep the original result if neither same-optimizer restart
                # reaches a verified first-order saddle.
                saddle_multistart_attempts = []
                target_mode_is_negative = getattr(
                    last_optimizer,
                    "_last_exact_target_mode_is_negative",
                    None,
                )
                if (
                    not last_optimizer.is_converged
                    and not getattr(last_optimizer, "is_stalled", False)
                    and reference_mode is not None
                    and (n_imag <= 1 or target_mode_is_negative is False)
                ):
                    baseline_coords = geometry.cart_coords.copy()
                    baseline_optimizer = last_optimizer
                    baseline_freqs = freqs_cm.copy()
                    baseline_modes = modes.detach().cpu().clone()
                    restart_modes = _path_restart_mode_candidates(
                        baseline_optimizer,
                        geometry,
                        reference_mode,
                    )
                    multistart_success = False
                    path_negative_candidate = None
                    for mode_source, restart_unit in restart_modes:
                        # If the explicit MEP tangent already recovered the
                        # requested negative mode, keep that chemically direct
                        # candidate for flattening rather than trying another
                        # identity. The soft-root shell is only a true fallback.
                        if (
                            mode_source == "initial-soft-root"
                            and path_negative_candidate is not None
                        ):
                            break
                        restart_kwargs = dict(rsirfo_kwargs)
                        restart_kwargs["reference_mode"] = restart_unit
                        for amplitude_ang in PATH_MODE_RESTART_AMPLITUDES_ANG:
                            geometry.cart_coords = (
                                initial_ts_cart_coords
                                + (amplitude_ang / BOHR2ANG) * restart_unit
                            )
                            _attach_rsirfo_calc()
                            optimizer = TSOPT_CLASS_MAP[kind](
                                geometry, **restart_kwargs
                            )
                            emit(
                                "\n====== TS optimization ("
                                + TSOPT_DISPLAY_NAME.get(kind, kind.upper())
                                + ", path-mode restart "
                                + f"{mode_source} {amplitude_ang:+.2f} Å) ======\n",
                                narrative=True,
                            )
                            optimizer.run()
                            emit_optimizer_terminal_status(
                                "tsopt",
                                converged=getattr(
                                    optimizer, "is_converged", None
                                ),
                                cycles=optimizer_cycle_count(optimizer),
                                max_cycles=(
                                    int(opt_cfg.get("max_cycles", 0)) or None
                                ),
                                stalled=getattr(optimizer, "is_stalled", False),
                                stop_reason=getattr(optimizer, "stop_reason", None) or None,
                            )
                            last_optimizer = optimizer
                            geometry.set_calculator(None)
                            freqs_cm, modes = _calc_freqs_and_modes()
                            neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                            n_imag = int(np.sum(neg_mask))
                            ims = [
                                float(x)
                                for x in freqs_cm
                                if x < -abs(neg_freq_thresh_cm)
                            ]
                            attempt = {
                                "mode_source": mode_source,
                                "displacement_ang": amplitude_ang,
                                "converged": bool(optimizer.is_converged),
                                "n_imaginary": n_imag,
                                "target_mode_negative": getattr(
                                    optimizer,
                                    "_last_exact_target_mode_is_negative",
                                    None,
                                ),
                            }
                            saddle_multistart_attempts.append(attempt)
                            emit(
                                "[path-mode restart] "
                                f"source={mode_source}, "
                                f"displacement={amplitude_ang:+.2f} Å, "
                                f"converged={attempt['converged']}, "
                                f"n_imag={n_imag}",
                                narrative=True,
                            )
                            if optimizer.is_converged and n_imag == 1:
                                multistart_success = True
                                break

                            if attempt["target_mode_negative"] is True:
                                force = (
                                    optimizer.forces[-1]
                                    if optimizer.forces
                                    else geometry.cart_forces
                                )
                                if isinstance(force, torch.Tensor):
                                    force = force.detach().cpu().numpy()
                                force = np.asarray(force, dtype=float)
                                score = (
                                    abs(int(n_imag) - 1),
                                    float(np.max(np.abs(force))),
                                )
                                if (
                                    path_negative_candidate is None
                                    or score
                                    < path_negative_candidate["score"]
                                ):
                                    path_negative_candidate = {
                                        "score": score,
                                        "coords": geometry.cart_coords.copy(),
                                        "optimizer": optimizer,
                                        "freqs": freqs_cm.copy(),
                                        "modes": modes.detach().cpu().clone(),
                                    }
                        if multistart_success:
                            break

                    if not multistart_success:
                        if path_negative_candidate is not None:
                            geometry.cart_coords = path_negative_candidate["coords"]
                            last_optimizer = path_negative_candidate["optimizer"]
                            freqs_cm = path_negative_candidate["freqs"]
                            modes = path_negative_candidate["modes"]
                            neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                            n_imag = int(np.sum(neg_mask))
                            ims = [
                                float(x)
                                for x in freqs_cm
                                if x < -abs(neg_freq_thresh_cm)
                            ]
                            emit(
                                "[path-mode restart] Retaining the best "
                                "path-negative higher-order candidate for "
                                "bounded extra-mode flattening.",
                                narrative=True,
                            )
                        else:
                            geometry.cart_coords = baseline_coords
                            last_optimizer = baseline_optimizer
                            freqs_cm = baseline_freqs
                            modes = baseline_modes
                            neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                            n_imag = int(np.sum(neg_mask))
                            ims = [
                                float(x)
                                for x in freqs_cm
                                if x < -abs(neg_freq_thresh_cm)
                            ]
                        geometry.set_calculator(None)
                        # Restart attempts overwrite the shared final file;
                        # restore it together with the selected coordinates.
                        (out_dir_path / "final_geometry.xyz").write_text(
                            geometry.as_xyz(), encoding="utf-8"
                        )

                flatten_max_iter = int(simple_cfg.get("flatten_max_iter", 0))
                target_mode_is_negative = getattr(
                    last_optimizer,
                    "_last_exact_target_mode_is_negative",
                    None,
                )
                flatten_max_iter, flatten_vetoed = _effective_flatten_iterations(
                    flatten_max_iter,
                    has_reference_mode=reference_mode is not None,
                    n_imag=n_imag,
                    target_mode_is_negative=target_mode_is_negative,
                )
                if flatten_vetoed:
                    # There is no path-correlated negative mode to preserve.
                    # Flattening unrelated negatives could manufacture the
                    # wrong n_imag=1 saddle, so leave this run non-converged.
                    click.echo(
                        "[flatten] Skipping extra-mode flattening because the "
                        "path-correlated mode is not negative.",
                        err=True,
                    )
                if (
                    flatten_max_iter > 0
                    and n_imag > 1
                    and not getattr(last_optimizer, "is_stalled", False)
                ):
                    click.echo("[flatten] Extra imaginary modes detected; starting RSIRFO flatten loop.")
                    masses_amu = _safe_masses_amu(geometry.atomic_numbers)
                    roots = rsirfo_kwargs.get("roots", [0])
                    main_root = int(roots[0]) if roots else 0

                    def _run_flatten_branch(
                        start_coords: np.ndarray,
                        branch_reference_mode: Optional[np.ndarray],
                        branch_label: str,
                    ) -> Dict[str, Any]:
                        """Optimize and exactly characterize one signed branch."""
                        geometry.cart_coords = np.asarray(
                            start_coords, dtype=float
                        ).copy()
                        _attach_rsirfo_calc()
                        retry_kwargs = dict(rsirfo_kwargs)
                        retry_kwargs["max_higher_order_checks"] = (
                            FLATTEN_RETRY_HIGHER_ORDER_CHECKS
                        )
                        if branch_reference_mode is not None:
                            retry_kwargs["reference_mode"] = (
                                branch_reference_mode
                            )
                        branch_optimizer = TSOPT_CLASS_MAP[kind](
                            geometry, **retry_kwargs
                        )
                        emit(
                            "\n====== TS optimization ("
                            + TSOPT_DISPLAY_NAME.get(kind, kind.upper())
                            + f", flatten retry {branch_label}) ======\n",
                            narrative=True,
                        )
                        branch_optimizer.run()
                        emit_optimizer_terminal_status(
                            "tsopt",
                            converged=getattr(
                                branch_optimizer, "is_converged", None
                            ),
                            cycles=optimizer_cycle_count(branch_optimizer),
                            max_cycles=(
                                int(opt_cfg.get("max_cycles", 0)) or None
                            ),
                            stalled=getattr(branch_optimizer, "is_stalled", False),
                            stop_reason=getattr(branch_optimizer, "stop_reason", None) or None,
                        )
                        geometry.set_calculator(None)
                        branch_freqs, branch_modes = _calc_freqs_and_modes()
                        branch_neg = branch_freqs < -abs(
                            neg_freq_thresh_cm
                        )
                        branch_n_imag = int(np.sum(branch_neg))
                        branch_ims = [
                            float(value)
                            for value in branch_freqs
                            if value < -abs(neg_freq_thresh_cm)
                        ]
                        emit(
                            f"[Imaginary modes:{branch_label}] "
                            f"n={branch_n_imag}  ({branch_ims})",
                            narrative=True,
                        )
                        return {
                            "label": branch_label,
                            "optimizer": branch_optimizer,
                            "coords": geometry.cart_coords.copy(),
                            "freqs": branch_freqs.copy(),
                            "modes": branch_modes.detach().cpu().clone(),
                            "n_imag": branch_n_imag,
                            "ims": branch_ims,
                        }

                    def _flatten_branch_score(
                        result: Dict[str, Any],
                    ) -> Tuple[int, int, float, float]:
                        """Prefer the requested n=1 saddle, then weak surplus."""
                        branch_optimizer = result["optimizer"]
                        target_negative = (
                            getattr(
                                branch_optimizer,
                                "_last_exact_target_mode_is_negative",
                                None,
                            )
                            is True
                        )
                        branch_freqs = np.asarray(result["freqs"], dtype=float)
                        negative_indices = np.flatnonzero(
                            branch_freqs < -abs(neg_freq_thresh_cm)
                        )
                        target_index = getattr(
                            branch_optimizer,
                            "_last_exact_target_mode_index",
                            None,
                        )
                        surplus_strength = 0.0
                        for mode_index in negative_indices:
                            if (
                                target_negative
                                and target_index is not None
                                and int(mode_index) == int(target_index)
                            ):
                                continue
                            surplus_strength += abs(
                                float(branch_freqs[int(mode_index)])
                            )
                        force = (
                            branch_optimizer.forces[-1]
                            if branch_optimizer.forces
                            else np.asarray([], dtype=float)
                        )
                        if isinstance(force, torch.Tensor):
                            force = force.detach().cpu().numpy()
                        force = np.asarray(force, dtype=float).reshape(-1)
                        max_force = (
                            float(np.max(np.abs(force)))
                            if force.size
                            else float("inf")
                        )
                        return (
                            0 if target_negative else 1,
                            abs(int(result["n_imag"]) - 1),
                            surplus_strength,
                            max_force,
                        )

                    for it in range(flatten_max_iter):
                        click.echo(f"[flatten] RSIRFO iteration {it + 1}/{flatten_max_iter}")
                        flatten_reference_mode = _transported_path_mode_full(
                            last_optimizer,
                            geometry,
                            reference_mode,
                        )
                        pre_flatten_coords = geometry.cart_coords.copy()
                        did_flatten = _flatten_once_with_modes_for_geom(
                            geometry,
                            masses_amu,
                            calc_kwargs_for_rsirfo,
                            freqs_cm,
                            modes,
                            neg_freq_thresh_cm,
                            float(simple_cfg.get("flatten_amp_ang", 0.10)),
                            float(simple_cfg.get("flatten_sep_cutoff", 0.0)),
                            int(simple_cfg.get("flatten_k", 10)),
                            main_root,
                            reference_mode=flatten_reference_mode,
                        )
                        if not did_flatten:
                            click.echo("[flatten] No eligible modes to flatten; stopping.")
                            break
                        # ``_flatten_once`` has already chosen its lower-energy
                        # signed displacement.  Mirror that displacement about
                        # the higher-order saddle to obtain the other physical
                        # descent branch.  Energy alone cannot decide which
                        # branch retains the requested reaction mode.
                        primary_start = geometry.cart_coords.copy()
                        primary_result = _run_flatten_branch(
                            primary_start,
                            flatten_reference_mode,
                            "primary",
                        )
                        selected_result = primary_result
                        if (
                            reference_mode is not None
                            and _flatten_branch_needs_alternate(primary_result)
                        ):
                            alternate_start = _mirrored_flatten_start(
                                pre_flatten_coords,
                                primary_start,
                            )
                            alternate_result = _run_flatten_branch(
                                alternate_start,
                                flatten_reference_mode,
                                "alternate",
                            )
                            primary_score = _flatten_branch_score(
                                primary_result
                            )
                            alternate_score = _flatten_branch_score(
                                alternate_result
                            )
                            if alternate_score < primary_score:
                                selected_result = alternate_result
                            emit(
                                "[flatten] Signed-branch probe selected "
                                f"{selected_result['label']} "
                                f"(primary={primary_score}, "
                                f"alternate={alternate_score}).",
                                narrative=True,
                            )

                        optimizer = selected_result["optimizer"]
                        geometry.cart_coords = selected_result["coords"].copy()
                        last_optimizer = optimizer
                        freqs_cm = selected_result["freqs"]
                        modes = selected_result["modes"]
                        n_imag = int(selected_result["n_imag"])
                        ims = list(selected_result["ims"])
                        # Both signed trials share the same final filename.
                        # Restore it to the selected coordinates as well as the
                        # in-memory geometry.
                        selected_final = getattr(
                            optimizer, "final_fn", out_dir_path / "final_geometry.xyz"
                        )
                        Path(selected_final).write_text(
                            geometry.as_xyz(), encoding="utf-8"
                        )
                        if (
                            reference_mode is not None
                            and getattr(
                                optimizer,
                                "_last_exact_target_mode_is_negative",
                                None,
                            )
                            is not True
                        ):
                            click.echo(
                                "[flatten] Path-correlated mode was lost; "
                                "stopping flatten retries without preserving "
                                "an unrelated negative mode.",
                                err=True,
                            )
                            break
                        if n_imag <= 1:
                            break

                needs_pdb = source_path.suffix.lower() == ".pdb"
                needs_gjf = prepared_input.is_gjf
                ref_pdb = source_path.resolve() if needs_pdb else None
                final_xyz_path = last_optimizer.final_fn if isinstance(last_optimizer.final_fn, Path) else Path(last_optimizer.final_fn)
                if convert_xyz_like_outputs(
                    final_xyz_path,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=out_dir_path / "final_geometry.pdb" if needs_pdb else None,
                    out_gjf_path=out_dir_path / "final_geometry.gjf" if needs_gjf else None,
                    context="final geometry",
                ):
                    click.echo("[convert] Wrote 'final_geometry' outputs.")

                if bool(opt_cfg.get("dump", False)) and needs_pdb:
                    trj_path = out_dir_path / "optimization_trj.xyz"
                    if trj_path.exists():
                        if convert_xyz_like_outputs(
                            trj_path,
                            prepared_input,
                            ref_pdb_path=ref_pdb,
                            out_pdb_path=out_dir_path / "optimization.pdb" if needs_pdb else None,
                            context="optimization trajectory",
                        ):
                            click.echo("[convert] Wrote 'optimization' outputs.")
                    else:
                        click.echo("[convert] WARNING: 'optimization_trj.xyz' not found; skipping conversion.", err=True)

                # --- RSIRFO: write all final imaginary modes ---
                neg_idx = _echo_imaginary_modes(freqs_cm, neg_freq_thresh_cm)
                _projection_certifiable = allows_saddle_certification(
                    geom_cfg.get("tr_projection", "constrained"),
                    geom_cfg.get("freeze_atoms", ()),
                )
                _hessian_saddle_verified = bool(
                    _projection_certifiable and len(neg_idx) == 1
                )
                if not _hessian_saddle_verified:
                    last_optimizer.is_converged = False
                    if len(neg_idx) == 0:
                        warning = "No imaginary mode found"
                    elif len(neg_idx) > 1:
                        warning = f"Found {len(neg_idx)} imaginary modes"
                    else:
                        warning = (
                            "legacy-active projection is comparison-only for "
                            "frozen systems"
                        )
                    click.echo(
                        f"[WARNING] {warning} in the final exact Hessian; "
                        "this geometry is not a first-order saddle and is "
                        "marked not_converged.",
                        err=True,
                    )
                if len(neg_idx) > 0:
                    vib_dir = out_dir_path / "vib"
                    ref_pdb_mode = source_path if source_path.suffix.lower() == ".pdb" else None
                    _write_all_imaginary_modes(
                        geometry,
                        freqs_cm,
                        modes,
                        neg_freq_thresh_cm,
                        _safe_masses_amu(geometry.atomic_numbers).tolist(),
                        vib_dir,
                        prepared_input=prepared_input,
                        ref_pdb=ref_pdb_mode,
                        amplitude_ang=0.8,
                        n_frames=20,
                    )

                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                # Collect result data for --out-json (hess-family: rsprfo/rsirfo/trim)
                if out_json:
                    _rsirfo_imag = [float(f) for f in freqs_cm if f < -abs(neg_freq_thresh_cm)]
                    # Frequency analysis and failed multistart restoration leave
                    # the shared Geometry intentionally detached.  Evaluate the
                    # final energy through the already-loaded calculator instead
                    # of dereferencing ``geometry.energy`` with calculator=None.
                    _rsirfo_energy = _calc_energy(
                        geometry, calc_cfg, calc=calc
                    )
                    _tsopt_result_data = {
                        "status": _tsopt_terminal_status(
                            last_optimizer,
                            saddle_verified=_hessian_saddle_verified,
                            projection_certifiable=_projection_certifiable,
                        ),
                        "energy_hartree": _rsirfo_energy,
                        "n_imaginary_modes": len(_rsirfo_imag),
                        "imaginary_frequencies_cm": _rsirfo_imag,
                        "opt_mode": kind,
                        "n_atoms": len(geometry.atoms),
                        "n_opt_cycles": last_optimizer.cur_cycle if hasattr(last_optimizer, "cur_cycle") else None,
                        "backend": calc_cfg.get("backend", backend),
                        "charge": calc_cfg["charge"],
                        "spin": calc_cfg["spin"],
                        "model": calc_cfg.get("model"),
                        "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                        "solvent": calc_cfg.get("solvent", "none"),
                        "thresh": opt_cfg.get("thresh", simple_cfg.get("thresh")),
                        "max_cycles": opt_cfg.get("max_cycles", simple_cfg.get("max_cycles")),
                        "input_file": str(input_path),
                        "reference_mode_file": (
                            None
                            if reference_mode_path is None
                            else str(reference_mode_path)
                        ),
                        "files": {"final_geometry_xyz": str(final_xyz_path.name)},
                        "rigid_projection": dict(rigid_projection_info),
                        "safeguards": {
                            "rejected_mode_loss_trials": int(
                                getattr(last_optimizer, "rejected_mode_loss_steps", 0)
                            ),
                            "exact_saddle_checks": int(
                                getattr(last_optimizer, "exact_saddle_checks", 0)
                            ),
                            "saddle_recovery_steps": int(
                                getattr(last_optimizer, "saddle_recovery_steps", 0)
                            ),
                            "saddle_recovery_check_interval": int(
                                getattr(
                                    last_optimizer,
                                    "saddle_recovery_check_interval",
                                    0,
                                )
                            ),
                            "saddle_recovery_max_cycles": int(
                                getattr(
                                    last_optimizer,
                                    "saddle_recovery_max_cycles",
                                    0,
                                )
                            ),
                            "last_exact_n_imaginary": getattr(
                                last_optimizer, "_last_exact_n_imaginary", None
                            ),
                            "last_exact_target_mode_index": getattr(
                                last_optimizer,
                                "_last_exact_target_mode_index",
                                None,
                            ),
                            "last_exact_target_mode_overlap": getattr(
                                last_optimizer,
                                "_last_exact_target_mode_overlap",
                                None,
                            ),
                            "last_exact_target_mode_is_negative": getattr(
                                last_optimizer,
                                "_last_exact_target_mode_is_negative",
                                None,
                            ),
                            "last_exact_target_mode_reanchored": bool(
                                getattr(
                                    last_optimizer,
                                    "_last_exact_target_mode_reanchored",
                                    False,
                                )
                            ),
                            "initial_reference_root_index": getattr(
                                last_optimizer,
                                "_initial_reference_root_index",
                                None,
                            ),
                            "initial_reference_root_overlap": getattr(
                                last_optimizer,
                                "_initial_reference_root_overlap",
                                None,
                            ),
                            "initial_reference_root_eigenvalue": getattr(
                                last_optimizer,
                                "_initial_reference_root_eigenvalue",
                                None,
                            ),
                            "last_recovery_mode_curvature": getattr(
                                last_optimizer,
                                "_last_recovery_mode_curvature",
                                None,
                            ),
                            "last_recovery_mode_frequency_cm": getattr(
                                last_optimizer,
                                "_last_recovery_mode_frequency_cm",
                                None,
                            ),
                            "stop_reason": str(
                                getattr(last_optimizer, "stop_reason", "")
                            ),
                            "path_mode_multistart_attempts": saddle_multistart_attempts,
                        },
                    }
                    # Convergence details from optimizer
                    if hasattr(last_optimizer, 'max_forces') and last_optimizer.max_forces:
                        _tsopt_result_data["final_max_force"] = float(last_optimizer.max_forces[-1])
                        _tsopt_result_data["final_rms_force"] = float(last_optimizer.rms_forces[-1])
                    if hasattr(last_optimizer, 'convergence') and last_optimizer.convergence:
                        _tsopt_result_data["convergence_thresholds"] = {k: float(v) for k, v in last_optimizer.convergence.items()}
                    if hasattr(last_optimizer, 'max_steps') and last_optimizer.max_steps:
                        _tsopt_result_data["final_max_step"] = float(last_optimizer.max_steps[-1])
                        _tsopt_result_data["final_rms_step"] = float(last_optimizer.rms_steps[-1])
                    for ext in (".pdb", ".cif", ".gjf"):
                        f = out_dir_path / f"final_geometry{ext}"
                        if f.exists():
                            _tsopt_result_data["files"][f"final_geometry_{ext[1:]}"] = f.name
                    # Add trajectory files if they exist
                    for _trj_name in ("optimization_trj.xyz", "optimization.pdb", "optimization.cif"):
                        _tf = out_dir_path / _trj_name
                        if _tf.exists():
                            _key = _trj_name.replace(".", "_").replace("-", "_")
                            _tsopt_result_data["files"][_key] = _trj_name
                    # List imaginary mode vib files
                    _vib_dir = out_dir_path / "vib"
                    if _vib_dir.exists():
                        _tsopt_result_data["files"]["imaginary_mode_files"] = sorted([
                            f"vib/{f.name}" for f in _vib_dir.iterdir() if f.suffix in ('.xyz', '.pdb', '.cif')
                        ])

            if _tsopt_result_data is not None:
                _tsopt_failed = _tsopt_result_data.get("status") != "converged"
            elif kind == "dimer":
                _tsopt_failed = not bool(getattr(runner, "is_converged", False))
            else:
                _tsopt_failed = not bool(
                    getattr(last_optimizer, "is_converged", False)
                )
            if _tsopt_failed:
                click.echo(
                    "[tsopt] TS optimization did not reach a validated "
                    "first-order saddle. Retry with --flatten when surplus "
                    "imaginary modes remain. In the all workflow, "
                    "--refine-path can improve a poor MEP before tsopt, "
                    "but recursive refinement may split the path into multiple "
                    "segments and substantially increase cost; it is off by "
                    "default.",
                    err=True,
                )

            # summary.md and key_* outputs are disabled.
            emit(format_elapsed("[time] Elapsed Time for TS Opt", time_start), narrative=True)

            # result.json (if --out-json)
            if out_json:
                from pdb2reaction.core.utils import calculator_provenance, write_result_json
                _tsopt_result_data.update(calculator_provenance(calc_cfg))
                write_result_json(
                    out_dir_path, _tsopt_result_data,
                    command="tsopt",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

        except ZeroStepLength as e:
            _write_error_json(out_dir_path, "tsopt", e, "ZeroStepLength", time_start)
            click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
            sys.exit(2)
        except OptimizationError as e:
            _write_error_json(out_dir_path, "tsopt", e, "OptimizationError", time_start)
            click.echo(f"ERROR: Optimization failed — {e}", err=True)
            sys.exit(3)
        except KeyboardInterrupt:
            click.echo("Interrupted by user.", err=True)
            sys.exit(130)
        except Exception as e:
            render_cli_exception(e, label="optimization", out_dir=out_dir_path, command="tsopt", time_start=time_start)
        finally:
            # Release GPU memory (model + Hessian) so subsequent stages don't OOM
            calc = geometry = optimizer = last_optimizer = None
            freqs_cm = modes = None
            gc.collect()  # break cyclic refs inside torch.nn.Module
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
