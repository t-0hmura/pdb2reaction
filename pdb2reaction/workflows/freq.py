# pdb2reaction/workflows/freq.py

"""
Vibrational frequency analysis with PHVA support and thermochemistry output.

Example:
    pdb2reaction freq -i a.pdb -q 0 -m 1

For detailed documentation, see: docs/freq.md
"""

from __future__ import annotations

import gc
import logging
import sys
import textwrap
from pathlib import Path
from typing import Any, Optional, Tuple, List, Sequence

import click
import numpy as np
import torch
from ase.data import atomic_masses
import ase.units as units
import yaml
import time

# ---------------- pysisyphus / pdb2reaction imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AMU2AU, AU2EV

# local helpers from pdb2reaction
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import GEOM_KW_DEFAULT, FREQ_CALC_KW, FREQ_KW, THERMO_KW, apply_backend_defaults
from pdb2reaction.core.utils import (
    apply_yaml_overrides,
    convert_xyz_like_outputs,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    resolve_freeze_atoms,
    cli_param_overridden,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
)
from pdb2reaction.cli.common_options import add_ml_charge_spin_options, add_precision_option, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option, add_coord_type_option
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, render_cli_exception

logger = logging.getLogger(__name__)


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


def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)



def _build_tr_basis(coords_bohr_t: torch.Tensor,
                    masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Mass-weighted translation/rotation basis (Tx, Ty, Tz, Rx, Ry, Rz), shape (3N, r<=6).
    """
    device, dtype = coords_bohr_t.device, coords_bohr_t.dtype
    N = coords_bohr_t.shape[0]
    m_au = masses_au_t.to(dtype=dtype, device=device)
    m_sqrt = torch.sqrt(m_au).reshape(-1, 1)

    com = (m_au.reshape(-1, 1) * coords_bohr_t).sum(0) / m_au.sum()
    x = coords_bohr_t - com

    eye3 = torch.eye(3, dtype=dtype, device=device)
    cols = []
    for i in range(3):
        cols.append((eye3[i].repeat(N, 1) * m_sqrt).reshape(-1, 1))
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


def _mw_projected_hessian(H: torch.Tensor,
                          coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor) -> torch.Tensor:
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

        Q, _ = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)
        Qt = Q.T

        QtH = Qt @ H
        H.addmm_(Q, QtH, beta=1.0, alpha=-1.0)

        HQ = QtH.T
        H.addmm_(HQ, Qt, beta=1.0, alpha=-1.0)

        QtHQ = QtH @ Q
        tmp = Q @ QtHQ
        H.addmm_(tmp, Qt, beta=1.0, alpha=1.0)

        # Bounded-peak symmetrization (helper writes both triangles; peak temp
        # <= chunk^2 instead of full N×N clone).
        from pdb2reaction.core.utils import symmetrize_inplace
        symmetrize_inplace(H)

        del masses_amu_t, m3, inv_sqrt_m, inv_sqrt_m_col, inv_sqrt_m_row
        del Q, Qt, QtH, HQ, QtHQ, tmp

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
                              freeze_idx: Optional[List[int]] = None) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize a (possibly PHVA/active-subspace) TR-projected mass-weighted Hessian
    to obtain frequencies (cm^-1) and mass-weighted eigenvectors (modes).

    If `freeze_idx` is provided (list of 0-based atom indices), perform
    Partial Hessian Vibrational Analysis (PHVA). Supports two cases:

      A) Full Hessian given (3N×3N):
         1) mass-weight the full Hessian
         2) take the active subspace by removing DOF of frozen atoms
         3) perform TR projection **only in the active subspace**
         4) diagonalize and embed eigenvectors back to 3N by zero-filling frozen DOF

      B) Already-reduced (active-block) Hessian given (3N_act×3N_act):
         1) mass-weight using only active masses
         2) TR projection in the active space
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
                freqs_cm = np.zeros((0,), dtype=float)
                modes = torch.zeros((0, 3 * N), dtype=H.dtype, device=H.device)
                return freqs_cm, modes

            expected_act_dim = 3 * n_active
            is_partial = (H.shape[0] == expected_act_dim and H.shape[1] == expected_act_dim)

            if is_partial:
                masses_act = masses_au_t[active_idx]
                coords_act = coords_bohr_t[active_idx, :]

                Hmw_act = _mass_weighted_hessian(H, masses_act)

                Q, _ = _tr_orthonormal_basis(coords_act, masses_act)
                Qt = Q.T
                QtH = Qt @ Hmw_act
                Hmw_act.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
                Hmw_act.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)
                QtHQ = QtH @ Q
                Hmw_act.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
                # Bounded-peak symmetrization (helper writes both triangles).
                from pdb2reaction.core.utils import symmetrize_inplace
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
                del Q, Qt, QtH, QtHQ, mask_dof

            else:
                H = _mass_weighted_hessian(H, masses_au_t)

                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=H.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False

                H = H[mask_dof][:, mask_dof]
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                coords_act = coords_bohr_t[active_idx, :]
                masses_act = masses_au_t[active_idx]
                Q, _ = _tr_orthonormal_basis(coords_act, masses_act)
                Qt = Q.T

                QtH = Qt @ H
                H.addmm_(Q, QtH, beta=1.0, alpha=-1.0)

                H.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)

                QtH = QtH @ Q
                H.addmm_(Q @ QtH, Qt, beta=1.0, alpha=1.0)
                # Bounded-peak symmetrization (helper writes both triangles).
                from pdb2reaction.core.utils import symmetrize_inplace
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
                del Vsub, mask_dof, Q, Qt, QtH

        else:
            H = _mw_projected_hessian(H, coords_bohr_t, masses_au_t)
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


def _calc_full_hessian_torch(
    geom,
    uma_kwargs: dict,
    device: torch.device,
    *,
    refresh_geom_meta: bool = False,
) -> torch.Tensor:
    """
    Return Hessian as torch.Tensor in Hartree/Bohr^2 (3N or 3N_act square).

    Reuses Geometry-level Hessian cache when available so repeated requests on the
    same coordinates do not trigger extra calculator evaluations.

    Parameters
    ----------
    refresh_geom_meta : bool, default False
        Accepted for signature parity with sibling toolchains; currently a no-op
        here. ``geom.set_results(results)`` is invoked unconditionally below and
        that call already refreshes ``geom.within_partial_hessian`` via
        ``Geometry.set_results`` (pysisyphus/Geometry.py:1422-1425), which is
        the only partial-Hessian metadata channel this package consumes.
    """
    def _to_torch(H_obj: Any, clone: bool = True) -> torch.Tensor:
        if isinstance(H_obj, torch.Tensor):
            t = H_obj.to(device=device)
            return t.clone() if clone else t
        return torch.as_tensor(H_obj, device=device)

    cached = getattr(geom, "_hessian", None)
    if cached is not None:
        return _to_torch(cached, clone=True)

    kw = dict(uma_kwargs or {})
    kw["out_hess_torch"] = True
    calc = create_calculator(**kw)
    echo_resolved_device()
    results = calc.get_hessian(geom.atoms, geom.cart_coords)

    # Keep Geometry cache in sync so optimizers/freq analysis can share one Hessian
    # evaluation on unchanged coordinates.
    try:
        geom.set_results(results)
    except Exception:
        import logging
        logging.getLogger(__name__).warning(
            "Failed to set Hessian results on Geometry cache", exc_info=True
        )

    # Clone so downstream in-place mass-weighting / TR projection does not
    # poison the Hessian cached on the Geometry object.
    return _to_torch(results["hessian"], clone=True)


def _calc_energy(geom, uma_kwargs: dict, calc=None) -> float:
    """
    Compute electronic energy (Hartree) from UMA calculator.
    """
    if calc is None:
        calc = create_calculator(**uma_kwargs)
    geom.set_calculator(calc)
    E = float(geom.energy)
    geom.set_calculator(None)
    return E


def _fmt_ha(x: float) -> str:
    return f"{float(x): .6f} Ha"


def _fmt_cal(x: float) -> str:
    return f"{float(x): .2f} cal/mol"


def _fmt_calK(x: float) -> str:
    return f"{float(x): .2f} cal/(mol*K)"


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            amplitude_ang: float = 0.8,
                            n_frames: int = 20,
                            comment: str = "mode",
                            ref_pdb: Optional[Path] = None,
                            write_pdb: bool = True,
                            prepared_input: Optional["PreparedInputStructure"] = None,
                            out_pdb: Optional[Path] = None) -> None:
    """
    Write a single mode animation as _trj.xyz (XYZ-like) and optionally .pdb.

    If `ref_pdb` is provided and is a .pdb file, the .pdb is generated by
    converting the _trj.xyz using the input PDB as the template.
    Set `write_pdb=False` to skip PDB generation.
    """
    ref_ang = geom.cart_coords.reshape(-1, 3) * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # _trj.xyz (concatenated XYZ-like trajectory)
    try:
        from pysisyphus.xyzloader import make_trj_str  # type: ignore
        amp_ang = amplitude_ang
        steps = np.sin(2.0 * np.pi * np.arange(n_frames) / n_frames)[:, None, None] * (amp_ang * mode[None, :, :])
        traj_ang = ref_ang[None, :, :] + steps  # (T,N,3) in Å
        comments = [f"{comment}  frame={i+1}/{n_frames}" for i in range(n_frames)]
        trj_str = make_trj_str(geom.atoms, traj_ang, comments=comments)
        out_trj.write_text(trj_str, encoding="utf-8")
    except Exception:
        with out_trj.open("w", encoding="utf-8") as f:
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                coords = ref_ang + phase * amplitude_ang * mode
                f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
                for sym, (x, y, z) in zip(geom.atoms, coords):
                    f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")

    needs_pdb = write_pdb and out_pdb is not None

    if not needs_pdb:
        return

    ref_for_conv = ref_pdb if (ref_pdb and ref_pdb.suffix.lower() == ".pdb") else None
    try:
        convert_xyz_like_outputs(
            out_trj,
            prepared_input,  # type: ignore[arg-type]
            ref_pdb_path=ref_for_conv,
            out_pdb_path=out_pdb if needs_pdb else None,
        )
    except Exception as e:
        click.echo(
            f"[convert] WARNING: Failed to convert mode trajectory '{out_trj.name}' to PDB: {e}",
            err=True,
        )


# Geometry defaults (local copy for CLI)
GEOM_KW = dict(GEOM_KW_DEFAULT)

# Calc defaults for freq (from defaults.py)
CALC_KW = FREQ_CALC_KW



@click.command(
    help="Vibrational frequency analysis and mode writer (+ default thermochemistry summary).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
# NOTE: the workers/workers_per_node help text below is intentionally explicit:
# freq DOES support workers>1 (FD Hessian path). Do not shorten the parenthetical
# when trimming docs — the verbose form is what disambiguates the supported case.
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, _trj.xyz, ...).",
)
@click.option(
    "--workers",
    type=int,
    default=CALC_KW["workers"],
    show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: analytical Hessian raises a RuntimeError when workers>1; pass --hessian-calc-mode FiniteDifference explicitly.",
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
    help="Freeze parent atoms of cap hydrogens (PDB input or XYZ/GJF with --ref-pdb).",
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
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB companions when a PDB template is available.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option("--max-write", type=int, default=FREQ_KW["max_write"], show_default=True,
              help="How many modes to export (after sorting per --sort).")
@click.option("--amplitude-ang", type=float, default=FREQ_KW["amplitude_ang"], show_default=True,
              help="Animation amplitude (Å) used for both _trj.xyz and .pdb.")
@click.option("--n-frames", type=int, default=FREQ_KW["n_frames"], show_default=True,
              help="Number of frames per mode animation.")
@click.option("--sort", type=click.Choice(["value", "abs"]), default="value", show_default=True,
              help="Sort modes by 'value' (cm^-1) or by absolute value.")
@click.option("-o", "--out-dir", type=str, default=FREQ_KW["out_dir"], show_default=True, help="Output directory.")
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Base YAML configuration file applied before explicit CLI options.",
)
# Thermochemistry options
@click.option("--temperature", type=float, default=THERMO_KW["temperature"], show_default=True,
              help="Temperature (K) for thermochemistry summary.")
@click.option("--pressure", "pressure_atm",
              type=float, default=THERMO_KW["pressure_atm"], show_default=True,
              help="Pressure (atm) for thermochemistry summary.")
@click.option(
    "--dump/--no-dump",
    "dump",
    default=THERMO_KW["dump"],
    show_default=True,
    help="When True, write 'thermoanalysis.yaml' under out-dir.",
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
    help="Validate options and print the execution plan without running frequency analysis.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
# Hessian calculation mode
@click.option("--hessian-calc-mode",
              type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
              default=None,
              help="How the ML backend computes Hessian. Defaults to 'FiniteDifference' (can also be set via YAML).")
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
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links: bool,
    freeze_atoms_text: Optional[str],
    convert_files: bool,
    ref_pdb: Optional[Path],
    max_write: int,
    amplitude_ang: float,
    n_frames: int,
    sort: str,
    out_dir: str,
    config_yaml: Optional[Path],
    # thermo
    temperature: float,
    pressure_atm: float,
    dump: bool,
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    # hessian
    hessian_calc_mode: Optional[str],
    # backend
    backend: str,
    solvent: str,
    solvent_model: str,
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: str,
    cli_coord_type: Optional[str],
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

    time_start = time.perf_counter()
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    try:
        apply_ref_pdb_override(prepared_input, ref_pdb)
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path
        charge, spin = resolve_charge_spin(
            prepared_input,
            charge,
            spin,
            ligand_charge=ligand_charge,
            prefix="[freq]",
        )
        validate_charge_spin_for_prepared(prepared_input, charge, spin)
    except BaseException:
        prepared_input.cleanup()
        raise

    geom_cfg = dict(GEOM_KW)
    calc_cfg = dict(CALC_KW)
    freq_cfg = dict(FREQ_KW)
    thermo_cfg = dict(THERMO_KW)

    apply_yaml_overrides(
        config_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",),)),
        ],
    )

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
    if precision is not None:
        from pdb2reaction.backends import apply_precision_to_calc_cfg
        apply_precision_to_calc_cfg(calc_cfg, precision)
    if backend_model is not None:
        from pdb2reaction.backends import apply_backend_model_to_calc_cfg
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
    from pdb2reaction.backends import apply_calc_file_to_calc_cfg
    apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
    apply_backend_defaults(calc_cfg)
    if cli_param_overridden(ctx, "hessian_calc_mode") and hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
    if cli_param_overridden(ctx, "max_write"):
        freq_cfg["max_write"] = int(max_write)
    if cli_param_overridden(ctx, "amplitude_ang"):
        freq_cfg["amplitude_ang"] = float(amplitude_ang)
    if cli_param_overridden(ctx, "n_frames"):
        freq_cfg["n_frames"] = int(n_frames)
    if cli_param_overridden(ctx, "sort"):
        freq_cfg["sort"] = str(sort)
    if cli_param_overridden(ctx, "out_dir"):
        freq_cfg["out_dir"] = str(out_dir)
    if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
        geom_cfg["coord_type"] = str(cli_coord_type).lower()
    if cli_param_overridden(ctx, "temperature"):
        thermo_cfg["temperature"] = float(temperature)
    if cli_param_overridden(ctx, "pressure_atm"):
        thermo_cfg["pressure_atm"] = float(pressure_atm)
    if cli_param_overridden(ctx, "dump"):
        thermo_cfg["dump"] = bool(dump)

    # `charge` / `spin` were already resolved (CLI -q/-m, gjf metadata,
    # --ligand-charge derivation) at the top of cli via resolve_charge_spin
    # and validated via validate_charge_spin_for_prepared. Assign directly;
    # an earlier .get("charge", charge) idiom silently returned the
    # CALC_KW default 0 when -q was not passed.
    calc_cfg["charge"] = int(charge)
    calc_cfg["spin"] = int(spin)

    apply_yaml_overrides(
        override_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",),)),
        ],
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
    resolve_freeze_atoms(geom_cfg, source_path, freeze_links, on_error="warn")

    # Ensure calc config reflects the geometry freeze list used in the run.
    calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
    calc_cfg["return_partial_hessian"] = True

    out_dir_path = Path(freq_cfg.get("out_dir", out_dir)).resolve()

    if show_config:
        click.echo(
            pretty_block(
                "yaml_layers",
                {
                    "config": None if config_yaml is None else str(config_yaml),
                    "override_yaml": None if override_yaml is None else str(override_yaml),
                    "merged_keys": sorted(merged_yaml_cfg.keys()),
                },
            )
        )

    if dry_run:
        click.echo(
            pretty_block(
                "dry_run_plan",
                {
                    "input_geometry": str(geom_input_path),
                    "output_dir": str(out_dir_path),
                    "freeze_links": bool(freeze_links),
                    "convert_files": bool(convert_files),
                    "will_run_frequency_analysis": True,
                    "will_write_modes": True,
                    "will_dump_thermo_yaml": bool(thermo_cfg.get("dump", False)),
                },
            )
        )
        click.echo("[dry-run] Validation complete. Frequency execution was skipped.")
        return

    out_dir_path.mkdir(parents=True, exist_ok=True)

    # Default-verbosity entry summary (skipped in child mode).
    from pdb2reaction.core.utils import echo_run_summary
    _model = calc_cfg.get("model")
    echo_run_summary({
        "input": str(input_path),
        "backend": (
            f"{calc_cfg.get('backend', '?')} ({_model}, {calc_cfg.get('precision', 'fp32')})"
            if _model else calc_cfg.get("backend", "?")
        ),
        "out": str(out_dir_path),
    })

    # Pretty-print config summary
    click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
    click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
    click.echo(pretty_block("freq", {**freq_cfg, "out_dir": str(out_dir_path)}))
    thermo_block = {
        "temperature": thermo_cfg["temperature"],
        "pressure_atm": thermo_cfg["pressure_atm"],
        "dump": thermo_cfg["dump"],
    }
    click.echo(pretty_block("thermo", thermo_block))

    coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
    coord_kwargs = dict(geom_cfg)
    coord_kwargs.pop("coord_type", None)
    geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

    # Masses (AU tensor for TR projection & MW->Cart conversion)
    masses_amu = _safe_masses_amu(geometry.atomic_numbers)
    masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=torch.float32)
    device = _torch_device(calc_cfg.get("device", "auto"))
    masses_au_t = masses_au_t.to(device=device)

    try:
        freeze_list = list(calc_cfg.get("freeze_atoms", []))
        from pdb2reaction.io.hessian_cache import load as _hess_load
        _cached_ts = _hess_load("ts")
        if _cached_ts is not None:
            click.echo("[freq] Reusing cached TS Hessian.", narrative=True)
            H = _cached_ts["hessian"]
            if isinstance(H, torch.Tensor):
                H = H.to(device=device)
            else:
                H = torch.as_tensor(H, device=device)
        else:
            H = _calc_full_hessian_torch(geometry, calc_cfg, device)
        coords_bohr = geometry.cart_coords.reshape(-1, 3)

        # PHVA: use the freeze list to carve out the active subspace and apply TR projection there.
        _n_atoms = len(geometry.atomic_numbers)
        _n_frozen = len(set(int(i) for i in freeze_list))
        _n_active = max(_n_atoms - _n_frozen, 0)
        click.echo(
            f"[freq] Hessian ready: shape={tuple(H.shape)}, active_atoms={_n_active}/{_n_atoms}, "
            f"frozen_atoms={_n_frozen}, active_dof={3 * _n_active}",
            detail=True,
        )
        freqs_cm, modes_mw = _frequencies_cm_and_modes(
            H,
            geometry.atomic_numbers,
            coords_bohr,
            device,
            freeze_idx=freeze_list if len(freeze_list) > 0 else None
        )

        del H
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        if freq_cfg["sort"] == "abs":
            order = np.argsort(np.abs(freqs_cm))
        else:
            order = np.argsort(freqs_cm)

        n_write = int(min(freq_cfg["max_write"], len(order)))
        _imag = [f for f in freqs_cm if f < 0]
        _preview_n = min(20, len(order))
        _freq_preview = ", ".join(f"{float(freqs_cm[j]):+.1f}" for j in order[:_preview_n])
        _suffix = ", ..." if len(order) > _preview_n else ""
        click.echo(
            f"[freq] {len(freqs_cm)} modes ({len(_imag)} imaginary); "
            f"first {_preview_n} by {freq_cfg['sort']}: [{_freq_preview}{_suffix}] cm⁻¹; "
            f"full list: {out_dir_path / 'frequencies_cm-1.txt'}",
            narrative=True,
        )
        if len(order) > _preview_n:
            _freq_str = ", ".join(f"{float(freqs_cm[j]):+.1f}" for j in order)
            click.echo(f"[freq:all] [{_freq_str}] cm⁻¹")
        click.echo(f"[INFO] Writing first {n_write} modes ({freq_cfg['sort']} ascending).", detail=True)

        # Reference PDB (only when input is PDB)
        ref_pdb = source_path if source_path.suffix.lower() == ".pdb" else None
        write_pdb = ref_pdb is not None

        # write modes
        for k, idx in enumerate(order[:n_write], start=1):
            freq = float(freqs_cm[idx])
            mode_cart_3N = _mw_mode_to_cart(modes_mw[idx], masses_au_t)
            out_trj = out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1_trj.xyz"
            _write_mode_trj_and_pdb(
                geometry,
                mode_cart_3N,
                out_trj,
                amplitude_ang=freq_cfg["amplitude_ang"],
                n_frames=freq_cfg["n_frames"],
                comment=f"mode {k}  {freq:+.2f} cm-1",
                ref_pdb=ref_pdb,
                write_pdb=write_pdb,
                prepared_input=prepared_input,
                out_pdb=out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1.pdb" if write_pdb else None,
            )
        (out_dir_path / "frequencies_cm-1.txt").write_text(
            "\n".join(f"{i+1:4d}  {float(freqs_cm[j]):+12.4f}" for i, j in enumerate(order)),
            encoding="utf-8"
        )

        del modes_mw
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        _thermo_data = None  # populated below if thermoanalysis succeeds
        try:
            from thermoanalysis.QCData import QCData
            from thermoanalysis.thermo import thermochemistry
            from thermoanalysis.constants import J2AU, NA, J2CAL

            qc_data = {
                "coords3d": geometry.cart_coords.reshape(-1, 3) * BOHR2ANG,  # Å
                "wavenumbers": freqs_cm,                                 # cm^-1
                "scf_energy": _calc_energy(geometry, calc_cfg),          # Hartree
                "masses": masses_amu,
                "mult": int(calc_cfg["spin"]),
            }
            qc = QCData(qc_data, point_group="c1", mult=int(calc_cfg["spin"]))

            T = float(thermo_cfg["temperature"])
            p_atm = float(thermo_cfg["pressure_atm"])
            p_pa = p_atm * 101325.0  # Pa

            tr = thermochemistry(qc, T, pressure=p_pa)  # default: QRRHO

            au2CalMol = (1.0 / J2AU) * NA * J2CAL
            to_cal_per_mol = lambda x: float(x) * au2CalMol
            J_per_Kmol_to_cal_per_Kmol = lambda j: float(j) * J2CAL

            n_imag = int(np.sum(freqs_cm < 0.0))

            EE = float(tr.U_el)
            ZPE = float(tr.ZPE)
            dE_therm = float(tr.U_therm)
            dH_therm = float(tr.H - tr.U_el)
            dG_therm = float(tr.dG)

            sum_EE_ZPE = EE + ZPE
            sum_EE_thermal_E = float(tr.U_tot)
            sum_EE_thermal_H = float(tr.H)
            sum_EE_thermal_G = float(tr.G)

            E_thermal_cal = to_cal_per_mol(tr.U_therm)
            Cv_cal_per_Kmol = J_per_Kmol_to_cal_per_Kmol(tr.c_tot)
            S_cal_per_Kmol  = to_cal_per_mol(tr.S_tot)

            click.echo("\n====== Thermochemistry summary ======\n", narrative=True)
            click.echo(f"Temperature (K)         = {T:.2f}")
            click.echo(f"Pressure    (atm)       = {p_atm:.4f}")
            if freeze_list:
                click.echo("[NOTE] Thermochemistry uses active DOF (PHVA) due to frozen atoms.", narrative=True)
            click.echo(f"Number of Imaginary Freq = {n_imag:d}")
            click.echo("")

            click.echo(f"Electronic Energy (EE)                 = {_fmt_ha(EE)}")
            click.echo(f"Zero-point Energy Correction           = {_fmt_ha(ZPE)}")
            click.echo(f"Thermal Correction to Energy           = {_fmt_ha(dE_therm)}")
            click.echo(f"Thermal Correction to Enthalpy         = {_fmt_ha(dH_therm)}")
            click.echo(f"Thermal Correction to Free Energy      = {_fmt_ha(dG_therm)}")
            click.echo(f"EE + Zero-point Energy                 = {_fmt_ha(sum_EE_ZPE)}")
            click.echo(f"EE + Thermal Energy Correction         = {_fmt_ha(sum_EE_thermal_E)}")
            click.echo(f"EE + Thermal Enthalpy Correction       = {_fmt_ha(sum_EE_thermal_H)}")
            click.echo(f"EE + Thermal Free Energy Correction    = {_fmt_ha(sum_EE_thermal_G)}")
            click.echo("")
            click.echo(f"E (Thermal)                            = {_fmt_cal(E_thermal_cal)}")
            click.echo(f"Heat Capacity (Cv)                     = {_fmt_calK(Cv_cal_per_Kmol)}")
            click.echo(f"Entropy (S)                            = {_fmt_calK(S_cal_per_Kmol)}")

            if bool(thermo_cfg["dump"]):
                out_yaml = out_dir_path / "thermoanalysis.yaml"
                payload = {
                    "temperature_K": T,
                    "pressure_atm": p_atm,
                    "num_imag_freq": n_imag,
                    "electronic_energy_ha": EE,
                    "zpe_correction_ha": ZPE,
                    "thermal_correction_energy_ha": dE_therm,
                    "thermal_correction_enthalpy_ha": dH_therm,
                    "thermal_correction_free_energy_ha": dG_therm,
                    "sum_EE_and_ZPE_ha": sum_EE_ZPE,
                    "sum_EE_and_thermal_energy_ha": sum_EE_thermal_E,
                    "sum_EE_and_thermal_enthalpy_ha": sum_EE_thermal_H,
                    "sum_EE_and_thermal_free_energy_ha": sum_EE_thermal_G,
                    "E_thermal_cal_per_mol": E_thermal_cal,
                    "Cv_cal_per_mol_K": Cv_cal_per_Kmol,
                    "S_cal_per_mol_K": S_cal_per_Kmol,
                }
                with out_yaml.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)
                click.echo(f"[dump] Wrote thermoanalysis summary → {out_yaml}", detail=True)

            _thermo_data = {
                "electronic_energy_ha": EE,
                "zpe_correction_ha": ZPE,
                "thermal_correction_energy_ha": dE_therm,
                "thermal_correction_enthalpy_ha": dH_therm,
                "thermal_correction_free_energy_ha": dG_therm,
                "sum_EE_and_ZPE_ha": sum_EE_ZPE,
                "sum_EE_and_thermal_energy_ha": sum_EE_thermal_E,
                "sum_EE_and_thermal_enthalpy_ha": sum_EE_thermal_H,
                "sum_EE_and_thermal_free_energy_ha": sum_EE_thermal_G,
                "E_thermal_cal_per_mol": E_thermal_cal,
                "Cv_cal_per_mol_K": Cv_cal_per_Kmol,
                "S_cal_per_mol_K": S_cal_per_Kmol,
            }

        except ImportError:
            click.echo("[thermo] WARNING: 'thermoanalysis' package not found; skipped thermochemistry summary.", err=True)
        except Exception as e:
            import traceback
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            click.echo("Unhandled error during thermochemistry summary:\n" + textwrap.indent(tb, "  "), err=True)

        # summary.md and key_* outputs are disabled.
        click.echo(f"[DONE] Wrote modes and list → {out_dir_path}", detail=True)

        click.echo(format_elapsed("[time] Elapsed Time for Freq", time_start), narrative=True)

        # result.json (if --out-json)
        if out_json:
            from pdb2reaction.core.utils import write_result_json
            _all_freqs = [float(f) for f in freqs_cm]
            _imag_freqs = [f for f in _all_freqs if f < 0.0]
            result_data = {
                "status": "completed",
                "n_modes": len(_all_freqs),
                "n_imaginary": len(_imag_freqs),
                "frequencies_cm": _all_freqs,
                "imaginary_frequencies_cm": _imag_freqs,
                "thermochemistry": _thermo_data,
                "backend": calc_cfg.get("backend", backend),
                "charge": calc_cfg["charge"],
                "spin": calc_cfg["spin"],
                "model": calc_cfg.get("model"),
                "n_atoms": len(geometry.atomic_numbers),
                "n_freeze_atoms": len(freeze_list) if 'freeze_list' in dir() else 0,
                "solvent": calc_cfg.get("solvent", "none"),
                "temperature_K": thermo_cfg.get("temperature", 298.15) if 'thermo_cfg' in dir() else 298.15,
                "pressure_atm": thermo_cfg.get("pressure_atm", 1.0) if 'thermo_cfg' in dir() else 1.0,
                "input_file": str(input_path),
                "files": {
                    "frequencies_txt": "frequencies_cm-1.txt",
                },
            }
            if _thermo_data is not None and bool(thermo_cfg.get("dump", False)):
                result_data["files"]["thermoanalysis_yaml"] = "thermoanalysis.yaml"
            write_result_json(
                out_dir_path, result_data,
                command="freq",
                elapsed_seconds=time.perf_counter() - time_start,
            )

    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        render_cli_exception(e, label="frequency analysis", out_dir=out_dir_path, command="freq", time_start=time_start)
    finally:
        prepared_input.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM:
        # rebind to None first (decrements heavy torch.nn.Module refs),
        # then del to drop the local names; gc.collect breaks cycles.
        calc = geometry = H = modes = None
        del calc, geometry, H, modes
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
