# pdb2reaction/freq.py

"""
freq — Vibrational frequency analysis, mode export, and thermochemistry
====================================================================

Usage (CLI)
-----------
    pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q <charge>] [--ligand-charge <number|'RES:Q,...'>] [-m <spin>] \
        [--workers <int>] [--workers-per-node <int>] \
        [--freeze-links {True|False}] [--max-write <int>] \
        [--amplitude-ang <float>] [--n-frames <int>] [--sort {value|abs}] \
        [--out-dir <dir>] [--args-yaml <file>] [--temperature <K>] \
        [--pressure <atm>] [--dump {True|False}] \
        [--hessian-calc-mode {Analytical|FiniteDifference}] \
        [--convert-files {True|False}] [--ref-pdb <file>]

Examples
--------
    # Minimal frequency run with explicit charge and spin
    pdb2reaction freq -i a.pdb -q 0 -m 1

    # PHVA with YAML overrides and a custom output directory
    pdb2reaction freq -i a.xyz -q -1 --args-yaml ./args.yaml --out-dir ./result_freq/

Description
-----------
- Computes vibrational frequencies and normal modes using the UMA calculator.
- Supports Partial Hessian Vibrational Analysis (PHVA) when atoms are frozen.
- Exports animated modes (.trj; and, for PDB inputs, optionally .pdb animations when ``--convert-files`` is enabled) and prints a Gaussian-style thermochemistry summary.
- For XYZ/GJF inputs, ``--ref-pdb`` supplies a reference PDB topology while keeping XYZ coordinates, enabling
  format-aware PDB/GJF output conversion.
- Configuration can be provided via YAML (sections: ``geom``, ``calc``, ``freq``, ``thermo``); YAML values override CLI.
- Thermochemistry uses the PHVA frequencies (respecting ``freeze_atoms``). CLI pressure (atm) is converted internally to Pa.
- The thermochemistry summary is printed when the optional ``thermoanalysis`` package is available; writing a YAML summary is controlled by ``--dump``.
- `-q/--charge` is required for non-`.gjf` inputs **unless** ``--ligand-charge`` is provided; when ``-q`` is omitted and
  ``--ligand-charge`` is set, the full complex is treated as an enzyme–substrate system and the total charge is inferred
  using ``extract.py``’s residue-aware logic. `.gjf` templates supply charge/spin when available and allow omitting
  `-q/--charge`. Explicit `-q` always overrides derived values.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_freq/)
  ├─ mode_XXXX_{±freq}cm-1.trj    # XYZ-like trajectory (Å) with sinusoidal motion per mode
  ├─ mode_XXXX_{±freq}cm-1.pdb    # Multi-MODEL PDB animation (PDB inputs only; attempted when --convert-files is enabled)
  ├─ frequencies_cm-1.txt         # All computed frequencies (cm^-1), sorted according to --sort
  └─ thermoanalysis.yaml          # Thermochemistry summary (written only when --dump True and thermoanalysis is available)

Notes
-----
- Input geometry: .pdb, .xyz, .trj (via pysisyphus ``geom_loader``). For PDB, ``--freeze-links`` (default: True)
  detects parent atoms of link hydrogens and merges them with any ``geom.freeze_atoms``; the merged list is echoed.
- UMA settings:
  - ``--hessian-calc-mode``: ``FiniteDifference`` (default) or ``Analytical``; may also be set in YAML.
  - ``return_partial_hessian`` defaults to True so PHVA can use a reduced active‑block Hessian when possible.
  - Device: ``auto`` selects CUDA if available; otherwise CPU.
  - Parallelism: ``--workers``/``--workers-per-node`` configure UMA predictors; when ``workers>1``,
    the parallel predictor disables analytical Hessians and always uses finite differences.
- PHVA and projections:
  - With frozen atoms, eigenanalysis is performed in the active DOF subspace and translation/rotation (TR) modes are projected in that subspace.
  - Both full Hessians (3N×3N) and pre‑reduced active blocks (3N_act×3N_act) are accepted.
  - Frequencies are reported in cm^-1 (negative values denote imaginary modes).
- Mode writing:
  - ``--max-write`` limits how many modes are exported (ascending by value, or by absolute value with ``--sort abs``).
  - ``--amplitude-ang`` (Å) and ``--n-frames`` control the sinusoidal animation.
- For PDB inputs and when ``--convert-files`` is enabled, the .trj is converted to a .pdb animation using the input PDB as a template.
- Thermochemistry:
  - Requires the optional ``thermoanalysis`` package; if absent, the summary is skipped with a warning.
  - Default model is QRRHO; the summary includes EE, ZPE, and thermal corrections to E/H/G. Values in cal·mol^-1 and cal·(mol·K)^-1 are also printed.
- Performance and numerical details:
  - GPU memory usage is minimized by keeping only one Hessian in memory and avoiding redundant allocations.
- Exit behavior:
  - On keyboard interrupt, exits with code 130; on errors during frequency analysis, prints a traceback and exits with code 1.
    Errors during the optional thermochemistry summary are reported but do not abort the run.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List

import click
import numpy as np
import torch
from ase.data import atomic_masses
import ase.units as units
import yaml
import time

# ---------------- pysisyphus / pdb2reaction imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AMU2AU, AU2EV

# local helpers from pdb2reaction
from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .utils import (
    load_yaml_dict,
    apply_yaml_overrides,
    convert_xyz_like_outputs,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
    resolve_freeze_atoms,
)


def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)


# ===================================================================
#          Mass-weighted TR projection & vibrational analysis
# ===================================================================

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

        H = 0.5 * (H + H.T)

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
        masses_amu = np.array([atomic_masses[z] for z in Z])
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
                Hmw_act = 0.5 * (Hmw_act + Hmw_act.T)

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
                H = 0.5 * (H + H.T)

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


def _calc_full_hessian_torch(geom, uma_kwargs: dict, device: torch.device) -> torch.Tensor:
    """
    UMA calculator producing Hessian as torch.Tensor in Hartree/Bohr^2 (3N or 3N_act square).
    """
    kw = dict(uma_kwargs or {})
    kw["out_hess_torch"] = True
    calc = uma_pysis(**kw)
    H = calc.get_hessian(geom.atoms, geom.cart_coords)["hessian"].to(device=device)
    return H


def _calc_energy(geom, uma_kwargs: dict) -> float:
    """
    Compute electronic energy (Hartree) from UMA calculator.
    """
    calc = uma_pysis(**uma_kwargs)
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
    Write a single mode animation as .trj (XYZ-like) and optionally .pdb.

    If `ref_pdb` is provided and is a .pdb file, the .pdb is generated by
    converting the .trj using the input PDB as the template.
    Set `write_pdb=False` to skip PDB generation.
    """
    ref_ang = geom.cart_coords.reshape(-1, 3) * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # .trj (concatenated XYZ-like trajectory)
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


# ===================================================================
#                         Defaults for CLI
# ===================================================================

# Geometry defaults
GEOM_KW = dict(GEOM_KW_DEFAULT)

CALC_KW = dict(_UMA_CALC_KW)
CALC_KW.update(
    {
        "return_partial_hessian": True,  # default to PHVA-friendly active-block Hessians when possible
    }
)

# Freq writer defaults
FREQ_KW = {
    "amplitude_ang": 0.8,     # animation amplitude (Å) applied to both .trj and .pdb outputs
    "n_frames": 20,           # number of frames per vibrational mode
    "max_write": 10,          # maximum number of modes to export
    "sort": "value",          # "value" (ascending by cm^-1) | "abs" (ascending by absolute value)
}

# Thermochemistry defaults
THERMO_KW = {
    "temperature": 298.15,    # temperature in Kelvin
    "pressure_atm": 1.0,      # pressure in atm (converted to Pa internally)
    "dump": False,            # write thermoanalysis.yaml when True
}


# ===================================================================
#                            CLI
# ===================================================================

@click.command(
    help="Vibrational frequency analysis and mode writer (+ default thermochemistry summary).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure (.pdb, .xyz, .trj, ...)",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=False,
    help=(
        "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
        "(PDB inputs or XYZ/GJF with --ref-pdb)."
    ),
)
@click.option(
    "--workers",
    type=int,
    default=CALC_KW["workers"],
    show_default=True,
    help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel UMA predictor (workers>1).",
)
@click.option(
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
        "when -q is omitted (requires PDB input or --ref-pdb)."
    ),
)
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1) for the ML region.")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze parent atoms of link hydrogens (PDB only).")
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
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
@click.option("--max-write", type=int, default=10, show_default=True,
              help="How many modes to export (after sorting per --sort).")
@click.option("--amplitude-ang", type=float, default=0.8, show_default=True,
              help="Animation amplitude (Å) used for both .trj and .pdb.")
@click.option("--n-frames", type=int, default=20, show_default=True,
              help="Number of frames per mode animation.")
@click.option("--sort", type=click.Choice(["value", "abs"]), default="value", show_default=True,
              help="Sort modes by 'value' (cm^-1) or by absolute value.")
@click.option("--out-dir", type=str, default="./result_freq/", show_default=True, help="Output directory")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, freq, thermo).",
)
# Thermochemistry options
@click.option("--temperature", type=float, default=THERMO_KW["temperature"], show_default=True,
              help="Temperature (K) for thermochemistry summary.")
@click.option("--pressure", "pressure_atm",
              type=float, default=THERMO_KW["pressure_atm"], show_default=True,
              help="Pressure (atm) for thermochemistry summary.")
@click.option("--dump", type=click.BOOL, default=THERMO_KW["dump"], show_default=True,
              help="When True, write 'thermoanalysis.yaml' under out-dir.")
# Hessian calculation mode
@click.option("--hessian-calc-mode",
              type=click.Choice(["FiniteDifference", "Analytical"]),
              default=None,
              help="How UMA computes Hessian. Defaults to 'FiniteDifference' (can also be set via YAML).")
def cli(
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    max_write: int,
    amplitude_ang: float,
    n_frames: int,
    sort: str,
    out_dir: str,
    args_yaml: Optional[Path],
    # thermo
    temperature: float,
    pressure_atm: float,
    dump: bool,
    # hessian
    hessian_calc_mode: Optional[str],
) -> None:
    time_start = time.perf_counter()
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    apply_ref_pdb_override(prepared_input, ref_pdb)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input,
        charge,
        spin,
        ligand_charge=ligand_charge,
        prefix="[freq]",
    )

    # --------------------------
    # 1) Assemble configuration (defaults ← CLI ← YAML)
    # --------------------------
    yaml_cfg = load_yaml_dict(args_yaml)
    geom_cfg = dict(GEOM_KW)
    calc_cfg = dict(CALC_KW)
    freq_cfg = dict(FREQ_KW)
    thermo_cfg = dict(THERMO_KW)

    # CLI overrides
    calc_cfg["charge"] = int(charge)
    calc_cfg["spin"]   = int(spin)
    calc_cfg["workers"] = int(workers)
    calc_cfg["workers_per_node"] = int(workers_per_node)
    if hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

    freq_cfg["max_write"]     = int(max_write)
    freq_cfg["amplitude_ang"] = float(amplitude_ang)
    freq_cfg["n_frames"]      = int(n_frames)
    freq_cfg["sort"]          = sort

    thermo_cfg["temperature"]   = float(temperature)
    thermo_cfg["pressure_atm"]  = float(pressure_atm)
    thermo_cfg["dump"]          = bool(dump)

    # YAML overrides (highest precedence)
    apply_yaml_overrides(
        yaml_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (freq_cfg, (("freq",),)),
        ],
    )
    thermo_yaml = yaml_cfg.get("thermo")
    thermo_yaml_dict = thermo_yaml if isinstance(thermo_yaml, dict) else None

    # Normalize freeze_atoms and optionally add link-parent indices for PDB inputs
    resolve_freeze_atoms(geom_cfg, source_path, freeze_links, on_error="warn")

    # Ensure calc config reflects the geometry freeze list used in the run.
    calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
    calc_cfg.setdefault("return_partial_hessian", True)

    out_dir_path = Path(out_dir).resolve()
    out_dir_path.mkdir(parents=True, exist_ok=True)

    # Pretty-print config summary
    click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
    click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
    click.echo(pretty_block("freq", {**freq_cfg, "out_dir": str(out_dir_path)}))
    thermo_block = {
        "temperature": thermo_cfg["temperature"],
        "pressure_atm": thermo_cfg["pressure_atm"],
        "dump": thermo_cfg["dump"],
    }
    if thermo_yaml_dict:
        thermo_block["note"] = "args-yaml thermo section is not applied by this command."
    click.echo(pretty_block("thermo", thermo_block))
    if thermo_yaml_dict:
        click.echo(pretty_block("thermo (args-yaml; ignored)", thermo_yaml_dict))

    # --------------------------
    # 2) Load geometry
    # --------------------------
    coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
    coord_kwargs = dict(geom_cfg)
    coord_kwargs.pop("coord_type", None)
    geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

    # Masses (AU tensor for TR projection & MW->Cart conversion)
    masses_amu = np.array([atomic_masses[z] for z in geometry.atomic_numbers])
    masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=torch.float32)
    device = _torch_device(calc_cfg.get("device", "auto"))
    masses_au_t = masses_au_t.to(device=device)

    # --------------------------
    # 3) Compute Hessian & modes
    # --------------------------
    try:
        freeze_list = list(calc_cfg.get("freeze_atoms", []))
        H = _calc_full_hessian_torch(geometry, calc_cfg, device)
        coords_bohr = geometry.cart_coords.reshape(-1, 3)

        # PHVA: use the freeze list to carve out the active subspace and apply TR projection there.
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
        click.echo(f"[INFO] Total modes: {len(freqs_cm)}  → write first {n_write} modes ({freq_cfg['sort']} ascending).")

        # Reference PDB (only when input is PDB)
        ref_pdb = source_path if source_path.suffix.lower() == ".pdb" else None
        write_pdb = ref_pdb is not None

        # write modes
        for k, idx in enumerate(order[:n_write], start=1):
            freq = float(freqs_cm[idx])
            mode_cart_3N = _mw_mode_to_cart(modes_mw[idx], masses_au_t)
            out_trj = out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1.trj"
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

        # --------------------------
        # 4) Thermochemistry summary
        # --------------------------
        try:
            from thermoanalysis.QCData import QCData
            from thermoanalysis.thermo import thermochemistry
            from thermoanalysis.constants import J2AU, NA, J2CAL

            qc_data = {
                "coords3d": geometry.cart_coords.reshape(-1, 3) * BOHR2ANG,  # Å
                "wavenumbers": freqs_cm,                                 # cm^-1
                "scf_energy": _calc_energy(geometry, calc_cfg),          # Hartree
                "masses": masses_amu,
                "mult": int(spin),
            }
            qc = QCData(qc_data, point_group="c1", mult=int(spin))

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

            click.echo("\nThermochemistry Summary")
            click.echo("------------------------")
            click.echo(f"Temperature (K)         = {T:.2f}")
            click.echo(f"Pressure    (atm)       = {p_atm:.4f}")
            if freeze_list:
                click.echo("[NOTE] Thermochemistry uses active DOF (PHVA) due to frozen atoms.")
            click.echo(f"Number of Imaginary Freq = {n_imag:d}\n")

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
            click.echo("")

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
                click.echo(f"[dump] Wrote thermoanalysis summary → {out_yaml}")

        except ImportError:
            click.echo("[thermo] WARNING: 'thermoanalysis' package not found; skipped thermochemistry summary.", err=True)
        except Exception as e:
            import traceback
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            click.echo("Unhandled error during thermochemistry summary:\n" + textwrap.indent(tb, "  "), err=True)

        click.echo(f"[DONE] Wrote modes and list → {out_dir_path}")

        click.echo(format_elapsed("[time] Elapsed Time for Freq", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        import traceback
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during frequency analysis:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()
