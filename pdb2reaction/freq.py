# pdb2reaction/freq.py

"""
Vibrational frequency analysis CLI for pdb2reaction.

Usage
-----
pdb2reaction freq \
  -i a.pdb -q 0 -s 1 \
  --max-write 20 \
  --out-dir ./result_freq/ \
  --args-yaml ./args.yaml \
  --temperature 298.15 \
  --pressure 1.0 \
  --dump False

Sections that can be overridden via YAML:
  geom, calc, freq

Output:
- In out_dir, write up to max-write modes in ascending order of frequency
  as 'mode_XXXX_{±freq}cm-1.trj' and '.pdb'.

Thermochemistry (default on):
- Prints a Gaussian-style summary at the end.
- When --dump True, writes 'thermoanalysis.yaml' under out_dir.

Notes:
- Thermochemistry uses frequencies computed here (PHVA respecting freeze_atoms).
- Pressure is specified in atm at CLI and converted to Pa internally.
"""

from __future__ import annotations

import sys
import textwrap
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

# ---------------- pysisyphus / pdb2reaction imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AMU2AU, AU2EV

# local helpers from pdb2reaction
from .uma_pysis import uma_pysis
from .utils import freeze_links as _freeze_links_from_utils, convert_xyz_to_pdb as _convert_xyz_to_pdb


# ===================================================================
#                         Generic helpers
# ===================================================================

def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively update dict *dst* with *src*, returning *dst*."""
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")
    return data


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if fa else ""
    return g


def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)


# ===================================================================
#          Mass-weighted TR projection & vibrational analysis
# ===================================================================

def _build_tr_basis(coords_bohr_t: torch.Tensor,
                    masses_au_t: torch.Tensor) -> torch.Tensor:
    """Mass-weighted translation/rotation basis (Tx, Ty, Tz, Rx, Ry, Rz), shape (3N, r<=6)."""
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
    """Orthonormalize TR basis in mass-weighted space by SVD. Returns (Q, rank)."""
    B = _build_tr_basis(coords_bohr_t, masses_au_t)
    U, S, Vh = torch.linalg.svd(B, full_matrices=False)
    r = int((S > rtol * S.max()).sum().item())
    Q = U[:, :r]
    del B, S, Vh, U
    return Q, r


def _mw_projected_hessian(H_t: torch.Tensor,
                          coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Project out translations/rotations in mass-weighted space:
    Hmw = M^{-1/2} H M^{-1/2};  P = I - QQ^T;  Hmw_proj = P Hmw P

    To save memory, update **H_t in-place** (no clone) and return it.
    No explicit symmetrization (0.5*(H+H^T)). The eigen-decomposition uses
    only the lower triangle (UPLO="L").
    """
    dtype, device = H_t.dtype, H_t.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)

        # In-place mass-weighting on input Hessian
        H_t.mul_(inv_sqrt_m_row)
        H_t.mul_(inv_sqrt_m_col)

        Q, _ = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)  # (3N, r)
        Qt = Q.T

        QtH = Qt @ H_t
        H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)

        HQ = QtH.T
        H_t.addmm_(HQ, Qt, beta=1.0, alpha=-1.0)

        QtHQ = QtH @ Q
        tmp = Q @ QtHQ
        H_t.addmm_(tmp, Qt, beta=1.0, alpha=1.0)

        del masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row
        del Q, Qt, QtH, HQ, QtHQ, tmp

        if torch.cuda.is_available() and device.type == "cuda":
            torch.cuda.empty_cache()
        # Return the in-place updated mass-weighted & TR-projected Hessian
        return H_t


# ---- PHVA helper: mass-weighted Hessian without TR projection (for active subspace) ----
def _mass_weighted_hessian(H_t: torch.Tensor,
                           masses_au_t: torch.Tensor) -> torch.Tensor:
    """Return Hmw = M^{-1/2} H M^{-1/2} (no symmetrization/TR projection; in-place)."""
    dtype, device = H_t.dtype, H_t.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)
        # In-place mass-weighting on input Hessian
        H_t.mul_(inv_sqrt_m_row)
        H_t.mul_(inv_sqrt_m_col)
        del masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row
        return H_t


def _frequencies_cm_and_modes(H_t: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              tol: float = 1e-6,
                              freeze_idx: Optional[List[int]] = None) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize a (possibly PHVA/active-subspace) TR-projected mass-weighted Hessian
    to obtain frequencies (cm^-1) and mass-weighted eigenvectors (modes).

    If `freeze_idx` is provided (list of 0-based atom indices), perform
    Partial Hessian Vibrational Analysis (PHVA):
      1) build Hmw = M^{-1/2} H M^{-1/2}
      2) take the active subspace by removing DOF of frozen atoms
      3) perform TR projection **only in the active subspace** (always applied; works for 0/1/2 frozen atoms)
      4) diagonalize and embed eigenvectors back to 3N by zero-filling frozen DOF

    Returns:
      freqs_cm : (nmode,) numpy, negatives are imaginary
      modes    : (nmode, 3N) torch (mass-weighted eigenvectors)
    """
    with torch.no_grad():
        Z = np.array(atomic_numbers, dtype=int)
        N = int(len(Z))
        masses_amu = np.array([atomic_masses[z] for z in Z])  # amu
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        # --------------------------------------------
        # PHVA path (active DOF subspace with TR-proj)
        # --------------------------------------------
        if freeze_idx is not None and len(freeze_idx) > 0:
            # Active atom indices
            frozen_set = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen_set]

            n_active = len(active_idx)
            if n_active == 0:
                # All atoms are frozen → no modes
                freqs_cm = np.zeros((0,), dtype=float)
                modes = torch.zeros((0, 3 * N), dtype=H_t.dtype, device=H_t.device)
                return freqs_cm, modes

            # Mass-weighted Hessian (full, in-place) → active-DOF submatrix
            H_t = _mass_weighted_hessian(H_t, masses_au_t)
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=H_t.device)
            for i in frozen_set:
                mask_dof[3 * i:3 * i + 3] = False
            H_t = H_t[mask_dof][:, mask_dof]
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            # Active-only TR basis and projection in mass-weighted space (always applied; in-place)
            coords_act = coords_bohr_t[active_idx, :]
            masses_act = masses_au_t[active_idx]
            Q, _ = _tr_orthonormal_basis(coords_act, masses_act)  # (3N_act, r)
            Qt = Q.T

            QtH = Qt @ H_t
            H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)

            H_t.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)

            QtH = QtH @ Q
            H_t.addmm_(Q @ QtH, Qt, beta=1.0, alpha=1.0)

            del Q, Qt, QtH

            # No symmetrization; diagonalize using lower triangle only
            omega2, Vsub = torch.linalg.eigh(H_t, UPLO="L")

            sel = torch.abs(omega2) > tol
            omega2 = omega2[sel]
            Vsub = Vsub[:, sel]  # (3N_act, nsel)

            modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=H_t.dtype, device=H_t.device)
            modes[:, mask_dof] = Vsub.T  # (nsel, 3N_act) → place into active DOF
            del Vsub, mask_dof

        else:
            # Legacy behavior: TR-projection in full DOF → diagonalization (both in-place; no symmetrization)
            H_t = _mw_projected_hessian(H_t, coords_bohr_t, masses_au_t)
            omega2, V = torch.linalg.eigh(H_t, UPLO="L")

            sel = torch.abs(omega2) > tol
            omega2 = omega2[sel]
            modes = V[:, sel].T
            del V

        # Convert to frequencies (cm^-1)
        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()

        del omega2, hnu, sel
        if torch.cuda.is_available() and H_t.is_cuda:
            torch.cuda.empty_cache()
        return freqs_cm, modes


def _mw_mode_to_cart(mode_mw_3N_t: torch.Tensor,
                     masses_au_t: torch.Tensor) -> np.ndarray:
    """Convert one mass-weighted eigenvector (3N,) to Cartesian (3N,) and L2-normalize."""
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=mode_mw_3N_t.dtype, device=mode_mw_3N_t.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        v_cart = torch.sqrt(1.0 / m3) * mode_mw_3N_t
        v_cart.div_(torch.linalg.norm(v_cart))
        arr = v_cart.detach().cpu().numpy()
        del masses_amu_t, m3, v_cart
        return arr


def _calc_full_hessian_torch(geom, uma_kwargs: dict, device: torch.device) -> torch.Tensor:
    """UMA calculator producing analytic Hessian as torch.Tensor in Hartree/Bohr^2 (3N,3N)."""
    kw = dict(uma_kwargs or {})
    kw["out_hess_torch"] = True
    calc = uma_pysis(**kw)
    H_t = calc.get_hessian(geom.atoms, geom.coords)["hessian"].to(device=device)
    return H_t


def _calc_energy(geom, uma_kwargs: dict) -> float:
    """Compute electronic energy (Hartree) from UMA calculator."""
    calc = uma_pysis(out_hess_torch=False, **uma_kwargs)
    geom.set_calculator(calc)
    E = float(geom.energy)
    geom.set_calculator(None)
    return E


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            out_pdb: Path,
                            amplitude_ang: float = 0.8,
                            n_frames: int = 20,
                            comment: str = "mode",
                            ref_pdb: Optional[Path] = None) -> None:
    """Write a single mode animation as .trj (XYZ-like) and .pdb.

    If `ref_pdb` is provided and is a .pdb file, the .pdb is generated by
    converting the .trj using the input PDB as the template (same as path_opt).
    """
    ref_ang = geom.coords.reshape(-1, 3) * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # .trj (concatenated XYZ-like trajectory)
    if ref_pdb is not None and ref_pdb.suffix.lower() == ".pdb":
        # Emit a simple XYZ-like trajectory in Å for the converter
        with out_trj.open("w", encoding="utf-8") as f:
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                coords = ref_ang + phase * amplitude_ang * mode  # Å
                f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
                for sym, (x, y, z) in zip(geom.atoms, coords):
                    f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")
        # Generate PDB using the input PDB as template
        try:
            _convert_xyz_to_pdb(out_trj, ref_pdb, out_pdb)
        except Exception:
            # Fallback: generate MODEL/ENDMDL using ASE
            atoms0 = Atoms(geom.atoms, positions=ref_ang, pbc=False)
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                ai = atoms0.copy()
                ai.set_positions(ref_ang + phase * amplitude_ang * mode)
                write(out_pdb, ai, append=(i != 0))
        return

    # If no ref_pdb is given, use the legacy behavior (use pysisyphus.make_trj_str if available)
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

    # .pdb (MODEL/ENDMDL via ASE)
    atoms0 = Atoms(geom.atoms, positions=ref_ang, pbc=False)
    for i in range(n_frames):
        phase = np.sin(2.0 * np.pi * i / n_frames)
        ai = atoms0.copy()
        ai.set_positions(ref_ang + phase * amplitude_ang * mode)
        write(out_pdb, ai, append=(i != 0))


# ===================================================================
#                         Defaults for CLI
# ===================================================================

# Geometry defaults
GEOM_KW = {
    "coord_type": "cart",     # str, coordinate representation for geom_loader
    "freeze_atoms": [],       # List[int], 0-based indices to freeze
}

# UMA calculator defaults
CALC_KW = {
    "charge": 0,              # int
    "spin": 1,                # int multiplicity (2S+1)
    "model": "uma-s-1p1",     # str, UMA pretrained model ID
    "task_name": "omol",      # str
    "device": "auto",         # "cuda" | "cpu" | "auto"
    "max_neigh": None,        # Optional[int]
    "radius": None,           # Optional[float] (Å)
    "r_edges": False,         # bool
    "out_hess_torch": True,   # Required: return Hessian as a torch tensor
}

# Freq writer defaults
FREQ_KW = {
    "amplitude_ang": 0.8,     # float, animation amplitude (Å) for both .trj and .pdb
    "n_frames": 20,           # int, number of frames
    "max_write": 20,          # int, number of modes to write
    "sort": "value",          # "value" (ascending by value) | "abs" (ascending by absolute value)
}

# Thermochemistry defaults (added)
THERMO_KW = {
    "temperature": 298.15,    # float, Kelvin
    "pressure_atm": 1.0,      # float, atm (converted to Pa internally)
    "dump": False,            # bool, write thermoanalysis.yaml when True
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
@click.option("-q", "--charge", type=int, required=True, help="Total charge")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1)")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze parent atoms of link hydrogens (PDB only).")
@click.option("--max-write", type=int, default=20, show_default=True,
              help="Number of modes to write (ascending by frequency value).")
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
    help="YAML with extra args (sections: geom, calc, freq).",
)
# ---- Thermochemistry options (added) ----
@click.option("--temperature", type=float, default=THERMO_KW["temperature"], show_default=True,
              help="Temperature (K) for thermochemistry summary.")
@click.option("--pressure", type=float, default=THERMO_KW["pressure_atm"], show_default=True,
              help="Pressure (atm) for thermochemistry summary.")
@click.option("--dump", type=click.BOOL, default=THERMO_KW["dump"], show_default=True,
              help="When True, write 'thermoanalysis.yaml' under out-dir.")
def cli(
    input_path: Path,
    charge: int,
    spin: int,
    freeze_links: bool,
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
) -> None:
    time_start = time.perf_counter()
    
    # --------------------------
    # 1) Assemble configuration
    # --------------------------
    yaml_cfg = _load_yaml(args_yaml)
    geom_cfg = dict(GEOM_KW)
    calc_cfg = dict(CALC_KW)
    freq_cfg = dict(FREQ_KW)
    thermo_cfg = dict(THERMO_KW)

    _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
    _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
    _deep_update(freq_cfg, yaml_cfg.get("freq", {}))

    # CLI overrides
    calc_cfg["charge"] = int(charge)
    calc_cfg["spin"]   = int(spin)
    freq_cfg["max_write"]     = int(max_write)
    freq_cfg["amplitude_ang"] = float(amplitude_ang)
    freq_cfg["n_frames"]      = int(n_frames)
    freq_cfg["sort"]          = sort

    thermo_cfg["temperature"]   = float(temperature)
    thermo_cfg["pressure_atm"]  = float(pressure_atm)
    thermo_cfg["dump"]          = bool(dump)

    # Freeze links (PDB only): merge with existing list
    if freeze_links and input_path.suffix.lower() == ".pdb":
        try:
            detected = _freeze_links_from_utils(input_path)
        except Exception as e:
            click.echo(f"[freeze-links] WARNING: Could not detect link parents: {e}", err=True)
            detected = []
        base_freeze = list(geom_cfg.get("freeze_atoms", []))
        merged = sorted(set(base_freeze).union(set(detected)))
        geom_cfg["freeze_atoms"] = merged
        if merged:
            click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, merged))}")

    out_dir_path = Path(out_dir).resolve()
    out_dir_path.mkdir(parents=True, exist_ok=True)

    # Pretty-print config summary
    click.echo(_pretty_block("geom", _format_geom_for_echo(geom_cfg)))
    click.echo(_pretty_block("calc", calc_cfg))
    click.echo(_pretty_block("freq", {**freq_cfg, "out_dir": str(out_dir_path)}))
    click.echo(_pretty_block("thermo", {
        "temperature": thermo_cfg["temperature"],
        "pressure_atm": thermo_cfg["pressure_atm"],
        "dump": thermo_cfg["dump"],
    }))

    # --------------------------
    # 2) Load geometry
    # --------------------------
    coord_type = geom_cfg.get("coord_type", "cart")
    coord_kwargs = dict(geom_cfg)
    coord_kwargs.pop("coord_type", None)
    geometry = geom_loader(input_path, coord_type=coord_type, **coord_kwargs)

    # Masses (AU tensor for TR projection & MW->Cart conversion)
    masses_amu = np.array([atomic_masses[z] for z in geometry.atomic_numbers])
    masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=torch.float32)
    device = _torch_device(calc_cfg.get("device", "auto"))
    masses_au_t = masses_au_t.to(device=device)

    # --------------------------
    # 3) Compute Hessian & modes
    # --------------------------
    try:
        H_t = _calc_full_hessian_torch(geometry, calc_cfg, device)
        coords_bohr = geometry.coords.reshape(-1, 3)

        # PHVA: use the freeze list to carve out the active subspace and apply TR projection there
        freeze_list = list(geom_cfg.get("freeze_atoms", []))
        freqs_cm, modes_mw = _frequencies_cm_and_modes(
            H_t,
            geometry.atomic_numbers,
            coords_bohr,
            device,
            freeze_idx=freeze_list if len(freeze_list) > 0 else None
        )

        # sort order
        if freq_cfg["sort"] == "abs":
            order = np.argsort(np.abs(freqs_cm))
        else:
            order = np.argsort(freqs_cm)

        n_write = int(min(freq_cfg["max_write"], len(order)))
        click.echo(f"[INFO] Total modes: {len(freqs_cm)}  → write first {n_write} modes ({freq_cfg['sort']} ascending).")

        # Reference PDB (only when input is PDB)
        ref_pdb = input_path if input_path.suffix.lower() == ".pdb" else None

        # write modes
        for k, idx in enumerate(order[:n_write], start=1):
            freq = float(freqs_cm[idx])
            mode_cart_3N = _mw_mode_to_cart(modes_mw[idx], masses_au_t)  # (3N,)
            out_trj = out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1.trj"
            out_pdb = out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1.pdb"
            _write_mode_trj_and_pdb(
                geometry,
                mode_cart_3N,
                out_trj,
                out_pdb,
                amplitude_ang=freq_cfg["amplitude_ang"],
                n_frames=freq_cfg["n_frames"],
                comment=f"mode {k}  {freq:+.2f} cm-1",
                ref_pdb=ref_pdb,
            )
        # also write a simple list
        (out_dir_path / "frequencies_cm-1.txt").write_text(
            "\n".join(f"{i+1:4d}  {float(freqs_cm[j]):+12.4f}" for i, j in enumerate(order)),
            encoding="utf-8"
        )

        if torch.cuda.is_available() and H_t.is_cuda:
            torch.cuda.empty_cache()

        # --------------------------
        # 4) Thermochemistry summary (default on)
        # --------------------------
        try:
            # Lazy import so that freq-only users can still run without thermoanalysis installed
            from thermoanalysis.QCData import QCData
            from thermoanalysis.thermo import thermochemistry
            from thermoanalysis.constants import J2AU, NA, J2CAL

            # Prepare QCData dict
            qc_data = {
                "coords3d": geometry.coords.reshape(-1, 3) * BOHR2ANG,  # Å
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

            # Converters
            au2CalMol = (1.0 / J2AU) * NA * J2CAL
            to_cal_per_mol = lambda x: float(x) * au2CalMol
            J_per_Kmol_to_cal_per_Kmol = lambda j: float(j) * J2CAL

            # Counts
            n_imag = int(np.sum(freqs_cm < 0.0))

            # Compose summary
            EE = float(tr.U_el)
            ZPE = float(tr.ZPE)
            dE_therm = float(tr.U_therm)               # Thermal correction to Energy (includes ZPE)
            dH_therm = float(tr.H - tr.U_el)           # Thermal correction to Enthalpy (= U_therm + kBT)
            dG_therm = float(tr.dG)                    # Thermal correction to Free Energy (= G - EE)

            sum_EE_ZPE = EE + ZPE
            sum_EE_thermal_E = float(tr.U_tot)         # = EE + U_therm
            sum_EE_thermal_H = float(tr.H)             # = H
            sum_EE_thermal_G = float(tr.G)             # = G

            E_thermal_cal = to_cal_per_mol(tr.U_therm)               # cal/mol
            Cv_cal_per_Kmol = J_per_Kmol_to_cal_per_Kmol(tr.c_tot)   # cal/(mol*K)
            S_cal_per_Kmol  = to_cal_per_mol(tr.S_tot)               # cal/(mol*K)

            # Echo summary (Gaussian-like)
            click.echo("\nThermochemistry Summary")
            click.echo("------------------------")
            click.echo(f"Temperature (K)         = {T:.2f}")
            click.echo(f"Pressure    (atm)       = {p_atm:.4f}")
            if freeze_list:
                click.echo("[NOTE] Thermochemistry uses active DOF (PHVA) due to frozen atoms.")
            click.echo(f"Number of Imaginary Freq = {n_imag:d}\n")

            def _ha(x): return f"{float(x): .6f} Ha"
            def _cal(x): return f"{float(x): .2f} cal/mol"
            def _calK(x): return f"{float(x): .2f} cal/(mol*K)"

            click.echo(f"Electronic Energy (EE)                 = {_ha(EE)}")
            click.echo(f"Zero-point Energy Correction           = {_ha(ZPE)}")
            click.echo(f"Thermal Correction to Energy           = {_ha(dE_therm)}")
            click.echo(f"Thermal Correction to Enthalpy         = {_ha(dH_therm)}")
            click.echo(f"Thermal Correction to Free Energy      = {_ha(dG_therm)}")
            click.echo(f"EE + Zero-point Energy                 = {_ha(sum_EE_ZPE)}")
            click.echo(f"EE + Thermal Energy Correction         = {_ha(sum_EE_thermal_E)}")
            click.echo(f"EE + Thermal Enthalpy Correction       = {_ha(sum_EE_thermal_H)}")
            click.echo(f"EE + Thermal Free Energy Correction    = {_ha(sum_EE_thermal_G)}")
            click.echo("")
            click.echo(f"E (Thermal)                            = {_cal(E_thermal_cal)}")
            click.echo(f"Heat Capacity (Cv)                     = {_calK(Cv_cal_per_Kmol)}")
            click.echo(f"Entropy (S)                            = {_calK(S_cal_per_Kmol)}")
            click.echo("")

            # Dump YAML when requested
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

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for Freq: {hh:02d}:{mm:02d}:{ss:06.3f}")

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        import traceback
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during frequency analysis:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


# Allow `python -m pdb2reaction.freq` direct execution
if __name__ == "__main__":
    cli()
