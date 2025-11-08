# pdb2reaction/path_search.py
"""
Multistep MEP search via recursive GSM segmentation.

Overview
--------
Starting from two endpoint structures (reactant, product):
(1) Run one GSM between the endpoints to obtain an initial MEP.
(2) Optimize the single images immediately left/right of the highest-energy image (HEI) to get nearby minima (End1, End2).
(3) Run a second GSM between End1–End2 (refinement) to finalize one "step" of the MEP.
(4) Detect covalent bond changes for pairs (A, End1) and (End2, B). Recursively apply the same procedure to the side(s)
    where bond changes are present.
(5) Concatenate all "step" MEPs into a continuous reactant→product MEP, removing duplicates and inserting bridge GSMs as
    needed. If the interface between steps shows covalent changes, insert a *new* recursive segment instead of a bridge.

**Support for multiple input structures (new)**
----------------------------------------------
- You can pass two or more structures to `-i/--input` **in reaction order**
  (e.g., `reac.pdb im1.pdb im2.pdb prod.pdb`). Due to Click constraints, specify by repeating `-i`:
  `-i reac.pdb -i im1.pdb -i im2.pdb -i prod.pdb`.
- The above procedure is applied to each adjacent pair (A→I1, I1→I2, I2→B, ...). All per-pair paths are then
  concatenated (with duplicate removal and automatic bridge/recursive segment insertion) to build one MEP.
- Pre-optimization and Kabsch pre-alignment are performed **per adjacent pair** sequentially.

Key specifications
------------------
- A single UMA calculator (uma_pysis) is shared across all stages (serial execution).
- GSM uses pysisyphus GrowingString + StringOptimizer.
- The HEI neighbors are optimized as single structures (LBFGS / RFO selectable).
- Covalent changes are detected via `bond_changes.compare_structures`.
- When concatenating "step" MEPs:
  - If endpoint mismatch exceeds a threshold, run a bridge GSM automatically.
  - **If there are covalent changes between the previous step’s tail and the next step’s head,
    insert a *new* recursive segment instead of a bridge.**
- **Two input endpoints are pre-aligned using Kabsch (do not use StringOptimizer’s align).**
- **For systems with `freeze_atoms`, the Kabsch alignment is solved using only the freeze atoms; the resulting rigid transform
  is applied to all atoms.**
- **Each input structure is single-point optimized *before* alignment (sopt-mode).**
  **However, when `--pre_opt True` / `--pre-opt True` is passed, this pre-optimization is skipped.**
- Priority of YAML (geom, calc, gs, opt, sopt, bond, search) and CLI is the same as existing tools (CLI > YAML > defaults).
- **Separate `max_nodes` can be set for segments and bridges
  (`search.max_nodes_segment` / `search.max_nodes_bridge`).**
  Bridge GSM enforces `climb=False` (except for the “connective” GSM at max-depth, which keeps `climb=True`).
- **Kink detection**: If the maximum per-atom displacement between the End1/End2 optimized structures is
  ≤ `search.kink_disp_thresh` (default 0.2 Å), treat the vicinity as a “kink” and **skip GSM**.
  Instead, generate `search.kink_max_nodes` (default 3) linear interpolation structures, then single-structure optimize
  them (using common `freeze_atoms`) and insert them as images.

Output
------
out_dir/
  ├─ summary.yaml                 : run summary
  ├─ final_geometries.trj         : concatenated MEP (energy on line 2 in each block)
  ├─ final_geometries.pdb         : produced when input is PDB
  └─ segments/
      ├─ seg_000_gsm/ ...         : initial GSM
      ├─ seg_000_left_opt/ ...    : single-structure optimization of HEI-1
      ├─ seg_000_right_opt/ ...   : single-structure optimization of HEI+1
      ├─ seg_000_refine_gsm/ ...  : refined GSM between End1–End2
      ├─ seg_001_...              : left-side substeps produced by recursion
      └─ seg_002_...              : right-side substeps produced by recursion
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Callable  # Using Callable for type hints

import sys
import traceback
import textwrap
import tempfile
import os
import time  # <<< elapsed time

import click
import numpy as np
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import AU2KCALPERMOL  # <<< for kcal conversion

from .uma_pysis import uma_pysis
from .utils import convert_xyz_to_pdb, freeze_links
from .trj2fig import run_trj2fig  # <<< auto-generate figure when .trj is produced

from .bond_changes import compare_structures, summarize_changes

# -----------------------------------------------
# Defaults
# -----------------------------------------------

# Geometry (input handling)
GEOM_KW: Dict[str, Any] = {
    "coord_type": "cart",   # GrowingString recommends Cartesian
    "freeze_atoms": [],     # 0-based indices
}

# UMA calculator settings
CALC_KW: Dict[str, Any] = {
    "charge": 0,
    "spin": 1,                  # multiplicity (=2S+1)
    "model": "uma-s-1p1",
    "task_name": "omol",
    "device": "auto",
    "max_neigh": None,
    "radius": None,
    "r_edges": False,
}

# GrowingString (path representation)
GS_KW: Dict[str, Any] = {
    "max_nodes": 30,            # including endpoints, string has max_nodes+2 images
    "perp_thresh": 5e-3,
    "reparam_check": "rms",     # "rms" | "norm"
    "reparam_every": 1,
    "reparam_every_full": 1,
    "param": "equi",            # "equi" | "energy"
    "max_micro_cycles": 10,
    "reset_dlc": True,
    "climb": True,              # allow True for the first segment
    "climb_rms": 5e-4,
    "climb_lanczos": True,
    "climb_lanczos_rms": 5e-4,
    "climb_fixed": False,
    "scheduler": None,          # serial computation (shared calculator)
}

# StringOptimizer (GSM optimization control)
OPT_KW: Dict[str, Any] = {
    "type": "string",
    "stop_in_when_full": 1000,  # tolerance after fully grown
    "align": False,             # do not use StringOptimizer's align
    "scale_step": "global",     # global | per_image
    "max_cycles": 1000,
    "dump": False,
    "dump_restart": False,
    "reparam_thresh": 1e-3,
    "coord_diff_thresh": 0.0,
    "out_dir": "./result_path_search/",
    "print_every": 1,
}

# Single-structure optimization (common settings for LBFGS/RFO)
SOPT_BASE_KW: Dict[str, Any] = {
    "thresh": "gau",            # convergence preset
    "max_cycles": 10000,
    "print_every": 1,
    "min_step_norm": 1e-8,
    "assert_min_step": True,
    "dump": False,
    "dump_restart": False,
    "prefix": "",
    "out_dir": "./result_path_search/",
}

# For LBFGS
LBFGS_KW: Dict[str, Any] = {
    **SOPT_BASE_KW,
    "keep_last": 7,
    "beta": 1.0,
    "gamma_mult": False,
    "max_step": 0.30,
    "control_step": True,
    "double_damp": True,
    "line_search": True,
    "mu_reg": None,
    "max_mu_reg_adaptions": 10,
}

# For RFO
RFO_KW: Dict[str, Any] = {
    **SOPT_BASE_KW,
    "trust_radius": 0.30,
    "trust_update": True,
    "trust_min": 0.01,
    "trust_max": 0.30,
    "max_energy_incr": None,
    "hessian_update": "bfgs",
    "hessian_init": "calc",
    "hessian_recalc": 100,
    "hessian_recalc_adapt": 2.0,
    "small_eigval_thresh": 1e-8,
    "line_search": True,
    "alpha0": 1.0,
    "max_micro_cycles": 25,
    "rfo_overlaps": False,
    "gediis": False,
    "gdiis": True,
    "gdiis_thresh": 2.5e-3,
    "gediis_thresh": 1.0e-2,
    "gdiis_test_direction": True,
    "adapt_step_func": False,
}

# Parameters for detecting covalent-bond changes
BOND_KW: Dict[str, Any] = {
    "device": "cuda",
    "bond_factor": 1.20,
    "margin_fraction": 0.05,
    "delta_fraction": 0.05,
}

# Global search control
SEARCH_KW: Dict[str, Any] = {
    "max_depth": 10,               # maximum recursion steps
    "stitch_rmsd_thresh": 1.0e-4,  # RMSD (Å) threshold for considering endpoints as duplicates during stitching
    "bridge_rmsd_thresh": 1.0e-4,  # if endpoint mismatch exceeds this, run a bridge GSM
    "rmsd_align": True,            # use Kabsch alignment when computing RMSD
    "max_nodes_segment": 20,
    "max_nodes_bridge": 5,
    # --- kink detection thresholds and interpolation count (new) ---
    "kink_disp_thresh": 0.25,      # Å: if max atomic displacement between End1–End2 ≤ this, skip GSM; use interpolation + sopt
    "kink_max_nodes": 3,           # number of linear interpolation internal nodes to create
}

# -----------------------------------------------
# Utilities
# -----------------------------------------------

def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively overwrite dict *dst* with *src*, and return *dst*."""
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    """Load YAML and return a dict (empty dict if *path* is None)."""
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")
    return data


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    """Produce a readable string dump for settings."""
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Pretty-format freeze_atoms for display."""
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple, np.ndarray)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if len(fa) else ""
    return g


def _freeze_links_for_pdb(pdb_path: Path) -> Sequence[int]:
    """Detect parent atoms of link hydrogens in PDB and return 0-based indices."""
    try:
        return freeze_links(pdb_path)
    except Exception as e:
        click.echo(f"[freeze-links] WARNING: Could not detect link parents for '{pdb_path.name}': {e}", err=True)
        return []


def _load_two_endpoints(
    paths: Sequence[Path],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> Sequence:
    """Load two or more input structures and set freeze_atoms when needed."""
    geoms = []
    for p in paths:
        g = geom_loader(p, coord_type=coord_type)
        freeze = list(base_freeze)
        if auto_freeze_links and p.suffix.lower() == ".pdb":
            detected = _freeze_links_for_pdb(p)
            if detected:
                freeze = sorted(set(freeze).union(detected))
                click.echo(f"[freeze-links] {p.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# --- For multiple-structure inputs (new) ---
def _load_structures(
    paths: Sequence[Path],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> List[Any]:
    """Load multiple structures and set freeze_atoms as needed; return a list of geometries."""
    geoms: List[Any] = []
    for p in paths:
        g = geom_loader(p, coord_type=coord_type)
        freeze = list(base_freeze)
        if auto_freeze_links and p.suffix.lower() == ".pdb":
            detected = _freeze_links_for_pdb(p)
            if detected:
                freeze = sorted(set(freeze).union(detected))
                click.echo(f"[freeze-links] {p.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms
# ---------------------------------


def _ensure_calc_on_geom(g, calc) -> None:
    """Attach a pysisyphus Calculator to Geometry (overwrite if necessary)."""
    try:
        g.set_calculator(calc)
    except Exception:
        g.set_calculator(calc)


def _write_xyz_trj_with_energy(images: Sequence, energies: Sequence[float], path: Path) -> None:
    """Write an XYZ .trj sequence with energy on the second line of each block."""
    blocks: List[str] = []
    E = np.array(energies, dtype=float)
    for geom, e in zip(images, E):
        s = geom.as_xyz()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{e:.12f}"
        s_mod = "\n".join(lines)
        if not s_mod.endswith("\n"):
            s_mod += "\n"
        blocks.append(s_mod)
    with open(path, "w") as f:
        f.write("".join(blocks))


def _maybe_convert_to_pdb(in_path: Path, ref_pdb_path: Optional[Path], out_path: Optional[Path] = None) -> Optional[Path]:
    """
    If inputs are PDB, convert the specified .xyz / .trj to PDB and save.
    Returns the output path on success; otherwise returns None.
    """
    try:
        if ref_pdb_path is None:
            return None
        if not in_path.exists():
            return None
        if in_path.suffix.lower() not in (".xyz", ".trj"):
            return None
        out_pdb = out_path if out_path is not None else in_path.with_suffix(".pdb")
        convert_xyz_to_pdb(in_path, ref_pdb_path, out_pdb)
        click.echo(f"[convert] Wrote '{out_pdb}'.")
        return out_pdb
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert '{in_path.name}' to PDB: {e}", err=True)
        return None


def _kabsch_rmsd(A: np.ndarray, B: np.ndarray, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """RMSD after Kabsch alignment (or simple RMSD if align=False). Supports subset selection via indices."""
    assert A.shape == B.shape and A.shape[1] == 3
    if indices is not None and len(indices) > 0:
        idx = np.array(sorted({int(i) for i in indices if 0 <= int(i) < A.shape[0]}), dtype=int)
        if idx.size == 0:
            idx = np.arange(A.shape[0], dtype=int)
        A = A[idx]
        B = B[idx]
    if not align:
        diff = A - B
        return float(np.sqrt((diff * diff).sum() / A.shape[0]))

    Ac = A - A.mean(axis=0, keepdims=True)
    Bc = B - B.mean(axis=0, keepdims=True)
    H = Ac.T @ Bc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # Handle improper rotation
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    Arot = Ac @ R
    diff = Arot - Bc
    return float(np.sqrt((diff * diff).sum() / A.shape[0]))


def _rmsd_between(ga, gb, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """RMSD between two Geometries (optional Kabsch alignment and subset)."""
    return _kabsch_rmsd(np.array(ga.coords3d), np.array(gb.coords3d), align=align, indices=indices)


def _orient_two_by_closeness(a, b, end1, end2, align: bool = True) -> Tuple:
    """Given endpoints a,b and candidates end1,end2, return (left,right) where left is closer to a and right closer to b."""
    r1a = _rmsd_between(end1, a, align)
    r2a = _rmsd_between(end2, a, align)
    if r1a <= r2a:
        left, right = end1, end2
    else:
        left, right = end2, end1
    return left, right


def _has_bond_change(x, y, bond_cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """Return (bool, summary) indicating if any covalent bonds form/break between x and y."""
    res = compare_structures(
        x, y,
        device=bond_cfg.get("device", "cuda"),
        bond_factor=float(bond_cfg.get("bond_factor", 1.20)),
        margin_fraction=float(bond_cfg.get("margin_fraction", 0.05)),
        delta_fraction=float(bond_cfg.get("delta_fraction", 0.05)),
    )
    formed = len(res.formed_covalent) > 0
    broken = len(res.broken_covalent) > 0
    summary = summarize_changes(x, res, one_based=True)
    return (formed or broken), summary


# ---- Kabsch alignment producing a new cloned geometry (subset-capable) ----

def _kabsch_aligned_clone(ref_geom, mob_geom, coord_type: str, use_indices: Optional[Sequence[int]] = None) -> Any:
    """
    Create a new Geometry by Kabsch-aligning *mob_geom* to *ref_geom* (mob is transformed).
    If *use_indices* is given, the optimal rigid transform is determined from the subset, but applied to all atoms.
    The original *mob_geom* is not modified.
    """
    A = np.asarray(ref_geom.coords3d, dtype=float)  # (N,3)
    B = np.asarray(mob_geom.coords3d, dtype=float)  # (N,3)
    assert A.shape == B.shape and A.shape[1] == 3, "Kabsch alignment requires same atom ordering."

    N = A.shape[0]
    if use_indices is not None and len(use_indices) > 0:
        idx = np.array(sorted({int(i) for i in use_indices if 0 <= int(i) < N}), dtype=int)
        if idx.size == 0:
            idx = np.arange(N, dtype=int)
    else:
        idx = np.arange(N, dtype=int)

    CA_sel = A[idx].mean(axis=0, keepdims=True)
    CB_sel = B[idx].mean(axis=0, keepdims=True)
    Ac_sel = A[idx] - CA_sel
    Bc_sel = B[idx] - CB_sel

    H = Ac_sel.T @ Bc_sel
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    B_aligned = (B - CB_sel) @ R + CA_sel  # (N,3)

    atoms = [a for a in mob_geom.atoms]
    lines = [str(len(atoms)), ""]
    for sym, (x, y, z) in zip(atoms, B_aligned):
        lines.append(f"{sym} {x:.15f} {y:.15f} {z:.15f}")
    xyz_str = "\n".join(lines) + "\n"

    tmp = None
    try:
        tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
        tmp.write(xyz_str)
        tmp.flush()
        tmp.close()
        g_new = geom_loader(Path(tmp.name), coord_type=coord_type)
        # propagate freeze atoms
        try:
            g_new.freeze_atoms = np.array(getattr(mob_geom, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
    finally:
        if tmp is not None:
            try:
                os.unlink(tmp.name)
            except Exception:
                pass

    return g_new


# ---------- Robust Kabsch R,t computation (minimal version) ----------

def _kabsch_R_t(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute optimal rotation R and translation t (Kabsch) that aligns Q to P in least squares sense.

    Parameters
    ----------
    P, Q : (N, 3) arrays

    Returns
    -------
    R : (3, 3) rotation matrix
    t : (3,) translation vector
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("Kabsch expects P, Q with shape (N, 3).")
    mu_P = P.mean(axis=0)
    mu_Q = Q.mean(axis=0)
    Pc = P - mu_P
    Qc = Q - mu_Q
    # NOTE: covariance order is Pc.T @ Qc (maps Q -> P). Do not swap!
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # ensure right-handed
    if np.linalg.det(R) < 0.0:
        Vt[-1, :] *= -1.0
        R = Vt.T @ U.T
    t = mu_P - mu_Q @ R
    return R, t


def _align_second_to_first_kabsch(geom_ref, geom_to_align) -> Tuple[float, float, int]:
    """
    Rigidly align *geom_to_align* to *geom_ref* with Kabsch (no scaling).
    If either endpoint has freeze_atoms, the transform is determined on that subset only,
    and applied to all atoms.

    Returns
    -------
    rmsd_before, rmsd_after, n_used  (RMSDs are evaluated on the selection)
    """
    P = np.array(geom_ref.coords3d, dtype=float)    # (N, 3)
    Q = np.array(geom_to_align.coords3d, dtype=float)
    if P.shape != Q.shape:
        raise ValueError(f"Different atom counts for endpoints: {P.shape[0]} vs {Q.shape[0]}")

    N = P.shape[0]
    # choose indices used to decide transform
    fa0 = getattr(geom_ref, "freeze_atoms", np.array([], dtype=int))
    fa1 = getattr(geom_to_align, "freeze_atoms", np.array([], dtype=int))
    freeze_union = sorted(set(map(int, fa0)) | set(map(int, fa1)))

    if len(freeze_union) > 0:
        use_mask = np.zeros(N, dtype=bool)
        # defensively ignore out-of-range indices
        valid_idx = [i for i in freeze_union if 0 <= i < N]
        use_mask[valid_idx] = True
    else:
        use_mask = np.ones(N, dtype=bool)

    P_sel = P[use_mask]
    Q_sel = Q[use_mask]
    n_used = int(P_sel.shape[0])

    # RMSD before alignment (on the selection; no extra alignment)
    def _rmsd(A: np.ndarray, B: np.ndarray) -> float:
        return float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1)))) if len(A) else float("nan")

    rmsd_before = _rmsd(P_sel, Q_sel)

    R, t = _kabsch_R_t(P_sel, Q_sel)
    Q_aligned = (Q @ R) + t  # apply to all atoms

    # set_coords expects 1D array (length 3N)
    geom_to_align.set_coords(Q_aligned.reshape(-1))

    rmsd_after = _rmsd(P[use_mask], Q_aligned[use_mask])

    return rmsd_before, rmsd_after, n_used
# ---------- end of additional robust alignment ----------


# ---------- GS configuration utility (minimal) ----------
def _gs_cfg_with_overrides(base: Dict[str, Any], **overrides: Any) -> Dict[str, Any]:
    """Return a shallow copy of GrowingString config with specified overrides."""
    cfg = dict(base)
    for k, v in overrides.items():
        cfg[k] = v
    return cfg
# -----------------------------------------------

# ---------- Kink detection & interpolation (new) ----------

def _max_displacement_between(ga, gb, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """
    Return maximum per-atom displacement (Å) between two structures.
    If *align* is True, Kabsch-align gb to ga first (transform decided by *indices* subset if provided;
    transform is applied to all atoms).
    """
    A = np.asarray(ga.coords3d, dtype=float)  # (N,3)
    B = np.asarray(gb.coords3d, dtype=float)
    if A.shape != B.shape or A.shape[1] != 3:
        raise ValueError("Geometries must have same number of atoms for displacement.")
    if align:
        if indices is None or len(indices) == 0:
            idx = np.arange(A.shape[0], dtype=int)
        else:
            idx = np.array(sorted({int(i) for i in indices if 0 <= int(i) < A.shape[0]}), dtype=int)
            if idx.size == 0:
                idx = np.arange(A.shape[0], dtype=int)
        R, t = _kabsch_R_t(A[idx], B[idx])
        B = (B @ R) + t
    disp = np.linalg.norm(A - B, axis=1)  # per-atom distances
    return float(np.max(disp))


def _new_geom_from_coords(atoms: Sequence[str], coords: np.ndarray, coord_type: str, freeze_atoms: Sequence[int]) -> Any:
    """Create a Geometry from coordinates via an XYZ string, attach freeze_atoms, and return it."""
    lines = [str(len(atoms)), ""]
    for sym, (x, y, z) in zip(atoms, coords):
        lines.append(f"{sym} {x:.15f} {y:.15f} {z:.15f}")
    s = "\n".join(lines) + "\n"
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()
        g = geom_loader(Path(tmp.name), coord_type=coord_type)
        g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _make_linear_interpolations(gL, gR, n_internal: int) -> List[Any]:
    """
    Return linear interpolation structures between gL→gR (n_internal internal points).
    Endpoints are not included. Atom order follows gL.
    """
    A = np.asarray(gL.coords3d, dtype=float)
    B = np.asarray(gR.coords3d, dtype=float)
    assert A.shape == B.shape and A.shape[1] == 3, "Atom counts must match for interpolation."
    atoms = [a for a in gL.atoms]
    coord_type = gL.coord_type
    # union of freeze_atoms
    faL = getattr(gL, "freeze_atoms", np.array([], dtype=int))
    faR = getattr(gR, "freeze_atoms", np.array([], dtype=int))
    freeze_union = sorted(set(map(int, faL)) | set(map(int, faR)))
    interps: List[Any] = []
    for k in range(1, n_internal + 1):
        t = k / (n_internal + 1.0)
        C = (1.0 - t) * A + t * B
        interps.append(_new_geom_from_coords(atoms, C, coord_type, freeze_union))
    return interps


def _energy_of(g) -> float:
    """Return the energy (Hartree) of a Geometry with an attached Calculator."""
    _ensure_calc_on_geom(g, getattr(g, "calculator", None))
    return float(g.energy)


@dataclass
class GSMResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int


# ---- Per-segment summary (additional) ----
@dataclass
class SegmentReport:
    tag: str
    barrier_kcal: float
    delta_kcal: float
    summary: str  # summarize_changes string (empty if no changes)


def _run_gsm_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # reference PDB for conversion
) -> GSMResult:
    """Run GSM between gA–gB and save segment outputs."""
    # Attach calculator to endpoints
    for g in (gA, gB):
        _ensure_calc_on_geom(g, shared_calc)

    def calc_getter():
        return shared_calc

    # Setup GrowingString
    gs = GrowingString(
        images=[gA, gB],
        calc_getter=calc_getter,
        **gs_cfg,
    )

    # Prepare optimizer
    _opt_args = dict(opt_cfg)
    seg_dir = out_dir / f"{tag}_gsm"
    seg_dir.mkdir(parents=True, exist_ok=True)
    _opt_args["out_dir"] = str(seg_dir)

    optimizer = StringOptimizer(
        geometry=gs,
        **{k: v for k, v in _opt_args.items() if k != "type"}
    )

    click.echo(f"\n=== [{tag}] GSM started ===\n")
    optimizer.run()
    click.echo(f"\n=== [{tag}] GSM finished ===\n")

    # energies & images
    energies = list(map(float, np.array(gs.energy, dtype=float)))
    images = list(gs.images)
    # write trajectory
    final_trj = seg_dir / "final_geometries.trj"
    wrote_with_energy = True
    try:
        _write_xyz_trj_with_energy(images, energies, final_trj)
        click.echo(f"[{tag}] Wrote '{final_trj}'.")
    except Exception:
        wrote_with_energy = False
        with open(final_trj, "w") as f:
            f.write(gs.as_xyz())
        click.echo(f"[{tag}] Wrote '{final_trj}'.")

    # >>> plot .trj if possible
    try:
        if wrote_with_energy:
            run_trj2fig(final_trj, [seg_dir / "energy.png"], unit="kcal", reference="init", reverse_x=False)
            click.echo(f"[{tag}] Plot  → '{seg_dir / 'energy.png'}'")
        else:
            click.echo(f"[{tag}] WARNING: energy comments missing; skip plotting.", err=True)
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    # If input is PDB, convert intermediate .trj to PDB
    _maybe_convert_to_pdb(final_trj, ref_pdb_path, seg_dir / "final_geometries.pdb")

    # write HEI
    try:
        E = np.array(energies, dtype=float)
        hei_idx = int(np.argmax(E))
        hei_geom = images[hei_idx]
        hei_E = float(E[hei_idx])

        hei_xyz = seg_dir / "gsm_hei.xyz"
        s = hei_geom.as_xyz()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{hei_E:.12f}"
            s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
        with open(hei_xyz, "w") as f:
            f.write(s)
        click.echo(f"[{tag}] Wrote '{hei_xyz}'.")
        # also convert HEI to PDB
        _maybe_convert_to_pdb(hei_xyz, ref_pdb_path, seg_dir / "gsm_hei.pdb")
    except Exception:
        hei_idx = int(np.argmax(np.array(energies, dtype=float)))

    return GSMResult(images=images, energies=energies, hei_idx=hei_idx)


def _optimize_single(
    g,
    shared_calc,
    sopt_kind: str,
    sopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
):
    """Run single-structure optimization (LBFGS/RFO) and return the final Geometry."""
    # Attach calculator
    _ensure_calc_on_geom(g, shared_calc)

    seg_dir = out_dir / f"{tag}_{sopt_kind}_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(sopt_cfg)
    args["out_dir"] = str(seg_dir)

    if sopt_kind == "lbfgs":
        opt = LBFGS(g, **args)
    else:
        opt = RFOptimizer(g, **args)

    click.echo(f"\n=== [{tag}] single-structure {sopt_kind.upper()} started ===\n")
    opt.run()
    click.echo(f"\n=== [{tag}] single-structure {sopt_kind.upper()} finished ===\n")

    # Optimizers write the final structure to file (final_fn)
    try:
        final_xyz = Path(opt.final_fn) if isinstance(opt.final_fn, (str, Path)) else Path(opt.final_fn)
        # if inputs are PDB, also convert final .xyz/.trj to PDB
        _maybe_convert_to_pdb(final_xyz, ref_pdb_path)
        g_final = geom_loader(final_xyz, coord_type=g.coord_type)
        # propagate freeze atoms
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        _ensure_calc_on_geom(g_final, shared_calc)
        return g_final
    except Exception:
        # fallback: return the current geometry
        return g


def _refine_between(
    gL,
    gR,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
) -> GSMResult:
    """Refine between End1–End2 via GSM (force climb=True)."""
    gs_refine_cfg = _gs_cfg_with_overrides(gs_cfg, climb=True, climb_lanczos=True)
    return _run_gsm_between(gL, gR, shared_calc, gs_refine_cfg, opt_cfg, out_dir, tag=f"{tag}_refine", ref_pdb_path=ref_pdb_path)


def _maybe_bridge_segments(
    tail_g,
    head_g,
    shared_calc,
    gs_cfg: Dict[str, Any],  # bridge-specific GS config
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    rmsd_thresh: float,
    ref_pdb_path: Optional[Path],  # for PDB conversion
) -> Optional[GSMResult]:
    """Run a bridge GSM if two segment endpoints are farther than the threshold."""
    rmsd = _rmsd_between(tail_g, head_g, align=True)
    if rmsd <= rmsd_thresh:
        return None
    click.echo(f"[{tag}] Gap detected between segments (RMSD={rmsd:.4e} Å) — bridging by GSM.")
    return _run_gsm_between(tail_g, head_g, shared_calc, gs_cfg, opt_cfg, out_dir, tag=f"{tag}_bridge", ref_pdb_path=ref_pdb_path)


def _stitch_paths(
    parts: List[Tuple[List[Any], List[float]]],
    stitch_rmsd_thresh: float,
    bridge_rmsd_thresh: float,
    shared_calc,
    gs_cfg,   # GrowingString config for bridges (climb=False, max_nodes=search.max_nodes_bridge)
    opt_cfg,
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
    bond_cfg: Optional[Dict[str, Any]] = None,  # to detect bond changes between adjacent segments
    segment_builder: Optional[Callable[[Any, Any, str], "CombinedPath"]] = None,  # <<< returns CombinedPath
    segments_out: Optional[List["SegmentReport"]] = None,  # <<< append inserted segment summaries in order
) -> Tuple[List[Any], List[float]]:
    """
    Concatenate path parts (images, energies). Insert bridge GSMs when needed.
    If covalent changes are detected between the tail of the previous part and the head of the next part,
    generate and insert a *new* segment using `segment_builder` (recursive GSM), instead of bridging.
    """
    all_imgs: List[Any] = []
    all_E: List[float] = []

    def append_part(imgs: List[Any], Es: List[float]) -> None:
        nonlocal all_imgs, all_E
        if not imgs:
            return
        if not all_imgs:
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return
        # previous tail and current head
        tail = all_imgs[-1]
        head = imgs[0]

        # check covalent changes between adjacent endpoints
        adj_changed, adj_summary = False, ""
        if segment_builder is not None and bond_cfg is not None:
            try:
                adj_changed, adj_summary = _has_bond_change(tail, head, bond_cfg)
            except Exception:
                adj_changed, adj_summary = False, ""

        if adj_changed and segment_builder is not None:
            click.echo(f"[{tag}] Bond-change detected between segment endpoints — inserting a new segment via recursive GSM.")
            if adj_summary:
                click.echo(textwrap.indent(adj_summary, prefix="  "))
            # build and insert a new segment (tail→head)
            sub = segment_builder(tail, head, f"{tag}_mid")
            seg_imgs, seg_E = sub.images, sub.energies
            if segments_out is not None and getattr(sub, "segments", None):
                segments_out.extend(sub.segments)  # <<< append in the appearance order
            if seg_imgs:
                # remove duplicate with existing tail
                if _rmsd_between(all_imgs[-1], seg_imgs[0], align=True) <= stitch_rmsd_thresh:
                    seg_imgs = seg_imgs[1:]
                    seg_E = seg_E[1:]
                all_imgs.extend(seg_imgs)
                all_E.extend(seg_E)
            # then connect the current part (deduplicate again)
            if _rmsd_between(all_imgs[-1], imgs[0], align=True) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return  # done for this part

        rmsd = _rmsd_between(tail, head, align=True)
        if rmsd <= stitch_rmsd_thresh:
            # perfect duplicate → drop head and concatenate
            all_imgs.extend(imgs[1:])
            all_E.extend(Es[1:])
        elif rmsd > bridge_rmsd_thresh:
            # larger mismatch → run bridge GSM
            br = _maybe_bridge_segments(
                tail, head, shared_calc, gs_cfg, opt_cfg, out_dir, tag=tag,
                rmsd_thresh=bridge_rmsd_thresh, ref_pdb_path=ref_pdb_path
            )
            if br is not None:
                # drop duplicate at bridge start if any
                b_imgs, b_E = br.images, br.energies
                if _rmsd_between(all_imgs[-1], b_imgs[0], align=True) <= stitch_rmsd_thresh:
                    b_imgs = b_imgs[1:]
                    b_E = b_E[1:]
                if b_imgs:
                    all_imgs.extend(b_imgs)
                    all_E.extend(b_E)
            # then connect the current part (re-check duplicates)
            if _rmsd_between(all_imgs[-1], imgs[0], align=True) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
        else:
            # minor mismatch → connect directly
            all_imgs.extend(imgs)
            all_E.extend(Es)

    for (imgs, Es) in parts:
        append_part(imgs, Es)

    return all_imgs, all_E


# -----------------------------------------------
# Recursive search (core)
# -----------------------------------------------

@dataclass
class CombinedPath:
    images: List[Any]
    energies: List[float]
    segments: List[SegmentReport]  # <<< segment summaries in output order


def _build_multistep_path(
    gA,
    gB,
    shared_calc,
    geom_cfg: Dict[str, Any],
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    sopt_kind: str,
    sopt_cfg: Dict[str, Any],
    bond_cfg: Dict[str, Any],
    search_cfg: Dict[str, Any],
    out_dir: Path,
    ref_pdb_path: Optional[Path],
    depth: int,
    seg_counter: List[int],
    branch_tag: str,
) -> CombinedPath:
    """
    Recursively construct a multistep MEP from A–B and return it (A→B order).
    """
    # GS config for *segment* (apply search.max_nodes_segment)
    seg_max_nodes = int(search_cfg.get("max_nodes_segment", gs_cfg.get("max_nodes", 30)))
    gs_seg_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=seg_max_nodes)

    if depth > int(search_cfg.get("max_depth", 10)):
        click.echo(f"[{branch_tag}] Reached max recursion depth. Returning current endpoints only.")
        # at upper bound: simply run GSM for A–B (climb=True)
        gsm = _run_gsm_between(gA, gB, shared_calc, gs_seg_cfg, opt_cfg, out_dir, tag=f"seg_{seg_counter[0]:03d}_maxdepth", ref_pdb_path=ref_pdb_path)
        seg_counter[0] += 1
        return CombinedPath(images=gsm.images, energies=gsm.energies, segments=[])

    seg_id = seg_counter[0]
    seg_counter[0] += 1
    tag0 = f"seg_{seg_id:03d}"

    # 1) Initial GSM A–B (segment settings)
    gs_seg_cfg_first = _gs_cfg_with_overrides(gs_seg_cfg, climb=True, climb_lanczos=True)
    gsm0 = _run_gsm_between(gA, gB, shared_calc, gs_seg_cfg_first, opt_cfg, out_dir, tag=tag0, ref_pdb_path=ref_pdb_path)

    # HEI and its neighbors
    hei = int(gsm0.hei_idx)
    if not (1 <= hei <= len(gsm0.images) - 2):
        # edge case: HEI at an endpoint
        click.echo(f"[{tag0}] WARNING: HEI is at an endpoint (idx={hei}). Returning raw GSM.")
        return CombinedPath(images=gsm0.images, energies=gsm0.energies, segments=[])

    left_img = gsm0.images[hei - 1]
    right_img = gsm0.images[hei + 1]

    # 2) Single-structure optimizations of both neighbors
    left_opt = _optimize_single(left_img, shared_calc, sopt_kind, sopt_cfg, out_dir, tag=f"{tag0}_left", ref_pdb_path=ref_pdb_path)
    right_opt = _optimize_single(right_img, shared_calc, sopt_kind, sopt_cfg, out_dir, tag=f"{tag0}_right", ref_pdb_path=ref_pdb_path)

    # 3) Orient ends so that left_end is closer to A and right_end closer to B
    left_end, right_end = _orient_two_by_closeness(gA, gB, left_opt, right_opt, align=bool(search_cfg.get("rmsd_align", True)))

    # --- 3.5) Kink detection (new) ---
    # If the max per-atom displacement between End1 and End2 is ≤ threshold, skip GSM and do interpolation + sopt
    try:
        fa_union = sorted(set(map(int, getattr(left_end, "freeze_atoms", []))) |
                          set(map(int, getattr(right_end, "freeze_atoms", []))))
        max_disp = _max_displacement_between(
            left_end, right_end,
            align=False,
            indices=fa_union if len(fa_union) > 0 else None
        )
    except Exception as e:
        click.echo(f"[{tag0}] WARNING: Failed to evaluate max displacement for kink detection: {e}", err=True)
        max_disp = float("inf")

    kink_thresh = float(search_cfg.get("kink_disp_thresh", 0.2))
    use_kink = (max_disp <= kink_thresh)

    # 4) Refined GSM between End1–End2 (or replace with kink interpolation path)
    if use_kink:
        n_inter = int(search_cfg.get("kink_max_nodes", 3))
        click.echo(f"[{tag0}] Kink detected (max disp = {max_disp:.4f} Å ≤ {kink_thresh:.4f} Å). "
                   f"Using {n_inter} linear interpolation nodes + single-structure optimization instead of GSM.")
        # generate linear interpolation structures
        inter_geoms = _make_linear_interpolations(left_end, right_end, n_inter)
        # sopt each interpolation
        opt_inters: List[Any] = []
        for i, g_int in enumerate(inter_geoms, 1):
            g_int.set_calculator(shared_calc)
            g_opt = _optimize_single(g_int, shared_calc, sopt_kind, sopt_cfg, out_dir, tag=f"{tag0}_kink_int{i}", ref_pdb_path=ref_pdb_path)
            opt_inters.append(g_opt)
        # images for the step: left_end, [opt inters], right_end
        step_imgs = [left_end] + opt_inters + [right_end]
        # energies
        step_E = [float(img.energy) for img in step_imgs]
        # wrap in a GSMResult-like form
        ref1 = GSMResult(images=step_imgs, energies=step_E, hei_idx=int(np.argmax(step_E)))
        step_tag_for_report = f"{tag0}_kink"
    else:
        ref1 = _refine_between(left_end, right_end, shared_calc, gs_seg_cfg, opt_cfg, out_dir, tag=tag0, ref_pdb_path=ref_pdb_path)
        step_tag_for_report = f"{tag0}_refine"

    step_imgs, step_E = ref1.images, ref1.energies

    # 5) Check for covalent changes (left: A–left_end, right: right_end–B)
    left_changed, left_summary = _has_bond_change(gA, left_end, bond_cfg)
    right_changed, right_summary = _has_bond_change(right_end, gB, bond_cfg)

    click.echo(f"[{tag0}] Bond-change(A vs left_end): {'Yes' if left_changed else 'No'}")
    if left_changed:
        click.echo(textwrap.indent(left_summary, prefix="  "))
    click.echo(f"[{tag0}] Bond-change(right_end vs B): {'Yes' if right_changed else 'No'}")
    if right_changed:
        click.echo(textwrap.indent(right_summary, prefix="  "))

    # Segment summary (this step)
    try:
        barrier_kcal = (max(step_E) - step_E[0]) * AU2KCALPERMOL
        delta_kcal = (step_E[-1] - step_E[0]) * AU2KCALPERMOL
    except Exception:
        barrier_kcal = float("nan")
        delta_kcal = float("nan")
    # covalent changes within this step (left_end→right_end)
    _changed, step_summary = _has_bond_change(step_imgs[0], step_imgs[-1], bond_cfg)
    seg_report = SegmentReport(
        tag=step_tag_for_report,
        barrier_kcal=float(barrier_kcal),
        delta_kcal=float(delta_kcal),
        summary=step_summary if _changed else "(no covalent-bond changes detected)"
    )

    # 6) Recurse (left/right)
    parts: List[Tuple[List[Any], List[float]]] = []
    seg_reports: List[SegmentReport] = []

    if left_changed:
        subL = _build_multistep_path(
            gA, left_end, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg, bond_cfg, search_cfg,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}L"
        )
        parts.append((subL.images, subL.energies))
        seg_reports.extend(subL.segments)

    # central step (this one)
    parts.append((step_imgs, step_E))
    seg_reports.append(seg_report)

    if right_changed:
        subR = _build_multistep_path(
            right_end, gB, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg, bond_cfg, search_cfg,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}R"
        )
        parts.append((subR.images, subR.energies))
        seg_reports.extend(subR.segments)

    # 7) Stitch (remove duplicates; optionally insert bridge/new segments)
    #    GS config for bridges (max_nodes_bridge, climb=False)
    bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 10))
    gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

    # Builder to create a new recursive segment when tail–head shows covalent changes
    def _segment_builder(tail_g, head_g, _tag: str) -> CombinedPath:
        sub = _build_multistep_path(
            tail_g, head_g,
            shared_calc,
            geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg,
            bond_cfg, search_cfg,
            out_dir=out_dir,
            ref_pdb_path=ref_pdb_path,
            depth=depth + 1,
            seg_counter=seg_counter,
            branch_tag=f"{branch_tag}B",
        )
        return sub

    stitched_imgs, stitched_E = _stitch_paths(
        parts,
        stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-3)),
        bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-2)),
        shared_calc=shared_calc,
        gs_cfg=gs_bridge_cfg,   # bridges: climb=False and dedicated max_nodes
        opt_cfg=opt_cfg,
        out_dir=out_dir,
        tag=tag0,
        ref_pdb_path=ref_pdb_path,
        bond_cfg=bond_cfg,
        segment_builder=_segment_builder,
        segments_out=seg_reports,   # <<< collect summaries of any inserted segments at the correct position
    )

    return CombinedPath(images=stitched_imgs, energies=stitched_E, segments=seg_reports)


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Multistep MEP search by recursive GSM segmentation.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,   # <<< FIX: In Click, options can't have nargs=-1; use multiple=True with repeated -i.
    required=True,
    help="Two or more structures in reaction order (repeat -i: e.g., -i reactant.pdb -i im1.pdb -i product.pdb).",
)
@click.option("-q", "--charge", type=int, required=True, help="Total charge")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1)")
@click.option("--freeze-links", "freeze_links_flag", type=click.BOOL, default=True, show_default=True,
              help="If PDB, freeze parent atoms of link hydrogens.")
@click.option("--max-nodes", type=int, default=30, show_default=True,
              help="Internal nodes (string has max_nodes+2 images incl. endpoints). "
                   "This value is used for *segment* GSM unless overridden by YAML search.max_nodes_segment.")
@click.option("--max-cycles", type=int, default=1000, show_default=True, help="Max GSM optimization cycles")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Search for transition state after path growth")
@click.option("--sopt-mode", type=click.Choice(["lbfgs", "rfo", "light", "heavy"], case_sensitive=False),
              default="lbfgs", show_default=True, help="Single-structure optimizer: lbfgs(light) or rfo(heavy)")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump GSM/single-optimization trajectories during run")
@click.option("--out-dir", "out_dir", type=str, default="./result_path_search/", show_default=True, help="Output directory")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt, sopt, bond, search).",
)
@click.option(
    "--pre-opt", "--pre_opt",
    "pre_opt",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="If True, SKIP initial single-structure optimizations of inputs (i.e., do not pre-opt). "
         "Default False keeps the original behavior (perform pre-optimization).",
)
def cli(
    input_paths: Sequence[Path],
    charge: int,
    spin: int,
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    sopt_mode: str,
    dump: bool,
    out_dir: str,
    args_yaml: Optional[Path],
    pre_opt: bool,
) -> None:
    time_start = time.perf_counter()  # <<< start timing
    try:
        # --------------------------
        # 0) Input validation (multi-structure)
        # --------------------------
        if len(input_paths) < 2:
            raise click.BadParameter("You must provide at least two structures for --input in reaction order (reactant [intermediates ...] product).")

        # --------------------------
        # 1) Resolve settings (defaults ← YAML ← CLI)
        # --------------------------
        yaml_cfg = _load_yaml(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg   = dict(GS_KW)
        opt_cfg  = dict(OPT_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)
        bond_cfg  = dict(BOND_KW)
        search_cfg = dict(SEARCH_KW)

        _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
        _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
        _deep_update(gs_cfg,   yaml_cfg.get("gs",   {}))
        _deep_update(opt_cfg,  yaml_cfg.get("opt",  {}))
        _deep_update(lbfgs_cfg, yaml_cfg.get("sopt", {}).get("lbfgs", yaml_cfg.get("lbfgs", {})))
        _deep_update(rfo_cfg,   yaml_cfg.get("sopt", {}).get("rfo",   yaml_cfg.get("rfo",   {})))
        _deep_update(bond_cfg,  yaml_cfg.get("bond", {}))
        _deep_update(search_cfg,yaml_cfg.get("search", {}))

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)

        gs_cfg["max_nodes"] = int(max_nodes)           # base (mainly for display/back-compat)
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)

        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir

        # Keep dump/out_dir consistent for single-structure optimizers
        lbfgs_cfg["dump"] = bool(dump)
        rfo_cfg["dump"]   = bool(dump)
        lbfgs_cfg["out_dir"] = out_dir
        rfo_cfg["out_dir"]   = out_dir

        # normalize sopt kind
        sopt_kind = sopt_mode.strip().lower()
        if sopt_kind in ("light", "lbfgs"):
            sopt_kind = "lbfgs"
            sopt_cfg = lbfgs_cfg
        elif sopt_kind in ("heavy", "rfo"):
            sopt_kind = "rfo"
            sopt_cfg = rfo_cfg
        else:
            raise click.BadParameter(f"Unknown --sopt-mode '{sopt_mode}'.")

        # finalize max_nodes for segment/bridge
        # - If YAML(search.max_nodes_segment) is absent, use CLI --max-nodes
        if "max_nodes_segment" not in yaml_cfg.get("search", {}):
            search_cfg["max_nodes_segment"] = int(max_nodes)
        # - For bridge, respect YAML/defaults (override via YAML search.max_nodes_bridge if needed)

        # For display: resolved settings
        out_dir_path = Path(out_dir).resolve()
        echo_geom = _format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_opt  = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        # pre_opt display
        click.echo(_pretty_block("geom", echo_geom))
        click.echo(_pretty_block("calc", echo_calc))
        click.echo(_pretty_block("gs",   echo_gs))
        click.echo(_pretty_block("opt",  echo_opt))
        click.echo(_pretty_block("sopt."+sopt_kind, sopt_cfg))
        click.echo(_pretty_block("bond", bond_cfg))
        click.echo(_pretty_block("search", search_cfg))
        click.echo(_pretty_block("run_flags", {"pre_opt_skip": bool(pre_opt)}))

        # --------------------------
        # 2) Prepare inputs
        #    (a) Load structures (multi)
        #    (b) Prepare calculator
        #    (c) Pre-optimize endpoints (optional; all inputs)
        #    (d) Freeze-aware Kabsch pre-alignment per adjacent pair
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        p_list = [Path(p) for p in input_paths]

        geoms = _load_structures(
            paths=p_list,
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )

        # (b) Shared UMA calculator (used across all optimizations and GSM)
        shared_calc = uma_pysis(**calc_cfg)
        for g in geoms:
            _ensure_calc_on_geom(g, shared_calc)

        # pick first PDB path as reference (if any)
        ref_pdb_for_segments: Optional[Path] = None
        for p in p_list:
            if p.suffix.lower() == ".pdb":
                ref_pdb_for_segments = p.resolve()
                break

        # (c) Single-structure optimization before alignment (unless pre_opt=True) — applied to all inputs
        if not pre_opt:
            new_geoms: List[Any] = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                g_opt = _optimize_single(g, shared_calc, sopt_kind, sopt_cfg, out_dir_path, tag=tag, ref_pdb_path=ref_pdb_for_segments)
                new_geoms.append(g_opt)
            geoms = new_geoms
        else:
            click.echo("[init] Skipping endpoint pre-optimization as requested by --pre_opt/--pre-opt True.")

        # (d) Per-pair freeze-aware Kabsch pre-alignment (transform decided on selection → applied to all atoms)
        for i in range(len(geoms) - 1):
            gi, gj = geoms[i], geoms[i + 1]
            N = int(np.asarray(gi.coords3d).shape[0])
            fa_i = getattr(gi, "freeze_atoms", np.array([], dtype=int))
            fa_j = getattr(gj, "freeze_atoms", np.array([], dtype=int))
            use_idx = sorted({int(x) for x in list(fa_i) + list(fa_j) if 0 <= int(x) < N})
            idx_for_msg = f"{len(use_idx)} (freeze_atoms)" if len(use_idx) > 0 else f"{N} (all atoms)"
            try:
                rmsd_before, rmsd_after, _n_used = _align_second_to_first_kabsch(gi, gj)
                click.echo(f"[align {i:02d}->{i+1:02d}] Kabsch pre-alignment (used {idx_for_msg}): "
                           f"RMSD_before={rmsd_before:.6f} Å → RMSD_after={rmsd_after:.6f} Å")
            except Exception as e:
                click.echo(f"[align {i:02d}->{i+1:02d}] WARNING: Kabsch pre-alignment skipped: {e}", err=True)

        # --------------------------
        # 3) Run recursive search for each adjacent pair and incrementally stitch
        # --------------------------
        click.echo("\n=== Multistep MEP search (multi-structure) started ===\n")
        seg_counter = [0]  # list for by-reference increment

        # GS config for bridges (max_nodes_bridge, climb=False)
        bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 10))
        gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

        # containers for incremental stitching
        combined_imgs: List[Any] = []
        combined_Es: List[float] = []
        seg_reports_all: List[SegmentReport] = []

        # builder for tail→head recursive segments during stitching
        def _segment_builder_for_pairs(tail_g, head_g, _tag: str) -> CombinedPath:
            return _build_multistep_path(
                tail_g, head_g,
                shared_calc,
                geom_cfg, gs_cfg, opt_cfg,
                sopt_kind, sopt_cfg,
                bond_cfg, search_cfg,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                depth=0,
                seg_counter=seg_counter,
                branch_tag="B",
            )

        for i in range(len(geoms) - 1):
            gA, gB = geoms[i], geoms[i + 1]
            pair_tag = f"pair_{i:02d}"
            click.echo(f"\n--- Processing pair {i:02d}: image {i} → {i+1} ---")
            pair_path = _build_multistep_path(
                gA, gB,
                shared_calc,
                geom_cfg, gs_cfg, opt_cfg,
                sopt_kind, sopt_cfg,
                bond_cfg, search_cfg,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                depth=0,
                seg_counter=seg_counter,
                branch_tag=pair_tag,
            )

            if i == 0:
                combined_imgs = list(pair_path.images)
                combined_Es = list(pair_path.energies)
                seg_reports_all.extend(pair_path.segments)
            else:
                # connect previous result + current pair as two parts (insert bridge/recursive segments if needed)
                parts = [(combined_imgs, combined_Es), (pair_path.images, pair_path.energies)]
                combined_imgs, combined_Es = _stitch_paths(
                    parts=parts,
                    stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-3)),
                    bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-2)),
                    shared_calc=shared_calc,
                    gs_cfg=gs_bridge_cfg,
                    opt_cfg=opt_cfg,
                    out_dir=out_dir_path,
                    tag=pair_tag,
                    ref_pdb_path=ref_pdb_for_segments,
                    bond_cfg=bond_cfg,
                    segment_builder=_segment_builder_for_pairs,
                    segments_out=seg_reports_all,  # add any new segment summaries immediately after insertion
                )
                # After stitching, append the current pair's own segment summaries (so they follow any inserted ones)
                seg_reports_all.extend(pair_path.segments)

        click.echo("\n=== Multistep MEP search (multi-structure) finished ===\n")

        combined_all = CombinedPath(images=combined_imgs, energies=combined_Es, segments=seg_reports_all)

        # --------------------------
        # 4) Outputs
        # --------------------------
        final_trj = out_dir_path / "final_geometries.trj"
        _write_xyz_trj_with_energy(combined_all.images, combined_all.energies, final_trj)
        click.echo(f"[write] Wrote '{final_trj}'.")

        # >>> final MEP energy figure
        try:
            run_trj2fig(final_trj, [out_dir_path / "energy.png"], unit="kcal", reference="init", reverse_x=False)
            click.echo(f"[plot] Figure → '{out_dir_path / 'energy.png'}'")
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot final energy: {e}", err=True)

        # If a reference PDB exists, also write a PDB for the final MEP
        if ref_pdb_for_segments is not None:
            try:
                final_pdb = out_dir_path / "final_geometries.pdb"
                convert_xyz_to_pdb(final_trj, ref_pdb_for_segments, final_pdb)
                click.echo(f"[convert] Wrote '{final_pdb}'.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final MEP to PDB: {e}", err=True)

        # Run summary YAML
        summary = {
            "n_segments": seg_counter[0],
            "n_images_final": len(combined_all.images),
            "out_dir": str(out_dir_path),
            "settings": {
                "geom": echo_geom,
                "calc": echo_calc,
                "gs": echo_gs,
                "opt": echo_opt,
                "sopt_kind": sopt_kind,
                "sopt": sopt_cfg,
                "bond": bond_cfg,
                "search": search_cfg,
                "pre_opt_skip": bool(pre_opt),
            },
            "n_inputs": len(p_list),
            "inputs": [str(p.resolve()) for p in p_list],
        }
        with open(out_dir_path / "summary.yaml", "w") as f:
            yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
        click.echo(f"[write] Wrote '{out_dir_path / 'summary.yaml'}'.")

        # --------------------------
        # 5) Additional console summary
        # --------------------------
        # Overall (first vs last image) covalent changes
        try:
            overall_changed, overall_summary = _has_bond_change(combined_all.images[0], combined_all.images[-1], bond_cfg)
        except Exception:
            overall_changed, overall_summary = False, ""

        # Segment-by-segment report (in output order)
        click.echo("\n=== MEP Summary ===\n")

        # elapsed time
        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed: {hh:02d}:{mm:02d}:{ss:06.3f}")

        # overall covalent changes
        click.echo("\n[overall] Covalent-bond changes between first and last image:")
        if overall_changed and overall_summary.strip():
            click.echo(textwrap.indent(overall_summary.strip(), prefix="  "))
        else:
            click.echo("  (no covalent-bond changes detected)")

        # per-segment summaries
        if combined_all.segments:
            click.echo("\n[segments] In the order they appear along the final MEP:")
            for i, seg in enumerate(combined_all.segments, 1):
                click.echo(f"  [{i:02d}] {seg.tag}  |  barrier = {seg.barrier_kcal:.2f} kcal/mol,  ΔE = {seg.delta_kcal:.2f} kcal/mol")
                if seg.summary.strip():
                    click.echo(textwrap.indent(seg.summary.strip(), prefix="      "))
        else:
            click.echo("\n[segments] (no segment reports)")

    except ZeroStepLength:
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Path search failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during path search:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
