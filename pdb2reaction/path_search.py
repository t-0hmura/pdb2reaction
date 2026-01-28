# pdb2reaction/path_search.py

"""
path_search — Recursive MEP segmentation (GSM/DMF) to build a continuous multistep MEP
=============================================================================

Usage (CLI)
-----------
    pdb2reaction path-search -i R.pdb [I.pdb ...] P.pdb [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [--multiplicity 2S+1]
                            [--workers N] [--workers-per-node N]
                            [--mep-mode {gsm|dmf}] [--freeze-links BOOL] [--thresh PRESET]
                            [--refine-mode peak|minima]
                            [--max-nodes N] [--max-cycles N] [--climb BOOL]
                            [--opt-mode light|heavy] [--dump BOOL]
                            [--convert-files {True|False}]
                            [--out-dir DIR] [--preopt BOOL]
                            [--align {True|False}] [--ref-full-pdb FILE ...]
                            [--ref-pdb FILE ...]
                            [--args-yaml FILE]

Core inputs (strongly recommended):
    -i/--input
        Two or more structures in reaction order (repeatable or space‑separated after a single -i).
    -q/--charge
        Total system charge for the ML region (required for non-`.gjf` inputs unless ``--ligand-charge``
        derivation succeeds; `.gjf` templates supply defaults when available). When ``-q`` is omitted but
        ``--ligand-charge`` is set, the full complex is treated as an enzyme–substrate system and the total
        charge is inferred using ``extract.py``’s residue-aware logic when inputs are PDBs. Explicit ``-q``
        overrides any derived charge; otherwise the run aborts for non-`.gjf` inputs.

Recommended/common:
    -m/--multiplicity
        Spin multiplicity (2S+1); defaults to a .gjf template value when available,
        otherwise 1 when omitted.
    --opt-mode
        Single-structure optimizer: light (=LBFGS) or heavy (=RFO); default light.
    --workers / --workers-per-node
        UMA predictor parallelism (workers > 1 disables analytical Hessians).
    --mep-mode
        Segment generator: GSM (string) or DMF (direct max flux); default gsm.
    --max-nodes
        Internal nodes per MEP segment (GSM string or DMF path); default 10.
    --max-cycles
        Max optimization cycles; default 300.
    --climb {True|False}
        Enable TS search (GSM climbing) for standard GSM segments; bridge segments always run with climbing disabled.
    --preopt {True|False}
        Preoptimize endpoints; default True.
    --align {True|False}
        Rigidly co‑align all inputs after pre‑opt; default on.
    --args-yaml PATH
        YAML with overrides (sections: geom, calc, gs, opt, sopt, bond, search, dmf).
    --thresh STR
        Convergence preset for GSM and single optimizations
        (gau_loose|gau|gau_tight|gau_vtight|baker|never). When omitted, the effective default is `gau`.
    --ref-full-pdb PATH [...]
        Full template PDB(s) for final merge (see Notes).
    --ref-pdb PATH [...]
        Pocket reference PDB(s) for the final merge when --input uses XYZ/GJF.
    --out-dir PATH
        Output directory; default ./result_path_search/
    --dump {True|False}
        Save optimizer dumps; default False.
    --convert-files {True|False}
        Convert outputs into format-aware companions; PDB companions are produced for trajectories
        when references exist, and GJF companions are produced for XYZ snapshots when templates exist;
        default on.
    --freeze-links {True|False}
        Freeze parents of link hydrogens for PDB input; default True.

Examples
--------
    # Minimal (pocket-only MEP; always writes mep.trj and emits mep.pdb when inputs are PDB)
    pdb2reaction path-search -i reactant.pdb product.pdb -q 0

    # Multistep with intermediates, YAML overrides, and PDB merge to a full system
    pdb2reaction path-search -i R.pdb IM1.pdb IM2.pdb P.pdb -q -1 \
        --args-yaml params.yaml --ref-full-pdb holo_template.pdb --out-dir ./run_ps

Description
-----------
Constructs a continuous minimum‑energy path (MEP) between two or more structures ordered along a reaction.
The method runs GSM or DMF segments to localize barriers and **recursively** refines only those regions that
exhibit covalent bond changes. Kinks (no bond change) are represented by linearly interpolated, individually
optimized images. Multi‑structure inputs are processed per adjacent pair and stitched into a single MEP. A
single UMA calculator (uma_pysis) is shared serially across all stages. Configuration precedence:
**YAML > CLI > defaults**.

Workflow
--------
1) Initial path (per adjacent pair A→B): run GSM or DMF to obtain a preliminary MEP.
2) Localize barrier: find the highest‑energy image (HEI); optimize HEI±1 as single structures → two nearby minima,
   End1 and End2.
3) Refine the step:
   - If no covalent bond change between End1–End2 (a “kink”): insert `search.kink_max_nodes` linearly interpolated
     nodes and optimize each; **skip** GSM.
   - Otherwise: run a refinement GSM/DMF segment between End1 and End2. With `refine-mode=minima`, End1/End2 are taken from the
     nearest local minima flanking the HEI instead of strictly HEI±1.
4) Recurse selectively: evaluate covalent changes for (A→End1) and (End2→B); recurse only on sides that change.
5) Stitch subpaths: concatenate sub‑MEPs with duplicate removal via RMSD. If endpoints mismatch beyond
   `search.bridge_rmsd_thresh`, insert a *bridge* MEP segment (GSM/DMF).
   - If the interface itself shows covalent changes, insert a **new recursive segment** instead of a bridge.
6) Optional alignment & merge: after pre‑opt, when `--align` (default), rigidly co‑align all inputs and
   refine `freeze_atoms` to match the first input. If `--ref-full-pdb` is supplied, merge pocket trajectories
   into full templates and annotate segments (requires PDB pocket inputs or `--ref-pdb`).

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_path_search/)
  ├─ summary.yaml                    # Run-level summary (no exhaustive settings dump)
  ├─ mep.trj                         # Final MEP as XYZ 
  ├─ mep.pdb                         # Final MEP as PDB (written when inputs were PDB)
  ├─ mep_w_ref.pdb                   # Full-system merged path (requires --ref-full-pdb and pocket PDBs)
  ├─ mep_w_ref_seg_XX.pdb            # Per-segment merged paths (bond-change segments; requires --ref-full-pdb)
  ├─ mep_seg_XX.trj / mep_seg_XX.pdb # Pocket-only segment paths (bond-change segments; format follows input)
  ├─ hei_seg_XX.xyz / hei_seg_XX.pdb # Highest-energy image snapshots; hei_seg_XX.gjf when a template is available
  ├─ hei_w_ref_seg_XX.pdb            # Merged HEI per bond-change segment (requires --ref-full-pdb)
  ├─ mep_plot.png                    # ΔE profile vs. image index (from trj2fig)
  ├─ energy_diagram_MEP.png          # PNG export of the diagram when kaleido is installed
  ├─ seg_000_mep/ ...                # Initial GSM/DMF for each primary segment
  ├─ seg_000_left_<lbfgs|rfo>_opt/   # HEI-1 single-structure optimization
  ├─ seg_000_right_<lbfgs|rfo>_opt/  # HEI+1 single-structure optimization
  ├─ seg_000_refine_mep/ ...         # Refinement GSM/DMF (bond-change segments)
  ├─ seg_000_kink_...                # Kink interpolation optimizations (when applicable)
  ├─ seg_000_seg_002_bridge_mep/ ... # Bridge GSM/DMF segments; path indicates bridged segments
  ├─ seg_001_...                     # Left-side recursive substeps (if any)
  └─ seg_002_...                     # Right-side recursive substeps (if any)

Notes
-----
- Inputs:
  - Provide ≥2 structures to `-i/--input` in reaction order. Either repeat `-i` or pass multiple paths after a single `-i`.
  - Endpoint pre‑optimization runs by default (`--preopt True`). With `--align True`, all inputs are co‑aligned and
    `freeze_atoms` are refined to the first input.
- Bond‑change detection: uses `bond_changes.compare_structures` with `bond_factor`, `margin_fraction`,
  and `delta_fraction` thresholds.
- Concatenation policy:
  - Endpoint duplicate removal when RMSD ≤ `search.stitch_rmsd_thresh` (default 1e‑4 Å).
  - Bridge GSM/DMF segment when gap RMSD > `search.bridge_rmsd_thresh` (default 1e‑4 Å).
  - If an interface shows covalent changes, insert a **new recursive segment** instead of a bridge.
- Nodes and recursion:
  - Segment vs bridge nodes can differ via `search.max_nodes_segment` and `search.max_nodes_bridge` (segment defaults to `--max-nodes`).
  - Kinks use `search.kink_max_nodes` (default 3) linear nodes, each optimized individually.
  - Recursion depth is capped by `search.max_depth` (default 10).
- Calculators & optimizers:
  - A single UMA calculator (`uma_pysis`, default model "uma-s-1p1") is shared serially across all stages.
  - GSM employs pysisyphus `GrowingString` + `StringOptimizer`; DMF uses the Direct Max Flux interpolator.
    Single‑structure optimization uses LBFGS or RFO.
- Final merge rule with `--align True`: when `--ref-full-pdb` is provided, the **first** reference PDB is used for *all* pairs
  in the final merge (passing one file is sufficient). Without `--align`, supply one reference PDB per input.
- Console output prints the linear state sequence (e.g., `R --> TS1 --> IM1_1 -|--> IM1_2 --> ... --> P`) and the exact
  labels/energies used to build the energy diagram.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Callable

import sys
import traceback
import textwrap
import tempfile
import os
import time
import re

import click
import numpy as np
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import AU2KCALPERMOL, BOHR2ANG

# Biopython for PDB I/O and parsing
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .path_opt import (
    _maybe_convert_to_pdb,
    _optimize_single,
    _run_dmf_mep,
    _write_ase_trj_with_energy,
    DMF_KW as _PATH_DMF_KW,
    GS_KW as _PATH_GS_KW,
    STOPT_KW as _PATH_STOPT_KW,
)
from .opt import (
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
)
from .utils import (
    detect_freeze_links_safe,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    build_energy_diagram,
    prepare_input_structure,
    fill_charge_spin_from_gjf,
    _derive_charge_from_ligand_charge,
    maybe_convert_xyz_to_gjf,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    PreparedInputStructure,
    GjfTemplate,
)
from .summary_log import write_summary_log
from .trj2fig import run_trj2fig
from .bond_changes import compare_structures, summarize_changes
from .align_freeze_atoms import align_and_refine_sequence_inplace, kabsch_R_t

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _close_matplotlib_figures() -> None:
    """Best-effort cleanup for matplotlib figures to avoid open-figure warnings."""

    try:
        import matplotlib.pyplot as plt  # type: ignore

        plt.close("all")
    except Exception:
        pass


# YAML helper to preserve multiline blocks for bond-change summaries in summary.yaml
class _LiteralStr(str):
    """String marker to force literal block style when dumping YAML."""


def _literal_str_representer(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")


yaml.add_representer(_LiteralStr, _literal_str_representer)
yaml.add_representer(_LiteralStr, _literal_str_representer, Dumper=yaml.SafeDumper)


def _bond_changes_block(text: Optional[str]):
    """
    Prepare bond-change summaries for YAML output, emitting structured lists instead
    of escaped ``\n`` strings when possible.
    """

    if text is None:
        return ""

    cleaned = str(text).strip()
    if not cleaned:
        return ""

    # Convert the ``summarize_changes`` output
    #   Bond formed (1):
    #     - C1-C2 : 3.0 Å --> 1.5 Å
    #   Bond broken: None
    # into a YAML-friendly nested list to avoid embedded newlines.
    lines = cleaned.splitlines()
    sections = []
    title: Optional[str] = None
    entries: List[str] = []

    def _flush() -> None:
        nonlocal title, entries
        if title is None:
            return
        payload = entries if entries else ["None"]
        sections.append({title: payload})
        title, entries = None, []

    for ln in lines:
        stripped = ln.strip()
        if stripped.startswith("Bond "):
            _flush()
            if stripped.endswith(": None"):
                sections.append({stripped[:-len(": None")]: ["None"]})
            else:
                title = stripped.rstrip(":")
        elif stripped.startswith("- "):
            entries.append(stripped[2:])

    _flush()

    # Fallback to literal block if the format is unexpected
    if sections:
        return sections
    if "\n" in cleaned:
        return _LiteralStr(cleaned)
    return cleaned


# -----------------------------------------------
# Configuration defaults
# -----------------------------------------------

# Geometry (input handling)
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)

CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

# DMF (Direct Max Flux + (C)FB-ENM)
DMF_KW: Dict[str, Any] = dict(_PATH_DMF_KW)

# GrowingString (path representation)
GS_KW: Dict[str, Any] = dict(_PATH_GS_KW)

# StringOptimizer (GSM optimization control)
STOPT_KW: Dict[str, Any] = dict(_PATH_STOPT_KW)
STOPT_KW.update({
    "out_dir": "./result_path_search/",  # output directory for string-optimizer artifacts
})

# LBFGS settings
LBFGS_KW: Dict[str, Any] = dict(_LBFGS_KW)
LBFGS_KW.update({
    "out_dir": "./result_path_search/",  # LBFGS output directory (restart, logs)
})

# RFO settings
RFO_KW: Dict[str, Any] = dict(_RFO_KW)
RFO_KW.update({
    "out_dir": "./result_path_search/",  # RFO output directory (restart, logs)
})

# Covalent‑bond change detection
BOND_KW: Dict[str, Any] = {
    "device": "cuda",               # compute UMA graph features on CUDA when available
    "bond_factor": 1.20,            # covalent cutoff multiplier (r_cov * bond_factor)
    "margin_fraction": 0.05,        # tolerance margin added to bond cutoff for stability
    "delta_fraction": 0.05,         # threshold fraction to flag bond formation/breaking
}

# Global search control
SEARCH_KW: Dict[str, Any] = {
    "max_depth": 10,                 # recursion depth for path-search branching
    "stitch_rmsd_thresh": 1.0e-4,    # RMSD cutoff when stitching partial segments
    "bridge_rmsd_thresh": 1.0e-4,    # RMSD cutoff when bridging paths
    "rmsd_align": True,              # retained for compatibility (ignored internally)
    "max_nodes_segment": 10,         # max nodes per segment during expansion
    "max_nodes_bridge": 5,           # max nodes per bridge construction
    "kink_max_nodes": 3,             # max nodes for kink resolution
    "max_seq_kink": 2,               # max sequential kinks allowed before aborting
    "refine_mode": None,             # optional refinement strategy tag
}


def _load_two_endpoints(
    paths: Sequence[Path],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> Sequence:
    """
    Load two or more geometries and assign `freeze_atoms`; return pysisyphus geometries.
    """
    geoms = []
    for p in paths:
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        if auto_freeze_links and p.suffix.lower() == ".pdb":
            detected = detect_freeze_links_safe(p)
            freeze = merge_freeze_atom_indices(cfg, detected)
            if detected and freeze:
                click.echo(f"[freeze-links] {p.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        else:
            freeze = merge_freeze_atom_indices(cfg)
        g = geom_loader(p, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# Multi‑structure loader
def _load_structures(
    inputs: Sequence[PreparedInputStructure],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> List[Any]:
    """
    Load multiple geometries and assign `freeze_atoms`; return a list of geometries.
    """
    geoms: List[Any] = []
    for prepared in inputs:
        geom_path = prepared.geom_path
        src_path = prepared.source_path
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        if auto_freeze_links and src_path.suffix.lower() == ".pdb":
            detected = detect_freeze_links_safe(src_path)
            freeze = merge_freeze_atom_indices(cfg, detected)
            if detected and freeze:
                click.echo(f"[freeze-links] {src_path.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        else:
            freeze = merge_freeze_atom_indices(cfg)
        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


def _write_xyz_trj_with_energy(images: Sequence, energies: Sequence[float], path: Path) -> None:
    """
    Write an XYZ `.trj` with the energy on line 2 of each block.
    """
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


def _maybe_convert_to_gjf(
    xyz_path: Path,
    template: Optional[GjfTemplate],
    out_path: Optional[Path] = None,
) -> Optional[Path]:
    """
    Convert XYZ to Gaussian input using a template, when available.
    """
    try:
        if template is None or (not xyz_path.exists()):
            return None
        target = out_path if out_path is not None else xyz_path.with_suffix(".gjf")
        maybe_convert_xyz_to_gjf(xyz_path, template, target)
        click.echo(f"[convert] Wrote '{target}'.")
        return target
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert '{xyz_path.name}' to GJF: {e}", err=True)
        return None


def _kabsch_rmsd(A: np.ndarray, B: np.ndarray, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """
    RMSD between A and B (no rigid alignment; `align` is ignored). Optional subset selection via `indices`.
    """
    assert A.shape == B.shape and A.shape[1] == 3
    if indices is not None and len(indices) > 0:
        idx = np.array(sorted({int(i) for i in indices if 0 <= int(i) < A.shape[0]}), dtype=int)
        if idx.size == 0:
            idx = np.arange(A.shape[0], dtype=int)
        A = A[idx]
        B = B[idx]
    diff = A - B
    return float(np.sqrt((diff * diff).sum() / A.shape[0]))


def _rmsd_between(ga, gb, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """
    RMSD between two pysisyphus Geometries (no alignment; optional subset selection).
    """
    return _kabsch_rmsd(np.array(ga.coords3d), np.array(gb.coords3d), align=False, indices=indices)


def _has_bond_change(x, y, bond_cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Determine if covalent bonds form/break between `x` and `y`. Returns (changed?, summary_text).
    """
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


def _gs_cfg_with_overrides(base: Dict[str, Any], **overrides: Any) -> Dict[str, Any]:
    """
    Shallow copy of a GS config with specified overrides.
    """
    cfg = dict(base)
    for k, v in overrides.items():
        cfg[k] = v
    return cfg


# -----------------------------------------------
# Kink detection & interpolation helpers
# -----------------------------------------------

def _max_displacement_between(ga, gb, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """
    Maximum per‑atom displacement (Å) between two structures (no alignment).
    """
    A = np.asarray(ga.coords3d, dtype=float)
    B = np.asarray(gb.coords3d, dtype=float)
    if A.shape != B.shape or A.shape[1] != 3:
        raise ValueError("Geometries must have the same number of atoms for displacement.")
    disp = np.linalg.norm(A - B, axis=1)
    return float(np.max(disp))


def _new_geom_from_coords(atoms: Sequence[str], coords: np.ndarray, coord_type: str, freeze_atoms: Sequence[int]) -> Any:
    """
    Create a pysisyphus Geometry from Bohr coords via temporary XYZ; attach `freeze_atoms`.
    """
    lines = [str(len(atoms)), ""]
    coords_ang = np.asarray(coords, dtype=float) * BOHR2ANG
    for sym, (x, y, z) in zip(atoms, coords_ang):
        lines.append(f"{sym} {x:.15f} {y:.15f} {z:.15f}")
    s = "\n".join(lines) + "\n"
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()
        g = geom_loader(
            Path(tmp.name), coord_type=coord_type, freeze_atoms=freeze_atoms
        )
        g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _make_linear_interpolations(gL, gR, n_internal: int) -> List[Any]:
    """
    Return `n_internal` linearly interpolated structures between gL → gR (excluding endpoints).
    Atom order follows `gL`.
    """
    A = np.asarray(gL.coords3d, dtype=float)
    B = np.asarray(gR.coords3d, dtype=float)
    assert A.shape == B.shape and A.shape[1] == 3, "Atom counts must match for interpolation."
    atoms = [a for a in gL.atoms]
    coord_type = gL.coord_type
    faL = getattr(gL, "freeze_atoms", None)
    faR = getattr(gR, "freeze_atoms", None)

    def _as_list(raw: Any) -> List[Any]:
        if raw is None:
            return []
        try:
            return list(raw)
        except Exception:
            return []

    freeze_union = sorted(
        set(map(int, _as_list(faL))) | set(map(int, _as_list(faR)))
    )
    interps: List[Any] = []
    for k in range(1, n_internal + 1):
        t = k / (n_internal + 1.0)
        C = (1.0 - t) * A + t * B
        interps.append(_new_geom_from_coords(atoms, C, coord_type, freeze_union))
    return interps


def _energy_of(g) -> float:
    """
    Return the energy (Hartree) of a Geometry (ensures calculator is attached).
    """
    g.set_calculator(getattr(g, "calculator", None))
    return float(g.energy)


# ---- Segment/bridge tagging helpers ----

def _tag_images(images: Sequence[Any], **attrs: Any) -> None:
    """
    Attach arbitrary attributes to Geometry images.
    """
    for im in images:
        for k, v in attrs.items():
            try:
                setattr(im, k, v)
            except Exception:
                pass


def _segment_base_id(tag: str) -> str:
    """
    Extract base id 'seg_XXX' from a tag like 'seg_000_refine'; fallback to `tag` or 'seg'.
    """
    m = re.search(r"(seg_\d{3})", tag or "")
    return m.group(1) if m else (tag or "seg")


def _is_local_minimum(idx: int, energies: Sequence[float]) -> bool:
    """
    Return True if the point at `idx` is a local minimum according to neighboring energies.

    Definition:
    - Interior point: both neighbors have higher energy.
    - Endpoint: the single neighbor has higher energy.
    """

    if idx < 0 or idx >= len(energies):
        return False
    if idx == 0:
        return len(energies) > 1 and energies[1] > energies[0]
    if idx == len(energies) - 1:
        return energies[-2] > energies[-1]
    return energies[idx - 1] > energies[idx] and energies[idx + 1] > energies[idx]


def _find_nearest_local_minimum(hei_idx: int, direction: int, energies: Sequence[float]) -> Optional[int]:
    """
    Starting from the HEI, walk `direction` (-1 for left, +1 for right) to find the nearest
    local minimum. Returns the index of the first minimum encountered, or None if absent.
    """

    i = hei_idx + direction
    while 0 <= i < len(energies):
        if _is_local_minimum(i, energies):
            return i
        i += direction
    return None


@dataclass
class GSMResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int


# ---- Per‑segment summary for the console report ----
@dataclass
class SegmentReport:
    tag: str
    barrier_kcal: float
    delta_kcal: float
    summary: str  # summarize_changes string (empty for bridges)
    kind: str = "seg"          # "seg" or "bridge"
    seg_index: int = 0         # 1‑based index along final MEP (assigned later)


def _run_mep_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    gjf_template: Optional[GjfTemplate] = None,
    prepared_input: Optional[PreparedInputStructure] = None,
) -> GSMResult:
    """
    Run GSM between `gA`–`gB`, save segment outputs, and return images/energies/HEI index.
    """
    if gjf_template is None:
        try:
            gjf_template = _PRIMARY_GJF_TEMPLATE  # type: ignore[name-defined]
        except NameError:
            pass
    for g in (gA, gB):
        g.set_calculator(shared_calc)

    def calc_getter():
        return shared_calc

    gs = GrowingString(
        images=[gA, gB],
        calc_getter=calc_getter,
        **gs_cfg,
    )

    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)
    _opt_args = dict(opt_cfg)
    _opt_args["out_dir"] = str(seg_dir)

    optimizer = StringOptimizer(
        geometry=gs,
        **{k: v for k, v in _opt_args.items() if k != "type"}
    )

    click.echo(f"\n=== [{tag}] GSM started ===\n")
    optimizer.run()
    click.echo(f"\n=== [{tag}] GSM finished ===\n")

    energies = list(map(float, np.array(gs.energy, dtype=float)))
    images = list(gs.images)

    # Choose HEI: prefer internal local maxima; fallback to highest internal node
    E = np.array(energies, dtype=float)
    nE = len(E)
    local_max_candidates = [i for i in range(1, nE - 1) if (E[i] > E[i - 1] and E[i] > E[i + 1])]
    if local_max_candidates:
        hei_idx = int(max(local_max_candidates, key=lambda i: E[i]))
    else:
        hei_idx = int(np.argmax(E[1:-1])) + 1 if nE >= 3 else int(np.argmax(E))

    # Write trajectory
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

    # Energy plot for the segment
    try:
        if wrote_with_energy:
            run_trj2fig(final_trj, [seg_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            _close_matplotlib_figures()
            click.echo(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'")
        else:
            click.echo(f"[{tag}] WARNING: Energies missing; skipping plot.", err=True)
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    # Convert trajectory and HEI outputs based on the input template
    prepared_for_outputs = prepared_input
    ref_for_conv = ref_pdb_path
    if prepared_for_outputs is None and gjf_template is not None:
        prepared_for_outputs = PreparedInputStructure(
            source_path=gjf_template.path,
            geom_path=gjf_template.path,
            gjf_template=gjf_template,
        )
    if prepared_for_outputs is None and ref_for_conv is not None:
        prepared_for_outputs = PreparedInputStructure(
            source_path=ref_for_conv,
            geom_path=ref_for_conv,
        )
    if prepared_for_outputs is not None and ref_for_conv is None:
        if prepared_for_outputs.source_path.suffix.lower() == ".pdb":
            ref_for_conv = prepared_for_outputs.source_path.resolve()

    needs_pdb = ref_for_conv is not None
    needs_gjf = bool(prepared_for_outputs and prepared_for_outputs.is_gjf)

    if prepared_for_outputs is not None and (needs_pdb or needs_gjf):
        try:
            convert_xyz_like_outputs(
                final_trj,
                prepared_for_outputs,
                ref_pdb_path=ref_for_conv,
                out_pdb_path=seg_dir / "final_geometries.pdb" if needs_pdb else None,
                out_gjf_path=seg_dir / "final_geometries.gjf" if needs_gjf else None,
            )
        except Exception as e:
            click.echo(
                f"[{tag}] WARNING: Failed to convert segment trajectory: {e}", err=True
            )

    # Write HEI structure
    try:
        hei_geom = images[hei_idx]
        hei_E = float(E[hei_idx])
        hei_xyz = seg_dir / "hei.xyz"
        s = hei_geom.as_xyz()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{hei_E:.12f}"
        s_out = "\n".join(lines)
        if not s_out.endswith("\n"):
            s_out += "\n"
        with open(hei_xyz, "w") as f:
            f.write(s_out)
        click.echo(f"[{tag}] Wrote '{hei_xyz}'.")

        if prepared_for_outputs is not None and (needs_pdb or needs_gjf):
            try:
                convert_xyz_like_outputs(
                    hei_xyz,
                    prepared_for_outputs,
                    ref_pdb_path=ref_for_conv,
                    out_pdb_path=seg_dir / "hei.pdb" if needs_pdb else None,
                    out_gjf_path=seg_dir / "hei.gjf" if needs_gjf else None,
                )
            except Exception as e:
                click.echo(
                    f"[{tag}] WARNING: Failed to convert HEI structure: {e}",
                    err=True,
                )
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to write HEI structure: {e}", err=True)

    return GSMResult(images=images, energies=energies, hei_idx=hei_idx)


def _ase_atoms_to_geom(atoms, coord_type: str, template_g=None, shared_calc=None):
    """
    Convert ASE Atoms → pysisyphus Geometry, preserving freeze_atoms when present.
    """

    from ase.io import write as ase_write

    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        ase_write(tmp, atoms, format="xyz")
        tmp.flush()
        tmp.close()
        g = geom_loader(
            Path(tmp.name),
            coord_type=coord_type,
            freeze_atoms=getattr(template_g, "freeze_atoms", []),
        )
        try:
            g.freeze_atoms = np.array(getattr(template_g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        if shared_calc is not None:
            g.set_calculator(shared_calc)
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _run_dmf_between(
    gA,
    gB,
    calc_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    max_nodes: int,
    prepared_inputs: Sequence[PreparedInputStructure],
    shared_calc,
    dmf_cfg: Dict[str, Any],
) -> GSMResult:
    """
    Run DMF for a segment and convert outputs to pysisyphus Geometries.
    """

    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)

    fix_atoms: List[int] = []
    try:
        fix_atoms = sorted(
            {int(i) for g in [gA,gB] for i in getattr(g, "freeze_atoms", [])}
        )
    except Exception:
        pass

    dmf_res = _run_dmf_mep(
        geoms=[gA, gB],
        calc_cfg=calc_cfg,
        out_dir_path=seg_dir,
        prepared_inputs=prepared_inputs,
        max_nodes=max_nodes,
        fix_atoms=fix_atoms,
        dmf_cfg=dmf_cfg,
    )

    energies = list(map(float, dmf_res.energies))

    final_trj = seg_dir / "final_geometries.trj"
    _write_ase_trj_with_energy(dmf_res.images, energies, final_trj)
    _maybe_convert_to_pdb(final_trj, ref_pdb_path)

    try:
        run_trj2fig(final_trj, [seg_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
        _close_matplotlib_figures()
        click.echo(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'")
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    imgs: List[Any] = []
    for atoms in dmf_res.images:
        imgs.append(
            _ase_atoms_to_geom(atoms, coord_type=gA.coord_type, template_g=gA, shared_calc=shared_calc)
        )

    return GSMResult(images=imgs, energies=energies, hei_idx=int(dmf_res.hei_idx))


def _refine_between(
    gL,
    gR,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    mep_mode_kind: str,
    calc_cfg: Dict[str, Any],
    max_nodes: int,
    prepared_inputs: Sequence[PreparedInputStructure],
    dmf_cfg: Dict[str, Any],
) -> GSMResult:
    """
    Refine End1–End2 via GSM or DMF depending on the selected mode.
    """

    if mep_mode_kind == "dmf":
        return _run_dmf_between(
            gL,
            gR,
            calc_cfg,
            out_dir,
            tag=f"{tag}_refine",
            ref_pdb_path=ref_pdb_path,
            max_nodes=max_nodes,
            prepared_inputs=prepared_inputs,
            shared_calc=shared_calc,
            dmf_cfg=dmf_cfg,
        )

    return _run_mep_between(gL, gR, shared_calc, gs_cfg, opt_cfg, out_dir, tag=f"{tag}_refine", ref_pdb_path=ref_pdb_path)


def _maybe_bridge_segments(
    tail_g,
    head_g,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    rmsd_thresh: float,
    ref_pdb_path: Optional[Path],
    mep_mode_kind: str,
    calc_cfg: Dict[str, Any],
    max_nodes: int,
    prepared_inputs: Sequence[PreparedInputStructure],
    dmf_cfg: Dict[str, Any],
) -> Optional[GSMResult]:
    """
    Run a bridge GSM if two segment endpoints are farther than the threshold.
    """
    rmsd = _rmsd_between(tail_g, head_g, align=False)
    if rmsd <= rmsd_thresh:
        return None
    click.echo(
        f"[{tag}] Gap detected between segments (RMSD={rmsd:.4e} Å) — bridging via {mep_mode_kind.upper()}."
    )
    if mep_mode_kind == "dmf":
        return _run_dmf_between(
            tail_g,
            head_g,
            calc_cfg,
            out_dir,
            tag=f"{tag}_bridge",
            ref_pdb_path=ref_pdb_path,
            max_nodes=max_nodes,
            prepared_inputs=prepared_inputs,
            shared_calc=shared_calc,
            dmf_cfg=dmf_cfg,
        )

    return _run_mep_between(tail_g, head_g, shared_calc, gs_cfg, opt_cfg, out_dir, tag=f"{tag}_bridge", ref_pdb_path=ref_pdb_path)


def _stitch_paths(
    parts: List[Tuple[List[Any], List[float]]],
    stitch_rmsd_thresh: float,
    bridge_rmsd_thresh: float,
    shared_calc,
    gs_cfg,
    opt_cfg,
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    bond_cfg: Optional[Dict[str, Any]] = None,
    segment_builder: Optional[Callable[[Any, Any, str], "CombinedPath"]] = None,
    segments_out: Optional[List["SegmentReport"]] = None,
    bridge_pair_index: Optional[int] = None,
    mep_mode_kind: str = "dmf",
    calc_cfg: Optional[Dict[str, Any]] = None,
    max_nodes: int = 10,
    prepared_inputs: Optional[Sequence[PreparedInputStructure]] = None,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[List[Any], List[float]]:
    """
    Concatenate path parts (images, energies). Insert bridge GSMs when needed.
    If covalent changes are detected across an interface, build and insert a *new* recursive segment
    using `segment_builder` instead of bridging. Update `segments_out` accordingly.
    """
    dmf_cfg = dict(dmf_cfg or DMF_KW)
    all_imgs: List[Any] = []
    all_E: List[float] = []

    def _last_known_seg_tag_from_images(imgs: List[Any]) -> Optional[str]:
        for im in reversed(imgs):
            t = getattr(im, "mep_seg_tag", None)
            if t:
                return t
        return None

    def _first_known_seg_tag_from_images(imgs: List[Any]) -> Optional[str]:
        for im in imgs:
            t = getattr(im, "mep_seg_tag", None)
            if t:
                return t
        return None

    def append_part(imgs: List[Any], Es: List[float]) -> None:
        nonlocal all_imgs, all_E
        if not imgs:
            return
        if not all_imgs:
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return
        tail = all_imgs[-1]
        head = imgs[0]

        adj_changed, adj_summary = False, ""
        if segment_builder is not None and bond_cfg is not None:
            try:
                adj_changed, adj_summary = _has_bond_change(tail, head, bond_cfg)
            except Exception:
                adj_changed, adj_summary = False, ""

        if adj_changed and segment_builder is not None:
            click.echo(f"[{tag}] Covalent changes detected at interface — inserting a new recursive segment.")
            if adj_summary:
                click.echo(textwrap.indent(adj_summary, prefix="  "))
            sub = segment_builder(tail, head, f"{tag}_mid")
            seg_imgs, seg_E = sub.images, sub.energies
            if segments_out is not None and getattr(sub, "segments", None):
                segments_out.extend(sub.segments)
            if seg_imgs:
                if _rmsd_between(all_imgs[-1], seg_imgs[0], align=False) <= stitch_rmsd_thresh:
                    seg_imgs = seg_imgs[1:]
                    seg_E = seg_E[1:]
                all_imgs.extend(seg_imgs)
                all_E.extend(seg_E)
            if _rmsd_between(all_imgs[-1], imgs[0], align=False) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return

        rmsd = _rmsd_between(tail, head, align=False)
        if rmsd <= stitch_rmsd_thresh:
            all_imgs.extend(imgs[1:])
            all_E.extend(Es[1:])
        elif rmsd > bridge_rmsd_thresh:
            left_tag_recent = _last_known_seg_tag_from_images(all_imgs) or "segL"
            right_tag_upcoming = _first_known_seg_tag_from_images(imgs) or "segR"
            left_base = _segment_base_id(left_tag_recent)
            right_base = _segment_base_id(right_tag_upcoming)
            bridge_name_base = f"{left_base}_{right_base}"

            br = _maybe_bridge_segments(
                tail, head, shared_calc, gs_cfg, opt_cfg, out_dir, tag=bridge_name_base,
                rmsd_thresh=bridge_rmsd_thresh, ref_pdb_path=ref_pdb_path,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg or {}, max_nodes=max_nodes,
                prepared_inputs=prepared_inputs or [],
                dmf_cfg=dmf_cfg,
            )
            if br is not None:
                _tag_images(br.images, mep_seg_tag=f"{bridge_name_base}_bridge", mep_seg_kind="bridge",
                            mep_has_bond_changes=False, pair_index=bridge_pair_index)
                b_imgs, b_E = br.images, br.energies
                if _rmsd_between(all_imgs[-1], b_imgs[0], align=False) <= stitch_rmsd_thresh:
                    b_imgs = b_imgs[1:]
                    b_E = b_E[1:]
                if b_imgs:
                    all_imgs.extend(b_imgs)
                    all_E.extend(b_E)

                if segments_out is not None:
                    try:
                        barrier_kcal = (max(br.energies) - br.energies[0]) * AU2KCALPERMOL
                        delta_kcal = (br.energies[-1] - br.energies[0]) * AU2KCALPERMOL
                    except Exception:
                        barrier_kcal = float("nan")
                        delta_kcal = float("nan")
                    bridge_report = SegmentReport(
                        tag=f"{bridge_name_base}_bridge",
                        barrier_kcal=float(barrier_kcal),
                        delta_kcal=float(delta_kcal),
                        summary="",
                        kind="bridge"
                    )
                    insert_pos: Optional[int] = None
                    try:
                        for j, sr in enumerate(segments_out):
                            if sr.tag == right_tag_upcoming:
                                insert_pos = j
                                break
                    except Exception:
                        insert_pos = None
                    if insert_pos is None:
                        segments_out.append(bridge_report)
                    else:
                        segments_out.insert(insert_pos, bridge_report)

            if _rmsd_between(all_imgs[-1], imgs[0], align=False) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
        else:
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
    segments: List[SegmentReport]


def _trailing_kink_count(segments: Sequence[SegmentReport]) -> int:
    """Return the number of consecutive kink segments at the end of ``segments``."""

    count = 0
    for seg in reversed(segments):
        if seg.tag and "kink" in seg.tag:
            count += 1
        else:
            break
    return count


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
    refine_mode_kind: str,
    mep_mode_kind: str,
    calc_cfg: Dict[str, Any],
    dmf_cfg: Dict[str, Any],
    prepared_inputs: Sequence[PreparedInputStructure],
    out_dir: Path,
    ref_pdb_path: Optional[Path],
    prepared_input: Optional[PreparedInputStructure],
    depth: int,
    seg_counter: List[int],
    branch_tag: str,
    pair_index: Optional[int] = None,
    kink_seq_count: int = 0,
) -> CombinedPath:
    """
    Recursively construct a multistep MEP from A–B and return it (A→B order).
    """
    seg_max_nodes = int(search_cfg.get("max_nodes_segment", gs_cfg.get("max_nodes", 10)))
    gs_seg_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=seg_max_nodes)
    max_seq_kink = int(search_cfg.get("max_seq_kink", 2))

    def _terminate_with_maxdepth(reason_msg: Optional[str] = None) -> CombinedPath:
        if reason_msg:
            click.echo(reason_msg)

        seg_tag = f"seg_{seg_counter[0]:03d}_maxdepth"
        gsm = (
            _run_dmf_between(
                gA,
                gB,
                calc_cfg,
                out_dir,
                tag=seg_tag,
                ref_pdb_path=ref_pdb_path,
                max_nodes=seg_max_nodes,
                prepared_inputs=prepared_inputs,
                shared_calc=shared_calc,
                dmf_cfg=dmf_cfg,
            )
            if mep_mode_kind == "dmf"
            else _run_mep_between(
                gA,
                gB,
                shared_calc,
                gs_seg_cfg,
                opt_cfg,
                out_dir,
                tag=seg_tag,
                ref_pdb_path=ref_pdb_path,
                prepared_input=prepared_input,
            )
        )
        seg_counter[0] += 1

        try:
            changed, step_summary = _has_bond_change(gsm.images[0], gsm.images[-1], bond_cfg)
        except Exception as e:
            click.echo(f"[{seg_tag}] WARNING: Failed to evaluate bond changes at max depth: {e}", err=True)
            changed, step_summary = True, ""

        try:
            barrier_kcal = (max(gsm.energies) - gsm.energies[0]) * AU2KCALPERMOL
            delta_kcal = (gsm.energies[-1] - gsm.energies[0]) * AU2KCALPERMOL
        except Exception:
            barrier_kcal = float("nan")
            delta_kcal = float("nan")

        seg_report = SegmentReport(
            tag=seg_tag,
            barrier_kcal=float(barrier_kcal),
            delta_kcal=float(delta_kcal),
            summary=step_summary if changed else "(no covalent changes detected)",
            kind="seg",
        )

        _tag_images(
            gsm.images,
            mep_seg_tag=seg_tag,
            mep_seg_kind="seg",
            mep_has_bond_changes=bool(changed),
            pair_index=pair_index,
        )

        return CombinedPath(images=gsm.images, energies=gsm.energies, segments=[seg_report])

    if depth > int(search_cfg.get("max_depth", 10)):
        click.echo(f"[{branch_tag}] Reached maximum recursion depth. Returning current endpoints only.")
        return _terminate_with_maxdepth()

    seg_id = seg_counter[0]
    seg_counter[0] += 1
    tag0 = f"seg_{seg_id:03d}"

    gsm0 = (
        _run_dmf_between(
            gA,
            gB,
            calc_cfg,
            out_dir,
            tag=tag0,
            ref_pdb_path=ref_pdb_path,
            max_nodes=seg_max_nodes,
            prepared_inputs=prepared_inputs,
            shared_calc=shared_calc,
            dmf_cfg=dmf_cfg,
        )
        if mep_mode_kind == "dmf"
        else _run_mep_between(
            gA,
            gB,
            shared_calc,
            gs_seg_cfg,
            opt_cfg,
            out_dir,
            tag=tag0,
            ref_pdb_path=ref_pdb_path,
            prepared_input=prepared_input,
        )
    )

    hei = int(gsm0.hei_idx)
    if not (1 <= hei <= len(gsm0.images) - 2):
        click.echo(f"[{tag0}] WARNING: HEI is at an endpoint (idx={hei}). Returning the raw GSM path.")
        _tag_images(gsm0.images, pair_index=pair_index)
        return CombinedPath(images=gsm0.images, energies=gsm0.energies, segments=[])

    if refine_mode_kind == "minima":
        left_idx = _find_nearest_local_minimum(hei_idx=hei, direction=-1, energies=gsm0.energies)
        right_idx = _find_nearest_local_minimum(hei_idx=hei, direction=1, energies=gsm0.energies)

        if left_idx is None:
            left_idx = hei - 1
        if right_idx is None:
            right_idx = hei + 1

        click.echo(
            f"[{tag0}] Using nearest local minima around HEI (left idx={left_idx}, right idx={right_idx})."
        )
        left_img = gsm0.images[left_idx]
        right_img = gsm0.images[right_idx]
    else:
        left_img = gsm0.images[hei - 1]
        right_img = gsm0.images[hei + 1]
        click.echo(f"[{tag0}] Refining HEI±1 (peak mode).")

    left_end = _optimize_single(
        left_img,
        shared_calc,
        sopt_kind,
        sopt_cfg,
        out_dir,
        tag=f"{tag0}_left",
        prepared_input=prepared_input,
    )
    right_end = _optimize_single(
        right_img,
        shared_calc,
        sopt_kind,
        sopt_cfg,
        out_dir,
        tag=f"{tag0}_right",
        prepared_input=prepared_input,
    )

    try:
        lr_changed, lr_summary = _has_bond_change(left_end, right_end, bond_cfg)
    except Exception as e:
        click.echo(f"[{tag0}] WARNING: Failed to evaluate bond changes for kink detection: {e}", err=True)
        lr_changed, lr_summary = True, ""
    use_kink = (not lr_changed)

    if use_kink:
        n_inter = int(search_cfg.get("kink_max_nodes", 3))
        click.echo(f"[{tag0}] Kink detected (no covalent changes between End1 and End2). "
                   f"Using {n_inter} linear interpolation nodes + single-structure optimizations instead of GSM.")
        inter_geoms = _make_linear_interpolations(left_end, right_end, n_inter)
        opt_inters: List[Any] = []
        for i, g_int in enumerate(inter_geoms, 1):
            g_int.set_calculator(shared_calc)
            g_opt = _optimize_single(
                g_int,
                shared_calc,
                sopt_kind,
                sopt_cfg,
                out_dir,
                tag=f"{tag0}_kink_int{i}",
                prepared_input=prepared_input,
            )
            opt_inters.append(g_opt)
        step_imgs = [left_end] + opt_inters + [right_end]
        step_E = [float(img.energy) for img in step_imgs]
        ref1 = GSMResult(images=step_imgs, energies=step_E, hei_idx=int(np.argmax(step_E)))
        step_tag_for_report = f"{tag0}_kink"
    else:
        click.echo(f"[{tag0}] Kink not detected (covalent changes present between End1 and End2).")
        if lr_summary:
            click.echo(textwrap.indent(lr_summary, prefix="  "))
        ref1 = _refine_between(
            left_end,
            right_end,
            shared_calc,
            gs_seg_cfg,
            opt_cfg,
            out_dir,
            tag=tag0,
            ref_pdb_path=ref_pdb_path,
            mep_mode_kind=mep_mode_kind,
            calc_cfg=calc_cfg,
            max_nodes=seg_max_nodes,
            prepared_inputs=prepared_inputs,
            dmf_cfg=dmf_cfg,
        )
        step_tag_for_report = f"{tag0}_refine"

    step_imgs, step_E = ref1.images, ref1.energies

    _changed, step_summary = _has_bond_change(step_imgs[0], step_imgs[-1], bond_cfg)
    _tag_images(step_imgs, mep_seg_tag=step_tag_for_report, mep_seg_kind="seg",
                mep_has_bond_changes=bool(_changed), pair_index=pair_index)

    left_changed, left_summary = _has_bond_change(gA, left_end, bond_cfg)
    right_changed, right_summary = _has_bond_change(right_end, gB, bond_cfg)

    click.echo(f"[{tag0}] Covalent changes (A vs left_end): {'Yes' if left_changed else 'No'}")
    if left_changed:
        click.echo(textwrap.indent(left_summary, prefix="  "))
    click.echo(f"[{tag0}] Covalent changes (right_end vs B): {'Yes' if right_changed else 'No'}")
    if right_changed:
        click.echo(textwrap.indent(right_summary, prefix="  "))

    try:
        barrier_kcal = (max(step_E) - step_E[0]) * AU2KCALPERMOL
        delta_kcal = (step_E[-1] - step_E[0]) * AU2KCALPERMOL
    except Exception:
        barrier_kcal = float("nan")
        delta_kcal = float("nan")

    seg_report = SegmentReport(
        tag=step_tag_for_report,
        barrier_kcal=float(barrier_kcal),
        delta_kcal=float(delta_kcal),
        summary=step_summary if _changed else "(no covalent changes detected)",
        kind="seg"
    )

    parts: List[Tuple[List[Any], List[float]]] = []
    seg_reports: List[SegmentReport] = []

    trailing_kink_run = kink_seq_count
    if left_changed:
        subL = _build_multistep_path(
            gA, left_end, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg, bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
            out_dir, ref_pdb_path, prepared_input, depth + 1, seg_counter, branch_tag=f"{branch_tag}L",
            pair_index=pair_index,
            kink_seq_count=kink_seq_count,
        )
        _tag_images(subL.images, pair_index=pair_index)
        parts.append((subL.images, subL.energies))
        seg_reports.extend(subL.segments)
        trailing_kink_run = _trailing_kink_count(seg_reports)

    current_kink_run = trailing_kink_run + 1 if use_kink else 0
    if use_kink and current_kink_run >= max_seq_kink:
        warning_msg = (
            f"[{tag0}] Consecutive kink segments were detected. Something seems wrong. "
            "Please check the initial structure and the generated intermediate structures. "
            "Alternatively, try switching the mep-mode. If that still fails, try including intermediate structures in the inputs."
        )
        return _terminate_with_maxdepth(reason_msg=warning_msg)

    parts.append((step_imgs, step_E))
    seg_reports.append(seg_report)

    if right_changed:
        subR = _build_multistep_path(
            right_end, gB, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg, bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
            out_dir, ref_pdb_path, prepared_input, depth + 1, seg_counter, branch_tag=f"{branch_tag}R",
            pair_index=pair_index,
            kink_seq_count=current_kink_run,
        )
        _tag_images(subR.images, pair_index=pair_index)
        parts.append((subR.images, subR.energies))
        seg_reports.extend(subR.segments)

    bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 5))
    gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

    def _segment_builder(tail_g, head_g, _tag: str) -> CombinedPath:
        sub = _build_multistep_path(
            tail_g, head_g,
            shared_calc,
            geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg,
            bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
            out_dir=out_dir,
            ref_pdb_path=ref_pdb_path,
            prepared_input=prepared_input,
            depth=depth + 1,
            seg_counter=seg_counter,
            branch_tag=f"{branch_tag}B",
            pair_index=pair_index,
            kink_seq_count=_trailing_kink_count(seg_reports),
        )
        _tag_images(sub.images, pair_index=pair_index)
        return sub

    stitched_imgs, stitched_E = _stitch_paths(
        parts,
        stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-4)),
        bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-4)),
        shared_calc=shared_calc,
        gs_cfg=gs_bridge_cfg,
        opt_cfg=opt_cfg,
        out_dir=out_dir,
        tag=tag0,
        ref_pdb_path=ref_pdb_path,
        bond_cfg=bond_cfg,
        segment_builder=_segment_builder,
        segments_out=seg_reports,
        bridge_pair_index=pair_index,
        mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg,
        max_nodes=bridge_max_nodes,
        dmf_cfg=dmf_cfg,
        prepared_inputs=prepared_inputs,
    )

    _tag_images(stitched_imgs, pair_index=pair_index)

    return CombinedPath(images=stitched_imgs, energies=stitched_E, segments=seg_reports)


# -----------------------------------------------
# Full‑system merge helpers (Biopython)
# -----------------------------------------------

def _atom_key_from_res_atom(res: PDB.Residue.Residue, atom: PDB.Atom.Atom) -> Tuple[str, str, str, str, str]:
    """
    Build a key for atom identity:
    (RESNAME, RESSEQ, ICODE, CHAIN, ATOMNAME) — uppercase where applicable.
    - RESSEQ is numeric (without insertion code).
    - ICODE is '' when blank (or ' ' in PDB).
    """
    resname = (res.get_resname() or "").strip().upper()
    het, resseq, icode = res.id
    icode_txt = "" if (icode == " " or icode is None) else str(icode).strip().upper()
    resseq_txt = str(int(resseq))
    chain_id = (res.get_parent().id or "").strip().upper()
    atname = atom.get_name().strip().upper()
    return (resname, resseq_txt, icode_txt, chain_id, atname)


def _structure_to_arrays(struct: PDB.Structure.Structure) -> Tuple[np.ndarray, List[PDB.Atom.Atom], List[Tuple[str, str, str, str, str]], Dict[Tuple[str,str,str,str,str], int]]:
    """
    Extract: coordinates (Å), atom list, key list, and key→index map from a Biopython Structure.
    Keys are (RESNAME, RESSEQ, ICODE, CHAIN, ATOMNAME).
    """
    atoms: List[PDB.Atom.Atom] = [a for a in struct.get_atoms()]
    coords = np.array([a.get_coord() for a in atoms], dtype=float)
    keys: List[Tuple[str, str, str, str, str]] = []
    key2idx: Dict[Tuple[str,str,str,str,str], int] = {}
    for idx, a in enumerate(atoms):
        res = a.get_parent()
        k = _atom_key_from_res_atom(res, a)
        keys.append(k)
        if k not in key2idx:
            key2idx[k] = idx
    return coords, atoms, keys, key2idx


def _load_structures_and_chain_align(ref_paths: Sequence[Path]) -> Tuple[List[PDB.Structure.Structure], List[np.ndarray], List[List[PDB.Atom.Atom]], List[Dict[Tuple[str,str,str,str,str], int]]]:
    """
    Load all full templates and rigidly chain‑align them into the coordinate frame of the first template.
    """
    parser = PDBParser(QUIET=True)
    structs: List[PDB.Structure.Structure] = [parser.get_structure(f"ref{i:02d}", str(p)) for i, p in enumerate(ref_paths)]
    N_expected = None
    coords_list: List[np.ndarray] = []
    atoms_list: List[List[PDB.Atom.Atom]] = []
    keymaps: List[Dict[Tuple[str,str,str,str,str], int]] = []

    for s in structs:
        coords, atoms, _keys, key2idx = _structure_to_arrays(s)
        if N_expected is None:
            N_expected = coords.shape[0]
        else:
            if coords.shape[0] != N_expected:
                raise click.BadParameter(f"[merge] Atom count mismatch among --ref-full-pdb templates: {N_expected} vs {coords.shape[0]}")
        coords_list.append(coords)
        atoms_list.append(atoms)
        keymaps.append(key2idx)

    aligned_coords: List[np.ndarray] = []
    aligned_coords.append(coords_list[0].copy())
    for j in range(1, len(coords_list)):
        P = aligned_coords[j-1]
        Q = coords_list[j]
        R, t = kabsch_R_t(P, Q)
        Qa = (Q @ R) + t
        aligned_coords.append(Qa)

    return structs, aligned_coords, atoms_list, keymaps


def _pocket_keys_from_pdb(pocket_pdb: Path) -> List[Tuple[str, str, str, str, str]]:
    """
    Return atom identity keys for a pocket PDB file.
    """
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("pocket", str(pocket_pdb))
    keys: List[Tuple[str, str, str, str, str]] = []
    for a in st.get_atoms():
        res = a.get_parent()
        k = _atom_key_from_res_atom(res, a)
        keys.append(k)
    return keys


def _write_model_block(structure: PDB.Structure.Structure,
                       remark_lines: List[str]) -> str:
    """
    Render a single MODEL block (without 'MODEL/ENDMDL') with provided REMARK lines.
    """
    io = PDBIO()
    io.set_structure(structure)
    from io import StringIO
    buf = StringIO()
    io.save(buf)
    body = "\n".join([ln for ln in buf.getvalue().splitlines() if ln.strip() != "END"])
    rem = ""
    for line in remark_lines:
        rem += f"REMARK   1 {line}\n"
    return rem + body + ("\n" if not body.endswith("\n") else "")


def _chunk_remark_indices(indices: List[int], width: int = 60) -> List[str]:
    """
    Wrap pocket atom indices into REMARK lines with limited width.
    """
    s = ",".join(map(str, indices))
    out: List[str] = []
    cur = ""
    for tok in s.split(","):
        add = (tok if not cur else "," + tok)
        if len(cur) + len(add) > width:
            out.append(f"POCKET_ATOM_INDICES {cur}")
            cur = tok
        else:
            cur += tok if not cur else "," + tok
    if cur:
        out.append(f"POCKET_ATOM_INDICES {cur}")
    return out


def _merge_pair_to_full(pair_images: List[Any],
                        pocket_ref_pdb: Path,
                        structA: PDB.Structure.Structure,
                        structB: PDB.Structure.Structure,
                        coordsA_aligned: np.ndarray,
                        coordsB_aligned: np.ndarray,
                        keymapA: Dict[Tuple[str,str,str,str,str], int],
                        keymapB: Dict[Tuple[str,str,str,str,str], int],
                        out_path: Optional[Path],
                        drop_first: bool = False,
                        seg_indices_for_frames: Optional[List[int]] = None,
                        seg_report_lookup: Optional[Dict[int, SegmentReport]] = None,
                        include_pocket_indices_for_first_model: bool = False) -> Tuple[List[str], List[int]]:
    """
    Merge a pocket‑only trajectory for a *pair* into the corresponding full templates (A,B),
    generating MODEL blocks and (optionally) writing a PDB. Returns (blocks, 1‑based active indices).
    """
    pocket_keys = _pocket_keys_from_pdb(pocket_ref_pdb)

    match_tpl_idx: List[int] = []
    for k in pocket_keys:
        ia = keymapA.get(k, None)
        ib = keymapB.get(k, None)
        if ia is None or ib is None:
            match_tpl_idx.append(-1)
        else:
            match_tpl_idx.append(int(ia))

    active_full_idx = sorted({i for i in match_tpl_idx if i >= 0})
    active_one_based = [i+1 for i in active_full_idx]

    Nfull = coordsA_aligned.shape[0]
    if Nfull != coordsB_aligned.shape[0]:
        raise click.BadParameter("[merge] Template A/B atom count mismatch.")

    atomsA: List[PDB.Atom.Atom] = [a for a in structA.get_atoms()]

    start_k = 1 if drop_first and len(pair_images) > 0 else 0

    model_blocks: List[str] = []

    M = len(pair_images)
    if M == 0:
        return model_blocks, active_one_based

    seg_idx_seq: List[int] = []
    if seg_indices_for_frames is not None and len(seg_indices_for_frames) == M:
        seg_idx_seq = seg_indices_for_frames[start_k:]
    else:
        seg_idx_seq = [0] * (M - start_k)

    for kk, k in enumerate(range(start_k, M)):
        tfrac = 0.0 if M == 1 else (k / (M - 1.0))

        C = (1.0 - tfrac) * coordsA_aligned + tfrac * coordsB_aligned

        idx_sel = [j for j in range(len(match_tpl_idx)) if match_tpl_idx[j] >= 0]
        if len(idx_sel) >= 3:
            Y = np.array([C[match_tpl_idx[j]] for j in idx_sel], dtype=float)
            P_bohr = np.array(pair_images[k].coords3d, dtype=float)
            P = P_bohr * BOHR2ANG
            P_sel = np.array([P[j] for j in idx_sel], dtype=float)
            R, t = kabsch_R_t(Y, P_sel)
            Paligned = (P @ R) + t
            for jj, pidx in enumerate(idx_sel):
                full_i = match_tpl_idx[pidx]
                if 0 <= full_i < Nfull:
                    C[full_i] = Paligned[pidx]
        else:
            P = np.array(pair_images[k].coords3d, dtype=float) * BOHR2ANG
            for j, full_i in enumerate(match_tpl_idx):
                if full_i >= 0 and 0 <= full_i < Nfull and j < P.shape[0]:
                    C[full_i] = P[j]

        for i, a in enumerate(atomsA):
            a.set_coord(C[i])
            a.set_bfactor(100.0 if i in active_full_idx else 0.0)

        remark_lines: List[str] = []
        remark_lines.append(f"PAIR_MERGE FRAC {tfrac:.6f}")

        if include_pocket_indices_for_first_model and kk == 0:
            remark_lines.extend(_chunk_remark_indices([i for i in active_one_based], width=60))

        seg_idx = seg_idx_seq[kk] if kk < len(seg_idx_seq) else 0
        if seg_idx and seg_report_lookup is not None:
            rep = seg_report_lookup.get(seg_idx, None)
            if rep is not None:
                remark_lines.append(
                    f"MEP_SEG_INDEX {int(seg_idx):02d} TAG {rep.tag} KIND {rep.kind} "
                    f"DELTAE_BARRIER_KCAL {rep.barrier_kcal:.6f} DELTAE_KCAL {rep.delta_kcal:.6f}"
                )
                if rep.kind != "bridge" and rep.summary and rep.summary.strip() and rep.summary.strip() != "(no covalent changes detected)":
                    for ln in rep.summary.strip().splitlines():
                        remark_lines.append(f"SEG_BONDS {ln.strip()}")
            else:
                remark_lines.append(f"MEP_SEG_INDEX {int(seg_idx):02d}")

        model_blocks.append(_write_model_block(structA, remark_lines))

    if out_path is not None:
        with open(out_path, "w") as f:
            for m, blk in enumerate(model_blocks, start=1):
                f.write(f"MODEL     {m}\n")
                f.write(blk)
                f.write("ENDMDL\n")
            f.write("END\n")
        click.echo(f"[merge] Wrote pair-merged PDB → '{out_path}'")

    return model_blocks, active_one_based


def _merge_final_and_write(final_images: List[Any],
                           pocket_inputs: Sequence[Path],
                           ref_pdbs: Sequence[Path],
                           segments: List[SegmentReport],
                           out_dir: Path,
                           pocket_ref_pdbs: Optional[Sequence[Path]] = None) -> None:
    """
    Merge the entire pocket MEP into full templates (for all pairs) and write outputs.
    """
    if len(ref_pdbs) != len(pocket_inputs):
        raise click.BadParameter("--ref-full-pdb must match the number of --input after preprocessing (caller should replicate the first ref for all pairs when --align True).")

    if pocket_ref_pdbs is None:
        pocket_ref_pdbs = pocket_inputs
    if len(pocket_ref_pdbs) != len(pocket_inputs):
        raise click.BadParameter("--ref-pdb must match the number of --input after preprocessing.")

    structs, aligned_coords, _atoms_list, keymaps = _load_structures_and_chain_align(ref_pdbs)

    seg_lookup: Dict[int, SegmentReport] = {int(s.seg_index): s for s in segments if int(s.seg_index) > 0}

    n_pairs = len(pocket_inputs) - 1
    groups: List[Tuple[int, List[Any]]] = []
    cur_idx = None
    cur_list: List[Any] = []
    for im in final_images:
        pi = getattr(im, "pair_index", None)
        if pi is None:
            pi = 0
        if cur_idx is None:
            cur_idx = int(pi)
            cur_list = [im]
        elif int(pi) == int(cur_idx):
            cur_list.append(im)
        else:
            groups.append((int(cur_idx), cur_list))
            cur_idx = int(pi)
            cur_list = [im]
    if cur_list:
        groups.append((int(cur_idx), cur_list))

    for (pi, _) in groups:
        if not (0 <= pi < n_pairs):
            raise click.BadParameter(f"[merge] Illegal pair_index {pi} (n_pairs={n_pairs}).")

    final_blocks: List[str] = []
    wrote_indices = False

    for gi, (pi, imgs) in enumerate(groups):
        pocket_ref = Path(pocket_ref_pdbs[pi])
        structA = structs[pi]
        structB = structs[pi+1]
        coordsA = aligned_coords[pi]
        coordsB = aligned_coords[pi+1]
        keymapA = keymaps[pi]
        keymapB = keymaps[pi+1]
        seg_indices_for_frames = [int(getattr(im, "mep_seg_index", 0) or 0) for im in imgs]

        blocks, active_one_based = _merge_pair_to_full(
            pair_images=imgs,
            pocket_ref_pdb=pocket_ref,
            structA=structA,
            structB=structB,
            coordsA_aligned=coordsA,
            coordsB_aligned=coordsB,
            keymapA=keymapA,
            keymapB=keymapB,
            out_path=None,
            drop_first=(gi > 0),
            seg_indices_for_frames=seg_indices_for_frames,
            seg_report_lookup=seg_lookup,
            include_pocket_indices_for_first_model=(not wrote_indices)
        )
        if blocks:
            wrote_indices = True
            final_blocks.extend(blocks)

    final_path = out_dir / "mep_w_ref.pdb"
    with open(final_path, "w") as f:
        for m, blk in enumerate(final_blocks, start=1):
            f.write(f"MODEL     {m}\n")
            f.write(blk)
            f.write("ENDMDL\n")
        f.write("END\n")
    click.echo(f"[merge] Wrote concatenated full-system trajectory → '{final_path}'")

    # Per‑segment merged MEPs (bond‑change segments only) + HEI merged only for bond‑change segments
    for s in segments:
        seg_idx = int(s.seg_index)
        seg_frames: List[Any] = [im for im in final_images if int(getattr(im, "mep_seg_index", 0) or 0) == seg_idx]
        if not seg_frames:
            continue

        # Determine pair index for this segment (assume consistent within the segment)
        pi_vals = sorted({int(getattr(im, "pair_index", 0)) for im in seg_frames})
        pi = pi_vals[0]
        pocket_ref = Path(pocket_ref_pdbs[pi])
        structA = structs[pi]
        structB = structs[pi+1]
        coordsA = aligned_coords[pi]
        coordsB = aligned_coords[pi+1]
        keymapA = keymaps[pi]
        keymapB = keymaps[pi+1]

        # Per‑segment merged MEP only when covalent changes are present
        if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
            seg_indices_for_frames = [seg_idx] * len(seg_frames)
            blocks, _ = _merge_pair_to_full(
                pair_images=seg_frames,
                pocket_ref_pdb=pocket_ref,
                structA=structA,
                structB=structB,
                coordsA_aligned=coordsA,
                coordsB_aligned=coordsB,
                keymapA=keymapA,
                keymapB=keymapB,
                out_path=None,
                drop_first=False,
                seg_indices_for_frames=seg_indices_for_frames,
                seg_report_lookup=seg_lookup,
                include_pocket_indices_for_first_model=True,
            )
            out_seg = out_dir / f"mep_w_ref_seg_{seg_idx:02d}.pdb"
            with open(out_seg, "w") as f:
                for m, blk in enumerate(blocks, start=1):
                    f.write(f"MODEL     {m}\n")
                    f.write(blk)
                    f.write("ENDMDL\n")
                f.write("END\n")
            click.echo(f"[merge] Wrote per-segment merged trajectory → '{out_seg}'")

        # Per‑segment HEI merged to reference (only for bond‑change segments)
        if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
            try:
                energies_au = [float(getattr(im, "energy")) for im in seg_frames]
                imax = int(np.argmax(np.array(energies_au, dtype=float)))
                hei_frame = seg_frames[imax]
                blocks_hei, _ = _merge_pair_to_full(
                    pair_images=[hei_frame],
                    pocket_ref_pdb=pocket_ref,
                    structA=structA,
                    structB=structB,
                    coordsA_aligned=coordsA,
                    coordsB_aligned=coordsB,
                    keymapA=keymapA,
                    keymapB=keymapB,
                    out_path=None,
                    drop_first=False,
                    seg_indices_for_frames=[seg_idx],
                    seg_report_lookup=seg_lookup,
                    include_pocket_indices_for_first_model=True,
                )
                out_hei = out_dir / f"hei_w_ref_seg_{seg_idx:02d}.pdb"
                with open(out_hei, "w") as f:
                    for m, blk in enumerate(blocks_hei, start=1):
                        f.write(f"MODEL     {m}\n")
                        f.write(blk)
                        f.write("ENDMDL\n")
                    f.write("END\n")
                click.echo(f"[merge] Wrote merged HEI for segment → '{out_hei}'")
            except Exception as e:
                click.echo(f"[merge] WARNING: Failed to write merged HEI for segment {seg_idx:02d}: {e}", err=True)


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Multistep MEP search via recursive GSM/DMF segmentation.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        # Allow a single '-i' followed by multiple paths (as extra args)
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    required=True,
    help=("Two or more structures in reaction order. "
          "Either repeat '-i' (e.g., '-i A -i B -i C') or use a single '-i' "
          "followed by multiple space-separated paths (e.g., '-i A B C').")
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
)
@click.option(
    "--refine-mode",
    type=click.Choice(["peak", "minima"], case_sensitive=False),
    default=None,
    show_default=True,
    help=(
        "Refinement seed selection around the highest-energy image: "
        "'peak' uses HEI±1, 'minima' uses the nearest local minima in each direction. "
        "Defaults to peak for gsm and minima for dmf when omitted."
    ),
)
@click.option(
    "-q",
    "--charge",
    type=int,
    default=None,
    show_default=False,
    help=(
        "Total charge. Required for non-.gjf inputs unless --ligand-charge derives it "
        "from PDB inputs."
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
        "when -q is omitted (PDB inputs only)."
    ),
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1) for the ML region (defaults from a .gjf template when available, otherwise 1).",
)
@click.option("--freeze-links", "freeze_links_flag", type=click.BOOL, default=True, show_default=True,
              help="For PDB input, freeze parent atoms of link hydrogens.")
@click.option("--max-nodes", type=int, default=10, show_default=True,
              help=("Number of internal nodes (string has max_nodes+2 images including endpoints). "
                    "Used for *segment* GSM unless overridden by YAML search.max_nodes_segment."))
@click.option("--max-cycles", type=int, default=300, show_default=True, help="Maximum GSM optimization cycles.")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Enable climbing image for standard GSM segments (bridge segments always disable climbing).")
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Single-structure optimizer: light (=LBFGS) or heavy (=RFO).",
)
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump GSM/single-optimization trajectories during the run.")
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option("--out-dir", "out_dir", type=str, default="./result_path_search/", show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Convergence preset for GSM and single optimizations "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
    ),
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt, sopt, bond, search, dmf)."
)
@click.option(
    "--preopt",
    "preopt",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If False, skip initial single-structure optimizations of inputs."
)
@click.option(
    "--align",
    "align",
    type=click.BOOL,
    default=True,
    show_default=True,
    help=("After preoptimization, align all inputs to the *first* input and match freeze_atoms "
          "using the align_freeze_atoms API. When --align is True and --ref-full-pdb is provided, "
          "the first reference PDB will be used for all pairs in the final merge.")
)
@click.option(
    "--ref-full-pdb",
    "ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Full-size template PDBs in the same reaction order as --input. "
          "With --align True, only the *first* provided reference PDB is used for all pairs "
          "in the final merge (you may pass just one).")
)
@click.option(
    "--ref-pdb",
    "pocket_ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Pocket reference PDBs used only for the final full-system merge. "
          "Useful when --input uses XYZ/GJF intermediates but PDB snapshots exist for merging. "
          "Must match the number and order of --input.")
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    mep_mode: str,
    refine_mode: Optional[str],
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    align: bool,
    ref_pdb_paths: Optional[Sequence[Path]],
    pocket_ref_pdb_paths: Optional[Sequence[Path]],
) -> None:
    set_convert_file_enabled(convert_files)
    prepared_inputs: List[PreparedInputStructure] = []
    global _PRIMARY_GJF_TEMPLATE
    _PRIMARY_GJF_TEMPLATE = None
    command_str = " ".join(sys.argv)

    # Robustly accept both styles for -i/--input, --ref-full-pdb, and --ref-pdb
    def _collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
        vals: List[str] = []
        i = 0
        while i < len(argv):
            tok = argv[i]
            if tok in names:
                j = i + 1
                while j < len(argv) and not argv[j].startswith("-"):
                    vals.append(argv[j])
                    j += 1
                i = j
            else:
                i += 1
        return vals

    argv_all = sys.argv[1:]
    i_vals = _collect_option_values(argv_all, ("-i", "--input"))
    if i_vals:
        i_parsed: List[Path] = []
        for tok in i_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Input path '{tok}' not found or is a directory. "
                    f"When using '-i', list only existing file paths (multiple paths may follow a single '-i')."
                )
            i_parsed.append(p)
        input_paths = tuple(i_parsed)

    ref_vals = _collect_option_values(argv_all, ("--ref-full-pdb",))
    if ref_vals:
        ref_parsed: List[Path] = []
        for tok in ref_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Reference PDB path '{tok}' not found or is a directory. "
                    f"When using '--ref-full-pdb', multiple files may follow a single option."
                )
            ref_parsed.append(p)
        ref_pdb_paths = tuple(ref_parsed)

    pocket_ref_vals = _collect_option_values(argv_all, ("--ref-pdb",))
    if pocket_ref_vals:
        pocket_ref_parsed: List[Path] = []
        for tok in pocket_ref_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Pocket reference PDB path '{tok}' not found or is a directory. "
                    f"When using '--ref-pdb', multiple files may follow a single option."
                )
            pocket_ref_parsed.append(p)
        pocket_ref_pdb_paths = tuple(pocket_ref_parsed)

    time_start = time.perf_counter()
    freeze_atoms_for_log: List[int] = []
    try:
        # --------------------------
        # 0) Input validation (multi‑structure)
        # --------------------------
        if len(input_paths) < 2:
            raise click.BadParameter(
                "Provide at least two structures for --input in reaction order (reactant [intermediates ...] product)."
            )

        mep_mode_kind = mep_mode.strip().lower()
        refine_mode_kind = refine_mode.strip().lower() if refine_mode else None

        do_merge = bool(ref_pdb_paths) and len(ref_pdb_paths) > 0
        if do_merge:
            if align:
                pass
            else:
                if len(ref_pdb_paths) != len(input_paths):
                    raise click.BadParameter("--ref-full-pdb must be given for each --input (same count and order). "
                                             "Alternatively, use --align to allow using only the first reference PDB for all pairs.")
            if pocket_ref_pdb_paths and len(pocket_ref_pdb_paths) != len(input_paths):
                raise click.BadParameter("--ref-pdb must be given for each --input (same count and order).")

        p_list = [Path(p) for p in input_paths]
        prepared_inputs = [prepare_input_structure(p) for p in p_list]
        any_non_gjf = any(not prep.is_gjf for prep in prepared_inputs)
        if _PRIMARY_GJF_TEMPLATE is None:
            _PRIMARY_GJF_TEMPLATE = next((prep.gjf_template for prep in prepared_inputs if prep.gjf_template), None)

        # --------------------------
        # 1) Resolve settings (defaults ← CLI ← YAML)
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        dmf_cfg  = dict(DMF_KW)
        gs_cfg   = dict(GS_KW)
        opt_cfg  = dict(STOPT_KW)
        lbfgs_cfg = dict(_LBFGS_KW)
        rfo_cfg   = dict(_RFO_KW)
        bond_cfg  = dict(BOND_KW)
        search_cfg = dict(SEARCH_KW)
        search_cfg["refine_mode"] = refine_mode_kind

        resolved_charge = charge
        resolved_spin = spin
        for prepared in prepared_inputs:
            resolved_charge, resolved_spin = fill_charge_spin_from_gjf(
                resolved_charge, resolved_spin, prepared.gjf_template
            )
        if resolved_charge is None and ligand_charge is not None:
            resolved_charge = _derive_charge_from_ligand_charge(
                prepared_inputs[0], ligand_charge, prefix="[path-search]"
            )
        if resolved_charge is None:
            if any_non_gjf:
                for prepared in prepared_inputs:
                    prepared.cleanup()
                raise click.ClickException(
                    "-q/--charge is required unless all inputs are .gjf templates with charge metadata."
                )
            resolved_charge = 0
        if resolved_spin is None:
            resolved_spin = 1
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)
        calc_cfg["workers"] = int(workers)
        calc_cfg["workers_per_node"] = int(workers_per_node)

        gs_cfg["max_nodes"] = int(max_nodes)
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)

        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["thresh"] = str(thresh)
            rfo_cfg["thresh"] = str(thresh)

        lbfgs_cfg["dump"] = bool(dump)
        rfo_cfg["dump"]   = bool(dump)
        lbfgs_cfg["out_dir"] = out_dir
        rfo_cfg["out_dir"]   = out_dir

        # YAML overrides (highest precedence)
        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("sopt", "lbfgs"), ("opt", "lbfgs"), ("lbfgs",))),
                (rfo_cfg, (("sopt", "rfo"), ("opt", "rfo"), ("rfo",))),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
            ],
        )

        refine_mode_kind = search_cfg.get("refine_mode")
        if refine_mode_kind is None:
            refine_mode_kind = "peak" if mep_mode_kind == "gsm" else "minima"
        else:
            refine_mode_kind = str(refine_mode_kind).strip().lower()
            if refine_mode_kind not in {"peak", "minima"}:
                raise click.BadParameter(f"Unknown --refine-mode '{refine_mode_kind}'.")
        search_cfg["refine_mode"] = refine_mode_kind

        opt_kind = opt_mode.strip().lower()
        if opt_kind == "light":
            sopt_kind = "lbfgs"
            sopt_cfg = lbfgs_cfg
        elif opt_kind == "heavy":
            sopt_kind = "rfo"
            sopt_cfg = rfo_cfg
        else:
            raise click.BadParameter(f"Unknown --opt-mode '{opt_mode}'.")

        if "max_nodes_segment" not in yaml_cfg.get("search", {}):
            search_cfg["max_nodes_segment"] = int(max_nodes)

        out_dir_path = Path(out_dir).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_freeze_atoms_for_echo(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_opt  = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)
        opt_yaml = yaml_cfg.get("opt")
        if isinstance(opt_yaml, dict) and "out_dir" in opt_yaml:
            yaml_out_dir = opt_yaml.get("out_dir")
            if yaml_out_dir is not None:
                try:
                    yaml_out_dir_resolved = str(Path(str(yaml_out_dir)).resolve())
                except Exception:
                    yaml_out_dir_resolved = str(yaml_out_dir)
                if yaml_out_dir_resolved != str(out_dir_path):
                    echo_opt["out_dir_yaml"] = str(yaml_out_dir)
                    echo_opt["out_dir_note"] = (
                        "CLI --out-dir takes precedence; YAML opt.out_dir is ignored."
                    )

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs",   echo_gs))
        click.echo(pretty_block("opt",  echo_opt))
        if mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg))
        echo_sopt = dict(sopt_cfg)
        echo_sopt["out_dir"] = str(out_dir_path)
        echo_sopt["out_dir_per_tag"] = f"{out_dir_path}/<tag>_{sopt_kind}_opt"
        click.echo(pretty_block("sopt."+sopt_kind, echo_sopt))
        click.echo(pretty_block("bond", bond_cfg))
        click.echo(pretty_block("search", search_cfg))
        click.echo(
            pretty_block(
                "run_flags",
                {"preopt": bool(preopt), "align": bool(align), "mep_mode": mep_mode_kind},
            )
        )

        # --------------------------
        # 2) Prepare inputs
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        geoms = _load_structures(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"]),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )
        if geoms:
            freeze_effective: Dict[str, List[int]] = {}
            for prepared, g in zip(prepared_inputs, geoms):
                try:
                    freeze_list = list(getattr(g, "freeze_atoms", []))
                except Exception:
                    freeze_list = list(geom_cfg.get("freeze_atoms", []))
                freeze_effective[prepared.source_path.name] = freeze_list
            click.echo(pretty_block("freeze_atoms (effective)", freeze_effective))

        main_prepared = prepared_inputs[0] if prepared_inputs else None

        shared_calc = uma_pysis(**calc_cfg)
        for g in geoms:
            g.set_calculator(shared_calc)

        # If any input is PDB, treat as "PDB input" for final output handling.
        ref_pdb_for_segments: Optional[Path] = None
        for p in p_list:
            if p.suffix.lower() == ".pdb":
                ref_pdb_for_segments = p.resolve()
                break

        if preopt:
            new_geoms: List[Any] = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                g_opt = _optimize_single(
                    g,
                    shared_calc,
                    sopt_kind,
                    sopt_cfg,
                    out_dir_path,
                    tag=tag,
                    prepared_input=prepared_inputs[i] if i < len(prepared_inputs) else main_prepared,
                )
                new_geoms.append(g_opt)
            geoms = new_geoms
        else:
            click.echo("[init] Skipping endpoint preoptimization as requested by --preopt False.")

        # Align all inputs to the first structure, guided by freeze constraints, when requested
        align_thresh = str(opt_cfg.get("thresh", "gau"))
        if align:
            try:
                click.echo("\n=== Aligning all inputs to the first structure (freeze-guided scan + relaxation) ===\n")
                _ = align_and_refine_sequence_inplace(
                    geoms,
                    thresh=align_thresh,
                    shared_calc=shared_calc,
                    out_dir=out_dir_path / "align_refine",
                    verbose=True,
                )
                click.echo("[align] Completed input alignment.")
            except Exception as e:
                click.echo(f"[align] WARNING: Alignment failed; continuing without alignment: {e}", err=True)
        else:
            click.echo("[align] Skipping input alignment as requested by --align False.")

        # --------------------------
        # 3) Run recursive search for each adjacent pair and stitch
        # --------------------------
        click.echo("\n=== Multistep MEP search (multi-structure) started ===\n")
        seg_counter = [0]

        bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 10))
        gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

        combined_imgs: List[Any] = []
        combined_Es: List[float] = []
        seg_reports_all: List[SegmentReport] = []

        def _segment_builder_for_pairs(tail_g, head_g, _tag: str) -> CombinedPath:
            sub = _build_multistep_path(
                tail_g, head_g,
                shared_calc,
                geom_cfg, gs_cfg, opt_cfg,
                sopt_kind, sopt_cfg,
                bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                prepared_input=main_prepared,
                depth=0,
                seg_counter=seg_counter,
                branch_tag="B",
                pair_index=None,
                kink_seq_count=_trailing_kink_count(seg_reports_all),
            )
            return sub

        for i in range(len(geoms) - 1):
            gA, gB = geoms[i], geoms[i + 1]
            pair_tag = f"pair_{i:02d}"
            click.echo(f"\n--- Processing pair {i:02d}: image {i} → {i+1} ---")
            pair_path = _build_multistep_path(
                gA, gB,
                shared_calc,
                geom_cfg, gs_cfg, opt_cfg,
                sopt_kind, sopt_cfg,
                bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                prepared_input=main_prepared,
                depth=0,
                seg_counter=seg_counter,
                branch_tag=pair_tag,
                pair_index=i,
                kink_seq_count=_trailing_kink_count(seg_reports_all),
            )

            if i == 0:
                combined_imgs = list(pair_path.images)
                combined_Es = list(pair_path.energies)
                seg_reports_all.extend(pair_path.segments)
            else:
                parts = [(combined_imgs, combined_Es), (pair_path.images, pair_path.energies)]
                combined_imgs, combined_Es = _stitch_paths(
                    parts=parts,
                    stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-4)),
                    bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-4)),
                    shared_calc=shared_calc,
                    gs_cfg=gs_bridge_cfg,
                    opt_cfg=opt_cfg,
                    out_dir=out_dir_path,
                    tag=pair_tag,
                    ref_pdb_path=ref_pdb_for_segments,
                    bond_cfg=bond_cfg,
                    segment_builder=_segment_builder_for_pairs,
                    segments_out=seg_reports_all,
                    bridge_pair_index=i,
                    mep_mode_kind=mep_mode_kind,
                    calc_cfg=calc_cfg,
                    max_nodes=bridge_max_nodes,
                    dmf_cfg=dmf_cfg,
                    prepared_inputs=prepared_inputs,
                )
                seg_reports_all.extend(pair_path.segments)

        click.echo("\n=== Multistep MEP search (multi-structure) finished ===\n")

        combined_all = CombinedPath(images=combined_imgs, energies=combined_Es, segments=seg_reports_all)

        # --------------------------
        # 4) Outputs
        # --------------------------
        for idx, srep in enumerate(combined_all.segments, 1):
            srep.seg_index = idx
        tag_to_index = {s.tag: int(s.seg_index) for s in combined_all.segments}
        for im in combined_all.images:
            tag = getattr(im, "mep_seg_tag", None)
            if tag and tag in tag_to_index:
                try:
                    setattr(im, "mep_seg_index", int(tag_to_index[tag]))
                except Exception:
                    pass

        # Final MEP output rule:
        # - Always write 'mep.trj' (XYZ) for intermediate handoff.
        # - If reference PDBs are available, also emit 'mep.pdb' (and GJF when applicable).
        main_prepared = prepared_inputs[0]
        needs_pdb = ref_pdb_for_segments is not None
        needs_gjf = main_prepared.is_gjf

        final_trj = out_dir_path / "mep.trj"
        _write_xyz_trj_with_energy(combined_all.images, combined_all.energies, final_trj)
        click.echo(f"[write] Wrote '{final_trj}'.")
        try:
            run_trj2fig(final_trj, [out_dir_path / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            _close_matplotlib_figures()
            click.echo(f"[plot] Saved energy plot → '{out_dir_path / 'mep_plot.png'}'")
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot final energy: {e}", err=True)

        if needs_pdb or needs_gjf:
            try:
                convert_xyz_like_outputs(
                    final_trj,
                    main_prepared,
                    ref_pdb_path=ref_pdb_for_segments,
                    out_pdb_path=out_dir_path / "mep.pdb" if needs_pdb else None,
                    out_gjf_path=out_dir_path / "mep.gjf" if needs_gjf else None,
                )
                click.echo("[convert] Wrote final MEP outputs.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final MEP outputs: {e}", err=True)

        # Pocket‑only per‑segment trajectories & HEIs (bond‑change segments only)
        try:
            frame_seg_indices: List[int] = [int(getattr(im, "mep_seg_index", 0) or 0) for im in combined_all.images]
            seg_to_frames: Dict[int, List[int]] = {}
            for ii, sidx in enumerate(frame_seg_indices):
                if sidx <= 0:
                    continue
                seg_to_frames.setdefault(int(sidx), []).append(ii)

            for s in combined_all.segments:
                seg_idx = int(s.seg_index)
                idxs = seg_to_frames.get(seg_idx, [])
                if not idxs:
                    continue

                if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    seg_imgs = [combined_all.images[j] for j in idxs]
                    seg_Es = [combined_all.energies[j] for j in idxs]
                    seg_trj = out_dir_path / f"mep_seg_{seg_idx:02d}.trj"
                    _write_xyz_trj_with_energy(seg_imgs, seg_Es, seg_trj)
                    click.echo(f"[write] Wrote per-segment pocket trajectory → '{seg_trj}'")
                    if needs_pdb or needs_gjf:
                        try:
                            convert_xyz_like_outputs(
                                seg_trj,
                                main_prepared,
                                ref_pdb_path=ref_pdb_for_segments,
                                out_pdb_path=out_dir_path / f"mep_seg_{seg_idx:02d}.pdb" if needs_pdb else None,
                                out_gjf_path=out_dir_path / f"mep_seg_{seg_idx:02d}.gjf" if needs_gjf else None,
                            )
                        except Exception as e:
                            click.echo(
                                f"[convert] WARNING: Failed to convert per-segment trajectory {seg_idx:02d}: {e}",
                                err=True,
                            )

                if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    energies_seg = [combined_all.energies[j] for j in idxs]
                    imax_rel = int(np.argmax(np.array(energies_seg, dtype=float)))
                    imax_abs = idxs[imax_rel]
                    hei_img = combined_all.images[imax_abs]
                    hei_E = [combined_all.energies[imax_abs]]
                    hei_trj = out_dir_path / f"hei_seg_{seg_idx:02d}.xyz"
                    _write_xyz_trj_with_energy([hei_img], hei_E, hei_trj)
                    click.echo(f"[write] Wrote segment HEI (pocket) → '{hei_trj}'")
                    if needs_pdb or needs_gjf:
                        try:
                            convert_xyz_like_outputs(
                                hei_trj,
                                main_prepared,
                                ref_pdb_path=ref_pdb_for_segments,
                                out_pdb_path=out_dir_path / f"hei_seg_{seg_idx:02d}.pdb" if needs_pdb else None,
                                out_gjf_path=out_dir_path / f"hei_seg_{seg_idx:02d}.gjf" if needs_gjf else None,
                            )
                        except Exception as e:
                            click.echo(
                                f"[convert] WARNING: Failed to convert HEI for segment {seg_idx:02d}: {e}",
                                err=True,
                            )
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to emit per-segment pocket outputs: {e}", err=True)

        if do_merge:
            click.echo("\n=== Full-system merge (pocket → templates) started ===\n")
            # With --align True, use only the first reference PDB for all pairs (replicate it).
            if align:
                if not ref_pdb_paths or len(ref_pdb_paths) < 1:
                    raise click.BadParameter("--ref-full-pdb must provide at least one file when performing final merge with --align True.")
                first_ref = Path(ref_pdb_paths[0])
                ref_list_for_merge = [first_ref for _ in input_paths]
            else:
                ref_list_for_merge = [Path(p) for p in ref_pdb_paths]

            _merge_final_and_write(
                final_images=list(combined_all.images),
                pocket_inputs=[Path(p) for p in input_paths],
                ref_pdbs=ref_list_for_merge,
                segments=combined_all.segments,
                out_dir=out_dir_path,
                pocket_ref_pdbs=[Path(p) for p in pocket_ref_pdb_paths] if pocket_ref_pdb_paths else None,
            )
            click.echo("\n=== Full-system merge finished ===\n")

        # --------------------------
        # 5) Console summary
        # --------------------------
        try:
            overall_changed, overall_summary = _has_bond_change(combined_all.images[0], combined_all.images[-1], bond_cfg)
        except Exception:
            overall_changed, overall_summary = False, ""

        click.echo("\n=== MEP Summary ===\n")

        click.echo("\n[overall] Covalent-bond changes between first and last image:")
        if overall_changed and overall_summary.strip():
            click.echo(textwrap.indent(overall_summary.strip(), prefix="  "))
        else:
            click.echo("  (no covalent changes detected)")

        if combined_all.segments:
            click.echo("\n[segments] Along the final MEP order (ΔE‡, ΔE). Bridges are shown between connected segments:")
            for i, seg in enumerate(combined_all.segments, 1):
                kind_label = "BRIDGE" if seg.kind == "bridge" else "SEG"
                click.echo(f"  [{i:02d}] ({kind_label}) {seg.tag}  |  ΔE‡ = {seg.barrier_kcal:.2f} kcal/mol,  ΔE = {seg.delta_kcal:.2f} kcal/mol")
                if seg.kind != "bridge" and seg.summary.strip():
                    click.echo(textwrap.indent(seg.summary.strip(), prefix="      "))
        else:
            click.echo("\n[segments] (no segment reports)")

        # --------------------------
        # 6) Energy diagram (compressed state sequence)
        # --------------------------
        diagram_payload: Optional[Dict[str, Any]] = None
        try:
            # Map frames to segment indices (for anchoring R/P energies in Hartree)
            frame_seg_indices: List[int] = [int(getattr(im, "mep_seg_index", 0) or 0) for im in combined_all.images]
            seg_to_frames: Dict[int, List[int]] = {}
            for ii, sidx in enumerate(frame_seg_indices):
                if sidx <= 0:
                    continue
                seg_to_frames.setdefault(int(sidx), []).append(ii)

            # Bond‑change segments in the final MEP order
            bc_segments_in_order: List[SegmentReport] = [
                s for s in combined_all.segments
                if (s.kind == "seg" and s.summary and s.summary.strip() != "(no covalent changes detected)")
            ]

            # Determine which frames to use as R and P for anchoring energies in au
            start_idx_for_diag = 0
            end_idx_for_diag = len(combined_all.energies) - 1
            if bc_segments_in_order:
                first_bc = bc_segments_in_order[0]
                last_bc = bc_segments_in_order[-1]
                idxs_first_bc = seg_to_frames.get(int(first_bc.seg_index), [])
                idxs_last_bc = seg_to_frames.get(int(last_bc.seg_index), [])
                if idxs_first_bc:
                    start_idx_for_diag = int(idxs_first_bc[0])
                if idxs_last_bc:
                    end_idx_for_diag = int(idxs_last_bc[-1])

            E0_au = float(combined_all.energies[start_idx_for_diag])
            EP_au = float(combined_all.energies[end_idx_for_diag])

            # Build TS groups and compressed state energies purely from segment-level ΔE / ΔE‡
            ts_groups: List[Dict[str, Any]] = []
            ts_count = 0
            current: Optional[Dict[str, Any]] = None
            # Energy of the current state relative to R (in kcal/mol)
            E_current_kcal = 0.0

            def _is_bond_change_seg(seg: SegmentReport) -> bool:
                return (
                    seg.kind == "seg"
                    and bool(seg.summary)
                    and seg.summary.strip() != "(no covalent changes detected)"
                )

            for s in combined_all.segments:
                if _is_bond_change_seg(s):
                    ts_count += 1
                    barrier_kcal = float(getattr(s, "barrier_kcal", float("nan")))
                    delta_kcal = float(getattr(s, "delta_kcal", float("nan")))
                    if not np.isfinite(barrier_kcal):
                        barrier_kcal = 0.0
                    if not np.isfinite(delta_kcal):
                        delta_kcal = 0.0

                    ts_e = E_current_kcal + barrier_kcal
                    first_im_e = E_current_kcal + delta_kcal

                    current = {
                        "ts_label": f"TS{ts_count}",
                        "ts_energy": ts_e,
                        "first_im_energy": first_im_e,
                        "tail_im_energy": first_im_e,
                        "has_extra": False,
                        "index": ts_count,
                        "bridge_peaks": [],
                    }
                    ts_groups.append(current)
                    E_current_kcal = first_im_e
                else:
                    # Segments without covalent changes (bridge or kinks) belong to
                    # the current TS block if one exists.
                    if current is None:
                        # Before the first bond‑change segment: accumulate net ΔE if available.
                        delta_kcal = float(getattr(s, "delta_kcal", float("nan")))
                        if np.isfinite(delta_kcal):
                            E_current_kcal += delta_kcal
                        continue

                    delta_kcal = float(getattr(s, "delta_kcal", float("nan")))
                    barrier_kcal = float(getattr(s, "barrier_kcal", float("nan")))

                    if s.kind == "bridge":
                        if np.isfinite(barrier_kcal) and barrier_kcal > 1.0e-3:
                            peak_e = E_current_kcal + barrier_kcal
                            peaks = current.setdefault("bridge_peaks", [])
                            suffix = "" if len(peaks) == 0 else f"_{len(peaks) + 1}"
                            peak_label = f"IM{current['index']}_TS{suffix}"
                            peaks.append({"label": peak_label, "energy": peak_e})
                            # Log the corresponding peak in au for debugging
                            peak_e_au = E0_au + peak_e / AU2KCALPERMOL
                            click.echo(
                                "    [bridge] Recorded diagram-only TS peak "
                                f"{peak_label} at {peak_e_au:.6f} au "
                                "(from segment-level ΔE‡; bridge segments skip tsopt/thermo/DFT)."
                            )

                    if np.isfinite(delta_kcal):
                        E_current_kcal += delta_kcal
                        current["tail_im_energy"] = E_current_kcal
                        current["has_extra"] = True

            # Assemble labels and energies in kcal/mol
            labels: List[str]
            energies_kcal: List[float]
            chain_tokens: List[str]

            if not ts_groups:
                # No bond‑change segments: simple R→P diagram
                labels = ["R", "P"]
                energies_kcal = [0.0, (EP_au - E0_au) * AU2KCALPERMOL]
                chain_tokens = ["R", "-->", "P"]
            else:
                labels = ["R"]
                energies_kcal = [0.0]
                chain_tokens = ["R"]

                for i, g in enumerate(ts_groups, start=1):
                    last_group = (i == len(ts_groups))

                    labels.append(g["ts_label"])
                    energies_kcal.append(float(g["ts_energy"]))
                    chain_tokens.extend(["-->", g["ts_label"]])

                    if last_group:
                        continue

                    labels.append(f"IM{i}_1")
                    energies_kcal.append(float(g["first_im_energy"]))
                    chain_tokens.extend(["-->", f"IM{i}_1"])

                    for bp in g.get("bridge_peaks", []):
                        labels.append(str(bp["label"]))
                        energies_kcal.append(float(bp["energy"]))
                        chain_tokens.extend(["-->", str(bp["label"])])

                    if g["has_extra"]:
                        labels.append(f"IM{i}_2")
                        energies_kcal.append(float(g["tail_im_energy"]))
                        chain_tokens.extend(["-|-->", f"IM{i}_2"])

                labels.append("P")
                energies_kcal.append(E_current_kcal)
                chain_tokens.extend(["-->", "P"])

            # Convert back to Hartree (au) for completeness in the summary
            energies_au: List[float] = [E0_au + ek / AU2KCALPERMOL for ek in energies_kcal]

            diagram_payload = {
                "name": "energy_diagram_MEP",
                "labels": labels,
                "energies_kcal": energies_kcal,
                "ylabel": "ΔE (kcal/mol)",
                "energies_au": energies_au,
                "image": str(out_dir_path / "energy_diagram_MEP.png"),
            }

            labels_repr = "[" + ", ".join(f'"{lab}"' for lab in labels) + "]"
            energies_repr = "[" + ", ".join(f"{val:.6f}" for val in energies_kcal) + "]"
            click.echo(f"[diagram] build_energy_diagram.labels = {labels_repr}")
            click.echo(f"[diagram] build_energy_diagram.energies_kcal = {energies_repr}")

            title_note = "(MEP; all segments)"
            fig = build_energy_diagram(
                energies=energies_kcal,
                labels=labels,
                ylabel="ΔE (kcal/mol)",
                baseline=True,
                showgrid=False,
            )
            fig.update_layout(title=title_note)

            try:
                png_path = out_dir_path / "energy_diagram_MEP.png"
                fig.write_image(str(png_path), scale=2)
                click.echo(f"[diagram] Wrote energy diagram (PNG) → '{png_path}'")
            except Exception as e:
                click.echo(f"[diagram] NOTE: PNG export skipped (install 'kaleido' to enable): {e}", err=True)

            chain_text = " ".join(chain_tokens)
            click.echo(f"[diagram] State label sequence: {chain_text}")

        except Exception as e:
            click.echo(f"[diagram] WARNING: Failed to build energy diagram: {e}", err=True)

        # --------------------------
        # 7) Summary (YAML)
        # --------------------------
        summary = {
            "out_dir": str(out_dir_path),
            "n_images": len(combined_all.images),
            "n_segments": len(combined_all.segments),
            "segments": [
                {
                    "index": int(s.seg_index),
                    "tag": s.tag,
                    "kind": s.kind,
                    "barrier_kcal": float(s.barrier_kcal),
                    "delta_kcal": float(s.delta_kcal),
                    "bond_changes": (
                        _bond_changes_block(s.summary) if (s.kind != "bridge") else ""
                    ),
                } for s in combined_all.segments
            ],
        }
        if diagram_payload is not None:
            summary["energy_diagrams"] = [diagram_payload]

        with open(out_dir_path / "summary.yaml", "w") as f:
            yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
        click.echo(f"[write] Wrote '{out_dir_path / 'summary.yaml'}'.")

        try:
            try:
                freeze_atoms_for_log = sorted(
                    {
                        int(i)
                        for g in getattr(combined_all, "images", [])
                        for i in getattr(g, "freeze_atoms", [])
                    }
                )
            except Exception:
                pass

            diag_for_log: Dict[str, Any] = {}
            if diagram_payload is not None:
                diag_for_log = dict(diagram_payload)
            mep_info = {
                "n_images": len(combined_all.images),
                "n_segments": len(combined_all.segments),
                "traj_pdb": str(out_dir_path / "mep.pdb") if (out_dir_path / "mep.pdb").exists() else None,
                "mep_plot": str(out_dir_path / "mep_plot.png") if (out_dir_path / "mep_plot.png").exists() else None,
                "diagram": diag_for_log,
            }
            summary_payload = {
                "root_out_dir": str(out_dir_path),
                "path_dir": str(out_dir_path),
                "path_module_dir": "path_search",
                "pipeline_mode": mep_mode,
                "refine_path": True,
                "tsopt": False,
                "thermo": False,
                "dft": False,
                "opt_mode": opt_mode,
                "mep_mode": mep_mode,
                "uma_model": calc_cfg.get("model"),
                "command": command_str,
                "charge": calc_cfg.get("charge"),
                "spin": calc_cfg.get("spin"),
                "freeze_atoms": freeze_atoms_for_log,
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {
                    "summary.yaml": "YAML-format summary",
                    "summary.log": "This summary",
                    "mep_plot.png": "UMA MEP energy vs image index",
                    "energy_diagram_MEP.png": "Compressed MEP diagram R–TS–IM–P",
                },
            }
            write_summary_log(out_dir_path / "summary.log", summary_payload)
            click.echo(f"[write] Wrote '{out_dir_path / 'summary.log'}'.")
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

        # --------------------------
        # 8) Elapsed time
        # --------------------------
        click.echo(format_elapsed("[time] Elapsed for Path Search", time_start))

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
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
        _PRIMARY_GJF_TEMPLATE = None
