# pdb2reaction/scan3d.py

"""
scan3d — Three-distance 3D scan with harmonic restraints 
===================================================================

Usage (CLI)
-----------
    pdb2reaction scan3d -i INPUT.{pdb,xyz,trj,...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] \
        [-m MULTIPLICITY] \
        --scan-list '[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2),(I3,J3,LOW3,HIGH3)]' \
        [--one-based {True|False}] \
        [--max-step-size FLOAT] \
        [--bias-k FLOAT] \
        [--relax-max-cycles INT] \
        [--opt-mode {light,heavy}] \
        [--freeze-links {True|False}] \
        [--dump {True|False}] \
        [--convert-files {True|False}] \
        [--out-dir PATH] \
        [--csv PATH] \
        [--args-yaml FILE] \
        [--preopt {True|False}] \
        [--baseline {first|min}] \
        [--thresh {gau_loose|gau|gau_tight|gau_vtight|baker|never}] \
        [--zmin FLOAT] [--zmax FLOAT]

Examples
--------
    # Minimal example (three distance ranges)
    pdb2reaction scan3d -i input.pdb -q 0 \
        --scan-list '[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]'

    # LBFGS with inner-path trajectory dumping and 3D energy isosurface plot
    pdb2reaction scan3d -i input.pdb -q 0 \
        --scan-list '[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]' \
        --max-step-size 0.20 --dump True --out-dir ./result_scan3d/ --opt-mode light \
        --preopt True --baseline min

    # Plot only from an existing surface.csv (skip new energy evaluation)
    pdb2reaction scan3d -i input.pdb -q 0 \
        --scan-list '[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]' \
        --csv ./result_scan3d/surface.csv --out-dir ./result_scan3d/

Description
-----------
- A 3D grid scan driven by harmonic restraints on three inter-atomic distances (d1, d2, d3).
- Provide exactly one Python-like list
      [(i1, j1, low1, high1), (i2, j2, low2, high2), (i3, j3, low3, high3)]
  via **--scan-list**.
  - Indices are **1-based by default**; pass **--one-based False** to interpret them as 0-based.
  - For PDB inputs, each atom entry can be an integer index or a selector string such as
    ``'TYR,285,CA'`` or ``'MMT,309,C10'`` (resname, resseq, atom).
- `-q/--charge` is required for non-`.gjf` inputs **unless** ``--ligand-charge`` is provided; `.gjf` templates supply
  charge/spin when available. When ``-q`` is omitted but ``--ligand-charge`` is set, the full complex is treated as an
  enzyme–substrate system and the total charge is inferred using ``extract.py``’s residue-aware logic. Explicit ``-q``
  always overrides any derived charge.
  `-m/--multiplicity` specifies the spin multiplicity (2S+1) and defaults to 1 if omitted.
- Step schedule (h = `--max-step-size` in Å):
  - `N1 = ceil(|high1 - low1| / h)`, `N2 = ceil(|high2 - low2| / h)`, `N3 = ceil(|high3 - low3| / h)`.
  - `d1_values = linspace(low1, high1, N1 + 1)` (or `[low1]` if the span is ~0)
    `d2_values = linspace(low2, high2, N2 + 1)` (or `[low2]` if the span is ~0)
    `d3_values = linspace(low3, high3, N3 + 1)` (or `[low3]` if the span is ~0).
  - Internally, these {d1,d2,d3} value lists are **reordered** so that the values
    closest to the (optionally pre-optimized) starting structure are scanned first
    along each axis.

- Nested scan procedure (outer d1, middle d2, inner d3):
  1) For each `d1[i]`, relax with **only the d1 restraint active** (d1 bias only),
     starting from the previously scanned structure whose d1 value is closest.
  2) Snapshot that minimum, then for each `d2[j]` relax with **d1 and d2 restraints**,
     starting from the previously scanned structure at that d1 whose d2 value is closest.
  3) Snapshot that (d1,d2) minimum, then for each `d3[k]` relax with **d1, d2, d3 restraints**
     starting from the previously scanned structure at that (d1,d2) whose d3 value is closest.
     The **unbiased** energy is recorded (harmonic bias removed for evaluation).

- Plot-only mode:
  - If **--csv PATH** is provided, the 3D scan and energy evaluation are skipped.
  - The script reads the precomputed grid from the given CSV (surface.csv format)
    and only performs the 3D RBF interpolation and HTML plot generation.
  - If the CSV already has an `energy_kcal` column, it is used as-is.
    Otherwise, `energy_kcal` is reconstructed from `energy_hartree` and `--baseline`.
  - In this mode, the existing CSV is not overwritten and no new geometries are written;
    only the HTML visualization is generated under `--out-dir`.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_scan3d/)
  ├─ surface.csv                      # Grid metadata:
  │                                   # i,j,k,d1_A,d2_A,d3_A,energy_hartree,energy_kcal,bias_converged;
  │                                   # may also contain a reference row for the starting structure
  │                                   # with i=j=k=-1 (pre-optimized when --preopt True).
  ├─ scan3d_density.html              # 3D energy landscape (isosurface mesh only)
  └─ grid/
      ├─ point_iXXX_jYYY_kZZZ.xyz     # Constrained, relaxed geometries for each grid point
      │                               # XXX,YYY,ZZZ = int(round(d(Å)*100)), e.g. d=1.25Å → "125"
      ├─ point_iXXX_jYYY_kZZZ.pdb     # PDB companion when the input was PDB and conversion is enabled
      ├─ point_iXXX_jYYY_kZZZ.gjf     # GJF companion when a template is available and conversion is enabled
      ├─ preopt_iXXX_jYYY_kZZZ.xyz    # Starting structure used for the scan:
      │                               # pre-optimized when --preopt True, otherwise the input structure;
      │                               # same naming convention as above.
      ├─ preopt_iXXX_jYYY_kZZZ.pdb    # PDB companion for the starting structure when conversion is enabled
      ├─ preopt_iXXX_jYYY_kZZZ.gjf    # GJF companion for the starting structure when conversion is enabled
      └─ inner_path_d1_###_d2_###.trj # When --dump True; captures inner d3 paths per (d1,d2)
         inner_path_d1_###_d2_###.pdb  # PDB conversion of the inner-path trajectory when input was PDB and conversion enabled

Notes
-----
- UMA only (`uma_pysis` calculator) and the same `HarmonicBiasCalculator` used in the 1D/2D scan.
- Convergence is controlled by LBFGS or RFO depending on `--opt-mode` (default: `light`).
  Ångström limits are converted to Bohr to cap LBFGS step and RFO trust radii.
- `--baseline min|first`:
  - `min`   : shift PES so that the global minimum is 0 kcal/mol (**default**)
  - `first` : shift so that the grid point with `(i,j,k) = (0,0,0)` is 0 kcal/mol;
              if that point is missing, the global minimum is used instead.
- Format-aware XYZ/TRJ → PDB/GJF conversions respect the global `--convert-files {True|False}` toggle (default: enabled).
- The 3D visualization:
  - 3D RBF interpolation on a **50×50×50 grid** in (d1,d2,d3)-space.
  - Several semi-transparent isosurfaces (mesh) at discrete energy levels with **step color bands**
    (color is piecewise-constant in energy between contour levels).
  - No XY/YZ/ZX planes are drawn; the view is isosurface-only for clarity.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import ast
import math
import sys
import textwrap
import traceback
import tempfile
import os
import time

import click
import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import plotly.graph_objects as go

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .opt import (
    HarmonicBiasCalculator,
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
)
from .utils import (
    detect_freeze_links_safe,
    pretty_block,
    format_geom_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    normalize_choice,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    resolve_atom_spec_index,
)

# Default keyword dictionaries for the 3D scan (override only the knobs we touch)
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

OPT_BASE_KW: Dict[str, Any] = dict(_OPT_BASE_KW)
OPT_BASE_KW.update({
    "out_dir": "./result_scan3d/",  # base output directory for 3D scan artifacts
    "dump": False,                   # disable trajectory dumping by default
    "max_cycles": 10000,             # safety cap on optimization cycles
})

LBFGS_KW: Dict[str, Any] = dict(_LBFGS_KW)
LBFGS_KW.update({"out_dir": "./result_scan3d/"})  # directory for LBFGS-specific files

RFO_KW: Dict[str, Any] = dict(_RFO_KW)
RFO_KW.update({"out_dir": "./result_scan3d/"})    # directory for RFO-specific files

BIAS_KW: Dict[str, Any] = {"k": 100.0}  # harmonic restraint strength (eV/Å^2)

_OPT_MODE_ALIASES = (
    (("light",), "lbfgs"),
    (("heavy",), "rfo"),
)

HARTREE_TO_KCAL_MOL = 627.50961
_VOLUME_GRID_N = 50  # 50×50×50 RBF interpolation grid


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _snapshot_geometry(g) -> Any:
    """Snapshot a Geometry via temporary XYZ to decouple state."""
    s = g.as_xyz()
    if not s.endswith("\n"):
        s += "\n"
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()
        snap = geom_loader(
            Path(tmp.name),
            coord_type=getattr(g, "coord_type", GEOM_KW_DEFAULT["coord_type"]),
            freeze_atoms=getattr(g, "freeze_atoms", []),
        )
        try:
            import numpy as _np  # noqa: PLC0415
            snap.freeze_atoms = _np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        return snap
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _format_distance_tag(d: float) -> str:
    """
    Format a distance in Å as an integer tag, e.g. 1.25 → '125'.

    We multiply by 100 and round to the nearest integer, then pad to at
    least three digits, so 1.25 Å → '125', 0.95 Å → '095'.
    """
    val = int(round(float(d) * 100.0))
    return f"{val:03d}"


def _measure_distances_A(
    geom,
    i1: int,
    j1: int,
    i2: int,
    j2: int,
    i3: int,
    j3: int,
) -> Tuple[float, float, float]:
    """Measure the three bias distances (Å) for a given geometry."""
    coords = np.asarray(getattr(geom, "coords"), dtype=float).reshape(-1, 3)

    def dist(i: int, j: int) -> float:
        v = coords[i] - coords[j]
        return float(np.linalg.norm(v) / ANG2BOHR)

    return dist(i1, j1), dist(i2, j2), dist(i3, j3)


def _parse_scan_list(
    raw: str,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
) -> Tuple[
    Tuple[int, int, float, float],
    Tuple[int, int, float, float],
    Tuple[int, int, float, float],
    List[Tuple[Any, Any, float, float]],
]:
    """
    Parse --scan-list into three quadruples and return them with 0-based indices.

    Expected format:
        [(i1,j1,low1,high1),(i2,j2,low2,high2),(i3,j3,low3,high3)]

    Parameters
    ----------
    raw
        String passed to --scan-list.
    one_based
        If True, the indices in `raw` are interpreted as 1-based and converted
        to 0-based. If False, they are assumed to be 0-based already.
    """
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for --scan-list: {e}")

    if not (isinstance(obj, (list, tuple)) and len(obj) == 3):
        raise click.BadParameter(
            "--scan-list must contain exactly three quadruples: "
            "[(i1,j1,low1,high1),(i2,j2,low2,high2),(i3,j3,low3,high3)]"
        )

    def _resolve_index(value: Any, entry_idx: int, side_label: str) -> int:
        if isinstance(value, (int, np.integer)):
            idx_val = int(value)
            if one_based:
                idx_val -= 1
            if idx_val < 0:
                raise click.BadParameter(
                    f"Negative atom index after base conversion: {idx_val} (0-based expected)."
                )
            return idx_val
        if isinstance(value, str):
            if not atom_meta:
                raise click.BadParameter(
                    f"--scan-list entry {entry_idx} ({side_label}) uses a string atom spec, "
                    "but no PDB metadata is available."
                )
            try:
                return resolve_atom_spec_index(value, atom_meta)
            except ValueError as exc:
                raise click.BadParameter(
                    f"--scan-list entry {entry_idx} ({side_label}) {exc}"
                )
        raise click.BadParameter(
            f"--scan-list entry {entry_idx} ({side_label}) must be an int index or atom spec string."
        )

    parsed: List[Tuple[int, int, float, float]] = []
    for q in obj:
        if not (
            isinstance(q, (list, tuple)) and len(q) == 4
            and isinstance(q[2], (int, float, np.floating))
            and isinstance(q[3], (int, float, np.floating))
        ):
            raise click.BadParameter(f"--scan-list entry must be (i,j,low,high): got {q}")

        i = _resolve_index(q[0], len(parsed) + 1, "i")
        j = _resolve_index(q[1], len(parsed) + 1, "j")
        low, high = float(q[2]), float(q[3])
        if low <= 0.0 or high <= 0.0:
            raise click.BadParameter(f"Distances must be positive: {(i, j, low, high)}")
        parsed.append((i, j, low, high))

    return parsed[0], parsed[1], parsed[2], list(obj)


def _values_from_bounds(low: float, high: float, h: float) -> np.ndarray:
    if h <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    delta = abs(high - low)
    if delta < 1e-12:
        return np.array([low], dtype=float)
    N = int(math.ceil(delta / h))
    return np.linspace(low, high, N + 1, dtype=float)


def _atom_label_from_meta(atom_meta: Sequence[Dict[str, Any]], index: int) -> str:
    if index < 0 or index >= len(atom_meta):
        return f"idx{index}"
    meta = atom_meta[index]
    resname = (meta.get("resname") or "?").strip() or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    atom = (meta.get("name") or "?").strip() or "?"
    return f"{resname}-{resseq_txt}-{atom}"


def _axis_label_csv(
    axis_name: str,
    i_idx: int,
    j_idx: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
    pair_raw: Optional[Tuple[Any, Any, float, float]] = None,
) -> str:
    """Return a CSV-safe axis label without commas."""
    if pair_raw and (isinstance(pair_raw[0], str) or isinstance(pair_raw[1], str)) and atom_meta:
        i_label = _atom_label_from_meta(atom_meta, i_idx)
        j_label = _atom_label_from_meta(atom_meta, j_idx)
        return f"{axis_name}_{i_label}_{j_label}_A"
    i_disp = i_idx + 1 if one_based else i_idx
    j_disp = j_idx + 1 if one_based else j_idx
    return f"{axis_name}_{i_disp}_{j_disp}_A"


def _axis_label_html(label: str) -> str:
    """Return a human-readable axis label for HTML output."""
    parts = label.split("_")
    if len(parts) >= 4 and parts[-1] == "A":
        axis = parts[0]
        i_disp = parts[1]
        j_disp = parts[2]
        return f"{axis} ({i_disp},{j_disp}) (Å)"
    return label


def _extract_axis_label(df: pd.DataFrame, column: str, fallback: Optional[str]) -> Optional[str]:
    if column not in df.columns:
        return fallback
    values = df[column].dropna()
    if values.empty:
        return fallback
    return str(values.iloc[0])


def _make_optimizer(
    geom,
    kind: str,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    max_step_bohr: float,
    relax_max_cycles: int,
    out_dir: Path,
    prefix: str,
):
    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    if kind == "lbfgs":
        args = {**lbfgs_cfg, **common}
        args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
        args["max_cycles"] = int(relax_max_cycles)
        return LBFGS(geom, **args)
    else:
        args = {**rfo_cfg, **common}
        tr = float(rfo_cfg.get("trust_radius", 0.30))
        args["trust_radius"] = min(tr, max_step_bohr)
        args["trust_max"] = min(float(rfo_cfg.get("trust_max", 0.30)), max_step_bohr)
        args["max_cycles"] = int(relax_max_cycles)
        return RFOptimizer(geom, **args)


def _unbiased_energy_hartree(geom, base_calc) -> float:
    """Evaluate UMA energy (Hartree) without harmonic bias."""
    coords_bohr = np.asarray(geom.coords)
    elems = getattr(geom, "atoms", None)

    if elems is None:
        return float("nan")

    try:
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception:
        return float("nan")


@click.command(
    help="3D distance scan with harmonic restraints.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="Input structure file (.pdb, .xyz, .trj, ...). Required unless --csv is provided.",
)
@click.option("-q", "--charge", type=int, required=False, help="Charge of the ML region.")
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
    help="Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) for unknown residues.",
)
@click.option(
    "-m", "--multiplicity", "spin",
    type=int,
    default=1,
    show_default=True,
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--scan-list", "scan_list_raw",
    type=str,
    required=False,
    help=(
        "Python-like list with three quadruples: "
        "'[(i1,j1,low1,high1),(i2,j2,low2,high2),(i3,j3,low3,high3)]'."
    ),
)
@click.option(
    "--one-based",
    "one_based",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Interpret (i,j) indices in --scan-list as 1-based (default) or 0-based.",
)
@click.option(
    "--max-step-size",
    type=float,
    default=0.20,
    show_default=True,
    help="Maximum step size in each distance [Å].",
)
@click.option(
    "--bias-k",
    type=float,
    default=100.0,
    show_default=True,
    help="Harmonic well strength k [eV/Å^2].",
)
@click.option(
    "--relax-max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum optimizer cycles per grid relaxation.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Relaxation mode: light (=LBFGS) or heavy (=RFO).",
)
@click.option(
    "--freeze-links",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If input is PDB, freeze parent atoms of link hydrogens.",
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Write inner d3 scan trajectories per (d1,d2) as TRJ under result_scan3d/grid/.",
)
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option(
    "--out-dir",
    type=str,
    default="./result_scan3d/",
    show_default=True,
    help="Base output directory.",
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "If provided, skip the 3D scan and read a precomputed surface.csv from this path. "
        "Only plotting is performed."
    ),
)
@click.option(
    "--thresh",
    type=str,
    default="baker",
    show_default=False,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file with extra args (sections: geom, calc, opt, lbfgs, rfo, bias).",
)
@click.option(
    "--preopt",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Pre-optimize the initial structure without bias before the scan.",
)
@click.option(
    "--baseline",
    type=click.Choice(["min", "first"]),
    default="min",
    show_default=True,
    help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0,k=0).",
)
@click.option(
    "--zmin",
    type=float,
    default=None,
    show_default=False,
    help="Lower bound of color scale for plots (kcal/mol).",
)
@click.option(
    "--zmax",
    type=float,
    default=None,
    show_default=False,
    help="Upper bound of color scale for plots (kcal/mol).",
)
def cli(
    input_path: Optional[Path],
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    scan_list_raw: Optional[str],
    one_based: bool,
    max_step_size: float,
    bias_k: float,
    relax_max_cycles: int,
    opt_mode: str,
    freeze_links: bool,
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    csv_path: Optional[Path],
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    baseline: str,
    zmin: Optional[float],
    zmax: Optional[float],
) -> None:
    from .utils import load_yaml_dict, apply_yaml_overrides

    set_convert_file_enabled(convert_files)
    prepared_input = None
    geom_input_path = None
    source_path = None
    if csv_path is None:
        if input_path is None:
            raise click.ClickException("-i/--input is required unless --csv is provided.")
        if scan_list_raw is None:
            raise click.ClickException("--scan-list is required unless --csv is provided.")
        prepared_input = prepare_input_structure(input_path)
        apply_ref_pdb_override(prepared_input, ref_pdb)
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path

        charge, spin = resolve_charge_spin_or_raise(
            prepared_input,
            charge,
            spin,
            ligand_charge=ligand_charge,
            prefix="[scan3d]",
        )

    try:
        time_start = time.perf_counter()

        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg = dict(RFO_KW)
        bias_cfg = dict(BIAS_KW)

        # CLI overrides (defaults ← CLI)
        if csv_path is None:
            calc_cfg["charge"] = int(charge)
            calc_cfg["spin"] = int(spin)
        calc_cfg["workers"] = int(workers)
        calc_cfg["workers_per_node"] = int(workers_per_node)
        opt_cfg["out_dir"] = out_dir
        opt_cfg["dump"] = False
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)

        # YAML overrides (highest precedence)
        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",),)),
                (rfo_cfg, (("rfo",),)),
                (bias_cfg, (("bias",),)),
            ],
        )

        if bias_k is not None:
            bias_cfg["k"] = float(bias_k)

        kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=_OPT_MODE_ALIASES,
            allowed_hint="light|heavy",
        )

        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        _ensure_dir(out_dir_path)
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_freeze_atoms_for_echo(calc_cfg)
        echo_opt = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)
        echo_bias = dict(bias_cfg)
        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("opt", echo_opt))
        click.echo(
            pretty_block("lbfgs" if kind == "lbfgs" else "rfo", (lbfgs_cfg if kind == "lbfgs" else rfo_cfg))
        )
        click.echo(pretty_block("bias", echo_bias))

        pdb_atom_meta: List[Dict[str, Any]] = []
        d1_label_csv = None
        d2_label_csv = None
        d3_label_csv = None
        if csv_path is None:
            if source_path and source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            (
                (i1, j1, low1, high1),
                (i2, j2, low2, high2),
                (i3, j3, low3, high3),
                raw_pairs,
            ) = _parse_scan_list(scan_list_raw, one_based=one_based, atom_meta=pdb_atom_meta)
            d1_label_csv = _axis_label_csv("d1", i1, j1, one_based, pdb_atom_meta, raw_pairs[0])
            d2_label_csv = _axis_label_csv("d2", i2, j2, one_based, pdb_atom_meta, raw_pairs[1])
            d3_label_csv = _axis_label_csv("d3", i3, j3, one_based, pdb_atom_meta, raw_pairs[2])
            click.echo(
                pretty_block(
                    "scan-list (0-based)",
                    {
                        "d1": (i1, j1, low1, high1),
                        "d2": (i2, j2, low2, high2),
                        "d3": (i3, j3, low3, high3),
                    },
                )
            )

            if pdb_atom_meta:
                click.echo("[scan3d] PDB atom details for scanned pairs:")
                legend = format_pdb_atom_metadata_header()
                click.echo(f"        legend: {legend}")
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}")
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}")
                click.echo(f"  d3 i: {format_pdb_atom_metadata(pdb_atom_meta, i3)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j3)}")

        final_dir = out_dir_path

        ref_pdb_path = None
        if csv_path is None and source_path and source_path.suffix.lower() == ".pdb":
            ref_pdb_path = source_path

        # ==== Either load existing surface.csv, or run the full 3D scan ====
        if csv_path is not None:
            csv_path = Path(csv_path).resolve()
            try:
                df = pd.read_csv(csv_path)
            except Exception as e:
                click.echo(f"[read] Failed to read CSV '{csv_path}': {e}", err=True)
                sys.exit(1)
            click.echo(f"[read] Loaded precomputed grid from '{csv_path}'.")
        else:
            tmp_root = Path(tempfile.mkdtemp(prefix="scan3d_tmp_"))
            grid_dir = out_dir_path / "grid"
            tmp_opt_dir = tmp_root / "opt"
            _ensure_dir(grid_dir)
            _ensure_dir(tmp_opt_dir)

            freeze = merge_freeze_atom_indices(geom_cfg)
            if freeze_links and source_path and source_path.suffix.lower() == ".pdb":
                detected = detect_freeze_links_safe(source_path)
                if detected:
                    freeze = merge_freeze_atom_indices(geom_cfg, detected)
                    if freeze:
                        click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, freeze))}")
            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            geom_outer = geom_loader(
                geom_input_path, coord_type=coord_type, freeze_atoms=freeze
            )
            if freeze:
                try:
                    import numpy as _np  # noqa: PLC0415
                    geom_outer.freeze_atoms = _np.array(freeze, dtype=int)
                except Exception:
                    pass

            base_calc = uma_pysis(**calc_cfg)
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

            # Optional pre-optimization of the starting structure
            if preopt:
                click.echo("[preopt] Unbiased relaxation of the initial structure ...")
                geom_outer.set_calculator(base_calc)
                max_step_bohr_local = float(max_step_size) * ANG2BOHR
                optimizer0 = _make_optimizer(
                    geom_outer,
                    kind,
                    lbfgs_cfg,
                    rfo_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr_local,
                    relax_max_cycles=relax_max_cycles,
                    out_dir=tmp_opt_dir,
                    prefix="preopt_",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.", err=True)
                except OptimizationError as e:
                    click.echo(f"[preopt] OptimizationError — {e}", err=True)

            # Measure the three bias distances on the starting structure
            # (pre-optimized when --preopt True, otherwise the input geometry)
            d1_ref, d2_ref, d3_ref = _measure_distances_A(
                geom_outer, i1, j1, i2, j2, i3, j3
            )
            click.echo(
                pretty_block(
                    "preopt distances (Å)",
                    {"d1_ref": d1_ref, "d2_ref": d2_ref, "d3_ref": d3_ref},
                )
            )

            # Save the starting structure (pre-optimized when requested) and prepare a record for plotting
            preopt_tag_i = _format_distance_tag(d1_ref)
            preopt_tag_j = _format_distance_tag(d2_ref)
            preopt_tag_k = _format_distance_tag(d3_ref)
            preopt_xyz_path = grid_dir / f"preopt_i{preopt_tag_i}_j{preopt_tag_j}_k{preopt_tag_k}.xyz"
            try:
                s_pre = geom_outer.as_xyz()
                if not s_pre.endswith("\n"):
                    s_pre += "\n"
                with open(preopt_xyz_path, "w") as f:
                    f.write(s_pre)
                click.echo(f"[preopt] Wrote '{preopt_xyz_path}'.")
            except Exception as e:
                click.echo(f"[preopt] WARNING: failed to write '{preopt_xyz_path.name}': {e}", err=True)

            try:
                convert_xyz_like_outputs(
                    preopt_xyz_path,
                    prepared_input,
                    ref_pdb_path=ref_pdb_path,
                    out_pdb_path=grid_dir / f"preopt_i{preopt_tag_i}_j{preopt_tag_j}_k{preopt_tag_k}.pdb",
                    out_gjf_path=grid_dir / f"preopt_i{preopt_tag_i}_j{preopt_tag_j}_k{preopt_tag_k}.gjf",
                )
            except Exception as e:
                click.echo(
                    f"[convert] WARNING: failed to convert '{preopt_xyz_path.name}' to PDB/GJF: {e}",
                    err=True,
                )

            E_pre_h = _unbiased_energy_hartree(geom_outer, base_calc)
            preopt_record = {
                "i": -1,
                "j": -1,
                "k": -1,
                "d1_A": float(d1_ref),
                "d2_A": float(d2_ref),
                "d3_A": float(d3_ref),
                "energy_hartree": float(E_pre_h),
                "bias_converged": True,
            }

            # Build and reorder the grids so that we scan from values closest to the reference distances
            d1_values = _values_from_bounds(low1, high1, float(max_step_size))
            d2_values = _values_from_bounds(low2, high2, float(max_step_size))
            d3_values = _values_from_bounds(low3, high3, float(max_step_size))

            d1_values = np.array(
                sorted(d1_values, key=lambda v: abs(v - d1_ref)),
                dtype=float,
            )
            d2_values = np.array(
                sorted(d2_values, key=lambda v: abs(v - d2_ref)),
                dtype=float,
            )
            d3_values = np.array(
                sorted(d3_values, key=lambda v: abs(v - d3_ref)),
                dtype=float,
            )

            N1, N2, N3 = len(d1_values), len(d2_values), len(d3_values)
            click.echo(f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x: f'{x:.3f}', d1_values))}")
            click.echo(f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x: f'{x:.3f}', d2_values))}")
            click.echo(f"[grid] d3 steps = {N3}  values(A)={list(map(lambda x: f'{x:.3f}', d3_values))}")
            click.echo(f"[grid] total grid points = {N1 * N2 * N3}")

            max_step_bohr = float(max_step_size) * ANG2BOHR

            records: List[Dict[str, Any]] = []

            # Store starting geometry snapshot for reuse
            geom_outer_initial = _snapshot_geometry(geom_outer)

            # Caches for nearest-neighbor starting geometries
            d1_geoms: Dict[int, Any] = {}
            d2_geoms: Dict[int, Dict[int, Any]] = {}
            d3_geoms: Dict[Tuple[int, int], Dict[int, Any]] = {}

            # ===== 3D nested scan: d1 (outer) → d2 (middle) → d3 (inner) =====
            for i_idx, d1_target in enumerate(d1_values):
                click.echo(f"\n=== d1 step {i_idx + 1}/{N1} : target = {d1_target:.3f} Å ===")

                # Choose initial geometry for this d1 from the previously scanned
                # structure with the closest d1 value (or the reference structure).
                if not d1_geoms:
                    geom_outer_i = _snapshot_geometry(geom_outer_initial)
                else:
                    nearest_i = min(
                        d1_geoms.keys(),
                        key=lambda p: abs(d1_values[p] - d1_target),
                    )
                    geom_outer_i = _snapshot_geometry(d1_geoms[nearest_i])

                biased.set_pairs([(i1, j1, float(d1_target))])
                geom_outer_i.set_calculator(biased)

                opt1 = _make_optimizer(
                    geom_outer_i,
                    kind,
                    lbfgs_cfg,
                    rfo_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr,
                    relax_max_cycles=relax_max_cycles,
                    out_dir=tmp_opt_dir,
                    prefix=f"d1_{i_idx:03d}_",
                )
                try:
                    opt1.run()
                except ZeroStepLength:
                    click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2/d3 scan.", err=True)
                except OptimizationError as e:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {e}", err=True)

                # Snapshot after d1 relaxation for inner loops and cache
                geom_after_d1 = _snapshot_geometry(geom_outer_i)
                d1_geoms[i_idx] = geom_after_d1

                if i_idx not in d2_geoms:
                    d2_geoms[i_idx] = {}

                for j_idx, d2_target in enumerate(d2_values):
                    click.echo(
                        f"\n--- (d1,d2) = ({i_idx + 1}/{N1}, {j_idx + 1}/{N2}) : "
                        f"targets = ({d1_target:.3f}, {d2_target:.3f}) Å ---"
                    )

                    # Choose initial geometry for this (d1,d2) from the previously
                    # scanned structure with the closest d2 value at this d1
                    # (or the d1-relaxed structure).
                    d2_store = d2_geoms[i_idx]
                    if not d2_store:
                        geom_mid = _snapshot_geometry(geom_after_d1)
                    else:
                        nearest_j = min(
                            d2_store.keys(),
                            key=lambda p: abs(d2_values[p] - d2_target),
                        )
                        geom_mid = _snapshot_geometry(d2_store[nearest_j])

                    biased.set_pairs(
                        [
                            (i1, j1, float(d1_target)),
                            (i2, j2, float(d2_target)),
                        ]
                    )
                    geom_mid.set_calculator(biased)

                    opt2 = _make_optimizer(
                        geom_mid,
                        kind,
                        lbfgs_cfg,
                        rfo_cfg,
                        opt_cfg,
                        max_step_bohr=max_step_bohr,
                        relax_max_cycles=relax_max_cycles,
                        out_dir=tmp_opt_dir,
                        prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}_",
                    )
                    try:
                        opt2.run()
                    except ZeroStepLength:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — continuing to d3 scan.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {e}", err=True)

                    geom_after_d2 = _snapshot_geometry(geom_mid)
                    d2_store[j_idx] = geom_after_d2

                    key_ij = (i_idx, j_idx)
                    if key_ij not in d3_geoms:
                        d3_geoms[key_ij] = {}
                    d3_store = d3_geoms[key_ij]

                    trj_blocks = [] if dump else None

                    for k_idx, d3_target in enumerate(d3_values):
                        # Choose initial geometry for this (d1,d2,d3) from the
                        # previously scanned structure with the closest d3 value
                        # at this (d1,d2), or from the d1/d2-relaxed structure.
                        if not d3_store:
                            geom_inner = _snapshot_geometry(geom_after_d2)
                        else:
                            nearest_k = min(
                                d3_store.keys(),
                                key=lambda p: abs(d3_values[p] - d3_target),
                            )
                            geom_inner = _snapshot_geometry(d3_store[nearest_k])

                        biased.set_pairs(
                            [
                                (i1, j1, float(d1_target)),
                                (i2, j2, float(d2_target)),
                                (i3, j3, float(d3_target)),
                            ]
                        )
                        geom_inner.set_calculator(biased)

                        opt3 = _make_optimizer(
                            geom_inner,
                            kind,
                            lbfgs_cfg,
                            rfo_cfg,
                            opt_cfg,
                            max_step_bohr=max_step_bohr,
                            relax_max_cycles=relax_max_cycles,
                            out_dir=tmp_opt_dir,
                            prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}_d3_{k_idx:03d}_",
                        )
                        try:
                            opt3.run()
                            converged = True
                        except ZeroStepLength:
                            click.echo(
                                f"[d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] ZeroStepLength — recorded anyway.",
                                err=True,
                            )
                            converged = False
                        except OptimizationError as e:
                            click.echo(
                                f"[d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] OptimizationError — {e}",
                                err=True,
                            )
                            converged = False

                        # Cache final geometry for nearest-neighbor reuse
                        d3_store[k_idx] = _snapshot_geometry(geom_inner)

                        E_h = _unbiased_energy_hartree(geom_inner, base_calc)

                        tag_i = _format_distance_tag(d1_target)
                        tag_j = _format_distance_tag(d2_target)
                        tag_k = _format_distance_tag(d3_target)
                        xyz_path = grid_dir / f"point_i{tag_i}_j{tag_j}_k{tag_k}.xyz"
                        try:
                            s = geom_inner.as_xyz()
                            if not s.endswith("\n"):
                                s += "\n"
                            with open(xyz_path, "w") as f:
                                f.write(s)
                        except Exception as e:
                            click.echo(f"[write] WARNING: failed to write {xyz_path.name}: {e}", err=True)
                        else:
                            try:
                                convert_xyz_like_outputs(
                                    xyz_path,
                                    prepared_input,
                                    ref_pdb_path=ref_pdb_path,
                                    out_pdb_path=grid_dir / f"point_i{tag_i}_j{tag_j}_k{tag_k}.pdb",
                                    out_gjf_path=grid_dir / f"point_i{tag_i}_j{tag_j}_k{tag_k}.gjf",
                                )
                            except Exception as e:
                                click.echo(
                                    f"[convert] WARNING: failed to convert '{xyz_path.name}' to PDB/GJF: {e}",
                                    err=True,
                                )

                        if dump and trj_blocks is not None:
                            sblock = geom_inner.as_xyz()
                            if not sblock.endswith("\n"):
                                sblock += "\n"
                            trj_blocks.append(sblock)

                        records.append(
                            {
                                "i": int(i_idx),
                                "j": int(j_idx),
                                "k": int(k_idx),
                                "d1_A": float(d1_target),
                                "d2_A": float(d2_target),
                                "d3_A": float(d3_target),
                                "energy_hartree": E_h,
                                "bias_converged": bool(converged),
                            }
                        )

                    if dump and trj_blocks:
                        trj_path = grid_dir / f"inner_path_d1_{i_idx:03d}_d2_{j_idx:03d}.trj"
                        try:
                            with open(trj_path, "w") as f:
                                f.write("".join(trj_blocks))
                            click.echo(f"[write] Wrote '{trj_path}'.")
                        except Exception as e:
                            click.echo(f"[write] WARNING: failed to write '{trj_path}': {e}", err=True)
                        else:
                            try:
                                convert_xyz_like_outputs(
                                    trj_path,
                                    prepared_input,
                                    ref_pdb_path=ref_pdb_path,
                                    out_pdb_path=grid_dir / f"inner_path_d1_{i_idx:03d}_d2_{j_idx:03d}.pdb",
                                )
                            except Exception as e:
                                click.echo(
                                    f"[convert] WARNING: failed to convert '{trj_path.name}' to PDB: {e}",
                                    err=True,
                                )

            # Add starting structure as an extra record for plotting
            records.append(preopt_record)

            # ===== surface.csv =====
            df = pd.DataFrame.from_records(records)

        # ===== surface.csv handling & baseline =====
        if df.empty:
            click.echo("No grid records produced; aborting.", err=True)
            sys.exit(1)

        d1_label_csv = _extract_axis_label(df, "d1_label", d1_label_csv)
        d2_label_csv = _extract_axis_label(df, "d2_label", d2_label_csv)
        d3_label_csv = _extract_axis_label(df, "d3_label", d3_label_csv)

        if d1_label_csv is None or d2_label_csv is None or d3_label_csv is None:
            click.echo(
                "[plot] WARNING: axis label metadata is missing in CSV; using generic labels.",
                err=True,
            )

        d1_label_html = _axis_label_html(d1_label_csv) if d1_label_csv else "d1 (Å)"
        d2_label_html = _axis_label_html(d2_label_csv) if d2_label_csv else "d2 (Å)"
        d3_label_html = _axis_label_html(d3_label_csv) if d3_label_csv else "d3 (Å)"

        # If energy_kcal is already present (e.g. loaded from existing CSV), reuse it.
        # Otherwise compute it from energy_hartree and baseline.
        if "energy_kcal" not in df.columns:
            if "energy_hartree" not in df.columns:
                click.echo(
                    "[baseline] energy_kcal is missing and energy_hartree is not available in CSV; aborting.",
                    err=True,
                )
                sys.exit(1)

            if baseline == "first":
                ref_mask = (df["i"] == 0) & (df["j"] == 0) & (df["k"] == 0)
                if not ref_mask.any():
                    click.echo(
                        "[baseline] 'first' requested but (i=0,j=0,k=0) missing; using global minimum instead.",
                        err=True,
                    )
                    ref = float(df["energy_hartree"].min())
                else:
                    ref = float(df.loc[ref_mask, "energy_hartree"].iloc[0])
            else:
                ref = float(df["energy_hartree"].min())

            df["energy_kcal"] = (df["energy_hartree"] - ref) * HARTREE_TO_KCAL_MOL

        # Only write surface.csv when we actually performed the scan in this run
        if csv_path is None:
            surface_csv = final_dir / "surface.csv"
            df["d1_label"] = d1_label_csv
            df["d2_label"] = d2_label_csv
            df["d3_label"] = d3_label_csv
            df.to_csv(surface_csv, index=False)
            click.echo(f"[write] Wrote '{surface_csv}'.")

        # ===== 3D RBF interpolation & visualization (isosurface only) =====
        d1_points = df["d1_A"].to_numpy(dtype=float)
        d2_points = df["d2_A"].to_numpy(dtype=float)
        d3_points = df["d3_A"].to_numpy(dtype=float)
        z_points = df["energy_kcal"].to_numpy(dtype=float)

        mask = (
            np.isfinite(d1_points)
            & np.isfinite(d2_points)
            & np.isfinite(d3_points)
            & np.isfinite(z_points)
        )
        if not np.any(mask):
            click.echo("[plot] No finite data for plotting.", err=True)
            sys.exit(1)

        x_min, x_max = float(np.min(d1_points[mask])), float(np.max(d1_points[mask]))
        y_min, y_max = float(np.min(d2_points[mask])), float(np.max(d2_points[mask]))
        z_min_val, z_max_val = float(np.min(d3_points[mask])), float(np.max(d3_points[mask]))

        xi = np.linspace(x_min, x_max, _VOLUME_GRID_N)
        yi = np.linspace(y_min, y_max, _VOLUME_GRID_N)
        zi = np.linspace(z_min_val, z_max_val, _VOLUME_GRID_N)

        click.echo("[plot] 3D RBF interpolation on a 50×50×50 grid ...")
        rbf3d = Rbf(
            d1_points[mask],
            d2_points[mask],
            d3_points[mask],
            z_points[mask],
            function="multiquadric",
        )

        XI, YI, ZI = np.meshgrid(xi, yi, zi, indexing="xy")
        X_flat = XI.flatten()
        Y_flat = YI.flatten()
        Z_flat = ZI.flatten()
        E_flat = rbf3d(X_flat, Y_flat, Z_flat)

        vmin = float(np.nanmin(E_flat)) if zmin is None else float(zmin)
        vmax = float(np.nanmax(E_flat)) if zmax is None else float(zmax)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = float(np.nanmin(E_flat)), float(np.nanmax(E_flat))

        # Discrete isosurfaces (mesh) with banded colors (no XY/YZ/ZX planes)
        n_levels = 8
        level_values = np.linspace(vmin, vmax, n_levels + 2)[1:-1]
        level_colors = [
            "#0d0887",
            "#5b02a3",
            "#9c179e",
            "#cb4679",
            "#ed7953",
            "#fb9f3a",
            "#fdca26",
            "#f0f921",
        ]

        # One opacity per isosurface level (outermost surfaces more transparent)
        level_opacity = [
            1.000,
            0.667,
            0.444,
            0.296,
            0.198,
            0.132,
            0.088,
            0.059,
        ]

        isosurfaces = []
        for lvl, color, opacity_lvl in zip(level_values, level_colors, level_opacity):
            trace = go.Isosurface(
                x=X_flat,
                y=Y_flat,
                z=Z_flat,
                value=E_flat,
                isomin=lvl,
                isomax=lvl,
                surface_count=1,
                opacity=opacity_lvl,
                showscale=False,
                colorscale=[[0.0, color], [1.0, color]],
                caps=dict(x_show=False, y_show=False, z_show=False),
                name=f"{lvl:.1f} kcal/mol",
            )
            isosurfaces.append(trace)

        # Add a dummy scatter trace to host a global colorbar
        colorbar_colorscale = [
            [idx / (len(level_colors) - 1), col]
            for idx, col in enumerate(level_colors)
        ]
        cb_tickvals = [float(v) for v in level_values]
        cb_ticktext = [f"{v:.1f}" for v in level_values]

        colorbar_trace = go.Scatter3d(
            x=[x_min],
            y=[y_min],
            z=[z_min_val],
            mode="markers",
            marker=dict(
                size=0,
                opacity=0.0,
                color=[vmin, vmax],
                colorscale=colorbar_colorscale,
                showscale=True,
                colorbar=dict(
                    title=dict(text="(kcal/mol)", side="top", font=dict(size=16, color="#1C1C1C")),
                    tickfont=dict(size=14, color="#1C1C1C"),
                    ticks="inside",
                    ticklen=10,
                    tickcolor="#1C1C1C",
                    outlinecolor="#1C1C1C",
                    outlinewidth=2,
                    lenmode="fraction",
                    len=1.11,
                    x=1.05,
                    y=0.53,
                    xanchor="left",
                    yanchor="middle",
                    tickvals=cb_tickvals,
                    ticktext=cb_ticktext,
                ),
            ),
            hoverinfo="none",
            showlegend=False,
        )

        fig3d = go.Figure(data=isosurfaces + [colorbar_trace])

        fig3d.update_layout(
            title="3D Energy Landscape (UMA)",
            width=900,
            height=800,
            scene=dict(
                bgcolor="rgba(0,0,0,0)",
                xaxis=dict(
                    title=d1_label_html,
                    range=[x_min, x_max],
                    showline=True,
                    linewidth=4,
                    linecolor="#1C1C1C",
                    mirror=True,
                    ticks="inside",
                    tickwidth=4,
                    tickcolor="#1C1C1C",
                    gridcolor="rgba(0,0,0,0.1)",
                    zerolinecolor="rgba(0,0,0,0.1)",
                    showbackground=False,
                ),
                yaxis=dict(
                    title=d2_label_html,
                    range=[y_min, y_max],
                    showline=True,
                    linewidth=4,
                    linecolor="#1C1C1C",
                    mirror=True,
                    ticks="inside",
                    tickwidth=4,
                    tickcolor="#1C1C1C",
                    gridcolor="rgba(0,0,0,0.1)",
                    zerolinecolor="rgba(0,0,0,0.1)",
                    showbackground=False,
                ),
                zaxis=dict(
                    title=d3_label_html,
                    range=[z_min_val, z_max_val],
                    showline=True,
                    linewidth=4,
                    linecolor="#1C1C1C",
                    mirror=True,
                    ticks="inside",
                    tickwidth=4,
                    tickcolor="#1C1C1C",
                    gridcolor="rgba(0,0,0,0.1)",
                    zerolinecolor="rgba(0,0,0,0.1)",
                    showbackground=False,
                ),
                aspectmode="cube",
            ),
            margin=dict(l=10, r=20, b=10, t=40),
            paper_bgcolor="white",
        )

        html3d = final_dir / "scan3d_density.html"
        fig3d.write_html(str(html3d))
        click.echo(f"[plot] Wrote '{html3d}'.")

        click.echo("\n=== 3D Scan finished ===\n")
        click.echo(format_elapsed("[time] Elapsed Time for 3D Scan", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled exception during 3D scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        if prepared_input:
            prepared_input.cleanup()
