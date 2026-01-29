# pdb2reaction/utils.py

"""
Utilities for configuration, plotting, coordinates, Gaussian templates, and link-freezing.

Categories:
    - Generic helpers: pretty_block, format_elapsed, deep_update, load_yaml_dict, etc.
    - Gaussian (.gjf): parse_gjf_template, prepare_input_structure, convert_xyz_to_gjf
    - Plotting: build_energy_diagram (Plotly-based energy diagrams)
    - Coordinate conversion: convert_xyz_to_pdb (XYZ to PDB with reference topology)
    - Link-freezing: detect_freeze_links, resolve_freeze_atoms (for PDB link hydrogens)

Dependencies: PyYAML, ASE, Plotly, Click
"""

import ast
import functools
import math
import os
import re
import tempfile
import time
from collections.abc import Iterable as _Iterable, Mapping, Sequence as _Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from numbers import Real, Integral
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, List, Tuple, Callable, Iterator

import click
from click.core import ParameterSource
import numpy as np
import yaml
from ase.data import chemical_symbols
from ase.io import read, write
import plotly.graph_objs as go

from .add_elem_info import guess_element
from pysisyphus.constants import AU2KCALPERMOL, ANG2BOHR
from pysisyphus.helpers import geom_loader

# =============================================================================
# YAML helpers (shared representers)
# =============================================================================


class YamlLiteralStr(str):
    """String marker to force literal block style when dumping YAML."""


class YamlFlowList(list):
    """List marker to force flow style when dumping YAML."""


def register_yaml_representers() -> None:
    """Register shared YAML representers (literal strings and flow lists)."""
    yaml.add_representer(
        YamlLiteralStr,
        lambda dumper, data: dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
    )
    yaml.add_representer(
        YamlLiteralStr,
        lambda dumper, data: dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|"),
        Dumper=yaml.SafeDumper
    )
    yaml.SafeDumper.add_representer(
        YamlFlowList,
        lambda dumper, data: dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)
    )


register_yaml_representers()

# =============================================================================
# Generic helpers
# =============================================================================


def pretty_block(title: str, content: Dict[str, Any]) -> str:
    """
    Return a YAML-formatted block with an underlined title.

    Empty mappings render as "{}".
    """
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + body + "\n"


def strip_inherited_keys(
    child_cfg: Dict[str, Any],
    base_cfg: Mapping[str, Any],
    *,
    mode: str = "present",
) -> Dict[str, Any]:
    """Return child_cfg without inherited keys (for concise logs)."""
    trimmed: Dict[str, Any] = {}
    if mode not in {"present", "same"}:
        raise ValueError(f"Unknown strip_inherited_keys mode: {mode}")
    for key, value in child_cfg.items():
        if key in base_cfg:
            if mode == "present":
                continue
            if base_cfg.get(key) == value:
                continue
        trimmed[key] = value
    return trimmed


def format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize geometry configuration for CLI echo output.
    """
    g = dict(geom_cfg)
    freeze_atoms = g.get("freeze_atoms")
    if freeze_atoms is None:
        return g

    if isinstance(freeze_atoms, str):
        return g

    try:
        items = list(freeze_atoms)
    except TypeError:
        return g

    joined = ",".join(map(str, items))
    g["freeze_atoms"] = f"[{joined}]" if items else "[]"
    return g


def format_elapsed(prefix: str, start_time: float, end_time: Optional[float] = None) -> str:
    """Return a formatted elapsed-time string with the provided ``prefix`` label."""
    elapsed = max(0.0, (end_time or time.perf_counter()) - start_time)
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{prefix}: {int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"


def xyz_string_with_energy(geom: Any, energy: Optional[float] = None) -> str:
    """Return an XYZ string, optionally overwriting the comment line with an energy value."""
    s = geom.as_xyz()
    lines = s.splitlines()
    if energy is not None and len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{energy:.12f}"
        s = "\n".join(lines)
    if not s.endswith("\n"):
        s += "\n"
    return s


def distance_A_from_coords(coords_bohr: "np.ndarray", i: int, j: int) -> float:
    """Return interatomic distance in Å given coords in Bohr."""
    diff = coords_bohr[i] - coords_bohr[j]
    return float(np.linalg.norm(diff) / ANG2BOHR)


def distance_tag(value_A: float, *, digits: int = 2, pad: int = 3) -> str:
    """Format a distance in Å as a zero-padded integer tag (default: ×10^2)."""
    scale = 10 ** digits
    return f"{int(round(value_A * scale)):0{pad}d}"


def as_list(raw: Any) -> List[Any]:
    """Return ``raw`` as a list, or [] when not iterable/None."""
    if raw is None:
        return []
    try:
        return list(raw)
    except Exception:
        return []


def ensure_dir(path: Path) -> None:
    """Create a directory (parents ok); noop if it already exists."""
    path.mkdir(parents=True, exist_ok=True)


def collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
    """
    Collect values following a flag that may appear once with multiple space-separated values,
    e.g., "-i A B C".
    """
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


def collect_single_option_values(
    argv: Sequence[str],
    names: Sequence[str],
    label: str,
) -> List[str]:
    """Collect values following a flag that must appear at most once."""
    vals: List[str] = []
    seen = 0
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            seen += 1
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    if seen > 1:
        raise click.BadParameter(
            f"Use a single {label} followed by multiple values; repeated flags are not accepted."
        )
    return vals


def geom_from_xyz_string(
    xyz_text: str,
    *,
    coord_type: str,
    freeze_atoms: Optional[Sequence[int]] = None,
) -> Any:
    """Load a pysisyphus Geometry from an XYZ text string (tempfile-backed)."""
    s = xyz_text if xyz_text.endswith("\n") else (xyz_text + "\n")
    freeze_atoms = list(freeze_atoms) if freeze_atoms is not None else []
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()

        g = geom_loader(
            Path(tmp.name),
            coord_type=coord_type,
            freeze_atoms=freeze_atoms,
        )
        try:
            g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        except Exception:
            click.echo(
                "[geom] WARNING: Failed to attach freeze_atoms to geometry.",
                err=True,
            )
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def snapshot_geometry(geom: Any, *, coord_type_default: str) -> Any:
    """Create an independent pysisyphus Geometry snapshot from the given Geometry."""
    s = geom.as_xyz()
    return geom_from_xyz_string(
        s,
        coord_type=getattr(geom, "coord_type", coord_type_default),
        freeze_atoms=getattr(geom, "freeze_atoms", []),
    )


def make_snapshot_geometry(coord_type_default: str) -> Callable[[Any], Any]:
    """Return a snapshot helper bound to a default coord_type (scan helpers)."""
    return functools.partial(snapshot_geometry, coord_type_default=coord_type_default)


def normalize_freeze_atoms(raw: Any) -> List[int]:
    """Normalize freeze_atoms values (string/list/iterable) into a list of integers."""
    if raw is None:
        return []
    if isinstance(raw, str):
        tokens = re.findall(r"-?\d+", raw)
        return [int(tok) for tok in tokens]
    try:
        return [int(i) for i in raw]
    except Exception:
        return []


def merge_freeze_atom_indices(
    geom_cfg: Dict[str, Any],
    *indices: _Iterable[int],
) -> List[int]:
    """Merge one or more iterables of indices into ``geom_cfg['freeze_atoms']``.

    Existing entries are preserved, duplicates removed, and the result sorted.
    The updated list is returned.
    """
    merged: set[int] = set()
    base = geom_cfg.get("freeze_atoms", None)
    merged.update(normalize_freeze_atoms(base))
    for seq in indices:
        merged.update(normalize_freeze_atoms(seq))
    result = sorted(merged)
    geom_cfg["freeze_atoms"] = result
    return result


def merge_freeze_atom_groups(*groups: Sequence[int]) -> List[int]:
    """Merge multiple freeze_atoms groups into a sorted list of ints."""
    merged: set[int] = set()
    for group in groups:
        merged.update(normalize_freeze_atoms(group))
    return sorted(merged)


def build_sopt_kwargs(
    kind: str,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    max_step_bohr: float,
    relax_max_cycles: int,
    relax_override_requested: bool,
    out_dir: Path,
    prefix: str,
) -> Dict[str, Any]:
    """Build LBFGS/RFO optimizer kwargs with a shared max-step cap."""
    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    if kind == "lbfgs":
        args = {**lbfgs_cfg, **common}
        args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
    else:
        args = {**rfo_cfg, **common}
        tr = float(rfo_cfg.get("trust_radius", 0.10))
        args["trust_radius"] = min(tr, max_step_bohr)
        args["trust_max"] = min(float(rfo_cfg.get("trust_max", 0.10)), max_step_bohr)
    if relax_override_requested:
        args["max_cycles"] = int(relax_max_cycles)
    return args


def make_sopt_optimizer(
    geom: Any,
    kind: str,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    max_step_bohr: float,
    relax_max_cycles: int,
    relax_override_requested: bool,
    out_dir: Path,
    prefix: str,
):
    """Construct an LBFGS/RFO optimizer based on shared settings."""
    args = build_sopt_kwargs(
        kind,
        lbfgs_cfg,
        rfo_cfg,
        opt_cfg,
        max_step_bohr,
        relax_max_cycles,
        relax_override_requested,
        out_dir,
        prefix,
    )
    from pysisyphus.optimizers.LBFGS import LBFGS
    from pysisyphus.optimizers.RFOptimizer import RFOptimizer

    if kind == "lbfgs":
        return LBFGS(geom, **args)
    return RFOptimizer(geom, **args)


def normalize_choice(
    value: str,
    *,
    param: str,
    alias_groups: _Sequence[tuple[_Sequence[str], str]],
    allowed_hint: str,
) -> str:
    """Normalize *value* using alias groups and raise ``click.BadParameter`` on failure."""
    key = (value or "").strip().lower()
    for aliases, canonical in alias_groups:
        if any(key == alias.lower() for alias in aliases):
            return canonical

    hint = allowed_hint.strip()
    detail = f" Allowed: {hint}." if hint else ""
    raise click.BadParameter(f"Unknown value for {param} '{value}'.{detail}")


def cli_param_overridden(ctx: click.Context, name: str) -> bool:
    """Return True when a CLI parameter value was explicitly provided."""
    try:
        source = ctx.get_parameter_source(name)
    except Exception:
        return True
    return source not in (None, ParameterSource.DEFAULT)


def deep_update(dst: Dict[str, Any], src: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Recursively update mapping *dst* with *src*, returning *dst*.
    """
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _get_mapping_section(cfg: Mapping[str, Any], path: _Sequence[str]) -> Optional[Dict[str, Any]]:
    cur: Any = cfg
    for key in path:
        if not isinstance(cur, Mapping):
            return None
        cur = cur.get(key)
        if cur is None:
            return None
    return cur if isinstance(cur, dict) else None


def apply_yaml_overrides(
    yaml_cfg: Mapping[str, Any],
    overrides: _Sequence[Tuple[Dict[str, Any], _Sequence[_Sequence[str]]]],
) -> None:
    """Apply YAML overrides to multiple target dictionaries.

    Parameters
    ----------
    yaml_cfg : Mapping[str, Any]
        Parsed YAML configuration (root-level mapping).
    overrides : Sequence[Tuple[Dict[str, Any], Sequence[Sequence[str]]]]
        Each entry consists of the target dictionary to update followed by one or
        more candidate key paths. The first existing path is used. For example::

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (lbfgs_cfg, (("sopt", "lbfgs"), ("lbfgs",))),
                ],
            )

        This mirrors the previous ``deep_update(..., yaml_cfg.get(...))`` pattern
        while centralizing the shared logic.
    """
    for target, paths in overrides:
        for path in paths:
            norm_path = tuple(path)
            section = _get_mapping_section(yaml_cfg, norm_path)
            if section is not None:
                deep_update(target, section)
                break


def load_yaml_dict(path: Optional[Path]) -> Dict[str, Any]:
    """
    Load a YAML file whose root must be a mapping. Return an empty dict if *path* is None.
    """
    if not path:
        return {}

    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")

    return data


def build_scan_configs(
    yaml_cfg: Mapping[str, Any],
    *,
    geom_kw: Dict[str, Any],
    calc_kw: Dict[str, Any],
    opt_kw: Dict[str, Any],
    lbfgs_kw: Dict[str, Any],
    rfo_kw: Dict[str, Any],
    bias_kw: Dict[str, Any],
    extra_overrides: Sequence[Tuple[Dict[str, Any], _Sequence[_Sequence[str]]]] = (),
    charge: Optional[int] = None,
    spin: Optional[int] = None,
    workers: int = 1,
    workers_per_node: int = 1,
    out_dir: str = ".",
    thresh: Optional[str] = None,
    bias_k: Optional[float] = None,
    set_charge_spin: bool = True,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Build common scan configs (defaults ← CLI ← YAML)."""
    geom_cfg = dict(geom_kw)
    calc_cfg = dict(calc_kw)
    opt_cfg = dict(opt_kw)
    lbfgs_cfg = dict(lbfgs_kw)
    rfo_cfg = dict(rfo_kw)
    bias_cfg = dict(bias_kw)

    if set_charge_spin:
        if charge is not None:
            calc_cfg["charge"] = int(charge)
        if spin is not None:
            calc_cfg["spin"] = int(spin)
    calc_cfg["workers"] = int(workers)
    calc_cfg["workers_per_node"] = int(workers_per_node)
    opt_cfg["out_dir"] = out_dir
    opt_cfg["dump"] = False
    if thresh is not None:
        opt_cfg["thresh"] = str(thresh)

    apply_yaml_overrides(
        yaml_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (opt_cfg, (("opt",),)),
            (lbfgs_cfg, (("lbfgs",),)),
            (rfo_cfg, (("rfo",),)),
            (bias_cfg, (("bias",),)),
            *list(extra_overrides),
        ],
    )

    if bias_k is not None:
        bias_cfg["k"] = float(bias_k)

    return geom_cfg, calc_cfg, opt_cfg, lbfgs_cfg, rfo_cfg, bias_cfg


def convert_xyz_to_gjf_if_enabled(
    xyz_path: Path,
    template: Optional["GjfTemplate"],
    *,
    out_path: Optional[Path] = None,
    context: str = "GJF",
    on_error: str = "raise",
) -> Optional[Path]:
    """Convert XYZ to GJF when enabled; return output path or None."""
    if not (_CONVERT_FILES_ENABLED and template is not None and xyz_path.exists()):
        return None
    target = out_path if out_path is not None else xyz_path.with_suffix(".gjf")
    try:
        convert_xyz_to_gjf(xyz_path, template, target)
        return target
    except Exception as e:
        if on_error == "warn":
            click.echo(
                f"[convert] WARNING: Failed to convert '{xyz_path.name}' to {context}: {e}",
                err=True,
            )
            return None
        raise click.ClickException(
            f"[convert] Failed to convert '{xyz_path.name}' to {context}: {e}"
        ) from e


# =============================================================================
# Plotly: Energy diagram builder
# =============================================================================
def build_energy_diagram(
    energies: Sequence[float],
    labels: Sequence[str],
    ylabel: str = "ΔE",
    baseline: bool = False,
    showgrid: bool = False,
) -> go.Figure:
    """
    Plot an energy diagram using Plotly.

    Parameters
    ----------
    energies : Sequence[float]
        Energies for each state (same unit). Values are plotted without conversion.
    labels : Sequence[str]
        Labels corresponding to each state (for example, ["R", "TS1", "IM1", "TS2", "P"]).
        Must be the same length as ``energies``.
    ylabel : str, optional
        Y-axis label (for example, "ΔE" or "ΔG"). Defaults to ``"ΔE"``.
    baseline : bool, optional
        If ``True``, draw a dotted baseline at the energy of the first state across the plot.
    showgrid : bool, optional
        If ``True``, show grid lines on both axes. Defaults to ``False``.

    Returns
    -------
    plotly.graph_objs.Figure
        Figure containing the energy diagram.

    Notes
    -----
    - Each state is rendered as a thick horizontal segment (width ``HLINE_WIDTH``).
    - Adjacent states are connected by dotted diagonal segments from the right end of
      the left state to the left end of the right state.
    - Segment length automatically shrinks with additional states so that gaps remain
      between neighbors.
    - X-axis ticks are centered on each state and labeled using ``labels``.
    """
    if len(energies) == 0:
        raise ValueError("`energies` must contain at least one value.")
    if len(energies) != len(labels):
        raise ValueError("`energies` and `labels` must have the same length.")

    n = len(energies)
    energies = [float(e) for e in energies]

    # -----------------------------
    # Layout/style constants
    # -----------------------------
    AXIS_WIDTH = 3
    FONT_SIZE = 18
    AXIS_TITLE_SIZE = 20
    HLINE_WIDTH = 6           # Width of the horizontal state segments
    CONNECTOR_WIDTH = 2       # Width of the dotted connectors
    LINE_COLOR = "#1C1C1C"
    GRID_COLOR = "lightgrey"

    # -----------------------------
    # Geometry along the X axis (centers and segment lengths)
    # -----------------------------
    # Place segment centers at 0.5, 1.5, 2.5, ... (equally spaced)
    centers = [i + 0.5 for i in range(n)]

    # Shorten the segment as n grows (min 0.35, max 0.85)
    # Examples: n=5 -> 0.7, n=10 -> 0.5, n>=20 -> 0.35
    seg_width = min(0.85, max(0.35, 0.90 - 0.04 * n))
    half = seg_width / 2.0

    lefts = [c - half for c in centers]
    rights = [c + half for c in centers]

    # -----------------------------
    # Assemble the figure
    # -----------------------------
    fig = go.Figure()

    # Baseline (dotted line at the first energy level)
    if baseline:
        fig.add_trace(
            go.Scatter(
                x=[lefts[0], rights[-1]],
                y=[energies[0], energies[0]],
                mode="lines",
                line=dict(color=GRID_COLOR, dash="dot", width=2),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Horizontal segments for each state
    for i, (e, lab) in enumerate(zip(energies, labels)):
        fig.add_trace(
            go.Scatter(
                x=[lefts[i], rights[i]],
                y=[e, e],
                mode="lines",
                line=dict(color=LINE_COLOR, width=HLINE_WIDTH),
                hovertemplate=f"{lab}: %{{y:.6f}}<extra></extra>",
                showlegend=False,
            )
        )

    # Dotted diagonals between adjacent states (right end -> left end)
    for i in range(n - 1):
        fig.add_trace(
            go.Scatter(
                x=[rights[i], lefts[i + 1]],
                y=[energies[i], energies[i + 1]],
                mode="lines",
                line=dict(color=LINE_COLOR, width=CONNECTOR_WIDTH, dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # -----------------------------
    # Axis ranges and styling
    # -----------------------------
    # Add a small margin beyond the first/last segments on X
    xpad = max(0.08, 0.15 * (1.0 - seg_width))
    x_min = lefts[0] - xpad
    x_max = rights[-1] + xpad

    # Add vertical padding above and below
    y_min = min(energies)
    y_max = max(energies)
    span = max(1e-6, y_max - y_min)  # Avoid zero span even if all values match
    ypad_low = 0.10 * span
    ypad_high = 0.20 * span
    y_range = [y_min - ypad_low, y_max + ypad_high]

    xaxis_config = dict(
        range=[x_min, x_max],
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        tickmode="array",
        tickvals=centers,
        ticktext=list(labels),
        title=dict(text="", font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    yaxis_config = dict(
        range=y_range,
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        title=dict(text=ylabel, font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    fig.update_layout(
        xaxis=xaxis_config,
        yaxis=yaxis_config,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=80, r=40, t=40, b=80),
    )

    return fig


# =============================================================================
# Coordinate conversion utilities
# =============================================================================
def convert_xyz_to_pdb(xyz_path: Path, ref_pdb_path: Path, out_pdb_path: Path) -> None:
    """Overlay coordinates from *xyz_path* onto the topology of *ref_pdb_path* and write to *out_pdb_path*.

    Notes
    -----
        - *xyz_path* may contain one or many frames. For multi‑frame trajectories,
          a MODEL/ENDMDL block is appended for each subsequent frame in the output PDB.
        - On the first frame the output file is created/overwritten; subsequent frames are appended.

    Args:
        xyz_path: Path to an XYZ file (single or multi-frame).
        ref_pdb_path: Path to a reference PDB providing atom ordering/topology.
        out_pdb_path: Destination PDB file to write.
    """
    ref_atoms = read(ref_pdb_path)  # Reference topology/ordering (single frame)
    traj = read(xyz_path, index=":", format="xyz")  # Load all frames from the XYZ
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")

    ref_symbols = ref_atoms.get_chemical_symbols()

    for step, frame in enumerate(traj):
        xyz_symbols = frame.get_chemical_symbols()
        xyz_positions = frame.get_positions()

        if xyz_symbols != ref_symbols:
            raise ValueError(
                "Atom ordering mismatch between XYZ and PDB; "
                "expected identical ordering when converting coordinates."
            )

        atoms = ref_atoms.copy()
        atoms.set_positions(xyz_positions)
        if step == 0:
            write(out_pdb_path, atoms)  # Create/overwrite on the first frame
        else:
            write(out_pdb_path, atoms, append=True)  # Append subsequent frames using MODEL/ENDMDL


# =============================================================================
# Gaussian input (.gjf) helpers
# =============================================================================

_FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
_GJF_COORD_RE = re.compile(
    rf"^(?P<prefix>.*?)(?P<sep0>\s*)(?P<x>{_FLOAT_PATTERN})(?P<sep1>\s+)"
    rf"(?P<y>{_FLOAT_PATTERN})(?P<sep2>\s+)(?P<z>{_FLOAT_PATTERN})(?P<suffix>\s*)$"
)


def _count_decimals(token: str) -> int:
    if "." not in token:
        return 6
    mantissa = token.split(".", 1)[1]
    mantissa = mantissa.split("e", 1)[0].split("E", 1)[0]
    mantissa = mantissa.split("d", 1)[0].split("D", 1)[0]
    return sum(ch.isdigit() for ch in mantissa)


def _format_like(template: str, value: float) -> str:
    stripped = template.strip()
    width = len(template)
    if not stripped:
        formatted = f"{value:.10f}"
    elif any(ch in stripped for ch in "eEdD"):
        mantissa = stripped.replace("d", "e").replace("D", "E")
        mantissa_only, _, _ = mantissa.partition("E")
        decimals = _count_decimals(mantissa_only)
        formatted = f"{value:.{decimals}E}"
        if "d" in template:
            formatted = formatted.replace("E", "d")
        elif "D" in template:
            formatted = formatted.replace("E", "D")
        elif "e" in template:
            formatted = formatted.replace("E", "e")
    else:
        decimals = _count_decimals(stripped)
        formatted = f"{value:.{decimals}f}"
    if len(formatted) < width:
        formatted = formatted.rjust(width)
    return formatted


def _token_to_symbol(text: str) -> str:
    stripped = text.strip()
    if not stripped:
        raise ValueError("Empty atom token in Gaussian coordinate line.")
    token = stripped.split()[0]
    m = re.match(r"([A-Za-z]{1,2})", token)
    if m:
        return m.group(1).title()
    if token.isdigit():
        z = int(token)
        if 0 <= z < len(chemical_symbols):
            return chemical_symbols[z]
    raise ValueError(f"Could not determine element symbol from '{token}'.")


def _convert_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


@dataclass
class GjfCoordinateLine:
    prefix: str
    sep0: str
    sep1: str
    sep2: str
    suffix: str
    x_template: str
    y_template: str
    z_template: str
    symbol: str
    x: float
    y: float
    z: float

    def render(self, coords: Tuple[float, float, float]) -> str:
        x_str = _format_like(self.x_template, coords[0])
        y_str = _format_like(self.y_template, coords[1])
        z_str = _format_like(self.z_template, coords[2])
        return f"{self.prefix}{self.sep0}{x_str}{self.sep1}{y_str}{self.sep2}{z_str}{self.suffix}"


@dataclass
class GjfTemplate:
    path: Path
    prefix_lines: List[str]
    suffix_lines: List[str]
    coord_lines: List[GjfCoordinateLine]
    charge: int
    spin: int

    @property
    def natoms(self) -> int:
        return len(self.coord_lines)

    def as_xyz_string(self) -> str:
        lines = [str(self.natoms), f"converted from {self.path.name}"]
        for atom in self.coord_lines:
            lines.append(f"{atom.symbol}  {atom.x:.10f}  {atom.y:.10f}  {atom.z:.10f}")
        return "\n".join(lines) + "\n"


def write_xyz_trj_with_energy(images: Sequence[Any], energies: Sequence[float], path: Path) -> None:
    """Write an XYZ `.trj` with the energy on line 2 of each block."""
    blocks: List[str] = []
    E = np.array(energies, dtype=float)
    for geom, e in zip(images, E):
        if hasattr(geom, "as_xyz"):
            blocks.append(xyz_string_with_energy(geom, energy=float(e)))
            continue
        # ASE Atoms fallback
        symbols = geom.get_chemical_symbols()
        coords = geom.get_positions()
        lines = [str(len(symbols)), f"{float(e):.12f}"]
        lines.extend(
            f"{sym} {x:.15f} {y:.15f} {z:.15f}"
            for sym, (x, y, z) in zip(symbols, coords)
        )
        blocks.append("\n".join(lines) + "\n")
    with open(path, "w") as f:
        f.write("".join(blocks))


def set_freeze_atoms_or_warn(
    geom: Any,
    freeze_atoms: Sequence[int],
    *,
    context: str,
) -> None:
    """Attach freeze_atoms to a geometry; warn once on failure."""
    if not freeze_atoms:
        return
    try:
        geom.freeze_atoms = np.array(sorted({int(i) for i in freeze_atoms}), dtype=int)
    except Exception:
        click.echo(f"[{context}] WARNING: Failed to attach freeze_atoms to geometry.", err=True)


def read_xyz_energies(path: Path | str) -> List[float]:
    """
    Extract energies from the second-line comment of each XYZ frame.
    The first numeric token found on the comment line is used.
    """
    energies: List[float] = []
    with open(path, encoding="utf-8") as fh:
        while (hdr := fh.readline()):
            try:
                nat = int(hdr.strip())
            except ValueError:
                break
            comment = fh.readline().strip()
            m = re.search(
                r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
                comment,
            )
            if not m:
                raise RuntimeError(f"Energy not found in comment: {comment}")
            energies.append(float(m.group(1)))
            for _ in range(nat):
                fh.readline()
    if not energies:
        raise RuntimeError(f"No energy data in {path}")
    return energies


def parse_xyz_block(
    block: Sequence[str],
    *,
    path: Path,
    frame_idx: int,
) -> Tuple[List[str], np.ndarray]:
    if not block:
        raise click.ClickException(f"[xyz] Empty XYZ frame in {path}")
    try:
        nat = int(block[0].strip().split()[0])
    except Exception:
        raise click.ClickException(
            f"[xyz] Malformed XYZ/TRJ header in frame {frame_idx} of {path}"
        )
    if len(block) < 2 + nat:
        raise click.ClickException(
            f"[xyz] Incomplete XYZ frame {frame_idx} in {path} (expected {nat} atoms)."
        )
    elems: List[str] = []
    coords: List[List[float]] = []
    for k in range(nat):
        parts = block[2 + k].split()
        if len(parts) < 4:
            raise click.ClickException(
                f"[xyz] Malformed atom line in frame {frame_idx} of {path}"
            )
        elems.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords, dtype=float)


def xyz_blocks_first_last(
    blocks: Sequence[Sequence[str]],
    *,
    path: Path,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    if not blocks:
        raise click.ClickException(f"[xyz] No frames found in {path}")
    first_elems, first_coords = parse_xyz_block(blocks[0], path=path, frame_idx=1)
    last_elems, last_coords = parse_xyz_block(blocks[-1], path=path, frame_idx=len(blocks))
    if first_elems != last_elems:
        raise click.ClickException(f"[xyz] Element list changed across frames in {path}")
    return first_elems, first_coords, last_coords


def read_xyz_first_last(trj_path: Path) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """
    Lightweight XYZ trajectory reader: return (elements, first_coords[Å], last_coords[Å]).
    Assumes standard multi-frame XYZ: natoms line, comment line, natoms atom lines.
    """
    blocks = read_xyz_as_blocks(trj_path, strict=True)
    return xyz_blocks_first_last(blocks, path=trj_path)


def read_xyz_as_blocks(trj_path: Path, *, strict: bool = False) -> List[List[str]]:
    """
    Read a multi-frame XYZ/TRJ file and return a list of frames, each as a list of lines.

    When *strict* is True, malformed headers or truncated frames raise a ClickException.
    """
    try:
        lines = trj_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as e:
        raise click.ClickException(f"[xyz] Failed to read XYZ/TRJ: {trj_path} ({e})")

    blocks: List[List[str]] = []
    i = 0
    n = len(lines)
    while i < n:
        while i < n and not lines[i].strip():
            i += 1
        if i >= n:
            break
        header = lines[i].strip()
        try:
            nat = int(header.split()[0])
        except Exception:
            if strict:
                raise click.ClickException(f"[xyz] Malformed header at line {i+1} in {trj_path}")
            break
        block = lines[i : i + 2 + nat]
        if len(block) < 2 + nat:
            if strict:
                raise click.ClickException(
                    f"[xyz] Incomplete frame at line {i+1} in {trj_path} (expected {nat} atoms)."
                )
            break
        blocks.append(block)
        i += 2 + nat
    return blocks


@dataclass
class PreparedInputStructure:
    source_path: Path
    geom_path: Path
    gjf_template: Optional[GjfTemplate] = None
    _tmp_geom_path: Optional[Path] = None

    @property
    def is_gjf(self) -> bool:
        return self.gjf_template is not None

    def cleanup(self) -> None:
        if self._tmp_geom_path and self._tmp_geom_path.exists():
            try:
                self._tmp_geom_path.unlink()
            except FileNotFoundError:
                pass

    def __enter__(self) -> "PreparedInputStructure":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.cleanup()


def _parse_coord_line(line: str) -> GjfCoordinateLine:
    match = _GJF_COORD_RE.match(line)
    if not match:
        raise ValueError(f"Could not parse Gaussian coordinate line: '{line}'.")
    prefix = match.group("prefix")
    sep0 = match.group("sep0")
    sep1 = match.group("sep1")
    sep2 = match.group("sep2")
    suffix = match.group("suffix")
    x_token = match.group("x")
    y_token = match.group("y")
    z_token = match.group("z")
    symbol = _token_to_symbol(f"{prefix}{sep0}")
    return GjfCoordinateLine(
        prefix=prefix,
        sep0=sep0,
        sep1=sep1,
        sep2=sep2,
        suffix=suffix,
        x_template=x_token,
        y_template=y_token,
        z_template=z_token,
        symbol=symbol,
        x=_convert_float(x_token),
        y=_convert_float(y_token),
        z=_convert_float(z_token),
    )


def parse_gjf_template(path: Path) -> GjfTemplate:
    lines = path.read_text().splitlines()
    section = 0
    charge_line_idx = None
    charge = None
    spin = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            if section < 2:
                section += 1
            continue
        if section < 2:
            continue
        tokens = stripped.split()
        if len(tokens) < 2:
            continue
        try:
            charge = int(tokens[0])
            spin = int(tokens[1])
        except ValueError:
            continue
        charge_line_idx = idx
        break
    if charge_line_idx is None or charge is None or spin is None:
        raise ValueError(f"Failed to locate charge/multiplicity line in '{path}'.")

    coord_start = charge_line_idx + 1
    while coord_start < len(lines) and not lines[coord_start].strip():
        coord_start += 1

    prefix_lines = lines[:coord_start]
    coord_lines_raw: List[str] = []
    idx = coord_start
    while idx < len(lines):
        line = lines[idx]
        if not line.strip():
            break
        coord_lines_raw.append(line)
        idx += 1
    if not coord_lines_raw:
        raise ValueError(f"No coordinates found in '{path}'.")
    suffix_lines = lines[idx:]

    coord_lines = [_parse_coord_line(line) for line in coord_lines_raw]
    return GjfTemplate(
        path=path,
        prefix_lines=prefix_lines,
        suffix_lines=suffix_lines,
        coord_lines=coord_lines,
        charge=charge,
        spin=spin,
    )


def prepare_input_structure(path: Path) -> PreparedInputStructure:
    if path.suffix.lower() != ".gjf":
        return PreparedInputStructure(source_path=path, geom_path=path)
    template = parse_gjf_template(path)
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(template.as_xyz_string())
        tmp.flush()
    finally:
        tmp.close()
    tmp_path = Path(tmp.name)
    return PreparedInputStructure(
        source_path=path,
        geom_path=tmp_path,
        gjf_template=template,
        _tmp_geom_path=tmp_path,
    )


def _read_atom_count(path: Path) -> int:
    try:
        atoms = read(path, index=0)
    except Exception as e:
        raise click.ClickException(f"Failed to read '{path}' for --ref-pdb validation: {e}")
    return len(atoms)


def _validate_ref_pdb_atom_count(geom_path: Path, ref_pdb_path: Path) -> None:
    geom_count = _read_atom_count(geom_path)
    ref_count = _read_atom_count(ref_pdb_path)
    if geom_count != ref_count:
        raise click.ClickException(
            f"--ref-pdb atom count ({ref_count}) does not match input ({geom_count})."
        )


def apply_ref_pdb_override(
    prepared_input: PreparedInputStructure,
    ref_pdb: Optional[Path],
) -> Optional[Path]:
    """Use a reference PDB topology while keeping XYZ coordinates for geometry loading."""
    if ref_pdb is None:
        return None
    ref_pdb = Path(ref_pdb).resolve()
    if ref_pdb.suffix.lower() != ".pdb":
        raise click.BadParameter("--ref-pdb must be a .pdb file.")
    _validate_ref_pdb_atom_count(prepared_input.geom_path, ref_pdb)
    prepared_input.source_path = ref_pdb
    return ref_pdb


def fill_charge_spin_from_gjf(
    charge: Optional[int],
    spin: Optional[int],
    template: Optional[GjfTemplate],
) -> Tuple[Optional[int], Optional[int]]:
    if template is not None:
        if charge is None:
            charge = template.charge
        if spin is None:
            spin = template.spin
    return charge, spin


def _round_charge_with_note(q: float, prefix: str) -> int:
    if not math.isfinite(q):
        raise click.BadParameter(f"{prefix} Computed total charge is non-finite: {q!r}")
    q_int = int(round(q))
    if not math.isclose(q_int, q):
        click.echo(
            f"{prefix} NOTE: extractor total charge = {q:+g} → rounded to integer {q_int:+d}."
        )
    return q_int


def _derive_charge_from_ligand_charge(
    prepared: PreparedInputStructure,
    ligand_charge: Optional[float | str | Dict[str, float]],
    *,
    prefix: str,
) -> Optional[int]:
    if ligand_charge is None:
        return None
    try:
        from Bio import PDB

        from .extract import compute_charge_summary, log_charge_summary

        parser = PDB.PDBParser(QUIET=True)
        complex_struct = parser.get_structure("complex", str(prepared.source_path))
        selected_ids = {res.get_full_id() for res in complex_struct.get_residues()}
        summary = compute_charge_summary(complex_struct, selected_ids, set(), ligand_charge)
        log_charge_summary(prefix, summary)
        q_total = float(summary.get("total_charge", 0.0))
        click.echo(
            f"{prefix} Charge summary from full complex (--ligand-charge without extraction):"
        )
        click.echo(
            f"  Protein: {summary.get('protein_charge', 0.0):+g},  "
            f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
            f"Ions: {summary.get('ion_total_charge', 0.0):+g},  "
            f"Total: {q_total:+g}"
        )
        return _round_charge_with_note(q_total, prefix)
    except Exception as e:
        click.echo(
            f"{prefix} NOTE: failed to derive total charge from --ligand-charge: {e}",
            err=True,
        )
        return None


def resolve_charge_spin(
    prepared_inputs: PreparedInputStructure | Sequence[PreparedInputStructure],
    charge: Optional[int],
    spin: Optional[int],
    *,
    spin_default: int = 1,
    charge_default: int = 0,
    ligand_charge: Optional[float | str | Dict[str, float]] = None,
    prefix: str = "[charge]",
    cleanup_on_error: Optional[Callable[[], None]] = None,
) -> Tuple[int, int]:
    """Resolve charge/spin from CLI args, GJF templates, and ligand metadata.

    Accepts either a single PreparedInputStructure or a sequence of them.
    """
    # Normalize to sequence
    if isinstance(prepared_inputs, PreparedInputStructure):
        inputs = [prepared_inputs]
        cleanup_on_error = cleanup_on_error or prepared_inputs.cleanup
    else:
        inputs = list(prepared_inputs)

    resolved_charge = charge
    resolved_spin = spin
    for prepared in inputs:
        resolved_charge, resolved_spin = fill_charge_spin_from_gjf(
            resolved_charge, resolved_spin, prepared.gjf_template
        )

    if ligand_charge is not None:
        for prepared in inputs:
            if prepared.source_path.suffix.lower() in {".xyz", ".gjf"}:
                if cleanup_on_error:
                    cleanup_on_error()
                raise click.ClickException(
                    "--ligand-charge is only supported for PDB inputs; it cannot be used with .xyz or .gjf files."
                )
        if resolved_charge is None:
            resolved_charge = _derive_charge_from_ligand_charge(
                inputs[0], ligand_charge, prefix=prefix
            )

    if resolved_charge is None:
        if any(not p.is_gjf for p in inputs):
            if cleanup_on_error:
                cleanup_on_error()
            raise click.ClickException(
                "-q/--charge is required unless the input is a .gjf template with charge metadata."
            )
        resolved_charge = charge_default

    if resolved_spin is None:
        resolved_spin = spin_default
    return int(resolved_charge), int(resolved_spin)


# Backwards compatibility aliases
resolve_charge_spin_or_raise = resolve_charge_spin
resolve_charge_spin_multi = resolve_charge_spin


@contextmanager
def prepared_cli_input(
    input_path: Path,
    *,
    ref_pdb: Optional[Path],
    charge: Optional[int],
    spin: Optional[int],
    ligand_charge: Optional[float | str | Dict[str, float]] = None,
    prefix: str = "[charge]",
) -> Iterator[Tuple[PreparedInputStructure, int, int]]:
    """Context-managed input preparation with charge/spin resolution."""
    with prepare_input_structure(input_path) as prepared:
        apply_ref_pdb_override(prepared, ref_pdb)
        charge_res, spin_res = resolve_charge_spin(
            prepared,
            charge,
            spin,
            ligand_charge=ligand_charge,
            prefix=prefix,
        )
        yield prepared, charge_res, spin_res


_CONVERT_FILES_ENABLED: bool = True


def set_convert_file_enabled(enabled: bool) -> None:
    """Globally enable or disable XYZ/TRJ conversions to PDB/GJF outputs.

    The toggle mirrors the ``--convert-files {True|False}`` CLI flag used
    by every workflow. When disabled, format-aware conversions are skipped even
    if reference templates are available.
    """

    global _CONVERT_FILES_ENABLED
    _CONVERT_FILES_ENABLED = bool(enabled)


def convert_xyz_to_gjf(xyz_path: Path, template: GjfTemplate, out_path: Path) -> None:
    """Render single- or multi-frame XYZ/TRJ coordinates into a Gaussian template.

    Multi-frame trajectories are emitted as blank-separated geometries suitable
    for QST-style Gaussian inputs.
    """
    traj = read(xyz_path, index=":", format="xyz")
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")
    new_lines: List[str] = list(template.prefix_lines)
    for frame_idx, atoms in enumerate(traj):
        if len(atoms) != template.natoms:
            raise ValueError(
                f"Atom count mismatch for '{xyz_path}': xyz has {len(atoms)} atoms, "
                f"but template has {template.natoms}."
            )
        coords = atoms.get_positions()
        if frame_idx > 0:
            new_lines.append("")  # Blank line between multiple geometries (QST-style)
        for idx, coord_line in enumerate(template.coord_lines):
            new_lines.append(coord_line.render(tuple(map(float, coords[idx]))))
    new_lines.extend(template.suffix_lines)
    text = "\n".join(new_lines)
    if not text.endswith("\n"):
        text += "\n"
    out_path.write_text(text)


def convert_xyz_like_outputs(
    xyz_path: Path,
    prepared_input: PreparedInputStructure,
    *,
    ref_pdb_path: Optional[Path],
    out_pdb_path: Optional[Path] = None,
    out_gjf_path: Optional[Path] = None,
    context: str = "outputs",
    on_error: str = "raise",
) -> bool:
    """Convert an XYZ/TRJ file to PDB outputs (and XYZ to GJF) based on the original input type.

    Parameters
    ----------
    xyz_path:
        Source XYZ-like file (single or multi-frame).
    prepared_input:
        Prepared input structure returned by :func:`prepare_input_structure`.
    ref_pdb_path:
        Reference PDB topology (required for PDB conversions).
    out_pdb_path / out_gjf_path:
        Targets for the converted files. Conversions are skipped when the
        corresponding output path is ``None`` or when the input type does not
        request that format.
    Returns True when at least one conversion was attempted and succeeded; False otherwise.
    """

    if not _CONVERT_FILES_ENABLED:
        return False

    source_suffix = prepared_input.source_path.suffix.lower()
    needs_pdb = source_suffix == ".pdb" and out_pdb_path is not None and ref_pdb_path is not None
    needs_gjf = (
        xyz_path.suffix.lower() == ".xyz"
        and prepared_input.is_gjf
        and prepared_input.gjf_template is not None
        and out_gjf_path is not None
    )

    if not (needs_pdb or needs_gjf):
        return False

    try:
        if needs_pdb:
            convert_xyz_to_pdb(xyz_path, ref_pdb_path, out_pdb_path)
        if needs_gjf:
            convert_xyz_to_gjf(xyz_path, prepared_input.gjf_template, out_gjf_path)
    except Exception as e:
        if on_error == "warn":
            click.echo(f"[convert] WARNING: Failed to convert {context}: {e}", err=True)
            return False
        raise click.ClickException(f"[convert] Failed to convert {context}: {e}") from e
    return True


def _convert_to_pdb_logged(
    src_path: Path, ref_pdb_path: Optional[Path], out_path: Optional[Path] = None
) -> Optional[Path]:
    """Convert an XYZ/TRJ to PDB when conversion is enabled; return path or None."""
    try:
        if ref_pdb_path is None or not _CONVERT_FILES_ENABLED:
            return None
        src_path = Path(src_path)
        if (not src_path.exists()) or src_path.suffix.lower() not in (".xyz", ".trj"):
            return None
        out_path = out_path if out_path is not None else src_path.with_suffix(".pdb")
        convert_xyz_to_pdb(src_path, ref_pdb_path, out_path)
        if out_path.exists():
            click.echo(f"[convert] Wrote '{out_path}'.")
            return out_path
        return None
    except Exception as e:
        click.echo(
            f"[convert] WARNING: Failed to convert '{src_path}' to PDB: {e}",
            err=True,
        )
        return None


# =============================================================================
# Link-freezing helpers
# =============================================================================
def parse_pdb_coords(pdb_path):
    """Parse ATOM/HETATM records from *pdb_path* and separate link hydrogen (HL) atoms.

    Returns:
        A tuple (others, lkhs) where:
            - others: list of tuples (index, x, y, z, line) for all atoms except the
              'HL' atom of residue 'LKH'. ``index`` is the 0-based position in the
              atom sequence as loaded from the *first* MODEL (or the full file if no
              MODEL records are present).
            - lkhs: list of tuples (x, y, z, line) for atoms where residue name is
              'LKH' and atom name is 'HL' in the same MODEL selection.

    Notes
    -----
        - Coordinates are read from standard PDB columns:
          X: columns 31–38, Y: 39–46, Z: 47–54 (1-based indexing).
        - If multiple MODEL blocks are present, only the first model is considered,
          matching typical geom_loader behavior.
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    others = []
    lkhs = []
    model_seen = False
    in_first_model = True
    atom_index = 0
    for line in lines:
        if line.startswith("MODEL"):
            if not model_seen:
                model_seen = True
                in_first_model = True
            else:
                in_first_model = False
            continue
        if line.startswith("ENDMDL"):
            if model_seen and in_first_model:
                break
            continue
        if model_seen and not in_first_model:
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        current_index = atom_index
        atom_index += 1

        name = line[12:16].strip()
        resname = line[17:20].strip()
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        if resname == "LKH" and name == "HL":
            lkhs.append((x, y, z, line))
        else:
            others.append((current_index, x, y, z, line))
    return others, lkhs


def nearest_index(point, pool):
    """Find the nearest point in *pool* to *point* using Euclidean distance.

    Args:
        point: Tuple (x, y, z) representing the query coordinate.
        pool: Iterable of tuples (index, x, y, z, line) to search.

    Returns:
        A tuple (index, distance) where:
            - index is the 0-based index of the nearest entry in *pool* (or -1 if *pool* is empty).
            - distance is the Euclidean distance to that entry (``inf`` if *pool* is empty).
    """
    x, y, z = point
    best_i = -1
    best_d2 = float("inf")
    for atom_index, a, b, c, _ in pool:
        d2 = (a - x) ** 2 + (b - y) ** 2 + (c - z) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_i = atom_index
    return best_i, math.sqrt(best_d2)


def load_pdb_atom_metadata(pdb_path: Path) -> List[Dict[str, Any]]:
    """Return per-atom metadata (serial, name, resname, resseq, element) in file order."""

    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            serial_txt = line[6:11].strip()
            resseq_txt = line[22:26].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            element_txt = line[76:78].strip()
            is_hetatm = line.startswith("HETATM")

            try:
                serial = int(serial_txt) if serial_txt else None
            except ValueError:
                serial = None
            try:
                resseq = int(resseq_txt) if resseq_txt else None
            except ValueError:
                resseq = None

            if not element_txt:
                inferred = guess_element(atom_name, res_name, is_hetatm)
                element_txt = inferred or ""

            atoms.append(
                {
                    "serial": serial,
                    "name": atom_name,
                    "resname": res_name,
                    "resseq": resseq,
                    "element": element_txt,
                }
            )

    return atoms


def _split_atom_spec_tokens(spec: str) -> List[str]:
    # Split an atom selector string into tokens using whitespace, comma, slash, backtick, or backslash.
    # Split the atom specification without parsing spaces by replacing spaces with commas before splitting.
    # Without replacing, it didn't work well for specs like "ALA 25 CA", somehow.
    tokens = [t for t in re.split(r"[\s/`,\\]+", spec.strip().replace(' ',',')) if t]
    return tokens


def resolve_atom_spec_index(spec: str, atom_meta: Sequence[Dict[str, Any]]) -> int:
    """
    Resolve an atom selector string into a 0-based atom index using PDB metadata.

    The selector must contain three fields (resname, resseq, atomname) separated by
    whitespace, comma, slash, backtick, or backslash. Field order is flexible; when
    unordered matching fails, the ordered interpretation (resname, resseq, atomname)
    is used as a fallback.
    """
    tokens = _split_atom_spec_tokens(spec)
    if len(tokens) != 3:
        raise ValueError(
            f"Atom spec '{spec}' must have exactly 3 fields (resname, resseq, atomname)."
        )

    tokens_upper = [t.upper() for t in tokens]
    matches: List[int] = []
    for idx, meta in enumerate(atom_meta):
        resname = (meta.get("resname") or "").strip().upper()
        resseq = meta.get("resseq")
        atom = (meta.get("name") or "").strip().upper()
        if resseq is None:
            continue
        fields = {resname, str(resseq), atom}
        if all(tok in fields for tok in tokens_upper):
            matches.append(idx)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Atom spec '{spec}' matches {len(matches)} atoms; use an explicit atom index."
        )

    resname, resseq_txt, atom = tokens_upper
    if not resseq_txt.isdigit():
        raise ValueError(
            f"Atom spec '{spec}' could not be resolved and residue number '{tokens[1]}' is not numeric."
        )
    resseq_int = int(resseq_txt)
    ordered_matches = [
        idx
        for idx, meta in enumerate(atom_meta)
        if (meta.get("resname") or "").strip().upper() == resname
        and meta.get("resseq") == resseq_int
        and (meta.get("name") or "").strip().upper() == atom
    ]
    if len(ordered_matches) == 1:
        return ordered_matches[0]
    if len(ordered_matches) > 1:
        raise ValueError(
            f"Atom spec '{spec}' matches {len(ordered_matches)} atoms after ordered fallback; "
            "use an explicit atom index."
        )

    raise ValueError(f"Atom spec '{spec}' did not match any atom.")


def values_from_bounds(low: float, high: float, h: float) -> "np.ndarray":
    """Return evenly spaced values from low→high with step cap h (inclusive)."""
    if h <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    delta = abs(high - low)
    if delta < 1e-12:
        return np.array([low], dtype=float)
    N = int(math.ceil(delta / h))
    return np.linspace(low, high, N + 1, dtype=float)


def atom_label_from_meta(atom_meta: Sequence[Dict[str, Any]], index: int) -> str:
    if index < 0 or index >= len(atom_meta):
        return f"idx{index}"
    meta = atom_meta[index]
    resname = (meta.get("resname") or "?").strip() or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    atom = (meta.get("name") or "?").strip() or "?"
    return f"{resname}-{resseq_txt}-{atom}"


def axis_label_csv(
    axis_name: str,
    i_idx: int,
    j_idx: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
    pair_raw: Optional[Tuple[Any, Any, float, float]] = None,
) -> str:
    if pair_raw and (isinstance(pair_raw[0], str) or isinstance(pair_raw[1], str)) and atom_meta:
        i_label = atom_label_from_meta(atom_meta, i_idx)
        j_label = atom_label_from_meta(atom_meta, j_idx)
        return f"{axis_name}_{i_label}_{j_label}_A"
    i_disp = i_idx + 1 if one_based else i_idx
    j_disp = j_idx + 1 if one_based else j_idx
    return f"{axis_name}_{i_disp}_{j_disp}_A"


def axis_label_html(label: str) -> str:
    parts = label.split("_")
    if len(parts) >= 4 and parts[-1] == "A":
        axis = parts[0]
        i_disp = parts[1]
        j_disp = parts[2]
        return f"{axis} ({i_disp},{j_disp}) (Å)"
    return label


def parse_scan_list_quads_checked(
    raw: str,
    *,
    expected_len: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str,
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]]]:
    parsed, raw_pairs = parse_scan_list_quads(
        raw,
        expected_len=expected_len,
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=option_name,
    )
    for i, j, low, high in parsed:
        if low <= 0.0 or high <= 0.0:
            raise click.BadParameter(f"Distances must be positive: {(i, j, low, high)}")
    return parsed, raw_pairs


def parse_scan_list_triples(
    raw: str,
    *,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str,
    return_one_based: bool = False,
) -> Tuple[List[Tuple[int, int, float]], List[Tuple[Any, Any, float]]]:
    """Parse --scan-list style triples into indices (0-based by default)."""
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not isinstance(obj, (list, tuple)):
        raise click.BadParameter(
            f"{option_name} must be a list/tuple of (i,j,target)."
        )

    parsed: List[Tuple[int, int, float]] = []
    for entry_idx, t in enumerate(obj, start=1):
        if not (
            isinstance(t, (list, tuple))
            and len(t) == 3
            and isinstance(t[2], Real)
        ):
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i,j,target): got {t}"
            )

        i = resolve_scan_index(
            t[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            t[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if return_one_based:
            i += 1
            j += 1
        parsed.append((i, j, float(t[2])))

    return parsed, list(obj)


def unbiased_energy_hartree(geom, base_calc) -> float:
    """Evaluate UMA energy (Hartree) without harmonic bias."""
    import numpy as np

    coords_bohr = np.asarray(geom.coords)
    elems = getattr(geom, "atoms", None)
    if elems is None:
        return float("nan")
    try:
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception:
        return float("nan")


def close_matplotlib_figures() -> None:
    """Best-effort cleanup for matplotlib figures to avoid open-figure warnings."""
    try:
        import matplotlib.pyplot as plt  # type: ignore
        plt.close("all")
    except Exception:
        pass


def resolve_scan_index(
    value: Any,
    *,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    context: str,
) -> int:
    """Resolve an index or atom-spec string for scan lists with consistent errors."""
    if isinstance(value, Integral):
        idx_val = int(value)
        if one_based:
            idx_val -= 1
        if idx_val < 0:
            raise click.BadParameter(
                f"Negative atom index after base conversion in {context}: {idx_val} (0-based expected)."
            )
        return idx_val
    if isinstance(value, str):
        if not atom_meta:
            raise click.BadParameter(
                f"{context} uses a string atom spec, but no PDB metadata is available."
            )
        try:
            return resolve_atom_spec_index(value, atom_meta)
        except ValueError as exc:
            raise click.BadParameter(f"{context} {exc}")
    raise click.BadParameter(f"{context} must be an int index or atom spec string.")


def parse_scan_list_quads(
    raw: str,
    *,
    expected_len: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str,
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]]]:
    """Parse --scan-list style quadruples into 0-based indices."""
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not (isinstance(obj, (list, tuple)) and len(obj) == expected_len):
        quads = ",".join([f"(i{n},j{n},low{n},high{n})" for n in range(1, expected_len + 1)])
        raise click.BadParameter(
            f"{option_name} must contain exactly {expected_len} quadruples: [{quads}]"
        )

    parsed: List[Tuple[int, int, float, float]] = []
    for entry_idx, q in enumerate(obj, start=1):
        if not (
            isinstance(q, (list, tuple))
            and len(q) == 4
            and isinstance(q[2], Real)
            and isinstance(q[3], Real)
        ):
            raise click.BadParameter(f"{option_name} entry must be (i,j,low,high): got {q}")

        i = resolve_scan_index(
            q[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            q[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        parsed.append((i, j, float(q[2]), float(q[3])))

    return parsed, list(obj)


def format_pdb_atom_metadata_header() -> str:
    """Column legend for :func:`format_pdb_atom_metadata`, aligned to match values."""

    return f"{'id':>5} {'atom':<4} {'res':<4} {'resid':>4} {'el':<2}"


def format_pdb_atom_metadata(atom_meta: Sequence[Dict[str, Any]], index: int) -> str:
    """Format metadata for atom *index* as aligned text: serial name resname resseq element."""

    fallback_serial = index + 1
    if index < 0 or index >= len(atom_meta):
        return f"{fallback_serial:>5} {'?':<4} {'?':<4} {'?':>4} {'?':<2}"

    meta = atom_meta[index]
    serial = meta.get("serial") or fallback_serial
    name = meta.get("name") or "?"
    resname = meta.get("resname") or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    element = (meta.get("element") or "?").strip() or "?"

    return f"{serial:>5} {name:<4} {resname:<4} {resseq_txt:>4} {element:<2}"


def detect_freeze_links(pdb_path):
    """Identify link-parent atom indices for 'LKH'/'HL' link hydrogens.

    For each 'HL' atom in residue 'LKH', find the nearest atom among all other
    ATOM/HETATM records and return the indices of those nearest neighbors in the
    same atom ordering used by geometry loading (first MODEL if present).

    Args:
        pdb_path: Path to the input PDB file.

    Returns:
        List of 0-based indices into the full atom sequence (including any link H atoms)
        corresponding to the nearest neighbors (link parents). Returns an empty list if
        no LKH/HL atoms are present or if link hydrogens exist without any other atoms.
    """
    others, lkhs = parse_pdb_coords(pdb_path)

    if not lkhs or not others:
        return []

    indices = []
    for (x, y, z, line) in lkhs:
        idx, dist = nearest_index((x, y, z), others)
        if idx >= 0:
            indices.append(idx)
    return indices


def detect_freeze_links_logged(pdb_path: Path) -> List[int]:
    """Return link-parent indices and raise a user-facing error on failure."""
    try:
        return list(detect_freeze_links(pdb_path))
    except Exception as e:  # pragma: no cover - defensive logging helper
        raise click.ClickException(
            f"[freeze-links] Failed to detect link parents for '{pdb_path.name}': {e}"
        ) from e


def merge_detected_freeze_links(
    geom_cfg: Dict[str, Any],
    pdb_path: Path,
    *,
    prefix: str = "[freeze-links]",
) -> List[int]:
    """Detect link-parent atoms and merge them into ``geom_cfg['freeze_atoms']``."""
    detected = detect_freeze_links_logged(pdb_path)
    merged = merge_freeze_atom_indices(geom_cfg, detected)
    if merged:
        click.echo(f"{prefix} Freeze atoms (0-based): {','.join(map(str, merged))}")
    return merged


def resolve_freeze_atoms(
    geom_cfg: Dict[str, Any],
    source_path: Optional[Path],
    freeze_links: bool,
    *,
    prefix: str = "[freeze-links]",
    on_error: str = "raise",
) -> List[int]:
    """Normalize freeze_atoms and optionally merge detected link-parent atoms."""
    merge_freeze_atom_indices(geom_cfg)
    if not freeze_links or source_path is None or source_path.suffix.lower() != ".pdb":
        return list(geom_cfg.get("freeze_atoms", []))
    try:
        return merge_detected_freeze_links(geom_cfg, source_path, prefix=prefix)
    except Exception as e:
        if on_error == "warn":
            click.echo(f"{prefix} WARNING: Could not detect link parents: {e}", err=True)
            return list(geom_cfg.get("freeze_atoms", []))
        raise


def load_prepared_geometries(
    prepared_inputs: Sequence["PreparedInputStructure"],
    *,
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
    prefix: str = "[freeze-links]",
) -> List[Any]:
    """Load multiple PreparedInputStructure geometries and apply freeze atom logic."""
    geoms: List[Any] = []

    for prepared in prepared_inputs:
        src_path = prepared.source_path
        geom_path = prepared.geom_path
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        freeze = resolve_freeze_atoms(
            cfg,
            src_path,
            auto_freeze_links,
            prefix=f"{prefix} {src_path.name}:",
        )

        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)

    return geoms
