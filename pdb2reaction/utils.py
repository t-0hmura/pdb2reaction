# pdb2reaction/utils.py

"""
utils — concise utilities for configuration, plotting, coordinates, and link-freezing
====================================================================

Usage (API)
-----
    from pdb2reaction.utils import (
        build_energy_diagram,
        convert_xyz_to_pdb,
        freeze_links,
        merge_freeze_atom_indices,
        normalize_choice,
        pretty_block,
    )

Examples::
    >>> from pathlib import Path
    >>> block = pretty_block("Geometry", {"freeze_atoms": [0, 1, 5]})
    >>> diagram = build_energy_diagram([0.0, 12.3, 5.4], ["R", "TS", "P"])
    >>> indices = freeze_links(Path("pocket.pdb"))

Description
-----
- **Generic helpers**
  - `pretty_block(title, content)`: Return a YAML-formatted block with an underlined title. Uses `yaml.safe_dump` with `allow_unicode=True`, `sort_keys=False`. Renders `(empty)` when `content` is empty.
  - `format_geom_for_echo(geom_cfg)`: Normalize geometry configuration for CLI echo. If `"freeze_atoms"` is an iterable (but not a string), convert it to a comma-separated string; `None`/string/other types are left unchanged. Empty iterables become `""`.
  - `format_elapsed(prefix, start_time, end_time=None)`: Format a wall-clock duration (HH:MM:SS.sss) given a start time and optional end time, using `time.perf_counter()` when the end time is omitted.
  - `merge_freeze_atom_indices(geom_cfg, *indices)`: Merge one or more iterables of atom indices into `geom_cfg["freeze_atoms"]`. Preserve existing entries, de-duplicate, sort numerically, and return the updated list (in place).
  - `normalize_choice(value, *, param, alias_groups, allowed_hint)`: Canonicalise CLI-style string options using alias groups. Returns the mapped value or raises `click.BadParameter` with the provided hint when no alias matches.
  - `deep_update(dst, src)`: Recursively update mapping `dst` with `src`. Nested dicts are merged, non-dicts overwrite; returns `dst`.
  - `_get_mapping_section(cfg, path)`: Internal helper to resolve a nested mapping section. Returns a `dict` or `None`.
  - `apply_yaml_overrides(yaml_cfg, overrides)`: For each target dictionary and its candidate key paths, find the first existing path in `yaml_cfg` and apply it via `deep_update`. Centralizes repeated `yaml_cfg.get(...)`-style merging.
  - `load_yaml_dict(path)`: Load a YAML file whose root must be a mapping. Returns `{}` when `path` is `None`. Raises `ValueError` if the YAML root is not a mapping.

- **Plotly: Energy diagram builder**
  - `build_energy_diagram(energies, labels, ylabel="ΔE", baseline=False, showgrid=False)`:
    Render an energy diagram where each state is a thick horizontal segment and adjacent states are connected by dotted diagonals (right end of left state → left end of right state). Segment length shrinks as the number of states grows to keep gaps readable. X ticks are centered on states and labeled by `labels`. Optional dotted baseline at the first state’s energy; optional grid. Energies are plotted as provided (no unit conversion). Returns a `plotly.graph_objs.Figure`. Validates equal lengths for `energies`/`labels` and non-empty input.

- **Coordinate conversion utilities**
  - `convert_xyz_to_pdb(xyz_path, ref_pdb_path, out_pdb_path)`:
    Overlay coordinates from an XYZ file (single or multi-frame) onto the atom ordering/topology of a reference PDB and write to `out_pdb_path`. The first frame creates/overwrites; subsequent frames append using `MODEL`/`ENDMDL`. Implemented with ASE (`ase.io.read`/`write`). Raises `ValueError` if no frames are found in the XYZ.

- **Link-freezing helpers**
  - `parse_pdb_coords(pdb_path)`: Parse `ATOM`/`HETATM` records, separating all atoms (as “others”) from link hydrogens `HL` in residue `LKH` (as “lkhs”). Coordinates are read from standard PDB columns (1‑based): X 31–38, Y 39–46, Z 47–54. Returns `(others, lkhs)` as lists of `(x, y, z, line)`.
  - `nearest_index(point, pool)`: Find the Euclidean nearest neighbor of a given `(x, y, z)` within `pool`; returns `(index, distance)` where `index` is 0‑based or `-1` if `pool` is empty (distance will be `inf` in that case).
  - `freeze_links(pdb_path)`: For each `LKH`/`HL` atom, find the nearest atom among all other `ATOM`/`HETATM` records and return the corresponding 0‑based indices into the sequence of non‑`LKH` atoms (“others”). Returns an empty list if no link hydrogens are present.

Outputs (& Directory Layout)
-----
- This module does not create directories.
- Functions primarily return Python objects or mutate dictionaries in place.
- On-disk output occurs only when explicitly requested by the caller:
  - `convert_xyz_to_pdb` writes a PDB file to `out_pdb_path` (first frame create/overwrite; subsequent frames append with `MODEL`/`ENDMDL` blocks).
  - `build_energy_diagram` returns a Plotly `Figure`; it does not write files unless the caller saves/exports the figure.

Notes:
-----
- Energy units in `build_energy_diagram` are passed through unchanged; ensure consistent units across states.
- Axis/line styling in `build_energy_diagram` is fixed-width with automatic padding; segment length adapts to the number of states.
- `load_yaml_dict` uses `yaml.safe_load` and enforces a mapping at the YAML root; empty files yield `{}`.
- `apply_yaml_overrides` tries candidate key paths in order and applies only the first existing mapping section per target.
- `parse_pdb_coords` skips unparseable coordinate fields and considers only `ATOM`/`HETATM` lines.
- Dependencies: PyYAML, ASE (`ase.io.read`/`write`), Plotly (graph objects), Click (for CLI error reporting). Ensure these are installed.
"""

import math
import time
from collections.abc import Iterable as _Iterable, Mapping, Sequence as _Sequence
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, List, Tuple

import click
import yaml
from ase.io import read, write
import plotly.graph_objs as go


# =============================================================================
# Generic helpers
# =============================================================================


def pretty_block(title: str, content: Dict[str, Any]) -> str:
    """
    Return a YAML-formatted block with an underlined title.
    """

    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalise geometry configuration for CLI echo output.
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

    g["freeze_atoms"] = ",".join(map(str, items)) if items else ""
    return g


def format_elapsed(prefix: str, start_time: float, end_time: Optional[float] = None) -> str:
    """Return a formatted elapsed-time string with the provided ``prefix`` label."""

    finish = end_time if end_time is not None else time.perf_counter()
    elapsed = max(0.0, finish - start_time)
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{prefix}: {int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"


def merge_freeze_atom_indices(
    geom_cfg: Dict[str, Any],
    *indices: _Iterable[int],
) -> List[int]:
    """Merge one or more iterables of indices into ``geom_cfg['freeze_atoms']``.

    Existing entries are preserved, duplicates removed, and the result sorted.
    The updated list is returned.
    """

    merged: set[int] = set()
    base = geom_cfg.get("freeze_atoms", [])
    if isinstance(base, _Iterable):
        merged.update(int(i) for i in base)
    for seq in indices:
        if seq is None:
            continue
        merged.update(int(i) for i in seq)
    result = sorted(merged)
    geom_cfg["freeze_atoms"] = result
    return result


def normalize_choice(
    value: str,
    *,
    param: str,
    alias_groups: _Sequence[tuple[_Sequence[str], str]],
    allowed_hint: str,
) -> str:
    """Normalise *value* using alias groups and raise ``click.BadParameter`` on failure."""

    key = (value or "").strip().lower()
    for aliases, canonical in alias_groups:
        if any(key == alias.lower() for alias in aliases):
            return canonical

    hint = allowed_hint.strip()
    detail = f" Allowed: {hint}." if hint else ""
    raise click.BadParameter(f"Unknown value for {param} '{value}'.{detail}")


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
        while centralising the shared logic.
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


# =============================================================================
# Plotly: Energy diagram builder
# =============================================================================
def build_energy_diagram(
    energies: Sequence[float],
    labels: Sequence[str],
    ylabel: str = "ΔE",
    baseline: bool = False,
    showgrid=False,
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

    Notes:
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

    for step, frame in enumerate(traj):
        atoms = ref_atoms.copy()
        atoms.set_positions(frame.get_positions())
        if step == 0:
            write(out_pdb_path, atoms)  # Create/overwrite on the first frame
        else:
            write(out_pdb_path, atoms, append=True)  # Append subsequent frames using MODEL/ENDMDL


# =============================================================================
# Link-freezing helpers
# =============================================================================
def parse_pdb_coords(pdb_path):
    """Parse ATOM/HETATM records from *pdb_path* and separate link hydrogen (HL) atoms.

    Returns:
        A tuple (others, lkhs) where:
            - others: list of tuples (x, y, z, line) for all atoms except the 'HL' atom
              of residue 'LKH'.
            - lkhs: list of tuples (x, y, z, line) for atoms where residue name is 'LKH'
              and atom name is 'HL'.

    Notes:
        - Coordinates are read from standard PDB columns:
          X: columns 31–38, Y: 39–46, Z: 47–54 (1-based indexing).
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    others = []
    lkhs = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        name    = line[12:16].strip()
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
            others.append((x, y, z, line))
    return others, lkhs


def nearest_index(point, pool):
    """Find the nearest point in *pool* to *point* using Euclidean distance.

    Args:
        point: Tuple (x, y, z) representing the query coordinate.
        pool: Iterable of tuples (x, y, z, line) to search.

    Returns:
        A tuple (index, distance) where:
            - index is the 0-based index of the nearest entry in *pool* (or -1 if *pool* is empty).
            - distance is the Euclidean distance to that entry.
    """
    x, y, z = point
    best_i = -1
    best_d2 = float("inf")
    for i, (a, b, c, _) in enumerate(pool):
        d2 = (a - x) ** 2 + (b - y) ** 2 + (c - z) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_i = i
    return best_i, math.sqrt(best_d2)


def freeze_links(pdb_path):
    """Identify link-parent atom indices for 'LKH'/'HL' link hydrogens.

    For each 'HL' atom in residue 'LKH', find the nearest atom among all other
    ATOM/HETATM records and return the indices of those nearest neighbors.

    Args:
        pdb_path: Path to the input PDB file.

    Returns:
        List of 0-based indices into the sequence of non-LKH atoms ("others") corresponding
        to the nearest neighbors (link parents). Returns an empty list if no LKH/HL atoms
        are present.
    """
    others, lkhs = parse_pdb_coords(pdb_path)

    if not lkhs:
        return []

    indices = []
    for (x, y, z, line) in lkhs:
        idx, dist = nearest_index((x, y, z), others)
        indices.append(idx)
    return indices
