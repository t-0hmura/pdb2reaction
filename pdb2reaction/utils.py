# pdb2reaction/utils.py

import sys
import math
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, List

import yaml
from ase.io import read, write
import plotly.graph_objs as go


# =============================================================================
# Generic helpers
# =============================================================================


def deep_update(dst: Dict[str, Any], src: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Recursively update mapping *dst* with *src*, returning *dst*."""

    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def load_yaml_dict(path: Optional[Path]) -> Dict[str, Any]:
    """Load a YAML file whose root must be a mapping. Return an empty dict if *path* is None."""

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
        各状態のエネルギー（同一単位）。ここでは値をそのまま描画し、変換は行いません。
    labels : Sequence[str]
        各状態に対応するラベル（例: ["R", "TS1", "IM1", "TS2", "P"]）。
        `energies` と同じ長さである必要があります。
    ylabel : str, optional
        縦軸ラベル（例: "ΔE", "ΔG" など）。既定は "ΔE"。
    baseline : bool, optional
        True の場合、最初のエネルギー高さに横点線のベースラインを図全体に描画します。

    Returns
    -------
    plotly.graph_objs.Figure
        エネルギーダイアグラムの Figure。

    描画仕様
    --------
    - 各状態は「太い横線」で描画します（線幅は HLINE_WIDTH）。
    - 隣接する状態は、左側の「右端」から右側の「左端」へ「斜めの点線」で接続します。
    - 横線の長さは状態数 n に応じて自動調整し、必ず隣接線との間に x 方向の隙間が生じます。
    - x 軸の目盛は各状態の中心に設定し、表示文字列は `labels` を使用します。
    """
    if len(energies) == 0:
        raise ValueError("`energies` は 1 つ以上の値を含む必要があります。")
    if len(energies) != len(labels):
        raise ValueError("`energies` と `labels` は同じ長さである必要があります。")

    n = len(energies)
    energies = [float(e) for e in energies]

    # -----------------------------
    # レイアウト/スタイル定数
    # -----------------------------
    AXIS_WIDTH = 3
    FONT_SIZE = 18
    AXIS_TITLE_SIZE = 20
    HLINE_WIDTH = 6           # 状態レベル（横線）の太さ
    CONNECTOR_WIDTH = 2       # 斜め点線の太さ
    LINE_COLOR = "#1C1C1C"
    GRID_COLOR = "lightgrey"

    # -----------------------------
    # X 方向の幾何（中心位置と横線長）
    # -----------------------------
    # 状態中心を 0.5, 1.5, 2.5, ... に配置（等間隔）
    centers = [i + 0.5 for i in range(n)]

    # 横線の長さは n が増えると短くする（最小 0.35, 最大 0.85）
    # 例: n=5 -> 0.7, n=10 -> 0.5, n>=20 -> 0.35
    seg_width = min(0.85, max(0.35, 0.90 - 0.04 * n))
    half = seg_width / 2.0

    lefts = [c - half for c in centers]
    rights = [c + half for c in centers]

    # -----------------------------
    # Figure 組み立て
    # -----------------------------
    fig = go.Figure()

    # baseline（最初のエネルギー高さに横点線）
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

    # 各状態の太い横線
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

    # 隣接状態間の斜めの点線（右端 -> 左端）
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
    # 軸レンジ・体裁
    # -----------------------------
    # X は最初と最後の横線から少しだけ余白をとる
    xpad = max(0.08, 0.15 * (1.0 - seg_width))
    x_min = lefts[0] - xpad
    x_max = rights[-1] + xpad

    # Y は上下に余白をとる
    y_min = min(energies)
    y_max = max(energies)
    span = max(1e-6, y_max - y_min)  # 全て同じ高さでもゼロ高にならないよう回避
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
# Coordinate conversion utilities (既存：変更なし)
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
# Link-freezing helpers (既存：変更なし)
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
