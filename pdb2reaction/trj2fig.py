# pdb2reaction/trj2fig.py

"""
Energy-profile utility for XYZ trajectories
======================================================

機能
----

• XYZ各フレームの2行目コメントから Hartree エネルギーを抽出
• 参照指定に応じて「ΔE」または「絶対E」を出力（kcal mol⁻¹ or hartree）
  - `-r init` : 初期フレームを基準（--reverse-x のときは最後のフレーム）
  - `-r None` : 絶対エネルギー（基準なし）
  - `-r <int>`: 指定フレームを基準
• 目盛太字・フォント・スプライン曲線など整えた Plotly 図を生成
  （PNG[default] / HTML / SVG / PDF）
• 表データを CSV で書き出し可能
• `--reverse-x` は x 軸を反転（最後のフレームが左側）

使い方例
--------

図（PNG高解像度）を生成（x軸を反転、基準は左端＝最後のフレーム）

    trj2fig -i traj.xyz --reverse-x

CSV と 図を同時出力（基準フレーム #5、単位 hartree）

    trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

複数出力（PNG, HTML, PDF を一括出力）

    trj2fig -i traj.xyz -o energy.png energy.html energy.pdf

備考
----

• 旧版の `--output-peak` 機能は削除しました。
• 既定の出力は `energy.png` です。`-o` に複数ファイル名を与えるか、
  `-o` を繰り返し指定すると複数出力できます。
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import click
import plotly.graph_objs as go
from pysisyphus.constants import AU2KCALPERMOL

AXIS_WIDTH = 3    # axis & tick thickness
FONT_SIZE = 18    # tick-label font size
AXIS_TITLE_SIZE = 20  # axis-label font size
LINE_WIDTH = 2    # curve width
MARKER_SIZE = 6   # marker size


# ---------------------------------------------------------------------
#  File helpers
# ---------------------------------------------------------------------
def read_energies_xyz(fname: Path | str) -> List[float]:
    """Return a list of Hartree energies extracted from 2-line comments."""
    energies: List[float] = []
    with open(fname, encoding="utf-8") as fh:
        while (hdr := fh.readline()):
            try:
                nat = int(hdr.strip())
            except ValueError:  # non-XYZ header reached
                break
            comment = fh.readline().strip()
            m = re.search(r"(-?\d+(?:\.\d+)?)", comment)
            if not m:
                raise RuntimeError(f"Energy not found in comment: {comment}")
            energies.append(float(m.group(1)))
            for _ in range(nat):  # skip coordinates
                fh.readline()
    if not energies:
        raise RuntimeError(f"No energy data in {fname}")
    return energies


# ---------------------------------------------------------------------
#  Transformation
# ---------------------------------------------------------------------
def _parse_reference_spec(spec: str | None) -> str | int | None:
    """
    Normalize reference spec:
      - "init" (case-insensitive) -> "init"
      - "none"/"null" -> None
      - integer-like string -> int
    """
    if spec is None:
        return "init"
    s = str(spec).strip()
    lower = s.lower()
    if lower in {"none", "null"}:
        return None
    if lower == "init":
        return "init"
    try:
        return int(s)
    except ValueError:
        raise ValueError(
            f'Invalid -r/--reference: {spec!r}. Use "init", "None", or an integer index.'
        )


def _resolve_reference_index(
    n_frames: int, ref_spec: str | int | None, reverse_x: bool
) -> Tuple[Optional[int], bool]:
    """
    Decide reference index and whether to produce delta.

    Returns (ref_index or None, is_delta)
    """
    if ref_spec is None:
        return None, False  # absolute energies
    if ref_spec == "init":
        idx = 0 if not reverse_x else n_frames - 1
        return idx, True
    # integer index
    idx = int(ref_spec)
    if idx < 0 or idx >= n_frames:
        raise IndexError(f"Reference index {idx} out of range (0..{n_frames-1}).")
    return idx, True


def transform_series(
    energies_hartree: Sequence[float],
    ref_spec_raw: str | None,
    unit: str,
    reverse_x: bool,
) -> Tuple[List[float], str, bool]:
    """
    Compute y-series and label.

    Returns (values, ylabel, is_delta)
    """
    ref_spec = _parse_reference_spec(ref_spec_raw)
    ref_idx, is_delta = _resolve_reference_index(len(energies_hartree), ref_spec, reverse_x)

    scale = AU2KCALPERMOL if unit == "kcal" else 1.0
    if is_delta:
        base = energies_hartree[ref_idx]  # type: ignore[index]
        values = [float((e - base) * scale) for e in energies_hartree]
        ylabel = f"ΔE ({'kcal/mol' if unit == 'kcal' else 'hartree'})"
    else:
        values = [float(e * scale) for e in energies_hartree]
        ylabel = f"E ({'kcal/mol' if unit == 'kcal' else 'hartree'})"

    return values, ylabel, is_delta


# ---------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------
def _axis_template() -> dict:
    return dict(
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor="#1C1C1C",
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor="#1C1C1C",
        tickfont=dict(size=FONT_SIZE, color="#1C1C1C"),
        gridcolor="lightgrey",
        gridwidth=0.5,
        zeroline=False,
    )


def build_figure(delta_or_abs: Sequence[float], ylabel: str, reverse_x: bool) -> go.Figure:
    """Build a title-less Plotly figure."""
    fig = go.Figure(
        go.Scatter(
            x=list(range(len(delta_or_abs)))),
    )
    fig.data[0].update(
        y=list(delta_or_abs),
        mode="lines+markers",
        marker=dict(size=MARKER_SIZE),
        line=dict(shape="spline", smoothing=1.0, width=LINE_WIDTH),
    )

    xaxis_conf = _axis_template() | {
        "title": dict(text="Frame", font=dict(size=AXIS_TITLE_SIZE, color="#1C1C1C"))
    }
    if reverse_x:
        xaxis_conf["autorange"] = "reversed"

    fig.update_layout(
        xaxis=xaxis_conf,
        yaxis=_axis_template() | {
            "title": dict(text=ylabel, font=dict(size=AXIS_TITLE_SIZE, color="#1C1C1C"))
        },
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=80, r=40, t=40, b=80),
    )
    return fig


def save_outputs(
    outs: Sequence[Path],
    fig: Optional[go.Figure],
    energies: Sequence[float],
    values: Sequence[float],
    unit: str,
    is_delta: bool,
) -> None:
    """Write all requested outputs."""
    for out in outs:
        ext = out.suffix.lower()
        if ext == ".csv":
            write_csv(out, energies, values, unit, is_delta)
        elif ext == ".html":
            assert fig is not None
            fig.write_html(out)
            print(f"[trj2fig] Figure → {out}")
        elif ext in {".png", ".jpg", ".jpeg", ".pdf", ".svg"}:
            assert fig is not None
            kw = {"engine": "kaleido"}
            if ext == ".png":
                kw["scale"] = 2  # hi-res PNG
            fig.write_image(out, **kw)
            print(f"[trj2fig] Figure → {out}")
        else:
            raise ValueError(f"Unsupported format: {ext}")


def write_csv(
    out: Path,
    energies_hartree: Sequence[float],
    series: Sequence[float],
    unit: str,
    is_delta: bool,
) -> None:
    """Save energies (hartree) and ΔE/E series to CSV."""
    colname = (f"delta_{unit}" if is_delta else f"energy_{unit}")
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "energy_hartree", colname])
        for i, (eh, y) in enumerate(zip(energies_hartree, series)):
            w.writerow([i, f"{eh:.8f}", f"{y:.6f}"])
    print(f"[trj2fig] CSV → {out}")


# ---------------------------------------------------------------------
#  CLI (argparse)
# ---------------------------------------------------------------------
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="trj2fig",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot ΔE/E from an XYZ trajectory, and export figure/CSV (no title).",
    )
    p.add_argument("-i", "--input", required=True, help="XYZ trajectory file")
    p.add_argument(
        "-o",
        "--out",
        nargs="+",
        default=["energy.png"],
        help="Output file(s) [.png/.html/.svg/.pdf/.csv]. Multiple names allowed.",
    )
    p.add_argument("--unit", choices=["kcal", "hartree"], default="kcal", help="Energy unit")
    p.add_argument(
        "-r",
        "--reference",
        default="init",
        help='Reference: "init" (or "None" for absolute E, or integer frame index)',
    )
    p.add_argument(
        "--reverse-x",
        action="store_true",
        help="Reverse the x-axis (last frame on the left).",
    )
    return p.parse_args()


def run_trj2fig(input_path: Path, outs: Sequence[Path], unit: str, reference: str, reverse_x: bool) -> None:
    traj = input_path.expanduser().resolve()
    if not traj.is_file():
        raise FileNotFoundError(traj)

    energies = read_energies_xyz(traj)
    values, ylabel, is_delta = transform_series(energies, reference, unit, reverse_x)

    need_plot = any(Path(o).suffix.lower() != ".csv" for o in outs)
    fig = build_figure(values, ylabel, reverse_x) if need_plot else None

    out_paths = [Path(o).expanduser().resolve() for o in outs]
    save_outputs(out_paths, fig, energies, values, unit, is_delta)


def main() -> None:
    args = parse_cli()
    run_trj2fig(Path(args.input), args.out, args.unit, args.reference, args.reverse_x)


# ---------------------------------------------------------------------
#  Click wrapper for package CLI integration
# ---------------------------------------------------------------------
@click.command(
    name="trj2fig",
    help="Plot ΔE/E from an XYZ trajectory and export figure/CSV.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",  # 内部引数名を明示
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="XYZ trajectory file",
)
@click.option(
    "-o",
    "--out",
    "outs",
    multiple=True,                      # -o を繰り返し指定できる
    default=(),                         # 既定は空（後で自前で既定値を補う）
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output file(s). You can repeat -o, and/or list extra filenames after options "
         "(.png/.html/.svg/.pdf/.csv). If nothing is given, defaults to energy.png.",
)
@click.argument(
    "extra_outs",                        # -o の後ろに続く余剰のファイル名も受け取る
    nargs=-1,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--unit",
    type=click.Choice(["kcal", "hartree"]),
    default="kcal",
    help="Energy unit.",
)
@click.option(
    "-r",
    "--reference",
    default="init",
    help='Reference: "init" (or "None" for absolute E, or integer frame index).',
)
@click.option(
    "--reverse-x",
    is_flag=True,
    help="Reverse the x-axis (last frame on the left).",
)
def cli(
    input_path: Path,
    outs: Tuple[Path, ...],
    extra_outs: Tuple[Path, ...],
    unit: str,
    reference: str,
    reverse_x: bool,
) -> None:
    # -o で与えられたもの + 位置引数で与えられたものを結合
    all_outs: List[Path] = list(outs) + list(extra_outs)
    if not all_outs:
        all_outs = [Path("energy.png")]
    run_trj2fig(input_path, all_outs, unit, reference, reverse_x)


if __name__ == "__main__":
    main()
