# pdb2reaction/trj2fig.py

"""
trj2fig — Energy-profile utility for XYZ trajectories
====================================================================

Usage (CLI)
-----------
    pdb2reaction trj2fig -i traj.xyz [-o OUTPUT ...] [-r {init|None|INT}] \
        [--unit {kcal|hartree}] [--reverse-x]

Examples
--------
    # High-resolution PNG with x-axis reversed (reference is the last frame)
    pdb2reaction trj2fig -i traj.xyz --reverse-x

    # CSV + figure (reference frame #5; output values in hartree)
    pdb2reaction trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

    # Produce multiple outputs in one run (PNG, HTML, PDF)
    pdb2reaction trj2fig -i traj.xyz -o energy.png energy.html energy.pdf

Description
-----------
- Extracts Hartree energies from the second-line comment of each XYZ frame
  (uses the first decimal number on that line; scientific notation/exponents are not parsed).
- When -q/--charge or -m/--multiplicity is supplied, energies are recomputed for every
  frame with `uma_pysis` using the provided charge/spin instead of reading the comments.
- Computes either ΔE (relative to a chosen reference) or absolute E; units: kcal/mol (default) or hartree.
  Reference specification:
    - -r init  : use the initial frame (or the last frame if --reverse-x is set).
    - -r None  : use absolute energies (no reference). Also accepts "none"/"null" (case-insensitive).
    - -r <int> : use the given 0-based frame index as the reference.
- Generates a polished Plotly figure (no title) with strong ticks, consistent fonts, markers,
  and a smoothed spline curve. Supported figure outputs: PNG (default), JPEG/JPG, HTML, SVG, PDF.
- Optionally writes a CSV table of the data.
- --reverse-x flips the x-axis so the last frame appears on the left
  (and makes -r init point to that last frame).

Outputs (& Directory Layout)
----------------------------
Current working directory (default output path)
  ├─ energy.png                    # Default PNG figure when no -o is supplied
  └─ <custom outputs from -o>      # Multiple filenames may be provided; Click also allows repeating -o
        • .png / .jpg / .jpeg  → Raster Plotly exports (PNG uses scale=2)
        • .html                → Standalone interactive Plotly figure
        • .svg / .pdf          → Vector exports
        • .csv                 → Tabular data with columns ``frame``, ``energy_hartree``, and the relevant ΔE/E column based on --unit/-r

Notes
-----
- The legacy --output-peak option has been removed.
- If a frame comment lacks a parseable decimal number, an error is raised; if no energies are found, the run fails.
- Unsupported file extensions in -o cause an error.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import click
import plotly.graph_objs as go
from ase import Atoms
from ase.io import read
from pysisyphus.constants import AU2KCALPERMOL, ANG2BOHR

from .uma_pysis import uma_pysis

AXIS_WIDTH = 3         # axis and tick thickness
FONT_SIZE = 18         # tick-label font size
AXIS_TITLE_SIZE = 20   # axis-title font size
LINE_WIDTH = 2         # curve width
MARKER_SIZE = 6        # marker size


# ---------------------------------------------------------------------
#  File helpers
# ---------------------------------------------------------------------
def read_energies_xyz(fname: Path | str) -> List[float]:
    """
    Extract Hartree energies from the second-line comment of each XYZ frame.

    The first decimal number found on the comment line is used
    (scientific notation/exponents are not parsed).
    """
    energies: List[float] = []
    with open(fname, encoding="utf-8") as fh:
        while (hdr := fh.readline()):
            try:
                nat = int(hdr.strip())
            except ValueError:  # reached a non-XYZ header
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


def recompute_energies(
    traj_path: Path, charge: Optional[int], multiplicity: Optional[int]
) -> List[float]:
    """
    Recalculate Hartree energies for every frame using uma_pysis.
    """

    frames_obj = read(traj_path, index=":", format="xyz")
    frames = [frames_obj] if isinstance(frames_obj, Atoms) else list(frames_obj)
    if not frames:
        raise RuntimeError(f"No frames found in {traj_path}")

    calc = uma_pysis(charge=charge or 0, spin=multiplicity or 1)
    energies: List[float] = []
    for atoms in frames:
        elems = atoms.get_chemical_symbols()
        coords_bohr = atoms.get_positions() * ANG2BOHR
        energies.append(float(calc.get_energy(elems, coords_bohr)["energy"]))

    return energies


# ---------------------------------------------------------------------
#  Transformation
# ---------------------------------------------------------------------
def _parse_reference_spec(spec: str | None) -> str | int | None:
    """
    Normalize the reference specification:
      - "init" (case-insensitive) -> "init"
      - "none"/"null"             -> None
      - integer-like string       -> int
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
    Decide the reference index and whether to compute a ΔE series.

    Returns (reference_index or None, is_delta).
    """
    if ref_spec is None:
        return None, False  # absolute energies
    if ref_spec == "init":
        idx = 0 if not reverse_x else n_frames - 1
        return idx, True
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
    Compute the y-series and its axis label.

    Returns (values, ylabel, is_delta).
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
    """
    Build a Plotly figure without a title.
    """
    fig = go.Figure(
        go.Scatter(
            x=list(range(len(delta_or_abs))))
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
    """
    Write all requested outputs.
    """
    for out in outs:
        ext = out.suffix.lower()
        if ext == ".csv":
            write_csv(out, energies, values, unit, is_delta)
        elif ext == ".html":
            assert fig is not None
            fig.write_html(out)
            print(f"[trj2fig] Saved figure -> {out}")
        elif ext in {".png", ".jpg", ".jpeg", ".pdf", ".svg"}:
            assert fig is not None
            kw = {"engine": "kaleido"}
            if ext == ".png":
                kw["scale"] = 2  # high-resolution PNG
            fig.write_image(out, **kw)
            print(f"[trj2fig] Saved figure -> {out}")
        else:
            raise ValueError(f"Unsupported format: {ext}")


def write_csv(
    out: Path,
    energies_hartree: Sequence[float],
    series: Sequence[float],
    unit: str,
    is_delta: bool,
) -> None:
    """
    Save energies (hartree) and ΔE/E series to CSV.
    """
    colname = (f"delta_{unit}" if is_delta else f"energy_{unit}")
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "energy_hartree", colname])
        for i, (eh, y) in enumerate(zip(energies_hartree, series)):
            w.writerow([i, f"{eh:.8f}", f"{y:.6f}"])
    print(f"[trj2fig] Saved CSV -> {out}")


# ---------------------------------------------------------------------
#  CLI (argparse)
# ---------------------------------------------------------------------
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="trj2fig",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot ΔE or E from an XYZ trajectory and export a figure and/or CSV (no title).",
    )
    p.add_argument("-i", "--input", required=True, help="XYZ trajectory file")
    p.add_argument(
        "-o",
        "--out",
        nargs="+",
        default=["energy.png"],
        help="Output file(s) [.png/.jpg/.jpeg/.html/.svg/.pdf/.csv]. Multiple names allowed.",
    )
    p.add_argument("--unit", choices=["kcal", "hartree"], default="kcal", help="Energy unit")
    p.add_argument(
        "-r",
        "--reference",
        default="init",
        help='Reference: "init" (initial frame; last frame if --reverse-x), "None" (absolute E), or an integer index.',
    )
    p.add_argument("-q", "--charge", type=int, required=False, help="Total charge (recompute energies when provided)")
    p.add_argument(
        "-m",
        "--multiplicity",
        type=int,
        required=False,
        help="Spin multiplicity (2S+1). Recompute energies when provided.",
    )
    p.add_argument(
        "--reverse-x",
        action="store_true",
        help="Reverse the x-axis (last frame on the left).",
    )
    return p.parse_args()


def run_trj2fig(
    input_path: Path,
    outs: Sequence[Path],
    unit: str,
    reference: str,
    reverse_x: bool,
    charge: Optional[int] = None,
    multiplicity: Optional[int] = None,
) -> None:
    traj = input_path.expanduser().resolve()
    if not traj.is_file():
        raise FileNotFoundError(traj)

    if charge is None and multiplicity is None:
        energies = read_energies_xyz(traj)
    else:
        energies = recompute_energies(traj, charge, multiplicity)
    values, ylabel, is_delta = transform_series(energies, reference, unit, reverse_x)

    need_plot = any(Path(o).suffix.lower() != ".csv" for o in outs)
    fig = build_figure(values, ylabel, reverse_x) if need_plot else None

    out_paths = [Path(o).expanduser().resolve() for o in outs]
    save_outputs(out_paths, fig, energies, values, unit, is_delta)


def main() -> None:
    args = parse_cli()
    run_trj2fig(
        Path(args.input),
        args.out,
        args.unit,
        args.reference,
        args.reverse_x,
        args.charge,
        args.multiplicity,
    )


# ---------------------------------------------------------------------
#  Click wrapper for package CLI integration
# ---------------------------------------------------------------------
@click.command(
    name="trj2fig",
    help="Plot ΔE or E from an XYZ trajectory and export figure/CSV.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="XYZ trajectory file",
)
@click.option(
    "-o",
    "--out",
    "outs",
    multiple=True,
    default=(),
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output file(s). You can repeat -o and/or list extra filenames after options "
         "(.png/.jpg/.jpeg/.html/.svg/.pdf/.csv). If nothing is given, defaults to energy.png.",
)
@click.argument(
    "extra_outs",
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
    help='Reference: "init" (initial frame; last frame if --reverse-x), "None" (absolute E), or an integer index.',
)
@click.option(
    "-q",
    "--charge",
    type=int,
    default=None,
    help="Total charge. Triggers energy recomputation when supplied.",
)
@click.option(
    "-m",
    "--multiplicity",
    type=int,
    default=None,
    help="Spin multiplicity (2S+1). Triggers energy recomputation when supplied.",
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
    charge: Optional[int],
    multiplicity: Optional[int],
    reverse_x: bool,
) -> None:
    # Combine outputs from -o with positional filenames that follow the options
    all_outs: List[Path] = list(outs) + list(extra_outs)
    if not all_outs:
        all_outs = [Path("energy.png")]
    run_trj2fig(input_path, all_outs, unit, reference, reverse_x, charge, multiplicity)
