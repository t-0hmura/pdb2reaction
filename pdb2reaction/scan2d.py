"""
scan2d — Two-distance 2D scan with harmonic restraints (UMA only)
==================================================================

Usage (CLI)
-----
    pdb2reaction scan2d -i INPUT.{pdb,xyz,trj,...} -q CHARGE \
        --scan-list "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2)]" \
        [--one-based|--zero-based] \
        --max-step-size FLOAT \
        --bias-k FLOAT \
        --relax-max-cycles INT \
        --opt-mode {light,lbfgs,heavy,rfo} \
        --freeze-links {True|False} \
        --dump {True|False} \
        --out-dir PATH \
        [--args-yaml FILE] \
        [--preopt {True|False}] \
        [--baseline {first|min}] \
        [--thresh {gau_loose|gau|gau_tight|gau_vtight|baker|never}] \
        [--zmin FLOAT] [--zmax FLOAT]

Examples::
    # 最小例（2つの距離レンジを与える）
    pdb2reaction scan2d -i input.pdb -q 0 \
        --scan-list "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

    # LBFGS・トラジェクトリ出力・正方形プロット（PNG + HTML）
    pdb2reaction scan2d -i input.pdb -q 0 \
        --scan-list "[(12,45,1.30,3.10),(10,55,1.20,3.20)]" \
        --max-step-size 0.20 --dump True --out-dir ./result_scan2d/ --opt-mode lbfgs \
        --preopt True --baseline min

Description
-----
- 2 本の原子間距離 (d1, d2) を同時に扱う 2D スキャンです。
- 入力は Python 風リスト `[(i1, j1, low1, high1), (i2, j2, low2, high2)]` を **--scan-list** で 1 つ渡します。
  - デフォルトは 1-based インデックスです（--zero-based で 0-based に切替可）。
- ステップ幅の決定は 1D スキャンと同様に、与えた `--max-step-size = h (Å)` を上限に
  - `N1 = ceil(|high1 - low1| / h)`, `N2 = ceil(|high2 - low2| / h)`
  - 値列: `d1_values = linspace(low1, high1, N1 + 1)`（差が 0 の場合は `[low1]`）
           `d2_values = linspace(low2, high2, N2 + 1)`（差が 0 の場合は `[low2]`）
- スキャンの流れ（ネスト走査）:
  1) d1 の各値 d1[i] に対して、d1 のみの拘束で最小化（バイアスは d1 のみ）。
  2) その最小化後構造をスナップショットし、**d1 を固定しつつ d2 を走査**（d2[j] へ順に）。
     各 (d1[i], d2[j]) に対して拘束付き最小化を行い、その座標で **バイアスなし（UMA 単発）** の
     エネルギーを評価・記録します（PES グリッド）。
- 出力:
  - `result_scan2d/surface.csv` … 列: i,j,d1_A,d2_A,energy_hartree,energy_kcal,bias_converged (True/False)
  - `result_scan2d/scan2d_contour.png` … 2D 等高線（PNG）
  - `result_scan2d/scan2d_surface.html` … 3D サーフェス + 底面投影（HTML）
- 中間構造や最適化の作業出力は **一時ディレクトリ** にのみ保存し、`result_scan2d` 直下には置きません。

Notes
-----
- UMA のみ（calculator は `uma_pysis`）。1D `scan.py` と同じ HarmonicBias を使用。
- 収束は LBFGS/RFO のいずれか（--opt-mode）。ステップ／トラスト半径は Å 指定上限を Bohr に変換して上限化。
- `--baseline min|first`:
  - `min`  : グリッド最小値を 0 kcal/mol として相対化（デフォルト）
  - `first`: (i=0, j=0) の点を 0 kcal/mol として相対化
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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
    format_elapsed,
    merge_freeze_atom_indices,
    normalize_choice,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    charge_option,
    spin_option,
)

# 2D スキャンで使う既定群（必要部分のみ上書き）
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

OPT_BASE_KW: Dict[str, Any] = dict(_OPT_BASE_KW)
OPT_BASE_KW.update({
    "out_dir": "./result_scan2d/",
    "dump": False,
    "max_cycles": 100,
})

LBFGS_KW: Dict[str, Any] = dict(_LBFGS_KW)
LBFGS_KW.update({"out_dir": "./result_scan2d/"})

RFO_KW: Dict[str, Any] = dict(_RFO_KW)
RFO_KW.update({"out_dir": "./result_scan2d/"})

BIAS_KW: Dict[str, Any] = {"k": 100.0}  # eV/Å^2

_OPT_MODE_ALIASES = (
    (("light", "lbfgs"), "lbfgs"),
    (("heavy", "rfo"),   "rfo"),
)

HARTREE_TO_KCAL_MOL = 627.50961


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _snapshot_geometry(g) -> Any:
    s = g.as_xyz()
    if not s.endswith("\n"):
        s += "\n"
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()
        snap = geom_loader(Path(tmp.name), coord_type=getattr(g, "coord_type", "cart"))
        try:
            import numpy as _np
            snap.freeze_atoms = _np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        return snap
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _parse_scan_list(raw: str, one_based: bool) -> Tuple[Tuple[int,int,float,float], Tuple[int,int,float,float]]:
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for --scan-list: {e}")

    if not (isinstance(obj, (list, tuple)) and len(obj) == 2):
        raise click.BadParameter("--scan-list must contain exactly two quadruples: [(i1,j1,low1,high1),(i2,j2,low2,high2)]")

    parsed: List[Tuple[int,int,float,float]] = []
    for q in obj:
        if not (
            isinstance(q, (list,tuple)) and len(q) == 4
            and isinstance(q[0], (int, np.integer))
            and isinstance(q[1], (int, np.integer))
            and isinstance(q[2], (int, float, np.floating))
            and isinstance(q[3], (int, float, np.floating))
        ):
            raise click.BadParameter(f"--scan-list entry must be (i,j,low,high): got {q}")

        i, j, low, high = int(q[0]), int(q[1]), float(q[2]), float(q[3])
        if one_based:
            i -= 1
            j -= 1
        if i < 0 or j < 0:
            raise click.BadParameter(f"Negative atom index after base conversion: {(i,j)} (0-based expected).")
        if low <= 0.0 or high <= 0.0:
            raise click.BadParameter(f"Distances must be positive: {(i,j,low,high)}")
        parsed.append((i, j, low, high))
    return parsed[0], parsed[1]


def _values_from_bounds(low: float, high: float, h: float) -> np.ndarray:
    if h <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    delta = abs(high - low)
    if delta < 1e-12:
        return np.array([low], dtype=float)
    N = int(math.ceil(delta / h))
    return np.linspace(low, high, N + 1, dtype=float)


def _make_optimizer(geom, kind: str, lbfgs_cfg: Dict[str,Any], rfo_cfg: Dict[str,Any],
                    opt_cfg: Dict[str,Any], max_step_bohr: float, relax_max_cycles: int,
                    out_dir: Path, prefix: str):
    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    if kind == "lbfgs":
        args = {**lbfgs_cfg, **common}
        args["max_step"]  = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
        args["max_cycles"] = int(relax_max_cycles)
        return LBFGS(geom, **args)
    else:
        args = {**rfo_cfg, **common}
        tr = float(rfo_cfg.get("trust_radius", 0.30))
        args["trust_radius"] = min(tr, max_step_bohr)
        args["trust_max"]    = min(float(rfo_cfg.get("trust_max", 0.30)), max_step_bohr)
        args["max_cycles"]   = int(relax_max_cycles)
        return RFOptimizer(geom, **args)


def _unbiased_energy_hartree(geom, base_calc) -> float:
    coords_bohr = np.asarray(geom.coords3d, dtype=float).reshape(-1, 3)
    for name in ("elem", "elements", "elems"):
        elems = getattr(geom, name, None)
        if elems is not None:
            try:
                return float(base_calc.get_energy(elems, coords_bohr)["energy"])
            except Exception:
                pass
    try:
        geom.set_calculator(base_calc)
        elems = getattr(geom, "elem", getattr(geom, "elements", getattr(geom, "elems", None)))
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception:
        return float("nan")


@click.command(
    help="2D distance scan with harmonic restraints (UMA only).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, ...).",
)
@charge_option()
@spin_option()
@click.option(
    "--scan-list", "scan_list_raw",
    type=str, required=True,
    help='Python-like list with two quadruples: "[(i1,j1,low1,high1),(i2,j2,low2,high2)]".',
)
@click.option("--one-based/--zero-based", "one_based", default=True, show_default=True,
              help="Interpret (i,j) indices in --scan-list as 1-based (default) or 0-based.")
@click.option("--max-step-size", type=float, default=0.20, show_default=True,
              help="Maximum step size in either distance [Å].")
@click.option("--bias-k", type=float, default=100.0, show_default=True,
              help="Harmonic well strength k [eV/Å^2].")
@click.option("--relax-max-cycles", type=int, default=10000, show_default=True,
              help="Maximum optimizer cycles per grid relaxation.")
@click.option("--opt-mode", type=str, default="light", show_default=True,
              help="Relaxation mode: light (=LBFGS) or heavy (=RFO).")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="If input is PDB, freeze parent atoms of link hydrogens.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="(Lightweight) write inner scan trajectories per d1-step as TRJ.")
@click.option("--out-dir", type=str, default="./result_scan2d/", show_default=True,
              help="Base output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file with extra args (sections: geom, calc, opt, lbfgs, rfo, bias).",
)
@click.option("--preopt", type=click.BOOL, default=True, show_default=True,
              help="preoptimize initial structure without bias before the scan.")
@click.option("--baseline", type=click.Choice(["min","first"]), default="min", show_default=True,
              help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0).")
@click.option("--zmin", type=float, default=None, show_default=False,
              help="Lower bound of color scale for plots (kcal/mol).")
@click.option("--zmax", type=float, default=None, show_default=False,
              help="Upper bound of color scale for plots (kcal/mol).")
def cli(
    input_path: Path,
    charge: Optional[int],
    spin: Optional[int],
    scan_list_raw: str,
    one_based: bool,
    max_step_size: float,
    bias_k: float,
    relax_max_cycles: int,
    opt_mode: str,
    freeze_links: bool,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    baseline: str,
    zmin: Optional[float],
    zmax: Optional[float],
) -> None:

    from .utils import load_yaml_dict, apply_yaml_overrides

    prepared_input = prepare_input_structure(input_path)
    geom_input_path = prepared_input.geom_path

    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

    try:
        time_start = time.perf_counter()

        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg   = dict(GEOM_KW)
        calc_cfg   = dict(CALC_KW)
        opt_cfg    = dict(OPT_BASE_KW)
        lbfgs_cfg  = dict(LBFGS_KW)
        rfo_cfg    = dict(RFO_KW)
        bias_cfg   = dict(BIAS_KW)

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg,  (("geom",),)),
                (calc_cfg,  (("calc",),)),
                (opt_cfg,   (("opt",),)),
                (lbfgs_cfg, (("lbfgs",),)),
                (rfo_cfg,   (("rfo",),)),
                (bias_cfg,  (("bias",),)),
            ],
        )

        calc_cfg["charge"]  = int(charge)
        calc_cfg["spin"]    = int(spin)
        opt_cfg["out_dir"]  = out_dir
        opt_cfg["dump"]     = False
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)
        if bias_k is not None:
            bias_cfg["k"] = float(bias_k)

        kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=_OPT_MODE_ALIASES,
            allowed_hint="light|lbfgs|heavy|rfo",
        )

        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        _ensure_dir(out_dir_path)
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_opt  = dict(opt_cfg); echo_opt["out_dir"] = str(out_dir_path)
        echo_bias = dict(bias_cfg)
        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("opt",  echo_opt))
        click.echo(pretty_block("lbfgs" if kind == "lbfgs" else "rfo", (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))
        click.echo(pretty_block("bias", echo_bias))

        (i1, j1, low1, high1), (i2, j2, low2, high2) = _parse_scan_list(scan_list_raw, one_based=one_based)
        click.echo(pretty_block("scan-list (0-based)",
            {"d1": (i1, j1, low1, high1), "d2": (i2, j2, low2, high2)}
        ))

        d1_values = _values_from_bounds(low1, high1, float(max_step_size))
        d2_values = _values_from_bounds(low2, high2, float(max_step_size))
        N1, N2 = len(d1_values), len(d2_values)
        click.echo(f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x:f'{x:.3f}', d1_values))}")
        click.echo(f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x:f'{x:.3f}', d2_values))}")
        click.echo(f"[grid] total grid points = {N1*N2}")

        tmp_root = Path(tempfile.mkdtemp(prefix="scan2d_tmp_"))
        grid_dir = out_dir_path / "grid"
        tmp_opt_dir  = tmp_root / "opt"
        _ensure_dir(grid_dir)
        _ensure_dir(tmp_opt_dir)

        final_dir = out_dir_path

        coord_type = geom_cfg.get("coord_type", "cart")
        geom_outer = geom_loader(geom_input_path, coord_type=coord_type)

        freeze = merge_freeze_atom_indices(geom_cfg)
        if freeze_links and input_path.suffix.lower() == ".pdb":
            detected = detect_freeze_links_safe(input_path)
            if detected:
                freeze = merge_freeze_atom_indices(geom_cfg, detected)
                if freeze:
                    click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, freeze))}")
        if freeze:
            try:
                import numpy as _np
                geom_outer.freeze_atoms = _np.array(freeze, dtype=int)
            except Exception:
                pass

        base_calc = uma_pysis(**calc_cfg)
        biased    = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

        if preopt:
            click.echo("[preopt] Unbiased relaxation of the initial structure ...")
            geom_outer.set_calculator(base_calc)
            max_step_bohr_local = float(max_step_size) * ANG2BOHR
            optimizer0 = _make_optimizer(
                geom_outer, kind, lbfgs_cfg, rfo_cfg, opt_cfg,
                max_step_bohr=max_step_bohr_local, relax_max_cycles=relax_max_cycles,
                out_dir=tmp_opt_dir, prefix="preopt_"
            )
            try:
                optimizer0.run()
            except ZeroStepLength:
                click.echo("[preopt] ZeroStepLength — continuing.", err=True)
            except OptimizationError as e:
                click.echo(f"[preopt] OptimizationError — {e}", err=True)

        max_step_bohr = float(max_step_size) * ANG2BOHR

        records: List[Dict[str, Any]] = []
        for i_idx, d1_target in enumerate(d1_values):
            click.echo(f"\n--- d1 step {i_idx+1}/{N1} : target = {d1_target:.3f} Å ---")
            geom_outer.set_calculator(biased)
            biased.set_pairs([(i1, j1, float(d1_target))])
            geom_outer.set_calculator(biased)

            opt1 = _make_optimizer(
                geom_outer, kind, lbfgs_cfg, rfo_cfg, opt_cfg,
                max_step_bohr=max_step_bohr, relax_max_cycles=relax_max_cycles,
                out_dir=tmp_opt_dir, prefix=f"d1_{i_idx:03d}_"
            )
            try:
                opt1.run()
            except ZeroStepLength:
                click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2 scan.", err=True)
            except OptimizationError as e:
                click.echo(f"[d1 {i_idx}] OptimizationError — {e}", err=True)

            geom_inner = _snapshot_geometry(geom_outer)
            geom_inner.set_calculator(biased)

            trj_blocks = [] if dump else None

            for j_idx, d2_target in enumerate(d2_values):
                biased.set_pairs([(i1, j1, float(d1_target)),
                                  (i2, j2, float(d2_target))])
                geom_inner.set_calculator(biased)

                opt2 = _make_optimizer(
                    geom_inner, kind, lbfgs_cfg, rfo_cfg, opt_cfg,
                    max_step_bohr=max_step_bohr, relax_max_cycles=relax_max_cycles,
                    out_dir=tmp_opt_dir, prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}_"
                )
                try:
                    opt2.run()
                    converged = True
                except ZeroStepLength:
                    click.echo(f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — recorded anyway.", err=True)
                    converged = False
                except OptimizationError as e:
                    click.echo(f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {e}", err=True)
                    converged = False

                E_h = _unbiased_energy_hartree(geom_inner, base_calc)

                # 中間XYZは一時ディレクトリにのみ保存
                xyz_path = grid_dir / f"point_i{i_idx:03d}_j{j_idx:03d}.xyz"
                try:
                    s = geom_inner.as_xyz()
                    if not s.endswith("\n"):
                        s += "\n"
                    with open(xyz_path, "w") as f:
                        f.write(s)
                except Exception as e:
                    click.echo(f"[write] WARNING: failed to write {xyz_path.name}: {e}", err=True)

                if dump and trj_blocks is not None:
                    sblock = geom_inner.as_xyz()
                    if not sblock.endswith("\n"):
                        sblock += "\n"
                    trj_blocks.append(sblock)

                records.append({
                    "i": int(i_idx),
                    "j": int(j_idx),
                    "d1_A": float(d1_target),
                    "d2_A": float(d2_target),
                    "energy_hartree": E_h,
                    "bias_converged": bool(converged),
                })

            if dump and trj_blocks:
                trj_path = grid_dir / f"inner_path_d1_{i_idx:03d}.trj"
                try:
                    with open(trj_path, "w") as f:
                        f.write("".join(trj_blocks))
                    click.echo(f"[write] Wrote '{trj_path}'.")
                except Exception as e:
                    click.echo(f"[write] WARNING: failed to write '{trj_path}': {e}", err=True)

        # ===== surface.csv（最終出力：result_scan2d直下）=====
        df = pd.DataFrame.from_records(records)
        if df.empty:
            click.echo("No grid records produced; aborting.", err=True)
            sys.exit(1)

        if baseline == "first":
            ref = float(df.loc[(df["i"]==0) & (df["j"]==0), "energy_hartree"].iloc[0])
        else:
            ref = float(df["energy_hartree"].min())
        df["energy_kcal"] = (df["energy_hartree"] - ref) * HARTREE_TO_KCAL_MOL

        surface_csv = final_dir / "surface.csv"
        df.to_csv(surface_csv, index=False)
        click.echo(f"[write] Wrote '{surface_csv}'.")

        # ===== プロット（RBF 200×200固定、フォーマット統一、plotはfinal_dir直下） =====
        d1_points = df["d1_A"].to_numpy(dtype=float)
        d2_points = df["d2_A"].to_numpy(dtype=float)
        z_points  = df["energy_kcal"].to_numpy(dtype=float)
        mask = np.isfinite(d1_points) & np.isfinite(d2_points) & np.isfinite(z_points)
        if not np.any(mask):
            click.echo("[plot] No finite data for plotting.", err=True)
            sys.exit(1)

        x_min, x_max = float(np.min(d1_points[mask])), float(np.max(d1_points[mask]))
        y_min, y_max = float(np.min(d2_points[mask])), float(np.max(d2_points[mask]))

        xi = np.linspace(x_min, x_max, 200)
        yi = np.linspace(y_min, y_max, 200)
        XI, YI = np.meshgrid(xi, yi)

        rbf = Rbf(d1_points[mask], d2_points[mask], z_points[mask], function="multiquadric")
        ZI = rbf(XI, YI)

        vmin = float(np.nanmin(ZI)) if zmin is None else float(zmin)
        vmax = float(np.nanmax(ZI)) if zmax is None else float(zmax)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = float(np.nanmin(ZI)), float(np.nanmax(ZI))

        # きりの良い等高線・tick用のステップ
        def _nice_step(span: float) -> float:
            if span <= 0:
                return 1.0
            raw = span / 6.0
            mag = 10 ** math.floor(math.log10(raw))
            for m in (0.5, 1, 2, 5, 10, 20):
                s = m * mag
                if raw <= s:
                    return s
            return 10.0 * mag

        c_step = _nice_step(vmax - vmin)
        c_start = math.floor(vmin / c_step) * c_step
        c_end   = math.ceil (vmax / c_step) * c_step

        # ---- 2D 等高線（PNG: サイズ明示）----
        fig2d = go.Figure(data=go.Contour(
            z=ZI,
            x=xi,
            y=yi,
            contours=dict(start=c_start, end=c_end, size=c_step),
            zmin=vmin,
            zmax=vmax,
            contours_coloring='heatmap',
            colorscale='plasma',
            colorbar=dict(
                title=dict(text="(kcal/mol)", side="top", font=dict(size=16, color='#1C1C1C')),
                tickfont=dict(size=14, color='#1C1C1C'),
                ticks="inside",
                ticklen=10,
                tickcolor='#1C1C1C',
                outlinecolor='#1C1C1C',
                outlinewidth=2,
                lenmode="fraction",
                len=1.11,
                x=1.05,
                y=0.53,
                xanchor="left",
                yanchor="middle"
            )
        ))
        fig2d.update_layout(
            width=640,
            height=600,
            xaxis_title=f"d1 ({i1+1 if one_based else i1},{j1+1 if one_based else j1}) distance (Å)",
            yaxis_title=f"d2 ({i2+1 if one_based else i2},{j2+1 if one_based else j2}) distance (Å)",
            plot_bgcolor='white',
            xaxis=dict(
                range=[x_min, x_max],
                showline=True,
                linewidth=3,
                linecolor='#1C1C1C',
                mirror=True,
                tickson='boundaries',
                ticks='inside',
                tickwidth=3,
                tickcolor='#1C1C1C',
                title_font=dict(size=18, color='#1C1C1C'),
                tickfont=dict(size=18, color='#1C1C1C'),
                tickvals=list(np.linspace(x_min, x_max, 6)),
                tickformat=".2f"
            ),
            yaxis=dict(
                range=[y_min, y_max],
                showline=True,
                linewidth=3,
                linecolor='#1C1C1C',
                mirror=True,
                tickson='boundaries',
                ticks='inside',
                tickwidth=3,
                tickcolor='#1C1C1C',
                title_font=dict(size=18, color='#1C1C1C'),
                tickfont=dict(size=18, color='#1C1C1C'),
                tickvals=list(np.linspace(y_min, y_max, 6)),
                tickformat=".2f"
            ),
            margin=dict(l=10, r=10, b=10, t=40)
        )
        png2d = final_dir / "scan2d_map.png"
        fig2d.write_image(str(png2d), scale=2, engine="kaleido", width=680, height=600)
        click.echo(f"[plot] Wrote '{png2d}'.")

        # ---- 3D サーフェス + 底面投影 ----
        # 地形は z の上半分に、平面は z 軸の完全な底面に
        spread = vmax - vmin if (vmax > vmin) else 1.0
        z_bottom = vmin - spread
        z_top    = vmax

        # z軸のtickは zmin(=vmin)未満を出さない + きりの良い数字
        z_step = _nice_step(vmax - vmin)
        z_start_tick = math.ceil(vmin / z_step) * z_step  # 最初のtickはvmin以上
        z_ticks = np.arange(z_start_tick, z_top + 0.5 * z_step, z_step).tolist()

        surface3d = go.Surface(
            x=XI, y=YI, z=ZI,
            colorscale='plasma',
            cmin=vmin, cmax=vmax,
            colorbar=dict(
                title=dict(text="(kcal/mol)", side="top", font=dict(size=16, color='#1C1C1C')),
                tickfont=dict(size=14, color='#1C1C1C'),
                ticks="inside",
                ticklen=10,
                tickcolor='#1C1C1C',
                outlinecolor='#1C1C1C',
                outlinewidth=2,
                lenmode="fraction",
                len=1.11,
                x=1.05,
                y=0.53,
                xanchor="left",
                yanchor="middle",
            ),
            contours={
                "z": {
                    "show": True,
                    "start": c_start,
                    "end": c_end,
                    "size": c_step,
                    "color": "black",
                    "project": {"z": True},
                }
            },
            name="3D Surface"
        )

        plane_proj = go.Surface(
            x=XI, y=YI, z=np.full_like(ZI, z_bottom),
            surfacecolor=ZI,
            colorscale='plasma',
            cmin=vmin, cmax=vmax,
            showscale=False,
            opacity=1.0,
            name="2D Contour Projection (Bottom)"
        )

        fig3d = go.Figure(data=[surface3d, plane_proj])
        fig3d.update_layout(
            title="Energy Landscape with 2D PES Scan",
            width=800, height=700,
            scene=dict(
                bgcolor="rgba(0,0,0,0)",
                xaxis=dict(
                    title=f"d1 ({i1+1 if one_based else i1},{j1+1 if one_based else j1}) (Å)",
                    range=[x_min, x_max],
                    showline=True, linewidth=4, linecolor='#1C1C1C', mirror=True,
                    ticks='inside', tickwidth=4, tickcolor='#1C1C1C',
                    gridcolor='rgba(0,0,0,0.1)', zerolinecolor='rgba(0,0,0,0.1)',
                    showbackground=False
                ),
                yaxis=dict(
                    title=f"d2 ({i2+1 if one_based else i2},{j2+1 if one_based else j2}) (Å)",
                    range=[y_min, y_max],
                    showline=True, linewidth=4, linecolor='#1C1C1C', mirror=True,
                    ticks='inside', tickwidth=4, tickcolor='#1C1C1C',
                    gridcolor='rgba(0,0,0,0.1)', zerolinecolor='rgba(0,0,0,0.1)',
                    showbackground=False
                ),
                zaxis=dict(
                    title="Potential Energy (kcal/mol)",
                    range=[z_bottom, z_top],
                    tickmode="array",
                    tickvals=z_ticks,  # vmin未満のtickは含めない
                    showline=True, linewidth=4, linecolor='#1C1C1C', mirror=True,
                    ticks='inside', tickwidth=4, tickcolor='#1C1C1C',
                    showgrid=True, gridcolor='rgba(0,0,0,0.1)',
                    zerolinecolor='rgba(0,0,0,0.1)',
                    showbackground=False
                ),
            ),
            margin=dict(l=10, r=20, b=10, t=40),
            paper_bgcolor="white",
        )

        html3d = final_dir / "scan2d_landscape.html"
        fig3d.write_html(str(html3d))
        click.echo(f"[plot] Wrote '{html3d}'.")

        click.echo("\n=== 2D Scan finished ===\n")
        click.echo(format_elapsed("[time] Elapsed Time for 2D Scan", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled exception during 2D scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()


if __name__ == "__main__":
    cli()
