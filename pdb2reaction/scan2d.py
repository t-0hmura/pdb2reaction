"""
2d-scan — Two-distance 2D scan with harmonic restraints (UMA only)
==================================================================

Usage (CLI)
-----
    pdb2reaction 2d-scan -i INPUT.{pdb,xyz,trj,...} -q CHARGE \
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
    pdb2reaction 2d-scan -i input.pdb -q 0 \
        --scan-list "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

    # LBFGS・トラジェクトリ出力・正方形プロット（PNG + HTML）
    pdb2reaction 2d-scan -i input.pdb -q 0 \
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
  - `out_dir/grid/point_i###_j###.xyz` … 各格子点の最終座標（バイアス付き緩和の結果）
  - `out_dir/surface.csv` … 列: i,j,d1_A,d2_A,energy_hartree,energy_kcal,bias_converged (True/False)
  - `out_dir/plots/scan2d_contour.png` … 2D 等高線（Plotly, kaleido があれば PNG。無ければ HTML）
  - `out_dir/plots/scan2d_surface.html` … 3D サーフェス + 底面投影
- 図は添付 Plotly スクリプトを参考にし、**正方形**になるように寸法・比率を整えています。
  （オプション `--zmin/--zmax` でカラースケールの上下限を指定可能。）

Notes
-----
- UMA のみ（calculator は `uma_pysis`）。1D `scan.py` と同じ HarmonicBias を使用。
- 収束は LBFGS/RFO のいずれか（--opt-mode）。ステップ／トラスト半径は Å 指定上限を Bohr に変換して上限化。
- `--baseline min|first`:
  - `min`  : グリッド最小値を 0 kcal/mol として相対化（デフォルト）
  - `first`: (i=0, j=0) の点を 0 kcal/mol として相対化
- 大量格子の場合の PDB 書き出しはファイル数が膨大になるため、既定では .xyz のみ。
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

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, BOHR2ANG

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
    "dump": False,            # per-grid での巨大化を避けるため、デフォルトは False
    "max_cycles": 100,        # 1 ステップの緩和を軽量に（必要に応じて --relax-max-cycles で上書き）
})

LBFGS_KW: Dict[str, Any] = dict(_LBFGS_KW)
LBFGS_KW.update({"out_dir": "./result_scan2d/"})

RFO_KW: Dict[str, Any] = dict(_RFO_KW)
RFO_KW.update({"out_dir": "./result_scan2d/"})

# ハーモニック拘束の既定（強め：グリッド点で距離を保つため）
BIAS_KW: Dict[str, Any] = {"k": 100.0}  # eV/Å^2

# opt-mode エイリアス
_OPT_MODE_ALIASES = (
    (("light", "lbfgs"), "lbfgs"),
    (("heavy", "rfo"),   "rfo"),
)

# 単位変換
HARTREE_TO_KCAL_MOL = 627.50961


# ---- 小物ユーティリティ ----

def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _snapshot_geometry(g) -> Any:
    """
    pysisyphus.Geometry をスナップショット。
    （既存 scan.py と同様：一時 XYZ 経由で独立ジオメトリを得る）
    """
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
    """
    --scan-list で渡された Python 風リスト:
      "[(i1,j1,low1,high1),(i2,j2,low2,high2)]"
    を 0-based に正規化して返す。
    """
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
    """
    [low, high] を最大ステップ h で均等分割。
    差が 0 の場合は長さ 1 の配列 [low] を返す。
    """
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


def _try_write_png(fig, path_png: Path, width: int, height: int) -> Optional[Path]:
    """
    kaleido があれば PNG を出力。無ければ None を返す。
    """
    try:
        fig.write_image(str(path_png), scale=2, engine="kaleido", width=width, height=height)
        return path_png
    except Exception as e:
        click.echo(f"[plot] PNG export failed (kaleido not available?): {e}", err=True)
        return None


# --- B 案の良点を取り込んだ: unbiased SP エネルギー取得の安全関数 ---
def _unbiased_energy_hartree(geom, base_calc) -> float:
    """
    現在座標で UMA の非拘束エネルギー（Hartree）を評価する（属性差異に頑健）。
    """
    coords_bohr = np.asarray(geom.coords3d, dtype=float).reshape(-1, 3)
    for name in ("elem", "elements", "elems"):
        elems = getattr(geom, name, None)
        if elems is not None:
            try:
                return float(base_calc.get_energy(elems, coords_bohr)["energy"])
            except Exception:
                pass
    # 最後の手段: Geometry 側の計算器を base に張替えて取得
    try:
        geom.set_calculator(base_calc)
        elems = getattr(geom, "elem", getattr(geom, "elements", getattr(geom, "elems", None)))
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception:
        return float("nan")


# ---- CLI 本体 ----

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
              help="Pre-optimize initial structure without bias before the scan.")
@click.option("--baseline", type=click.Choice(["min","first"]), default="min", show_default=True,
              help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0).")
# --- B 案から取り込み：プロットのカラースケールを CLI から制御 ---
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

    # charge/spin の解決（CLI > ファイル）
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

    try:
        time_start = time.perf_counter()

        # -------------------------
        # 1) 構成の組み立て
        # -------------------------
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

        # CLI override
        calc_cfg["charge"]  = int(charge)
        calc_cfg["spin"]    = int(spin)
        opt_cfg["out_dir"]  = out_dir
        opt_cfg["dump"]     = False  # per-grid の巨大化防止（--dump は内側のみ）
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

        # Echo config
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_opt  = dict(opt_cfg); echo_opt["out_dir"] = str(out_dir_path)
        echo_bias = dict(bias_cfg)
        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("opt",  echo_opt))
        click.echo(pretty_block("lbfgs" if kind == "lbfgs" else "rfo", (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))
        click.echo(pretty_block("bias", echo_bias))

        # -------------------------
        # 2) scan-list の解釈
        # -------------------------
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

        # -------------------------
        # 3) Geometry 読み込み & 計算器
        # -------------------------
        _ensure_dir(out_dir_path)
        grid_dir  = out_dir_path / "grid"
        plots_dir = out_dir_path / "plots"
        _ensure_dir(grid_dir)
        _ensure_dir(plots_dir)

        coord_type = geom_cfg.get("coord_type", "cart")
        geom_outer = geom_loader(geom_input_path, coord_type=coord_type)

        # freeze-atoms（PDB のリンク親の凍結をマージ）
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

        # UMA 計算器
        base_calc = uma_pysis(**calc_cfg)
        biased    = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

        # 任意の preopt（無拘束）
        if preopt:
            click.echo("[preopt] Unbiased relaxation of the initial structure ...")
            geom_outer.set_calculator(base_calc)
            max_step_bohr_local = float(max_step_size) * ANG2BOHR
            optimizer0 = _make_optimizer(
                geom_outer, kind, lbfgs_cfg, rfo_cfg, opt_cfg,
                max_step_bohr=max_step_bohr_local, relax_max_cycles=relax_max_cycles,
                out_dir=out_dir_path, prefix="preopt_"
            )
            try:
                optimizer0.run()
            except ZeroStepLength:
                click.echo("[preopt] ZeroStepLength — continuing.", err=True)
            except OptimizationError as e:
                click.echo(f"[preopt] OptimizationError — {e}", err=True)

        # ステップ上限（Bohr）
        max_step_bohr = float(max_step_size) * ANG2BOHR

        # -------------------------
        # 4) ネスト走査 & エネルギー収集
        # -------------------------
        records: List[Dict[str, Any]] = []
        # 外側：d1（まず d1 の拘束のみで位置合わせ）
        for i_idx, d1_target in enumerate(d1_values):
            click.echo(f"\n--- d1 step {i_idx+1}/{N1} : target = {d1_target:.3f} Å ---")
            # d1 のみ拘束
            geom_outer.set_calculator(biased)
            biased.set_pairs([(i1, j1, float(d1_target))])
            geom_outer.set_calculator(biased)  # flush caches

            opt1 = _make_optimizer(
                geom_outer, kind, lbfgs_cfg, rfo_cfg, opt_cfg,
                max_step_bohr=max_step_bohr, relax_max_cycles=relax_max_cycles,
                out_dir=out_dir_path, prefix=f"d1_{i_idx:03d}_"
            )
            try:
                opt1.run()
            except ZeroStepLength:
                click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2 scan.", err=True)
            except OptimizationError as e:
                click.echo(f"[d1 {i_idx}] OptimizationError — {e}", err=True)

            # 内側：d2（d1 を固定し、d2 を走査）
            geom_inner = _snapshot_geometry(geom_outer)
            geom_inner.set_calculator(biased)

            trj_blocks = [] if dump else None

            for j_idx, d2_target in enumerate(d2_values):
                biased.set_pairs([(i1, j1, float(d1_target)),
                                  (i2, j2, float(d2_target))])
                geom_inner.set_calculator(biased)  # flush

                opt2 = _make_optimizer(
                    geom_inner, kind, lbfgs_cfg, rfo_cfg, opt_cfg,
                    max_step_bohr=max_step_bohr, relax_max_cycles=relax_max_cycles,
                    out_dir=out_dir_path, prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}_"
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

                # その点での UMA 単発エネルギー（無拘束）を評価（B 案の堅牢版関数を使用）
                E_h = _unbiased_energy_hartree(geom_inner, base_calc)

                # 書き出し
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
                    # 軽量：ステップの最終スナップのみを積む
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

        # -------------------------
        # 5) surface.csv に保存 & 相対化
        # -------------------------
        df = pd.DataFrame.from_records(records)
        if df.empty:
            click.echo("No grid records produced; aborting.", err=True)
            sys.exit(1)

        # 相対化（kcal/mol）
        if baseline == "first":
            ref = float(df.loc[(df["i"]==0) & (df["j"]==0), "energy_hartree"].iloc[0])
        else:
            ref = float(df["energy_hartree"].min())
        df["energy_kcal"] = (df["energy_hartree"] - ref) * HARTREE_TO_KCAL_MOL

        surface_csv = out_dir_path / "surface.csv"
        df.to_csv(surface_csv, index=False)
        click.echo(f"[write] Wrote '{surface_csv}'.")

        # -------------------------
        # 6) プロット（2D 等高線 & 3D + 投影）
        # -------------------------
        try:
            # グリッド形状にピボット（不足があれば RBF 補間にフォールバック）
            piv = df.pivot(index="d2_A", columns="d1_A", values="energy_kcal")
            x_vals = np.array(piv.columns, dtype=float)
            y_vals = np.array(piv.index, dtype=float)
            Z = piv.values  # shape (len(y), len(x))

            # 欠損があれば RBF 補間（任意・フォールバック）
            if np.isnan(Z).any():
                click.echo("[plot] Missing grid points detected; filling with RBF interpolation.")
                try:
                    from scipy.interpolate import Rbf
                    # 有効点のみ
                    dfe = df[np.isfinite(df["energy_kcal"])]
                    rbf = Rbf(dfe["d1_A"], dfe["d2_A"], dfe["energy_kcal"], function="multiquadric")
                    XI, YI = np.meshgrid(x_vals, y_vals)
                    Z = rbf(XI, YI)
                except Exception as e:
                    click.echo(f"[plot] RBF interpolation failed: {e}", err=True)
                    # NaN を単純埋め
                    Z = np.nan_to_num(Z, nan=np.nanmean(Z))

            # --- B 案の良点：zmin/zmax を CLI で上書き可能に ---
            zmin_eff = float(np.nanmin(Z)) if zmin is None else float(zmin)
            zmax_eff = float(np.nanmax(Z)) if zmax is None else float(zmax)

            # 2D 等高線（正方形 + 軸スケール固定）
            import plotly.graph_objects as go

            fig2d = go.Figure(data=go.Contour(
                z=Z,
                x=x_vals,  # d1
                y=y_vals,  # d2
                contours=dict(
                    start=float(np.floor(zmin_eff/10.0)*10.0),
                    end=float(np.ceil(zmax_eff/10.0)*10.0),
                    size=10.0
                ),
                zmin=zmin_eff, zmax=zmax_eff,
                contours_coloring="heatmap",
                colorscale="plasma",
                colorbar=dict(
                    title=dict(text="(kcal/mol)", side="top"),
                    ticks="inside",
                ),
            ))
            # 正方形に（B 案取り込み：scaleanchor）
            fig2d.update_layout(
                width=720, height=720,
                xaxis_title=f"d1 ({i1+1 if one_based else i1},{j1+1 if one_based else j1}) distance (Å)",
                yaxis_title=f"d2 ({i2+1 if one_based else i2},{j2+1 if one_based else j2}) distance (Å)",
                plot_bgcolor="white",
                xaxis=dict(
                    range=[min(x_vals), max(x_vals)],
                    showline=True, linewidth=3, linecolor="#1C1C1C", mirror=True,
                    ticks="inside", tickwidth=3, tickcolor="#1C1C1C",
                ),
                yaxis=dict(
                    range=[min(y_vals), max(y_vals)],
                    showline=True, linewidth=3, linecolor="#1C1C1C", mirror=True,
                    ticks="inside", tickwidth=3, tickcolor="#1C1C1C",
                ),
                margin=dict(l=10, r=10, t=40, b=10),
            )
            # 軸スケール固定で真の正方
            fig2d.update_yaxes(scaleanchor="x", scaleratio=1)

            png2d = _try_write_png(fig2d, plots_dir / "scan2d_contour.png", width=720, height=720)
            if png2d:
                click.echo(f"[plot] Wrote '{png2d}'.")
            else:
                html2d = plots_dir / "scan2d_contour.html"
                fig2d.write_html(str(html2d))
                click.echo(f"[plot] Wrote '{html2d}' (HTML fallback).")

            # 3D サーフェス + 底面投影（立方体比）
            XI, YI = np.meshgrid(x_vals, y_vals)
            surface = go.Surface(
                x=XI, y=YI, z=Z,
                colorscale="plasma",
                cmin=zmin_eff, cmax=zmax_eff,
                colorbar=dict(
                    title=dict(text="(kcal/mol)", side="top"),
                    ticks="inside",
                ),
                contours={
                    "z": {
                        "show": True,
                        "start": float(np.floor(zmin_eff/10.0)*10.0),
                        "end": float(np.ceil(zmax_eff/10.0)*10.0),
                        "size": 10.0,
                        "color": "black",
                        "project": {"z": True},
                    }
                },
                name="3D Surface",
            )

            contour_proj = go.Surface(
                x=XI, y=YI, z=np.full_like(Z, zmin_eff - 1.0),  # 少し下に投影
                surfacecolor=Z, showscale=False, colorscale="plasma",
                cmin=zmin_eff, cmax=zmax_eff, opacity=1.0, name="2D Projection"
            )

            fig3d = go.Figure(data=[surface, contour_proj])
            fig3d.update_layout(
                title="2D PES with bottom projection",
                width=800, height=800,
                scene=dict(
                    bgcolor="rgba(0,0,0,0)",
                    xaxis=dict(title=f"d1 ({i1+1 if one_based else i1},{j1+1 if one_based else j1}) (Å)",
                               showline=True, linewidth=4, linecolor="#1C1C1C", mirror=True,
                               ticks="inside", tickwidth=4, tickcolor="#1C1C1C",
                               showbackground=False),
                    yaxis=dict(title=f"d2 ({i2+1 if one_based else i2},{j2+1 if one_based else j2}) (Å)",
                               showline=True, linewidth=4, linecolor="#1C1C1C", mirror=True,
                               ticks="inside", tickwidth=4, tickcolor="#1C1C1C",
                               showbackground=False),
                    zaxis=dict(title="Energy (kcal/mol)",
                               showline=True, linewidth=4, linecolor="#1C1C1C", mirror=True,
                               ticks="inside", tickwidth=4, tickcolor="#1C1C1C",
                               showbackground=False),
                    aspectmode="cube",   # 立方体（x,y,z を等尺に）
                ),
                margin=dict(l=10, r=10, t=40, b=10),
                paper_bgcolor="white",
            )

            html3d = plots_dir / "scan2d_surface.html"
            fig3d.write_html(str(html3d))
            click.echo(f"[plot] Wrote '{html3d}'.")

        except Exception as e:
            click.echo(f"[plot] WARNING: plotting failed: {e}", err=True)

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
