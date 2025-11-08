# pdb2reaction/path_search.py
"""
Multistep MEP search by recursive GSM segmentation.

概要
----
2 つの端点構造（reactant, product）から開始し，
(1) まず端点間に GSM を 1 回走らせて MEP を取得。
(2) 最高エネルギー画像（HEI）の左右 1 画像ずつを単独最適化して近接極小（End1, End2）を得る。
(3) End1–End2 間で再度 GSM（精緻化）を行い，1「段階（step）」の MEP を確定。
(4) A（reactant）と End1，End2 と B（product）の各対で共有結合の形成/解離を判定し，
    結合変化がある側に対して同手順を再帰的に適用。
(5) 得られた各「段階」MEP を，重複排除・必要に応じたブリッジ GSM の自動挿入を行いながら連結し，
    reactant→product の連続 MEP を構築・書き出す。

主な仕様
--------
- UMA 計算器（uma_pysis）を全工程で共有（直列実行）。
- GSM は pysisyphus の GrowingString + StringOptimizer を使用。
- HEI の左右画像は単独最適化（LBFGS / RFO を選択）。
- 共有結合変化の検出に `bond_changes.compare_structures` を使用。
- 「段階」MEP の連結時：
  - 端点のズレが閾値を超える場合はブリッジ GSM を自動実行。
  - **前ステップ末端–次ステップ先頭に共有結合変化がある場合は，
    ブリッジではなく新規セグメントを自動生成して挿入（再帰探索）。**
- **入力 2 構造は既定で Kabsch により事前アライン（StringOptimizer の align は使用しない）。**
- **freeze_atoms がある系では，Kabsch の剛体整列は freeze_atoms の原子座標のみで最適化し，
  得た剛体変換を全原子に適用。**
- **入力 2 構造はアライン前に単独最適化（sopt-mode）を実施。**
  **ただし `--pre_opt True` / `--pre-opt True` 指定時はこの事前最適化をスキップ。**
- YAML（geom, calc, gs, opt, sopt, bond, search）と CLI の優先順位は既存ツールと同じ（CLI > YAML > defaults）。
- **セグメント用とブリッジ用で max_nodes を分離設定可能
  （search.max_nodes_segment / search.max_nodes_bridge）。**
  ブリッジ GSM は強制的に `climb=False`（ただし max_depth 到達時の「つなぎ」GSM は従来通り `climb=True`）。
- **初回のセグメント GSM は `climb=False`，End1–End2 の再 GSM（refine）は `climb=True`。**

出力
----
out_dir/
  ├─ summary.yaml                 : 実行サマリ
  ├─ final_geometries.trj         : 連結 MEP（各ブロック2行目にエネルギー）
  ├─ final_geometries.pdb         : 入力が PDB の場合に生成
  └─ segments/
      ├─ seg_000_gsm/ ...         : 初回 GSM
      ├─ seg_000_left_opt/ ...    : HEI-1 の単独最適化
      ├─ seg_000_right_opt/ ...   : HEI+1 の単独最適化
      ├─ seg_000_refine_gsm/ ...  : End1–End2 の再 GSM（精緻化）
      ├─ seg_001_...              : 再帰で発生する左側のサブステップ
      └─ seg_002_...              : 再帰で発生する右側のサブステップ
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Callable  # 型ヒントで Callable を使用

import sys
import traceback
import textwrap
import tempfile
import os

import click
import numpy as np
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength

from .uma_pysis import uma_pysis
from .utils import convert_xyz_to_pdb, freeze_links

from .bond_changes import compare_structures, summarize_changes

# -----------------------------------------------
# Defaults
# -----------------------------------------------

# 幾何（入力の扱い）
GEOM_KW: Dict[str, Any] = {
    "coord_type": "cart",   # GrowingString は Cartesian 推奨
    "freeze_atoms": [],     # 0-based indices
}

# UMA 計算条件
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

# GrowingString（経路表現）
GS_KW: Dict[str, Any] = {
    "max_nodes": 30,            # 端点含めて max_nodes+2 画像
    "perp_thresh": 5e-3,
    "reparam_check": "rms",     # "rms" | "norm"
    "reparam_every": 1,
    "reparam_every_full": 1,
    "param": "equi",            # "equi" | "energy"
    "max_micro_cycles": 10,
    "reset_dlc": True,
    "climb": True,              # 初回セグメントでは後で強制的に False にする
    "climb_rms": 5e-4,
    "climb_lanczos": True,
    "climb_lanczos_rms": 5e-4,
    "climb_fixed": False,
    "scheduler": None,          # 直列計算（同一 calculator を共有）
}

# StringOptimizer（GSM の最適化制御）
OPT_KW: Dict[str, Any] = {
    "type": "string",
    "stop_in_when_full": 1000,  # fully grown 後の停止猶予
    "align": False,             # StringOptimizer の align は使用しない
    "scale_step": "global",     # global | per_image
    "max_cycles": 1000,
    "dump": False,
    "dump_restart": False,
    "reparam_thresh": 1e-3,
    "coord_diff_thresh": 0.0,
    "out_dir": "./result_path_search/",
    "print_every": 1,
}

# 単独最適化（LBFGS/RFO）の共通
SOPT_BASE_KW: Dict[str, Any] = {
    "thresh": "gau",            # 収束閾値プリセット
    "max_cycles": 10000,
    "print_every": 1,
    "min_step_norm": 1e-8,
    "assert_min_step": True,
    "dump": False,
    "dump_restart": False,
    "prefix": "",
    "out_dir": "./result_path_search/",
}

# LBFGS 用
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

# RFO 用
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

# 共有結合変化の判定パラメタ
BOND_KW: Dict[str, Any] = {
    "device": "cuda",
    "bond_factor": 1.20,
    "margin_fraction": 0.05,
    "delta_fraction": 0.05,
}

# 探索全体の制御
SEARCH_KW: Dict[str, Any] = {
    "max_depth": 10,               # 再帰の最大段数
    "stitch_rmsd_thresh": 1.0e-4,  # 連結時に重複とみなす RMSD（Å）
    "bridge_rmsd_thresh": 1.0e-4,  # 端点ズレがこれを超えたらブリッジ GSM 実行
    "rmsd_align": True,            # RMSD 計算で Kabsch アラインを使用
    "max_nodes_segment": 30,
    "max_nodes_bridge": 10,
}

# -----------------------------------------------
# Utilities
# -----------------------------------------------

def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """dict *dst* を *src* で再帰的に上書きして返す。"""
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    """YAML を読み込んで dict を返す（未指定なら空 dict）。"""
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")
    return data


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    """設定の見やすいダンプ文字列を生成。"""
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """freeze_atoms を表示用に整形。"""
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple, np.ndarray)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if len(fa) else ""
    return g


def _freeze_links_for_pdb(pdb_path: Path) -> Sequence[int]:
    """PDB の link 水素の親原子を検出し 0-based index を返す。"""
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
    """2 つの端点構造を読み込み，必要に応じて freeze_atoms を設定。"""
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


def _ensure_calc_on_geom(g, calc) -> None:
    """pysisyphus Geometry に Calculator を付与（必要なら上書き）。"""
    try:
        g.set_calculator(calc)
    except Exception:
        g.set_calculator(calc)


def _write_xyz_trj_with_energy(images: Sequence, energies: Sequence[float], path: Path) -> None:
    """画像列 + エネルギーから XYZ 形式の .trj を出力（各ブロック 2 行目に E）。"""
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
    入力が PDB の場合のみ，指定された .xyz / .trj を PDB に変換して保存。
    成功時は出力パスを返す（不要・失敗時は None）。
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
    """Kabsch アライン後 RMSD（align=False なら単純 RMSD）。indices 指定で部分集合 RMSD。"""
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
    # 反転の取り扱い
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    Arot = Ac @ R
    diff = Arot - Bc
    return float(np.sqrt((diff * diff).sum() / A.shape[0]))


def _rmsd_between(ga, gb, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """Geometry 間の RMSD（必要に応じ Kabsch アライン・部分集合）。"""
    return _kabsch_rmsd(np.array(ga.coords3d), np.array(gb.coords3d), align=align, indices=indices)


def _orient_two_by_closeness(a, b, end1, end2, align: bool = True) -> Tuple:
    """端点 a,b に対して，a に近い方を left，b に近い方を right に並べ替える。"""
    r1a = _rmsd_between(end1, a, align)
    r2a = _rmsd_between(end2, a, align)
    if r1a <= r2a:
        left, right = end1, end2
    else:
        left, right = end2, end1
    return left, right


def _has_bond_change(x, y, bond_cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """共有結合の形成/解離の有無を判定（True/False, 要約文字列）。"""
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


# ---- Kabsch によるアライン（mob を ref に合わせるクローンを作成；部分集合対応） ----

def _kabsch_aligned_clone(ref_geom, mob_geom, coord_type: str, use_indices: Optional[Sequence[int]] = None) -> Any:
    """
    ref_geom に対して mob_geom を Kabsch で最適回転・並進し，
    新しい Geometry を生成して返す（元の mob_geom は不変）。
    use_indices 指定時はその部分集合のみで剛体変換を決定（変換自体は全原子に適用）。
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
        # freeze を引き継ぐ
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


# ---------- path_opt.py と同等の堅牢な Kabsch 実装（最小限） ----------

def _kabsch_R_t(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Q を P に最小二乗で合わせる最適回転 R と並進 t（Kabsch）を求める。

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
    # 右手系を保証
    if np.linalg.det(R) < 0.0:
        Vt[-1, :] *= -1.0
        R = Vt.T @ U.T
    t = mu_P - mu_Q @ R
    return R, t


def _align_second_to_first_kabsch(geom_ref, geom_to_align) -> Tuple[float, float, int]:
    """
    geom_to_align を geom_ref に Kabsch で剛体整列（スケーリング無し）。
    いずれかの端点が freeze_atoms を持つ場合は，その **み** を用いて変換を決定し，
    変換は全原子に適用する。

    Returns
    -------
    rmsd_before, rmsd_after, n_used  （RMSD は選択集合上で評価）
    """
    P = np.array(geom_ref.coords3d, dtype=float)    # (N, 3)
    Q = np.array(geom_to_align.coords3d, dtype=float)
    if P.shape != Q.shape:
        raise ValueError(f"Different atom counts for endpoints: {P.shape[0]} vs {Q.shape[0]}")

    N = P.shape[0]
    # 変換決定に用いるインデックス
    fa0 = getattr(geom_ref, "freeze_atoms", np.array([], dtype=int))
    fa1 = getattr(geom_to_align, "freeze_atoms", np.array([], dtype=int))
    freeze_union = sorted(set(map(int, fa0)) | set(map(int, fa1)))

    if len(freeze_union) > 0:
        use_mask = np.zeros(N, dtype=bool)
        # 範囲外を防御的に除外
        valid_idx = [i for i in freeze_union if 0 <= i < N]
        use_mask[valid_idx] = True
    else:
        use_mask = np.ones(N, dtype=bool)

    P_sel = P[use_mask]
    Q_sel = Q[use_mask]
    n_used = int(P_sel.shape[0])

    # 整列前 RMSD（選択集合上；追加アラインなし）
    def _rmsd(A: np.ndarray, B: np.ndarray) -> float:
        return float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1)))) if len(A) else float("nan")

    rmsd_before = _rmsd(P_sel, Q_sel)

    R, t = _kabsch_R_t(P_sel, Q_sel)
    Q_aligned = (Q @ R) + t  # 全原子に適用

    # set_coords は 1 次元（長さ 3N）を受け取る
    geom_to_align.set_coords(Q_aligned.reshape(-1))

    rmsd_after = _rmsd(P[use_mask], Q_aligned[use_mask])

    return rmsd_before, rmsd_after, n_used
# ---------- 追加ここまで ----------


# ---------- GS 設定ユーティリティ（最小限） ----------
def _gs_cfg_with_overrides(base: Dict[str, Any], **overrides: Any) -> Dict[str, Any]:
    """GrowingString 用設定の浅いコピーに上書きを適用して返す。"""
    cfg = dict(base)
    for k, v in overrides.items():
        cfg[k] = v
    return cfg
# -----------------------------------------------


@dataclass
class GSMResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int


def _run_gsm_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # PDB 変換用の参照 PDB
) -> GSMResult:
    """gA–gB 間で GSM を実行し，セグメント出力として保存。"""
    # Calculator を端点に付与
    for g in (gA, gB):
        _ensure_calc_on_geom(g, shared_calc)

    def calc_getter():
        return shared_calc

    # GrowingString のセットアップ
    gs = GrowingString(
        images=[gA, gB],
        calc_getter=calc_getter,
        **gs_cfg,
    )

    # Optimizer 準備
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
    # 書き出し
    final_trj = seg_dir / "final_geometries.trj"
    try:
        _write_xyz_trj_with_energy(images, energies, final_trj)
        click.echo(f"[{tag}] Wrote '{final_trj}'.")
    except Exception:
        with open(final_trj, "w") as f:
            f.write(gs.as_xyz())
        click.echo(f"[{tag}] Wrote '{final_trj}'.")

    # 入力が PDB の場合，中間 .trj を PDB に変換
    _maybe_convert_to_pdb(final_trj, ref_pdb_path, seg_dir / "final_geometries.pdb")

    # HEI 書き出し
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
        # HEI も PDB へ
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
    ref_pdb_path: Optional[Path],  # PDB 変換用
):
    """単独構造最適化（LBFGS/RFO）。最終 Geometry を返す。"""
    # Calculator を付与
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

    # Optimizer は最終構造をファイル出力（final_fn）
    try:
        final_xyz = Path(opt.final_fn) if isinstance(opt.final_fn, (str, Path)) else Path(opt.final_fn)
        # 入力が PDB の場合，最終 .xyz/.trj を PDB へ変換
        _maybe_convert_to_pdb(final_xyz, ref_pdb_path)
        g_final = geom_loader(final_xyz, coord_type=g.coord_type)
        # freeze を引き継ぐ
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        _ensure_calc_on_geom(g_final, shared_calc)
        return g_final
    except Exception:
        # 失敗時は現在の g を返す
        return g


def _refine_between(
    gL,
    gR,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # PDB 変換用
) -> GSMResult:
    """End1–End2 の再 GSM（精緻化）。ここでは climb を強制的に有効化。"""
    gs_refine_cfg = _gs_cfg_with_overrides(gs_cfg, climb=True, climb_lanczos=True)
    return _run_gsm_between(gL, gR, shared_calc, gs_refine_cfg, opt_cfg, out_dir, tag=f"{tag}_refine", ref_pdb_path=ref_pdb_path)


def _maybe_bridge_segments(
    tail_g,
    head_g,
    shared_calc,
    gs_cfg: Dict[str, Any],  # ブリッジ用設定を受け取る
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    rmsd_thresh: float,
    ref_pdb_path: Optional[Path],  # PDB 変換用
) -> Optional[GSMResult]:
    """2 セグメント端点のズレが閾値超ならブリッジ GSM を実行。"""
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
    gs_cfg,   # ブリッジ用 GrowingString 設定（climb=False, max_nodes=search.max_nodes_bridge）
    opt_cfg,
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # PDB 変換用
    bond_cfg: Optional[Dict[str, Any]] = None,  # 隣接セグメント間の結合変化検出に使用
    segment_builder: Optional[Callable[[Any, Any, str], Tuple[List[Any], List[float]]]] = None,  # 結合変化時の新規セグメント生成
) -> Tuple[List[Any], List[float]]:
    """
    パート列（images, energies）を連結する。必要に応じてブリッジ GSM を挿入。
    さらに，直前末尾と今回先頭の間で共有結合変化が検出された場合は，
    `segment_builder` により新規セグメントを生成・挿入する。
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
        # 直前末尾と今回先頭
        tail = all_imgs[-1]
        head = imgs[0]

        # 隣接端点間の共有結合変化チェック
        adj_changed, adj_summary = False, ""
        if segment_builder is not None and bond_cfg is not None:
            try:
                adj_changed, adj_summary = _has_bond_change(tail, head, bond_cfg)
            except Exception:
                adj_changed, adj_summary = False, ""

        if adj_changed:
            click.echo(f"[{tag}] Bond-change detected between segment endpoints — inserting a new segment via recursive GSM.")
            if adj_summary:
                click.echo(textwrap.indent(adj_summary, prefix="  "))
            # 新規セグメント（tail→head）を生成・挿入
            seg_imgs, seg_E = segment_builder(tail, head, f"{tag}_mid")
            if seg_imgs:
                # 既存末尾との重複を除去
                if _rmsd_between(all_imgs[-1], seg_imgs[0], align=True) <= stitch_rmsd_thresh:
                    seg_imgs = seg_imgs[1:]
                    seg_E = seg_E[1:]
                all_imgs.extend(seg_imgs)
                all_E.extend(seg_E)
            # 続けて今回パートを結合（重複処理）
            if _rmsd_between(all_imgs[-1], imgs[0], align=True) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return  # このパートの処理は完了

        rmsd = _rmsd_between(tail, head, align=True)
        if rmsd <= stitch_rmsd_thresh:
            # 完全重複 → 先頭を落として結合
            all_imgs.extend(imgs[1:])
            all_E.extend(Es[1:])
        elif rmsd > bridge_rmsd_thresh:
            # ズレが大 → ブリッジ GSM
            br = _maybe_bridge_segments(
                tail, head, shared_calc, gs_cfg, opt_cfg, out_dir, tag=tag,
                rmsd_thresh=bridge_rmsd_thresh, ref_pdb_path=ref_pdb_path
            )
            if br is not None:
                # ブリッジ先頭の重複を除去
                b_imgs, b_E = br.images, br.energies
                if _rmsd_between(all_imgs[-1], b_imgs[0], align=True) <= stitch_rmsd_thresh:
                    b_imgs = b_imgs[1:]
                    b_E = b_E[1:]
                if b_imgs:
                    all_imgs.extend(b_imgs)
                    all_E.extend(b_E)
            # その後，今回パートを連結（改めて重複処理）
            if _rmsd_between(all_imgs[-1], imgs[0], align=True) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
        else:
            # 軽微なズレ → そのまま接続
            all_imgs.extend(imgs)
            all_E.extend(Es)

    for (imgs, Es) in parts:
        append_part(imgs, Es)

    return all_imgs, all_E


# -----------------------------------------------
# 再帰探索（本体）
# -----------------------------------------------

@dataclass
class CombinedPath:
    images: List[Any]
    energies: List[float]


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
    A–B の多段階 MEP を再帰的に構築して返す（A→B の順）。
    """
    # セグメント用 GS 設定（search.max_nodes_segment を適用）
    seg_max_nodes = int(search_cfg.get("max_nodes_segment", gs_cfg.get("max_nodes", 30)))
    gs_seg_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=seg_max_nodes)

    if depth > int(search_cfg.get("max_depth", 10)):
        click.echo(f"[{branch_tag}] Reached max recursion depth. Returning current endpoints only.")
        # 上限到達時は，A–B を単純に GSM（ここは climb=True のまま）
        gsm = _run_gsm_between(gA, gB, shared_calc, gs_seg_cfg, opt_cfg, out_dir, tag=f"seg_{seg_counter[0]:03d}_maxdepth", ref_pdb_path=ref_pdb_path)
        seg_counter[0] += 1
        return CombinedPath(images=gsm.images, energies=gsm.energies)

    seg_id = seg_counter[0]
    seg_counter[0] += 1
    tag0 = f"seg_{seg_id:03d}"

    # 1) A–B 初回 GSM（セグメント設定）
    #    初回は kink 抑制のため climb を強制 OFF
    gs_seg_cfg_first = _gs_cfg_with_overrides(gs_seg_cfg, climb=False, climb_lanczos=False)
    gsm0 = _run_gsm_between(gA, gB, shared_calc, gs_seg_cfg_first, opt_cfg, out_dir, tag=tag0, ref_pdb_path=ref_pdb_path)

    # HEI とその両隣
    hei = int(gsm0.hei_idx)
    if not (1 <= hei <= len(gsm0.images) - 2):
        # 端に寄った場合はそのまま返す
        click.echo(f"[{tag0}] WARNING: HEI is at an endpoint (idx={hei}). Returning raw GSM.")
        return CombinedPath(images=gsm0.images, energies=gsm0.energies)

    left_img = gsm0.images[hei - 1]
    right_img = gsm0.images[hei + 1]

    # 2) 両隣画像の単独最適化
    left_opt = _optimize_single(left_img, shared_calc, sopt_kind, sopt_cfg, out_dir, tag=f"{tag0}_left", ref_pdb_path=ref_pdb_path)
    right_opt = _optimize_single(right_img, shared_calc, sopt_kind, sopt_cfg, out_dir, tag=f"{tag0}_right", ref_pdb_path=ref_pdb_path)

    # 3) 端点の向き（A に近い方を left，B に近い方を right）
    left_end, right_end = _orient_two_by_closeness(gA, gB, left_opt, right_opt, align=bool(search_cfg.get("rmsd_align", True)))

    # 4) End1–End2 の再 GSM（精緻化；セグメント設定）
    ref1 = _refine_between(left_end, right_end, shared_calc, gs_seg_cfg, opt_cfg, out_dir, tag=tag0, ref_pdb_path=ref_pdb_path)
    step_imgs, step_E = ref1.images, ref1.energies

    # 5) 共有結合変化の有無（左：A–left_end，右：right_end–B）
    left_changed, left_summary = _has_bond_change(gA, left_end, bond_cfg)
    right_changed, right_summary = _has_bond_change(right_end, gB, bond_cfg)

    click.echo(f"[{tag0}] Bond-change(A vs left_end): {'Yes' if left_changed else 'No'}")
    if left_changed:
        click.echo(textwrap.indent(left_summary, prefix="  "))
    click.echo(f"[{tag0}] Bond-change(right_end vs B): {'Yes' if right_changed else 'No'}")
    if right_changed:
        click.echo(textwrap.indent(right_summary, prefix="  "))

    # 6) 再帰（左・右）
    parts: List[Tuple[List[Any], List[float]]] = []

    if left_changed:
        subL = _build_multistep_path(
            gA, left_end, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg, bond_cfg, search_cfg,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}L"
        )
        parts.append((subL.images, subL.energies))

    # 中央（今回の段階）
    parts.append((step_imgs, step_E))

    if right_changed:
        subR = _build_multistep_path(
            right_end, gB, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg, bond_cfg, search_cfg,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}R"
        )
        parts.append((subR.images, subR.energies))

    # 7) 連結（重複除去・必要に応じブリッジ/新規セグメント挿入）
    #    ブリッジ用の GS 設定（max_nodes_bridge, climb=False）
    bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 10))
    gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

    # tail–head 間に結合変化がある場合は，再帰で新規セグメントを構築するビルダ
    def _segment_builder(tail_g, head_g, _tag: str) -> Tuple[List[Any], List[float]]:
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
        return sub.images, sub.energies

    stitched_imgs, stitched_E = _stitch_paths(
        parts,
        stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-3)),
        bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-2)),
        shared_calc=shared_calc,
        gs_cfg=gs_bridge_cfg,   # ブリッジは climb=False・専用 max_nodes
        opt_cfg=opt_cfg,
        out_dir=out_dir,
        tag=tag0,
        ref_pdb_path=ref_pdb_path,
        bond_cfg=bond_cfg,
        segment_builder=_segment_builder
    )

    return CombinedPath(images=stitched_imgs, energies=stitched_E)


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
    nargs=2,
    required=True,
    help="Two endpoint structures (reactant and product), e.g. .pdb / .xyz",
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
    help="If True, SKIP initial single-structure optimizations of endpoints (i.e., do not pre-opt). "
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
    try:
        # --------------------------
        # 1) 設定の統合（defaults ← YAML ← CLI）
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

        gs_cfg["max_nodes"] = int(max_nodes)           # ベース（主に表示用／後方互換）
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)

        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir

        # 単独最適化の dump/out_dir も揃える
        lbfgs_cfg["dump"] = bool(dump)
        rfo_cfg["dump"]   = bool(dump)
        lbfgs_cfg["out_dir"] = out_dir
        rfo_cfg["out_dir"]   = out_dir

        # sopt kind 正規化
        sopt_kind = sopt_mode.strip().lower()
        if sopt_kind in ("light", "lbfgs"):
            sopt_kind = "lbfgs"
            sopt_cfg = lbfgs_cfg
        elif sopt_kind in ("heavy", "rfo"):
            sopt_kind = "rfo"
            sopt_cfg = rfo_cfg
        else:
            raise click.BadParameter(f"Unknown --sopt-mode '{sopt_mode}'.")

        # segment/bridge の max_nodes を確定
        # - YAML(search.max_nodes_segment) が無ければ CLI --max-nodes を採用
        if "max_nodes_segment" not in yaml_cfg.get("search", {}):
            search_cfg["max_nodes_segment"] = int(max_nodes)
        # - bridge は YAML/既定を尊重（必要なら YAML で search.max_nodes_bridge を指定）

        # 表示用：解決後の設定
        out_dir_path = Path(out_dir).resolve()
        echo_geom = _format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_opt  = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        # pre_opt の表示
        click.echo(_pretty_block("geom", echo_geom))
        click.echo(_pretty_block("calc", echo_calc))
        click.echo(_pretty_block("gs",   echo_gs))
        click.echo(_pretty_block("opt",  echo_opt))
        click.echo(_pretty_block("sopt."+sopt_kind, sopt_cfg))
        click.echo(_pretty_block("bond", bond_cfg))
        click.echo(_pretty_block("search", search_cfg))
        click.echo(_pretty_block("run_flags", {"pre_opt_skip": bool(pre_opt)}))

        # --------------------------
        # 2) 入力準備
        #    (a) 端点読み込み
        #    (b) Calculator 用意
        #    (c) 端点の単独最適化（オプション）
        #    (d) freeze-aware Kabsch 事前アライン
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        p0, p1 = Path(input_paths[0]), Path(input_paths[1])

        geoms = _load_two_endpoints(
            paths=[p0, p1],
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )
        gA, gB = geoms[0], geoms[1]

        # (b) 共通 UMA 計算器（以後の最適化・GSM で共有）
        shared_calc = uma_pysis(**calc_cfg)
        for g in (gA, gB):
            _ensure_calc_on_geom(g, shared_calc)

        # (c) 入力端点の単独最適化（アライン前；pre_opt=True ならスキップ）
        if not pre_opt:
            gA = _optimize_single(gA, shared_calc, sopt_kind, sopt_cfg, out_dir_path, tag="initA", ref_pdb_path=(p0 if p0.suffix.lower() == ".pdb" else (p1 if p1.suffix.lower() == ".pdb" else None)))
            gB = _optimize_single(gB, shared_calc, sopt_kind, sopt_cfg, out_dir_path, tag="initB", ref_pdb_path=(p0 if p0.suffix.lower() == ".pdb" else (p1 if p1.suffix.lower() == ".pdb" else None)))
        else:
            click.echo("[init] Skipping endpoint pre-optimization as requested by --pre_opt/--pre-opt True.")

        # (d) freeze-aware Kabsch pre-alignment（選択集合で変換決定→全原子適用）
        N = int(np.asarray(gA.coords3d).shape[0])
        faA = getattr(gA, "freeze_atoms", np.array([], dtype=int))
        faB = getattr(gB, "freeze_atoms", np.array([], dtype=int))
        use_idx = sorted({int(i) for i in list(faA) + list(faB) if 0 <= int(i) < N})
        idx_for_msg = f"{len(use_idx)} (freeze_atoms)" if len(use_idx) > 0 else f"{N} (all atoms)"

        try:
            # path_opt.py と同等の実装で 2nd → 1st をアライン
            rmsd_before, rmsd_after, _n_used = _align_second_to_first_kabsch(gA, gB)
            click.echo(f"[align] Kabsch pre-alignment (used {idx_for_msg}): "
                       f"RMSD_before={rmsd_before:.6f} Å → RMSD_after={rmsd_after:.6f} Å")
        except Exception as e:
            click.echo(f"[align] WARNING: Kabsch pre-alignment skipped: {e}", err=True)

        # --------------------------
        # 3) 再帰探索の実行
        # --------------------------
        click.echo("\n=== Multistep MEP search started ===\n")
        seg_counter = [0]  # 参照で更新するためリスト

        # 入力が PDB の場合の参照 PDB（以降の中間変換に使用）
        ref_pdb_for_segments = p0 if p0.suffix.lower() == ".pdb" else (p1 if p1.suffix.lower() == ".pdb" else None)

        combined = _build_multistep_path(
            gA, gB,
            shared_calc,
            geom_cfg, gs_cfg, opt_cfg,
            sopt_kind, sopt_cfg,
            bond_cfg, search_cfg,
            out_dir=out_dir_path,
            ref_pdb_path=ref_pdb_for_segments,
            depth=0,
            seg_counter=seg_counter,
            branch_tag="",
        )
        click.echo("\n=== Multistep MEP search finished ===\n")

        # --------------------------
        # 4) 出力
        # --------------------------
        final_trj = out_dir_path / "final_geometries.trj"
        _write_xyz_trj_with_energy(combined.images, combined.energies, final_trj)
        click.echo(f"[write] Wrote '{final_trj}'.")

        # 参照 PDB があれば PDB も出力
        ref_pdb = None
        if p0.suffix.lower() == ".pdb":
            ref_pdb = p0.resolve()
        elif p1.suffix.lower() == ".pdb":
            ref_pdb = p1.resolve()
        if ref_pdb is not None:
            try:
                final_pdb = out_dir_path / "final_geometries.pdb"
                convert_xyz_to_pdb(final_trj, ref_pdb, final_pdb)
                click.echo(f"[convert] Wrote '{final_pdb}'.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final MEP to PDB: {e}", err=True)

        # 実行サマリ
        summary = {
            "n_segments": seg_counter[0],
            "n_images_final": len(combined.images),
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
        }
        with open(out_dir_path / "summary.yaml", "w") as f:
            yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
        click.echo(f"[write] Wrote '{out_dir_path / 'summary.yaml'}'.")

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
