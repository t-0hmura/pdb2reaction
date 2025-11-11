# pdb2reaction/irc.py
# -*- coding: utf-8 -*-
"""
EulerPC による IRC 計算用 CLI

例:
  pdb2reaction irc \
    -i a.pdb -q 0 -s 1 \
    --max-cycles 20 \
    --step-size 0.5 \
    --root 0 \
    --forward True \
    --backward False \
    --out-dir "./result_irc/" \
    --args-yaml args.yaml

方針
-----
- CLI で受け付けるのは上記例にある引数のみ。
- それ以外のパラメータ（幾何, UMA 設定, IRC/EulerPC の詳細設定など）は YAML から読み込む。
- 入力が .pdb の場合は、生成された IRC 全体のトラジェクトリ (finished_irc.trj)、
  forward/backward の各トラジェクトリ (forward_irc.trj, backward_irc.trj) を PDB へ変換する。

YAML セクション例（任意）
------------------------
geom:
  coord_type: cart
  freeze_atoms: []

calc:
  charge: 0
  spin: 1
  model: "uma-s-1p1"
  task_name: "omol"
  device: "auto"
  max_neigh: null
  radius: null
  r_edges: false
  out_hess_torch: false

irc:
  # IRC(基底)の設定
  downhill: false
  forward: true
  backward: true
  hessian_init: "calc"
  displ: "energy"
  displ_energy: 1.0e-3
  displ_length: 0.1
  rms_grad_thresh: 1.0e-3
  hard_rms_grad_thresh: null
  energy_thresh: 1.0e-6
  imag_below: 0.0
  force_inflection: true
  check_bonds: false
  prefix: ""
  dump_fn: "irc_data.h5"
  dump_every: 5

  # EulerPC 固有の設定
  hessian_update: "bofill"
  hessian_recalc: null
  max_pred_steps: 500
  loose_cycles: 3
  corr_func: "mbs"

備考
----
- CLI で与えられた値は YAML より優先します（charge/spin, step-size→step_length,
  max-cycles→max_cycles, root, forward/backward, out-dir）。
- UMA 設定は `pdb2reaction.uma_pysis.uma_pysis` にそのまま渡します。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import sys
import textwrap

import click
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC  # EulerPC(IRC) 実装
from pdb2reaction.uma_pysis import uma_pysis
from pdb2reaction.utils import convert_xyz_to_pdb


# --------------------------
# デフォルト設定
# --------------------------

GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",
    "freeze_atoms": [],  # 0-based indices
}

CALC_KW_DEFAULT: Dict[str, Any] = {
    "charge": 0,
    "spin": 1,
    "model": "uma-s-1p1",
    "task_name": "omol",
    "device": "auto",         # "cuda" | "cpu" | "auto"
    "max_neigh": None,
    "radius": None,           # Å
    "r_edges": False,
    "out_hess_torch": False,  # EulerPCでは numpy 側で良い
}

IRC_KW_DEFAULT: Dict[str, Any] = {
    # IRC.__init__ の引数（EulerPC に **kwargs で渡る）
    "step_length": 0.10,         # CLI --step-size で上書き
    "max_cycles": 125,           # CLI --max-cycles で上書き
    "downhill": False,
    "forward": True,             # CLI --forward で上書き
    "backward": True,            # CLI --backward で上書き
    "root": 0,                   # CLI --root で上書き
    "hessian_init": "calc",      # TS からの通常 IRC 想定
    "displ": "energy",
    "displ_energy": 1.0e-3,
    "displ_length": 0.10,
    "rms_grad_thresh": 1.0e-3,
    "hard_rms_grad_thresh": None,
    "energy_thresh": 1.0e-6,
    "imag_below": 0.0,
    "force_inflection": True,
    "check_bonds": False,
    "out_dir": "./result_irc/",
    "prefix": "",
    "dump_fn": "irc_data.h5",
    "dump_every": 5,

    # EulerPC 固有
    "hessian_update": "bofill",
    "hessian_recalc": None,
    "max_pred_steps": 500,
    "loose_cycles": 3,
    "corr_func": "mbs",
}


# --------------------------
# ヘルパ
# --------------------------

def _deep_update(dst: Dict[str, Any], src: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """辞書を再帰的にマージ（src が優先）。"""
    if not src:
        return dst
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML のトップレベルは mapping である必要があります (got: {type(data)})")
    return data


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if fa else ""
    return g


def _echo_convert_trj_to_pdb_if_exists(trj_path: Path, ref_pdb: Path, out_path: Path) -> None:
    if trj_path.exists():
        try:
            convert_xyz_to_pdb(trj_path, ref_pdb, out_path)
            click.echo(f"[convert] Wrote '{out_path}'.")
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert '{trj_path.name}' to PDB: {e}", err=True)


# --------------------------
# CLI
# --------------------------

@click.command(
    help="EulerPC による IRC 計算。CLI では例にある引数のみ受け付け、それ以外は YAML から読み込みます。",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="入力構造ファイル (.pdb, .xyz, .trj, など)。",
)
@click.option("-q", "--charge", type=int, required=True, help="全電荷。calc セクションをこの値で上書き。")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="多重度 (2S+1)。calc セクションをこの値で上書き。")
@click.option("--max-cycles", type=int, default=None, help="IRC の最大ステップ数（YAML の irc.max_cycles を上書き）。")
@click.option("--step-size", type=float, default=None, help="IRC の非質量重み付け座標でのステップ長（YAML の irc.step_length を上書き）。")
@click.option("--root", type=int, default=None, help="初期変位に用いる虚振動のモード番号（YAML の irc.root を上書き）。")
@click.option("--forward", type=bool, default=None, help="正方向 IRC を実行（YAML の irc.forward を上書き）。True/False を明示。")
@click.option("--backward", type=bool, default=None, help="逆方向 IRC を実行（YAML の irc.backward を上書き）。True/False を明示。")
@click.option("--out-dir", type=str, default="./result_irc/", show_default=True, help="出力ディレクトリ。YAML の irc.out_dir を上書き。")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="追加パラメータを含む YAML（セクション: geom, calc, irc）。",
)
def cli(
    input_path: Path,
    charge: int,
    spin: int,
    max_cycles: Optional[int],
    step_size: Optional[float],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) 設定の構築（YAML → デフォルト, その後 CLI で上書き）
        # --------------------------
        yaml_cfg = _load_yaml(args_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg:  Dict[str, Any] = dict(IRC_KW_DEFAULT)

        _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
        _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
        _deep_update(irc_cfg,  yaml_cfg.get("irc",  {}))

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)

        if max_cycles is not None:
            irc_cfg["max_cycles"] = int(max_cycles)
        if step_size is not None:
            irc_cfg["step_length"] = float(step_size)
        if root is not None:
            irc_cfg["root"] = int(root)
        if forward is not None:
            irc_cfg["forward"] = bool(forward)
        if backward is not None:
            irc_cfg["backward"] = bool(backward)
        if out_dir:
            irc_cfg["out_dir"] = str(out_dir)

        out_dir_path = Path(irc_cfg["out_dir"]).resolve()
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # 画面表示用（freeze_atoms を見やすく）
        click.echo(_pretty_block("geom", _format_geom_for_echo(geom_cfg)))
        click.echo(_pretty_block("calc", calc_cfg))
        click.echo(_pretty_block("irc",  {**irc_cfg, "out_dir": str(out_dir_path)}))

        # --------------------------
        # 2) 幾何のロード & UMA 設定
        # --------------------------
        coord_type = geom_cfg.get("coord_type", "cart")
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)

        geometry = geom_loader(input_path, coord_type=coord_type, **coord_kwargs)

        calc_builder_or_instance = uma_pysis(**calc_cfg)
        try:
            # 工場関数が返る場合
            geometry.set_calculator(calc_builder_or_instance())
        except TypeError:
            # すでにインスタンスが返ってきた場合
            geometry.set_calculator(calc_builder_or_instance)

        # --------------------------
        # 3) EulerPC を構築して実行
        # --------------------------
        # EulerPC.__init__ は IRC.__init__ に **kwargs をそのまま渡せる
        eulerpc = EulerPC(geometry, **irc_cfg)

        click.echo("\n=== IRC (EulerPC) started ===\n")
        eulerpc.run()
        click.echo("\n=== IRC (EulerPC) finished ===\n")

        # --------------------------
        # 4) PDB 変換（入力が PDB の場合）
        # --------------------------
        if input_path.suffix.lower() == ".pdb":
            ref_pdb = input_path.resolve()

            # IRC 全体
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.trj'}",
                ref_pdb,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.pdb'}",
            )
            # forward / backward
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.trj'}",
                ref_pdb,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.pdb'}",
            )
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.trj'}",
                ref_pdb,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.pdb'}",
            )

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for IRC: {hh:02d}:{mm:02d}:{ss:06.3f}")

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = textwrap.indent("".join(__import__("traceback").format_exception(type(e), e, e.__traceback__)), "  ")
        click.echo("Unhandled exception during IRC:\n" + tb, err=True)
        sys.exit(1)


# 直接実行
if __name__ == "__main__":
    cli()
