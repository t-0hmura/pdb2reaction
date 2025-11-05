# pdb2reaction/path_opt.py
"""
MEP optimization by Growing String.

- Endpoints: two structures (e.g. reactant/product) given by -i/--input.
- Calculator: UMA (via uma_pysis wrapper).
- Path optimizer: GrowingString + StringOptimizer (pysisyphus).
- YAML config with sub-sections: geom, calc, gs, opt  （CLI > YAML > defaults）
- If inputs are PDB and --freeze-links=True, freeze parent atoms of link hydrogens.
- Dump the highest-energy image (HEI) as gsm_hei.xyz and gsm_hei.pdb (if a PDB reference is available).
Examples
--------
pdb2reaction path_opt -i reac.pdb prod.pdb -q 0 -s 1 \
  --freeze-links True --max-nodes 10 --max-cycles 100 \
  --dump True --dump_dir ./result_path_opt/ --args-yaml ./args.yaml
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import sys
import traceback
import textwrap

import click
import numpy as np
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError

from .uma_pysis import uma_pysis
from .utils import convert_xyz_to_pdb, freeze_links


# -----------------------------------------------
# Defaults (YAML/CLI で上書き)
# -----------------------------------------------

# 幾何（入力の扱い）
GEOM_KW: Dict[str, Any] = {
    "coord_type": "cart",   # GrowingString は cart 推奨
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
    "max_nodes": 10,            # 端点含めて max_nodes+2 画像
    "perp_thresh": 5e-3,        # フロンティア成長判定（⊥力のRMS/NORM）
    "reparam_check": "rms",     # "rms" | "norm"
    "reparam_every": 2,
    "reparam_every_full": 3,
    "param": "equi",
    "max_micro_cycles": 5,
    "reset_dlc": True,
    "climb": True,              # climbing image 有効
    "climb_rms": 5e-3,          # CI 開始の閾値（RMS(force)）
    "climb_lanczos": True,      # HEI 接線を Lanczos で
    "climb_lanczos_rms": 5e-3,
    "scheduler": None,         # 直列計算（同一 calculator 使い回し前提）
}

# StringOptimizer（最適化制御）
OPT_KW: Dict[str, Any] = {
    "type": "string",           # 記録用タグ
    "stop_in_when_full": 20,    # fully grown 後に N サイクルで停止
    "align": False,             # Kabsch アライン（cart では通常 True 推奨だが既定 False）
    "scale_step": "global",
    "max_cycles": 10000,
    "dump": False,
    "dump_restart": False,
    "reparam_thresh": 0.0,      # 連続 reparam 後のRMS座標差チェック無効化
    "coord_diff_thresh": 0.0,   # 近接画像チェック（0で無効）
    "out_dir": "./result_path_opt/",
    "print_every": 1,
}


# -----------------------------------------------
# Utilities
# -----------------------------------------------

def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively update dict *dst* with *src*, returning *dst*."""
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")
    return data


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
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
    """2 つの端点構造を読み込み、必要に応じて freeze_atoms を設定。"""
    geoms = []
    for p in paths:
        g = geom_loader(p, coord_type=coord_type)
        # freeze_atoms は Geometry に後付けで設定
        freeze = list(base_freeze)
        if auto_freeze_links and p.suffix.lower() == ".pdb":
            detected = _freeze_links_for_pdb(p)
            if detected:
                freeze = sorted(set(freeze).union(detected))
                click.echo(f"[freeze-links] {p.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="MEP optimization by Growing String.",
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
@click.option("--max-nodes", type=int, default=10, show_default=True, help="Internal nodes (string has max_nodes+2 images incl. endpoints)")
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Max optimization cycles")
@click.option("--dump", type=click.BOOL, default=False, show_default=True, help="Dump optimizer trajectory/restarts during run")
@click.option("--out-dir", "out_dir", type=str, default="./result_path_opt/", show_default=True, help="Output directory")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt).",
)
def cli(
    input_paths: Sequence[Path],
    charge: int,
    spin: int,
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    dump: bool,
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
    try:
        # --------------------------
        # 1) 設定の組み立て（defaults ← YAML ← CLI）
        # --------------------------
        yaml_cfg = _load_yaml(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg   = dict(GS_KW)
        opt_cfg  = dict(OPT_KW)

        _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
        _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
        _deep_update(gs_cfg,   yaml_cfg.get("gs",   {}))
        _deep_update(opt_cfg,  yaml_cfg.get("opt",  {}))

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)

        gs_cfg["max_nodes"] = int(max_nodes)

        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir  # CLI では --out-dir を受け取り Optimizer には out_dir を渡す

        # 表示用：解決後の設定
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = _format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_opt  = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        click.echo(_pretty_block("geom", echo_geom))
        click.echo(_pretty_block("calc", echo_calc))
        click.echo(_pretty_block("gs",   echo_gs))
        click.echo(_pretty_block("opt",  echo_opt))

        # --------------------------
        # 2) 構造の準備（端点2つの読み込みと freeze）
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # 入力パス
        p0, p1 = Path(input_paths[0]), Path(input_paths[1])

        # 端点の読み込み（PDB の場合は link 親凍結を合流）
        geoms = _load_two_endpoints(
            paths=[p0, p1],
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )

        # 共通 UMA 計算器（同一インスタンスを全画像で共有）
        shared_calc = uma_pysis(**calc_cfg)
        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # GrowingString が新規ノードを生成する際に使う
            return shared_calc

        # --------------------------
        # 3) 経路オブジェクトとオプティマイザの構築
        # --------------------------
        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        # StringOptimizer は out_dir を "out_dir" で受け取る
        opt_args = dict(opt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"}  # type は単なるタグ
        )

        # --------------------------
        # 4) 最適化実行
        # --------------------------
        click.echo("\n=== Growing String optimization started ===\n")
        optimizer.run()
        click.echo("\n=== Growing String optimization finished ===\n")

        # --------------------------
        # 5) 最終経路の書き出し（final_geometries.trj）
        # --------------------------
        final_trj = out_dir_path / "final_geometries.trj"
        try:
            with open(final_trj, "w") as f:
                f.write(gs.as_xyz())
            click.echo(f"[write] Wrote '{final_trj}'.")

            if input_paths[0].suffix.lower() == ".pdb":
                ref_pdb = input_paths[0].resolve()

                try:
                    out_pdb = out_dir_path / "final_geometries.pdb"
                    convert_xyz_to_pdb(final_trj, ref_pdb, out_pdb)
                    click.echo(f"[convert] Wrote '{out_pdb}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert MEP path trajectory to PDB: {e}", err=True)

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        try:
            energies = np.array(gs.energy, dtype=float)
            hei_idx = int(np.argmax(energies))
            hei_geom = gs.images[hei_idx]
            hei_E = float(energies[hei_idx])

            hei_xyz = out_dir_path / "gsm_hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"HEI idx={hei_idx}; E={hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            ref_pdb = None
            if input_paths[0].suffix.lower() == ".pdb":
                ref_pdb = input_paths[0].resolve()
            if ref_pdb is not None:
                hei_pdb = out_dir_path / "gsm_hei.pdb"
                convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
                click.echo(f"[convert] Wrote '{hei_pdb}'.")
            else:
                click.echo("[convert] Skipped 'gsm_hei.pdb' (no PDB reference among inputs).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)


    except OptimizationError as e:
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during path optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


# PDB link 親検出の補助（opt.py と同様に名前衝突回避用）
def freeze_links_helper(pdb_path: Path):
    return freeze_links(pdb_path)


# 直接実行対応
if __name__ == "__main__":
    cli()
