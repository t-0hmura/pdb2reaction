# pdb2reaction/opt.py
from __future__ import annotations

from pathlib import Path
import click


@click.command(help="単一構造の最適化（プレースホルダ）")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="入力構造ファイル（PDB/XYZ 等）"
)
def cli(input_path: Path) -> None:
    print(f"opt {input_path.resolve()}")
