# pdb2reaction/commands/opt.py
from pathlib import Path
import click

@click.command(help="(例) 構造最適化ステップのみを実行する。")
@click.option("-i", "--input", "input_pdb", type=click.Path(path_type=Path), required=True,
              help="入力PDBファイル")
@click.option("--max-steps", type=int, default=200, show_default=True, help="最大ステップ数")
def cli(input_pdb: Path, max_steps: int) -> None:
    click.echo(f"[opt] input={input_pdb} max_steps={max_steps}")
    # TODO: 実処理 run_opt(input_pdb, max_steps)
