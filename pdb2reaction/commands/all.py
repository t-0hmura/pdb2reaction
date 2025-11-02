# pdb2reaction/commands/all.py
from pathlib import Path
import click

@click.command(help="パイプラインの全ステージを実行する。")
@click.option("-i", "--input", "input_pdb", type=click.Path(path_type=Path), required=True,
              help="入力PDBファイル")
@click.option("--keep-intermediates/--no-keep-intermediates", default=False, show_default=True,
              help="中間生成物を残すか")
def cli(input_pdb: Path, keep_intermediates: bool) -> None:
    click.echo(f"[all] input={input_pdb} keep={keep_intermediates}")
    # TODO: 実処理 run_all(input_pdb, keep_intermediates)
