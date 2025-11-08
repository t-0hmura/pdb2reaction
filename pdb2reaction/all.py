# pdb2reaction/all.py

from pathlib import Path
import click


@click.command(help="Execute whole pipeline")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structures (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"all {input_path.resolve()}")
