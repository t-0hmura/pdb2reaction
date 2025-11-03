# pdb2reaction/path_opt.py
from pathlib import Path
import click


@click.command(help="MEP optimization by Growing String")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structures (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"path_opt {input_path.resolve()}")
