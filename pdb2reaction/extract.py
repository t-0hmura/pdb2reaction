# pdb2reaction/extract.py
from pathlib import Path
import click


@click.command(help="Extract Active Site structures")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structures (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"extract {input_path.resolve()}")
