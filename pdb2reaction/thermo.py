# pdb2reaction/thermo.py
from pathlib import Path
import click


@click.command(help="Thermochemistry Analysis")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structure (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"thermo {input_path.resolve()}")
