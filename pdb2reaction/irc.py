# pdb2reaction/irc.py
from pathlib import Path
import click


@click.command(help="IRC calculation")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structure (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"irc {input_path.resolve()}")
