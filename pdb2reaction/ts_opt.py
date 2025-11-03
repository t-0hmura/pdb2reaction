# pdb2reaction/ts_opt.py
from pathlib import Path
import click


@click.command(help="TS optimization")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structure (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"ts_opt {input_path.resolve()}")
