# pdb2reaction/cli.py
import click

class DefaultGroup(click.Group):
    def __init__(self, *args, default: str | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self._default_cmd = default

    def parse_args(self, ctx, args):
        if any(a in ("-h", "--help") for a in args):
            return super().parse_args(ctx, args)

        if self._default_cmd is not None:
            if not args or args[0].startswith("-"):
                args.insert(0, self._default_cmd)
        return super().parse_args(ctx, args)


from .all import cli as all_cmd
from .scan import cli as scan_cmd
from .opt import cli as opt_cmd
from .path_opt import cli as path_opt_cmd
from .ts_opt import cli as ts_opt_cmd
from .irc import cli as irc_cmd
from .thermo import cli as thermo_cmd
from .extract import cli as extract_cmd
from .trj2fig import cli as trj2fig_cmd


@click.group(
    cls=DefaultGroup,
    default="all",
    help="pdb2reaction: Execute each step by subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
def cli() -> None:
    pass


cli.add_command(all_cmd, name="all")
cli.add_command(scan_cmd, name="scan")
cli.add_command(opt_cmd, name="opt")
cli.add_command(path_opt_cmd, name="path_opt")
cli.add_command(ts_opt_cmd, name="ts_opt")
cli.add_command(irc_cmd, name="irc")
cli.add_command(thermo_cmd, name="thermo")
cli.add_command(extract_cmd, name="extract")
cli.add_command(trj2fig_cmd, name="trj2fig")


if __name__ == "__main__":
    cli()
