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
from .path_search import cli as path_search_cmd
from .tsopt import cli as tsopt_cmd
from .freq import cli as freq_cmd
from .irc import cli as irc_cmd
from .trj2fig import cli as trj2fig_cmd
from .add_elem_info import cli as add_elem_info_cmd
from .dft import cli as dft_cmd
from .scan2d import cli as scan2d_cmd
from .scan3d import cli as scan3d_cmd


@click.group(
    cls=DefaultGroup,
    default="all",
    help="pdb2reaction: Execute each step by subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
def cli() -> None:
    pass


@click.command(
    name="extract",
    help="Extract a binding pocket.",
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "help_option_names": [],
    },
)
@click.pass_context
def extract_cmd(ctx: click.Context) -> None:
    import sys
    import os
    from . import extract as _extract_mod

    argv_backup = sys.argv[:]
    try:
        prog_base = (ctx.find_root().info_name or os.path.basename(sys.argv[0]))
        sys.argv = [f"{prog_base} extract"] + list(ctx.args)
        _extract_mod.extract()
    finally:
        sys.argv = argv_backup


cli.add_command(all_cmd, name="all")
cli.add_command(scan_cmd, name="scan")
cli.add_command(opt_cmd, name="opt")
cli.add_command(path_opt_cmd, name="path-opt")
cli.add_command(path_search_cmd, name="path-search")
cli.add_command(tsopt_cmd, name="tsopt")
cli.add_command(freq_cmd, name="freq")
cli.add_command(irc_cmd, name="irc")
cli.add_command(extract_cmd, name="extract")
cli.add_command(trj2fig_cmd, name="trj2fig")
cli.add_command(add_elem_info_cmd, name="add-elem-info")
cli.add_command(dft_cmd, name="dft")
cli.add_command(scan2d_cmd, name="scan2d")
cli.add_command(scan3d_cmd, name="scan3d")

# Disable pysisyphus logging
import logging
logging.disable(logging.CRITICAL)