# pdb2reaction/cli.py

import click


class DefaultGroup(click.Group):
    def __init__(self, *args, default: str | None = None, command_loaders=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._default_cmd = default
        self._command_loaders = command_loaders or {}

    def parse_args(self, ctx, args):
        if any(a in ("-h", "--help") for a in args):
            return super().parse_args(ctx, args)

        if self._default_cmd is not None:
            if not args or args[0].startswith("-"):
                args.insert(0, self._default_cmd)
        return super().parse_args(ctx, args)

    def get_command(self, ctx, cmd_name):
        if cmd_name in self._command_loaders:
            cmd = self._command_loaders[cmd_name]()
            self.add_command(cmd, name=cmd_name)
            return cmd
        return super().get_command(ctx, cmd_name)


def _load_all():
    from .all import cli as cmd
    return cmd


def _load_scan():
    from .scan import cli as cmd
    return cmd


def _load_opt():
    from .opt import cli as cmd
    return cmd


def _load_path_opt():
    from .path_opt import cli as cmd
    return cmd


def _load_path_search():
    from .path_search import cli as cmd
    return cmd


def _load_tsopt():
    from .tsopt import cli as cmd
    return cmd


def _load_freq():
    from .freq import cli as cmd
    return cmd


def _load_irc():
    from .irc import cli as cmd
    return cmd


def _load_trj2fig():
    from .trj2fig import cli as cmd
    return cmd


def _load_add_elem_info():
    from .add_elem_info import cli as cmd
    return cmd


def _load_dft():
    from .dft import cli as cmd
    return cmd


def _load_scan2d():
    from .scan2d import cli as cmd
    return cmd


def _load_scan3d():
    from .scan3d import cli as cmd
    return cmd


COMMAND_LOADERS = {
    "all": _load_all,
    "scan": _load_scan,
    "opt": _load_opt,
    "path-opt": _load_path_opt,
    "path-search": _load_path_search,
    "tsopt": _load_tsopt,
    "freq": _load_freq,
    "irc": _load_irc,
    "trj2fig": _load_trj2fig,
    "add-elem-info": _load_add_elem_info,
    "dft": _load_dft,
    "scan2d": _load_scan2d,
    "scan3d": _load_scan3d,
}


@click.group(
    cls=DefaultGroup,
    default="all",
    help="pdb2reaction: Root command to execute each subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
    command_loaders=COMMAND_LOADERS,
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


cli.add_command(extract_cmd, name="extract")

# Disable pysisyphus logging
import logging
logging.disable(logging.CRITICAL)

# Filter noisy UMA/torch_dmf warnings that clutter CLI output
import warnings

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"var\(\): degrees of freedom is <= 0\. Correction should be strictly less than the reduction factor.*",
    module=r"fairchem\.core\.models\.uma\.escn_moe"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"Sparse CSR tensor support is in beta state.*",
    module=r"torch_dmf"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"t_eval update skipped due to insufficient candidates",
    module=r"torch_dmf"
)
