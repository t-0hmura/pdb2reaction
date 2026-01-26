# pdb2reaction/cli.py

import importlib
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


class LazyGroup(DefaultGroup):
    def __init__(self, *args, lazy_commands: dict[str, tuple[str, str]] | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self._lazy_commands = lazy_commands or {}

    def get_command(self, ctx, cmd_name):
        cmd = super().get_command(ctx, cmd_name)
        if cmd is not None:
            return cmd
        target = self._lazy_commands.get(cmd_name)
        if target is None:
            return None
        module_name, attr = target
        module = importlib.import_module(module_name)
        cmd = getattr(module, attr)
        self.commands[cmd_name] = cmd
        return cmd

    def list_commands(self, ctx):
        return sorted(set(self.commands) | set(self._lazy_commands))


LAZY_COMMANDS = {
    "all": ("pdb2reaction.all", "cli"),
    "scan": ("pdb2reaction.scan", "cli"),
    "opt": ("pdb2reaction.opt", "cli"),
    "path-opt": ("pdb2reaction.path_opt", "cli"),
    "path-search": ("pdb2reaction.path_search", "cli"),
    "tsopt": ("pdb2reaction.tsopt", "cli"),
    "freq": ("pdb2reaction.freq", "cli"),
    "irc": ("pdb2reaction.irc", "cli"),
    "trj2fig": ("pdb2reaction.trj2fig", "cli"),
    "add-elem-info": ("pdb2reaction.add_elem_info", "cli"),
    "dft": ("pdb2reaction.dft", "cli"),
    "scan2d": ("pdb2reaction.scan2d", "cli"),
    "scan3d": ("pdb2reaction.scan3d", "cli"),
}


@click.group(
    cls=LazyGroup,
    default="all",
    help="pdb2reaction: Root command to execute each subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
    lazy_commands=LAZY_COMMANDS,
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
