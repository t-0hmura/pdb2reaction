# pdb2reaction/cli.py
from __future__ import annotations
from typing import List
import click

# Subcommands
from .commands import all as all_cmd
from .commands import opt as opt_cmd

CONTEXT_SETTINGS = dict(
    help_option_names=["-h", "--help"],
    ignore_unknown_options=True,   # accept unknown options at the group level
    allow_extra_args=True,         # pass remaining args to the subcommand
)

def _normalize_args(args: List[str]) -> List[str]:
    """
    Make the CLI forgiving:
    - Convert single-dash long options to double-dash (e.g., '-charge' -> '--charge')
    - Convert underscore style to dash style ('--opt_mode' -> '--opt-mode')
    - Leave short options like '-i' intact
    - Do not touch values (tokens that don't start with '-')
    """
    out: List[str] = []
    for tok in args:
        if tok.startswith('-') and not tok.startswith('--') and len(tok) > 2:
            # Treat '-long' as '--long'
            tok = '--' + tok.lstrip('-')
        if tok.startswith('--'):
            if '=' in tok:
                name, val = tok.split('=', 1)
                name = name.replace('_', '-')
                tok = f"{name}={val}"
            else:
                tok = tok.replace('_', '-')
        out.append(tok)
    return out

@click.group(context_settings=CONTEXT_SETTINGS, invoke_without_command=True)
@click.option(
    "-m", "--mode",
    type=click.Choice(["all", "opt"], case_sensitive=False),
    help="Explicit subcommand. Defaults to 'all' when omitted.",
)
@click.pass_context
def cli(ctx: click.Context, mode: str | None) -> None:
    """pdb2reaction: a CLI toolbox for PDB-driven workflows.

    Examples:
      - pdb2reaction opt -i a.pdb ...
      - pdb2reaction -m opt -i a.pdb ...
      - pdb2reaction -i a.pdb            # no subcommand -> defaults to 'all'
    """
    # If a subcommand was explicitly invoked, do nothing here
    if ctx.invoked_subcommand:
        return

    # If user asked for help/version at the group level, show it and exit
    if any(x in ctx.args for x in ("-h", "--help", "-V", "--version")) and not mode:
        return

    # Decide subcommand: explicit --mode or default 'all'
    subcmd = (mode or "all").lower()
    if subcmd not in ("all", "opt"):
        raise click.UsageError(f"Unknown mode: {subcmd}")

    # Normalize remaining args and dispatch to the chosen subcommand
    norm_args = _normalize_args(list(ctx.args))
    # Forward to the subcommand, preserving exit status
    ctx.exit(cli.main(args=[subcmd] + norm_args, standalone_mode=False))

# Register subcommands
cli.add_command(all_cmd.cli, name="all")
cli.add_command(opt_cmd.cli, name="opt")
