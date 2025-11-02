# pdb2reaction/cli.py
from __future__ import annotations

import click

class DefaultGroup(click.Group):
    """引数が空、または先頭がオプションのとき既定サブコマンドを挿入する Group"""
    def __init__(self, *args, default: str | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self._default_cmd = default

    def parse_args(self, ctx, args):
        if any(a in ("-h", "--help") for a in args):
            return super().parse_args(ctx, args)

        # 引数が空 もしくは 先頭引数がオプション（- 始まり）の場合は既定サブコマンドを先頭に挿入
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


@click.group(
    cls=DefaultGroup,
    default="all",
    help="pdb2reaction: サブコマンドで各処理を起動します。例: pdb2reaction opt -i a.pdb",
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


if __name__ == "__main__":
    cli()
