# pdb2reaction/cli.py
from __future__ import annotations
import click
from .commands import all as all_cmd
from .commands import opt as opt_cmd

CONTEXT_SETTINGS = dict(
    help_option_names=["-h", "--help"],
    ignore_unknown_options=True,  # 親で未定義のオプションを許容
    allow_extra_args=True,        # 残りの引数を保持（後でサブコマンドへ委譲）
)

@click.group(context_settings=CONTEXT_SETTINGS, invoke_without_command=True)
@click.option(
    "-m", "--mode",
    type=click.Choice(["all", "opt"], case_sensitive=False),
    help="処理モード。サブコマンド未指定時はここで指定、未指定なら 'all' を既定にします。",
)
@click.pass_context
def cli(ctx: click.Context, mode: str | None) -> None:
    """pdb2reaction: PDB から反応解析パイプラインを起動するCLI.

    使い方:
      pdb2reaction opt -i a.pdb ...
      pdb2reaction -m opt -i a.pdb ...
      pdb2reaction -i a.pdb ...         # サブコマンド未指定 → 既定の 'all' を実行
    """
    # サブコマンドが明示された場合は親では何もしない
    if ctx.invoked_subcommand:
        if mode:
            raise click.UsageError(
                "サブコマンドと -m/--mode は同時指定できません。どちらか一方にしてください。"
            )
        return

    # サブコマンド未指定 → mode か 'all' にフォールバック
    subcmd = (mode or "all").lower()
    if subcmd not in ("all", "opt"):
        raise click.UsageError(f"未知のモード: {subcmd}")

    # 親で受けた未解釈の引数（-i 等）を、そのまま該当サブコマンドへ委譲
    args = [subcmd] + list(ctx.args)
    ctx.exit(cli.main(args=args, standalone_mode=False))

# サブコマンドを登録
cli.add_command(all_cmd.cli, name="all")
cli.add_command(opt_cmd.cli, name="opt")
