import click
from click.testing import CliRunner

from pdb2reaction.cli.app import cli as root_cli


def _option(command: str, name: str) -> click.Option:
    cmd = root_cli.get_command(click.Context(root_cli), command)
    assert cmd is not None
    return next(param for param in cmd.params if param.name == name)


def test_all_exposes_canonical_names_with_compatibility_aliases() -> None:
    assert _option("all", "workers").opts == ["--uma-workers", "--workers"]
    assert _option("all", "dft_func_basis").opts == ["--func-basis", "--dft-func-basis"]
    assert _option("all", "max_cycles_dmf").opts == ["--dmf-max-iterations", "--max-cycles-dmf"]
    assert _option("all", "thresh_dmf").opts == ["--dmf-tol", "--thresh-dmf"]
    assert _option("all", "resume_segment").opts == ["--resume-segment"]


def test_all_public_booleans_are_toggle_options() -> None:
    for name in ("include_h2o", "do_tsopt", "do_thermo", "do_dft", "scan_preopt_override"):
        option = _option("all", name)
        assert option.is_bool_flag
        assert option.secondary_opts


def test_dft_canonical_names_keep_old_aliases() -> None:
    assert _option("dft", "max_cycle").opts == ["--scf-max-cycles", "--max-cycle"]
    assert _option("dft", "conv_tol").opts == ["--scf-tol", "--conv-tol"]
    assert _option("dft", "engine").opts == ["--dft-engine", "--engine"]
    assert _option("dft", "memory").opts == ["--dft-memory", "--dft-mem"]


def test_conflicting_alias_values_are_rejected_before_execution() -> None:
    result = CliRunner().invoke(
        root_cli,
        ["dft", "--scf-tol", "1e-8", "--conv-tol", "1e-7"],
    )
    assert result.exit_code == 2
    assert "Conflicting values were supplied through aliases" in result.output


def test_every_multiplicity_option_requires_a_positive_integer() -> None:
    context = click.Context(root_cli)
    checked = []
    for command_name in root_cli.list_commands(context):
        command = root_cli.get_command(context, command_name)
        for parameter in getattr(command, "params", ()):
            if isinstance(parameter, click.Option) and "--multiplicity" in parameter.opts:
                checked.append(command_name)
                assert isinstance(parameter.type, click.IntRange), command_name
                assert parameter.type.min == 1, command_name
    assert checked


def test_shared_choice_options_use_one_case_rule() -> None:
    expected = {"--backend": True, "--solvent-model": True, "--dft-solvent-model": True, "--dft-engine": False}
    context = click.Context(root_cli)
    checked = set()
    for command_name in root_cli.list_commands(context):
        command = root_cli.get_command(context, command_name)
        for parameter in getattr(command, "params", ()):
            if not isinstance(parameter, click.Option) or not isinstance(parameter.type, click.Choice):
                continue
            for flag, case_sensitive in expected.items():
                if flag in parameter.opts:
                    checked.add(flag)
                    assert parameter.type.case_sensitive is case_sensitive, (command_name, flag)
    assert checked == set(expected)
