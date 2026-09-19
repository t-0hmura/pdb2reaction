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
