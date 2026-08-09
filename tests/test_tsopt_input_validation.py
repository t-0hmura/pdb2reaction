"""TS-optimization configuration validation."""

from __future__ import annotations

from copy import deepcopy
from inspect import signature

import click
import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli


def test_subspace_line_search_kwargs_exist_only_for_rsprfo(tmp_path) -> None:
    from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer
    from pdb2reaction.workflows.tsopt import _build_rsirfo_kwargs

    assert signature(TSHessianOptimizer).parameters["min_line_search"].default is False
    assert signature(TSHessianOptimizer).parameters["max_line_search"].default is False

    defaults = _build_rsirfo_kwargs({}, {}, tmp_path, kind="rsprfo")
    assert defaults["min_line_search"] is False
    assert defaults["max_line_search"] is False

    configured = {"min_line_search": True, "max_line_search": True}

    rsprfo = _build_rsirfo_kwargs({}, configured, tmp_path, kind="rsprfo")
    assert rsprfo["min_line_search"] is True
    assert rsprfo["max_line_search"] is True

    for kind in ("rsirfo", "trim"):
        kwargs = _build_rsirfo_kwargs({}, configured, tmp_path, kind=kind)
        assert "min_line_search" not in kwargs
        assert "max_line_search" not in kwargs


@pytest.mark.parametrize("kind", ("rsprfo", "rsirfo", "trim"))
def test_ts_rfo_omits_rfo_overlaps_and_requires_one_root(kind, tmp_path) -> None:
    from pdb2reaction.workflows.tsopt import _build_rsirfo_kwargs

    kwargs = _build_rsirfo_kwargs(
        {"rfo_overlaps": True},
        {"rfo_overlaps": True, "roots": [0]},
        tmp_path,
        kind=kind,
    )
    assert "rfo_overlaps" not in kwargs
    assert kwargs["roots"] == [0]

    for roots in ([], [0, 1]):
        with pytest.raises(click.BadParameter, match="exactly one root index"):
            _build_rsirfo_kwargs({}, {"roots": roots}, tmp_path, kind=kind)


def test_dimer_has_independent_yaml_line_search_setting() -> None:
    from pdb2reaction.core.defaults import HESSIAN_DIMER_CLI_KW
    from pdb2reaction.core.utils import apply_yaml_overrides

    configured = deepcopy(HESSIAN_DIMER_CLI_KW)
    assert configured["lbfgs"]["line_search"] is True

    apply_yaml_overrides(
        {"hessian_dimer": {"lbfgs": {"line_search": False}}},
        [(configured, (("hessian_dimer",),))],
    )
    assert configured["lbfgs"]["line_search"] is False


def test_rsprfo_yaml_line_search_values_reach_constructor(
    monkeypatch, tmp_path,
) -> None:
    """An explicit YAML choice must not be replaced by workflow policy."""
    from pdb2reaction.workflows import tsopt

    class ConstructorReached(RuntimeError):
        pass

    captured = {}

    class CaptureOptimizer:
        def __init__(self, geometry, **kwargs):
            captured.update(kwargs)
            raise ConstructorReached

    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", CaptureOptimizer)
    monkeypatch.setattr(tsopt, "create_calculator", lambda **kwargs: object())

    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "rsirfo:\n"
        "  min_line_search: true\n"
        "  max_line_search: true\n"
        "  rfo_overlaps: true\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt", "-i", str(source), "-q", "0", "--opt-mode", "rsprfo",
            "-o", str(tmp_path / "out"), "--config", str(config),
        ],
    )

    assert result.exit_code == 1
    assert "ConstructorReached" in result.output
    assert captured["min_line_search"] is True
    assert captured["max_line_search"] is True
    assert "rfo_overlaps" not in captured


def test_dimer_hessian_update_interval_must_advance(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "hessian_dimer:\n  update_interval_hessian: 0\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt",
            "-i",
            str(source),
            "-q",
            "0",
            "--opt-mode",
            "dimer",
            "--dry-run",
            "--config",
            str(config),
        ],
    )

    assert result.exit_code != 0
    assert "update_interval_hessian must be a positive integer" in result.output


def test_tsopt_dry_run_rejects_invalid_configured_hessian_mode(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("calc:\n  hessian_calc_mode: typo\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt", "-i", str(source), "-q", "0", "--dry-run",
            "--config", str(config),
        ],
    )

    assert result.exit_code == 1
    assert "Unsupported hessian_calc_mode 'typo'" in result.output


def test_tsopt_cli_print_every_overrides_nested_dimer_config(
    monkeypatch, tmp_path,
) -> None:
    from pdb2reaction.workflows import tsopt

    stripped = []
    original_strip = tsopt.strip_inherited_keys

    def capture_strip(values, *args, **kwargs):
        stripped.append(dict(values))
        return original_strip(values, *args, **kwargs)

    monkeypatch.setattr(tsopt, "strip_inherited_keys", capture_strip)
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "hessian_dimer:\n  lbfgs:\n    print_every: 99\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt", "-i", str(source), "-q", "0", "--opt-mode", "dimer",
            "--print-every", "7", "--show-config", "--dry-run",
            "--config", str(config),
        ],
    )

    assert result.exit_code == 0, result.output
    assert any(values.get("print_every") == 7 for values in stripped)


def test_tsopt_rejects_input_aliasing_final_geometry(tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    source = out / "final_geometry.xyz"
    original = "1\natom\nHe 0 0 0\n"
    source.write_text(original, encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        ["tsopt", "-i", str(source), "-q", "0", "-o", str(out)],
    )

    assert result.exit_code == 2
    assert "collides with a reserved TSOPT output" in result.output
    assert source.read_text(encoding="utf-8") == original
