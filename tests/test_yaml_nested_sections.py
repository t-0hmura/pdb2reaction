"""Nested YAML spellings: opt.lbfgs = lbfgs, opt.rfo = rfo, freq.thermo = thermo."""

from __future__ import annotations

import json
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.utils import apply_yaml_overrides, build_scan_configs
from pdb2reaction.workflows import freq as freq_workflow
from pdb2reaction.workflows import opt as opt_workflow
from pdb2reaction.workflows import tsopt

H2_XYZ = "2\nH2\nH 0 0 0\nH 0 0 0.74\n"


def test_both_spellings_are_combined_and_conflicts_rejected():
    lbfgs_cfg = {"max_step": 0.3, "memory": 5}
    yaml_cfg = {"lbfgs": {"memory": 7}, "opt": {"lbfgs": {"max_step": 0.2, "memory": 7}}}
    apply_yaml_overrides(yaml_cfg, [(lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs")))])
    assert lbfgs_cfg == {"max_step": 0.2, "memory": 7}

    yaml_cfg["lbfgs"]["memory"] = 9
    with pytest.raises(click.BadParameter, match="lbfgs.memory and opt.lbfgs.memory conflict"):
        apply_yaml_overrides(yaml_cfg, [({}, (("lbfgs",), ("opt", "lbfgs")))])


def test_nested_sections_never_reach_the_parent_target():
    opt_cfg, freq_cfg, thermo_cfg = {"max_cycles": 1}, {"max_write": 1}, {"temperature": 298.15}
    yaml_cfg = {
        "opt": {"max_cycles": 5, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}},
        "freq": {"max_write": 3, "thermo": {"temperature": 310.0}},
    }
    # No lbfgs/rfo target is registered, as in tsopt.
    apply_yaml_overrides(
        yaml_cfg,
        [
            (opt_cfg, (("opt",),)),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",), ("freq", "thermo"))),
        ],
    )
    assert opt_cfg == {"max_cycles": 5}
    assert freq_cfg == {"max_write": 3}
    assert thermo_cfg == {"temperature": 310.0}


def test_dimer_line_search_section_is_not_a_nested_spelling():
    simple_cfg, lbfgs_cfg = {}, {"max_step": 0.3}
    apply_yaml_overrides(
        {"hessian_dimer": {"lbfgs": {"max_step": 0.1}}},
        [(simple_cfg, (("hessian_dimer",),)), (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs")))],
    )
    assert simple_cfg == {"lbfgs": {"max_step": 0.1}}
    assert lbfgs_cfg == {"max_step": 0.3}


def _opt_dry_run(tmp_path: Path, monkeypatch, config: dict, *extra: str):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text(H2_XYZ)
    path = tmp_path / "config.yaml"
    path.write_text(json.dumps(config))
    blocks = {}

    def capture(title, content, **_kwargs):
        blocks[title] = dict(content)
        return ""

    monkeypatch.setattr(opt_workflow, "pretty_block", capture)
    result = CliRunner().invoke(
        opt_workflow.cli,
        ["-i", str(xyz), "-q", "0", "--config", str(path), "--dry-run",
         "-o", str(tmp_path / "out"), *extra],
    )
    return result, blocks


@pytest.mark.parametrize(("kind", "mode"), [("lbfgs", "grad"), ("rfo", "hess")])
def test_opt_reads_nested_optimizer_section(tmp_path, monkeypatch, kind, mode):
    result, blocks = _opt_dry_run(
        tmp_path, monkeypatch, {"opt": {kind: {"max_cycles": 50}}}, "--opt-mode", mode,
    )
    assert result.exit_code == 0, result.output
    assert blocks["opt"]["max_cycles"] == 50
    assert kind not in blocks["opt"]


def test_opt_rejects_conflicting_spellings(tmp_path, monkeypatch):
    result, _ = _opt_dry_run(
        tmp_path, monkeypatch,
        {"lbfgs": {"max_cycles": 50}, "opt": {"lbfgs": {"max_cycles": 60}}},
        "--opt-mode", "grad",
    )
    assert result.exit_code == 2, result.output
    assert "lbfgs.max_cycles and opt.lbfgs.max_cycles conflict" in result.output


def test_opt_nested_value_counts_as_explicit(tmp_path, monkeypatch):
    result, _ = _opt_dry_run(
        tmp_path, monkeypatch,
        {"opt": {"max_cycles": 10, "lbfgs": {"max_cycles": 50}}},
        "--opt-mode", "grad",
    )
    assert result.exit_code == 2, result.output
    assert "opt.max_cycles and lbfgs.max_cycles conflict" in result.output


@pytest.mark.parametrize("kind", ["lbfgs", "rfo"])
def test_scan_reads_nested_optimizer_section(kind):
    _, _, opt_cfg, lbfgs_cfg, rfo_cfg, _ = build_scan_configs(
        {"opt": {kind: {"max_step": 0.05}}},
        kind=kind,
        geom_kw={},
        calc_kw={"workers": 1, "workers_per_node": 1},
        opt_kw={"thresh": "default"},
        lbfgs_kw={"max_step": 0.3},
        rfo_kw={"max_step": 0.3},
        bias_kw={"k": 100.0},
    )
    assert (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)["max_step"] == 0.05
    assert kind not in opt_cfg


def test_scan_rejects_explicit_opt_and_nested_conflict():
    with pytest.raises(click.BadParameter, match="opt.max_cycles and rfo.max_cycles conflict"):
        build_scan_configs(
            {"opt": {"max_cycles": 10, "rfo": {"max_cycles": 50}}},
            kind="rfo",
            geom_kw={},
            calc_kw={"workers": 1, "workers_per_node": 1},
            opt_kw={"max_cycles": 100},
            lbfgs_kw={"max_cycles": 100},
            rfo_kw={"max_cycles": 100},
            bias_kw={"k": 100.0},
        )


def test_tsopt_constructor_gets_no_nested_optimizer_section(monkeypatch, tmp_path):
    class ConstructorReached(RuntimeError):
        pass

    captured = {}

    class CaptureOptimizer:
        def __init__(self, geometry, **kwargs):
            captured.update(kwargs)
            raise ConstructorReached("ConstructorReached: no model evaluation requested.")

    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", CaptureOptimizer)
    monkeypatch.setattr(tsopt, "create_calculator", lambda **_kwargs: object())
    source = tmp_path / "input.xyz"
    source.write_text("1\nconstructor-only input\nHe 0 0 0\n")
    config = tmp_path / "config.yaml"
    config.write_text(json.dumps(
        {"opt": {"max_cycles": 7, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}}}
    ))
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(root_cli, [
        "tsopt", "-i", str(source), "-q", "0", "-m", "1", "-o", str(tmp_path / "out"),
        "--no-flatten", "--no-dump", "--config", str(config),
    ])
    assert result.exit_code == 1 and "ConstructorReached" in result.output, result.output
    assert captured["max_cycles"] == 7
    assert "lbfgs" not in captured and "rfo" not in captured


def _freq_dry_run(tmp_path: Path, monkeypatch, config: dict):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text(H2_XYZ)
    path = tmp_path / "config.yaml"
    path.write_text(json.dumps(config))
    blocks = {}

    def capture(title, content, **_kwargs):
        blocks[title] = dict(content)
        return ""

    monkeypatch.setattr(freq_workflow, "pretty_block", capture)
    result = CliRunner().invoke(root_cli, [
        "freq", "-i", str(xyz), "-q", "0", "--config", str(path), "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    return result, blocks


def test_freq_reads_nested_thermo_section(tmp_path, monkeypatch):
    result, blocks = _freq_dry_run(tmp_path, monkeypatch, {"freq": {"thermo": {"dump": True}}})
    assert result.exit_code == 0, result.output
    assert blocks["dry_run_plan"]["will_dump_thermo_yaml"] is True

    result, _ = _freq_dry_run(
        tmp_path, monkeypatch, {"freq": {"thermo": {"temperature": -5.0}}},
    )
    assert result.exit_code != 0
    assert "thermo.temperature must be a finite number greater than zero" in result.output


def test_freq_rejects_conflicting_thermo_spellings(tmp_path, monkeypatch):
    result, _ = _freq_dry_run(
        tmp_path, monkeypatch,
        {"thermo": {"temperature": 300.0}, "freq": {"thermo": {"temperature": 310.0}}},
    )
    assert result.exit_code != 0
    assert "thermo.temperature and freq.thermo.temperature conflict" in result.output


def test_all_thermo_preflight_reads_nested_thermo_section(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text(H2_XYZ)
    config = tmp_path / "config.yaml"
    config.write_text(json.dumps({"freq": {"thermo": {"temperature": -5.0}}}))
    result = CliRunner().invoke(root_cli, [
        "all", "-i", str(xyz), "-q", "0", "--tsopt", "--thermo",
        "--config", str(config), "--dry-run",
    ])
    assert result.exit_code != 0
    assert "thermo.temperature must be a finite number greater than zero" in result.output
