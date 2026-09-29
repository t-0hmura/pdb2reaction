"""YAML routing of single-structure optimizer settings in path-opt and path-search."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli

SMOKE = Path(__file__).resolve().parent / "smoke"
CALC_YAML = "calc:\n  charge: -1\n  spin: 1\n"


def _invoke(tmp_path: Path, command: str, yaml_text: str, *extra: str):
    config = tmp_path / "config.yaml"
    config.write_text(CALC_YAML + yaml_text, encoding="utf-8")
    structure = SMOKE / "r.pdb"
    return CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(structure), str(structure),
            "--config", str(config), "--dry-run",
            "--out-dir", str(tmp_path / command), *extra,
        ],
    )


def _block(output: str, title: str) -> str:
    head = f"\n{title}\n{'-' * len(title)}\n"
    start = output.index(head) + len(head)
    end = output.find("\n\n", start)
    return output[start:] if end < 0 else output[start:end]


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
@pytest.mark.parametrize(
    ("yaml_text", "message"),
    [
        (
            "opt:\n  thresh: gau\nlbfgs:\n  thresh: baker\n",
            "opt.thresh and lbfgs.thresh conflict",
        ),
        (
            "opt:\n  thresh: gau\n  lbfgs:\n    thresh: baker\n",
            "opt.thresh and opt.lbfgs.thresh conflict",
        ),
    ],
)
def test_conflicting_opt_and_lbfgs_values_are_rejected(
    tmp_path: Path, command: str, yaml_text: str, message: str
) -> None:
    result = _invoke(tmp_path, command, yaml_text)

    assert result.exit_code == 2, result.output
    assert message in result.output


def test_stopt_lbfgs_configures_preoptimization_only(tmp_path: Path) -> None:
    result = _invoke(
        tmp_path, "path-opt", "stopt:\n  lbfgs:\n    max_cycles: 7\n",
        "--show-config", "-v", "3",
    )

    assert result.exit_code == 0, result.output
    assert "preopt_max_cycles: 7" in _block(result.output, "dry_run_plan")
    stopt_keys = [line.split(":", 1)[0] for line in _block(result.output, "stopt").splitlines()]
    assert "lbfgs" not in stopt_keys


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_conflicting_lbfgs_and_stopt_lbfgs_values_are_rejected(
    tmp_path: Path, command: str
) -> None:
    result = _invoke(
        tmp_path, command,
        "lbfgs:\n  max_cycles: 5\nstopt:\n  lbfgs:\n    max_cycles: 7\n",
    )

    assert result.exit_code == 2, result.output
    assert "lbfgs.max_cycles and stopt.lbfgs.max_cycles conflict" in result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_non_mapping_opt_section_is_rejected(tmp_path: Path, command: str) -> None:
    result = _invoke(tmp_path, command, "opt: 5\n")

    assert result.exit_code == 2, result.output
    assert "YAML section 'opt' must be a mapping" in result.output


def test_apply_single_opt_yaml_layer_keeps_optimizer_keys_out_of_stopt() -> None:
    from pdb2reaction.core.utils import deep_update
    from pdb2reaction.workflows._path_yaml_helpers import apply_single_opt_yaml_layer

    layer = {
        "opt": {"thresh": "gau_tight"},
        "stopt": {"max_cycles": 40, "lbfgs": {"max_cycles": 7}},
    }
    stopt_cfg = {"max_cycles": 40, "lbfgs": {"max_cycles": 7}}
    lbfgs_cfg = {"thresh": "gau", "max_cycles": 100}
    rfo_cfg = {"thresh": "gau", "max_cycles": 100}

    apply_single_opt_yaml_layer(
        layer,
        lbfgs_cfg=lbfgs_cfg,
        rfo_cfg=rfo_cfg,
        stopt_cfg=stopt_cfg,
        opt_base_kw={"thresh": "gau", "max_cycles": 100},
        deep_update=deep_update,
    )

    assert stopt_cfg == {"max_cycles": 40}
    assert lbfgs_cfg == {"thresh": "gau_tight", "max_cycles": 7}
    assert rfo_cfg == {"thresh": "gau_tight", "max_cycles": 100}
