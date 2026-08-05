"""Resolved calculator provenance used by the actual opt result payload."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from click.testing import CliRunner

from pdb2reaction.backends import apply_calc_file_to_calc_cfg
from pdb2reaction.workflows import opt as opt_workflow
from pdb2reaction.workflows.opt import _opt_result_provenance


def test_yaml_resolved_orb_fields_keep_legacy_and_canonical_parity() -> None:
    # This is the final calc_cfg shape after a YAML-only ORB selection; no raw
    # Click callback value is passed to the result seam.
    calc_cfg = {
        "backend": "orb",
        "model": "orb_v3_conservative_omol",
        "precision": "fp64",
    }
    original = deepcopy(calc_cfg)

    result = _opt_result_provenance(calc_cfg)

    assert result == {
        "backend": "orb",
        "model": "orb_v3_conservative_omol",
        "mlip_backend": "orb",
        "mlip_model": "orb_v3_conservative_omol",
        "mlip_precision": "fp64",
    }
    assert calc_cfg == original


def test_custom_fixture_fields_keep_legacy_and_canonical_parity(
    tmp_path: Path,
) -> None:
    calc_file = tmp_path / "toy.py"
    calc_file.write_text("def get_calculator(**kwargs):\n    return None\n")
    calc_cfg = {"backend": "uma", "model": "uma-s-1p1"}
    apply_calc_file_to_calc_cfg(calc_cfg, str(calc_file), "get_calculator")
    original = deepcopy(calc_cfg)

    result = _opt_result_provenance(calc_cfg)

    assert result["backend"] == result["mlip_backend"] == "custom"
    assert result["model"] == result["mlip_model"] == "toy.py:get_calculator"
    assert result["mlip_precision"] is None
    assert calc_cfg == original


def test_yaml_opt_mode_selects_optimizer_unless_cli_is_explicit(tmp_path: Path) -> None:
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n")
    config = tmp_path / "config.yaml"
    config.write_text("opt:\n  opt_mode: lbfgs\n")

    base_args = [
        "-i", str(xyz), "-q", "0", "--config", str(config), "--dry-run",
    ]
    yaml_only = CliRunner().invoke(
        opt_workflow.cli,
        [*base_args, "-o", str(tmp_path / "yaml")],
    )
    explicit = CliRunner().invoke(
        opt_workflow.cli,
        [*base_args, "--opt-mode", "rfo", "-o", str(tmp_path / "cli")],
    )

    assert yaml_only.exit_code == 0
    assert "[opt] lbfgs," in yaml_only.output
    assert explicit.exit_code == 0
    assert "[opt] rfo," in explicit.output
