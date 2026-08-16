import json

import click
import pytest
from click.testing import CliRunner

from pdb2reaction.cli.app import cli


def test_def2_ecp_is_only_auto_enabled_for_covered_elements() -> None:
    from pdb2reaction.workflows.dft import _def2_ecp_required

    atomic_numbers = {"H": 1, "C": 6, "O": 8, "Rb": 37, "I": 53}
    charge = atomic_numbers.__getitem__

    assert not _def2_ecp_required(
        [("H", (0, 0, 0)), ("C", (0, 0, 1)), ("O", (0, 1, 0))],
        charge,
    )
    assert _def2_ecp_required([("Rb", (0, 0, 0)), ("I", (0, 0, 3))], charge)


def test_dft_rejects_unknown_yaml_engine_before_execution(tmp_path) -> None:
    xyz = tmp_path / "h2.xyz"
    xyz.write_text(
        "2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n",
        encoding="utf-8",
    )
    config = tmp_path / "config.yaml"
    config.write_text("dft:\n  engine: quantum-potato\n", encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "dft",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
        ],
    )

    assert result.exit_code == 2
    assert "dft.engine must be either 'cpu' or 'gpu'" in result.output


def test_dft_nonconvergence_commits_json_before_exit(tmp_path) -> None:
    from pdb2reaction.workflows.dft import _finalize_dft_result

    with pytest.raises(SystemExit) as caught:
        _finalize_dft_result(
            out_json=True,
            out_dir=tmp_path,
            payload={"status": "not_converged", "converged": False},
            elapsed_seconds=1.0,
        )

    assert caught.value.code == 3
    payload = json.loads((tmp_path / "result.json").read_text())
    assert payload["status"] == "not_converged"
    assert payload["converged"] is False


def test_prepare_dft_output_dir_invalidates_prior_public_results(tmp_path) -> None:
    from pdb2reaction.workflows.dft import _prepare_dft_output_dir

    for name in ("result.yaml", "result.json", "summary.json"):
        (tmp_path / name).write_text("stale\n", encoding="utf-8")

    _prepare_dft_output_dir(tmp_path)

    assert all(
        not (tmp_path / name).exists()
        for name in ("result.yaml", "result.json", "summary.json")
    )


def test_prepare_dft_output_dir_rejects_input_alias_before_mutation(tmp_path) -> None:
    from pdb2reaction.workflows.dft import _prepare_dft_output_dir

    result = tmp_path / "result.yaml"
    result.write_text("input: retained\n", encoding="utf-8")
    alias = tmp_path / "input.yaml"
    alias.hardlink_to(result)

    with pytest.raises(click.UsageError, match="collides with reserved DFT output"):
        _prepare_dft_output_dir(tmp_path, protected_inputs=(alias,))

    assert result.read_text(encoding="utf-8") == "input: retained\n"


def test_dft_unexpected_config_failure_uses_yaml_effective_output(
    monkeypatch,
    tmp_path,
) -> None:
    from pdb2reaction.workflows import dft

    xyz = tmp_path / "h2.xyz"
    xyz.write_text(
        "2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n",
        encoding="utf-8",
    )
    effective_out = tmp_path / "configured-output"
    config = tmp_path / "config.yaml"
    config.write_text(
        "dft:\n"
        f"  out_dir: {effective_out}\n"
        "  func_basis: wb97m-v/def2-tzvpd\n",
        encoding="utf-8",
    )

    def fail_after_output_resolution(_value):
        raise RuntimeError("config probe failed")

    monkeypatch.setattr(dft, "_parse_func_basis", fail_after_output_resolution)
    result = CliRunner().invoke(
        cli,
        [
            "dft",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
        ],
    )

    assert result.exit_code == 1
    payload = json.loads((effective_out / "result.json").read_text())
    assert payload["status"] == "error"
    assert payload["command"] == "dft"
    assert payload["error"] == "config probe failed"


@pytest.mark.parametrize(
    ("label", "config_text", "expected"),
    [
        ("syntax", "calc: [1,\n", "invalid YAML"),
        ("root", "- 1\n- 2\n", "must be a mapping"),
        ("section", "geom: 5\n", "YAML section 'geom' must be a mapping"),
    ],
)
def test_malformed_config_yaml_is_an_input_error(
    tmp_path, label, config_text, expected
) -> None:
    """A malformed --config file exits 2 with a message, not a raw traceback.

    The YAML layer is loaded before each command's own exception rendering, so
    the loader and the section resolver own these diagnostics.
    """
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", encoding="utf-8")
    config = tmp_path / f"{label}.yaml"
    config.write_text(config_text, encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "sp",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
            "-o",
            str(tmp_path / f"out_{label}"),
        ],
    )

    assert result.exit_code == 2
    assert expected in result.output
    assert "Traceback" not in result.output
