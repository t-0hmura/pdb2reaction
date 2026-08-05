"""Unit tests for the structured error envelope produced by _write_error_json."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.cli.decorators import _write_error_json


class _CustomOptError(RuntimeError):
    pass


def test_error_envelope_includes_class_chain() -> None:
    with tempfile.TemporaryDirectory() as d:
        _write_error_json(
            Path(d),
            "opt",
            _CustomOptError("optimization diverged"),
            "opt-stage",
            time_start=None,
        )
        r = json.loads((Path(d) / "result.json").read_text())
        assert r["status"] == "error"
        assert r["error"] == "optimization diverged"
        assert r["error_type"] == "_CustomOptError"
        chain = r["error_class_chain"]
        assert "RuntimeError" in chain
        assert "Exception" in chain
        assert "BaseException" in chain
        assert chain[0] == "_CustomOptError"
        assert r["error_module"] == _CustomOptError.__module__
        assert r["error_label"] == "opt-stage"


def test_error_envelope_mirrors_to_summary_json() -> None:
    with tempfile.TemporaryDirectory() as d:
        _write_error_json(
            Path(d),
            "tsopt",
            ValueError("malformed input"),
            "tsopt-stage",
        )
        r = json.loads((Path(d) / "result.json").read_text())
        s = json.loads((Path(d) / "summary.json").read_text())
        assert r == s


def test_path_search_bad_parameter_uses_public_command_and_exit_two(
    tmp_path: Path,
) -> None:
    reactant = tmp_path / "r.xyz"
    product = tmp_path / "p.xyz"
    reactant.write_text("1\nr\nHe 0 0 0\n", encoding="utf-8")
    product.write_text("1\np\nHe 0 0 0.1\n", encoding="utf-8")
    out = tmp_path / "out"

    result = CliRunner().invoke(
        root_cli,
        [
            "path-search", "-i", str(reactant), "-i", str(product),
            "-q", "0", "--freeze-atoms", "bad", "-o", str(out),
        ],
    )

    assert result.exit_code == 2
    payload = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert payload["command"] == "path-search"
