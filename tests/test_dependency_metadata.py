"""Dependency constraints backed by current runtime/API requirements."""

from __future__ import annotations

from pathlib import Path
import tomllib


def _project() -> dict:
    path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    return tomllib.loads(path.read_text(encoding="utf-8"))["project"]


def test_runtime_dependency_floors_match_consumed_apis() -> None:
    project = _project()
    dependencies = set(project["dependencies"])
    extras = project["optional-dependencies"]

    assert "pydmf>=1.2" in dependencies
    assert "plotly>=6.1.1" in dependencies
    assert extras["orb"] == ["orb-models>=0.7.0"]
    assert extras["aimnet"] == ["aimnet>=0.2.0"]
    assert extras["mcp"] == ["mcp[cli]>=1.29,<2"]
