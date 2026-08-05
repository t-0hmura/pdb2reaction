"""Contract tests for release metadata, landing substitution, and licensing."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import tomllib
from packaging.markers import default_environment
from packaging.requirements import Requirement


def _load(name: str):
    path = Path(__file__).parents[1] / ".github" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


rel = _load("check_release_versions")
REPO_ROOT = Path(__file__).parents[1]


def test_dft_gpu_dependencies_are_linux_only() -> None:
    data = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    requirements = {
        req.name: req
        for req in map(Requirement, data["project"]["optional-dependencies"]["dft"])
    }
    for package in ("gpu4pyscf-cuda12x", "cupy-cuda12x"):
        marker = requirements[package].marker
        assert marker is not None
        for platform, expected in (("linux", True), ("darwin", False), ("win32", False)):
            env = default_environment()
            env.update({"sys_platform": platform, "platform_machine": "x86_64"})
            assert marker.evaluate(env) is expected


def test_mcp_environment_override_sample_is_posix_only() -> None:
    sample = json.loads(
        (REPO_ROOT / "examples" / "mcp_client_config.json").read_text(
            encoding="utf-8"
        )
    )
    assert sample["_comment"].startswith("POSIX sample")
    assert "%APPDATA%" not in sample["_comment"]
    assert "/bin:" in sample["mcpServers"]["pdb2reaction"]["env"]["PATH"]


def test_public_single_model_examples_have_balanced_terminators() -> None:
    for filename in ("1.R.pdb", "2.IM.pdb"):
        lines = (REPO_ROOT / "examples" / filename).read_text(
            encoding="utf-8"
        ).splitlines()
        assert lines[-1] == "END"
        assert not any(line.startswith(("MODEL", "ENDMDL")) for line in lines)


# --- landing-page release substitution ------------------------------- #
def test_real_landing_pages_use_release_substitution() -> None:
    assert rel.check_landing_pages() == []
    for page in rel.LANDING_PAGES:
        assert rel.LANDING_SUBSTITUTION in page.read_text(encoding="utf-8")


def test_hardcoded_landing_literal_is_rejected(tmp_path, monkeypatch) -> None:
    page = tmp_path / "index.md"
    page.write_text("# Title\n\n*Version: v0.4.12*\n", encoding="utf-8")
    monkeypatch.setattr(rel, "LANDING_PAGES", (page,))
    errors = rel.check_landing_pages()
    assert errors
    assert any("hardcoded" in e or "must use" in e for e in errors)


def test_misspelled_substitution_is_rejected(tmp_path, monkeypatch) -> None:
    page = tmp_path / "index.md"
    page.write_text("# Title\n\n*Version: v{{ relese }}*\n", encoding="utf-8")
    monkeypatch.setattr(rel, "LANDING_PAGES", (page,))
    errors = rel.check_landing_pages()
    assert errors and any("must use" in e for e in errors)


# --- CFF <-> pyproject license identity ---------------------------------- #
def test_real_license_metadata_matches_exact_spdx() -> None:
    assert rel.check_license_metadata() == []
    assert rel._cff_license() == "GPL-3.0-only"
    assert rel._pyproject_license() == "GPL-3.0-only"


def test_license_mismatch_is_rejected(monkeypatch) -> None:
    monkeypatch.setattr(rel, "_cff_license", lambda: "GPL-3.0")
    errors = rel.check_license_metadata()
    assert errors and "GPL-3.0" in errors[0]


def test_cffconvert_validates_citation_when_available() -> None:
    # cffconvert is optional. When available, validate the CFF against its
    # schema; exact SPDX identity above is enforced independently.
    try:
        from cffconvert.cli.create_citation import create_citation
    except ImportError:
        return
    create_citation(str(rel.REPO_ROOT / "CITATION.cff"), None).validate()


# --- whole-checker positive control on the real repo ---------------------- #
def test_release_checker_passes_on_real_repo(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["check_release_versions.py"])
    assert rel.main() == 0
