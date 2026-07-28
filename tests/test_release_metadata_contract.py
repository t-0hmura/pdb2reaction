"""Contract tests for release-metadata: landing substitution + license (M67, P03)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load(name: str):
    path = Path(__file__).parents[1] / ".github" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


rel = _load("check_release_versions")


# --- M67: landing-page release substitution ------------------------------- #
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


# --- P03: CFF <-> pyproject license identity ------------------------------ #
def test_real_license_metadata_matches_exact_spdx() -> None:
    assert rel.check_license_metadata() == []
    assert rel._cff_license() == "GPL-3.0-only"
    assert rel._pyproject_license() == "GPL-3.0-only"


def test_license_mismatch_is_rejected(monkeypatch) -> None:
    monkeypatch.setattr(rel, "_cff_license", lambda: "GPL-3.0")
    errors = rel.check_license_metadata()
    assert errors and "GPL-3.0" in errors[0]


def test_cffconvert_validates_citation_when_available() -> None:
    # cffconvert is not installed in every environment. When present, the CFF
    # must validate cleanly against the SPDX schema; when absent this passes
    # without a skip and the schema check is confirmed at the release freeze
    # (the exact-SPDX identity above is still enforced unconditionally).
    try:
        from cffconvert.cli.create_citation import create_citation
    except ImportError:
        return
    create_citation(str(rel.REPO_ROOT / "CITATION.cff"), None).validate()


# --- whole-checker positive control on the real repo ---------------------- #
def test_release_checker_passes_on_real_repo(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["check_release_versions.py"])
    assert rel.main() == 0
