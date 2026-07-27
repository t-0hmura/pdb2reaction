#!/usr/bin/env python3
"""Fail when release metadata and the bundled Colab ref disagree."""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = REPO_ROOT / "examples" / "pdb2reaction_colab.ipynb"
NOTEBOOK_REF = "pdb2reaction_version"
# The notebook's debug-install branch pins the same release a second time, via
# setuptools-scm. Gate it too, or it goes stale silently at the next release.
_NOTEBOOK_PRETEND_RE = re.compile(
    r"SETUPTOOLS_SCM_PRETEND_VERSION[\"']\]\s*=\s*[\"']v?([^\"']+)[\"']"
)

# EN/JA landing pages must render their version header via the MyST `release`
# substitution rather than a hardcoded literal (M67).
LANDING_PAGES = (REPO_ROOT / "docs" / "index.md", REPO_ROOT / "docs" / "ja" / "index.md")
LANDING_SUBSTITUTION = "v{{ release }}"
# A version header line (EN `*Version: ...*` / JA `*バージョン: ...*`) carrying a
# hardcoded `v<digits>` literal. The BibTeX `version = {...}` code block uses a
# different shape and is intentionally not matched here.
_LANDING_LITERAL_RE = re.compile(
    r"(?im)^\s*\*\s*(?:version|バージョン)\s*[:：].*v\d",
)


def _cff_version() -> str:
    text = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^version:\s*[\"']?([^\"'#\s]+)", text, re.MULTILINE)
    if match is None:
        raise ValueError("CITATION.cff has no top-level version field")
    return match.group(1)


def _cff_license() -> str:
    text = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^license:\s*[\"']?([^\"'#\s]+)", text, re.MULTILINE)
    if match is None:
        raise ValueError("CITATION.cff has no top-level license field")
    return match.group(1)


def _pyproject_license() -> str:
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r"^license\s*=\s*[\"']([^\"']+)[\"']", text, re.MULTILINE)
    if match is None:
        raise ValueError("pyproject.toml has no license expression")
    return match.group(1)


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def check_landing_pages() -> list[str]:
    """Require the EN/JA landing headers to use the MyST release substitution."""
    errors: list[str] = []
    for page in LANDING_PAGES:
        text = page.read_text(encoding="utf-8")
        if LANDING_SUBSTITUTION not in text:
            errors.append(
                f"{_rel(page)}: landing header must use "
                f"'{LANDING_SUBSTITUTION}' (MyST release substitution)"
            )
        if _LANDING_LITERAL_RE.search(text):
            errors.append(
                f"{_rel(page)}: landing header has a hardcoded "
                f"version literal; use '{LANDING_SUBSTITUTION}'"
            )
    return errors


def check_license_metadata() -> list[str]:
    """Require the CFF license to equal pyproject's exact SPDX expression (P03)."""
    cff = _cff_license()
    pyproject = _pyproject_license()
    if cff != pyproject:
        return [
            f"CITATION.cff license '{cff}' != pyproject.toml license '{pyproject}' "
            "(use the exact SPDX identifier in both)"
        ]
    return []


def _docs_release() -> str:
    tree = ast.parse((REPO_ROOT / "docs" / "conf.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "release"
            for target in node.targets
        ):
            value = ast.literal_eval(node.value)
            return str(value)
    raise ValueError("docs/conf.py has no literal release assignment")


def _notebook_version() -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    pattern = re.compile(rf"^{NOTEBOOK_REF}\s*=\s*[\"']v?([^\"']+)[\"']", re.MULTILINE)
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        match = pattern.search(str(cell.get("source", "")))
        if match is not None:
            return match.group(1)
    raise ValueError(f"{NOTEBOOK.name} has no {NOTEBOOK_REF} assignment")


def _notebook_pretend_version() -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        match = _NOTEBOOK_PRETEND_RE.search(str(cell.get("source", "")))
        if match is not None:
            return match.group(1)
    raise ValueError(
        f"{NOTEBOOK.name} has no SETUPTOOLS_SCM_PRETEND_VERSION assignment"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-version")
    args = parser.parse_args()

    values = {
        "CITATION.cff": _cff_version(),
        "docs/conf.py": _docs_release(),
        NOTEBOOK.name: _notebook_version(),
        f"{NOTEBOOK.name} (setuptools-scm)": _notebook_pretend_version(),
    }
    expected = str(args.expected_version or values["CITATION.cff"]).removeprefix("v")
    errors: list[str] = []
    mismatches = {name: value for name, value in values.items() if value != expected}
    if mismatches:
        details = ", ".join(f"{name}={value}" for name, value in mismatches.items())
        errors.append(f"expected {expected}; mismatched: {details}")
    errors.extend(check_landing_pages())
    errors.extend(check_license_metadata())
    if errors:
        raise SystemExit("[release-version] failed:\n  " + "\n  ".join(errors))
    print(
        f"[release-version] OK: {expected} ({', '.join(values)}); "
        f"landing headers use {LANDING_SUBSTITUTION}; "
        f"license={_pyproject_license()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
