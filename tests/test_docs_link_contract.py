"""Contract tests for the public-markdown link checker (M66, M68)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load(name: str):
    path = Path(__file__).parents[1] / ".github" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # register before exec_module so the loaded script resolves its own module name
    spec.loader.exec_module(module)
    return module


links = _load("check_markdown_links")


def test_public_markdown_paths_includes_root_and_tree_pages() -> None:
    paths = links.public_markdown_paths()
    names = {p.name for p in paths}
    assert "README.md" in names
    assert "CONTRIBUTING.md" in names
    assert any(links.DOCS_ROOT in p.parents for p in paths)
    assert any((links.REPO_ROOT / "skills") in p.parents for p in paths)


def test_real_repo_links_resolve() -> None:
    assert links.main() == 0


def test_readme_example_target_resolves() -> None:
    # M68: the README advertises examples/run.sh, which must resolve.
    readme = links.REPO_ROOT / "README.md"
    errors: list[str] = []
    links._check_path(readme, errors)
    assert errors == []
    assert (links.REPO_ROOT / "examples" / "run.sh").exists()


def test_broken_link_in_readme_like_page_is_reported(tmp_path) -> None:
    page = tmp_path / "README.md"
    page.write_text("See [missing](does_not_exist.md).\n", encoding="utf-8")
    errors: list[str] = []
    links._check_path(page, errors)
    assert len(errors) == 1
    assert "broken local link" in errors[0]


def test_broken_links_in_skill_and_example_pages_reported_independently(tmp_path) -> None:
    for name in ("SKILL.md", "index.md"):
        page = tmp_path / name
        page.write_text("[x](./nope.md)\n", encoding="utf-8")
        errors: list[str] = []
        links._check_path(page, errors)
        assert errors, f"expected a broken-link error for {name}"
