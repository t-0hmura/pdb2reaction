"""Enforceable AST import-graph gate with negative controls.

The gate asserts three real invariants of the product graph:
 - no bundled ``pysisyphus/**`` file imports ``pdb2reaction``;
 - the product package has no strongly connected component of its own modules;
 - the product layer back-edges ``core -> workflows`` / ``domain -> workflows``
   are absent.

The negative controls prove each check actually fires on a known-bad graph
(so silently dropping a check breaks a focused test) and that the relative-import
resolver does not misclassify a package<->submodule reference as a cycle.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_CHECKER_PATH = _REPO_ROOT / ".github" / "scripts" / "check_import_graph.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_import_graph", _CHECKER_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


checker = _load_checker()


def _write(base: Path, rel: str, text: str = "") -> None:
    path = base / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# The real repository must pass the gate.
# ---------------------------------------------------------------------------


def test_real_repository_passes_import_graph_gate():
    messages = checker.check_repository(_REPO_ROOT)
    assert messages == [], "\n".join(messages)


# ---------------------------------------------------------------------------
# Known-bad control: every check must fire.  Removing any single sub-check
# would leave one of these categories empty and fail this focused test.
# ---------------------------------------------------------------------------


def test_gate_fires_on_cycle_forbidden_edge_and_engine_import(tmp_path):
    prod = tmp_path / "fakeproduct"
    eng = tmp_path / "fakeengine"

    _write(prod, "__init__.py")
    _write(prod, "core/__init__.py")
    # core -> workflows (forbidden back-edge) and a product cycle a <-> b
    _write(prod, "core/utils.py",
           "from fakeproduct.workflows.run import go\n"
           "from fakeproduct.core.a import A\n")
    _write(prod, "core/a.py", "from fakeproduct.core.b import B\n")
    _write(prod, "core/b.py", "from fakeproduct.core.a import A\n")
    _write(prod, "workflows/__init__.py")
    _write(prod, "workflows/run.py", "go = 1\n")

    _write(eng, "__init__.py")
    _write(eng, "mod.py", "from fakeproduct.core.utils import x\n")  # engine -> product

    v = checker.collect_violations(
        [prod, eng],
        product_prefix="fakeproduct",
        engine_prefix="fakeengine",
        forbidden_layer_edges=[("core", "workflows")],
    )
    assert v["sccs"], "SCC check did not fire on the a<->b cycle"
    assert any("fakeproduct.core.a" in s and "fakeproduct.core.b" in s for s in v["sccs"])
    assert v["engine_to_product"], "engine->product check did not fire"
    assert any("fakeengine.mod -> fakeproduct.core.utils" == e for e in v["engine_to_product"])
    assert v["forbidden_edges"], "forbidden layer back-edge check did not fire"
    assert any("core.utils -> fakeproduct.workflows.run" in e for e in v["forbidden_edges"])


# ---------------------------------------------------------------------------
# Resolver correctness controls.
# ---------------------------------------------------------------------------


def test_package_to_submodule_relative_import_is_not_a_cycle(tmp_path):
    """`from . import child` (package re-export) + `from . import sibling`
    (sibling reference) must NOT read as a package<->submodule cycle."""
    pkg = tmp_path / "good"
    _write(pkg, "__init__.py", "from . import child\n")
    _write(pkg, "child.py", "from . import sibling\n")
    _write(pkg, "sibling.py", "VALUE = 1\n")

    v = checker.collect_violations(
        [pkg],
        product_prefix="good",
        engine_prefix="unused",
        forbidden_layer_edges=[],
    )
    assert v["sccs"] == [], f"false cycle reported: {v['sccs']}"


def test_genuine_relative_submodule_cycle_is_detected(tmp_path):
    """A real two-submodule cycle expressed with relative imports IS caught."""
    pkg = tmp_path / "good"
    _write(pkg, "__init__.py")
    _write(pkg, "child.py", "from . import sibling\n")
    _write(pkg, "sibling.py", "from . import child\n")

    v = checker.collect_violations(
        [pkg],
        product_prefix="good",
        engine_prefix="unused",
        forbidden_layer_edges=[],
    )
    assert any("good.child" in s and "good.sibling" in s for s in v["sccs"]), v["sccs"]


def test_clean_layered_graph_reports_nothing(tmp_path):
    """A well-layered product (workflows -> core, no back-edge) is clean."""
    prod = tmp_path / "fakeproduct"
    _write(prod, "__init__.py")
    _write(prod, "core/__init__.py")
    _write(prod, "core/utils.py", "VALUE = 1\n")
    _write(prod, "workflows/__init__.py")
    _write(prod, "workflows/run.py", "from fakeproduct.core.utils import VALUE\n")

    v = checker.collect_violations(
        [prod],
        product_prefix="fakeproduct",
        engine_prefix="unused",
        forbidden_layer_edges=[("core", "workflows"), ("domain", "workflows")],
    )
    assert v == {"sccs": [], "engine_to_product": [], "forbidden_edges": []}
