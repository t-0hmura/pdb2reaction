#!/usr/bin/env python3
"""Static AST import-graph gate for pdb2reaction.

Enforces the product's real (post-C12) module-dependency invariants:

1. No file under ``pysisyphus/**`` (the bundled engine fork) imports the product
   package ``pdb2reaction``.
2. The product package has no strongly connected component among its own modules
   (no product import cycle).
3. Product layer back-edges ``core -> workflows`` and ``domain -> workflows`` are
   absent.

The graph is built from static imports without importing any optional MLIP
stack. It resolves absolute imports, relative imports, package ``__init__.py``,
aliases, and ``from . import child`` correctly: a ``from . import <submodule>``
is an edge to the SUBMODULE, never to the package ``__init__`` — otherwise every
sibling reference through a package would masquerade as a package<->submodule
cycle. Dynamic factory strings are out of scope (covered by runtime smoke tests).

Run as a script to check the repository (exit non-zero on any violation); import
``collect_violations`` for unit tests and negative controls.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------


def discover_modules(roots: Iterable[Path]) -> Dict[str, Path]:
    """Map dotted module name -> file path for every .py under each root package."""
    modules: Dict[str, Path] = {}
    for root in roots:
        root = Path(root)
        for path in root.rglob("*.py"):
            rel = path.relative_to(root.parent)
            parts = list(rel.with_suffix("").parts)
            if parts and parts[-1] == "__init__":
                parts = parts[:-1]
            if not parts:
                continue
            modules[".".join(parts)] = path
    return modules


def _package_of(modname: str, path: Path) -> str:
    """Package that contains ``modname`` (base for a relative import)."""
    if path.name == "__init__.py":
        return modname  # a package: `.` refers to itself
    return modname.rsplit(".", 1)[0] if "." in modname else ""


def _resolve_relative(pkg: str, level: int, module: str) -> str:
    base = pkg
    for _ in range(level - 1):
        base = base.rsplit(".", 1)[0] if "." in base else ""
    if module:
        return f"{base}.{module}" if base else module
    return base  # `from . import X` -> base is the package itself


def edges_for(modname: str, path: Path, modules: Dict[str, Path]) -> Set[str]:
    """Return the set of in-repo modules imported by this file."""
    tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"), filename=str(path))
    pkg = _package_of(modname, path)
    targets: Set[str] = set()

    def _add_longest_prefix(name: str) -> None:
        if name in modules:
            targets.add(name)
            return
        parts = name.split(".")
        for i in range(len(parts) - 1, 0, -1):
            prefix = ".".join(parts[:i])
            if prefix in modules:
                targets.add(prefix)
                return

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                _add_longest_prefix(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level and node.level > 0:
                base = _resolve_relative(pkg, node.level, node.module or "")
            else:
                base = node.module or ""
            if not base:
                continue
            # `from . import submodule` / `from .pkg import submodule` resolves to
            # the submodule; only genuine attribute imports pull in the package.
            base_is_package = base in modules and modules[base].name == "__init__.py"
            need_base_edge = not base_is_package
            for alias in node.names:
                child = f"{base}.{alias.name}"
                if child in modules:
                    targets.add(child)
                else:
                    need_base_edge = True
            if need_base_edge:
                _add_longest_prefix(base)

    targets.discard(modname)
    return targets


def build_graph(roots: Iterable[Path]) -> Tuple[Dict[str, Path], Dict[str, Set[str]]]:
    modules = discover_modules(roots)
    graph = {m: edges_for(m, p, modules) for m, p in modules.items()}
    return modules, graph


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def tarjan_scc(graph: Dict[str, Set[str]]) -> List[List[str]]:
    index_counter = [0]
    stack: List[str] = []
    lowlink: Dict[str, int] = {}
    index: Dict[str, int] = {}
    on_stack: Dict[str, bool] = {}
    result: List[List[str]] = []

    def strongconnect(v: str) -> None:
        work = [(v, 0)]
        while work:
            node, pi = work[-1]
            if pi == 0:
                index[node] = index_counter[0]
                lowlink[node] = index_counter[0]
                index_counter[0] += 1
                stack.append(node)
                on_stack[node] = True
            recurse = False
            succs = sorted(graph.get(node, ()))
            for j in range(pi, len(succs)):
                w = succs[j]
                if w not in index:
                    work[-1] = (node, j + 1)
                    work.append((w, 0))
                    recurse = True
                    break
                if on_stack.get(w):
                    lowlink[node] = min(lowlink[node], index[w])
            if recurse:
                continue
            if lowlink[node] == index[node]:
                comp: List[str] = []
                while True:
                    w = stack.pop()
                    on_stack[w] = False
                    comp.append(w)
                    if w == node:
                        break
                result.append(comp)
            work.pop()
            if work:
                parent = work[-1][0]
                lowlink[parent] = min(lowlink[parent], lowlink[node])

    for v in graph:
        if v not in index:
            strongconnect(v)
    return result


def _top(name: str) -> str:
    return name.split(".", 1)[0]


def _layer(name: str) -> str:
    parts = name.split(".")
    return parts[1] if len(parts) > 1 else ""


def collect_violations(
    roots: Iterable[Path],
    *,
    product_prefix: str,
    engine_prefix: str,
    forbidden_layer_edges: Iterable[Tuple[str, str]],
) -> Dict[str, List[str]]:
    """Return the three violation classes for the given package roots.

    Keys: ``sccs`` (product multi-module cycles), ``engine_to_product`` (engine
    fork importing the product), ``forbidden_edges`` (product layer back-edges).
    An all-empty result means the graph is clean.
    """
    modules, graph = build_graph(roots)
    forbidden_pairs = set(forbidden_layer_edges)

    product = {m for m in modules if _top(m) == product_prefix}
    product_graph = {m: {t for t in graph[m] if t in product} for m in product}
    sccs = [
        " <-> ".join(sorted(c))
        for c in tarjan_scc(product_graph)
        if len(c) > 1
    ]

    engine_to_product = [
        f"{m} -> {t}"
        for m in modules
        if _top(m) == engine_prefix
        for t in sorted(graph[m])
        if _top(t) == product_prefix
    ]

    forbidden_edges: List[str] = []
    for m in sorted(product):
        for t in sorted(graph[m]):
            if _top(t) != product_prefix:
                continue
            if (_layer(m), _layer(t)) in forbidden_pairs:
                forbidden_edges.append(f"{m} -> {t}")

    return {
        "sccs": sccs,
        "engine_to_product": engine_to_product,
        "forbidden_edges": forbidden_edges,
    }


# p2r-specific configuration
PRODUCT_PREFIX = "pdb2reaction"
ENGINE_PREFIX = "pysisyphus"
FORBIDDEN_LAYER_EDGES = (("core", "workflows"), ("domain", "workflows"))


def check_repository(repo_root: Path) -> List[str]:
    """Return a flat list of violation messages for the real repository."""
    roots = [repo_root / PRODUCT_PREFIX, repo_root / ENGINE_PREFIX]
    violations = collect_violations(
        roots,
        product_prefix=PRODUCT_PREFIX,
        engine_prefix=ENGINE_PREFIX,
        forbidden_layer_edges=FORBIDDEN_LAYER_EDGES,
    )
    messages: List[str] = []
    for cycle in violations["sccs"]:
        messages.append(f"product import cycle: {cycle}")
    for edge in violations["engine_to_product"]:
        messages.append(f"bundled engine imports product: {edge}")
    for edge in violations["forbidden_edges"]:
        messages.append(f"forbidden layer back-edge: {edge}")
    return messages


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    messages = check_repository(repo_root)
    if messages:
        print("Import-graph gate FAILED:")
        for m in messages:
            print(f"  - {m}")
        return 1
    print("Import-graph gate OK: no product cycle, layer back-edge, or engine->product import.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
