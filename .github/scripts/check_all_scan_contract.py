#!/usr/bin/env python3
"""Verify that options forwarded from `pdb2reaction all` to `pdb2reaction scan` exist."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Set


REPO_ROOT = Path(__file__).resolve().parents[2]
ALL_PY = REPO_ROOT / "pdb2reaction" / "workflows" / "all.py"
SCAN_PY = REPO_ROOT / "pdb2reaction" / "workflows" / "scan.py"
SCAN_COMMON_PY = REPO_ROOT / "pdb2reaction" / "workflows" / "scan_common.py"
COMMON_OPTIONS_PY = REPO_ROOT / "pdb2reaction" / "cli" / "common_options.py"


def _is_click_option_call(node: ast.Call) -> bool:
    func = node.func
    return (
        isinstance(func, ast.Attribute)
        and func.attr == "option"
        and isinstance(func.value, ast.Name)
        and func.value.id == "click"
    )


def _const_flag(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        value = node.value.strip()
        if value.startswith("--"):
            return value
    return None


def _const_flags(node: ast.AST) -> Set[str]:
    """Return every literal long option contained in an expression."""
    flags: Set[str] = set()
    for child in ast.walk(node):
        flag = _const_flag(child)
        if flag is not None:
            flags.add(flag)
    return flags


def _expand_flag_spec(spec: str) -> Set[str]:
    parts = [p.strip() for p in spec.split("/") if p.strip()]
    flags = {p for p in parts if p.startswith("--")}
    if flags:
        return flags
    return {spec}


def _collect_scan_declared_flags(scan_tree: ast.AST) -> Set[str]:
    flags: Set[str] = set()
    for node in ast.walk(scan_tree):
        if not isinstance(node, ast.Call) or not _is_click_option_call(node):
            continue
        for arg in node.args:
            flag = _const_flag(arg)
            if flag is not None:
                flags.update(_expand_flag_spec(flag))
    return flags


class _ForwardedScanFlags(ast.NodeVisitor):
    def __init__(self) -> None:
        self.targets: Set[str] = set()
        self.flags: Set[str] = set()

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        previous = self.targets
        self.targets = {
            "cli": {"scan_args"},
            "_forward_calc_file_argv": {"child_args"},
        }.get(node.name, set())
        self.generic_visit(node)
        self.targets = previous

    def visit_Assign(self, node: ast.Assign) -> None:
        if self.targets:
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in self.targets:
                    self.flags.update(_const_flags(node.value))
        self.generic_visit(node)

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        if (
            self.targets
            and isinstance(node.target, ast.Name)
            and node.target.id in self.targets
            and node.value is not None
        ):
            self.flags.update(_const_flags(node.value))
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        if self.targets:
            if (
                isinstance(node.func, ast.Name)
                and node.func.id == "_append_explicit_child_runtime_argv"
                and node.args
                and isinstance(node.args[0], ast.Name)
                and node.args[0].id in self.targets
            ):
                runtime_flags = {
                    "workers": "--workers",
                    "workers_per_node": "--workers-per-node",
                    "dmf_backend": "--dmf-backend",
                    "backend": "--backend",
                    "solvent": "--solvent",
                    "solvent_model": "--solvent-model",
                }
                for keyword in node.keywords:
                    if keyword.arg not in runtime_flags:
                        continue
                    if (
                        isinstance(keyword.value, ast.Constant)
                        and keyword.value.value is None
                    ):
                        continue
                    self.flags.add(runtime_flags[keyword.arg])

            if isinstance(node.func, ast.Name) and node.func.id == "_append_cli_arg":
                if (
                    len(node.args) >= 2
                    and isinstance(node.args[0], ast.Name)
                    and node.args[0].id in self.targets
                ):
                    self.flags.update(_const_flags(node.args[1]))

            if (
                isinstance(node.func, ast.Attribute)
                and node.func.attr in {"append", "extend"}
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id in self.targets
            ):
                for arg in node.args:
                    self.flags.update(_const_flags(arg))
        self.generic_visit(node)


def main() -> int:
    all_tree = ast.parse(ALL_PY.read_text(encoding="utf-8"), filename=str(ALL_PY))
    scan_tree = ast.parse(SCAN_PY.read_text(encoding="utf-8"), filename=str(SCAN_PY))
    scan_common_tree = ast.parse(
        SCAN_COMMON_PY.read_text(encoding="utf-8"),
        filename=str(SCAN_COMMON_PY),
    )
    common_options_tree = ast.parse(
        COMMON_OPTIONS_PY.read_text(encoding="utf-8"),
        filename=str(COMMON_OPTIONS_PY),
    )

    declared = _collect_scan_declared_flags(scan_tree)
    # `scan` receives many flags through @add_scan_common_options(...).
    declared.update(_collect_scan_declared_flags(scan_common_tree))
    # Calculator-file options are attached by @add_calc_file_option.
    declared.update(_collect_scan_declared_flags(common_options_tree))
    collector = _ForwardedScanFlags()
    collector.visit(all_tree)
    forwarded = collector.flags

    missing = sorted(flag for flag in forwarded if flag not in declared)
    if missing:
        print("Detected all->scan option contract drift.")
        for flag in missing:
            print(f"  missing in scan CLI: {flag}")
        return 1

    print(f"[scan-contract] all forwards {len(forwarded)} scan options; all are declared.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
