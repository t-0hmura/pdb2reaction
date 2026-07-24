from __future__ import annotations

import ast
import importlib.util
from pathlib import Path


def _load_checker():
    path = (
        Path(__file__).parents[1]
        / ".github"
        / "scripts"
        / "check_all_scan_contract.py"
    )
    spec = importlib.util.spec_from_file_location("p2r_scan_contract", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_checker_collects_only_runtime_keywords_forwarded_to_scan() -> None:
    module = _load_checker()
    tree = ast.parse(
        """
def cli():
    scan_args: list[str] = ["--positive" if enabled else "--negative"]
    scan_args.extend(["--extended", value])
    _append_explicit_child_runtime_argv(
        scan_args,
        workers=workers,
        workers_per_node=workers_per_node,
        dmf_backend=None,
    )

def _forward_calc_file_argv(child_args, cfg):
    child_args.extend(["--from-helper", cfg])

def _append_explicit_child_runtime_argv(
    child_args, *, workers=None, workers_per_node=None, dmf_backend=None
):
    if workers is not None:
        child_args.extend(["--workers", str(workers)])
    if workers_per_node is not None:
        child_args.extend(["--workers-per-node", str(workers_per_node)])
    if dmf_backend is not None:
        child_args.extend(["--dmf-backend", dmf_backend])
"""
    )
    collector = module._ForwardedScanFlags()
    collector.visit(tree)
    assert collector.flags == {
        "--extended",
        "--from-helper",
        "--negative",
        "--positive",
        "--workers",
        "--workers-per-node",
    }
