from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import pytest


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


class ReachedScan(BaseException):
    """Escape the parent pipeline at the scan dispatch."""


@pytest.mark.parametrize(
    "extra",
    [
        ["--dump", "--scan-zero-based", "--scan-preopt", "--scan-endopt",
         "--scan-max-step-size", "0.2", "--scan-bias-k", "100",
         "--scan-relax-max-cycles", "50", "--thresh", "gau", "--print-every", "5",
         "--backend", "uma", "--workers", "2", "--workers-per-node", "1",
         "--solvent", "water", "--solvent-model", "alpb"],
        ["--no-dump", "--scan-one-based", "--no-scan-preopt", "--no-scan-endopt"],
    ],
    ids=["set", "unset"],
)
def test_all_scan_argv_is_accepted_by_the_scan_command(tmp_path, monkeypatch, extra) -> None:
    """Every option ``all`` forwards to ``scan`` parses with the real scan command."""
    from click.testing import CliRunner

    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.workflows import all as all_workflow
    from pdb2reaction.workflows import scan as scan_workflow
    from test_all_root_summary_golden import pdb_frame, write

    captured = []

    def first_child(name, _cli, args, **_kwargs):
        captured.append((name, list(args)))
        raise ReachedScan

    monkeypatch.setattr(all_workflow, "_run_cli_main", first_child)
    source = write(tmp_path / "input0.pdb", pdb_frame(0))
    scan_lists = "[(0,1,1.97)]" if "--scan-zero-based" in extra else "[(1,2,1.97)]"
    args = ["all", "-i", str(source), "-q", "0", "-m", "1",
            "--out-dir", str(tmp_path / "out"), "--no-preopt", "--no-freeze-links",
            "--scan-lists", scan_lists, *extra]
    with pytest.raises(ReachedScan):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    ((name, scan_args),) = captured
    assert name == "scan"
    context = scan_workflow.cli.make_context("scan", list(scan_args))
    # The scan command keeps unknown tokens in ctx.args; none may be an option.
    assert [token for token in context.args if token.startswith("-")] == []
    if "--thresh" in extra:
        assert scan_args[scan_args.index("--thresh") + 1] == "gau"
