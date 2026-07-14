"""MCP tool-schema and argv-construction regressions.

These tests do not launch calculators or the MCP transport.  A tiny fake MCP
registry captures the nested tool functions, while ``run_subcmd`` is replaced
with a recorder so we can assert the exact CLI contract agents will execute.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from pdb2reaction.mcp import _tools


class _FakeMCP:
    def __init__(self) -> None:
        self.tools: dict[str, object] = {}

    def tool(self):
        def decorator(func):
            self.tools[func.__name__] = func
            return func

        return decorator


class _FakeResult:
    def __init__(self, argv: list[str], out_dir: Path | None) -> None:
        self.argv = argv
        self.out_dir = out_dir

    def to_dict(self) -> dict:
        return {
            "status": "ok",
            "argv": self.argv,
            "out_dir": str(self.out_dir) if self.out_dir else None,
        }


@pytest.fixture()
def registered_tools(monkeypatch):
    calls: list[dict] = []

    def fake_run_subcmd(argv, *, out_dir=None, timeout=None, **_kwargs):
        record = {
            "argv": list(argv),
            "out_dir": out_dir,
            "timeout": timeout,
        }
        calls.append(record)
        return _FakeResult(record["argv"], out_dir)

    monkeypatch.setattr(_tools, "run_subcmd", fake_run_subcmd)
    mcp = _FakeMCP()
    _tools.register_all(mcp)
    return mcp.tools, calls


def test_registers_one_tool_per_cli_subcommand(registered_tools) -> None:
    tools, _calls = registered_tools
    assert len(tools) == 18


def test_charge_is_optional_when_ligand_mapping_derives_it(
    registered_tools, tmp_path: Path
) -> None:
    tools, calls = registered_tools

    tools["optimize_geometry"](
        "model.pdb",
        ligand_charge="LIG:-1",
        out_dir=str(tmp_path / "opt"),
    )

    argv = calls[-1]["argv"]
    assert "-q" not in argv
    assert argv[argv.index("-l") + 1] == "LIG:-1"


def test_search_paths_always_passes_two_ordered_endpoints(
    registered_tools, tmp_path: Path
) -> None:
    tools, calls = registered_tools
    signature = inspect.signature(tools["search_paths"])
    assert signature.parameters["product_pdb"].default is inspect.Parameter.empty

    tools["search_paths"](
        "R.pdb",
        product_pdb="P.pdb",
        intermediate_pdbs=["IM1.pdb", "IM2.pdb"],
        ligand_charge="LIG:-1",
        out_dir=str(tmp_path / "path-search"),
    )

    argv = calls[-1]["argv"]
    input_at = argv.index("-i")
    assert argv[input_at + 1 : input_at + 5] == [
        "R.pdb",
        "IM1.pdb",
        "IM2.pdb",
        "P.pdb",
    ]
    assert "-q" not in argv


def test_scan_1d_uses_one_scan_flag_for_all_stages(
    registered_tools, tmp_path: Path
) -> None:
    tools, calls = registered_tools
    first = "[(1,2,1.4)]"
    second = "[(3,4,2.0)]"

    tools["scan_1d"](
        "R.pdb",
        first,
        additional_scan_stages=[second],
        charge=0,
        out_dir=str(tmp_path / "scan"),
    )

    argv = calls[-1]["argv"]
    assert argv.count("--scan-lists") == 1
    scan_at = argv.index("--scan-lists")
    assert argv[scan_at + 1 : scan_at + 3] == [first, second]


def test_full_pipeline_uses_explicit_all_and_one_scan_flag(
    registered_tools, tmp_path: Path
) -> None:
    tools, calls = registered_tools
    first = "[(1,2,1.4)]"
    second = "[(3,4,2.0)]"

    tools["run_full_pipeline"](
        "R.pdb",
        scan_lists=first,
        additional_scan_stages=[second],
        ligand_charge="LIG:-1",
        do_tsopt=True,
        out_dir=str(tmp_path / "all"),
    )

    argv = calls[-1]["argv"]
    assert argv[:3] == ["pdb2reaction", "all", "-i"]
    assert argv.count("--scan-lists") == 1
    scan_at = argv.index("--scan-lists")
    assert argv[scan_at + 1 : scan_at + 3] == [first, second]


def test_helper_extra_args_are_appended_once(
    registered_tools,
) -> None:
    tools, calls = registered_tools

    tools["add_element_info"](
        "raw.pdb",
        "fixed.pdb",
        extra_args=["--verbose", "1"],
    )

    argv = calls[-1]["argv"]
    assert argv.count("--verbose") == 1
