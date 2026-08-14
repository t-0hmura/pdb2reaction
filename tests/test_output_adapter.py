"""Direct/library output must not depend on CLI echo bootstrap."""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import click

from pdb2reaction.core.output import _TAG_AWARE_MARKER, emit


PRIVATE_TAGS = {"narrative", "detail", "force", "raw_path"}


def test_fresh_process_direct_output_adapter_and_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "profile.csv"
    code = f"""
from pathlib import Path
from pdb2reaction.core.utils import emit, echo_run_summary, emit_optimizer_terminal_status
from pdb2reaction.io.trj2fig import write_csv

emit("narrative", narrative=True)
emit("detail", narrative=False, detail=True)
emit("forced", narrative=False, force=True)
emit("raw", narrative=False, raw_path=True)
echo_run_summary({{"backend": "custom"}})
emit_optimizer_terminal_status("opt", converged=True, cycles=2, max_cycles=4)
write_csv(Path({str(csv_path)!r}), [-1.0, -0.9], [0.0, 62.75], "kcal", True)
"""

    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[1],
    )

    assert proc.returncode == 0, proc.stderr
    assert "narrative" in proc.stdout
    assert "[opt] Converged!" in proc.stdout
    assert csv_path.exists()
    assert "delta_kcal" in csv_path.read_text(encoding="utf-8")


def test_only_output_adapter_passes_private_tags_to_native_click_echo() -> None:
    package_root = Path(__file__).resolve().parents[1] / "pdb2reaction"
    violations: list[str] = []

    class Visitor(ast.NodeVisitor):
        def __init__(self, path: Path) -> None:
            self.path = path
            self.functions: list[str] = []

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            self.functions.append(node.name)
            self.generic_visit(node)
            self.functions.pop()

        visit_AsyncFunctionDef = visit_FunctionDef

        def visit_Call(self, node: ast.Call) -> None:
            tagged = PRIVATE_TAGS & {kw.arg for kw in node.keywords if kw.arg}
            native_echo = (
                isinstance(node.func, ast.Attribute)
                and node.func.attr == "echo"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id in {"click", "_click"}
            )
            allowed_adapter_call = (
                self.path.name == "output.py"
                and self.path.parent.name == "core"
                and self.functions
                and self.functions[-1] == "emit"
            )
            if tagged and native_echo and not allowed_adapter_call:
                violations.append(
                    f"{self.path.relative_to(package_root)}:{node.lineno}"
                )
            self.generic_visit(node)

    for path in sorted(package_root.rglob("*.py")):
        Visitor(path).visit(ast.parse(path.read_text(encoding="utf-8")))

    assert violations == []


def test_adapter_rechecks_current_echo_after_post_bootstrap_monkeypatch(
    monkeypatch,
) -> None:
    tagged_calls: list[dict] = []
    native_calls: list[dict] = []

    def tagged_echo(_message=None, **kwargs):
        tagged_calls.append(kwargs)

    setattr(tagged_echo, _TAG_AWARE_MARKER, True)
    monkeypatch.setattr(click, "echo", tagged_echo)
    emit("first", narrative=False, detail=True, force=True, raw_path=True)
    assert tagged_calls == [
        {"narrative": False, "detail": True, "force": True, "raw_path": True}
    ]

    def replacement_native_echo(_message=None, **kwargs):
        native_calls.append(kwargs)

    monkeypatch.setattr(click, "echo", replacement_native_echo)
    emit(
        "second",
        narrative=False,
        detail=True,
        force=True,
        raw_path=True,
        err=True,
    )
    assert native_calls == [{"err": True}]


def test_hessian_status_does_not_add_blank_lines() -> None:
    code = """
from pdb2reaction.core.utils import (
    _patch_click_echo,
    emit,
    set_console_gating,
    set_pipeline_mode,
    set_verbose_level,
)
set_verbose_level(2)
set_console_gating(True)
set_pipeline_mode(True)
_patch_click_echo()
print('cycle 400')
emit('[hessian] Completed FiniteDifference Hessian: 15.44 s', detail=True)
print('separator')
print('cycle 500')
"""
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[1],
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == [
        "cycle 400",
        "[hessian] Completed FiniteDifference Hessian: 15.44 s",
        "separator",
        "cycle 500",
    ]
