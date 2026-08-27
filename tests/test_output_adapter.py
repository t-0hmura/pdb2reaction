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


def test_stdout_spacing_respects_an_intervening_stderr_line() -> None:
    code = """
import click
from pdb2reaction.core.output import emit
from pdb2reaction.core.utils import _patch_click_echo

_patch_click_echo()
emit('[mode]')
emit('')
click.echo('[warning]', err=True)
emit('\\n[backend] Preparing MLIP model (uma / uma-s-1p2)...')
"""
    proc = subprocess.run(
        [sys.executable, "-u", "-c", code],
        cwd=Path(__file__).resolve().parents[1],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    assert proc.returncode == 0, proc.stdout
    assert proc.stdout.splitlines() == [
        "[mode]",
        "",
        "[warning]",
        "",
        "[backend] Preparing MLIP model (uma / uma-s-1p2)...",
    ]


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
print('Convergence thresholds:')
emit('[hessian] Completed FiniteDifference Hessian: 15.30 s', detail=True)
emit('[Imaginary modes] n=1', detail=True)
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
        "Convergence thresholds:",
        "[hessian] Completed FiniteDifference Hessian: 15.30 s",
        "[Imaginary modes] n=1",
        "cycle 500",
    ]


def test_mep_summary_body_is_visible_at_default_pipeline_verbosity() -> None:
    from pdb2reaction.workflows import path_search

    source = Path(path_search.__file__).read_text(encoding="utf-8")
    lines = source.splitlines()
    start = next(
        i for i, line in enumerate(lines, 1) if "MEP summary started" in line
    )
    finish = next(
        i
        for i, line in enumerate(lines[start:], start + 1)
        if "mep summary finished" in line.lower()
    )
    tree = ast.parse(source)
    output_calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not start <= node.lineno <= finish:
            continue
        is_emit = isinstance(node.func, ast.Name) and node.func.id == "emit"
        is_click_echo = (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "echo"
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "click"
        )
        if is_emit or is_click_echo:
            output_calls.append((node, is_emit))

    assert len(output_calls) >= 8
    for node, is_emit in output_calls:
        keywords = {kw.arg: kw.value for kw in node.keywords if kw.arg}
        if is_emit:
            assert isinstance(keywords.get("narrative"), ast.Constant)
            assert keywords["narrative"].value is True
        else:
            assert isinstance(keywords.get("err"), ast.Constant)
            assert keywords["err"].value is True
