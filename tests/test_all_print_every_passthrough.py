"""``all --print-every`` reaches optimizing children as a CLI flag, not via YAML.

path-opt and path-search accept a hidden ``--print-every`` that only sets the
single-structure optimizers (LBFGS/RFO) and wins over YAML without a conflict.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.result_commit import RUN_ID_ENV
from pdb2reaction.workflows import all as all_workflow

SMOKE = Path(__file__).resolve().parent / "smoke"


class ReachedChild(BaseException):
    """Escape the parent pipeline at the first child dispatch."""


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _pdb(index: int) -> str:
    coords = [(index * .125, .25, -.5), (1.25 + index * .25, .5, .75)]
    return "".join(
        f"HETATM{serial:5d} {name:<4s} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n"
        for serial, (name, element, (x, y, z)) in enumerate(
            zip(("C1", "O1"), ("C", "O"), coords), 1
        )
    ) + "END\n"


_CASES = {
    "path-opt": (2, []),
    "path_search": (2, ["--refine-path"]),
    "scan": (1, ["--scan-lists", "[(1,2,1.50)]"]),
    "tsopt": (1, ["--tsopt"]),
}


@pytest.mark.parametrize("print_every", [7, None], ids=["explicit", "omitted"])
@pytest.mark.parametrize("child_name", list(_CASES))
def test_all_forwards_print_every_as_child_cli_flag(
    tmp_path, monkeypatch, child_name, print_every
):
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    n_inputs, extra = _CASES[child_name]
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in range(n_inputs)]
    config = _write(tmp_path / "config.yaml", "lbfgs:\n  print_every: 50\n")
    captured = {}

    def child(name, _cli, args, **kwargs):
        captured["name"] = name
        captured["args"] = list(args)
        if "--config" in args:
            captured["config"] = yaml.safe_load(
                Path(args[args.index("--config") + 1]).read_text(encoding="utf-8")
            )
        raise ReachedChild

    def no_calculator(*args, **kwargs):
        raise AssertionError("No calculator may be built before the child dispatch")

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "create_calculator", no_calculator)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
             "--no-preopt", "--no-freeze-links", "--config", str(config), *extra]
    if print_every is not None:
        args += ["--print-every", str(print_every)]

    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    assert captured["name"] == child_name
    child_args = captured["args"]
    if print_every is None:
        assert "--print-every" not in child_args
    else:
        assert child_args.count("--print-every") == 1
        assert child_args[child_args.index("--print-every") + 1] == str(print_every)
    forwarded = captured.get("config") or {}
    assert "print_every" not in (forwarded.get("opt") or {})
    if forwarded:
        assert forwarded["lbfgs"]["print_every"] == 50


def _path_cfg(tmp_path: Path, command: str, yaml_text: str, *extra: str):
    config = _write(tmp_path / "config.yaml", "calc:\n  charge: -1\n  spin: 1\n" + yaml_text)
    structure = SMOKE / "r.pdb"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(structure), str(structure),
            "--config", str(config), "--dry-run", "--show-config", "-v", "3",
            "--out-dir", str(tmp_path / command), *extra,
        ],
    )
    return result


def _block(output: str, title: str) -> dict:
    head = f"\n{title}\n{'-' * len(title)}\n"
    start = output.index(head) + len(head)
    end = output.find("\n\n", start)
    return yaml.safe_load(output[start:] if end < 0 else output[start:end])


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
@pytest.mark.parametrize(
    ("opt_mode", "kind"), [("grad", "lbfgs"), ("hess", "rfo")]
)
@pytest.mark.parametrize(
    "yaml_text",
    ["", "opt:\n  print_every: 50\n", "lbfgs:\n  print_every: 50\nrfo:\n  print_every: 50\n"],
    ids=["no-yaml", "opt-yaml", "optimizer-yaml"],
)
def test_path_hidden_print_every_sets_single_structure_optimizers(
    tmp_path, command, opt_mode, kind, yaml_text
):
    result = _path_cfg(
        tmp_path, command, yaml_text, "--opt-mode", opt_mode, "--print-every", "7"
    )

    assert result.exit_code == 0, result.output
    assert _block(result.output, f"opt.{kind}")["print_every"] == 7
    assert _block(result.output, "stopt")["print_every"] == 10


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_path_print_every_without_cli_keeps_yaml_value(tmp_path, command):
    result = _path_cfg(tmp_path, command, "lbfgs:\n  print_every: 50\n")

    assert result.exit_code == 0, result.output
    assert _block(result.output, "opt.lbfgs")["print_every"] == 50


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_path_print_every_is_hidden_from_help(command):
    for flag in ("--help", "--help-advanced"):
        result = CliRunner().invoke(root_cli, [command, flag])
        assert result.exit_code == 0, result.output
        assert "--print-every" not in result.output
