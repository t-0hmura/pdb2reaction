"""Real parent CLI/publication regression; no computational child is run."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from pdb2reaction.core.result_commit import RUN_ID_ENV, apply_current_run_id
from pdb2reaction.workflows import all as all_workflow


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _coords(index: int):
    return [(index * 0.125, 0.25, -0.5), (1.25 + index * 0.25, 0.5, 0.75)]


def _frame(index: int) -> str:
    atoms = "".join(
        f"{element} {x:.3f} {y:.3f} {z:.3f}\n"
        for element, (x, y, z) in zip(("C", "O"), _coords(index))
    )
    return f"2\nE={index * 0.01:.3f} unit=hartree frame{index}\n{atoms}"


def _pdb(residue: str) -> str:
    # Reference coordinates intentionally differ from the child trajectory.
    return "".join(
        f"HETATM{serial:5d} {name:<4s} {residue:3s} A   1    "
        f"{x:8.3f}{8.0:8.3f}{9.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        for serial, name, element, x in ((1, "C1", "C", 8.0), (2, "O1", "O", 9.0))
    ) + "END\n"


def _assert_pdb(path: Path, frames, residue: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    assert sum(line.startswith("MODEL") for line in lines) == len(frames)
    assert sum(line.startswith("ENDMDL") for line in lines) == len(frames)
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    assert len(atoms) == 2 * len(frames)
    assert [line[12:16].strip() for line in atoms] == ["C1", "O1"] * len(frames)
    assert [line[76:78].strip() for line in atoms] == ["C", "O"] * len(frames)
    assert {line[17:20] for line in atoms} == {residue}
    actual = [tuple(float(line[start:start + 8]) for start in (30, 38, 46)) for line in atoms]
    assert actual == [xyz for index in frames for xyz in _coords(index)]


@pytest.mark.parametrize("case", ["scan_preopt", "direct_pdb", "reference_free_xyz", "no_convert"])
def test_all_path_opt_mep_pdb_publication(tmp_path: Path, monkeypatch, case: str) -> None:
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    out = tmp_path / "out"
    has_reference = case != "reference_free_xyz"
    convert = case != "no_convert"
    scan_preopt = case in {"scan_preopt", "no_convert"}
    references = []
    child_calls = []
    if case == "direct_pdb":
        inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(label)) for i, label in enumerate(("RAW", "MID", "END"))]
        extra = ["--no-preopt"]
    else:
        inputs = [_write(tmp_path / "input.xyz", _frame(0))]
        extra = ["--scan-lists", "[(1,2,1.5)]", "--scan-lists", "[(1,2,1.75)]"]
        extra += ["--preopt", "true"] if scan_preopt else ["--no-preopt"]

    def fake_child(name, _cli, args, **_kwargs):
        child_calls.append(name)
        child_out = Path(args[args.index("--out-dir") + 1])
        if name == "scan":
            if scan_preopt:
                _write(child_out / "preopt/result.xyz", _frame(0))
                _write(child_out / "preopt/result.pdb", _pdb("PRE"))
            for index in (1, 2):
                _write(child_out / f"stage_{index:02d}/result.xyz", _frame(index))
                if has_reference:
                    _write(child_out / f"stage_{index:02d}/result.pdb", _pdb(f"S{index:02d}"))
            payload = {"scientific_status": "success", "preopt_converged": True if scan_preopt else None, "stages": []}
        elif name == "path-opt":
            pair = len(references)
            refs = [Path(args[i + 1]) for i, token in enumerate(args) if token == "--ref-pdb"]
            references.append(refs)
            left = Path(args[args.index("-i") + 1])
            assert left.suffix == (".pdb" if case == "direct_pdb" else ".xyz")
            _write(child_out / "final_geometries_trj.xyz", _frame(pair) + _frame(pair + 1))
            payload = {"stage_outcomes": [{"converged": True}], "preopt_converged": True}
        else:
            raise AssertionError(f"Unexpected computational child: {name}")
        _write(child_out / "result.json", json.dumps(apply_current_run_id(payload)))

    def no_calculator(*_args, **_kwargs):
        raise AssertionError("This publication test must not construct a calculator")

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_child)
    monkeypatch.setattr(all_workflow, "create_calculator", no_calculator)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow._path_search, "has_bond_change", lambda *_a, **_k: (False, ""))
    args = [arg for path in inputs for arg in ("-i", str(path))]
    args += ["-q", "0", "-m", "1", "--out-dir", str(out)]
    args += ["--convert-files", "true"] if convert else ["--no-convert-files"]
    args += extra
    result = CliRunner().invoke(all_workflow.cli, args)
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"
    assert child_calls == (["scan"] if case != "direct_pdb" else []) + ["path-opt", "path-opt"]
    summary = json.loads((out / "summary.json").read_text())
    manifest = json.loads((out / "_work/_run_manifest.json").read_text())
    assert summary["run_id"] == manifest["run_id"]
    assert summary["n_images"] == 3 and summary["n_segments"] == 2
    assert (out / "mep_trj.xyz").read_text() == "".join(_frame(i) for i in range(3))
    assert not (out / "mep_trj.pdb").exists()  # Current release filename is mep.pdb.
    expected_pdb = convert and has_reference
    assert (out / "mep.pdb").exists() is expected_pdb
    assert ("mep.pdb" in summary["key_output_files"]) is expected_pdb
    assert ("output.public.mep.pdb" in manifest["produced"]) is expected_pdb
    path_dir = out / "_work/path_opt"
    if has_reference:
        labels = ["RAW", "MID"] if case == "direct_pdb" else ["PRE", "S01"]
        assert [len(refs) for refs in references] == [1, 1]
        for refs, label in zip(references, labels):
            assert next(line[17:20] for line in refs[0].read_text().splitlines() if line.startswith("HETATM")) == label
    else:
        assert references == [[], []]
    for index in (1, 2):
        segment_pdb = path_dir / f"mep_seg_{index:02d}.pdb"
        assert segment_pdb.exists() is expected_pdb
        if expected_pdb:
            _assert_pdb(segment_pdb, [index - 1, index], labels[index - 1])
    if expected_pdb:
        _assert_pdb(out / "mep.pdb", [0, 1, 2], labels[0])
