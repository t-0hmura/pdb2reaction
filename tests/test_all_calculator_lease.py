"""Current-run TS/IRC artifacts and parent calculator lifetime."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pdb2reaction.workflows import all as all_workflow
from pdb2reaction.workflows._run_session import (
    ArtifactClaimError,
    InvocationManifest,
    RunSession,
)


class _Prepared:
    def __init__(self, path: Path) -> None:
        self.source_path = Path(path)
        self.is_gjf = False
        self.cleanup_count = 0

    def cleanup(self) -> None:
        self.cleanup_count += 1


class _Geometry:
    def __init__(self) -> None:
        self.calculator = None

    def set_calculator(self, calculator) -> None:
        self.calculator = calculator


class _Calculator:
    def __init__(self, counters: dict[str, int]) -> None:
        self.counters = counters
        self.closed = False
        counters["live"] += 1
        counters["max_live"] = max(counters["max_live"], counters["live"])

    def close(self) -> None:
        if not self.closed:
            self.closed = True
            self.counters["close"] += 1
            self.counters["live"] -= 1


def _xyz_frames() -> str:
    return (
        "1\nleft\nH 0.0 0.0 0.0\n"
        "1\nright\nH 0.0 0.0 0.1\n"
    )


def _install_lightweight_pipeline(monkeypatch, events, calculators, geometries):
    monkeypatch.setattr(
        all_workflow,
        "prepare_input_structure",
        lambda path: _Prepared(Path(path)),
    )
    monkeypatch.setattr(all_workflow, "apply_ref_pdb_override", lambda *_a: None)
    monkeypatch.setattr(all_workflow, "convert_xyz_like_outputs", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow, "_get_freeze_atoms", lambda *_a: [])
    monkeypatch.setattr(all_workflow, "_freeze_atoms_for_log", lambda: [])
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)

    def new_geometry(*_args, **_kwargs):
        geometry = _Geometry()
        geometries.append(geometry)
        return geometry

    monkeypatch.setattr(all_workflow, "geom_loader", new_geometry)
    monkeypatch.setattr(all_workflow, "_geom_from_angstrom", new_geometry)

    def create_calculator(**kwargs):
        events.append(("factory", dict(kwargs)))
        calculator = _Calculator(calculators)
        return calculator

    monkeypatch.setattr(all_workflow, "create_calculator", create_calculator)

    def run_child(name, _cli, args, **_kwargs):
        events.append(("child", name))
        out_dir = Path(args[args.index("--out-dir") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        if name == "tsopt":
            (out_dir / "result.json").write_text(
                json.dumps({"status": "converged", "n_imaginary_modes": 1}),
                encoding="utf-8",
            )
            (out_dir / "final_geometry.xyz").write_text(
                "1\nts\nH 0 0 0\n", encoding="utf-8"
            )
        elif name == "irc":
            (out_dir / "finished_irc_trj.xyz").write_text(
                _xyz_frames(), encoding="utf-8"
            )

    monkeypatch.setattr(all_workflow, "_run_cli_main", run_child)


def test_ts_then_irc_constructs_one_complete_shared_calculator_after_child(
    tmp_path: Path, monkeypatch,
) -> None:
    events = []
    counters = {"live": 0, "max_live": 0, "close": 0}
    geometries = []
    _install_lightweight_pipeline(monkeypatch, events, counters, geometries)
    hei = tmp_path / "hei.xyz"
    hei.write_text("1\nhei\nH 0 0 0\n", encoding="utf-8")
    manifest = InvocationManifest()
    session = RunSession(manifest=manifest)
    calc_cfg = {
        "backend": "custom",
        "calc_file": "factory.py",
        "calc_factory": "build",
        "workers": 3,
        "solvent": "water",
        "precision": "fp64",
    }
    segment_root = tmp_path / "segments" / "seg_01"

    ts_path, g_ts = all_workflow._run_tsopt_on_hei(
        hei, 0, 1, calc_cfg, None, segment_root, False, "hess", None, False,
        manifest=manifest,
        artifact_prefix="post.01.tsopt",
        public_root=tmp_path,
    )
    assert ts_path.name == "final_geometry.xyz"
    assert [event for event in events if event[0] == "factory"] == []

    result = all_workflow._irc_and_match(
        1,
        segment_root,
        ts_path,
        hei,
        None,
        g_ts,
        0,
        1,
        False,
        calc_cfg,
        None,
        False,
        session=session,
        manifest=manifest,
        artifact_prefix="post.01.irc",
        public_root=tmp_path,
    )

    assert events.index(("child", "irc")) < next(
        index for index, event in enumerate(events) if event[0] == "factory"
    )
    assert [event for event in events if event[0] == "factory"] == [
        ("factory", calc_cfg)
    ]
    assert manifest.paths("output.public.") == [
        segment_root / "ts" / "result.json",
        segment_root / "ts" / "final_geometry.xyz",
        segment_root / "irc" / "finished_irc_trj.xyz",
    ]
    calculator = result["ts_geom"].calculator
    assert calculator is result["left_min_geom"].calculator
    assert calculator is result["right_min_geom"].calculator
    result["calculator_lease"].release()
    session.close()
    assert all(geometry.calculator is None for geometry in geometries)
    assert counters == {"live": 0, "max_live": 1, "close": 1}


def test_irc_exception_releases_and_closes_lease_once(tmp_path: Path, monkeypatch) -> None:
    events = []
    counters = {"live": 0, "max_live": 0, "close": 0}
    geometries = []
    _install_lightweight_pipeline(monkeypatch, events, counters, geometries)
    monkeypatch.setattr(
        all_workflow,
        "_orient_irc_endpoints",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("mapping failed")),
    )
    source = tmp_path / "ts.xyz"
    source.write_text("1\nts\nH 0 0 0\n", encoding="utf-8")
    g_ts = _Geometry()
    geometries.append(g_ts)
    manifest = InvocationManifest()
    session = RunSession(manifest=manifest)

    with pytest.raises(RuntimeError, match="mapping failed"):
        all_workflow._irc_and_match(
            1,
            tmp_path / "seg",
            source,
            source,
            None,
            g_ts,
            0,
            1,
            False,
            {"backend": "uma"},
            None,
            False,
            session=session,
            manifest=manifest,
            artifact_prefix="post.01.irc",
        )
    session.close()

    assert all(geometry.calculator is None for geometry in geometries)
    assert counters == {"live": 0, "max_live": 1, "close": 1}


def test_two_irc_segments_never_hold_two_parent_calculators(
    tmp_path: Path, monkeypatch,
) -> None:
    events = []
    counters = {"live": 0, "max_live": 0, "close": 0}
    geometries = []
    _install_lightweight_pipeline(monkeypatch, events, counters, geometries)
    source = tmp_path / "ts.xyz"
    source.write_text("1\nts\nH 0 0 0\n", encoding="utf-8")
    manifest = InvocationManifest()
    session = RunSession(manifest=manifest)

    for index in (1, 2):
        result = all_workflow._irc_and_match(
            index,
            tmp_path / f"seg_{index}",
            source,
            source,
            None,
            _Geometry(),
            0,
            1,
            False,
            {"backend": "uma"},
            None,
            False,
            session=session,
            manifest=manifest,
            artifact_prefix=f"post.{index:02d}.irc",
        )
        result["calculator_lease"].release()
    session.close()

    assert counters == {"live": 0, "max_live": 1, "close": 2}


def test_tsopt_rejects_unchanged_prior_result(tmp_path: Path, monkeypatch) -> None:
    events = []
    counters = {"live": 0, "max_live": 0, "close": 0}
    geometries = []
    _install_lightweight_pipeline(monkeypatch, events, counters, geometries)
    hei = tmp_path / "hei.xyz"
    hei.write_text("1\nhei\nH 0 0 0\n", encoding="utf-8")
    ts_dir = tmp_path / "seg" / "ts"
    ts_dir.mkdir(parents=True)
    (ts_dir / "result.json").write_text(
        json.dumps({"status": "converged", "n_imaginary_modes": 1}),
        encoding="utf-8",
    )
    (ts_dir / "final_geometry.xyz").write_text(
        "1\nold\nH 0 0 0\n", encoding="utf-8"
    )
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *_a, **_k: None)

    with pytest.raises(ArtifactClaimError, match="did not produce"):
        all_workflow._run_tsopt_on_hei(
            hei,
            0,
            1,
            {"backend": "uma"},
            None,
            tmp_path / "seg",
            False,
            "hess",
            None,
            False,
            manifest=InvocationManifest(),
            artifact_prefix="post.01.tsopt",
        )


def test_tsopt_never_converts_an_unclaimed_stale_xyz(
    tmp_path: Path, monkeypatch,
) -> None:
    events = []
    counters = {"live": 0, "max_live": 0, "close": 0}
    geometries = []
    _install_lightweight_pipeline(monkeypatch, events, counters, geometries)
    hei = tmp_path / "hei.xyz"
    hei.write_text("1\nhei\nH 0 0 0\n", encoding="utf-8")
    ts_dir = tmp_path / "seg" / "ts"
    ts_dir.mkdir(parents=True)
    (ts_dir / "final_geometry.xyz").write_text(
        "1\nstale\nH 9 9 9\n", encoding="utf-8"
    )
    conversion_sources: list[Path] = []
    monkeypatch.setattr(
        all_workflow,
        "convert_xyz_like_outputs",
        lambda source, *_a, **_k: conversion_sources.append(Path(source)),
    )

    def write_current_pdb(name, _cli, args, **_kwargs):
        assert name == "tsopt"
        out_dir = Path(args[args.index("--out-dir") + 1])
        (out_dir / "result.json").write_text(
            json.dumps({"status": "converged", "n_imaginary_modes": 1}),
            encoding="utf-8",
        )
        (out_dir / "final_geometry.pdb").write_text("CURRENT\n", encoding="utf-8")

    monkeypatch.setattr(all_workflow, "_run_cli_main", write_current_pdb)

    ts_path, _ = all_workflow._run_tsopt_on_hei(
        hei,
        0,
        1,
        {"backend": "uma"},
        None,
        tmp_path / "seg",
        False,
        "hess",
        None,
        False,
        manifest=InvocationManifest(),
        artifact_prefix="post.01.tsopt",
    )

    assert ts_path.name == "final_geometry.pdb"
    assert conversion_sources == []
