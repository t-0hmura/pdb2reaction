"""Lightweight non-computational integration for ``all`` run ownership."""

from __future__ import annotations

import json
import inspect
import os
from pathlib import Path

from click.testing import CliRunner

from pdb2reaction.workflows import all as all_workflow
from pdb2reaction.core.result_commit import RUN_ID_ENV, apply_current_run_id


def test_refine_path_publishes_each_segment_mep_trajectory() -> None:
    source = inspect.getsource(all_workflow.cli.callback)
    assert "segment_mep_publications" in source
    assert 'f"seg_{segment_index:02d}"' in source
    assert '"mep_trj.xyz"' in source
    assert "Failed to publish segment" in source


def test_stopped_ts_summary_reuses_the_final_citation_payload() -> None:
    source = inspect.getsource(all_workflow.cli.callback)
    stopped_branch = source.split(
        'if not bool(_tsopt_decision.get("continue_irc")):', 1
    )[1]
    stopped_writer = stopped_branch.split("summary_payload = {", 1)[1].split(
        "write_summary_log", 1
    )[0]

    assert "citation_post_segments = [_stop_log]" in stopped_branch
    assert "summary_payload.update(_all_method_citation_payload())" in stopped_writer


def test_all_extraction_remaps_original_freeze_indices_before_path_search(
    tmp_path: Path, monkeypatch,
) -> None:
    examples = Path(__file__).parents[1] / "examples"
    reactant = examples / "1.R.pdb"
    product = examples / "3.P.pdb"
    caller_config = tmp_path / "config.yaml"
    caller_config.write_text(
        "geom:\n  freeze_atoms: [4345]\n",
        encoding="utf-8",
    )
    captured: dict[str, object] = {}

    def capture_child(name, _cli, args, **_kwargs) -> None:
        assert name == "path_search"
        config = Path(args[args.index("--config") + 1])
        captured["config"] = all_workflow.yaml.safe_load(
            config.read_text(encoding="utf-8")
        )
        captured["inputs"] = [
            Path(args[index + 1])
            for index, token in enumerate(args)
            if token == "-i"
        ]
        raise RuntimeError("stop after mapped child config capture")

    monkeypatch.setattr(all_workflow, "_run_cli_main", capture_child)
    result = CliRunner().invoke(
        all_workflow.cli,
        [
            "-i", str(reactant), "-i", str(product),
            "-b", "uma",
            "-o", str(tmp_path / "result"),
            "-c", "320,321,322",
            "-l", "GPP:-3,MG:2,SAM:1",
            "-r", "1.6",
            "--selected-resn", "186",
            "--config", str(caller_config),
            "--freeze-atoms", "4385,4399,4422",
            "--refine-path", "true",
        ],
    )

    assert result.exit_code != 0
    assert isinstance(result.exception, RuntimeError)
    assert str(result.exception) == "stop after mapped child config capture"
    child_freezes = captured["config"]["geom"]["freeze_atoms"]
    assert {47, 87, 101, 124} <= set(child_freezes)
    assert len(captured["inputs"]) == 2
    assert all(path.name.startswith("model_") for path in captured["inputs"])
    model_atom_count = sum(
        line.startswith(("ATOM", "HETATM"))
        for line in captured["inputs"][0].read_text(encoding="utf-8").splitlines()
    )
    assert all(1 <= index <= model_atom_count for index in child_freezes)


def _replace_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.next"
    temporary.write_bytes(payload)
    os.replace(temporary, path)


def _trajectory() -> bytes:
    return (
        b"2\nE=0.0 unit=hartree left\nH 0 0 0\nH 0 0 0.74\n"
        b"2\nE=0.1 unit=hartree right\nH 0 0 0\nH 0 0 0.80\n"
    )


def test_all_manifest_ignores_stale_scan_stages_and_records_distinct_runs(
    tmp_path: Path, monkeypatch,
) -> None:
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    source = tmp_path / "input.xyz"
    source.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    scan_dir = out_dir / "_work" / "scan"
    for index in (1, 2, 3):
        stale = scan_dir / f"stage_{index:02d}" / "result.xyz"
        _replace_bytes(stale, _trajectory())
    stale_hei = out_dir / "_work" / "path_opt" / "hei_seg_01.xyz"
    _replace_bytes(stale_hei, b"2\nold hei\nH 0 0 0\nH 0 0 0.74\n")
    stale_summary = out_dir / "_work" / "path_opt" / "summary.json"
    _replace_bytes(stale_summary, b'{"status": "STALE_SENTINEL"}')
    unrelated_public = out_dir / "unrelated.txt"
    _replace_bytes(unrelated_public, b"not a pipeline producer\n")

    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(
        all_workflow,
        "_write_segment_energy_diagram",
        lambda *_a, **_k: None,
    )
    monkeypatch.setattr(
        all_workflow._path_search,
        "has_bond_change",
        lambda *_a, **_k: (False, ""),
    )

    emit_hei = {"value": True}
    child_run_ids: list[tuple[str | None, str | None]] = []

    def fake_child(name, _cli, args, **_kwargs) -> None:
        child_run_ids.append(
            (
                os.environ.get(RUN_ID_ENV),
                apply_current_run_id({}).get("run_id"),
            )
        )
        child_out = Path(args[args.index("--out-dir") + 1])
        if name == "scan":
            _replace_bytes(child_out / "stage_01" / "result.xyz", _trajectory())
            _replace_bytes(
                child_out / "result.json",
                json.dumps(
                    apply_current_run_id(
                        {
                            "scientific_status": "success",
                            "preopt_converged": None,
                        }
                    )
                ).encode("utf-8"),
            )
        elif name == "path-opt":
            _replace_bytes(child_out / "final_geometries_trj.xyz", _trajectory())
            _replace_bytes(
                child_out / "result.json",
                json.dumps(
                    apply_current_run_id({
                        "stage_outcomes": [{"converged": True}],
                    })
                ).encode("utf-8"),
            )
            if emit_hei["value"]:
                _replace_bytes(
                    child_out / "hei.xyz",
                    b"2\nE=0.1 unit=hartree hei\nH 0 0 0\nH 0 0 0.77\n",
                )
        else:  # pragma: no cover - this test requests no post-processing
            raise AssertionError(name)

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_child)

    run_ids = []
    for _ in range(2):
        child_start = len(child_run_ids)
        result = CliRunner().invoke(
            all_workflow.cli,
            [
                "-i",
                str(source),
                "-q",
                "0",
                "--scan-lists",
                "[(1,2,0.80)]",
                "--no-preopt",
                "--out-dir",
                str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        internal = json.loads(
            (out_dir / "_work" / "_run_manifest.json").read_text(encoding="utf-8")
        )
        produced = internal["produced"]
        assert "scan.stage.01" in produced
        assert "scan.stage.02" not in produced
        assert "scan.stage.03" not in produced
        assert "path.hei.01" in produced
        assert produced["path.hei.01"]["stamp"]["sha256"]
        public_summary = json.loads(
            (out_dir / "summary.json").read_text(encoding="utf-8")
        )
        assert public_summary["run_id"] == internal["run_id"]
        assert "STALE_SENTINEL" not in json.dumps(public_summary)
        assert "summary.log" in public_summary["key_output_files"]
        assert "unrelated.txt" not in public_summary["key_output_files"]
        assert "output.public.summary.log" in produced
        assert "output.public.unrelated.txt" not in produced
        assert child_run_ids[child_start:]
        assert all(
            environment_id == internal["run_id"]
            and payload_id == internal["run_id"]
            for environment_id, payload_id in child_run_ids[child_start:]
        )
        assert RUN_ID_ENV not in os.environ
        run_ids.append(internal["run_id"])

    assert run_ids[0] != run_ids[1]

    emit_hei["value"] = False
    stale_hei_run = CliRunner().invoke(
        all_workflow.cli,
        [
            "-i",
            str(source),
            "-q",
            "0",
            "--scan-lists",
            "[(1,2,0.80)]",
            "--no-preopt",
            "--out-dir",
            str(out_dir),
        ],
    )
    assert stale_hei_run.exit_code == 0, stale_hei_run.output
    stale_manifest = json.loads(
        (out_dir / "_work" / "_run_manifest.json").read_text(encoding="utf-8")
    )
    assert "path.hei.01" not in stale_manifest["produced"]


def test_path_search_cannot_consume_an_unchanged_prior_summary(
    tmp_path: Path, monkeypatch,
) -> None:
    source = tmp_path / "input.xyz"
    source.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    stale_summary = out_dir / "_work" / "path_search" / "summary.json"
    _replace_bytes(
        stale_summary,
        b'{"status": "STALE_SENTINEL", "segments": []}',
    )

    def fake_child(name, _cli, args, **_kwargs) -> None:
        child_out = Path(args[args.index("--out-dir") + 1])
        if name == "scan":
            _replace_bytes(child_out / "stage_01" / "result.xyz", _trajectory())
            _replace_bytes(
                child_out / "result.json",
                json.dumps(
                    apply_current_run_id(
                        {
                            "scientific_status": "success",
                            "preopt_converged": None,
                        }
                    )
                ).encode("utf-8"),
            )
        elif name == "path_search":
            return
        else:  # pragma: no cover
            raise AssertionError(name)

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_child)

    result = CliRunner().invoke(
        all_workflow.cli,
        [
            "-i",
            str(source),
            "-q",
            "0",
            "--scan-lists",
            "[(1,2,0.80)]",
            "--no-preopt",
            "--refine-path",
            "true",
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code != 0
    assert "did not produce 'path.summary'" in result.output
    assert not (out_dir / "summary.json").exists()


def test_all_stops_before_path_search_when_scan_outcome_failed(
    tmp_path: Path, monkeypatch,
) -> None:
    source = tmp_path / "input.xyz"
    source.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    dispatched = []

    def fake_child(name, _cli, args, **_kwargs) -> None:
        dispatched.append(name)
        child_out = Path(args[args.index("--out-dir") + 1])
        if name != "scan":  # pragma: no cover
            raise AssertionError(name)
        _replace_bytes(child_out / "stage_01" / "result.xyz", _trajectory())
        _replace_bytes(
            child_out / "result.json",
            json.dumps(
                apply_current_run_id(
                    {
                        "scientific_status": "failed",
                        "scientific_status_reasons": ["stage_1: not converged"],
                        "preopt_converged": None,
                    }
                )
            ).encode("utf-8"),
        )

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_child)

    result = CliRunner().invoke(
        all_workflow.cli,
        [
            "-i",
            str(source),
            "-q",
            "0",
            "--scan-lists",
            "[(1,2,0.80)]",
            "--no-preopt",
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code != 0
    assert "did not produce a scientifically usable path" in result.output
    assert dispatched == ["scan"]
