"""``all --precision-post`` sets the precision of the stages after the MEP.

TSOPT, IRC, the post-IRC endpoint optimizations (their calculator) and freq
use ``--precision-post`` when it is given and ``--precision`` when it is not;
scan, path-opt and path-search always use ``--precision``. Children are stubbed
at the dispatcher and the parent calculator is a recording stub, so no backend
runs.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.result_commit import RUN_ID_ENV, apply_current_run_id
from pdb2reaction.workflows import all as all_workflow


class ReachedChild(BaseException):
    """Escape the parent pipeline at a child dispatch."""


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _coords(index):
    return [(index * .125, .25, -.5), (1.25 + index * .25, .5, .75)]


def _frame(index):
    atoms = "".join(
        f"{element} {x:.3f} {y:.3f} {z:.3f}\n"
        for element, (x, y, z) in zip(("C", "O"), _coords(index))
    )
    return f"2\nE={(0., .1, .02)[index]:.3f} unit=hartree frame{index}\n{atoms}"


def _pdb(index):
    return "".join(
        f"HETATM{serial:5d} {name:<4s} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n"
        for serial, (name, element, (x, y, z)) in enumerate(
            zip(("C1", "O1"), ("C", "O"), _coords(index)), 1
        )
    ) + "END\n"


def _child_precision(args):
    """Precision a child resolves: its --precision, else calc.precision of --config."""
    if "--precision" in args:
        return args[args.index("--precision") + 1]
    config = yaml.safe_load(
        Path(args[args.index("--config") + 1]).read_text(encoding="utf-8")
    )
    return config["calc"]["precision"]


_TS_RESULT = {
    "optimization_status": "converged", "converged": True,
    "execution_status": "completed", "scientific_status": "success",
    "hessian_status": "completed", "saddle_validation": "first_order",
    "n_imaginary_modes": 1, "imaginary_frequencies_cm": [-500.0],
    "reaction_mode_index": 0, "reaction_mode_frequency_cm": -500.0,
    "energy_hartree": 0.1,
}


class _EnergyOnly:
    """Parent-side calculator stub: a finite energy, nothing else."""

    def get_energy(self, atoms, coords, **kwargs):
        return {"energy": -1.0}


def _run(tmp_path, monkeypatch, n_inputs, extra, precision_post):
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)[:n_inputs]]
    calls, calculators, identities = [], [], []

    def child(name, _cli, args, **kwargs):
        args = list(args)
        # Read now: the child YAML is removed when the run unwinds.
        calls.append((name, args, _child_precision(args)))
        out = Path(args[args.index("--out-dir") + 1])
        if name == "path-opt":
            _write(out / "final_geometries_trj.xyz", "".join(_frame(i) for i in range(3)))
            _write(out / "hei.xyz", _frame(1))
            _write(out / "hei.pdb", _pdb(1))
            _write(out / "result.json", json.dumps(apply_current_run_id({
                "optimization_status": "converged", "converged": True,
                "preopt_requested": False, "preopt_converged": None,
                "path_optimizers": [],
                "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                    "executed": True, "converged": True,
                                    "usable": True, "required": True}],
            })))
        elif name == "tsopt":
            _write(out / "final_geometry.xyz", _frame(1))
            _write(out / "result.json", json.dumps(apply_current_run_id(_TS_RESULT)))
        elif name == "irc":
            trajectory = _write(out / "finished_irc_trj.xyz", _frame(0) + _frame(2))
            _write(out / "result.json", json.dumps(apply_current_run_id(
                {"files": {"finished_irc": str(trajectory)}}
            )))
        else:
            raise ReachedChild
        return 0

    def calculator(**cfg):
        calculators.append(dict(cfg))
        return _EnergyOnly()

    def endpoint_opt(geom, _mode, out_dir, tag, *, calc_identity_cfg=None, outcome=None, **kwargs):
        identities.append(dict(calc_identity_cfg))
        if outcome is not None:
            outcome.update({"optimization_status": "converged", "converged": True})
        return geom, Path(out_dir) / f"{tag}.xyz", True

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "create_calculator", calculator)
    monkeypatch.setattr(all_workflow, "_optimize_endpoint_geom", endpoint_opt)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(
        all_workflow, "_write_segment_energy_diagram",
        lambda prefix, *, labels, energies_au, **kwargs: {
            "name": prefix.stem, "labels": labels, "energies_au": list(energies_au),
        },
    )
    monkeypatch.setattr(
        all_workflow._path_search, "has_bond_change",
        lambda *args, **kwargs: (True, "Bond formed"),
    )
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
             "--no-preopt", "--no-freeze-links", "--no-tsopt-from-mep-tan",
             "--precision", "tf32", *extra]
    if precision_post is not None:
        args += ["--precision-post", precision_post]
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)
    return calls, calculators, identities


_PRE_MEP = {
    "path-opt": (2, []),
    "path_search": (2, ["--refine-path"]),
    "scan": (1, ["--scan-lists", "[(1,2,1.50)]"]),
}


@pytest.mark.parametrize("precision_post", [None, "fp64"], ids=["post-unset", "post-fp64"])
@pytest.mark.parametrize("child_name", list(_PRE_MEP))
def test_pre_mep_children_keep_precision(tmp_path, monkeypatch, child_name, precision_post):
    n_inputs, extra = _PRE_MEP[child_name]

    def first_child(name, _cli, args, **kwargs):
        captured.append((name, list(args), _child_precision(list(args))))
        raise ReachedChild

    captured = []
    monkeypatch.setattr(all_workflow, "_run_cli_main", first_child)
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)[:n_inputs]]
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
             "--no-preopt", "--no-freeze-links", "--precision", "tf32", *extra]
    if precision_post is not None:
        args += ["--precision-post", precision_post]
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    (name, child_args, resolved), = captured
    assert name == child_name
    assert "--precision" not in child_args
    assert resolved == "tf32"


@pytest.mark.parametrize("precision_post", [None, "fp64"], ids=["post-unset", "post-fp64"])
def test_tsopt_only_post_stages_follow_precision_post(tmp_path, monkeypatch, precision_post):
    expected = precision_post or "tf32"
    calls, calculators, identities = _run(
        tmp_path, monkeypatch, 1, ["--tsopt", "--thermo"], precision_post
    )

    assert [name for name, _, _ in calls] == ["tsopt", "irc", "freq"]
    for _name, child_args, resolved in calls:
        assert resolved == expected
        assert ("--precision" in child_args) is (precision_post is not None)
    endpoint_calc, = calculators
    assert endpoint_calc["precision"] == expected
    assert [cfg["precision"] for cfg in identities] == [expected, expected]


@pytest.mark.parametrize("precision_post", [None, "fp64"], ids=["post-unset", "post-fp64"])
def test_mep_post_stages_follow_precision_post(tmp_path, monkeypatch, precision_post):
    expected = precision_post or "tf32"
    calls, calculators, identities = _run(
        tmp_path, monkeypatch, 2, ["--tsopt", "--thermo"], precision_post
    )

    assert [name for name, _, _ in calls] == ["path-opt", "tsopt", "irc", "freq"]
    (_, mep_args, mep_resolved), *post_calls = calls
    assert "--precision" not in mep_args
    assert mep_resolved == "tf32"
    for _name, child_args, resolved in post_calls:
        assert resolved == expected
        assert ("--precision" in child_args) is (precision_post is not None)
    endpoint_calc, = calculators
    assert endpoint_calc["precision"] == expected
    assert [cfg["precision"] for cfg in identities] == [expected, expected]


def test_summary_records_post_precision(tmp_path):
    from pdb2reaction.io.summary import write_summary_log

    for post, recorded in ((None, "tf32"), ("fp64", "fp64")):
        summary = {"segments": []}
        all_workflow._enrich_summary(
            summary, pipeline_mode="path-opt", out_dir=tmp_path,
            mlip_backend="uma", mlip_model="uma-s-1p2p1", mlip_precision="tf32",
            mlip_precision_post=post, charge=0, spin=1,
            config={"tsopt": True, "thermo": False, "dft": False},
        )
        assert summary["mlip_precision"] == "tf32"
        assert summary["mlip_precision_post"] == recorded

        dest = tmp_path / f"summary_{recorded}.log"
        write_summary_log(dest, {"root_out_dir": str(tmp_path), **summary})
        text = dest.read_text(encoding="utf-8")
        assert "MLIP precision     : tf32" in text
        assert ("MLIP post precision: fp64" in text) is (post is not None)
