"""Golden snapshots of what ``pdb2reaction all`` publishes at its output root.

Every scenario runs the real ``all`` callback through the root CLI. The child
commands are stubbed at ``_run_cli_main`` and write the files the parent reads,
DFT is stubbed at ``_run_dft_for_state``, the parent calculator returns fixed
energies, endpoint optimization returns its IRC geometry, and images are
placeholders, so no backend runs.

Each snapshot holds the root ``summary.json`` bytes, the root ``summary.log``,
the child calls, the output files and the manifest entries. Only run-specific
values are replaced (the temporary directory, the generated child YAML path,
the package version and the compute environment); key order and every other
value are pinned. Set ``P2R_UPDATE_GOLDEN=1`` to rewrite the snapshots after an
intended change, and review the diff.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.result_commit import (
    RUN_ID_ENV,
    apply_current_run_id,
    serialize_json_bytes,
)
from pdb2reaction.workflows import all as all_workflow
from pdb2reaction.workflows import path_search as path_search_workflow
from pdb2reaction.workflows._outcomes import aggregate_workflow_truth, attach_outcomes

GOLDEN_DIR = Path(__file__).parent / "golden" / "all_root_summary"
UPDATE_GOLDEN = os.environ.get("P2R_UPDATE_GOLDEN") == "1"
RUN_ID = "golden-run"

ELEMENTS = ("C", "O")
# Frame k: minima at even k, saddles at odd k (hartree).
FRAME_ENERGIES = (-100.0, -99.95, -100.01, -99.955, -100.02)
# Thermal corrections to G (hartree) and the DFT energy offset.
FRAME_G_CORRECTIONS = (0.010, 0.008, 0.011, 0.009, 0.012)
DFT_OFFSET = -0.25


def coords(index: int) -> List[tuple]:
    return [(index * .125, .25, -.5), (1.25 + index * .25, .5, .75)]


def frame_of(cart_ang: Sequence[Sequence[float]]) -> int:
    """Recover the frame index from the first atom's x coordinate."""
    index = int(round(float(cart_ang[0][0]) / .125))
    if not 0 <= index < len(FRAME_ENERGIES):
        raise AssertionError(f"Coordinates are not a fixture frame: {cart_ang}")
    return index


def fixture_frame(cart_ang) -> Optional[int]:
    """The fixture frame index of ``cart_ang``, or None for other geometries."""
    cart = np.asarray(cart_ang, dtype=float).reshape(-1, 3)
    if not np.isfinite(cart).all():
        return None
    index = int(round(float(cart[0][0]) / .125))
    if 0 <= index < len(FRAME_ENERGIES) and np.allclose(cart, coords(index), atol=1e-4):
        return index
    return None


def energy_of(cart_ang) -> float:
    index = fixture_frame(cart_ang)
    if index is not None:
        return FRAME_ENERGIES[index]
    # Geometries off the fixture frames (diagnostic iterates) get a finite energy.
    cart = np.asarray(cart_ang, dtype=float).reshape(-1, 3)
    return -100.0 + 1e-3 * float(np.linalg.norm(cart[1] - cart[0]))


def xyz_frame(index: int) -> str:
    atoms = "".join(
        f"{element} {x:.6f} {y:.6f} {z:.6f}\n"
        for element, (x, y, z) in zip(ELEMENTS, coords(index))
    )
    return f"{len(ELEMENTS)}\nE={FRAME_ENERGIES[index]:.6f} unit=hartree frame{index}\n{atoms}"


def pdb_frame(index: int) -> str:
    return "".join(
        f"HETATM{serial:5d} {name:<4s} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n"
        for serial, (name, element, (x, y, z)) in enumerate(
            zip(("C1", "O1"), ELEMENTS, coords(index)), 1
        )
    ) + "END\n"


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def read_first_coords(path: Path) -> np.ndarray:
    """First-frame Cartesian coordinates (Å) of an XYZ or PDB fixture file."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if Path(path).suffix.lower() == ".pdb":
        rows = [
            (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            for line in lines
            if line.startswith(("ATOM", "HETATM"))
        ]
    else:
        count = int(lines[0].split()[0])
        rows = [tuple(float(v) for v in line.split()[1:4]) for line in lines[2:2 + count]]
    return np.asarray(rows, dtype=float)


def option_value(args: Sequence[str], flag: str) -> Optional[str]:
    args = list(args)
    return args[args.index(flag) + 1] if flag in args else None


def ts_result(energy: float, *, converged: bool = True, files: Optional[dict] = None) -> dict:
    if not converged:
        return {
            "optimization_status": "not_converged", "converged": False,
            "execution_status": "completed", "scientific_status": "failed",
            "hessian_status": "skipped", "n_imaginary_modes": None,
            "energy_hartree": energy, "n_opt_cycles": 40, "max_cycles": 40,
            "stop_reason": "max_cycles_reached", "files": files or {},
        }
    return {
        "optimization_status": "converged", "converged": True,
        "execution_status": "completed", "scientific_status": "success",
        "hessian_status": "completed", "saddle_validation": "first_order",
        "n_imaginary_modes": 1, "n_negative_modes": 1,
        "imaginary_frequencies_cm": [-500.0], "frequency_zero_cutoff_cm": 5.0,
        "imaginary_mode_criterion": "exact_phva",
        "imaginary_frequency_threshold_cm": 5.0,
        "reaction_mode_index": 0, "reaction_mode_frequency_cm": -500.0,
        "reaction_mode_overlap": 0.9, "reaction_mode_source": "mep-reference-overlap",
        "energy_hartree": energy, "n_opt_cycles": 7, "max_cycles": 100000,
        "stop_reason": "converged", "files": files or {},
    }


class EnergyCalculator:
    """Parent calculator stub: fixed energies by frame, nothing else."""

    def __init__(self, closed: List[str]):
        self._closed = closed

    def get_energy(self, atoms, coords_bohr, **kwargs):
        cart = np.asarray(coords_bohr, dtype=float).reshape(-1, 3) * all_workflow.BOHR2ANG
        return {"energy": energy_of(cart)}

    def close(self):
        self._closed.append("close")


@dataclass
class Scenario:
    name: str
    input_frames: tuple
    args: List[str]
    path_search_segments: Optional[List[dict]] = None
    tsopt_converged: bool = True
    endpoint_failure: Optional[str] = None
    endpoint_not_converged: Optional[str] = None
    reversed_irc_segments: tuple = ()
    dft_failed_states: tuple = ()
    resume_from: Optional[List[str]] = None


@dataclass
class Run:
    out: Path
    exit_code: int
    output: str
    children: List[List[str]] = field(default_factory=list)
    endpoint_calls: List[str] = field(default_factory=list)
    dft_calls: List[str] = field(default_factory=list)
    closed: List[str] = field(default_factory=list)


class PipelineStubs:
    """Install the child, DFT, calculator and image stubs for one invocation."""

    def __init__(self, tmp_path: Path, scenario: Scenario, endpoint_opt=None):
        self.tmp_path = tmp_path
        self.scenario = scenario
        self.children: List[List[str]] = []
        self.endpoint_calls: List[str] = []
        self.dft_calls: List[str] = []
        self.closed: List[str] = []
        self.endpoint_opt = endpoint_opt

    # ----- child commands -----
    def child(self, name, _cli, args, **kwargs):
        args = [str(arg) for arg in args]
        self.children.append([name, *args])
        out = Path(option_value(args, "--out-dir"))
        handler = {
            "path-opt": self._path_opt,
            "path_search": self._path_search,
            "scan": self._scan,
            "tsopt": self._tsopt,
            "irc": self._irc,
            "freq": self._freq,
        }[name]
        handler(out, args)
        return 0

    def _path_opt(self, out: Path, args):
        # The parent reads the child's MEP convergence from this result.json.
        assert "--out-json" in args
        index = int(out.name.split("_")[1])
        first = 2 * (index - 1)
        write(out / "final_geometries_trj.xyz", "".join(xyz_frame(k) for k in (first, first + 1, first + 2)))
        write(out / "hei.xyz", xyz_frame(first + 1))
        write(out / "hei.pdb", pdb_frame(first + 1))
        write(out / "result.json", json.dumps(apply_current_run_id({
            "optimization_status": "converged", "converged": True,
            "preopt_requested": False, "preopt_converged": None,
            "path_optimizers": [],
            "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        })))

    def _path_search(self, out: Path, args):
        segments = self.scenario.path_search_segments
        all_frames: List[int] = []
        for segment in segments:
            index = int(segment["index"])
            first = 2 * (index - 1)
            frames = (first, first + 1, first + 2)
            trajectory = "".join(xyz_frame(k) for k in frames)
            write(out / f"mep_seg_{index:02d}_trj.xyz", trajectory)
            write(out / f"seg_{index:03d}_mep" / "final_geometries_trj.xyz", trajectory)
            write(out / f"hei_seg_{index:02d}.xyz", xyz_frame(first + 1))
            write(out / f"hei_seg_{index:02d}.pdb", pdb_frame(first + 1))
            all_frames.extend(frames if not all_frames else frames[1:])
        if all_frames:
            write(out / "mep_trj.xyz", "".join(xyz_frame(k) for k in all_frames))
            (out / "energy_diagram_MEP.png").write_bytes(b"PNG")
        summary = {
            "out_dir": str(out), "n_images": len(all_frames),
            "n_segments": len(segments), "segments": segments,
            "path_optimizers": ["lbfgs"],
        }
        # The verdict comes from the real path-search outcome helpers, so a run
        # without a reactive segment is not a success.
        reports = [
            path_search_workflow.SegmentReport(
                tag=segment["tag"], barrier_kcal=segment["barrier_kcal"],
                delta_kcal=segment["delta_kcal"], summary=segment["bond_changes"],
                kind=segment["kind"], seg_index=int(segment["index"]),
                converged=segment["converged"],
            )
            for segment in segments
        ]
        leaves, expected = path_search_workflow._path_leaves_and_expected(reports)
        attach_outcomes(summary, truth=aggregate_workflow_truth(leaves, expected),
                        stage_outcomes=leaves)
        write(out / "summary.json", json.dumps(apply_current_run_id(summary)))

    def _scan(self, out: Path, args):
        stage = out / "stage_01"
        write(stage / "result.xyz", xyz_frame(2))
        write(stage / "result.pdb", pdb_frame(2))
        write(out / "result.json", json.dumps(apply_current_run_id({
            "execution_status": "completed", "scientific_status": "success",
            "preopt_converged": None,
            "stages": [{"index": 1, "n_steps": 4, "converged": True,
                        "optimizer_status": "converged",
                        "final_energy_hartree": FRAME_ENERGIES[2]}],
            "stage_outcomes": [{"stage": "scan", "item_id": "stage_1",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        })))

    def _tsopt(self, out: Path, args):
        frame = frame_of(read_first_coords(Path(option_value(args, "-i"))))
        geometry = write(out / "final_geometry.xyz", xyz_frame(frame))
        write(out / "result.json", json.dumps(apply_current_run_id(ts_result(
            FRAME_ENERGIES[frame], converged=self.scenario.tsopt_converged,
            files={"final_geometry": str(geometry)},
        ))))

    def _irc(self, out: Path, args):
        ts = frame_of(read_first_coords(Path(option_value(args, "-i"))))
        frames = [ts - 1, ts, ts + 1]
        segment = (ts + 1) // 2
        if segment in self.scenario.reversed_irc_segments:
            frames.reverse()
        trajectory = write(out / "finished_irc_trj.xyz", "".join(xyz_frame(k) for k in frames))
        write(out / "result.json", json.dumps(apply_current_run_id({
            "files": {"finished_irc": str(trajectory)},
            "forward_requested": True, "forward_integration_converged": True,
            "forward_integration_stop_reason": "converged",
            "backward_requested": True, "backward_integration_converged": True,
            "backward_integration_stop_reason": "converged",
            "n_frames_forward": 1, "n_frames_backward": 1,
        })))

    def _freq(self, out: Path, args):
        cart = read_first_coords(Path(option_value(args, "-i")))
        frame = fixture_frame(cart)
        correction = 0.01 if frame is None else FRAME_G_CORRECTIONS[frame]
        out.mkdir(parents=True, exist_ok=True)
        (out / "thermoanalysis.yaml").write_text(
            "sum_EE_and_thermal_free_energy_ha: "
            f"{energy_of(cart) + correction!r}\n"
            f"thermal_correction_free_energy_ha: {correction!r}\n"
            "symmetry_number: 1\nsymmetry_number_source: fixture\n"
            "point_group: C*v\npoint_group_source: fixture\n",
            encoding="utf-8",
        )
        saddle = frame is not None and frame % 2 == 1
        frequencies = ("-500.0", "120.0", "980.0") if saddle else ("15.0", "120.0", "980.0")
        (out / "frequencies_cm-1.txt").write_text(
            "".join(f"{i} {value}\n" for i, value in enumerate(frequencies)),
            encoding="utf-8",
        )

    # ----- in-process stand-ins -----
    def dft(self, pdb_path, q_int, spin, out_dir, args_yaml, **kwargs):
        state = Path(out_dir).name
        self.dft_calls.append(state)
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cart = read_first_coords(Path(pdb_path))
        if state in self.scenario.dft_failed_states:
            return {"_dft_converged": False, "_dft_failed": True, "_dft_returncode": 1}
        energy = energy_of(cart) + DFT_OFFSET
        (out_dir / "result.yaml").write_text(
            f"energy:\n  hartree: {energy!r}\n  converged: true\n", encoding="utf-8"
        )
        return {
            "energy": {"hartree": energy, "converged": True},
            "_dft_converged": True, "_dft_failed": False, "_dft_returncode": 0,
        }

    def calculator(self, **cfg):
        return EnergyCalculator(self.closed)

    def optimize_endpoint(self, geom, _mode, out_dir, tag, *, dump=False, outcome=None, **kwargs):
        if self.endpoint_opt is not None:
            return self.endpoint_opt(self, geom, _mode, out_dir, tag, dump=dump, outcome=outcome)
        self.endpoint_calls.append(tag)
        work = Path(out_dir) / f"{tag}_rfo_opt"
        work.mkdir(parents=True, exist_ok=True)
        if self.scenario.endpoint_failure and self.scenario.endpoint_failure in tag:
            raise ValueError(f"endpoint optimization failed for {tag}")
        converged = not (
            self.scenario.endpoint_not_converged
            and self.scenario.endpoint_not_converged in tag
        )
        if outcome is not None:
            outcome.update({
                "optimization_status": "converged" if converged else "not_converged",
                "converged": converged,
                "n_opt_cycles": 3 + len(self.endpoint_calls),
                "max_cycles": 100000,
            })
        final = work / "final.xyz"
        final.write_text(geom.as_xyz() + "\n", encoding="utf-8")
        return geom, final, converged

    def install(self, monkeypatch):
        monkeypatch.setenv(RUN_ID_ENV, RUN_ID)
        monkeypatch.setattr(all_workflow, "_run_cli_main", self.child)
        monkeypatch.setattr(all_workflow, "_run_dft_for_state", self.dft)
        monkeypatch.setattr(all_workflow, "create_calculator", self.calculator)
        monkeypatch.setattr(all_workflow, "_optimize_endpoint_geom", self.optimize_endpoint)
        monkeypatch.setattr(
            all_workflow, "write_plotly_image",
            lambda _fig, png, **_kw: Path(png).write_bytes(b"PNG"),
        )

        def trj2fig(_trajectory, outputs, **_kwargs):
            for output in outputs:
                Path(output).parent.mkdir(parents=True, exist_ok=True)
                Path(output).write_bytes(b"PNG")

        monkeypatch.setattr(all_workflow, "run_trj2fig", trj2fig)
        monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)

        def bond_change(first, second, _cfg):
            d1 = np.linalg.norm(np.diff(np.asarray(first.cart_coords).reshape(-1, 3), axis=0))
            d2 = np.linalg.norm(np.diff(np.asarray(second.cart_coords).reshape(-1, 3), axis=0))
            changed = abs(d1 - d2) * all_workflow.BOHR2ANG > 0.15
            return changed, "C1-O1 distance change" if changed else ""

        monkeypatch.setattr(all_workflow._path_search, "has_bond_change", bond_change)


def inputs_for(tmp_path: Path, frames: Sequence[int]) -> List[Path]:
    return [write(tmp_path / f"input{k}.pdb", pdb_frame(k)) for k in frames]


def invoke(tmp_path: Path, monkeypatch, scenario: Scenario, extra: Sequence[str] = (),
           endpoint_opt=None) -> Run:
    """Run ``all`` once for ``scenario`` and return what it did."""
    inputs = inputs_for(tmp_path, scenario.input_frames)
    out = tmp_path / "out"
    base = ["all", *[arg for path in inputs for arg in ("-i", str(path))],
            "-q", "0", "-m", "1", "--out-dir", str(out),
            "--no-preopt", "--no-freeze-links", "--no-tsopt-from-mep-tan"]
    if scenario.resume_from is not None:
        first = PipelineStubs(tmp_path, scenario)
        first.install(monkeypatch)
        warmup = CliRunner().invoke(root_cli, [*base, *scenario.resume_from])
        assert warmup.exit_code == 0, warmup.output + repr(warmup.exception)
    stubs = PipelineStubs(tmp_path, scenario, endpoint_opt=endpoint_opt)
    stubs.install(monkeypatch)
    result = CliRunner().invoke(root_cli, [*base, *scenario.args, *extra])
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    return Run(out=out, exit_code=result.exit_code, output=result.output,
               children=stubs.children, endpoint_calls=stubs.endpoint_calls,
               dft_calls=stubs.dft_calls, closed=stubs.closed)


# ----- normalization -----
def _replace_roots(value: Any, roots: Sequence[str]) -> Any:
    if isinstance(value, str):
        for root in roots:
            value = value.replace(root, "<TMP>")
        return value
    if isinstance(value, list):
        return [_replace_roots(item, roots) for item in value]
    if isinstance(value, dict):
        return {key: _replace_roots(item, roots) for key, item in value.items()}
    return value


def _roots(tmp_path: Path) -> List[str]:
    return sorted({str(tmp_path), str(tmp_path.resolve())}, key=len, reverse=True)


def normalized_summary_json(raw: bytes, tmp_path: Path) -> str:
    """The published bytes with only the run-specific values replaced."""
    payload = json.loads(raw.decode("utf-8"))
    assert serialize_json_bytes(payload) == raw, "summary.json is not in the public format"
    for key, placeholder in (
        ("pdb2reaction_version", "<VERSION>"),
        ("environment", "<ENVIRONMENT>"),
    ):
        if key in payload:
            payload[key] = placeholder
    payload = _replace_roots(payload, _roots(tmp_path))
    return serialize_json_bytes(payload).decode("utf-8")


def normalized_text(text: str, tmp_path: Path) -> str:
    from pdb2reaction import __version__

    for root in _roots(tmp_path):
        text = text.replace(root, "<TMP>")
    return text.replace(f"pdb2reaction {__version__}", "pdb2reaction <VERSION>")


def normalized_argv(argv: Sequence[str], tmp_path: Path) -> List[str]:
    argv = _replace_roots(list(argv), _roots(tmp_path))
    if "--config" in argv:
        argv[argv.index("--config") + 1] = "<CHILD_YAML>"
    return argv


def snapshot(run: Run, tmp_path: Path) -> Dict[str, Any]:
    out = run.out
    manifest = json.loads((out / "_work" / "_run_manifest.json").read_text(encoding="utf-8"))
    produced = [
        [key, _replace_roots(entry["path"], _roots(tmp_path))]
        for key, entry in manifest["produced"].items()
    ]
    log_path = out / "summary.log"
    return {
        "exit_code": run.exit_code,
        "children": [normalized_argv(call, tmp_path) for call in run.children],
        "endpoint_calls": list(run.endpoint_calls),
        "dft_calls": list(run.dft_calls),
        "files": sorted(
            path.relative_to(out).as_posix()
            for path in out.rglob("*")
            if path.is_file() and path.name != "run.log"
        ),
        "manifest_produced": produced,
        "summary_json": normalized_summary_json((out / "summary.json").read_bytes(), tmp_path),
        "summary_log": (
            normalized_text(log_path.read_text(encoding="utf-8"), tmp_path)
            if log_path.is_file() else None
        ),
    }


SEGMENT_1 = {"index": 1, "tag": "seg_001", "kind": "seg", "converged": True,
             "barrier_kcal": 31.38, "delta_kcal": -6.28,
             "bond_changes": "C1-O1 distance change"}
SEGMENT_2 = {"index": 2, "tag": "seg_002", "kind": "seg", "converged": True,
             "barrier_kcal": 34.51, "delta_kcal": -6.28,
             "bond_changes": "C1-O1 distance change"}
BRIDGE_1 = {"index": 1, "tag": "seg_001", "kind": "bridge", "converged": True,
            "barrier_kcal": 0.5, "delta_kcal": -0.1, "bond_changes": ""}

POST = ["--tsopt", "--thermo", "--dft"]

SCENARIOS = [
    # MEP via path-opt
    Scenario("mep_no_post", (0, 2), []),
    Scenario("mep_post_two_segments", (0, 2, 4), POST, reversed_irc_segments=(2,)),
    Scenario("mep_tsopt_not_converged", (0, 2), ["--tsopt"], tsopt_converged=False),
    Scenario("mep_endpoint_failure", (0, 2), ["--tsopt", "--thermo"], endpoint_failure="product"),
    Scenario("mep_endpoint_not_converged", (0, 2), ["--tsopt"], endpoint_not_converged="reactant"),
    Scenario("mep_dft_failure", (0, 2), POST, dft_failed_states=("TS",)),
    Scenario("scan_post", (0,), ["--scan-lists", "[(1,2,1.97)]", "--tsopt", "--thermo"]),
    # MEP via path-search (--refine-path)
    Scenario("refine_no_post", (0, 2), ["--refine-path"], path_search_segments=[SEGMENT_1]),
    Scenario("refine_post", (0, 2), ["--refine-path", *POST],
             path_search_segments=[SEGMENT_1]),
    Scenario("refine_no_segments", (0, 2), ["--refine-path", "--tsopt"], path_search_segments=[]),
    Scenario("refine_no_reactive_segment", (0, 2), ["--refine-path", "--tsopt"],
             path_search_segments=[BRIDGE_1]),
    # TS-only (one structure with --tsopt)
    Scenario("ts_only_post", (1,), POST),
    Scenario("ts_only_tsopt_not_converged", (1,), ["--tsopt", "--thermo"], tsopt_converged=False),
    Scenario("ts_only_endpoint_failure", (1,), ["--tsopt", "--thermo"], endpoint_failure="product"),
    Scenario("ts_only_dft_failure", (1,), POST, dft_failed_states=("P",)),
    Scenario("ts_only_reversed_irc", (1,), POST, reversed_irc_segments=(1,)),
    # Resume of the post-processing on a saved MEP
    Scenario("mep_resume_segment", (0, 2), ["--resume-segment", "1", "--tsopt", "--thermo"],
             resume_from=["--tsopt"]),
]


@pytest.mark.parametrize("scenario", SCENARIOS, ids=[s.name for s in SCENARIOS])
def test_all_root_outputs_match_golden(tmp_path, monkeypatch, scenario):
    run = invoke(tmp_path, monkeypatch, scenario)
    current = snapshot(run, tmp_path)
    golden_json = GOLDEN_DIR / f"{scenario.name}.summary.json"
    golden_meta = GOLDEN_DIR / f"{scenario.name}.meta.json"
    golden_log = GOLDEN_DIR / f"{scenario.name}.summary.log"
    meta = {key: value for key, value in current.items()
            if key not in {"summary_json", "summary_log"}}
    if UPDATE_GOLDEN:
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        golden_json.write_text(current["summary_json"], encoding="utf-8")
        golden_meta.write_text(json.dumps(meta, indent=1) + "\n", encoding="utf-8")
        if current["summary_log"] is None:
            golden_log.unlink(missing_ok=True)
        else:
            golden_log.write_text(current["summary_log"], encoding="utf-8")
    assert current["summary_json"] == golden_json.read_text(encoding="utf-8")
    assert meta == json.loads(golden_meta.read_text(encoding="utf-8"))
    expected_log = golden_log.read_text(encoding="utf-8") if golden_log.is_file() else None
    assert current["summary_log"] == expected_log
