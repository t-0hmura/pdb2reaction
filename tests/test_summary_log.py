"""Tests for pdb2reaction.io.summary formatting helpers."""

import json
import time
from pathlib import Path

import pytest

from pdb2reaction.io.summary import (
    _fmt_bool,
    _shorten_path,
    _format_energy_rows,
    _format_bond_changes,
    emit_method_citations,
    format_result_warning,
    format_method_citations,
    method_references,
    write_summary_log,
)


class TestFmtBool:
    def test_true(self):
        assert _fmt_bool(True) == "True"

    def test_false(self):
        assert _fmt_bool(False) == "False"

    def test_none(self):
        assert _fmt_bool(None) == "-"

    def test_truthy_int(self):
        assert _fmt_bool(1) == "True"

    def test_falsy_zero(self):
        assert _fmt_bool(0) == "False"


class TestShortenPath:
    def test_none_path(self):
        assert _shorten_path(None, None) == "(not available)"

    def test_empty_path(self):
        # Path("") is truthy in Python (resolves to "."), so _shorten_path returns "."
        result = _shorten_path(Path(""), None)
        assert isinstance(result, str)

    def test_relative_to_root(self):
        result = _shorten_path(Path("/a/b/c/file.txt"), Path("/a/b/c"))
        assert result == "file.txt"

    def test_relative_to_parent(self):
        result = _shorten_path(Path("/a/b/c/file.txt"), Path("/a/b/c/sub"))
        assert result == "file.txt"

    def test_no_root(self):
        p = Path("/some/absolute/path.txt")
        assert _shorten_path(p, None) == str(p)


class TestFormatEnergyRows:
    def test_basic_formatting(self):
        rows = _format_energy_rows(
            labels=["R", "TS", "P"],
            energies_au=[-100.0, -99.5, -100.1],
            energies_kcal=[0.0, 313.75, -62.75],
        )
        assert len(rows) == 3
        assert "R" in rows[0]
        assert "TS" in rows[1]
        assert "P" in rows[2]

    def test_none_energies(self):
        rows = _format_energy_rows(
            labels=["R"],
            energies_au=None,
            energies_kcal=None,
        )
        assert len(rows) == 1
        assert "n/a" in rows[0]

    def test_auto_relative_energy(self):
        """When energies_kcal is None, relative is computed from AU."""
        rows = _format_energy_rows(
            labels=["R", "TS"],
            energies_au=[-100.0, -99.9],
            energies_kcal=None,
        )
        assert len(rows) == 2
        # R should have 0.0 relative (base - base = 0)
        # TS should have a positive relative energy


class TestFormatBondChanges:
    def test_empty_text(self):
        lines = _format_bond_changes("")
        assert len(lines) == 1
        assert "no covalent changes detected" in lines[0]

    def test_nonempty_text(self):
        lines = _format_bond_changes("C1-O2 formed\nN3-H4 broken")
        assert len(lines) == 2
        assert "C1-O2 formed" in lines[0]
        assert "N3-H4 broken" in lines[1]


def test_summary_log_uses_mlip_provenance_and_labels(tmp_path):
    dest = tmp_path / "summary.log"
    payload = {
        "root_out_dir": str(tmp_path),
        "pipeline_mode": "path-opt",
        "mlip_backend": "mace",
        "mlip_model": "mace-off23-small",
        "segments": [{"index": 1, "tag": "seg_01", "kind": "seg"}],
        "post_segments": [
            {
                "index": 1,
                "tag": "seg_01",
                "kind": "seg",
                "mlip": {"barrier_kcal": 12.3, "delta_kcal": -1.2},
                "gibbs_mlip": {"barrier_kcal": 13.4, "delta_kcal": -0.7},
                "gibbs_dft_mlip": {"barrier_kcal": 14.5, "delta_kcal": -0.2},
            }
        ],
        "energy_diagrams": [
            {
                "name": "energy_diagram_MLIP_all",
                "labels": ["R", "TS", "P"],
                "energies_kcal": [0.0, 12.3, -1.2],
                "ylabel": "ΔE (kcal/mol)",
            }
        ],
    }

    write_summary_log(dest, payload)
    text = dest.read_text(encoding="utf-8")

    assert "MLIP backend       : mace" in text
    assert "MLIP model         : mace-off23-small" in text
    assert "MLIP energies (TSOPT+IRC)" in text
    assert "DFT//MLIP Gibbs" in text
    assert "UMA model" not in text
    assert "DFT//UMA" not in text


def test_summary_log_does_not_invent_default_backend(tmp_path):
    dest = tmp_path / "summary.log"

    write_summary_log(dest, {})

    text = dest.read_text(encoding="utf-8")
    assert "MLIP backend       : -" in text
    assert "MLIP backend       : uma" not in text


def test_write_summary_log_marks_non_successful_results(tmp_path):
    dest = tmp_path / "summary.log"

    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "path_module_dir": "path_search",
            "pipeline_mode": "path-search",
            "segments": [{"index": 1, "barrier_kcal": 12.3}],
            "energy_diagrams": [],
            "execution_status": "completed",
            "scientific_status": "partial",
            "scientific_status_reasons": ["segment 1 did not converge"],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "Execution status    : completed" in text
    assert "Scientific status   : partial" in text
    assert "RESULT WARNING      : Segment 1 did not converge." in text
    assert "Status reason" not in text


@pytest.mark.parametrize(
    ("reason", "expected"),
    [
        (
            "all:segment_2:mep_not_converged",
            "Segment 2: MEP optimization did not converge. Review the MEP trajectory and convergence log.",
        ),
        (
            "all:segment_2:endpoint_opt:product_converged",
            "Segment 2: the product endpoint optimization did not converge or could not be confirmed. "
            "Review the endpoint structure and optimizer log.",
        ),
        (
            "all:segment_3:irc:irc:forward:not_converged;irc:backward:energy_invalid",
            "Segment 3: Forward IRC did not converge. Review its trajectory and IRC log. "
            "Backward IRC did not produce a valid energy profile. Review its trajectory and IRC log.",
        ),
        (
            "path-search:endpoint_hei;engine_nonconverged",
            "No reactive segment was identified and the path-search engine did not converge. "
            "Review the path and path-search log.",
        ),
        (
            "missing:segment_4",
            "Segment 4: the expected segment result is missing. Review the workflow outputs.",
        ),
    ],
)
def test_format_result_warning_explains_priority_status_codes(reason, expected):
    assert format_result_warning(reason) == expected


def test_summary_log_tree_lists_only_current_run_paths(tmp_path):
    current = tmp_path / "segments" / "seg_01" / "structures" / "ts.pdb"
    stale = tmp_path / "segments" / "seg_02" / "structures" / "old.pdb"
    stale_diagram = tmp_path / "irc_plot_all.png"
    for path in (current, stale, stale_diagram):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    dest = tmp_path / "summary.log"

    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "current_output_paths": [
                "segments/seg_01/structures/ts.pdb",
            ],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "seg_01/" in text
    assert "ts.pdb" in text
    assert "seg_02" not in text
    assert "irc_plot_all.png" not in text


def test_summary_log_ts_only_uses_refined_provenance_and_top_level_counts(tmp_path):
    dest = tmp_path / "summary.log"
    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "pipeline_mode": "tsopt-only",
            "n_images": 5,
            "n_segments": 1,
            "mep": {"n_images": 0, "n_segments": 0},
            "segments": [{
                "index": 1,
                "tag": "seg_01",
                "kind": "tsopt",
                "barrier_kcal": 8.0,
                "delta_kcal": -1.0,
            }],
            "post_segments": [{
                "index": 1,
                "tag": "seg_01",
                "kind": "tsopt",
                "dft": {"barrier_kcal": 7.8, "delta_kcal": -1.2},
            }],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "Number of IRC frames : 5" in text
    assert "Number of segments   : 1" in text
    assert "refined TS − assigned endpoint" in text
    assert "MEP ΔE" not in text
    assert "DFT ΔE‡" in text
    assert "DFT//MLIP ΔE" not in text


def test_method_citations_follow_resolved_methods_and_match_stdout(
    tmp_path, capsys
):
    payload = {
        "root_out_dir": str(tmp_path),
        "pipeline_mode": "path-search",
        "mep_mode": "dmf",
        "dmf_correlated": True,
        "path_opt_mode": "hess",
        "post_opt_mode": "hess",
        "ts_opt_mode": "hess",
        "endpoint_opt_mode": "hess",
        "post_segments": [
            {
                "endpoint_opt": {"reactant_converged": True},
                "ts_imag": {"n_imag": 1},
                "thermo_symmetry": {"R": {"symmetry_number": 1}},
            }
        ],
    }
    dest = tmp_path / "summary.log"

    write_summary_log(dest, payload)
    lines = format_method_citations(payload)
    references = method_references(payload)
    emit_method_citations(payload)

    text = dest.read_text(encoding="utf-8")
    stdout = capsys.readouterr().out
    block = "\n".join(lines)
    assert block in text
    assert text.rstrip().endswith(block)
    # Same citations, destination-appropriate header: the log keeps its section
    # index, stdout heads the block like every other console section.
    assert stdout.splitlines()[0] == "====== Citations & References ======"
    assert stdout.splitlines()[1:] == lines[1:]
    assert lines[0] == "[6] Methods and citations"
    assert "pdb2reaction:" in block
    assert "Direct Max Flux (DMF)" in block
    assert "FB-ENM initialization for DMF" in block
    assert "Correlated FB-ENM (CFB-ENM)" in block
    assert "RFO / P-RFO" in block
    assert "RS-P-RFO" in block
    assert "Euler predictor-corrector IRC" in block
    assert "quasi-RRHO thermochemistry" in block
    assert "Growing String Method" not in block
    assert all(set(ref) == {"method", "citation", "doi"} for ref in references)
    assert len({ref["doi"] for ref in references}) == len(references)
    assert lines[1] == "Please cite the software and methods used:"
    for index, reference in enumerate(references, start=1):
        offset = 2 * index
        assert lines[offset] == f"- {reference['method']}:"
        assert lines[offset + 1] == f"[{index}] {reference['citation']}"


def test_method_citations_use_actual_path_and_post_stages() -> None:
    path_only = {
        "pipeline_mode": "path-search",
        "mep_mode": "gsm",
        "path_opt_mode": "grad",
        "post_opt_mode": "hess",
        "post_segments": [],
    }
    mixed = {
        **path_only,
        "post_segments": [
            {"endpoint_opt": {}, "ts_imag": {"n_imag": 1}}
        ],
    }

    path_text = "\n".join(format_method_citations(path_only))
    mixed_text = "\n".join(format_method_citations(mixed))

    assert "Limited-memory BFGS (L-BFGS)" in path_text
    assert "RFO / P-RFO" not in path_text
    assert "quasi-RRHO thermochemistry" not in path_text
    assert "Limited-memory BFGS (L-BFGS)" in mixed_text
    assert "RFO / P-RFO" in mixed_text
    assert "quasi-RRHO thermochemistry" not in mixed_text
    assert "Keil, F. J." not in mixed_text
    assert "Chakraborty, A." in mixed_text


def test_dmf_and_post_references_follow_effective_stage_settings() -> None:
    base = {
        "pipeline_mode": "path-search",
        "mep_mode": "dmf",
        "path_opt_mode": "grad",
        "post_segments": [{"endpoint_opt": {}}],
        "ts_opt_mode": "hess",
        "endpoint_opt_mode": "grad",
    }

    uncorrelated = "\n".join(
        format_method_citations({**base, "dmf_correlated": False})
    )
    correlated = "\n".join(
        format_method_citations({**base, "dmf_correlated": True})
    )

    assert "Correlated FB-ENM (CFB-ENM)" not in uncorrelated
    assert "Correlated FB-ENM (CFB-ENM)" in correlated
    assert "RS-P-RFO" in correlated
    assert "Limited-memory BFGS (L-BFGS)" in correlated


def test_final_stdout_places_citations_immediately_before_elapsed(capsys) -> None:
    from pdb2reaction.workflows.all import _emit_final_summary

    _emit_final_summary(
        None,
        time.time(),
        citation_payload={
            "pipeline_mode": "path-search",
            "mep_mode": "gsm",
            "path_opt_mode": "grad",
            "post_segments": [],
        },
    )

    lines = [line for line in capsys.readouterr().out.splitlines() if line]
    assert lines[-1].startswith("[time] Elapsed Time for Whole Pipeline")
    # stdout heads its blocks like every other section; the numbered form
    # belongs to summary.log alone.
    assert "====== Citations & References ======" in lines[:-1]
    assert not any(line.startswith("[6] ") for line in lines)


def test_final_stdout_explains_non_success_scientific_status(
    tmp_path, capsys
) -> None:
    from pdb2reaction.workflows.all import _emit_final_summary

    (tmp_path / "summary.json").write_text(
        json.dumps(
            {
                "execution_status": "completed",
                "scientific_status": "partial",
                "scientific_status_reasons": ["IRC endpoint was not validated"],
            }
        ),
        encoding="utf-8",
    )

    _emit_final_summary(tmp_path, time.time())

    output = capsys.readouterr().out
    assert "Scientific status: partial" in output
    assert "RESULT WARNING: IRC endpoint was not validated." in output
    assert "Status reason:" not in output
    assert output.rstrip().splitlines()[-1].startswith(
        "[time] Elapsed Time for Whole Pipeline"
    )


def test_citation_block_headers_match_their_destination() -> None:
    """`summary.log` numbers its sections, stdout uses `====== ... ======`. The
    two shared one renderer, so the log's section index leaked to the console.
    """
    from pdb2reaction.io.summary import format_method_citations

    payload = {"pipeline_mode": "all", "ts_opt_mode": "rsirfo", "post_segments": []}

    log_block = format_method_citations(payload)
    assert log_block[0] == "[6] Methods and citations"

    stdout_block = format_method_citations(
        payload, header="====== Citations & References ======"
    )
    assert stdout_block[0] == "====== Citations & References ======"
    # Only the header differs; the citations themselves are one source.
    assert log_block[1:] == stdout_block[1:]
