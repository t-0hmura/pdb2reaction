"""Every bond judgment of ``all`` reads the YAML ``bond`` section.

The stubs are those of test_all_root_summary_golden.py; the stubbed bond check
also records the bond configuration each call received.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import test_all_root_summary_golden as golden
from pdb2reaction.workflows import all as all_workflow

SEEN: list = []


class RecordingStubs(golden.PipelineStubs):
    def install(self, monkeypatch):
        super().install(monkeypatch)
        stub = all_workflow._path_search.has_bond_change

        def record(first, second, cfg):
            SEEN.append((cfg["bond_factor"], cfg["metal_delta_fraction"]))
            return stub(first, second, cfg)

        monkeypatch.setattr(all_workflow._path_search, "has_bond_change", record)


@pytest.mark.parametrize(
    "frames,calls",
    [
        # MEP: the segment's bond_changes, the IRC endpoint matching (four
        # pairings) and the check of the optimized endpoints (four pairings).
        pytest.param((0, 2), 9, id="mep"),
        # TS-only: the optimized R/P comparison.
        pytest.param((1,), 1, id="ts-only"),
    ],
)
@pytest.mark.parametrize(
    "bond_yaml,expected",
    [
        pytest.param(None, (1.20, 0.40), id="defaults"),
        pytest.param("bond:\n  bond_factor: 1.30\n  metal_delta_fraction: 0.25\n",
                     (1.30, 0.25), id="yaml"),
    ],
)
def test_all_bond_judgments_read_the_yaml_bond_section(
    tmp_path: Path, monkeypatch, frames, calls, bond_yaml, expected
) -> None:
    SEEN.clear()
    monkeypatch.setattr(golden, "PipelineStubs", RecordingStubs)
    extra = []
    if bond_yaml is not None:
        config = tmp_path / "bond.yaml"
        config.write_text(bond_yaml, encoding="utf-8")
        extra = ["--config", str(config)]

    run = golden.invoke(
        tmp_path, monkeypatch, golden.Scenario("bond", frames, ["--tsopt"]), extra=extra
    )

    assert run.exit_code == 0, run.output
    assert SEEN == [expected] * calls
