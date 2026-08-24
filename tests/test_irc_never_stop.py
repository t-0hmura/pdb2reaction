"""Regression tests for IRC physical-stop bypass and its workflow output."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from pdb2reaction.core.defaults import IRC_KW
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.IRC import IRC


def _irc_stop_probe(*, never_stop: bool, increased: bool, converged: bool) -> IRC:
    irc = object.__new__(IRC)
    irc.never_stop = never_stop
    irc.energy_increased = increased
    irc.energy_converged = converged
    irc.never_stop_energy_increase_bypasses = 0
    irc.never_stop_energy_convergence_bypasses = 0
    return irc


def test_never_stop_is_opt_in() -> None:
    assert IRC_KW["never_stop"] is False
    assert IRC_KW["energy_increase_thresh"] == pytest.approx(0.0)


def test_default_irc_stops_on_energy_increase_or_one_step_energy_change() -> None:
    assert _irc_stop_probe(
        never_stop=False, increased=True, converged=False
    )._energy_stop_message() == "Energy increased!"
    assert _irc_stop_probe(
        never_stop=False, increased=False, converged=True
    )._energy_stop_message() == "Energy converged!"


def test_never_stop_ignores_energy_only_stops() -> None:
    assert _irc_stop_probe(
        never_stop=True, increased=True, converged=False
    )._energy_stop_message() == ""


def test_never_stop_ignores_gradient_stops_without_claiming_convergence() -> None:
    irc = _irc_stop_probe(
        never_stop=True, increased=False, converged=False
    )
    irc.past_inflection = True
    irc.rms_grad_thresh = 1.0e-3
    irc.hard_rms_grad_thresh = 2.0e-3
    irc.converged = False

    assert not irc._gradient_converged(0.0)
    irc.past_inflection = False
    assert not irc._hard_gradient_stop(0.0)
    assert irc.converged is False


def test_never_stop_reports_physical_thresholds_as_disabled(capsys) -> None:
    irc = _irc_stop_probe(
        never_stop=True, increased=False, converged=False
    )

    irc.report_conv_thresholds()

    output = capsys.readouterr().out
    assert "stop thresholds are disabled" in output
    assert "rms(|gradient|) <=" not in output


def test_active_gradient_accepts_list_indices_from_refined_partial_hessian() -> None:
    irc = object.__new__(IRC)
    irc._act_dofs = [0, 2]

    np.testing.assert_array_equal(
        irc.active_gradient(np.array([1.0, 2.0, 3.0])),
        np.array([1.0, 3.0]),
    )


@pytest.mark.parametrize("threshold", [-1.0, np.nan, np.inf])
def test_energy_increase_threshold_must_be_finite_and_nonnegative(
    tmp_path, threshold
) -> None:
    geometry = Geometry(
        ("H",), np.zeros(3), coord_type="cart"
    )

    with pytest.raises(
        ValueError,
        match="must be finite and non-negative",
    ):
        IRC(
            geometry,
            out_dir=tmp_path,
            energy_increase_thresh=threshold,
        )


def test_default_irc_stops_on_any_energy_increase() -> None:
    irc = _irc_stop_probe(
        never_stop=False, increased=False, converged=False
    )
    irc.energy_increase_thresh = IRC_KW["energy_increase_thresh"]

    assert not irc._energy_increase_exceeds_tolerance(-100.0, -100.0)
    assert irc._energy_increase_exceeds_tolerance(-100.0, -99.999999)
    assert _irc_stop_probe(
        never_stop=True, increased=False, converged=True
    )._energy_stop_message() == ""


def test_directional_endpoint_energy_fields_keep_legacy_aliases() -> None:
    from pdb2reaction.workflows.irc import _directional_endpoint_energy_fields

    fields = _directional_endpoint_energy_fields([-10.0, -9.0, -11.0], -8.5)

    assert fields["energy_first_hartree"] == -10.0
    assert fields["energy_last_hartree"] == -11.0
    assert fields["energy_ts_hartree"] == -8.5
    assert fields["endpoint_energy_orientation"] == "finished_first_to_finished_last"
    assert fields["energy_reactant_hartree"] == fields["energy_first_hartree"]
    assert fields["energy_product_hartree"] == fields["energy_last_hartree"]


def test_workflow_uses_engine_normalized_prefix(tmp_path) -> None:
    from pdb2reaction.workflows.irc import _irc_output_path

    irc = object.__new__(IRC)
    irc.out_dir = tmp_path
    irc.prefix = "segment_"

    assert _irc_output_path(irc, "finished_irc_trj.xyz") == (
        tmp_path / "segment_finished_irc_trj.xyz"
    )


def test_engine_writes_first_last_and_stitched_endpoints(tmp_path) -> None:
    irc = object.__new__(IRC)
    irc.out_dir = tmp_path
    irc.prefix = ""
    irc.atoms = ("H",)
    irc._m_sqrt = np.ones(3)
    irc.all_energies = np.array([-1.0, -0.9])
    coords = np.array(
        [
            [[0.0, 0.0, 0.0]],
            [[0.0, 0.0, 1.0]],
        ]
    )

    irc.dump_ends(tmp_path, "finished", coords=coords, trj=True)

    first = (tmp_path / "finished_first.xyz").read_text(encoding="utf-8")
    last = (tmp_path / "finished_last.xyz").read_text(encoding="utf-8")
    trajectory = (tmp_path / "finished_irc_trj.xyz").read_text(encoding="utf-8")
    assert "0.00000000" in first
    assert "0.52917721" in last
    assert trajectory.count("\n1\n") + trajectory.startswith("1\n") >= 2


def test_real_irc_generation_invalidates_prefixed_directional_outputs(
    tmp_path,
) -> None:
    from pdb2reaction.workflows.irc import _prepare_irc_output_dir

    stale = [
        tmp_path / "segment_forward_irc_trj.xyz",
        tmp_path / "segment_backward_irc.pdb",
        tmp_path / "segment_finished_first.xyz",
        tmp_path / "result.json",
        tmp_path / "summary.json",
    ]
    for path in stale:
        path.write_text("stale\n", encoding="utf-8")
    unrelated = tmp_path / "notes.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    assert _prepare_irc_output_dir(tmp_path, prefix="segment") == tmp_path.resolve()

    assert all(not path.exists() for path in stale)
    assert unrelated.read_text(encoding="utf-8") == "keep\n"


def test_irc_output_prefix_cannot_escape_output_directory(tmp_path) -> None:
    import click
    from pdb2reaction.workflows.irc import _prepare_irc_output_dir

    sibling = tmp_path / "escape_finished_irc.pdb"
    sibling.write_text("keep\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    with pytest.raises(click.BadParameter, match="without path components"):
        _prepare_irc_output_dir(out_dir, prefix="../escape")
    assert sibling.read_text(encoding="utf-8") == "keep\n"


def test_irc_original_input_is_protected_from_invalidation(tmp_path) -> None:
    import click
    from pdb2reaction.workflows.irc import _prepare_irc_output_dir

    source = tmp_path / "finished_irc.cif"
    source.write_text("data_input\n", encoding="utf-8")

    with pytest.raises(click.UsageError, match="collides with a reserved"):
        _prepare_irc_output_dir(
            tmp_path,
            protected_inputs=(source,),
        )
    assert source.read_text(encoding="utf-8") == "data_input\n"


@pytest.mark.parametrize(
    "device",
    ["cpu"] + (["cuda"] if torch.cuda.is_available() else []),
)
def test_terminal_hessian_conversion_consumes_tensor_in_place(device) -> None:
    from pdb2reaction.workflows.irc import (
        _consume_mw_hessian_to_cartesian_active,
    )

    hessian = torch.tensor(
        [[2.0, 0.5], [0.5, 3.0]],
        dtype=torch.float64,
        device=device,
    )
    original = hessian.detach().cpu().numpy().copy()
    data_ptr = hessian.data_ptr()
    masses = np.array([2.0, 3.0])

    result = _consume_mw_hessian_to_cartesian_active(hessian, masses)

    assert hessian.data_ptr() == data_ptr
    np.testing.assert_allclose(result, masses[:, None] * original * masses[None, :])
    np.testing.assert_allclose(hessian.detach().cpu().numpy(), result)


def test_terminal_hessian_conversion_consumes_numpy_in_place() -> None:
    from pdb2reaction.workflows.irc import (
        _consume_mw_hessian_to_cartesian_active,
    )

    hessian = np.array([[2.0, 0.5], [0.5, 3.0]])
    original = hessian.copy()
    masses = np.array([2.0, 3.0])

    result = _consume_mw_hessian_to_cartesian_active(hessian, masses)

    assert result is hessian
    np.testing.assert_allclose(result, masses[:, None] * original * masses[None, :])
