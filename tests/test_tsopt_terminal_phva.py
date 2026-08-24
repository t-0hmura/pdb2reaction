"""Terminal PHVA runs only after numerical convergence."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch

from pdb2reaction.workflows import tsopt


class _Geometry:
    atomic_numbers = np.array([1])
    atoms = ["H"]
    cart_coords = np.zeros(3)


def test_optimizer_terminal_phva_carries_the_exact_hessian_for_irc_cache():
    exact_hessian = torch.diag(torch.tensor([-1.0, 2.0, 3.0]))
    optimizer = SimpleNamespace(
        _last_exact_cart_coords=np.zeros(3),
        _last_exact_frequencies_cm=np.array([-100.0, 20.0, 30.0]),
        _last_exact_modes=torch.eye(3),
        _last_rigid_projection_info={},
        cur_H=exact_hessian,
    )

    reused = tsopt._optimizer_exact_frequency_data(optimizer, _Geometry())

    assert reused is not None
    assert len(reused) == 4
    cached_hessian = reused[3]
    assert torch.equal(cached_hessian, exact_hessian)
    assert cached_hessian.data_ptr() != exact_hessian.data_ptr()


def _runner(tmp_path, monkeypatch, *, stalled):
    runner = object.__new__(tsopt.HessianDimer)
    runner.geom = _Geometry()
    runner.dump = False
    runner.optim_all_path = tmp_path / "optimization_all_trj.xyz"
    runner.mode_path = tmp_path / ".dimer_mode.dat"
    runner.out_dir = tmp_path
    runner.vib_dir = tmp_path / "vib"
    runner.freeze_atoms = []
    runner.masses_au_t = torch.ones(1)
    runner.masses_amu = np.ones(1)
    runner.device = torch.device("cpu")
    runner.root = 0
    runner.thresh_loose = "baker"
    runner.thresh = "baker"
    runner.flatten_max_iter = 1
    runner.flatten_loop_bofill = False
    runner.neg_freq_thresh_cm = 5.0
    runner.tr_projection = "none"
    runner.rigid_projection_info = {}
    runner.max_total_cycles = 1
    runner._cycles_spent = 0
    runner.is_converged = False
    runner.is_stalled = False
    runner.stop_reason = ""
    runner.n_imaginary_modes = None
    runner.imaginary_frequencies_cm = []
    runner.saddle_order_verified = False
    runner.prepared_input = None
    runner.ref_pdb = None

    hessian_calls = []
    mode_exports = []

    def fake_hessian(*, allow_reuse):
        hessian_calls.append(allow_reuse)
        return torch.eye(3)

    def fake_loop(_threshold):
        runner._cycles_spent = runner.max_total_cycles
        if stalled:
            runner.is_stalled = True
            runner.stop_reason = "energy plateau"
        return 1, False, False

    runner._calc_full_hessian_cached = fake_hessian
    runner._dimer_loop = fake_loop
    monkeypatch.setattr(
        tsopt,
        "_mode_direction_by_root",
        lambda *args, **kwargs: (np.ones((1, 3)), -100.0),
    )
    monkeypatch.setattr(
        tsopt,
        "_frequencies_cm_and_modes",
        lambda *args, **kwargs: (
            np.array([-100.0, 20.0, 30.0]),
            torch.eye(3),
        ),
    )
    monkeypatch.setattr(
        tsopt,
        "_write_all_imaginary_modes",
        lambda *args, **kwargs: mode_exports.append(True),
    )
    monkeypatch.setattr(
        tsopt,
        "write",
        lambda path, _atoms: Path(path).write_text("final\n", encoding="utf-8"),
    )
    return runner, hessian_calls, mode_exports


def test_dimer_max_cycles_saves_final_structure_and_skips_phva(
    monkeypatch, tmp_path, capsys,
):
    runner, hessian_calls, mode_exports = _runner(
        tmp_path, monkeypatch, stalled=False
    )

    runner.run()

    assert len(hessian_calls) == 1
    assert mode_exports == []
    assert runner.n_imaginary_modes is None
    assert runner.hessian_status == "skipped"
    assert (tmp_path / "final_geometry.xyz").is_file()
    assert "ERROR: Not converged." not in capsys.readouterr().err


def test_dimer_plateau_saves_final_structure_and_skips_phva(
    monkeypatch, tmp_path, capsys,
):
    runner, hessian_calls, mode_exports = _runner(
        tmp_path, monkeypatch, stalled=True
    )

    runner.run()

    assert len(hessian_calls) == 1
    assert mode_exports == []
    assert runner.n_imaginary_modes is None
    assert runner.hessian_status == "skipped"
    assert (tmp_path / "final_geometry.xyz").is_file()
    assert "ERROR: Not converged." not in capsys.readouterr().err
    assert runner.is_stalled is True
