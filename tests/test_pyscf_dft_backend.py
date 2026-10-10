"""Small CPU contracts for the stateful PySCF backend."""

from __future__ import annotations

import json
import shutil
import numpy as np
import pytest


pytest.importorskip("pyscf")


def _settings(**overrides):
    values = {
        "func_basis": "hf/sto-3g",
        "engine": "cpu",
        "charge": 0,
        "multiplicity": 1,
        "density_fit": False,
        "save_scf_checkpoint": False,
    }
    values.update(overrides)
    return values


def test_he_energy_force_and_analytical_hessian() -> None:
    from pdb2reaction.backends import create_calculator

    calc = create_calculator(
        backend="dft",
        dft_settings=_settings(),
        hessian_calc_mode="Analytical",
        out_hess_torch=False,
        print_timing=False,
    )
    energy = calc.get_energy(["He"], np.zeros(3))["energy"]
    forces = calc.get_forces(["He"], np.zeros(3))["forces"]
    hessian = calc.get_hessian(["He"], np.zeros(3))["hessian"]

    assert energy == pytest.approx(-2.80778395754, abs=1.0e-9)
    assert np.linalg.norm(forces) < 1.0e-10
    assert np.asarray(hessian).shape == (3, 3)
    assert len(calc.session.metrics) == 1
    assert calc.session._last_good is None


def test_exact_cache_and_next_geometry_reuse_previous_density() -> None:
    from pdb2reaction.backends import create_calculator

    calc = create_calculator(
        backend="dft",
        dft_settings=_settings(func_basis="lda/sto-3g"),
        print_timing=False,
    )
    first = np.array([0.0, 0.0, -0.7, 0.0, 0.0, 0.7])
    second = np.array([0.0, 0.0, -0.72, 0.0, 0.0, 0.72])
    calc.get_forces(["H", "H"], first)
    calc.get_forces(["H", "H"], first)
    calc.get_forces(["H", "H"], second)

    assert len(calc.session.metrics) == 2
    assert calc.session.metrics[0]["guess_source"] == "fresh"
    assert calc.session.metrics[1]["guess_source"] == "previous_density"


@pytest.mark.parametrize("from_checkpoint", [False, True])
def test_gpu_lowmem_passes_reused_density_to_rebuilt_method(
    monkeypatch, from_checkpoint
) -> None:
    from pdb2reaction.backends.pyscf_dft import PySCFDFTSession
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    settings = resolve_dft_settings({
        "backend": "dft",
        "dft": {"func_basis": "lda/sto-3g", "engine": "gpu"},
    })
    density = object()

    class FakeMethod:
        converged = True

        def __init__(self):
            self.dm0 = None

        def make_rdm1(self):
            return density

        def kernel(self, dm0=None):
            self.dm0 = dm0
            return -1.0

    session = PySCFDFTSession(settings)
    rebuilt = FakeMethod()
    if from_checkpoint:
        session._pending_checkpoint = True
        session._last_good = {"loaded": True}
        monkeypatch.setattr(
            session, "_restore_last_good_into_scanner", lambda mol: None
        )
    else:
        session._scanner = FakeMethod()
        session._using_rks_lowmem = True

    def build_method(mol, mm_coords, mm_charges):
        session._scanner = rebuilt
        session._using_rks_lowmem = True

    monkeypatch.setattr(session, "_build_method", build_method)
    monkeypatch.setattr(session, "_update_mm_mol", lambda *args: None)

    _, guess_source = session._run_scf(object(), None, None)

    assert rebuilt.dm0 is density
    assert guess_source == (
        "checkpoint" if from_checkpoint else "previous_density"
    )


def test_checkpoint_is_opt_in_and_round_trips(tmp_path) -> None:
    from pdb2reaction.backends import create_calculator
    from pdb2reaction.backends.base import BackendError

    atoms = (["He"], np.zeros((1, 3)))
    disabled = create_calculator(
        backend="dft", dft_settings=_settings(), print_timing=False
    )
    disabled.get_energy(atoms[0], atoms[1].reshape(-1))
    with pytest.raises(BackendError, match="disabled"):
        disabled.save_scf_checkpoint(tmp_path / "disabled.chk", atoms)

    checkpoint = tmp_path / "state.chk"
    enabled_settings = _settings(
        save_scf_checkpoint=True, checkpoint_path=str(checkpoint)
    )
    enabled = create_calculator(
        backend="dft", dft_settings=enabled_settings, print_timing=False
    )
    enabled.get_energy(atoms[0], atoms[1].reshape(-1))
    enabled.close()
    restored = create_calculator(
        backend="dft", dft_settings=enabled_settings, print_timing=False
    )
    restored.get_energy(atoms[0], atoms[1].reshape(-1))

    assert checkpoint.is_file()
    assert checkpoint.with_suffix(".chk.json").is_file()
    assert restored.session.metrics[0]["guess_source"] == "checkpoint"

    mismatched = create_calculator(
        backend="dft", dft_settings=enabled_settings, print_timing=False
    )
    assert not mismatched.load_scf_checkpoint(
        checkpoint, (["He"], np.array([[0.1, 0.0, 0.0]]))
    )
    assert mismatched.session.checkpoint_status["reason"] == "coordinates_mismatch"


def test_last_good_checkpoint_survives_lost_scanner(tmp_path) -> None:
    from pdb2reaction.backends import create_calculator

    checkpoint = tmp_path / "last-good.chk"
    atoms = (["He"], np.zeros((1, 3)))
    calc = create_calculator(
        backend="dft",
        dft_settings=_settings(save_scf_checkpoint=True),
        print_timing=False,
    )
    calc.get_energy(atoms[0], atoms[1].reshape(-1))
    calc.session._scanner = None

    calc.save_scf_checkpoint(checkpoint, atoms)

    assert checkpoint.is_file()
    assert checkpoint.with_suffix(".chk.json").is_file()


def test_checkpoint_rejects_mixed_binary_and_metadata_generations(tmp_path) -> None:
    from pdb2reaction.backends import create_calculator

    atoms = (["He"], np.zeros((1, 3)))
    first_path = tmp_path / "first.chk"
    second_path = tmp_path / "second.chk"
    for path in (first_path, second_path):
        calc = create_calculator(
            backend="dft",
            dft_settings=_settings(save_scf_checkpoint=True),
            print_timing=False,
        )
        calc.get_energy(atoms[0], atoms[1].reshape(-1))
        calc.save_scf_checkpoint(path, atoms)
        calc.close()

    shutil.copyfile(
        second_path.with_suffix(".chk.json"),
        first_path.with_suffix(".chk.json"),
    )
    metadata_path = first_path.with_suffix(".chk.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["schema"] = 1
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    restored = create_calculator(
        backend="dft",
        dft_settings=_settings(save_scf_checkpoint=True),
        print_timing=False,
    )

    assert not restored.load_scf_checkpoint(first_path, atoms)
    assert restored.session.checkpoint_status["reason"] == "schema_mismatch"
    metadata["schema"] = 2
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    assert not restored.load_scf_checkpoint(first_path, atoms)
    assert restored.session.checkpoint_status["reason"] == "generation_mismatch"


def test_explicit_checkpoint_load_does_not_make_it_an_automatic_save_target(
    tmp_path,
) -> None:
    from pdb2reaction.backends import create_calculator

    checkpoint = tmp_path / "ts.chk"
    atoms = (["He"], np.zeros((1, 3)))
    producer = create_calculator(
        backend="dft",
        dft_settings=_settings(
            save_scf_checkpoint=True, checkpoint_path=str(checkpoint)
        ),
        print_timing=False,
    )
    producer.get_energy(atoms[0], atoms[1].reshape(-1))
    producer.close()
    original_checkpoint = checkpoint.read_bytes()
    original_metadata = checkpoint.with_suffix(".chk.json").read_bytes()

    consumer = create_calculator(
        backend="dft",
        dft_settings=_settings(save_scf_checkpoint=True, checkpoint_path=None),
        print_timing=False,
    )
    assert consumer.load_scf_checkpoint(checkpoint, atoms)
    consumer.get_energy(["He"], np.array([0.1, 0.0, 0.0]))
    consumer.close()

    assert checkpoint.read_bytes() == original_checkpoint
    assert checkpoint.with_suffix(".chk.json").read_bytes() == original_metadata


def test_pcm_force_path_uses_native_pyscf_solvent() -> None:
    from pdb2reaction.backends import create_calculator

    calc = create_calculator(
        backend="dft",
        dft_settings=_settings(solvent="water", solvent_model="pcm"),
        print_timing=False,
    )
    result = calc.get_forces(["He"], np.zeros(3))
    assert np.isfinite(result["energy"])
    assert np.all(np.isfinite(result["forces"]))


def test_settings_forward_pyscf_object_config() -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    settings = resolve_dft_settings(
        {
            "backend": "dft",
            "dft": {
                "func_basis": "hf/sto-3g",
                "engine": "cpu",
                "pyscf": {
                    "mf": {"max_cycle": 42},
                    "grids": {"level": 1},
                },
            },
        }
    )
    assert settings.max_cycle == 42
    assert settings.grid_level == 1


def test_checkpoint_path_is_assigned_by_the_workflow_not_the_settings() -> None:
    from pdb2reaction.core.dft_settings import resolve_dft_settings

    settings = resolve_dft_settings(
        {"backend": "dft", "dft": {"save_scf_checkpoint": True}}
    )
    assert settings.checkpoint_path is None


def test_stepwise_grid_density_converges_a_coarse_stage_first() -> None:
    from pyscf import dft, gto

    from pdb2reaction.backends.pyscf_dft import (
        SCF_STEPWISE_CONV_TOL,
        SCF_STEPWISE_GRID_LEVEL,
        stepwise_grid_density,
    )

    mol = gto.M(
        atom="O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59", basis="sto-3g", verbose=0
    )
    built = []

    def make_method():
        mf = dft.RKS(mol)
        mf.xc = "lda"
        mf.grids.level = 3
        built.append(mf)
        return mf

    density = stepwise_grid_density(make_method)

    assert density is not None
    assert len(built) == 1
    assert built[0].grids.level == SCF_STEPWISE_GRID_LEVEL
    assert built[0].nlcgrids.level == SCF_STEPWISE_GRID_LEVEL
    assert built[0].conv_tol == SCF_STEPWISE_CONV_TOL


def test_stepwise_grid_density_falls_back_without_grid_or_convergence() -> None:
    from pdb2reaction.backends.pyscf_dft import stepwise_grid_density

    class Grids:
        level = 3

    class Unconverged:
        converged = False
        grids = Grids()

        def kernel(self):
            return 0.0

    class NoGrid:
        def kernel(self):
            raise AssertionError("a method without a grid is not run")

    assert stepwise_grid_density(Unconverged) is None
    assert stepwise_grid_density(NoGrid) is None


def test_stepwise_grid_applies_only_to_the_first_scf() -> None:
    from pdb2reaction.backends import create_calculator

    first = np.array([0.0, 0.0, -0.7, 0.0, 0.0, 0.7])
    second = np.array([0.0, 0.0, -0.72, 0.0, 0.0, 0.72])
    normal = create_calculator(
        backend="dft",
        dft_settings=_settings(func_basis="lda/sto-3g", scf_stepwise_grid=False),
        print_timing=False,
    )
    staged = create_calculator(
        backend="dft",
        dft_settings=_settings(func_basis="lda/sto-3g", scf_stepwise_grid=True),
        print_timing=False,
    )
    e_normal = normal.get_energy(["H", "H"], first)["energy"]
    e_staged = staged.get_energy(["H", "H"], first)["energy"]
    staged.get_energy(["H", "H"], second)

    assert e_staged == pytest.approx(e_normal, abs=1.0e-6)
    assert staged.session.metrics[0]["cycles"] < normal.session.metrics[0]["cycles"]
    assert staged.session.metrics[0]["guess_source"] == "fresh"
    assert staged.session.metrics[1]["guess_source"] == "previous_density"


def test_stepwise_grid_is_skipped_for_hartree_fock() -> None:
    from pdb2reaction.backends import create_calculator

    plain = create_calculator(
        backend="dft", dft_settings=_settings(scf_stepwise_grid=False), print_timing=False
    )
    calc = create_calculator(
        backend="dft",
        dft_settings=_settings(scf_stepwise_grid=True),
        print_timing=False,
    )
    plain.get_energy(["He"], np.zeros(3))
    calc.get_energy(["He"], np.zeros(3))

    assert calc.session.metrics[0]["cycles"] == plain.session.metrics[0]["cycles"]
