"""DMF consumes the same resolved calculator contract as other workflows."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import click
import numpy as np
import pytest

import pdb2reaction.backends as backends
from pdb2reaction.workflows._outcomes import make_leaf
from pdb2reaction.workflows import path_opt


@pytest.mark.parametrize(
    "config",
    [
        {
            "backend": "orb",
            "model": "orb-sentinel",
            "device": "cpu",
            "precision": "float64",
            "compile_model": True,
            "freeze_atoms": [0],
        },
        {
            "backend": "mace",
            "model": "mace-sentinel",
            "device": "cuda",
            "default_dtype": "float64",
            "hessian_calc_mode": "FiniteDifference",
        },
        {
            "backend": "custom",
            "calc_file": "/tmp/calculator.py",
            "calc_factory": "build",
            "charge": -1,
            "spin": 2,
            "device": "cpu",
            "return_partial_hessian": True,
        },
        {
            "backend": "uma",
            "model": "uma-sentinel",
            "task_name": "omol",
            "workers": 3,
            "workers_per_node": 2,
            "precision": "fp64",
            "solvent": "none",
        },
    ],
)
def test_dmf_forwards_complete_resolved_mapping_without_mutation(
    monkeypatch, config,
) -> None:
    original = deepcopy(config)
    captured = []
    sentinel = object()

    def fake_factory(**kwargs):
        captured.append(kwargs)
        return sentinel

    monkeypatch.setattr(path_opt, "create_ase_calculator", fake_factory)

    result = path_opt._create_dmf_ase_calculator(config)

    assert result is sentinel
    assert captured == [original]
    assert config == original


@pytest.mark.parametrize(
    ("backend", "kwargs", "expected"),
    [
        (
            "orb",
            {
                "model": "orb-sentinel",
                "precision": "float64",
                "compile_model": True,
                "workers": 9,
            },
            {
                "model": "orb-sentinel",
                "precision": "float64",
                "compile_model": True,
            },
        ),
        (
            "mace",
            {
                "model": "mace-sentinel",
                "default_dtype": "float64",
                "precision": "ignored",
            },
            {"model": "mace-sentinel", "default_dtype": "float64"},
        ),
        (
            "custom",
            {
                "calc_file": "/tmp/calculator.py",
                "calc_factory": "build",
                "charge": -1,
                "spin": 2,
                "device": "cpu",
                "model": "ignored",
            },
            {
                "calc_file": "/tmp/calculator.py",
                "calc_factory": "build",
                "charge": -1,
                "spin": 2,
                "device": "cpu",
            },
        ),
        (
            "uma",
            {
                "model": "uma-sentinel",
                "task_name": "omol",
                "workers": 3,
                "workers_per_node": 2,
                "precision": "fp64",
                "freeze_atoms": [0],
            },
            {
                "model": "uma-sentinel",
                "task_name": "omol",
                "workers": 3,
                "workers_per_node": 2,
                "precision": "fp64",
            },
        ),
    ],
)
def test_central_ase_factory_filters_only_for_selected_backend(
    monkeypatch, backend, kwargs, expected,
) -> None:
    captured = []

    def constructor(**accepted):
        captured.append(accepted)
        return accepted

    monkeypatch.setattr(backends, "resolve_backend", lambda value: value)
    monkeypatch.setattr(backends, "_import_cls", lambda *_args: constructor)

    result = backends.create_ase_calculator(backend=backend, **kwargs)

    assert captured == [expected]
    assert result == expected


def test_dmf_solvent_rejection_precedes_calculator_construction(
    monkeypatch, tmp_path: Path,
) -> None:
    constructed = False

    def forbidden(_config):
        nonlocal constructed
        constructed = True
        raise AssertionError("calculator construction must not run")

    monkeypatch.setattr(path_opt, "_create_dmf_ase_calculator", forbidden)

    with pytest.raises(click.ClickException, match="not compatible with --solvent"):
        path_opt._run_dmf_mep(
            [],
            {"backend": "orb", "solvent": "water"},
            tmp_path,
            [],
            4,
            [],
        )

    assert constructed is False


@pytest.mark.parametrize("backend", ["gup", "", None])
def test_dmf_unknown_backend_is_rejected_before_import(
    monkeypatch, tmp_path: Path, backend,
) -> None:
    constructed = False

    def forbidden(_config):
        nonlocal constructed
        constructed = True
        raise AssertionError("calculator construction must not run")

    monkeypatch.setattr(path_opt, "_create_dmf_ase_calculator", forbidden)

    with pytest.raises(click.ClickException, match="either 'cpu' or 'gpu'"):
        path_opt._run_dmf_mep(
            [],
            {"backend": "orb", "solvent": "none"},
            tmp_path,
            [],
            4,
            [],
            dmf_cfg={"backend": backend},
        )

    assert constructed is False


def test_main_dmf_ipopt_options_reach_primary_solve(tmp_path: Path) -> None:
    options = {"tol": 1.0e-8, "print_level": 3}

    resolved = path_opt._main_dmf_ipopt_options(options, tmp_path, 321)

    assert resolved == {
        "tol": 1.0e-8,
        "print_level": 3,
        "output_file": str(tmp_path / "dmf_ipopt.out"),
        "max_iter": 321,
    }
    assert options == {"tol": 1.0e-8, "print_level": 3}


@pytest.mark.parametrize("invalid", [0, -1, 1.5, True, float("nan")])
def test_main_dmf_ipopt_options_reject_invalid_explicit_cycle_budget(
    tmp_path: Path, invalid,
) -> None:
    with pytest.raises(click.BadParameter, match="dmf.max_cycles"):
        path_opt._main_dmf_ipopt_options({}, tmp_path, invalid)


@pytest.mark.parametrize("invalid", [-1.0, float("nan"), float("inf")])
def test_harmonic_fix_atoms_rejects_invalid_force_constant(invalid) -> None:
    from pdb2reaction.workflows.restraints import HarmonicFixAtoms

    with pytest.raises(ValueError, match="finite and non-negative"):
        HarmonicFixAtoms([0], np.zeros((1, 3)), k_fix=invalid)


def test_dmf_interpolation_cache_is_released_without_emptying_mid_phase() -> None:
    calls = []

    class Stage:
        def release_device_cache(self, *, empty_cache):
            calls.append(empty_cache)

    path_opt._release_dmf_interpolation_cache(Stage())
    path_opt._release_dmf_interpolation_cache(object())  # pydmf 1.2 fallback

    assert calls == [False]


def test_torch_dmf_runtime_options_disable_unused_history_and_preserve_precision(
    monkeypatch,
) -> None:
    monkeypatch.setattr(path_opt.torch.cuda, "is_available", lambda: True)
    assert path_opt._torch_dmf_runtime_kwargs(
        "cpu", {"keep_history": True}, {}, {}
    ) == {}
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu",
        {"device": "cpu", "dtype": "float64"},
        {"device": "cuda", "dtype": "float32"},
        {"dtype": "float64"},
    ) == {"keep_history": False, "device": "cuda", "dtype": "float64"}
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu", {"keep_history": True}, {}, {}
    ) == {"keep_history": True, "device": "cuda"}
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu", {"keep_history": True}, {}, {}, supports_keep_history=False
    ) == {"device": "cuda"}


def test_default_gpu_dmf_backend_requires_cuda(monkeypatch) -> None:
    """The public gpu backend must not silently compute on the CPU."""
    monkeypatch.setattr(path_opt.torch.cuda, "is_available", lambda: False)

    with pytest.raises(click.ClickException, match="requires CUDA"):
        path_opt._torch_dmf_runtime_kwargs("gpu", {}, {}, {})

    # An explicit expert device stays in charge, and the cpu backend is unaffected.
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu", {"device": "cpu"}, {}, {}
    ) == {"keep_history": False, "device": "cpu"}
    assert path_opt._torch_dmf_runtime_kwargs("cpu", {}, {}, {}) == {}


def test_requested_preopt_failure_prevents_path_opt_success() -> None:
    preopt = make_leaf(
        "path-opt",
        "preopt_endpoint_0",
        executed=True,
        converged=False,
    )
    mep = make_leaf(
        "path-opt",
        "gsm_mep",
        executed=True,
        converged=True,
    )

    truth, outcomes = path_opt._combine_path_opt_outcomes([preopt], mep)

    assert [outcome.item_id for outcome in outcomes] == [
        "preopt_endpoint_0",
        "gsm_mep",
    ]
    assert truth.scientific_status == "partial"
