"""DMF consumes the same resolved calculator contract as other workflows."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import click
import pytest

import pdb2reaction.backends as backends
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
