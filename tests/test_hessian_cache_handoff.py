"""Regression tests for geometry-safe in-process Hessian handoff."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from pdb2reaction.io import hessian_cache


def setup_function() -> None:
    hessian_cache.clear()


def teardown_function() -> None:
    hessian_cache.clear()


def test_endpoint_cache_matches_xyz_round_trip_coordinates() -> None:
    coords = np.array([0.0, 1.234567890, -2.345678901])
    hessian_cache.store(
        "irc_endpoint",
        np.eye(3),
        meta={"cart_coords": coords, "irc_direction": "forward"},
    )
    entry = hessian_cache.load("irc_endpoint")

    # XYZ serialization changes coordinates slightly; the endpoint cache is
    # still valid within the explicit bohr tolerance.
    round_tripped = coords + np.array([0.0, 2.0e-6, -2.0e-6])
    assert entry is not None
    assert hessian_cache.matches_cart_coords(entry, round_tripped)


def test_endpoint_cache_rejects_swapped_or_stale_geometry() -> None:
    hessian_cache.store(
        "irc_endpoint",
        np.eye(3),
        meta={"cart_coords": np.array([0.0, 0.0, 0.0])},
    )
    entry = hessian_cache.load("irc_endpoint")

    assert entry is not None
    assert not hessian_cache.matches_cart_coords(
        entry, np.array([0.0, 0.0, 1.0e-3])
    )
    assert not hessian_cache.matches_cart_coords(entry, np.zeros(6))


def test_cache_without_coordinate_identity_is_not_reused() -> None:
    hessian_cache.store("ts", np.eye(3))
    entry = hessian_cache.load("ts")

    assert entry is not None
    assert not hessian_cache.matches_cart_coords(entry, np.zeros(3))


def test_store_copies_tensor_and_coordinate_metadata() -> None:
    H = torch.eye(3, dtype=torch.float64)
    coords = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)
    hessian_cache.store(
        "irc_endpoint",
        H,
        active_dofs=[0, 1, 2],
        meta={"cart_coords": coords},
    )
    H.add_(5.0)
    coords.add_(5.0)

    entry = hessian_cache.load("irc_endpoint")
    assert entry is not None
    torch.testing.assert_close(entry["hessian"], torch.eye(3, dtype=torch.float64))
    np.testing.assert_allclose(entry["meta"]["cart_coords"], [1.0, 2.0, 3.0])


def test_discard_prevents_missing_endpoint_from_reusing_previous_seed() -> None:
    hessian_cache.store("irc_endpoint", np.eye(3))
    hessian_cache.discard("irc_endpoint")
    assert hessian_cache.load("irc_endpoint") is None


# ---------------------------------------------------------------------------
# Defensive in-process ownership
# ---------------------------------------------------------------------------
def test_load_returns_independent_tensor_snapshots() -> None:
    hessian_cache.store(
        "ts",
        torch.eye(3, dtype=torch.float64),
        meta={"cart_coords": np.zeros(3)},
    )
    first = hessian_cache.load("ts")
    second = hessian_cache.load("ts")

    assert first is not None and second is not None
    assert first["hessian"].data_ptr() != second["hessian"].data_ptr()

    # Mutating a loaded snapshot must not corrupt the retained raw artifact.
    first["hessian"].add_(5.0)
    again = hessian_cache.load("ts")
    torch.testing.assert_close(again["hessian"], torch.eye(3, dtype=torch.float64))


def test_load_snapshot_isolates_numpy_meta_and_active_dofs() -> None:
    hessian_cache.store(
        "ts",
        np.eye(3),
        active_dofs=[0, 1, 2],
        meta={"cart_coords": np.zeros(3), "source": "tsopt_exact"},
    )
    snap = hessian_cache.load("ts")
    snap["meta"]["cart_coords"][0] = 99.0
    snap["active_dofs"].append(999)
    # Mutating a loaded numpy Hessian through torch.as_tensor must not persist.
    torch.as_tensor(snap["hessian"]).mul_(0.0)

    fresh = hessian_cache.load("ts")
    assert fresh["meta"]["cart_coords"][0] == 0.0
    assert fresh["active_dofs"] == [0, 1, 2]
    np.testing.assert_array_equal(fresh["hessian"], np.eye(3))


# ---------------------------------------------------------------------------
# Complete reuse identity
# ---------------------------------------------------------------------------
def _identity(
    *,
    run="run-A",
    backend="uma",
    model="m",
    precision="fp64",
    charge=0,
    spin=1,
    active_atoms=(0, 1),
    active_dofs=(0, 1, 2, 3, 4, 5),
    potential=None,
    constraints=None,
    source="tsopt_exact",
    coords=None,
):
    if coords is None:
        coords = np.zeros(6)
    return hessian_cache.build_identity(
        atoms=[1, 1],
        cart_coords=coords,
        run_id=run,
        backend=backend,
        model=model,
        precision=precision,
        charge=charge,
        spin=spin,
        potential=potential or {},
        active_atoms=list(active_atoms),
        active_dofs=list(active_dofs),
        constraints=constraints or {"freeze_atoms": []},
        source=source,
        method=source,
    )


def test_load_matching_accepts_only_full_identity() -> None:
    hessian_cache.store("ts", np.eye(6), identity=_identity())

    # Exact identity, exact coordinates.
    assert hessian_cache.load_matching("ts", _identity()) is not None
    # Coordinates within the bohr round-trip tolerance.
    assert hessian_cache.load_matching("ts", _identity(coords=np.full(6, 2.0e-6))) is not None
    # Every single-field change rejects.
    assert hessian_cache.load_matching("ts", _identity(coords=np.full(6, 1.0e-3))) is None
    assert hessian_cache.load_matching("ts", _identity(backend="orb")) is None
    assert hessian_cache.load_matching("ts", _identity(model="other")) is None
    assert hessian_cache.load_matching("ts", _identity(precision="fp32")) is None
    assert hessian_cache.load_matching("ts", _identity(charge=-1)) is None
    assert hessian_cache.load_matching("ts", _identity(spin=3)) is None
    assert hessian_cache.load_matching("ts", _identity(active_dofs=(0, 1, 2))) is None
    assert hessian_cache.load_matching("ts", _identity(active_atoms=(0,))) is None
    assert hessian_cache.load_matching(
        "ts", _identity(constraints={"freeze_atoms": [0]})
    ) is None
    assert hessian_cache.load_matching(
        "ts", _identity(potential={"solvent": "water"})
    ) is None
    assert hessian_cache.load_matching("ts", _identity(source="irc_endpoint_quasi_newton")) is None
    assert hessian_cache.load_matching("ts", _identity(run="run-B")) is None


def test_load_matching_returns_defensive_snapshot() -> None:
    hessian_cache.store("ts", np.eye(6), active_dofs=[0, 1, 2, 3, 4, 5], identity=_identity())
    snap = hessian_cache.load_matching("ts", _identity())
    assert snap is not None
    snap["hessian"][0, 0] = 42.0
    snap["active_dofs"].append(7)
    again = hessian_cache.load_matching("ts", _identity())
    np.testing.assert_array_equal(again["hessian"], np.eye(6))
    assert again["active_dofs"] == [0, 1, 2, 3, 4, 5]


def test_legacy_coordinate_only_entry_is_never_reused_by_load_matching() -> None:
    # A legacy entry carries no identity token.
    hessian_cache.store("ts", np.eye(6), meta={"cart_coords": np.zeros(6)})
    assert hessian_cache.load_matching("ts", _identity()) is None
    # The coordinate-only load path still returns a snapshot for legacy callers.
    assert hessian_cache.load("ts") is not None


def test_missing_run_id_never_reuses() -> None:
    hessian_cache.store("ts", np.eye(6), identity=_identity(run=None))
    assert hessian_cache.load_matching("ts", _identity(run=None)) is None
    assert hessian_cache.load_matching("ts", _identity(run="run-A")) is None


def test_identity_from_context_round_trips_through_cache(monkeypatch) -> None:
    from pdb2reaction.core.result_commit import RUN_ID_ENV

    monkeypatch.setenv(RUN_ID_ENV, "run-X")

    class _Geom:
        atomic_numbers = np.array([1, 1, 8])
        cart_coords = np.arange(9, dtype=float)
        freeze_atoms = np.array([0])

    calc_cfg = {
        "backend": "uma",
        "model": "m",
        "precision": "fp64",
        "charge": 0,
        "spin": 1,
        "freeze_atoms": [0],
    }
    hessian_cache.store(
        "ts",
        np.eye(6),
        identity=hessian_cache.identity_from_context(_Geom(), calc_cfg, role="ts"),
    )
    # Same evaluator context reuses.
    assert hessian_cache.load_matching(
        "ts", hessian_cache.identity_from_context(_Geom(), calc_cfg, role="ts")
    ) is not None
    # A different evaluator does not.
    other = dict(calc_cfg, backend="orb")
    assert hessian_cache.load_matching(
        "ts", hessian_cache.identity_from_context(_Geom(), other, role="ts")
    ) is None


@pytest.mark.parametrize(
    "change",
    [
        {"hessian_calc_mode": "Analytical"},
        {"hessian_double": True},
        {"return_partial_hessian": False},
    ],
)
def test_identity_rejects_a_different_hessian_construction_mode(
    change,
    monkeypatch,
) -> None:
    from pdb2reaction.core.result_commit import RUN_ID_ENV

    monkeypatch.setenv(RUN_ID_ENV, "run-hessian-policy")

    class _Geom:
        atoms = ["H", "H"]
        atomic_numbers = [1, 1]
        cart_coords = np.zeros(6)
        freeze_atoms = []

    base = {
        "backend": "uma",
        "model": "uma-s-1p2",
        "precision": "fp32",
        "charge": 0,
        "spin": 1,
        "hessian_calc_mode": "FiniteDifference",
        "hessian_double": False,
        "return_partial_hessian": True,
    }

    base_identity = hessian_cache.identity_from_context(_Geom(), base, role="ts")
    other_identity = hessian_cache.identity_from_context(
        _Geom(), {**base, **change}, role="ts"
    )
    assert (
        base_identity["evaluator"]["potential"]
        != other_identity["evaluator"]["potential"]
    )
    hessian_cache.store("ts", np.eye(6), identity=base_identity)
    assert hessian_cache.load_matching("ts", other_identity) is None
