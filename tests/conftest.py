"""Shared pytest fixtures for pdb2reaction test suite."""

import pytest


@pytest.fixture(autouse=True)
def _reset_allow_charge_mult_mismatch():
    """Reset the process-global ``--allow-charge-mult-mismatch`` toggle before each test.

    The flag is applied as a process-global side effect via an eager CLI callback (like
    ``--deterministic``). The bool-compat CLI tests invoke many commands with that option's
    negative/value forms, which fires the callback and leaves the global set; without this reset
    it would leak into later tests (e.g. ``validate_charge_spin``)."""
    from pdb2reaction.core.utils import set_allow_charge_mult_mismatch

    set_allow_charge_mult_mismatch(False)
    yield
