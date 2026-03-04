"""Unit tests for pdb2reaction.bond_changes pure helper functions."""

from __future__ import annotations

import pytest


class TestBondStr:
    def test_one_based(self):
        from pdb2reaction.bond_changes import _bond_str

        result = _bond_str(0, 1, ["C", "O"], one_based=True)
        assert result == "C1-O2"

    def test_zero_based(self):
        from pdb2reaction.bond_changes import _bond_str

        result = _bond_str(0, 1, ["C", "O"], one_based=False)
        assert result == "C0-O1"

    def test_same_element(self):
        from pdb2reaction.bond_changes import _bond_str

        result = _bond_str(2, 5, ["H", "H", "C", "N", "O", "C"], one_based=True)
        assert result == "C3-C6"


class TestElementArrays:
    def test_basic(self):
        from pdb2reaction.bond_changes import _element_arrays

        elems, radii = _element_arrays(["C", "O", "N"])
        assert elems == ["C", "O", "N"]
        assert len(radii) == 3
        # Carbon covalent radius should be positive
        assert radii[0] > 0

    def test_lowercase_normalization(self):
        from pdb2reaction.bond_changes import _element_arrays

        elems, radii = _element_arrays(["c", "o"])
        assert elems == ["C", "O"]

    def test_empty(self):
        from pdb2reaction.bond_changes import _element_arrays

        elems, radii = _element_arrays([])
        assert elems == []
        assert len(radii) == 0


class TestResolveDevice:
    def test_cpu(self):
        import torch
        from pdb2reaction.bond_changes import _resolve_device

        dev = _resolve_device("cpu")
        assert dev == torch.device("cpu")

    def test_auto_returns_device(self):
        import torch
        from pdb2reaction.bond_changes import _resolve_device

        dev = _resolve_device("auto")
        assert isinstance(dev, torch.device)
