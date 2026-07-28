"""Extended unit tests for pdb2reaction.core.utils pure functions."""

from __future__ import annotations



class TestDistanceTag:
    def test_basic(self):
        from pdb2reaction.core.utils import distance_tag

        assert distance_tag(1.50) == "150"

    def test_zero(self):
        from pdb2reaction.core.utils import distance_tag

        assert distance_tag(0.0) == "000"

    def test_custom_digits(self):
        from pdb2reaction.core.utils import distance_tag

        assert distance_tag(1.5, digits=3) == "1500"

    def test_custom_pad(self):
        from pdb2reaction.core.utils import distance_tag

        assert distance_tag(1.5, pad=5) == "00150"

    def test_rounding(self):
        from pdb2reaction.core.utils import distance_tag

        assert distance_tag(1.555) == "156"
        assert distance_tag(1.554) == "155"


class TestNormalizeFreezeAtoms:
    def test_none(self):
        from pdb2reaction.core.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms(None) == []

    def test_string_csv(self):
        from pdb2reaction.core.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms("1,3,5") == [1, 3, 5]

    def test_string_spaces(self):
        from pdb2reaction.core.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms("1 3 5") == [1, 3, 5]

    def test_list(self):
        from pdb2reaction.core.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms([0, 2, 4]) == [0, 2, 4]

    def test_empty_string(self):
        from pdb2reaction.core.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms("") == []

    def test_non_iterable(self):
        from pdb2reaction.core.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms(42) == []


class TestMergeFreezeAtomGroups:
    def test_merge_overlapping(self):
        from pdb2reaction.core.utils import merge_freeze_atom_groups

        assert merge_freeze_atom_groups([1, 3], [2, 3]) == [1, 2, 3]

    def test_empty(self):
        from pdb2reaction.core.utils import merge_freeze_atom_groups

        assert merge_freeze_atom_groups() == []

    def test_single_group(self):
        from pdb2reaction.core.utils import merge_freeze_atom_groups

        assert merge_freeze_atom_groups([3, 1, 2]) == [1, 2, 3]


class TestAsList:
    def test_none(self):
        from pdb2reaction.core.utils import as_list

        assert as_list(None) == []

    def test_list(self):
        from pdb2reaction.core.utils import as_list

        assert as_list([1, 2]) == [1, 2]

    def test_tuple(self):
        from pdb2reaction.core.utils import as_list

        assert as_list((1, 2)) == [1, 2]

    def test_non_iterable(self):
        from pdb2reaction.core.utils import as_list

        assert as_list(42) == []

