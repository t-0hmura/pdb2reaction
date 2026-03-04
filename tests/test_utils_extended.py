"""Extended unit tests for pdb2reaction.utils pure functions."""

from __future__ import annotations

import pytest


class TestDistanceTag:
    def test_basic(self):
        from pdb2reaction.utils import distance_tag

        assert distance_tag(1.50) == "150"

    def test_zero(self):
        from pdb2reaction.utils import distance_tag

        assert distance_tag(0.0) == "000"

    def test_custom_digits(self):
        from pdb2reaction.utils import distance_tag

        assert distance_tag(1.5, digits=3) == "1500"

    def test_custom_pad(self):
        from pdb2reaction.utils import distance_tag

        assert distance_tag(1.5, pad=5) == "00150"

    def test_rounding(self):
        from pdb2reaction.utils import distance_tag

        assert distance_tag(1.555) == "156"
        assert distance_tag(1.554) == "155"


class TestNormalizeFreezeAtoms:
    def test_none(self):
        from pdb2reaction.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms(None) == []

    def test_string_csv(self):
        from pdb2reaction.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms("1,3,5") == [1, 3, 5]

    def test_string_spaces(self):
        from pdb2reaction.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms("1 3 5") == [1, 3, 5]

    def test_list(self):
        from pdb2reaction.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms([0, 2, 4]) == [0, 2, 4]

    def test_empty_string(self):
        from pdb2reaction.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms("") == []

    def test_non_iterable(self):
        from pdb2reaction.utils import normalize_freeze_atoms

        assert normalize_freeze_atoms(42) == []


class TestMergeFreezeAtomGroups:
    def test_merge_overlapping(self):
        from pdb2reaction.utils import merge_freeze_atom_groups

        assert merge_freeze_atom_groups([1, 3], [2, 3]) == [1, 2, 3]

    def test_empty(self):
        from pdb2reaction.utils import merge_freeze_atom_groups

        assert merge_freeze_atom_groups() == []

    def test_single_group(self):
        from pdb2reaction.utils import merge_freeze_atom_groups

        assert merge_freeze_atom_groups([3, 1, 2]) == [1, 2, 3]


class TestAsList:
    def test_none(self):
        from pdb2reaction.utils import as_list

        assert as_list(None) == []

    def test_list(self):
        from pdb2reaction.utils import as_list

        assert as_list([1, 2]) == [1, 2]

    def test_tuple(self):
        from pdb2reaction.utils import as_list

        assert as_list((1, 2)) == [1, 2]

    def test_non_iterable(self):
        from pdb2reaction.utils import as_list

        assert as_list(42) == []


class TestYamlValid:
    def test_valid_with_keys(self, tmp_path):
        import yaml
        from pdb2reaction.all import _yaml_valid

        f = tmp_path / "test.yaml"
        f.write_text(yaml.dump({"segments": [1], "out_dir": "/tmp"}))
        assert _yaml_valid(f, required_keys=("segments",)) is True

    def test_missing_key(self, tmp_path):
        import yaml
        from pdb2reaction.all import _yaml_valid

        f = tmp_path / "test.yaml"
        f.write_text(yaml.dump({"out_dir": "/tmp"}))
        assert _yaml_valid(f, required_keys=("segments",)) is False

    def test_corrupted(self, tmp_path):
        from pdb2reaction.all import _yaml_valid

        f = tmp_path / "bad.yaml"
        f.write_text("{{{{invalid yaml")
        assert _yaml_valid(f) is False

    def test_nonexistent(self, tmp_path):
        from pdb2reaction.all import _yaml_valid

        assert _yaml_valid(tmp_path / "nope.yaml") is False

    def test_no_required_keys(self, tmp_path):
        import yaml
        from pdb2reaction.all import _yaml_valid

        f = tmp_path / "test.yaml"
        f.write_text(yaml.dump({"any": "thing"}))
        assert _yaml_valid(f) is True


class TestStageDone:
    def test_all_files_exist(self, tmp_path):
        from pdb2reaction.all import _stage_done

        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_text("x")
        f2.write_text("y")
        assert _stage_done([f1, f2]) is True

    def test_missing_file(self, tmp_path):
        from pdb2reaction.all import _stage_done

        f1 = tmp_path / "a.txt"
        f1.write_text("x")
        f2 = tmp_path / "missing.txt"
        assert _stage_done([f1, f2]) is False

    def test_empty_dir(self, tmp_path):
        from pdb2reaction.all import _stage_done

        d = tmp_path / "empty_dir"
        d.mkdir()
        assert _stage_done([d]) is False

    def test_nonempty_dir(self, tmp_path):
        from pdb2reaction.all import _stage_done

        d = tmp_path / "full_dir"
        d.mkdir()
        (d / "file.txt").write_text("content")
        assert _stage_done([d]) is True

    def test_empty_list(self):
        from pdb2reaction.all import _stage_done

        assert _stage_done([]) is True
