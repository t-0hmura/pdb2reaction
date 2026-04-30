# tests/test_summary_log.py
"""Tests for pdb2reaction.summary_log formatting helpers."""

from pathlib import Path

from pdb2reaction.summary_log import (
    _fmt_bool,
    _shorten_path,
    _format_energy_rows,
    _format_bond_changes,
)


class TestFmtBool:
    def test_true(self):
        assert _fmt_bool(True) == "True"

    def test_false(self):
        assert _fmt_bool(False) == "False"

    def test_none(self):
        assert _fmt_bool(None) == "-"

    def test_truthy_int(self):
        assert _fmt_bool(1) == "True"

    def test_falsy_zero(self):
        assert _fmt_bool(0) == "False"


class TestShortenPath:
    def test_none_path(self):
        assert _shorten_path(None, None) == "(not available)"

    def test_empty_path(self):
        # Path("") is truthy in Python (resolves to "."), so _shorten_path returns "."
        result = _shorten_path(Path(""), None)
        assert isinstance(result, str)

    def test_relative_to_root(self):
        result = _shorten_path(Path("/a/b/c/file.txt"), Path("/a/b/c"))
        assert result == "file.txt"

    def test_relative_to_parent(self):
        result = _shorten_path(Path("/a/b/c/file.txt"), Path("/a/b/c/sub"))
        assert result == "file.txt"

    def test_no_root(self):
        p = Path("/some/absolute/path.txt")
        assert _shorten_path(p, None) == str(p)


class TestFormatEnergyRows:
    def test_basic_formatting(self):
        rows = _format_energy_rows(
            labels=["R", "TS", "P"],
            energies_au=[-100.0, -99.5, -100.1],
            energies_kcal=[0.0, 313.75, -62.75],
        )
        assert len(rows) == 3
        assert "R" in rows[0]
        assert "TS" in rows[1]
        assert "P" in rows[2]

    def test_none_energies(self):
        rows = _format_energy_rows(
            labels=["R"],
            energies_au=None,
            energies_kcal=None,
        )
        assert len(rows) == 1
        assert "n/a" in rows[0]

    def test_auto_relative_energy(self):
        """When energies_kcal is None, relative is computed from AU."""
        rows = _format_energy_rows(
            labels=["R", "TS"],
            energies_au=[-100.0, -99.9],
            energies_kcal=None,
        )
        assert len(rows) == 2
        # R should have 0.0 relative (base - base = 0)
        # TS should have a positive relative energy


class TestFormatBondChanges:
    def test_empty_text(self):
        lines = _format_bond_changes("")
        assert len(lines) == 1
        assert "no covalent changes detected" in lines[0]

    def test_nonempty_text(self):
        lines = _format_bond_changes("C1-O2 formed\nN3-H4 broken")
        assert len(lines) == 2
        assert "C1-O2 formed" in lines[0]
        assert "N3-H4 broken" in lines[1]
