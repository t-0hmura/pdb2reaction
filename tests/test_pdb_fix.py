from pathlib import Path

import pytest

from pdb2reaction.io.pdb_fix import clean_pdb_file, has_altloc
from pdb2reaction.io.structure_formats import read_pdb_atom_sites


def _atom_line(
    serial: int,
    atom_name: str,
    altloc: str,
    resname: str,
    occupancy: float,
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:^4s}{altloc:1s}{resname:>3s} A   1    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"{occupancy:6.2f}{10.0:6.2f}           C\n"
    )


def test_clean_pdb_rejects_same_input_and_output(tmp_path: Path) -> None:
    source = tmp_path / "same.pdb"
    original = _atom_line(1, "CA", "A", "ALA", 1.0) + "END\n"
    source.write_text(original, encoding="utf-8")

    with pytest.raises(ValueError, match="different files"):
        clean_pdb_file(source, source)

    assert source.read_text(encoding="utf-8") == original


def test_clean_pdb_rejects_hardlink_output(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"
    output = tmp_path / "alias.pdb"
    original = _atom_line(1, "CA", "A", "ALA", 1.0) + "END\n"
    source.write_text(original, encoding="utf-8")
    output.hardlink_to(source)

    with pytest.raises(ValueError, match="different files"):
        clean_pdb_file(source, output)

    assert source.read_text(encoding="utf-8") == original


def test_has_altloc_propagates_read_failure(tmp_path: Path) -> None:
    with pytest.raises(IsADirectoryError):
        has_altloc(tmp_path)


def test_blank_duplicate_atoms_survive_altloc_selection(tmp_path: Path) -> None:
    source = tmp_path / "duplicate-shared.pdb"
    source.write_text(
        _atom_line(1, "H", " ", "LIG", 1.0)
        + _atom_line(2, "H", " ", "LIG", 1.0)
        + _atom_line(3, "CA", "A", "LIG", 0.6)
        + _atom_line(4, "CA", "B", "LIG", 0.4)
        + "END\n",
        encoding="utf-8",
    )
    output = tmp_path / "clean.pdb"

    clean_pdb_file(source, output)

    atom_lines = [
        line
        for line in output.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]
    assert [int(line[6:11]) for line in atom_lines] == [1, 2, 3]


def test_microheterogeneous_altlocs_choose_one_residue(tmp_path: Path) -> None:
    source = tmp_path / "microheterogeneity.pdb"
    source.write_text(
        _atom_line(1, "CA", "A", "SER", 0.4)
        + _atom_line(2, "CA", "B", "CYS", 0.6)
        + "END\n",
        encoding="utf-8",
    )
    output = tmp_path / "clean.pdb"

    clean_pdb_file(source, output)

    atom_lines = [
        line
        for line in output.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]
    assert len(atom_lines) == 1
    assert int(atom_lines[0][6:11]) == 2
    assert atom_lines[0][17:20] == "CYS"
    records, _ = read_pdb_atom_sites(source)
    assert len(records) == 1
    assert records[0].resname == "CYS"


def test_segment_id_keeps_colliding_pdb_residues_distinct(tmp_path: Path) -> None:
    def with_segid(line: str, segid: str) -> str:
        core = line.rstrip("\n").ljust(80)
        return core[:72] + f"{segid:<4s}" + core[76:] + "\n"

    source = tmp_path / "segments.pdb"
    source.write_text(
        with_segid(_atom_line(1, "CA", "A", "ALA", 0.4), "S1")
        + with_segid(_atom_line(2, "CA", "B", "ALA", 0.6), "S2")
        + "END\n",
        encoding="utf-8",
    )

    records, _ = read_pdb_atom_sites(source)

    assert len(records) == 2
    assert [record.segid for record in records] == ["S1", "S2"]
