"""Tests for mmCIF and PDB-size-limit-safe structure I/O."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def _write_minimal_cif(path: Path) -> None:
    path.write_text(
        """data_input
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . SAM LONG_CHAIN 1 . ? 0.0 1.0 2.0 1.00 12.0 . 10001 SAM LONG_CHAIN C1 1
HETATM 2 O O1 . SAM LONG_CHAIN 1 . ? 1.0 1.0 2.0 1.00 13.0 . 10001 SAM LONG_CHAIN O1 1
#
""",
        encoding="utf-8",
    )


def test_cif_is_normalized_and_output_restores_auth_ids(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.core.utils import (
        convert_xyz_to_pdb,
        load_pdb_atom_metadata,
        prepare_input_structure,
    )

    source = tmp_path / "large-id.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "optimized.xyz"
    xyz.write_text("2\nframe\nC 3.0 4.0 5.0\nO 6.0 7.0 8.0\n", encoding="utf-8")
    out_pdb = tmp_path / "optimized.pdb"

    prepared = prepare_input_structure(source)
    try:
        assert prepared.is_cif
        assert prepared.source_path.suffix == ".pdb"
        metadata = load_pdb_atom_metadata(prepared.source_path)
        assert metadata[0]["chain"] == "LONG_CHAIN"
        assert metadata[0]["resseq"] == 10001

        convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)
        out_cif = out_pdb.with_suffix(".cif")
        assert out_cif.exists()
        data = MMCIF2Dict(str(out_cif))
        assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN", "LONG_CHAIN"]
        assert data["_atom_site.auth_seq_id"] == ["10001", "10001"]
        assert np.allclose([float(value) for value in data["_atom_site.Cartn_x"]], [3.0, 6.0])
    finally:
        from pdb2reaction.io.structure_formats import unregister_coordinate_template

        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_pdb_output_rejects_coordinate_field_overflow(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import convert_xyz_to_pdb

    ref = tmp_path / "ref.pdb"
    ref.write_text(
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    xyz = tmp_path / "far.xyz"
    xyz.write_text("1\nframe\nC 10000.0 0.0 0.0\n", encoding="utf-8")

    with pytest.raises(ValueError, match="fixed-column PDB range"):
        convert_xyz_to_pdb(xyz, ref, tmp_path / "out.pdb")


@pytest.mark.parametrize(
    ("invalid_frame", "message"),
    [
        ("1\nlate-count\nC 2.0 0.0 0.0\n", "Atom count mismatch"),
        (
            "2\npermuted\nO 2.0 0.0 0.0\nC 3.0 0.0 0.0\n",
            "Ordered elements differ",
        ),
        (
            "2\nnon-finite\nC nan 0.0 0.0\nO 3.0 0.0 0.0\n",
            "non-finite",
        ),
        (
            "2\noverflow\nC 9999.9999 0.0 0.0\nO 3.0 0.0 0.0\n",
            "fixed-column PDB range",
        ),
    ],
)
def test_xyz_overlay_validates_all_frames_before_any_output_mutation(
    tmp_path: Path,
    invalid_frame: str,
    message: str,
) -> None:
    from pdb2reaction.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from pdb2reaction.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text(
        "2\nvalid-first\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\n" + invalid_frame,
        encoding="utf-8",
    )
    out_pdb = tmp_path / "existing.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"existing-pdb-generation\n")
    out_cif.write_bytes(b"existing-cif-generation\n")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        with pytest.raises(ValueError, match=message):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)
        assert out_pdb.read_bytes() == b"existing-pdb-generation\n"
        assert out_cif.read_bytes() == b"existing-cif-generation\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        from pdb2reaction.io.structure_formats import unregister_coordinate_template

        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_xyz_overlay_valid_multiframe_publishes_pdb_cif_and_registry(
    tmp_path: Path,
) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from pdb2reaction.io.structure_formats import coordinate_template_for

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text(
        "2\nframe-1\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\n"
        "2\nframe-2\nC 2.0 0.0 0.0\nO 3.0 0.0 0.0\n",
        encoding="utf-8",
    )
    out_pdb = tmp_path / "trajectory.pdb"

    prepared = prepare_input_structure(source)
    try:
        convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)
        text = out_pdb.read_text(encoding="utf-8")
        assert text.count("MODEL") == 2
        assert text.count("ENDMDL") == 2
        data = MMCIF2Dict(str(out_pdb.with_suffix(".cif")))
        assert data["_atom_site.pdbx_PDB_model_num"] == ["1", "1", "2", "2"]
        assert coordinate_template_for(out_pdb) is prepared.structure_template
    finally:
        from pdb2reaction.io.structure_formats import unregister_coordinate_template

        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_xyz_overlay_companion_publish_failure_keeps_existing_generation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdb2reaction.core import result_commit
    from pdb2reaction.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from pdb2reaction.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "frame.xyz"
    xyz.write_text("2\nframe\nC 2 0 0\nO 3 0 0\n", encoding="utf-8")
    out_pdb = tmp_path / "existing.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old-pdb\n")
    out_cif.write_bytes(b"old-cif\n")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        replace = result_commit._replace_exact

        def fail_companion(staged: Path, destination: Path) -> None:
            if destination == out_cif:
                raise OSError("injected companion publication failure")
            replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", fail_companion)
        with pytest.raises(result_commit.ResultCommitError, match="companion publication"):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)
        assert out_pdb.read_bytes() == b"old-pdb\n"
        assert out_cif.read_bytes() == b"old-cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        from pdb2reaction.io.structure_formats import unregister_coordinate_template

        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_xyz_overlay_primary_publish_failure_rolls_back_companion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdb2reaction.core import result_commit
    from pdb2reaction.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from pdb2reaction.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "frame.xyz"
    xyz.write_text("2\nframe\nC 2 0 0\nO 3 0 0\n", encoding="utf-8")
    out_pdb = tmp_path / "existing.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old-pdb\n")
    out_cif.write_bytes(b"old-cif\n")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        replace = result_commit._replace_exact
        failed = False

        def fail_primary_once(staged: Path, destination: Path) -> None:
            nonlocal failed
            if destination == out_pdb and not failed:
                failed = True
                raise OSError("injected primary publication failure")
            replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", fail_primary_once)
        with pytest.raises(result_commit.ResultCommitError, match="primary publication"):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)
        assert out_pdb.read_bytes() == b"old-pdb\n"
        assert out_cif.read_bytes() == b"old-cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        from pdb2reaction.io.structure_formats import unregister_coordinate_template

        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_pdb_with_more_than_ten_thousand_residues_uses_safe_bridge(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.core.utils import load_pdb_atom_metadata, prepare_input_structure
    from pdb2reaction.io.structure_formats import write_pdb_as_mmcif
    from pysisyphus.io.pdb import parse_pdb

    source = tmp_path / "ten-thousand.pdb"
    lines = []
    for index in range(10_001):
        chain = "A"
        resseq = index + 1
        serial = index + 1
        lines.append(
            f"ATOM  {serial:5d}  CA  GLY {chain}{resseq:4d}    "
            f"{float(index % 10):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}           C\n"
        )
    lines.append("END\n")
    source.write_text("".join(lines), encoding="utf-8")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.source_path != source
        assert prepared.structure_template is not None
        assert prepared.structure_template.natoms == 10_001
        atoms, _coords, fragments, _atom_map = parse_pdb(str(prepared.source_path))
        assert len(atoms) == 10_001
        assert len(fragments) == 10_001
        metadata = load_pdb_atom_metadata(prepared.source_path)
        assert metadata[0]["chain"] == "A"
        assert metadata[-1]["chain"] == "A"
        assert metadata[-1]["resseq"] == 10_001
        out_cif = tmp_path / "large-output.cif"
        write_pdb_as_mmcif(
            prepared.source_path, prepared.structure_template, out_cif
        )
        data = MMCIF2Dict(str(out_cif))
        assert data["_atom_site.auth_asym_id"][-1] == "A"
        assert data["_atom_site.auth_seq_id"][-1] == "10001"
    finally:
        prepared.cleanup()


def test_bundled_pdb_parser_distinguishes_two_letter_atoms_from_hydrogens(
    tmp_path: Path,
) -> None:
    from pysisyphus.io.pdb import parse_pdb

    records = [
        ("HG  ", "HG"),
        ("HE  ", "HE"),
        (" HG ", "H"),
        (" HE ", "H"),
        (" CA ", "C"),
        (" CD ", "C"),
        (" CE ", "C"),
        (" NE ", "N"),
        (" PA ", "P"),
        (" PB ", "P"),
        (" SG ", "S"),
        (" OG ", "O"),
        (" NH1", "NH"),
        ("ZN  ", "N"),
    ]
    lines = []
    for serial, (name, element) in enumerate(records, start=1):
        lines.append(
            f"HETATM{serial:5d} {name:4s} RES A{serial:4d}    "
            f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
    path = tmp_path / "elements.pdb"
    path.write_text("".join(lines) + "END\n", encoding="utf-8")

    atoms, *_ = parse_pdb(str(path))
    assert atoms == [
        "Hg",
        "He",
        "H",
        "H",
        "C",
        "C",
        "C",
        "N",
        "P",
        "P",
        "S",
        "O",
        "N",
        "Zn",
    ]


def test_chain_qualified_atom_selector_disambiguates_repeated_ids() -> None:
    from pdb2reaction.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "name": "C1"},
        {"chain": "B", "resname": "SAM", "resseq": 12, "name": "C1"},
    ]
    assert resolve_atom_spec_index("B:SAM:12:C1", metadata) == 1


def test_chain_qualified_atom_selector_keeps_duplicate_field_roles_distinct() -> None:
    from pdb2reaction.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "name": "A"},
        {"chain": "B", "resname": "SAM", "resseq": 12, "name": "A"},
    ]
    assert resolve_atom_spec_index("A:SAM:12:A", metadata) == 0


def test_chain_qualified_atom_selector_accepts_insertion_code() -> None:
    from pdb2reaction.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "icode": "A", "name": "C1"},
        {"chain": "A", "resname": "SAM", "resseq": 12, "icode": "B", "name": "C1"},
    ]
    assert resolve_atom_spec_index("A:SAM:12B:C1", metadata) == 1


def test_chain_qualified_atom_selector_keeps_chain_case_distinct() -> None:
    from pdb2reaction.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "name": "C1"},
        {"chain": "a", "resname": "SAM", "resseq": 12, "name": "C1"},
    ]
    assert resolve_atom_spec_index("a:SAM:12:C1", metadata) == 1


def test_atom_label_includes_chain_for_repeated_residues() -> None:
    from pdb2reaction.core.utils import atom_label_from_meta

    metadata = [
        {"chain": "LONG_CHAIN", "resname": "SAM", "resseq": 10001, "name": "C1"}
    ]
    assert atom_label_from_meta(metadata, 0) == "LONG_CHAIN:SAM:10001:C1"


def test_atom_label_includes_insertion_code() -> None:
    from pdb2reaction.core.utils import atom_label_from_meta

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "icode": "B", "name": "C1"}
    ]
    assert atom_label_from_meta(metadata, 0) == "A:SAM:12B:C1"


def test_hybrid36_upper_and_lower_ranges_are_decoded() -> None:
    from pdb2reaction.io.structure_formats import _hy36decode

    assert _hy36decode(4, "A000") == 10_000
    assert _hy36decode(4, "ZZZZ") == 1_223_055
    assert _hy36decode(4, "a000") == 1_223_056
    assert _hy36decode(5, "A0000") == 100_000


def test_internal_atom_serial_wraps_without_six_digit_pdb_field() -> None:
    from pdb2reaction.io.structure_formats import _internal_atom_serial

    assert _internal_atom_serial(0) == 1
    assert _internal_atom_serial(99_998) == 99_999
    assert _internal_atom_serial(99_999) == 1


def test_decimal_overflow_pdb_fields_are_normalized(tmp_path: Path) -> None:
    from pdb2reaction.io.structure_formats import read_pdb_atom_sites

    source = tmp_path / "overflow.pdb"
    source.write_text(
        f"ATOM  {100000:5d}  CA  GLY A{100000:4d}    "
        f"{1.25:8.3f}{2.5:8.3f}{3.75:8.3f}{1.0:6.2f}{4.0:6.2f}           C1+\nEND\n",
        encoding="utf-8",
    )

    records, nonstandard = read_pdb_atom_sites(source)
    assert nonstandard
    assert len(records) == 1
    assert records[0].chain_id == "A"
    assert records[0].resseq == "100000"
    assert records[0].formal_charge == "1"
    assert np.allclose([records[0].x, records[0].y, records[0].z], [1.25, 2.5, 3.75])


def test_duplicate_atom_names_without_altloc_are_preserved(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import prepare_input_structure
    from pdb2reaction.io.structure_formats import (
        pdb_requires_normalization,
        read_pdb_atom_sites,
    )

    source = tmp_path / "small-molecule.pdb"
    source.write_text(
        "ATOM      1  C   MOL     1       0.000   0.000   0.000                       C\n"
        "ATOM      2  H   MOL     1       1.000   0.000   0.000                       H\n"
        "ATOM      3  H   MOL     1      -1.000   0.000   0.000                       H\n"
        "END\n",
        encoding="utf-8",
    )

    records, nonstandard = read_pdb_atom_sites(source)
    assert [record.atom_name for record in records] == ["C", "H", "H"]
    assert not nonstandard
    assert pdb_requires_normalization(source)
    prepared = prepare_input_structure(source)
    try:
        normalized, _ = read_pdb_atom_sites(prepared.source_path)
        assert len(normalized) == 3
        assert [record.element for record in normalized] == ["C", "H", "H"]
    finally:
        prepared.cleanup()


def test_duplicate_atom_bridge_preserves_elements_for_geometry_loader(
    tmp_path: Path,
) -> None:
    from pdb2reaction.core.utils import prepare_input_structure
    from pysisyphus.io.pdb import parse_pdb

    elements = [
        "C", "C", "C", "C", "C", "C", "H", "H", "H", "H",
        "H", "C", "H", "O", "C", "H", "C", "H", "H", "H",
    ]
    source = tmp_path / "repeated-element-names.pdb"
    source.write_text(
        "".join(
            f"ATOM  {serial:5d} {element:>2s}   MOL     1"
            f"{float(serial):12.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
            for serial, element in enumerate(elements, start=1)
        )
        + "END\n",
        encoding="utf-8",
    )

    prepared = prepare_input_structure(source)
    try:
        parsed, *_ = parse_pdb(str(prepared.source_path))
        assert parsed == elements
    finally:
        prepared.cleanup()


def _altloc_atom_line(
    serial: int,
    label: str,
    x: float,
    occupancy: float | str | None,
) -> str:
    if occupancy is None:
        occupancy_field = "      "
    elif isinstance(occupancy, str):
        occupancy_field = f"{occupancy:>6s}"
    else:
        occupancy_field = f"{occupancy:6.2f}"
    return (
        f"ATOM  {serial:5d}  CA {label}ALA A   1    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy_field}{10.0:6.2f}           C\n"
    )


@pytest.mark.parametrize(
    ("occupancy_a", "occupancy_b", "expected_serial"),
    [
        (None, 0.10, 2),
        (0.00, None, 1),
        (0.50, 0.50, 1),
        ("nan", 0.10, 2),
        ("inf", 0.10, 2),
    ],
)
def test_altloc_policy_matches_typed_and_text_preserving_pdb_paths(
    tmp_path: Path,
    occupancy_a: float | str | None,
    occupancy_b: float | str | None,
    expected_serial: int,
) -> None:
    from pdb2reaction.io.pdb_fix import clean_pdb_file
    from pdb2reaction.io.structure_formats import read_pdb_atom_sites

    source = tmp_path / "altloc-policy.pdb"
    source.write_text(
        _altloc_atom_line(1, "A", 1.0, occupancy_a)
        + _altloc_atom_line(2, "B", 2.0, occupancy_b)
        + "END\n",
        encoding="utf-8",
    )
    cleaned = tmp_path / "cleaned.pdb"

    records, _ = read_pdb_atom_sites(source)
    clean_pdb_file(source, cleaned)
    retained = next(
        line
        for line in cleaned.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    )

    assert len(records) == 1
    assert records[0].x == pytest.approx(float(expected_serial))
    assert records[0].altloc == ""
    assert int(retained[6:11]) == expected_serial
    assert retained[16] == " "
    assert float(retained[30:38]) == pytest.approx(records[0].x)


def test_selected_altloc_supersedes_same_atom_blank_record_in_both_paths(
    tmp_path: Path,
) -> None:
    from pdb2reaction.io.pdb_fix import clean_pdb_file
    from pdb2reaction.io.structure_formats import read_pdb_atom_sites

    source = tmp_path / "altloc-shared.pdb"
    source.write_text(
        _altloc_atom_line(1, " ", 9.0, 1.00)
        + _altloc_atom_line(2, "A", 2.0, 0.60)
        + _altloc_atom_line(3, "B", 3.0, 0.40)
        + "END\n",
        encoding="utf-8",
    )
    cleaned = tmp_path / "cleaned.pdb"

    records, _ = read_pdb_atom_sites(source)
    clean_pdb_file(source, cleaned)
    retained = next(
        line
        for line in cleaned.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    )

    assert [record.x for record in records] == [2.0]
    assert int(retained[6:11]) == 2
    assert float(retained[30:38]) == pytest.approx(2.0)


def test_cif_altloc_treats_unknown_occupancy_as_missing(tmp_path: Path) -> None:
    from pdb2reaction.io.structure_formats import read_mmcif_atom_sites

    source = tmp_path / "missing-altloc.cif"
    source.write_text(
        """data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA A ALA A 1 1 ? 1 0 0 ? 0 1 ALA A CA 1
ATOM 2 C CA B ALA A 1 1 ? 2 0 0 0.1 0 1 ALA A CA 1
#
""",
        encoding="utf-8",
    )

    records = read_mmcif_atom_sites(source)
    assert len(records) == 1
    assert records[0].x == pytest.approx(2.0)
    assert records[0].altloc == ""


def test_failed_internal_pdb_write_removes_temporary_bridge(
    tmp_path: Path, monkeypatch
) -> None:
    from pdb2reaction.io import structure_formats

    source = tmp_path / "out-of-range.cif"
    _write_minimal_cif(source)
    source.write_text(
        source.read_text(encoding="utf-8").replace("0.0 1.0 2.0", "10000.0 1.0 2.0"),
        encoding="utf-8",
    )
    bridge_dir = tmp_path / "bridge"

    def controlled_mkdtemp(*, prefix):
        assert prefix == "p2r_structure_"
        bridge_dir.mkdir()
        return str(bridge_dir)

    monkeypatch.setattr(structure_formats.tempfile, "mkdtemp", controlled_mkdtemp)
    with np.testing.assert_raises_regex(ValueError, "fixed-column PDB range"):
        structure_formats.normalize_structure_to_pdb(source)
    assert not bridge_dir.exists()


def test_ref_cif_atom_count_error_cleans_temporary_bridge(
    tmp_path: Path, monkeypatch
) -> None:
    from click import ClickException
    from pdb2reaction.core import utils

    xyz = tmp_path / "one.xyz"
    xyz.write_text("1\nframe\nC 0 0 0\n", encoding="utf-8")
    ref = tmp_path / "two.cif"
    _write_minimal_cif(ref)
    bridge_dirs = []
    normalize = utils.normalize_structure_to_pdb

    def record_bridge(path):
        result = normalize(path)
        bridge_dirs.append(result[2])
        return result

    monkeypatch.setattr(utils, "normalize_structure_to_pdb", record_bridge)
    prepared = utils.prepare_input_structure(xyz)
    try:
        with pytest.raises(ClickException, match="atom count"):
            utils.apply_ref_pdb_override(prepared, ref)
    finally:
        prepared.cleanup()

    assert bridge_dirs
    assert all(not path.exists() for path in bridge_dirs)


def test_sp_dry_run_cleans_temporary_cif_bridge(tmp_path: Path, monkeypatch) -> None:
    from click.testing import CliRunner
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.core import utils

    source = tmp_path / "input.cif"
    _write_minimal_cif(source)
    bridge_dirs = []
    normalize = utils.normalize_structure_to_pdb

    def record_bridge(path):
        result = normalize(path)
        bridge_dirs.append(result[2])
        return result

    monkeypatch.setattr(utils, "normalize_structure_to_pdb", record_bridge)
    result = CliRunner().invoke(
        root_cli,
        ["sp", "-i", str(source), "-q", "0", "--dry-run", "-o", str(tmp_path / "out")],
    )

    assert result.exit_code == 0, result.output
    assert bridge_dirs
    assert all(not path.exists() for path in bridge_dirs)


@pytest.mark.parametrize("extra", [(), ("--freeze-atoms", "not-an-index")])
def test_freq_owns_prepared_cif_through_dry_run_and_preexecution_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    extra: tuple[str, ...],
) -> None:
    from click.testing import CliRunner
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.core import utils

    source = tmp_path / "input.cif"
    _write_minimal_cif(source)
    bridge_dirs = []
    normalize = utils.normalize_structure_to_pdb

    def record_bridge(path):
        result = normalize(path)
        bridge_dirs.append(result[2])
        return result

    monkeypatch.setattr(utils, "normalize_structure_to_pdb", record_bridge)
    result = CliRunner().invoke(
        root_cli,
        [
            "freq", "-i", str(source), "-q", "0", "--dry-run",
            "-o", str(tmp_path / "out"), *extra,
        ],
    )

    assert result.exit_code == (0 if not extra else 2), result.output
    assert bridge_dirs
    assert all(not path.exists() for path in bridge_dirs)


def test_pdb_to_cif_preserves_output_occupancy_and_bfactor(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        write_pdb_as_mmcif,
    )

    pdb = tmp_path / "marked.pdb"
    pdb.write_text(
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000  0.50 99.00           C\nEND\n",
        encoding="utf-8",
    )
    record = AtomSiteRecord(
        group_pdb="HETATM",
        element="C",
        atom_name="C1",
        altloc="",
        resname="SAM",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=12.0,
    )
    template = CoordinateTemplate((record,), pdb, "mmcif", "test")
    out = tmp_path / "marked.cif"
    write_pdb_as_mmcif(pdb, template, out)

    data = MMCIF2Dict(str(out))
    assert data["_atom_site.occupancy"] == ["0.50"]
    assert data["_atom_site.B_iso_or_equiv"] == ["99.00"]


def test_multimodel_pdb_is_written_as_multimodel_cif(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        write_pdb_as_mmcif,
    )

    pdb = tmp_path / "trajectory.pdb"
    pdb.write_text(
        "MODEL        1\n"
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000  1.00 10.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "HETATM    1  C1  SAM A   1       4.000   5.000   6.000  1.00 20.00           C\n"
        "ENDMDL\nEND\n",
        encoding="utf-8",
    )
    record = AtomSiteRecord(
        group_pdb="HETATM",
        element="C",
        atom_name="C1",
        altloc="",
        resname="SAM",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
    )
    template = CoordinateTemplate((record,), pdb, "mmcif", "test")
    out = tmp_path / "trajectory.cif"
    write_pdb_as_mmcif(pdb, template, out)

    data = MMCIF2Dict(str(out))
    assert data["_atom_site.pdbx_PDB_model_num"] == ["1", "2"]
    assert data["_atom_site.Cartn_x"] == ["1.000000", "4.000000"]
    assert data["_atom_site.B_iso_or_equiv"] == ["10.00", "20.00"]


def test_cif_altloc_selection_is_coherent_per_residue(tmp_path: Path) -> None:
    from pdb2reaction.io.structure_formats import read_mmcif_atom_sites

    source = tmp_path / "altloc.cif"
    source.write_text(
        """data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0 0 0 1.0 0 1 ALA A CA 1
ATOM 2 C CB A ALA A 1 1 ? 1 0 0 0.9 0 1 ALA A CB 1
ATOM 3 C CG A ALA A 1 1 ? 2 0 0 0.1 0 1 ALA A CG 1
ATOM 4 C CB B ALA A 1 1 ? 3 0 0 0.2 0 1 ALA A CB 1
ATOM 5 C CD B ALA A 1 1 ? 4 0 0 0.8 0 1 ALA A CD 1
#
""",
        encoding="utf-8",
    )

    records = read_mmcif_atom_sites(source)
    assert [record.atom_name for record in records] == ["CA", "CB", "CG"]
    assert all(record.altloc == "" for record in records)


def test_altloc_selection_preserves_repeated_blank_atom_names(tmp_path: Path) -> None:
    from pdb2reaction.io.structure_formats import read_mmcif_atom_sites

    source = tmp_path / "blank-duplicates.cif"
    source.write_text(
        """data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 H H . MOL A 1 1 ? 0 0 0 1.0 0 1 MOL A H 1
HETATM 2 H H . MOL A 1 1 ? 1 0 0 1.0 0 1 MOL A H 1
HETATM 3 C C1 A MOL A 1 1 ? 2 0 0 0.8 0 1 MOL A C1 1
HETATM 4 C C1 B MOL A 1 1 ? 3 0 0 0.2 0 1 MOL A C1 1
#
""",
        encoding="utf-8",
    )

    records = read_mmcif_atom_sites(source)
    assert [record.atom_name for record in records] == ["H", "H", "C1"]


def test_cif_duplicate_atom_names_survive_internal_pdb_and_structure_load(
    tmp_path: Path,
) -> None:
    from pdb2reaction.core.utils import prepare_input_structure
    from pdb2reaction.workflows.extract import load_structure

    source = tmp_path / "duplicate-names.cif"
    source.write_text(
        """data_duplicate
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 H H . MOL A 1 1 ? 0 0 0 1.0 0 1 MOL A H 1
HETATM 2 H H . MOL A 1 1 ? 1 0 0 1.0 0 1 MOL A H 1
HETATM 3 C C1 . MOL A 1 1 ? 2 0 0 1.0 0 1 MOL A C1 1
#
""",
        encoding="utf-8",
    )

    prepared = prepare_input_structure(source)
    try:
        structure = load_structure(str(prepared.source_path), "duplicate")
    finally:
        prepared.cleanup()

    atoms = list(structure.get_atoms())
    assert len(atoms) == 3
    assert len({atom.get_name() for atom in atoms}) == 3
    assert [atom.xtra["p2r_atom_site"].atom_name for atom in atoms] == ["H", "H", "C1"]


def test_extract_cif_accepts_chain_resname_and_emits_cif(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.workflows.extract import extract_api

    source = tmp_path / "complex.cif"
    _write_minimal_cif(source)
    out_pdb = tmp_path / "new" / "nested" / "model.pdb"
    result = extract_api(
        complex_pdb=[str(source)],
        center="LONG_CHAIN:SAM",
        output=[str(out_pdb)],
        radius=0.1,
        include_h2o=False,
        add_linkh=False,
    )

    out_cif = out_pdb.with_suffix(".cif")
    assert out_pdb.exists()
    assert out_cif.exists()
    assert str(out_cif) in result["outputs"]
    data = MMCIF2Dict(str(out_cif))
    assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN", "LONG_CHAIN"]
    assert data["_atom_site.auth_seq_id"] == ["10001", "10001"]


def test_extract_multi_input_creates_nested_output_directories(tmp_path: Path) -> None:
    from pdb2reaction.workflows.extract import extract_api

    source = tmp_path / "complex.cif"
    _write_minimal_cif(source)
    first = tmp_path / "per_file" / "a" / "model_a.pdb"
    second = tmp_path / "per_file" / "b" / "model_b.pdb"
    extract_api(
        complex_pdb=[str(source), str(source)],
        center="LONG_CHAIN:SAM",
        output=[str(first), str(second)],
        radius=0.1,
        include_h2o=False,
        add_linkh=False,
    )
    assert first.exists()
    assert second.exists()

    combined = tmp_path / "combined" / "nested" / "models.pdb"
    extract_api(
        complex_pdb=[str(source), str(source)],
        center="LONG_CHAIN:SAM",
        output=[str(combined)],
        radius=0.1,
        include_h2o=False,
        add_linkh=False,
    )
    assert combined.exists()
    assert combined.read_text().count("MODEL") == 2


def test_all_dry_run_accepts_cif_input(tmp_path: Path, monkeypatch) -> None:
    from click.testing import CliRunner
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.core import utils

    source = tmp_path / "complex.cif"
    _write_minimal_cif(source)
    bridge_dirs = []
    normalize = utils.normalize_structure_to_pdb

    def record_bridge(path):
        result = normalize(path)
        bridge_dirs.append(result[2])
        return result

    monkeypatch.setattr(utils, "normalize_structure_to_pdb", record_bridge)
    result = CliRunner().invoke(
        root_cli,
        [
            "all",
            "-i",
            str(source),
            "-q",
            "0",
            "--tsopt",
            "True",
            "--dry-run",
            "True",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Structure bridge" in result.output
    assert "Dry-run mode" in result.output
    assert bridge_dirs
    assert all(not path.exists() for path in bridge_dirs)


def test_path_merge_keys_use_original_cif_ids_for_subset_models(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import prepare_input_structure
    from pdb2reaction.workflows.path_search import (
        _load_structures_and_chain_align,
        _model_keys_from_pdb,
    )

    columns = """loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""
    full = tmp_path / "full.cif"
    full.write_text(
        "data_full\n" + columns
        + "ATOM 1 C CA . ALA PROTEIN 1 1 ? 0 0 0 1 0 1 ALA PROTEIN CA 1\n"
        + "HETATM 2 C C1 . SAM LIGAND 2 . ? 1 0 0 1 0 10001 SAM LIGAND C1 1\n#\n",
        encoding="utf-8",
    )
    model = tmp_path / "model.cif"
    model.write_text(
        "data_model\n" + columns
        + "HETATM 1 C C1 . SAM LIGAND 2 . ? 1 0 0 1 0 10001 SAM LIGAND C1 1\n#\n",
        encoding="utf-8",
    )

    full_prepared = prepare_input_structure(full)
    model_prepared = prepare_input_structure(model)
    try:
        _structs, _coords, _atoms, keymaps = _load_structures_and_chain_align(
            [full_prepared.source_path]
        )
        model_keys = _model_keys_from_pdb(model_prepared.source_path)
        assert model_keys == [("SAM", "10001", "", "LIGAND", "C1")]
        assert model_keys[0] in keymaps[0]
    finally:
        full_prepared.cleanup()
        model_prepared.cleanup()


def test_all_segment_copy_emits_cif_for_bridge_input(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from pdb2reaction.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        register_coordinate_template,
        unregister_coordinate_template,
    )
    from pdb2reaction.workflows.all import _copy_structures_to_seg_dir
    from pdb2reaction.workflows._run_session import InvocationManifest

    xyz = tmp_path / "reactant.xyz"
    xyz.write_text("1\nframe\nC 1.0 2.0 3.0\n", encoding="utf-8")
    pdb = tmp_path / "reactant.pdb"
    pdb.write_text(
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000  0.75 42.00           C\nEND\n",
        encoding="utf-8",
    )
    record = AtomSiteRecord(
        group_pdb="HETATM",
        element="C",
        atom_name="C1",
        altloc="",
        resname="SAM",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
    )
    template = CoordinateTemplate((record,), tmp_path / "input.cif", "mmcif", "test")
    register_coordinate_template(pdb, template)
    manifest = InvocationManifest()
    try:
        seg_dir = _copy_structures_to_seg_dir(
            {"R": xyz},
            tmp_path / "result",
            1,
            ".cif",
            manifest=manifest,
        )
        out_pdb = seg_dir / "reactant.pdb"
        out_cif = seg_dir / "reactant.cif"
        assert out_pdb.exists()
        assert out_cif.exists()
        data = MMCIF2Dict(str(out_cif))
        assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN"]
        assert data["_atom_site.auth_seq_id"] == ["10001"]
        assert data["_atom_site.occupancy"] == ["0.75"]
        assert data["_atom_site.B_iso_or_equiv"] == ["42.00"]
        assert manifest.paths("output.public.") == [out_pdb, out_cif]
    finally:
        unregister_coordinate_template(pdb)


def test_path_merge_template_falls_back_to_any_bridged_reference(tmp_path: Path) -> None:
    from pdb2reaction.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        register_coordinate_template,
        unregister_coordinate_template,
    )
    from pdb2reaction.workflows.path_search import _coordinate_template_for_refs

    ordinary = tmp_path / "ordinary.pdb"
    bridged = tmp_path / "bridged.pdb"
    record = AtomSiteRecord(
        group_pdb="ATOM",
        element="C",
        atom_name="CA",
        altloc="",
        resname="GLY",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
    )
    template = CoordinateTemplate((record,), tmp_path / "endpoint.cif", "mmcif", "test")
    register_coordinate_template(bridged, template)
    try:
        assert _coordinate_template_for_refs([ordinary, bridged]) is template
    finally:
        unregister_coordinate_template(bridged)


def test_template_free_generation_removes_prior_cif_companion(tmp_path: Path) -> None:
    from pdb2reaction.io.structure_formats import (
        register_output_template_and_write_cif,
    )

    out_pdb = tmp_path / "current.pdb"
    out_pdb.write_text("END\n", encoding="utf-8")
    old_companion = out_pdb.with_suffix(".cif")
    old_companion.write_text("data_stale\n", encoding="utf-8")

    assert register_output_template_and_write_cif(out_pdb, None) is None
    assert not old_companion.exists()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("occupancy", float("nan"), "non-finite occupancy"),
        ("bfactor", float("inf"), "non-finite B-factor"),
    ],
)
def test_mmcif_render_rejects_nonfinite_atom_metadata(
    tmp_path: Path, field: str, value: float, message: str
) -> None:
    from dataclasses import replace
    from pdb2reaction.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        render_mmcif_frames,
    )

    record = AtomSiteRecord(
        group_pdb="ATOM",
        element="C",
        atom_name="CA",
        altloc="",
        resname="GLY",
        chain_id="A",
        resseq="1",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
    )
    record = replace(record, **{field: value})
    template = CoordinateTemplate((record,), tmp_path / "input.pdb", "pdb", "test")

    with pytest.raises(ValueError, match=message):
        render_mmcif_frames([np.zeros((1, 3))], template)


@pytest.mark.parametrize(
    ("keyword", "values", "message"),
    [
        ("occupancy_frames", [np.array([np.nan])], "Occupancy frame"),
        ("bfactor_frames", [np.array([np.inf])], "B-factor frame"),
    ],
)
def test_mmcif_render_rejects_nonfinite_metadata_frames(
    tmp_path: Path, keyword: str, values, message: str
) -> None:
    from pdb2reaction.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        render_mmcif_frames,
    )

    record = AtomSiteRecord(
        group_pdb="ATOM",
        element="C",
        atom_name="CA",
        altloc="",
        resname="GLY",
        chain_id="A",
        resseq="1",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
    )
    template = CoordinateTemplate((record,), tmp_path / "input.pdb", "pdb", "test")

    with pytest.raises(ValueError, match=message):
        render_mmcif_frames(
            [np.zeros((1, 3))],
            template,
            **{keyword: values},
        )


def test_multimodel_pdb_is_reduced_to_one_geometry(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import prepare_input_structure

    source = tmp_path / "ensemble.pdb"
    source.write_text(
        "MODEL        1\n"
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      1  C   MOL A   1       9.000   0.000   0.000  1.00  0.00           C\n"
        "ENDMDL\nEND\n",
        encoding="utf-8",
    )

    prepared = prepare_input_structure(source)
    try:
        geometry = prepared.geom_path.read_text(encoding="utf-8")
        assert geometry.count("ATOM  ") == 1
        assert "   9.000" not in geometry
    finally:
        prepared.cleanup()


@pytest.mark.parametrize("coordinate", ["nan", "inf", "-inf"])
def test_pdb_reader_rejects_nonfinite_coordinates(
    tmp_path: Path, coordinate: str
) -> None:
    from pdb2reaction.io.structure_formats import read_pdb_atom_sites

    source = tmp_path / "nonfinite.pdb"
    source.write_text(
        "ATOM      1  C   MOL A   1    "
        f"{coordinate:>8}   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Non-finite coordinates"):
        read_pdb_atom_sites(source)


def test_long_mmcif_resname_matches_ligand_charge_mapping(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import prepare_input_structure, resolve_charge_spin

    source = tmp_path / "long-resname.cif"
    _write_minimal_cif(source)
    source.write_text(
        source.read_text(encoding="utf-8").replace(" SAM ", " ABCD "),
        encoding="utf-8",
    )

    prepared = prepare_input_structure(source)
    try:
        charge, spin = resolve_charge_spin(
            prepared,
            charge=None,
            spin=1,
            ligand_charge="ABCD:-1",
        )
        assert (charge, spin) == (-1, 1)
    finally:
        prepared.cleanup()
