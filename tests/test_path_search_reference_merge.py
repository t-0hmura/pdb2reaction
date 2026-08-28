"""Contracts for full-system coordinate composites used for inspection."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from Bio.PDB import PDBParser
from pysisyphus.constants import BOHR2ANG
from pdb2reaction.workflows.path_search import (
    _INSPECTION_MERGE_REMARKS,
    cli,
    _full_system_merge_enabled,
    _merge_pair_to_full,
    _structure_to_arrays,
)


def test_reference_merge_requires_explicit_request_and_alignment() -> None:
    refs = (Path("template.pdb"),)

    assert _full_system_merge_enabled(
        refs, align=True, write_ref_merge=True
    ) is True
    assert _full_system_merge_enabled(
        refs, align=True, write_ref_merge=False
    ) is False
    assert _full_system_merge_enabled(
        refs, align=False, write_ref_merge=True
    ) is False
    assert _full_system_merge_enabled(
        (), align=True, write_ref_merge=True
    ) is False


def test_reference_merge_dependency_is_documented_on_requesting_option() -> None:
    options = {param.name: param for param in cli.params}

    assert "write-ref-merge" not in (options["align"].help or "")
    assert "write-ref-merge" not in (options["ref_pdb_paths"].help or "")
    assert "first -i input" in (options["ref_pdb_paths"].help or "")
    request_help = options["write_ref_merge"].help or ""
    assert "--align" in request_help
    assert "--ref-full-pdb" in request_help


def test_reference_merge_pdb_remark_is_compact_and_unambiguous() -> None:
    assert _INSPECTION_MERGE_REMARKS == (
        "FOR_INSPECTION ACTIVE_SITE_PATH_IN_STATIC_FULL_SYSTEM_TEMPLATE",
    )
    assert all(len(f"REMARK   1 {line}") <= 80 for line in _INSPECTION_MERGE_REMARKS)


def test_reference_merge_writes_inspection_remark(tmp_path: Path) -> None:
    template = tmp_path / "template.pdb"
    template.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       1.500   1.500   0.000  1.00  0.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    parser = PDBParser(QUIET=True)
    struct_a = parser.get_structure("a", str(template))
    struct_b = parser.get_structure("b", str(template))
    coords_a, _, _, keymap_a = _structure_to_arrays(struct_a)
    coords_b, _, _, keymap_b = _structure_to_arrays(struct_b)
    image = SimpleNamespace(coords3d=coords_a / BOHR2ANG)

    blocks, _ = _merge_pair_to_full(
        pair_images=[image],
        model_ref_pdb=template,
        structA=struct_a,
        structB=struct_b,
        coordsA_aligned=coords_a,
        coordsB_aligned=coords_b,
        keymapA=keymap_a,
        keymapB=keymap_b,
        out_path=None,
    )

    assert len(blocks) == 1
    assert (
        "REMARK   1 FOR_INSPECTION ACTIVE_SITE_PATH_IN_STATIC_FULL_SYSTEM_TEMPLATE"
        in blocks[0]
    )
