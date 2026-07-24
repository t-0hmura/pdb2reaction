"""Unit tests for pdb2reaction.utils helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def test_pretty_block_with_numpy_scalars() -> None:
    """pretty_block should normalize NumPy scalars for YAML-safe dumping."""
    from pdb2reaction.core.utils import pretty_block, set_verbose_level

    set_verbose_level(3)  # pretty_block (config dump) returns "" below -v 3
    try:
        payload = {
            "freeze_atoms": [np.int64(0), np.int64(3)],
            "ratio": np.float64(1.25),
        }
        text = pretty_block("freeze_atoms (effective)", payload)

        assert "freeze_atoms" in text
        assert "- 0" in text
        assert "- 3" in text
        assert "ratio: 1.25" in text
    finally:
        set_verbose_level(0)


def _prepared_with_template(
    tmp_path: Path,
    *,
    source_name: str = "input.pdb",
    template_charge: int = 7,
    template_spin: int = 3,
):
    from pdb2reaction.core.utils import GjfTemplate, PreparedInputStructure

    source_path = tmp_path / source_name
    source_path.write_text("")
    geom_path = tmp_path / "geom.xyz"
    geom_path.write_text("")
    template = GjfTemplate(
        path=tmp_path / "template.gjf",
        prefix_lines=[],
        suffix_lines=[],
        coord_lines=[],
        charge=template_charge,
        spin=template_spin,
    )
    return PreparedInputStructure(
        source_path=source_path,
        geom_path=geom_path,
        gjf_template=template,
    )


def test_resolve_charge_spin_prefers_ligand_derivation_over_gjf(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from pdb2reaction.core import utils as u

    prepared = _prepared_with_template(tmp_path)
    monkeypatch.setattr(
        u,
        "_derive_charge_from_ligand_charge",
        lambda *_args, **_kwargs: -2,
    )

    charge, spin = u.resolve_charge_spin(
        prepared,
        charge=None,
        spin=None,
        ligand_charge="LIG:-2",
    )

    assert charge == -2
    assert spin == 3


def test_resolve_charge_spin_falls_back_to_gjf_after_ligand_derivation(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from pdb2reaction.core import utils as u

    prepared = _prepared_with_template(tmp_path, template_charge=5, template_spin=2)
    monkeypatch.setattr(
        u,
        "_derive_charge_from_ligand_charge",
        lambda *_args, **_kwargs: None,
    )

    charge, spin = u.resolve_charge_spin(
        prepared,
        charge=None,
        spin=None,
        ligand_charge="LIG:-2",
    )

    assert charge == 5
    assert spin == 2


def test_resolve_charge_spin_skips_ligand_validation_when_charge_is_explicit(
    tmp_path: Path,
) -> None:
    from pdb2reaction.core import utils as u

    prepared = _prepared_with_template(
        tmp_path,
        source_name="input.gjf",
        template_charge=9,
        template_spin=2,
    )

    charge, spin = u.resolve_charge_spin(
        prepared,
        charge=1,
        spin=None,
        ligand_charge="LIG:-2",
    )

    assert charge == 1
    assert spin == 2


@pytest.mark.parametrize(
    ("headers", "message"),
    [
        (((0, 1), (2, 1)), "Conflicting GJF charge headers"),
        (((0, 1), (0, 3)), "Conflicting GJF multiplicity headers"),
    ],
)
def test_resolve_charge_spin_rejects_conflicting_gjf_headers(
    tmp_path: Path,
    headers: tuple[tuple[int, int], tuple[int, int]],
    message: str,
) -> None:
    from click import ClickException
    from pdb2reaction.core.utils import resolve_charge_spin

    prepared = []
    for index, (charge, spin) in enumerate(headers):
        directory = tmp_path / str(index)
        directory.mkdir()
        item = _prepared_with_template(
            directory,
            source_name=f"image-{index}.gjf",
            template_charge=charge,
            template_spin=spin,
        )
        item.gjf_template.path = directory / f"image-{index}.gjf"
        prepared.append(item)

    with pytest.raises(ClickException, match=message):
        resolve_charge_spin(prepared, charge=None, spin=None)


def test_explicit_gjf_fields_override_only_their_matching_conflicts(
    tmp_path: Path,
) -> None:
    from pdb2reaction.core.utils import resolve_charge_spin

    prepared = []
    for index, charge in enumerate((0, 2)):
        directory = tmp_path / str(index)
        directory.mkdir()
        prepared.append(
            _prepared_with_template(
                directory,
                source_name=f"image-{index}.gjf",
                template_charge=charge,
                template_spin=1,
            )
        )

    assert resolve_charge_spin(prepared, charge=-1, spin=None) == (-1, 1)


def test_resolve_configured_charge_spin_precedence() -> None:
    from pdb2reaction.core.utils import resolve_configured_charge_spin

    cfg = {"calc": {"charge": -2, "spin": 3}}
    assert resolve_configured_charge_spin(
        cfg, charge=None, spin=None,
    ) == (-2, 3)
    assert resolve_configured_charge_spin(
        cfg, charge=1, spin=2,
    ) == (1, 2)
    # An explicit residue-derived charge request outranks calc.charge while
    # the independent configured multiplicity remains effective.
    assert resolve_configured_charge_spin(
        cfg, charge=None, spin=None, ligand_charge="LIG:-1",
    ) == (None, 3)


def test_prepared_cli_input_uses_ref_pdb_for_xyz_ligand_charge(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """XYZ coordinates may derive -l charge from an atom-matched --ref-pdb."""
    from pdb2reaction.core import utils as u

    xyz = tmp_path / "geom.xyz"
    xyz.write_text("1\nframe\nC 0.0 0.0 0.0\n")
    ref = tmp_path / "model.pdb"
    ref.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "END\n"
    )
    seen: dict[str, Path] = {}

    def fake_derive(prepared, ligand_charge, *, prefix):
        seen["source_path"] = prepared.source_path
        assert ligand_charge == "LIG:-2"
        return -2

    monkeypatch.setattr(u, "_derive_charge_from_ligand_charge", fake_derive)

    with u.prepared_cli_input(
        xyz,
        ref_pdb=ref,
        charge=None,
        spin=None,
        ligand_charge="LIG:-2",
    ) as (prepared, charge, spin):
        assert prepared.geom_path == xyz
        assert prepared.source_path == ref.resolve()
        assert charge == -2
        assert spin == 1

    assert seen["source_path"] == ref.resolve()


def test_validate_charge_spin_water_neutral_singlet_passes() -> None:
    from pdb2reaction.core.utils import validate_charge_spin

    validate_charge_spin(["O", "H", "H"], charge=0, multiplicity=1)


def test_validate_charge_spin_water_neutral_doublet_raises() -> None:
    import pytest

    from pdb2reaction.core.utils import validate_charge_spin

    with pytest.raises(ValueError, match="electron count inconsistent"):
        validate_charge_spin(["O", "H", "H"], charge=0, multiplicity=2)


@pytest.mark.parametrize("multiplicity", [0, -2])
def test_validate_charge_spin_rejects_nonpositive_multiplicity(
    multiplicity: int,
) -> None:
    from pdb2reaction.core.utils import (
        set_allow_charge_mult_mismatch,
        validate_charge_spin,
    )

    set_allow_charge_mult_mismatch(True)
    with pytest.raises(ValueError, match="multiplicity must be an integer >= 1"):
        validate_charge_spin(["H"], charge=0, multiplicity=multiplicity)


def test_multiframe_gjf_conversion_warns_that_output_is_coordinate_archive(
    caplog,
    tmp_path: Path,
) -> None:
    """Do not let agents/users mistake the historical multi-frame layout for QST."""
    import logging

    from pdb2reaction.core.utils import convert_xyz_to_gjf, parse_gjf_template

    gjf = tmp_path / "template.gjf"
    gjf.write_text(
        "%chk=test.chk\n"
        "#p hf/sto-3g\n\n"
        "title\n\n"
        "0 1\n"
        "H 0.000000 0.000000 0.000000\n\n"
    )
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text(
        "1\nframe 1\nH 0.0 0.0 0.0\n"
        "1\nframe 2\nH 0.1 0.0 0.0\n"
    )
    out = tmp_path / "trajectory.gjf"

    with caplog.at_level(logging.WARNING):
        convert_xyz_to_gjf(xyz, parse_gjf_template(gjf), out)

    assert "not directly executable as a Gaussian QST/Link1 job" in caplog.text
    assert out.exists()


def test_read_xyz_energies_rejects_bare_integer_frame_comment(
    tmp_path: Path,
) -> None:
    from pdb2reaction.core.utils import read_xyz_energies

    trajectory = tmp_path / "frames.xyz"
    trajectory.write_text("1\n0\nH 0.0 0.0 0.0\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Ambiguous integer-only"):
        read_xyz_energies(trajectory)


def test_read_xyz_energies_accepts_keyed_integer_and_decimal_values(
    tmp_path: Path,
) -> None:
    from pdb2reaction.core.utils import read_xyz_energies

    trajectory = tmp_path / "energies.xyz"
    trajectory.write_text(
        "1\nE=-100\nH 0.0 0.0 0.0\n"
        "1\nEnergy=-99.500 unit=hartree frame 1\nH 0.1 0.0 0.0\n",
        encoding="utf-8",
    )

    assert read_xyz_energies(trajectory) == [-100.0, -99.5]
