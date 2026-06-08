"""Unit tests for pdb2reaction.utils helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np


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


def test_validate_charge_spin_water_neutral_singlet_passes() -> None:
    from pdb2reaction.core.utils import validate_charge_spin

    validate_charge_spin(["O", "H", "H"], charge=0, multiplicity=1)


def test_validate_charge_spin_water_neutral_doublet_raises() -> None:
    import pytest

    from pdb2reaction.core.utils import validate_charge_spin

    with pytest.raises(ValueError, match="electron count inconsistent"):
        validate_charge_spin(["O", "H", "H"], charge=0, multiplicity=2)
