"""Unit tests for pdb2reaction.domain.bond_changes helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner
import numpy as np
import pytest



class TestBondStr:
    def test_one_based(self):
        from pdb2reaction.domain.bond_changes import _bond_str

        result = _bond_str(0, 1, ["C", "O"], one_based=True)
        assert result == "C1-O2"

    def test_zero_based(self):
        from pdb2reaction.domain.bond_changes import _bond_str

        result = _bond_str(0, 1, ["C", "O"], one_based=False)
        assert result == "C0-O1"

    def test_same_element(self):
        from pdb2reaction.domain.bond_changes import _bond_str

        result = _bond_str(2, 5, ["H", "H", "C", "N", "O", "C"], one_based=True)
        assert result == "C3-C6"


class TestElementArrays:
    def test_basic(self):
        from pdb2reaction.domain.bond_changes import _element_arrays

        elems, radii = _element_arrays(["C", "O", "N"])
        assert elems == ["C", "O", "N"]
        assert len(radii) == 3
        # Carbon covalent radius should be positive
        assert radii[0] > 0

    def test_lowercase_normalization(self):
        from pdb2reaction.domain.bond_changes import _element_arrays

        elems, radii = _element_arrays(["c", "o"])
        assert elems == ["C", "O"]

    def test_empty(self):
        from pdb2reaction.domain.bond_changes import _element_arrays

        elems, radii = _element_arrays([])
        assert elems == []
        assert len(radii) == 0


class TestResolveDevice:
    def test_cpu(self):
        import torch
        from pdb2reaction.domain.bond_changes import _resolve_device

        dev = _resolve_device("cpu")
        assert dev == torch.device("cpu")

    def test_auto_returns_device(self):
        import torch
        from pdb2reaction.domain.bond_changes import _resolve_device

        dev = _resolve_device("auto")
        assert isinstance(dev, torch.device)


# Bohr; well inside the C-C covalent cutoff (~3.27 Bohr after the default margin)
_BOND_BOHR = 2.8


def _line_geom(coords):
    """Minimal geom stub: all-carbon chain with explicit Bohr coordinates."""
    coords = np.asarray(coords, dtype=float)
    return SimpleNamespace(atoms=["C"] * len(coords), coords3d=coords)


def _as_int_pairs(pairs):
    # compare_structures returns numpy-int tuples; normalize for comparison
    return {tuple(int(x) for x in p) for p in pairs}


class TestCompareStructures:
    @pytest.mark.parametrize("bond_factor", [0.0, -1.0, np.nan, np.inf])
    def test_rejects_invalid_bond_factor(self, bond_factor):
        from pdb2reaction.domain.bond_changes import compare_structures

        geom = _line_geom([[0.0, 0.0, 0.0]])
        with pytest.raises(ValueError, match="bond_factor"):
            compare_structures(
                geom,
                geom,
                device="cpu",
                bond_factor=bond_factor,
            )

    @pytest.mark.parametrize("nonfinite", [np.nan, np.inf, -np.inf])
    def test_rejects_nonfinite_coordinates(self, nonfinite):
        from pdb2reaction.domain.bond_changes import compare_structures

        finite = _line_geom([[0.0, 0.0, 0.0]])
        invalid = _line_geom([[nonfinite, 0.0, 0.0]])
        with pytest.raises(ValueError, match="coordinates must be finite"):
            compare_structures(finite, invalid, device="cpu")

    def test_chunk_boundary_formed_and_broken(self):
        """Row-chunked detection must be correct across the 1024-row block
        boundary: one bond forms in the first block, one breaks in the second
        (exercises the per-block upper-triangle offset and the ``+= i0`` global
        index correction)."""
        from pdb2reaction.domain.bond_changes import compare_structures

        N = 1026
        base = np.zeros((N, 3), dtype=float)
        base[:, 0] = np.arange(N) * 20.0  # 20 Bohr spacing -> no bonds anywhere

        r1 = base.copy()
        r2 = base.copy()
        # formed (2, 3): apart in r1, bonded in r2  -> first row block
        r2[3, 0] = base[2, 0] + _BOND_BOHR
        # broken (1024, 1025): bonded in r1, apart in r2  -> second row block
        r1[1025, 0] = base[1024, 0] + _BOND_BOHR

        res = compare_structures(_line_geom(r1), _line_geom(r2), device="cpu")

        assert _as_int_pairs(res.formed_covalent) == {(2, 3)}
        assert _as_int_pairs(res.broken_covalent) == {(1024, 1025)}
        # memory-bounded path: full N x N matrices dropped, sparse lengths kept
        assert res.distances_1 is None and res.distances_2 is None
        lengths = {tuple(int(x) for x in k): v for k, v in res.changed_lengths.items()}
        assert lengths[(2, 3)][1] == pytest.approx(_BOND_BOHR, abs=1e-6)      # d2 (formed)
        assert lengths[(1024, 1025)][0] == pytest.approx(_BOND_BOHR, abs=1e-6)  # d1 (broken)

    def test_summarize_uses_sparse_lengths(self):
        """summarize_changes still prints the ``D1 Å --> D2 Å`` lengths from the
        sparse map even though the full distance matrices are gone."""
        from pdb2reaction.domain.bond_changes import compare_structures, summarize_changes

        r1 = [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]            # apart
        r2 = [[0.0, 0.0, 0.0], [_BOND_BOHR, 0.0, 0.0]]      # bonded -> formed
        res = compare_structures(_line_geom(r1), _line_geom(r2), device="cpu")

        assert _as_int_pairs(res.formed_covalent) == {(0, 1)}
        text = summarize_changes(_line_geom(r2), res, one_based=True)
        assert "C1-C2" in text
        assert "Å -->" in text


def _write_ordinary_gjf(path: Path, *, product: bool) -> None:
    fluorine_x = 1.35 if product else 4.0
    chlorine_x = -4.0 if product else -1.8
    path.write_text(
        "%chk=test.chk\n"
        "# hf/sto-3g\n\n"
        "ordinary Gaussian input\n\n"
        "0 1\n"
        "C   0.000000  0.000000  0.000000\n"
        "H   0.000000  1.000000  0.000000\n"
        "H   0.000000 -1.000000  0.000000\n"
        "H   0.000000  0.000000  1.000000\n"
        f"F   {fluorine_x:.6f}  0.000000  0.000000\n"
        f"Cl  {chlorine_x:.6f}  0.000000  0.000000\n\n",
        encoding="utf-8",
    )


def test_bond_summary_loads_ordinary_gjf_through_p2r_preparation(
    tmp_path: Path,
) -> None:
    from pysisyphus.constants import BOHR2ANG
    from pdb2reaction.domain.bond_summary import _load_geom

    reactant = tmp_path / "reactant.gjf"
    _write_ordinary_gjf(reactant, product=False)

    geom = _load_geom(str(reactant))
    assert tuple(geom.atoms) == ("C", "H", "H", "H", "F", "Cl")
    np.testing.assert_allclose(
        np.asarray(geom.coords3d)[:, 0] * BOHR2ANG,
        [0.0, 0.0, 0.0, 0.0, 4.0, -1.8],
        atol=1.0e-10,
    )


def test_real_bond_summary_cli_reports_ordinary_gjf_bond_changes(
    tmp_path: Path,
) -> None:
    from pdb2reaction.cli import cli as root_cli

    reactant = tmp_path / "reactant.gjf"
    product = tmp_path / "product.gjf"
    _write_ordinary_gjf(reactant, product=False)
    _write_ordinary_gjf(product, product=True)

    result = CliRunner().invoke(
        root_cli,
        ["bond-summary", str(reactant), str(product)],
    )

    assert result.exit_code == 0, result.output
    assert "C1-F5" in result.output
    assert "C1-Cl6" in result.output


def test_metal_pairs_need_a_larger_change():
    """Metal coordination fluctuations are not bond changes; large moves are."""
    from pysisyphus.constants import BOHR2ANG
    from pdb2reaction.core.defaults import BOND_KW
    from pdb2reaction.domain.bond_changes import compare_structures

    def _pair(elems, d_ang):
        coords = np.array([[0.0, 0.0, 0.0], [d_ang / BOHR2ANG, 0.0, 0.0]])
        return SimpleNamespace(atoms=list(elems), coords3d=coords)

    mgo = ("Mg", "O")
    assert BOND_KW["metal_delta_fraction"] == 0.40
    # The Mg-O boundary is 2.36 A; 2.1 -> 2.6 A crosses it with a 0.5 A change.
    res = compare_structures(_pair(mgo, 2.1), _pair(mgo, 2.6), device="cpu")
    assert not res.formed_covalent and not res.broken_covalent
    res = compare_structures(
        _pair(mgo, 2.1), _pair(mgo, 2.6), device="cpu", metal_delta_fraction=0.05
    )
    assert _as_int_pairs(res.broken_covalent) == {(0, 1)}
    res = compare_structures(_pair(mgo, 2.1), _pair(mgo, 3.3), device="cpu")
    assert _as_int_pairs(res.broken_covalent) == {(0, 1)}
    # Pairs without a metal keep delta_fraction: C-O 1.43 -> 1.80 A still breaks.
    co = ("C", "O")
    res = compare_structures(_pair(co, 1.43), _pair(co, 1.80), device="cpu")
    assert _as_int_pairs(res.broken_covalent) == {(0, 1)}
    from pdb2reaction.domain.bond_changes import has_bond_change

    changed, _ = has_bond_change(_pair(mgo, 2.1), _pair(mgo, 2.6), dict(BOND_KW, device="cpu"))
    assert changed is False



def _write_mgo_xyz(path: Path, distance: float) -> Path:
    path.write_text(f"2\nMgO\nMg 0 0 0\nO {distance} 0 0\n", encoding="utf-8")
    return path


@pytest.mark.parametrize(
    "extra,broken",
    [
        pytest.param([], False, id="default-0.40"),
        pytest.param(["--metal-delta-fraction", "0.05"], True, id="metal-delta-0.05"),
    ],
)
def test_bond_summary_metal_delta_fraction_option(tmp_path: Path, extra, broken) -> None:
    from pdb2reaction.cli import cli as root_cli

    # Mg-O 2.1 -> 2.6 A crosses the bond threshold with a 0.5 A change.
    first = _write_mgo_xyz(tmp_path / "a.xyz", 2.1)
    second = _write_mgo_xyz(tmp_path / "b.xyz", 2.6)

    result = CliRunner().invoke(root_cli, ["bond-summary", str(first), str(second), *extra])

    assert result.exit_code == 0, result.output
    assert ("Mg1-O2" in result.output) is broken


@pytest.mark.parametrize(
    "bond_yaml,broken",
    [
        pytest.param(None, [], id="default-0.40"),
        pytest.param("bond:\n  metal_delta_fraction: 0.05\n", ["Mg1-O2"], id="yaml-0.05"),
    ],
)
def test_irc_endpoint_bond_changes_read_the_yaml_bond_section(
    tmp_path: Path, monkeypatch, bond_yaml, broken
) -> None:
    import json

    import pdb2reaction.io.hessian_cache as hessian_cache
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.workflows import irc as irc_workflow

    class EndpointsEulerPC:
        """Writes Mg-O 2.1 A and 2.6 A as the two IRC endpoints."""

        def __init__(self, geometry, **kwargs):
            self.out_dir = Path(kwargs["out_dir"])
            self._act_atoms = [0, 1]
            self.max_cycles = 50
            self.never_stop = False
            self.all_energies = [-1.3, -1.0, -1.3]
            self.ts_energy = -1.0
            self.init_hessian_shape = (6, 6)
            for direction in ("forward", "backward"):
                setattr(self, direction, True)
                setattr(self, f"{direction}_energies", [-1.0, -1.1, -1.2, -1.3])
                setattr(self, f"{direction}_coords", [np.zeros(6)] * 4)
                setattr(self, f"{direction}_cycle", 3)

        def get_path_for_fn(self, name):
            return self.out_dir / name

        def run(self):
            self.out_dir.mkdir(parents=True, exist_ok=True)
            _write_mgo_xyz(self.out_dir / "finished_first.xyz", 2.1)
            _write_mgo_xyz(self.out_dir / "finished_last.xyz", 2.6)

    class Calculator:
        def get_hessian(self, atoms, coords):
            return {"energy": -1.0, "hessian": np.eye(len(coords))}

    monkeypatch.setattr(irc_workflow, "create_calculator", lambda **_kw: Calculator())
    monkeypatch.setattr(irc_workflow, "EulerPC", EndpointsEulerPC)
    monkeypatch.setattr(hessian_cache, "load_matching", lambda *args, **kwargs: None)
    monkeypatch.setattr(hessian_cache, "identity_from_context", lambda *args, **kwargs: {})
    args = ["irc", "-i", str(_write_mgo_xyz(tmp_path / "ts.xyz", 2.3)), "-q", "0", "-m", "1",
            "-o", str(tmp_path / "out"), "--out-json"]
    if bond_yaml is not None:
        config = tmp_path / "bond.yaml"
        config.write_text(bond_yaml, encoding="utf-8")
        args += ["--config", str(config)]

    result = CliRunner().invoke(root_cli, args)

    assert result.exit_code == 0, result.output
    data = json.loads((tmp_path / "out" / "result.json").read_text(encoding="utf-8"))
    assert data["bond_changes"] == {"formed": [], "broken": broken}
