from __future__ import annotations

import numpy as np
import pytest

from pdb2reaction.backends import xtb_alpb_correction as xtb


def test_xtb_engrad_stdout_and_unit_conversions_accept_d_exponents(tmp_path) -> None:
    engrad = tmp_path / "sample.engrad"
    engrad.write_text(
        "# The current total energy in Eh\n"
        "-1.250000000000D+00\n"
        "# The current gradient in Eh/bohr\n"
        "1.0D-02\n-2.0D-02\n3.0D-02\n"
        "-4.0D-02\n5.0D-02\n-6.0D-02\n",
        encoding="utf-8",
    )
    energy, gradient = xtb._parse_engrad(engrad, 2)
    assert energy == pytest.approx(-1.25)
    np.testing.assert_allclose(
        gradient,
        [[0.01, -0.02, 0.03], [-0.04, 0.05, -0.06]],
    )
    assert xtb._parse_energy_from_stdout(
        "cycle\n TOTAL ENERGY       -2.500000D+00 Eh\n",
    ) == pytest.approx(-2.5)

    energy_ev, forces, hessian = xtb.convert_units_xtb_to_mlip(
        energy_ha=energy,
        gradient_ha_bohr=gradient,
        hessian_ha_bohr2=np.eye(6),
    )
    assert energy_ev == pytest.approx(energy * xtb.EV_PER_HARTREE)
    np.testing.assert_allclose(
        forces, -gradient * xtb.FORCE_EV_ANG_PER_HA_BOHR,
    )
    np.testing.assert_allclose(
        hessian, np.eye(6) * xtb.HESS_EV_ANG2_PER_HA_BOHR2,
    )


def test_xtb_engrad_rejects_truncated_gradient(tmp_path) -> None:
    engrad = tmp_path / "truncated.engrad"
    engrad.write_text(
        "# The current total energy in Eh\n-1.0\n"
        "# The current gradient in Eh/bohr\n0.1\n0.2\n",
        encoding="utf-8",
    )
    with pytest.raises(xtb.XTBError, match="expected 3, got 2"):
        xtb._parse_engrad(engrad, 1)


@pytest.mark.parametrize("blocked", [False, True])
def test_xtb_hessian_parses_dense_and_blocked_literals(tmp_path, blocked) -> None:
    path = tmp_path / "hessian"
    matrix = np.array(
        [[1.0, 0.2, 0.3], [0.2, 2.0, 0.4], [0.3, 0.4, 3.0]],
    )
    if blocked:
        body = (
            "$hessian\n3\n1 2 3\n"
            "1 1.0D+00 2.0D-01 3.0D-01\n"
            "2 2.0D-01 2.0D+00 4.0D-01\n"
            "3 3.0D-01 4.0D-01 3.0D+00\n$end\n"
        )
    else:
        body = (
            "$hessian\n"
            "1.0D+00 2.0D-01 3.0D-01\n"
            "2.0D-01 2.0D+00 4.0D-01\n"
            "3.0D-01 4.0D-01 3.0D+00\n$end\n"
        )
    path.write_text(body, encoding="utf-8")
    np.testing.assert_allclose(xtb._parse_xtb_hessian(path, 1), matrix)


def test_xtb_dense_hessian_rejects_truncation(tmp_path) -> None:
    path = tmp_path / "hessian"
    path.write_text("$hessian\n1.0 0.0\n$end\n", encoding="utf-8")
    with pytest.raises(xtb.XTBError, match="expected at least 9 values, got 2"):
        xtb._parse_xtb_hessian(path, 1)


def test_xtb_raw_hessian_requires_exact_finite_square(tmp_path) -> None:
    path = tmp_path / "hessian.out"
    matrix = np.arange(9, dtype=float).reshape(3, 3)
    path.write_text("\n".join(" ".join(map(str, row)) for row in matrix))
    np.testing.assert_allclose(
        xtb._parse_xtb_hessian(path, 1),
        0.5 * (matrix + matrix.T),
    )

    path.write_text("0 1 2 3 4 5 6 7")
    with pytest.raises(xtb.XTBError, match="expected exactly 9 finite values"):
        xtb._parse_xtb_hessian(path, 1)


def test_xtb_command_absolutizes_only_path_like_executable(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    assert xtb._xtb_cmd_tokens("xtb --copy") == ["xtb", "--copy"]
    assert xtb._xtb_cmd_tokens("./bin/xtb --copy")[0] == str(
        (tmp_path / "bin" / "xtb").resolve()
    )

    cmd = xtb._build_xtb_cmd(
        xtb_cmd="xtb",
        xyz_filename="input.xyz",
        charge=0,
        multiplicity=1,
        solvent="none",
        solvent_model="alpb",
        xtb_acc=0.2,
        mode="hess",
        xcontrol_filename="raw_hessian.inp",
    )
    assert cmd[cmd.index("--input") + 1] == "raw_hessian.inp"
