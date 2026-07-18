from __future__ import annotations

import json

import pandas as pd
from click.testing import CliRunner


def test_plot_only_csv_does_not_claim_calculator_provenance(tmp_path, monkeypatch) -> None:
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.workflows import scan3d

    rows = []
    for i, d1 in enumerate((1.0, 2.0)):
        for j, d2 in enumerate((1.5, 2.5)):
            for k, d3 in enumerate((2.0, 3.0)):
                rows.append(
                    {
                        "i": i,
                        "j": j,
                        "k": k,
                        "d1_A": d1,
                        "d2_A": d2,
                        "d3_A": d3,
                        "energy_hartree": float(i + 2 * j + 4 * k) / 1000.0,
                    }
                )
    csv_path = tmp_path / "surface.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    monkeypatch.setattr(scan3d, "_VOLUME_GRID_N", 3)
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(
        root_cli,
        [
            "scan3d",
            "--csv",
            str(csv_path),
            "--out-json",
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["backend"] is None
    assert payload["model"] is None
    assert payload["mlip_backend"] is None
    assert payload["mlip_model"] is None
    assert payload["mlip_precision"] is None
    assert payload["charge"] is None
    assert payload["spin"] is None
    assert payload["solvent"] is None


def test_fresh_scan3d_preserves_calculator_and_solvent_provenance() -> None:
    from pdb2reaction.workflows.scan3d import _result_calculator_fields

    fields = _result_calculator_fields(
        {
            "backend": "orb",
            "model": "orb_v3_conservative_omol",
            "precision": "float64",
            "solvent": "water",
        },
        backend="uma",
        charge=-1,
        spin=2,
        plot_only=False,
    )

    assert fields == {
        "charge": -1,
        "spin": 2,
        "backend": "orb",
        "model": "orb_v3_conservative_omol",
        "mlip_backend": "orb",
        "mlip_model": "orb_v3_conservative_omol",
        "mlip_precision": "fp64",
        "solvent": "water",
    }
