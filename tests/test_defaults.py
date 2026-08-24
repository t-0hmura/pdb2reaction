"""Tests for pdb2reaction.core.defaults configuration constants."""

from copy import deepcopy
from inspect import signature

from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    CALC_KW_DEFAULT,
    OPT_BASE_KW,
    LBFGS_KW,
    LBFGS_TS_KW,
    RFO_KW,
    RSIRFO_KW,
    GS_KW,
    SEARCH_KW,
    IRC_KW,
    BOND_KW,
    BIAS_KW,
    THERMO_KW,
    OPT_MODE_ALIASES,
    TSOPT_MODE_ALIASES,
    DMF_KW,
    fresh_dmf_config,
)


class TestDefaultsStructure:
    """Verify that default dictionaries have expected keys and types."""

    def test_geom_kw_has_coord_type(self):
        assert "coord_type" in GEOM_KW_DEFAULT
        assert isinstance(GEOM_KW_DEFAULT["coord_type"], str)

    def test_calc_kw_has_backend(self):
        assert "backend" in CALC_KW_DEFAULT
        assert CALC_KW_DEFAULT["backend"] == "uma"

    def test_calc_kw_has_charge_spin(self):
        assert "charge" in CALC_KW_DEFAULT
        assert "spin" in CALC_KW_DEFAULT

    def test_opt_base_kw_has_thresh(self):
        assert "thresh" in OPT_BASE_KW
        assert isinstance(OPT_BASE_KW["thresh"], str)

    def test_lbfgs_has_max_cycles(self):
        assert "max_cycles" in LBFGS_KW
        assert LBFGS_KW["max_cycles"] == 100000

    def test_rfo_has_max_cycles(self):
        assert "max_cycles" in RFO_KW
        assert RFO_KW["max_cycles"] == 100000

    def test_minimizer_trial_rejection_defaults(self):
        assert LBFGS_KW["reject_uphill"] is False
        assert RFO_KW["reject_uphill"] is False
        assert LBFGS_KW["uphill_tolerance"] == 1e-4
        assert RFO_KW["uphill_tolerance"] == 1e-4
        assert LBFGS_TS_KW["reject_uphill"] is False

        from pysisyphus.optimizers.LBFGS import LBFGS
        from pysisyphus.optimizers.RFOptimizer import RFOptimizer

        assert signature(LBFGS).parameters["reject_uphill"].default is False
        assert signature(RFOptimizer).parameters["reject_uphill"].default is False
        assert signature(LBFGS).parameters["uphill_tolerance"].default == 1e-4
        assert signature(RFOptimizer).parameters["uphill_tolerance"].default == 1e-4

    def test_ts_search_uses_terminal_exact_validation_without_automatic_recovery(self):
        assert RSIRFO_KW["check_eigval_structure"] is False
        assert RSIRFO_KW["reject_mode_loss"] is False
        assert RSIRFO_KW["verify_saddle"] is True
        assert RSIRFO_KW["saddle_imaginary_threshold_cm"] == 5.0
        assert RSIRFO_KW["saddle_recovery_check_interval"] == 50
        assert RSIRFO_KW["saddle_recovery_max_cycles"] == 0

    def test_path_workflow_max_nodes_release_defaults(self):
        from pdb2reaction.workflows.path_search import _stitch_paths

        assert GS_KW["max_nodes"] == 20
        assert SEARCH_KW["max_nodes_segment"] == 20
        assert signature(_stitch_paths).parameters["max_nodes"].default == 20

    def test_irc_kw_has_step_length(self):
        assert "step_length" in IRC_KW
        assert IRC_KW["step_length"] > 0

    def test_bond_kw_has_factor(self):
        assert "bond_factor" in BOND_KW
        assert BOND_KW["bond_factor"] > 0

    def test_bias_kw_has_k(self):
        assert "k" in BIAS_KW
        assert float(BIAS_KW["k"]) > 0

    def test_thermo_kw_structure(self):
        assert "temperature" in THERMO_KW
        assert THERMO_KW["temperature"] > 0

    def test_opt_mode_aliases_cover_expected(self):
        # OPT_MODE_ALIASES is a tuple of ((aliases,...), canonical) pairs
        canonical_modes = {pair[1] for pair in OPT_MODE_ALIASES}
        assert "lbfgs" in canonical_modes
        assert "rfo" in canonical_modes

    def test_tsopt_mode_aliases_cover_expected(self):
        canonical_modes = {pair[1] for pair in TSOPT_MODE_ALIASES}
        assert "dimer" in canonical_modes
        assert "rsirfo" in canonical_modes


class TestDefaultsConsistency:
    """Cross-check defaults for internal consistency."""

    def test_no_negative_cycles(self):
        for name, d in [("LBFGS_KW", LBFGS_KW), ("RFO_KW", RFO_KW), ("OPT_BASE_KW", OPT_BASE_KW)]:
            if "max_cycles" in d:
                assert d["max_cycles"] is None or d["max_cycles"] > 0, (
                    f"{name} has an invalid max_cycles"
                )

    def test_solvent_default_is_none(self):
        assert CALC_KW_DEFAULT.get("solvent", "none") == "none"

    def test_fresh_dmf_config_isolates_nested_request_overrides(self):
        canonical = deepcopy(DMF_KW)

        first = fresh_dmf_config(
            {
                "fbenm_options": {"delta_scale": 9.0},
                "cfbenm_options": {"eps": 8.0},
                "dmf_options": {"beta": 7.0},
            }
        )
        first["fbenm_options"]["bond_scale"] = 6.0

        second = fresh_dmf_config()
        assert first["fbenm_options"]["delta_scale"] == 9.0
        assert first["cfbenm_options"]["eps"] == 8.0
        assert first["dmf_options"]["beta"] == 7.0
        assert second == canonical
        assert DMF_KW == canonical
        assert second["fbenm_options"] is not DMF_KW["fbenm_options"]
        assert second["cfbenm_options"] is not DMF_KW["cfbenm_options"]
        assert second["dmf_options"] is not DMF_KW["dmf_options"]


def test_fresh_dmf_config_rejects_unknown_top_level_keys() -> None:
    """A retained but unconsumed top-level dmf key would keep the default."""
    import click
    import pytest

    from pdb2reaction.core.defaults import DMF_TOP_LEVEL_KEYS, fresh_dmf_config

    # Every consumed top-level key is accepted.
    assert fresh_dmf_config({"max_cycles": 5})["max_cycles"] == 5
    assert fresh_dmf_config({"tol": "baker"})["tol"] == "baker"
    assert fresh_dmf_config({"ipopt_options": {"tol": 1e-6}})["ipopt_options"] == {
        "tol": 1e-6
    }
    # Nested solver maps keep their pass-through contract.
    nested = fresh_dmf_config({"dmf_options": {"beta": 3.0, "future_knob": 1}})
    assert nested["dmf_options"]["beta"] == 3.0
    assert nested["dmf_options"]["future_knob"] == 1

    with pytest.raises(click.BadParameter, match="dmf.max_cycle"):
        fresh_dmf_config({"max_cycle": 5})
    assert "max_cycle" not in DMF_TOP_LEVEL_KEYS
