# tests/test_defaults.py
"""Tests for pdb2reaction.defaults configuration constants."""

from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    CALC_KW_DEFAULT,
    OPT_BASE_KW,
    LBFGS_KW,
    LBFGS_TS_KW,
    RFO_KW,
    RSIRFO_KW,
    IRC_KW,
    BOND_KW,
    BIAS_KW,
    THERMO_KW,
    OPT_MODE_ALIASES,
    TSOPT_MODE_ALIASES,
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
        assert LBFGS_KW["max_cycles"] > 0

    def test_rfo_has_max_cycles(self):
        assert "max_cycles" in RFO_KW
        assert RFO_KW["max_cycles"] > 0

    def test_minimizer_trial_rejection_defaults(self):
        assert LBFGS_KW["reject_uphill"] is True
        assert RFO_KW["reject_uphill"] is True
        assert LBFGS_TS_KW["reject_uphill"] is False

    def test_ts_saddle_safeguards_are_enabled(self):
        assert RSIRFO_KW["check_eigval_structure"] is True
        assert RSIRFO_KW["reject_mode_loss"] is True
        assert RSIRFO_KW["verify_saddle"] is True
        assert RSIRFO_KW["saddle_imaginary_threshold_cm"] == 5.0

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
                assert d["max_cycles"] > 0, f"{name} has non-positive max_cycles"

    def test_solvent_default_is_none(self):
        assert CALC_KW_DEFAULT.get("solvent", "none") == "none"
