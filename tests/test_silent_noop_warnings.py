"""Flags and settings that do nothing must say so.

Requested settings that cannot be honored must emit a diagnostic instead of
appearing to have been applied. These output-only guards do not change
behavior, thresholds, or defaults.
"""

import warnings

import click
import pytest

import pdb2reaction.backends as backends_mod
import pdb2reaction.backends.orb as orb_mod
import pdb2reaction.cli.common_options as common_options_mod
import pdb2reaction.core.utils as utils_mod


def test_auto_backend_says_which_one_it_picked(monkeypatch, capsys) -> None:
    # Provenance records the resolved name, so without this line an automatic
    # choice is indistinguishable from an explicit one.
    def fake_import(name):
        if name in {"fairchem.core", "orb_models.forcefield.pretrained"}:
            raise ImportError(name)
        return object()

    monkeypatch.setattr(backends_mod.importlib, "import_module", fake_import)
    monkeypatch.setattr(backends_mod, "_import_cls", lambda backend, key: object())

    assert backends_mod.resolve_backend("auto") == "mace"

    assert "backend='auto' resolved to 'mace'" in capsys.readouterr().err


def test_dropped_calc_keys_are_reported(monkeypatch, capsys) -> None:
    # `_filter_kwargs` silently drops calc: keys a backend does not accept.
    built = []

    class _Calculator:
        def __init__(self, **kwargs):
            built.append(kwargs)

    monkeypatch.setattr(backends_mod, "_import_cls", lambda backend, key: _Calculator)

    backends_mod.create_calculator("orb", nonexistent_knob=1, radius=6.0)

    assert "nonexistent_knob" not in built[0]
    err = capsys.readouterr().err
    assert "orb ignored calc setting(s) nonexistent_knob;" in err


def test_unusable_calc_keys_are_reported_but_seeded_defaults_are_not() -> None:
    """Every workflow seeds `calc_cfg` from the UMA-flavoured block, so a key UMA
    accepts and this backend does not is an untouched default, not a choice.
    Reporting those buried the line under seven or eight settings per run. A key
    outside UMA's own accepted set is genuinely unusable and still reported.
    """
    accepted_map = backends_mod._BACKEND_ACCEPTED_KEYS
    seeded = accepted_map["uma"]

    def dropped(backend: str, kwargs: dict) -> list:
        filtered = backends_mod._filter_kwargs(kwargs, accepted_map[backend])
        return sorted(k for k in kwargs if k not in filtered and k not in seeded)

    # The measured seeded noise: UMA accepts these, the others do not.
    noise = {
        "max_neigh": 30,
        "print_vram": True,
        "r_edges": True,
        "radius": 6.0,
        "task_name": "omol",
        "workers": 1,
        "workers_per_node": 1,
    }
    for backend in ("orb", "mace", "aimnet2"):
        assert set(noise) - accepted_map[backend], f"{backend} accepts all of {noise}"
        assert dropped(backend, dict(noise)) == []
        # A key no backend knows is still surfaced.
        assert dropped(backend, {**noise, "nonexistent_knob": 1}) == ["nonexistent_knob"]


def test_non_uma_worker_parallelism_is_reported_and_removed() -> None:
    from pdb2reaction.core.defaults import apply_backend_defaults

    cfg = {"backend": "orb", "workers": 4, "workers_per_node": 2}
    with pytest.warns(UserWarning, match="does not use UMA worker parallelism"):
        apply_backend_defaults(cfg)

    assert "workers" not in cfg
    assert "workers_per_node" not in cfg


def test_non_uma_default_worker_values_are_silent() -> None:
    from pdb2reaction.core.defaults import apply_backend_defaults

    cfg = {"backend": "orb", "workers": 1, "workers_per_node": 1}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        apply_backend_defaults(cfg)

    assert "workers" not in cfg
    assert "workers_per_node" not in cfg


@pytest.mark.parametrize(("compile_model", "reported"), [(True, True), (None, False)])
def test_orb_reports_an_unhonoured_compile_request(capsys, compile_model, reported) -> None:
    # The second loader base drops `compile`; reaching it means
    # calc.compile_model=true was not honoured.
    def loader(**kwargs):
        if "compile" in kwargs:
            raise TypeError("unexpected keyword argument 'compile'")
        return "model", "adapter"

    calculator = object.__new__(orb_mod.OrbCalculator)
    calculator.device_str = "cpu"
    calculator.precision = "float32-highest"
    calculator.compile_model = compile_model
    calculator._loader = loader

    assert calculator._load_model() == ("model", "adapter")

    err = capsys.readouterr().err
    assert ("does not accept `compile`" in err) is reported
    assert ("calc.compile_model=true" in err) is reported


def test_freeze_links_reports_having_no_cap_hydrogen(tmp_path, capsys) -> None:
    # --freeze-links defaults on; with no LKH/HL the run freezes nothing.
    structure = tmp_path / "plain.pdb"
    structure.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    assert utils_mod.detect_freeze_links(structure) == []

    err = capsys.readouterr().err
    assert "no LKH/HL cap hydrogen found" in err
    assert "no cap parent is frozen despite --freeze-links" in err


def test_no_deterministic_cannot_override_the_environment(capsys) -> None:
    """`--no-deterministic` does not switch off what the env var turned on."""
    import os

    from pdb2reaction.cli.common_options import _deterministic_callback

    class _Ctx:
        resilient_parsing = False

        @staticmethod
        def get_parameter_source(_name):
            return click.core.ParameterSource.COMMANDLINE

    class _Param:
        name = "deterministic"

    prev = os.environ.get("PDB2REACTION_STRICT_DETERMINISTIC")
    os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = "1"
    try:
        _deterministic_callback(_Ctx(), _Param(), False)
        err = capsys.readouterr().err
        assert "stays strictly deterministic despite" in err
    finally:
        if prev is None:
            os.environ.pop("PDB2REACTION_STRICT_DETERMINISTIC", None)
        else:
            os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = prev


def test_no_deterministic_is_silent_without_the_environment(capsys) -> None:
    import os

    from pdb2reaction.cli.common_options import _deterministic_callback

    class _Ctx:
        resilient_parsing = False

    prev = os.environ.pop("PDB2REACTION_STRICT_DETERMINISTIC", None)
    try:
        _deterministic_callback(_Ctx(), None, False)
        assert capsys.readouterr().err == ""
    finally:
        if prev is not None:
            os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = prev
