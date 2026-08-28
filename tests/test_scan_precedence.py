"""Regression tests for ``defaults < --config < CLI`` scan precedence."""

from click.testing import CliRunner

from pdb2reaction.core.utils import build_scan_configs, build_sopt_kwargs
from pdb2reaction.workflows import scan3d as scan3d_workflow


def _base_kw():
    return dict(
        geom_kw={},
        calc_kw={"workers": 1, "workers_per_node": 1},
        opt_kw={"thresh": "default"},
        lbfgs_kw={},
        rfo_kw={},
        bias_kw={"k": 100.0},
    )


def test_cli_thresh_overrides_yaml():
    """Explicit CLI --thresh wins over a conflicting YAML opt.thresh."""
    _, _, opt_cfg, _, _, _ = build_scan_configs(
        {"opt": {"thresh": "yamlval"}},
        **_base_kw(),
        thresh="cli_baker",  # CLI explicit (option default is None)
    )
    assert opt_cfg["thresh"] == "cli_baker"


def test_yaml_thresh_wins_when_no_cli():
    """YAML opt.thresh applies over the built-in default when CLI is unset."""
    _, _, opt_cfg, _, _, _ = build_scan_configs(
        {"opt": {"thresh": "yamlval"}},
        **_base_kw(),
        thresh=None,  # CLI not given
    )
    assert opt_cfg["thresh"] == "yamlval"


def test_yaml_workers_wins_when_cli_not_overridden():
    """YAML calc.workers must not be clobbered by the non-explicit CLI default."""
    _, calc_cfg, _, _, _, _ = build_scan_configs(
        {"calc": {"workers": 4}},
        **_base_kw(),
        workers=1,
        workers_overridden=False,  # CLI default, not explicitly set
    )
    assert calc_cfg["workers"] == 4


def test_cli_workers_overrides_yaml_when_explicit():
    """Explicit CLI --workers wins over YAML calc.workers."""
    _, calc_cfg, _, _, _, _ = build_scan_configs(
        {"calc": {"workers": 4}},
        **_base_kw(),
        workers=8,
        workers_overridden=True,  # cli_param_overridden(ctx, "workers")
    )
    assert calc_cfg["workers"] == 8


def test_cli_bias_k_overrides_yaml():
    """Explicit CLI --bias-k wins over YAML bias.k."""
    *_, bias_cfg = build_scan_configs(
        {"bias": {"k": 250.0}},
        **_base_kw(),
        bias_k=999.0,  # CLI explicit (option default is None)
    )
    assert bias_cfg["k"] == 999.0


def test_yaml_relax_cycles_survive_when_cli_value_is_omitted():
    _, _, opt_cfg, _, _, _ = build_scan_configs(
        {"opt": {"max_cycles": 321}},
        **_base_kw(),
        relax_max_cycles=10000,
        relax_max_cycles_overridden=False,
    )
    assert opt_cfg["max_cycles"] == 321


def test_explicit_relax_cycles_override_yaml():
    _, _, opt_cfg, _, _, _ = build_scan_configs(
        {"opt": {"max_cycles": 321}},
        **_base_kw(),
        relax_max_cycles=17,
        relax_max_cycles_overridden=True,
    )
    assert opt_cfg["max_cycles"] == 17


def test_scan_optimizer_builder_preserves_legacy_cycle_override_signature(tmp_path):
    legacy = build_sopt_kwargs(
        "lbfgs",
        {},
        {},
        {"max_cycles": 321},
        1.0,
        17,
        True,
        tmp_path,
        "legacy",
    )
    resolved = build_sopt_kwargs(
        "lbfgs",
        {},
        {},
        {"max_cycles": 321},
        1.0,
        out_dir=tmp_path,
        prefix="resolved",
    )

    assert legacy["max_cycles"] == 17
    assert resolved["max_cycles"] == 321


def test_grid_scan_default_is_baker_without_masking_yaml():
    grid_kw = _base_kw()
    grid_kw["opt_kw"] = {"thresh": "baker"}

    _, _, default_opt, _, _, _ = build_scan_configs({}, **grid_kw, thresh=None)
    _, _, yaml_opt, _, _, _ = build_scan_configs(
        {"opt": {"thresh": "gau_tight"}}, **grid_kw, thresh=None
    )

    assert default_opt["thresh"] == "baker"
    assert yaml_opt["thresh"] == "gau_tight"


def test_scan3d_forwards_actual_click_parameter_sources(monkeypatch, tmp_path):
    csv_path = tmp_path / "surface.csv"
    csv_path.write_text("energy_hartree\n0.0\n", encoding="utf-8")
    captured = []

    def _capture_and_stop(**kwargs):
        captured.append(kwargs)
        raise RuntimeError("stop after effective-config capture")

    monkeypatch.setattr(scan3d_workflow, "_build_scan_context", _capture_and_stop)
    runner = CliRunner()

    omitted = runner.invoke(
        scan3d_workflow.cli,
        ["--csv", str(csv_path), "--out-dir", str(tmp_path / "omitted")],
    )
    assert omitted.exit_code == 1
    assert len(captured) == 1
    assert captured[0]["workers_overridden"] is False
    assert captured[0]["workers_per_node_overridden"] is False
    assert captured[0]["relax_max_cycles_overridden"] is False
    assert captured[0]["thresh"] is None

    captured.clear()
    explicit = runner.invoke(
        scan3d_workflow.cli,
        [
            "--csv",
            str(csv_path),
            "--out-dir",
            str(tmp_path / "explicit"),
            "--workers",
            "3",
            "--workers-per-node",
            "2",
            "--relax-max-cycles",
            "19",
            "--thresh",
            "gau_tight",
        ],
    )
    assert explicit.exit_code == 1
    assert len(captured) == 1
    assert captured[0]["workers_overridden"] is True
    assert captured[0]["workers_per_node_overridden"] is True
    assert captured[0]["relax_max_cycles_overridden"] is True
    assert captured[0]["relax_max_cycles"] == 19
    assert captured[0]["thresh"] == "gau_tight"


def _path_search_endpoints(tmp_path):
    reactant = tmp_path / "r.xyz"
    reactant.write_text("2\nR\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", encoding="utf-8")
    product = tmp_path / "p.xyz"
    product.write_text("2\nP\nH 0.0 0.0 0.0\nH 0.0 0.0 1.10\n", encoding="utf-8")
    return reactant, product


def _path_search_dry_run(tmp_path, label, extra):
    from pdb2reaction.workflows import path_search as path_search_workflow

    reactant, product = _path_search_endpoints(tmp_path)
    config = tmp_path / "config.yaml"
    config.write_text("search:\n  refine_mode: minima\n", encoding="utf-8")
    return CliRunner().invoke(
        path_search_workflow.cli,
        [
            "-i", str(reactant),
            "-i", str(product),
            "-q", "0",
            "--config", str(config),
            "--mep-mode", "dmf",
            "--dry-run",
            "-o", str(tmp_path / label),
        ]
        + extra,
    )


def test_explicit_cli_refine_mode_outranks_yaml_search_section(tmp_path):
    """An explicit --refine-mode must not lose to YAML search.refine_mode."""
    omitted = _path_search_dry_run(tmp_path, "omitted", [])
    assert omitted.exit_code == 0
    assert "refine_mode: minima" in omitted.output

    explicit = _path_search_dry_run(tmp_path, "explicit", ["--refine-mode", "peak"])
    assert explicit.exit_code == 0
    assert "refine_mode: peak" in explicit.output
    assert "refine_mode: minima" not in explicit.output


def test_path_search_dry_run_shows_effective_calculator_settings(tmp_path):
    """Requested blocks render at default verbosity with resolved calc values."""
    result = _path_search_dry_run(tmp_path, "display", ["--precision", "fp64"])

    assert result.exit_code == 0
    # The block itself is printed without raising the verbosity ...
    assert "dry_run_plan" in result.output
    assert "calc" in result.output
    # ... and it carries the effective CLI calculator setting.
    assert "precision: fp64" in result.output


def test_path_search_omitted_reference_merge_is_silent(tmp_path):
    template = tmp_path / "template.pdb"
    template.write_text("END\n", encoding="utf-8")

    result = _path_search_dry_run(
        tmp_path,
        "no-reference-merge",
        ["--ref-full-pdb", str(template)],
    )

    assert result.exit_code == 0, result.output
    assert "[merge]" not in result.output


def test_path_search_explicit_reference_merge_missing_template_warns(tmp_path):
    result = _path_search_dry_run(
        tmp_path,
        "missing-reference-template",
        ["--write-ref-merge"],
    )

    assert result.exit_code == 0, result.output
    assert (
        "[merge] WARNING: --write-ref-merge requires --ref-full-pdb; skipped."
        in result.output
    )


def test_path_search_dmf_solvent_is_rejected_before_dry_run(tmp_path):
    result = _path_search_dry_run(
        tmp_path,
        "dmf-solvent",
        ["--solvent", "water"],
    )

    assert result.exit_code != 0
    assert "not compatible with --solvent 'water'" in result.output


def test_path_search_yaml_preserves_stopt_cycle_pair_and_output(tmp_path):
    from pdb2reaction.workflows import path_search as path_search_workflow

    reactant, product = _path_search_endpoints(tmp_path)
    configured_out = tmp_path / "yaml-output"
    config = tmp_path / "stopt.yaml"
    config.write_text(
        "stopt:\n"
        f"  out_dir: {configured_out}\n"
        "  max_cycles: 31\n"
        "  stop_in_when_full: 7\n"
    )

    result = CliRunner().invoke(
        path_search_workflow.cli,
        [
            "-i", str(reactant), "-i", str(product), "-q", "0",
            "--config", str(config), "--mep-mode", "dmf", "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "stop_in_when_full: 7" in result.output
    assert "max_cycles: 31" in result.output
    assert f"out_dir: {configured_out}" in result.output
