"""Regression test for scan config precedence.

Guards the C1a fix: ``build_scan_configs`` must apply YAML (--config) over
built-in defaults and then let explicit CLI-derived values override YAML, i.e.
the uniform ``defaults < --config < CLI`` ordering used by the other
subcommands.  Before the fix, ``apply_yaml_overrides`` ran *after* the
CLI-derived assignments, so YAML silently overrode explicit CLI options for
scan / scan2d.
"""

from pdb2reaction.utils import build_scan_configs


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
