"""Shared YAML-overlay helpers for `path-opt` and `path-search` cli().

Both subcommands accept the same single-structure optimiser YAML schema:

  opt:
    <common keys>     # OPT_BASE_KW keys overlay both LBFGS + RFO
    lbfgs: {...}      # LBFGS-only overrides
    rfo: {...}        # RFO-only overrides
"""

from __future__ import annotations

from typing import Any, Dict, Mapping


def apply_single_opt_yaml_layer(
    layer_cfg: Dict[str, Any],
    *,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    opt_base_kw: Mapping[str, Any],
    deep_update,
    apply_yaml_overrides,
) -> None:
    """Apply single-structure optimiser overrides from one YAML layer.

    Mutates ``lbfgs_cfg`` and ``rfo_cfg`` in place. The ``deep_update`` and
    ``apply_yaml_overrides`` callables are taken as parameters so this
    module stays free of upward imports from ``pdb2reaction.core.utils``.
    """
    if not isinstance(layer_cfg, dict):
        return
    opt_section = layer_cfg.get("opt")
    if isinstance(opt_section, dict):
        common_updates = {k: v for k, v in opt_section.items() if k in opt_base_kw}
        if common_updates:
            deep_update(lbfgs_cfg, common_updates)
            deep_update(rfo_cfg, common_updates)

        lbfgs_section = opt_section.get("lbfgs")
        if isinstance(lbfgs_section, dict):
            deep_update(lbfgs_cfg, lbfgs_section)
        else:
            apply_yaml_overrides(
                layer_cfg,
                [(lbfgs_cfg, (("lbfgs",),))],
            )

        rfo_section = opt_section.get("rfo")
        if isinstance(rfo_section, dict):
            deep_update(rfo_cfg, rfo_section)
        else:
            apply_yaml_overrides(
                layer_cfg,
                [(rfo_cfg, (("rfo",),))],
            )
        return

    apply_yaml_overrides(
        layer_cfg,
        [
            (lbfgs_cfg, (("lbfgs",),)),
            (rfo_cfg, (("rfo",),)),
        ],
    )


__all__ = ["apply_single_opt_yaml_layer"]
