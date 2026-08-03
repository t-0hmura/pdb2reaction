"""Helpers extracted from `pdb2reaction/workflows/all.py:cli()`.

The `all` subcommand's `cli()` function is ~4,900 LOC; progressive
helper extracts land in this module so the cli() body shrinks one piece
at a time while behavior is preserved.

Anything imported here must be safe to use from ``cli()`` body callers
without changing observable behavior.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple


def build_energy_level_dict(
    *,
    labels: Sequence[str],
    energies_au: Sequence[float],
    ref_energy: float,
    au_to_kcal: float,
    diagram_path: str,
    structures: Dict[str, Any],
) -> Dict[str, Any]:
    """Assemble a segment_log sub-dict for one energy level (MLIP / Gibbs / DFT).

    Extracted from the 4-way repeated R/TS/P pattern in
    ``workflows/all.py:cli()`` TSOPT segment_log construction. Pure
    function: kcal/mol projection lives here, three call sites in cli()
    feed in their level-specific energies + reference.
    """
    kcal = [(e - ref_energy) * au_to_kcal for e in energies_au]
    return {
        "labels": list(labels),
        "energies_au": list(energies_au),
        "energies_kcal": kcal,
        "diagram": diagram_path,
        "structures": dict(structures),
        "barrier_kcal": kcal[1] if len(kcal) > 1 else 0.0,
        "delta_kcal": kcal[-1] if kcal else 0.0,
    }


def build_pipeline_summary_payload(
    *,
    out_dir: Path,
    path_dir: Path,
    summary: Dict[str, Any],
    refine_path: bool,
    do_tsopt: bool,
    do_thermo: bool,
    do_dft: bool,
    dft_func_basis_use: Optional[str],
    opt_mode: Optional[str],
    opt_mode_post: Optional[str],
    path_opt_mode: Optional[str],
    post_opt_mode: Optional[str],
    ts_opt_mode: Optional[str],
    endpoint_opt_mode: Optional[str],
    mep_mode_kind: Optional[str],
    dmf_correlated: bool,
    mlip_backend: str,
    mlip_model: Optional[str],
    mlip_precision: Optional[str],
    command_str: str,
    q_int: int,
    spin: int,
    freeze_atoms: Sequence[int],
    post_segment_logs: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    """Assemble the summary_log payload for the `all` pipeline.

    Extracted from the inner body of the nested
    ``_write_pipeline_summary_log`` in ``workflows/all.py:cli()`` so the
    dict-construction is unit-testable separately from the I/O wrapper.
    """
    diag_for_log: Dict[str, Any] = {}
    for diag in summary.get("energy_diagrams", []) or []:
        if isinstance(diag, dict) and str(diag.get("name", "")).lower().endswith("mep"):
            diag_for_log = diag
            break
    # MEP products are MOVED to the pipeline root (out_dir); they no longer
    # live under _work/path_search (path_dir).  Point the existence checks at
    # the root, and reference the canonical energy diagram name (the legacy
    # mep_plot.png is engine scratch, not a deliverable).
    mep_info = {
        "n_images": summary.get("n_images"),
        "n_segments": summary.get("n_segments"),
        "traj_pdb": str(out_dir / "mep.pdb") if (out_dir / "mep.pdb").exists() else None,
        "mep_plot": str(out_dir / "energy_diagram_MEP.png") if (out_dir / "energy_diagram_MEP.png").exists() else None,
        "diagram": diag_for_log,
    }
    return {
        "root_out_dir": str(out_dir),
        "path_dir": str(path_dir),
        "path_module_dir": path_dir.name,
        "pipeline_mode": "path-search" if refine_path else "path-opt",
        "refine_path": refine_path,
        "tsopt": do_tsopt,
        "thermo": do_thermo,
        "dft": do_dft,
        "dft_func_basis": dft_func_basis_use if do_dft else None,
        "opt_mode": opt_mode.lower() if opt_mode else None,
        "opt_mode_post": opt_mode_post.lower() if opt_mode_post else None,
        "path_opt_mode": (
            path_opt_mode.lower() if path_opt_mode else None
        ),
        "post_opt_mode": (
            post_opt_mode.lower() if post_opt_mode else None
        ),
        "ts_opt_mode": ts_opt_mode.lower() if ts_opt_mode else None,
        "endpoint_opt_mode": (
            endpoint_opt_mode.lower() if endpoint_opt_mode else None
        ),
        "mep_mode": mep_mode_kind,
        "dmf_correlated": bool(dmf_correlated),
        "mlip_backend": mlip_backend,
        "mlip_model": mlip_model,
        "mlip_precision": mlip_precision,
        "status": summary.get("status"),
        "status_reasons": summary.get("status_reasons", []),
        "execution_status": summary.get("execution_status"),
        "scientific_status": summary.get("scientific_status"),
        "scientific_status_reasons": summary.get(
            "scientific_status_reasons", []
        ),
        "command": command_str,
        "charge": q_int,
        "spin": spin,
        "freeze_atoms": list(freeze_atoms),
        "mep": mep_info,
        "segments": summary.get("segments", []),
        "energy_diagrams": summary.get("energy_diagrams", []),
        "post_segments": list(post_segment_logs),
        "key_files": {},
    }


def build_thermo_symmetry_provenance(
    thermo_payloads: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Copy complete, child-reported R/TS/P symmetry provenance."""

    provenance: Dict[str, Dict[str, Any]] = {}
    for label in ("R", "TS", "P"):
        payload = thermo_payloads.get(label)
        if not isinstance(payload, Mapping):
            continue
        value = payload.get("symmetry_number")
        source = payload.get("symmetry_number_source")
        if (
            isinstance(value, int)
            and not isinstance(value, bool)
            and value > 0
            and isinstance(source, str)
            and source.strip()
        ):
            state_provenance = {
                "symmetry_number": value,
                "symmetry_number_source": source,
            }
            point_group = payload.get("point_group")
            point_group_source = payload.get("point_group_source")
            if (
                isinstance(point_group, str)
                and point_group.strip()
                and isinstance(point_group_source, str)
                and point_group_source.strip()
            ):
                state_provenance.update(
                    point_group=point_group,
                    point_group_source=point_group_source,
                )
            provenance[label] = state_provenance
    return provenance


__all__ = [
    "build_energy_level_dict",
    "build_pipeline_summary_payload",
    "build_thermo_symmetry_provenance",
]
