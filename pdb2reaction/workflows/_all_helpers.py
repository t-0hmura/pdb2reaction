"""Helpers extracted from `pdb2reaction/workflows/all.py:cli()`.

The `all` subcommand's `cli()` function is ~4,900 LOC; progressive
helper extracts land in this module so the cli() body shrinks one piece
at a time while behavior is preserved.

Anything imported here must be safe to use from ``cli()`` body callers
without changing observable behavior.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Sequence


@dataclass(frozen=True)
class AllContext:
    """Frozen bundle of the `pdb2reaction all` CLI parameters.

    Mirrors the ``cli()`` signature in `pdb2reaction/workflows/all.py`
    so helper functions can accept a single argument instead of
    re-listing 66 keyword parameters.

    Status: transitional foundation. The dataclass holds the canonical
    parameter set, and ``tests/test_all_helpers.py``'s drift guard
    enforces that the field names stay in lockstep with the
    ``cli.callback`` signature. cli() itself still binds parameters
    individually for back-compat with its existing nested closures,
    so AllContext is not yet exercised at runtime. Incremental
    decomposition extracts can migrate one helper at a time to take
    ``ctx: AllContext`` instead of an exploding kwargs list, knowing
    the drift guard prevents the dataclass from silently going out of
    sync while that migration proceeds.

    Field order matches the CLI option declaration order in
    `workflows/all.py` so the two read side-by-side cleanly.
    """

    # Inputs / output
    input_paths: Sequence[Path]
    center_spec: Optional[str]
    out_dir: Path
    # Extract knobs
    radius: float
    radius_het2het: float
    include_h2o: bool
    exclude_backbone: bool
    add_linkh: bool
    selected_resn: str
    modified_residue: str
    ligand_charge: Optional[str]
    charge_override: Optional[int]
    # Workers / backend
    workers: int
    workers_per_node: int
    backend: Optional[str]
    solvent: Optional[str]
    solvent_model: Optional[str]
    # Charges
    spin: int
    # Freeze / MEP
    freeze_links_flag: Optional[bool]
    mep_mode: str
    dmf_backend: str
    max_nodes: int
    max_cycles: int
    climb: bool
    opt_mode: str
    opt_mode_post: Optional[str]
    dump: bool
    convert_files: bool
    refine_path: bool
    thresh: Optional[str]
    thresh_post: str
    config_yaml: Optional[Path]
    show_config: bool
    dry_run: bool
    preopt: bool
    hessian_calc_mode: Optional[str]
    # Stage toggles
    do_tsopt: bool
    do_thermo: bool
    do_dft: bool
    # Scan
    scan_lists_raw: Sequence[str]
    scan_out_dir: Optional[Path]
    scan_one_based: Optional[bool]
    scan_max_step_size: Optional[float]
    scan_bias_k: Optional[float]
    scan_relax_max_cycles: Optional[int]
    scan_preopt_override: Optional[bool]
    scan_endopt_override: Optional[bool]
    ref_pdb_cli: Optional[Path]
    # TSOPT / FREQ / DFT
    tsopt_max_cycles: Optional[int]
    tsopt_out_dir: Optional[Path]
    flatten: bool
    freq_out_dir: Optional[Path]
    freq_max_write: Optional[int]
    freq_amplitude_ang: Optional[float]
    freq_n_frames: Optional[int]
    freq_sort: Optional[str]
    freq_temperature: Optional[float]
    freq_pressure: Optional[float]
    dft_out_dir: Optional[Path]
    dft_func_basis: Optional[str]
    dft_max_cycle: Optional[int]
    dft_conv_tol: Optional[float]
    dft_grid_level: Optional[int]
    dft_engine: Optional[str]
    cli_coord_type: Optional[str]
    precision: Optional[str]
    backend_model: Optional[str]
    calc_file: Optional[str] = None
    calc_factory: str = "get_calculator"


def build_energy_level_dict(
    *,
    labels: Sequence[str],
    energies_au: Sequence[float],
    ref_energy: float,
    au_to_kcal: float,
    diagram_path: str,
    structures: Dict[str, Any],
) -> Dict[str, Any]:
    """Assemble a segment_log sub-dict for one energy level (UMA / Gibbs / DFT).

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
    mep_mode_kind: Optional[str],
    uma_model: Optional[str],
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
        "mep_mode": mep_mode_kind,
        "uma_model": uma_model,
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


__all__ = [
    "AllContext",
    "build_energy_level_dict",
    "build_pipeline_summary_payload",
]
