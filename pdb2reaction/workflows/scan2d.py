# pdb2reaction/workflows/scan2d.py

"""
2D grid scan with harmonic restraints on two inter-atomic distances.

Example:
    pdb2reaction scan2d -i input.pdb -q 0 --scan-lists '[(12,45,1.30,3.10),(10,55,1.20,3.20)]'

For detailed documentation, see: docs/scan2d.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import logging
import math
import sys
import textwrap
import traceback
import shutil
import tempfile
import time

import click
import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import plotly.graph_objects as go

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, AU2KCALPERMOL

from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    BIAS_KW,
    OPT_MODE_ALIASES,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    UMA_CALC_KW,
    OUT_DIR_SCAN2D,
    apply_backend_defaults,
)
from pdb2reaction.backends import create_calculator
from pdb2reaction.workflows.restraints import HarmonicBiasCalculator
from pdb2reaction.core.utils import (
    axis_label_csv,
    axis_label_html,
    build_sopt_kwargs,
    is_scan_spec_file,
    make_sopt_optimizer,
    parse_scan_list_quads_checked,
    parse_scan_spec_quads,
    unbiased_energy_hartree,
    values_from_bounds,
    pretty_block,
    strip_inherited_keys,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    ensure_dir,
    make_snapshot_geometry,
    convert_xyz_like_outputs,
    build_scan_configs,
    cli_param_overridden,
    resolve_freeze_atoms,
    distance_A_from_coords,
    distance_tag,
    set_freeze_atoms_or_warn,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
)
from pdb2reaction.workflows.scan_common import (
    add_scan_common_options,
    load_merged_yaml_cfg,
    resolve_yaml_sources,
)
from pdb2reaction.cli.decorators import _write_error_json
from pdb2reaction.cli.common_options import (
    add_coord_type_option,
    add_print_every_option,
    add_precision_option,
    add_deterministic_option,
)

logger = logging.getLogger(__name__)

_snapshot_geometry = make_snapshot_geometry(GEOM_KW_DEFAULT["coord_type"])


def _sort_values_by_reference(values: np.ndarray, ref: Optional[float]) -> np.ndarray:
    """Sort scan values so that those closest to ref come first."""
    if ref is None or not np.isfinite(ref):
        return values
    order = np.argsort(np.abs(values - ref))
    return values[order]


def _build_scan_context(
    *,
    yaml_cfg: Dict[str, Any],
    geom_kw: Dict[str, Any],
    calc_kw: Dict[str, Any],
    opt_kw: Dict[str, Any],
    lbfgs_kw: Dict[str, Any],
    rfo_kw: Dict[str, Any],
    bias_kw: Dict[str, Any],
    charge: Optional[int],
    spin: Optional[int],
    workers: int,
    workers_per_node: int,
    out_dir: str,
    thresh: Optional[str],
    bias_k: float,
    opt_mode: str,
    relax_max_cycles: int,
    relax_override_requested: bool,
    max_step_size: float,
    source_path: Optional[Path],
    freeze_links: bool,
    set_charge_spin: bool = True,
    workers_overridden: bool = True,
    workers_per_node_overridden: bool = True,
) -> Tuple[
    Dict[str, Any],
    Dict[str, Any],
    Dict[str, Any],
    Dict[str, Any],
    Dict[str, Any],
    Dict[str, Any],
    str,
    List[int],
    Path,
]:
    geom_cfg, calc_cfg, opt_cfg, lbfgs_cfg, rfo_cfg, bias_cfg = build_scan_configs(
        yaml_cfg,
        geom_kw=geom_kw,
        calc_kw=calc_kw,
        opt_kw=opt_kw,
        lbfgs_kw=lbfgs_kw,
        rfo_kw=rfo_kw,
        bias_kw=bias_kw,
        charge=charge,
        spin=spin,
        workers=workers,
        workers_per_node=workers_per_node,
        out_dir=out_dir,
        thresh=thresh,
        bias_k=bias_k,
        set_charge_spin=set_charge_spin,
        workers_overridden=workers_overridden,
        workers_per_node_overridden=workers_per_node_overridden,
    )

    kind = normalize_choice(
        opt_mode,
        param="--opt-mode",
        alias_groups=OPT_MODE_ALIASES,
        allowed_hint="grad|hess",
    )

    # Convert 1-based YAML freeze_atoms to 0-based internal
    if geom_cfg.get("freeze_atoms"):
        geom_cfg["freeze_atoms"] = yaml_freeze_to_internal(geom_cfg["freeze_atoms"])
    freeze: List[int] = []
    if source_path is not None:
        freeze = resolve_freeze_atoms(geom_cfg, source_path, freeze_links)
    calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
    calc_cfg["return_partial_hessian"] = True

    out_dir_path = Path(opt_cfg["out_dir"]).resolve()
    ensure_dir(out_dir_path)
    echo_geom = format_geom_for_echo(geom_cfg)
    echo_calc = format_geom_for_echo(calc_cfg)
    echo_opt = dict(opt_cfg)
    if relax_override_requested:
        echo_opt["max_cycles"] = int(relax_max_cycles)
    echo_opt["out_dir"] = str(out_dir_path)
    echo_bias = dict(bias_cfg)
    click.echo(pretty_block("geom", echo_geom))
    click.echo(pretty_block("calc", echo_calc))
    click.echo(pretty_block("opt", echo_opt))
    max_step_bohr_for_log = float(max_step_size) * ANG2BOHR
    echo_sopt = build_sopt_kwargs(
        kind,
        lbfgs_cfg,
        rfo_cfg,
        opt_cfg,
        max_step_bohr_for_log,
        relax_max_cycles,
        relax_override_requested,
        out_dir_path,
        str(opt_cfg.get("prefix", "")),
    )
    echo_sopt = strip_inherited_keys(echo_sopt, opt_cfg)
    click.echo(
        pretty_block("lbfgs" if kind == "lbfgs" else "rfo", echo_sopt)
    )
    click.echo(pretty_block("bias", echo_bias))

    return (
        geom_cfg,
        calc_cfg,
        opt_cfg,
        lbfgs_cfg,
        rfo_cfg,
        bias_cfg,
        kind,
        freeze,
        out_dir_path,
    )


@click.command(
    help="2D distance scan with harmonic restraints.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, _trj.xyz, ...).",
)
@click.option(
    "-s", "--scan-lists",
    "scan_list_raw",
    type=str,
    required=True,
    help=(
        "Scan targets: inline Python literal or a YAML/JSON spec file path. "
        "scan2d expects EXACTLY 2 quadruples (i, j, low, high) — one per "
        "scanned bond axis — e.g. '[(12,45,1.30,3.10),(10,55,1.20,3.20)]'. "
        "Atom indices may also be PDB-style strings like 'CE  SAM   216'. "
        "Step count per axis is set via --max-step-size, NOT inside the tuple "
        "(scan2d does not accept a 5th element; the scan command's 3-tuple "
        "(i, j, target) form is rejected here too)."
    ),
)
@add_scan_common_options(
    workers_default=UMA_CALC_KW["workers"],
    workers_per_node_default=UMA_CALC_KW["workers_per_node"],
    out_dir_default=OUT_DIR_SCAN2D,
    baseline_help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0).",
    dump_help="Write inner scan trajectories per d1-step as TRJ under result_scan2d/grid/.",
)
@click.option(
    "--print-parsed/--no-print-parsed",
    "print_parsed",
    default=False,
    show_default=True,
    help="Print parsed scan targets after resolving --scan-lists.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help=(
        "Resolve and validate options (input, charge/spin parity, "
        "--scan-lists parse) and print the planned scan, then exit "
        "without running any optimization."
    ),
)
@add_coord_type_option()
@add_print_every_option()
@add_precision_option()
@add_deterministic_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    scan_list_raw: Optional[str],
    one_based: bool,
    max_step_size: float,
    bias_k: float,
    relax_max_cycles: int,
    opt_mode: str,
    freeze_links: bool,
    freeze_atoms_text: Optional[str],
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    preopt: bool,
    print_parsed: bool,
    baseline: str,
    zmin: Optional[float],
    zmax: Optional[float],
    out_json: bool,
    dry_run: bool,
    backend: str,
    solvent: str,
    solvent_model: str,
    cli_coord_type: Optional[str],
    print_every: Optional[int],
    precision: Optional[str],
) -> None:

    set_convert_file_enabled(convert_files)
    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )

    cycles_overridden = cli_param_overridden(ctx, "relax_max_cycles")

    with prepared_cli_input(
        input_path,
        ref_pdb=ref_pdb,
        charge=charge,
        spin=spin,
        ligand_charge=ligand_charge,
        prefix="[scan2d]",
    ) as (prepared_input, resolved_charge, resolved_spin):
        validate_charge_spin_for_prepared(prepared_input, resolved_charge, resolved_spin)
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path

        tmp_root = None
        out_dir_path = Path(out_dir).resolve()
        if dry_run:
            # Stop here AFTER parsing --scan-lists (the help text promises that
            # --dry-run validates the scan-lists spec). Input has been resolved,
            # charge/spin parity has been validated against the prepared
            # geometry, and --scan-lists is parsed below before exit.
            _dry_parsed, _ = parse_scan_list_quads_checked(
                scan_list_raw,
                expected_len=2,
                one_based=bool(one_based),
                atom_meta=[],
                option_name="--scan-lists",
            )
            click.echo("[scan2d] --dry-run: input, charge/spin parity, and --scan-lists parse OK.")
            click.echo(f"[scan2d] input geometry  : {geom_input_path}")
            click.echo(f"[scan2d] resolved charge : {int(resolved_charge):+d}")
            click.echo(f"[scan2d] resolved spin   : {int(resolved_spin)} (multiplicity)")
            click.echo(f"[scan2d] out_dir         : {out_dir_path}")
            click.echo(f"[scan2d] --scan-lists    : {scan_list_raw} → {len(_dry_parsed)} axis tuples")
            click.echo(f"[scan2d] preopt={bool(preopt)}  freeze_links={bool(freeze_links)}")
            click.echo("[scan2d] No 2D scan was executed.")
            return
        try:
            time_start = time.perf_counter()

            yaml_cfg = load_merged_yaml_cfg(
                config_yaml=config_yaml,
                override_yaml=None,
            )
            yaml_opt = yaml_cfg.get("opt") if isinstance(yaml_cfg, dict) else None
            relax_override_requested = cycles_overridden and not (
                isinstance(yaml_opt, dict) and "max_cycles" in yaml_opt
            )

            (
                geom_cfg,
                calc_cfg,
                opt_cfg,
                lbfgs_cfg,
                rfo_cfg,
                bias_cfg,
                kind,
                freeze,
                out_dir_path,
            ) = _build_scan_context(
                yaml_cfg=yaml_cfg,
                geom_kw=dict(GEOM_KW_DEFAULT),
                calc_kw=dict(UMA_CALC_KW),
                opt_kw=dict(OPT_BASE_KW),
                lbfgs_kw=dict(LBFGS_KW),
                rfo_kw=dict(RFO_KW),
                bias_kw=dict(BIAS_KW),
                charge=resolved_charge,
                spin=resolved_spin,
                workers=workers,
                workers_per_node=workers_per_node,
                out_dir=out_dir,
                thresh=thresh,
                # bias_k is None when neither CLI --bias-k nor YAML bias.k set
                # (the common-decorator default flipped to None to enable YAML
                # override). `build_scan_configs` handles None via
                # `if bias_k is not None: bias_cfg["k"] = float(bias_k)`.
                bias_k=bias_k,
                opt_mode=opt_mode,
                relax_max_cycles=relax_max_cycles,
                relax_override_requested=relax_override_requested,
                max_step_size=max_step_size,
                source_path=source_path,
                freeze_links=freeze_links,
                workers_overridden=cli_param_overridden(ctx, "workers"),
                workers_per_node_overridden=cli_param_overridden(ctx, "workers_per_node"),
            )

            # Merge CLI --freeze-atoms (already 0-based)
            try:
                freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            if freeze_atoms_cli:
                merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)
                freeze = list(geom_cfg.get("freeze_atoms", []))
                calc_cfg["freeze_atoms"] = freeze

            if cli_param_overridden(ctx, "backend"):
                calc_cfg["backend"] = backend
            if cli_param_overridden(ctx, "solvent"):
                calc_cfg["solvent"] = solvent
            if cli_param_overridden(ctx, "solvent_model"):
                calc_cfg["solvent_model"] = solvent_model
            if precision is not None:
                from pdb2reaction.backends import apply_precision_to_calc_cfg
                apply_precision_to_calc_cfg(calc_cfg, precision)
            if cli_param_overridden(ctx, "print_every") and print_every is not None:
                opt_cfg["print_every"] = int(print_every)
            if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
                geom_cfg["coord_type"] = str(cli_coord_type).lower()
            apply_backend_defaults(calc_cfg)

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            if scan_list_raw is None:
                raise click.BadParameter("--scan-lists is required.")
            scan_one_based = bool(one_based)
            scan_source = "--scan-lists"
            if is_scan_spec_file(scan_list_raw):
                spec_path = Path(scan_list_raw)
                parsed, raw_pairs, scan_one_based = parse_scan_spec_quads(
                    spec_path,
                    expected_len=2,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
                scan_source = f"--scan-lists ({spec_path})"
            else:
                parsed, raw_pairs = parse_scan_list_quads_checked(
                    scan_list_raw,
                    expected_len=2,
                    one_based=scan_one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
            (i1, j1, low1, high1), (i2, j2, low2, high2) = parsed
            d1_label_csv = axis_label_csv("d1", i1, j1, scan_one_based, pdb_atom_meta, raw_pairs[0])
            d2_label_csv = axis_label_csv("d2", i2, j2, scan_one_based, pdb_atom_meta, raw_pairs[1])
            d1_label_html = axis_label_html(d1_label_csv)
            d2_label_html = axis_label_html(d2_label_csv)
            if print_parsed:
                click.echo(
                    pretty_block(
                        "scan-parsed",
                        {
                            "source": scan_source,
                            "one_based": bool(scan_one_based),
                            "pairs": parsed,
                        },
                    )
                )
            click.echo(
                pretty_block(
                    "scan-list (1-based)",
                    {"d1": (i1+1, j1+1, low1, high1), "d2": (i2+1, j2+1, low2, high2)},
                )
            )

            if pdb_atom_meta:
                click.echo("[scan2d] PDB atom details for scanned pairs:", detail=True)
                legend = format_pdb_atom_metadata_header()
                click.echo(f"        legend: {legend}", detail=True)
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}", detail=True)
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}", detail=True)

            # Temporary and grid directories
            tmp_root = Path(tempfile.mkdtemp(prefix="scan2d_tmp_"))
            grid_dir = out_dir_path / "grid"
            tmp_opt_dir = tmp_root / "opt"
            ensure_dir(grid_dir)
            ensure_dir(tmp_opt_dir)

            final_dir = out_dir_path

            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            geom_outer = geom_loader(
                geom_input_path, coord_type=coord_type, freeze_atoms=freeze
            )
            set_freeze_atoms_or_warn(geom_outer, freeze, context="scan2d")

            base_calc = create_calculator(**calc_cfg)
            echo_resolved_device()
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

            # Records (including preopt) will be accumulated here
            records: List[Dict[str, Any]] = []
            ref_pdb_path = source_path if source_path.suffix.lower() == ".pdb" else None

            # Reference distances from the (pre)optimized structure, used for scan ordering
            d1_ref: Optional[float] = None
            d2_ref: Optional[float] = None

            # Cache of previously converged geometries for nearest-start logic:
            # each entry is (d1_A, d2_A, geometry_snapshot)
            visited_geoms: List[Tuple[float, float, Any]] = []

            if preopt:
                click.echo("[preopt] Unbiased relaxation of the initial structure ...")
                geom_outer.set_calculator(base_calc)
                max_step_bohr_local = float(max_step_size) * ANG2BOHR
                optimizer0 = make_sopt_optimizer(
                    geom_outer,
                    kind,
                    lbfgs_cfg,
                    rfo_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr_local,
                    relax_max_cycles=relax_max_cycles,
                    relax_override_requested=relax_override_requested,
                    out_dir=tmp_opt_dir,
                    prefix="preopt",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.")
                except OptimizationError as e:
                    click.echo(f"[preopt] OptimizationError — {e}")

                # Measure optimized distances and record preopt structure
                try:
                    coords_outer = np.asarray(getattr(geom_outer, "coords3d"), dtype=float)
                    d1_ref = distance_A_from_coords(coords_outer, i1, j1)
                    d2_ref = distance_A_from_coords(coords_outer, i2, j2)

                    d1_tag = distance_tag(d1_ref)
                    d2_tag = distance_tag(d2_ref)

                    preopt_xyz_path = grid_dir / f"preopt_i{d1_tag}_j{d2_tag}.xyz"
                    s = geom_outer.as_xyz()
                    if not s.endswith("\n"):
                        s += "\n"
                    with open(preopt_xyz_path, "w") as f:
                        f.write(s)

                    convert_xyz_like_outputs(
                        preopt_xyz_path,
                        prepared_input,
                        ref_pdb_path=ref_pdb_path,
                        out_pdb_path=grid_dir / f"preopt_i{d1_tag}_j{d2_tag}.pdb",
                        out_gjf_path=grid_dir / f"preopt_i{d1_tag}_j{d2_tag}.gjf",
                        context=f"'{preopt_xyz_path.name}' to PDB/GJF",
                    )

                    E_pre_h = unbiased_energy_hartree(geom_outer, base_calc)
                    records.append(
                        {
                            "i": int(-1),
                            "j": int(-1),
                            "d1_A": float(d1_ref),
                            "d2_A": float(d2_ref),
                            "energy_hartree": E_pre_h,
                            "bias_converged": True,
                        }
                    )
                    # Store preoptimized geometry as a candidate for nearest-start
                    visited_geoms.append(
                        (float(d1_ref), float(d2_ref), _snapshot_geometry(geom_outer))
                    )

                    click.echo(
                        f"[preopt] Recorded preoptimized structure at d1={d1_ref:.3f} Å, d2={d2_ref:.3f} Å."
                    )
                except Exception as e:
                    click.echo(
                        f"[preopt] WARNING: failed to record preoptimized structure: {e}",
                        err=True,
                    )

            max_step_bohr = float(max_step_size) * ANG2BOHR

            # Construct scan grids and reorder so that points near the preopt geometry are visited first
            d1_values = values_from_bounds(low1, high1, float(max_step_size))
            d2_values = values_from_bounds(low2, high2, float(max_step_size))

            d1_values = _sort_values_by_reference(d1_values, d1_ref)
            d2_values = _sort_values_by_reference(d2_values, d2_ref)

            N1, N2 = len(d1_values), len(d2_values)
            click.echo(
                f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x:f'{x:.3f}', d1_values))}"
            )
            click.echo(
                f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x:f'{x:.3f}', d2_values))}"
            )
            click.echo(f"[grid] total grid points = {N1*N2}", narrative=True)

            # Track the start index of the previous row's entries in visited_geoms
            # so we can prune old rows and bound memory usage.
            _prev_row_start: int = 0

            for i_idx, d1_target in enumerate(d1_values):
                # Prune visited_geoms: keep only the previous row's entries
                # (current row hasn't started yet).  Entries before _prev_row_start
                # belong to even older rows and are unlikely nearest neighbors.
                if _prev_row_start > 0:
                    visited_geoms = visited_geoms[_prev_row_start:]
                    _prev_row_start = 0
                # Mark where the current row starts; this becomes _prev_row_start
                # at the beginning of the *next* outer iteration.
                _cur_row_start = len(visited_geoms)

                click.echo(
                    f"[stage] d1 step {i_idx+1}/{N1}: target = {d1_target:.3f} Å"
                )
                biased.set_pairs([(i1, j1, float(d1_target))])
                geom_outer.set_calculator(biased)

                opt1 = make_sopt_optimizer(
                    geom_outer,
                    kind,
                    lbfgs_cfg,
                    rfo_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr,
                    relax_max_cycles=relax_max_cycles,
                    relax_override_requested=relax_override_requested,
                    out_dir=tmp_opt_dir,
                    prefix=f"d1_{i_idx:03d}",
                )
                try:
                    opt1.run()
                except ZeroStepLength:
                    click.echo(
                        f"[d1 {i_idx}] ZeroStepLength — continuing to d2 scan."
                    )
                except OptimizationError as e:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {e}")

                geom_inner = _snapshot_geometry(geom_outer)
                geom_inner.set_calculator(biased)

                # Store the d1-relaxed structure as a candidate for nearest-start
                try:
                    coords_inner = np.asarray(getattr(geom_inner, "coords3d"), dtype=float)
                    d1_cur = distance_A_from_coords(coords_inner, i1, j1)
                    d2_cur = distance_A_from_coords(coords_inner, i2, j2)
                    visited_geoms.append(
                        (float(d1_cur), float(d2_cur), _snapshot_geometry(geom_inner))
                    )
                except Exception as e:
                    click.echo(
                        f"[nearest-start] WARNING: failed to store d1-relaxed structure for d1={d1_target:.3f} Å: {e}",
                        err=True,
                    )

                trj_blocks = [] if dump else None

                for j_idx, d2_target in enumerate(d2_values):
                    # Choose initial structure: nearest previously converged (d1,d2) point
                    if visited_geoms:
                        try:
                            target_vec = np.array(
                                [float(d1_target), float(d2_target)], dtype=float
                            )
                            prev_coords = np.array(
                                [(g[0], g[1]) for g in visited_geoms],
                                dtype=float,
                            )
                            dists2 = np.sum((prev_coords - target_vec) ** 2, axis=1)
                            best_idx = int(np.argmin(dists2))
                            _, _, best_geom = visited_geoms[best_idx]
                            # Reset geom_inner coordinates to the best previous geometry
                            try:
                                geom_inner.coords[:] = np.array(
                                    best_geom.coords, copy=True
                                )
                            except Exception:
                                geom_inner.coords = np.array(
                                    best_geom.coords, copy=True
                                )
                        except Exception as e:
                            click.echo(
                                f"[nearest-start] WARNING: failed to select nearest previous structure for d1={d1_target:.3f}, d2={d2_target:.3f}: {e}",
                                err=True,
                            )

                    biased.set_pairs(
                        [
                            (i1, j1, float(d1_target)),
                            (i2, j2, float(d2_target)),
                        ]
                    )
                    geom_inner.set_calculator(biased)

                    opt2 = make_sopt_optimizer(
                        geom_inner,
                        kind,
                        lbfgs_cfg,
                        rfo_cfg,
                        opt_cfg,
                        max_step_bohr=max_step_bohr,
                        relax_max_cycles=relax_max_cycles,
                        relax_override_requested=relax_override_requested,
                        out_dir=tmp_opt_dir,
                        prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}",
                    )
                    try:
                        opt2.run()
                        converged = True
                    except ZeroStepLength:
                        click.echo(
                            f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — recorded anyway.",
                            err=True,
                        )
                        converged = False
                    except OptimizationError as e:
                        click.echo(
                            f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {e}"
                        )
                        converged = False

                    E_h = unbiased_energy_hartree(geom_inner, base_calc)

                    # Write per-grid XYZ snapshots under result_scan2d/grid/
                    d1_tag = distance_tag(d1_target)
                    d2_tag = distance_tag(d2_target)
                    xyz_path = grid_dir / f"point_i{d1_tag}_j{d2_tag}.xyz"
                    try:
                        s = geom_inner.as_xyz()
                        if not s.endswith("\n"):
                            s += "\n"
                        with open(xyz_path, "w") as f:
                            f.write(s)
                        convert_xyz_like_outputs(
                            xyz_path,
                            prepared_input,
                            ref_pdb_path=ref_pdb_path,
                            out_pdb_path=grid_dir / f"point_i{d1_tag}_j{d2_tag}.pdb",
                            out_gjf_path=grid_dir / f"point_i{d1_tag}_j{d2_tag}.gjf",
                            context=f"'{xyz_path.name}' to PDB/GJF",
                        )
                    except Exception as e:
                        click.echo(
                            f"[write] WARNING: failed to write {xyz_path.name}: {e}",
                            err=True,
                        )

                    # Store this converged grid point for nearest-start initialization
                    try:
                        coords_inner = np.asarray(getattr(geom_inner, "coords3d"), dtype=float)
                        d1_cur = distance_A_from_coords(coords_inner, i1, j1)
                        d2_cur = distance_A_from_coords(coords_inner, i2, j2)
                        visited_geoms.append(
                            (float(d1_cur), float(d2_cur), _snapshot_geometry(geom_inner))
                        )
                    except Exception as e:
                        click.echo(
                            f"[nearest-start] WARNING: failed to store geometry for d1={d1_target:.3f}, d2={d2_target:.3f}: {e}",
                            err=True,
                        )

                    if dump and trj_blocks is not None:
                        sblock = geom_inner.as_xyz()
                        if not sblock.endswith("\n"):
                            sblock += "\n"
                        trj_blocks.append(sblock)

                    records.append(
                        {
                            "i": int(i_idx),
                            "j": int(j_idx),
                            "d1_A": float(d1_target),
                            "d2_A": float(d2_target),
                            "energy_hartree": E_h,
                            "bias_converged": bool(converged),
                        }
                    )

                # Update row tracking: current row entries start at _cur_row_start
                _prev_row_start = _cur_row_start

                if dump and trj_blocks:
                    trj_path = grid_dir / f"inner_path_d1_{i_idx:03d}_trj.xyz"
                    try:
                        with open(trj_path, "w") as f:
                            f.write("".join(trj_blocks))
                        click.echo(f"[write] Wrote '{trj_path}'.")
                        convert_xyz_like_outputs(
                            trj_path,
                            prepared_input,
                            ref_pdb_path=ref_pdb_path,
                            out_pdb_path=grid_dir / f"inner_path_d1_{i_idx:03d}.pdb",
                            context=f"'{trj_path.name}' to PDB",
                        )
                    except Exception as e:
                        click.echo(
                            f"[write] WARNING: failed to write '{trj_path}': {e}"
                        )

            # ===== surface.csv (final output directly under result_scan2d) =====
            df = pd.DataFrame.from_records(records)
            if df.empty:
                click.echo("No grid records produced; aborting.")
                sys.exit(1)

            if baseline == "first":
                mask = (df["i"] == 0) & (df["j"] == 0)
                if mask.sum() == 0:
                    click.echo("WARNING: baseline='first' but grid point (0,0) not found; falling back to min.", err=True)
                    ref = float(df["energy_hartree"].min())
                else:
                    ref = float(df.loc[mask, "energy_hartree"].iloc[0])
            else:
                ref = float(df["energy_hartree"].min())
            df["energy_kcal"] = (df["energy_hartree"] - ref) * AU2KCALPERMOL
            df["d1_label"] = d1_label_csv
            df["d2_label"] = d2_label_csv

            surface_csv = final_dir / "surface.csv"
            df.to_csv(surface_csv, index=False)
            click.echo(f"[write] Wrote '{surface_csv}'.")

            # ===== Plots (RBF on a fixed 50×50 grid, unified layout, placed under final_dir) =====
            d1_points = df["d1_A"].to_numpy(dtype=float)
            d2_points = df["d2_A"].to_numpy(dtype=float)
            z_points = df["energy_kcal"].to_numpy(dtype=float)
            mask = (
                np.isfinite(d1_points)
                & np.isfinite(d2_points)
                & np.isfinite(z_points)
            )
            if not np.any(mask):
                click.echo("[plot] No finite data for plotting.")
                sys.exit(1)

            x_min, x_max = float(np.min(d1_points[mask])), float(
                np.max(d1_points[mask])
            )
            y_min, y_max = float(np.min(d2_points[mask])), float(
                np.max(d2_points[mask])
            )

            xi = np.linspace(x_min, x_max, 50)
            yi = np.linspace(y_min, y_max, 50)
            XI, YI = np.meshgrid(xi, yi)

            rbf = Rbf(
                d1_points[mask], d2_points[mask], z_points[mask], function="multiquadric"
            )
            ZI = rbf(XI, YI)

            vmin = float(np.nanmin(ZI)) if zmin is None else float(zmin)
            vmax = float(np.nanmax(ZI)) if zmax is None else float(zmax)
            if (
                not np.isfinite(vmin)
                or not np.isfinite(vmax)
                or vmax <= vmin
            ):
                vmin, vmax = float(np.nanmin(ZI)), float(np.nanmax(ZI))

            # Choose neat contour/tick steps
            def _nice_step(span: float) -> float:
                if span <= 0:
                    return 1.0
                raw = span / 6.0
                mag = 10 ** math.floor(math.log10(raw))
                candidates = (0.5, 1, 2, 5, 10, 20)
                best = candidates[0] * mag
                best_err = abs(best - raw)
                for m in candidates[1:]:
                    s = m * mag
                    err = abs(s - raw)
                    if err < best_err:
                        best, best_err = s, err
                return best

            c_step = _nice_step(vmax - vmin)
            c_start = math.floor(vmin / c_step) * c_step
            c_end = math.ceil(vmax / c_step) * c_step

            # ---- 2D contour plot (PNG with explicit size) ----
            fig2d = go.Figure(
                data=go.Contour(
                    z=ZI,
                    x=xi,
                    y=yi,
                    contours=dict(start=c_start, end=c_end, size=c_step),
                    zmin=vmin,
                    zmax=vmax,
                    contours_coloring="heatmap",
                    colorscale="plasma",
                    colorbar=dict(
                        title=dict(
                            text="(kcal/mol)", side="top", font=dict(size=16, color="#1C1C1C")
                        ),
                        tickfont=dict(size=14, color="#1C1C1C"),
                        ticks="inside",
                        ticklen=10,
                        tickcolor="#1C1C1C",
                        outlinecolor="#1C1C1C",
                        outlinewidth=2,
                        lenmode="fraction",
                        len=1.11,
                        x=1.05,
                        y=0.53,
                        xanchor="left",
                        yanchor="middle",
                    ),
                )
            )
            fig2d.update_layout(
                width=640,
                height=600,
                xaxis_title=d1_label_html,
                yaxis_title=d2_label_html,
                plot_bgcolor="white",
                xaxis=dict(
                    range=[x_min, x_max],
                    showline=True,
                    linewidth=3,
                    linecolor="#1C1C1C",
                    mirror=True,
                    tickson="boundaries",
                    ticks="inside",
                    tickwidth=3,
                    tickcolor="#1C1C1C",
                    title_font=dict(size=18, color="#1C1C1C"),
                    tickfont=dict(size=18, color="#1C1C1C"),
                    tickvals=list(np.linspace(x_min, x_max, 6)),
                    tickformat=".2f",
                ),
                yaxis=dict(
                    range=[y_min, y_max],
                    showline=True,
                    linewidth=3,
                    linecolor="#1C1C1C",
                    mirror=True,
                    tickson="boundaries",
                    ticks="inside",
                    tickwidth=3,
                    tickcolor="#1C1C1C",
                    title_font=dict(size=18, color="#1C1C1C"),
                    tickfont=dict(size=18, color="#1C1C1C"),
                    tickvals=list(np.linspace(y_min, y_max, 6)),
                    tickformat=".2f",
                ),
                margin=dict(l=10, r=10, b=10, t=40),
            )
            png2d = final_dir / "scan2d_map.png"
            fig2d.write_image(str(png2d), scale=2, engine="kaleido", width=680, height=600)
            click.echo(f"[plot] Wrote '{png2d}'.")

            # ---- 3D surface plus base-plane projection ----
            spread = vmax - vmin if (vmax > vmin) else 1.0
            z_bottom = vmin - spread
            z_top = vmax

            # Avoid ticks below zmin (= vmin) and snap to sensible values
            z_step = _nice_step(vmax - vmin)
            z_start_tick = math.ceil(vmin / z_step) * z_step  # First tick must be ≥ vmin
            z_ticks = np.arange(z_start_tick, z_top + 0.5 * z_step, z_step).tolist()

            surface3d = go.Surface(
                x=XI,
                y=YI,
                z=ZI,
                colorscale="plasma",
                cmin=vmin,
                cmax=vmax,
                colorbar=dict(
                    title=dict(
                        text="(kcal/mol)", side="top", font=dict(size=16, color="#1C1C1C")
                    ),
                    tickfont=dict(size=14, color="#1C1C1C"),
                    ticks="inside",
                    ticklen=10,
                    tickcolor="#1C1C1C",
                    outlinecolor="#1C1C1C",
                    outlinewidth=2,
                    lenmode="fraction",
                    len=1.11,
                    x=1.05,
                    y=0.53,
                    xanchor="left",
                    yanchor="middle",
                ),
                contours={
                    "z": {
                        "show": True,
                        "start": c_start,
                        "end": c_end,
                        "size": c_step,
                        "color": "black",
                        "project": {"z": True},
                    }
                },
                name="3D Surface",
            )

            plane_proj = go.Surface(
                x=XI,
                y=YI,
                z=np.full_like(ZI, z_bottom),
                surfacecolor=ZI,
                colorscale="plasma",
                cmin=vmin,
                cmax=vmax,
                showscale=False,
                opacity=1.0,
                name="2D Contour Projection (Bottom)",
            )

            fig3d = go.Figure(data=[surface3d, plane_proj])
            fig3d.update_layout(
                title="Energy Landscape with 2D PES Scan",
                width=800,
                height=700,
                scene=dict(
                    bgcolor="rgba(0,0,0,0)",
                    xaxis=dict(
                        title=d1_label_html,
                        range=[x_min, x_max],
                        showline=True,
                        linewidth=4,
                        linecolor="#1C1C1C",
                        mirror=True,
                        ticks="inside",
                        tickwidth=4,
                        tickcolor="#1C1C1C",
                        gridcolor="rgba(0,0,0,0.1)",
                        zerolinecolor="rgba(0,0,0,0.1)",
                        showbackground=False,
                    ),
                    yaxis=dict(
                        title=d2_label_html,
                        range=[y_min, y_max],
                        showline=True,
                        linewidth=4,
                        linecolor="#1C1C1C",
                        mirror=True,
                        ticks="inside",
                        tickwidth=4,
                        tickcolor="#1C1C1C",
                        gridcolor="rgba(0,0,0,0.1)",
                        zerolinecolor="rgba(0,0,0,0.1)",
                        showbackground=False,
                    ),
                    zaxis=dict(
                        title="Potential Energy (kcal/mol)",
                        range=[z_bottom, z_top],
                        tickmode="array",
                        tickvals=z_ticks,
                        showline=True,
                        linewidth=4,
                        linecolor="#1C1C1C",
                        mirror=True,
                        ticks="inside",
                        tickwidth=4,
                        tickcolor="#1C1C1C",
                        showgrid=True,
                        gridcolor="rgba(0,0,0,0.1)",
                        zerolinecolor="rgba(0,0,0,0.1)",
                        showbackground=False,
                    ),
                ),
                margin=dict(l=10, r=20, b=10, t=40),
                paper_bgcolor="white",
            )

            html3d = final_dir / "scan2d_landscape.html"
            fig3d.write_html(str(html3d))
            click.echo(f"[plot] Wrote '{html3d}'.")

            click.echo("\n====== 2D Scan finished ======\n", narrative=True)
            click.echo(format_elapsed("[time] Elapsed Time for 2D Scan", time_start), narrative=True)

            # result.json (if --out-json)
            if out_json:
                from pdb2reaction.core.utils import write_result_json, atom_label_from_meta
                min_energy = float(df["energy_hartree"].min()) if not df.empty else None
                _pair1: Dict[str, Any] = {"i": int(i1 + 1), "j": int(j1 + 1), "low": float(low1), "high": float(high1)}
                _pair2: Dict[str, Any] = {"i": int(i2 + 1), "j": int(j2 + 1), "low": float(low2), "high": float(high2)}
                if pdb_atom_meta:
                    _pair1["label_i"] = atom_label_from_meta(pdb_atom_meta, i1)
                    _pair1["label_j"] = atom_label_from_meta(pdb_atom_meta, j1)
                    _pair2["label_i"] = atom_label_from_meta(pdb_atom_meta, i2)
                    _pair2["label_j"] = atom_label_from_meta(pdb_atom_meta, j2)
                result_data: Dict[str, Any] = {
                    "status": "completed",
                    "charge": resolved_charge,
                    "spin": resolved_spin,
                    "backend": calc_cfg.get("backend", backend),
                    "model": calc_cfg.get("model"),
                    "solvent": calc_cfg.get("solvent", "none"),
                    "max_step_size_angstrom": float(max_step_size),
                    "n_grid_points": len(df),
                    "grid_shape": [len(d1_values), len(d2_values)],
                    "pair1": _pair1,
                    "pair2": _pair2,
                    "min_energy_hartree": min_energy,
                    "files": {
                        "surface_csv": "surface.csv",
                        "scan2d_map_png": "scan2d_map.png",
                        "scan2d_landscape_html": "scan2d_landscape.html",
                    },
                }
                write_result_json(
                    final_dir, result_data,
                    command="scan2d",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

        except KeyboardInterrupt:
            click.echo("Interrupted by user.", err=True)
            sys.exit(130)
        except Exception as e:
            _write_error_json(out_dir_path, "scan2d", e, "UnhandledError", time_start)
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            click.echo(
                "Unhandled exception during 2D scan:\n"
                + textwrap.indent(tb, "  "),
                err=True,
            )
            sys.exit(1)
        finally:
            if tmp_root is not None:
                shutil.rmtree(tmp_root, ignore_errors=True)
