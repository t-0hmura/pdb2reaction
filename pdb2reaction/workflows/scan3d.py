# pdb2reaction/workflows/scan3d.py

"""
3D grid scan with harmonic restraints on three inter-atomic distances.

Example:
    pdb2reaction scan3d -i input.pdb -q 0 --scan-lists '[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]'

For detailed documentation, see: docs/scan3d.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import logging
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
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    UMA_CALC_KW,
    OUT_DIR_SCAN3D,
    apply_backend_defaults,
)
from pdb2reaction.backends import create_calculator
from pdb2reaction.workflows.restraints import HarmonicBiasCalculator
from pdb2reaction.core.utils import (
    axis_label_csv,
    axis_label_html,
    is_scan_spec_file,
    make_sopt_optimizer,
    parse_scan_list_quads_checked,
    parse_scan_spec_quads,
    unbiased_energy_hartree,
    values_from_bounds,
    pretty_block,
    format_elapsed,
    prepared_cli_input,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    ensure_dir,
    make_snapshot_geometry,
    cli_param_overridden,
    distance_A_from_coords,
    distance_tag,
    set_freeze_atoms_or_warn,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
)
from pdb2reaction.workflows.scan_common import (
    add_scan_common_options,
    load_merged_yaml_cfg,
    resolve_yaml_sources,
)
from pdb2reaction.workflows.scan2d import _build_scan_context
from pdb2reaction.cli.decorators import _write_error_json
from pdb2reaction.cli.common_options import (
    add_coord_type_option,
    add_print_every_option,
    add_precision_option, add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)

logger = logging.getLogger(__name__)

_VOLUME_GRID_N = 50  # 50×50×50 RBF interpolation grid


_snapshot_geometry = make_snapshot_geometry(GEOM_KW_DEFAULT["coord_type"])


def _extract_axis_label(df: pd.DataFrame, column: str, fallback: Optional[str]) -> Optional[str]:
    if column not in df.columns:
        return fallback
    values = df[column].dropna()
    if values.empty:
        return fallback
    return str(values.iloc[0])


@click.command(
    help="3D distance scan with harmonic restraints.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="Input structure file (.pdb, .xyz, _trj.xyz, ...). Required unless --csv is provided.",
)
@click.option(
    "-s", "--scan-lists",
    "scan_list_raw",
    type=str,
    required=False,
    help=(
        "Scan targets: inline Python literal or a YAML/JSON spec file path. "
        "scan3d expects EXACTLY 3 quadruples (i, j, low, high) — one per "
        "scanned bond axis — e.g. "
        "'[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]'. "
        "Atom indices may also be PDB-style strings like 'CE  SAM   216'. "
        "Step count per axis is set via --max-step-size, NOT inside the tuple "
        "(scan3d does not accept a 5th element)."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "If provided, skip the 3D scan and read a precomputed surface.csv from this path. "
        "When used, -i/--input and --scan-lists are optional."
    ),
)
@add_scan_common_options(
    workers_default=UMA_CALC_KW["workers"],
    workers_per_node_default=UMA_CALC_KW["workers_per_node"],
    out_dir_default=OUT_DIR_SCAN3D,
    baseline_help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0,k=0).",
    dump_help="Write inner d3 scan trajectories per (d1,d2) as TRJ under result_scan3d/grid/.",
    max_step_help="Maximum step size per scanned distance [Å].",
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
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Optional[Path],
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
    csv_path: Optional[Path],
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
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: str,
) -> None:
    set_convert_file_enabled(convert_files)
    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )

    cycles_overridden = cli_param_overridden(ctx, "relax_max_cycles")
    _scan3d_tmp_root: Optional[Path] = None

    def _run_scan3d(
        prepared_input: Optional["PreparedInputStructure"],
        charge_val: Optional[int],
        spin_val: Optional[int],
        geom_input: Optional[Path],
        source: Optional[Path],
    ) -> None:
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
            charge=charge_val,
            spin=spin_val,
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
            source_path=source,
            freeze_links=freeze_links,
            set_charge_spin=(csv_path is None),
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
        from pdb2reaction.backends import apply_effective_precision
        apply_effective_precision(calc_cfg, precision)
        if backend_model is not None:
            from pdb2reaction.backends import apply_backend_model_to_calc_cfg
            apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        from pdb2reaction.backends import apply_calc_file_to_calc_cfg
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        if cli_param_overridden(ctx, "print_every") and print_every is not None:
            opt_cfg["print_every"] = int(print_every)
        if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
            geom_cfg["coord_type"] = str(cli_coord_type).lower()
        apply_backend_defaults(calc_cfg)

        pdb_atom_meta: List[Dict[str, Any]] = []
        d1_label_csv = None
        d2_label_csv = None
        d3_label_csv = None
        if csv_path is None:
            if source and source.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source)

            if scan_list_raw is None:
                raise click.BadParameter(
                    "--scan-lists is required.; "
                    'recover: pass --scan-lists "[(i,j,low,high),(i,j,low,high),(i,j,low,high)]" '
                    "— exactly 3 quadruples (1-based atom indices) for scan3d."
                )
            scan_one_based = bool(one_based)
            scan_source = "--scan-lists"
            if is_scan_spec_file(scan_list_raw):
                spec_path = Path(scan_list_raw)
                parsed, raw_pairs, scan_one_based = parse_scan_spec_quads(
                    spec_path,
                    expected_len=3,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
                scan_source = f"--scan-lists ({spec_path})"
            else:
                parsed, raw_pairs = parse_scan_list_quads_checked(
                    scan_list_raw,
                    expected_len=3,
                    one_based=scan_one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
            (i1, j1, low1, high1), (i2, j2, low2, high2), (i3, j3, low3, high3) = parsed
            d1_label_csv = axis_label_csv("d1", i1, j1, scan_one_based, pdb_atom_meta, raw_pairs[0])
            d2_label_csv = axis_label_csv("d2", i2, j2, scan_one_based, pdb_atom_meta, raw_pairs[1])
            d3_label_csv = axis_label_csv("d3", i3, j3, scan_one_based, pdb_atom_meta, raw_pairs[2])
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
                    {
                        "d1": (i1+1, j1+1, low1, high1),
                        "d2": (i2+1, j2+1, low2, high2),
                        "d3": (i3+1, j3+1, low3, high3),
                    },
                )
            )

            if pdb_atom_meta:
                click.echo("[scan3d] PDB atom details for scanned pairs:", detail=True)
                legend = format_pdb_atom_metadata_header()
                click.echo(f"        legend: {legend}", detail=True)
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}", detail=True)
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}", detail=True)
                click.echo(f"  d3 i: {format_pdb_atom_metadata(pdb_atom_meta, i3)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j3)}", detail=True)

        final_dir = out_dir_path

        ref_pdb_path = None
        if csv_path is None and source and source.suffix.lower() == ".pdb":
            ref_pdb_path = source

        # ==== Either load existing surface.csv, or run the full 3D scan ====
        if csv_path is not None:
            resolved_csv_path = Path(csv_path).resolve()
            try:
                df = pd.read_csv(resolved_csv_path)
            except Exception as e:
                click.echo(f"[read] Failed to read CSV '{resolved_csv_path}': {e}")
                sys.exit(1)
            click.echo(f"[read] Loaded precomputed grid from '{resolved_csv_path}'.")
        else:
            nonlocal _scan3d_tmp_root
            tmp_root = Path(tempfile.mkdtemp(prefix="scan3d_tmp_"))
            _scan3d_tmp_root = tmp_root
            grid_dir = out_dir_path / "grid"
            tmp_opt_dir = tmp_root / "opt"
            ensure_dir(grid_dir)
            ensure_dir(tmp_opt_dir)

            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            geom_outer = geom_loader(
                geom_input, coord_type=coord_type, freeze_atoms=freeze
            )
            set_freeze_atoms_or_warn(geom_outer, freeze, context="scan3d")

            base_calc = create_calculator(**calc_cfg)
            echo_resolved_device()
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

            # Optional pre-optimization of the starting structure
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

            # Measure the three bias distances on the starting structure
            # (pre-optimized when --preopt True, otherwise the input geometry)
            coords_outer = np.asarray(getattr(geom_outer, "coords3d"), dtype=float)
            d1_ref = distance_A_from_coords(coords_outer, i1, j1)
            d2_ref = distance_A_from_coords(coords_outer, i2, j2)
            d3_ref = distance_A_from_coords(coords_outer, i3, j3)
            click.echo(
                pretty_block(
                    "preopt distances (Å)",
                    {"d1_ref": d1_ref, "d2_ref": d2_ref, "d3_ref": d3_ref},
                )
            )

            # Save the starting structure (pre-optimized when requested) and prepare a record for plotting
            preopt_tag_i = distance_tag(d1_ref)
            preopt_tag_j = distance_tag(d2_ref)
            preopt_tag_k = distance_tag(d3_ref)
            preopt_xyz_path = grid_dir / f"preopt_i{preopt_tag_i}_j{preopt_tag_j}_k{preopt_tag_k}.xyz"
            try:
                s_pre = geom_outer.as_xyz()
                if not s_pre.endswith("\n"):
                    s_pre += "\n"
                with open(preopt_xyz_path, "w") as f:
                    f.write(s_pre)
                click.echo(f"[preopt] Wrote '{preopt_xyz_path}'.")
            except Exception as e:
                click.echo(f"[preopt] WARNING: failed to write '{preopt_xyz_path.name}': {e}", err=True)

            convert_xyz_like_outputs(
                preopt_xyz_path,
                prepared_input,
                ref_pdb_path=ref_pdb_path,
                out_pdb_path=grid_dir / f"preopt_i{preopt_tag_i}_j{preopt_tag_j}_k{preopt_tag_k}.pdb",
                out_gjf_path=grid_dir / f"preopt_i{preopt_tag_i}_j{preopt_tag_j}_k{preopt_tag_k}.gjf",
                context=f"'{preopt_xyz_path.name}' to PDB/GJF",
            )

            E_pre_h = unbiased_energy_hartree(geom_outer, base_calc)
            preopt_record = {
                "i": -1,
                "j": -1,
                "k": -1,
                "d1_A": float(d1_ref),
                "d2_A": float(d2_ref),
                "d3_A": float(d3_ref),
                "energy_hartree": float(E_pre_h),
                "bias_converged": True,
            }

            # Build and reorder the grids so that we scan from values closest to the reference distances
            d1_values = values_from_bounds(low1, high1, float(max_step_size))
            d2_values = values_from_bounds(low2, high2, float(max_step_size))
            d3_values = values_from_bounds(low3, high3, float(max_step_size))

            d1_values = np.array(
                sorted(d1_values, key=lambda v: abs(v - d1_ref)),
                dtype=float,
            )
            d2_values = np.array(
                sorted(d2_values, key=lambda v: abs(v - d2_ref)),
                dtype=float,
            )
            d3_values = np.array(
                sorted(d3_values, key=lambda v: abs(v - d3_ref)),
                dtype=float,
            )

            N1, N2, N3 = len(d1_values), len(d2_values), len(d3_values)
            click.echo(f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x: f'{x:.3f}', d1_values))}", narrative=True)
            click.echo(f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x: f'{x:.3f}', d2_values))}", narrative=True)
            click.echo(f"[grid] d3 steps = {N3}  values(A)={list(map(lambda x: f'{x:.3f}', d3_values))}", narrative=True)
            click.echo(f"[grid] total grid points = {N1 * N2 * N3}", narrative=True)

            max_step_bohr = float(max_step_size) * ANG2BOHR

            records: List[Dict[str, Any]] = []

            # Store starting geometry snapshot for reuse
            geom_outer_initial = _snapshot_geometry(geom_outer)

            # Caches for nearest-neighbor starting geometries
            d1_geoms: Dict[int, Any] = {}
            d2_geoms: Dict[int, Dict[int, Any]] = {}
            d3_geoms: Dict[Tuple[int, int], Dict[int, Any]] = {}

            # ===== 3D nested scan: d1 (outer) → d2 (middle) → d3 (inner) =====
            for i_idx, d1_target in enumerate(d1_values):
                click.echo(f"[stage] d1 step {i_idx + 1}/{N1}: target = {d1_target:.3f} Å", narrative=True)

                # Choose initial geometry for this d1 from the previously scanned
                # structure with the closest d1 value (or the reference structure).
                if not d1_geoms:
                    geom_outer_i = _snapshot_geometry(geom_outer_initial)
                else:
                    nearest_i = min(
                        d1_geoms.keys(),
                        key=lambda p: abs(d1_values[p] - d1_target),
                    )
                    geom_outer_i = _snapshot_geometry(d1_geoms[nearest_i])

                biased.set_pairs([(i1, j1, float(d1_target))])
                geom_outer_i.set_calculator(biased)

                opt1 = make_sopt_optimizer(
                    geom_outer_i,
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
                    click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2/d3 scan.")
                except OptimizationError as e:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {e}")

                # Snapshot after d1 relaxation for inner loops and cache
                geom_after_d1 = _snapshot_geometry(geom_outer_i)
                d1_geoms[i_idx] = geom_after_d1

                if i_idx not in d2_geoms:
                    d2_geoms[i_idx] = {}

                for j_idx, d2_target in enumerate(d2_values):
                    click.echo(
                        f"[stage] d1/d2 step ({i_idx + 1}/{N1}, {j_idx + 1}/{N2}): "
                        f"targets = ({d1_target:.3f}, {d2_target:.3f}) Å"
                    )

                    # Choose initial geometry for this (d1,d2) from the previously
                    # scanned structure with the closest d2 value at this d1
                    # (or the d1-relaxed structure).
                    d2_store = d2_geoms[i_idx]
                    if not d2_store:
                        geom_mid = _snapshot_geometry(geom_after_d1)
                    else:
                        nearest_j = min(
                            d2_store.keys(),
                            key=lambda p: abs(d2_values[p] - d2_target),
                        )
                        geom_mid = _snapshot_geometry(d2_store[nearest_j])

                    biased.set_pairs(
                        [
                            (i1, j1, float(d1_target)),
                            (i2, j2, float(d2_target)),
                        ]
                    )
                    geom_mid.set_calculator(biased)

                    opt2 = make_sopt_optimizer(
                        geom_mid,
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
                    except ZeroStepLength:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — continuing to d3 scan.")
                    except OptimizationError as e:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {e}")

                    geom_after_d2 = _snapshot_geometry(geom_mid)
                    d2_store[j_idx] = geom_after_d2

                    key_ij = (i_idx, j_idx)
                    if key_ij not in d3_geoms:
                        d3_geoms[key_ij] = {}
                    d3_store = d3_geoms[key_ij]

                    trj_blocks = [] if dump else None

                    for k_idx, d3_target in enumerate(d3_values):
                        # Choose initial geometry for this (d1,d2,d3) from the
                        # previously scanned structure with the closest d3 value
                        # at this (d1,d2), or from the d1/d2-relaxed structure.
                        if not d3_store:
                            geom_inner = _snapshot_geometry(geom_after_d2)
                        else:
                            nearest_k = min(
                                d3_store.keys(),
                                key=lambda p: abs(d3_values[p] - d3_target),
                            )
                            geom_inner = _snapshot_geometry(d3_store[nearest_k])

                        biased.set_pairs(
                            [
                                (i1, j1, float(d1_target)),
                                (i2, j2, float(d2_target)),
                                (i3, j3, float(d3_target)),
                            ]
                        )
                        geom_inner.set_calculator(biased)

                        opt3 = make_sopt_optimizer(
                            geom_inner,
                            kind,
                            lbfgs_cfg,
                            rfo_cfg,
                            opt_cfg,
                            max_step_bohr=max_step_bohr,
                            relax_max_cycles=relax_max_cycles,
                            relax_override_requested=relax_override_requested,
                            out_dir=tmp_opt_dir,
                            prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}_d3_{k_idx:03d}",
                        )
                        try:
                            opt3.run()
                            converged = True
                        except ZeroStepLength:
                            click.echo(
                                f"[d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] ZeroStepLength — recorded anyway.",
                                err=True,
                            )
                            converged = False
                        except OptimizationError as e:
                            click.echo(
                                f"[d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] OptimizationError — {e}",
                                err=True,
                            )
                            converged = False

                        # Cache final geometry for nearest-neighbor reuse
                        d3_store[k_idx] = _snapshot_geometry(geom_inner)

                        E_h = unbiased_energy_hartree(geom_inner, base_calc)

                        tag_i = distance_tag(d1_target)
                        tag_j = distance_tag(d2_target)
                        tag_k = distance_tag(d3_target)
                        xyz_path = grid_dir / f"point_i{tag_i}_j{tag_j}_k{tag_k}.xyz"
                        try:
                            s = geom_inner.as_xyz()
                            if not s.endswith("\n"):
                                s += "\n"
                            with open(xyz_path, "w") as f:
                                f.write(s)
                        except Exception as e:
                            click.echo(f"[write] WARNING: failed to write {xyz_path.name}: {e}", err=True)
                        else:
                            convert_xyz_like_outputs(
                                xyz_path,
                                prepared_input,
                                ref_pdb_path=ref_pdb_path,
                                out_pdb_path=grid_dir / f"point_i{tag_i}_j{tag_j}_k{tag_k}.pdb",
                                out_gjf_path=grid_dir / f"point_i{tag_i}_j{tag_j}_k{tag_k}.gjf",
                                context=f"'{xyz_path.name}' to PDB/GJF",
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
                                "k": int(k_idx),
                                "d1_A": float(d1_target),
                                "d2_A": float(d2_target),
                                "d3_A": float(d3_target),
                                "energy_hartree": E_h,
                                "bias_converged": bool(converged),
                            }
                        )

                    # Free d3 geometry snapshots for this (d1,d2) pair; they are
                    # only used within the d3 inner loop and won't be accessed again.
                    d3_geoms.pop(key_ij, None)

                    if dump and trj_blocks:
                        trj_path = grid_dir / f"inner_path_d1_{i_idx:03d}_d2_{j_idx:03d}_trj.xyz"
                        try:
                            with open(trj_path, "w") as f:
                                f.write("".join(trj_blocks))
                            click.echo(f"[write] Wrote '{trj_path}'.")
                        except Exception as e:
                            click.echo(f"[write] WARNING: failed to write '{trj_path}': {e}", err=True)
                        else:
                            convert_xyz_like_outputs(
                                trj_path,
                                prepared_input,
                                ref_pdb_path=ref_pdb_path,
                                out_pdb_path=grid_dir / f"inner_path_d1_{i_idx:03d}_d2_{j_idx:03d}.pdb",
                                context=f"'{trj_path.name}' to PDB",
                            )

            # Add starting structure as an extra record for plotting
            records.append(preopt_record)

            # ===== surface.csv =====
            df = pd.DataFrame.from_records(records)

        # ===== surface.csv handling & baseline =====
        if df.empty:
            click.echo("No grid records produced; aborting.")
            sys.exit(1)

        d1_label_csv = _extract_axis_label(df, "d1_label", d1_label_csv)
        d2_label_csv = _extract_axis_label(df, "d2_label", d2_label_csv)
        d3_label_csv = _extract_axis_label(df, "d3_label", d3_label_csv)

        if d1_label_csv is None or d2_label_csv is None or d3_label_csv is None:
            click.echo(
                "[plot] WARNING: axis label metadata is missing in CSV; using generic labels.",
                err=True,
            )

        d1_label_html = axis_label_html(d1_label_csv) if d1_label_csv else "d1 (Å)"
        d2_label_html = axis_label_html(d2_label_csv) if d2_label_csv else "d2 (Å)"
        d3_label_html = axis_label_html(d3_label_csv) if d3_label_csv else "d3 (Å)"

        # If energy_kcal is already present (e.g. loaded from existing CSV), reuse it.
        # Otherwise compute it from energy_hartree and baseline.
        if "energy_kcal" not in df.columns:
            if "energy_hartree" not in df.columns:
                click.echo(
                    "[baseline] energy_kcal is missing and energy_hartree is not available in CSV; aborting.",
                    err=True,
                )
                sys.exit(1)

            if baseline == "first":
                ref_mask = (df["i"] == 0) & (df["j"] == 0) & (df["k"] == 0)
                if not ref_mask.any():
                    click.echo(
                        "[baseline] 'first' requested but (i=0,j=0,k=0) missing; using global minimum instead.",
                        err=True,
                    )
                    ref = float(df["energy_hartree"].min())
                else:
                    ref = float(df.loc[ref_mask, "energy_hartree"].iloc[0])
            else:
                ref = float(df["energy_hartree"].min())

            df["energy_kcal"] = (df["energy_hartree"] - ref) * AU2KCALPERMOL

        # Only write surface.csv when we actually performed the scan in this run
        if csv_path is None:
            surface_csv = final_dir / "surface.csv"
            df["d1_label"] = d1_label_csv
            df["d2_label"] = d2_label_csv
            df["d3_label"] = d3_label_csv
            df.to_csv(surface_csv, index=False)
            click.echo(f"[write] Wrote '{surface_csv}'.")

        # ===== 3D RBF interpolation & visualization (isosurface only) =====
        d1_points = df["d1_A"].to_numpy(dtype=float)
        d2_points = df["d2_A"].to_numpy(dtype=float)
        d3_points = df["d3_A"].to_numpy(dtype=float)
        z_points = df["energy_kcal"].to_numpy(dtype=float)

        mask = (
            np.isfinite(d1_points)
            & np.isfinite(d2_points)
            & np.isfinite(d3_points)
            & np.isfinite(z_points)
        )
        if not np.any(mask):
            click.echo("[plot] No finite data for plotting.")
            sys.exit(1)

        x_min, x_max = float(np.min(d1_points[mask])), float(np.max(d1_points[mask]))
        y_min, y_max = float(np.min(d2_points[mask])), float(np.max(d2_points[mask]))
        z_min_val, z_max_val = float(np.min(d3_points[mask])), float(np.max(d3_points[mask]))

        xi = np.linspace(x_min, x_max, _VOLUME_GRID_N)
        yi = np.linspace(y_min, y_max, _VOLUME_GRID_N)
        zi = np.linspace(z_min_val, z_max_val, _VOLUME_GRID_N)

        click.echo("[plot] 3D RBF interpolation on a 50×50×50 grid ...")
        rbf3d = Rbf(
            d1_points[mask],
            d2_points[mask],
            d3_points[mask],
            z_points[mask],
            function="multiquadric",
        )

        XI, YI, ZI = np.meshgrid(xi, yi, zi, indexing="xy")
        X_flat = XI.flatten()
        Y_flat = YI.flatten()
        Z_flat = ZI.flatten()
        E_flat = rbf3d(X_flat, Y_flat, Z_flat)

        vmin = float(np.nanmin(E_flat)) if zmin is None else float(zmin)
        vmax = float(np.nanmax(E_flat)) if zmax is None else float(zmax)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = float(np.nanmin(E_flat)), float(np.nanmax(E_flat))

        # Discrete isosurfaces (mesh) with banded colors (no XY/YZ/ZX planes)
        n_levels = 8
        level_values = np.linspace(vmin, vmax, n_levels + 2)[1:-1]
        level_colors = [
            "#0d0887",
            "#5b02a3",
            "#9c179e",
            "#cb4679",
            "#ed7953",
            "#fb9f3a",
            "#fdca26",
            "#f0f921",
        ]

        # One opacity per isosurface level (outermost surfaces more transparent)
        level_opacity = [
            1.000,
            0.667,
            0.444,
            0.296,
            0.198,
            0.132,
            0.088,
            0.059,
        ]

        isosurfaces = []
        for lvl, color, opacity_lvl in zip(level_values, level_colors, level_opacity):
            trace = go.Isosurface(
                x=X_flat,
                y=Y_flat,
                z=Z_flat,
                value=E_flat,
                isomin=lvl,
                isomax=lvl,
                surface_count=1,
                opacity=opacity_lvl,
                showscale=False,
                colorscale=[[0.0, color], [1.0, color]],
                caps=dict(x_show=False, y_show=False, z_show=False),
                name=f"{lvl:.1f} kcal/mol",
            )
            isosurfaces.append(trace)

        # Add a dummy scatter trace to host a global colorbar
        colorbar_colorscale = [
            [idx / (len(level_colors) - 1), col]
            for idx, col in enumerate(level_colors)
        ]
        cb_tickvals = [float(v) for v in level_values]
        cb_ticktext = [f"{v:.1f}" for v in level_values]

        colorbar_trace = go.Scatter3d(
            x=[x_min],
            y=[y_min],
            z=[z_min_val],
            mode="markers",
            marker=dict(
                size=0,
                opacity=0.0,
                color=[vmin, vmax],
                colorscale=colorbar_colorscale,
                showscale=True,
                colorbar=dict(
                    title=dict(text="(kcal/mol)", side="top", font=dict(size=16, color="#1C1C1C")),
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
                    tickvals=cb_tickvals,
                    ticktext=cb_ticktext,
                ),
            ),
            hoverinfo="none",
            showlegend=False,
        )

        fig3d = go.Figure(data=isosurfaces + [colorbar_trace])

        fig3d.update_layout(
            title="3D Energy Landscape",
            width=900,
            height=800,
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
                    title=d3_label_html,
                    range=[z_min_val, z_max_val],
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
                aspectmode="cube",
            ),
            margin=dict(l=10, r=20, b=10, t=40),
            paper_bgcolor="white",
        )

        html3d = final_dir / "scan3d_density.html"
        fig3d.write_html(str(html3d))
        click.echo(f"[plot] Wrote '{html3d}'.")

        click.echo("\n====== 3D Scan finished ======\n", narrative=True)
        click.echo(format_elapsed("[time] Elapsed Time for 3D Scan", time_start), narrative=True)

        # result.json (if --out-json)
        if out_json:
            from pdb2reaction.core.utils import write_result_json, atom_label_from_meta
            min_energy = float(df["energy_hartree"].min()) if (not df.empty and "energy_hartree" in df.columns) else None
            result_data: Dict[str, Any] = {
                "status": "completed",
                "charge": charge_val,
                "spin": spin_val,
                "backend": calc_cfg.get("backend", backend),
                "model": calc_cfg.get("model"),
                "solvent": calc_cfg.get("solvent", "none"),
                "n_grid_points": len(df),
                "min_energy_hartree": min_energy,
                "files": {
                    "scan3d_density_html": "scan3d_density.html",
                },
            }
            if csv_path is None:
                _pair1: Dict[str, Any] = {"i": int(i1 + 1), "j": int(j1 + 1), "low": float(low1), "high": float(high1)}
                _pair2: Dict[str, Any] = {"i": int(i2 + 1), "j": int(j2 + 1), "low": float(low2), "high": float(high2)}
                _pair3: Dict[str, Any] = {"i": int(i3 + 1), "j": int(j3 + 1), "low": float(low3), "high": float(high3)}
                if pdb_atom_meta:
                    _pair1["label_i"] = atom_label_from_meta(pdb_atom_meta, i1)
                    _pair1["label_j"] = atom_label_from_meta(pdb_atom_meta, j1)
                    _pair2["label_i"] = atom_label_from_meta(pdb_atom_meta, i2)
                    _pair2["label_j"] = atom_label_from_meta(pdb_atom_meta, j2)
                    _pair3["label_i"] = atom_label_from_meta(pdb_atom_meta, i3)
                    _pair3["label_j"] = atom_label_from_meta(pdb_atom_meta, j3)
                result_data["pair1"] = _pair1
                result_data["pair2"] = _pair2
                result_data["pair3"] = _pair3
                result_data["grid_shape"] = [len(d1_values), len(d2_values), len(d3_values)]
                result_data["files"]["surface_csv"] = "surface.csv"
            write_result_json(
                final_dir, result_data,
                command="scan3d",
                elapsed_seconds=time.perf_counter() - time_start,
            )

    try:
        if csv_path is None:
            if input_path is None:
                raise click.ClickException("-i/--input is required unless --csv is provided.")
            if scan_list_raw is None:
                raise click.ClickException("--scan-lists is required unless --csv is provided.")
            with prepared_cli_input(
                input_path,
                ref_pdb=ref_pdb,
                charge=charge,
                spin=spin,
                ligand_charge=ligand_charge,
                prefix="[scan3d]",
            ) as (prepared_input, charge_val, spin_val):
                validate_charge_spin_for_prepared(prepared_input, charge_val, spin_val)
                if dry_run:
                    # Stop here. Input has been resolved, charge/spin parity
                    # has been validated against the prepared geometry, and
                    # the scan-lists literal has been parsed by Click. Print
                    # the plan and exit without burning GPU/queue time.
                    click.echo("[scan3d] --dry-run: input, charge/spin parity, and option parsing OK.")
                    click.echo(f"[scan3d] input geometry  : {prepared_input.geom_path}")
                    click.echo(f"[scan3d] resolved charge : {int(charge_val):+d}")
                    click.echo(f"[scan3d] resolved spin   : {int(spin_val)} (multiplicity)")
                    click.echo(f"[scan3d] out_dir         : {Path(out_dir).resolve()}")
                    click.echo(f"[scan3d] --scan-lists    : {scan_list_raw if scan_list_raw else '(none)'}")
                    click.echo(f"[scan3d] preopt={bool(preopt)}  freeze_links={bool(freeze_links)}")
                    click.echo("[scan3d] No 3D scan was executed.")
                    return
                _run_scan3d(
                    prepared_input,
                    charge_val,
                    spin_val,
                    prepared_input.geom_path,
                    prepared_input.source_path,
                )
        else:
            if dry_run:
                click.echo("[scan3d] --dry-run with --csv: option parsing OK.")
                click.echo(f"[scan3d] csv input  : {csv_path}")
                click.echo(f"[scan3d] out_dir    : {Path(out_dir).resolve()}")
                click.echo("[scan3d] No 3D scan was executed.")
                return
            _run_scan3d(None, charge, spin, None, None)
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        _write_error_json(Path(out_dir).resolve(), "scan3d", e, "UnhandledError", None)
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled exception during 3D scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        if _scan3d_tmp_root is not None:
            shutil.rmtree(_scan3d_tmp_root, ignore_errors=True)
