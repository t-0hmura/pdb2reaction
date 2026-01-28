# pdb2reaction/scan2d.py

"""
scan2d — Two-distance 2D scan with harmonic restraints 
==================================================================

Usage (CLI)
-----------
    pdb2reaction scan2d -i INPUT.{pdb,xyz,trj,...} [-q <charge>] [--ligand-charge <number|'RES:Q,...'>] [-m <multiplicity>] \
        --scan-list(s) '[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2)]' \
        [--one-based {True|False}] \
        [--max-step-size FLOAT] \
        [--bias-k FLOAT] \
        [--relax-max-cycles INT] \
        [--opt-mode {light,heavy}] \
        [--freeze-links {True|False}] \
        [--dump {True|False}] \
        [--convert-files {True|False}] [--ref-pdb FILE] \
        [--out-dir PATH] \
        [--args-yaml FILE] \
        [--preopt {True|False}] \
        [--baseline {first|min}] \
        [--thresh {gau_loose|gau|gau_tight|gau_vtight|baker|never}] \
        [--zmin FLOAT] [--zmax FLOAT]

Examples
--------
    # Minimal example (two distance ranges)
    pdb2reaction scan2d -i input.pdb -q 0 \
        --scan-list(s) '[(12,45,1.30,3.10),(10,55,1.20,3.20)]'

    # LBFGS with trajectory dumping and PNG + HTML plots
    pdb2reaction scan2d -i input.pdb -q 0 \
        --scan-list(s) '[(12,45,1.30,3.10),(10,55,1.20,3.20)]' \
        --max-step-size 0.20 --dump True --out-dir ./result_scan2d/ --opt-mode light \
        --preopt True --baseline min

Description
-----------
- A 2D grid scan driven by harmonic restraints on two inter-atomic distances (d1, d2).
- Provide exactly one Python-like list `[(i1, j1, low1, high1), (i2, j2, low2, high2)]` via **--scan-list(s)**.
  - Indices are **1-based by default**; pass **--one-based False** to interpret them as 0-based.
  - For PDB inputs, each atom entry can be an integer index or a selector string such as
    ``'TYR,285,CA'`` or ``'MMT,309,C10'`` (resname, resseq, atom).
  - For XYZ/GJF inputs, ``--ref-pdb`` supplies a reference PDB topology while keeping the XYZ coordinates,
    enabling format-aware PDB/GJF output conversion.
- Step schedule (h = `--max-step-size` in Å):
  - `N1 = ceil(|high1 - low1| / h)`, `N2 = ceil(|high2 - low2| / h)`.
  - `d1_values = linspace(low1, high1, N1 + 1)` (or `[low1]` if the span is ~0)
    `d2_values = linspace(low2, high2, N2 + 1)` (or `[low2]` if the span is ~0).
- Optional pre-optimization:
  - If `--preopt True` (default), the initial structure is first relaxed with the unbiased UMA calculator.
    The resulting geometry is written as `preopt_i***_j***.*` in `grid/` and its energy is recorded.
- Nested scan procedure (outer d1, inner d2):
  1) For each `d1[i]`, relax with **only the d1 restraint active** (d1 bias only) and store this minimum.
  2) Each grid point `(d1[i], d2[j])` is then relaxed with **both d1 and d2 restraints active**, starting
     from the previously converged structure whose `(d1, d2)` distances are closest to the current targets.
     The **unbiased** single-point energy (harmonic bias removed for evaluation) is recorded for the surface.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_scan2d/)
  ├─ surface.csv                  # Grid metadata: i,j,d1_A,d2_A,energy_hartree,energy_kcal,bias_converged
  ├─ scan2d_map.png               # 2D contour/heatmap of the PES
  ├─ scan2d_landscape.html        # 3D surface with a base-plane projection
  └─ grid/
      ├─ point_i125_j324.xyz      # Constrained, relaxed geometries; i/j tags are 100×distance in Å (rounded)
      ├─ point_i125_j324.pdb      # Same, when the input was a PDB and conversion is enabled
      ├─ point_i125_j324.gjf      # Same, when the input provided a GJF template and conversion is enabled
      ├─ preopt_i126_j326.xyz     # Preoptimized unbiased structure; tags from its d1/d2 distances
      ├─ preopt_i126_j326.pdb     # Optional PDB version (for PDB input and conversion enabled)
      ├─ preopt_i126_j326.gjf     # Optional GJF version (for GJF input and conversion enabled)
      ├─ inner_path_d1_###.trj    # Written only when --dump True; captures inner d2 trajectories per outer step
      └─ inner_path_d1_###.pdb    # PDB conversion of inner_path_d1_###.trj when the input was a PDB and conversion is enabled

Optimizer scratch artifacts live in temporary directories; only the files above persist under ``out_dir``.

Notes
-----
- UMA only (`uma_pysis` calculator) and the same `HarmonicBiasCalculator` used in the 1D scan.
- Convergence is controlled by LBFGS or RFO depending on `--opt-mode` (`light` = LBFGS, `heavy` = RFO; default: light).
- Ångström limits are converted to Bohr to cap LBFGS step and RFO trust radii.
- The `-m/--multiplicity` option sets the spin multiplicity (2S+1) for the ML region.
- `-q/--charge` is required for non-`.gjf` inputs **unless** ``--ligand-charge`` is provided; `.gjf` templates supply
  charge/spin when available. When ``-q`` is omitted but ``--ligand-charge`` is set, the full complex is treated as an
  enzyme–substrate system and the total charge is inferred using ``extract.py``’s residue-aware logic. Explicit ``-q``
  always overrides any derived charge.
- Format-aware XYZ/TRJ → PDB/GJF conversions respect the global `--convert-files {True|False}` toggle (default: enabled).
- Convergence preset: `--thresh` defaults to `baker` for this command (override with `--thresh` or YAML).
- Cycle cap: `opt.max_cycles` defaults to `10000`; `--relax-max-cycles` overrides it only when explicitly set.
- `--baseline min|first`:
  - `min`   : shift PES so that the global minimum is 0 kcal/mol (**default**)
  - `first` : shift so that the first grid point (i=0, j=0) is 0 kcal/mol
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import math
import sys
import textwrap
import traceback
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

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .opt import (
    HarmonicBiasCalculator,
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
)
from .utils import (
    axis_label_csv,
    axis_label_html,
    build_sopt_kwargs,
    make_sopt_optimizer,
    parse_scan_list_quads_checked,
    unbiased_energy_hartree,
    values_from_bounds,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    set_convert_file_enabled,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    ensure_dir,
    make_snapshot_geometry,
    convert_xyz_like_outputs,
    build_scan_configs,
    cli_param_overridden,
    load_yaml_dict,
    resolve_freeze_atoms,
    distance_A_from_coords,
    distance_tag,
)
from .scan_common import add_scan_common_options, build_scan_kw_defaults

# Default keyword dictionaries for the 2D scan (override only the knobs we touch)
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

OPT_BASE_KW, LBFGS_KW, RFO_KW = build_scan_kw_defaults(
    opt_base_kw=_OPT_BASE_KW,
    lbfgs_kw=_LBFGS_KW,
    rfo_kw=_RFO_KW,
    out_dir="./result_scan2d/",
    thresh="baker",
)

BIAS_KW: Dict[str, Any] = {"k": 100.0}  # harmonic restraint strength (eV/Å^2)

_OPT_MODE_ALIASES = (
    (("light",), "lbfgs"),
    (("heavy",), "rfo"),
)

_snapshot_geometry = make_snapshot_geometry(GEOM_KW_DEFAULT["coord_type"])


def _sort_values_by_reference(values: np.ndarray, ref: Optional[float]) -> np.ndarray:
    """Sort scan values so that those closest to ref come first."""
    if ref is None or not np.isfinite(ref):
        return values
    order = np.argsort(np.abs(values - ref))
    return values[order]


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
    help="Input structure file (.pdb, .xyz, .trj, ...).",
)
@click.option(
    "--scan-list",
    "--scan-lists",
    "scan_list_raw",
    type=str,
    required=True,
    help="Python-like list with two quadruples: '[(i1,j1,low1,high1),(i2,j2,low2,high2)]'.",
)
@add_scan_common_options(
    workers_default=CALC_KW["workers"],
    workers_per_node_default=CALC_KW["workers_per_node"],
    out_dir_default="./result_scan2d/",
    baseline_help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0).",
    dump_help="Write inner scan trajectories per d1-step as TRJ under result_scan2d/grid/.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    scan_list_raw: str,
    one_based: bool,
    max_step_size: float,
    bias_k: float,
    relax_max_cycles: int,
    opt_mode: str,
    freeze_links: bool,
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    baseline: str,
    zmin: Optional[float],
    zmax: Optional[float],
) -> None:

    set_convert_file_enabled(convert_files)

    relax_max_cycles_override_requested = cli_param_overridden(ctx, "relax_max_cycles")

    with prepared_cli_input(
        input_path,
        ref_pdb=ref_pdb,
        charge=charge,
        spin=spin,
        ligand_charge=ligand_charge,
        prefix="[scan2d]",
    ) as (prepared_input, charge, spin):
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path

        try:
            time_start = time.perf_counter()

            yaml_cfg = load_yaml_dict(args_yaml)

            geom_cfg = dict(GEOM_KW)
            calc_cfg = dict(CALC_KW)
            opt_cfg = dict(OPT_BASE_KW)
            lbfgs_cfg = dict(LBFGS_KW)
            rfo_cfg = dict(RFO_KW)
            bias_cfg = dict(BIAS_KW)

            geom_cfg, calc_cfg, opt_cfg, lbfgs_cfg, rfo_cfg, bias_cfg = build_scan_configs(
                yaml_cfg,
                geom_kw=geom_cfg,
                calc_kw=calc_cfg,
                opt_kw=opt_cfg,
                lbfgs_kw=lbfgs_cfg,
                rfo_kw=rfo_cfg,
                bias_kw=bias_cfg,
                charge=charge,
                spin=spin,
                workers=workers,
                workers_per_node=workers_per_node,
                out_dir=out_dir,
                thresh=thresh,
                bias_k=bias_k,
            )

            kind = normalize_choice(
                opt_mode,
                param="--opt-mode",
                alias_groups=_OPT_MODE_ALIASES,
                allowed_hint="light|heavy",
            )

            # Resolve freeze list before logging so printed config matches runtime.
            freeze = resolve_freeze_atoms(geom_cfg, source_path, freeze_links)

            out_dir_path = Path(opt_cfg["out_dir"]).resolve()
            ensure_dir(out_dir_path)
            echo_geom = format_geom_for_echo(geom_cfg)
            echo_calc = format_geom_for_echo(calc_cfg)
            echo_opt = dict(opt_cfg)
            if relax_max_cycles_override_requested:
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
                relax_max_cycles_override_requested,
                out_dir_path,
                str(opt_cfg.get("prefix", "")),
            )
            click.echo(
                pretty_block(
                    "lbfgs" if kind == "lbfgs" else "rfo",
                    echo_sopt,
                )
            )
            click.echo(pretty_block("bias", echo_bias))

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            parsed, raw_pairs = parse_scan_list_quads_checked(
                scan_list_raw,
                expected_len=2,
                one_based=one_based,
                atom_meta=pdb_atom_meta,
                option_name="--scan-list",
            )
            (i1, j1, low1, high1), (i2, j2, low2, high2) = parsed
            d1_label_csv = axis_label_csv("d1", i1, j1, one_based, pdb_atom_meta, raw_pairs[0])
            d2_label_csv = axis_label_csv("d2", i2, j2, one_based, pdb_atom_meta, raw_pairs[1])
            d1_label_html = axis_label_html(d1_label_csv)
            d2_label_html = axis_label_html(d2_label_csv)
            click.echo(
                pretty_block(
                    "scan-list (0-based)",
                    {"d1": (i1, j1, low1, high1), "d2": (i2, j2, low2, high2)},
                )
            )

            if pdb_atom_meta:
                click.echo("[scan2d] PDB atom details for scanned pairs:")
                legend = format_pdb_atom_metadata_header()
                click.echo(f"        legend: {legend}")
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}")
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}")

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
            if freeze:
                try:
                    import numpy as _np

                    geom_outer.freeze_atoms = _np.array(freeze, dtype=int)
                except Exception:
                    pass

            base_calc = uma_pysis(**calc_cfg)
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
                    relax_override_requested=relax_max_cycles_override_requested,
                    out_dir=tmp_opt_dir,
                    prefix="preopt_",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.", err=True)
                except OptimizationError as e:
                    click.echo(f"[preopt] OptimizationError — {e}", err=True)

                # Measure optimized distances and record preopt structure
                try:
                    coords_outer = np.asarray(getattr(geom_outer, "coords"), dtype=float).reshape(-1, 3)
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
            click.echo(f"[grid] total grid points = {N1*N2}")

            for i_idx, d1_target in enumerate(d1_values):
                click.echo(
                    f"\n--- d1 step {i_idx+1}/{N1} : target = {d1_target:.3f} Å ---"
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
                    relax_override_requested=relax_max_cycles_override_requested,
                    out_dir=tmp_opt_dir,
                    prefix=f"d1_{i_idx:03d}_",
                )
                try:
                    opt1.run()
                except ZeroStepLength:
                    click.echo(
                        f"[d1 {i_idx}] ZeroStepLength — continuing to d2 scan.", err=True
                    )
                except OptimizationError as e:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {e}", err=True)

                geom_inner = _snapshot_geometry(geom_outer)
                geom_inner.set_calculator(biased)

                # Store the d1-relaxed structure as a candidate for nearest-start
                try:
                    coords_inner = np.asarray(getattr(geom_inner, "coords"), dtype=float).reshape(-1, 3)
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
                        relax_override_requested=relax_max_cycles_override_requested,
                        out_dir=tmp_opt_dir,
                        prefix=f"d1_{i_idx:03d}_d2_{j_idx:03d}_",
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
                            f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {e}", err=True
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
                        coords_inner = np.asarray(getattr(geom_inner, "coords"), dtype=float).reshape(-1, 3)
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

                if dump and trj_blocks:
                    trj_path = grid_dir / f"inner_path_d1_{i_idx:03d}.trj"
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
                            f"[write] WARNING: failed to write '{trj_path}': {e}", err=True
                        )

            # ===== surface.csv (final output directly under result_scan2d) =====
            df = pd.DataFrame.from_records(records)
            if df.empty:
                click.echo("No grid records produced; aborting.", err=True)
                sys.exit(1)

            if baseline == "first":
                ref = float(
                    df.loc[(df["i"] == 0) & (df["j"] == 0), "energy_hartree"].iloc[0]
                )
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
                click.echo("[plot] No finite data for plotting.", err=True)
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

            click.echo("\n=== 2D Scan finished ===\n")
            click.echo(format_elapsed("[time] Elapsed Time for 2D Scan", time_start))

        except KeyboardInterrupt:
            click.echo("\nInterrupted by user.", err=True)
            sys.exit(130)
        except Exception as e:
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            click.echo(
                "Unhandled exception during 2D scan:\n"
                + textwrap.indent(tb, "  "),
                err=True,
            )
            sys.exit(1)
