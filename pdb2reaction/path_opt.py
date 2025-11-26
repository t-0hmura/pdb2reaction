# pdb2reaction/path_opt.py

"""
path_opt — Minimum-energy path (MEP) optimization via the Growing String method (pysisyphus) with a UMA calculator
===================================================================================================================

Usage (CLI)
-----------
    pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} \
        [-q <charge>] [-m <multiplicity>] [--freeze-links {True|False}] \
        [--max-nodes <int>] [--max-cycles <int>] [--climb {True|False}] \
        [--opt-mode {light|heavy}] [--dump {True|False}] [--out-dir <dir>] \
        [--thresh <preset>] [--args-yaml <file>] [--preopt {True|False}] \
        [--preopt-max-cycles <int>] [--fix-ends {True|False}]

Examples
--------
    # Minimal: two endpoints, neutral singlet
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1

    # Typical full run with YAML overrides, dumps, and a convergence preset
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1 \
      --freeze-links True --max-nodes 10 --max-cycles 100 \
      --climb True --dump True --out-dir ./result_path_opt/ \
      --thresh gau_tight --args-yaml ./args.yaml

Description
-----------
- Optimizes a minimum-energy path between two endpoints using pysisyphus `GrowingString` and `StringOptimizer`, with UMA as the calculator (via `uma_pysis`).
- Inputs: two structures (.pdb or .xyz). If a PDB is provided and `--freeze-links=True` (default), parent atoms of link hydrogens are added to `freeze_atoms` (0-based indices).
- Configuration via YAML with sections `geom`, `calc`, `gs`, `opt`, and single-structure optimizer sections such as `sopt.lbfgs` / `sopt.rfo` (also accepting `opt.lbfgs` / `opt.rfo` and top-level `lbfgs` / `rfo`). Precedence: CLI > YAML > built-in defaults.
- Optional endpoint pre-optimization: with `--preopt=True` (default False), each endpoint is relaxed individually via single-structure LBFGS ("light") or RFO ("heavy") before alignment and GSM. The iteration limit for this pre-optimization is controlled independently by `--preopt-max-cycles` (default: 10000) for whichever optimizer is selected.
- Alignment: before optimization, all inputs after the first are rigidly Kabsch-aligned to the first structure using an external routine with a short relaxation. `StringOptimizer.align` is disabled. If either endpoint specifies `freeze_atoms`, the RMSD fit uses only those atoms and the resulting rigid transform is applied to all atoms.
- With `--climb=True` (default), a climbing-image step refines the highest-energy image. Lanczos-based tangent estimation (`gs.climb_lanczos`) is available via YAML but is disabled by default; the CLI does not toggle it.
- `--thresh` sets the convergence preset used by the string optimizer, the optional endpoint pre-optimization, and the pre-alignment refinement (e.g., `gau_loose|gau|gau_tight|gau_vtight|baker|never`).
- `--fix-ends=True` fixes both endpoint geometries during GSM (`fix_first=True`, `fix_last=True`).
- After optimization, the highest-energy image (HEI) is identified as the highest-energy internal local maximum (preferring internal nodes). If none exist, the maximum among internal nodes is used; if there are no internal nodes, the global maximum is used. The selected HEI is exported.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_path_opt/)
  ├─ final_geometries.trj        # XYZ trajectory of the optimized path; comment line carries per-image energy when available
  ├─ final_geometries.pdb        # Converted from .trj when the *first* endpoint is a PDB
  ├─ gsm_hei.xyz                 # Highest-energy image with energy on the comment line
  ├─ gsm_hei.pdb                 # HEI converted to PDB when the *first* endpoint is a PDB
  ├─ gsm_hei.gjf                 # HEI written using a detected .gjf template, when available
  ├─ align_refine/               # Files from external alignment/refinement
  └─ <optimizer dumps / restarts>  # Emitted when dumping is enabled (e.g., via `--dump` and/or YAML `dump_restart` settings)

Notes
-----
- Charge/spin: `-q/--charge` (required for non-`.gjf` inputs) and `-m/--multiplicity` are reconciled with any `.gjf` template values: explicit CLI options win; otherwise template values are used, and omitting `-q/--charge` for non-`.gjf` inputs causes the CLI to abort. Override explicitly to avoid unphysical conditions.
- Coordinates are Cartesian; `freeze_atoms` use 0-based indices. With `--freeze-links=True` and PDB inputs, link-hydrogen parents are added automatically.
- `--max-nodes` sets the number of internal nodes; the string has (max_nodes + 2) images including endpoints.
- `--max-cycles` limits optimization; after full growth, the same bound applies to additional refinement.
- `--preopt-max-cycles` limits only the optional endpoint single-structure preoptimization (LBFGS or RFO, selected via `--opt-mode`) and does not affect `--max-cycles`.
- Exit codes: 0 (success); 3 (optimization failed); 4 (final trajectory write error); 5 (HEI dump error); 130 (keyboard interrupt); 1 (unhandled error).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import sys
import traceback
import textwrap
import time

import click
import numpy as np

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .utils import (
    convert_xyz_to_pdb,
    detect_freeze_links,
    detect_freeze_links_safe,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    fill_charge_spin_from_gjf,
    maybe_convert_xyz_to_gjf,
    PreparedInputStructure,
)
from .opt import (
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
)
from .align_freeze_atoms import align_and_refine_sequence_inplace


# -----------------------------------------------
# Defaults (overridden by YAML/CLI)
# -----------------------------------------------

# Geometry (input handling)
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)

# UMA calculator settings
CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

# GrowingString (path representation)
GS_KW: Dict[str, Any] = {
    "fix_first": False,
    "fix_last": False,
    "max_nodes": 10,            # int, internal nodes; total images = max_nodes + 2 including endpoints
    "perp_thresh": 5e-3,        # float, frontier growth criterion (RMS/NORM of perpendicular force)
    "reparam_check": "rms",     # str, "rms" | "norm"; convergence check metric after reparam
    "reparam_every": 1,         # int, reparametrize every N steps while growing
    "reparam_every_full": 1,    # int, reparametrize every N steps after fully grown
    "param": "equi",            # str, "equi" (even spacing) | "energy" (weight by energy)
    "max_micro_cycles": 10,     # int, micro-optimization cycles per macro iteration
    "reset_dlc": True,          # bool, reset DLC coordinates when appropriate
    "climb": True,              # bool, enable climbing image
    "climb_rms": 5e-4,          # float, RMS force threshold to start climbing image
    "climb_lanczos": True,     # bool, use Lanczos to estimate the HEI tangent (disabled by default)
    "climb_lanczos_rms": 5e-4,  # float, RMS force threshold for Lanczos tangent
    "climb_fixed": False,       # bool, fix the HEI image index instead of adapting it
    "scheduler": None,          # Optional[str], execution scheduler; None = serial (shared calculator)
}

# StringOptimizer (optimization control)
STOPT_KW: Dict[str, Any] = {
    "type": "string",           # str, tag for bookkeeping / output labelling
    "stop_in_when_full": 300,   # int, allow N extra cycles after the string is fully grown
    "align": False,             # bool, keep internal align disabled; use external Kabsch alignment instead
    "scale_step": "global",     # str, "global" | "per_image" scaling policy
    "max_cycles": 300,          # int, maximum macro cycles for the optimizer
    "dump": False,              # bool, write optimizer trajectory to disk
    "dump_restart": False,      # bool | int, write restart YAML every N cycles (False disables)
    "reparam_thresh": 0.0,      # float, convergence threshold for reparametrization
    "coord_diff_thresh": 0.0,   # float, tolerance for coordinate difference before pruning
    "out_dir": "./result_path_opt/",  # str, output directory for optimizer artifacts
    "print_every": 10,          # int, status print frequency (cycles)
}


def _load_two_endpoints(
    inputs: Sequence[PreparedInputStructure],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> Sequence:
    """
    Load the two endpoint structures and set `freeze_atoms` as needed.
    """
    geoms = []
    for prepared in inputs:
        geom_path = prepared.geom_path
        src_path = prepared.source_path
        g = geom_loader(geom_path, coord_type=coord_type)
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        if auto_freeze_links and src_path.suffix.lower() == ".pdb":
            detected = detect_freeze_links_safe(src_path)
            freeze = merge_freeze_atom_indices(cfg, detected)
            if detected and freeze:
                click.echo(
                    f"[freeze-links] {src_path.name}: Freeze atoms (0-based): "
                    f"{','.join(map(str, freeze))}"
                )
        else:
            freeze = merge_freeze_atom_indices(cfg)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


def _maybe_convert_to_pdb(
    src_path: Path, ref_pdb_path: Optional[Path], out_path: Optional[Path] = None
) -> Optional[Path]:
    """
    Convert an XYZ/TRJ path to PDB using a reference topology when available.

    Returns the written path on success; otherwise None. Conversion is skipped when
    the source does not exist, lacks an XYZ/TRJ suffix, or no reference PDB is
    provided.
    """
    try:
        if ref_pdb_path is None:
            return None
        src_path = Path(src_path)
        if (not src_path.exists()) or src_path.suffix.lower() not in (".xyz", ".trj"):
            return None

        out_path = out_path if out_path is not None else src_path.with_suffix(".pdb")
        convert_xyz_to_pdb(src_path, ref_pdb_path, out_path)
        click.echo(f"[convert] Wrote '{out_path}'.")
        return out_path
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert '{src_path}' to PDB: {e}", err=True)
        return None


def _optimize_single(
    g,
    shared_calc,
    sopt_kind: str,
    sopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
):
    """
    Single-structure optimization (LBFGS or RFO) shared by path_opt and path_search.
    """
    g.set_calculator(shared_calc)

    seg_dir = out_dir / f"{tag}_{sopt_kind}_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(sopt_cfg)
    args["out_dir"] = str(seg_dir)

    if sopt_kind == "lbfgs":
        opt = LBFGS(g, **args)
    else:
        opt = RFOptimizer(g, **args)

    click.echo(f"\n=== [{tag}] Single-structure {sopt_kind.upper()} started ===\n")
    opt.run()
    click.echo(f"\n=== [{tag}] Single-structure {sopt_kind.upper()} finished ===\n")

    try:
        final_xyz = Path(opt.final_fn)
        _maybe_convert_to_pdb(final_xyz, ref_pdb_path)
        g_final = geom_loader(final_xyz, coord_type=g.coord_type)
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        g_final.set_calculator(shared_calc)
        return g_final
    except Exception:
        return g


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="MEP optimization via the Growing String method.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    nargs=2,
    required=True,
    help="Two endpoint structures (reactant and product); accepts .pdb or .xyz.",
)
@click.option("-q", "--charge", type=int, required=False, help="Charge of the ML region.")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=1,
    show_default=True,
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--freeze-links",
    "freeze_links_flag",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If a PDB is provided, freeze the parent atoms of link hydrogens.",
)
@click.option(
    "--max-nodes",
    type=int,
    default=10,
    show_default=True,
    help="Number of internal nodes (string has max_nodes+2 images including endpoints).",
)
@click.option(
    "--max-cycles",
    type=int,
    default=300,
    show_default=True,
    help="Maximum optimization cycles.",
)
@click.option(
    "--climb",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Search for a transition state (climbing image) after path growth.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Single-structure optimizer for endpoint preoptimization: light (=LBFGS) or heavy (=RFO).",
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Dump optimizer trajectory/restarts during the run.",
)
@click.option(
    "--out-dir",
    "out_dir",
    type=str,
    default="./result_path_opt/",
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset for the string optimizer, pre-alignment refinement, and endpoint preoptimization (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt, sopt.lbfgs, sopt.rfo).",
)
@click.option(
    "--preopt",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="If True, preoptimize each endpoint via the selected single-structure optimizer (LBFGS/RFO) before alignment and GSM.",
)
@click.option(
    "--preopt-max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum cycles for endpoint preoptimization (applies to the chosen optimizer; only used when --preopt True).",
)
@click.option(
    "--fix-ends",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Fix structures of input endpoints during GSM.",
)
def cli(
    input_paths: Sequence[Path],
    charge: Optional[int],
    spin: Optional[int],
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    preopt_max_cycles: int,
    fix_ends: bool,
) -> None:
    input_paths = tuple(Path(p) for p in input_paths)
    prepared_inputs = [prepare_input_structure(p) for p in input_paths]
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble final config (defaults ← YAML ← CLI)
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg = dict(GS_KW)
        opt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(_LBFGS_KW)
        rfo_cfg = dict(_RFO_KW)

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (gs_cfg, (("gs",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("sopt", "lbfgs"), ("opt", "lbfgs"), ("lbfgs",))),
                (rfo_cfg, (("sopt", "rfo"), ("opt", "rfo"), ("rfo",))),
            ],
        )

        # CLI overrides (highest precedence)
        resolved_charge = charge
        resolved_spin = spin
        for prepared in prepared_inputs:
            resolved_charge, resolved_spin = fill_charge_spin_from_gjf(
                resolved_charge, resolved_spin, prepared.gjf_template
            )
        if resolved_charge is None:
            resolved_charge = 0
        if resolved_spin is None:
            resolved_spin = 1
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        gs_cfg["max_nodes"] = int(max_nodes)
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)
        gs_cfg["fix_first"] = bool(fix_ends)
        gs_cfg["fix_last"] = bool(fix_ends)

        # Lanczos tangent estimation can be enabled via YAML (`gs.climb_lanczos`);
        # the CLI does not modify this setting.
        opt_cfg["dump"] = bool(dump)
        opt_cfg["out_dir"] = out_dir
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["thresh"] = str(thresh)
            rfo_cfg["thresh"] = str(thresh)

        # Use external Kabsch alignment; keep internal align disabled.
        opt_cfg["align"] = False

        lbfgs_cfg["dump"] = bool(dump)
        rfo_cfg["dump"] = bool(dump)
        lbfgs_cfg["out_dir"] = out_dir
        rfo_cfg["out_dir"] = out_dir

        opt_kind = opt_mode.strip().lower()
        if opt_kind == "light":
            sopt_kind = "lbfgs"
            sopt_cfg = lbfgs_cfg
        elif opt_kind == "heavy":
            sopt_kind = "rfo"
            sopt_cfg = rfo_cfg
        else:
            raise click.BadParameter(f"Unknown --opt-mode '{opt_mode}'.")

        sopt_cfg = dict(sopt_cfg)
        sopt_cfg["max_cycles"] = int(preopt_max_cycles)

        # For display: resolved configuration
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_gs = dict(gs_cfg)
        echo_opt = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs", echo_gs))
        click.echo(pretty_block("opt", echo_opt))
        click.echo(pretty_block("sopt." + sopt_kind, sopt_cfg))
        click.echo(
            pretty_block(
                "run_flags",
                {"preopt": bool(preopt), "preopt_max_cycles": int(preopt_max_cycles)},
            )
        )

        # --------------------------
        # 2) Prepare structures (load two endpoints and apply freezing)
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        source_paths = [prep.source_path for prep in prepared_inputs]

        geoms = _load_two_endpoints(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )

        # Shared UMA calculator (reuse the same instance for all images)
        shared_calc = uma_pysis(**calc_cfg)

        # Optional endpoint pre-optimization (LBFGS/RFO) before alignment/GSM
        if preopt:
            click.echo(
                "\n=== Preoptimizing endpoints via single-structure optimizer ===\n"
            )
            ref_pdb_for_preopt: Optional[Path] = None
            for p in source_paths:
                if p.suffix.lower() == ".pdb":
                    ref_pdb_for_preopt = p.resolve()
                    break

            new_geoms = []
            for i, g in enumerate(geoms):
                tag = f"preopt{i:02d}"
                try:
                    g_opt = _optimize_single(
                        g,
                        shared_calc,
                        sopt_kind,
                        sopt_cfg,
                        out_dir_path,
                        tag=tag,
                        ref_pdb_path=ref_pdb_for_preopt,
                    )
                    new_geoms.append(g_opt)
                except Exception as e:
                    click.echo(
                        f"[preopt] WARNING: Failed to preoptimize endpoint {i}: {e}",
                        err=True,
                    )
                    new_geoms.append(g)
            geoms = new_geoms
        else:
            click.echo(
                "[preopt] Skipping endpoint preoptimization (use --preopt True to enable)."
            )

        # External Kabsch alignment (if freeze_atoms exist, use only them)
        align_thresh = str(opt_cfg.get("thresh", "gau"))
        try:
            click.echo(
                "\n=== Aligning all inputs to the first structure "
                "(freeze-guided scan + relaxation) ===\n"
            )
            _ = align_and_refine_sequence_inplace(
                geoms,
                thresh=align_thresh,
                shared_calc=shared_calc,
                out_dir=out_dir_path / "align_refine",
                verbose=True,
            )
            click.echo("[align] Completed input alignment.")
        except Exception as e:
            click.echo(f"[align] WARNING: alignment skipped: {e}", err=True)

        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # Used when GrowingString generates new nodes
            return shared_calc

        # --------------------------
        # 3) Build path object and optimizer
        # --------------------------
        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        opt_args = dict(opt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"},
        )

        # --------------------------
        # 4) Run optimization
        # --------------------------
        click.echo("\n=== Growing String optimization started ===\n")
        optimizer.run()
        click.echo("\n=== Growing String optimization finished ===\n")

        # --------------------------
        # 5) Write final path (final_geometries.trj)
        # --------------------------
        final_trj = out_dir_path / "final_geometries.trj"
        try:
            try:
                energies = np.array(gs.energy, dtype=float)
                blocks = []
                for idx, (geom, E) in enumerate(zip(gs.images, energies)):
                    s = geom.as_xyz()
                    lines = s.splitlines()
                    if len(lines) >= 2 and lines[0].strip().isdigit():
                        lines[1] = f"{E:.12f}"
                    s_mod = "\n".join(lines)
                    if not s_mod.endswith("\n"):
                        s_mod += "\n"
                    blocks.append(s_mod)
                annotated = "".join(blocks)
                with open(final_trj, "w") as f:
                    f.write(annotated)
                click.echo(f"[write] Wrote '{final_trj}' with energy.")
            except Exception:
                with open(final_trj, "w") as f:
                    f.write(gs.as_xyz())
                click.echo(f"[write] Wrote '{final_trj}'.")

            if input_paths[0].suffix.lower() == ".pdb":
                ref_pdb = input_paths[0].resolve()

                try:
                    out_pdb = out_dir_path / "final_geometries.pdb"
                    convert_xyz_to_pdb(final_trj, ref_pdb, out_pdb)
                    click.echo(f"[convert] Wrote '{out_pdb}'.")
                except Exception as e:
                    click.echo(
                        f"[convert] WARNING: Failed to convert MEP path trajectory to PDB: {e}",
                        err=True,
                    )

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        # --------------------------
        # 6) Identify and write HEI (gsm_hei.xyz[/pdb/.gjf])
        # --------------------------
        try:
            energies = np.array(gs.energy, dtype=float)
            nE = int(len(energies))
            hei_idx = None
            if nE >= 3:
                candidates = [
                    i
                    for i in range(1, nE - 1)
                    if energies[i] > energies[i - 1] and energies[i] > energies[i + 1]
                ]
                if candidates:
                    cand_es = energies[candidates]
                    rel = int(np.argmax(cand_es))
                    hei_idx = int(candidates[rel])
                else:
                    rel = int(np.argmax(energies[1:-1]))
                    hei_idx = 1 + rel
            if hei_idx is None:
                hei_idx = int(np.argmax(energies))

            hei_geom = gs.images[hei_idx]
            hei_E = float(energies[hei_idx])

            hei_xyz = out_dir_path / "gsm_hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            ref_pdb = None
            if source_paths[0].suffix.lower() == ".pdb":
                ref_pdb = source_paths[0].resolve()
            if ref_pdb is not None:
                hei_pdb = out_dir_path / "gsm_hei.pdb"
                convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
                click.echo(f"[convert] Wrote '{hei_pdb}'.")
            else:
                click.echo(
                    "[convert] Skipped 'gsm_hei.pdb' (no PDB reference among inputs)."
                )

            template = next(
                (prep.gjf_template for prep in prepared_inputs if prep.gjf_template),
                None,
            )
            if template is not None:
                try:
                    hei_gjf = out_dir_path / "gsm_hei.gjf"
                    maybe_convert_xyz_to_gjf(hei_xyz, template, hei_gjf)
                    click.echo(f"[convert] Wrote '{hei_gjf}'.")
                except Exception as e:
                    click.echo(
                        f"[convert] WARNING: Failed to convert HEI to GJF: {e}", err=True
                    )

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

        click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start))

    except OptimizationError as e:
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo(
            "Unhandled error during path optimization:\n"
            + textwrap.indent(tb, "  "),
            err=True,
        )
        sys.exit(1)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()


def freeze_links_helper(pdb_path: Path):
    """
    Expose detect_freeze_links for external callers/tests.
    """
    return detect_freeze_links(pdb_path)
