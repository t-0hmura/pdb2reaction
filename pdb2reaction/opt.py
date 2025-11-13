# pdb2reaction/opt.py

"""
opt — Single-structure geometry optimization (LBFGS or RFO)
====================================================================

Usage (CLI)
-----
    pdb2reaction opt -i INPUT -q CHARGE [-s SPIN] [--opt-mode {light|lbfgs|heavy|rfo}] [--freeze-links/--no-freeze-links] [--dump] [--out-dir DIR] [--max-cycles N] [--args-yaml FILE]

Examples::
    pdb2reaction opt -i input.pdb -q 0
    pdb2reaction opt -i input.pdb -q 0 -s 1 --opt-mode rfo --dump --out-dir ./result_opt/ --args-yaml ./args.yaml

Description
-----
- Single-structure geometry optimization using pysisyphus with a UMA calculator.
- Input formats: .pdb, .xyz, .trj, etc., via pysisyphus `geom_loader`.
- Optimizers: LBFGS ("light") or RFOptimizer ("heavy"); aliases: light|lbfgs|heavy|rfo.
- Configuration via YAML sections `geom`, `calc`, `opt`, `lbfgs`, `rfo`. Precedence: CLI > YAML > built-in defaults.
- PDB-aware post-processing: if the input is a PDB, convert `final_geometry.xyz` → `final_geometry.pdb` and, when `--dump` is set, `optimization.trj` → `optimization.pdb` using the input PDB as the topology reference.
- Optional link-atom handling for PDBs: `--freeze-links` (default: enabled) detects link hydrogen parents and freezes those (0-based indices), merged with any `geom.freeze_atoms`.

Key options (YAML keys → meaning; defaults)
- Geometry (`geom`):
  - `coord_type`: "cart" (default) | "dlc" (often better for small molecules).
  - `freeze_atoms`: list[int], 0-based indices to freeze (default: []).
- Calculator (`calc`, UMA via `uma_pysis`):
  - `charge` (required via `-q`), `spin` (multiplicity; default 1—set explicitly for the correct state).
  - `model`: "uma-s-1p1" (default) | "uma-m-1p1"; `task_name`: "omol".
  - `device`: "auto" (GPU if available) | "cuda" | "cpu".
  - `max_neigh`: Optional[int]; `radius`: Optional[float] (Å); `r_edges`: bool.
  - `out_hess_torch`: bool; when True, provide a torch.Tensor Hessian (CUDA); else numpy on CPU.
- Optimizer base (`opt`):
  - `thresh` presets (forces in Hartree/bohr, steps in bohr):
    - `gau_loose`: max|F| 2.5e-3, RMS(F) 1.7e-3, max|step| 1.0e-2, RMS(step) 6.7e-3.
    - `gau` (default): max|F| 4.5e-4, RMS(F) 3.0e-4, max|step| 1.8e-3, RMS(step) 1.2e-3.
    - `gau_tight`: max|F| 1.5e-5, RMS(F) 1.0e-5, max|step| 6.0e-5, RMS(step) 4.0e-5.
    - `gau_vtight`: max|F| 2.0e-6, RMS(F) 1.0e-6, max|step| 6.0e-6, RMS(step) 4.0e-6.
    - `baker`: converged if (max|F| < 3.0e-4) AND (|ΔE| < 1.0e-6 OR max|step| < 3.0e-4).
    - `never`: disable built-in convergence (for external stopping).
  - `max_cycles` 10000; `print_every` 1; `min_step_norm` 1e-8 with `assert_min_step` True.
  - Convergence toggles: `rms_force`, `rms_force_only`, `max_force_only`, `force_only`.
  - Extras: `converge_to_geom_rms_thresh` (RMSD target), `overachieve_factor`, `check_eigval_structure` (TS mode checks).
  - Bookkeeping: `dump` (write `optimization.trj`), `dump_restart` (YAML every N cycles), `prefix`, `out_dir` (default `./result_opt/`).

- LBFGS-specific (`lbfgs`, used when `--opt-mode light|lbfgs`):
  - Memory: `keep_last` 7.
  - Scaling: `beta` 1.0; `gamma_mult` False.
  - Step control: `max_step` 0.30; `control_step` True; `double_damp` True; `line_search` True.
  - Regularized L‑BFGS: `mu_reg` (disables double_damp/control_step/line_search); `max_mu_reg_adaptions` 10.

- RFO-specific (`rfo`, used when `--opt-mode heavy|rfo`):
  - Trust region: `trust_radius` 0.30; `trust_update` True; bounds: `trust_min` 0.01, `trust_max` 0.30; `max_energy_incr` Optional[float].
  - Hessian: `hessian_update` "bfgs" | "bofill"; `hessian_init` "calc"; `hessian_recalc` Optional[int] (e.g., 100); `hessian_recalc_adapt` 2.0.
  - Numerics: `small_eigval_thresh` 1e-8; `line_search` True; `alpha0` 1.0; `max_micro_cycles` 25; `rfo_overlaps` False.
  - DIIS helpers: `gdiis` True; `gediis` False; `gdiis_thresh` 2.5e-3 (vs RMS(step)); `gediis_thresh` 1.0e-2 (vs RMS(force)); `gdiis_test_direction` True.
  - Step model: `adapt_step_func` False.

Outputs (& Directory Layout)
-----
- `out_dir/` (default: `./result_opt/`)
  - `final_geometry.xyz` — final optimized geometry (always).
  - `final_geometry.pdb` — converted from XYZ when the input was a PDB.
  - `optimization.trj` — trajectory (written when `--dump` or `opt.dump: true`).
  - `optimization.pdb` — converted from TRJ when input was a PDB and dumping is enabled.
  - `restart*.yml` — optional restart files every N cycles when `opt.dump_restart` is set (file names depend on optimizer).
- The tool prints resolved configuration blocks (`geom`, `calc`, `opt`, and `lbfgs` or `rfo`), progress every `print_every` cycles, and a final wall-clock time summary.

Notes:
-----
- **Required physics:** `-q/--charge` is required. Always set `-s/--spin` for the correct multiplicity (default is 1).
- **Input handling:** Supports .pdb/.xyz/.trj and other formats accepted by pysisyphus `geom_loader`. `geom.coord_type="dlc"` can improve stability for small molecules.
- **Freeze links (PDB only):** With `--freeze-links` (default), parent atoms of link hydrogens are detected and frozen; indices are 0-based and merged with `geom.freeze_atoms`.
- **PDB conversion caveat:** Conversions reuse the input PDB topology; ensure atom order/topology match the optimized coordinates.
- **Devices:** `calc.device="auto"` selects CUDA when available; otherwise CPU.
- **Hessian form:** Set `calc.out_hess_torch=True` to receive a torch/CUDA Hessian; otherwise numpy/CPU.
- **Stopping safeguards:** A `ZeroStepLength` triggers an error (min step < 1e-8). `max_energy_incr` (RFO) aborts on large uphill steps.
- **Exit codes:** 0 (success); 2 (`ZeroStepLength`); 3 (`OptimizationError`); 130 (keyboard interrupt); 1 (unhandled error).
- **Precedence:** Settings are applied with the precedence **CLI > YAML > internal defaults**.
"""

from pathlib import Path
from typing import Any, Dict, Optional

import inspect
import sys
import textwrap
import traceback

import click
import numpy as np
import yaml
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .utils import (
    convert_xyz_to_pdb,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    merge_detected_link_parents,
    merge_freeze_atom_indices,
)

# -----------------------------------------------
# Default settings (overridable via YAML/CLI)
# -----------------------------------------------

# Geometry options (YAML key: geom)
GEOM_KW = dict(GEOM_KW_DEFAULT)

# Calculator: UMA / uma_pysis  (YAML key: calc)
CALC_KW = dict(_UMA_CALC_KW)

# Optimizer base (common to LBFGS & RFO)  (YAML key: opt)
OPT_BASE_KW = {
    # Convergence threshold preset
    "thresh": "gau",             # "gau_loose" | "gau" | "gau_tight" | "gau_vtight" | "baker" | "never"

    # Convergence presets (forces in Hartree/bohr, steps in bohr)
    # +------------+------------------------------------------------------------+---------+--------+-----------+-------------+
    # |  Preset    | Purpose                                                    | max|F|  | RMS(F) | max|step| | RMS(step)   |
    # +------------+------------------------------------------------------------+---------+--------+-----------+-------------+
    # | gau_loose  | Loose/quick pre-optimization; rough path searches          | 2.5e-3  | 1.7e-3 | 1.0e-2    | 6.7e-3      |
    # | gau        | Standard "Gaussian-like" tightness for routine work        | 4.5e-4  | 3.0e-4 | 1.8e-3    | 1.2e-3      |
    # | gau_tight  | Tighter; better structures / freq / TS refinement          | 1.5e-5  | 1.0e-5 | 6.0e-5    | 4.0e-5      |
    # | gau_vtight | Very tight; benchmarking/high-precision final structures    | 2.0e-6  | 1.0e-6 | 6.0e-6    | 4.0e-6      |
    # | baker      | Baker-style rule (special condition; see below)            | 3.0e-4* | 2.0e-4 | 3.0e-4*   | 2.0e-4      |
    # | never      | Disable built-in convergence (debug/external stopping)     |   —     |   —    |    —      |    —        |
    # +------------+------------------------------------------------------------+---------+--------+-----------+-------------+
    # * Baker rule in this implementation: Converged if (max|F| < 3.0e-4) AND (|ΔE| < 1.0e-6 OR max|step| < 3.0e-4).

    "max_cycles": 10000,         # int, hard cap on optimization cycles (tool default)
    "print_every": 1,            # int, progress print frequency in cycles

    # Step-size safeguarding
    "min_step_norm": 1e-8,       # float, minimum ||step|| before raising ZeroStepLength
    "assert_min_step": True,     # bool, enforce the min_step_norm check

    # Convergence criteria toggles
    "rms_force": None,           # Optional[float], if set, derive thresholds from this RMS(F)
    "rms_force_only": False,     # bool, only check RMS(force)
    "max_force_only": False,     # bool, only check max(|force|)
    "force_only": False,         # bool, check RMS(force) and max(|force|) only

    # Extra convergence mechanisms
    "converge_to_geom_rms_thresh": 0.05,  # float, RMSD to reference geometry (Growing-NT)
    "overachieve_factor": 0.0,            # float, consider converged if forces < thresh/this_factor
    "check_eigval_structure": False,      # bool, TS search: require expected negative modes

    # Line Search
    "line_search": True,         # bool, enable polynomial line search

    # Dumping / restart / bookkeeping
    "dump": False,               # bool, write optimization trajectory
    "dump_restart": False,       # False | int, write restart YAML every N cycles (False disables)
    "prefix": "",                # str, file name prefix
    "out_dir": "./result_opt/",  # str, output directory (tool default)
}

# LBFGS-specific (YAML key: lbfgs)
LBFGS_KW = {
    **OPT_BASE_KW,

    # History / memory
    "keep_last": 7,              # int, number of (s, y) pairs to retain

    # Preconditioner / initial scaling
    "beta": 1.0,                 # float, β in -(H + βI)^{-1} g
    "gamma_mult": False,         # bool, estimate β from previous cycle (Nocedal Eq. 7.20)

    # Step-size control
    "max_step": 0.30,            # float, maximum allowed component-wise step (structural change limiter)
    "control_step": True,        # bool, scale step to satisfy |max component| <= max_step

    # Safeguards & line search
    "double_damp": True,         # bool, double-damping to enforce s·y > 0

    # Regularized L-BFGS (μ_reg)
    "mu_reg": None,              # Optional[float], initial regularization; if set: disables double_damp, control_step, line_search
    "max_mu_reg_adaptions": 10,  # int, maximum trial steps for μ adaptation
}

# RFO-specific (YAML key: rfo)
RFO_KW = {
    **OPT_BASE_KW,

    # Trust-region (step-size) control
    "trust_radius": 0.30,        # float, initial trust radius (in working coordinates)
    "trust_update": True,        # bool, adapt the trust radius based on step quality
    "trust_min": 0.01,           # float, lower bound for trust radius
    "trust_max": 0.30,           # float, upper bound for trust radius
    "max_energy_incr": None,     # Optional[float], abort if ΔE exceeds this after a bad step

    # Hessian model / refresh
    "hessian_update": "bfgs",    # "bfgs" (faster convergence) | "bofill" (more robust)
    "hessian_init": "calc",      # Initial Hessian calculation (do not change)
    "hessian_recalc": 100,       # Optional[int], recompute exact Hessian every N cycles
    "hessian_recalc_adapt": 2.0, # Heuristic: trigger exact Hessian recompute based on force norm thresholding

    # Numerical hygiene & mode filtering
    "small_eigval_thresh": 1e-8, # float, treat |λ| < threshold as zero / remove corresponding modes

    # RFO/RS micro-iterations
    "alpha0": 1.0,               # float, initial α for restricted-step RFO
    "max_micro_cycles": 25,      # int, max inner iterations to hit the trust radius
    "rfo_overlaps": False,       # bool, mode following via eigenvector overlap across cycles

    # Inter/Extrapolation helpers
    "gediis": False,             # bool, enable GEDIIS (energy-based DIIS)
    "gdiis": True,               # bool, enable GDIIS (gradient-based DIIS)

    # Thresholds for enabling DIIS (semantics matter)
    "gdiis_thresh": 2.5e-3,      # float, compared to RMS(step)  → enable GDIIS when small
    "gediis_thresh": 1.0e-2,     # float, compared to RMS(force) → enable GEDIIS when small

    "gdiis_test_direction": True,# bool, compare DIIS step direction to the RFO step

    # Choice of step model
    "adapt_step_func": False,    # bool, switch to shifted-Newton on trust when PD & gradient is small
}


# -----------------------------------------------
# Utilities
# -----------------------------------------------

def _norm_opt_mode(mode: str) -> str:
    m = (mode or "").strip().lower()
    if m in ("light", "lbfgs"):
        return "lbfgs"
    if m in ("heavy", "rfo"):
        return "rfo"
    raise click.BadParameter(f"Unknown value for --opt-mode '{mode}'. Allowed: light|lbfgs|heavy|rfo")


def _maybe_convert_outputs_to_pdb(
    input_path: Path,
    out_dir: Path,
    dump: bool,
    get_trj_fn,
    final_xyz_path: Path,
) -> None:
    """
    If the input is a PDB, convert outputs (final_geometry.xyz and, if dump, optimization.trj) to PDB.
    """
    if input_path.suffix.lower() != ".pdb":
        return

    ref_pdb = input_path.resolve()
    # final_geometry.xyz → final_geometry.pdb
    final_pdb = out_dir / "final_geometry.pdb"
    try:
        convert_xyz_to_pdb(final_xyz_path, ref_pdb, final_pdb)
        click.echo(f"[convert] Wrote '{final_pdb}'.")
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert final geometry to PDB: {e}", err=True)

    # optimization.trj → optimization.pdb (if dump)
    if dump:
        try:
            trj_path = get_trj_fn("optimization.trj")
            if trj_path.exists():
                opt_pdb = out_dir / "optimization.pdb"
                convert_xyz_to_pdb(trj_path, ref_pdb, opt_pdb)
                click.echo(f"[convert] Wrote '{opt_pdb}'.")
            else:
                click.echo("[convert] WARNING: 'optimization.trj' not found; skipping PDB conversion.", err=True)
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)
# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Single-structure geometry optimization using LBFGS or RFO.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, ...).",
)
@click.option("-q", "--charge", type=int, required=True, help="Total charge.")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1).")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze the parent atoms of link hydrogens (PDB only).")
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Maximum number of optimization cycles.")
@click.option("--opt-mode", type=str, default="light", show_default=True,
              help="Optimization mode: 'light' (=LBFGS) or 'heavy' (=RFO). Aliases: light|lbfgs|heavy|rfo.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Write optimization trajectory to 'optimization.trj'.")
@click.option("--out-dir", type=str, default="./result_opt/", show_default=True, help="Output directory.")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra arguments (sections: geom, calc, opt, lbfgs, rfo).",
)
def cli(
    input_path: Path,
    charge: int,
    spin: int,
    freeze_links: bool,
    max_cycles: int,
    opt_mode: str,
    dump: bool,
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
    time_start = time.perf_counter()

    # --------------------------
    # 1) Assemble configuration
    # --------------------------
    try:
        yaml_cfg = load_yaml_dict(args_yaml)
        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg  = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)

        # Merge YAML → defaults
        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",),)),
                (rfo_cfg, (("rfo",),)),
            ],
        )

        # CLI overrides (CLI > YAML)
        calc_cfg["charge"] = charge
        calc_cfg["spin"]   = spin
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir

        freeze = merge_freeze_atom_indices(geom_cfg)
        # Optionally infer "freeze_atoms" from link hydrogens in PDB
        if freeze_links and input_path.suffix.lower() == ".pdb":
            freeze = merge_detected_link_parents(
                geom_cfg,
                input_path,
                on_warning=lambda exc: click.echo(
                    f"[freeze-links] WARNING: Could not detect link parents: {exc}",
                    err=True,
                ),
            )
            if freeze:
                click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, freeze))}")

        # Normalize and select optimizer kind
        kind = _norm_opt_mode(opt_mode)

        # Pretty-print the resolved configuration
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
        click.echo(pretty_block("calc", calc_cfg))
        click.echo(pretty_block("opt",  {**opt_cfg, "out_dir": str(out_dir_path)}))
        click.echo(pretty_block(kind, (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))

        # --------------------------
        # 2) Prepare geometry
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        coord_type = geom_cfg.get("coord_type", "cart")
        # Pass all geometry kwargs except coord_type as coord_kwargs
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geometry = geom_loader(
            input_path,
            coord_type=coord_type,
            **coord_kwargs,
        )

        # Attach UMA calculator
        calc_builder_or_instance = uma_pysis(**calc_cfg)
        try:
            # Some wrappers return a factory; others return a ready instance
            geometry.set_calculator(calc_builder_or_instance())
        except TypeError:
            geometry.set_calculator(calc_builder_or_instance)

        # --------------------------
        # 3) Build optimizer
        # --------------------------
        common_kwargs = dict(opt_cfg)
        # Ensure paths (strings) are OK; Optimizer expects str, not Path
        common_kwargs["out_dir"] = str(out_dir_path)

        if kind == "lbfgs":
            lbfgs_args = {**lbfgs_cfg, **common_kwargs}
            # Keep only supported keys for LBFGS.__init__
            optimizer = LBFGS(geometry, **lbfgs_args)
        else:
            rfo_args = {**rfo_cfg, **common_kwargs}
            optimizer = RFOptimizer(geometry, **rfo_args)

        # --------------------------
        # 4) Run optimization
        # --------------------------
        click.echo("\n=== Optimization started ===\n")
        optimizer.run()
        click.echo("\n=== Optimization finished ===\n")

        # --------------------------
        # 5) Post-processing: PDB conversions (if input is PDB)
        # --------------------------
        # Final geometry location (Optimizer sets final_fn during run)
        final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
        _maybe_convert_outputs_to_pdb(
            input_path=input_path,
            out_dir=out_dir_path,
            dump=bool(opt_cfg["dump"]),
            get_trj_fn=optimizer.get_path_for_fn,
            final_xyz_path=final_xyz_path,
        )

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for Opt: {hh:02d}:{mm:02d}:{ss:06.3f}")

    except ZeroStepLength:
        click.echo("ERROR: Step length fell below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Optimization failed - {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled exception during optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
# Allow `python -m pdb2reaction.commands.opt` direct execution
if __name__ == "__main__":
    cli()
