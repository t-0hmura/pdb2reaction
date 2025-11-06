# pdb2reaction/opt.py

"""
Single-structure optimization entrypoint (LBFGS or RFO).

- Supports .pdb, .xyz, etc. via pysisyphus' geom_loader.
- Calculator: UMA (via uma_pysis wrapper).
- Optimizers: LBFGS ("light") or RFOptimizer ("heavy").
- YAML config with sub-sections: geom, calc, opt, lbfgs, rfo (CLI > YAML > defaults).
- If the input is a PDB, convert outputs (final_geometry.xyz and, if dump=True,
  optimization.trj) to PDB using the reference input topology.

Example
-------
pdb2reaction opt -i input.pdb -q 0
pdb2reaction opt -i input.pdb -q 0 -s 1 --opt-mode rfo --dump --out-dir ./result_opt/ --args-yaml ./args.yaml
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

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength

from .uma_pysis import uma_pysis
from .utils import convert_xyz_to_pdb, freeze_links

# -----------------------------------------------
# Defaults variables (overridable via YAML/CLI)
# -----------------------------------------------

# Freeze atoms except indicies selected by freeze_links (yaml: geom)
GEOM_KW = {
    "coord_type"  : "cart", # Coordinate type: "cart" | "dlc". dlc is sometimies better for small molecules. 
    'freeze_atoms': []      # 0-based freeze atom indiceis
}

# Calculator: UMA / uma_pysis  (yaml: calc)
CALC_KW = {
    # Charge and multiplicity
    "charge": 0,                 # int, total charge
    "spin": 1,                   # int, multiplicity (= 2S+1); 1 for singlet

    # Model selection
    "model": "uma-s-1p1",        # str, UMA pretrained model name. "uma-s-1p1" | "uma-m-1p1"
    "task_name": "omol",         # str, UMA dataset/task tag carried into AtomicData. (Use "omol" for molecular systems.)

    # Device & graph construction
    "device": "auto",            # "cuda" | "cpu" | "auto"
    "max_neigh": None,           # Optional[int], override model's max neighbors
    "radius": None,              # Optional[float], cutoff radius (Å); default = model cutoff
    "r_edges": False,            # bool, store edge vectors in graph (UMA option)

    # Hessian output form
    "out_hess_torch": False,     # bool, if True return torch.Tensor Hessian on CUDA device; otherwise numpy on CPU
}

# Optimizer base (common to LBFGS & RFO,  (yaml: opt)
OPT_BASE_KW = {
    # Convergence threshold set
    "thresh": "gau",             # "gau_loose"|"gau"|"gau_tight"|"gau_vtight"|"baker"|"never"

    # Convergence Presets (forces: Hartree/bohr, steps : bohr)
    # +------------+--------------------------------------------------------------+---------+--------+-----------+-------------+
    # |  Preset    | Purpose                                                      | max|F|  | RMS(F) | max|step| | RMS(step)   |
    # +------------+--------------------------------------------------------------+---------+--------+-----------+-------------+
    # | gau_loose  | Loose/quick pre-optimization; rough path searches            | 2.5e-3  | 1.7e-3 | 1.0e-2    | 6.7e-3      |
    # | gau        | Standard “Gaussian-like” tightness for routine work          | 4.5e-4  | 3.0e-4 | 1.8e-3    | 1.2e-3      |
    # | gau_tight  | Tighter; for higher-quality structures / freq / TS refine    | 1.5e-5  | 1.0e-5 | 6.0e-5    | 4.0e-5      |
    # | gau_vtight | Very tight; benchmarking/high-precision final structures     | 2.0e-6  | 1.0e-6 | 6.0e-6    | 4.0e-6      |
    # | baker      | Baker-style rule (special stopping condition, see below)     | 3.0e-4* | 2.0e-4 | 3.0e-4*   | 2.0e-4      |
    # | never      | Disable auto convergence (debug/external stopping)           | —       | —      | —         | —           |
    # +------------+--------------------------------------------------------------+---------+--------+-----------+-------------+
    # * Baker rule in this implementation: Converged if ( max|F| < 3.0e-4 ) AND ( |ΔE| < 1.0e-6  OR  max|step| < 3.0e-4 ).

    "max_cycles": 10000,         # int, hard cap on optimization cycles (tool default)
    "print_every": 1,            # int, reporting cadence (cycles)

    # Step-size safeguarding
    "min_step_norm": 1e-8,       # float, min ||step|| before raising ZeroStepLength
    "assert_min_step": True,     # bool, enforce the min_step_norm check

    # Convergence criteria toggles
    "rms_force": None,           # Optional[float], if set, derive thresholds from this RMS(F)
    "rms_force_only": False,     # bool, only check RMS(force)
    "max_force_only": False,     # bool, only check max(|force|)
    "force_only": False,         # bool, check both RMS(force) and max(|force|) only

    # Extra convergence mechanisms
    "converge_to_geom_rms_thresh": 0.05,  # float, RMSD to reference geometry (Growing-NT)
    "overachieve_factor": 0.0,            # float, signal conv. if forces < thresh/this factor
    "check_eigval_structure": False,      # bool, TS-search: require desired negative modes

    # Dumping / restart / bookkeeping
    "dump": False,               # bool, write optimization trajectory
    "dump_restart": False,       # False | int, every N cycles dump restart YAML (False disables)
    "prefix": "",                # str, file name prefix
    "out_dir": "./result_opt/",  # str, output directory (tool default)
}

# LBFGS-specific (yaml: lbfgs)
LBFGS_KW = {    
    **OPT_BASE_KW,

    # History / memory
    "keep_last": 7,              # int, number of (s, y) pairs to keep

    # Preconditioner / initial scaling
    "beta": 1.0,                 # float, β in -(H + βI)^{-1} g
    "gamma_mult": False,         # bool, estimate β from prev cycle (Nocedal Eq. 7.20)

    # Step-size control
    "max_step": 0.30,            # float, max step for structual changes
    "control_step": True,        # bool, scale down step to satisfy |max component| <= max_step

    # Safeguards & line search
    "double_damp": True,         # bool, double-damping to ensure s·y > 0
    "line_search": True,         # bool, polynomial line search on last step

    # Regularized L-BFGS (μ_reg)
    "mu_reg": None,              # Optional[float], initial regularization; if set: disables double_damp, control_step, and line_search
    "max_mu_reg_adaptions": 10,  # int, max trial steps for μ adaptation
}

# RFOptimizer-specific (yaml: rfo)
RFO_KW = {
    **OPT_BASE_KW,

    # Trust-region (Step-size) control.
    "trust_radius": 0.30,         # float, initial trust radius (in working coords)
    "trust_update": True,        # bool, enable adaptive trust radius update
    "trust_min": 0.01,           # float, lower bound for trust radius
    "trust_max": 0.30,           # float, upper bound for trust radius
    "max_energy_incr": None,     # Optional[float], abort if ΔE > this after a bad step

    # Hessian model / refresh
    "hessian_update": "bfgs",    # "bfgs" (fast convergence) | "bofill" (robust)
    "hessian_init": "calc",      # Initial Hessian calculation, don't change.
    "hessian_recalc": 100,       # Optional[int], recalc exact Hessian every N cycles
    "hessian_recalc_adapt": 2.0, # If norm(force) become 1/hessian_recalc_adapt, recalc hessian

    # Numerical hygiene & mode filtering
    "small_eigval_thresh": 1e-8, # float, eigenvalues |λ| < thresh are treated as zero / removed

    # RFO/RS micro-iterations
    "line_search": True,         # bool, enable polynomial line search fallback (RFOptimizer)
    "alpha0": 1.0,               # float, initial α for restricted-step RFO
    "max_micro_cycles": 25,      # int, max inner cycles to hit trust radius
    "rfo_overlaps": False,       # bool, mode-following via eigenvector overlap across cycles

    # Inter/Extrapolation helpers
    "gediis": False,             # bool, enable GEDIIS (energy-based DIIS)
    "gdiis": True,               # bool, enable GDIIS (gradient-based DIIS)

    # Thresholds for enabling DIIS (note the semantics!)
    "gdiis_thresh": 2.5e-3,      # float, compared to RMS(step)  → enable GDIIS when small
    "gediis_thresh": 1.0e-2,     # float, compared to RMS(forces) → enable GEDIIS when small

    "gdiis_test_direction": True,# bool, compare DIIS step direction to RFO step

    # Choice of step model
    "adapt_step_func": False,    # bool, switch to shifted-Newton on trust when PD & small grad
}


# -----------------------------------------------
# Utilities
# -----------------------------------------------

def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively update dict *dst* with *src*, returning *dst*."""
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")
    return data


def _norm_opt_mode(mode: str) -> str:
    m = (mode or "").strip().lower()
    if m in ("light", "lbfgs"):
        return "lbfgs"
    if m in ("heavy", "rfo"):
        return "rfo"
    raise click.BadParameter(f"Unknown --opt-mode '{mode}'. Use: light|lbfgs|heavy|rfo")


def _maybe_convert_outputs_to_pdb(
    input_path: Path,
    out_dir: Path,
    dump: bool,
    get_trj_fn,
    final_xyz_path: Path,
) -> None:
    """If input is PDB, convert outputs (final_geometry.xyz and, if dump, optimization.trj) to PDB."""
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
                click.echo("[convert] WARNING: 'optimization.trj' not found; skip PDB conversion.", err=True)
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"

def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if fa else ""
    return g

# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Single-structure optimization (LBFGS or RFO).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, ...)",
)
@click.option("-q", "--charge", type=int, required=True, help="Total charge")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1)")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze parent atoms of link hydrogens (PDB only).")
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Max cycles")
@click.option("--opt-mode", type=str, default="light", show_default=True,
              help="light (=LBFGS) or heavy (=RFO). Aliases: light|lbfgs|heavy|rfo")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump optimization trajectory (optimization.trj)")
@click.option("--out-dir", type=str, default="./result_opt/", show_default=True, help="Output directory")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, opt, lbfgs, rfo).",
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
    # --------------------------
    # 1) Assemble configuration
    # --------------------------
    try:
        yaml_cfg = _load_yaml(args_yaml)
        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg  = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)

        # Merge YAML → defaults
        _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
        _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
        _deep_update(opt_cfg,  yaml_cfg.get("opt",  {}))
        _deep_update(lbfgs_cfg, yaml_cfg.get("lbfgs", {}))
        _deep_update(rfo_cfg,   yaml_cfg.get("rfo",   {}))

        # CLI overrides (CLI > YAML)
        calc_cfg["charge"] = charge
        calc_cfg["spin"]   = spin
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir

        # Optionally infer "freeze_atoms" from link hydrogens in PDB
        if freeze_links and input_path.suffix.lower() == ".pdb":
            try:
                detected = freeze_links_indices = freeze_links_helper(input_path)
            except NameError:
                # The helper is named freeze_links in utils.py; avoid shadowing option name
                detected = freeze_links(input_path)
            except Exception as e:
                click.echo(f"[freeze-links] WARNING: Could not detect link parents: {e}", err=True)
                detected = []
            base_freeze = list(geom_cfg.get("freeze_atoms", []))
            merged = sorted(set(base_freeze).union(set(detected)))
            geom_cfg["freeze_atoms"] = merged
            if merged:
                click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, merged))}")

        # Normalize and pick optimizer kind
        kind = _norm_opt_mode(opt_mode)

        # Pretty-print config summary (after precedence resolution)
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        click.echo(_pretty_block("geom", _format_geom_for_echo(geom_cfg)))
        click.echo(_pretty_block("calc", calc_cfg))
        click.echo(_pretty_block("opt",  {**opt_cfg, "out_dir": str(out_dir_path)}))
        click.echo(_pretty_block(kind, (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))

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
        # Ensure paths (strings) are ok, Optimizer takes str, not Path
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
        # 5) Post: PDB conversions (if input is PDB)
        # --------------------------
        # final geometry location (Optimizer final_fn is set during run)
        final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
        _maybe_convert_outputs_to_pdb(
            input_path=input_path,
            out_dir=out_dir_path,
            dump=bool(opt_cfg["dump"]),
            get_trj_fn=optimizer.get_path_for_fn,
            final_xyz_path=final_xyz_path,
        )

    except ZeroStepLength:
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


# Avoid shadowing the click option name `freeze_links` above
def freeze_links_helper(pdb_path: Path):
    """Small shim to keep the code intent readable."""
    return freeze_links(pdb_path)


# Allow `python -m pdb2reaction.commands.opt` direct execution
if __name__ == "__main__":
    cli()