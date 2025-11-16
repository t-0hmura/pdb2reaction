# pdb2reaction/opt.py

"""
opt — Single-structure geometry optimization (LBFGS or RFO)
====================================================================

Usage (CLI)
-----
    pdb2reaction opt -i INPUT -q CHARGE [-s SPIN]
        [--opt-mode {light|lbfgs|heavy|rfo}] [--freeze-links {True|False}]
        [--dist-freeze "[(I,J,TARGET_A), ...]"] [--one-based|--zero-based] [--bias-k FLOAT]
        [--dump {True|False}] [--out-dir DIR] [--max-cycles N] [--args-yaml FILE]

Examples::
    pdb2reaction opt -i input.pdb -q 0
    pdb2reaction opt -i input.pdb -q 0 -s 1 --opt-mode rfo --dump True --out-dir ./result_opt/ --args-yaml ./args.yaml

Description
-----
- Single-structure geometry optimization using pysisyphus with a UMA calculator.
- Input formats: .pdb, .xyz, .trj, etc., via pysisyphus `geom_loader`.
- Optimizers: LBFGS ("light") or RFOptimizer ("heavy"); aliases: light|lbfgs|heavy|rfo.
- Configuration via YAML sections `geom`, `calc`, `opt`, `lbfgs`, `rfo`. Precedence: CLI > YAML > built-in defaults.
- PDB-aware post-processing: if the input is a PDB, convert `final_geometry.xyz` → `final_geometry.pdb` and, when `--dump True`, `optimization.trj` → `optimization.pdb` using the input PDB as the topology reference.
- Optional link-atom handling for PDBs: `--freeze-links True` (default) detects link hydrogen parents and freezes those (0-based indices), merged with any `geom.freeze_atoms`.
- Harmonic restraints: `--dist-freeze` accepts (i,j,target Å) tuples to apply harmonic wells during the optimization. Omit the target to freeze the initial distance; the strength is set via `--bias-k` (eV/Å²).

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
  - `max_cycles` 10000; `print_every` 100; `min_step_norm` 1e-8 with `assert_min_step` True.
  - Convergence toggles: `rms_force`, `rms_force_only`, `max_force_only`, `force_only`.
  - Extras: `converge_to_geom_rms_thresh` (RMSD target), `overachieve_factor`, `check_eigval_structure` (TS mode checks).
  - Bookkeeping: `dump` (`--dump True` writes `optimization.trj`), `dump_restart` (YAML every N cycles), `prefix`, `out_dir` (default `./result_opt/`).

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
from typing import Any, Dict, List, Optional, Sequence, Tuple

import ast
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
from pysisyphus.constants import ANG2BOHR, BOHR2ANG, AU2EV

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .utils import (
    convert_xyz_to_pdb,
    detect_freeze_links,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    normalize_choice,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    maybe_convert_xyz_to_gjf,
    GjfTemplate,
)

EV2AU = 1.0 / AU2EV  # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)

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
    "print_every": 100,            # int, progress print frequency in cycles

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

# Normalization helpers
_OPT_MODE_ALIASES = (
    (("light", "lbfgs"), "lbfgs"),
    (("heavy", "rfo"), "rfo"),
)


class HarmonicBiasCalculator:
    """Wrap a base UMA calculator with harmonic distance restraints."""

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])

    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        coords = np.array(coords_bohr, dtype=float).reshape(-1, 3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        k = self.k_au_bohr2
        for (i, j, target_ang) in self._pairs:
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]
            rij = float(np.linalg.norm(rij_vec))
            if rij < 1e-14:
                continue
            target_bohr = float(target_ang) * ANG2BOHR
            diff_bohr = rij - target_bohr
            E_bias += 0.5 * k * diff_bohr * diff_bohr
            u = rij_vec / max(rij, 1e-14)
            Fi = -k * diff_bohr * u
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_forces(elem, coords_bohr)
        E0 = float(base["energy"])
        F0 = np.asarray(base["forces"], dtype=float).reshape(-1)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias, "forces": F0 + Fbias}

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        E0 = float(self.base.get_energy(elem, coords_bohr)["energy"])
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias}

    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

    def __getattr__(self, name: str):
        return getattr(self.base, name)


def _parse_dist_freeze(
    args: Sequence[str],
    one_based: bool,
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse --dist-freeze arguments into 0-based pairs with optional targets."""

    parsed: List[Tuple[int, int, Optional[float]]] = []
    for idx, raw in enumerate(args, start=1):
        try:
            obj = ast.literal_eval(raw)
        except Exception as e:
            raise click.BadParameter(f"Invalid literal for --dist-freeze #{idx}: {e}")
        if isinstance(obj, (list, tuple)) and not obj:
            iterable = []
        elif isinstance(obj, (list, tuple)) and isinstance(obj[0], (list, tuple)):
            iterable = obj
        else:
            iterable = [obj]
        for entry in iterable:
            if not (isinstance(entry, (list, tuple)) and len(entry) in (2, 3)):
                raise click.BadParameter(
                    f"--dist-freeze #{idx} entries must be (i,j) or (i,j,target_A): {entry}"
                )
            if not (
                isinstance(entry[0], (int, np.integer))
                and isinstance(entry[1], (int, np.integer))
            ):
                raise click.BadParameter(f"Atom indices in --dist-freeze #{idx} must be integers: {entry}")
            i = int(entry[0])
            j = int(entry[1])
            target = None
            if len(entry) == 3:
                if not isinstance(entry[2], (int, float, np.floating)):
                    raise click.BadParameter(f"Target distance must be numeric in --dist-freeze #{idx}: {entry}")
                target = float(entry[2])
                if target <= 0.0:
                    raise click.BadParameter(f"Target distance must be > 0 in --dist-freeze #{idx}: {entry}")
            if one_based:
                i -= 1
                j -= 1
            if i < 0 or j < 0:
                raise click.BadParameter(f"--dist-freeze #{idx} produced negative index after conversion: {entry}")
            parsed.append((i, j, target))
    return parsed


def _resolve_dist_freeze_targets(
    geometry,
    tuples: List[Tuple[int, int, Optional[float]]],
) -> List[Tuple[int, int, float]]:
    coords_bohr = np.array(geometry.coords3d, dtype=float).reshape(-1, 3)
    coords_ang = coords_bohr * BOHR2ANG
    n = coords_ang.shape[0]
    resolved: List[Tuple[int, int, float]] = []
    for (i, j, target) in tuples:
        if not (0 <= i < n and 0 <= j < n):
            raise click.BadParameter(
                f"--dist-freeze indices {(i, j)} are out of bounds for the loaded geometry (N={n})."
            )
        if target is None:
            vec = coords_ang[i] - coords_ang[j]
            dist = float(np.linalg.norm(vec))
        else:
            dist = float(target)
        resolved.append((i, j, dist))
    return resolved


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


def _maybe_write_final_gjf(
    template: Optional[GjfTemplate],
    final_xyz_path: Path,
    out_dir: Path,
) -> None:
    if template is None:
        return
    final_gjf = out_dir / "final_geometry.gjf"
    try:
        maybe_convert_xyz_to_gjf(final_xyz_path, template, final_gjf)
        click.echo(f"[convert] Wrote '{final_gjf}'.")
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert final geometry to GJF: {e}", err=True)
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
@click.option("-q", "--charge", type=int, default=None, show_default=False, help="Total charge.")
@click.option(
    "-s",
    "--spin",
    type=int,
    default=None,
    show_default=False,
    help="Multiplicity (2S+1). Defaults to 1 when not provided.",
)
@click.option(
    "--dist-freeze",
    "dist_freeze_raw",
    type=str,
    multiple=True,
    default=(),
    show_default=False,
    help="Python-like list(s) of (i,j,target_A) to restrain distances (target optional).",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret --dist-freeze indices as 1-based (default) or 0-based.",
)
@click.option(
    "--bias-k",
    type=float,
    default=10.0,
    show_default=True,
    help="Harmonic restraint strength k [eV/Å^2] for --dist-freeze.",
)
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="Freeze the parent atoms of link hydrogens (PDB only).")
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Maximum number of optimization cycles.")
@click.option("--opt-mode", type=str, default="light", show_default=True,
              help="Optimization mode: 'light' (=LBFGS) or 'heavy' (=RFO). Aliases: light|lbfgs|heavy|rfo.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Write optimization trajectory to 'optimization.trj'.")
@click.option("--out-dir", type=str, default="./result_opt/", show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra arguments (sections: geom, calc, opt, lbfgs, rfo).",
)
def cli(
    input_path: Path,
    charge: Optional[int],
    spin: Optional[int],
    dist_freeze_raw: Sequence[str],
    one_based: bool,
    bias_k: float,
    freeze_links: bool,
    max_cycles: int,
    opt_mode: str,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
) -> None:
    time_start = time.perf_counter()
    prepared_input = prepare_input_structure(input_path)
    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

    try:
        dist_freeze = _parse_dist_freeze(dist_freeze_raw, one_based=bool(one_based))
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

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
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)

        # Optionally infer "freeze_atoms" from link hydrogens in PDB
        if freeze_links and input_path.suffix.lower() == ".pdb":
            try:
                detected = detect_freeze_links(input_path)
            except Exception as e:
                click.echo(f"[freeze-links] WARNING: Could not detect link parents: {e}", err=True)
                detected = []
            merged = merge_freeze_atom_indices(geom_cfg, detected)
            if merged:
                click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, merged))}")

        # Normalize and select optimizer kind
        kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=_OPT_MODE_ALIASES,
            allowed_hint="light|lbfgs|heavy|rfo",
        )

        # Pretty-print the resolved configuration
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
        click.echo(pretty_block("calc", calc_cfg))
        click.echo(pretty_block("opt",  {**opt_cfg, "out_dir": str(out_dir_path)}))
        click.echo(pretty_block(kind, (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))
        if dist_freeze:
            display_pairs = []
            for (i, j, target) in dist_freeze:
                label = (f"{target:.4f}" if target is not None else "<current>")
                display_pairs.append((int(i) + 1, int(j) + 1, label))
            click.echo(
                pretty_block(
                    "dist_freeze (input)",
                    {
                        "k (eV/Å^2)": float(bias_k),
                        "pairs_1based": display_pairs,
                    },
                )
            )

        # --------------------------
        # 2) Prepare geometry
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        coord_type = geom_cfg.get("coord_type", "cart")
        # Pass all geometry kwargs except coord_type as coord_kwargs
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geometry = geom_loader(
            geom_input_path,
            coord_type=coord_type,
            **coord_kwargs,
        )

        # Attach UMA calculator
        calc_builder_or_instance = uma_pysis(**calc_cfg)
        try:
            base_calc = calc_builder_or_instance()
        except TypeError:
            base_calc = calc_builder_or_instance
        geometry.set_calculator(base_calc)

        resolved_dist_freeze: List[Tuple[int, int, float]] = []
        if dist_freeze:
            try:
                resolved_dist_freeze = _resolve_dist_freeze_targets(geometry, dist_freeze)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            click.echo(
                pretty_block(
                    "dist_freeze (active)",
                    {
                        "k (eV/Å^2)": float(bias_k),
                        "pairs_1based": [
                            (int(i) + 1, int(j) + 1, float(f"{t:.4f}"))
                            for (i, j, t) in resolved_dist_freeze
                        ],
                    },
                )
            )
            bias_calc = HarmonicBiasCalculator(base_calc, k=float(bias_k))
            bias_calc.set_pairs(resolved_dist_freeze)
            geometry.set_calculator(bias_calc)

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
        _maybe_write_final_gjf(prepared_input.gjf_template, final_xyz_path, out_dir_path)

        click.echo(format_elapsed("[time] Elapsed Time for Opt", time_start))

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
    finally:
        prepared_input.cleanup()


# Avoid shadowing the click option name `freeze_links` above
def freeze_links_helper(pdb_path: Path):
    """
    Small shim to keep the intent readable.
    """
    return detect_freeze_links(pdb_path)


# Allow `python -m pdb2reaction.commands.opt` direct execution
if __name__ == "__main__":
    cli()
