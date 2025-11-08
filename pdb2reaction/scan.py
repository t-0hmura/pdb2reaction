# pdb2reaction/scan.py

"""
Bond-length driven staged scan with harmonic distance restraints and full relaxation.

Overview
--------
- Input geometry: .pdb / .xyz / ... via pysisyphus.geom_loader (Cartesian recommended).
- Calculator: UMA (via uma_pysis) wrapped by HarmonicBiasCalculator to add per-step
  harmonic distance wells toward evolving targets for specified atom pairs.
- Optimizers: LBFGS ("light") or RFOptimizer ("heavy").

Per‑stage scheduling
--------------------
For a given list of scan tuples [(i, j, target), ...]:

1) Compute each pair's displacement Δ = (target − current_distance).
2) Let d_max = max(|Δ|). With --max-step-size = h, set N = ceil(d_max / h).
3) Per-pair step width is δ_k = Δ / N.
4) At step s (1..N), the temporary target becomes r_k(s) = r_k(0) + s * δ_k.
5) Relax the full structure under the harmonic wells.

Outputs per stage (k = 1..K)
----------------------------
  stage_{k:02d}/result.xyz
  (if the input was PDB) stage_{k:02d}/result.pdb
  If --dump:
    stage_{k:02d}/scan.trj
    (if the input was PDB) stage_{k:02d}/scan.pdb

Example
-------
pdb2reaction scan -i input.pdb -q 0 --scan-lists "[(12,45,1.35)]" \
  --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
  --max-step-size 0.3 --dump True --out-dir ./result_scan/ --opt-mode rfo

Notes
-----
- Indices are 1-based by default. Use --zero-based if your tuples are 0-based.
- Units: distances in Å. The bias strength 'k' is in the calculator's
  native energy units per Å^2. Start with k=10.0 and adjust as needed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import ast
import math
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
from .utils import convert_xyz_to_pdb, freeze_links as _freeze_links_util


# --------------------------------------------------------------------------------------
# Defaults (merge order: defaults ← YAML ← CLI)
# --------------------------------------------------------------------------------------

# Geometry handling (Cartesian recommended for scans)
GEOM_KW: Dict[str, Any] = {
    "coord_type": "cart",  # "cart" | "dlc" (Cartesians recommended for scans)
    "freeze_atoms": [],    # 0-based indices to freeze (optional)
}

# UMA calculator defaults (aligned with opt.py)
CALC_KW: Dict[str, Any] = {
    "charge": 0,
    "spin": 1,                # multiplicity (= 2S+1)
    "model": "uma-s-1p1",
    "task_name": "omol",
    "device": "auto",
    "max_neigh": None,
    "radius": None,
    "r_edges": False,
    "out_hess_torch": False,
}

# Optimizer base (convergence, dumping, etc.)
OPT_BASE_KW: Dict[str, Any] = {
    "thresh": "gau",          # "gau_loose"|"gau"|"gau_tight"|"gau_vtight"|"baker"|"never"
    "max_cycles": 10000,
    "print_every": 1,
    "min_step_norm": 1e-8,
    "assert_min_step": True,
    "rms_force": None,
    "rms_force_only": False,
    "max_force_only": False,
    "force_only": False,
    "converge_to_geom_rms_thresh": 0.05,
    "overachieve_factor": 0.0,
    "check_eigval_structure": False,
    "dump": False,
    "dump_restart": False,
    "prefix": "",
    "out_dir": "./result_scan/",
}

# LBFGS specifics
LBFGS_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "keep_last": 7,
    "beta": 1.0,
    "gamma_mult": False,
    "max_step": 0.30,      # Cartesian step cap
    "control_step": True,
    "double_damp": True,
    "line_search": True,
    "mu_reg": None,
    "max_mu_reg_adaptions": 10,
}

# RFO specifics
RFO_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "trust_radius": 0.30,
    "trust_update": True,
    "trust_min": 0.01,
    "trust_max": 0.30,
    "max_energy_incr": None,
    "hessian_update": "bfgs",
    "hessian_init": "calc",
    "hessian_recalc": 100,
    "hessian_recalc_adapt": 2.0,
    "small_eigval_thresh": 1e-8,
    "line_search": True,
    "alpha0": 1.0,
    "max_micro_cycles": 25,
    "rfo_overlaps": False,
    "gediis": False,
    "gdiis": True,
    "gdiis_thresh": 2.5e-3,
    "gediis_thresh": 1.0e-2,
    "gdiis_test_direction": True,
    "adapt_step_func": False,
}

# Bias (harmonic well) defaults; can be overridden via YAML: section "bias"
BIAS_KW: Dict[str, Any] = {
    "k": 10.0,  # energy units per Å^2
}


# --------------------------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------------------------

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


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    import yaml as _yaml
    body = _yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple, np.ndarray)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if fa else ""
    return g


def _freeze_links_for_pdb(pdb_path: Path) -> List[int]:
    try:
        return list(_freeze_links_util(pdb_path))
    except Exception as e:
        click.echo(f"[freeze-links] WARNING: Could not detect link parents for '{pdb_path.name}': {e}", err=True)
        return []


def _ensure_stage_dir(base: Path, k: int) -> Path:
    d = base / f"stage_{k:02d}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _coords3d_to_xyz_string(geom, energy: Optional[float] = None) -> str:
    s = geom.as_xyz()
    lines = s.splitlines()
    if energy is not None and len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{energy:.12f}"
        s = "\n".join(lines)
    if not s.endswith("\n"):
        s += "\n"
    return s


def _parse_scan_lists(args: Sequence[str], one_based: bool) -> List[List[Tuple[int, int, float]]]:
    """
    Parse multiple Python-like list strings:
      ["[(0,1,1.5), (2,3,2.0)]", "[(5,7,1.2)]", ...]
    Returns: [[(i,j,t), ...], [(i,j,t), ...], ...] with 0-based indices.
    """
    if not args:
        raise click.BadParameter("--scan-lists must be provided at least once.")
    stages: List[List[Tuple[int, int, float]]] = []
    for idx, s in enumerate(args, start=1):
        try:
            obj = ast.literal_eval(s)
        except Exception as e:
            raise click.BadParameter(f"Invalid literal for --scan-lists #{idx}: {e}")
        if not isinstance(obj, (list, tuple)):
            raise click.BadParameter(f"--scan-lists #{idx} must be a list/tuple of (i,j,target).")
        tuples: List[Tuple[int, int, float]] = []
        for t in obj:
            if (
                isinstance(t, (list, tuple)) and len(t) == 3
                and isinstance(t[0], (int, np.integer))
                and isinstance(t[1], (int, np.integer))
                and isinstance(t[2], (int, float, np.floating))
            ):
                i, j, r = int(t[0]), int(t[1]), float(t[2])
                if one_based:
                    i -= 1
                    j -= 1
                if i < 0 or j < 0:
                    raise click.BadParameter(f"Negative atom index in --scan-lists #{idx}: {(i,j,r)} (0-based expected).")
                if r <= 0.0:
                    raise click.BadParameter(f"Non-positive target length in --scan-lists #{idx}: {(i,j,r)}.")
                tuples.append((i, j, r))
            else:
                raise click.BadParameter(f"--scan-lists #{idx} contains an invalid triple: {t}")
        stages.append(tuples)
    return stages


def _pair_distances(coords: np.ndarray, pairs: Iterable[Tuple[int, int]]) -> List[float]:
    """coords: (N,3) in Å; returns a list of distances for the given pairs."""
    dists: List[float] = []
    for i, j in pairs:
        v = coords[i] - coords[j]
        d = float(np.linalg.norm(v))
        dists.append(d)
    return dists


def _schedule_for_stage(
    coords: np.ndarray,
    tuples: List[Tuple[int, int, float]],
    max_step_size: float,
) -> Tuple[int, List[float], List[float], List[float]]:
    """
    Given current coords and stage tuples, compute:
      N: number of steps
      r0: initial distances per tuple
      rT: target distances per tuple
      step_widths: δ_k per tuple (signed)
    """
    pairs = [(i, j) for (i, j, _) in tuples]
    r0 = _pair_distances(coords, pairs)
    rT = [t for (_, _, t) in tuples]
    deltas = [RT - R0 for (R0, RT) in zip(r0, rT)]
    d_max = max((abs(d) for d in deltas), default=0.0)
    if d_max <= 0.0:
        return 0, r0, rT, [0.0] * len(tuples)
    if max_step_size <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    N = int(math.ceil(d_max / max_step_size))
    step_widths = [d / N for d in deltas]
    return N, r0, rT, step_widths


# --------------------------------------------------------------------------------------
# Harmonic bias (well) calculator wrapper
# --------------------------------------------------------------------------------------

class HarmonicBiasCalculator:
    """
    Wrap a base calculator and add harmonic distance wells:
        E_bias = sum_k 0.5 * k * (|r_i − r_j| − target_k)^2

    Supports several method names to be compatible with different calculator APIs:
      - get_energy_and_forces(geometry) -> (E, forces_flat or (N,3))
      - get_energy_and_gradient(geometry) -> (E, grad_flat)
      - get_energy(geometry) -> E
      - get_forces(geometry) -> forces_flat
      - get_gradient(geometry) -> grad_flat
    Unknown attributes are forwarded to the base calculator.
    """

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k = float(k)
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])

    # ---- control API ----
    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        # pairs: list of (i, j, target_length) with 0-based indices
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    # ---- utilities ----
    @staticmethod
    def _flatten_forces(F: np.ndarray, n_atoms: Optional[int] = None) -> np.ndarray:
        arr = np.array(F, dtype=float)
        if arr.ndim == 2 and arr.shape[1] == 3:
            return arr.reshape(-1)
        if arr.ndim == 1:
            return arr
        if n_atoms is not None:
            return arr.reshape(n_atoms * 3)
        return arr.reshape(-1)

    def _bias_energy_forces(self, geometry) -> Tuple[float, np.ndarray]:
        coords = np.array(geometry.coords3d, dtype=float)  # (N,3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        for (i, j, target) in self._pairs:
            # defensive index check
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]
            rij = float(np.linalg.norm(rij_vec))
            if rij < 1e-14:
                continue
            diff = rij - float(target)
            E_bias += 0.5 * self.k * diff * diff
            u = rij_vec / rij
            Fi = -self.k * diff * u     # -dE/dr_i
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    # ---- base calls ----
    def _base_energy_and_forces(self, geometry) -> Tuple[float, np.ndarray]:
        # Try combined methods first
        if hasattr(self.base, "get_energy_and_forces"):
            e, F = self.base.get_energy_and_forces(geometry)
            Ff = self._flatten_forces(F, getattr(geometry, "natoms", None))
            return float(e), Ff
        if hasattr(self.base, "get_energy_and_gradient"):
            e, g = self.base.get_energy_and_gradient(geometry)
            gf = np.array(g, dtype=float).reshape(-1)
            Ff = -gf
            return float(e), Ff
        # Fall back to separate
        e = None
        if hasattr(self.base, "get_energy"):
            e = float(self.base.get_energy(geometry))
        Ff = None
        if hasattr(self.base, "get_forces"):
            Ff = self._flatten_forces(self.base.get_forces(geometry), getattr(geometry, "natoms", None))
        elif hasattr(self.base, "get_gradient"):
            g = np.array(self.base.get_gradient(geometry), dtype=float).reshape(-1)
            Ff = -g
        if e is None or Ff is None:
            raise RuntimeError("Base calculator lacks required energy/forces API.")
        return e, Ff

    # ---- public API expected by Geometry/Optimizers ----
    def get_energy_and_forces(self, geometry):
        e_base, F_base = self._base_energy_and_forces(geometry)
        e_bias, F_bias = self._bias_energy_forces(geometry)
        return e_base + e_bias, F_base + F_bias

    def get_energy_and_gradient(self, geometry):
        e, F = self.get_energy_and_forces(geometry)
        grad = -np.array(F, dtype=float).reshape(-1)
        return e, grad

    def get_energy(self, geometry):
        e_base, _ = self._base_energy_and_forces(geometry)
        e_bias, _ = self._bias_energy_forces(geometry)
        return e_base + e_bias

    def get_forces(self, geometry):
        _, F = self.get_energy_and_forces(geometry)
        return np.array(F, dtype=float).reshape(-1)

    def get_gradient(self, geometry):
        _, g = self.get_energy_and_gradient(geometry)
        return np.array(g, dtype=float).reshape(-1)

    # Delegate unknown attributes to base calculator
    def __getattr__(self, name: str):
        return getattr(self.base, name)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def _norm_opt_mode(mode: str) -> str:
    m = (mode or "").strip().lower()
    if m in ("light", "lbfgs"):
        return "lbfgs"
    if m in ("heavy", "rfo"):
        return "rfo"
    raise click.BadParameter(f"Unknown --opt-mode '{mode}'. Use: light|lbfgs|heavy|rfo")


@click.command(
    help="Bond-length driven scan with staged harmonic restraints and relaxation.",
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
@click.option(
    "--scan-lists", "scan_lists_raw",
    type=str, multiple=True, required=True,
    help='Python-like list of (i,j,target) per stage. Repeatable. Example: '
         '"[(0,1,1.50),(2,3,2.00)]" "[(5,7,1.20)]".',
)
@click.option("--one-based/--zero-based", "one_based", default=True, show_default=True,
              help="Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.")
@click.option("--max-step-size", type=float, default=0.30, show_default=True,
              help="Maximum change in any scanned bond length per step [Å].")
@click.option("--bias-k", type=float, default=10.0, show_default=True,
              help="Harmonic well strength k [energy/Å^2].")
@click.option("--relax-max-cycles", type=int, default=10000, show_default=True,
              help="Maximum optimizer cycles per step.")
@click.option("--opt-mode", type=str, default="light", show_default=True,
              help="Per-step relaxation mode: light (=LBFGS) or heavy (=RFO).")
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="If input is PDB, freeze parent atoms of link hydrogens.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Write stage trajectory as scan.trj (and scan.pdb for PDB input).")
@click.option("--out-dir", type=str, default="./result_scan/", show_default=True,
              help="Base output directory.")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file with extra args (sections: geom, calc, opt, lbfgs, rfo, bias).",
)
def cli(
    input_path: Path,
    charge: int,
    spin: int,
    scan_lists_raw: Sequence[str],
    one_based: bool,
    max_step_size: float,
    bias_k: Optional[float],
    relax_max_cycles: int,
    opt_mode: str,
    freeze_links: bool,
    dump: bool,
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
    try:
        # ------------------------------------------------------------------
        # 1) Assemble configuration (defaults ← YAML ← CLI)
        # ------------------------------------------------------------------
        yaml_cfg = _load_yaml(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg  = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)
        bias_cfg  = dict(BIAS_KW)

        _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
        _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
        _deep_update(opt_cfg,  yaml_cfg.get("opt",  {}))
        _deep_update(lbfgs_cfg, yaml_cfg.get("lbfgs", {}))
        _deep_update(rfo_cfg,   yaml_cfg.get("rfo",   {}))
        _deep_update(bias_cfg,  yaml_cfg.get("bias",  {}))

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)
        opt_cfg["out_dir"] = out_dir
        # Do not use the optimizer's own dump per step; stage dumping is controlled separately.
        opt_cfg["dump"]    = False
        kind = _norm_opt_mode(opt_mode)

        # Bias strength override
        if bias_k is not None:
            bias_cfg["k"] = float(bias_k)

        # Present final config
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = _format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_opt  = dict(opt_cfg); echo_opt["out_dir"] = str(out_dir_path)
        echo_bias = dict(bias_cfg)
        click.echo(_pretty_block("geom", echo_geom))
        click.echo(_pretty_block("calc", echo_calc))
        click.echo(_pretty_block("opt",  echo_opt))
        click.echo(_pretty_block("lbfgs" if kind == "lbfgs" else "rfo", (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))
        click.echo(_pretty_block("bias", echo_bias))

        # ------------------------------------------------------------------
        # 2) Parse scan lists
        # ------------------------------------------------------------------
        stages = _parse_scan_lists(scan_lists_raw, one_based=one_based)
        K = len(stages)
        click.echo(f"[scan] Received {K} stage(s).")

        # ------------------------------------------------------------------
        # 3) Load geometry (Cartesian) and set calculator (UMA → harmonic-bias wrapper)
        # ------------------------------------------------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Load
        coord_type = geom_cfg.get("coord_type", "cart")
        # Pass only kwargs geom_loader understands; set freezes later.
        geom = geom_loader(input_path, coord_type=coord_type)

        # Merge freeze_atoms with link parents (PDB)
        freeze = list(geom_cfg.get("freeze_atoms", []))
        if freeze_links and input_path.suffix.lower() == ".pdb":
            detected = _freeze_links_for_pdb(input_path)
            if detected:
                freeze = sorted(set(freeze).union(detected))
                if freeze:
                    click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, freeze))}")
        # Attach freeze indices to Geometry for optimizer awareness
        if freeze:
            try:
                import numpy as _np
                geom.freeze_atoms = _np.array(freeze, dtype=int)
            except Exception:
                # Not all Geometry types expose freeze_atoms; ignore if absent
                pass

        # Build UMA calculator
        calc_builder_or_instance = uma_pysis(**calc_cfg)
        try:
            base_calc = calc_builder_or_instance()
        except TypeError:
            base_calc = calc_builder_or_instance

        # Wrap with bias calculator
        biased = HarmonicBiasCalculator(base_calc=base_calc, k=float(bias_cfg["k"]))
        geom.set_calculator(biased)

        # ------------------------------------------------------------------
        # 4) Stage-by-stage scan
        # ------------------------------------------------------------------
        is_pdb_input = (input_path.suffix.lower() == ".pdb")

        # Optimizer factory (per step, fresh instance)
        def _make_optimizer(kind: str, _out_dir: Path, _prefix: str):
            common = dict(opt_cfg)
            common["out_dir"] = str(_out_dir)
            common["prefix"] = _prefix
            # Tighten step control to the global --max-step-size
            if kind == "lbfgs":
                args = {**lbfgs_cfg, **common}
                # Keep the LBFGS step size conservative.
                args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), float(max_step_size))
                args["max_cycles"] = int(relax_max_cycles)
                return LBFGS(geom, **args)
            else:
                args = {**rfo_cfg, **common}
                # Trust radius no larger than --max-step-size
                tr = float(rfo_cfg.get("trust_radius", 0.30))
                args["trust_radius"] = min(tr, float(max_step_size))
                args["trust_max"] = min(float(rfo_cfg.get("trust_max", 0.30)), float(max_step_size))
                args["max_cycles"] = int(relax_max_cycles)
                return RFOptimizer(geom, **args)

        # Iterate stages
        for k, tuples in enumerate(stages, start=1):
            stage_dir = _ensure_stage_dir(out_dir_path, k)
            click.echo(f"\n--- Stage {k}/{K} ---")
            click.echo(f"Targets (i,j,target Å): {tuples}")

            # Current coordinates and schedule
            R = np.array(geom.coords3d, dtype=float)  # (N,3)
            Nsteps, r0, rT, step_widths = _schedule_for_stage(R, tuples, float(max_step_size))
            click.echo(f"[stage {k}] initial distances = {['{:.3f}'.format(x) for x in r0]}")
            click.echo(f"[stage {k}] target distances  = {['{:.3f}'.format(x) for x in rT]}")
            click.echo(f"[stage {k}] steps N = {Nsteps}")

            trj_blocks: List[str] = []

            if Nsteps == 0:
                # Nothing to do: just write current geometry as the stage result
                final_xyz = stage_dir / "result.xyz"
                with open(final_xyz, "w") as f:
                    f.write(_coords3d_to_xyz_string(geom))
                click.echo(f"[write] Wrote '{final_xyz}'.")
                if is_pdb_input:
                    try:
                        convert_xyz_to_pdb(final_xyz, input_path.resolve(), stage_dir / "result.pdb")
                        click.echo(f"[convert] Wrote '{stage_dir / 'result.pdb'}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert stage result to PDB: {e}", err=True)
                continue

            # Run N step(s)
            # Precompute r0 per tuple for linear schedule
            pairs = [(i, j) for (i, j, _) in tuples]
            for s in range(1, Nsteps + 1):
                # Compute per-pair step target for this step
                step_targets = [r0_i + s * dw for (r0_i, dw) in zip(r0, step_widths)]
                # Update bias well targets
                biased.set_pairs([(i, j, t) for ((i, j), t) in zip(pairs, step_targets)])

                # Build optimizer and relax
                prefix = f"scan_s{s:04d}_"
                optimizer = _make_optimizer(kind, stage_dir, prefix)
                click.echo(f"[stage {k}] step {s}/{Nsteps}: relaxation ({kind}) ...")
                try:
                    optimizer.run()
                except ZeroStepLength:
                    click.echo(f"[stage {k}] step {s}: ZeroStepLength — continuing to next step.", err=True)
                except OptimizationError as e:
                    click.echo(f"[stage {k}] step {s}: OptimizationError — {e}", err=True)
                    # Keep proceeding; structure is still updated to latest coords

                # Record trajectory block (optional energy annotation if available)
                try:
                    e_now = float(biased.get_energy(geom))
                except Exception:
                    e_now = None
                trj_blocks.append(_coords3d_to_xyz_string(geom, energy=e_now))

            # Stage outputs
            if dump and len(trj_blocks) > 0:
                trj_path = stage_dir / "scan.trj"
                with open(trj_path, "w") as f:
                    f.write("".join(trj_blocks))
                click.echo(f"[write] Wrote '{trj_path}'.")
                if is_pdb_input:
                    try:
                        convert_xyz_to_pdb(trj_path, input_path.resolve(), stage_dir / "scan.pdb")
                        click.echo(f"[convert] Wrote '{stage_dir / 'scan.pdb'}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert stage trajectory to PDB: {e}", err=True)

            final_xyz = stage_dir / "result.xyz"
            with open(final_xyz, "w") as f:
                f.write(_coords3d_to_xyz_string(geom))
            click.echo(f"[write] Wrote '{final_xyz}'.")

            if is_pdb_input:
                try:
                    convert_xyz_to_pdb(final_xyz, input_path.resolve(), stage_dir / "result.pdb")
                    click.echo(f"[convert] Wrote '{stage_dir / 'result.pdb'}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert stage result to PDB: {e}", err=True)

        click.echo("\n=== Scan finished ===\n")

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
