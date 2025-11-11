# pdb2reaction/scan.py

"""
Bond-length driven staged scan with harmonic distance restraints and full relaxation.

**This version is simplified to support only `uma_pysis`** and removes extra
general-purpose handling to reduce overhead and speed up scans.

Overview
--------
- Input geometry: .pdb / .xyz / ... via pysisyphus.geom_loader (Cartesian recommended).
- Calculator: UMA (via uma_pysis) wrapped by a *lean* HarmonicBiasCalculator that
  adds per-step harmonic distance wells toward evolving targets for specified atom pairs.
- Optimizers: LBFGS ("light") or RFOptimizer ("heavy").

Key simplifications for speed
-----------------------------
- HarmonicBiasCalculator now implements only the (elem, coords) 2-argument API that
  pysisyphus.Geometry uses with `uma_pysis`. Geometry-style 1-arg fallbacks and broad
  dict/array shape inference were removed.
- Bias is computed in a.u. directly (Hartree, Bohr) using a pre‑converted k to
  avoid repeated unit conversions.
- Trajectory blocks are only accumulated when --dump is True.
- No attempt to re-query energy for per-frame annotation during scan (saves an extra call).

Per‑stage scheduling
--------------------
For a given list of scan tuples [(i, j, target), ...] (targets in Å):

1) Compute each pair's *Å-space* displacement Δ = (target − current_distance_Å).
2) Let d_max = max(|Δ|). With --max-step-size = h (Å), set N = ceil(d_max / h).
3) Per-pair step width is δ_k = Δ / N (Å).
4) At step s (1..N), the temporary target becomes r_k(s) = r_k(0) + s * δ_k (Å).
5) Relax the full structure under the harmonic wells.

Outputs per stage (k = 1..K)
----------------------------
  stage_{k:02d}/result.xyz
  (if the input was PDB) stage_{k:02d}/result.pdb
  If --dump:
    stage_{k:02d}/scan.trj
    (if the input was PDB) stage_{k:02d}/scan.pdb

Additional optional optimizations
---------------------------------
- With --preopt True: pre-optimize the initial structure **without bias** and continue the scan
  from that geometry. Results written to `preopt/result.xyz` (and `.pdb` if input was PDB).
- With --endopt True: after **each stage** completes its biased stepping, perform an **additional
  unbiased** geometry optimization of that stage's final structure before writing outputs.
- For each stage, echo covalent-bond **formation/breaking** between the stage's *first* structure
  and its *final* structure (the latter is the end-of-stage optimized structure when `--endopt True`).

Example
-------
pdb2reaction scan -i input.pdb -q 0 --scan-lists "[(12,45,1.35)]" \
  --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
  --max-step-size 0.2 --dump True --out-dir ./result_scan/ --opt-mode lbfgs \
  --preopt True --endopt True

Notes
-----
- Indices are 1-based by default. Use --zero-based if your tuples are 0-based.
- Units: distances in Å in the CLI/YAML. Internally, the bias is applied in a.u.
  (Hartree/Bohr) using a pre-converted k.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import ast
import math
import sys
import textwrap
import traceback
import tempfile
import os

import click
import numpy as np
import yaml
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV

EV2AU = 1.0 / AU2EV  # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)

from .uma_pysis import uma_pysis
from .utils import (
    convert_xyz_to_pdb,
    freeze_links as _freeze_links_util,
    load_yaml_dict,
    apply_yaml_overrides,
)
from .bond_changes import compare_structures, summarize_changes


# --------------------------------------------------------------------------------------
# Defaults (merge order: defaults ← YAML ← CLI)
# --------------------------------------------------------------------------------------

# Geometry handling (Cartesian recommended for scans)
GEOM_KW: Dict[str, Any] = {
    "coord_type": "cart",  # "cart" | "dlc" (Cartesians recommended for scans)
    "freeze_atoms": [],    # 0-based indices to freeze (optional)
}

# UMA calculator defaults
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
    "max_cycles": 100,
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
    "max_step": 0.30,      # Cartesian step cap (Bohr, for cart coords)
    "control_step": True,
    "double_damp": True,
    "line_search": True,
    "mu_reg": None,
    "max_mu_reg_adaptions": 10,
}

# RFO specifics
RFO_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "trust_radius": 0.30,  # Bohr
    "trust_update": True,
    "trust_min": 0.01,
    "trust_max": 0.30,     # Bohr
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
    "k": 100,  # eV / Å^2
}

# Bond-change detection (as in path_search)
BOND_KW: Dict[str, Any] = {
    "device": "cuda",
    "bond_factor": 1.20,
    "margin_fraction": 0.05,
    "delta_fraction": 0.05,
}

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


def _pair_distances(coords_ang: np.ndarray, pairs: Iterable[Tuple[int, int]]) -> List[float]:
    """coords_ang: (N,3) in Å; returns a list of distances (Å) for the given pairs."""
    dists: List[float] = []
    for i, j in pairs:
        v = coords_ang[i] - coords_ang[j]
        d = float(np.linalg.norm(v))
        dists.append(d)
    return dists


def _schedule_for_stage(
    coords_ang: np.ndarray,
    tuples: List[Tuple[int, int, float]],
    max_step_size_ang: float,
) -> Tuple[int, List[float], List[float], List[float]]:
    """
    Given current *Å* coords and stage tuples, compute:
      N: number of steps
      r0: initial distances per tuple (Å)
      rT: target distances per tuple (Å)
      step_widths: δ_k per tuple (Å, signed)
    """
    pairs = [(i, j) for (i, j, _) in tuples]
    r0 = _pair_distances(coords_ang, pairs)
    rT = [t for (_, _, t) in tuples]
    deltas = [RT - R0 for (R0, RT) in zip(r0, rT)]
    d_max = max((abs(d) for d in deltas), default=0.0)
    if d_max <= 0.0:
        return 0, r0, rT, [0.0] * len(tuples)
    if max_step_size_ang <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    N = int(math.ceil(d_max / max_step_size_ang))
    step_widths = [d / N for d in deltas]
    return N, r0, rT, step_widths


# --------------------------------------------------------------------------------------
# Bond‑change helpers
# --------------------------------------------------------------------------------------

def _has_bond_change(x, y, bond_cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """Return for covalent bonds forming/breaking between `x` and `y`."""
    res = compare_structures(
        x, y,
        device=bond_cfg.get("device", "cuda"),
        bond_factor=float(bond_cfg.get("bond_factor", 1.20)),
        margin_fraction=float(bond_cfg.get("margin_fraction", 0.05)),
        delta_fraction=float(bond_cfg.get("delta_fraction", 0.05)),
    )
    formed = len(getattr(res, "formed_covalent", [])) > 0
    broken = len(getattr(res, "broken_covalent", [])) > 0
    summary = summarize_changes(x, res, one_based=True)
    return (formed or broken), summary


def _snapshot_geometry(g) -> Any:
    """
    Create an independent pysisyphus Geometry snapshot from the given Geometry.
    Implemented via temporary XYZ serialization to avoid mutating the original.
    """
    s = g.as_xyz()
    if not s.endswith("\n"):
        s += "\n"
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()
        snap = geom_loader(Path(tmp.name), coord_type=getattr(g, "coord_type", "cart"))
        try:
            snap.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        return snap
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


# --------------------------------------------------------------------------------------
# Harmonic bias (well) calculator wrapper — optimized for uma_pysis
# --------------------------------------------------------------------------------------

class HarmonicBiasCalculator:
    """
    Wrap a base *uma_pysis* calculator and add harmonic distance wells:
        E_bias = sum_k 0.5 * k * (|r_i − r_j| − target_k)^2

    API (only the (elem, coords) call pattern used by pysisyphus.Geometry):
        get_energy(elem, coords) -> {"energy": E_total}
        get_forces(elem, coords) -> {"energy": E_total, "forces": F_total_flat}

    Units
    -----
    - Incoming coords are Bohr; base `uma_pysis` returns energy [Hartree] and forces [Hartree/Bohr].
    - `k` is provided in eV/Å^2 (CLI). Internally we convert once to Hartree/Bohr^2.
    - Targets are specified in Å and converted on the fly to Bohr.
    """

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)  # eV / Å^2
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU  # Hartree / Bohr^2
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])  # targets in Å

    # ---- control API ----
    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        # pairs: list of (i, j, target_length_in_Å) with 0-based indices
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    # ---- bias core (coords in Bohr) ----
    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        """coords_bohr: (N,3) in Bohr; returns (E_bias [Hartree], F_bias_flat [Hartree/Bohr])."""
        coords = np.array(coords_bohr, dtype=float).reshape(-1, 3)  # (N,3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        k = self.k_au_bohr2
        for (i, j, target_ang) in self._pairs:
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]               # Bohr
            rij = float(np.linalg.norm(rij_vec))          # Bohr
            if rij < 1e-14:
                continue
            target_bohr = float(target_ang) * ANG2BOHR    # Bohr
            diff_bohr = rij - target_bohr                 # Bohr
            E_bias += 0.5 * k * diff_bohr * diff_bohr     # Hartree
            u = rij_vec / max(rij, 1e-14)                 # unit vector (Bohr cancels)
            Fi = -k * diff_bohr * u                       # Hartree/Bohr
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    # ---- public API expected by Geometry/Optimizers ----
    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        # One base call provides both energy and forces in a.u.
        base = self.base.get_forces(elem, coords_bohr)
        E0 = float(base["energy"])
        F0 = np.asarray(base["forces"], dtype=float).reshape(-1)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias, "forces": F0 + Fbias}

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        # Use base energy + bias energy (no extra base call for forces)
        E0 = float(self.base.get_energy(elem, coords_bohr)["energy"])
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias}

    # Optional convenience for components that might use combined calls
    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

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
    help="Bond-length driven scan with staged harmonic restraints and relaxation (UMA only).",
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
@click.option("--max-step-size", type=float, default=0.20, show_default=True,
              help="Maximum change in any scanned bond length per step [Å].")
@click.option("--bias-k", type=float, default=100, show_default=True,
              help="Harmonic well strength k [eV/Å^2].")
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
    help="YAML file with extra args (sections: geom, calc, opt, lbfgs, rfo, bias, bond).",
)
@click.option("--preopt", type=click.BOOL, default=True, show_default=True,
              help="Pre-optimize initial structure without bias before the scan.")
@click.option("--endopt", type=click.BOOL, default=True, show_default=True,
              help="After each stage, run an additional unbiased optimization of the stage result.")
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
    preopt: bool,
    endopt: bool,
) -> None:
    try:
        time_start = time.perf_counter()  # start timing
        
        # ------------------------------------------------------------------
        # 1) Assemble configuration (defaults ← YAML ← CLI)
        # ------------------------------------------------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg  = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)
        bias_cfg  = dict(BIAS_KW)
        bond_cfg  = dict(BOND_KW)  # <-- added

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",),)),
                (rfo_cfg, (("rfo",),)),
                (bias_cfg, (("bias",),)),
                (bond_cfg, (("bond",),)),
            ],
        )

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
        echo_bond = dict(bond_cfg)  # <-- added
        click.echo(_pretty_block("geom", echo_geom))
        click.echo(_pretty_block("calc", echo_calc))
        click.echo(_pretty_block("opt",  echo_opt))
        click.echo(_pretty_block("lbfgs" if kind == "lbfgs" else "rfo", (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))
        click.echo(_pretty_block("bias", echo_bias))
        click.echo(_pretty_block("bond", echo_bond))  # <-- added

        # ------------------------------------------------------------------
        # 2) Parse scan lists
        # ------------------------------------------------------------------
        stages = _parse_scan_lists(scan_lists_raw, one_based=one_based)
        K = len(stages)
        click.echo(f"[scan] Received {K} stage(s).")

        # Prepare end-of-run summary collector (minimal additions)
        stages_summary: List[Dict[str, Any]] = []

        # ------------------------------------------------------------------
        # 3) Load geometry (Cartesian) and set calculator (UMA → harmonic-bias wrapper)
        # ------------------------------------------------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Load
        coord_type = geom_cfg.get("coord_type", "cart")
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
                pass

        # Build UMA calculator (only uma_pysis is supported)
        base_calc = uma_pysis(**calc_cfg)

        # Input type flag (used for optional PDB conversion)
        is_pdb_input = (input_path.suffix.lower() == ".pdb")

        # ------------------------------------------------------------------
        # Optional pre-optimization WITHOUT bias
        # ------------------------------------------------------------------
        if preopt:
            pre_dir = out_dir_path / "preopt"
            pre_dir.mkdir(parents=True, exist_ok=True)
            geom.set_calculator(base_calc)
            click.echo(f"[preopt] Unbiased relaxation ({kind}) ...")
            # Local optimizer factory for preopt
            max_step_bohr_local = float(max_step_size) * ANG2BOHR  # ensure step controls are in Bohr
            def _make_optimizer_pre(kind_local: str, _out_dir: Path, _prefix: str):
                common = dict(opt_cfg)
                common["out_dir"] = str(_out_dir)
                common["prefix"] = _prefix
                if kind_local == "lbfgs":
                    args = {**lbfgs_cfg, **common}
                    args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr_local)
                    args["max_cycles"] = int(relax_max_cycles)
                    return LBFGS(geom, **args)
                else:
                    args = {**rfo_cfg, **common}
                    tr = float(rfo_cfg.get("trust_radius", 0.30))
                    args["trust_radius"] = min(tr, max_step_bohr_local)
                    args["trust_max"] = min(float(rfo_cfg.get("trust_max", 0.30)), max_step_bohr_local)
                    args["max_cycles"] = int(relax_max_cycles)
                    return RFOptimizer(geom, **args)
            optimizer0 = _make_optimizer_pre(kind, pre_dir, "preopt_")
            try:
                optimizer0.run()
            except ZeroStepLength:
                click.echo(f"[preopt] ZeroStepLength — continuing.", err=True)
            except OptimizationError as e:
                click.echo(f"[preopt] OptimizationError — {e}", err=True)

            # Write preopt result
            pre_xyz = pre_dir / "result.xyz"
            with open(pre_xyz, "w") as f:
                f.write(_coords3d_to_xyz_string(geom))
            click.echo(f"[write] Wrote '{pre_xyz}'.")
            if is_pdb_input:
                try:
                    convert_xyz_to_pdb(pre_xyz, input_path.resolve(), pre_dir / "result.pdb")
                    click.echo(f"[convert] Wrote '{pre_dir / 'result.pdb'}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert preopt result to PDB: {e}", err=True)

        # Wrap with bias calculator for the scan
        biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))
        geom.set_calculator(biased)

        # ------------------------------------------------------------------
        # 4) Stage-by-stage scan
        # ------------------------------------------------------------------

        # Optimizer factory (per step or end-of-stage, fresh instance)
        max_step_bohr = float(max_step_size) * ANG2BOHR  # ensure step controls are in Bohr
        def _make_optimizer(kind: str, _out_dir: Path, _prefix: str):
            common = dict(opt_cfg)
            common["out_dir"] = str(_out_dir)
            common["prefix"] = _prefix
            if kind == "lbfgs":
                args = {**lbfgs_cfg, **common}
                # Cap LBFGS step size in Bohr by the global Å limit converted to Bohr.
                args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
                args["max_cycles"] = int(relax_max_cycles)
                return LBFGS(geom, **args)
            else:
                args = {**rfo_cfg, **common}
                tr = float(rfo_cfg.get("trust_radius", 0.30))
                args["trust_radius"] = min(tr, max_step_bohr)
                args["trust_max"] = min(float(rfo_cfg.get("trust_max", 0.30)), max_step_bohr)
                args["max_cycles"] = int(relax_max_cycles)
                return RFOptimizer(geom, **args)

        # Iterate stages
        for k, tuples in enumerate(stages, start=1):
            stage_dir = _ensure_stage_dir(out_dir_path, k)
            click.echo(f"\n--- Stage {k}/{K} ---")
            click.echo(f"Targets (i,j,target Å): {tuples}")

            # Snapshot **beginning** geometry of this stage for bond-change comparison
            start_geom_for_stage = _snapshot_geometry(geom)

            # Current coordinates (Bohr) and schedule computed in Å
            R_bohr = np.array(geom.coords3d, dtype=float)      # (N,3) Bohr
            R_ang  = R_bohr * BOHR2ANG                         # (N,3) Å
            Nsteps, r0, rT, step_widths = _schedule_for_stage(R_ang, tuples, float(max_step_size))
            click.echo(f"[stage {k}] initial distances (Å) = {['{:.3f}'.format(x) for x in r0]}")
            click.echo(f"[stage {k}] target distances  (Å) = {['{:.3f}'.format(x) for x in rT]}")
            click.echo(f"[stage {k}] steps N = {Nsteps}")

            # ---- record per-stage summary (for final echo) ----
            srec: Dict[str, Any] = {
                "index": int(k),
                "pairs_1based": [(int(i)+1, int(j)+1) for (i, j, _) in tuples],
                "initial_distances_A": [float(f"{x:.3f}") for x in r0],
                "target_distances_A": [float(f"{x:.3f}") for x in rT],
                "per_pair_step_A": [float(f"{x:.3f}") for x in step_widths],
                "num_steps": int(Nsteps),
                "bond_change": {"changed": None, "summary": ""},
            }
            stages_summary.append(srec)
            # ---------------------------------------------------

            trj_blocks: List[str] = [] if dump else None

            pairs = [(i, j) for (i, j, _) in tuples]

            if Nsteps == 0:
                # No stepping; optionally perform end-of-stage unbiased optimization
                if endopt:
                    geom.set_calculator(base_calc)
                    click.echo(f"[stage {k}] endopt (unbiased) ...")
                    try:
                        end_optimizer = _make_optimizer(kind, stage_dir, "endopt_")
                        end_optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)

                # ---- Echo bond changes: start vs final (possibly endopt) ----
                try:
                    changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                    click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                    if changed and summary and summary.strip():
                        click.echo(textwrap.indent(summary.strip(), prefix="  "))
                    if not changed:
                        click.echo("  (no covalent changes detected)")
                    # record to summary
                    try:
                        srec["bond_change"]["changed"] = bool(changed)
                        srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                    except Exception:
                        pass
                except Exception as e:
                    click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)
                # -------------------------------------------------------------

                # Write current (possibly endopted) geometry as the stage result
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

            # Run N step(s) with bias
            for s in range(1, Nsteps + 1):
                # Compute per-pair step target (Å) for this step
                step_targets = [r0_i + s * dw for (r0_i, dw) in zip(r0, step_widths)]

                # Update bias well targets (still in Å; wrapper converts internally)
                biased.set_pairs([(i, j, t) for ((i, j), t) in zip(pairs, step_targets)])
                # IMPORTANT: only the calculator's internal state changed -> flush Geometry caches
                # Re-attaching the same calculator marks the geometry "dirty" so the next call recomputes E/F
                geom.set_calculator(biased)

                # Build optimizer and relax (with bias)
                prefix = f"scan_s{s:04d}_"
                optimizer = _make_optimizer(kind, stage_dir, prefix)
                click.echo(f"[stage {k}] step {s}/{Nsteps}: relaxation ({kind}) ...")
                try:
                    optimizer.run()
                except ZeroStepLength:
                    click.echo(f"[stage {k}] step {s}: ZeroStepLength — continuing to next step.", err=True)
                except OptimizationError as e:
                    click.echo(f"[stage {k}] step {s}: OptimizationError — {e}", err=True)

                # Record trajectory block only when requested (biased result)
                if dump and trj_blocks is not None:
                    trj_blocks.append(_coords3d_to_xyz_string(geom))

            # Optional end-of-stage UNBIASED optimization
            if endopt:
                geom.set_calculator(base_calc)
                click.echo(f"[stage {k}] endopt (unbiased) ...")
                try:
                    end_optimizer = _make_optimizer(kind, stage_dir, "endopt_")
                    end_optimizer.run()
                except ZeroStepLength:
                    click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                except OptimizationError as e:
                    click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)

            # ---- Echo bond changes: start vs final (possibly endopt) ----
            try:
                changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                if changed and summary and summary.strip():
                    click.echo(textwrap.indent(summary.strip(), prefix="  "))
                if not changed:
                    click.echo("  (no covalent changes detected)")
                # record to summary
                try:
                    srec["bond_change"]["changed"] = bool(changed)
                    srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                except Exception:
                    pass
            except Exception as e:
                click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)
            # -------------------------------------------------------------

            # Stage outputs
            if dump and trj_blocks:
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

        # ------------------------------------------------------------------
        # 5) Final summary echo (human‑friendly)
        # ------------------------------------------------------------------
        def _echo_human_summary(_stages: List[Dict[str, Any]], _max_step_size: float) -> None:
            """Print a readable end-of-run summary like the requested example."""
            def _fmt_target_value(x: float) -> str:
                # 2.600 -> "2.6", 1.500 -> "1.5"
                s = f"{x:.3f}".rstrip("0").rstrip(".")
                return s

            def _targets_triplet_str(pairs_1based: List[Tuple[int, int]], targets: List[float]) -> str:
                triples = [f"({i}, {j}, {_fmt_target_value(t)})" for (i, j), t in zip(pairs_1based, targets)]
                return "[" + ", ".join(triples) + "]"

            def _list_of_str_3f(values: List[float]) -> str:
                return "[" + ", ".join(f"'{v:.3f}'" for v in values) + "]"

            click.echo("\nSummary")
            click.echo("------------------")
            for s in _stages:
                idx = int(s.get("index", 0))
                pairs_1b = list(s.get("pairs_1based", []))
                r0 = list(s.get("initial_distances_A", []))
                rT = list(s.get("target_distances_A", []))
                dA = list(s.get("per_pair_step_A", []))
                N = int(s.get("num_steps", 0))
                bchg = s.get("bond_change", {}) or {}
                changed = bool(bchg.get("changed"))
                summary_txt = (bchg.get("summary") or "").strip()

                click.echo(f"[stage {idx}] Targets (i,j,target Å): { _targets_triplet_str(pairs_1b, rT) }")
                click.echo(f"[stage {idx}] initial distances (Å) = { _list_of_str_3f(r0) }")
                click.echo(f"[stage {idx}] target distances  (Å) = { _list_of_str_3f(rT) }")
                click.echo(f"[stage {idx}] per_pair_step     (Å) = { _list_of_str_3f(dA) }")
                click.echo(f"[stage {idx}] steps N = {N}")
                click.echo(f"[stage {idx}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                if changed and summary_txt:
                    click.echo(textwrap.indent(summary_txt, prefix="  "))
                if not changed:
                    click.echo("  (no covalent changes detected)")
                click.echo("")  # blank line between stages

        _echo_human_summary(stages_summary, float(max_step_size))
        # ------------------------------------------------------------------

        click.echo("\n=== Scan finished ===\n")

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for Scan: {hh:02d}:{mm:02d}:{ss:06.3f}")

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
