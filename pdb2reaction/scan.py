# pdb2reaction/scan.py

"""
scan — Bond‑length–driven staged scan with harmonic distance restraints and full relaxation
=====================================================================================================

Usage (CLI)
-----------
    pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} [-q <charge>] [--ligand-charge <number|'RES:Q,...'>] \
        [--scan-list(s) '[(I,J,TARGET_ANG), ...]' ...] [-m <spin>] \
        [--one-based {True|False}] [--max-step-size <float>] \
        [--bias-k <float>] [--relax-max-cycles <int>] \
        [--opt-mode {light|heavy}] [--freeze-links {True|False}] \
        [--dump {True|False}] [--convert-files {True|False}] [--ref-pdb <file>] \
        [--out-dir <dir>] [--thresh <preset>] [--args-yaml <file>] \
        [--preopt {True|False}] [--endopt {True|False}]

Examples
--------
    # Single-stage, minimal inputs (PDB)
    pdb2reaction scan -i input.pdb -q 0 --scan-lists '[(12,45,1.35)]'

    # Two stages, LBFGS, dumping trajectories
    pdb2reaction scan -i input.pdb -q 0 --scan-lists \
        '[(12,45,1.35)]' '[(10,55,2.20),(23,34,1.80)]' \
        --max-step-size 0.2 --dump True --out-dir ./result_scan/ --opt-mode light \
        --preopt True --endopt True


Description
-----------
Runs a staged, bond‑length–driven scan. At each step, harmonic distance wells are applied to
specified atom pairs and the full structure is relaxed. This implementation supports only
the UMA calculator via `uma_pysis` and removes general-purpose handling to reduce overhead.
For PDB inputs, scan tuples can use integer indices or selector strings such as
``'TYR,285,CA'`` and ``'MMT,309,C10'`` to reference atoms (resname, resseq, atom).
For non-PDB inputs, only integer indices are supported.
If you pass one ``--scan-list(s)`` literal, the scan runs in a single stage; multiple
literals are executed as sequential stages, each starting from the previous stage’s
relaxed final structure.
Use ``--ref-pdb`` with XYZ/GJF inputs to load a reference PDB topology while keeping
the XYZ coordinates, enabling format-aware PDB/GJF output conversions.

Scheduling
  - For scan tuples [(i, j, target_Å)], compute the Å‑space displacement Δ = target − current_distance_Å.
  - Let d_max = max(|Δ|). With --max-step-size = h (Å), set N = ceil(d_max / h).
  - Per‑pair step width δ_k = Δ / N (Å).
  - At step s (1..N), the temporary target is r_k(s) = r_k(0) + s · δ_k (Å).
  - Relax the full structure under the harmonic wells for that step.

Harmonic bias model
  - E_bias = Σ ½ · k · (|r_i − r_j| − target_k)².
  - k is provided in eV/Å² (CLI/YAML) and is converted once to Hartree/Bohr² by the bias wrapper.
  - Coordinates are in Bohr; the UMA base returns energies in Hartree and forces in Hartree/Bohr.

Optional optimizations
  - --preopt: preoptimize the initial structure **without bias** before the scan; writes
    `preopt/result.xyz` (and `.pdb` if the input was PDB), then continues from that geometry.
  - --endopt: after **each stage** completes its biased stepping, perform an additional **unbiased**
    optimization of that stage’s final structure before writing outputs.
  - For each stage, covalent‑bond **formation/breaking** is reported between the stage’s first
    structure and its final structure (the latter is the end‑of‑stage optimized result when
    `--endopt True`).
  - By default, both `--preopt` and `--endopt` are enabled.

Charge/spin resolution
  - `-q/--charge` is required unless the input is a `.gjf` template **or** ``--ligand-charge`` is provided.
    When ``-q`` is omitted but ``--ligand-charge`` is set, the full complex is treated as an enzyme–substrate
    system and the total charge is inferred using ``extract.py``’s residue-aware logic. Multiplicity defaults to 1
    when omitted, and an explicit `-q` overrides any derived charge.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_scan/)
  ├─ preopt/                      # Created when --preopt True; holds the unbiased starting optimization
  │   ├─ result.xyz
  │   ├─ result.gjf               # When a GJF template is available **and** conversion is enabled
  │   └─ result.pdb               # When the input was PDB **and** conversion is enabled
  └─ stage_{k:02d}/               # One directory per stage (k = 1..K)
      ├─ result.xyz               # Final structure for the stage (after optional endopt)
      ├─ result.gjf               # When a GJF template is available **and** conversion is enabled
      ├─ result.pdb               # When the input was PDB **and** conversion is enabled
      ├─ scan.trj                 # Written only when --dump True (concatenated biased frames)
      └─ scan.pdb                 # Written only when --dump True, the input was PDB, and conversion is enabled

Notes
-----
- Optimizers: `--opt-mode light` (default) selects LBFGS; `--opt-mode heavy` selects RFOptimizer.
  Step/trust radii are capped in Bohr based on `--max-step-size` (Å).
- Format-aware XYZ/TRJ → PDB/GJF conversions honor the global
  `--convert-files {True|False}` toggle (default: enabled). For stage trajectories (`scan.trj`),
  only PDB companions are produced; GJF companions are written for final `result.xyz` (and preopt outputs).
- Indexing: (i, j) are 1‑based by default; use `--one-based False` if your tuples are 0‑based.
- Provide multiple literals after a single ``--scan-list(s)`` to define sequential stages.
- Units: Distances in CLI/YAML are Å; the bias is applied internally in a.u. (Hartree/Bohr) with
  k converted from eV/Å² to Hartree/Bohr².
- Performance simplifications:
  - Trajectories are accumulated only when `--dump` is True.
  - Energy is not re‑queried for per‑frame annotation during the scan to avoid extra calls.
- PDB convenience: With `--freeze-links` (default True), parent atoms of link hydrogens
  are detected and frozen for PDB inputs, and merged with any user‑specified frozen atoms.
- YAML: Additional arguments can be supplied via `--args-yaml` under sections:
  `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, and `bond`.
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
from pysisyphus.constants import BOHR2ANG, ANG2BOHR

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from .opt import (
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
    HarmonicBiasCalculator,
)
from .utils import (
    convert_xyz_like_outputs,
    detect_freeze_links_safe,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    normalize_choice,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    resolve_atom_spec_index,
)
from .bond_changes import compare_structures, summarize_changes


# --------------------------------------------------------------------------------------
# Defaults (merge order: defaults ← CLI ← YAML)
# --------------------------------------------------------------------------------------

# Geometry handling (Cartesian recommended for scans)
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)

CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

# Optimizer base (convergence, dumping, etc.)
OPT_BASE_KW: Dict[str, Any] = dict(_OPT_BASE_KW)
OPT_BASE_KW.update({
    "max_cycles": 100,          # hard cap on optimization cycles per biased segment
    "out_dir": "./result_scan/",  # output directory for scan artifacts
})

# LBFGS specifics
LBFGS_KW: Dict[str, Any] = dict(_LBFGS_KW)
LBFGS_KW.update({
    "out_dir": "./result_scan/",  # location for LBFGS-specific outputs (restart, etc.)
})

# RFO specifics
RFO_KW: Dict[str, Any] = dict(_RFO_KW)
RFO_KW.update({
    "out_dir": "./result_scan/",  # location for RFO-specific outputs (restart, etc.)
})

# Bias (harmonic well) defaults; can be overridden via YAML: section "bias"
BIAS_KW: Dict[str, Any] = {
    "k": 100,  # float, harmonic bias strength in eV/Å^2
}

# Bond-change detection (as in path_search)
BOND_KW: Dict[str, Any] = {
    "device": "cuda",            # str, device for UMA graph analysis during bond detection
    "bond_factor": 1.20,         # float, scaling of covalent radii for bond cutoff
    "margin_fraction": 0.05,     # float, fractional margin to tolerate small deviations
    "delta_fraction": 0.05,      # float, change threshold to flag bond formation/breaking
}

# Normalization helper
_OPT_MODE_ALIASES = (
    (("light",), "lbfgs"),
    (("heavy",), "rfo"),
)


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


def _parse_scan_lists(
    args: Sequence[str],
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
) -> List[List[Tuple[int, int, float]]]:
    """
    Parse multiple Python-like list strings:
      ['[(0,1,1.5), (2,3,2.0)]', '[(5,7,1.2)]', ...]
    Returns: [[(i,j,t), ...], [(i,j,t), ...], ...] with 0-based indices.
    """
    if not args:
        raise click.BadParameter("--scan-lists must be provided at least once.")
    stages: List[List[Tuple[int, int, float]]] = []
    def _resolve_index(value: Any, stage_idx: int, side_label: str) -> int:
        if isinstance(value, (int, np.integer)):
            idx_val = int(value)
            if one_based:
                idx_val -= 1
            if idx_val < 0:
                raise click.BadParameter(
                    f"Negative atom index in --scan-lists #{stage_idx}: {idx_val} (0-based expected)."
                )
            return idx_val
        if isinstance(value, str):
            if not atom_meta:
                raise click.BadParameter(
                    f"--scan-lists #{stage_idx} ({side_label}) uses a string atom spec, "
                    "but no PDB metadata is available (non-PDB inputs require integer indices)."
                )
            try:
                return resolve_atom_spec_index(value, atom_meta)
            except ValueError as exc:
                raise click.BadParameter(
                    f"--scan-lists #{stage_idx} ({side_label}) {exc}"
                )
        raise click.BadParameter(
            f"--scan-lists #{stage_idx} ({side_label}) must be an int index or atom spec string."
        )

    for idx, s in enumerate(args, start=1):
        try:
            obj = ast.literal_eval(s)
        except Exception as e:
            raise click.BadParameter(f"Invalid literal for --scan-lists #{idx}: {e}")
        if not isinstance(obj, (list, tuple)):
            raise click.BadParameter(f"--scan-lists #{idx} must be a list/tuple of (i,j,target).")
        tuples: List[Tuple[int, int, float]] = []
        for t in obj:
            if not (
                isinstance(t, (list, tuple)) and len(t) == 3
                and isinstance(t[2], (int, float, np.floating))
            ):
                raise click.BadParameter(f"--scan-lists #{idx} contains an invalid triple: {t}")
            i = _resolve_index(t[0], idx, "i")
            j = _resolve_index(t[1], idx, "j")
            r = float(t[2])
            if r <= 0.0:
                raise click.BadParameter(f"Non-positive target length in --scan-lists #{idx}: {(i,j,r)}.")
            tuples.append((i, j, r))
        stages.append(tuples)
    return stages


def _collect_scan_list_values(argv: Sequence[str], names: Sequence[str]) -> Tuple[List[str], int]:
    """Return scan-list literals following a single flag and the number of flag occurrences."""
    values: List[str] = []
    flag_count = 0
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            flag_count += 1
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                values.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    return values, flag_count


def _pair_distances(coords_ang: np.ndarray, pairs: Iterable[Tuple[int, int]]) -> List[float]:
    """
    coords_ang: (N,3) in Å; returns a list of distances (Å) for the given pairs.
    """
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
    """
    Return whether covalent bonds formed or broke between `x` and `y`,
    and a one-based human-readable summary.
    """
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
        snap = geom_loader(
            Path(tmp.name),
            coord_type=getattr(g, "coord_type", GEOM_KW_DEFAULT["coord_type"]),
            freeze_atoms=getattr(g, "freeze_atoms", []),
        )
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


@click.command(
    help="Bond-length driven scan with staged harmonic restraints and relaxation.",
    context_settings={"help_option_names": ["-h", "--help"], "allow_extra_args": True},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, ...).",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=False,
    help=(
        "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
        "(PDB inputs or XYZ/GJF with --ref-pdb)."
    ),
)
@click.option(
    "--workers",
    type=int,
    default=CALC_KW["workers"],
    show_default=True,
    help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel UMA predictor (workers>1).",
)
@click.option(
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
        "when -q is omitted (requires PDB input or --ref-pdb)."
    ),
)
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1) for the ML region.")
@click.option(
    "--scan-lists",
    "--scan-list",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help="Python-like list of (i,j,target) per stage. Pass a single --scan-list(s) followed by "
         "multiple literals to run sequential stages, e.g. --scan-lists '[(0,1,1.50)]' '[(5,7,1.20)]'.",
)
@click.option("--one-based", "one_based", type=click.BOOL, default=True, show_default=True,
              help="Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.")
@click.option("--max-step-size", type=float, default=0.20, show_default=True,
              help="Maximum change in any scanned bond length per step [Å].")
@click.option("--bias-k", type=float, default=100, show_default=True,
              help="Harmonic well strength k [eV/Å^2].")
@click.option("--relax-max-cycles", type=int, default=10000, show_default=True,
              help="Maximum optimizer cycles per relaxation (preopt, per-step, and end-of-stage).")
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Relaxation mode: light (=LBFGS) or heavy (=RFO).",
)
@click.option("--freeze-links", type=click.BOOL, default=True, show_default=True,
              help="If input is PDB, freeze parent atoms of link hydrogens.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Write stage trajectory as scan.trj (and scan.pdb for PDB input).")
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option("--out-dir", type=str, default="./result_scan/", show_default=True,
              help="Base output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help="Convergence preset for relaxations (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file with extra args (sections: geom, calc, opt, lbfgs, rfo, bias, bond).",
)
@click.option("--preopt", type=click.BOOL, default=True, show_default=True,
              help="Preoptimize initial structure without bias before the scan.")
@click.option("--endopt", type=click.BOOL, default=True, show_default=True,
              help="After each stage, run an additional unbiased optimization of the stage result.")
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    scan_lists_raw: Sequence[str],
    one_based: bool,
    max_step_size: float,
    bias_k: Optional[float],
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
    endopt: bool,
) -> None:
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    apply_ref_pdb_override(prepared_input, ref_pdb)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input,
        charge,
        spin,
        ligand_charge=ligand_charge,
        prefix="[scan]",
    )
    needs_pdb = source_path.suffix.lower() == ".pdb"
    needs_gjf = prepared_input.is_gjf
    ref_pdb = source_path.resolve() if needs_pdb else None
    try:
        time_start = time.perf_counter()

        # ------------------------------------------------------------------
        # 1) Assemble configuration (defaults ← CLI ← YAML)
        # ------------------------------------------------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg  = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)
        bias_cfg  = dict(BIAS_KW)
        bond_cfg  = dict(BOND_KW)

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)
        calc_cfg["workers"] = int(workers)
        calc_cfg["workers_per_node"] = int(workers_per_node)
        opt_cfg["out_dir"] = out_dir
        # Do not use the optimizer's own dump per step; stage dumping is controlled separately.
        opt_cfg["dump"]    = False
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)
        kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=_OPT_MODE_ALIASES,
            allowed_hint="light|heavy",
        )

        # YAML overrides (highest precedence)
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

        # Bias strength override
        if bias_k is not None:
            bias_cfg["k"] = float(bias_k)

        # Present final config
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_freeze_atoms_for_echo(calc_cfg)
        echo_opt  = dict(opt_cfg); echo_opt["out_dir"] = str(out_dir_path)
        echo_bias = dict(bias_cfg)
        echo_bond = dict(bond_cfg)
        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("opt",  echo_opt))
        click.echo(pretty_block("lbfgs" if kind == "lbfgs" else "rfo", (lbfgs_cfg if kind == "lbfgs" else rfo_cfg)))
        click.echo(pretty_block("bias", echo_bias))
        click.echo(pretty_block("bond", echo_bond))

        pdb_atom_meta: List[Dict[str, Any]] = []
        if source_path.suffix.lower() == ".pdb":
            pdb_atom_meta = load_pdb_atom_metadata(source_path)

        # ------------------------------------------------------------------
        # 2) Parse scan lists
        # ------------------------------------------------------------------
        argv_scan = sys.argv[1:]
        scan_lists_raw_final, scan_flag_count = _collect_scan_list_values(
            argv_scan, ("--scan-lists", "--scan-list")
        )
        if scan_flag_count > 1:
            raise click.BadParameter(
                "Use a single --scan-list/--scan-lists followed by multiple stage literals; "
                "repeated flags are not accepted."
            )
        if not scan_lists_raw_final:
            raise click.BadParameter("--scan-list(s) must be provided at least once.")
        stages = _parse_scan_lists(
            scan_lists_raw_final, one_based=one_based, atom_meta=pdb_atom_meta
        )
        K = len(stages)
        click.echo(f"[scan] Received {K} stage(s).")

        if pdb_atom_meta:
            click.echo("[scan] PDB atom details for scanned pairs:")
            legend = format_pdb_atom_metadata_header()
            click.echo(f"        legend: {legend}")
            for stage_idx, tuples in enumerate(stages, start=1):
                click.echo(f"  Stage {stage_idx}:")
                for pair_idx, (i, j, _) in enumerate(tuples, start=1):
                    click.echo(
                        f"    pair {pair_idx} i: {format_pdb_atom_metadata(pdb_atom_meta, i)}"
                    )
                    click.echo(
                        f"           j: {format_pdb_atom_metadata(pdb_atom_meta, j)}"
                    )

        # Prepare end-of-run summary collector
        stages_summary: List[Dict[str, Any]] = []

        # ------------------------------------------------------------------
        # 3) Load geometry (Cartesian) and set calculator (UMA → harmonic-bias wrapper)
        # ------------------------------------------------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Load
        freeze = merge_freeze_atom_indices(geom_cfg)
        if freeze_links and source_path.suffix.lower() == ".pdb":
            detected = detect_freeze_links_safe(source_path)
            if detected:
                freeze = merge_freeze_atom_indices(geom_cfg, detected)
                if freeze:
                    click.echo(f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, freeze))}")

        coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
        geom = geom_loader(geom_input_path, coord_type=coord_type, freeze_atoms=freeze)

        max_step_bohr = float(max_step_size) * ANG2BOHR  # shared cap for LBFGS step / RFO trust radii

        def _make_optimizer(kind_local: str, _out_dir: Path, _prefix: str):
            common = dict(opt_cfg)
            common["out_dir"] = str(_out_dir)
            common["prefix"] = _prefix
            common["max_cycles"] = int(relax_max_cycles)
            if kind_local == "lbfgs":
                args = {**lbfgs_cfg, **common}
                args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
                return LBFGS(geom, **args)
            args = {**rfo_cfg, **common}
            tr = float(rfo_cfg.get("trust_radius", 0.30))
            args["trust_radius"] = min(tr, max_step_bohr)
            args["trust_max"] = min(float(rfo_cfg.get("trust_max", 0.30)), max_step_bohr)
            return RFOptimizer(geom, **args)

        # Merge freeze_atoms with link parents (PDB)
        # Attach freeze indices to Geometry for optimizer awareness
        if freeze:
            try:
                geom.freeze_atoms = np.array(freeze, dtype=int)
            except Exception:
                pass

        # Build UMA calculator (only uma_pysis is supported)
        base_calc = uma_pysis(**calc_cfg)

        # ------------------------------------------------------------------
        # Optional preoptimization WITHOUT bias
        # ------------------------------------------------------------------
        if preopt:
            pre_dir = out_dir_path / "preopt"
            pre_dir.mkdir(parents=True, exist_ok=True)
            geom.set_calculator(base_calc)
            click.echo(f"[preopt] Unbiased relaxation ({kind}) ...")
            optimizer0 = _make_optimizer(kind, pre_dir, "preopt_")
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
            try:
                convert_xyz_like_outputs(
                    pre_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=pre_dir / "result.pdb" if needs_pdb else None,
                    out_gjf_path=pre_dir / "result.gjf" if needs_gjf else None,
                )
                if needs_pdb or needs_gjf:
                    written = []
                    if needs_pdb:
                        written.append("'result.pdb'")
                    if needs_gjf:
                        written.append("'result.gjf'")
                    click.echo(f"[convert] Wrote {', '.join(written)}.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert preopt result: {e}", err=True)

        # Wrap with bias calculator for the scan
        biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))
        geom.set_calculator(biased)

        # ------------------------------------------------------------------
        # 4) Stage-by-stage scan
        # ------------------------------------------------------------------

        # Iterate stages
        for k, tuples in enumerate(stages, start=1):
            stage_dir = _ensure_stage_dir(out_dir_path, k)
            click.echo(f"\n--- Stage {k}/{K} ---")
            click.echo(f"Targets (i,j,target Å): {tuples}")

            # Snapshot beginning geometry of this stage for bond-change comparison
            start_geom_for_stage = _snapshot_geometry(geom)

            # Current coordinates (Bohr) and schedule computed in Å
            R_bohr = np.array(geom.coords3d, dtype=float)      # (N,3) Bohr
            R_ang  = R_bohr * BOHR2ANG                         # (N,3) Å
            Nsteps, r0, rT, step_widths = _schedule_for_stage(R_ang, tuples, float(max_step_size))
            click.echo(f"[stage {k}] initial distances (Å) = {['{:.3f}'.format(x) for x in r0]}")
            click.echo(f"[stage {k}] target distances  (Å) = {['{:.3f}'.format(x) for x in rT]}")
            click.echo(f"[stage {k}] steps N = {Nsteps}")

            # Record per-stage summary
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

                # Bond changes: start vs final (possibly endopt)
                try:
                    changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                    click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                    if changed and summary and summary.strip():
                        click.echo(textwrap.indent(summary.strip(), prefix="  "))
                    if not changed:
                        click.echo("  (no covalent changes detected)")
                    srec["bond_change"]["changed"] = bool(changed)
                    srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                except Exception as e:
                    click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                # Write current (possibly endopted) geometry as the stage result
                final_xyz = stage_dir / "result.xyz"
                with open(final_xyz, "w") as f:
                    f.write(_coords3d_to_xyz_string(geom))
                click.echo(f"[write] Wrote '{final_xyz}'.")
                try:
                    convert_xyz_like_outputs(
                        final_xyz,
                        prepared_input,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=stage_dir / "result.pdb" if needs_pdb else None,
                        out_gjf_path=stage_dir / "result.gjf" if needs_gjf else None,
                    )
                    if needs_pdb or needs_gjf:
                        written = []
                        if needs_pdb:
                            written.append("'result.pdb'")
                        if needs_gjf:
                            written.append("'result.gjf'")
                        click.echo(f"[convert] Wrote {', '.join(written)}.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert stage result: {e}", err=True)
                continue

            # Run N step(s) with bias
            for s in range(1, Nsteps + 1):
                # Compute per-pair step target (Å) for this step
                step_targets = [r0_i + s * dw for (r0_i, dw) in zip(r0, step_widths)]

                # Update bias well targets (still in Å; wrapper converts internally)
                biased.set_pairs([(i, j, t) for ((i, j), t) in zip(pairs, step_targets)])
                # Flushing Geometry caches by re-attaching the calculator
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

            # Bond changes: start vs final (possibly endopt)
            try:
                changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                if changed and summary and summary.strip():
                    click.echo(textwrap.indent(summary.strip(), prefix="  "))
                if not changed:
                    click.echo("  (no covalent changes detected)")
                srec["bond_change"]["changed"] = bool(changed)
                srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
            except Exception as e:
                click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

            # Stage outputs
            if dump and trj_blocks:
                trj_path = stage_dir / "scan.trj"
                with open(trj_path, "w") as f:
                    f.write("".join(trj_blocks))
                click.echo(f"[write] Wrote '{trj_path}'.")
                try:
                    convert_xyz_like_outputs(
                        trj_path,
                        prepared_input,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=stage_dir / "scan.pdb" if needs_pdb else None,
                        out_gjf_path=stage_dir / "scan.gjf" if needs_gjf else None,
                    )
                    if needs_pdb or needs_gjf:
                        written = []
                        if needs_pdb:
                            written.append("'scan.pdb'")
                        if needs_gjf:
                            written.append("'scan.gjf'")
                        click.echo(f"[convert] Wrote {', '.join(written)}.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert stage trajectory: {e}", err=True)

            final_xyz = stage_dir / "result.xyz"
            with open(final_xyz, "w") as f:
                f.write(_coords3d_to_xyz_string(geom))
            click.echo(f"[write] Wrote '{final_xyz}'.")
            try:
                convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=stage_dir / "result.pdb" if needs_pdb else None,
                    out_gjf_path=stage_dir / "result.gjf" if needs_gjf else None,
                )
                if needs_pdb or needs_gjf:
                    written = []
                    if needs_pdb:
                        written.append("'result.pdb'")
                    if needs_gjf:
                        written.append("'result.gjf'")
                    click.echo(f"[convert] Wrote {', '.join(written)}.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert stage result: {e}", err=True)

        # ------------------------------------------------------------------
        # 5) Final summary echo (human‑friendly)
        # ------------------------------------------------------------------
        def _echo_summary(_stages: List[Dict[str, Any]]) -> None:
            """
            Print a readable end-of-run summary.
            """
            def _fmt_target_value(x: float) -> str:
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

        _echo_summary(stages_summary)
        # ------------------------------------------------------------------

        click.echo("\n=== Scan finished ===\n")

        click.echo(format_elapsed("[time] Elapsed Time for Scan", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()
