# pdb2reaction/all.py

"""
all — SINGLE command to execute an end-to-end enzymatic reaction workflow:
(Optional PDB element-field fix) → (optional) extract active-site pockets → (optional) staged scan on a single
structure → MEP search (recursive ``path_search`` by default; pairwise ``path-opt`` with ``--refine-path False``)
→ (optional) merge back into full-system PDB templates → optional TS optimization, IRC (EulerPC),
thermochemistry (freq), DFT single-points, and DFT//UMA diagrams.
================================================================================================================

Usage (CLI)
-----------
    pdb2reaction all -i INPUT1 [INPUT2 ...] [-c <substrate-spec>] \
        [--ligand-charge <number|"RES:Q,...">] [-q/--charge <forced_net_charge>] [-m/--mult <2S+1>] \
        [--freeze-links {True|False}] [--mep-mode {gsm|dmf}] [--max-nodes <int>] [--max-cycles <int>] \
        [--climb {True|False}] [--opt-mode {light|heavy}] [--dump {True|False}] \
        [--convert-files/--no-convert-files] [--refine-path {True|False}] [--thresh <preset>] \
        [--args-yaml <file>] [--preopt {True|False}] \
        [--hessian-calc-mode {Analytical|FiniteDifference}] [--out-dir <dir>] \
        [--tsopt {True|False}] [--thermo {True|False}] [--dft {True|False}] \
        [--tsopt-max-cycles <int>] [--freq-* overrides] [--dft-* overrides] \
        [--dft-engine {gpu|cpu|auto}] [--scan-lists "[(...)]" ...]

    # Single-structure scan builder (repeat --scan-lists to define sequential stages; works with or without extraction)
    pdb2reaction all -i INPUT.pdb [-c <substrate-spec>] --scan-lists "[(...)]" [...]

    # Single-structure TSOPT-only mode (no MEP search) when --scan-lists is omitted and --tsopt True
    pdb2reaction all -i INPUT.pdb [-c <substrate-spec>] --tsopt True [other options]


Examples
--------
    # Minimal end-to-end run with explicit substrate and ligand charges (multi-structure)
    pdb2reaction all -i reactant.pdb product.pdb -c "GPP,MMT" \
        --ligand-charge "GPP:-3,MMT:-1"

    # Full ensemble with an intermediate, residue-ID substrate spec, and full post-processing
    pdb2reaction all -i A.pdb B.pdb C.pdb -c "308,309" --ligand-charge "-1" \
        --mult 1 --freeze-links True --max-nodes 10 --max-cycles 100 --climb True \
        --opt-mode light --dump False --args-yaml params.yaml --preopt True \
        --out-dir result_all --tsopt True --thermo True --dft True

    # Single-structure + staged scan to build an ordered series before the MEP step + post-processing
    pdb2reaction all -i A.pdb -c "308,309" \
        --scan-lists "[(10,55,2.20),(23,34,1.80)]" --mult 1 --out-dir result_scan_all \
        --tsopt True --thermo True --dft True

    # Single-structure TSOPT-only mode (skips the MEP step entirely)
    pdb2reaction all -i A.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
        --tsopt True --thermo True --dft True --out-dir result_tsopt_only


Description
-----------
Runs a one-shot pipeline centered on pocket models. The command is intentionally permissive about how
``-i/--input`` is given: it accepts repeated ``-i`` flags and the sloppy form ``-i A.pdb B.pdb C.pdb``.

Pipeline overview
-----------------
(0) **Preflight: PDB element-field repair (automatic, best-effort)**
    - If any input PDB contains ATOM/HETATM records with an empty element column (PDB cols 77–78),
      ``add_elem_info`` is invoked to fill elements before downstream steps.
    - Fixed copies are written to ``<out-dir>/add_elem_info/`` and used for extraction/MEP.

(1) **Active-site pocket extraction** *(enabled when ``-c/--center`` is provided)*
    - Define the substrate/center by:
        • a PDB path,
        • residue IDs like ``"123,124"`` or ``"A:123,B:456"`` (insertion codes allowed: ``"123A"`` / ``"A:123A"``),
        • residue names like ``"GPP,MMT"``.
    - The extractor writes per-input pocket PDBs under ``<out-dir>/pockets/`` as ``pocket_<stem>.pdb``.
    - The extractor’s **model #1 total pocket charge** is used as the workflow charge for later stages,
      cast to the nearest integer; if rounding occurs a NOTE is printed.
    - Common extractor knobs: ``--radius``, ``--radius-het2het``, ``--include-H2O``,
      ``--exclude-backbone``, ``--add-linkH``, ``--selected_resn``, ``--verbose``.
    - ``--ligand-charge`` can be either:
        • a numeric total charge (distributed across unknown residues), or
        • a residue-name mapping like ``"GPP:-3,MMT:-1"``.

(1′) **Extraction skipped** *(when ``-c/--center`` is omitted)*
    - No pocket is built; the **full input structure(s)** are used directly as the “pocket inputs”.
    - Charge/spin are resolved as described in **Charge handling** below.
    - Full-system merge is skipped (there is no separate “template” step).

(1b) **Optional staged scan (single-structure only; requires ``--scan-lists``)**
    - Triggered only when **exactly one** input is given **and** ``--scan-lists`` is provided.
    - The scan is a staged, bond-length–driven scan using the UMA calculator:
        • repeat ``--scan-lists`` to define sequential stages,
        • each stage starts from the previous stage’s relaxed final structure (stages are chained).
    - When extraction is enabled, indices in ``--scan-lists`` refer to the **original full input PDB**
      (**ATOM/HETATM file order**, 1-based) and are **auto-mapped** onto the extracted pocket using structural atom identity
      (chain/residue/atom name, etc.) with occurrence counting for duplicates.
    - Stage endpoints are collected as ``scan/stage_XX/result.(pdb|xyz|gjf)`` and appended to the ordered
      input series for the MEP step:
      ``[initial pocket (or full input), stage_01/result.*, stage_02/result.*, ...]``.

(2) **MEP search on the ordered pocket inputs**
    Two modes are supported:

    **(2a) Recursive MEP search (default):** ``--refine-path True``
      - Runs the recursive ``path_search`` workflow on the ordered series.
      - ``--mep-mode`` selects the optimizer: Growing String Method (``gsm``) or Direct Max Flux (``dmf``).
      - When the original inputs are PDBs and extraction is enabled, the original **full-system** PDBs are
        passed as merge templates (``--ref-pdb``) automatically:
          • multi-input: one template per pocket input in reaction order,
          • single+scan: the same original template is reused for all scan-derived structures.

    **(2b) Pairwise GSM concatenation:** ``--refine-path False``
      - Runs ``path-opt`` once for each adjacent pair in the ordered series and concatenates the resulting
        trajectories into ``<out-dir>/path_opt/mep.trj`` (no recursive search).
      - Full-system merge is not performed in this mode.

    Shared knobs for the MEP step:
      - ``--mult`` (multiplicity), ``--freeze-links``, ``--max-nodes``, ``--max-cycles``, ``--climb``,
        ``--opt-mode``, ``--dump``, ``--thresh``, ``--preopt``, ``--args-yaml``, ``--out-dir``.
      - Note: ``--dump`` is always passed to the MEP step; for scan/tsopt/freq it is forwarded only when
        explicitly set on this command (otherwise each subcommand’s defaults apply; freq is run with dump=True by default).
      - ``--convert-files/--no-convert-files`` controls whether XYZ/TRJ outputs are converted into
        PDB/GJF companions when possible.

(3) **Merge to full systems (PDB templates only)**
    - Only applicable in the recursive ``path_search`` branch (``--refine-path True``) and only when
      full-system PDB templates are available (see (2a)).
    - In that case, merged trajectories (e.g., ``mep_w_ref.pdb`` and per-segment merged files such as ``*_seg_XX.pdb``)
      are produced by ``path_search`` under ``<out-dir>/path_search/``.
      This wrapper mirrors only the main merged trajectory and key summary/plot files to ``<out-dir>/``;
      per-segment merged files remain under ``<out-dir>/path_search/``.
    - When templates are unavailable (non-PDB inputs or extraction skipped), only pocket-level outputs exist.

(4) **Optional per-segment post-processing (reactive segments only)**
    Post-processing is applied only to segments that are flagged as having **covalent/bond changes** in
    ``summary.yaml``. For each such segment:

    - ``--tsopt True``:
        • Optimize a TS starting from the segment HEI snapshot,
        • Run **IRC (EulerPC)** from the optimized TS,
        • Map IRC endpoints onto the segment’s left/right MEP endpoints (bond-state matching first,
          RMSD fallback) to assign backward/forward consistently across segments,
        • Optimize both IRC endpoints to minima (LBFGS/RFO depending on ``--opt-mode``),
        • Use these optimized minima as the final R and P for downstream analyses.
        • Write per-segment UMA energy diagram: ``post_seg_XX/energy_diagram_UMA.png``.
        • Copy the TS structure to ``<out-dir>/ts_seg_XX.pdb`` (when a PDB representation is available)
          or ``<out-dir>/ts_seg_XX.xyz`` otherwise.

    - ``--thermo True`` (with or without ``--tsopt``):
        • Run ``freq`` on (R, TS, P) and build a UMA Gibbs diagram
          ``post_seg_XX/energy_diagram_G_UMA.png``.
        • If ``--tsopt`` is **off**, “TS” is the raw HEI snapshot from the MEP; R/P come from the segment
          endpoints (no IRC-based reassignment).

    - ``--dft True`` (with or without ``--tsopt``):
        • Run DFT single-point on (R, TS, P) and build ``post_seg_XX/energy_diagram_DFT.png``.
        • With ``--thermo True``, also build **DFT//UMA** Gibbs diagram
          ``post_seg_XX/energy_diagram_G_DFT_plus_UMA.png``.

    Across all reactive segments, aggregated plots are written at ``<out-dir>/`` when the required data exist:
      - ``energy_diagram_UMA_all.png`` (UMA R–TS–P; requires ``--tsopt``),
      - ``energy_diagram_G_UMA_all.png`` (UMA Gibbs; requires ``--thermo``),
      - ``energy_diagram_DFT_all.png`` (DFT; requires ``--dft``),
      - ``energy_diagram_G_DFT_plus_UMA_all.png`` (DFT//UMA Gibbs; requires ``--dft`` and ``--thermo``),
      - ``irc_plot_all.png`` (aggregated IRC plot; requires ``--tsopt``).

(Alt) **Single-structure TSOPT-only mode**
    - Triggered when **exactly one** input is given, **no** ``--scan-lists`` is provided, and ``--tsopt True``.
    - The tool skips the MEP step entirely and instead:
        • Runs ``tsopt`` on the pocket (or full input when extraction is skipped),
        • Runs IRC (EulerPC) and obtains both endpoints,
        • Optimizes both endpoints to minima and builds UMA diagrams for **R–TS–P**,
        • Optionally adds UMA Gibbs, DFT, and **DFT//UMA** diagrams.
    - In this specific mode only, the higher-energy IRC endpoint (as returned by IRC, before endpoint minimization)
      is treated as the reactant (R).
    - Outputs live under ``<out-dir>/tsopt_single/`` and key plots are mirrored to ``<out-dir>/`` as
      ``*_all.png``; the IRC plot is mirrored as ``irc_plot_all.png``; the TS is copied as
      ``ts_seg_01.(pdb|xyz)``.

Charge handling
---------------
  - ``-q/--charge`` **forces** the total system charge with a console **WARNING**, overriding
    extractor/GJF/``--ligand-charge``-derived values.
  - With extraction enabled (``-c/--center``):
      • the extractor’s model #1 **total pocket charge** is used (rounded to an integer) unless overridden.
  - With extraction skipped (no ``--center``):
      1) a numeric ``--ligand-charge`` value is interpreted as the **TOTAL system charge** (rounded),
      2) else, if the first input is ``.gjf``, its charge is used,
      3) else default is **0**.
    Non-numeric ``--ligand-charge`` mappings are ignored in this mode.
  - Spin precedence when extraction is skipped:
      explicit ``--mult`` > GJF multiplicity (if available) > default.

Inputs
------
  - ``-i/--input`` accepts:
      • Two or more structures in reaction order (reactant [intermediates ...] product), or
      • A single structure when using ``--scan-lists`` (staged scan), or
      • A single structure when using ``--tsopt True`` (TSOPT-only mode).
  - When using extraction (``-c/--center``), inputs must be **PDB**.
  - When using ``--scan-lists`` with extraction skipped, the single input must be a **PDB**.
  - When extraction is skipped, multi-structure runs accept **PDB/XYZ/GJF** inputs.

Forwarded / relevant options
----------------------------
  - MEP search: ``--mult``, ``--freeze-links``, ``--mep-mode``, ``--max-nodes``, ``--max-cycles``,
    ``--climb``, ``--opt-mode``, ``--dump``, ``--thresh``, ``--preopt``, ``--args-yaml``, ``--out-dir``.
  - File conversion: ``--convert-files/--no-convert-files`` toggles conversion of XYZ/TRJ outputs into
    PDB/GJF companions when possible.
  - Scan (single-structure): inherits charge/spin and shared knobs; per-stage/scan overrides include
    ``--scan-out-dir``, ``--scan-one-based``, ``--scan-max-step-size``, ``--scan-bias-k``,
    ``--scan-relax-max-cycles``, ``--scan-preopt``, ``--scan-endopt``.
  - TS optimization / IRC: ``--tsopt``, ``--tsopt-max-cycles``, ``--tsopt-out-dir`` and shared knobs above.
    ``--hessian-calc-mode`` is forwarded to TSOPT and freq.
  - Frequency analysis: ``--freq-out-dir``, ``--freq-max-write``, ``--freq-amplitude-ang``,
    ``--freq-n-frames``, ``--freq-sort``, ``--freq-temperature``, ``--freq-pressure`` plus shared knobs.
  - DFT single-points: ``--dft-out-dir``, ``--dft-func-basis``, ``--dft-max-cycle``, ``--dft-conv-tol``,
    ``--dft-grid-level``, ``--dft-engine``.
  - YAML forwarding: ``--args-yaml`` is passed unchanged to ``path_search``, ``path-opt``, ``scan``,
    ``tsopt``, ``freq``, and ``dft`` so a single file can host per-module sections
    (see each subcommand for accepted keys).

Outputs (& Directory Layout)
----------------------------
<out-dir>/ (default: ./result_all/)
  ├─ add_elem_info/                      # created only when element fields are missing in inputs
  │   ├─ <input1>.pdb
  │   └─ ...
  ├─ pockets/                            # created when extraction (-c/--center) is used
  │   ├─ pocket_<input1_stem>.pdb
  │   ├─ pocket_<input2_stem>.pdb
  │   └─ ...
  ├─ scan/                               # present only in single-structure + scan mode
  │   ├─ stage_01/result.(pdb|xyz|gjf)
  │   ├─ stage_02/result.(pdb|xyz|gjf)
  │   └─ ...
  ├─ path_search/                        # when --refine-path True (recursive path_search)
  │   ├─ mep.trj                         # final pocket MEP as XYZ/TRJ trajectory
  │   ├─ summary.yaml
  │   ├─ summary.log                     # human-readable summary (also mirrored to <out-dir>/)
  │   ├─ mep_seg_XX.(trj|pdb)            # pocket-only segment trajectories
  │   ├─ hei_seg_XX.(xyz|pdb|gjf)        # HEI snapshots per reactive segment
  │   ├─ hei_w_ref_seg_XX.pdb            # merged HEI per segment (when --ref-pdb / PDB input)
  │   ├─ mep_w_ref_seg_XX.pdb            # merged per-segment trajectories (when --ref-pdb / PDB input; not mirrored)
  │   ├─ seg_XXX_~~~/ ...                # GSM internals / recursion tree
  │   └─ post_seg_XX/                    # created when downstream post-processing runs
  │       ├─ ts/ ...
  │       ├─ irc/ ...
  │       ├─ structures/
  │       ├─ freq/ ...                   # with --thermo True
  │       ├─ dft/  ...                   # with --dft True
  │       ├─ energy_diagram_UMA.png
  │       ├─ energy_diagram_G_UMA.png
  │       ├─ energy_diagram_DFT.png
  │       └─ energy_diagram_G_DFT_plus_UMA.png
  ├─ path_opt/                           # when --refine-path False (pairwise path-opt + concatenation)
  │   ├─ mep.trj, summary.yaml, summary.log
  │   ├─ mep_plot.png
  │   ├─ energy_diagram_mep.png          # compressed diagram for concatenated MEP (when available; not mirrored)
  │   ├─ mep_seg_XX.* , hei_seg_XX.*     # per-segment outputs copied from each pairwise run
  │   ├─ seg_XXX_mep/ ...                # per-pair path-opt run directories
  │   └─ post_seg_XX/ ...                # same layout as path_search/ when post-processing runs
  ├─ mep_plot.png                        # copied from path_* / (MEP energy profile)
  ├─ mep.trj / mep.xyz                   # copied from path_*/ when available (pocket MEP trajectory)
  ├─ mep_w_ref.trj / mep_w_ref.xyz       # copied when produced (merged XYZ/TRJ companions)
  ├─ energy_diagram_MEP.png              # produced by path_search and copied to <out-dir>/ when available
  │                                      # (in --refine-path False branch, the compressed diagram is
  │                                      #  written as path_opt/energy_diagram_mep.png and is not copied by default)
  ├─ mep.pdb                             # copied from path_* when a PDB representation is available
  ├─ mep_w_ref.pdb                       # copied from path_search/ when full-system merge is available
  ├─ summary.yaml                        # mirrored from path_*/summary.yaml
  ├─ summary.log                         # mirrored from path_*/summary.log
  ├─ energy_diagram_UMA_all.png
  ├─ energy_diagram_G_UMA_all.png
  ├─ energy_diagram_DFT_all.png
  ├─ energy_diagram_G_DFT_plus_UMA_all.png
  ├─ irc_plot_all.png
  ├─ ts_seg_XX.pdb / ts_seg_XX.xyz       # TS structure per reactive segment (format depends on inputs)
  │
  └─ tsopt_single/                       # present only in single-structure TSOPT-only mode
      ├─ ts/ ...
      ├─ irc/ ...
      ├─ structures/
      ├─ freq/ ...                       # with --thermo True
      ├─ dft/  ...                       # with --dft True
      ├─ energy_diagram_UMA.png
      ├─ energy_diagram_G_UMA.png
      ├─ energy_diagram_DFT.png
      ├─ energy_diagram_G_DFT_plus_UMA.png
      ├─ summary.yaml
      ├─ summary.log
      └─ (mirrored diagrams as *_all.png at <out-dir>/)

Notes
-----
- ``-c/--center`` is strongly recommended for meaningful pocket models; without it the full structure is used.
- A **single-structure** run requires either ``--scan-lists`` (staged scan) or ``--tsopt True`` (TSOPT-only).
- Energies in diagrams are plotted relative to the first state in kcal/mol (converted from Hartree).
- ``--no-convert-files`` may suppress generation of PDB/GJF companion files from XYZ/TRJ outputs; the core
  numeric results and YAML summaries are unaffected.
"""

from __future__ import annotations

import ast
from pathlib import Path
from collections import defaultdict
from typing import List, Sequence, Optional, Tuple, Dict, Any

import sys, os
import math
import click
from click.core import ParameterSource
import time
import yaml
import numpy as np
import shutil

# Biopython for PDB parsing (post-processing helpers)
from Bio import PDB

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength

AtomKey = Tuple[str, str, str, str, str, str]

# Local imports from the package
from .extract import extract_api
from . import path_search as _path_search
from . import path_opt as _path_opt
from . import tsopt as _tsopt
from . import freq as _freq_cli
from . import dft as _dft_cli
from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
DEFAULT_COORD_TYPE = GEOM_KW_DEFAULT["coord_type"]
from .trj2fig import run_trj2fig
from .summary_log import write_summary_log
from .utils import (
    build_energy_diagram,
    convert_xyz_like_outputs,
    detect_freeze_links_safe,
    format_elapsed,
    prepare_input_structure,
    maybe_convert_xyz_to_gjf,
    set_convert_file_enabled,
    resolve_charge_spin_or_raise,
    load_yaml_dict,
    apply_yaml_overrides,
)
from . import scan as _scan_cli
from .add_elem_info import assign_elements as _assign_elem_info
from . import irc as _irc_cli


# -----------------------------
# Helpers
# -----------------------------


def _close_matplotlib_figures() -> None:
    """Best-effort cleanup for matplotlib figures to avoid open-figure warnings."""

    try:
        import matplotlib.pyplot as plt  # type: ignore

        plt.close("all")
    except Exception:
        pass


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
    """
    Robustly collect values following a flag that may appear **once** followed by multiple space-separated values,
    e.g., "-i A B C". This mirrors the bauavior implemented in `path_search.cli`.
    """
    vals: List[str] = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    return vals


def _append_cli_arg(args: List[str], flag: str, value: Any | None) -> None:
    """Append ``flag`` and ``value`` (converted to string) to ``args`` when ``value`` is not ``None``."""
    if value is None:
        return
    if isinstance(value, bool):
        args.extend([flag, "True" if value else "False"])
    else:
        args.extend([flag, str(value)])


def _resolve_override_dir(default: Path, override: Path | None) -> Path:
    """Return ``override`` when provided (respecting absolute paths); otherwise ``default``."""
    if override is None:
        return default
    if override.is_absolute():
        return override
    return default.parent / override


CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)


def _build_calc_cfg(charge: int, spin: int, yaml_cfg: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Return a UMA calculator configuration honoring YAML overrides when provided."""
    cfg: Dict[str, Any] = dict(CALC_KW)
    cfg["charge"] = int(charge)
    cfg["spin"] = int(spin)
    if yaml_cfg:
        apply_yaml_overrides(
            yaml_cfg,
            [
                (cfg, (("calc",),)),
            ],
        )
    return cfg


def _parse_atom_key_from_line(line: str) -> Optional[AtomKey]:
    """
    Extract a structural identity key from a PDB ATOM/HETATM record.

    Returns:
        (chainID, resName, resSeq, iCode, atomName, altLoc), with blanks normalized to ''.
    """
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    atomname = line[12:16].strip()
    altloc = (line[16] if len(line) > 16 else " ").strip()
    resname = line[17:20].strip()
    chain = (line[21] if len(line) > 21 else " ").strip()
    resseq = line[22:26].strip()
    icode = (line[26] if len(line) > 26 else " ").strip()
    return (chain, resname, resseq, icode, atomname, altloc)


def _key_variants(key: AtomKey) -> List[AtomKey]:
    """Return key variants with progressively relaxed identity fields (deduplicated)."""
    chain, resn, resseq, icode, atom, alt = key
    raw_variants = [
        (chain, resn, resseq, icode, atom, alt),
        (chain, resn, resseq, icode, atom, ""),
        (chain, resn, resseq, "", atom, alt),
        (chain, resn, resseq, "", atom, ""),
    ]
    seen: set[AtomKey] = set()
    variants: List[AtomKey] = []
    for variant in raw_variants:
        if variant in seen:
            continue
        seen.add(variant)
        variants.append(variant)
    return variants


def _build_variant_occurrence_table(keys: Sequence[AtomKey]) -> List[Dict[AtomKey, int]]:
    """
    Track how many times each relaxed key variant has appeared up to each atom index.
    Returns a per-atom list of dicts: variant -> occurrence count (1-based).
    """
    counts: Dict[AtomKey, int] = defaultdict(int)
    per_atom: List[Dict[AtomKey, int]] = []
    for key in keys:
        current: Dict[AtomKey, int] = {}
        for variant in _key_variants(key):
            counts[variant] += 1
            current[variant] = counts[variant]
        per_atom.append(current)
    return per_atom


def _pocket_key_to_index(pocket_pdb: Path) -> Dict[AtomKey, List[int]]:
    """
    Build mapping: structural atom key -> list of pocket indices (1-based by file order).
    """
    key2idx: Dict[AtomKey, List[int]] = defaultdict(list)
    idx = 0
    try:
        with open(pocket_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    key = _parse_atom_key_from_line(line)
                    if key is None:
                        continue
                    idx += 1
                    for variant in _key_variants(key):
                        key2idx[variant].append(idx)
    except FileNotFoundError:
        raise click.ClickException(f"[all] Pocket PDB not found: {pocket_pdb}")
    if not key2idx:
        raise click.ClickException(f"[all] Pocket PDB {pocket_pdb} has no ATOM/HETATM records.")
    return dict(key2idx)


def _read_full_atom_keys_in_file_order(full_pdb: Path) -> List[AtomKey]:
    """
    Read ATOM/HETATM lines and return keys in the original file order.
    """
    keys: List[AtomKey] = []
    try:
        with open(full_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    key = _parse_atom_key_from_line(line)
                    if key is not None:
                        keys.append(key)
    except FileNotFoundError:
        raise click.ClickException(f"[all] File not found while parsing PDB: {full_pdb}")
    if not keys:
        raise click.ClickException(f"[all] No ATOM/HETATM records detected in {full_pdb}.")
    return keys


def _format_atom_key_for_msg(key: AtomKey) -> str:
    """Pretty string for diagnostics."""
    chain, resn, resseq, icode, atom, alt = key
    res = f"{chain}:{resn}{resseq}{(icode if icode else '')}"
    alt_sfx = f",alt={alt}" if alt else ""
    return f"{res}:{atom}{alt_sfx}"


def _parse_scan_lists_literals(scan_lists_raw: Sequence[str]) -> List[List[Tuple[int, int, float]]]:
    """Parse ``--scan-lists`` literals without re-basing atom indices."""
    stages: List[List[Tuple[int, int, float]]] = []
    for idx_stage, literal in enumerate(scan_lists_raw, start=1):
        try:
            obj = ast.literal_eval(literal)
        except Exception as exc:
            raise click.BadParameter(f"Invalid literal for --scan-lists #{idx_stage}: {exc}")
        if not isinstance(obj, (list, tuple)):
            raise click.BadParameter(
                f"--scan-lists #{idx_stage} must be a list/tuple of (i,j,target)."
            )
        tuples: List[Tuple[int, int, float]] = []
        for t in obj:
            if (
                isinstance(t, (list, tuple))
                and len(t) == 3
                and isinstance(t[0], (int, np.integer))
                and isinstance(t[1], (int, np.integer))
                and isinstance(t[2], (int, float, np.floating))
            ):
                tuples.append((int(t[0]), int(t[1]), float(t[2])))
            else:
                raise click.BadParameter(
                    f"--scan-lists #{idx_stage} contains an invalid triple: {t}"
                )
        if not tuples:
            raise click.BadParameter(
                f"--scan-lists #{idx_stage} must contain at least one (i,j,target) triple."
            )
        stages.append(tuples)
    return stages


def _format_scan_stage(stage: List[Tuple[int, int, float]]) -> str:
    """Serialize a scan stage back into a Python-like literal string."""
    return "[" + ", ".join(f"({i},{j},{target})" for (i, j, target) in stage) + "]"


def _convert_scan_lists_to_pocket_indices(
    scan_lists_raw: Sequence[str],
    full_input_pdb: Path,
    pocket_pdb: Path,
) -> List[List[Tuple[int, int, float]]]:
    """
    Convert user-provided atom indices (based on the full input PDB) to pocket indices.
    Returns the converted stages as lists of (i,j,target) with 1-based pocket indices.

    Structural keys (chainID, resName, resSeq, iCode, atomName, altLoc) are used instead of serial numbers,
    with per-variant occurrence counts to distinguish atoms that otherwise share the same key.
    """
    if not scan_lists_raw:
        return []

    stages = _parse_scan_lists_literals(scan_lists_raw)

    orig_keys_in_order = _read_full_atom_keys_in_file_order(full_input_pdb)
    key_to_pocket_idx = _pocket_key_to_index(pocket_pdb)
    variant_occ_table = _build_variant_occurrence_table(orig_keys_in_order)

    n_atoms_full = len(orig_keys_in_order)

    def _map_full_index_to_pocket(idx_one_based: int, stage_idx: int, tuple_idx: int, side_label: str) -> int:
        """
        Convert a 1-based index from the full PDB into the pocket's 1-based index.
        Fall back in the order: strict match → ignore altloc → ignore iCode → ignore both,
        and use the atom index (occurrence count) when multiple atoms share a structural key.
        """
        key = orig_keys_in_order[idx_one_based - 1]

        variant_occ = variant_occ_table[idx_one_based - 1]
        for variant in _key_variants(key):
            occurrence = variant_occ.get(variant)
            indices = key_to_pocket_idx.get(variant)
            if occurrence is None or not indices:
                continue
            if occurrence <= len(indices):
                return indices[occurrence - 1]

        msg_key = _format_atom_key_for_msg(key)
        raise click.BadParameter(
            f"--scan-lists #{stage_idx} tuple #{tuple_idx} ({side_label}) references atom index {idx_one_based} "
            f"(key {msg_key}) which is not present in the pocket after extraction. "
            "Increase extraction coverage (e.g., --radius/--radius-het2het, --selected_resn, or set --exclude-backbone False), "
            "or choose atoms that survive in the pocket."
        )

    converted: List[List[Tuple[int, int, float]]] = []
    for stage_idx, stage in enumerate(stages, start=1):
        stage_converted: List[Tuple[int, int, float]] = []
        for tuple_idx, (idx_i, idx_j, target) in enumerate(stage, start=1):
            if idx_i <= 0 or idx_j <= 0:
                raise click.BadParameter(
                    f"--scan-lists #{stage_idx} tuple #{tuple_idx} must use 1-based atom indices."
                )
            if idx_i > n_atoms_full or idx_j > n_atoms_full:
                raise click.BadParameter(
                    f"--scan-lists #{stage_idx} tuple #{tuple_idx} references an atom index "
                    f"beyond the input PDB atom count ({n_atoms_full})."
                )

            pi = _map_full_index_to_pocket(idx_i, stage_idx, tuple_idx, "i")
            pj = _map_full_index_to_pocket(idx_j, stage_idx, tuple_idx, "j")

            stage_converted.append((pi, pj, target))
        converted.append(stage_converted)
    return converted


def _round_charge_with_note(q: float) -> int:
    """
    Cast the extractor's total charge (float) to an integer suitable for the path search.
    If it is not already an integer within 1e-6, round to the nearest integer with a console note.
    """
    q_rounded = int(round(float(q)))
    if not math.isfinite(q):
        raise click.BadParameter(f"Computed total charge is non-finite: {q!r}")
    if abs(float(q) - q_rounded) > 1e-6:
        click.echo(f"[all] NOTE: extractor total charge = {q:g} → rounded to integer {q_rounded} for the path search.")
    return q_rounded


def _pdb_needs_elem_fix(p: Path) -> bool:
    """
    Return True if the file contains ATOM/HETATM records and at least one has an empty element field (cols 77–78).
    This is a light-weight check to decide whether to run add_elem_info.
    """
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if len(line) < 78 or not line[76:78].strip():
                        return True
        return False
    except Exception:
        return False


_FREEZE_ATOMS_GLOBAL: Optional[List[int]] = None


def _get_freeze_atoms(pdb_path: Optional[Path], freeze_links_flag: bool) -> List[int]:
    """
    Determine freeze atom indices once and reuse them globally.

    The first time this is called with a PDB path and freeze_links_flag=True,
    link-parent atoms are detected from that PDB. The resulting indices are
    reused for subsequent calls (even if a non-PDB path is provided), under
    the assumption that atom indexing is consistent across the trajectory.
    """
    global _FREEZE_ATOMS_GLOBAL
    if not freeze_links_flag:
        return []
    if _FREEZE_ATOMS_GLOBAL is not None:
        return _FREEZE_ATOMS_GLOBAL
    if pdb_path is None or pdb_path.suffix.lower() != ".pdb":
        # No suitable PDB available yet to determine freeze atoms.
        return []
    try:
        fa = detect_freeze_links_safe(pdb_path)
        _FREEZE_ATOMS_GLOBAL = [int(i) for i in fa]
    except Exception as e:
        click.echo(
            f"[all] WARNING: detect_freeze_links_safe failed for {pdb_path}: {e}; no atoms will be frozen.",
            err=True,
        )
        _FREEZE_ATOMS_GLOBAL = []
    return _FREEZE_ATOMS_GLOBAL


def _freeze_atoms_for_log() -> List[int]:
    """Return a sorted freeze_atoms list for summary logs (may be empty)."""

    try:
        return sorted({int(i) for i in (_FREEZE_ATOMS_GLOBAL or [])})
    except Exception:
        return []


# ---------- Post-processing helpers ----------


def _read_summary(summary_yaml: Path) -> List[Dict[str, Any]]:
    """
    Read path_search/summary.yaml and return segments list (empty if not found).
    """
    try:
        if not summary_yaml.exists():
            return []
        data = yaml.safe_load(summary_yaml.read_text(encoding="utf-8")) or {}
        segs = data.get("segments", []) or []
        if not isinstance(segs, list):
            return []
        return segs
    except Exception:
        return []


def _pdb_models_to_coords_and_elems(pdb_path: Path) -> Tuple[List[np.ndarray], List[str]]:
    """
    Return ([coords_model1, coords_model2, ...] in Å), [elements] from a multi-model PDB.
    """
    parser = PDB.PDBParser(QUIET=True)
    st = parser.get_structure("seg", str(pdb_path))
    models = list(st.get_models())
    if not models:
        raise click.ClickException(f"[post] No MODEL found in PDB: {pdb_path}")
    atoms0 = [a for a in models[0].get_atoms()]
    elems: List[str] = []
    for a in atoms0:
        el = (a.element or "").strip()
        if not el:
            nm = a.get_name().strip()
            el = "".join([c for c in nm if c.isalpha()])[:2].title() or "C"
        elems.append(el)
    coords_list: List[np.ndarray] = []
    for m in models:
        atoms = [a for a in m.get_atoms()]
        if len(atoms) != len(atoms0):
            raise click.ClickException(f"[post] Atom count mismatch across models in {pdb_path}")
        coords = np.array([a.get_coord() for a in atoms], dtype=float)
        coords_list.append(coords)
    return coords_list, elems


def _geom_from_angstrom(
    elems: Sequence[str],
    coords_ang: np.ndarray,
    freeze_atoms: Sequence[int],
) -> Any:
    """
    Create a Geometry from Å coordinates using _path_search._new_geom_from_coords (expects Bohr).
    """
    coords_bohr = np.asarray(coords_ang, dtype=float) / BOHR2ANG
    return _path_search._new_geom_from_coords(
        elems,
        coords_bohr,
        coord_type=DEFAULT_COORD_TYPE,
        freeze_atoms=freeze_atoms,
    )


def _load_segment_endpoints(
    path_dir: Path,
    seg_tag: str,
    freeze_atoms: Sequence[int],
) -> Optional[Tuple[Any, Any]]:
    """
    Load left/right endpoints for a segment from
    ``<path_dir>/<seg_tag>_refine_mep/final_geometries.trj``.
    If it does not exist, load structures from
    ``<path_dir>/<seg_tag>_mep/final_geometries.trj``.

    Uses seg_tag (e.g. 'seg_000') and returns (gL_ref, gR_ref).
    """
    base_tag = _path_search._segment_base_id(seg_tag)
    refine_trj = path_dir / f"{base_tag}_refine_mep" / "final_geometries.trj"
    gsm_trj = path_dir / f"{base_tag}_mep" / "final_geometries.trj"

    if refine_trj.exists():
        trj_path = refine_trj
    elif gsm_trj.exists():
        trj_path = gsm_trj
    else:
        return None

    base = str(trj_path)
    gL_ref = geom_loader(
        base + "[0]", coord_type=DEFAULT_COORD_TYPE, freeze_atoms=freeze_atoms
    )
    gR_ref = geom_loader(
        base + "[-1]", coord_type=DEFAULT_COORD_TYPE, freeze_atoms=freeze_atoms
    )

    try:
        fa = np.array(freeze_atoms, dtype=int)
        gL_ref.freeze_atoms = fa
        gR_ref.freeze_atoms = fa
    except Exception:
        pass

    return gL_ref, gR_ref


def _save_single_geom_as_pdb_for_tools(
    g: Any,
    ref_pdb: Path,
    out_dir: Path,
    name: str,
) -> Path:
    """
    Write a single-geometry XYZ trajectory with energy and convert to PDB using the pocket reference
    when possible, for downstream CLI tools. Returns a path to the written structure (PDB when possible,
    otherwise XYZ).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_trj = out_dir / f"{name}.xyz"
    _path_search._write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)

    if ref_pdb.suffix.lower() == ".pdb":
        pdb_out = out_dir / f"{name}.pdb"
        try:
            _path_search._maybe_convert_to_pdb(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
            if pdb_out.exists():
                return pdb_out
        except Exception:
            pass

    return xyz_trj


def _read_xyz_first_last(trj_path: Path) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """
    Lightweight XYZ trajectory reader: return (elements, first_coords[Å], last_coords[Å]).
    Assumes standard multi-frame XYZ: natoms line, comment line, natoms atom lines.
    """
    try:
        lines = trj_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as e:
        raise click.ClickException(f"[irc] Failed to read XYZ/TRJ: {trj_path} ({e})")

    i = 0
    n = len(lines)
    first_elems: Optional[List[str]] = None
    first_coords: Optional[np.ndarray] = None
    last_elems: Optional[List[str]] = None
    last_coords: Optional[np.ndarray] = None

    while i < n:
        try:
            nat = int(lines[i].strip().split()[0])
        except Exception:
            raise click.ClickException(f"[irc] Malformed XYZ/TRJ at line {i + 1}: {trj_path}")
        i += 1
        if i < n:
            i += 1
        if i + nat > n:
            raise click.ClickException(f"[irc] Unexpected EOF while reading frame in {trj_path}")
        elems: List[str] = []
        coords: List[List[float]] = []
        for k in range(nat):
            parts = lines[i + k].split()
            if len(parts) < 4:
                raise click.ClickException(f"[irc] Malformed atom line at {i + k + 1} in {trj_path}")
            elems.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        i += nat
        arr = np.array(coords, dtype=float)
        if first_elems is None or first_coords is None:
            first_elems = elems
            first_coords = arr
        last_elems = elems
        last_coords = arr

    if first_elems is None or first_coords is None or last_elems is None or last_coords is None:
        raise click.ClickException(f"[irc] No frames found in {trj_path}")

    if first_elems != last_elems:
        raise click.ClickException(f"[irc] Element list changed across frames in {trj_path}")

    return first_elems, first_coords, last_coords


def _read_xyz_as_blocks(trj_path: Path) -> List[List[str]]:
    """
    Read a multi-frame XYZ/TRJ file and return a list of frames, each as a list of lines.

    Used for building a concatenated IRC trajectory without parsing coordinates/energies.
    """
    try:
        lines = trj_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as e:
        raise click.ClickException(f"[irc_all] Failed to read XYZ/TRJ: {trj_path} ({e})")

    blocks: List[List[str]] = []
    i = 0
    n = len(lines)
    while i < n:
        # Skip empty lines
        while i < n and not lines[i].strip():
            i += 1
        if i >= n:
            break
        try:
            nat = int(lines[i].split()[0])
        except Exception:
            click.echo(
                f"[irc_all] WARNING: Malformed XYZ/TRJ header at line {i + 1} in {trj_path}; stopping here.",
                err=True,
            )
            break
        start = i
        i += 1  # comment line
        if i < n:
            i += 1
        if i + nat > n:
            click.echo(
                f"[irc_all] WARNING: Unexpected EOF while reading frame from {trj_path}; last frame skipped.",
                err=True,
            )
            break
        end = i + nat
        blocks.append(lines[start:end])
        i = end
    return blocks


def _find_with_suffixes(base_no_ext: Path, suffixes: Sequence[str]) -> Optional[Path]:
    """
    Given a base path without extension, return the first existing file among base.suffix for suffixes.
    """
    for s in suffixes:
        p = base_no_ext.with_suffix(s)
        if p.exists():
            return p
    return None


def _write_segment_energy_diagram(
    prefix: Path,
    labels: List[str],
    energies_au: List[float],
    title_note: str,
    ylabel: str = "ΔE (kcal/mol)",
) -> Optional[Dict[str, Any]]:
    """
    Write energy diagram (PNG) using utils.build_energy_diagram, optionally annotating the title.
    """
    if not energies_au:
        return None
    e0 = energies_au[0]
    energies_kcal = [(e - e0) * AU2KCALPERMOL for e in energies_au]
    fig = build_energy_diagram(
        energies=energies_kcal,
        labels=labels,
        ylabel=ylabel,
        baseline=True,
        showgrid=False,
    )
    if title_note:
        fig.update_layout(title=title_note)
    png = prefix.with_suffix(".png")
    try:
        fig.write_image(str(png), scale=2)
        click.echo(f"[diagram] Wrote energy diagram → {png.name}")
    except Exception as e:
        click.echo(f"[diagram] WARNING: Failed to write energy diagram {png.name}: {e}", err=True)

    payload: Dict[str, Any] = {
        "name": prefix.stem,
        "labels": labels,
        "energies_kcal": energies_kcal,
        "ylabel": ylabel,
    }
    if title_note:
        payload["title"] = title_note
    payload["energies_au"] = list(energies_au)
    payload["image"] = str(png)
    return payload


def _build_global_segment_labels(n_segments: int) -> List[str]:
    """
    Build GSM-like labels for aggregated R/TS/P diagrams over multiple segments.

    Pattern:
      - n = 1: ["R", "TS1", "P"]
      - n ≥ 2: R, TS1, IM1_1, IM1_2, TS2, IM2_1, IM2_2, ..., TSN, P
    """
    if n_segments <= 0:
        return []
    if n_segments == 1:
        return ["R", "TS1", "P"]

    labels: List[str] = []
    for seg_idx in range(1, n_segments + 1):
        if seg_idx == 1:
            labels.extend(["R", "TS1", "IM1_1"])
        elif seg_idx == n_segments:
            labels.extend([f"IM{seg_idx - 1}_2", f"TS{seg_idx}", "P"])
        else:
            labels.extend(
                [f"IM{seg_idx - 1}_2", f"TS{seg_idx}", f"IM{seg_idx}_1"]
            )
    return labels


def _concat_images_horizontally(
    image_paths: Sequence[Path],
    out_path: Path,
    gap: int = 20,
) -> None:
    """
    Concatenate multiple PNG images horizontally into a single PNG.
    Gaps between images are left blank; no interpolation between segments.
    """
    existing = [p for p in image_paths if p is not None and p.exists()]
    if not existing:
        return
    try:
        from PIL import Image
    except Exception:
        click.echo(f"[irc_all] Pillow not available; skipping '{out_path.name}'.", err=True)
        return

    images = [Image.open(str(p)) for p in existing]
    widths = [im.width for im in images]
    heights = [im.height for im in images]
    max_height = max(heights)
    total_width = sum(widths) + gap * (len(images) - 1)
    canvas = Image.new("RGB", (total_width, max_height), "white")
    x = 0
    for im in images:
        y = (max_height - im.height) // 2
        canvas.paste(im, (x, y))
        x += im.width + gap
    _ensure_dir(out_path.parent)
    canvas.save(str(out_path))
    click.echo(f"[irc_all] Wrote aggregated IRC plot → {out_path}")


def _merge_irc_trajectories_to_single_plot(
    trj_and_flags: Sequence[Tuple[Path, bool]],
    out_png: Path,
) -> None:
    """
    Build a single IRC plot over all reactive segments using trj2fig.

    Parameters
    ----------
    trj_and_flags : Sequence[Tuple[Path, bool]]
        For each segment: (finished_irc.trj path, reverse_flag). When reverse_flag is True,
        the frame order of that segment is reversed before concatenation.
    out_png : Path
        Output PNG path for the aggregated plot.
    """
    # Collect blocks from each segment
    all_blocks: List[str] = []
    for trj_path, reverse in trj_and_flags:
        if not isinstance(trj_path, Path) or not trj_path.exists():
            continue
        try:
            blocks = _read_xyz_as_blocks(trj_path)
        except click.ClickException as e:
            click.echo(str(e), err=True)
            continue
        if not blocks:
            continue
        if reverse:
            blocks = list(reversed(blocks))
        all_blocks.extend("\n".join(b) for b in blocks)

    if not all_blocks:
        return

    tmp_trj = out_png.with_suffix(".trj")
    _ensure_dir(tmp_trj.parent)
    try:
        tmp_trj.write_text("\n".join(all_blocks) + "\n", encoding="utf-8")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: Failed to write concatenated IRC trajectory: {e}", err=True)
        return

    try:
        run_trj2fig(tmp_trj, [out_png], unit="kcal", reference="init", reverse_x=False)
        _close_matplotlib_figures()
        click.echo(f"[irc_all] Wrote aggregated IRC plot → {out_png}")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: failed to plot concatenated IRC trajectory: {e}", err=True)
    finally:
        try:
            tmp_trj.unlink()
        except Exception:
            pass


def _optimize_endpoint_geom(
    geom: Any,
    opt_mode_default: str,
    out_dir: Path,
    tag: str,
    dump: bool,
) -> Tuple[Any, Path]:
    """
    Optimize an endpoint geometry using LBFGS/RFO with settings mirroring path_search defaults.

    Args:
        geom: pysisyphus Geometry with calculator attached.
        opt_mode_default: "lbfgs"/"light" or "rfo"/"heavy".
        out_dir: base directory for the optimization outputs.
        tag: tag prefix for the subdirectory.
        dump: whether to dump optimizer trajectory.

    Returns:
        (optimized_geometry, final_xyz_path)
    """
    mode = (opt_mode_default or "light").lower()
    if mode == "light":
        sopt_kind = "lbfgs"
        base_cfg = dict(_path_search.LBFGS_KW)
        OptClass = LBFGS
    else:
        sopt_kind = "rfo"
        base_cfg = dict(_path_search.RFO_KW)
        OptClass = RFOptimizer

    cfg = dict(base_cfg)
    opt_dir = out_dir / f"{tag}_{sopt_kind}_opt"
    _ensure_dir(opt_dir)
    cfg["out_dir"] = str(opt_dir)
    cfg["dump"] = bool(dump)
    max_cycles = int(cfg.get("max_cycles", 300))
    cfg["max_cycles"] = max_cycles

    geom.set_calculator(getattr(geom, "calculator", None))

    click.echo(f"[endpoint-opt] Optimizing '{tag}' with {sopt_kind.upper()} → {opt_dir}")
    opt = OptClass(geom, **cfg)
    try:
        opt.run()
    except (OptimizationError, ZeroStepLength) as e:
        click.echo(
            f"[endpoint-opt] WARNING: optimization for '{tag}' terminated early ({e}); using last geometry.",
            err=True,
        )

    final_xyz = Path(opt.final_fn) if isinstance(opt.final_fn, (str, Path)) else opt.final_fn
    g_final = geom_loader(
        final_xyz,
        coord_type=DEFAULT_COORD_TYPE,
        freeze_atoms=getattr(geom, "freeze_atoms", []),
    )
    try:
        g_final.freeze_atoms = np.array(getattr(geom, "freeze_atoms", []), dtype=int)
    except Exception:
        pass
    g_final.set_calculator(getattr(geom, "calculator", None))
    return g_final, final_xyz


def _run_freq_for_state(
    pdb_path: Path,
    q_int: int,
    spin: int,
    out_dir: Path,
    args_yaml: Optional[Path],
    freeze_links: bool,
    overrides: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Run freq CLI; return parsed thermo dict (may be empty).
    """
    fdir = out_dir
    _ensure_dir(fdir)
    overrides = overrides or {}

    freeze_use = overrides.get("freeze_links")
    if freeze_use is None:
        freeze_use = freeze_links

    dump_use = overrides.get("dump")
    if dump_use is None:
        dump_use = True

    args = [
        "-i",
        str(pdb_path),
        "-q",
        str(int(q_int)),
        "-m",
        str(int(spin)),
        "--freeze-links",
        "True" if freeze_use and (pdb_path.suffix.lower() == ".pdb") else "False",
        "--out-dir",
        str(fdir),
    ]

    _append_cli_arg(args, "--max-write", overrides.get("max_write"))
    _append_cli_arg(args, "--amplitude-ang", overrides.get("amplitude_ang"))
    _append_cli_arg(args, "--n-frames", overrides.get("n_frames"))
    if overrides.get("sort") is not None:
        args.extend(["--sort", str(overrides.get("sort"))])
    _append_cli_arg(args, "--temperature", overrides.get("temperature"))
    _append_cli_arg(args, "--pressure", overrides.get("pressure"))
    _append_cli_arg(args, "--dump", dump_use)

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        args.extend(["--hessian-calc-mode", str(hess_mode)])

    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _saved = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "freq"] + args
        _freq_cli.cli.main(args=args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            click.echo(f"[freq] WARNING: freq exited with code {code}", err=True)
    finally:
        sys.argv = _saved
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _read_imaginary_frequency_info(freq_dir: Path) -> Optional[Dict[str, Any]]:
    """Return diagnostic info about imaginary frequencies if present."""

    freq_file = freq_dir / "frequencies_cm-1.txt"
    if not freq_file.exists():
        return None
    try:
        vals: List[float] = []
        for line in freq_file.read_text(encoding="utf-8").splitlines():
            try:
                tok = line.strip().split()[1]
                vals.append(float(tok))
            except Exception:
                continue
        if not vals:
            return None
        negatives = [v for v in vals if v < 0.0]
        nu_imag = min(negatives) if negatives else None
        min_abs_imag = min((abs(v) for v in negatives), default=None)
        return {
            "n_imag": len(negatives),
            "nu_imag_max_cm": nu_imag,
            "min_abs_imag_cm": min_abs_imag,
            "min_freq_cm": min(vals),
        }
    except Exception:
        return None


def _read_imaginary_frequency(freq_dir: Path) -> Optional[float]:
    """Return the most negative (imaginary) frequency from frequencies_cm-1.txt if present."""

    diag = _read_imaginary_frequency_info(freq_dir)
    if not diag:
        return None
    return diag.get("nu_imag_max_cm") or diag.get("min_freq_cm")


def _run_dft_for_state(
    pdb_path: Path,
    q_int: int,
    spin: int,
    out_dir: Path,
    args_yaml: Optional[Path],
    func_basis: str = "wb97m-v/def2-tzvpd",
    overrides: Optional[Dict[str, Any]] = None,
    engine: str = "gpu",
) -> Dict[str, Any]:
    """
    Run dft CLI; return parsed result.yaml dict (may be empty).
    """
    ddir = out_dir
    _ensure_dir(ddir)
    overrides = overrides or {}

    func_basis_use = overrides.get("func_basis", func_basis)

    args = [
        "-i",
        str(pdb_path),
        "-q",
        str(int(q_int)),
        "-m",
        str(int(spin)),
        "--func-basis",
        str(func_basis_use),
        "--out-dir",
        str(ddir),
    ]
    if engine:
        args.extend(["--engine", str(engine)])

    _append_cli_arg(args, "--max-cycle", overrides.get("max_cycle"))
    _append_cli_arg(args, "--conv-tol", overrides.get("conv_tol"))
    _append_cli_arg(args, "--grid-level", overrides.get("grid_level"))

    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _saved = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "dft"] + args
        _dft_cli.cli.main(args=args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            click.echo(f"[dft] WARNING: dft exited with code {code}", err=True)
    finally:
        sys.argv = _saved
    y = out_dir / "result.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _run_dft_sequence(
    state_jobs: Sequence[Tuple[str, Optional[Path], Path]],
    q_int: int,
    spin: int,
    args_yaml: Optional[Path],
    func_basis: str,
    overrides: Optional[Dict[str, Any]],
    engine: str,
) -> Dict[str, Dict[str, Any]]:
    """Run DFT on a sequence of states."""
    results: Dict[str, Dict[str, Any]] = {}
    for label, pdb_path, out_dir in state_jobs:
        res = _run_dft_for_state(
            pdb_path,
            q_int,
            spin,
            out_dir,
            args_yaml,
            func_basis=func_basis,
            overrides=overrides,
            engine=engine,
        )
        results[label] = res
    return results


def _run_tsopt_on_hei(
    hei_pdb: Path,
    charge: int,
    spin: int,
    calc_cfg: Dict[str, Any],
    args_yaml: Optional[Path],
    out_dir: Path,
    freeze_links: bool,
    opt_mode_default: str,
    overrides: Optional[Dict[str, Any]] = None,
) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI pocket structure; return (final_geom_path, ts_geom).

    Honor the input extension (.pdb/.xyz/.gjf) and pick the first available final_geometry.{pdb|xyz|gjf} output.
    """
    overrides = overrides or {}
    prepared_input = prepare_input_structure(hei_pdb)
    needs_pdb = prepared_input.source_path.suffix.lower() == ".pdb"
    needs_gjf = prepared_input.is_gjf
    ref_pdb = prepared_input.source_path if needs_pdb else None
    ts_dir = _resolve_override_dir(out_dir / "ts", overrides.get("out_dir"))
    _ensure_dir(ts_dir)

    freeze_use = overrides.get("freeze_links")
    if freeze_use is None:
        freeze_use = freeze_links

    opt_mode = overrides.get("opt_mode", opt_mode_default)

    ts_args: List[str] = [
        "-i",
        str(hei_pdb),
        "-q",
        str(int(charge)),
        "-m",
        str(int(spin)),
        "--freeze-links",
        "True" if freeze_use else "False",
        "--out-dir",
        str(ts_dir),
    ]

    if opt_mode is not None:
        ts_args.extend(["--opt-mode", str(opt_mode)])

    _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
    _append_cli_arg(ts_args, "--dump", overrides.get("dump"))
    _append_cli_arg(ts_args, "--thresh", overrides.get("thresh"))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

    if args_yaml is not None:
        ts_args.extend(["--args-yaml", str(args_yaml)])

    click.echo(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
    _saved = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "tsopt"] + ts_args
        _tsopt.cli.main(args=ts_args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            raise click.ClickException(f"[tsopt] tsopt exit code {code}.")
    finally:
        sys.argv = _saved

    ts_pdb = ts_dir / "final_geometry.pdb"
    ts_xyz = ts_dir / "final_geometry.xyz"
    ts_gjf = ts_dir / "final_geometry.gjf"

    ts_geom_path: Optional[Path] = None

    if ts_xyz.exists():
        try:
            convert_xyz_like_outputs(
                ts_xyz,
                prepared_input,
                ref_pdb_path=ref_pdb,
                out_pdb_path=ts_pdb if needs_pdb else None,
                out_gjf_path=ts_gjf if needs_gjf else None,
            )
        except Exception as e:
            click.echo(f"[tsopt] WARNING: Failed to convert TS geometry: {e}", err=True)

    if needs_pdb and ts_pdb.exists():
        ts_geom_path = ts_pdb
    elif ts_xyz.exists():
        ts_geom_path = ts_xyz
    elif ts_pdb.exists():
        ts_geom_path = ts_pdb
    elif needs_gjf and ts_gjf.exists():
        ts_geom_path = ts_gjf
    elif ts_gjf.exists():
        ts_geom_path = ts_gjf
    else:
        raise click.ClickException("[tsopt] TS outputs not found.")

    g_ts = geom_loader(
        ts_geom_path,
        coord_type=DEFAULT_COORD_TYPE,
        freeze_atoms=_freeze_atoms_for_log(),
    )

    calc_args = dict(calc_cfg)
    calc = uma_pysis(**calc_args)
    g_ts.set_calculator(calc)

    prepared_input.cleanup()
    return ts_geom_path, g_ts


def _irc_and_match(
    seg_idx: int,
    seg_dir: Path,
    ref_pdb_for_seg: Path,
    seg_pocket_pdb: Path,
    g_ts: Any,
    q_int: int,
    spin: int,
    freeze_links_flag: bool,
    calc_cfg: Dict[str, Any],
    args_yaml: Optional[Path],
    seg_tag: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run IRC via the irc CLI (EulerPC), then map the IRC endpoints to (left, right).

    - Run irc on the TS structure into ``seg_dir/irc``.
    - Read endpoints from ``finished_irc.trj``.
    - When ``seg_tag`` is provided, map the two IRC endpoints to the corresponding
      GSM segment endpoints loaded from
      ``<path_dir>/<seg_tag>_refine_mep/final_geometries.trj`` (or
      ``<path_dir>/<seg_tag>_mep/final_geometries.trj`` as a fallback). Bond-state
      matching is attempted first; if that fails, the assignment that minimizes the
      sum of RMSDs to the GSM endpoints is used.
    - For TSOPT-only mode (``seg_tag`` is ``None``), the original (first, last)
      IRC endpoints are kept as (left, right).
    - Returns the endpoint geometries, tags, and paths to the per-segment IRC plot
      and ``finished_irc.trj``, together with a flag indicating whether the IRC
      trajectory should be reversed when constructing the global IRC plot.
    """
    freeze_atoms: List[int] = _get_freeze_atoms(seg_pocket_pdb, freeze_links_flag)

    irc_dir = seg_dir / "irc"
    _ensure_dir(irc_dir)

    irc_args: List[str] = [
        "-i",
        str(ref_pdb_for_seg),
        "-q",
        str(int(q_int)),
        "-m",
        str(int(spin)),
        "--out-dir",
        str(irc_dir),
        "--freeze-links",
        "True" if (freeze_links_flag and ref_pdb_for_seg.suffix.lower() == ".pdb") else "False",
    ]

    if args_yaml is not None:
        irc_args.extend(["--args-yaml", str(args_yaml)])
    click.echo(f"[irc] Running EulerPC IRC → out={irc_dir}")
    _saved_argv = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "irc"] + irc_args
        _irc_cli.cli.main(args=irc_args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            raise click.ClickException(f"[irc] irc terminated with exit code {code}.")
    finally:
        sys.argv = _saved_argv

    finished_pdb = irc_dir / "finished_irc.pdb"
    finished_trj = irc_dir / "finished_irc.trj"
    irc_plot = irc_dir / "irc_plot.png"

    # Ensure we have a PDB for visualization if possible
    try:
        if finished_trj.exists() and (not finished_pdb.exists()):
            ref_for_conv: Optional[Path] = None
            if seg_pocket_pdb.suffix.lower() == ".pdb":
                ref_for_conv = seg_pocket_pdb
            elif ref_pdb_for_seg.suffix.lower() == ".pdb":
                ref_for_conv = ref_pdb_for_seg
            if ref_for_conv is not None:
                _path_search._maybe_convert_to_pdb(finished_trj, ref_pdb_path=ref_for_conv, out_path=finished_pdb)
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to convert finished_irc.trj to PDB: {e}", err=True)

    elems, c_first, c_last = _read_xyz_first_last(finished_trj)

    calc_args = dict(calc_cfg)
    shared_calc = uma_pysis(**calc_args)
    g_left = _geom_from_angstrom(elems, c_first, freeze_atoms)
    g_right = _geom_from_angstrom(elems, c_last, freeze_atoms)
    g_left.set_calculator(shared_calc)
    g_right.set_calculator(shared_calc)

    left_tag = "backward"
    right_tag = "forward"
    reverse_irc = False

    path_root = seg_dir.parent

    # Preferred mapping: use endpoints from path_search/<seg_tag>_~~~_mep
    if seg_tag is not None:
        try:
            endpoints = _load_segment_endpoints(path_root, seg_tag, freeze_atoms)
            if endpoints is not None:
                gL_end, gR_end = endpoints
                bond_cfg = dict(_path_search.BOND_KW)

                def _matches(x, y) -> bool:
                    try:
                        changed, _ = _path_search._has_bond_change(x, y, bond_cfg)
                        return not changed
                    except Exception:
                        return False

                L_L = _matches(g_left, gL_end)
                L_R = _matches(g_left, gR_end)
                R_L = _matches(g_right, gL_end)
                R_R = _matches(g_right, gR_end)

                matched_bond = False
                if L_L and R_R:
                    matched_bond = True
                    # orientation already consistent
                elif L_R and R_L:
                    matched_bond = True
                    g_left, g_right = g_right, g_left
                    left_tag, right_tag = right_tag, left_tag
                    reverse_irc = True

                if not matched_bond:
                    # Fallback: minimize total RMSD between (left,right) and (L_end,R_end)
                    try:
                        d_LL = _path_search._rmsd_between(g_left, gL_end, align=True)
                        d_LR = _path_search._rmsd_between(g_left, gR_end, align=True)
                        d_RL = _path_search._rmsd_between(g_right, gL_end, align=True)
                        d_RR = _path_search._rmsd_between(g_right, gR_end, align=True)
                        opt1 = d_LL + d_RR  # left->L_end, right->R_end
                        opt2 = d_LR + d_RL  # left->R_end, right->L_end
                        if opt2 < opt1:
                            g_left, g_right = g_right, g_left
                            left_tag, right_tag = right_tag, left_tag
                            reverse_irc = True
                    except Exception as e:
                        click.echo(
                            f"[irc] WARNING: segment endpoint mapping via RMSD failed: {e}",
                            err=True,
                        )
            else:
                click.echo(
                    f"[irc] WARNING: LBFGS endpoints not found for segment tag '{seg_tag}' under 'segments/'; "
                    "using raw IRC orientation.",
                    err=True,
                )
        except Exception as e:
            click.echo(f"[irc] WARNING: segment endpoint mapping failed: {e}", err=True)
    else:
        # TSOPT-only mode: use raw IRC orientation.
        click.echo(f"[irc] TSOPT-only mode: Use raw irc orientation.")

    # Per-segment IRC plot
    try:
        if finished_trj.exists():
            run_trj2fig(finished_trj, [irc_plot], unit="kcal", reference="init", reverse_x=False)
            _close_matplotlib_figures()
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to plot finished IRC trajectory: {e}", err=True)

    return {
        "left_min_geom": g_left,
        "right_min_geom": g_right,
        "ts_geom": g_ts,
        "left_tag": left_tag,
        "right_tag": right_tag,
        "freeze_atoms": freeze_atoms,
        "irc_plot_path": irc_plot if irc_plot.exists() else None,
        "irc_trj_path": finished_trj if finished_trj.exists() else None,
        "reverse_irc": reverse_irc,
    }


# -----------------------------
# CLI
# -----------------------------


@click.command(
    help=(
        "Run pocket extraction → (optional single-structure staged scan) → MEP search → merge to full PDBs in one shot.\n"
        "If exactly one input is provided: (a) with --scan-lists, run staged scan on the pocket (or full structure "
        "when extraction is skipped) and use stage results as inputs for path_search; "
        "(b) with --tsopt True and no --scan-lists, run TSOPT-only mode."
    ),
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
# ===== Inputs =====
@click.option(
    "-i",
    "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    required=True,
    help=(
        "Two or more **full structures** (PDB/XYZ/GJF) in reaction order (reactant [intermediates ...] product), "
        "or a single **full structure** (with --scan-lists or with --tsopt True). "
        "Extraction (-c/--center) requires PDB inputs. When using --scan-lists without extraction, "
        "the input must be a PDB. You may pass a single '-i' followed by multiple space-separated files "
        "(e.g., '-i A.pdb B.pdb C.pdb')."
    ),
)
@click.option(
    "-c",
    "--center",
    "center_spec",
    type=str,
    required=False,
    default=None,
    help=(
        "Substrate specification for the extractor: "
        "a PDB path, a residue-ID list like '123,124' or 'A:123,B:456' "
        "(insertion codes OK: '123A' / 'A:123A'), "
        "or a residue-name list like 'GPP,MMT'. "
        "When omitted, extraction is skipped and the **full input structure(s)** are used directly as pockets."
    ),
)
@click.option(
    "--out-dir",
    "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("./result_all/"),
    show_default=True,
    help="Top-level output directory for the pipeline.",
)
# ===== Extractor knobs =====
@click.option(
    "-r",
    "--radius",
    type=float,
    default=2.6,
    show_default=True,
    help="Inclusion cutoff (Å) around substrate atoms.",
)
@click.option(
    "--radius-het2het",
    type=float,
    default=0.0,
    show_default=True,
    help="Independent hetero–hetero cutoff (Å) for non‑C/H pairs.",
)
@click.option(
    "--include-H2O",
    "--include-h2o",
    "include_h2o",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Include waters (HOH/WAT/TIP3/SOL) in the pocket.",
)
@click.option(
    "--exclude-backbone",
    "exclude_backbone",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).",
)
@click.option(
    "--add-linkH",
    "add_linkh",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Add link hydrogens for severed bonds (carbon-only) in pockets.",
)
@click.option(
    "--selected_resn",
    type=str,
    default="",
    show_default=True,
    help="Force-include residues (comma/space separated; chain/insertion codes allowed).",
)
@click.option(
    "--ligand-charge",
    type=str,
    default=None,
    help=(
        "Either a total charge (number) to distribute across unknown residues "
        "or a mapping like 'GPP:-3,MMT:-1'."
    ),
)
@click.option(
    "-q",
    "--charge",
    "charge_override",
    type=int,
    default=None,
    help=(
        "Force the total system charge (overrides extractor/GJF/--ligand-charge-derived values; emits a warning when used)."
    ),
)
@click.option(
    "--verbose",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Enable INFO-level logging inside extractor.",
)
# ===== Path search knobs =====
@click.option(
    "-m",
    "--mult",
    "spin",
    type=int,
    default=1,
    show_default=True,
    help="Multiplicity (2S+1).",
)
@click.option(
    "--freeze-links",
    "freeze_links_flag",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="For pocket PDB input, freeze parent atoms of link hydrogens.",
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
)
@click.option(
    "--max-nodes",
    type=int,
    default=10,
    show_default=True,
    help="Max internal nodes for **segment** GSM (String has max_nodes+2 images including endpoints).",
)
@click.option(
    "--max-cycles",
    type=int,
    default=300,
    show_default=True,
    help="Maximum GSM optimization cycles.",
)
@click.option(
    "--climb",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Enable transition-state climbing after growth for the **first** segment in each pair.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help=(
        "Optimizer mode forwarded to scan/tsopt and used for single optimizations: "
        "light (=LBFGS/Dimer) or heavy (=RFO/RSIRFO)."
    ),
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help=(
        "Dump GSM / single-structure trajectories during the run, forwarding the same flag "
        "to scan/tsopt/freq."
    ),
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option(
    "--refine-path",
    "refine_path",
    type=click.BOOL,
    default=True,
    show_default=True,
    help=(
        "If True, run recursive path_search on the full ordered series; if False, run a single-pass "
        "path-opt GSM between each adjacent pair and concatenate the segments (no path_search)."
    ),
)
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "YAML with extra args for path_search (sections: geom, calc, gs, opt, sopt, bond, search)."
    ),
)
@click.option(
    "--preopt",
    "preopt",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If False, skip initial single-structure optimizations of the pocket inputs.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="Common UMA Hessian calculation mode forwarded to tsopt and freq. Defaults to 'FiniteDifference'.",
)
# ===== Post-processing toggles =====
@click.option(
    "--tsopt",
    "do_tsopt",
    type=click.BOOL,
    default=False,
    show_default=True,
    help=(
        "TS optimization + IRC per reactive segment (or TSOPT-only mode for single-structure), "
        "and build energy diagrams."
    ),
)
@click.option(
    "--thermo",
    "do_thermo",
    type=click.BOOL,
    default=False,
    show_default=True,
    help=(
        "Run freq on (R, TS, P) per reactive segment (or TSOPT-only mode) "
        "and build Gibbs free-energy diagram (UMA)."
    ),
)
@click.option(
    "--dft",
    "do_dft",
    type=click.BOOL,
    default=False,
    show_default=True,
    help=(
        "Run DFT single-point on (R, TS, P) and build DFT energy diagram. "
        "With --thermo True, also generate a DFT//UMA Gibbs diagram."
    ),
)
@click.option(
    "--tsopt-max-cycles",
    type=int,
    default=None,
    help="Override tsopt --max-cycles value.",
)
@click.option(
    "--tsopt-out-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help="Override tsopt output subdirectory (relative paths are resolved against the default).",
)
@click.option(
    "--freq-out-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help=(
        "Override freq output base directory (relative paths resolved against the default)."
    ),
)
@click.option(
    "--freq-max-write",
    type=int,
    default=None,
    help="Override freq --max-write value.",
)
@click.option(
    "--freq-amplitude-ang",
    type=float,
    default=None,
    help="Override freq --amplitude-ang (Å).",
)
@click.option(
    "--freq-n-frames",
    type=int,
    default=None,
    help="Override freq --n-frames value.",
)
@click.option(
    "--freq-sort",
    type=click.Choice(["value", "abs"], case_sensitive=False),
    default=None,
    help="Override freq mode sorting.",
)
@click.option(
    "--freq-temperature",
    type=float,
    default=None,
    help="Override freq thermochemistry temperature (K).",
)
@click.option(
    "--freq-pressure",
    type=float,
    default=None,
    help="Override freq thermochemistry pressure (atm).",
)
@click.option(
    "--dft-out-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help=(
        "Override dft output base directory (relative paths resolved against the default)."
    ),
)
@click.option(
    "--dft-func-basis",
    type=str,
    default=None,
    help="Override dft --func-basis value.",
)
@click.option(
    "--dft-max-cycle",
    type=int,
    default=None,
    help="Override dft --max-cycle value.",
)
@click.option(
    "--dft-conv-tol",
    type=float,
    default=None,
    help="Override dft --conv-tol value.",
)
@click.option(
    "--dft-grid-level",
    type=int,
    default=None,
    help="Override dft --grid-level value.",
)
@click.option(
    "--dft-engine",
    type=click.Choice(["gpu", "cpu", "auto"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="Preferred DFT backend: GPU (default), CPU, or auto (try GPU then CPU).",
)
@click.option(
    "--scan-lists",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help=(
        "Python-like list of (i,j,target_Å) per stage for **single-structure** scan. Repeatable; "
        "when repeated, stages run **sequentially**, each starting from the prior stage's relaxed "
        "structure. "
        'Example: "[(12,45,1.35)]" "--scan-lists \'[(10,55,2.20),(23,34,1.80)]\'". '
        "Indices refer to the original full input PDB (1-based). When extraction is used, they are "
        "auto-mapped to the pocket after extraction. Stage results feed into the MEP step (path_search or path_opt)."
    ),
)
@click.option(
    "--scan-out-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help=(
        "Override the scan output directory (default: <out-dir>/scan/). Relative paths are resolved "
        "against the default parent."
    ),
)
@click.option(
    "--scan-one-based",
    type=click.BOOL,
    default=None,
    help=(
        "Override the scan subcommand indexing interpretation (True = 1-based, False = 0-based)."
    ),
)
@click.option(
    "--scan-max-step-size",
    type=float,
    default=None,
    help="Override scan --max-step-size (Å).",
)
@click.option(
    "--scan-bias-k",
    type=float,
    default=None,
    help="Override scan harmonic bias strength k (eV/Å^2).",
)
@click.option(
    "--scan-relax-max-cycles",
    type=int,
    default=None,
    help="Override scan relaxation max cycles per step.",
)
@click.option(
    "--scan-preopt",
    "scan_preopt_override",
    type=click.BOOL,
    default=None,
    help="Override scan --preopt flag.",
)
@click.option(
    "--scan-endopt",
    "scan_endopt_override",
    type=click.BOOL,
    default=None,
    help="Override scan --endopt flag.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    center_spec: Optional[str],
    out_dir: Path,
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    ligand_charge: Optional[str],
    charge_override: Optional[int],
    verbose: bool,
    spin: int,
    freeze_links_flag: bool,
    mep_mode: str,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    refine_path: bool,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    hessian_calc_mode: Optional[str],
    do_tsopt: bool,
    do_thermo: bool,
    do_dft: bool,
    scan_lists_raw: Sequence[str],
    scan_out_dir: Optional[Path],
    scan_one_based: Optional[bool],
    scan_max_step_size: Optional[float],
    scan_bias_k: Optional[float],
    scan_relax_max_cycles: Optional[int],
    scan_preopt_override: Optional[bool],
    scan_endopt_override: Optional[bool],
    tsopt_max_cycles: Optional[int],
    tsopt_out_dir: Optional[Path],
    freq_out_dir: Optional[Path],
    freq_max_write: Optional[int],
    freq_amplitude_ang: Optional[float],
    freq_n_frames: Optional[int],
    freq_sort: Optional[str],
    freq_temperature: Optional[float],
    freq_pressure: Optional[float],
    dft_out_dir: Optional[Path],
    dft_func_basis: Optional[str],
    dft_max_cycle: Optional[int],
    dft_conv_tol: Optional[float],
    dft_grid_level: Optional[int],
    dft_engine: str,
) -> None:
    """
    The **all** command composes `extract` → (optional `scan` on pocket or full input) → MEP search
    (`path_search` with ``--refine-path True`` or concatenated `path-opt` otherwise) and hides ref-template
    bookkeeping. It also accepts the sloppy `-i A B C` style like `path_search` does.
    With single input:
      - with --scan-lists: run staged scan and use stage results as inputs for path_search,
      - with --tsopt True and no --scan-lists: run TSOPT-only mode (no path_search).

    With ``--refine-path True``, the recursive ``path_search`` workflow is used. When ``False`` (default),
    a single-pass ``path-opt`` GSM is run between each adjacent pair of inputs and the segments are
    concatenated into the final MEP without invoking ``path_search``.
    """
    set_convert_file_enabled(convert_files)
    command_str = " ".join(sys.argv)
    time_start = time.perf_counter()
    energy_diagrams: List[Dict[str, Any]] = []

    dump_override_requested = False
    try:
        dump_source = ctx.get_parameter_source("dump")
        dump_override_requested = dump_source not in (None, ParameterSource.DEFAULT)
    except Exception:
        dump_override_requested = False

    argv_all = sys.argv[1:]
    i_vals = _collect_option_values(argv_all, ("-i", "--input"))
    if i_vals:
        i_parsed: List[Path] = []
        for tok in i_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Input path '{tok}' not found or is a directory. "
                    "When using '-i', list only existing file paths (multiple paths may follow a single '-i')."
                )
            i_parsed.append(p)
        input_paths = tuple(i_parsed)

    is_single = len(input_paths) == 1
    has_scan = bool(scan_lists_raw)
    single_tsopt_mode = is_single and (not has_scan) and do_tsopt

    if (len(input_paths) < 2) and not (is_single and (has_scan or do_tsopt)):
        raise click.BadParameter(
            "Provide at least two structures with -i/--input in reaction order, "
            "or use a single structure with --scan-lists, or a single structure with --tsopt True."
        )

    tsopt_opt_mode_default = opt_mode.lower()
    tsopt_overrides: Dict[str, Any] = {}
    if tsopt_max_cycles is not None:
        tsopt_overrides["max_cycles"] = int(tsopt_max_cycles)
    if dump_override_requested:
        tsopt_overrides["dump"] = bool(dump)
    if tsopt_out_dir is not None:
        tsopt_overrides["out_dir"] = tsopt_out_dir
    if hessian_calc_mode is not None:
        tsopt_overrides["hessian_calc_mode"] = hessian_calc_mode
    if thresh is not None:
        tsopt_overrides["thresh"] = str(thresh)

    freq_overrides: Dict[str, Any] = {}
    if freq_max_write is not None:
        freq_overrides["max_write"] = int(freq_max_write)
    if freq_amplitude_ang is not None:
        freq_overrides["amplitude_ang"] = float(freq_amplitude_ang)
    if freq_n_frames is not None:
        freq_overrides["n_frames"] = int(freq_n_frames)
    if freq_sort is not None:
        freq_overrides["sort"] = freq_sort.lower()
    if freq_temperature is not None:
        freq_overrides["temperature"] = float(freq_temperature)
    if freq_pressure is not None:
        freq_overrides["pressure"] = float(freq_pressure)
    if dump_override_requested:
        freq_overrides["dump"] = bool(dump)
    if hessian_calc_mode is not None:
        freq_overrides["hessian_calc_mode"] = hessian_calc_mode

    dft_overrides: Dict[str, Any] = {}
    if dft_max_cycle is not None:
        dft_overrides["max_cycle"] = int(dft_max_cycle)
    if dft_conv_tol is not None:
        dft_overrides["conv_tol"] = float(dft_conv_tol)
    if dft_grid_level is not None:
        dft_overrides["grid_level"] = int(dft_grid_level)

    dft_func_basis_use = dft_func_basis or "wb97m-v/def2-tzvpd"

    yaml_cfg = load_yaml_dict(args_yaml)

    skip_extract = center_spec is None or str(center_spec).strip() == ""

    out_dir = out_dir.resolve()
    pockets_dir = out_dir / "pockets"
    path_dir = out_dir / ("path_search" if refine_path else "path_opt")
    scan_dir = _resolve_override_dir(out_dir / "scan", scan_out_dir)
    stage_total = 3
    _ensure_dir(out_dir)
    if not skip_extract:
        _ensure_dir(pockets_dir)
    if not single_tsopt_mode:
        _ensure_dir(path_dir)

    elem_tmp_dir = out_dir / "add_elem_info"
    inputs_for_extract: List[Path] = []
    elem_fix_echo = False
    for p in input_paths:
        if _pdb_needs_elem_fix(p):
            if not elem_fix_echo:
                click.echo(
                    "\n=== [all] Preflight — add_elem_info (only when element fields are missing) ===\n"
                )
                elem_fix_echo = True
            _ensure_dir(elem_tmp_dir)
            out_p = (elem_tmp_dir / p.name).resolve()
            try:
                _assign_elem_info(str(p), str(out_p), overwrite=False)
                click.echo(f"[all] add_elem_info: fixed elements → {out_p}")
                inputs_for_extract.append(out_p)
            except SystemExit as e:
                code = getattr(e, "code", 1)
                click.echo(
                    f"[all] WARNING: add_elem_info exited with code {code} for {p}; using original.",
                    err=True,
                )
                inputs_for_extract.append(p.resolve())
            except Exception as e:
                click.echo(
                    f"[all] WARNING: add_elem_info failed for {p}: {e} — using original file.",
                    err=True,
                )
                inputs_for_extract.append(p.resolve())
        else:
            inputs_for_extract.append(p.resolve())

    extract_inputs = tuple(inputs_for_extract)

    pocket_outputs: List[Path] = []
    if not skip_extract:
        for p in extract_inputs:
            pocket_outputs.append((pockets_dir / f"pocket_{p.stem}.pdb").resolve())

    q_from_flow: Optional[int] = None

    if not skip_extract:
        click.echo(
            f"\n=== [all] Stage 1/{stage_total} — Active-site pocket extraction (multi-structure union when applicable) ===\n"
        )
        try:
            ex_res = extract_api(
                complex_pdb=[str(p) for p in extract_inputs],
                center=center_spec,
                output=[str(p) for p in pocket_outputs],
                radius=float(radius),
                radius_het2het=float(radius_het2het),
                include_H2O=bool(include_h2o),
                exclude_backbone=bool(exclude_backbone),
                add_linkH=bool(add_linkh),
                selected_resn=selected_resn or "",
                ligand_charge=ligand_charge,
                verbose=bool(verbose),
            )
        except Exception as e:
            raise click.ClickException(f"[all] Extractor failed: {e}")

        click.echo("[all] Pocket files:")
        for op in pocket_outputs:
            click.echo(f"  - {op}")

        try:
            cs = ex_res.get("charge_summary", {})
            q_total = float(cs.get("total_charge", 0.0))
            q_prot = float(cs.get("protein_charge", 0.0))
            q_lig = float(cs.get("ligand_total_charge", 0.0))
            q_ion = float(cs.get("ion_total_charge", 0.0))
            click.echo("\n[all] Charge summary from extractor (model #1):")
            click.echo(
                f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}"
            )
            q_from_flow = _round_charge_with_note(q_total)
        except Exception as e:
            raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")
    else:
        click.echo(
            f"\n=== [all] Stage 1/{stage_total} — Extraction skipped (no -c/--center); using FULL structures as pockets ===\n"
        )
        first_input = input_paths[0].resolve()
        gjf_charge: Optional[int] = None
        gjf_spin: Optional[int] = None
        if first_input.suffix.lower() == ".gjf":
            try:
                with prepare_input_structure(first_input) as prepared:
                    gjf_charge, gjf_spin = resolve_charge_spin_or_raise(
                        prepared, charge=None, spin=None
                    )
                click.echo(
                    f"[all] Detected from GJF (first input): charge={gjf_charge:+d}, spin={gjf_spin}"
                )
            except Exception as e:
                click.echo(
                    f"[all] NOTE: failed to parse charge/spin from GJF '{first_input.name}': {e}",
                    err=True,
                )

        user_provided_spin = True
        try:
            spin_source = ctx.get_parameter_source("spin")
            user_provided_spin = spin_source not in (None, ParameterSource.DEFAULT)
        except Exception:
            user_provided_spin = True

        q_total_fallback: float
        numeric_ligand_charge: Optional[float] = None
        if ligand_charge is not None:
            try:
                numeric_ligand_charge = float(ligand_charge)
            except Exception:
                numeric_ligand_charge = None

        if numeric_ligand_charge is not None:
            q_total_fallback = numeric_ligand_charge
            click.echo(
                f"[all] Using --ligand-charge as TOTAL system charge: {q_total_fallback:+g}"
            )
        elif gjf_charge is not None:
            q_total_fallback = float(gjf_charge)
            click.echo(f"[all] Using total charge from first GJF: {q_total_fallback:+g}")
        else:
            q_total_fallback = 0.0
            if ligand_charge is not None:
                click.echo(
                    "[all] NOTE: non-numeric --ligand-charge is ignored without --center; "
                    "defaulting total charge to 0 (no GJF charge available).",
                    err=True,
                )
            else:
                if charge_override is None:
                    click.echo(
                        "[all] NOTE: No total charge provided; defaulting to 0. "
                        "Supply '--ligand-charge <number>' to override."
                    )
        q_from_flow = _round_charge_with_note(q_total_fallback)

        if (not user_provided_spin) and (gjf_spin is not None):
            spin = int(gjf_spin)
            click.echo(f"[all] Spin multiplicity set from GJF: {spin}")

    if charge_override is not None:
        q_int = int(charge_override)
        override_msg = (
            f"[all] WARNING: -q/--charge override supplied; forcing TOTAL system charge to {q_int:+d}"
        )
        if q_from_flow is not None:
            override_msg += f" (would otherwise use {int(q_from_flow):+d} from workflow)"
        click.echo(override_msg, err=True)
    else:
        q_int = int(q_from_flow) if q_from_flow is not None else 0

    calc_cfg_shared = _build_calc_cfg(q_int, spin, yaml_cfg)

    # -------------------------------------------------------------------------
    # TSOPT-only single-structure mode
    # -------------------------------------------------------------------------
    if single_tsopt_mode:
        click.echo("\n=== [all] TSOPT-only single-structure mode ===\n")
        tsroot = out_dir / "tsopt_single"
        _ensure_dir(tsroot)

        # In TSOPT-only mode, no MEP search is performed. Use a placeholder for
        # MEP-related fields in downstream summaries.
        mep_mode_kind = "---"

        ts_initial_pdb = pocket_outputs[0] if not skip_extract else input_paths[0].resolve()

        ts_pdb, g_ts = _run_tsopt_on_hei(
            ts_initial_pdb,
            q_int,
            spin,
            calc_cfg_shared,
            args_yaml,
            tsroot,
            freeze_links_flag,
            tsopt_opt_mode_default,
            overrides=tsopt_overrides,
        )

        irc_res = _irc_and_match(
            seg_idx=1,
            seg_dir=tsroot,
            ref_pdb_for_seg=ts_pdb,
            seg_pocket_pdb=ts_initial_pdb,
            g_ts=g_ts,
            q_int=q_int,
            spin=spin,
            freeze_links_flag=freeze_links_flag,
            calc_cfg=calc_cfg_shared,
            args_yaml=args_yaml,
            seg_tag=None,
        )
        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        irc_plot_path = irc_res.get("irc_plot_path")

        eL = float(gL.energy)
        eR_raw = float(gR.energy)
        eT = float(gT.energy)

        if eL >= eR_raw:
            g_react_irc, e_react_irc = gL, eL
            g_prod_irc, e_prod_irc = gR, eR_raw
        else:
            g_react_irc, e_react_irc = gR, eR_raw
            g_prod_irc, e_prod_irc = gL, eL

        struct_dir = tsroot / "structures"
        _ensure_dir(struct_dir)
        pocket_ref = ts_initial_pdb
        pR_irc = _save_single_geom_as_pdb_for_tools(
            g_react_irc, pocket_ref, struct_dir, "reactant_irc"
        )
        pP_irc = _save_single_geom_as_pdb_for_tools(
            g_prod_irc, pocket_ref, struct_dir, "product_irc"
        )
        pT = _save_single_geom_as_pdb_for_tools(gT, pocket_ref, struct_dir, "ts")

        endpoint_opt_dir = tsroot / "endpoint_opt"
        _ensure_dir(endpoint_opt_dir)
        try:
            g_react_opt, _ = _optimize_endpoint_geom(
                g_react_irc,
                tsopt_opt_mode_default,
                endpoint_opt_dir,
                "reactant",
                dump=dump,
            )
        except Exception as e:
            click.echo(
                f"[post] WARNING: Reactant endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            g_react_opt = g_react_irc
        try:
            g_prod_opt, _ = _optimize_endpoint_geom(
                g_prod_irc,
                tsopt_opt_mode_default,
                endpoint_opt_dir,
                "product",
                dump=dump,
            )
        except Exception as e:
            click.echo(
                f"[post] WARNING: Product endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            g_prod_opt = g_prod_irc

        # Clean up endpoint_opt as a temporary working directory
        shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
        click.echo(f"[endpoint-opt] Clean endpoint-opt working dir.") 

        pR = _save_single_geom_as_pdb_for_tools(
            g_react_opt, pocket_ref, struct_dir, "reactant"
        )
        pP = _save_single_geom_as_pdb_for_tools(
            g_prod_opt, pocket_ref, struct_dir, "product"
        )

        e_react = float(g_react_opt.energy)
        e_prod = float(g_prod_opt.energy)

        diag_payload = _write_segment_energy_diagram(
            tsroot / "energy_diagram_UMA",
            labels=["R", "TS", "P"],
            energies_au=[e_react, eT, e_prod],
            title_note="(UMA, TSOPT + IRC)",
        )
        if diag_payload:
            energy_diagrams.append(diag_payload)

        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        freq_root = _resolve_override_dir(tsroot / "freq", freq_out_dir)
        dft_root = _resolve_override_dir(tsroot / "dft", dft_out_dir)

        if do_thermo:
            click.echo("[thermo] Single TSOPT: freq on R/TS/P")
            tR = _run_freq_for_state(
                pR, q_int, spin, freq_root / "R", args_yaml, freeze_links_flag, overrides=freq_overrides
            )
            tT = _run_freq_for_state(
                pT, q_int, spin, freq_root / "TS", args_yaml, freeze_links_flag, overrides=freq_overrides
            )
            tP = _run_freq_for_state(
                pP, q_int, spin, freq_root / "P", args_yaml, freeze_links_flag, overrides=freq_overrides
            )
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            try:
                GR = float(
                    tR.get("sum_EE_and_thermal_free_energy_ha", e_react)
                )
                GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", eT))
                GP = float(
                    tP.get("sum_EE_and_thermal_free_energy_ha", e_prod)
                )
                diag_payload = _write_segment_energy_diagram(
                    tsroot / "energy_diagram_G_UMA",
                    labels=["R", "TS", "P"],
                    energies_au=[GR, GT, GP],
                    title_note="(UMA + Thermal Correction)",
                    ylabel="ΔG (kcal/mol)",
                )
                if diag_payload:
                    energy_diagrams.append(diag_payload)
            except Exception as e:
                click.echo(
                    f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True
                )

        if do_dft:
            click.echo("[dft] Single TSOPT: DFT on R/TS/P")
            dft_jobs = [
                ("R", pR, dft_root / "R"),
                ("TS", pT, dft_root / "TS"),
                ("P", pP, dft_root / "P"),
            ]
            dft_payloads = _run_dft_sequence(
                dft_jobs,
                q_int,
                spin,
                args_yaml,
                dft_func_basis_use,
                dft_overrides,
                dft_engine,
            )
            dR = dft_payloads.get("R")
            dT = dft_payloads.get("TS")
            dP = dft_payloads.get("P")
            try:
                eR_dft = float(
                    ((dR or {}).get("energy", {}) or {}).get("hartree", e_react)
                )
                eT_dft = float(
                    ((dT or {}).get("energy", {}) or {}).get("hartree", eT)
                )
                eP_dft = float(
                    ((dP or {}).get("energy", {}) or {}).get("hartree", e_prod)
                )
                diag_payload = _write_segment_energy_diagram(
                    tsroot / "energy_diagram_DFT",
                    labels=["R", "TS", "P"],
                    energies_au=[eR_dft, eT_dft, eP_dft],
                    title_note=f"({dft_func_basis_use} // UMA)",
                )
                if diag_payload:
                    energy_diagrams.append(diag_payload)
            except Exception as e:
                click.echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            if do_thermo:
                try:
                    dG_R = float(
                        (thermo_payloads.get("R", {}) or {}).get(
                            "thermal_correction_free_energy_ha", 0.0
                        )
                    )
                    dG_T = float(
                        (thermo_payloads.get("TS", {}) or {}).get(
                            "thermal_correction_free_energy_ha", 0.0
                        )
                    )
                    dG_P = float(
                        (thermo_payloads.get("P", {}) or {}).get(
                            "thermal_correction_free_energy_ha", 0.0
                        )
                    )
                    GR_dftUMA = eR_dft + dG_R
                    GT_dftUMA = eT_dft + dG_T
                    GP_dftUMA = eP_dft + dG_P
                    diag_payload = _write_segment_energy_diagram(
                        tsroot / "energy_diagram_G_DFT_plus_UMA",
                        labels=["R", "TS", "P"],
                        energies_au=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                        title_note=f"({dft_func_basis_use} // UMA + Thermal Correction)",
                        ylabel="ΔG (kcal/mol)",
                    )
                    if diag_payload:
                        energy_diagrams.append(diag_payload)
                except Exception as e:
                    click.echo(
                        f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}",
                        err=True,
                    )

        # Build summary.yaml for TSOPT-only runs, mirroring path_search/path_opt outputs
        bond_cfg = dict(_path_search.BOND_KW)
        bond_summary = ""
        try:
            changed, bond_summary = _path_search._has_bond_change(
                g_react_opt, g_prod_opt, bond_cfg
            )
            if not changed:
                bond_summary = "(no covalent changes detected)"
        except Exception as e:
            click.echo(
                f"[post] WARNING: Failed to detect bond changes for TSOPT-only endpoints: {e}",
                err=True,
            )
            bond_summary = "(no covalent changes detected)"

        barrier = (eT - e_react) * AU2KCALPERMOL
        delta = (e_prod - e_react) * AU2KCALPERMOL

        n_images = 0
        try:
            irc_trj_path = irc_res.get("irc_trj_path")
            if isinstance(irc_trj_path, Path) and irc_trj_path.exists():
                n_images = len(_read_xyz_as_blocks(irc_trj_path))
        except Exception:
            n_images = 0

        summary = {
            "out_dir": str(tsroot),
            "n_images": n_images,
            "n_segments": 1,
            "segments": [
                {
                    "index": 1,
                    "tag": "seg_01",
                    "kind": "tsopt",
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": _path_search._bond_changes_block(bond_summary),
                }
            ],
        }
        if energy_diagrams:
            summary["energy_diagrams"] = list(energy_diagrams)
        try:
            with open(tsroot / "summary.yaml", "w") as f:
                yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
            click.echo(f"[write] Wrote '{tsroot / 'summary.yaml'}'.")
            try:
                dst_summary = out_dir / "summary.yaml"
                shutil.copy2(tsroot / "summary.yaml", dst_summary)
                click.echo(f"[all] Copied summary.yaml → {dst_summary}")
            except Exception as e:
                click.echo(
                    f"[all] WARNING: Failed to mirror summary.yaml in TSOPT-only mode: {e}",
                    err=True,
                )
            try:
                ts_freq_info = (
                    _read_imaginary_frequency_info(freq_root / "TS") if do_thermo else None
                )
                segment_log = {
                    "index": 1,
                    "tag": "seg_01",
                    "kind": "tsopt",
                    "bond_changes": summary["segments"][0].get("bond_changes"),
                    "mep_barrier_kcal": barrier,
                    "mep_delta_kcal": delta,
                    "post_dir": str(tsroot),
                    "irc_plot": str(irc_plot_path) if isinstance(irc_plot_path, Path) else None,
                    "irc_traj": str(irc_trj_path) if isinstance(irc_trj_path, Path) else None,
                }
                if ts_freq_info is not None:
                    segment_log["ts_imag"] = ts_freq_info
                    if ts_freq_info.get("nu_imag_max_cm") is not None:
                        segment_log["ts_imag_freq_cm"] = ts_freq_info["nu_imag_max_cm"]
                segment_log["uma"] = {
                    "labels": ["R", "TS", "P"],
                    "energies_au": [e_react, eT, e_prod],
                    "energies_kcal": [0.0, (eT - e_react) * AU2KCALPERMOL, (e_prod - e_react) * AU2KCALPERMOL],
                    "diagram": str(tsroot / "energy_diagram_UMA.png"),
                    "structures": {"R": pR, "TS": pT, "P": pP},
                    "barrier_kcal": barrier,
                    "delta_kcal": delta,
                }
                if do_thermo and thermo_payloads:
                    GR = float(thermo_payloads.get("R", {}).get("sum_EE_and_thermal_free_energy_ha", e_react))
                    GT = float(thermo_payloads.get("TS", {}).get("sum_EE_and_thermal_free_energy_ha", eT))
                    GP = float(thermo_payloads.get("P", {}).get("sum_EE_and_thermal_free_energy_ha", e_prod))
                    segment_log["gibbs_uma"] = {
                        "labels": ["R", "TS", "P"],
                        "energies_au": [GR, GT, GP],
                        "energies_kcal": [
                            0.0,
                            (GT - GR) * AU2KCALPERMOL,
                            (GP - GR) * AU2KCALPERMOL,
                        ],
                        "diagram": str(tsroot / "energy_diagram_G_UMA.png"),
                        "structures": {"R": pR, "TS": pT, "P": pP},
                        "barrier_kcal": (GT - GR) * AU2KCALPERMOL,
                        "delta_kcal": (GP - GR) * AU2KCALPERMOL,
                    }
                if do_dft:
                    segment_log["dft"] = {
                        "labels": ["R", "TS", "P"],
                        "energies_au": [eR_dft, eT_dft, eP_dft],
                        "energies_kcal": [
                            0.0,
                            (eT_dft - eR_dft) * AU2KCALPERMOL,
                            (eP_dft - eR_dft) * AU2KCALPERMOL,
                        ],
                        "diagram": str(tsroot / "energy_diagram_DFT.png"),
                        "structures": {"R": pR, "TS": pT, "P": pP},
                        "barrier_kcal": (eT_dft - eR_dft) * AU2KCALPERMOL,
                        "delta_kcal": (eP_dft - eR_dft) * AU2KCALPERMOL,
                    }
                    if do_thermo:
                        segment_log["gibbs_dft_uma"] = {
                            "labels": ["R", "TS", "P"],
                            "energies_au": [GR_dftUMA, GT_dftUMA, GP_dftUMA],
                            "energies_kcal": [
                                0.0,
                                (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                                (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                            ],
                            "diagram": str(tsroot / "energy_diagram_G_DFT_plus_UMA.png"),
                            "structures": {"R": pR, "TS": pT, "P": pP},
                            "barrier_kcal": (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                            "delta_kcal": (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                        }

                summary_payload = {
                    "root_out_dir": str(out_dir),
                    "path_dir": str(tsroot),
                    "path_module_dir": "tsopt_single",
                    "pipeline_mode": "tsopt-only",
                    "refine_path": refine_path,
                    "tsopt": do_tsopt,
                    "thermo": do_thermo,
                    "dft": do_dft,
                    "opt_mode": tsopt_opt_mode_default,
                    "mep_mode": mep_mode_kind,
                    "uma_model": calc_cfg_shared.get("model"),
                    "command": command_str,
                    "charge": q_int,
                    "spin": spin,
                    "freeze_atoms": _freeze_atoms_for_log(),
                    "mep": {"n_images": n_images, "n_segments": 1},
                    "segments": summary.get("segments", []),
                    "energy_diagrams": summary.get("energy_diagrams", []),
                    "post_segments": [segment_log],
                    "key_files": {
                        "summary.yaml": "YAML-format summary",
                        "summary.log": "This summary",
                        "energy_diagram_UMA_all.png": "UMA R–TS–P energies",
                    },
                }
                write_summary_log(tsroot / "summary.log", summary_payload)
                try:
                    shutil.copy2(tsroot / "summary.log", out_dir / "summary.log")
                except Exception:
                    pass
            except Exception as e:
                click.echo(f"[write] WARNING: Failed to write summary.log in TSOPT-only mode: {e}", err=True)
        except Exception as e:
            click.echo(
                f"[write] WARNING: Failed to write summary.yaml for TSOPT-only run: {e}",
                err=True,
            )

        try:
            for stem in (
                "energy_diagram_UMA",
                "energy_diagram_G_UMA",
                "energy_diagram_DFT",
                "energy_diagram_G_DFT_plus_UMA",
            ):
                src = tsroot / f"{stem}.png"
                if src.exists():
                    dst = out_dir / f"{stem}_all.png"
                    shutil.copy2(src, dst)
                    click.echo(f"[all] Copied {src.name} → {dst}")
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to mirror *_all diagrams in TSOPT-only mode: {e}",
                err=True,
            )

        try:
            if isinstance(irc_plot_path, Path) and irc_plot_path.exists():
                dst = out_dir / "irc_plot_all.png"
                shutil.copy2(irc_plot_path, dst)
                click.echo(f"[all] Copied IRC plot → {dst}")
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to mirror IRC plot in TSOPT-only mode: {e}",
                err=True,
            )

        try:
            if pT.suffix.lower() == ".pdb":
                ts_copy = out_dir / "ts_seg_01.pdb"
                shutil.copy2(pT, ts_copy)
            else:
                ts_copy = out_dir / "ts_seg_01.xyz"
                _path_search._write_xyz_trj_with_energy(
                    [gT], [float(gT.energy)], ts_copy
                )
            click.echo(
                f"[all] Copied TS structure for TSOPT-only run → {ts_copy}"
            )
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to write TS structure in TSOPT-only mode: {e}",
                err=True,
            )

        click.echo(
            "\n=== [all] TSOPT-only pipeline finished successfully ===\n"
        )
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # -------------------------------------------------------------------------
    # Stage 1b: optional staged scan (single-structure)
    # -------------------------------------------------------------------------
    pockets_for_path: List[Path]
    if is_single and has_scan and skip_extract and input_paths[0].suffix.lower() != ".pdb":
        raise click.ClickException(
            "[all] When using --scan-lists while skipping extraction, the input must be a PDB file. "
            "Specify -c/--center to build a pocket PDB before running the scan."
        )

    if is_single and has_scan:
        click.echo("\n=== [all] Stage 1b — Staged scan on input ===\n")
        _ensure_dir(scan_dir)

        if skip_extract:
            scan_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _parse_scan_lists_literals(scan_lists_raw)
        else:
            scan_input_pdb = Path(pocket_outputs[0]).resolve()
            full_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _convert_scan_lists_to_pocket_indices(
                scan_lists_raw, full_input_pdb, scan_input_pdb
            )
            click.echo(
                "[all] Remapped --scan-lists indices from the full PDB to the pocket ordering."
            )

        scan_one_based_effective = True if scan_one_based is None else bool(
            scan_one_based
        )
        scan_stage_literals: List[str] = []
        for stage in converted_scan_stages:
            if scan_one_based_effective:
                stage_use = stage
            else:
                stage_use = [(i - 1, j - 1, target) for (i, j, target) in stage]
            scan_stage_literals.append(_format_scan_stage(stage_use))

        scan_preopt_use = preopt if scan_preopt_override is None else bool(
            scan_preopt_override
        )
        scan_endopt_use = True if scan_endopt_override is None else bool(
            scan_endopt_override
        )
        scan_opt_mode_use = opt_mode.lower()

        scan_args: List[str] = [
            "-i",
            str(scan_input_pdb),
            "-q",
            str(int(q_int)),
            "-m",
            str(int(spin)),
            "--out-dir",
            str(scan_dir),
            "--freeze-links",
            "True" if freeze_links_flag else "False",
            "--preopt",
            "True" if scan_preopt_use else "False",
            "--endopt",
            "True" if scan_endopt_use else "False",
            "--opt-mode",
            str(scan_opt_mode_use),
        ]

        if dump_override_requested:
            scan_args.extend(
                ["--dump", "True" if dump else "False"]
            )

        if scan_one_based is not None:
            scan_args.append("--one-based" if scan_one_based else "--zero-based")

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        if args_yaml is not None:
            scan_args.extend(["--args-yaml", str(args_yaml)])
        for literal in scan_stage_literals:
            scan_args.extend(["--scan-lists", literal])

        click.echo("[all] Invoking scan with arguments:")
        click.echo("  " + " ".join(scan_args))

        _saved_argv = list(sys.argv)
        try:
            sys.argv = ["pdb2reaction", "scan"] + scan_args
            _scan_cli.cli.main(args=scan_args, standalone_mode=False)
        except SystemExit as e:
            code = getattr(e, "code", 1)
            if code not in (None, 0):
                raise click.ClickException(
                    f"[all] scan terminated with exit code {code}."
                )
        except Exception as e:
            raise click.ClickException(f"[all] scan failed: {e}")
        finally:
            sys.argv = _saved_argv

        stage_results: List[Path] = []
        for st in sorted(scan_dir.glob("stage_*")):
            res = _find_with_suffixes(st / "result", [".pdb", ".xyz", ".gjf"])
            if res:
                stage_results.append(res.resolve())
        if not stage_results:
            raise click.ClickException(
                "[all] No stage result structures found under scan/ "
                "(looked for result.[pdb|xyz|gjf])."
            )
        click.echo("[all] Collected scan stage pocket files:")
        for p in stage_results:
            click.echo(f"  - {p}")

        pockets_for_path = [scan_input_pdb] + stage_results
    else:
        if skip_extract:
            pockets_for_path = [p.resolve() for p in inputs_for_extract]
        else:
            pockets_for_path = list(pocket_outputs)

    # Determine availability of full-system templates for downstream merge/copies
    def _is_pdb(path: Path) -> bool:
        return path.suffix.lower() == ".pdb"

    gave_ref_pdb = False

    mep_mode_kind = mep_mode.strip().lower()

    if skip_extract:
        click.echo(
            "[all] NOTE: skipping --ref-pdb (no --center; inputs already represent full structures)."
        )
    elif is_single and has_scan:
        if _is_pdb(input_paths[0]):
            gave_ref_pdb = True
        else:
            click.echo(
                "[all] NOTE: skipping --ref-pdb (single+scan: original input is not a PDB)."
            )
    else:
        if all(_is_pdb(p) for p in input_paths):
            gave_ref_pdb = True
        else:
            click.echo(
                "[all] NOTE: skipping --ref-pdb (one or more original inputs are not PDB)."
            )

    # -------------------------------------------------------------------------
    # Stage 2: MEP search
    # -------------------------------------------------------------------------
    if not refine_path:
        click.echo(
            f"\n=== [all] Stage 2/{stage_total} — Pairwise MEP search via path-opt (no recursive path_search) ===\n"
        )

        if len(pockets_for_path) < 2:
            raise click.ClickException("[all] Need at least two structures for path-opt MEP concatenation.")

        combined_blocks: List[str] = []
        path_opt_segments: List[Dict[str, Any]] = []
        for idx, (pL, pR) in enumerate(zip(pockets_for_path, pockets_for_path[1:]), start=1):
            seg_dir = (path_dir / f"seg_{idx:03d}_mep").resolve()
            seg_tag = f"seg_{idx:03d}"
            po_args: List[str] = [
                "-i",
                str(pL),
                str(pR),
                "-q",
                str(q_int),
                "-m",
                str(int(spin)),
                "--freeze-links",
                "True" if freeze_links_flag else "False",
                "--mep-mode",
                mep_mode_kind,
                "--max-nodes",
                str(int(max_nodes)),
                "--max-cycles",
                str(int(max_cycles)),
                "--climb",
                "True" if climb else "False",
                "--opt-mode",
                str(opt_mode),
                "--dump",
                "True" if dump else "False",
                "--out-dir",
                str(seg_dir),
                "--preopt",
                "True" if preopt else "False",
            ]
            if thresh is not None:
                po_args.extend(["--thresh", str(thresh)])
            if args_yaml is not None:
                po_args.extend(["--args-yaml", str(args_yaml)])

            click.echo(f"[all] Invoking path-opt for segment {idx}:")
            click.echo("  " + " ".join(po_args))

            _saved_argv = list(sys.argv)
            try:
                sys.argv = ["pdb2reaction", "path-opt"] + po_args
                _path_opt.cli.main(args=po_args, standalone_mode=False)
            except SystemExit as e:
                code = getattr(e, "code", 1)
                if code not in (None, 0):
                    raise click.ClickException(
                        f"[all] path-opt terminated with exit code {code} (segment {idx})."
                    )
            except Exception as e:
                raise click.ClickException(f"[all] path-opt failed for segment {idx}: {e}")
            finally:
                sys.argv = _saved_argv

            seg_trj = seg_dir / "final_geometries.trj"
            if not seg_trj.exists():
                raise click.ClickException(
                    f"[all] path-opt segment {idx} did not produce final_geometries.trj"
                )

            try:
                mirror_dir = path_dir / f"{seg_tag}_mep"
                mirror_trj = mirror_dir / "final_geometries.trj"

                _ensure_dir(mirror_dir)
                if seg_trj.resolve() != mirror_trj.resolve():
                    shutil.copy2(seg_trj, mirror_trj)
            except Exception as e:
                click.echo(
                    f"[all] WARNING: failed to mirror path-opt trajectory for segment {idx:02d}: {e}",
                    err=True,
                )

            try:
                seg_mep_trj = path_dir / f"mep_seg_{idx:02d}.trj"
                shutil.copy2(seg_trj, seg_mep_trj)
                if pockets_for_path[0].suffix.lower() == ".pdb":
                    _path_search._maybe_convert_to_pdb(
                        seg_mep_trj,
                        ref_pdb_path=pockets_for_path[0],
                        out_path=path_dir / f"mep_seg_{idx:02d}.pdb",
                    )
            except Exception as e:
                click.echo(
                    f"[all] WARNING: failed to emit per-segment trajectory copies for segment {idx:02d}: {e}",
                    err=True,
                )

            hei_src = seg_dir / "hei.xyz"
            if hei_src.exists():
                try:
                    shutil.copy2(hei_src, path_dir / f"hei_seg_{idx:02d}.xyz")
                    hei_pdb_src = seg_dir / "hei.pdb"
                    if hei_pdb_src.exists():
                        shutil.copy2(hei_pdb_src, path_dir / f"hei_seg_{idx:02d}.pdb")
                    hei_gjf_src = seg_dir / "hei.gjf"
                    if hei_gjf_src.exists():
                        shutil.copy2(hei_gjf_src, path_dir / f"hei_seg_{idx:02d}.gjf")
                except Exception as e:
                    click.echo(
                        f"[all] WARNING: failed to prepare HEI artifacts for segment {idx:02d}: {e}",
                        err=True,
                    )

            blocks = ["\n".join(b) + "\n" for b in _read_xyz_as_blocks(seg_trj)]
            if not blocks:
                raise click.ClickException(
                    f"[all] No frames read from path-opt segment {idx} trajectory: {seg_trj}"
                )
            if idx > 1:
                blocks = blocks[1:]
            combined_blocks.extend(blocks)

            energies_seg: List[float] = []
            for blk in _read_xyz_as_blocks(seg_trj):
                E = np.nan
                if len(blk) >= 2:
                    try:
                        E = float(blk[1].split()[0])
                    except Exception:
                        E = np.nan
                energies_seg.append(E)

            path_opt_segments.append(
                {
                    "tag": seg_tag,
                    "energies": energies_seg,
                    "traj": seg_trj,
                    "inputs": (pL, pR),
                }
            )

        final_trj = path_dir / "mep.trj"
        try:
            final_trj.write_text("".join(combined_blocks), encoding="utf-8")
            click.echo(f"[all] Wrote concatenated MEP trajectory: {final_trj}")
        except Exception as e:
            raise click.ClickException(f"[all] Failed to write concatenated MEP: {e}")

        try:
            run_trj2fig(final_trj, [path_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            _close_matplotlib_figures()
            click.echo(f"[plot] Saved energy plot → '{path_dir / 'mep_plot.png'}'")
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot concatenated MEP: {e}", err=True)

        try:
            if pockets_for_path[0].suffix.lower() == ".pdb":
                mep_pdb = _path_search._maybe_convert_to_pdb(
                    final_trj, ref_pdb_path=pockets_for_path[0], out_path=path_dir / "mep.pdb"
                )
                if mep_pdb and mep_pdb.exists():
                    dst = out_dir / mep_pdb.name
                    shutil.copy2(mep_pdb, dst)
                    click.echo(f"[all] Copied concatenated MEP PDB → {dst}")
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to convert/copy concatenated MEP to PDB: {e}", err=True
            )

        try:
            labels = _build_global_segment_labels(len(path_opt_segments))
            energies_chain: List[float] = []
            for si, seg_info in enumerate(path_opt_segments):
                Es = [float(x) for x in seg_info.get("energies", [])]
                if not Es:
                    continue
                if si == 0:
                    energies_chain.append(Es[0])
                energies_chain.append(float(np.nanmax(Es)))
                energies_chain.append(Es[-1])
            if labels and energies_chain and len(labels) == len(energies_chain):
                title_note = "(GSM; all segments)" if len(path_opt_segments) > 1 else "(GSM)"
                diag_payload = _write_segment_energy_diagram(
                    path_dir / "energy_diagram_mep",
                    labels=labels,
                    energies_au=energies_chain,
                    title_note=title_note,
                )
                if diag_payload:
                    energy_diagrams.append(diag_payload)
        except Exception as e:
            click.echo(f"[diagram] WARNING: Failed to build GSM diagram for path-opt branch: {e}", err=True)

        segments_summary: List[Dict[str, Any]] = []
        bond_cfg = dict(_path_search.BOND_KW)
        for seg_idx, info in enumerate(path_opt_segments, start=1):
            Es = [float(x) for x in info.get("energies", []) if np.isfinite(x)]
            if not Es:
                continue
            barrier = (max(Es) - Es[0]) * AU2KCALPERMOL
            delta = (Es[-1] - Es[0]) * AU2KCALPERMOL
            bond_summary = ""
            try:
                elems, c_first, c_last = _read_xyz_first_last(Path(info["traj"]))
                freeze_atoms = _get_freeze_atoms(info["inputs"][0], freeze_links_flag)
                gL = _geom_from_angstrom(elems, c_first, freeze_atoms)
                gR = _geom_from_angstrom(elems, c_last, freeze_atoms)
                changed, bond_summary = _path_search._has_bond_change(gL, gR, bond_cfg)
                if not changed:
                    bond_summary = "(no covalent changes detected)"
            except Exception as e:
                click.echo(
                    f"[all] WARNING: Failed to detect bond changes for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                bond_summary = "(no covalent changes detected)"

            segments_summary.append(
                {
                    "index": seg_idx,
                    "tag": info.get("tag", f"seg_{seg_idx:03d}"),
                    "kind": "seg",
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": _path_search._bond_changes_block(bond_summary),
                }
            )

        summary = {
            "out_dir": str(path_dir),
            "n_images": len(_read_xyz_as_blocks(final_trj)),
            "n_segments": len(segments_summary),
            "segments": segments_summary,
        }
        if energy_diagrams:
            summary["energy_diagrams"] = list(energy_diagrams)
        try:
            with open(path_dir / "summary.yaml", "w") as f:
                yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
            click.echo(f"[write] Wrote '{path_dir / 'summary.yaml'}'.")
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to write summary.yaml for path-opt branch: {e}", err=True)

        try:
            for name in (
                "mep_plot.png",
                "energy_diagram_MEP.png",
                "mep.pdb",
                "mep_w_ref.pdb",
                "summary.yaml",
                "summary.log",
            ):
                src = path_dir / name
                if src.exists():
                    dst = out_dir / name
                    try:
                        shutil.copy2(src, dst)
                        click.echo(f"[all] Copied {name} → {dst}")
                    except Exception as e:
                        click.echo(
                            f"[all] WARNING: Failed to copy {src} to {dst}: {e}",
                            err=True,
                        )

            for stem in ("mep", "mep_w_ref"):
                for ext in (".trj", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        dst = out_dir / src.name
                        try:
                            shutil.copy2(src, dst)
                            click.echo(f"[all] Copied {src.name} → {dst}")
                        except Exception as e:
                            click.echo(
                                f"[all] WARNING: Failed to copy {src} to {dst}: {e}",
                                err=True,
                            )
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to relocate path-opt summary files: {e}", err=True
            )
        try:
            diag_for_log: Dict[str, Any] = {}
            for diag in summary.get("energy_diagrams", []) or []:
                if isinstance(diag, dict) and str(diag.get("name", "")).lower().endswith("mep"):
                    diag_for_log = diag
                    break
            mep_info = {
                "n_images": summary.get("n_images"),
                "n_segments": summary.get("n_segments"),
                "traj_pdb": str(path_dir / "mep.pdb") if (path_dir / "mep.pdb").exists() else None,
                "mep_plot": str(path_dir / "mep_plot.png") if (path_dir / "mep_plot.png").exists() else None,
                "diagram": diag_for_log,
            }
            summary_payload = {
                "root_out_dir": str(out_dir),
                "path_dir": str(path_dir),
                "path_module_dir": path_dir.name,
                "pipeline_mode": "path-opt",
                "refine_path": refine_path,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": opt_mode.lower() if opt_mode else None,
                "mep_mode": mep_mode_kind,
                "uma_model": calc_cfg_shared.get("model"),
                "command": command_str,
                "charge": q_int,
                "spin": spin,
                "freeze_atoms": _freeze_atoms_for_log(),
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {
                    "summary.yaml": "YAML-format summary",
                    "summary.log": "This summary",
                    "mep_plot.png": "UMA MEP energy vs image index",
                    "energy_diagram_MEP.png": "Compressed MEP diagram R–TS–IM–P",
                },
            }
            write_summary_log(path_dir / "summary.log", summary_payload)
            try:
                shutil.copy2(path_dir / "summary.log", out_dir / "summary.log")
                click.echo(f"[all] Copied summary.log → {out_dir / 'summary.log'}")
            except Exception:
                pass
        except Exception as e:
            click.echo(
                f"[write] WARNING: Failed to write summary.log for path-opt branch: {e}",
                err=True,
            )
    if refine_path:
        # --- recursive GSM path_search branch ---
        click.echo(
            f"\n=== [all] Stage 2/{stage_total} — MEP search on input structures (recursive GSM) ===\n"
        )

        ps_args: List[str] = []

        for p in pockets_for_path:
            ps_args.extend(["-i", str(p)])

        ps_args.extend(["-q", str(q_int)])
        ps_args.extend(["-m", str(int(spin))])

        ps_args.extend(["--freeze-links", "True" if freeze_links_flag else "False"])
        ps_args.extend(["--mep-mode", mep_mode_kind])
        ps_args.extend(["--max-nodes", str(int(max_nodes))])
        ps_args.extend(["--max-cycles", str(int(max_cycles))])
        ps_args.extend(["--climb", "True" if climb else "False"])
        ps_args.extend(["--opt-mode", str(opt_mode.lower())])
        ps_args.extend(["--dump", "True" if dump else "False"])
        if thresh is not None:
            ps_args.extend(["--thresh", str(thresh)])
        ps_args.extend(["--out-dir", str(path_dir)])
        ps_args.extend(["--preopt", "True" if preopt else "False"])
        if args_yaml is not None:
            ps_args.extend(["--args-yaml", str(args_yaml)])

        if gave_ref_pdb:
            for p in (input_paths if not (is_single and has_scan) else (input_paths[:1] * len(pockets_for_path))):
                ps_args.extend(["--ref-pdb", str(p)])

        click.echo("[all] Invoking path_search with arguments:")
        click.echo("  " + " ".join(ps_args))

        _saved_argv = list(sys.argv)
        try:
            sys.argv = ["pdb2reaction", "path_search"] + ps_args
            _path_search.cli.main(args=ps_args, standalone_mode=False)
        except SystemExit as e:
            code = getattr(e, "code", 1)
            if code not in (None, 0):
                raise click.ClickException(
                    f"[all] path_search terminated with exit code {code}."
                )
        except Exception as e:
            raise click.ClickException(f"[all] path_search failed: {e}")
        finally:
            sys.argv = _saved_argv

        try:
            for name in (
                "mep_plot.png",
                "energy_diagram_MEP.png",
                "mep.pdb",
                "mep_w_ref.pdb",
                "summary.yaml",
                "summary.log",
            ):
                src = path_dir / name
                if src.exists():
                    dst = out_dir / name
                    try:
                        shutil.copy2(src, dst)
                        click.echo(f"[all] Copied {name} → {dst}")
                    except Exception as e:
                        click.echo(
                            f"[all] WARNING: Failed to copy {src} to {dst}: {e}",
                            err=True,
                        )

            for stem in ("mep", "mep_w_ref"):
                for ext in (".trj", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        dst = out_dir / src.name
                        try:
                            shutil.copy2(src, dst)
                            click.echo(f"[all] Copied {src.name} → {dst}")
                        except Exception as e:
                            click.echo(
                                f"[all] WARNING: Failed to copy {src} to {dst}: {e}",
                                err=True,
                            )
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to relocate path_search summary files: {e}", err=True
            )

    # -------------------------------------------------------------------------
    # Stage 3: merge to full systems (performed by path_search when enabled)
    # -------------------------------------------------------------------------
    click.echo(f"\n=== [all] Stage 3/{stage_total} — Merge into full-system templates ===\n")
    if refine_path and gave_ref_pdb:
        click.echo(
            "[all] Merging was carried out by path_search using the original inputs as templates."
        )
        click.echo(f"[all] Final products can be found under: {path_dir}")
        click.echo(
            "  - mep_w_ref.pdb            (full-system merged trajectory; also copied to <out-dir>/)"
        )
        click.echo(
            "  - mep_w_ref_seg_XX.pdb     (per-segment merged trajectories for covalent-change segments)"
        )
    elif refine_path:
        click.echo(
            "[all] --ref-pdb was not provided; full-system merged trajectories are not produced."
        )
        click.echo(f"[all] Pocket-only outputs are under: {path_dir}")
    else:
        click.echo(
            "[all] path-opt mode produces pocket-level outputs only; full-system merge is not performed."
        )
        click.echo(f"[all] Aggregated products are under: {path_dir}")
    click.echo("  - summary.yaml             (segment barriers, ΔE, labels)")
    click.echo(
        "  - energy_diagram_MEP.png / energy_diagram.* (copied summary at <out-dir>/)"
    )
    click.echo("\n=== [all] Pipeline finished successfully (core path) ===\n")

    summary_yaml = path_dir / "summary.yaml"
    summary_loaded = load_yaml_dict(summary_yaml) if summary_yaml.exists() else {}
    summary: Dict[str, Any] = summary_loaded if isinstance(summary_loaded, dict) else {}
    segments = _read_summary(summary_yaml)
    if not energy_diagrams:
        existing_diagrams = summary.get("energy_diagrams", [])
        if isinstance(existing_diagrams, list):
            energy_diagrams.extend(existing_diagrams)

    def _write_pipeline_summary_log(post_segment_logs: Sequence[Dict[str, Any]]) -> None:
        try:
            diag_for_log: Dict[str, Any] = {}
            for diag in summary.get("energy_diagrams", []) or []:
                if isinstance(diag, dict) and str(diag.get("name", "")).lower().endswith("mep"):
                    diag_for_log = diag
                    break
            mep_info = {
                "n_images": summary.get("n_images"),
                "n_segments": summary.get("n_segments"),
                "traj_pdb": str(path_dir / "mep.pdb") if (path_dir / "mep.pdb").exists() else None,
                "mep_plot": str(path_dir / "mep_plot.png") if (path_dir / "mep_plot.png").exists() else None,
                "diagram": diag_for_log,
            }
            summary_payload = {
                "root_out_dir": str(out_dir),
                "path_dir": str(path_dir),
                "path_module_dir": path_dir.name,
                "pipeline_mode": "path-search" if refine_path else "path-opt",
                "refine_path": refine_path,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": opt_mode.lower() if opt_mode else None,
                "mep_mode": mep_mode_kind,
                "uma_model": calc_cfg_shared.get("model"),
                "command": command_str,
                "charge": q_int,
                "spin": spin,
                "freeze_atoms": _freeze_atoms_for_log(),
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "post_segments": list(post_segment_logs),
                "key_files": {
                    "summary.yaml": "YAML-format summary",
                    "summary.log": "This summary",
                    "mep_plot.png": "UMA MEP energy vs image index (copied from path_*/)",
                    "energy_diagram_MEP.png": "Compressed MEP diagram R–TS–IM–P (copied from path_*/)",
                },
            }
            write_summary_log(path_dir / "summary.log", summary_payload)
            try:
                shutil.copy2(path_dir / "summary.log", out_dir / "summary.log")
                click.echo(f"[all] Copied summary.log → {out_dir / 'summary.log'}")
            except Exception:
                pass
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

    if not (do_tsopt or do_thermo or do_dft):
        if energy_diagrams:
            summary["energy_diagrams"] = list(energy_diagrams)
        _write_pipeline_summary_log([])
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # -------------------------------------------------------------------------
    # Stage 4: post-processing per reactive segment
    # -------------------------------------------------------------------------
    click.echo(
        "\n=== [all] Stage 4 — Post-processing per reactive segment ===\n"
    )

    if not segments:
        click.echo("[post] No segments found in summary; nothing to do.")
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    reactive = [
        s
        for s in segments
        if (
            s.get("kind", "seg") == "seg"
            and str(s.get("bond_changes", "")).strip()
            and str(s.get("bond_changes", "")).strip()
            != "(no covalent changes detected)"
        )
    ]
    if not reactive:
        click.echo("[post] No bond-change segments. Skipping TS/thermo/DFT.")
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # Per-category per-segment energies
    tsopt_seg_energies: List[Tuple[float, float, float]] = []
    g_uma_seg_energies: List[Tuple[float, float, float]] = []
    dft_seg_energies: List[Tuple[float, float, float]] = []
    g_dftuma_seg_energies: List[Tuple[float, float, float]] = []
    irc_trj_for_all: List[Tuple[Path, bool]] = []
    post_segment_logs: List[Dict[str, Any]] = []

    for s in reactive:
        seg_idx = int(s.get("index", 0) or 0)
        seg_tag = s.get("tag", f"seg_{seg_idx:02d}")
        click.echo(f"\n--- [post] Segment {seg_idx:02d} ({seg_tag}) ---")

        seg_root = path_dir
        seg_dir = seg_root / f"post_seg_{seg_idx:02d}"
        _ensure_dir(seg_dir)

        segment_log = {
            "index": seg_idx,
            "tag": seg_tag,
            "kind": s.get("kind", "seg"),
            "bond_changes": s.get("bond_changes", ""),
            "mep_barrier_kcal": s.get("barrier_kcal"),
            "mep_delta_kcal": s.get("delta_kcal"),
            "post_dir": str(seg_dir),
        }
        post_segment_logs.append(segment_log)

        hei_base = seg_root / f"hei_seg_{seg_idx:02d}"
        hei_pocket_path = _find_with_suffixes(hei_base, [".pdb", ".xyz", ".gjf"])
        if hei_pocket_path is None:
            click.echo(
                f"[post] WARNING: HEI pocket file not found for segment {seg_idx:02d} (searched .pdb/.xyz/.gjf); skipping TSOPT.",
                err=True,
            )
            continue

        struct_dir = seg_dir / "structures"
        _ensure_dir(struct_dir)
        state_structs: Dict[str, Path] = {}
        uma_ref_energies: Dict[str, float] = {}

        if do_tsopt:
            ts_pdb, g_ts = _run_tsopt_on_hei(
                hei_pocket_path,
                q_int,
                spin,
                calc_cfg_shared,
                args_yaml,
                seg_dir,
                freeze_links_flag,
                tsopt_opt_mode_default,
                overrides=tsopt_overrides,
            )

            irc_res = _irc_and_match(
                seg_idx=seg_idx,
                seg_dir=seg_dir,
                ref_pdb_for_seg=ts_pdb,
                seg_pocket_pdb=hei_pocket_path,
                g_ts=g_ts,
                q_int=q_int,
                spin=spin,
                freeze_links_flag=freeze_links_flag,
                calc_cfg=calc_cfg_shared,
                args_yaml=args_yaml,
                seg_tag=str(seg_tag),
            )

            gL = irc_res["left_min_geom"]
            gR = irc_res["right_min_geom"]
            gT = irc_res["ts_geom"]
            irc_plot_path = irc_res.get("irc_plot_path")
            irc_trj_path = irc_res.get("irc_trj_path")
            reverse_irc = bool(irc_res.get("reverse_irc", False))

            if isinstance(irc_plot_path, Path) and irc_plot_path.exists():
                segment_log["irc_plot"] = str(irc_plot_path)
            if isinstance(irc_trj_path, Path) and irc_trj_path.exists():
                segment_log["irc_traj"] = str(irc_trj_path)

            if isinstance(irc_trj_path, Path) and irc_trj_path.exists():
                irc_trj_for_all.append((irc_trj_path, reverse_irc))

            pL_irc = _save_single_geom_as_pdb_for_tools(
                gL, hei_pocket_path, struct_dir, "reactant_irc"
            )
            pT = _save_single_geom_as_pdb_for_tools(
                gT, hei_pocket_path, struct_dir, "ts"
            )
            pR_irc = _save_single_geom_as_pdb_for_tools(
                gR, hei_pocket_path, struct_dir, "product_irc"
            )

            endpoint_opt_dir = seg_dir / "endpoint_opt"
            _ensure_dir(endpoint_opt_dir)
            try:
                g_react_opt, _ = _optimize_endpoint_geom(
                    gL,
                    tsopt_opt_mode_default,
                    endpoint_opt_dir,
                    f"seg_{seg_idx:02d}_reactant",
                    dump=dump,
                )
            except Exception as e:
                click.echo(
                    f"[post] WARNING: Reactant endpoint optimization failed for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                g_react_opt = gL
            try:
                g_prod_opt, _ = _optimize_endpoint_geom(
                    gR,
                    tsopt_opt_mode_default,
                    endpoint_opt_dir,
                    f"seg_{seg_idx:02d}_product",
                    dump=dump,
                )
            except Exception as e:
                click.echo(
                    f"[post] WARNING: Product endpoint optimization failed for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                g_prod_opt = gR

            shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
            click.echo(f"[endpoint-opt] Clean endpoint-opt working dir.") 

            pL = _save_single_geom_as_pdb_for_tools(
                g_react_opt, hei_pocket_path, struct_dir, "reactant"
            )
            pR = _save_single_geom_as_pdb_for_tools(
                g_prod_opt, hei_pocket_path, struct_dir, "product"
            )
            state_structs = {"R": pL, "TS": pT, "P": pR}

            eR = float(g_react_opt.energy)
            eT = float(gT.energy)
            eP = float(g_prod_opt.energy)
            uma_ref_energies = {"R": eR, "TS": eT, "P": eP}
            diag_payload = _write_segment_energy_diagram(
                seg_dir / "energy_diagram_UMA",
                labels=["R", f"TS{seg_idx}", "P"],
                energies_au=[eR, eT, eP],
                title_note="(UMA, TSOPT + IRC)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

            tsopt_seg_energies.append((eR, eT, eP))

            segment_log["uma"] = {
                "labels": ["R", f"TS{seg_idx}", "P"],
                "energies_au": [eR, eT, eP],
                "energies_kcal": [0.0, (eT - eR) * AU2KCALPERMOL, (eP - eR) * AU2KCALPERMOL],
                "diagram": str(seg_dir / "energy_diagram_UMA.png"),
                "structures": state_structs,
                "barrier_kcal": (eT - eR) * AU2KCALPERMOL,
                "delta_kcal": (eP - eR) * AU2KCALPERMOL,
            }

            try:
                if pT.suffix.lower() == ".pdb":
                    ts_copy = out_dir / f"ts_seg_{seg_idx:02d}.pdb"
                    shutil.copy2(pT, ts_copy)
                else:
                    ts_copy = out_dir / f"ts_seg_{seg_idx:02d}.xyz"
                    _path_search._write_xyz_trj_with_energy(
                        [gT], [float(gT.energy)], ts_copy
                    )
                click.echo(
                    f"[all] Copied TS structure for segment {seg_idx:02d} → {ts_copy}"
                )
            except Exception as e:
                click.echo(
                    f"[all] WARNING: Failed to write TS structure for segment {seg_idx:02d}: {e}",
                    err=True,
                )

        elif do_thermo or do_dft:
            seg_pocket_path = _find_with_suffixes(
                seg_root / f"mep_seg_{seg_idx:02d}", [".pdb"]
            )

            # Decide reference PDB (if any) for freeze-atoms detection / PDB conversion
            freeze_ref: Optional[Path] = None
            if seg_pocket_path is not None and seg_pocket_path.suffix.lower() == ".pdb":
                freeze_ref = seg_pocket_path
            elif hei_pocket_path.suffix.lower() == ".pdb":
                freeze_ref = hei_pocket_path

            freeze_atoms: List[int] = _get_freeze_atoms(freeze_ref, freeze_links_flag)

            try:
                endpoints = _load_segment_endpoints(seg_root, str(seg_tag), freeze_atoms)
                if endpoints is None:
                    click.echo(
                        f"[post] WARNING: final_geometries.trj not found for segment {seg_idx:02d}; cannot run thermo/DFT without --tsopt. Skipping segment.",
                        err=True,
                    )
                    continue
                gL, gR = endpoints
            except Exception as e:
                click.echo(
                    f"[post] WARNING: failed to load segment endpoints from final_geometries.trj for segment {seg_idx:02d}: {e}. Skipping segment.",
                    err=True,
                )
                continue

            try:
                g_ts = geom_loader(
                    hei_pocket_path,
                    coord_type=DEFAULT_COORD_TYPE,
                    freeze_atoms=freeze_atoms,
                )
                if freeze_atoms:
                    fa = np.array(freeze_atoms, dtype=int)
                    g_ts.freeze_atoms = fa
            except Exception as e:
                click.echo(
                    f"[post] WARNING: failed to load HEI geometry for segment {seg_idx:02d}: {e}. Skipping segment.",
                    err=True,
                )
                continue

            calc_args = dict(calc_cfg_shared)
            calc = uma_pysis(**calc_args)
            gL.set_calculator(calc)
            gR.set_calculator(calc)
            g_ts.set_calculator(calc)

            ref_for_structs = seg_pocket_path if seg_pocket_path is not None else hei_pocket_path
            pL = _save_single_geom_as_pdb_for_tools(
                gL, ref_for_structs, struct_dir, "reactant_mep"
            )
            pR = _save_single_geom_as_pdb_for_tools(
                gR, ref_for_structs, struct_dir, "product_mep"
            )
            pT = _save_single_geom_as_pdb_for_tools(
                g_ts, ref_for_structs, struct_dir, "ts_from_hei"
            )
            state_structs = {"R": pL, "TS": pT, "P": pR}

        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        freq_seg_root = _resolve_override_dir(seg_dir / "freq", freq_out_dir)
        dft_seg_root = _resolve_override_dir(seg_dir / "dft", dft_out_dir)

        if (do_thermo or do_dft) and not state_structs:
            click.echo(
                f"[post] WARNING: No segment structures prepared for segment {seg_idx:02d}; skipping thermo/DFT.",
                err=True,
            )
            continue

        p_react = state_structs.get("R")
        p_ts = state_structs.get("TS")
        p_prod = state_structs.get("P")

        if do_thermo:
            if not (p_react and p_ts and p_prod):
                click.echo(
                    f"[thermo] WARNING: Missing R/TS/P structures for segment {seg_idx:02d}; skipping thermo.",
                    err=True,
                )
            else:
                click.echo(
                    f"[thermo] Segment {seg_idx:02d}: freq on R/TS/P"
                )
                tR = _run_freq_for_state(
                    p_react,
                    q_int,
                    spin,
                    freq_seg_root / "R",
                    args_yaml,
                    freeze_links_flag,
                    overrides=freq_overrides,
                )
                tT = _run_freq_for_state(
                    p_ts,
                    q_int,
                    spin,
                    freq_seg_root / "TS",
                    args_yaml,
                    freeze_links_flag,
                    overrides=freq_overrides,
                )
                tP = _run_freq_for_state(
                    p_prod,
                    q_int,
                    spin,
                    freq_seg_root / "P",
                    args_yaml,
                    freeze_links_flag,
                    overrides=freq_overrides,
                )
                thermo_payloads = {"R": tR, "TS": tT, "P": tP}
                ts_freq_info = _read_imaginary_frequency_info(freq_seg_root / "TS")
                if ts_freq_info is not None:
                    segment_log["ts_imag"] = ts_freq_info
                    if ts_freq_info.get("nu_imag_max_cm") is not None:
                        segment_log["ts_imag_freq_cm"] = ts_freq_info["nu_imag_max_cm"]
                try:
                    GR = float(
                        tR.get(
                            "sum_EE_and_thermal_free_energy_ha",
                            uma_ref_energies.get("R", np.nan),
                        )
                    )
                    GT = float(
                        tT.get(
                            "sum_EE_and_thermal_free_energy_ha",
                            uma_ref_energies.get("TS", np.nan),
                        )
                    )
                    GP = float(
                        tP.get(
                            "sum_EE_and_thermal_free_energy_ha",
                            uma_ref_energies.get("P", np.nan),
                        )
                    )
                    gibbs_vals = [GR, GT, GP]
                    if all(np.isfinite(gibbs_vals)):
                        diag_payload = _write_segment_energy_diagram(
                            seg_dir / "energy_diagram_G_UMA",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_au=gibbs_vals,
                            title_note="(UMA + Thermal Correction)",
                            ylabel="ΔG (kcal/mol)",
                        )
                        if diag_payload:
                            energy_diagrams.append(diag_payload)
                        g_uma_seg_energies.append((GR, GT, GP))
                        segment_log["gibbs_uma"] = {
                            "labels": ["R", f"TS{seg_idx}", "P"],
                            "energies_au": gibbs_vals,
                            "energies_kcal": [
                                (GR - GR) * AU2KCALPERMOL,
                                (GT - GR) * AU2KCALPERMOL,
                                (GP - GR) * AU2KCALPERMOL,
                            ],
                            "diagram": str(seg_dir / "energy_diagram_G_UMA.png"),
                            "structures": state_structs,
                            "barrier_kcal": (GT - GR) * AU2KCALPERMOL,
                            "delta_kcal": (GP - GR) * AU2KCALPERMOL,
                        }
                    else:
                        click.echo(
                            "[thermo] NOTE: Gibbs energies non-finite; diagram skipped."
                        )
                except Exception as e:
                    click.echo(
                        f"[thermo] WARNING: failed to build Gibbs diagram: {e}",
                        err=True,
                    )

        if do_dft:
            if not (p_react and p_ts and p_prod):
                click.echo(
                    f"[dft] WARNING: Missing R/TS/P structures for segment {seg_idx:02d}; skipping DFT.",
                    err=True,
                )
            else:
                click.echo(f"[dft] Segment {seg_idx:02d}: DFT on R/TS/P")
                dft_jobs = [
                    ("R", p_react, dft_seg_root / "R"),
                    ("TS", p_ts, dft_seg_root / "TS"),
                    ("P", p_prod, dft_seg_root / "P"),
                ]
                dft_payloads = _run_dft_sequence(
                    dft_jobs,
                    q_int,
                    spin,
                    args_yaml,
                    dft_func_basis_use,
                    dft_overrides,
                    dft_engine,
                )
                dR = dft_payloads.get("R")
                dT = dft_payloads.get("TS")
                dP = dft_payloads.get("P")
                try:
                    eR_dft = float(
                        ((dR or {}).get("energy", {}) or {}).get(
                            "hartree", uma_ref_energies.get("R", np.nan)
                        )
                    )
                    eT_dft = float(
                        ((dT or {}).get("energy", {}) or {}).get(
                            "hartree", uma_ref_energies.get("TS", np.nan)
                        )
                    )
                    eP_dft = float(
                        ((dP or {}).get("energy", {}) or {}).get(
                            "hartree", uma_ref_energies.get("P", np.nan)
                        )
                    )
                    if all(map(np.isfinite, [eR_dft, eT_dft, eP_dft])):
                        diag_payload = _write_segment_energy_diagram(
                            seg_dir / "energy_diagram_DFT",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_au=[eR_dft, eT_dft, eP_dft],
                            title_note=f"({dft_func_basis_use})",
                        )
                        if diag_payload:
                            energy_diagrams.append(diag_payload)
                        dft_seg_energies.append((eR_dft, eT_dft, eP_dft))
                        segment_log["dft"] = {
                            "labels": ["R", f"TS{seg_idx}", "P"],
                            "energies_au": [eR_dft, eT_dft, eP_dft],
                            "energies_kcal": [
                                0.0,
                                (eT_dft - eR_dft) * AU2KCALPERMOL,
                                (eP_dft - eR_dft) * AU2KCALPERMOL,
                            ],
                            "diagram": str(seg_dir / "energy_diagram_DFT.png"),
                            "structures": state_structs,
                            "barrier_kcal": (eT_dft - eR_dft) * AU2KCALPERMOL,
                            "delta_kcal": (eP_dft - eR_dft) * AU2KCALPERMOL,
                        }
                    else:
                        click.echo(
                            "[dft] WARNING: some DFT energies missing; diagram skipped.",
                            err=True,
                        )
                except Exception as e:
                    click.echo(
                        f"[dft] WARNING: failed to build DFT diagram: {e}", err=True
                    )

                if do_thermo:
                    try:
                        dG_R = float(
                            (thermo_payloads.get("R", {}) or {}).get(
                                "thermal_correction_free_energy_ha", 0.0
                            )
                        )
                        dG_T = float(
                            (thermo_payloads.get("TS", {}) or {}).get(
                                "thermal_correction_free_energy_ha", 0.0
                            )
                        )
                        dG_P = float(
                            (thermo_payloads.get("P", {}) or {}).get(
                                "thermal_correction_free_energy_ha", 0.0
                            )
                        )
                        GR_dftUMA = eR_dft + dG_R
                        GT_dftUMA = eT_dft + dG_T
                        GP_dftUMA = eP_dft + dG_P
                        if all(
                            np.isfinite([GR_dftUMA, GT_dftUMA, GP_dftUMA])
                        ):
                            diag_payload = _write_segment_energy_diagram(
                                seg_dir / "energy_diagram_G_DFT_plus_UMA",
                                labels=["R", f"TS{seg_idx}", "P"],
                                energies_au=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                                title_note=f"({dft_func_basis_use} // UMA  + Thermal Correction)",
                                ylabel="ΔG (kcal/mol)",
                            )
                            if diag_payload:
                                energy_diagrams.append(diag_payload)
                            g_dftuma_seg_energies.append(
                                (GR_dftUMA, GT_dftUMA, GP_dftUMA)
                            )
                            segment_log["gibbs_dft_uma"] = {
                                "labels": ["R", f"TS{seg_idx}", "P"],
                                "energies_au": [GR_dftUMA, GT_dftUMA, GP_dftUMA],
                                "energies_kcal": [
                                    (GR_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                                    (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                                    (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                                ],
                                "diagram": str(seg_dir / "energy_diagram_G_DFT_plus_UMA.png"),
                                "structures": state_structs,
                                "barrier_kcal": (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                                "delta_kcal": (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                            }
                        else:
                            click.echo(
                                "[dft//uma] WARNING: DFT//UMA Gibbs energies non-finite; diagram skipped.",
                                err=True,
                            )
                    except Exception as e:
                        click.echo(
                            f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}",
                            err=True,
                        )

    # -------------------------------------------------------------------------
    # Aggregate diagrams over all reactive segments, with GSM-style labels
    # -------------------------------------------------------------------------
    if tsopt_seg_energies:
        tsopt_all_energies = [e for triple in tsopt_seg_energies for e in triple]
        tsopt_all_labels = _build_global_segment_labels(len(tsopt_seg_energies))
        if tsopt_all_labels and len(tsopt_all_labels) == len(tsopt_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_UMA_all",
                labels=tsopt_all_labels,
                energies_au=tsopt_all_energies,
                title_note="(UMA, TSOPT + IRC; all segments)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_thermo and g_uma_seg_energies:
        g_uma_all_energies = [e for triple in g_uma_seg_energies for e in triple]
        g_uma_all_labels = _build_global_segment_labels(len(g_uma_seg_energies))
        if g_uma_all_labels and len(g_uma_all_labels) == len(g_uma_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_G_UMA_all",
                labels=g_uma_all_labels,
                energies_au=g_uma_all_energies,
                title_note="(UMA + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_dft and dft_seg_energies:
        dft_all_energies = [e for triple in dft_seg_energies for e in triple]
        dft_all_labels = _build_global_segment_labels(len(dft_seg_energies))
        if dft_all_labels and len(dft_all_labels) == len(dft_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_DFT_all",
                labels=dft_all_labels,
                energies_au=dft_all_energies,
                title_note=f"({dft_func_basis_use}; all segments)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_dft and do_thermo and g_dftuma_seg_energies:
        g_dftuma_all_energies = [e for triple in g_dftuma_seg_energies for e in triple]
        g_dftuma_all_labels = _build_global_segment_labels(len(g_dftuma_seg_energies))
        if g_dftuma_all_labels and len(g_dftuma_all_labels) == len(g_dftuma_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_G_DFT_plus_UMA_all",
                labels=g_dftuma_all_labels,
                energies_au=g_dftuma_all_energies,
                title_note=f"({dft_func_basis_use} // UMA  + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    # -------------------------------------------------------------------------
    # Aggregated IRC plot over all reactive segments (single trj + trj2fig)
    # -------------------------------------------------------------------------
    if irc_trj_for_all:
        _merge_irc_trajectories_to_single_plot(
            irc_trj_for_all, out_dir / "irc_plot_all.png"
        )

    # Refresh summary.yaml with final energy diagram metadata (including aggregated diagrams)
    try:
        summary["energy_diagrams"] = list(energy_diagrams)
        with open(path_dir / "summary.yaml", "w") as f:
            yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
        click.echo(f"[write] Updated '{path_dir / 'summary.yaml'}' with energy diagrams.")
        try:
            dst_summary = out_dir / "summary.yaml"
            shutil.copy2(path_dir / "summary.yaml", dst_summary)
            click.echo(f"[all] Copied summary.yaml → {dst_summary}")
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to mirror summary.yaml to {out_dir}: {e}",
                err=True,
            )
        _write_pipeline_summary_log(post_segment_logs)
    except Exception as e:
        click.echo(
            f"[write] WARNING: Failed to refresh summary.yaml with energy diagram metadata: {e}",
            err=True,
        )

    click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
