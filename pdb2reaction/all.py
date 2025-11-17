# pdb2reaction/all.py

"""
all — SINGLE command to execute an end-to-end enzymatic reaction workflow:
Extract pockets → (optional) staged scan on a single structure → MEP (recursive GSM) → merge to full systems,
with optional TS optimization, pseudo‑IRC, thermochemistry, DFT, and DFT//UMA diagrams
================================================================================================================

Usage (CLI)
-----
    # Standard (multi-structure ensemble in reaction order)
    pdb2reaction all -i R.pdb [I1.pdb ...] P.pdb -c <substrate-spec> [--ligand-charge <map-or-number>]
                      [--spin <2S+1>] [--freeze-links True|False] [--max-nodes N] [--max-cycles N]
                      [--climb True|False] [--sopt-mode lbfgs|rfo|light|heavy]
                      [--opt-mode light|lbfgs|heavy|rfo]
                      [--dump True|False] [--args-yaml params.yaml]
                      [--preopt True|False] [--hessian-calc-mode Analytical|FiniteDifference] [--out-dir DIR]
                      [--tsopt True|False] [--thermo True|False] [--dft True|False]
                      [--tsopt-max-cycles N] [--freq-* overrides] [--dft-* overrides]

    # Override examples (repeatable; use only what you need)
      ... --scan-lists "[(12,45,1.35)]" --scan-one-based True --scan-bias-k 0.05 --scan-relax-max-cycles 150 \
          --tsopt-max-cycles 250 --freq-temperature 298.15 --freq-max-write 15 --dft-func-basis "wb97x-v/def2-tzvp"

    # Single-structure + staged scan (the scan creates intermediate/product candidates after extraction)
    pdb2reaction all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" [--scan-lists "..."] \
                      --spin 1 --freeze-links True --sopt-mode lbfgs --preopt True \
                      --out-dir result_all --tsopt True --thermo True --dft True

    # Single-structure TSOPT-only mode (no path search):
    # if exactly one input is given, --scan-lists is NOT provided, and --tsopt True:
    # run TS optimization on the extracted pocket, do pseudo‑IRC, minimize both ends, and make diagrams.
    pdb2reaction all -i single.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
                      --tsopt True --thermo True --dft True --out-dir result_tsopt_single


Examples
-----
    # Minimal end-to-end run with explicit substrate and ligand charges (multi-structure)
    pdb2reaction all -i reactant.pdb product.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1"

    # Full ensemble with an intermediate, residue-ID substrate spec, and full post-processing
    pdb2reaction all -i A.pdb B.pdb C.pdb -c "308,309" --ligand-charge "-1" \
      --spin 1 --freeze-links True --max-nodes 10 --max-cycles 100 --climb True \
      --sopt-mode lbfgs --dump False --args-yaml params.yaml --preopt True \
      --out-dir result_all --tsopt True --thermo True --dft True

    # Single-structure + scan to build an ordered series (initial + stage results) → path_search + post
    pdb2reaction all -i A.pdb -c "308,309" --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
      --spin 1 --out-dir result_scan_all --tsopt True --thermo True --dft True

    # Single-structure TSOPT-only mode (no path_search)
    pdb2reaction all -i A.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
      --tsopt True --thermo True --dft True --out-dir result_tsopt_only


Description
-----
Runs a one-shot pipeline centered on pocket models:

(1) **Active-site pocket extraction** (multi-structure union when multiple inputs)
    - Define the substrate (`-c/--center`, by PDB, residue IDs, or residue names).
    - Optionally provide `--ligand-charge` as a total number (distributed) or a mapping (e.g., `GPP:-3,MMT:-1`).
    - The extractor writes per-input pocket PDBs under `<out-dir>/pockets/`.
    - The extractor’s **first-model total pocket charge** is used as the total charge in later steps,
      cast to the nearest integer with a console note if rounding occurs.
    - Additional extractor toggles: `--radius`, `--radius-het2het`, `--include-H2O True|False`,
      `--exclude-backbone True|False`, `--add-linkH True|False`, `--selected_resn`, `--verbose True|False`.

(1b) **Optional staged scan (single-structure only)**
    - If **exactly one** full input PDB is provided and `--scan-lists` is given, the tool performs a
      **staged, bond-length–driven scan** *on the extracted pocket PDB* using the UMA calculator.
    - For each stage, the final relaxed structure (`stage_XX/result.pdb`) is collected as an
      **intermediate/product candidate**.
    - The ordered input series for the path search becomes:
      `[initial pocket, stage_01/result.pdb, stage_02/result.pdb, ...]`.

(2) **MEP search (recursive GSM) on pocket inputs**
    - Runs `path_search` with options forwarded from this command.
    - For multi-input runs, the original **full** PDBs are supplied as **merge references** automatically.
      In the scan-derived series (single-structure case), the single original full PDB is reused (repeated)
      as the reference template for all pocket inputs.

(3) **Merge to full systems**
    - The pocket MEP is merged back into the original full-system template(s) within `<out-dir>/path_search/`.
    - Pocket-only and full-system trajectories, per-segment merged PDBs, and a summary are written.

(4) **Optional per-segment post-processing** (for segments with covalent changes)
    - `--tsopt True`: Optimize TS on the HEI pocket; perform a pseudo‑IRC by displacing along the
      imaginary mode; assign forward/backward correspondence by bond-state matching; render a segment diagram.
    - `--thermo True`: Compute UMA thermochemistry on (R, TS, P) and add a Gibbs diagram.
    - `--dft True`: Do DFT single-point on (R, TS, P) and add a DFT diagram.
      With `--thermo True`, also generate a **DFT//UMA** Gibbs diagram.

(Alt) **Single-structure TSOPT-only mode**
    - If **exactly one** input is given, **no** `--scan-lists` is provided, and `--tsopt True`,
      the tool skips (2)-(3) and:
        • Runs `tsopt` on the **pocket** of that structure,  
        • Does a pseudo‑IRC and minimizes both ends,  
        • Builds UMA energy diagrams for **R–TS–P**,  
        • Optionally adds UMA Gibbs, DFT, and **DFT//UMA** diagrams.  
      ※ このモードに**限り**、IRCで得た両端のうち**エネルギーが高い方を反応物 (R)** として採用します。

**Charge handling**
  - The extractor’s **first-model total pocket charge** is used as the path/scan/TSOPT total charge (rounded to int).

**Inputs**
  - `-i/--input` accepts two or more **full** PDBs in reaction order (reactant [intermediates ...] product), or
    **a single** full PDB (with `--scan-lists` *or* `--tsopt True`).

**Forwarded / relevant options**
  - Path search: `--spin`, `--freeze-links`, `--max-nodes`, `--max-cycles`, `--climb`, `--sopt-mode`,
    `--dump`, `--preopt`, `--args-yaml`, `--out-dir`. (`--freeze-links` / `--dump` propagate to scan/tsopt/freq as shared flags.)
  - Scan (single-structure, pocket input): inherits charge/spin, `--freeze-links`, `--opt-mode`, `--dump`,
    `--args-yaml`, `--preopt`, and per-stage overrides (`--scan-out-dir`, `--scan-one-based` True|False
    (omit to keep scan's default 1-based indexing),
    `--scan-max-step-size`, `--scan-bias-k`, `--scan-relax-max-cycles`, `--scan-preopt`, `--scan-endopt`).
  - Shared knobs: `--opt-mode light|lbfgs|heavy|rfo` applies to both scan and tsopt; when omitted, scan defaults to
    LBFGS or RFO based on `--sopt-mode`, and tsopt falls back to `light`. `--hessian-calc-mode` applies to tsopt and freq.
  - TS optimization / pseudo-IRC: `--tsopt-max-cycles`, `--tsopt-out-dir`, and the shared knobs above tune downstream tsopt.
  - Frequency analysis: `--freq-out-dir`, `--freq-max-write`, `--freq-amplitude-ang`, `--freq-n-frames`, `--freq-sort`,
    `--freq-temperature`, `--freq-pressure`, plus shared `--freeze-links`, `--dump`, `--hessian-calc-mode`.
  - DFT single-points: `--dft-out-dir`, `--dft-func-basis`, `--dft-max-cycle`, `--dft-conv-tol`, `--dft-grid-level`.
  - Post-processing toggles: `--tsopt`, `--thermo`, `--dft`.
  - YAML forwarding: `--args-yaml` is passed unchanged to `path_search`, `scan`, `tsopt`, `freq`, and `dft` so a single file can
    host per-module sections (see the respective subcommand docs for accepted keys).

Outputs (& Directory Layout)
-----
<out-dir>/
  pockets/
    pocket_<input1_basename>.pdb
    pocket_<input2_basename>.pdb
    ...
  scan/                                  # present only in single-structure+scan mode
    stage_01/result.pdb
    stage_02/result.pdb
    ...
  path_search/                           # present when path_search is executed
    mep.trj
    energy.png
    mep.pdb
    mep_w_ref.pdb
    mep_w_ref_seg_XX.pdb
    summary.yaml
    tsopt_seg_XX/                        # when post-processing is enabled
      ts/ ...
      irc/ ...
      freq/ ...                          # with --thermo True
      dft/  ...                          # with --dft True
      energy_diagram_tsopt.(html|png)
      energy_diagram_G_UMA.(html|png)
      energy_diagram_DFT.(html|png)
      energy_diagram_G_DFT_plus_UMA.(html|png)
  tsopt_single/                          # present only in single-structure TSOPT-only mode
    ts/ ...
    irc/ ...
    structures/
      reactant.pdb
      ts.pdb
      product.pdb
    freq/ ...                            # with --thermo True
    dft/  ...                            # with --dft True
    energy_diagram_tsopt.(html|png)
    energy_diagram_G_UMA.(html|png)
    energy_diagram_DFT.(html|png)
    energy_diagram_G_DFT_plus_UMA.(html|png)

Notes
-----
- **Python ≥ 3.10** is required.
- **Substrate (`-c/--center`) と（必要に応じて）`--ligand-charge` の指定が実用上ほぼ必須**です。
- **Single-structure** で `--scan-lists` **または** `--tsopt True` のいずれかが与えられていれば実行可能。
  それ以外は従来通り、**少なくとも2構造**の入力が必要です。
- Energies in diagrams are plotted relative to the first state in kcal/mol (converted from Hartree).
- **NEW**: `-c/--center` を省略した場合は抽出をスキップし、**入力の全体構造**をそのままポケット入力として用います。
  この場合の全体電荷は既定で 0 とし、`--ligand-charge <number>` を与えたときはその数値を**全体系の電荷**として用います（四捨五入で整数化）。
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import List, Sequence, Optional, Tuple, Dict, Any

import sys
import math
import click
from click.core import ParameterSource
import time  # timing
import yaml
import numpy as np
import torch

# Biopython for PDB parsing (post-processing helpers)
from Bio import PDB

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL

# Local imports from the package
from .extract import extract_api
from . import path_search as _path_search
from . import tsopt as _tsopt
from . import freq as _freq_cli
from . import dft as _dft_cli
from .uma_pysis import uma_pysis
from .trj2fig import run_trj2fig
from .utils import (
    build_energy_diagram,
    detect_freeze_links_safe,
    format_elapsed,
    prepare_input_structure,
    maybe_convert_xyz_to_gjf,
    # ↓↓↓ 追加（GJF から電荷・多重度を拾うため）
    resolve_charge_spin_or_raise,
)
from . import scan as _scan_cli
from .add_elem_info import assign_elements as _assign_elem_info

# -----------------------------
# Helpers
# -----------------------------

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
    """
    Robustly collect values following a flag that may appear **once** followed by multiple space-separated values,
    e.g., "-i A B C". This mirrors the behavior implemented in `path_search.cli`.
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


# [FIX] --- begin: robust atom key helpers (serial 依存の廃止) ---------------------------

def _parse_atom_key_from_line(line: str) -> Optional[Tuple[str, str, str, str, str, str]]:
    """
    Extract a structural identity key from a PDB ATOM/HETATM record.

    Returns:
        (chainID, resName, resSeq, iCode, atomName, altLoc), with blanks normalized to ''.
    """
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    atomname = line[12:16].strip()
    altloc   = (line[16] if len(line) > 16 else " ").strip()
    resname  = line[17:20].strip()
    chain    = (line[21] if len(line) > 21 else " ").strip()
    resseq   = line[22:26].strip()
    icode    = (line[26] if len(line) > 26 else " ").strip()
    return (chain, resname, resseq, icode, atomname, altloc)


def _pocket_key_to_index(pocket_pdb: Path) -> Dict[Tuple[str, str, str, str, str, str], int]:
    """
    Build mapping: structural atom key -> pocket index (1-based by file order).
    """
    key2idx: Dict[Tuple[str, str, str, str, str, str], int] = {}
    idx = 0
    try:
        with open(pocket_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    key = _parse_atom_key_from_line(line)
                    if key is None:
                        continue
                    idx += 1
                    key2idx[key] = idx
    except FileNotFoundError:
        raise click.ClickException(f"[all] Pocket PDB not found: {pocket_pdb}")
    if not key2idx:
        raise click.ClickException(f"[all] Pocket PDB {pocket_pdb} has no ATOM/HETATM records.")
    return key2idx


def _read_full_atom_keys_in_file_order(full_pdb: Path) -> List[Tuple[str, str, str, str, str, str]]:
    """
    Read ATOM/HETATM lines and return keys in the original file order.
    """
    keys: List[Tuple[str, str, str, str, str, str]] = []
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


def _format_atom_key_for_msg(key: Tuple[str, str, str, str, str, str]) -> str:
    """Pretty string for diagnostics."""
    chain, resn, resseq, icode, atom, alt = key
    res = f"{chain}:{resn}{resseq}{(icode if icode else '')}"
    alt_sfx = f",alt={alt}" if alt else ""
    return f"{res}:{atom}{alt_sfx}"

# [FIX] --- end: robust atom key helpers (serial 依存の廃止) -----------------------------


def _parse_scan_lists_literals(scan_lists_raw: Sequence[str]) -> List[List[Tuple[int, int, float]]]:
    """Parse ``--scan-lists`` literals without re-basing atom indices."""
    stages: List[List[Tuple[int, int, float]]] = []
    for idx_stage, literal in enumerate(scan_lists_raw, start=1):
        try:
            obj = ast.literal_eval(literal)
        except Exception as exc:  # pragma: no cover - defensive
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

    [FIX] serial には依存せず、(chainID, resName, resSeq, iCode, atomName, altLoc) の「構造的キー」で対応付ける。
    """
    if not scan_lists_raw:
        return []

    stages = _parse_scan_lists_literals(scan_lists_raw)

    # [FIX] 元構造のキー配列（1-based indexと整合）と、ポケット側のキー→index辞書を用意
    orig_keys_in_order = _read_full_atom_keys_in_file_order(full_input_pdb)
    key_to_pocket_idx  = _pocket_key_to_index(pocket_pdb)

    n_atoms_full = len(orig_keys_in_order)

    def _map_full_index_to_pocket(idx_one_based: int, stage_idx: int, tuple_idx: int, side_label: str) -> int:
        """
        1-based のフル PDB index を、ポケット PDB の 1-based index に変換。
        厳密一致→(altloc無視)→(icode無視)→(両方無視)の順でフォールバック。
        """
        key = orig_keys_in_order[idx_one_based - 1]

        # 候補キー列を順に試す
        chain, resn, resseq, icode, atom, alt = key
        candidates = [
            key,
            (chain, resn, resseq, icode, atom, ''),   # altloc 無視
            (chain, resn, resseq, '',    atom, alt), # iCode 無視
            (chain, resn, resseq, '',    atom, ''),  # 両方無視
        ]
        for k in candidates:
            hit = key_to_pocket_idx.get(k, None)
            if hit is not None:
                return hit

        # 見つからない → 抽出で除外された等。親切なメッセージを出す。
        msg_key = _format_atom_key_for_msg(key)
        raise click.BadParameter(
            f"--scan-lists #{stage_idx} tuple #{tuple_idx} ({side_label}) references atom index {idx_one_based} "
            f"(key {msg_key}) which is not present in the pocket after extraction. "
            f"Increase extraction coverage (e.g., --radius/--radius-het2het, --selected_resn, or set --exclude-backbone False), "
            f"or choose atoms that survive in the pocket."
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
    Return True if the PDB has at least one ATOM/HETATM record whose element field (cols 77–78) is empty.
    This is a light-weight check to decide whether to run add_elem_info.
    """
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            saw_atom = False
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    saw_atom = True
                    if len(line) < 78 or not line[76:78].strip():
                        return True
        # If no ATOM/HETATM was seen, fall back to "no fix"
        return False
    except Exception:
        # On I/O errors, skip fixing (use original)
        return False


# ---------- Post-processing helpers (minimal, reuse internals) ----------

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
    # atom order taken from first model
    atoms0 = [a for a in models[0].get_atoms()]
    elems: List[str] = []
    for a in atoms0:
        el = (a.element or "").strip()
        if not el:
            # fall back: derive from atom name
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


def _geom_from_angstrom(elems: Sequence[str],
                        coords_ang: np.ndarray,
                        freeze_atoms: Sequence[int]) -> Any:
    """
    Create a Geometry from Å coordinates using _path_search._new_geom_from_coords (expects Bohr).
    """
    coords_bohr = np.asarray(coords_ang, dtype=float) / BOHR2ANG
    return _path_search._new_geom_from_coords(elems, coords_bohr, coord_type="cart", freeze_atoms=freeze_atoms)


def _load_segment_end_geoms(seg_pdb: Path, freeze_atoms: Sequence[int]) -> Tuple[Any, Any]:
    """
    Load first/last model as Geometries from a per-segment pocket PDB.
    """
    coords_list, elems = _pdb_models_to_coords_and_elems(seg_pdb)
    gL = _geom_from_angstrom(elems, coords_list[0], freeze_atoms)
    gR = _geom_from_angstrom(elems, coords_list[-1], freeze_atoms)
    return gL, gR


def _compute_imag_mode_direction(ts_geom: Any,
                                 uma_kwargs: Dict[str, Any],
                                 freeze_atoms: Sequence[int]) -> np.ndarray:
    """
    Compute imaginary mode direction (N×3, unit vector in Cartesian space) at TS geometry.
    Uses tsopt internal helpers to minimize new code.
    """
    # ---- FIX: resolve torch.device safely (avoid "auto") ----
    dev_arg = uma_kwargs.get("device", None)
    try:
        if dev_arg is None or (isinstance(dev_arg, str) and dev_arg.strip().lower() in ("", "auto")):
            device_obj = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            device_obj = torch.device(dev_arg)
    except Exception:
        device_obj = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # full analytic Hessian (torch tensor)
    H_t = _tsopt._calc_full_hessian_torch(
        ts_geom,
        uma_kwargs=uma_kwargs,
        device=device_obj,
    )
    coords_bohr_t = torch.as_tensor(ts_geom.coords.reshape(-1, 3), dtype=H_t.dtype, device=H_t.device)
    # masses in a.u.
    from ase.data import atomic_masses
    masses_amu = np.array([atomic_masses[z] for z in ts_geom.atomic_numbers])
    masses_au_t = torch.as_tensor(masses_amu * _tsopt.AMU2AU, dtype=H_t.dtype, device=H_t.device)
    mode = _tsopt._mode_direction_by_root(
        H_t,
        coords_bohr_t,
        masses_au_t,
        root=0,
        freeze_idx=list(freeze_atoms) if len(freeze_atoms) > 0 else None,
    )
    # ensure unit length
    norm = float(np.linalg.norm(mode.reshape(-1)))
    if norm <= 0:
        raise click.ClickException("[post] Imaginary mode direction has zero norm.")
    return (mode / norm)


def _displaced_geometry_along_mode(geom: Any,
                                   mode_xyz: np.ndarray,
                                   amplitude_ang: float,
                                   freeze_atoms: Sequence[int]) -> Any:
    """
    Displace geometry along mode by ± amplitude (Å). Returns new Geometry.
    """
    coords_bohr = np.asarray(geom.coords3d, dtype=float)  # Bohr
    disp_bohr = (amplitude_ang / BOHR2ANG) * np.asarray(mode_xyz, dtype=float)  # (N,3)
    new_coords_bohr = coords_bohr + disp_bohr
    return _path_search._new_geom_from_coords(geom.atoms, new_coords_bohr, coord_type=geom.coord_type, freeze_atoms=freeze_atoms)


def _save_single_geom_as_pdb_for_tools(g: Any, ref_pdb: Path, out_dir: Path, name: str) -> Path:
    """
    Write a single-geometry XYZ/TRJ with energy and convert to PDB using the pocket ref (for downstream CLI tools).
    Returns PDB path.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_trj = out_dir / f"{name}.trj"
    _path_search._write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)
    pdb_out = out_dir / f"{name}.pdb"
    _path_search._maybe_convert_to_pdb(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
    return pdb_out


# ==== 追加: 汎用ファイル探索（最小限） ===============================================

def _find_with_suffixes(base_no_ext: Path, suffixes: Sequence[str]) -> Optional[Path]:
    """
    Given a base path without extension, return the first existing file among base.suffix for suffixes.
    """
    for s in suffixes:
        p = base_no_ext.with_suffix(s)
        if p.exists():
            return p
    return None


# -------------------------------------------------------------------------------------

def _run_tsopt_on_hei(hei_pdb: Path,
                      charge: int,
                      spin: int,
                      args_yaml: Optional[Path],
                      out_dir: Path,
                      freeze_links: bool,
                      opt_mode_default: str,
                      overrides: Optional[Dict[str, Any]] = None) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI pocket structure; return (final_geom_path, ts_geom).

    【変更点（最小限）】
      - 入力ファイルは拡張子を尊重してそのまま渡す（.pdb/.xyz/.gjf を想定）
      - 出力は final_geometry.{pdb|xyz|gjf} の存在を優先順で検出し、存在するものを返す
      - 非PDB入力時は、PDB変換はリファレンスPDBが無い場合は行わず、そのまま読み込む
    """
    overrides = overrides or {}
    prepared_input = prepare_input_structure(hei_pdb)
    template = prepared_input.gjf_template
    ts_dir = _resolve_override_dir(out_dir / "ts", overrides.get("out_dir"))
    _ensure_dir(ts_dir)

    freeze_use = overrides.get("freeze_links")
    if freeze_use is None:
        freeze_use = freeze_links

    opt_mode = overrides.get("opt_mode", opt_mode_default)

    ts_args: List[str] = [
        "-i", str(hei_pdb),  # ★ 修正: prepared_input.source_path ではなく元のパスをそのまま渡す
        "-q", str(int(charge)),
        "-s", str(int(spin)),
        "--freeze-links", "True" if freeze_use else "False",
        "--out-dir", str(ts_dir),
    ]

    if opt_mode is not None:
        ts_args.extend(["--opt-mode", str(opt_mode)])

    _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
    _append_cli_arg(ts_args, "--dump", overrides.get("dump"))

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

    # 出力候補を優先順に検出
    ts_pdb = ts_dir / "final_geometry.pdb"
    ts_xyz = ts_dir / "final_geometry.xyz"
    ts_gjf = ts_dir / "final_geometry.gjf"

    ts_geom_path: Optional[Path] = None

    if ts_pdb.exists():
        ts_geom_path = ts_pdb
    elif ts_xyz.exists():
        # 入力がPDBのときだけ、参照PDBを用いてPDB変換を試みる
        if hei_pdb.suffix.lower() == ".pdb":
            try:
                _path_search._maybe_convert_to_pdb(ts_xyz, hei_pdb, ts_pdb)
                if ts_pdb.exists():
                    ts_geom_path = ts_pdb
                else:
                    ts_geom_path = ts_xyz
            except Exception:
                ts_geom_path = ts_xyz
        else:
            ts_geom_path = ts_xyz
    elif ts_gjf.exists():
        ts_geom_path = ts_gjf
    else:
        # 旧ロジック: 非PDB出力でも XYZ が無ければエラー
        raise click.ClickException("[tsopt] TS outputs not found.")

    # 読み込み
    g_ts = geom_loader(ts_geom_path, coord_type="cart")

    # 可能なら GJF も書き出し（テンプレートがあり、XYZ が存在する場合）
    if (template is not None) and ts_xyz.exists():
        try:
            final_gjf = ts_dir / "final_geometry.gjf"
            maybe_convert_xyz_to_gjf(ts_xyz, template, final_gjf)
            click.echo(f"[tsopt] Wrote '{final_gjf}'.")
        except Exception as e:
            click.echo(f"[tsopt] WARNING: Failed to convert TS geometry to GJF: {e}", err=True)

    # Ensure calculator to have energy on g_ts
    calc = uma_pysis(charge=int(charge), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
    g_ts.set_calculator(calc)
    _ = float(g_ts.energy)

    prepared_input.cleanup()
    return ts_geom_path, g_ts


def _pseudo_irc_and_match(seg_idx: int,
                          seg_dir: Path,
                          ref_pdb_for_seg: Path,
                          seg_pocket_pdb: Path,
                          g_ts: Any,
                          q_int: int,
                          spin: int,
                          freeze_links_flag: bool) -> Dict[str, Any]:
    """
    From a TS pocket geometry, perform pseudo-IRC:
      - compute imag. mode
      - displace ± (0.25 Å) and optimize both to minima (LBFGS)
      - map each min to left/right segment endpoint by bond-change check (if segment endpoints exist)
        *If no segment endpoints are available (single-structure TSOPT-only mode), fall back to (plus, minus).*
    Returns dict with paths/energies/geoms for {left, ts, right}, and small IRC plots.
    """
    # Freeze parents of link-H if requested
    freeze_atoms: List[int] = []
    if freeze_links_flag and seg_pocket_pdb.suffix.lower() == ".pdb":
        freeze_atoms = detect_freeze_links_safe(seg_pocket_pdb)

    # Mode direction
    uma_kwargs = dict(charge=int(q_int), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
    mode_xyz = _compute_imag_mode_direction(g_ts, uma_kwargs=uma_kwargs, freeze_atoms=freeze_atoms)

    # Displace ± and optimize
    irc_dir = seg_dir / "irc"
    _ensure_dir(irc_dir)
    amp = 0.25  # Å; small stable displacement
    g_plus0 = _displaced_geometry_along_mode(g_ts,  mode_xyz, +amp, freeze_atoms)
    g_minus0 = _displaced_geometry_along_mode(g_ts, mode_xyz, -amp, freeze_atoms)

    # Shared UMA calc
    shared_calc = uma_pysis(charge=int(q_int), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
    # LBFGS settings (reuse defaults)
    sopt_cfg = dict(_path_search.LBFGS_KW)
    sopt_cfg["dump"] = True
    sopt_cfg["out_dir"] = str(irc_dir)

    # Optimize
    g_plus  = _path_search._optimize_single(g_plus0, shared_calc, "lbfgs", sopt_cfg, irc_dir, tag=f"seg_{seg_idx:02d}_irc_plus",  ref_pdb_path=seg_pocket_pdb)
    g_minus = _path_search._optimize_single(g_minus0, shared_calc, "lbfgs", sopt_cfg, irc_dir, tag=f"seg_{seg_idx:02d}_irc_minus", ref_pdb_path=seg_pocket_pdb)

    # IRC mini plots (TS→min)
    try:
        trj_plus  = irc_dir / f"seg_{seg_idx:02d}_irc_plus_opt/optimization.trj"
        trj_minus = irc_dir / f"seg_{seg_idx:02d}_irc_minus_opt/optimization.trj"
        if trj_plus.exists():
            run_trj2fig(trj_plus, [irc_dir / f"irc_plus_plot.png"], unit="kcal", reference="init", reverse_x=False)
        if trj_minus.exists():
            run_trj2fig(trj_minus, [irc_dir / f"irc_minus_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to plot IRC mini plots: {e}", err=True)

    # Try to load segment endpoints (pocket-only) — available only after path_search
    gL_end = None
    gR_end = None
    seg_pocket_path = seg_dir.parent / f"mep_seg_{seg_idx:02d}.pdb"
    if seg_pocket_path.exists():
        try:
            gL_end, gR_end = _load_segment_end_geoms(seg_pocket_path, freeze_atoms)
        except Exception as e:
            click.echo(f"[post] WARNING: failed to load segment endpoints: {e}", err=True)

    # Decide mapping
    bond_cfg = dict(_path_search.BOND_KW)

    def _matches(x, y) -> bool:
        try:
            chg, _ = _path_search._has_bond_change(x, y, bond_cfg)
            return (not chg)
        except Exception:
            # fallback: small RMSD threshold
            return (_path_search._rmsd_between(x, y, align=True) < 1e-3)

    candidates = [("plus", g_plus), ("minus", g_minus)]
    mapping: Dict[str, Any] = {"left": None, "right": None}

    if (gL_end is not None) and (gR_end is not None):
        # First pass: exact match on bond changes
        for tag, g in candidates:
            if _matches(g, gL_end) and not _matches(g, gR_end):
                mapping["left"] = (tag, g)
            elif _matches(g, gR_end) and not _matches(g, gL_end):
                mapping["right"] = (tag, g)
        # Second pass: fill missing by RMSD
        for side, g_end in (("left", gL_end), ("right", gR_end)):
            if mapping[side] is None:
                remain = [(t, gg) for (t, gg) in candidates if mapping.get("left", (None, None))[0] != t and mapping.get("right", (None, None))[0] != t]
                if not remain:
                    remain = candidates
                best = min(remain, key=lambda p: _path_search._rmsd_between(p[1], g_end, align=True))
                mapping[side] = best
    else:
        # Fallback (single-structure TSOPT-only mode): keep a deterministic assignment
        mapping["left"] = ("minus", g_minus)
        mapping["right"] = ("plus", g_plus)

    # Energies (ensure calculator)
    for _, g in candidates:
        _path_search._ensure_calc_on_geom(g, shared_calc)
        _ = float(g.energy)
    _path_search._ensure_calc_on_geom(g_ts, shared_calc); _ = float(g_ts.energy)

    # Dump tiny TS↔min trj for each direction
    try:
        for side in ("left", "right"):
            tag, gmin = mapping[side]
            trj = irc_dir / f"irc_{side}.trj"
            _path_search._write_xyz_trj_with_energy([g_ts, gmin], [float(g_ts.energy), float(gmin.energy)], trj)
            run_trj2fig(trj, [irc_dir / f"irc_{side}_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception:
        pass

    return {
        "left_min_geom": mapping["left"][1],
        "right_min_geom": mapping["right"][1],
        "ts_geom": g_ts,
        "left_tag": mapping["left"][0],
        "right_tag": mapping["right"][0],
        "freeze_atoms": freeze_atoms,
    }


def _write_segment_energy_diagram(prefix: Path,
                                  labels: List[str],
                                  energies_eh: List[float],
                                  title_note: str) -> None:
    """
    Write energy diagram (HTML + PNG) using utils.build_energy_diagram.
    """
    if not energies_eh:
        return
    e0 = energies_eh[0]
    energies_kcal = [(e - e0) * AU2KCALPERMOL for e in energies_eh]
    fig = build_energy_diagram(
        energies=energies_kcal,
        labels=labels,
        ylabel="ΔE (kcal/mol)",
        baseline=True,
        showgrid=False,
        title=f"{prefix.name} {title_note}",
    )
    html = prefix.with_suffix(".html")
    png = prefix.with_suffix(".png")
    fig.write_html(str(html))
    try:
        fig.write_image(str(png), scale=2)
    except Exception:
        pass
    click.echo(f"[diagram] Wrote energy diagram → {html.name} / {png.name}")


def _run_freq_for_state(pdb_path: Path,
                        q_int: int,
                        spin: int,
                        out_dir: Path,
                        args_yaml: Optional[Path],
                        freeze_links: bool,
                        overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
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
        "-i", str(pdb_path),
        "-q", str(int(q_int)),
        "-s", str(int(spin)),
        "--freeze-links", "True" if freeze_use else "False",
        "--out-dir", str(fdir),
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
    # parse thermoanalysis.yaml if any
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _run_dft_for_state(pdb_path: Path,
                       q_int: int,
                       spin: int,
                       out_dir: Path,
                       args_yaml: Optional[Path],
                       func_basis: str = "wb97x-v/def2-tzvp",
                       overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Run dft CLI; return parsed result.yaml dict (may be empty).
    """
    ddir = out_dir
    _ensure_dir(ddir)
    overrides = overrides or {}

    func_basis_use = overrides.get("func_basis", func_basis)

    args = [
        "-i", str(pdb_path),
        "-q", str(int(q_int)),
        "-s", str(int(spin)),
        "--func-basis", str(func_basis_use),
        "--out-dir", str(ddir),
    ]

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


# -----------------------------
# CLI
# -----------------------------

@click.command(
    help="Run pocket extraction → (optional single-structure staged scan) → MEP search → merge to full PDBs in one shot.\n"
         "If exactly one input is provided: (a) with --scan-lists, stage results feed into path_search; "
         "(b) with --tsopt True and no --scan-lists, run TSOPT-only mode.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
# ===== Inputs =====
@click.option(
    "-i", "--input", "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True, required=True,
    help=("Two or more **full** PDBs in reaction order (reactant [intermediates ...] product), "
          "or a single **full** PDB (with --scan-lists or with --tsopt True). "
          "You may pass a single '-i' followed by multiple space-separated files (e.g., '-i A.pdb B.pdb C.pdb').")
)
@click.option(
    "-c", "--center", "center_spec",
    type=str, required=False, default=None,
    help=("Substrate specification for the extractor: "
          "a PDB path, a residue-ID list like '123,124' or 'A:123,B:456' "
          "(insertion codes OK: '123A' / 'A:123A'), "
          "or a residue-name list like 'GPP,MMT'. "
          "When omitted, extraction is skipped and the **full input structure(s)** are used directly as pockets.")
)
@click.option(
    "--out-dir", "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("./result_all/"), show_default=True,
    help="Top-level output directory for the pipeline."
)
# ===== Extractor knobs (subset of extract.parse_args) =====
@click.option("-r", "--radius", type=float, default=2.6, show_default=True,
              help="Inclusion cutoff (Å) around substrate atoms.")
@click.option("--radius-het2het", type=float, default=0.0, show_default=True,
              help="Independent hetero–hetero cutoff (Å) for non‑C/H pairs.")
@click.option("--include-H2O", "--include-h2o", "include_h2o", type=click.BOOL, default=True, show_default=True,
              help="Include waters (HOH/WAT/TIP3/SOL) in the pocket.")
@click.option("--exclude-backbone", "exclude_backbone", type=click.BOOL, default=True, show_default=True,
              help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).")
@click.option("--add-linkH", "add_linkh", type=click.BOOL, default=True, show_default=True,
              help="Add link hydrogens for severed bonds (carbon-only) in pockets.")
@click.option("--selected_resn", type=str, default="", show_default=True,
              help="Force-include residues (comma/space separated; chain/insertion codes allowed).")
@click.option("--ligand-charge", type=str, default=None,
              help=("Either a total charge (number) to distribute across unknown residues "
                    "or a mapping like 'GPP:-3,MMT:-1'."))
@click.option("--verbose", type=click.BOOL, default=True, show_default=True, help="Enable INFO-level logging inside extractor.")
# ===== Path search knobs (subset of path_search.cli) =====
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1).")
@click.option("--freeze-links", "freeze_links_flag", type=click.BOOL, default=True, show_default=True,
              help="For pocket PDB input, freeze parent atoms of link hydrogens.")
@click.option("--max-nodes", type=int, default=10, show_default=True,
              help="Max internal nodes for **segment** GSM (String has max_nodes+2 images including endpoints).")
@click.option("--max-cycles", type=int, default=100, show_default=True, help="Maximum GSM optimization cycles.")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Enable transition-state climbing after growth for the **first** segment in each pair.")
@click.option("--sopt-mode", type=click.Choice(["lbfgs", "rfo", "light", "heavy"], case_sensitive=False),
              default="lbfgs", show_default=True,
              help="Single-structure optimizer kind for HEI±1 and kink nodes.")
@click.option("--opt-mode", type=str, default=None,
              help="Common optimizer mode forwarded to scan/tsopt (--opt-mode). When unset, tools use their defaults.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump GSM / single-structure trajectories during the run, forwarding the same flag to scan/tsopt/freq.")
@click.option("--args-yaml", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="YAML with extra args for path_search (sections: geom, calc, gs, opt, sopt, bond, search).")
@click.option("--preopt", "preopt", type=click.BOOL, default=True, show_default=True,
              help="If False, skip initial single-structure optimizations of the pocket inputs.")
@click.option("--hessian-calc-mode",
              type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
              default=None,
              help="Common UMA Hessian calculation mode forwarded to tsopt and freq.")
# ===== Post-processing toggles =====
@click.option("--tsopt", "do_tsopt", type=click.BOOL, default=False, show_default=True,
              help="TS optimization + pseudo-IRC per reactive segment (or TSOPT-only mode for single-structure), and build energy diagrams.")
@click.option("--thermo", "do_thermo", type=click.BOOL, default=False, show_default=True,
              help="Run freq on (R,TS,P) per reactive segment (or TSOPT-only mode) and build Gibbs free-energy diagram (UMA).")
@click.option("--dft", "do_dft", type=click.BOOL, default=False, show_default=True,
              help="Run DFT single-point on (R,TS,P) and build DFT energy diagram. With --thermo True, also generate a DFT//UMA Gibbs diagram.")
              
@click.option("--tsopt-max-cycles", type=int, default=None,
              help="Override tsopt --max-cycles value.")
@click.option("--tsopt-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override tsopt output subdirectory (relative paths are resolved against the default).")
@click.option("--freq-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override freq output base directory (relative paths resolved against the default).")
@click.option("--freq-max-write", type=int, default=None,
              help="Override freq --max-write value.")
@click.option("--freq-amplitude-ang", type=float, default=None,
              help="Override freq --amplitude-ang (Å).")
@click.option("--freq-n-frames", type=int, default=None,
              help="Override freq --n-frames value.")
@click.option("--freq-sort", type=click.Choice(["value", "abs"], case_sensitive=False), default=None,
              help="Override freq mode sorting.")
@click.option("--freq-temperature", type=float, default=None,
              help="Override freq thermochemistry temperature (K).")
@click.option("--freq-pressure", type=float, default=None,
              help="Override freq thermochemistry pressure (atm).")
@click.option("--dft-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override dft output base directory (relative paths resolved against the default).")
@click.option("--dft-func-basis", type=str, default=None,
              help="Override dft --func-basis value.")
@click.option("--dft-max-cycle", type=int, default=None,
              help="Override dft --max-cycle value.")
@click.option("--dft-conv-tol", type=float, default=None,
              help="Override dft --conv-tol value.")
@click.option("--dft-grid-level", type=int, default=None,
              help="Override dft --grid-level value.")
@click.option(
    "--scan-lists", "scan_lists_raw",
    type=str, multiple=True, required=False,
    help='Python-like list of (i,j,target_Å) per stage for **single-structure** scan. Repeatable. '
         'Example: "[(12,45,1.35)]" "--scan-lists \'[(10,55,2.20),(23,34,1.80)]\'". '
         'Indices refer to the original all-input PDB (1-based); they are auto-mapped to the pocket after extraction. '
         'Stage results feed into path_search.',
)
@click.option("--scan-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override the scan output directory (default: <out-dir>/scan/). Relative paths are resolved against the default parent.")
@click.option("--scan-one-based", type=click.BOOL, default=None,
              help="Override scan indexing interpretation (True = 1-based, False = 0-based).")
@click.option("--scan-max-step-size", type=float, default=None,
              help="Override scan --max-step-size (Å).")
@click.option("--scan-bias-k", type=float, default=None,
              help="Override scan harmonic bias strength k (eV/Å^2).")
@click.option("--scan-relax-max-cycles", type=int, default=None,
              help="Override scan relaxation max cycles per step.")
@click.option("--scan-preopt", "scan_preopt_override", type=click.BOOL, default=None,
              help="Override scan --preopt flag.")
@click.option("--scan-endopt", "scan_endopt_override", type=click.BOOL, default=None,
              help="Override scan --endopt flag.")
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
    verbose: bool,
    spin: int,
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    sopt_mode: str,
    opt_mode: Optional[str],
    dump: bool,
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
) -> None:
    """
    The **all** command composes `extract` → (optional `scan` on pocket) → `path_search` and hides ref-template bookkeeping.
    It also accepts the sloppy `-i A B C` style like `path_search` does. With single input:
      - with --scan-lists: run staged scan on the pocket and use stage results as inputs for path_search,
      - with --tsopt True and no --scan-lists: run TSOPT-only mode (no path_search).
    """
    time_start = time.perf_counter()

    dump_override_requested = False
    try:
        dump_source = ctx.get_parameter_source("dump")
        dump_override_requested = dump_source not in (None, ParameterSource.DEFAULT)
    except Exception:
        dump_override_requested = False

    # --- Robustly accept a single "-i" followed by multiple paths (like path_search.cli) ---
    argv_all = sys.argv[1:]
    i_vals = _collect_option_values(argv_all, ("-i", "--input"))
    if i_vals:
        i_parsed: List[Path] = []
        for tok in i_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Input path '{tok}' not found or is a directory. "
                    f"When using '-i', list only existing file paths (multiple paths may follow a single '-i')."
                )
            i_parsed.append(p)
        input_paths = tuple(i_parsed)

    # --------------------------
    # Validate input count / single-structure modes
    # --------------------------
    is_single = (len(input_paths) == 1)
    has_scan = bool(scan_lists_raw)
    single_tsopt_mode = (is_single and (not has_scan) and do_tsopt)

    if (len(input_paths) < 2) and (not (is_single and (has_scan or do_tsopt))):
        raise click.BadParameter(
            "Provide at least two PDBs with -i/--input in reaction order, "
            "or use a single PDB with --scan-lists, or a single PDB with --tsopt True."
        )

    tsopt_opt_mode_default = opt_mode.lower() if opt_mode else "light"
    tsopt_overrides: Dict[str, Any] = {}
    if tsopt_max_cycles is not None:
        tsopt_overrides["max_cycles"] = int(tsopt_max_cycles)
    if dump_override_requested:
        tsopt_overrides["dump"] = bool(dump)
    if tsopt_out_dir is not None:
        tsopt_overrides["out_dir"] = tsopt_out_dir
    if hessian_calc_mode is not None:
        tsopt_overrides["hessian_calc_mode"] = hessian_calc_mode
    if opt_mode is not None:
        tsopt_overrides["opt_mode"] = tsopt_opt_mode_default

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

    dft_func_basis_use = dft_func_basis or "wb97x-v/def2-tzvp"

    # --------------------------
    # Prepare directories
    # --------------------------
    out_dir = out_dir.resolve()
    pockets_dir = out_dir / "pockets"
    path_dir = out_dir / "path_search"
    scan_dir = _resolve_override_dir(out_dir / "scan", scan_out_dir)  # for single-structure scan mode
    _ensure_dir(out_dir)
    _ensure_dir(pockets_dir)
    if not single_tsopt_mode:
        _ensure_dir(path_dir)  # path_search might be skipped only in tsopt-only mode

    # --------------------------
    # Preflight: add_elem_info only for inputs lacking element fields
    # → Create fixed copies under a temporary folder inside out_dir (usable for extraction or direct use)
    # --------------------------
    elem_tmp_dir = out_dir / "add_elem_info"
    inputs_for_extract: List[Path] = []
    elem_fix_echo=False
    for p in input_paths:
        if _pdb_needs_elem_fix(p):
            if elem_fix_echo==False:
                click.echo("\n=== [all] Preflight — add_elem_info (only when element fields are missing) ===\n")
                elem_fix_echo=True
            _ensure_dir(elem_tmp_dir)
            out_p = (elem_tmp_dir / p.name).resolve()
            try:
                _assign_elem_info(str(p), str(out_p), overwrite=False)
                click.echo(f"[all] add_elem_info: fixed elements → {out_p}")
                inputs_for_extract.append(out_p)
            except SystemExit as e:
                code = getattr(e, "code", 1)
                click.echo(f"[all] WARNING: add_elem_info exited with code {code} for {p}; using original.", err=True)
                inputs_for_extract.append(p.resolve())
            except Exception as e:
                click.echo(f"[all] WARNING: add_elem_info failed for {p}: {e} — using original file.", err=True)
                inputs_for_extract.append(p.resolve())
        else:
            inputs_for_extract.append(p.resolve())

    extract_inputs = tuple(inputs_for_extract)

    # Determine if extraction is skipped
    skip_extract = (center_spec is None) or (str(center_spec).strip() == "")

    # --------------------------
    # Stage 1: Extraction (or skip)
    # --------------------------
    # Build per-structure pocket output file list (one per input full PDB for extraction)
    pocket_outputs: List[Path] = []
    for p in extract_inputs:
        pocket_outputs.append((pockets_dir / f"pocket_{p.stem}.pdb").resolve())

    if not skip_extract:
        click.echo("\n=== [all] Stage 1/3 — Active-site pocket extraction (multi-structure union when applicable) ===\n")
        # Run extractor via its public API (multi-structure union mode)
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

        # Report extractor outputs and charge breakdown
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
            click.echo(f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}")
            q_int = _round_charge_with_note(q_total)
        except Exception as e:
            raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")
    else:
        click.echo("\n=== [all] Stage 1 — Extraction skipped (no -c/--center); using FULL structures as pockets ===\n")
        # ---- New: Prefer GJF charge/spin (first input) when available ----
        first_input = input_paths[0].resolve()
        gjf_charge: Optional[int] = None
        gjf_spin: Optional[int] = None
        if first_input.suffix.lower() == ".gjf":
            try:
                with prepare_input_structure(first_input) as prepared:
                    gjf_charge, gjf_spin = resolve_charge_spin_or_raise(prepared, charge=None, spin=None)
                click.echo(f"[all] Detected from GJF (first input): charge={gjf_charge:+d}, spin={gjf_spin}")
            except Exception as e:
                click.echo(f"[all] NOTE: failed to parse charge/spin from GJF '{first_input.name}': {e}", err=True)

        # Determine if user explicitly provided spin
        user_provided_spin = True
        try:
            spin_source = ctx.get_parameter_source("spin")
            user_provided_spin = (spin_source not in (None, ParameterSource.DEFAULT))
        except Exception:
            user_provided_spin = True  # be conservative

        # Charge precedence: numeric --ligand-charge > GJF (if any) > 0
        q_total_fallback: float
        numeric_ligand_charge: Optional[float] = None
        if ligand_charge is not None:
            try:
                numeric_ligand_charge = float(ligand_charge)
            except Exception:
                numeric_ligand_charge = None

        if numeric_ligand_charge is not None:
            q_total_fallback = numeric_ligand_charge
            click.echo(f"[all] Using --ligand-charge as TOTAL system charge: {q_total_fallback:+g}")
        elif gjf_charge is not None:
            q_total_fallback = float(gjf_charge)
            click.echo(f"[all] Using total charge from first GJF: {q_total_fallback:+g}")
        else:
            q_total_fallback = 0.0
            if ligand_charge is not None:
                click.echo("[all] NOTE: non-numeric --ligand-charge is ignored without --center; defaulting total charge to 0 (no GJF charge available).", err=True)
            else:
                click.echo("[all] NOTE: No total charge provided; defaulting to 0. Supply '--ligand-charge <number>' to override.")
        q_int = _round_charge_with_note(q_total_fallback)

        # Spin precedence: explicit --spin > GJF (if any) > existing default
        if (not user_provided_spin) and (gjf_spin is not None):
            spin = int(gjf_spin)
            click.echo(f"[all] Spin multiplicity set from GJF: {spin}")

    # --------------------------
    # Other path: single-structure + --tsopt True (and NO scan-lists) → TSOPT-only mode
    # --------------------------
    if single_tsopt_mode:
        click.echo("\n=== [all] TSOPT-only single-structure mode ===\n")
        tsroot = out_dir / "tsopt_single"
        _ensure_dir(tsroot)

        # Use the single pocket PDB as TS initial guess (or full PDB when extraction is skipped)
        ts_initial_pdb = (pocket_outputs[0] if not skip_extract else input_paths[0].resolve())
        # TS optimization
        ts_pdb, g_ts = _run_tsopt_on_hei(
            ts_initial_pdb,
            q_int,
            spin,
            args_yaml,
            tsroot,
            freeze_links_flag,
            tsopt_opt_mode_default,
            overrides=tsopt_overrides,
        )

        # Pseudo-IRC & minimize both ends (no segment endpoints exist → fallback mapping in helper)
        irc_res = _pseudo_irc_and_match(
            seg_idx=1,
            seg_dir=tsroot,
            ref_pdb_for_seg=ts_pdb,
            seg_pocket_pdb=ts_initial_pdb,
            g_ts=g_ts,
            q_int=q_int,
            spin=spin,
            freeze_links_flag=freeze_links_flag,
        )
        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]

        # Ensure UMA energies
        eL = float(gL.energy)
        eT = float(gT.energy)
        eR = float(gR.energy)

        # In this mode ONLY: assign Reactant/Product so that higher-energy end is the Reactant
        if eL >= eR:
            g_react, e_react = gL, eL
            g_prod,  e_prod  = gR, eR
        else:
            g_react, e_react = gR, eR
            g_prod,  e_prod  = gL, eL

        # Save standardized PDBs
        struct_dir = tsroot / "structures"
        _ensure_dir(struct_dir)
        pocket_ref = ts_initial_pdb
        pR = _save_single_geom_as_pdb_for_tools(g_react, pocket_ref, struct_dir, "reactant")
        pT = _save_single_geom_as_pdb_for_tools(gT,       pocket_ref, struct_dir, "ts")
        pP = _save_single_geom_as_pdb_for_tools(g_prod,   pocket_ref, struct_dir, "product")

        # UMA energy diagram (R, TS, P)
        _write_segment_energy_diagram(
            tsroot / "energy_diagram_tsopt",
            labels=["R", "TS", "P"],
            energies_eh=[e_react, eT, e_prod],
            title_note="(UMA, TSOPT/IRC)",
        )

        # Thermochemistry (UMA) Gibbs
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        freq_root = _resolve_override_dir(tsroot / "freq", freq_out_dir)
        dft_root = _resolve_override_dir(tsroot / "dft", dft_out_dir)

        if do_thermo:
            click.echo(f"[thermo] Single TSOPT: freq on R/TS/P")
            tR = _run_freq_for_state(pR, q_int, spin, freq_root / "R",  args_yaml, freeze_links_flag, overrides=freq_overrides)
            tT = _run_freq_for_state(pT, q_int, spin, freq_root / "TS", args_yaml, freeze_links_flag, overrides=freq_overrides)
            tP = _run_freq_for_state(pP, q_int, spin, freq_root / "P",  args_yaml, freeze_links_flag, overrides=freq_overrides)
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            try:
                GR = float(tR.get("sum_EE_and_thermal_free_energy_ha", e_react))
                GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", eT))
                GP = float(tP.get("sum_EE_and_thermal_free_energy_ha", e_prod))
                _write_segment_energy_diagram(
                    tsroot / "energy_diagram_G_UMA",
                    labels=["R", "TS", "P"],
                    energies_eh=[GR, GT, GP],
                    title_note="(Gibbs, UMA)",
                )
            except Exception as e:
                click.echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)

        # DFT & DFT//UMA
        if do_dft:
            click.echo(f"[dft] Single TSOPT: DFT on R/TS/P")
            dR = _run_dft_for_state(pR, q_int, spin, dft_root / "R",  args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            dT = _run_dft_for_state(pT, q_int, spin, dft_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            dP = _run_dft_for_state(pP, q_int, spin, dft_root / "P",  args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            try:
                eR_dft = float(((dR or {}).get("energy", {}) or {}).get("hartree", e_react))
                eT_dft = float(((dT or {}).get("energy", {}) or {}).get("hartree", eT))
                eP_dft = float(((dP or {}).get("energy", {}) or {}).get("hartree", e_prod))
                _write_segment_energy_diagram(
                    tsroot / "energy_diagram_DFT",
                    labels=["R", "TS", "P"],
                    energies_eh=[eR_dft, eT_dft, eP_dft],
                    title_note="(DFT wb97x-v/def2-tzvp)",
                )
            except Exception as e:
                click.echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            if do_thermo:
                try:
                    dG_R = float((thermo_payloads.get("R",  {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_T = float((thermo_payloads.get("TS", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_P = float((thermo_payloads.get("P",  {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                    GR_dftUMA = eR_dft + dG_R
                    GT_dftUMA = eT_dft + dG_T
                    GP_dftUMA = eP_dft + dG_P
                    _write_segment_energy_diagram(
                        tsroot / "energy_diagram_G_DFT_plus_UMA",
                        labels=["R", "TS", "P"],
                        energies_eh=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                        title_note="(Gibbs, DFT//UMA)",
                    )
                except Exception as e:
                    click.echo(f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}", err=True)

        click.echo("\n=== [all] TSOPT-only pipeline finished successfully ===\n")
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # --------------------------
    # Stage 1b: Optional scan (single-structure only) to build ordered pocket inputs
    # --------------------------
    pockets_for_path: List[Path]
    if is_single and has_scan:
        click.echo("\n=== [all] Stage 1b — Staged scan on input ===\n")
        _ensure_dir(scan_dir)

        if skip_extract:
            # Use the FULL input directly; indices already refer to the full PDB.
            scan_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _parse_scan_lists_literals(scan_lists_raw)
        else:
            # Use the pocket and re-map indices from the full PDB to the pocket order.
            scan_input_pdb = Path(pocket_outputs[0]).resolve()
            full_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _convert_scan_lists_to_pocket_indices(
                scan_lists_raw, full_input_pdb, scan_input_pdb
            )
            click.echo("[all] Remapped --scan-lists indices from the full PDB to the pocket ordering.")

        scan_one_based_effective = True if scan_one_based is None else bool(scan_one_based)
        scan_stage_literals: List[str] = []
        for stage in converted_scan_stages:
            if scan_one_based_effective:
                stage_use = stage
            else:
                stage_use = [(i - 1, j - 1, target) for (i, j, target) in stage]
            scan_stage_literals.append(_format_scan_stage(stage_use))

        scan_preopt_use = preopt if scan_preopt_override is None else bool(scan_preopt_override)
        scan_endopt_use = True if scan_endopt_override is None else bool(scan_endopt_override)
        scan_opt_mode_use = (opt_mode.lower() if opt_mode else
                              ("lbfgs" if sopt_mode.lower() in ("lbfgs", "light") else "rfo"))

        scan_args: List[str] = [
            "-i", str(scan_input_pdb),
            "-q", str(int(q_int)),
            "-s", str(int(spin)),
            "--out-dir", str(scan_dir),
            "--freeze-links", "True" if freeze_links_flag else "False",
            "--preopt", "True" if scan_preopt_use else "False",
            "--endopt", "True" if scan_endopt_use else "False",
            "--opt-mode", str(scan_opt_mode_use),
        ]

        if dump_override_requested:
            scan_args.extend(["--dump", "True" if dump else "False"])

        if scan_one_based is not None:
            scan_args.append("--one-based" if scan_one_based else "--zero-based")

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        if args_yaml is not None:
            scan_args.extend(["--args-yaml", str(args_yaml)])
        # Forward all converted --scan-lists
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
                raise click.ClickException(f"[all] scan terminated with exit code {code}.")
        except Exception as e:
            raise click.ClickException(f"[all] scan failed: {e}")
        finally:
            sys.argv = _saved_argv

        # Collect stage results as pocket inputs
        stage_pdbs = sorted((scan_dir).glob("stage_*/*"))
        stage_results: List[Path] = []
        for p in stage_pdbs:
            if p.name == "result.pdb":
                stage_results.append(p.resolve())
        if not stage_results:
            raise click.ClickException("[all] No stage result PDBs found under scan/.")
        click.echo("[all] Collected scan stage pocket files:")
        for p in stage_results:
            click.echo(f"  - {p}")

        # Input series to path_search: [initial (pocket or full), stage_01/result.pdb, stage_02/result.pdb, ...]
        pockets_for_path = [scan_input_pdb] + stage_results
    else:
        # Multi-structure standard route
        if skip_extract:
            # Use FULL inputs (possibly element-fixed) directly
            pockets_for_path = [p.resolve() for p in inputs_for_extract]
        else:
            # Use per-input pocket PDBs
            pockets_for_path = list(pocket_outputs)

    # --------------------------
    # Stage 2: Path search on pockets (or full structures) (auto-supplying ref templates = original full PDBs)
    # --------------------------
    click.echo("\n=== [all] Stage 2/3 — MEP search on input structures (recursive GSM) ===\n")

    # Build path_search CLI args using *repeated* options (robust for Click)
    ps_args: List[str] = []

    # Inputs: repeat "-i" per pocket to satisfy Click even without their argv aggregator
    for p in pockets_for_path:
        ps_args.extend(["-i", str(p)])

    # Charge & spin
    ps_args.extend(["-q", str(q_int)])
    ps_args.extend(["-s", str(int(spin))])

    # Freeze-links, nodes, cycles, climb, optimizer, dump, out-dir, preopt, args-yaml
    ps_args.extend(["--freeze-links", "True" if freeze_links_flag else "False"])
    ps_args.extend(["--max-nodes", str(int(max_nodes))])
    ps_args.extend(["--max-cycles", str(int(max_cycles))])
    ps_args.extend(["--climb", "True" if climb else "False"])
    ps_args.extend(["--sopt-mode", str(sopt_mode)])
    ps_args.extend(["--dump", "True" if dump else "False"])
    ps_args.extend(["--out-dir", str(path_dir)])
    ps_args.extend(["--preopt", "True" if preopt else "False"])
    if args_yaml is not None:
        ps_args.extend(["--args-yaml", str(args_yaml)])

    # Auto-provide ref templates (original full PDBs) for full-system merge,
    # but only when the ALL inputs are PDB files.
    def _is_pdb(path: Path) -> bool:
        return path.suffix.lower() == ".pdb"

    gave_ref_pdb = False  # optional: used below to print accurate merge info

    # Multi-structure: one ref per original input; single+scan: reuse the single template for all pockets.
    if is_single and has_scan:
        if _is_pdb(input_paths[0]):
            for _ in pockets_for_path:
                ps_args.extend(["--ref-pdb", str(input_paths[0])])
            gave_ref_pdb = True
        else:
            click.echo("[all] NOTE: skipping --ref-pdb (single+scan: original input is not a PDB).")
    else:
        if all(_is_pdb(p) for p in input_paths):
            for p in input_paths:
                ps_args.extend(["--ref-pdb", str(p)])
            gave_ref_pdb = True
        else:
            click.echo("[all] NOTE: skipping --ref-pdb (one or more original inputs are not PDB).")

    click.echo("[all] Invoking path_search with arguments:")
    click.echo("  " + " ".join(ps_args))

    _saved_argv = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "path_search"] + ps_args
        _path_search.cli.main(args=ps_args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            raise click.ClickException(f"[all] path_search terminated with exit code {code}.")
    except Exception as e:
        raise click.ClickException(f"[all] path_search failed: {e}")
    finally:
        sys.argv = _saved_argv

    # --------------------------
    # Stage 3: Merge (performed by path_search when --ref-pdb was supplied)
    # --------------------------
    click.echo("\n=== [all] Stage 3/3 — Merge into full-system templates ===\n")
    if gave_ref_pdb:
        click.echo("[all] Merging was carried out by path_search using the original inputs as templates.")
        click.echo(f"[all] Final products can be found under: {path_dir}")
        click.echo("  - mep_w_ref.pdb            (full-system merged trajectory)")
        click.echo("  - mep_w_ref_seg_XX.pdb     (per-segment merged trajectories for covalent-change segments)")
    else:
        click.echo("[all] --ref-pdb was not provided; full-system merged trajectories are not produced.")
        click.echo(f"[all] Pocket-only outputs are under: {path_dir}")
    click.echo("  - summary.yaml             (segment barriers, ΔE, labels)")
    click.echo("  - energy.png / energy_diagram.*")
    click.echo("\n=== [all] Pipeline finished successfully (core path) ===\n")

    # --------------------------
    # Optional Stage 4: TSOPT / THERMO / DFT (per reactive segment)
    # --------------------------
    if not (do_tsopt or do_thermo or do_dft):
        # Elapsed time
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    click.echo("\n=== [all] Stage 4 — Post-processing per reactive segment ===\n")

    # Load segment summary
    summary_yaml = path_dir / "summary.yaml"
    segments = _read_summary(summary_yaml)
    if not segments:
        click.echo("[post] No segments found in summary; nothing to do.")
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # Iterate only bond-change segments (kind='seg' and bond_changes not empty and not '(no covalent...)')
    reactive = [s for s in segments if (s.get("kind", "seg") == "seg" and str(s.get("bond_changes", "")).strip() and str(s.get("bond_changes", "")).strip() != "(no covalent changes detected)")]
    if not reactive:
        click.echo("[post] No bond-change segments. Skipping TS/thermo/DFT.")
        click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # For each reactive segment
    for s in reactive:
        seg_idx = int(s.get("index", 0) or 0)
        seg_tag = s.get("tag", f"seg_{seg_idx:02d}")
        click.echo(f"\n--- [post] Segment {seg_idx:02d} ({seg_tag}) ---")

        seg_root = path_dir  # base
        seg_dir = seg_root / f"tsopt_seg_{seg_idx:02d}"
        _ensure_dir(seg_dir)

        # HEI pocket file prepared by path_search (only for bond-change segments)
        hei_base = seg_root / f"hei_seg_{seg_idx:02d}"
        hei_pocket_path = _find_with_suffixes(hei_base, [".pdb", ".xyz", ".gjf"])
        if hei_pocket_path is None:
            click.echo(f"[post] WARNING: HEI pocket file not found for segment {seg_idx:02d} (searched .pdb/.xyz/.gjf); skipping TSOPT.", err=True)
            continue

        struct_dir = seg_dir / "structures"
        _ensure_dir(struct_dir)
        state_structs: Dict[str, Path] = {}
        uma_ref_energies: Dict[str, float] = {}

        # 4.1 TS optimization (optional; still needed to drive IRC & diagrams)
        if do_tsopt:
            ts_pdb, g_ts = _run_tsopt_on_hei(
                hei_pocket_path,
                q_int,
                spin,
                args_yaml,
                seg_dir,
                freeze_links_flag,
                tsopt_opt_mode_default,
                overrides=tsopt_overrides,
            )

            # 4.2 Pseudo-IRC & mapping to (left,right)
            irc_res = _pseudo_irc_and_match(
                seg_idx=seg_idx,
                seg_dir=seg_dir,
                ref_pdb_for_seg=ts_pdb,
                seg_pocket_pdb=hei_pocket_path,
                g_ts=g_ts,
                q_int=q_int,
                spin=spin,
                freeze_links_flag=freeze_links_flag,
            )

            gL = irc_res["left_min_geom"]
            gR = irc_res["right_min_geom"]
            gT = irc_res["ts_geom"]
            # Save standardized PDBs for tools
            pL = _save_single_geom_as_pdb_for_tools(gL, hei_pocket_path, struct_dir, "reactant_like")
            pT = _save_single_geom_as_pdb_for_tools(gT, hei_pocket_path, struct_dir, "ts")
            pR = _save_single_geom_as_pdb_for_tools(gR, hei_pocket_path, struct_dir, "product_like")
            state_structs = {"R": pL, "TS": pT, "P": pR}

            # 4.3 Segment-level energy diagram from UMA (R,TS,P)
            eR = float(gL.energy)
            eT = float(gT.energy)
            eP = float(gR.energy)
            uma_ref_energies = {"R": eR, "TS": eT, "P": eP}
            _write_segment_energy_diagram(
                seg_dir / "energy_diagram_tsopt",
                labels=["R", f"TS{seg_idx}", "P"],
                energies_eh=[eR, eT, eP],
                title_note="(UMA, TSOPT/IRC)",
            )
        elif (do_thermo or do_dft):
            seg_pocket_path = _find_with_suffixes(seg_root / f"mep_seg_{seg_idx:02d}", [".pdb"])
            if seg_pocket_path is None:
                click.echo(
                    f"[post] WARNING: mep_seg_{seg_idx:02d}.pdb not found; cannot run thermo/DFT without --tsopt. Skipping segment.",
                    err=True,
                )
                continue

            freeze_atoms: List[int] = []
            if freeze_links_flag and seg_pocket_path.suffix.lower() == ".pdb":
                freeze_atoms = detect_freeze_links_safe(seg_pocket_path)

            try:
                gL, gR = _load_segment_end_geoms(seg_pocket_path, freeze_atoms)
            except Exception as e:
                click.echo(
                    f"[post] WARNING: failed to load segment endpoints from {seg_pocket_path.name}: {e}. Skipping segment.",
                    err=True,
                )
                continue

            try:
                g_ts = geom_loader(hei_pocket_path, coord_type="cart")
            except Exception as e:
                click.echo(
                    f"[post] WARNING: failed to load HEI geometry for segment {seg_idx:02d}: {e}. Skipping segment.",
                    err=True,
                )
                continue

            ref_for_ts = seg_pocket_path if seg_pocket_path.suffix.lower() == ".pdb" else hei_pocket_path
            pL = _save_single_geom_as_pdb_for_tools(gL, seg_pocket_path, struct_dir, "reactant_like")
            pR = _save_single_geom_as_pdb_for_tools(gR, seg_pocket_path, struct_dir, "product_like")
            pT = _save_single_geom_as_pdb_for_tools(g_ts, ref_for_ts, struct_dir, "ts_from_hei")
            state_structs = {"R": pL, "TS": pT, "P": pR}

        # 4.4 Thermochemistry (UMA freq) and Gibbs diagram
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        freq_seg_root = _resolve_override_dir(seg_dir / "freq", freq_out_dir)
        dft_seg_root = _resolve_override_dir(seg_dir / "dft", dft_out_dir)

        if (do_thermo or do_dft) and (not state_structs):
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
                click.echo(f"[thermo] Segment {seg_idx:02d}: freq on R/TS/P")
                tR = _run_freq_for_state(p_react, q_int, spin, freq_seg_root / "R", args_yaml, freeze_links_flag, overrides=freq_overrides)
                tT = _run_freq_for_state(p_ts, q_int, spin, freq_seg_root / "TS", args_yaml, freeze_links_flag, overrides=freq_overrides)
                tP = _run_freq_for_state(p_prod, q_int, spin, freq_seg_root / "P", args_yaml, freeze_links_flag, overrides=freq_overrides)
                thermo_payloads = {"R": tR, "TS": tT, "P": tP}
                try:
                    GR = float(tR.get("sum_EE_and_thermal_free_energy_ha", uma_ref_energies.get("R", np.nan)))
                    GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", uma_ref_energies.get("TS", np.nan)))
                    GP = float(tP.get("sum_EE_and_thermal_free_energy_ha", uma_ref_energies.get("P", np.nan)))
                    gibbs_vals = [GR, GT, GP]
                    if all(np.isfinite(gibbs_vals)):
                        _write_segment_energy_diagram(
                            seg_dir / "energy_diagram_G_UMA",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_eh=gibbs_vals,
                            title_note="(Gibbs, UMA)",
                        )
                    else:
                        click.echo("[thermo] NOTE: Gibbs energies non-finite; diagram skipped.")
                except Exception as e:
                    click.echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)

        # 4.5 DFT single-point and (optionally) DFT//UMA Gibbs
        if do_dft:
            if not (p_react and p_ts and p_prod):
                click.echo(
                    f"[dft] WARNING: Missing R/TS/P structures for segment {seg_idx:02d}; skipping DFT.",
                    err=True,
                )
            else:
                click.echo(f"[dft] Segment {seg_idx:02d}: DFT on R/TS/P")
                dR = _run_dft_for_state(p_react, q_int, spin, dft_seg_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
                dT = _run_dft_for_state(p_ts, q_int, spin, dft_seg_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
                dP = _run_dft_for_state(p_prod, q_int, spin, dft_seg_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
                try:
                    eR_dft = float(((dR or {}).get("energy", {}) or {}).get("hartree", np.nan))
                    eT_dft = float(((dT or {}).get("energy", {}) or {}).get("hartree", np.nan))
                    eP_dft = float(((dP or {}).get("energy", {}) or {}).get("hartree", np.nan))
                    if all(map(np.isfinite, [eR_dft, eT_dft, eP_dft])):
                        _write_segment_energy_diagram(
                            seg_dir / "energy_diagram_DFT",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_eh=[eR_dft, eT_dft, eP_dft],
                            title_note="(DFT wb97x-v/def2-tzvp)",
                        )
                    else:
                        click.echo("[dft] WARNING: some DFT energies missing; diagram skipped.", err=True)
                except Exception as e:
                    click.echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

                # DFT//UMA thermal Gibbs (E_DFT + ΔG_therm(UMA))
                if do_thermo:
                    try:
                        dG_R = float((thermo_payloads.get("R", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                        dG_T = float((thermo_payloads.get("TS", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                        dG_P = float((thermo_payloads.get("P", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                        eR_dft = float(((dR or {}).get("energy", {}) or {}).get("hartree", uma_ref_energies.get("R", np.nan)))
                        eT_dft = float(((dT or {}).get("energy", {}) or {}).get("hartree", uma_ref_energies.get("TS", np.nan)))
                        eP_dft = float(((dP or {}).get("energy", {}) or {}).get("hartree", uma_ref_energies.get("P", np.nan)))
                        GR_dftUMA = eR_dft + dG_R
                        GT_dftUMA = eT_dft + dG_T
                        GP_dftUMA = eP_dft + dG_P
                        if all(np.isfinite([GR_dftUMA, GT_dftUMA, GP_dftUMA])):
                            _write_segment_energy_diagram(
                                seg_dir / "energy_diagram_G_DFT_plus_UMA",
                                labels=["R", f"TS{seg_idx}", "P"],
                                energies_eh=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                                title_note="(Gibbs, DFT//UMA)",
                            )
                        else:
                            click.echo("[dft//uma] WARNING: DFT//UMA Gibbs energies non-finite; diagram skipped.", err=True)
                    except Exception as e:
                        click.echo(f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}", err=True)

    click.echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))


if __name__ == "__main__":
    cli()
