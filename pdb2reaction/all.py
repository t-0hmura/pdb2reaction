# pdb2reaction/all.py

"""
End-to-end enzymatic reaction workflow: extraction, MEP search, TS optimization, IRC, and post-processing.

Example:
    pdb2reaction all -i reactant.pdb product.pdb -c 'GPP,MMT' --ligand-charge 'GPP:-3,MMT:-1'

For detailed documentation, see: docs/all.md
"""

from __future__ import annotations

import ast
from pathlib import Path
from collections import defaultdict
from typing import List, Sequence, Optional, Tuple, Dict, Any

import sys, os
import math
import tempfile
import re
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
    collect_option_values,
    collect_single_option_values,
    convert_xyz_like_outputs,
    detect_freeze_links_logged,
    format_elapsed,
    merge_freeze_atom_groups,
    prepare_input_structure,
    normalize_freeze_atoms,
    set_convert_file_enabled,
    resolve_charge_spin_or_raise,
    load_yaml_dict,
    apply_yaml_overrides,
    load_pdb_atom_metadata,
    merge_freeze_atom_indices,
    _round_charge_with_note,
    apply_ref_pdb_override,
    ensure_dir,
    parse_scan_list_triples,
    close_matplotlib_figures,
    _derive_charge_from_ligand_charge,
    write_xyz_trj_with_energy,
    read_xyz_as_blocks,
    read_xyz_first_last,
    xyz_blocks_first_last,
    set_freeze_atoms_or_warn,
)
from . import scan as _scan_cli
from .add_elem_info import assign_elements as _assign_elem_info
from . import irc as _irc_cli


# -----------------------------
# Helpers
# -----------------------------


def _copy_logged(src: Path, dst: Path, *, label: Optional[str] = None, echo: bool = True) -> bool:
    """Copy files with consistent warning messages; return success."""
    try:
        shutil.copy2(src, dst)
        if echo:
            shown = label or src.name
            click.echo(f"[all] Copied {shown} → {dst}")
        return True
    except Exception as e:
        shown = label or src
        click.echo(f"[all] WARNING: Failed to copy {shown} to {dst}: {e}", err=True)
        return False


def _run_cli_main(
    cmd_name: str,
    cli_obj,
    args: Sequence[str],
    *,
    on_nonzero: str = "warn",
    on_exception: str = "raise",
    prefix: Optional[str] = None,
) -> None:
    """Run a Click command with a temporary argv and consistent error handling."""
    saved = list(sys.argv)
    label = prefix or cmd_name
    try:
        sys.argv = ["pdb2reaction", cmd_name] + list(args)
        print()
        cli_obj.main(args=list(args), standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            if on_nonzero == "raise":
                raise click.ClickException(f"[{label}] {cmd_name} exit code {code}.")
            click.echo(f"[{label}] WARNING: {cmd_name} exited with code {code}", err=True)
    except Exception as e:
        if on_exception == "raise":
            raise click.ClickException(f"[{label}] {cmd_name} failed: {e}")
        click.echo(f"[{label}] WARNING: {cmd_name} failed: {e}", err=True)
    finally:
        sys.argv = saved


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


def _build_calc_cfg(
    charge: int,
    spin: int,
    workers: Optional[int] = None,
    workers_per_node: Optional[int] = None,
    yaml_cfg: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return a UMA calculator configuration honoring YAML overrides when provided."""
    cfg: Dict[str, Any] = dict(CALC_KW)
    cfg["charge"] = int(charge)
    cfg["spin"] = int(spin)
    if workers is not None:
        cfg["workers"] = int(workers)
    if workers_per_node is not None:
        cfg["workers_per_node"] = int(workers_per_node)
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


def _parse_scan_lists_literals(
    scan_lists_raw: Sequence[str],
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
) -> List[List[Tuple[int, int, float]]]:
    """Parse ``--scan-lists`` literals without re-basing atom indices."""
    stages: List[List[Tuple[int, int, float]]] = []
    for idx_stage, literal in enumerate(scan_lists_raw, start=1):
        tuples, _ = parse_scan_list_triples(
            literal,
            one_based=True,
            atom_meta=atom_meta,
            option_name=f"--scan-lists #{idx_stage}",
            return_one_based=True,
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

    full_atom_meta = load_pdb_atom_metadata(full_input_pdb)
    stages = _parse_scan_lists_literals(scan_lists_raw, atom_meta=full_atom_meta)

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


def _pdb_needs_elem_fix(p: Path) -> bool:
    """
    Return True if the file contains ATOM/HETATM records and at least one has an empty element field (cols 77–78).
    This is a light-weight check to decide whether to run add_elem_info.
    Raises a ClickException if the file cannot be inspected.
    """
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if len(line) < 78 or not line[76:78].strip():
                        return True
        return False
    except Exception as e:
        raise click.ClickException(
            f"[all] Failed to inspect PDB element fields for '{p}': {e}"
        ) from e


_FREEZE_ATOMS_GLOBAL: Optional[List[int]] = None
_FREEZE_ATOMS_YAML: Optional[List[int]] = None


def _set_yaml_freeze_atoms(yaml_cfg: Optional[Dict[str, Any]]) -> None:
    """Cache freeze_atoms from args-yaml for merging with freeze-links."""
    global _FREEZE_ATOMS_YAML
    if not isinstance(yaml_cfg, dict):
        _FREEZE_ATOMS_YAML = []
        return
    geom_cfg = yaml_cfg.get("geom")
    if not isinstance(geom_cfg, dict):
        _FREEZE_ATOMS_YAML = []
        return
    _FREEZE_ATOMS_YAML = normalize_freeze_atoms(geom_cfg.get("freeze_atoms"))


def _get_freeze_atoms(pdb_path: Optional[Path], freeze_links_flag: bool) -> List[int]:
    """
    Determine freeze atom indices once and reuse them globally.

    The first time this is called with a PDB path and freeze_links_flag=True,
    link-parent atoms are detected from that PDB. The resulting indices are
    reused for subsequent calls (even if a non-PDB path is provided), under
    the assumption that atom indexing is consistent across the trajectory.
    """
    global _FREEZE_ATOMS_GLOBAL
    if freeze_links_flag:
        if _FREEZE_ATOMS_GLOBAL is not None:
            return merge_freeze_atom_groups(_FREEZE_ATOMS_GLOBAL, _FREEZE_ATOMS_YAML or [])
        if pdb_path is None or pdb_path.suffix.lower() != ".pdb":
            # No suitable PDB available yet to determine freeze atoms.
            return merge_freeze_atom_groups(_FREEZE_ATOMS_YAML or [])
        fa = detect_freeze_links_logged(pdb_path)
        _FREEZE_ATOMS_GLOBAL = [int(i) for i in fa]
        return merge_freeze_atom_groups(_FREEZE_ATOMS_GLOBAL, _FREEZE_ATOMS_YAML or [])
    return merge_freeze_atom_groups(_FREEZE_ATOMS_YAML or [])


def _freeze_atoms_for_log() -> List[int]:
    """Return a sorted freeze_atoms list for summary logs (may be empty)."""

    try:
        return merge_freeze_atom_groups(_FREEZE_ATOMS_GLOBAL or [], _FREEZE_ATOMS_YAML or [])
    except Exception:
        return []


def _write_args_yaml_with_freeze_atoms(
    args_yaml: Optional[Path],
    freeze_atoms: Sequence[int],
) -> Optional[Path]:
    """
    Merge ``freeze_atoms`` into a YAML config under ``geom`` and write a temporary YAML file.
    Returns the new YAML path, or the original ``args_yaml`` when no freeze atoms are provided.
    """
    if not freeze_atoms:
        return args_yaml

    cfg = {} if args_yaml is None else load_yaml_dict(args_yaml)
    if not isinstance(cfg, dict):
        cfg = {}

    geom_cfg = cfg.get("geom")
    if not isinstance(geom_cfg, dict):
        geom_cfg = {}
    geom_cfg = dict(geom_cfg)

    merge_freeze_atom_indices(geom_cfg, freeze_atoms)
    cfg["geom"] = geom_cfg

    tmp_dir = Path(tempfile.mkdtemp(prefix="tmp_path_search_"))
    out_path = tmp_dir / "args_freeze_atoms.yaml"
    with out_path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)
    return out_path


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

    set_freeze_atoms_or_warn(gL_ref, freeze_atoms, context="all")
    set_freeze_atoms_or_warn(gR_ref, freeze_atoms, context="all")

    return gL_ref, gR_ref


def _save_single_geom_as_pdb_for_tools(
    g: Any,
    ref_pdb: Path,
    out_dir: Path,
    name: str,
) -> Path:
    """
    Write a single-geometry XYZ trajectory with energy for downstream CLI tools.
    If a PDB reference is available, also write a PDB companion for visualization.
    Returns the XYZ path to avoid passing rounded PDB coordinates between steps.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_trj = out_dir / f"{name}.xyz"
    write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)

    if ref_pdb.suffix.lower() == ".pdb":
        pdb_out = out_dir / f"{name}.pdb"
        try:
            _path_search._convert_to_pdb_logged(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
        except Exception as e:
            click.echo(
                f"[all] WARNING: failed to convert '{xyz_trj.name}' to PDB for {name}: {e}",
                err=True,
            )

    return xyz_trj



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
    ensure_dir(out_path.parent)
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
            blocks = read_xyz_as_blocks(trj_path)
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
    ensure_dir(tmp_trj.parent)
    try:
        tmp_trj.write_text("\n".join(all_blocks) + "\n", encoding="utf-8")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: Failed to write concatenated IRC trajectory: {e}", err=True)
        return

    try:
        run_trj2fig(tmp_trj, [out_png], unit="kcal", reference="init", reverse_x=False)
        close_matplotlib_figures()
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
    thresh: Optional[str],
) -> Tuple[Any, Path]:
    """
    Optimize an endpoint geometry using LBFGS/RFO with settings mirroring path_search defaults.

    Args:
        geom: pysisyphus Geometry with calculator attached.
        opt_mode_default: "lbfgs"/"light" or "rfo"/"heavy".
        out_dir: base directory for the optimization outputs.
        tag: tag prefix for the subdirectory.
        dump: whether to dump optimizer trajectory.
        thresh: optional convergence preset to override defaults.

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
    ensure_dir(opt_dir)
    cfg["out_dir"] = str(opt_dir)
    cfg["dump"] = bool(dump)
    max_cycles = int(cfg.get("max_cycles", 300))
    cfg["max_cycles"] = max_cycles

    geom.set_calculator(getattr(geom, "calculator", None))

    if thresh is not None:
        cfg["thresh"] = str(thresh)

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
    ref_pdb: Optional[Path],
    convert_files: bool,
    overrides: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Run freq CLI; return parsed thermo dict (may be empty).
    """
    fdir = out_dir
    ensure_dir(fdir)
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
        "True"
        if freeze_use
        and (
            pdb_path.suffix.lower() == ".pdb"
            or (ref_pdb is not None and ref_pdb.suffix.lower() == ".pdb")
        )
        else "False",
        "--convert-files",
        "True" if convert_files else "False",
        "--out-dir",
        str(fdir),
    ]
    if ref_pdb is not None:
        args.extend(["--ref-pdb", str(ref_pdb)])

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
    _run_cli_main("freq", _freq_cli.cli, args, on_nonzero="warn", prefix="freq")
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _read_imaginary_frequency(freq_dir: Path) -> Optional[Dict[str, Any]]:
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




def _run_dft_for_state(
    pdb_path: Path,
    q_int: int,
    spin: int,
    out_dir: Path,
    args_yaml: Optional[Path],
    func_basis: str = "wb97m-v/def2-tzvpd",
    overrides: Optional[Dict[str, Any]] = None,
    engine: str = "gpu",
    ref_pdb: Optional[Path] = None,
    convert_files: bool = True,
) -> Dict[str, Any]:
    """
    Run dft CLI; return parsed result.yaml dict (may be empty).
    """
    ddir = out_dir
    ensure_dir(ddir)
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
        "--convert-files",
        "True" if convert_files else "False",
        "--out-dir",
        str(ddir),
    ]
    if ref_pdb is not None:
        args.extend(["--ref-pdb", str(ref_pdb)])
    if engine:
        args.extend(["--engine", str(engine)])

    _append_cli_arg(args, "--max-cycle", overrides.get("max_cycle"))
    _append_cli_arg(args, "--conv-tol", overrides.get("conv_tol"))
    _append_cli_arg(args, "--grid-level", overrides.get("grid_level"))

    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _run_cli_main("dft", _dft_cli.cli, args, on_nonzero="warn", prefix="dft")
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
    ref_pdb: Optional[Path],
    convert_files: bool,
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
            ref_pdb=ref_pdb,
            convert_files=convert_files,
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
    opt_mode_default: Optional[str],
    ref_pdb: Optional[Path],
    convert_files: bool,
    overrides: Optional[Dict[str, Any]] = None,
) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI pocket structure; return (final_geom_path, ts_geom).

    Prefer the XYZ output to preserve coordinate precision between workflow steps, while still writing
    PDB/GJF companions when requested by the original input type.
    """
    overrides = overrides or {}
    prepared_input = prepare_input_structure(hei_pdb)
    apply_ref_pdb_override(prepared_input, ref_pdb)
    needs_pdb = prepared_input.source_path.suffix.lower() == ".pdb"
    needs_gjf = prepared_input.is_gjf
    ref_pdb = prepared_input.source_path if needs_pdb else None
    ts_dir = _resolve_override_dir(out_dir / "ts", overrides.get("out_dir"))
    ensure_dir(ts_dir)

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
        "--convert-files",
        "True" if convert_files else "False",
        "--out-dir",
        str(ts_dir),
    ]

    if opt_mode is not None:
        ts_args.extend(["--opt-mode", str(opt_mode)])

    _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
    _append_cli_arg(ts_args, "--dump", overrides.get("dump"))
    _append_cli_arg(ts_args, "--thresh", overrides.get("thresh"))
    _append_cli_arg(ts_args, "--flatten-imag-mode", overrides.get("flatten_imag_mode"))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

    if args_yaml is not None:
        ts_args.extend(["--args-yaml", str(args_yaml)])
    if ref_pdb is not None:
        ts_args.extend(["--ref-pdb", str(ref_pdb)])

    click.echo(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
    _run_cli_main("tsopt", _tsopt.cli, ts_args, on_nonzero="raise", prefix="tsopt")

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

    if ts_xyz.exists():
        ts_geom_path = ts_xyz
    elif needs_pdb and ts_pdb.exists():
        ts_geom_path = ts_pdb
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
    ref_pdb_template: Optional[Path],
    g_ts: Any,
    q_int: int,
    spin: int,
    freeze_links_flag: bool,
    calc_cfg: Dict[str, Any],
    args_yaml: Optional[Path],
    convert_files: bool,
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
    ensure_dir(irc_dir)

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
        "True"
        if freeze_links_flag
        and (
            ref_pdb_for_seg.suffix.lower() == ".pdb"
            or (ref_pdb_template is not None and ref_pdb_template.suffix.lower() == ".pdb")
        )
        else "False",
        "--convert-files",
        "True" if convert_files else "False",
    ]
    if ref_pdb_template is not None:
        irc_args.extend(["--ref-pdb", str(ref_pdb_template)])

    if args_yaml is not None:
        irc_args.extend(["--args-yaml", str(args_yaml)])
    click.echo(f"[irc] Running EulerPC IRC → out={irc_dir}")
    _run_cli_main("irc", _irc_cli.cli, irc_args, on_nonzero="raise", prefix="irc")

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
                _path_search._convert_to_pdb_logged(finished_trj, ref_pdb_path=ref_for_conv, out_path=finished_pdb)
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to convert finished_irc.trj to PDB: {e}", err=True)

    elems, c_first, c_last = read_xyz_first_last(finished_trj)

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
                        changed, _ = _path_search.has_bond_change(x, y, bond_cfg)
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
                        d_LL = _path_search._rmsd_between(g_left, gL_end)
                        d_LR = _path_search._rmsd_between(g_left, gR_end)
                        d_RL = _path_search._rmsd_between(g_right, gL_end)
                        d_RR = _path_search._rmsd_between(g_right, gR_end)
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
                    f"[irc] WARNING: LBFGS endpoints not found for segment tag '{seg_tag}' in the segment directories; "
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
            close_matplotlib_figures()
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
        "the input may be PDB/XYZ/GJF (integer indices only for non-PDB inputs). You may pass a single '-i' "
        "followed by multiple space-separated files "
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
    "--selected-resn",
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
        "Total charge (number) or per-resname mapping like 'GPP:-3,MMT:-1'. "
        "Used for extractor charge summaries; when extraction is skipped, PDB inputs "
        "derive the total charge and numeric values act as a total-charge fallback."
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
    "--opt-mode-post",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default=None,
    show_default=False,
    help=(
        "Optimizer mode override for TSOPT/post-IRC endpoint optimizations. "
        "If unset, uses --opt-mode when explicitly provided; otherwise falls back to tsopt defaults."
    ),
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help=(
        "Dump GSM/MEP trajectories. Always forwarded to path_search/path-opt; "
        "scan/tsopt receive it only when explicitly set here. "
        "The freq stage uses dump=True by default; set --dump False explicitly to disable it."
    ),
)
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
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
    help=(
        "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
    ),
)
@click.option(
    "--thresh-post",
    type=str,
    default="baker",
    show_default=True,
    help=(
        "Convergence preset for post-IRC endpoint optimizations "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "YAML forwarded unchanged to downstream subcommands (path_search/path-opt/scan/tsopt/freq/dft); "
        "include the sections those commands accept (e.g., geom, calc, gs, opt, sopt, bond, search, dmf, bias, freq, thermo, dft)."
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
    help="Override tsopt --max-cycles value. Defaults to 10000 when not provided.",
)
@click.option(
    "--tsopt-out-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help="Override tsopt output subdirectory (relative paths are resolved against the default).",
)
@click.option(
    "--flatten-imag-mode",
    "flatten_imag_mode",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Enable the extra-imaginary-mode flattening loop in tsopt (light: dimer loop, heavy: post-RSIRFO); False forces flatten_max_iter=0.",
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
    help="Override freq --max-write value. Defaults to 10.",
)
@click.option(
    "--freq-amplitude-ang",
    type=float,
    default=None,
    help="Override freq --amplitude-ang (Å). Defaults to 0.8.",
)
@click.option(
    "--freq-n-frames",
    type=int,
    default=None,
    help="Override freq --n-frames value. Defaults to 20.",
)
@click.option(
    "--freq-sort",
    type=click.Choice(["value", "abs"], case_sensitive=False),
    default=None,
    help="Override freq mode sorting. Defaults to 'value'.",
)
@click.option(
    "--freq-temperature",
    type=float,
    default=None,
    help="Override freq thermochemistry temperature (K). Defaults to 298.15 K.",
)
@click.option(
    "--freq-pressure",
    type=float,
    default=None,
    help="Override freq thermochemistry pressure (atm). Defaults to 1.0 atm.",
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
    help="Override dft --func-basis value. Defaults to 'wb97m-v/def2-tzvpd'.",
)
@click.option(
    "--dft-max-cycle",
    type=int,
    default=None,
    help="Override dft --max-cycle value. Defaults to 100.",
)
@click.option(
    "--dft-conv-tol",
    type=float,
    default=None,
    help="Override dft --conv-tol value. Defaults to 1e-9.",
)
@click.option(
    "--dft-grid-level",
    type=int,
    default=None,
    help="Override dft --grid-level value. Defaults to 3.",
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
    "--scan-list",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help=(
        "Python-like list of (i,j,target_Å) per stage for **single-structure** scan. A single "
        "literal runs one stage; multiple literals run **sequentially**, each starting from the "
        "prior stage's relaxed structure. "
        "Example: --scan-lists '[(12,45,1.35)]' '[(10,55,2.20),(23,34,1.80)]'. "
        "Pass a single --scan-list(s) followed by multiple values to define multiple stages "
        "(repeated flags are not accepted). "
        "Indices refer to the original full input PDB (1-based). When extraction is used, they are "
        "auto-mapped to the pocket after extraction. For non-PDB single-structure scans, only integer "
        "indices are supported (1-based by default). Stage results feed into the MEP step (path_search or path_opt)."
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
        "Override the scan subcommand indexing interpretation (True = 1-based, False = 0-based). "
        "Defaults to 1-based."
    ),
)
@click.option(
    "--scan-max-step-size",
    type=float,
    default=None,
    help="Override scan --max-step-size (Å). Defaults to 0.20 Å.",
)
@click.option(
    "--scan-bias-k",
    type=float,
    default=None,
    help="Override scan harmonic bias strength k (eV/Å^2). Defaults to 100.",
)
@click.option(
    "--scan-relax-max-cycles",
    type=int,
    default=None,
    help="Override scan relaxation max cycles per step. Defaults to 10000.",
)
@click.option(
    "--scan-preopt",
    "scan_preopt_override",
    type=click.BOOL,
    default=None,
    help="Override scan --preopt flag. When omitted, this follows --preopt (default True).",
)
@click.option(
    "--scan-endopt",
    "scan_endopt_override",
    type=click.BOOL,
    default=None,
    help="Override scan --endopt flag. Defaults to True.",
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
    workers: int,
    workers_per_node: int,
    verbose: bool,
    spin: int,
    freeze_links_flag: bool,
    mep_mode: str,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    opt_mode_post: Optional[str],
    dump: bool,
    convert_files: bool,
    refine_path: bool,
    thresh: Optional[str],
    thresh_post: str,
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
    flatten_imag_mode: bool,
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

    With ``--refine-path True`` (default), the recursive ``path_search`` workflow is used. When ``False``,
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

    opt_mode_set = False
    opt_mode_post_set = False
    try:
        opt_mode_source = ctx.get_parameter_source("opt_mode")
        opt_mode_set = opt_mode_source not in (None, ParameterSource.DEFAULT)
        opt_mode_post_source = ctx.get_parameter_source("opt_mode_post")
        opt_mode_post_set = opt_mode_post_source not in (None, ParameterSource.DEFAULT)
    except Exception:
        opt_mode_set = False
        opt_mode_post_set = False

    argv_all = sys.argv[1:]
    i_vals = collect_option_values(argv_all, ("-i", "--input"))
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

    scan_vals = collect_single_option_values(argv_all, ("--scan-lists", "--scan-list"), "--scan-list(s)")
    if scan_vals:
        scan_lists_raw = tuple(scan_vals)

    is_single = len(input_paths) == 1
    has_scan = bool(scan_lists_raw)
    single_tsopt_mode = is_single and (not has_scan) and do_tsopt

    if (len(input_paths) < 2) and not (is_single and (has_scan or do_tsopt)):
        raise click.BadParameter(
            "Provide at least two structures with -i/--input in reaction order, "
            "or use a single structure with --scan-lists, or a single structure with --tsopt True."
        )

    tsopt_opt_mode_default: Optional[str] = None
    if opt_mode_post_set and opt_mode_post is not None:
        tsopt_opt_mode_default = opt_mode_post.lower()
    elif opt_mode_set:
        tsopt_opt_mode_default = opt_mode.lower()
    else:
        tsopt_opt_mode_default = "heavy"
    tsopt_overrides: Dict[str, Any] = {}
    if tsopt_max_cycles is not None:
        tsopt_overrides["max_cycles"] = int(tsopt_max_cycles)
    if dump_override_requested:
        tsopt_overrides["dump"] = bool(dump)
    if tsopt_out_dir is not None:
        tsopt_overrides["out_dir"] = tsopt_out_dir
    if hessian_calc_mode is not None:
        tsopt_overrides["hessian_calc_mode"] = hessian_calc_mode
    if thresh_post is not None:
        tsopt_overrides["thresh"] = str(thresh_post)
    tsopt_overrides["flatten_imag_mode"] = bool(flatten_imag_mode)

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
    _set_yaml_freeze_atoms(yaml_cfg)

    skip_extract = center_spec is None or str(center_spec).strip() == ""

    out_dir = out_dir.resolve()
    pockets_dir = out_dir / "pockets"
    path_dir = out_dir / ("path_search" if refine_path else "path_opt")
    scan_dir = _resolve_override_dir(out_dir / "scan", scan_out_dir)
    stage_total = 3
    ensure_dir(out_dir)
    if not skip_extract:
        ensure_dir(pockets_dir)
    if not single_tsopt_mode:
        ensure_dir(path_dir)

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
            ensure_dir(elem_tmp_dir)
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

    resolved_charge: Optional[int] = None

    if not skip_extract:
        click.echo(
            f"=== [all] Stage 1/{stage_total} — Active-site pocket extraction (multi-structure union when applicable) ===\n"
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
            resolved_charge = _round_charge_with_note(q_total, prefix="[all]")
        except Exception as e:
            raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")
    else:
        click.echo(
            f"=== [all] Stage 1/{stage_total} — Extraction skipped (no -c/--center); using FULL structures as pockets ===\n"
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

        charge_total: float
        ligand_charge_numeric: Optional[float] = None
        if ligand_charge is not None:
            try:
                ligand_charge_numeric = float(ligand_charge)
            except Exception:
                ligand_charge_numeric = None

            if first_input.suffix.lower() == ".pdb":
                try:
                    with prepare_input_structure(first_input) as prepared:
                        resolved_charge = _derive_charge_from_ligand_charge(
                            prepared, ligand_charge, prefix="[all]"
                        )
                except Exception as e:
                    click.echo(
                        f"[all] NOTE: failed to derive total charge from full complex: {e}; "
                        "falling back to legacy handling.",
                        err=True,
                    )
            else:
                click.echo(
                    "[all] NOTE: --ligand-charge derivation requires a PDB input; skipping full-complex derivation.",
                    err=True,
                )

        if resolved_charge is None:
            if ligand_charge_numeric is not None:
                charge_total = ligand_charge_numeric
                click.echo(
                    f"[all] Using --ligand-charge as TOTAL system charge: {charge_total:+g}"
                )
                resolved_charge = _round_charge_with_note(charge_total, prefix="[all]")
            elif gjf_charge is not None:
                charge_total = float(gjf_charge)
                click.echo(f"[all] Using total charge from first GJF: {charge_total:+g}")
                resolved_charge = _round_charge_with_note(charge_total, prefix="[all]")
            else:
                if charge_override is None:
                    raise click.ClickException(
                        "[all] Total charge could not be resolved. Provide -q/--charge, "
                        "--ligand-charge, or a .gjf input with charge metadata."
                    )

        if (not user_provided_spin) and (gjf_spin is not None):
            spin = int(gjf_spin)
            click.echo(f"[all] Spin multiplicity set from GJF: {spin}")

    if charge_override is not None:
        q_int = int(charge_override)
        override_msg = (
            f"[all] WARNING: -q/--charge override supplied; forcing TOTAL system charge to {q_int:+d}"
        )
        if resolved_charge is not None:
            override_msg += f" (would otherwise use {int(resolved_charge):+d} from workflow)"
        click.echo(override_msg, err=True)
    else:
        q_int = int(resolved_charge) if resolved_charge is not None else 0

    freeze_ref: Optional[Path] = None
    if freeze_links_flag:
        if not skip_extract and pocket_outputs:
            freeze_ref = pocket_outputs[0]
        else:
            for p in input_paths:
                if p.suffix.lower() == ".pdb":
                    freeze_ref = p.resolve()
                    break
        if freeze_ref is not None:
            _get_freeze_atoms(freeze_ref, freeze_links_flag)

    args_yaml = _write_args_yaml_with_freeze_atoms(
        args_yaml,
        _freeze_atoms_for_log(),
    )

    calc_cfg_shared = _build_calc_cfg(
        q_int,
        spin,
        workers=workers,
        workers_per_node=workers_per_node,
        yaml_cfg=yaml_cfg,
    )

    # -------------------------------------------------------------------------
    # TSOPT-only single-structure mode
    # -------------------------------------------------------------------------
    if single_tsopt_mode:
        click.echo("\n=== [all] TSOPT-only single-structure mode ===\n")
        tsroot = out_dir / "tsopt_single"
        ensure_dir(tsroot)

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
            ts_initial_pdb if ts_initial_pdb.suffix.lower() == ".pdb" else None,
            convert_files,
            overrides=tsopt_overrides,
        )

        irc_res = _irc_and_match(
            seg_idx=1,
            seg_dir=tsroot,
            ref_pdb_for_seg=ts_pdb,
            seg_pocket_pdb=ts_initial_pdb,
            ref_pdb_template=ts_initial_pdb if ts_initial_pdb.suffix.lower() == ".pdb" else None,
            g_ts=g_ts,
            q_int=q_int,
            spin=spin,
            freeze_links_flag=freeze_links_flag,
            calc_cfg=calc_cfg_shared,
            args_yaml=args_yaml,
            convert_files=convert_files,
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
        ensure_dir(struct_dir)
        pocket_ref = ts_initial_pdb
        pR_irc = _save_single_geom_as_pdb_for_tools(
            g_react_irc, pocket_ref, struct_dir, "reactant_irc"
        )
        pP_irc = _save_single_geom_as_pdb_for_tools(
            g_prod_irc, pocket_ref, struct_dir, "product_irc"
        )
        pT = _save_single_geom_as_pdb_for_tools(gT, pocket_ref, struct_dir, "ts")

        endpoint_opt_dir = tsroot / "endpoint_opt"
        ensure_dir(endpoint_opt_dir)
        try:
            g_react_opt, _ = _optimize_endpoint_geom(
                g_react_irc,
                tsopt_opt_mode_default,
                endpoint_opt_dir,
                "reactant",
                dump=dump,
                thresh=thresh_post,
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
                thresh=thresh_post,
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
            ref_pdb_for_tsopt_only = (
                ts_initial_pdb if ts_initial_pdb.suffix.lower() == ".pdb" else None
            )
            tR = _run_freq_for_state(
                pR,
                q_int,
                spin,
                freq_root / "R",
                args_yaml,
                freeze_links_flag,
                ref_pdb_for_tsopt_only,
                convert_files,
                overrides=freq_overrides,
            )
            tT = _run_freq_for_state(
                pT,
                q_int,
                spin,
                freq_root / "TS",
                args_yaml,
                freeze_links_flag,
                ref_pdb_for_tsopt_only,
                convert_files,
                overrides=freq_overrides,
            )
            tP = _run_freq_for_state(
                pP,
                q_int,
                spin,
                freq_root / "P",
                args_yaml,
                freeze_links_flag,
                ref_pdb_for_tsopt_only,
                convert_files,
                overrides=freq_overrides,
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
                ref_pdb_for_tsopt_only,
                convert_files,
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
            changed, bond_summary = _path_search.has_bond_change(
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
                n_images = len(read_xyz_as_blocks(irc_trj_path))
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
            _copy_logged(tsroot / "summary.yaml", out_dir / "summary.yaml", label="summary.yaml")
            try:
                ts_freq_info = (
                    _read_imaginary_frequency(freq_root / "TS") if do_thermo else None
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
                _copy_logged(tsroot / "summary.log", out_dir / "summary.log", label="summary.log", echo=False)
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
                    _copy_logged(src, dst, label=src.name)
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to mirror *_all diagrams in TSOPT-only mode: {e}",
                err=True,
            )

        try:
            if isinstance(irc_plot_path, Path) and irc_plot_path.exists():
                dst = out_dir / "irc_plot_all.png"
                _copy_logged(irc_plot_path, dst, label="irc_plot_all.png")
        except Exception as e:
            click.echo(
                f"[all] WARNING: Failed to mirror IRC plot in TSOPT-only mode: {e}",
                err=True,
            )

        try:
            ts_pdb = pT.with_suffix(".pdb")
            if ts_pdb.exists():
                ts_copy = out_dir / "ts_seg_01.pdb"
                shutil.copy2(ts_pdb, ts_copy)
            else:
                ts_copy = out_dir / "ts_seg_01.xyz"
                write_xyz_trj_with_energy([gT], [float(gT.energy)], ts_copy)
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
    pocket_ref_pdbs: Optional[List[Path]] = None
    if is_single and has_scan:
        click.echo("=== [all] Stage 1b — Staged scan on input ===\n")
        ensure_dir(scan_dir)

        if skip_extract:
            scan_input_pdb = Path(input_paths[0]).resolve()
            scan_atom_meta = None
            if scan_input_pdb.suffix.lower() == ".pdb":
                scan_atom_meta = load_pdb_atom_metadata(scan_input_pdb)
            converted_scan_stages = _parse_scan_lists_literals(
                scan_lists_raw, atom_meta=scan_atom_meta
            )
        else:
            scan_input_pdb = Path(pocket_outputs[0]).resolve()
            full_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _convert_scan_lists_to_pocket_indices(
                scan_lists_raw, full_input_pdb, scan_input_pdb
            )
            click.echo(
                "[all] Remapped --scan-lists indices from the full PDB to the pocket ordering."
            )

        # For extracted pockets, scan stages are already 1-based and can be rebased
        # to 0-based if the user explicitly requests it. For non-extracted inputs,
        # keep indices as provided and let scan.py interpret them via --one-based.
        if skip_extract:
            scan_one_based_effective = True
        else:
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
            "--convert-files",
            "True" if convert_files else "False",
            "--preopt",
            "True" if scan_preopt_use else "False",
            "--endopt",
            "True" if scan_endopt_use else "False",
            "--opt-mode",
            str(scan_opt_mode_use),
        ]
        if scan_input_pdb.suffix.lower() == ".pdb":
            scan_args.extend(["--ref-pdb", str(scan_input_pdb)])

        if dump_override_requested:
            scan_args.extend(
                ["--dump", "True" if dump else "False"]
            )

        if scan_one_based is not None:
            scan_args.extend(["--one-based", "True" if scan_one_based else "False"])

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        if args_yaml is not None:
            scan_args.extend(["--args-yaml", str(args_yaml)])
        if scan_stage_literals:
            scan_args.append("--scan-lists")
            scan_args.extend(scan_stage_literals)

        click.echo("[all] Invoking scan with arguments:")
        click.echo("  " + " ".join(scan_args))

        _run_cli_main("scan", _scan_cli.cli, scan_args, on_nonzero="raise", prefix="all")

        stage_results: List[Path] = []
        for st in sorted(scan_dir.glob("stage_*")):
            res = _find_with_suffixes(st / "result", [".xyz", ".pdb", ".gjf"])
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

        if scan_input_pdb.suffix.lower() == ".pdb":
            candidate_pdbs: List[Path] = [scan_input_pdb]
            missing_pdb = False
            for stage_path in stage_results:
                if stage_path.suffix.lower() == ".pdb":
                    candidate_pdbs.append(stage_path)
                else:
                    pdb_candidate = stage_path.with_suffix(".pdb")
                    if pdb_candidate.exists():
                        candidate_pdbs.append(pdb_candidate)
                    else:
                        missing_pdb = True
                        break
            if not missing_pdb:
                pocket_ref_pdbs = candidate_pdbs
            else:
                click.echo(
                    "[all] WARNING: pocket PDB snapshots for staged scan were not found; "
                    "full-system merge will use input paths instead.",
                    err=True,
                )
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
            "[all] NOTE: skipping --ref-full-pdb (no --center; inputs already represent full structures)."
        )
    elif is_single and has_scan:
        if _is_pdb(input_paths[0]):
            gave_ref_pdb = True
        else:
            click.echo(
                "[all] NOTE: skipping --ref-full-pdb (single+scan: original input is not a PDB)."
            )
    else:
        if all(_is_pdb(p) for p in input_paths):
            gave_ref_pdb = True
        else:
            click.echo(
                "[all] NOTE: skipping --ref-full-pdb (one or more original inputs are not PDB)."
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
                "--convert-files",
                "True" if convert_files else "False",
                "--out-dir",
                str(seg_dir),
                "--preopt",
                "True" if preopt else "False",
            ]
            ref_pdb_for_seg: Optional[Path] = None
            if pocket_ref_pdbs and len(pocket_ref_pdbs) >= idx:
                ref_pdb_for_seg = pocket_ref_pdbs[idx - 1]
            elif pL.suffix.lower() == ".pdb":
                ref_pdb_for_seg = pL
            elif pR.suffix.lower() == ".pdb":
                ref_pdb_for_seg = pR
            elif is_single and has_scan and input_paths[0].suffix.lower() == ".pdb":
                ref_pdb_for_seg = input_paths[0]
            if ref_pdb_for_seg is not None:
                po_args.extend(["--ref-pdb", str(ref_pdb_for_seg)])
            if thresh is not None:
                po_args.extend(["--thresh", str(thresh)])
            if args_yaml is not None:
                po_args.extend(["--args-yaml", str(args_yaml)])

            click.echo(f"[all] Invoking path-opt for segment {idx}:")
            click.echo("  " + " ".join(po_args))

            _run_cli_main(
                "path-opt",
                _path_opt.cli,
                po_args,
                on_nonzero="raise",
                prefix=f"all seg {idx:02d}",
            )

            seg_trj = seg_dir / "final_geometries.trj"
            if not seg_trj.exists():
                raise click.ClickException(
                    f"[all] path-opt segment {idx} did not produce final_geometries.trj"
                )

            try:
                mirror_dir = path_dir / f"{seg_tag}_mep"
                mirror_trj = mirror_dir / "final_geometries.trj"

                ensure_dir(mirror_dir)
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
                    _path_search._convert_to_pdb_logged(
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

            raw_blocks = read_xyz_as_blocks(seg_trj, strict=True)
            blocks = ["\n".join(b) + "\n" for b in raw_blocks]
            if not blocks:
                raise click.ClickException(
                    f"[all] No frames read from path-opt segment {idx} trajectory: {seg_trj}"
                )
            if idx > 1:
                blocks = blocks[1:]
            combined_blocks.extend(blocks)

            energies_seg: List[float] = []
            for blk in raw_blocks:
                E = np.nan
                if len(blk) >= 2:
                    try:
                        E = float(blk[1].split()[0])
                    except Exception:
                        E = np.nan
                energies_seg.append(E)

            first_last = None
            try:
                first_last = xyz_blocks_first_last(raw_blocks, path=seg_trj)
            except Exception as e:
                click.echo(
                    f"[all] WARNING: failed to parse first/last frames for segment {idx:02d}: {e}",
                    err=True,
                )

            path_opt_segments.append(
                {
                    "tag": seg_tag,
                    "energies": energies_seg,
                    "traj": seg_trj,
                    "inputs": (pL, pR),
                    "first_last": first_last,
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
            close_matplotlib_figures()
            click.echo(f"[plot] Saved energy plot → '{path_dir / 'mep_plot.png'}'")
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot concatenated MEP: {e}", err=True)

        try:
            if pockets_for_path[0].suffix.lower() == ".pdb":
                mep_pdb = _path_search._convert_to_pdb_logged(
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
                first_last = info.get("first_last")
                if first_last:
                    elems, c_first, c_last = first_last
                else:
                    elems, c_first, c_last = read_xyz_first_last(Path(info["traj"]))
                freeze_atoms = _get_freeze_atoms(info["inputs"][0], freeze_links_flag)
                gL = _geom_from_angstrom(elems, c_first, freeze_atoms)
                gR = _geom_from_angstrom(elems, c_last, freeze_atoms)
                changed, bond_summary = _path_search.has_bond_change(gL, gR, bond_cfg)
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
            "n_images": len(read_xyz_as_blocks(final_trj)),
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
                    _copy_logged(src, dst, label=name)

            for stem in ("mep", "mep_w_ref"):
                for ext in (".trj", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        dst = out_dir / src.name
                        _copy_logged(src, dst, label=src.name)
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
            _copy_logged(path_dir / "summary.log", out_dir / "summary.log", label="summary.log")
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
        ps_args.extend(["--convert-files", "True" if convert_files else "False"])
        if args_yaml is not None:
            ps_args.extend(["--args-yaml", str(args_yaml)])

        if gave_ref_pdb:
            for p in (input_paths if not (is_single and has_scan) else (input_paths[:1] * len(pockets_for_path))):
                ps_args.extend(["--ref-full-pdb", str(p)])
            if pocket_ref_pdbs:
                for p in pocket_ref_pdbs:
                    ps_args.extend(["--ref-pdb", str(p)])

        click.echo("[all] Invoking path_search with arguments:")
        click.echo("  " + " ".join(ps_args))

        _run_cli_main(
            "path_search",
            _path_search.cli,
            ps_args,
            on_nonzero="raise",
            prefix="all",
        )

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
                    _copy_logged(src, dst, label=name)

            for stem in ("mep", "mep_w_ref"):
                for ext in (".trj", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        dst = out_dir / src.name
                        _copy_logged(src, dst, label=src.name)
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
            "[all] --ref-full-pdb was not provided; full-system merged trajectories are not produced."
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
        ensure_dir(seg_dir)

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
        hei_pocket_path = _find_with_suffixes(hei_base, [".xyz", ".pdb", ".gjf"])
        if hei_pocket_path is None:
            click.echo(
                f"[post] WARNING: HEI pocket file not found for segment {seg_idx:02d} (searched .pdb/.xyz/.gjf); skipping TSOPT.",
                err=True,
            )
            continue
        ref_pdb_for_seg: Optional[Path] = None
        if hei_pocket_path.suffix.lower() == ".pdb":
            ref_pdb_for_seg = hei_pocket_path
        else:
            candidate_ref = hei_base.with_suffix(".pdb")
            if candidate_ref.exists():
                ref_pdb_for_seg = candidate_ref

        struct_dir = seg_dir / "structures"
        ensure_dir(struct_dir)
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
                ref_pdb_for_seg,
                convert_files,
                overrides=tsopt_overrides,
            )

            irc_res = _irc_and_match(
                seg_idx=seg_idx,
                seg_dir=seg_dir,
                ref_pdb_for_seg=ts_pdb,
                seg_pocket_pdb=hei_pocket_path,
                ref_pdb_template=ref_pdb_for_seg,
                g_ts=g_ts,
                q_int=q_int,
                spin=spin,
                freeze_links_flag=freeze_links_flag,
                calc_cfg=calc_cfg_shared,
                args_yaml=args_yaml,
                convert_files=convert_files,
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

            ref_struct_template = ref_pdb_for_seg or hei_pocket_path
            pL_irc = _save_single_geom_as_pdb_for_tools(
                gL, ref_struct_template, struct_dir, "reactant_irc"
            )
            pT = _save_single_geom_as_pdb_for_tools(
                gT, ref_struct_template, struct_dir, "ts"
            )
            pR_irc = _save_single_geom_as_pdb_for_tools(
                gR, ref_struct_template, struct_dir, "product_irc"
            )

            endpoint_opt_dir = seg_dir / "endpoint_opt"
            ensure_dir(endpoint_opt_dir)
            try:
                g_react_opt, _ = _optimize_endpoint_geom(
                    gL,
                    tsopt_opt_mode_default,
                    endpoint_opt_dir,
                    f"seg_{seg_idx:02d}_reactant",
                    dump=dump,
                    thresh=thresh_post,
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
                    thresh=thresh_post,
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
                g_react_opt, ref_struct_template, struct_dir, "reactant"
            )
            pR = _save_single_geom_as_pdb_for_tools(
                g_prod_opt, ref_struct_template, struct_dir, "product"
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
                ts_pdb = pT.with_suffix(".pdb")
                if ts_pdb.exists():
                    ts_copy = out_dir / f"ts_seg_{seg_idx:02d}.pdb"
                    shutil.copy2(ts_pdb, ts_copy)
                else:
                    ts_copy = out_dir / f"ts_seg_{seg_idx:02d}.xyz"
                    write_xyz_trj_with_energy([gT], [float(gT.energy)], ts_copy)
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
                seg_root / f"mep_seg_{seg_idx:02d}", [".xyz", ".pdb"]
            )

            # Decide reference PDB (if any) for freeze-atoms detection / PDB conversion
            freeze_ref: Optional[Path] = ref_pdb_for_seg
            if freeze_ref is None and seg_pocket_path is not None and seg_pocket_path.suffix.lower() == ".pdb":
                freeze_ref = seg_pocket_path
            elif freeze_ref is None and hei_pocket_path.suffix.lower() == ".pdb":
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

            ref_for_structs = ref_pdb_for_seg or (seg_pocket_path if seg_pocket_path is not None else hei_pocket_path)
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
                    ref_pdb_for_seg,
                    convert_files,
                    overrides=freq_overrides,
                )
                tT = _run_freq_for_state(
                    p_ts,
                    q_int,
                    spin,
                    freq_seg_root / "TS",
                    args_yaml,
                    freeze_links_flag,
                    ref_pdb_for_seg,
                    convert_files,
                    overrides=freq_overrides,
                )
                tP = _run_freq_for_state(
                    p_prod,
                    q_int,
                    spin,
                    freq_seg_root / "P",
                    args_yaml,
                    freeze_links_flag,
                    ref_pdb_for_seg,
                    convert_files,
                    overrides=freq_overrides,
                )
                thermo_payloads = {"R": tR, "TS": tT, "P": tP}
                ts_freq_info = _read_imaginary_frequency(freq_seg_root / "TS")
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
                    ref_pdb_for_seg,
                    convert_files,
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
