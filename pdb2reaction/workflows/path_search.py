"""Recursive GSM/DMF paths with bond-change detection and refinement."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Callable

import logging
import shlex
import sys
import textwrap
import tempfile
import os
import time
import re

import click
from pdb2reaction.core.output import emit
import numpy as np

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import AU2KCALPERMOL, BOHR2ANG

# Biopython for PDB I/O and parsing
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO

from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    THRESH_CHOICES,
    UMA_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    OPT_MODE_ALIASES,
    DMF_KW,
    fresh_dmf_config,
    GS_KW,
    STOPT_KW,
    BOND_KW,
    SEARCH_KW,
    OUT_DIR_PATH_SEARCH,
    apply_backend_defaults,
)
from pdb2reaction.workflows.path_opt import (
    _validate_dmf_solvent_compatibility,
    _optimize_single,
    _run_dmf_mep,
    resolve_dmf_solve_tol,
)
from pdb2reaction.workflows._path_yaml_helpers import apply_single_opt_yaml_layer
from pdb2reaction.core.utils import (
    as_list,
    collect_option_values,
    deep_update,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    emit_dry_run_complete,
    format_elapsed,
    build_energy_diagram,
    prepare_input_structure,
    fill_charge_spin_from_gjf_inputs,
    _derive_charge_from_ligand_charge,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    _convert_to_pdb_logged,
    PreparedInputStructure,
    validate_charge_spin_for_prepared,
    GjfTemplate,
    geom_from_xyz_string,
    close_matplotlib_figures,
    write_xyz_trj_with_energy,
    set_freeze_atoms_or_warn,
    YamlLiteralStr,
    load_prepared_geometries,
    normalize_choice,
    cli_param_overridden,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    lossless_int,
    optional_positive_int,
)
from pdb2reaction.core.result_commit import apply_current_run_id, commit_json
from pdb2reaction.io.summary import (
    emit_method_citations,
    method_references,
    write_summary_log,
)
from pdb2reaction.io.trj2fig import run_trj2fig
from pdb2reaction.io.path_mode_cache import write_path_mode_cache
from pdb2reaction.io.structure_formats import (
    attach_template_metadata,
    atom_site_from_biopython_atom,
    coordinate_template_for,
    register_output_template_and_write_cif,
    residue_auth_identity,
)
from pdb2reaction.domain.bond_changes import has_bond_change
from pdb2reaction.workflows.align_freeze import (
    align_and_refine_sequence_inplace,
    alignment_failed_pair_indices,
    kabsch_R_t,
)
from pdb2reaction.cli.common_options import (
    add_allow_charge_mult_mismatch_option,
    add_backend_model_option,
    add_calc_file_option,
    add_coord_type_option,
    add_deterministic_option,
    add_precision_option,
)
from pdb2reaction.cli.decorators import run_cli, resolve_yaml_sources, load_merged_yaml_cfg

logger = logging.getLogger(__name__)

def _bond_changes_block(text: Optional[str]):
    """
    Prepare bond-change summaries for YAML output, emitting structured lists instead
    of escaped ``\n`` strings when possible.
    """

    if text is None:
        return ""

    cleaned = str(text).strip()
    if not cleaned:
        return ""

    # Convert the ``summarize_changes`` output
    #   Bond formed (1):
    #     - C1-C2 : 3.0 Å --> 1.5 Å
    #   Bond broken: None
    # into a YAML-friendly nested list to avoid embedded newlines.
    lines = cleaned.splitlines()
    sections = []
    title: Optional[str] = None
    entries: List[str] = []

    def _flush() -> None:
        nonlocal title, entries
        if title is None:
            return
        payload = entries if entries else ["None"]
        sections.append({title: payload})
        title, entries = None, []

    for ln in lines:
        stripped = ln.strip()
        if stripped.startswith("Bond "):
            _flush()
            if stripped.endswith(": None"):
                sections.append({stripped[:-len(": None")]: ["None"]})
            else:
                title = stripped.rstrip(":")
        elif stripped.startswith("- "):
            entries.append(stripped[2:])

    _flush()

    # Fallback to literal block if the format is unexpected
    if sections:
        return sections
    if "\n" in cleaned:
        return YamlLiteralStr(cleaned)
    return cleaned




def _calc_rmsd(A: np.ndarray, B: np.ndarray, indices: Optional[Sequence[int]] = None) -> float:
    """RMSD between A and B (no rigid alignment). Optional subset selection via `indices`."""
    assert A.shape == B.shape and A.shape[1] == 3
    if indices is not None and len(indices) > 0:
        idx = np.array(sorted({int(i) for i in indices if 0 <= int(i) < A.shape[0]}), dtype=int)
        if idx.size == 0:
            idx = np.arange(A.shape[0], dtype=int)
        A = A[idx]
        B = B[idx]
    diff = A - B
    return float(np.sqrt((diff * diff).sum() / A.shape[0]))


def _rmsd_between(ga, gb, indices: Optional[Sequence[int]] = None) -> float:
    """RMSD between two pysisyphus Geometries (no alignment; optional subset selection)."""
    return _calc_rmsd(np.array(ga.coords3d), np.array(gb.coords3d), indices=indices)


def _gs_cfg_with_overrides(base: Dict[str, Any], **overrides: Any) -> Dict[str, Any]:
    """
    Shallow copy of a GS config with specified overrides.
    """
    cfg = dict(base)
    for k, v in overrides.items():
        cfg[k] = v
    return cfg



def _new_geom_from_coords(atoms: Sequence[str], coords: np.ndarray, coord_type: str, freeze_atoms: Sequence[int]) -> Any:
    """
    Create a pysisyphus Geometry from Bohr coords via temporary XYZ; attach `freeze_atoms`.
    """
    lines = [str(len(atoms)), ""]
    coords_ang = np.asarray(coords, dtype=float) * BOHR2ANG
    for sym, (x, y, z) in zip(atoms, coords_ang):
        lines.append(f"{sym} {x:.15f} {y:.15f} {z:.15f}")
    s = "\n".join(lines) + "\n"
    return geom_from_xyz_string(s, coord_type=coord_type, freeze_atoms=freeze_atoms)


def _make_linear_interpolations(gL, gR, n_internal: int) -> List[Any]:
    """
    Return `n_internal` linearly interpolated structures between gL → gR (excluding endpoints).
    Atom order follows `gL`.
    """
    A = np.asarray(gL.coords3d, dtype=float)
    B = np.asarray(gR.coords3d, dtype=float)
    assert A.shape == B.shape and A.shape[1] == 3, "Atom counts must match for interpolation."
    atoms = [a for a in gL.atoms]
    coord_type = gL.coord_type
    faL = getattr(gL, "freeze_atoms", None)
    faR = getattr(gR, "freeze_atoms", None)

    freeze_union = sorted(
        set(map(int, as_list(faL))) | set(map(int, as_list(faR)))
    )
    interps: List[Any] = []
    for k in range(1, n_internal + 1):
        t = k / (n_internal + 1.0)
        C = (1.0 - t) * A + t * B
        interps.append(_new_geom_from_coords(atoms, C, coord_type, freeze_union))
    return interps


# ---- Segment/bridge tagging helpers ----

def _tag_images(images: Sequence[Any], **attrs: Any) -> None:
    """
    Attach arbitrary attributes to Geometry images.
    """
    warned = False
    for im in images:
        for k, v in attrs.items():
            try:
                setattr(im, k, v)
            except Exception:
                if not warned:
                    click.echo(
                        f"[tag] WARNING: Failed to set attribute '{k}' on an image.",
                        err=True,
                    )
                    warned = True


def _frame_ranges_by_segment(images: Sequence[Any]) -> Dict[int, Dict[str, Any]]:
    """Return additive half-open frame ranges for each tagged MEP segment."""

    result: Dict[int, Dict[str, Any]] = {}
    run_index: Optional[int] = None
    run_kind = ""
    run_start = 0

    def close_run(stop: int) -> None:
        if run_index is None or run_index <= 0:
            return
        entry = result.setdefault(
            run_index,
            {"kind": run_kind or "seg", "frame_ranges": []},
        )
        entry["frame_ranges"].append([run_start, stop])

    for frame_index, image in enumerate(images):
        segment_index = int(getattr(image, "mep_seg_index", 0) or 0)
        segment_kind = str(getattr(image, "mep_seg_kind", "") or "seg")
        if (segment_index, segment_kind) != (run_index, run_kind):
            close_run(frame_index)
            run_index = segment_index
            run_kind = segment_kind
            run_start = frame_index
    close_run(len(images))

    for entry in result.values():
        ranges = entry["frame_ranges"]
        if len(ranges) == 1:
            entry["frame_start"], entry["frame_stop"] = ranges[0]
    return result


def _segment_base_id(tag: str) -> str:
    """
    Extract base id 'seg_XXX' from a tag like 'seg_000_refine'; fallback to `tag` or 'seg'.
    """
    m = re.search(r"(seg_\d{3})", tag or "")
    return m.group(1) if m else (tag or "seg")


def _select_hei_index(energies: Sequence[float]) -> int:
    """Return the global highest-energy-image index."""
    E = np.asarray(energies, dtype=float)
    if E.size == 0:
        raise ValueError("Cannot select an HEI from an empty energy profile.")
    if not np.all(np.isfinite(E)):
        raise ValueError("Cannot select an HEI from non-finite energies.")
    return int(np.argmax(E))


def _is_local_minimum(idx: int, energies: Sequence[float]) -> bool:
    """
    Return True if the point at `idx` is a local minimum according to neighboring energies.

    Definition:
    - Interior point: both neighbors have higher energy.
    - Endpoint: the single neighbor has higher energy.
    """

    if idx < 0 or idx >= len(energies):
        return False
    if idx == 0:
        return len(energies) > 1 and energies[1] > energies[0]
    if idx == len(energies) - 1:
        return energies[-2] > energies[-1]
    return energies[idx - 1] > energies[idx] and energies[idx + 1] > energies[idx]


def _find_nearest_local_minimum(hei_idx: int, direction: int, energies: Sequence[float]) -> Optional[int]:
    """
    Starting from the HEI, walk `direction` (-1 for left, +1 for right) to find the nearest
    local minimum. Returns the index of the first minimum encountered, or None if absent.
    """

    i = hei_idx + direction
    while 0 <= i < len(energies):
        if _is_local_minimum(i, energies):
            return i
        i += direction
    return None


@dataclass
class GSMResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int
    # reported convergence of the string/DMF optimizer that produced this
    # MEP. ``None`` means no readable convergence signal (fail-closed: never
    # promoted to a usable segment by artifact existence alone).
    is_converged: Optional[bool] = None


# ---- Per‑segment summary for the console report ----
@dataclass
class SegmentReport:
    tag: str
    barrier_kcal: Optional[float]
    delta_kcal: float
    summary: str  # summarize_changes string (empty for bridges)
    kind: str = "seg"          # "seg", "kink", or "bridge"
    seg_index: int = 0         # 1‑based index along final MEP (assigned later)
    # the segment's optimizer convergence, threaded from GSMResult. A
    # reactive segment whose optimizer did not explicitly converge is unusable
    # and cannot make the path aggregate a scientific success.
    converged: Optional[bool] = None


_PRIMARY_GJF_TEMPLATE: Optional[GjfTemplate] = None


def _run_mep_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    gjf_template: Optional[GjfTemplate] = None,
    prepared_input: Optional[PreparedInputStructure] = None,
) -> GSMResult:
    """
    Run GSM between `gA`–`gB`, save segment outputs, and return images/energies/HEI index.
    """
    if gjf_template is None:
        gjf_template = _PRIMARY_GJF_TEMPLATE
    for g in (gA, gB):
        g.set_calculator(shared_calc)

    def calc_getter():
        return shared_calc

    gs = GrowingString(
        images=[gA, gB],
        calc_getter=calc_getter,
        **gs_cfg,
    )

    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)
    _opt_args = dict(stopt_cfg)
    _opt_args["out_dir"] = str(seg_dir)

    optimizer = StringOptimizer(
        geometry=gs,
        **{k: v for k, v in _opt_args.items() if k != "type"}
    )

    emit(f"\n====== [{tag}] GSM ======\n", narrative=True)
    optimizer.run()

    from pdb2reaction.workflows._outcomes import optimizer_converged_bit
    _seg_converged = optimizer_converged_bit(optimizer)

    energies = list(map(float, np.array(gs.energy, dtype=float)))
    images = list(gs.images)

    try:
        hei_idx = _select_hei_index(energies)
    except ValueError as exc:
        raise click.ClickException(f"{tag}: {exc}") from exc

    # Write trajectory
    final_trj = seg_dir / "final_geometries_trj.xyz"
    wrote_with_energy = True
    try:
        write_xyz_trj_with_energy(images, energies, final_trj)
        click.echo(f"[{tag}] Wrote '{final_trj}'.")
    except Exception:
        wrote_with_energy = False
        with open(final_trj, "w") as f:
            f.write(gs.as_xyz())
        click.echo(f"[{tag}] Wrote '{final_trj}'.")

    # Energy plot for the segment
    try:
        if wrote_with_energy:
            run_trj2fig(final_trj, [seg_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            close_matplotlib_figures()
            emit(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'", detail=True)
        else:
            click.echo(f"[{tag}] WARNING: Energies missing; skipping plot.", err=True)
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    # Convert trajectory and HEI outputs based on the input template
    prepared_for_outputs = prepared_input
    ref_for_conv = ref_pdb_path
    if prepared_for_outputs is None and gjf_template is not None:
        prepared_for_outputs = PreparedInputStructure(
            source_path=gjf_template.path,
            geom_path=gjf_template.path,
            gjf_template=gjf_template,
        )
    if prepared_for_outputs is None and ref_for_conv is not None:
        prepared_for_outputs = PreparedInputStructure(
            source_path=ref_for_conv,
            geom_path=ref_for_conv,
        )
    if prepared_for_outputs is not None and ref_for_conv is None:
        if prepared_for_outputs.source_path.suffix.lower() == ".pdb":
            ref_for_conv = prepared_for_outputs.source_path.resolve()

    needs_pdb = ref_for_conv is not None
    needs_gjf = bool(prepared_for_outputs and prepared_for_outputs.is_gjf)

    if prepared_for_outputs is not None and (needs_pdb or needs_gjf):
        try:
            convert_xyz_like_outputs(
                final_trj,
                prepared_for_outputs,
                ref_pdb_path=ref_for_conv,
                out_pdb_path=seg_dir / "final_geometries.pdb" if needs_pdb else None,
                out_gjf_path=seg_dir / "final_geometries.gjf" if needs_gjf else None,
            )
        except Exception as e:
            click.echo(
                f"[{tag}] WARNING: Failed to convert segment trajectory: {e}"
            )

    # Write HEI structure
    try:
        hei_geom = images[hei_idx]
        hei_E = float(energies[hei_idx])
        hei_xyz = seg_dir / "hei.xyz"
        s = hei_geom.as_xyz()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{hei_E:.12f}"
        s_out = "\n".join(lines)
        if not s_out.endswith("\n"):
            s_out += "\n"
        with open(hei_xyz, "w") as f:
            f.write(s_out)
        click.echo(f"[{tag}] Wrote '{hei_xyz}'.")

        if prepared_for_outputs is not None and (needs_pdb or needs_gjf):
            try:
                convert_xyz_like_outputs(
                    hei_xyz,
                    prepared_for_outputs,
                    ref_pdb_path=ref_for_conv,
                    out_pdb_path=seg_dir / "hei.pdb" if needs_pdb else None,
                    out_gjf_path=seg_dir / "hei.gjf" if needs_gjf else None,
                )
            except Exception as e:
                click.echo(
                    f"[{tag}] WARNING: Failed to convert HEI structure: {e}",
                    err=True,
                )
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to write HEI structure: {e}", err=True)

    return GSMResult(images=images, energies=energies, hei_idx=hei_idx,
                     is_converged=_seg_converged)


def _ase_atoms_to_geom(atoms, coord_type: str, template_g=None, shared_calc=None):
    """
    Convert ASE Atoms → pysisyphus Geometry, preserving freeze_atoms when present.
    """

    from ase.io import write as ase_write

    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        ase_write(tmp, atoms, format="xyz")
        tmp.flush()
        tmp.close()
        g = geom_loader(
            Path(tmp.name),
            coord_type=coord_type,
            freeze_atoms=getattr(template_g, "freeze_atoms", []),
        )
        set_freeze_atoms_or_warn(
            g,
            getattr(template_g, "freeze_atoms", []),
            context="path_search",
        )
        if shared_calc is not None:
            g.set_calculator(shared_calc)
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _run_dmf_between(
    gA,
    gB,
    calc_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    max_nodes: int,
    prepared_inputs: Sequence[PreparedInputStructure],
    shared_calc,
    dmf_cfg: Dict[str, Any],
) -> GSMResult:
    """
    Run DMF for a segment and convert outputs to pysisyphus Geometries.
    """

    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)

    fix_atoms: List[int] = []
    try:
        fix_atoms = sorted(
            {int(i) for g in [gA,gB] for i in getattr(g, "freeze_atoms", [])}
        )
    except Exception as exc:
        logger.debug("Failed to collect freeze_atoms for DMF segment %s: %s", tag, exc)

    dmf_res = _run_dmf_mep(
        geoms=[gA, gB],
        calc_cfg=calc_cfg,
        out_dir_path=seg_dir,
        prepared_inputs=prepared_inputs,
        max_nodes=max_nodes,
        fix_atoms=fix_atoms,
        dmf_cfg=dmf_cfg,
    )

    energies = list(map(float, dmf_res.energies))

    final_trj = seg_dir / "final_geometries_trj.xyz"
    write_xyz_trj_with_energy(dmf_res.images, energies, final_trj)
    _convert_to_pdb_logged(final_trj, ref_pdb_path)

    try:
        run_trj2fig(final_trj, [seg_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
        close_matplotlib_figures()
        emit(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'", detail=True)
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    imgs: List[Any] = []
    for atoms in dmf_res.images:
        imgs.append(
            _ase_atoms_to_geom(atoms, coord_type=gA.coord_type, template_g=gA, shared_calc=shared_calc)
        )

    return GSMResult(images=imgs, energies=energies, hei_idx=int(dmf_res.hei_idx),
                     is_converged=getattr(dmf_res, "is_converged", None))


def _refine_between(
    gL,
    gR,
    shared_calc,
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    mep_mode_kind: str,
    calc_cfg: Dict[str, Any],
    max_nodes: int,
    prepared_inputs: Sequence[PreparedInputStructure],
    dmf_cfg: Dict[str, Any],
) -> GSMResult:
    """
    Refine End1–End2 via GSM or DMF depending on the selected mode.
    """

    if mep_mode_kind == "dmf":
        return _run_dmf_between(
            gL,
            gR,
            calc_cfg,
            out_dir,
            tag=f"{tag}_refine",
            ref_pdb_path=ref_pdb_path,
            max_nodes=max_nodes,
            prepared_inputs=prepared_inputs,
            shared_calc=shared_calc,
            dmf_cfg=dmf_cfg,
        )

    return _run_mep_between(gL, gR, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=f"{tag}_refine", ref_pdb_path=ref_pdb_path)


def _bridge_segments(
    tail_g,
    head_g,
    shared_calc,
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    rmsd_thresh: float,
    ref_pdb_path: Optional[Path],
    mep_mode_kind: str,
    calc_cfg: Dict[str, Any],
    max_nodes: int,
    prepared_inputs: Sequence[PreparedInputStructure],
    dmf_cfg: Dict[str, Any],
) -> Optional[GSMResult]:
    """
    Run a bridge GSM if two segment endpoints are farther than the threshold.
    """
    rmsd = _rmsd_between(tail_g, head_g)
    if rmsd <= rmsd_thresh:
        return None
    click.echo(
        f"[{tag}] Gap detected between segments (RMSD={rmsd:.4e} bohr) — bridging via {mep_mode_kind.upper()}."
    )
    if mep_mode_kind == "dmf":
        return _run_dmf_between(
            tail_g,
            head_g,
            calc_cfg,
            out_dir,
            tag=f"{tag}_bridge",
            ref_pdb_path=ref_pdb_path,
            max_nodes=max_nodes,
            prepared_inputs=prepared_inputs,
            shared_calc=shared_calc,
            dmf_cfg=dmf_cfg,
        )

    return _run_mep_between(tail_g, head_g, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=f"{tag}_bridge", ref_pdb_path=ref_pdb_path)


def _stitch_paths(
    parts: List[Tuple[List[Any], List[float]]],
    stitch_rmsd_thresh: float,
    bridge_rmsd_thresh: float,
    shared_calc,
    gs_cfg,
    stopt_cfg,
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    bond_cfg: Optional[Dict[str, Any]] = None,
    segment_builder: Optional[
        Callable[[Any, Any, str, Optional[int]], "CombinedPath"]
    ] = None,
    segments_out: Optional[List["SegmentReport"]] = None,
    bridge_pair_index: Optional[int] = None,
    mep_mode_kind: str = "dmf",
    calc_cfg: Optional[Dict[str, Any]] = None,
    max_nodes: int = 20,
    prepared_inputs: Optional[Sequence[PreparedInputStructure]] = None,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[List[Any], List[float]]:
    """
    Concatenate path parts (images, energies). Insert bridge GSMs when needed.
    If covalent changes are detected across an interface, build and insert a *new* recursive segment
    using `segment_builder` instead of bridging. Update `segments_out` accordingly.
    """
    dmf_cfg = fresh_dmf_config(dmf_cfg)
    all_imgs: List[Any] = []
    all_E: List[float] = []

    def _last_known_seg_tag_from_images(imgs: List[Any]) -> Optional[str]:
        for im in reversed(imgs):
            t = getattr(im, "mep_seg_tag", None)
            if t:
                return t
        return None

    def _first_known_seg_tag_from_images(imgs: List[Any]) -> Optional[str]:
        for im in imgs:
            t = getattr(im, "mep_seg_tag", None)
            if t:
                return t
        return None

    def append_part(imgs: List[Any], Es: List[float]) -> None:
        nonlocal all_imgs, all_E
        if not imgs:
            return
        if not all_imgs:
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return
        tail = all_imgs[-1]
        head = imgs[0]

        adj_changed, adj_summary = False, ""
        if segment_builder is not None and bond_cfg is not None:
            adj_changed, adj_summary = has_bond_change(tail, head, bond_cfg)

        if adj_changed and segment_builder is not None:
            emit(f"[{tag}] Covalent changes detected at interface — inserting a new recursive segment.", narrative=True)
            if adj_summary:
                click.echo(textwrap.indent(adj_summary, prefix="  "))
            sub = segment_builder(
                tail,
                head,
                f"{tag}_mid",
                bridge_pair_index,
            )
            seg_imgs, seg_E = sub.images, sub.energies
            if segments_out is not None and getattr(sub, "segments", None):
                segments_out.extend(sub.segments)
            if seg_imgs:
                if _rmsd_between(all_imgs[-1], seg_imgs[0]) <= stitch_rmsd_thresh:
                    seg_imgs = seg_imgs[1:]
                    seg_E = seg_E[1:]
                all_imgs.extend(seg_imgs)
                all_E.extend(seg_E)
            if _rmsd_between(all_imgs[-1], imgs[0]) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return

        rmsd = _rmsd_between(tail, head)
        if rmsd <= stitch_rmsd_thresh:
            all_imgs.extend(imgs[1:])
            all_E.extend(Es[1:])
        elif rmsd > bridge_rmsd_thresh:
            left_tag_recent = _last_known_seg_tag_from_images(all_imgs) or "segL"
            right_tag_upcoming = _first_known_seg_tag_from_images(imgs) or "segR"
            left_base = _segment_base_id(left_tag_recent)
            right_base = _segment_base_id(right_tag_upcoming)
            bridge_name_base = f"{left_base}_{right_base}"

            br = _bridge_segments(
                tail, head, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=bridge_name_base,
                rmsd_thresh=bridge_rmsd_thresh, ref_pdb_path=ref_pdb_path,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg or {}, max_nodes=max_nodes,
                prepared_inputs=prepared_inputs or [],
                dmf_cfg=dmf_cfg,
            )
            if br is not None:
                _tag_images(br.images, mep_seg_tag=f"{bridge_name_base}_bridge", mep_seg_kind="bridge",
                            mep_has_bond_changes=False, pair_index=bridge_pair_index)
                b_imgs, b_E = br.images, br.energies
                if _rmsd_between(all_imgs[-1], b_imgs[0]) <= stitch_rmsd_thresh:
                    b_imgs = b_imgs[1:]
                    b_E = b_E[1:]
                if b_imgs:
                    all_imgs.extend(b_imgs)
                    all_E.extend(b_E)

                if segments_out is not None:
                    try:
                        barrier_kcal = (max(br.energies) - br.energies[0]) * AU2KCALPERMOL
                        delta_kcal = (br.energies[-1] - br.energies[0]) * AU2KCALPERMOL
                    except Exception as exc:
                        logger.debug("Failed to compute bridge energy barrier/delta: %s", exc)
                        barrier_kcal = float("nan")
                        delta_kcal = float("nan")
                    bridge_report = SegmentReport(
                        tag=f"{bridge_name_base}_bridge",
                        barrier_kcal=float(barrier_kcal),
                        delta_kcal=float(delta_kcal),
                        summary="",
                        kind="bridge",
                        converged=getattr(br, "is_converged", None),
                    )
                    insert_pos: Optional[int] = None
                    try:
                        for j, sr in enumerate(segments_out):
                            if sr.tag == right_tag_upcoming:
                                insert_pos = j
                                break
                    except Exception as exc:
                        logger.debug("Failed to find insert position for bridge segment: %s", exc)
                        insert_pos = None
                    if insert_pos is None:
                        segments_out.append(bridge_report)
                    else:
                        segments_out.insert(insert_pos, bridge_report)

            if _rmsd_between(all_imgs[-1], imgs[0]) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
        else:
            all_imgs.extend(imgs)
            all_E.extend(Es)

    for (imgs, Es) in parts:
        append_part(imgs, Es)

    return all_imgs, all_E


# Recursive search (core)

@dataclass
class CombinedPath:
    images: List[Any]
    energies: List[float]
    segments: List[SegmentReport]
    required_outcomes: List[Any] = field(default_factory=list)


def _raw_path_outcome(
    item_id: str,
    *,
    engine_converged: Optional[bool],
    artifacts: Sequence[str] = (),
):
    """Return an unusable diagnostic path outcome with engine provenance."""
    from pdb2reaction.workflows._outcomes import LeafOutcome

    converged = (
        engine_converged if isinstance(engine_converged, bool) else None
    )
    reason = "endpoint_hei"
    if converged is False:
        reason = "endpoint_hei;engine_nonconverged"
    return LeafOutcome(
        stage="path",
        item_id=item_id,
        required=True,
        executed=True,
        converged=converged,
        usable=False,
        reason=reason,
        artifacts=tuple(str(artifact) for artifact in artifacts),
    )


def _path_leaves_and_expected(
    segments: Sequence[SegmentReport],
    *,
    raw_artifacts: Sequence[str] = (),
    engine_converged: Optional[bool] = True,
    required_outcomes: Sequence[Any] = (),
):
    """Build path :class:`LeafOutcome` list + expected required IDs.

    Every segment used in the final path is required, including connector
    bridges.  When there is no reactive segment at all — the
    endpoint-HEI branch returns ``segments=[]`` even though an R/P energy diagram
    can still be drawn — an unusable ``raw_path`` leaf is emitted so the aggregate
    mapper cannot promote the diagnostic diagram to success. The raw
    trajectory/diagram remain reportable as artifacts.
    """

    from pdb2reaction.workflows._outcomes import make_leaf

    leaves: List[Any] = list(required_outcomes)
    reactive = [s for s in segments if getattr(s, "kind", "seg") == "seg"]
    for s in segments:
        # a reactive segment is usable only when its optimizer explicitly
        # converged. A nonconverged (max-cycle) StringOptimizer segment retains
        # its trajectory artifact but must not count toward completeness.
        _seg_conv = getattr(s, "converged", None)
        leaves.append(
            make_leaf(
                "path",
                f"segment_{int(s.seg_index)}",
                required=True,
                executed=True,
                converged=_seg_conv,
            )
        )
    expected = [
        outcome.item_id
        for outcome in required_outcomes
        if getattr(outcome, "required", True)
    ]
    expected.extend(f"segment_{int(s.seg_index)}" for s in segments)
    if not reactive:
        if not required_outcomes:
            leaves.append(
                _raw_path_outcome(
                    "raw_path",
                    engine_converged=engine_converged,
                    artifacts=raw_artifacts,
                )
            )
            expected.append("raw_path")
        expected.append("reactive_segment_1")
    return leaves, expected


def _trailing_kink_count(segments: Sequence[SegmentReport]) -> int:
    """Return the number of consecutive kink segments at the end of ``segments``."""

    count = 0
    for seg in reversed(segments):
        if seg.tag and "kink" in seg.tag:
            count += 1
        else:
            break
    return count


def _build_multistep_path(
    gA,
    gB,
    shared_calc,
    geom_cfg: Dict[str, Any],
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
    single_opt_kind: str,
    single_opt_cfg: Dict[str, Any],
    bond_cfg: Dict[str, Any],
    search_cfg: Dict[str, Any],
    refine_mode_kind: str,
    mep_mode_kind: str,
    calc_cfg: Dict[str, Any],
    dmf_cfg: Dict[str, Any],
    prepared_inputs: Sequence[PreparedInputStructure],
    out_dir: Path,
    ref_pdb_path: Optional[Path],
    prepared_input: Optional[PreparedInputStructure],
    depth: int,
    seg_counter: List[int],
    branch_tag: str,
    pair_index: Optional[int] = None,
    kink_seq_count: int = 0,
) -> CombinedPath:
    """
    Recursively construct a multistep MEP from A–B and return it (A→B order).
    """
    seg_max_nodes = int(
        search_cfg.get(
            "max_nodes_segment",
            gs_cfg.get("max_nodes", GS_KW["max_nodes"]),
        )
    )
    gs_seg_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=seg_max_nodes)
    max_seq_kink = int(search_cfg.get("max_seq_kink", 2))

    def _terminate_with_maxdepth(reason_msg: Optional[str] = None) -> CombinedPath:
        if reason_msg:
            click.echo(reason_msg)

        seg_tag = f"seg_{seg_counter[0]:03d}_maxdepth"
        gsm = (
            _run_dmf_between(
                gA,
                gB,
                calc_cfg,
                out_dir,
                tag=seg_tag,
                ref_pdb_path=ref_pdb_path,
                max_nodes=seg_max_nodes,
                prepared_inputs=prepared_inputs,
                shared_calc=shared_calc,
                dmf_cfg=dmf_cfg,
            )
            if mep_mode_kind == "dmf"
            else _run_mep_between(
                gA,
                gB,
                shared_calc,
                gs_seg_cfg,
                stopt_cfg,
                out_dir,
                tag=seg_tag,
                ref_pdb_path=ref_pdb_path,
                prepared_input=prepared_input,
            )
        )
        seg_counter[0] += 1

        if not (1 <= int(gsm.hei_idx) <= len(gsm.images) - 2):
            click.echo(
                f"[{seg_tag}] WARNING: HEI is at an endpoint. Returning the raw path.",
                err=True,
            )
            _tag_images(gsm.images, pair_index=pair_index)
            return CombinedPath(
                images=gsm.images,
                energies=gsm.energies,
                segments=[],
                required_outcomes=[
                    _raw_path_outcome(
                        f"raw_{seg_tag}",
                        engine_converged=getattr(gsm, "is_converged", None),
                    )
                ],
            )

        try:
            changed, step_summary = has_bond_change(gsm.images[0], gsm.images[-1], bond_cfg)
        except Exception as e:
            click.echo(f"[{seg_tag}] WARNING: Failed to evaluate bond changes at max depth: {e}", err=True)
            changed, step_summary = True, ""

        try:
            barrier_kcal = (max(gsm.energies) - gsm.energies[0]) * AU2KCALPERMOL
            delta_kcal = (gsm.energies[-1] - gsm.energies[0]) * AU2KCALPERMOL
        except Exception as exc:
            logger.debug("Failed to compute segment energy barrier/delta: %s", exc)
            barrier_kcal = float("nan")
            delta_kcal = float("nan")

        seg_report = SegmentReport(
            tag=seg_tag,
            barrier_kcal=float(barrier_kcal),
            delta_kcal=float(delta_kcal),
            summary=step_summary if changed else "(no covalent changes detected)",
            kind="seg",
            converged=getattr(gsm, "is_converged", None),
        )

        _tag_images(
            gsm.images,
            mep_seg_tag=seg_tag,
            mep_seg_kind="seg",
            mep_has_bond_changes=bool(changed),
            pair_index=pair_index,
        )

        return CombinedPath(images=gsm.images, energies=gsm.energies, segments=[seg_report])

    if depth > int(search_cfg.get("max_depth", 10)):
        click.echo(f"[{branch_tag}] Reached maximum recursion depth. Returning current endpoints only.")
        return _terminate_with_maxdepth()

    seg_id = seg_counter[0]
    seg_counter[0] += 1
    tag0 = f"seg_{seg_id:03d}"

    gsm0 = (
        _run_dmf_between(
            gA,
            gB,
            calc_cfg,
            out_dir,
            tag=tag0,
            ref_pdb_path=ref_pdb_path,
            max_nodes=seg_max_nodes,
            prepared_inputs=prepared_inputs,
            shared_calc=shared_calc,
            dmf_cfg=dmf_cfg,
        )
        if mep_mode_kind == "dmf"
        else _run_mep_between(
            gA,
            gB,
            shared_calc,
            gs_seg_cfg,
            stopt_cfg,
            out_dir,
            tag=tag0,
            ref_pdb_path=ref_pdb_path,
            prepared_input=prepared_input,
        )
    )

    hei = int(gsm0.hei_idx)
    if not (1 <= hei <= len(gsm0.images) - 2):
        click.echo(f"[{tag0}] WARNING: HEI is at an endpoint (idx={hei}). Returning the raw GSM path.", err=True)
        _tag_images(gsm0.images, pair_index=pair_index)
        return CombinedPath(
            images=gsm0.images,
            energies=gsm0.energies,
            segments=[],
            required_outcomes=[
                _raw_path_outcome(
                    f"raw_{tag0}",
                    engine_converged=getattr(gsm0, "is_converged", None),
                )
            ],
        )

    if refine_mode_kind == "minima":
        left_idx = _find_nearest_local_minimum(hei_idx=hei, direction=-1, energies=gsm0.energies)
        right_idx = _find_nearest_local_minimum(hei_idx=hei, direction=1, energies=gsm0.energies)

        if left_idx is None:
            left_idx = hei - 1
        if right_idx is None:
            right_idx = hei + 1

        click.echo(
            f"[{tag0}] Using nearest local minima around HEI (left idx={left_idx}, right idx={right_idx})."
        )
        left_img = gsm0.images[left_idx]
        right_img = gsm0.images[right_idx]
    else:
        left_img = gsm0.images[hei - 1]
        right_img = gsm0.images[hei + 1]
        emit(f"[{tag0}] Refining HEI±1 (peak mode).", narrative=True)

    left_end, left_conv = _optimize_single(
        left_img,
        shared_calc,
        single_opt_kind,
        single_opt_cfg,
        out_dir,
        tag=f"{tag0}_left",
        prepared_input=prepared_input,
        ref_pdb=ref_pdb_path,
    )
    right_end, right_conv = _optimize_single(
        right_img,
        shared_calc,
        single_opt_kind,
        single_opt_cfg,
        out_dir,
        tag=f"{tag0}_right",
        prepared_input=prepared_input,
        ref_pdb=ref_pdb_path,
    )

    try:
        lr_changed, _ = has_bond_change(left_end, right_end, bond_cfg)
    except Exception as e:
        click.echo(f"[{tag0}] WARNING: Failed to evaluate bond changes for kink detection: {e}", err=True)
        lr_changed, _ = True, ""
    use_kink = (not lr_changed)
    from pdb2reaction.workflows._outcomes import combine_step_convergence

    if use_kink:
        n_inter = int(search_cfg.get("kink_max_nodes", 3))
        emit(f"[{tag0}] Kink detected (no covalent changes between End1 and End2). "
                   f"Using {n_inter} linear interpolation nodes + single-structure optimizations instead of GSM.",
                   narrative=True)
        inter_geoms = _make_linear_interpolations(left_end, right_end, n_inter)
        opt_inters: List[Any] = []
        inter_convs: List[Optional[bool]] = []
        for i, g_int in enumerate(inter_geoms, 1):
            g_int.set_calculator(shared_calc)
            g_opt, _inter_conv = _optimize_single(
                g_int,
                shared_calc,
                single_opt_kind,
                single_opt_cfg,
                out_dir,
                tag=f"{tag0}_kink_int{i}",
                prepared_input=prepared_input,
                ref_pdb=ref_pdb_path,
            )
            opt_inters.append(g_opt)
            inter_convs.append(_inter_conv)
        step_imgs = [left_end] + opt_inters + [right_end]
        step_E = [float(img.energy) for img in step_imgs]
        _kink_hei = int(np.argmax(step_E[1:-1])) + 1 if len(step_E) > 2 else int(np.argmax(step_E))
        # A kink segment is assembled from single-structure optimizations (not a
        # StringOptimizer). It is usable only when EVERY endpoint/intermediate
        # optimization explicitly converged; fold their convergence rather than
        # hardcode True, so a nonconverged (max-cycle) single-structure opt cannot
        # silently become a usable reactive leaf (fail-closed).
        _kink_converged = combine_step_convergence(
            [left_conv] + inter_convs + [right_conv]
        )
        ref1 = GSMResult(images=step_imgs, energies=step_E, hei_idx=_kink_hei,
                         is_converged=_kink_converged)
        step_tag_for_report = f"{tag0}_kink"
    else:
        ref1 = _refine_between(
            left_end,
            right_end,
            shared_calc,
            gs_seg_cfg,
            stopt_cfg,
            out_dir,
            tag=tag0,
            ref_pdb_path=ref_pdb_path,
            mep_mode_kind=mep_mode_kind,
            calc_cfg=calc_cfg,
            max_nodes=seg_max_nodes,
            prepared_inputs=prepared_inputs,
            dmf_cfg=dmf_cfg,
        )
        step_tag_for_report = f"{tag0}_refine"

    step_imgs, step_E = ref1.images, ref1.energies

    if not (1 <= int(ref1.hei_idx) <= len(step_imgs) - 2):
        click.echo(
            f"[{step_tag_for_report}] WARNING: HEI is at an endpoint. "
            "Returning the raw refined path.",
            err=True,
        )
        _tag_images(step_imgs, pair_index=pair_index)
        return CombinedPath(
            images=step_imgs,
            energies=step_E,
            segments=[],
            required_outcomes=[
                _raw_path_outcome(
                    f"raw_{step_tag_for_report}",
                    engine_converged=getattr(ref1, "is_converged", None),
                )
            ],
        )

    _changed, step_summary = has_bond_change(step_imgs[0], step_imgs[-1], bond_cfg)
    step_kind = "kink" if use_kink else "seg"
    _tag_images(step_imgs, mep_seg_tag=step_tag_for_report, mep_seg_kind=step_kind,
                mep_has_bond_changes=bool(_changed), pair_index=pair_index)

    left_changed, left_summary = has_bond_change(gA, left_end, bond_cfg)
    right_changed, right_summary = has_bond_change(right_end, gB, bond_cfg)

    emit(f"[{tag0}] Covalent changes (A vs left_end): {'Yes' if left_changed else 'No'}", narrative=True)
    if left_changed:
        click.echo(textwrap.indent(left_summary, prefix="  "))
    emit(f"[{tag0}] Covalent changes (right_end vs B): {'Yes' if right_changed else 'No'}", narrative=True)
    if right_changed:
        click.echo(textwrap.indent(right_summary, prefix="  "))

    try:
        barrier_kcal = (max(step_E) - step_E[0]) * AU2KCALPERMOL
        delta_kcal = (step_E[-1] - step_E[0]) * AU2KCALPERMOL
    except Exception as exc:
        logger.debug("Failed to compute step energy barrier/delta: %s", exc)
        barrier_kcal = float("nan")
        delta_kcal = float("nan")

    segment_converged = combine_step_convergence(
        [
            getattr(gsm0, "is_converged", None),
            left_conv,
            right_conv,
            getattr(ref1, "is_converged", None),
        ]
    )
    seg_report = SegmentReport(
        tag=step_tag_for_report,
        barrier_kcal=None if use_kink else float(barrier_kcal),
        delta_kcal=float(delta_kcal),
        summary=step_summary if _changed else "(no covalent changes detected)",
        kind=step_kind,
        converged=segment_converged,
    )

    parts: List[Tuple[List[Any], List[float]]] = []
    seg_reports: List[SegmentReport] = []
    required_outcomes: List[Any] = []

    trailing_kink_run = kink_seq_count
    if left_changed:
        subL = _build_multistep_path(
            gA, left_end, shared_calc, geom_cfg, gs_cfg, stopt_cfg,
            single_opt_kind, single_opt_cfg, bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
            out_dir, ref_pdb_path, prepared_input, depth + 1, seg_counter, branch_tag=f"{branch_tag}L",
            pair_index=pair_index,
            kink_seq_count=kink_seq_count,
        )
        _tag_images(subL.images, pair_index=pair_index)
        parts.append((subL.images, subL.energies))
        seg_reports.extend(subL.segments)
        required_outcomes.extend(subL.required_outcomes)
        trailing_kink_run = _trailing_kink_count(seg_reports)

    current_kink_run = trailing_kink_run + 1 if use_kink else 0
    if use_kink and current_kink_run >= max_seq_kink:
        warning_msg = (
            f"[{tag0}] Consecutive kink segments were detected. Something seems wrong. "
            "Please check the initial structure and the generated intermediate structures. "
            "Alternatively, try switching the mep-mode. If that still fails, try including intermediate structures in the inputs."
        )
        return _terminate_with_maxdepth(reason_msg=warning_msg)

    parts.append((step_imgs, step_E))
    seg_reports.append(seg_report)

    if right_changed:
        subR = _build_multistep_path(
            right_end, gB, shared_calc, geom_cfg, gs_cfg, stopt_cfg,
            single_opt_kind, single_opt_cfg, bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
            out_dir, ref_pdb_path, prepared_input, depth + 1, seg_counter, branch_tag=f"{branch_tag}R",
            pair_index=pair_index,
            kink_seq_count=current_kink_run,
        )
        _tag_images(subR.images, pair_index=pair_index)
        parts.append((subR.images, subR.energies))
        seg_reports.extend(subR.segments)
        required_outcomes.extend(subR.required_outcomes)

    bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 5))
    gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

    def _segment_builder(
        tail_g,
        head_g,
        _tag: str,
        interface_pair_index: Optional[int],
    ) -> CombinedPath:
        sub = _build_multistep_path(
            tail_g, head_g,
            shared_calc,
            geom_cfg, gs_cfg, stopt_cfg,
            single_opt_kind, single_opt_cfg,
            bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
            out_dir=out_dir,
            ref_pdb_path=ref_pdb_path,
            prepared_input=prepared_input,
            depth=depth + 1,
            seg_counter=seg_counter,
            branch_tag=f"{branch_tag}B",
            pair_index=interface_pair_index,
            kink_seq_count=_trailing_kink_count(seg_reports),
        )
        _tag_images(sub.images, pair_index=interface_pair_index)
        required_outcomes.extend(sub.required_outcomes)
        return sub

    stitched_imgs, stitched_E = _stitch_paths(
        parts,
        stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-4)),
        bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-4)),
        shared_calc=shared_calc,
        gs_cfg=gs_bridge_cfg,
        stopt_cfg=stopt_cfg,
        out_dir=out_dir,
        tag=tag0,
        ref_pdb_path=ref_pdb_path,
        bond_cfg=bond_cfg,
        segment_builder=_segment_builder,
        segments_out=seg_reports,
        bridge_pair_index=pair_index,
        mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg,
        max_nodes=bridge_max_nodes,
        dmf_cfg=dmf_cfg,
        prepared_inputs=prepared_inputs,
    )

    _tag_images(stitched_imgs, pair_index=pair_index)

    return CombinedPath(
        images=stitched_imgs,
        energies=stitched_E,
        segments=seg_reports,
        required_outcomes=required_outcomes,
    )


# Full‑system merge helpers (Biopython)

def _atom_key_from_res_atom(res: PDB.Residue.Residue, atom: PDB.Atom.Atom) -> Tuple[str, str, str, str, str]:
    """
    Build a key for atom identity:
    (RESNAME, RESSEQ, ICODE, CHAIN, ATOMNAME) — uppercase where applicable.
    - RESSEQ is numeric (without insertion code).
    - ICODE is '' when blank (or ' ' in PDB).
    """
    chain_id, resseq_txt, icode_txt, resname = residue_auth_identity(res)
    retained_atom = atom_site_from_biopython_atom(atom)
    atname = (retained_atom.atom_name if retained_atom is not None else atom.get_name()).strip().upper()
    resname = resname.strip().upper()
    chain_id = chain_id.strip().upper()
    icode_txt = icode_txt.strip().upper()
    return (resname, resseq_txt, icode_txt, chain_id, atname)


def _structure_to_arrays(struct: PDB.Structure.Structure) -> Tuple[np.ndarray, List[PDB.Atom.Atom], List[Tuple[str, str, str, str, str]], Dict[Tuple[str,str,str,str,str], int]]:
    """
    Extract: coordinates (Å), atom list, key list, and key→index map from a Biopython Structure.
    Keys are (RESNAME, RESSEQ, ICODE, CHAIN, ATOMNAME).
    """
    atoms: List[PDB.Atom.Atom] = [a for a in struct.get_atoms()]
    coords = np.array([a.get_coord() for a in atoms], dtype=float)
    keys: List[Tuple[str, str, str, str, str]] = []
    key2idx: Dict[Tuple[str,str,str,str,str], int] = {}
    for idx, a in enumerate(atoms):
        res = a.get_parent()
        k = _atom_key_from_res_atom(res, a)
        keys.append(k)
        if k not in key2idx:
            key2idx[k] = idx
    return coords, atoms, keys, key2idx


def _load_structures_and_chain_align(ref_paths: Sequence[Path]) -> Tuple[List[PDB.Structure.Structure], List[np.ndarray], List[List[PDB.Atom.Atom]], List[Dict[Tuple[str,str,str,str,str], int]]]:
    """
    Load all full templates and rigidly chain‑align them into the coordinate frame of the first template.
    """
    parser = PDBParser(QUIET=True)
    structs: List[PDB.Structure.Structure] = []
    for i, path in enumerate(ref_paths):
        structure = parser.get_structure(f"ref{i:02d}", str(path))
        template = coordinate_template_for(path)
        if template is not None:
            attach_template_metadata(structure, template)
        structs.append(structure)
    N_expected = None
    coords_list: List[np.ndarray] = []
    atoms_list: List[List[PDB.Atom.Atom]] = []
    keymaps: List[Dict[Tuple[str,str,str,str,str], int]] = []

    for s in structs:
        coords, atoms, _keys, key2idx = _structure_to_arrays(s)
        if N_expected is None:
            N_expected = coords.shape[0]
        else:
            if coords.shape[0] != N_expected:
                raise click.BadParameter(f"[merge] Atom count mismatch among --ref-full-pdb templates: {N_expected} vs {coords.shape[0]}")
        coords_list.append(coords)
        atoms_list.append(atoms)
        keymaps.append(key2idx)

    aligned_coords: List[np.ndarray] = []
    aligned_coords.append(coords_list[0].copy())
    for j in range(1, len(coords_list)):
        P = aligned_coords[j-1]
        Q = coords_list[j]
        R, t = kabsch_R_t(P, Q)
        Qa = (Q @ R) + t
        aligned_coords.append(Qa)

    return structs, aligned_coords, atoms_list, keymaps


def _model_keys_from_pdb(model_pdb: Path) -> List[Tuple[str, str, str, str, str]]:
    """
    Return atom identity keys for a model PDB file.
    """
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("model", str(model_pdb))
    template = coordinate_template_for(model_pdb)
    if template is not None:
        attach_template_metadata(st, template)
    keys: List[Tuple[str, str, str, str, str]] = []
    for a in st.get_atoms():
        res = a.get_parent()
        k = _atom_key_from_res_atom(res, a)
        keys.append(k)
    return keys


def _coordinate_template_for_refs(
    ref_pdbs: Sequence[Path],
    preferred_index: int = 0,
):
    """Choose retained public IDs from a preferred or any bridged reference."""

    if 0 <= preferred_index < len(ref_pdbs):
        preferred = coordinate_template_for(ref_pdbs[preferred_index])
        if preferred is not None:
            return preferred
    for ref_path in ref_pdbs:
        template = coordinate_template_for(ref_path)
        if template is not None:
            return template
    return None


def _write_model_block(structure: PDB.Structure.Structure,
                       remark_lines: List[str]) -> str:
    """
    Render a single MODEL block (without 'MODEL/ENDMDL') with provided REMARK lines.
    """
    io = PDBIO()
    io.set_structure(structure)
    from io import StringIO
    buf = StringIO()
    # Internal bridges wrap atom serials at 99,999. Preserve those safe serials
    # instead of asking Biopython to generate an unrepresentable six-digit PDB
    # serial for very large structures; the public companion is mmCIF.
    io.save(buf, preserve_atom_numbering=True)
    body = "\n".join([ln for ln in buf.getvalue().splitlines() if ln.strip() != "END"])
    rem = ""
    for line in remark_lines:
        rem += f"REMARK   1 {line}\n"
    return rem + body + ("\n" if not body.endswith("\n") else "")


def _chunk_remark_indices(indices: List[int], width: int = 60) -> List[str]:
    """
    Wrap model atom indices into REMARK lines with limited width.
    """
    s = ",".join(map(str, indices))
    out: List[str] = []
    cur = ""
    for tok in s.split(","):
        add = (tok if not cur else "," + tok)
        if len(cur) + len(add) > width:
            out.append(f"POCKET_ATOM_INDICES {cur}")
            cur = tok
        else:
            cur += tok if not cur else "," + tok
    if cur:
        out.append(f"POCKET_ATOM_INDICES {cur}")
    return out


def _merge_pair_to_full(pair_images: List[Any],
                        model_ref_pdb: Path,
                        structA: PDB.Structure.Structure,
                        structB: PDB.Structure.Structure,
                        coordsA_aligned: np.ndarray,
                        coordsB_aligned: np.ndarray,
                        keymapA: Dict[Tuple[str,str,str,str,str], int],
                        keymapB: Dict[Tuple[str,str,str,str,str], int],
                        out_path: Optional[Path],
                        drop_first: bool = False,
                        seg_indices_for_frames: Optional[List[int]] = None,
                        seg_report_lookup: Optional[Dict[int, SegmentReport]] = None,
                        include_model_indices_for_first_image: bool = False) -> Tuple[List[str], List[int]]:
    """
    Merge a model‑only trajectory for a *pair* into the corresponding full templates (A,B),
    generating MODEL blocks and (optionally) writing a PDB. Returns (blocks, 1‑based active indices).
    """
    model_keys = _model_keys_from_pdb(model_ref_pdb)

    match_tpl_idx: List[int] = []
    for k in model_keys:
        ia = keymapA.get(k, None)
        ib = keymapB.get(k, None)
        if ia is None or ib is None:
            match_tpl_idx.append(-1)
        else:
            match_tpl_idx.append(int(ia))

    active_full_idx = sorted({i for i in match_tpl_idx if i >= 0})
    active_one_based = [i+1 for i in active_full_idx]

    Nfull = coordsA_aligned.shape[0]
    if Nfull != coordsB_aligned.shape[0]:
        raise click.BadParameter("[merge] Template A/B atom count mismatch.")

    atomsA: List[PDB.Atom.Atom] = [a for a in structA.get_atoms()]

    start_k = 1 if drop_first and len(pair_images) > 0 else 0

    model_blocks: List[str] = []

    M = len(pair_images)
    if M == 0:
        return model_blocks, active_one_based

    seg_idx_seq: List[int] = []
    if seg_indices_for_frames is not None and len(seg_indices_for_frames) == M:
        seg_idx_seq = seg_indices_for_frames[start_k:]
    else:
        seg_idx_seq = [0] * (M - start_k)

    for kk, k in enumerate(range(start_k, M)):
        tfrac = 0.0 if M == 1 else (k / (M - 1.0))

        C = (1.0 - tfrac) * coordsA_aligned + tfrac * coordsB_aligned

        idx_sel = [j for j in range(len(match_tpl_idx)) if match_tpl_idx[j] >= 0]
        if len(idx_sel) >= 3:
            Y = np.array([C[match_tpl_idx[j]] for j in idx_sel], dtype=float)
            P_bohr = np.array(pair_images[k].coords3d, dtype=float)
            P = P_bohr * BOHR2ANG
            P_sel = np.array([P[j] for j in idx_sel], dtype=float)
            R, t = kabsch_R_t(Y, P_sel)
            Paligned = (P @ R) + t
            for jj, pidx in enumerate(idx_sel):
                full_i = match_tpl_idx[pidx]
                if 0 <= full_i < Nfull:
                    C[full_i] = Paligned[jj]
        else:
            P = np.array(pair_images[k].coords3d, dtype=float) * BOHR2ANG
            for j, full_i in enumerate(match_tpl_idx):
                if full_i >= 0 and 0 <= full_i < Nfull and j < P.shape[0]:
                    C[full_i] = P[j]

        for i, a in enumerate(atomsA):
            a.set_coord(C[i])
            a.set_bfactor(100.0 if i in active_full_idx else 0.0)

        remark_lines: List[str] = []
        remark_lines.append(f"PAIR_MERGE FRAC {tfrac:.6f}")

        if include_model_indices_for_first_image and kk == 0:
            remark_lines.extend(_chunk_remark_indices([i for i in active_one_based], width=60))

        seg_idx = seg_idx_seq[kk] if kk < len(seg_idx_seq) else 0
        if seg_idx and seg_report_lookup is not None:
            rep = seg_report_lookup.get(seg_idx, None)
            if rep is not None:
                remark_lines.append(
                    f"MEP_SEG_INDEX {int(seg_idx):02d} TAG {rep.tag} KIND {rep.kind} "
                    f"DELTAE_BARRIER_KCAL {rep.barrier_kcal:.6f} DELTAE_KCAL {rep.delta_kcal:.6f}"
                )
                if rep.kind != "bridge" and rep.summary and rep.summary.strip() and rep.summary.strip() != "(no covalent changes detected)":
                    for ln in rep.summary.strip().splitlines():
                        remark_lines.append(f"SEG_BONDS {ln.strip()}")
            else:
                remark_lines.append(f"MEP_SEG_INDEX {int(seg_idx):02d}")

        model_blocks.append(_write_model_block(structA, remark_lines))

    if out_path is not None:
        with open(out_path, "w") as f:
            for m, blk in enumerate(model_blocks, start=1):
                f.write(f"MODEL     {m}\n")
                f.write(blk)
                f.write("ENDMDL\n")
            f.write("END\n")
        click.echo(f"[merge] Wrote pair-merged PDB → '{out_path}'")

    return model_blocks, active_one_based


def _merge_final_and_write(final_images: List[Any],
                           model_inputs: Sequence[Path],
                           ref_pdbs: Sequence[Path],
                           segments: List[SegmentReport],
                           out_dir: Path,
                           model_ref_pdbs: Optional[Sequence[Path]] = None) -> None:
    """
    Merge the entire model MEP into full templates (for all pairs) and write outputs.
    """
    if len(ref_pdbs) != len(model_inputs):
        raise click.BadParameter(
            "--ref-full-pdb must match the number of --input after preprocessing "
            "(caller should replicate the first ref for all pairs when --align).; "
            f"recover: provide --ref-full-pdb {len(model_inputs)} times "
            f"(currently {len(ref_pdbs)}), or omit it to default to --input itself."
        )

    if model_ref_pdbs is None:
        model_ref_pdbs = model_inputs
    if len(model_ref_pdbs) != len(model_inputs):
        raise click.BadParameter(
            "--ref-pdb must match the number of --input after preprocessing.; "
            f"recover: provide --ref-pdb {len(model_inputs)} times "
            f"(currently {len(model_ref_pdbs)}), in the same order as --input."
        )

    structs, aligned_coords, _atoms_list, keymaps = _load_structures_and_chain_align(ref_pdbs)

    seg_lookup: Dict[int, SegmentReport] = {int(s.seg_index): s for s in segments if int(s.seg_index) > 0}

    n_pairs = len(model_inputs) - 1
    groups: List[Tuple[int, List[Any]]] = []
    cur_idx = None
    cur_list: List[Any] = []
    for im in final_images:
        pi = getattr(im, "pair_index", None)
        if pi is None:
            pi = 0
        if cur_idx is None:
            cur_idx = int(pi)
            cur_list = [im]
        elif int(pi) == int(cur_idx):
            cur_list.append(im)
        else:
            groups.append((int(cur_idx), cur_list))
            cur_idx = int(pi)
            cur_list = [im]
    if cur_list:
        groups.append((int(cur_idx), cur_list))

    for (pi, _) in groups:
        if not (0 <= pi < n_pairs):
            raise click.BadParameter(f"[merge] Illegal pair_index {pi} (n_pairs={n_pairs}).")

    final_blocks: List[str] = []
    wrote_indices = False

    for gi, (pi, imgs) in enumerate(groups):
        model_ref = Path(model_ref_pdbs[pi])
        structA = structs[pi]
        structB = structs[pi+1]
        coordsA = aligned_coords[pi]
        coordsB = aligned_coords[pi+1]
        keymapA = keymaps[pi]
        keymapB = keymaps[pi+1]
        seg_indices_for_frames = [int(getattr(im, "mep_seg_index", 0) or 0) for im in imgs]

        blocks, active_one_based = _merge_pair_to_full(
            pair_images=imgs,
            model_ref_pdb=model_ref,
            structA=structA,
            structB=structB,
            coordsA_aligned=coordsA,
            coordsB_aligned=coordsB,
            keymapA=keymapA,
            keymapB=keymapB,
            out_path=None,
            drop_first=(gi > 0),
            seg_indices_for_frames=seg_indices_for_frames,
            seg_report_lookup=seg_lookup,
            include_model_indices_for_first_image=(not wrote_indices)
        )
        if blocks:
            wrote_indices = True
            final_blocks.extend(blocks)

    final_path = out_dir / "mep_w_ref.pdb"
    with open(final_path, "w") as f:
        for m, blk in enumerate(final_blocks, start=1):
            f.write(f"MODEL     {m}\n")
            f.write(blk)
            f.write("ENDMDL\n")
        f.write("END\n")
    click.echo(f"[merge] Wrote concatenated full-system trajectory → '{final_path}'")
    final_cif = register_output_template_and_write_cif(
        final_path,
        _coordinate_template_for_refs(ref_pdbs),
    )
    if final_cif is not None:
        click.echo(f"[merge] Wrote concatenated full-system mmCIF → '{final_cif}'")

    # Per‑segment merged MEPs (bond‑change segments only) + HEI merged only for bond‑change segments
    for s in segments:
        seg_idx = int(s.seg_index)
        seg_frames: List[Any] = [im for im in final_images if int(getattr(im, "mep_seg_index", 0) or 0) == seg_idx]
        if not seg_frames:
            continue

        # Determine pair index for this segment (assume consistent within the segment)
        pi_vals = sorted({int(getattr(im, "pair_index", 0)) for im in seg_frames})
        pi = pi_vals[0]
        model_ref = Path(model_ref_pdbs[pi])
        structA = structs[pi]
        structB = structs[pi+1]
        coordsA = aligned_coords[pi]
        coordsB = aligned_coords[pi+1]
        keymapA = keymaps[pi]
        keymapB = keymaps[pi+1]

        # Per‑segment merged MEP only when covalent changes are present
        if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
            seg_indices_for_frames = [seg_idx] * len(seg_frames)
            blocks, _ = _merge_pair_to_full(
                pair_images=seg_frames,
                model_ref_pdb=model_ref,
                structA=structA,
                structB=structB,
                coordsA_aligned=coordsA,
                coordsB_aligned=coordsB,
                keymapA=keymapA,
                keymapB=keymapB,
                out_path=None,
                drop_first=False,
                seg_indices_for_frames=seg_indices_for_frames,
                seg_report_lookup=seg_lookup,
                include_model_indices_for_first_image=True,
            )
            out_seg = out_dir / f"mep_w_ref_seg_{seg_idx:02d}.pdb"
            with open(out_seg, "w") as f:
                for m, blk in enumerate(blocks, start=1):
                    f.write(f"MODEL     {m}\n")
                    f.write(blk)
                    f.write("ENDMDL\n")
                f.write("END\n")
            click.echo(f"[merge] Wrote per-segment merged trajectory → '{out_seg}'")
            seg_cif = register_output_template_and_write_cif(
                out_seg,
                _coordinate_template_for_refs(ref_pdbs, pi),
            )
            if seg_cif is not None:
                click.echo(f"[merge] Wrote per-segment merged mmCIF → '{seg_cif}'")

        # Per‑segment HEI merged to reference (only for bond‑change segments)
        if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
            try:
                energies_au = [float(getattr(im, "energy")) for im in seg_frames]
                imax = int(np.argmax(np.array(energies_au, dtype=float)))
                hei_frame = seg_frames[imax]
                blocks_hei, _ = _merge_pair_to_full(
                    pair_images=[hei_frame],
                    model_ref_pdb=model_ref,
                    structA=structA,
                    structB=structB,
                    coordsA_aligned=coordsA,
                    coordsB_aligned=coordsB,
                    keymapA=keymapA,
                    keymapB=keymapB,
                    out_path=None,
                    drop_first=False,
                    seg_indices_for_frames=[seg_idx],
                    seg_report_lookup=seg_lookup,
                    include_model_indices_for_first_image=True,
                )
                out_hei = out_dir / f"hei_w_ref_seg_{seg_idx:02d}.pdb"
                with open(out_hei, "w") as f:
                    for m, blk in enumerate(blocks_hei, start=1):
                        f.write(f"MODEL     {m}\n")
                        f.write(blk)
                        f.write("ENDMDL\n")
                    f.write("END\n")
                click.echo(f"[merge] Wrote merged HEI for segment → '{out_hei}'")
                hei_cif = register_output_template_and_write_cif(
                    out_hei,
                    _coordinate_template_for_refs(ref_pdbs, pi),
                )
                if hei_cif is not None:
                    click.echo(f"[merge] Wrote merged HEI mmCIF → '{hei_cif}'")
            except Exception as e:
                click.echo(f"[merge] WARNING: Failed to write merged HEI for segment {seg_idx:02d}: {e}", err=True)



@click.command(
    help="Multistep MEP search via recursive GSM/DMF segmentation.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    required=True,
    help="Two or more PDB, mmCIF, XYZ, or GJF structures in reaction order. Repeat -i/--input for each path.",
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
)
@click.option(
    "--dmf-backend",
    type=click.Choice(["cpu", "gpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="DMF compute backend (--mep-mode dmf only): gpu (dmf.torch / CUDA) or cpu (dmf / NumPy). "
    "On a GPU out-of-memory error, retry with cpu.",
)
@click.option(
    "--refine-mode",
    type=click.Choice(["peak", "minima"], case_sensitive=False),
    default=None,
    show_default="peak for gsm, minima for dmf",
    help=(
        "Refinement seed selection around the highest-energy image: "
        "'peak' uses HEI±1, 'minima' uses the nearest local minima in each direction. "
        "Defaults to peak for gsm and minima for dmf when omitted."
    ),
)
@click.option(
    "-q",
    "--charge",
    type=int,
    default=None,
    show_default=False,
    help=(
        "Total charge. Required for non-.gjf inputs unless --ligand-charge derives it "
        "from PDB/mmCIF inputs."
    ),
)
@click.option(
    "--workers",
    type=int,
    default=UMA_CALC_KW["workers"],
    show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=UMA_CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel MLIP predictor (workers>1).",
)
@click.option(
    "-l",
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
        "when -q is omitted (PDB/mmCIF inputs only)."
    ),
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default="1 (or the .gjf template value)",
    help="Spin multiplicity (2S+1; inherits from a .gjf template when available).",
)
@click.option(
    "--freeze-links/--no-freeze-links",
    "freeze_links_flag",
    default=True,
    show_default=True,
    help="Freeze parent atoms of cap hydrogens (PDB/mmCIF input only).",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
)
@click.option("--max-nodes", type=int, default=GS_KW["max_nodes"], show_default=True,
              help=("Number of movable internal images per GSM/DMF segment; the complete segment "
                    "has max_nodes+2 images including endpoints. When not given, YAML "
                    "search.max_nodes_segment applies."))
@click.option(
    "--gsm-param",
    type=click.Choice(["equi", "energy"], case_sensitive=False),
    default=None,
    show_default="equi",
    help=(
        "GSM node parameterization after string growth. The energy scheme "
        "concentrates nodes in high-energy regions and may be tried when an "
        "equidistant path skips the reaction-coordinate region near the HEI."
    ),
)
@click.option("--max-cycles-gsm", type=click.IntRange(min=1), default=None, show_default="300",
              help="Maximum GSM string-optimizer cycles for the MEP stage.")
@click.option("--max-cycles-dmf", type=click.IntRange(min=1), default=None, show_default="300",
              help=("Maximum IPOPT iterations for the DMF MEP stage. This is a solver "
                    "iteration count, not a string-optimizer cycle count."))
@click.option(
    "--climb/--no-climb",
    default=True,
    show_default=True,
    help="Enable climbing image for standard GSM segments (bridge segments always disable climbing).",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="grad",
    show_default=True,
    help="Single-structure optimizer: grad (=LBFGS) or hess (=RFO).",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write GSM/single-optimization trajectories during the run.",
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/CIF/GJF companions based on the input format.",
)
@click.option(
    "--write-hei-mode-cache/--no-write-hei-mode-cache",
    default=True,
    show_default="on",
    hidden=True,
    help="Internal all-workflow switch for HEI path-mode cache emission.",
)
@click.option("-o", "--out-dir", "out_dir", type=str, default=OUT_DIR_PATH_SEARCH, show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau",
    help=(
        "Convergence preset for single-structure optimizations only "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option(
    "--thresh-gsm",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau_loose",
    help=(
        "Convergence preset for the GSM string optimizer "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option(
    "--thresh-dmf",
    type=str,
    default=None,
    show_default="tight",
    help=(
        "IPOPT dual-infeasibility tolerance for the DMF path optimizer: "
        "tight (0.04) | middle (0.10) | loose (0.20) or a positive float. "
        "This is not a Gaussian preset."
    ),
)
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Base YAML configuration file applied before explicit CLI options.",
)
@click.option(
    "--show-config/--no-show-config",
    "show_config",
    default=False,
    show_default=True,
    help="Print resolved configuration and continue execution.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running path search.",
)
@click.option(
    "--preopt/--no-preopt",
    "preopt",
    default=True,
    show_default=True,
    help="If True, run initial single-structure optimizations of inputs."
)
@click.option(
    "--align/--no-align",
    "align",
    default=True,
    show_default=True,
    help=("After preoptimization, align all inputs to the *first* input and match freeze_atoms "
          "using the align_freeze_atoms API. When --align is True and --ref-full-pdb is provided, "
          "the first reference PDB will be used for all pairs in the final merge.")
)
@click.option(
    "--ref-full-pdb",
    "ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Full-size template PDB/mmCIF files in the same reaction order as --input. "
          "With --align, only the *first* provided reference structure is used for all pairs "
          "in the final merge (you may pass just one).")
)
@click.option(
    "--ref-pdb",
    "model_ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Pocket reference PDB/mmCIF files used only for the final full-system merge. "
          "Useful when --input uses XYZ/GJF intermediates but PDB snapshots exist for merging. "
          "Must match the number and order of --input.")
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Experimental, computationally expensive xTB solvent delta correction. Examples: water, methanol, acetonitrile, dmso, thf, toluene. 'none' disables it.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              show_default=True, help="xTB solvent model.")
@add_coord_type_option(choices=("cart", "dlc"))
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    mep_mode: str,
    dmf_backend: str,
    refine_mode: Optional[str],
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links_flag: bool,
    freeze_atoms_text: Optional[str],
    max_nodes: int,
    gsm_param: Optional[str],
    max_cycles_gsm: Optional[int],
    max_cycles_dmf: Optional[int],
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    write_hei_mode_cache: bool,
    out_dir: str,
    thresh: Optional[str],
    thresh_gsm: Optional[str],
    thresh_dmf: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    preopt: bool,
    align: bool,
    ref_pdb_paths: Optional[Sequence[Path]],
    model_ref_pdb_paths: Optional[Sequence[Path]],
    backend: str,
    solvent: str,
    solvent_model: str,
    cli_coord_type: Optional[str],
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    from pdb2reaction.core.utils import current_cli_args, reject_option_like_extra_args

    argv_all = current_cli_args(ctx)
    _claimed_values = collect_option_values(
        argv_all,
        ("-i", "--input", "--ref-full-pdb", "--ref-pdb"),
    )
    reject_option_like_extra_args(
        ctx.args,
        allowed_values=_claimed_values,
        consumed_values=[
            *input_paths,
            *(ref_pdb_paths or ()),
            *(model_ref_pdb_paths or ()),
        ],
    )
    set_convert_file_enabled(convert_files)
    prepared_inputs: List[PreparedInputStructure] = []
    prepared_auxiliary: List[PreparedInputStructure] = []
    global _PRIMARY_GJF_TEMPLATE
    _PRIMARY_GJF_TEMPLATE = None
    command_str = shlex.join(["pdb2reaction", *map(str, argv_all)])
    # Robustly accept both styles for -i/--input, --ref-full-pdb, and --ref-pdb
    i_vals = collect_option_values(argv_all, ("-i", "--input"))
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

    ref_vals = collect_option_values(argv_all, ("--ref-full-pdb",))
    if ref_vals:
        ref_parsed: List[Path] = []
        for tok in ref_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Reference PDB/mmCIF path '{tok}' not found or is a directory. "
                    f"When using '--ref-full-pdb', multiple files may follow a single option."
                )
            ref_parsed.append(p)
        ref_pdb_paths = tuple(ref_parsed)

    model_ref_vals = collect_option_values(argv_all, ("--ref-pdb",))
    if model_ref_vals:
        model_ref_parsed: List[Path] = []
        for tok in model_ref_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Pocket reference PDB/mmCIF path '{tok}' not found or is a directory. "
                    f"When using '--ref-pdb', multiple files may follow a single option."
                )
            model_ref_parsed.append(p)
        model_ref_pdb_paths = tuple(model_ref_parsed)

    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )
    from pdb2reaction.core.utils import resolve_configured_charge_spin
    charge, spin = resolve_configured_charge_spin(
        merged_yaml_cfg, charge=charge, spin=spin, ligand_charge=ligand_charge,
    )

    time_start = time.perf_counter()

    def _prepare_many(paths: Sequence[Path]) -> List[PreparedInputStructure]:
        prepared_batch: List[PreparedInputStructure] = []
        try:
            for path in paths:
                prepared_batch.append(prepare_input_structure(Path(path)))
        except BaseException:
            for prepared in reversed(prepared_batch):
                prepared.cleanup()
            raise
        return prepared_batch

    def _run() -> None:
        nonlocal prepared_inputs, prepared_auxiliary, ref_pdb_paths, model_ref_pdb_paths
        global _PRIMARY_GJF_TEMPLATE
        if len(input_paths) < 2:
            raise click.BadParameter(
                "Provide at least two structures for --input in reaction order (reactant [intermediates ...] product)."
            )

        mep_mode_kind = mep_mode.strip().lower()
        refine_mode_kind = refine_mode.strip().lower() if refine_mode else None

        do_merge = bool(ref_pdb_paths) and len(ref_pdb_paths) > 0
        if do_merge:
            if align:
                pass
            else:
                if len(ref_pdb_paths) != len(input_paths):
                    raise click.BadParameter("--ref-full-pdb must be given for each --input (same count and order). "
                                             "Alternatively, use --align to allow using only the first reference PDB for all pairs.")
            if model_ref_pdb_paths and len(model_ref_pdb_paths) != len(input_paths):
                raise click.BadParameter("--ref-pdb must be given for each --input (same count and order).")

        p_list = [Path(p) for p in input_paths]
        prepared_inputs = _prepare_many(p_list)
        if ref_pdb_paths:
            prepared_full_refs = _prepare_many(ref_pdb_paths)
            prepared_auxiliary.extend(prepared_full_refs)
            ref_pdb_paths = tuple(prepared.source_path for prepared in prepared_full_refs)
        if model_ref_pdb_paths:
            prepared_model_refs = _prepare_many(model_ref_pdb_paths)
            prepared_auxiliary.extend(prepared_model_refs)
            model_ref_pdb_paths = tuple(
                prepared.source_path for prepared in prepared_model_refs
            )
        any_non_gjf = any(not prep.is_gjf for prep in prepared_inputs)
        if _PRIMARY_GJF_TEMPLATE is None:
            _PRIMARY_GJF_TEMPLATE = next((prep.gjf_template for prep in prepared_inputs if prep.gjf_template), None)

        geom_cfg = dict(GEOM_KW_DEFAULT)
        calc_cfg = dict(UMA_CALC_KW)
        dmf_cfg  = fresh_dmf_config()
        gs_cfg   = dict(GS_KW)
        stopt_cfg = dict(STOPT_KW)
        stopt_cfg["out_dir"] = out_dir
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg   = dict(RFO_KW)
        bond_cfg  = dict(BOND_KW)
        search_cfg = dict(SEARCH_KW)
        search_cfg["refine_mode"] = refine_mode_kind

        def _apply_single_opt_yaml_layer(layer_cfg: Dict[str, Any]) -> None:
            apply_single_opt_yaml_layer(
                layer_cfg,
                lbfgs_cfg=lbfgs_cfg,
                rfo_cfg=rfo_cfg,
                opt_base_kw=OPT_BASE_KW,
                deep_update=deep_update,
                apply_yaml_overrides=apply_yaml_overrides,
            )

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",),)),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
            ],
        )
        _apply_single_opt_yaml_layer(config_layer_cfg)

        resolved_charge = charge
        resolved_spin = spin
        resolved_charge, resolved_spin = fill_charge_spin_from_gjf_inputs(
            resolved_charge,
            resolved_spin,
            [prepared.gjf_template for prepared in prepared_inputs],
        )
        if resolved_charge is None and ligand_charge is not None:
            resolved_charge = _derive_charge_from_ligand_charge(
                prepared_inputs[0], ligand_charge, prefix="[path-search]"
            )
        if resolved_charge is None:
            if any_non_gjf:
                for prepared in prepared_inputs:
                    prepared.cleanup()
                raise click.ClickException(
                    "-q/--charge is required unless all inputs are .gjf templates with charge metadata."
                )
            resolved_charge = 0
        if resolved_spin is None:
            resolved_spin = 1

        # resolved_charge / resolved_spin already incorporate explicit CLI
        # -q/-m (set at the top of the block), gjf metadata fill, and
        # --ligand-charge derivation. Assign them directly. An earlier
        # `calc_cfg.get("charge", resolved_charge)` idiom silently returned
        # the UMA_CALC_KW default 0 when -q was not passed, dropping the
        # derived value and bypassing validate_charge_spin_for_prepared.
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        validate_charge_spin_for_prepared(prepared_inputs, calc_cfg["charge"], calc_cfg["spin"])

        if cli_param_overridden(ctx, "workers"):
            calc_cfg["workers"] = int(workers)
        if cli_param_overridden(ctx, "workers_per_node"):
            calc_cfg["workers_per_node"] = int(workers_per_node)
        if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
            geom_cfg["coord_type"] = str(cli_coord_type).lower()
        if cli_param_overridden(ctx, "max_nodes"):
            gs_cfg["max_nodes"] = int(max_nodes)
            search_cfg["max_nodes_segment"] = int(max_nodes)
        if cli_param_overridden(ctx, "gsm_param") and gsm_param is not None:
            gs_cfg["param"] = str(gsm_param).lower()
        # The GSM cycle budget also bounds the fully-grown string; DMF's budget
        # is a separate IPOPT iteration count.
        if cli_param_overridden(ctx, "max_cycles_gsm") and max_cycles_gsm is not None:
            stopt_cfg["max_cycles"] = int(max_cycles_gsm)
            stopt_cfg["stop_in_when_full"] = int(max_cycles_gsm)
        if cli_param_overridden(ctx, "max_cycles_dmf") and max_cycles_dmf is not None:
            dmf_cfg["max_cycles"] = int(max_cycles_dmf)
        if cli_param_overridden(ctx, "dmf_backend"):
            dmf_cfg["backend"] = str(dmf_backend).lower()
        if cli_param_overridden(ctx, "climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if cli_param_overridden(ctx, "dump"):
            stopt_cfg["dump"] = bool(dump)
            lbfgs_cfg["dump"] = bool(dump)
            rfo_cfg["dump"] = bool(dump)
        if cli_param_overridden(ctx, "out_dir"):
            stopt_cfg["out_dir"] = out_dir
            lbfgs_cfg["out_dir"] = out_dir
            rfo_cfg["out_dir"] = out_dir
        if cli_param_overridden(ctx, "thresh") and thresh is not None:
            lbfgs_cfg["thresh"] = str(thresh)
            rfo_cfg["thresh"] = str(thresh)
        if cli_param_overridden(ctx, "thresh_gsm") and thresh_gsm is not None:
            stopt_cfg["thresh"] = str(thresh_gsm)
        if cli_param_overridden(ctx, "thresh_dmf") and thresh_dmf is not None:
            dmf_cfg["tol"] = str(thresh_dmf)

        # Final YAML overrides (highest precedence)
        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",),)),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
            ],
        )
        _apply_single_opt_yaml_layer(override_layer_cfg)

        if mep_mode_kind == "dmf":
            dmf_cycles = optional_positive_int(dmf_cfg.get("max_cycles"), "dmf.max_cycles")
            dmf_cfg["max_cycles"] = dmf_cycles
            try:
                k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))
            except (TypeError, ValueError, OverflowError) as exc:
                raise click.BadParameter(
                    "dmf.k_fix must be finite and non-negative."
                ) from exc
            if not np.isfinite(k_fix) or k_fix < 0.0:
                raise click.BadParameter(
                    "dmf.k_fix must be finite and non-negative."
                )
            dmf_cfg["k_fix"] = k_fix
        else:
            string_cycles = optional_positive_int(
                stopt_cfg.get("max_cycles"), "stopt.max_cycles"
            )
            stopt_cfg["max_cycles"] = string_cycles
            stopt_cfg["stop_in_when_full"] = string_cycles

        # A dormant YAML DMF section does not affect GSM. An explicit CLI
        # tolerance is still validated as user input, regardless of MEP mode.
        if mep_mode_kind == "dmf" or cli_param_overridden(ctx, "thresh_dmf"):
            resolve_dmf_solve_tol(dmf_cfg, prefix="[path-search]")

        # Convert 1-based YAML freeze_atoms to 0-based internal
        if geom_cfg.get("freeze_atoms"):
            geom_cfg["freeze_atoms"] = yaml_freeze_to_internal(geom_cfg["freeze_atoms"])
        # Merge CLI --freeze-atoms (already 0-based)
        try:
            freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
        except click.BadParameter:
            raise
        if freeze_atoms_cli:
            merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)

        # An explicit CLI --refine-mode outranks YAML search.refine_mode, which
        # the merge above applied on top of the pre-seeded CLI value.
        if cli_param_overridden(ctx, "refine_mode") and refine_mode is not None:
            search_cfg["refine_mode"] = refine_mode.strip().lower()
        refine_mode_kind = search_cfg.get("refine_mode")
        if refine_mode_kind is None:
            refine_mode_kind = "peak" if mep_mode_kind == "gsm" else "minima"
        else:
            refine_mode_kind = str(refine_mode_kind).strip().lower()
            if refine_mode_kind not in {"peak", "minima"}:
                raise click.BadParameter(f"Unknown --refine-mode '{refine_mode_kind}'.")
        search_cfg["refine_mode"] = refine_mode_kind

        opt_kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=OPT_MODE_ALIASES,
            allowed_hint="grad|hess",
        )
        if opt_kind == "lbfgs":
            single_opt_kind = "lbfgs"
            single_opt_cfg = lbfgs_cfg
        else:
            single_opt_kind = "rfo"
            single_opt_cfg = rfo_cfg
        if preopt:
            preopt_cycles = optional_positive_int(
                single_opt_cfg.get("max_cycles"), "preopt max_cycles"
            )
            single_opt_cfg["max_cycles"] = preopt_cycles

        out_dir_path = Path(stopt_cfg["out_dir"]).resolve()
        # Resolve the effective calculator settings before any display: these
        # calls only normalize the config mapping and instantiate no calculator,
        # so --show-config/--dry-run must not print pre-CLI values.
        if cli_param_overridden(ctx, "backend"):
            calc_cfg["backend"] = backend
        if cli_param_overridden(ctx, "solvent"):
            calc_cfg["solvent"] = solvent
        if cli_param_overridden(ctx, "solvent_model"):
            calc_cfg["solvent_model"] = solvent_model
        from pdb2reaction.backends import apply_backend_model_to_calc_cfg
        # Unconditional: also pops a raw backend_model token from a --config YAML
        # (the helper no-ops when neither the CLI arg nor the YAML names one).
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        from pdb2reaction.backends import apply_calc_file_to_calc_cfg
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        from pdb2reaction.backends import apply_effective_precision
        apply_effective_precision(calc_cfg, precision)
        apply_backend_defaults(calc_cfg)

        if mep_mode_kind == "dmf":
            _validate_dmf_solvent_compatibility(calc_cfg)

        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_geom_for_echo(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_stopt = dict(stopt_cfg)
        echo_stopt["out_dir"] = str(out_dir_path)

        # --show-config/--dry-run exist to print these blocks, so they must not
        # be suppressed by the default verbosity gate.
        requested = bool(show_config or dry_run)
        click.echo(pretty_block("geom", echo_geom, force=requested))
        click.echo(pretty_block("calc", echo_calc, force=requested))
        click.echo(pretty_block("gs",   echo_gs, force=requested))
        click.echo(pretty_block("stopt", echo_stopt, force=requested))
        if mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg, force=requested))
        echo_opt = dict(single_opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)
        echo_opt["out_dir_per_tag"] = f"{out_dir_path}/<tag>_{single_opt_kind}_opt"
        click.echo(pretty_block("opt." + single_opt_kind, echo_opt, force=requested))
        click.echo(pretty_block("bond", bond_cfg, force=requested))
        click.echo(pretty_block("search", search_cfg, force=requested))
        click.echo(
            pretty_block(
                "run_flags",
                {"preopt": bool(preopt), "align": bool(align), "mep_mode": mep_mode_kind},
                force=requested,
            )
        )

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                force=True)
            )

        if dry_run:
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_count": len(p_list),
                        "input_first": str(p_list[0]) if p_list else None,
                        "input_last": str(p_list[-1]) if p_list else None,
                        "output_dir": str(out_dir_path),
                        "mep_mode": mep_mode_kind,
                        "gsm_param": str(gs_cfg.get("param", GS_KW["param"])),
                        "refine_mode": refine_mode_kind,
                        "opt_mode": ("grad" if opt_kind == "lbfgs" else "hess"),
                        "preopt": bool(preopt),
                        "align": bool(align),
                        "freeze_links": bool(freeze_links_flag),
                        "convert_files": bool(convert_files),
                        "max_depth": int(search_cfg.get("max_depth", SEARCH_KW["max_depth"])),
                        "max_nodes_segment": int(search_cfg.get("max_nodes_segment", gs_cfg.get("max_nodes", 0))),
                        "will_run_path_search": True,
                        "will_write_summary": True,
                    },
                    force=True,
                )
            )
            emit_dry_run_complete()
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)
        for name in ("mep.pdb", "mep.cif", "mep_plot.png"):
            (out_dir_path / name).unlink(missing_ok=True)

        geoms = load_prepared_geometries(
            prepared_inputs,
            coord_type=geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"]),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )
        main_prepared = prepared_inputs[0] if prepared_inputs else None

        if geoms:
            freeze_union = sorted(
                {int(i) for g in geoms for i in getattr(g, "freeze_atoms", [])}
            )
            calc_cfg["freeze_atoms"] = freeze_union

        shared_calc = create_calculator(**calc_cfg)
        for g in geoms:
            g.set_calculator(shared_calc)

        # If any input is PDB, treat as "PDB input" for final output handling.
        # Fall back to --ref-pdb (model_ref_pdb_paths) when inputs are XYZ.
        ref_pdb_for_segments: Optional[Path] = None
        for prepared in prepared_inputs:
            if prepared.source_path.suffix.lower() == ".pdb":
                ref_pdb_for_segments = prepared.source_path.resolve()
                break
        if ref_pdb_for_segments is None and model_ref_pdb_paths:
            for p in model_ref_pdb_paths:
                if Path(p).suffix.lower() == ".pdb":
                    ref_pdb_for_segments = Path(p).resolve()
                    break

        preopt_outcomes: List[Any] = []
        if preopt:
            from pdb2reaction.workflows._outcomes import make_leaf as _mk_leaf

            new_geoms: List[Any] = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                g_opt, preopt_converged = _optimize_single(
                    g,
                    shared_calc,
                    single_opt_kind,
                    single_opt_cfg,
                    out_dir_path,
                    tag=tag,
                    prepared_input=prepared_inputs[i] if i < len(prepared_inputs) else main_prepared,
                    ref_pdb=ref_pdb_for_segments,
                )
                new_geoms.append(g_opt)
                preopt_outcomes.append(
                    _mk_leaf(
                        "path",
                        f"preopt_endpoint_{i}",
                        executed=True,
                        converged=preopt_converged,
                    )
                )
            geoms = new_geoms
        else:
            click.echo("[init] Skipping endpoint preoptimization as requested by --no-preopt.")

        # Align all inputs to the first structure, guided by freeze constraints, when requested
        align_thresh = str(single_opt_cfg.get("thresh", "gau"))
        if align:
            try:
                emit("\n====== Aligning all inputs to the first structure (freeze-guided scan + relaxation) ======\n", narrative=True)
                alignment_results = align_and_refine_sequence_inplace(
                    geoms,
                    thresh=align_thresh,
                    shared_calc=shared_calc,
                    out_dir=out_dir_path / "align_refine",
                    verbose=True,
                )
                failed_pairs = alignment_failed_pair_indices(alignment_results)
                if failed_pairs:
                    raise click.ClickException(
                        "Input alignment did not converge for pair(s): "
                        + ", ".join(str(index) for index in failed_pairs)
                    )
                click.echo("[align] Completed input alignment.")
            except Exception as e:
                raise click.ClickException(f"Input alignment failed: {e}") from e
        else:
            click.echo("[align] Skipping input alignment as requested by --no-align.")

        _mep_search_start = time.perf_counter()
        emit("\n====== Multistep MEP search (multi-structure) started ======\n", narrative=True)
        seg_counter = [0]

        bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 5))
        gs_bridge_cfg = _gs_cfg_with_overrides(gs_cfg, max_nodes=bridge_max_nodes, climb=False, climb_lanczos=False)

        combined_imgs: List[Any] = []
        combined_Es: List[float] = []
        seg_reports_all: List[SegmentReport] = []
        required_outcomes_all: List[Any] = list(preopt_outcomes)

        for i in range(len(geoms) - 1):
            gA, gB = geoms[i], geoms[i + 1]
            pair_tag = f"pair_{i:02d}"
            emit(f"[stage] Processing pair {i:02d}: image {i} → {i+1}", narrative=True)
            pair_path = _build_multistep_path(
                gA, gB,
                shared_calc,
                geom_cfg, gs_cfg, stopt_cfg,
                single_opt_kind, single_opt_cfg,
                bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                prepared_input=main_prepared,
                depth=0,
                seg_counter=seg_counter,
                branch_tag=pair_tag,
                pair_index=i,
                kink_seq_count=_trailing_kink_count(seg_reports_all),
            )

            if i == 0:
                combined_imgs = list(pair_path.images)
                combined_Es = list(pair_path.energies)
                seg_reports_all.extend(pair_path.segments)
                required_outcomes_all.extend(pair_path.required_outcomes)
            else:
                def _segment_builder_for_pairs(
                    tail_g,
                    head_g,
                    _tag: str,
                    interface_pair_index: Optional[int],
                ) -> CombinedPath:
                    sub = _build_multistep_path(
                        tail_g, head_g,
                        shared_calc,
                        geom_cfg, gs_cfg, stopt_cfg,
                        single_opt_kind, single_opt_cfg,
                        bond_cfg, search_cfg, refine_mode_kind, mep_mode_kind, calc_cfg, dmf_cfg, prepared_inputs,
                        out_dir=out_dir_path,
                        ref_pdb_path=ref_pdb_for_segments,
                        prepared_input=main_prepared,
                        depth=0,
                        seg_counter=seg_counter,
                        branch_tag="B",
                        pair_index=interface_pair_index,
                        kink_seq_count=_trailing_kink_count(seg_reports_all),
                    )
                    required_outcomes_all.extend(sub.required_outcomes)
                    return sub

                parts = [(combined_imgs, combined_Es), (pair_path.images, pair_path.energies)]
                combined_imgs, combined_Es = _stitch_paths(
                    parts=parts,
                    stitch_rmsd_thresh=float(search_cfg.get("stitch_rmsd_thresh", 1e-4)),
                    bridge_rmsd_thresh=float(search_cfg.get("bridge_rmsd_thresh", 1e-4)),
                    shared_calc=shared_calc,
                    gs_cfg=gs_bridge_cfg,
                    stopt_cfg=stopt_cfg,
                    out_dir=out_dir_path,
                    tag=pair_tag,
                    ref_pdb_path=ref_pdb_for_segments,
                    bond_cfg=bond_cfg,
                    segment_builder=_segment_builder_for_pairs,
                    segments_out=seg_reports_all,
                    bridge_pair_index=i,
                    mep_mode_kind=mep_mode_kind,
                    calc_cfg=calc_cfg,
                    max_nodes=bridge_max_nodes,
                    dmf_cfg=dmf_cfg,
                    prepared_inputs=prepared_inputs,
                )
                seg_reports_all.extend(pair_path.segments)
                required_outcomes_all.extend(pair_path.required_outcomes)
            emit(
                f"[stage] Pair {i:02d} done: images={len(pair_path.images)}, "
                f"segments={len(pair_path.segments)}",
                detail=True,
            )

        emit(
            "====== Multistep MEP search (multi-structure) finished "
            f"(pairs={max(len(geoms) - 1, 0)}, segments={len(seg_reports_all)}, "
            f"elapsed={time.perf_counter() - _mep_search_start:.1f}s) ======\n",
            narrative=True,
        )

        combined_all = CombinedPath(
            images=combined_imgs,
            energies=combined_Es,
            segments=seg_reports_all,
            required_outcomes=required_outcomes_all,
        )

        for idx, srep in enumerate(combined_all.segments, 1):
            srep.seg_index = idx
        tag_to_index = {s.tag: int(s.seg_index) for s in combined_all.segments}
        for im in combined_all.images:
            tag = getattr(im, "mep_seg_tag", None)
            if tag and tag in tag_to_index:
                try:
                    setattr(im, "mep_seg_index", int(tag_to_index[tag]))
                except Exception:
                    pass

        # Final MEP output rule:
        # - Always write 'mep_trj.xyz' (XYZ) for intermediate handoff.
        # - If reference topologies are available, also emit mep.pdb and a
        #   public CIF/GJF companion when applicable.
        main_prepared = prepared_inputs[0]
        needs_pdb = ref_pdb_for_segments is not None
        needs_gjf = main_prepared.is_gjf

        final_trj = out_dir_path / "mep_trj.xyz"
        write_xyz_trj_with_energy(combined_all.images, combined_all.energies, final_trj)
        emit(f"[write] Wrote '{final_trj}'.", detail=True)
        try:
            run_trj2fig(final_trj, [out_dir_path / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            close_matplotlib_figures()
            emit(f"[plot] Saved energy plot → '{out_dir_path / 'mep_plot.png'}'", detail=True)
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot final energy: {e}", err=True)

        if needs_pdb or needs_gjf:
            try:
                did_convert = convert_xyz_like_outputs(
                    final_trj,
                    main_prepared,
                    ref_pdb_path=ref_pdb_for_segments,
                    out_pdb_path=out_dir_path / "mep.pdb" if needs_pdb else None,
                    out_gjf_path=out_dir_path / "mep.gjf" if needs_gjf else None,
                )
                if did_convert:
                    emit("[convert] Wrote final MEP outputs.", detail=True)
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final MEP outputs: {e}", err=True)

        # Pocket‑only per‑segment trajectories & HEIs (bond‑change segments only)
        try:
            frame_seg_indices: List[int] = [int(getattr(im, "mep_seg_index", 0) or 0) for im in combined_all.images]
            seg_to_frames: Dict[int, List[int]] = {}
            for ii, sidx in enumerate(frame_seg_indices):
                if sidx <= 0:
                    continue
                seg_to_frames.setdefault(int(sidx), []).append(ii)

            for s in combined_all.segments:
                seg_idx = int(s.seg_index)
                idxs = seg_to_frames.get(seg_idx, [])
                if not idxs:
                    continue

                if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    seg_imgs = [combined_all.images[j] for j in idxs]
                    seg_Es = [combined_all.energies[j] for j in idxs]
                    seg_trj = out_dir_path / f"mep_seg_{seg_idx:02d}_trj.xyz"
                    write_xyz_trj_with_energy(seg_imgs, seg_Es, seg_trj)
                    emit(f"[write] Wrote per-segment active site model trajectory → '{seg_trj}'", detail=True)
                    if needs_pdb or needs_gjf:
                        try:
                            convert_xyz_like_outputs(
                                seg_trj,
                                main_prepared,
                                ref_pdb_path=ref_pdb_for_segments,
                                out_pdb_path=out_dir_path / f"mep_seg_{seg_idx:02d}.pdb" if needs_pdb else None,
                                out_gjf_path=out_dir_path / f"mep_seg_{seg_idx:02d}.gjf" if needs_gjf else None,
                            )
                        except Exception as e:
                            click.echo(
                                f"[convert] WARNING: Failed to convert per-segment trajectory {seg_idx:02d}: {e}",
                                err=True,
                            )

                if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    energies_seg = [combined_all.energies[j] for j in idxs]
                    imax_rel = int(np.argmax(np.array(energies_seg, dtype=float)))
                    imax_abs = idxs[imax_rel]
                    hei_img = combined_all.images[imax_abs]
                    hei_E = [combined_all.energies[imax_abs]]
                    hei_trj = out_dir_path / f"hei_seg_{seg_idx:02d}.xyz"
                    write_xyz_trj_with_energy([hei_img], hei_E, hei_trj)
                    emit(f"[write] Wrote segment HEI (active site model) → '{hei_trj}'", detail=True)
                    if write_hei_mode_cache and len(idxs) >= 2:
                        cache_path = out_dir_path / f"hei_mode_seg_{seg_idx:02d}.npz"
                        legacy_mode_path = out_dir_path / f"hei_mode_seg_{seg_idx:02d}.txt"
                        cache = write_path_mode_cache(
                            cache_path,
                            [image.cart_coords for image in seg_imgs],
                            imax_rel,
                            energies=energies_seg,
                            trajectory_path=seg_trj,
                            hei_path=hei_trj,
                            atom_numbers=getattr(hei_img, "atomic_numbers", None),
                            primary_text_path=legacy_mode_path,
                            source="path-search",
                        )
                        if cache is not None:
                            emit(
                                "[write] Wrote HEI path-mode candidates "
                                f"({len(cache.labels)}; CPU/file cache) → '{cache.path}'",
                                detail=True,
                            )
                    if needs_pdb or needs_gjf:
                        try:
                            convert_xyz_like_outputs(
                                hei_trj,
                                main_prepared,
                                ref_pdb_path=ref_pdb_for_segments,
                                out_pdb_path=out_dir_path / f"hei_seg_{seg_idx:02d}.pdb" if needs_pdb else None,
                                out_gjf_path=out_dir_path / f"hei_seg_{seg_idx:02d}.gjf" if needs_gjf else None,
                            )
                        except Exception as e:
                            click.echo(
                                f"[convert] WARNING: Failed to convert HEI for segment {seg_idx:02d}: {e}",
                                err=True,
                            )
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to emit per-segment active site model outputs: {e}", err=True)

        if do_merge:
            _merge_start = time.perf_counter()
            emit("\n====== Full-system merge (active site model → templates) started ======\n", narrative=True)
            # With --align, use only the first reference PDB for all pairs (replicate it).
            if align:
                if not ref_pdb_paths or len(ref_pdb_paths) < 1:
                    raise click.BadParameter("--ref-full-pdb must provide at least one file when performing final merge with --align.")
                first_ref = Path(ref_pdb_paths[0])
                ref_list_for_merge = [first_ref for _ in input_paths]
            else:
                ref_list_for_merge = [Path(p) for p in ref_pdb_paths]

            _merge_final_and_write(
                final_images=list(combined_all.images),
                model_inputs=[Path(p) for p in input_paths],
                ref_pdbs=ref_list_for_merge,
                segments=combined_all.segments,
                out_dir=out_dir_path,
                model_ref_pdbs=[Path(p) for p in model_ref_pdb_paths] if model_ref_pdb_paths else None,
            )
            emit(
                f"====== Full-system merge finished (elapsed={time.perf_counter() - _merge_start:.1f}s) ======\n",
                narrative=True,
            )

        overall_changed: Optional[bool] = None
        overall_summary = ""
        try:
            overall_changed, overall_summary = has_bond_change(combined_all.images[0], combined_all.images[-1], bond_cfg)
        except Exception as exc:
            logger.debug("Failed to evaluate overall bond changes: %s", exc)
            click.echo(
                f"[overall] WARNING: Failed to evaluate covalent-bond changes: {exc}",
                err=True,
            )

        emit("\n====== MEP summary started ======\n", narrative=True)

        emit("[overall] Covalent-bond changes between first and last image:", narrative=True)
        if overall_changed is None:
            emit("  (covalent-bond change detection unavailable)", narrative=True)
        elif overall_changed and overall_summary.strip():
            emit(textwrap.indent(overall_summary.strip(), prefix="  "), narrative=True)
        else:
            emit("  (no covalent changes detected)", narrative=True)

        if combined_all.segments:
            emit("[segments] Along the final MEP order (ΔE‡, ΔE). Bridges are shown between connected segments:", narrative=True)
            for i, seg in enumerate(combined_all.segments, 1):
                kind_label = seg.kind.upper()
                if seg.kind == "kink":
                    emit(
                        f"  [{i:02d}] ({kind_label}) {seg.tag}  |  "
                        f"non-MEP connector, ΔE = {seg.delta_kcal:.2f} kcal/mol",
                        narrative=True,
                    )
                else:
                    emit(
                        f"  [{i:02d}] ({kind_label}) {seg.tag}  |  "
                        f"ΔE‡ = {seg.barrier_kcal:.2f} kcal/mol,  "
                        f"ΔE = {seg.delta_kcal:.2f} kcal/mol",
                        narrative=True,
                    )
                if seg.kind != "bridge" and seg.summary.strip():
                    emit(
                        textwrap.indent(seg.summary.strip(), prefix="      "),
                        narrative=True,
                    )
        else:
            emit("[segments] (no segment reports)", narrative=True)

        emit("====== MEP summary finished ======\n", narrative=True)

        diagram_payload: Optional[Dict[str, Any]] = None
        png_path = out_dir_path / "energy_diagram_MEP.png"
        png_path.unlink(missing_ok=True)
        try:
            # Map frames to segment indices (for anchoring R/P energies in Hartree)
            frame_seg_indices: List[int] = [int(getattr(im, "mep_seg_index", 0) or 0) for im in combined_all.images]
            seg_to_frames: Dict[int, List[int]] = {}
            for ii, sidx in enumerate(frame_seg_indices):
                if sidx <= 0:
                    continue
                seg_to_frames.setdefault(int(sidx), []).append(ii)

            # Bond‑change segments in the final MEP order
            bc_segments_in_order: List[SegmentReport] = [
                s for s in combined_all.segments
                if (s.kind == "seg" and s.summary and s.summary.strip() != "(no covalent changes detected)")
            ]

            # Determine which frames to use as R and P for anchoring energies in au
            start_idx_for_diag = 0
            end_idx_for_diag = len(combined_all.energies) - 1
            if bc_segments_in_order:
                first_bc = bc_segments_in_order[0]
                last_bc = bc_segments_in_order[-1]
                idxs_first_bc = seg_to_frames.get(int(first_bc.seg_index), [])
                idxs_last_bc = seg_to_frames.get(int(last_bc.seg_index), [])
                if idxs_first_bc:
                    start_idx_for_diag = int(idxs_first_bc[0])
                if idxs_last_bc:
                    end_idx_for_diag = int(idxs_last_bc[-1])

            E0_au = float(combined_all.energies[start_idx_for_diag])
            EP_au = float(combined_all.energies[end_idx_for_diag])

            # Build TS groups and compressed state energies purely from segment-level ΔE / ΔE‡
            ts_groups: List[Dict[str, Any]] = []
            ts_count = 0
            current: Optional[Dict[str, Any]] = None
            # Energy of the current state relative to R (in kcal/mol)
            E_current_kcal = 0.0

            def _is_bond_change_seg(seg: SegmentReport) -> bool:
                return (
                    seg.kind == "seg"
                    and bool(seg.summary)
                    and seg.summary.strip() != "(no covalent changes detected)"
                )

            for s in combined_all.segments:
                if _is_bond_change_seg(s):
                    ts_count += 1
                    _barrier = getattr(s, "barrier_kcal", None)
                    barrier_kcal = float(_barrier) if _barrier is not None else float("nan")
                    delta_kcal = float(getattr(s, "delta_kcal", float("nan")))
                    if not np.isfinite(barrier_kcal):
                        barrier_kcal = 0.0
                    if not np.isfinite(delta_kcal):
                        delta_kcal = 0.0

                    ts_e = E_current_kcal + barrier_kcal
                    first_im_e = E_current_kcal + delta_kcal

                    current = {
                        "ts_label": f"TS{ts_count}",
                        "ts_energy": ts_e,
                        "first_im_energy": first_im_e,
                        "tail_im_energy": first_im_e,
                        "has_extra": False,
                        "index": ts_count,
                        "bridge_peaks": [],
                    }
                    ts_groups.append(current)
                    E_current_kcal = first_im_e
                else:
                    # Segments without covalent changes (bridge or kinks) belong to
                    # the current TS block if one exists.
                    if current is None:
                        # Before the first bond‑change segment: accumulate net ΔE if available.
                        delta_kcal = float(getattr(s, "delta_kcal", float("nan")))
                        if np.isfinite(delta_kcal):
                            E_current_kcal += delta_kcal
                        continue

                    delta_kcal = float(getattr(s, "delta_kcal", float("nan")))
                    _barrier = getattr(s, "barrier_kcal", None)
                    barrier_kcal = float(_barrier) if _barrier is not None else float("nan")

                    if s.kind == "bridge":
                        if np.isfinite(barrier_kcal) and barrier_kcal > 1.0e-3:
                            peak_e = E_current_kcal + barrier_kcal
                            peaks = current.setdefault("bridge_peaks", [])
                            suffix = "" if len(peaks) == 0 else f"_{len(peaks) + 1}"
                            peak_label = f"IM{current['index']}_TS{suffix}"
                            peaks.append({"label": peak_label, "energy": peak_e})
                            # Log the corresponding peak in au for debugging
                            peak_e_au = E0_au + peak_e / AU2KCALPERMOL
                            click.echo(
                                "    [bridge] Recorded diagram-only TS peak "
                                f"{peak_label} at {peak_e_au:.6f} au "
                                "(from segment-level ΔE‡; bridge segments skip tsopt/thermo/DFT)."
                            )

                    if np.isfinite(delta_kcal):
                        E_current_kcal += delta_kcal
                        current["tail_im_energy"] = E_current_kcal
                        current["has_extra"] = True

            # Assemble labels and energies in kcal/mol
            labels: List[str]
            energies_kcal: List[float]
            chain_tokens: List[str]

            if not ts_groups:
                # No bond‑change segments: simple R→P diagram
                labels = ["R", "P"]
                energies_kcal = [0.0, (EP_au - E0_au) * AU2KCALPERMOL]
                chain_tokens = ["R", "-->", "P"]
            else:
                labels = ["R"]
                energies_kcal = [0.0]
                chain_tokens = ["R"]

                for i, g in enumerate(ts_groups, start=1):
                    last_group = (i == len(ts_groups))

                    labels.append(g["ts_label"])
                    energies_kcal.append(float(g["ts_energy"]))
                    chain_tokens.extend(["-->", g["ts_label"]])

                    if last_group:
                        continue

                    labels.append(f"IM{i}_1")
                    energies_kcal.append(float(g["first_im_energy"]))
                    chain_tokens.extend(["-->", f"IM{i}_1"])

                    for bp in g.get("bridge_peaks", []):
                        labels.append(str(bp["label"]))
                        energies_kcal.append(float(bp["energy"]))
                        chain_tokens.extend(["-->", str(bp["label"])])

                    if g["has_extra"]:
                        labels.append(f"IM{i}_2")
                        energies_kcal.append(float(g["tail_im_energy"]))
                        chain_tokens.extend(["-|-->", f"IM{i}_2"])

                labels.append("P")
                energies_kcal.append(E_current_kcal)
                chain_tokens.extend(["-->", "P"])

            # Convert back to Hartree (au) for completeness in the summary
            energies_au: List[float] = [E0_au + ek / AU2KCALPERMOL for ek in energies_kcal]

            diagram_payload = {
                "name": "energy_diagram_MEP",
                "labels": labels,
                "energies_kcal": energies_kcal,
                "ylabel": "ΔE (kcal/mol)",
                "energies_au": energies_au,
            }

            labels_repr = "[" + ", ".join(f'"{lab}"' for lab in labels) + "]"
            energies_repr = "[" + ", ".join(f"{val:.6f}" for val in energies_kcal) + "]"
            click.echo(f"[diagram] build_energy_diagram.labels = {labels_repr}")
            click.echo(f"[diagram] build_energy_diagram.energies_kcal = {energies_repr}")

            title_note = "(MEP; all segments)"
            fig = build_energy_diagram(
                energies=energies_kcal,
                labels=labels,
                ylabel="ΔE (kcal/mol)",
                baseline=True,
                showgrid=False,
            )
            fig.update_layout(title=title_note)

            try:
                fig.write_image(str(png_path), scale=2)
                diagram_payload["image"] = str(png_path)
                emit(f"[diagram] Wrote energy diagram (PNG) → '{png_path}'", detail=True)
            except Exception as e:
                click.echo(f"[diagram] NOTE: PNG export skipped (install 'kaleido' to enable): {e}")

            chain_text = " ".join(chain_tokens)
            emit(f"[diagram] State label sequence: {chain_text}", detail=True)

        except Exception as e:
            click.echo(f"[diagram] WARNING: Failed to build energy diagram: {e}", err=True)

        frame_ranges = _frame_ranges_by_segment(combined_all.images)
        summary = {
            "out_dir": str(out_dir_path),
            "n_images": len(combined_all.images),
            "n_segments": len(combined_all.segments),
            "segments": [
                {
                    "index": int(s.seg_index),
                    "tag": s.tag,
                    "kind": s.kind,
                    # Record the segment's own StringOptimizer convergence
                    # (tri-state) so a downstream consumer / the `all`-pipeline
                    # no-tsopt aggregate can gate on real per-segment convergence
                    # instead of assuming a segment converged.
                    "converged": s.converged,
                    "barrier_kcal": (
                        float(s.barrier_kcal)
                        if s.barrier_kcal is not None
                        else None
                    ),
                    "delta_kcal": float(s.delta_kcal),
                    "bond_changes": (
                        _bond_changes_block(s.summary) if (s.kind != "bridge") else ""
                    ),
                    **frame_ranges.get(int(s.seg_index), {}),
                } for s in combined_all.segments
            ],
        }
        if diagram_payload is not None:
            summary["energy_diagrams"] = [diagram_payload]

        # Enrich summary with metadata (inline to avoid circular import from all.py)
        try:
            from pdb2reaction._version import __version__
        except Exception:
            __version__ = "unknown"
        summary["pdb2reaction_version"] = __version__
        from pdb2reaction.core.utils import RESULT_JSON_SCHEMA_VERSION
        summary["schema_version"] = RESULT_JSON_SCHEMA_VERSION
        summary["pipeline_mode"] = "path-search"
        from pdb2reaction.workflows._outcomes import (
            aggregate_workflow_truth as _agg_truth,
            attach_outcomes as _attach_outcomes,
        )
        _raw_arts = [
            _f for _f in ("mep.pdb", "mep.cif", "mep_plot.png", "energy_diagram_MEP.png")
            if (out_dir_path / _f).exists()
        ]
        _path_leaves, _path_expected = _path_leaves_and_expected(
            combined_all.segments,
            raw_artifacts=_raw_arts,
            required_outcomes=combined_all.required_outcomes,
        )
        _path_truth = _agg_truth(_path_leaves, _path_expected)
        # Legacy byte-compat: `status` stays the diagram-based value it always
        # had — "success" when an energy diagram was produced, else "partial".
        # The endpoint-HEI demotion applies when there
        # is no reactive segment (segments=[]), a raw R/P diagram is a diagnostic
        # artifact only and must not read "success". All other convergence truth
        # is carried by the additive scientific_status, not this legacy field.
        _reactive_segs = [
            s for s in combined_all.segments if getattr(s, "kind", "seg") != "bridge"
        ]
        _legacy_status = "success" if summary.get("energy_diagrams") else "partial"
        if _legacy_status == "success" and not _reactive_segs:
            _legacy_status = "partial"  # endpoint-HEI demotion
        summary["status"] = _legacy_status
        _attach_outcomes(summary, truth=_path_truth, stage_outcomes=_path_leaves)
        from pdb2reaction.core.utils import calculator_provenance

        summary.update(calculator_provenance(calc_cfg))
        summary["charge"] = calc_cfg.get("charge")
        summary["spin"] = calc_cfg.get("spin")
        summary["command"] = command_str
        summary["references"] = method_references(
            {
                "pipeline_mode": "path-search",
                "mep_mode": mep_mode,
                "path_opt_mode": opt_mode,
                "dmf_correlated": bool(dmf_cfg.get("correlated", False)),
                "mlip_backend": summary.get("mlip_backend"),
                "mlip_model": summary.get("mlip_model"),
                "mlip_task": summary.get("mlip_task"),
            }
        )
        try:
            from pdb2reaction.core.utils import _collect_environment_info
            summary["environment"] = _collect_environment_info()
        except Exception:
            pass

        summary = apply_current_run_id(summary)
        commit_json(out_dir_path / "summary.json", summary)
        emit(f"[write] Wrote '{out_dir_path / 'summary.json'}'.", detail=True)

        summary_payload_for_citations: Dict[str, Any] = {}
        try:
            # `freeze_atoms_for_log` is also assigned later in this function
            # (so Python treats it as local to _run); if the inner sorted()
            # raised before its result was bound, the L2755 read below hit
            # an UnboundLocalError instead of falling back to the empty list
            # from the cli() scope. Pre-initialize to ensure the name exists.
            freeze_atoms_for_log = []
            try:
                freeze_atoms_for_log = sorted(
                    {
                        int(i)
                        for g in getattr(combined_all, "images", [])
                        for i in getattr(g, "freeze_atoms", [])
                    }
                )
            except (AttributeError, TypeError, ValueError) as exc:
                logger.debug("Failed to collect freeze_atoms for summary log: %s", exc)

            diag_for_log: Dict[str, Any] = {}
            if diagram_payload is not None:
                diag_for_log = dict(diagram_payload)
            mep_info = {
                "n_images": len(combined_all.images),
                "n_segments": len(combined_all.segments),
                "traj_pdb": str(out_dir_path / "mep.pdb") if (out_dir_path / "mep.pdb").exists() else None,
                "traj_cif": str(out_dir_path / "mep.cif") if (out_dir_path / "mep.cif").exists() else None,
                "mep_plot": str(out_dir_path / "mep_plot.png") if (out_dir_path / "mep_plot.png").exists() else None,
                "diagram": diag_for_log,
            }
            _provenance = calculator_provenance(calc_cfg)
            mlip_backend = str(_provenance["mlip_backend"])
            mlip_model = _provenance["mlip_model"]
            summary_payload = {
                "root_out_dir": str(out_dir_path),
                "path_dir": str(out_dir_path),
                "path_module_dir": "path_search",
                "pipeline_mode": "path-search",
                "refine_path": True,
                "tsopt": False,
                "thermo": False,
                "dft": False,
                "opt_mode": opt_mode,
                "mep_mode": mep_mode,
                "dmf_correlated": bool(dmf_cfg.get("correlated", False)),
                "mlip_backend": mlip_backend,
                "mlip_model": mlip_model,
                "mlip_model_label": _provenance["mlip_model_label"],
                "mlip_task": _provenance["mlip_task"],
                "mlip_precision": _provenance["mlip_precision"],
                "status": summary.get("status"),
                "status_reasons": summary.get("status_reasons", []),
                "execution_status": summary.get("execution_status"),
                "scientific_status": summary.get("scientific_status"),
                "scientific_status_reasons": summary.get(
                    "scientific_status_reasons", []
                ),
                "command": command_str,
                "charge": calc_cfg.get("charge"),
                "spin": calc_cfg.get("spin"),
                "freeze_atoms": freeze_atoms_for_log,
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {},
            }
            summary_payload_for_citations = summary_payload
            write_summary_log(out_dir_path / "summary.log", summary_payload)
            emit(f"[write] Wrote '{out_dir_path / 'summary.log'}'.", detail=True)
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

        # summary.md and key_* outputs are disabled.
        from pdb2reaction.core.utils import is_child_mode

        if summary_payload_for_citations and not is_child_mode():
            emit_method_citations(summary_payload_for_citations)
        emit(
            format_elapsed("[time] Elapsed Time for Path Search", time_start),
            narrative=True,
        )

    out_dir_path = Path(out_dir).resolve()
    try:
        run_cli(
            _run,
            label="path search",
            zero_step_exc=ZeroStepLength,
            zero_step_msg="ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).",
            opt_exc=OptimizationError,
            opt_msg="ERROR: Path search failed — {exc}",
            out_dir=out_dir_path,
            command="path-search",
            time_start=time_start,
        )
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
        for prepared in prepared_auxiliary:
            prepared.cleanup()
        _PRIMARY_GJF_TEMPLATE = None
