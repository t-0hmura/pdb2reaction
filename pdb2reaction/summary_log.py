"""
Human-friendly summary log writer used by ``path_search`` and ``all``.

The goal is to provide a compact, readable ``summary.log`` alongside the
machine-readable ``summary.yaml``. The log aggregates MEP details, segment
barriers, post-processing energies, and key output paths in a single place.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from pysisyphus.constants import AU2KCALPERMOL


def _as_path_str(path: Optional[Path]) -> str:
    return str(path) if path else "(not available)"


def _format_energy_rows(
    labels: Sequence[str],
    energies_au: Optional[Sequence[Optional[float]]],
    energies_kcal: Optional[Sequence[Optional[float]]],
) -> List[str]:
    rows: List[str] = []
    energies_au = list(energies_au or [])
    energies_kcal = list(energies_kcal or [])
    base_e = energies_au[0] if energies_au else None

    for i, lab in enumerate(labels):
        abs_e = energies_au[i] if i < len(energies_au) else None
        rel_e = energies_kcal[i] if i < len(energies_kcal) else None
        if rel_e is None and abs_e is not None and base_e is not None:
            rel_e = (abs_e - base_e) * AU2KCALPERMOL

        # Abs in Eh, and Rel column width tuned to match sample layout
        abs_txt = f"{abs_e:.6f} Eh" if abs_e is not None else "n/a"
        rel_txt = f"{rel_e:14.4f}" if rel_e is not None else "   n/a"
        rows.append(f"        {lab:<8}{rel_txt}    ({abs_txt})")
    return rows


def _format_bond_changes(text: str, indent: int = 6) -> List[str]:
    if not text:
        return ["".rjust(indent) + "(no covalent changes detected)"]
    blocks = [ln.rstrip() for ln in textwrap.dedent(text).splitlines() if ln.strip()]
    return ["".rjust(indent) + ln for ln in blocks]


def _emit_energy_block(
    lines: List[str],
    title: str,
    payload: Optional[Dict[str, Any]],
) -> None:
    if not payload:
        return
    labels: Sequence[str] = payload.get("labels") or ["R", "TS", "P"]
    energies_au = payload.get("energies_au")
    energies_kcal = payload.get("energies_kcal")
    lines.append(f"    -- {title} --")
    # Header uses Eh to match abs_txt
    lines.append("       State    Rel [kcal/mol]            Abs [Eh]")
    lines.extend(_format_energy_rows(labels, energies_au, energies_kcal))

    diagram = payload.get("diagram") or payload.get("image")
    if diagram:
        lines.append(f"       Diagram  : {diagram}")
    structs: Dict[str, Any] = payload.get("structures", {})
    if structs:
        lines.append("       Structures:")
        for key in ("R", "TS", "P"):
            if key in structs:
                lines.append(f"         {key}: {_as_path_str(structs.get(key))}")


def write_summary_log(dest: Path, payload: Dict[str, Any]) -> None:
    """Write a human-friendly summary.log at ``dest`` from a pre-collected payload."""

    root_out = payload.get("root_out_dir") or "-"
    path_module = payload.get("path_module_dir") or "-"
    pipeline_mode = payload.get("pipeline_mode") or "-"
    charge = payload.get("charge")
    spin = payload.get("spin")
    command = payload.get("command") or payload.get("cli_command")

    lines: List[str] = []
    lines.append("========================================================================")
    lines.append("pdb2reaction summary.log")
    lines.append("========================================================================")
    if command:
        lines.append(f"Input            : {command}")
    lines.append(f"Root out_dir       : {root_out}")
    lines.append(f"Path module dir    : {path_module}")
    lines.append(f"Pipeline mode      : {pipeline_mode}")
    lines.append(f"Total charge (ML)  : {charge if charge is not None else '-'}")
    lines.append(f"Multiplicity (2S+1): {spin if spin is not None else '-'}")
    lines.append("")

    mep = payload.get("mep", {}) or {}
    diag = mep.get("diagram") or {}
    lines.append("[1] Global MEP overview")
    lines.append(f"  Number of MEP images : {mep.get('n_images', '-')}")
    lines.append(f"  Number of segments   : {mep.get('n_segments', '-')}")
    if mep.get("traj_pdb"):
        lines.append(f"  MEP trajectory (PDB) : {mep.get('traj_pdb')}")
    if mep.get("mep_plot"):
        lines.append(f"  MEP energy plot      : {mep.get('mep_plot')}")
    lines.append("")
    lines.append("  MEP energy diagram (ΔE, kcal/mol)")
    if diag:
        if diag.get("image"):
            lines.append(f"    Image : {diag.get('image')}")
        lines.append("    State    ΔE [kcal/mol]")
        labels = diag.get("labels", [])
        energies = diag.get("energies_kcal", [])
        for i, lab in enumerate(labels):
            rel = energies[i] if i < len(energies) else None
            rel_txt = f"{rel:9.4f}" if rel is not None else "   n/a"
            lines.append(f"        {lab:<8}{rel_txt}")
    else:
        lines.append("    (no diagram available)")

    segments: Iterable[Dict[str, Any]] = payload.get("segments", []) or []
    lines.append("")
    lines.append("[2] Segment-level MEP summary (UMA path)")
    if segments:
        for seg in segments:
            idx = int(seg.get("index", 0) or 0)
            tag = seg.get("tag", f"seg_{idx:03d}")
            kind = seg.get("kind", "seg")
            lines.append(f"  - Segment {idx:02d} [{kind}]  tag={tag}")
            barrier = seg.get("barrier_kcal")
            delta = seg.get("delta_kcal")
            b_txt = f"{barrier:7.2f}" if barrier is not None else "   n/a"
            d_txt = f"{delta:7.2f}" if delta is not None else "   n/a"
            lines.append(f"      ΔE‡ = {b_txt} kcal/mol,  ΔE = {d_txt} kcal/mol")
            lines.append("      Bond changes:")
            lines.extend(_format_bond_changes(str(seg.get("bond_changes", ""))))
    else:
        lines.append("  (no segment reports)")

    post_segments: Iterable[Dict[str, Any]] = payload.get("post_segments", []) or []
    lines.append("")
    lines.append("[3] Per-segment post-processing (TSOPT / Thermo / DFT)")
    if post_segments:
        for seg in post_segments:
            idx = int(seg.get("index", 0) or 0)
            tag = seg.get("tag", f"seg_{idx:02d}")
            kind = seg.get("kind", "seg")
            lines.append(f"  === Segment {idx:02d} ({kind}) tag={tag} ===")
            if seg.get("post_dir"):
                lines.append(f"    Post-process dir : {seg.get('post_dir')}")
            if seg.get("ts_imag_freq_cm") is not None:
                lines.append(f"    TS imaginary freq: {seg.get('ts_imag_freq_cm'):.1f} cm^-1")
            if seg.get("irc_plot"):
                lines.append(f"    IRC plot         : {seg.get('irc_plot')}")
            if seg.get("irc_traj"):
                lines.append(f"    IRC trajectory   : {seg.get('irc_traj')}")
            _emit_energy_block(lines, "UMA energies (TSOPT+IRC)", seg.get("uma"))
            _emit_energy_block(lines, "UMA Gibbs (thermo)", seg.get("gibbs_uma"))
            _emit_energy_block(lines, "DFT single-point", seg.get("dft"))
            _emit_energy_block(lines, "DFT//UMA Gibbs", seg.get("gibbs_dft_uma"))
    else:
        lines.append("  (no post-processing results)")

    lines.append("")
    lines.append("[4] Energy diagrams (overview)")
    diagrams: Iterable[Dict[str, Any]] = payload.get("energy_diagrams", []) or []
    if diagrams:
        for diag_payload in diagrams:
            name = diag_payload.get("name", "diagram")
            ylabel = diag_payload.get("ylabel", "ΔE (kcal/mol)")
            lines.append(f"  {name}  (ylabel: {ylabel})")
            labels = diag_payload.get("labels", [])
            energies = diag_payload.get("energies_kcal", [])
            lines.append("    State   ΔE [kcal/mol]")
            for i, lab in enumerate(labels):
                rel = energies[i] if i < len(energies) else None
                rel_txt = f"{rel:7.3f}" if rel is not None else "   n/a"
                lines.append(f"        {lab:<8}{rel_txt}")
            if diag_payload.get("image"):
                lines.append(f"    Image : {diag_payload.get('image')}")
    else:
        lines.append("  (no energy diagrams recorded)")

    lines.append("")
    lines.append("[5] Output directories / key files (cheat sheet)")
    lines.append(f"  Root out_dir : {root_out}")
    if payload.get("path_dir"):
        lines.append(f"  Path outputs : {payload.get('path_dir')}")
    key_files: Dict[str, Any] = payload.get("key_files", {}) or {}
    if key_files:
        lines.append("  Key files (root):")
        for name, desc in key_files.items():
            lines.append(f"    {name:<18}: {desc}")

    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
