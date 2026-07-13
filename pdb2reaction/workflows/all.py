"""
End-to-end enzymatic reaction workflow: extraction, MEP search, TS optimization, IRC, and post-processing.

Example:
    pdb2reaction all -i reactant.pdb product.pdb -c 'GPP,SAM' -l 'GPP:-3,SAM:1'

For detailed documentation, see: docs/all.md

Table of contents (top-level definitions; refresh manually after structural edits):
    class _EchoState
    def _echo
    def _echo_section
    def _copy_logged
    def _run_cli_main
    def _append_cli_arg
    def _append_toggle_arg
    def _resolve_override_dir
    def _build_calc_cfg
    def _parse_atom_key_from_line
    def _key_variants
    def _build_variant_occurrence_table
    def _model_key_to_index
    def _read_full_atom_keys_in_file_order
    def _format_atom_key_for_msg
    def _parse_scan_lists_literals
    def _format_scan_stage
    def _convert_scan_lists_to_model_indices
    def _pdb_needs_elem_fix
    def _set_yaml_freeze_atoms
    def _get_freeze_atoms
    def _freeze_atoms_for_log
    def _write_args_yaml_with_freeze_atoms
    def _read_summary
    def _geom_from_angstrom
    def _load_segment_endpoints
    def _save_single_geom_as_pdb_for_tools
    def _copy_structures_to_seg_dir
    def _enrich_summary
    def _json_safe
    def _find_with_suffixes
    def _write_segment_energy_diagram
    def _build_global_segment_labels
    def _merge_irc_trajectories_to_single_plot
    def _optimize_endpoint_geom
    def _run_freq_for_state
    def _read_imaginary_frequency
    def _dft_succeeded
    def _dft_energy_ha
    def _run_dft_for_state
    def _run_dft_sequence
    def _run_tsopt_on_hei
    def _irc_and_match
    def _show_advanced_help
    def _configure_all_help_visibility
    def cli
"""

from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from typing import List, Sequence, Optional, Tuple, Dict, Any

import gc
import json
import torch
import logging
import signal
import sys
import tempfile
import click
from click.core import ParameterSource
import time
import yaml
import numpy as np
import shutil

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength

AtomKey = Tuple[str, str, str, str, str, str]

# Local imports from the package
from pdb2reaction.workflows.extract import extract_api
from pdb2reaction.workflows import path_search as _path_search
from pdb2reaction.workflows import path_opt as _path_opt
from pdb2reaction.workflows import tsopt as _tsopt
from pdb2reaction.workflows import freq as _freq_cli
from pdb2reaction.workflows.align_freeze import align_and_refine_sequence_inplace
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import GEOM_KW_DEFAULT, OUT_DIR_ALL, SEGMENTS_DIRNAME, THRESH_CHOICES, WORK_DIRNAME, UMA_CALC_KW as _UMA_CALC_KW
DEFAULT_COORD_TYPE = GEOM_KW_DEFAULT["coord_type"]
from pdb2reaction.io.trj2fig import run_trj2fig
from pdb2reaction.io.summary import write_summary_log
from pdb2reaction.core.utils import (
    build_energy_diagram,
    collect_option_values,
    collect_single_option_values,
    convert_xyz_like_outputs,
    convert_xyz_to_gjf,
    detect_freeze_links_logged,
    format_elapsed,
    merge_freeze_atom_groups,
    prepare_input_structure,
    normalize_freeze_atoms,
    yaml_freeze_to_internal,
    set_convert_file_enabled,
    resolve_charge_spin,
    validate_charge_spin_at_path,
    load_yaml_dict,
    apply_yaml_overrides,
    load_pdb_atom_metadata,
    _round_charge_with_note,
    apply_ref_pdb_override,
    ensure_dir,
    parse_scan_list_triples,
    close_matplotlib_figures,
    convert_xyz_to_pdb,
    _derive_charge_from_ligand_charge,
    write_xyz_trj_with_energy,
    read_xyz_as_blocks,
    read_xyz_first_last,
    xyz_blocks_first_last,
    set_freeze_atoms_or_warn,
    cli_param_overridden,
    verbose_level,
)
from pdb2reaction.cli.common_options import add_coord_type_option, add_precision_option, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option

logger = logging.getLogger(__name__)

class _EchoState:
    """Encapsulate CLI output state for section-spacing logic.

    Replaces the previous module-level ``_log_started`` global so that state
    is scoped to a single CLI invocation and is easier to test.
    """

    def __init__(self) -> None:
        self._started = False

    def reset(self) -> None:
        self._started = False

    def echo(self, *args, **kwargs) -> None:
        click.echo(*args, **kwargs)
        self._started = True

    def section(self, message: str, **kwargs) -> None:
        # Section banners are the narrative spine of the pipeline log: keep
        # them visible at default verbosity (caller may pass narrative=False
        # for a detail-only sub-banner). The leading blank inherits the same
        # tag so the spacing around a shown banner survives the gate.
        narrative = kwargs.setdefault("narrative", True)
        if self._started:
            click.echo(narrative=narrative)
        click.echo(message, **kwargs)
        self._started = True


_echo_state = _EchoState()


def _echo(*args, **kwargs) -> None:
    """Echo a line with output tracking for section spacing.

    Untagged by default (visible at ``-v 3`` inside ``all``). Use
    ``_echo_detail`` for default ``-v 2`` stage details and ``narrative=True``
    for milestone lines.
    """
    _echo_state.echo(*args, **kwargs)


def _echo_detail(*args, **kwargs) -> None:
    """Echo a level-2 detail line with local output tracking."""
    kwargs.setdefault("detail", True)
    _echo_state.echo(*args, **kwargs)


def _echo_section(message: str, **kwargs) -> None:
    """Echo a section header (narrative) with a leading blank line unless first."""
    _echo_state.section(message, **kwargs)


def _emit_final_summary(out_dir: Path | None, time_start: float) -> None:
    """Print a visual `====== Pipeline summary ======` block + Elapsed line.

    Reads ``summary.json`` if present and lifts the most-asked-for numbers
    (status, rate-limiting barrier, reactive-segment count, output dir) so
    the user sees them at the bottom of the log without scrolling back
    through `[diagram] Wrote ...` / `[time] Elapsed Time for X:` clutter.
    Falls back to just the Elapsed line when summary.json is absent
    (dry-run, early failure, TSOPT-only without aggregation).
    """
    summary: Dict[str, Any] = {}
    if out_dir is not None:
        summary_path = Path(out_dir) / "summary.json"
        if summary_path.exists():
            try:
                _loaded = json.loads(summary_path.read_text(encoding="utf-8"))
                if isinstance(_loaded, dict):
                    summary = _loaded
            except (OSError, json.JSONDecodeError):
                summary = {}
    if summary:
        _echo_section("====== Pipeline summary ======")
        status = summary.get("status")
        if status:
            _echo(f"Status: {status}", narrative=True)
        rls = summary.get("rate_limiting_step")
        if isinstance(rls, dict):
            barrier = rls.get("barrier_kcal")
            seg_idx = rls.get("segment")
            method = rls.get("method", "?")
            if barrier is not None:
                _echo(
                    f"Rate-limiting barrier: {float(barrier):.2f} kcal/mol (segment {seg_idx}, method {method})",
                    narrative=True,
                )
        n_reactive = summary.get("n_segments_reactive")
        if n_reactive is not None:
            _echo(f"Reactive segments: {n_reactive}", narrative=True)
        # Report the pipeline out-dir the user passed (-o), not the stage
        # sub-dir that summary.json happens to record (e.g. <out>/path_search).
        out_dir_show = str(out_dir) if out_dir is not None else summary.get("out_dir")
        if out_dir_show:
            _echo(f"Output dir: {out_dir_show}", narrative=True)
        _echo(narrative=True)
    _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start), narrative=True)


from pdb2reaction.workflows import scan as _scan_cli
from pdb2reaction.domain.add_elem_info import assign_elements as _assign_elem_info
from pdb2reaction.io.pdb_fix import has_altloc as _has_altloc, fix_altloc_file as _fix_altloc
from pdb2reaction.workflows import irc as _irc_cli




def _copy_logged(src: Path, dst: Path, *, label: Optional[str] = None, echo: bool = True) -> bool:
    """Copy files with consistent warning messages; return success."""
    try:
        shutil.copy2(src, dst)
        if echo:
            shown = label or src.name
            _echo_detail(f"[all] Copied {shown} → {dst}")
        return True
    except Exception as e:
        shown = label or src
        _echo(f"[all] WARNING: Failed to copy {shown} to {dst}: {e}", err=True)
        return False


def _move_logged(src: Path, dst: Path, *, label: Optional[str] = None, echo: bool = True) -> bool:
    """Move files with consistent warning messages; return success.

    Used to promote MEP deliverables from the ``_work`` scratch tree up to the
    pipeline root without leaving a duplicate behind. ``shutil.move`` replaces
    an existing destination, so re-runs land cleanly.
    """
    try:
        dst = Path(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        # shutil.move treats an existing destination *directory* as a target to
        # move into; clear any existing destination *file* first so the result
        # is always the intended path (idempotent across re-runs).
        if dst.is_file():
            dst.unlink()
        shutil.move(str(src), str(dst))
        if echo:
            shown = label or src.name
            _echo_detail(f"[all] Moved {shown} → {dst}")
        return True
    except Exception as e:
        shown = label or src
        _echo(f"[all] WARNING: Failed to move {shown} to {dst}: {e}", err=True)
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
    # In-proc subcommand dispatch — flag the child's banner / device echo to
    # stay silent so a 4-stage `all` pipeline doesn't reprint the same lines
    # `pdb2reaction ver. X` / `[calc] Resolved device: cuda` once per stage.
    from pdb2reaction.core.utils import set_child_mode
    set_child_mode(True)
    try:
        sys.argv = ["pdb2reaction", cmd_name] + list(args)
        _echo("")
        cli_obj.main(args=list(args), standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            if on_nonzero == "raise":
                raise click.ClickException(f"[{label}] {cmd_name} exit code {code}.")
            _echo(f"[{label}] WARNING: {cmd_name} exited with code {code}", err=True)
    except Exception as e:
        if on_exception == "raise":
            raise click.ClickException(f"[{label}] {cmd_name} failed: {e}")
        _echo(f"[{label}] WARNING: {cmd_name} failed: {e}", err=True)
    finally:
        sys.argv = saved
        set_child_mode(False)
        # DO NOT INLINE: required between every staged subcommand call (cyclic-ref break before allocator-cache reclaim).
        # Release GPU memory between pipeline stages to prevent OOM.
        # gc.collect() breaks cyclic refs inside torch.nn.Module.
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        _echo("")


def _append_cli_arg(args: List[str], flag: str, value: Any | None) -> None:
    """Append ``flag`` and ``value`` (converted to string) to ``args`` when ``value`` is not ``None``."""
    if value is None:
        return
    if isinstance(value, bool):
        args.extend([flag, "True" if value else "False"])
    else:
        args.extend([flag, str(value)])


def _append_toggle_arg(args: List[str], flag: str, value: bool | None) -> None:
    """Append Click bool-toggle option as ``--flag`` / ``--no-flag`` when value is not ``None``."""
    if value is None:
        return
    if not flag.startswith("--"):
        raise ValueError(f"Toggle flag must start with '--': {flag}")
    base = flag if not flag.startswith("--no-") else f"--{flag[5:]}"
    neg = f"--no-{base[2:]}"
    args.append(base if bool(value) else neg)


def _resolve_override_dir(default: Path, override: Path | None) -> Path:
    """Return ``override`` when provided (respecting absolute paths); otherwise ``default``."""
    if override is None:
        return default
    if override.is_absolute():
        return override
    return default.parent / override


CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)


def _forward_calc_file_argv(child_args: List[str], calc_cfg: Dict[str, Any]) -> bool:
    """Forward --calc-file/--calc-factory to a child stage's argv when the custom
    ML backend is active. Children validate -b against the MLIP choices (which
    exclude 'custom'), so the calc-file is forwarded instead of '-b custom'.
    Returns True when it forwarded, so callers skip the --backend forward."""
    cf = calc_cfg.get("calc_file")
    if not cf:
        return False
    child_args.extend(["--calc-file", str(cf)])
    fac = calc_cfg.get("calc_factory")
    if fac and str(fac) != "get_calculator":
        child_args.extend(["--calc-factory", str(fac)])
    return True


def _build_calc_cfg(
    charge: int,
    spin: int,
    workers: Optional[int] = None,
    workers_per_node: Optional[int] = None,
    yaml_cfg: Optional[Dict[str, Any]] = None,
    backend: Optional[str] = None,
    solvent: Optional[str] = None,
    solvent_model: Optional[str] = None,
) -> Dict[str, Any]:
    """Return a calculator configuration honoring YAML overrides when provided.

    Precedence: defaults → YAML → CLI explicit.
    Pass ``None`` for backend/solvent/solvent_model to skip CLI override.
    """
    cfg: Dict[str, Any] = dict(CALC_KW)
    cfg["charge"] = int(charge)
    cfg["spin"] = int(spin)
    if workers is not None:
        cfg["workers"] = int(workers)
    if workers_per_node is not None:
        cfg["workers_per_node"] = int(workers_per_node)
    # Apply YAML first (overrides defaults)
    if yaml_cfg:
        apply_yaml_overrides(
            yaml_cfg,
            [
                (cfg, (("calc",),)),
            ],
        )
    # CLI explicit values override YAML
    if backend is not None:
        cfg["backend"] = backend
    if solvent is not None:
        cfg["solvent"] = solvent
    if solvent_model is not None:
        cfg["solvent_model"] = solvent_model
    # Apply backend-specific defaults (model, precision, etc.) when switching
    # away from UMA.  Only overwrites keys that still hold UMA default values.
    from pdb2reaction.core.defaults import apply_backend_defaults
    apply_backend_defaults(cfg)
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


def _model_key_to_index(model_pdb: Path) -> Dict[AtomKey, List[int]]:
    """
    Build mapping: structural atom key -> list of model indices (1-based by file order).
    """
    key2idx: Dict[AtomKey, List[int]] = defaultdict(list)
    idx = 0
    try:
        with open(model_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    key = _parse_atom_key_from_line(line)
                    if key is None:
                        continue
                    idx += 1
                    for variant in _key_variants(key):
                        key2idx[variant].append(idx)
    except FileNotFoundError:
        raise click.ClickException(f"[all] Active site model PDB not found: {model_pdb}")
    if not key2idx:
        raise click.ClickException(f"[all] Active site model PDB {model_pdb} has no ATOM/HETATM records.")
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


def _convert_scan_lists_to_model_indices(
    scan_lists_raw: Sequence[str],
    full_input_pdb: Path,
    model_pdb: Path,
) -> List[List[Tuple[int, int, float]]]:
    """
    Convert user-provided atom indices (based on the full input PDB) to model indices.
    Returns the converted stages as lists of (i,j,target) with 1-based model indices.

    Structural keys (chainID, resName, resSeq, iCode, atomName, altLoc) are used instead of serial numbers,
    with per-variant occurrence counts to distinguish atoms that otherwise share the same key.
    """
    if not scan_lists_raw:
        return []

    full_atom_meta = load_pdb_atom_metadata(full_input_pdb)
    stages = _parse_scan_lists_literals(scan_lists_raw, atom_meta=full_atom_meta)

    orig_keys_in_order = _read_full_atom_keys_in_file_order(full_input_pdb)
    key_to_model_idx = _model_key_to_index(model_pdb)
    variant_occ_table = _build_variant_occurrence_table(orig_keys_in_order)

    n_atoms_full = len(orig_keys_in_order)

    def _map_full_index_to_model(idx_one_based: int, stage_idx: int, tuple_idx: int, side_label: str) -> int:
        """
        Convert a 1-based index from the full PDB into the model's 1-based index.
        Fall back in the order: strict match → ignore altloc → ignore iCode → ignore both,
        and use the atom index (occurrence count) when multiple atoms share a structural key.
        """
        key = orig_keys_in_order[idx_one_based - 1]

        variant_occ = variant_occ_table[idx_one_based - 1]
        for variant in _key_variants(key):
            occurrence = variant_occ.get(variant)
            indices = key_to_model_idx.get(variant)
            if occurrence is None or not indices:
                continue
            if occurrence <= len(indices):
                return indices[occurrence - 1]

        msg_key = _format_atom_key_for_msg(key)
        raise click.BadParameter(
            f"--scan-lists #{stage_idx} tuple #{tuple_idx} ({side_label}) references atom index {idx_one_based} "
            f"(key {msg_key}) which is not present in the active site model after extraction. "
            "Increase extraction coverage (e.g., --radius/--radius-het2het, --selected-resn, or set --exclude-backbone False), "
            "or choose atoms that survive in the active site model."
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

            pi = _map_full_index_to_model(idx_i, stage_idx, tuple_idx, "i")
            pj = _map_full_index_to_model(idx_j, stage_idx, tuple_idx, "j")

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
    raw = normalize_freeze_atoms(geom_cfg.get("freeze_atoms"))
    _FREEZE_ATOMS_YAML = yaml_freeze_to_internal(raw)


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
    except Exception as exc:
        logger.debug("Failed to merge freeze atom groups for log: %s", exc)
        return []


def _write_args_yaml_with_freeze_atoms(
    args_yaml: Optional[Path],
    freeze_atoms: Sequence[int],
    coord_type: Optional[str] = None,
    precision: Optional[str] = None,
    backend_model: Optional[str] = None,
) -> Optional[Path]:
    """
    Write ``freeze_atoms`` and (optionally) ``coord_type`` / ``precision`` into a
    YAML config under ``geom`` / ``calc`` and produce a temporary YAML file.
    Returns the new YAML path, or the original ``args_yaml`` when nothing was
    provided.

    ``freeze_atoms`` are 0-based internal indices. They are converted to 1-based
    before writing so that downstream consumers (path_search, etc.) can apply
    ``yaml_freeze_to_internal()`` consistently.

    ``coord_type`` (when set) overrides ``geom.coord_type`` for every child
    stage that reads this YAML — the path is how ``all --coord-type`` is
    propagated uniformly to opt / tsopt / freq / scan / path-opt / path-search
    without per-call argv plumbing. ``precision`` is propagated the same way
    via ``calc.precision`` so ``all --precision fp64`` reaches every child.
    """
    if not freeze_atoms and coord_type is None and precision is None and backend_model is None:
        return args_yaml

    cfg = {} if args_yaml is None else load_yaml_dict(args_yaml)
    if not isinstance(cfg, dict):
        cfg = {}

    geom_cfg = cfg.get("geom")
    if not isinstance(geom_cfg, dict):
        geom_cfg = {}
    geom_cfg = dict(geom_cfg)

    if freeze_atoms:
        # Overwrite (not merge) freeze_atoms: the caller already merged YAML +
        # freeze-links atoms into a single 0-based list.  Convert to 1-based for
        # the YAML file so that yaml_freeze_to_internal() in the consumer works.
        geom_cfg["freeze_atoms"] = sorted(int(i) + 1 for i in freeze_atoms)
    if coord_type is not None:
        geom_cfg["coord_type"] = coord_type
    cfg["geom"] = geom_cfg

    if precision is not None or backend_model is not None:
        calc_cfg = cfg.get("calc")
        if not isinstance(calc_cfg, dict):
            calc_cfg = {}
        calc_cfg = dict(calc_cfg)
        if precision is not None:
            calc_cfg["precision"] = precision
        if backend_model is not None:
            calc_cfg["model"] = backend_model
        cfg["calc"] = calc_cfg

    tmp_dir = Path(tempfile.mkdtemp(prefix="tmp_path_search_"))
    out_path = tmp_dir / "args_freeze_atoms.yaml"
    with out_path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)

    # Register cleanup so the temp directory is removed when the process exits.
    import atexit, shutil as _shutil
    atexit.register(lambda d=tmp_dir: _shutil.rmtree(d, ignore_errors=True))

    return out_path


# ---------- Post-processing helpers ----------


def _read_summary(summary_path: Path) -> List[Dict[str, Any]]:
    """
    Read path_search/summary.json and return segments list (empty if not found).
    """
    try:
        if not summary_path.exists():
            return []
        data = json.loads(summary_path.read_text(encoding="utf-8")) or {}
        segs = data.get("segments", []) or []
        if not isinstance(segs, list):
            return []
        return segs
    except Exception as exc:
        logger.debug("Failed to read summary JSON %s: %s", summary_path, exc)
        return []


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
    ``<path_dir>/<seg_tag>_refine_mep/final_geometries_trj.xyz``.
    If it does not exist, load structures from
    ``<path_dir>/<seg_tag>_mep/final_geometries_trj.xyz``.

    Uses seg_tag (e.g. 'seg_000') and returns (gL_ref, gR_ref).
    """
    base_tag = _path_search._segment_base_id(seg_tag)
    refine_trj = path_dir / f"{base_tag}_refine_mep" / "final_geometries_trj.xyz"
    gsm_trj = path_dir / f"{base_tag}_mep" / "final_geometries_trj.xyz"

    if refine_trj.exists():
        trj_path = refine_trj
    elif gsm_trj.exists():
        trj_path = gsm_trj
    else:
        return None

    # Read first/last frames explicitly to support `_trj.xyz` trajectories.
    elems, c_first, c_last = read_xyz_first_last(trj_path)
    gL_ref = _geom_from_angstrom(elems, c_first, freeze_atoms)
    gR_ref = _geom_from_angstrom(elems, c_last, freeze_atoms)

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
            _echo(
                f"[all] WARNING: failed to convert '{xyz_trj.name}' to PDB for {name}: {e}",
                err=True,
            )

    return xyz_trj


def _copy_structures_to_seg_dir(
    state_structs: Dict[str, Any],
    out_dir: Path,
    seg_idx: int,
    input_suffix: str,
    prepared_input: Any = None,
    ref_pdb_path: Optional[Path] = None,
) -> Path:
    """Copy R/TS/P structures to ``out_dir/seg_XX/`` in the input format.

    Parameters
    ----------
    state_structs : dict
        ``{"R": xyz_path, "TS": xyz_path, "P": xyz_path}``
    input_suffix : str
        Original input file suffix (e.g. ``".pdb"``, ``".xyz"``, ``".gjf"``).
    prepared_input :
        ``PreparedInputStructure`` holding a GJF template (optional).
    ref_pdb_path :
        Reference PDB for coordinate conversion (optional).

    Returns
    -------
    Path to the created ``seg_XX/`` directory.
    """
    seg_dir = out_dir / SEGMENTS_DIRNAME / f"seg_{seg_idx:02d}"
    seg_dir.mkdir(parents=True, exist_ok=True)

    name_map = {"R": "reactant", "TS": "ts", "P": "product"}

    for key, src_xyz in state_structs.items():
        src = Path(src_xyz)
        if not src.exists():
            continue
        dst_name = name_map.get(key, key.lower())

        if input_suffix == ".pdb":
            src_pdb = src.with_suffix(".pdb")
            if src_pdb.exists():
                shutil.copy2(src_pdb, seg_dir / f"{dst_name}.pdb")
            else:
                shutil.copy2(src, seg_dir / f"{dst_name}.xyz")
        elif input_suffix == ".gjf":
            if (
                prepared_input is not None
                and getattr(prepared_input, "gjf_template", None) is not None
            ):
                try:
                    convert_xyz_to_gjf(
                        src,
                        prepared_input.gjf_template,
                        seg_dir / f"{dst_name}.gjf",
                    )
                except Exception:
                    shutil.copy2(src, seg_dir / f"{dst_name}.xyz")
            else:
                shutil.copy2(src, seg_dir / f"{dst_name}.xyz")
        else:
            shutil.copy2(src, seg_dir / f"{dst_name}.xyz")

    return seg_dir


def _enrich_summary(
    summary: dict,
    *,
    version: str,
    pipeline_mode: str,
    mlip_backend: str,
    mlip_model: Optional[str],
    charge: int,
    spin: int,
    command: str = "",
    post_segments: Optional[list] = None,
    config: Optional[dict] = None,
    freeze_atoms: Optional[str] = None,
    out_dir: Optional[Path] = None,
) -> dict:
    """Add machine-readable metadata to summary dict for AI agent consumption.

    The resulting dict is written as summary.json and is intended to be the
    single machine-readable output consumed by MCP tools and AI agents.
    It should contain ALL information present in summary.log.
    """
    from pdb2reaction._version import __version__

    segments = summary.get("segments", [])
    reactive = [s for s in segments if s.get("kind") != "bridge"]
    n_reactive = len(reactive)

    # Determine status
    has_diagrams = bool(summary.get("energy_diagrams"))
    status = "success" if has_diagrams else ("partial" if segments else "failed")

    # Rate-limiting step (highest barrier among reactive segments)
    best_method = None
    rls = None
    if reactive:
        for diag in reversed(summary.get("energy_diagrams", [])):
            name = diag.get("name", "")
            if "G_MLIP" in name and "all" in name:
                best_method = "MLIP_Gibbs"
                break
            elif "MLIP" in name and "all" in name:
                best_method = "MLIP"
                break
        if best_method is None:
            best_method = "MEP"

        # Prefer the refined TSOPT+IRC barrier (post_segments' uma/gibbs_uma) matching
        # best_method; segments[*].barrier_kcal is only the un-refined MEP band and must
        # NOT be reported under an MLIP/Gibbs label (the two can differ by several kcal/mol).
        # An agent reading this field gets the TS/IRC-refined ΔE‡/ΔG‡, with the raw MEP
        # value kept alongside as `mep_barrier_kcal` for transparency.
        method_key = {"MLIP_Gibbs": "gibbs_uma", "MLIP": "uma"}.get(best_method)
        post_by_idx = {
            ps.get("index"): ps
            for ps in (post_segments or [])
            if isinstance(ps, dict)
        }
        max_barrier = -1e9
        for s in reactive:
            idx = s.get("index")
            refined = (post_by_idx.get(idx) or {}).get(method_key) if method_key else None
            if isinstance(refined, dict) and refined.get("barrier_kcal") is not None:
                b = refined.get("barrier_kcal") or 0
                cur_method = best_method
            else:
                # No refined TS/IRC value for this segment: fall back to the MEP band,
                # but label it MEP honestly rather than claiming an MLIP/Gibbs number.
                b = s.get("barrier_kcal", 0) or 0
                cur_method = "MEP"
            if b > max_barrier:
                max_barrier = b
                rls = {
                    "segment": idx,
                    "barrier_kcal": round(b, 2),
                    "method": cur_method,
                    "mep_barrier_kcal": round(s.get("barrier_kcal", 0) or 0, 2),
                }

    # Overall reaction energy from the best all-segment diagram
    overall_rxn_e = None
    for diag in reversed(summary.get("energy_diagrams", [])):
        name = diag.get("name", "")
        if "all" in name:
            energies = diag.get("energies_kcal", [])
            if len(energies) >= 2:
                overall_rxn_e = round(energies[-1] - energies[0], 2)
                break

    # Core metadata
    summary["pdb2reaction_version"] = __version__
    summary["pipeline_mode"] = pipeline_mode
    summary["status"] = status
    summary["mlip_backend"] = mlip_backend
    summary["mlip_model"] = mlip_model
    summary["charge"] = charge
    summary["spin"] = spin
    summary["n_segments_reactive"] = n_reactive
    if rls:
        summary["rate_limiting_step"] = rls
    if overall_rxn_e is not None:
        summary["overall_reaction_energy_kcal"] = overall_rxn_e
    if command:
        summary["command"] = command

    # Pipeline config (mirrors summary.log header)
    if config:
        summary["config"] = config
    if freeze_atoms:
        summary["freeze_atoms"] = freeze_atoms

    # Per-segment post-processing results. Legacy UMA-named machine keys are
    # retained for Supporting Data compatibility; provenance lives in the
    # top-level mlip_backend/mlip_model fields.
    if post_segments:
        summary["post_segments"] = _json_safe(post_segments)

    # Key output file paths for AI agent consumption
    if "out_dir" in summary:
        # Real pipeline root. Fallback to the legacy module_dir.parent for any
        # caller that does not pass out_dir explicitly.
        root = Path(out_dir) if out_dir is not None else Path(summary["out_dir"]).parent
        key_files: Dict[str, Any] = {}
        # Root-level deliverables (MEP products + authored/mirrored summaries live at root)
        for name, desc in [
            ("summary.log", "Human-readable results summary"),
            ("summary.json", "Machine-readable results summary"),
            ("mep_trj.xyz", "Full MEP trajectory"),
            ("mep.pdb", "Full MEP as PDB"),
            ("energy_diagram_MEP.png", "MEP energy plot"),
            ("irc_plot_all.png", "Aggregated IRC plot"),
        ]:
            if (root / name).exists():
                key_files.setdefault(name, desc)
        # Per-segment deliverables under segments/seg_NN/
        seg_parent = root / SEGMENTS_DIRNAME
        if seg_parent.exists():
            for child in sorted(seg_parent.iterdir()):
                if child.is_dir() and child.name.startswith("seg_"):
                    seg_files = [f.name for f in sorted(child.iterdir()) if f.is_file()]
                    key_files[child.name] = {
                        "description": f"Per-segment results for {child.name}",
                        "files": seg_files,
                    }
        if key_files:
            summary["key_output_files"] = key_files

    # Environment info
    try:
        from pdb2reaction.core.utils import _collect_environment_info
        summary.setdefault("environment", _collect_environment_info())
    except Exception:
        pass

    return summary


def _json_safe(obj):
    """Recursively convert Path objects to strings for JSON serialization."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(item) for item in obj]
    return obj


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
        _echo(f"[diagram] Wrote energy diagram → {png.name}")
    except Exception as e:
        _echo(f"[diagram] WARNING: Failed to write energy diagram {png.name}: {e}", err=True)

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


def _merge_irc_trajectories_to_single_plot(
    trj_and_flags: Sequence[Tuple[Path, bool]],
    out_png: Path,
) -> None:
    """
    Build a single IRC plot over all reactive segments using trj2fig.

    Parameters
    ----------
    trj_and_flags : Sequence[Tuple[Path, bool]]
        For each segment: (finished_irc_trj.xyz path, reverse_flag). When reverse_flag is True,
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
            _echo(str(e))
            continue
        if not blocks:
            continue
        if reverse:
            blocks = list(reversed(blocks))
        all_blocks.extend("\n".join(b) for b in blocks)

    if not all_blocks:
        return

    tmp_trj = out_png.with_name(f"{out_png.stem}_trj.xyz")
    ensure_dir(tmp_trj.parent)
    try:
        tmp_trj.write_text("\n".join(all_blocks) + "\n", encoding="utf-8")
    except Exception as e:
        _echo(f"[irc_all] WARNING: Failed to write concatenated IRC trajectory: {e}", err=True)
        return

    try:
        run_trj2fig(tmp_trj, [out_png], unit="kcal", reference="init", reverse_x=False)
        close_matplotlib_figures()
        _echo(f"[irc_all] Wrote aggregated IRC plot → {out_png}")
    except Exception as e:
        _echo(f"[irc_all] WARNING: failed to plot concatenated IRC trajectory: {e}", err=True)
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
        opt_mode_default: "grad"/"lbfgs"/"dimer" or "hess"/"rfo"/"rsprfo".
        out_dir: base directory for the optimization outputs.
        tag: tag prefix for the subdirectory.
        dump: whether to dump optimizer trajectory.
        thresh: optional convergence preset to override defaults.

    Returns:
        (optimized_geometry, final_xyz_path)
    """
    geom.set_calculator(getattr(geom, "calculator", None))
    mode = (opt_mode_default or "hess").lower()
    if mode in ("grad", "lbfgs", "dimer"):
        run_sequence = ("lbfgs",)
    elif mode in ("hess", "rfo", "rsprfo"):
        run_sequence = ("rfo",)
    else:
        run_sequence = ("rfo",)

    final_xyz: Optional[Path] = None
    for sopt_kind in run_sequence:
        if sopt_kind == "lbfgs":
            base_cfg = dict(_path_search.LBFGS_KW)
            OptClass = LBFGS
        else:
            base_cfg = dict(_path_search.RFO_KW)
            OptClass = RFOptimizer

        cfg = dict(base_cfg)
        opt_dir = out_dir / f"{tag}_{sopt_kind}_opt"
        label = sopt_kind.upper()

        ensure_dir(opt_dir)
        cfg["out_dir"] = str(opt_dir)
        cfg["dump"] = bool(dump)
        cfg["max_cycles"] = int(cfg.get("max_cycles", 300))
        if thresh is not None:
            cfg["thresh"] = str(thresh)

        # Seed cached IRC endpoint Hessian for RFO when available
        if sopt_kind == "rfo":
            from pdb2reaction.io.hessian_cache import load as _hess_load
            _cached = _hess_load("irc_endpoint")
            if _cached is not None:
                _echo_detail("[endpoint-opt] Reusing IRC endpoint Hessian for RFO seeding.")
                _active_dofs = _cached.get("active_dofs")
                _h_raw = _cached["hessian"]
                if isinstance(_h_raw, torch.Tensor):
                    _h = _h_raw.clone()
                else:
                    import torch as _torch
                    _h = _torch.as_tensor(_h_raw, dtype=_torch.float64)
                if _active_dofs is not None:
                    geom.within_partial_hessian = {
                        "active_n_dof": len(_active_dofs),
                        "full_n_dof": geom.cart_coords.size,
                        "active_dofs": _active_dofs,
                        "active_atoms": sorted(set(d // 3 for d in _active_dofs)),
                    }
                geom.cart_hessian = _h
                del _h

        _echo_detail(f"[endpoint-opt] Optimizing '{tag}' with {label} → {opt_dir}")
        opt = OptClass(geom, **cfg)
        try:
            opt.run()
        except (OptimizationError, ZeroStepLength) as e:
            _echo(
                f"[endpoint-opt] WARNING: optimization for '{tag}' terminated early ({e}); using last geometry.",
                err=True,
            )

        final_xyz = Path(opt.final_fn) if isinstance(opt.final_fn, (str, Path)) else opt.final_fn

    if final_xyz is None:
        raise click.ClickException(f"[endpoint-opt] No optimized geometry was produced for '{tag}'.")

    g_final = geom_loader(
        final_xyz,
        coord_type=DEFAULT_COORD_TYPE,
        freeze_atoms=getattr(geom, "freeze_atoms", []),
    )
    try:
        g_final.freeze_atoms = np.array(getattr(geom, "freeze_atoms", []), dtype=int)
    except Exception as exc:
        logger.debug("Failed to propagate freeze_atoms to optimized endpoint geometry: %s", exc)
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
        "--out-dir",
        str(fdir),
    ]
    _append_toggle_arg(
        args,
        "--freeze-links",
        bool(
            freeze_use
            and (
                pdb_path.suffix.lower() == ".pdb"
                or (ref_pdb is not None and ref_pdb.suffix.lower() == ".pdb")
            )
        ),
    )
    _append_toggle_arg(args, "--convert-files", bool(convert_files))
    if ref_pdb is not None:
        args.extend(["--ref-pdb", str(ref_pdb)])

    _append_cli_arg(args, "--max-write", overrides.get("max_write"))
    _append_cli_arg(args, "--amplitude-ang", overrides.get("amplitude_ang"))
    _append_cli_arg(args, "--n-frames", overrides.get("n_frames"))
    if overrides.get("sort") is not None:
        args.extend(["--sort", str(overrides.get("sort"))])
    _append_cli_arg(args, "--temperature", overrides.get("temperature"))
    _append_cli_arg(args, "--pressure", overrides.get("pressure"))
    _append_toggle_arg(args, "--dump", bool(dump_use))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        args.extend(["--hessian-calc-mode", str(hess_mode)])

    # Pass backend so freq uses the same MLIP (or the custom calc-file).
    if not _forward_calc_file_argv(args, overrides):
        _freq_backend = overrides.get("backend", "uma")
        if _freq_backend != "uma":
            args.extend(["--backend", _freq_backend])

    # Forward MLIP runtime knobs when the parent CLI explicitly set them; values are
    # injected into overrides by the cli() body before this helper is called.
    _freq_workers = overrides.get("workers")
    if _freq_workers is not None and int(_freq_workers) != 1:
        args.extend(["--workers", str(int(_freq_workers))])
    _freq_workers_per_node = overrides.get("workers_per_node")
    if _freq_workers_per_node is not None and int(_freq_workers_per_node) != 1:
        args.extend(["--workers-per-node", str(int(_freq_workers_per_node))])
    _freq_solvent = overrides.get("solvent")
    if _freq_solvent is not None and str(_freq_solvent).lower() != "none":
        args.extend(["--solvent", str(_freq_solvent)])
    _freq_solvent_model = overrides.get("solvent_model")
    if _freq_solvent_model is not None and str(_freq_solvent_model).lower() != "alpb":
        args.extend(["--solvent-model", str(_freq_solvent_model)])

    if args_yaml is not None:
        args.extend(["--config", str(args_yaml)])
    _run_cli_main("freq", _freq_cli.cli, args, on_nonzero="warn", prefix="freq")
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception as exc:
            logger.debug("Failed to parse thermoanalysis YAML %s: %s", y, exc)
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
            except Exception as exc:
                logger.debug("Skipping unparseable frequency line %r: %s", line.strip(), exc)
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
    except Exception as exc:
        click.echo(
            f"[freq] WARNING: Failed to parse imaginary frequencies from {freq_dir}: {exc}",
            err=True,
        )
        return None




def _dft_succeeded(result: Dict[str, Any]) -> bool:
    """Return True only if DFT converged and produced a valid energy."""
    return bool(result) and not result.get("_dft_failed", True)


def _dft_energy_ha(result: Dict[str, Any]) -> Optional[float]:
    """Extract DFT energy in hartree, or None if DFT failed."""
    if not _dft_succeeded(result):
        return None
    try:
        return float((result.get("energy") or {}).get("hartree"))
    except (TypeError, ValueError):
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
        "--out-dir",
        str(ddir),
    ]
    _append_toggle_arg(args, "--convert-files", bool(convert_files))
    if ref_pdb is not None:
        args.extend(["--ref-pdb", str(ref_pdb)])
    if engine:
        args.extend(["--engine", str(engine)])

    _append_cli_arg(args, "--max-cycle", overrides.get("max_cycle"))
    _append_cli_arg(args, "--conv-tol", overrides.get("conv_tol"))
    _append_cli_arg(args, "--grid-level", overrides.get("grid_level"))

    if args_yaml is not None:
        args.extend(["--config", str(args_yaml)])
    # DO NOT INLINE: torch (UMA/ORB) and gpu4pyscf both link libcusolver at different pinned versions; in-process DFT triggers dynamic-loader clash. Fresh interpreter isolation is mandatory.
    # Run DFT as a real subprocess to avoid libcusolver conflict with torch.
    # Free GPU memory before spawning subprocess.
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    import subprocess as _sp
    cmd = [sys.executable, "-m", "pdb2reaction", "dft"] + list(args)
    _echo(f"\n[dft] subprocess: {' '.join(cmd)}")
    proc = _sp.run(cmd, capture_output=True, text=True)
    if proc.stdout:
        _echo(proc.stdout.rstrip())
    if proc.returncode != 0:
        _echo(f"[dft] WARNING: dft exited with code {proc.returncode}", err=True)
        if proc.stderr:
            _echo(proc.stderr.rstrip(), err=True)
    y = out_dir / "result.yaml"
    if y.exists():
        try:
            data = yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception as exc:
            logger.debug("Failed to parse DFT result YAML %s: %s", y, exc)
            data = {}
    else:
        data = {}
    # Mark DFT success/failure based on convergence
    converged = (data.get("energy") or {}).get("converged", False)
    data["_dft_converged"] = bool(converged)
    data["_dft_failed"] = not bool(converged) or proc.returncode != 0
    return data


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
    Run tsopt CLI on a HEI model structure; return (final_geom_path, ts_geom).

    Prefer the XYZ output to preserve coordinate precision between workflow steps, while still writing
    PDB/GJF companions when requested by the original input type.
    """
    overrides = overrides or {}
    prepared_input = prepare_input_structure(hei_pdb)
    try:
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
        tsopt_mode = None if opt_mode is None else str(opt_mode).strip().lower()
        if tsopt_mode in ("grad", "lbfgs", "dimer"):
            tsopt_mode = "grad"
        elif tsopt_mode in ("hess", "rfo", "rsprfo"):
            tsopt_mode = "hess"

        ts_args: List[str] = [
            "-i",
            str(hei_pdb),
            "-q",
            str(int(charge)),
            "-m",
            str(int(spin)),
            "--out-dir",
            str(ts_dir),
        ]
        _append_toggle_arg(ts_args, "--freeze-links", bool(freeze_use))
        _append_toggle_arg(ts_args, "--convert-files", bool(convert_files))

        if tsopt_mode is not None:
            ts_args.extend(["--opt-mode", tsopt_mode])

        _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
        _append_toggle_arg(ts_args, "--dump", overrides.get("dump"))
        _append_cli_arg(ts_args, "--thresh", overrides.get("thresh"))
        _append_toggle_arg(ts_args, "--flatten", overrides.get("flatten"))

        hess_mode = overrides.get("hessian_calc_mode")
        if hess_mode:
            ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

        # Pass backend from calc_cfg so tsopt uses the same MLIP (or custom calc-file).
        if not _forward_calc_file_argv(ts_args, calc_cfg):
            backend_name = calc_cfg.get("backend", "uma")
            if backend_name != "uma":
                ts_args.extend(["--backend", backend_name])

        # Forward MLIP runtime knobs when CLI-explicit (calc_cfg only contains keys
        # supplied via --workers / --workers-per-node / --solvent / --solvent-model;
        # YAML-only values stay in --config args_yaml below). Without these explicit
        # CLI argv entries the subprocess silently falls back to defaults.
        _workers_cli = calc_cfg.get("workers")
        if _workers_cli is not None and int(_workers_cli) != 1:
            ts_args.extend(["--workers", str(int(_workers_cli))])
        _workers_per_node_cli = calc_cfg.get("workers_per_node")
        if _workers_per_node_cli is not None and int(_workers_per_node_cli) != 1:
            ts_args.extend(["--workers-per-node", str(int(_workers_per_node_cli))])
        _solvent_cli = calc_cfg.get("solvent")
        if _solvent_cli is not None and str(_solvent_cli).lower() != "none":
            ts_args.extend(["--solvent", str(_solvent_cli)])
        _solvent_model_cli = calc_cfg.get("solvent_model")
        if _solvent_model_cli is not None and str(_solvent_model_cli).lower() != "alpb":
            ts_args.extend(["--solvent-model", str(_solvent_model_cli)])

        if args_yaml is not None:
            ts_args.extend(["--config", str(args_yaml)])
        if ref_pdb is not None:
            ts_args.extend(["--ref-pdb", str(ref_pdb)])

        _echo()
        _echo_detail(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
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
                _echo(f"[tsopt] WARNING: Failed to convert TS geometry: {e}", err=True)

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
        calc = create_calculator(**calc_args)
        g_ts.set_calculator(calc)

        return ts_geom_path, g_ts
    finally:
        prepared_input.cleanup()


def _irc_and_match(
    seg_idx: int,
    seg_dir: Path,
    ref_pdb_for_seg: Path,
    seg_model_pdb: Path,
    ref_pdb_template: Optional[Path],
    g_ts: Any,
    q_int: int,
    spin: int,
    freeze_links_flag: bool,
    calc_cfg: Dict[str, Any],
    args_yaml: Optional[Path],
    convert_files: bool,
    seg_tag: Optional[str] = None,
    mep_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Run IRC via the irc CLI (EulerPC), then map the IRC endpoints to (left, right).

    - Run irc on the TS structure into ``seg_dir/irc``.
    - Read endpoints from ``finished_irc_trj.xyz``.
    - When ``seg_tag`` is provided, map the two IRC endpoints to the corresponding
      GSM segment endpoints loaded from
      ``<path_dir>/<seg_tag>_refine_mep/final_geometries_trj.xyz`` (or
      ``<path_dir>/<seg_tag>_mep/final_geometries_trj.xyz`` as a fallback). Bond-state
      matching is attempted first; if that fails, the assignment that minimizes the
      sum of RMSDs to the GSM endpoints is used.
    - For TSOPT-only mode (``seg_tag`` is ``None``), the original (first, last)
      IRC endpoints are kept as (left, right).
    - Returns the endpoint geometries, tags, and paths to the per-segment IRC plot
      and ``finished_irc_trj.xyz``, together with a flag indicating whether the IRC
      trajectory should be reversed when constructing the global IRC plot.
    """
    freeze_atoms: List[int] = _get_freeze_atoms(seg_model_pdb, freeze_links_flag)

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
    ]
    _append_toggle_arg(
        irc_args,
        "--freeze-links",
        bool(
            freeze_links_flag
            and (
                ref_pdb_for_seg.suffix.lower() == ".pdb"
                or (ref_pdb_template is not None and ref_pdb_template.suffix.lower() == ".pdb")
            )
        ),
    )
    _append_toggle_arg(irc_args, "--convert-files", bool(convert_files))
    if ref_pdb_template is not None:
        irc_args.extend(["--ref-pdb", str(ref_pdb_template)])

    # Pass backend from calc_cfg so IRC uses the same MLIP (or custom calc-file).
    if not _forward_calc_file_argv(irc_args, calc_cfg):
        backend_name = calc_cfg.get("backend", "uma")
        if backend_name != "uma":
            irc_args.extend(["--backend", backend_name])

    # Forward MLIP runtime knobs from calc_cfg when CLI-explicit (only set in the
    # shared calc_cfg when user passed them; otherwise downstream irc CLI would
    # fall back to its own defaults regardless of parent --workers / --solvent etc.).
    _workers_cli = calc_cfg.get("workers")
    if _workers_cli is not None and int(_workers_cli) != 1:
        irc_args.extend(["--workers", str(int(_workers_cli))])
    _workers_per_node_cli = calc_cfg.get("workers_per_node")
    if _workers_per_node_cli is not None and int(_workers_per_node_cli) != 1:
        irc_args.extend(["--workers-per-node", str(int(_workers_per_node_cli))])
    _solvent_cli = calc_cfg.get("solvent")
    if _solvent_cli is not None and str(_solvent_cli).lower() != "none":
        irc_args.extend(["--solvent", str(_solvent_cli)])
    _solvent_model_cli = calc_cfg.get("solvent_model")
    if _solvent_model_cli is not None and str(_solvent_model_cli).lower() != "alpb":
        irc_args.extend(["--solvent-model", str(_solvent_model_cli)])

    if args_yaml is not None:
        irc_args.extend(["--config", str(args_yaml)])
    _echo()
    _echo_detail(f"[irc] Running EulerPC IRC → out={irc_dir}")
    _run_cli_main("irc", _irc_cli.cli, irc_args, on_nonzero="raise", prefix="irc")

    finished_pdb = irc_dir / "finished_irc.pdb"
    finished_trj = irc_dir / "finished_irc_trj.xyz"
    irc_plot = irc_dir / "irc_plot.png"

    # Ensure we have a PDB for visualization if possible
    try:
        if finished_trj.exists() and (not finished_pdb.exists()):
            ref_for_conv: Optional[Path] = None
            if seg_model_pdb.suffix.lower() == ".pdb":
                ref_for_conv = seg_model_pdb
            elif ref_pdb_for_seg.suffix.lower() == ".pdb":
                ref_for_conv = ref_pdb_for_seg
            if ref_for_conv is not None:
                _path_search._convert_to_pdb_logged(finished_trj, ref_pdb_path=ref_for_conv, out_path=finished_pdb)
    except Exception as e:
        _echo(f"[irc] WARNING: failed to convert finished_irc_trj.xyz to PDB: {e}", err=True)

    elems, c_first, c_last = read_xyz_first_last(finished_trj)

    calc_args = dict(calc_cfg)
    shared_calc = create_calculator(**calc_args)
    g_left = _geom_from_angstrom(elems, c_first, freeze_atoms)
    g_right = _geom_from_angstrom(elems, c_last, freeze_atoms)
    g_left.set_calculator(shared_calc)
    g_right.set_calculator(shared_calc)

    left_tag = "backward"
    right_tag = "forward"
    reverse_irc = False

    path_root = mep_dir if mep_dir is not None else seg_dir.parent

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
                    except Exception as exc:
                        logger.debug("Bond-change check failed during IRC endpoint matching: %s", exc)
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
                        _echo(
                            f"[irc] WARNING: segment endpoint mapping via RMSD failed: {e}",
                            err=True,
                        )
            else:
                _echo(
                    f"[irc] WARNING: LBFGS endpoints not found for segment tag '{seg_tag}' in the segment directories; "
                    "using raw IRC orientation.",
                    err=True,
                )
        except Exception as e:
            _echo(f"[irc] WARNING: segment endpoint mapping failed: {e}", err=True)
    else:
        # TSOPT-only mode: use raw IRC orientation.
        _echo_detail("[irc] TSOPT-only mode: Use raw irc orientation.")

    # Per-segment IRC plot
    try:
        if finished_trj.exists():
            run_trj2fig(finished_trj, [irc_plot], unit="kcal", reference="init", reverse_x=False)
            close_matplotlib_figures()
    except Exception as e:
        _echo(f"[irc] WARNING: failed to plot finished IRC trajectory: {e}", err=True)

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



_ALL_PRIMARY_HELP_OPTIONS = frozenset(
    {
        "-i",
        "--input",
        "-c",
        "--center",
        "-l",
        "--ligand-charge",
        "-q",
        "--charge",
        "--out-dir",
        "--tsopt",
        "--thermo",
        "--dft",
        "--dft-func-basis",
        "--config",
        "--dry-run",
        "-s",
        "--scan-lists",
        "-b",
        "--backend",
        "--refine-path",
        "-o",
        "--help-advanced",
    }
)


def _show_advanced_help(
    ctx: click.Context, _param: click.Parameter, value: bool
) -> None:
    """Print full option help (including hidden advanced options) and exit."""
    if not value or ctx.resilient_parsing:
        return

    hidden = getattr(ctx.command, "_advanced_hidden_options", ())
    restored: list[click.Option] = []
    for opt in hidden:
        if opt.hidden:
            opt.hidden = False
            restored.append(opt)
    try:
        click.echo(ctx.command.get_help(ctx))
    finally:
        for opt in restored:
            opt.hidden = True
    ctx.exit()


def _configure_all_help_visibility(command: click.Command) -> None:
    """Hide advanced options from default --help while keeping them functional."""
    hidden_options: list[click.Option] = []
    for param in command.params:
        if not isinstance(param, click.Option):
            continue
        names = set(param.opts + param.secondary_opts)
        if names & _ALL_PRIMARY_HELP_OPTIONS:
            continue
        if param.hidden:
            continue
        param.hidden = True
        hidden_options.append(param)
    setattr(command, "_advanced_hidden_options", tuple(hidden_options))


@click.command(
    help=(
        "Run active site model extraction → (optional single-structure staged scan) → MEP search → merge to full PDBs in a single run.\n"
        "If exactly one input is provided: (a) with --scan-lists, run staged scan on the active site model (or full structure "
        "when extraction is skipped) and use stage results as inputs for path-opt (path_search with --refine-path True); "
        "(b) with --tsopt True and no --scan-lists, run TSOPT-only mode."
    ),
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "--help-advanced",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_show_advanced_help,
    help="Show all options (including advanced settings) and exit.",
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
        "the input may be PDB/XYZ/GJF (integer indices only for non-PDB inputs). "
        "Repeat -i/--input for each file."
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
        "or a residue-name list like 'GPP,SAM'. "
        "When omitted, extraction is skipped and the **full input structure(s)** are used directly as active site models."
    ),
)
@click.option(
    "-o", "--out-dir",
    "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path(OUT_DIR_ALL),
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
    "--include-h2o",
    "include_h2o",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Include waters (HOH/WAT/TIP3/SOL) in the active site model.",
)
@click.option(
    "--exclude-backbone",
    "exclude_backbone",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).",
)
@click.option(
    "--add-linkh",
    "add_linkh",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Add cap hydrogens for severed bonds (carbon boundaries only) in active site models.",
)
@click.option(
    "--selected-resn",
    type=str,
    default="",
    show_default=True,
    help=(
        "Force-include residues by residue ID (not name; e.g. '123', 'A:123A', 'B:456'); "
        "comma/space separated. Use '-c/--center GPP,SAM' for residue-name selection."
    ),
)
@click.option(
    "--modified-residue",
    type=str,
    default="",
    show_default=True,
    help=(
        "Comma-separated residue names (with optional charge) to treat as amino acids "
        "for backbone truncation and charge assignment. "
        "Examples: 'HD1,HD2,HD3' or 'HD1:0,SEP:-2'."
    ),
)
@click.option(
    "-l",
    "--ligand-charge",
    type=str,
    default=None,
    help=(
        "Total charge (number) or per-resname mapping like 'GPP:-3,SAM:1'. "
        "The per-resname mapping is applied whether or not extraction runs: with "
        "-c/--center it feeds the extractor charge summary; with -c omitted "
        "(extraction skipped) the same mapping is applied to the full input PDB to "
        "derive the total system charge. A bare number sets the total directly. "
        "PDB inputs only. To force a total charge regardless of residues, use "
        "-q/--charge (emits a warning)."
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
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel MLIP predictor (workers>1).",
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              show_default=True, help="xTB solvent model.")
# ===== Path search knobs =====
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=1,
    show_default=True,
    help="Spin multiplicity (2S+1).",
)
@click.option(
    "--freeze-links",
    "freeze_links_flag",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Freeze parent atoms of cap hydrogens (PDB input or XYZ/GJF with --ref-pdb).",
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
    "--max-nodes",
    type=int,
    default=20,
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
    help="Enable climbing image for standard GSM segments (bridge segments always disable climbing).",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="grad",
    show_default=True,
    help=(
        "Optimizer mode forwarded to scan/tsopt and used for single optimizations: "
        "grad (=LBFGS/Dimer) or hess (=RFO for scan/opt; RS-P-RFO for tsopt)."
    ),
)
@click.option(
    "--opt-mode-post",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="hess",
    show_default=True,
    help=(
        "Optimizer mode override for TSOPT/post-IRC endpoint optimizations. "
        "If unset, uses --opt-mode when explicitly provided; otherwise falls back to the default ('hess' = RS-P-RFO)."
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
    default=False,
    show_default=True,
    help=(
        "Run a single-pass path-opt GSM between each adjacent pair and concatenate the "
        "segments (default; no path_search). Use --refine-path True to run recursive "
        "path_search on the full ordered series for automatic multistep discovery."
    ),
)
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default=False,
    help=(
        "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
    ),
)
@click.option(
    "--thresh-post",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default="baker",
    show_default=True,
    help=(
        "Convergence preset for post-IRC endpoint optimizations "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
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
    help=(
        "Validate options and print the execution plan. With --center, runs extraction "
        "in a temporary directory to validate derived charge and electron parity; no "
        "computational stage or persistent output is produced."
    ),
)
@click.option(
    "--preopt",
    "preopt",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If True, run initial single-structure optimizations of the active site model inputs.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="Common MLIP Hessian calculation mode forwarded to tsopt and freq. Defaults to 'FiniteDifference'.",
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
        "and build Gibbs free-energy diagram (MLIP)."
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
        "With --thermo True, also generate a DFT//MLIP Gibbs diagram."
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
    "--flatten/--no-flatten",
    "flatten",
    default=False,
    show_default=True,
    help="Enable the extra-imaginary-mode flattening loop in tsopt (grad: dimer loop, hess: post-RS-P-RFO); --no-flatten forces flatten_max_iter=0.",
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
    type=click.Choice(["gpu", "cpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="DFT backend: gpu (GPU4PySCF, raises error if unavailable) or cpu (PySCF).",
)
@click.option(
    "-s", "--scan-lists",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help=(
        "Scan targets: inline Python literal or a YAML/JSON spec file path. "
        "Multiple inline literals define sequential stages, e.g. "
        "'[(12,45,1.35)]' '[(10,55,2.20),(23,34,1.80)]'. "
        "Indices refer to the original full input PDB (1-based). When extraction is used, they are "
        "auto-mapped to the active site model after extraction."
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
    help="Override scan harmonic bias strength k (eV/Å^2). Defaults to 300.",
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
    help="Override scan --preopt flag. Inherits from --preopt when omitted.",
)
@click.option(
    "--scan-endopt",
    "scan_endopt_override",
    type=click.BOOL,
    default=None,
    help="Override scan --endopt flag. Defaults to False.",
)
@click.option(
    "--ref-pdb",
    "ref_pdb_cli",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "Reference PDB for topology when -i provides XYZ inputs. "
        "Enables PDB output conversion in TSOPT-only, scan, and path_search pipelines."
    ),
)
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
    center_spec: Optional[str],
    out_dir: Path,
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    modified_residue: str,
    ligand_charge: Optional[str],
    charge_override: Optional[int],
    workers: int,
    workers_per_node: int,
    backend: str,
    solvent: str,
    solvent_model: str,
    spin: int,
    freeze_links_flag: bool,
    mep_mode: str,
    dmf_backend: str,
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
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
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
    ref_pdb_cli: Optional[Path],
    tsopt_max_cycles: Optional[int],
    tsopt_out_dir: Optional[Path],
    flatten: bool,
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
    cli_coord_type: Optional[str],
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: str,
) -> None:
    """
    The **all** command composes `extract` → (optional `scan` on model or full input) → MEP search
    (single-pass `path-opt` by default, or recursive `path_search` with ``--refine-path True``) and hides
    ref-template bookkeeping.
    With single input:
      - with --scan-lists: run staged scan and use stage results as inputs for path-opt (or path_search),
      - with --tsopt True and no --scan-lists: run TSOPT-only mode (no MEP search).

    By default, a single-pass ``path-opt`` GSM is run between each adjacent pair of inputs and the segments
    are concatenated into the final MEP.  With ``--refine-path True``, the recursive ``path_search`` workflow
    is used instead for automatic multistep discovery.
    """
    # Engage pipeline-scoped default-verbosity suppression for the duration of
    # this `all` run (reset per-invocation by DefaultGroup.parse_args). Standalone
    # leaf/report commands are unaffected and keep full output at default.
    from pdb2reaction.core.utils import set_pipeline_mode
    set_pipeline_mode(True)
    argv_all = sys.argv[1:]

    def _sigint_handler(signum, frame):
        _echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    signal.signal(signal.SIGINT, _sigint_handler)

    _echo_state.reset()
    global _FREEZE_ATOMS_GLOBAL, _FREEZE_ATOMS_YAML
    _FREEZE_ATOMS_GLOBAL = None
    _FREEZE_ATOMS_YAML = None
    set_convert_file_enabled(convert_files)
    command_str = " ".join(sys.argv)
    time_start = time.perf_counter()
    energy_diagrams: List[Dict[str, Any]] = []

    dump_override_requested = False
    try:
        dump_source = ctx.get_parameter_source("dump")
        dump_override_requested = dump_source not in (None, ParameterSource.DEFAULT)
    except Exception as exc:
        logger.debug("Failed to check dump parameter source: %s", exc)
        dump_override_requested = False

    args_yaml = config_yaml
    merged_yaml_cfg = load_yaml_dict(config_yaml) if config_yaml is not None else {}

    opt_mode_set = False
    opt_mode_post_set = False
    try:
        opt_mode_source = ctx.get_parameter_source("opt_mode")
        opt_mode_set = opt_mode_source not in (None, ParameterSource.DEFAULT)
        opt_mode_post_source = ctx.get_parameter_source("opt_mode_post")
        opt_mode_post_set = opt_mode_post_source not in (None, ParameterSource.DEFAULT)
    except Exception as exc:
        logger.debug("Failed to check opt_mode parameter sources: %s", exc)
        opt_mode_set = False
        opt_mode_post_set = False

    _mode_alias = {
        "grad": "grad",
        "hess": "hess",
        "lbfgs": "grad",
        "rfo": "hess",
        "dimer": "grad",
        "rsprfo": "hess",
    }
    opt_mode_norm = _mode_alias.get(str(opt_mode).strip().lower(), "hess")
    opt_mode_post_norm = (
        None
        if opt_mode_post is None
        else _mode_alias.get(str(opt_mode_post).strip().lower(), "hess")
    )

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

    scan_vals = collect_single_option_values(argv_all, ("-s", "--scan-lists"), "--scan-lists")
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

    if single_tsopt_mode:
        all_mode = "tsopt-only"
    elif has_scan:
        all_mode = "scan-to-path-search" if refine_path else "scan-to-path-opt"
    else:
        all_mode = "path-search" if refine_path else "path-opt"
    all_mode_label = "ts-only" if single_tsopt_mode else ("scan-lists" if has_scan else "mep")
    if verbose_level() >= 2:
        _echo(
            f"[mode] all ({all_mode_label}) inputs={len(input_paths)} "
            f"extract={'yes' if center_spec is not None and str(center_spec).strip() else 'no'} "
            f"scan={'yes' if has_scan else 'no'} tsopt={'yes' if do_tsopt else 'no'} "
            f"thermo={'yes' if do_thermo else 'no'} dft={'yes' if do_dft else 'no'} "
            f"dry_run={'yes' if dry_run else 'no'} internal={all_mode}",
            narrative=True,
        )

    tsopt_opt_mode_default: Optional[str] = None
    if opt_mode_post_set and opt_mode_post_norm is not None:
        tsopt_opt_mode_default = opt_mode_post_norm
    elif opt_mode_set:
        tsopt_opt_mode_default = opt_mode_norm
    else:
        tsopt_opt_mode_default = "hess"
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
    tsopt_overrides["flatten"] = bool(flatten)

    freq_overrides: Dict[str, Any] = {}
    # backend will be injected after calc_cfg_shared is built (see below)
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
    # all.py reads thermochemistry EXCLUSIVELY from freq's thermoanalysis.yaml
    # (its in-memory return is discarded via _run_cli_main). Always let freq
    # write that yaml so `--dump False` cannot silently zero out tR/tT/tP and
    # relabel electronic energies as Gibbs.
    freq_overrides["dump"] = True
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

    if show_config or (dry_run and verbose_level() >= 3):
        config_payload: Dict[str, Any] = {
            "yaml": {
                "config": str(config_yaml) if config_yaml else None,
                "effective_args_yaml": str(args_yaml) if args_yaml else None,
            },
            "all": {
                "inputs": [str(p) for p in input_paths],
                "center": (None if center_spec is None else str(center_spec)),
                "out_dir": str(out_dir),
                "charge_override": charge_override,
                "spin": int(spin),
                "mep_mode": str(mep_mode),
                "max_nodes": int(max_nodes),
                "max_cycles": int(max_cycles),
                "climb": bool(climb),
                "opt_mode": str(opt_mode),
                "opt_mode_post": (None if opt_mode_post is None else str(opt_mode_post)),
                "dump": bool(dump),
                "convert_files": bool(convert_files),
                "refine_path": bool(refine_path),
                "preopt": bool(preopt),
                "tsopt": bool(do_tsopt),
                "thermo": bool(do_thermo),
                "dft": bool(do_dft),
                "dft_engine": str(dft_engine),
            },
            "overrides": {
                "tsopt": tsopt_overrides,
                "freq": freq_overrides,
                "dft": dft_overrides,
            },
        }
        if merged_yaml_cfg:
            config_payload["effective_yaml"] = merged_yaml_cfg
        _echo_section("====== [all] Effective configuration ======")
        # `--show-config` is an explicit output request; dry-run's automatic
        # config dump is level-3 debug context so -v 1/2 stay compact.
        click.echo(
            yaml.safe_dump(config_payload, sort_keys=False, allow_unicode=True).rstrip(),
            narrative=show_config,
        )

    if dry_run:
        # Domain validation — run extract pre-pass so silent charge/parity
        # drift (extract says -2, run.sh hardcodes -1 → only fails 30 min
        # into the actual job) is caught here. We invoke extract_api against
        # a throwaway out_dir, then compare the derived total_charge with
        # -q (charge_override) and run the electron parity check
        # (validate_charge_spin_at_path) on the extracted geometry. No
        # tsopt/path_search/freq is started.
        _skip_extract_dry = center_spec is None or str(center_spec).strip() == ""
        if not _skip_extract_dry:
            _dry_tmp = Path(tempfile.mkdtemp(prefix="pdb2reaction_dry_extract_"))
            _first_in = input_paths[0].resolve() if input_paths else None
            if _first_in is None:
                raise click.BadParameter("[all] --dry-run requires at least one -i input.")
            _dry_out = (_dry_tmp / f"model_{_first_in.stem}.pdb").resolve()
            try:
                _ex = extract_api(
                    complex_pdb=[str(_first_in)],
                    center=str(center_spec),
                    output=[str(_dry_out)],
                    radius=float(radius),
                    radius_het2het=float(radius_het2het),
                    include_h2o=bool(include_h2o),
                    exclude_backbone=bool(exclude_backbone),
                    add_linkh=bool(add_linkh),
                    selected_resn=selected_resn or "",
                    modified_residue=modified_residue or "",
                    ligand_charge=ligand_charge,
                    verbose=False,
                )
            except click.ClickException:
                raise
            except Exception as e:
                raise click.ClickException(f"[all] --dry-run extract pre-check failed: {e}")
            _cs = _ex.get("charge_summary", {}) if isinstance(_ex, dict) else {}
            _q_total_raw = _cs.get("total_charge", None)
            if _q_total_raw is None:
                _echo("[all] --dry-run: extractor produced no total_charge; skipping parity check.", err=True)
            else:
                _q_extracted = _round_charge_with_note(float(_q_total_raw), prefix="[all] dry-run")
                _echo(f"[all] --dry-run extract: model total_charge = {_q_extracted:+d}", narrative=True)
                if charge_override is not None and int(charge_override) != int(_q_extracted):
                    raise click.BadParameter(
                        f"[all] -q {int(charge_override)} does not match extract-derived "
                        f"charge {int(_q_extracted):+d}. Either omit -q to auto-derive, "
                        f"pass --ligand-charge for non-standard residues, or set -q to {int(_q_extracted):+d}."
                    )
                _q_check = int(charge_override) if charge_override is not None else int(_q_extracted)
                # Use the module-level import (line 68). A nested
                # `from pdb2reaction.core.utils import validate_charge_spin_at_path` inside
                # this dry_run block would mark the name as a local in
                # cli()'s scope, breaking the non-dry-run reference at the
                # end of cli() with UnboundLocalError when --dry-run is not
                # passed.
                try:
                    validate_charge_spin_at_path(_dry_out, _q_check, int(spin))
                except click.ClickException:
                    raise
                except ValueError as e:
                    raise click.BadParameter(f"[all] --dry-run parity check failed: {e}")
                _echo(f"[all] --dry-run parity check OK: charge={_q_check:+d}, spin(multiplicity)={int(spin)}", narrative=True)
            shutil.rmtree(_dry_tmp, ignore_errors=True)
        _echo("[all] Dry-run mode: no search/post-processing was executed (extract pre-check ran).", narrative=True)
        _echo(
            "[all] Planned stages: extract -> optional scan -> path_opt/path_search -> optional tsopt/irc/freq/dft.",
            narrative=True,
        )
        _emit_final_summary(out_dir, time_start)
        return

    yaml_cfg = load_yaml_dict(args_yaml)
    _set_yaml_freeze_atoms(yaml_cfg)

    skip_extract = center_spec is None or str(center_spec).strip() == ""
    first_input = input_paths[0].resolve() if input_paths else None

    out_dir = out_dir.resolve()
    models_dir = out_dir / WORK_DIRNAME / "models"
    path_dir = out_dir / WORK_DIRNAME / ("path_search" if refine_path else "path_opt")
    scan_dir = _resolve_override_dir(out_dir / WORK_DIRNAME / "scan", scan_out_dir)
    stage_total = 4 if (do_tsopt or do_thermo or do_dft) else 3
    ensure_dir(out_dir)
    if not skip_extract:
        ensure_dir(models_dir)
    if not single_tsopt_mode:
        ensure_dir(path_dir)

    elem_tmp_dir = out_dir / WORK_DIRNAME / "add_elem_info"
    inputs_for_extract: List[Path] = []
    elem_fix_echo = False
    for p in input_paths:
        if _pdb_needs_elem_fix(p):
            if not elem_fix_echo:
                _echo_section(
                    "====== [all] Preflight — add_elem_info (only when element fields are missing) ======"
                )
                elem_fix_echo = True
            ensure_dir(elem_tmp_dir)
            out_p = (elem_tmp_dir / p.name).resolve()
            try:
                _assign_elem_info(str(p), str(out_p), overwrite=False)
                _echo(f"[all] add_elem_info: fixed elements → {out_p}")
                inputs_for_extract.append(out_p)
            except SystemExit as e:
                code = getattr(e, "code", 1)
                _echo(
                    f"[all] WARNING: add_elem_info exited with code {code} for {p}; using original.",
                    err=True,
                )
                inputs_for_extract.append(p.resolve())
            except Exception as e:
                _echo(
                    f"[all] WARNING: add_elem_info failed for {p}: {e} — using original file.",
                    err=True,
                )
                inputs_for_extract.append(p.resolve())
        else:
            inputs_for_extract.append(p.resolve())

    # --- fix_altloc: drop alternate locations (only when altLoc is detected) ---
    altloc_tmp_dir = out_dir / WORK_DIRNAME / "fix_altloc"
    final_inputs: List[Path] = []
    altloc_fix_echo = False
    for p in inputs_for_extract:
        if _has_altloc(p):
            if not altloc_fix_echo:
                _echo_section(
                    "====== [all] Preflight — fix_altloc (only when altLoc is detected) ======"
                )
                altloc_fix_echo = True
            ensure_dir(altloc_tmp_dir)
            out_p = (altloc_tmp_dir / p.name).resolve()
            try:
                _fix_altloc(str(p), str(out_p), overwrite=True, skip_if_no_altloc=False)
                _echo(f"[all] fix_altloc: fixed altLoc → {out_p}")
                final_inputs.append(out_p)
            except Exception as e:
                _echo(
                    f"[all] WARNING: fix_altloc failed for {p}: {e} — using original file.",
                    err=True,
                )
                final_inputs.append(p)
        else:
            final_inputs.append(p)

    extract_inputs = tuple(final_inputs)

    model_outputs: List[Path] = []
    if not skip_extract:
        for p in extract_inputs:
            model_outputs.append((models_dir / f"model_{p.stem}.pdb").resolve())

    resolved_charge: Optional[int] = None

    if not skip_extract:
        _echo_section(
            f"====== [all] Stage 1/{stage_total} — Active site model extraction ======"
        )
        try:
            ex_res = extract_api(
                complex_pdb=[str(p) for p in extract_inputs],
                center=center_spec,
                output=[str(p) for p in model_outputs],
                radius=float(radius),
                radius_het2het=float(radius_het2het),
                include_h2o=bool(include_h2o),
                exclude_backbone=bool(exclude_backbone),
                add_linkh=bool(add_linkh),
                selected_resn=selected_resn or "",
                modified_residue=modified_residue or "",
                ligand_charge=ligand_charge,
                verbose=True,  # extractor INFO now gated by the unified global -v level
            )
        except Exception as e:
            raise click.ClickException(f"[all] Extractor failed: {e}")

        _echo("[all] Active site model files:")
        for op in model_outputs:
            _echo(f"  - {op}")

        try:
            cs = ex_res.get("charge_summary", {})
            q_total = float(cs.get("total_charge", 0.0))
            q_prot = float(cs.get("protein_charge", 0.0))
            q_lig = float(cs.get("ligand_total_charge", 0.0))
            q_ion = float(cs.get("ion_total_charge", 0.0))
            _echo("[all] Charge summary from extractor (model #1):")
            _echo(
                f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}"
            )
            resolved_charge = _round_charge_with_note(q_total, prefix="[all]")
        except Exception as e:
            raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")
    else:
        _echo_section(
            f"====== [all] Stage 1/{stage_total} — Extraction skipped (no -c/--center); using FULL structures as active site models ======"
        )
        first_input = input_paths[0].resolve()
        gjf_charge: Optional[int] = None
        gjf_spin: Optional[int] = None

        user_provided_spin = True
        try:
            spin_source = ctx.get_parameter_source("spin")
            user_provided_spin = spin_source not in (None, ParameterSource.DEFAULT)
        except Exception as exc:
            logger.debug("Failed to check spin parameter source: %s", exc)
            user_provided_spin = True

        ligand_charge_numeric: Optional[float] = None
        if ligand_charge is not None:
            try:
                ligand_charge_numeric = float(ligand_charge)
            except Exception as exc:
                logger.debug("Failed to parse ligand_charge as float: %s", exc)
                ligand_charge_numeric = None

            if first_input.suffix.lower() == ".pdb":
                try:
                    with prepare_input_structure(first_input) as prepared:
                        resolved_charge = _derive_charge_from_ligand_charge(
                            prepared, ligand_charge, prefix="[all]"
                        )
                except Exception as e:
                    _echo(
                        f"[all] NOTE: failed to derive total charge from full complex: {e}; "
                        "falling back to legacy handling.",
                        err=True,
                    )
            else:
                _echo(
                    "[all] NOTE: --ligand-charge derivation requires a PDB input; skipping full-complex derivation.",
                    err=True,
                )

        if resolved_charge is None:
            if ligand_charge_numeric is not None:
                _echo(
                    f"[all] Using --ligand-charge as TOTAL system charge: {ligand_charge_numeric:+g}"
                )
                resolved_charge = _round_charge_with_note(ligand_charge_numeric, prefix="[all]")

        if first_input.suffix.lower() == ".gjf" and (resolved_charge is None or not user_provided_spin):
            try:
                with prepare_input_structure(first_input) as prepared:
                    gjf_charge, gjf_spin = resolve_charge_spin(
                        prepared, charge=None, spin=None
                    )
                _echo(
                    f"[all] Parsed GJF metadata (first input): charge={gjf_charge:+d}, spin={gjf_spin}"
                )
            except Exception as e:
                _echo(
                    f"[all] NOTE: failed to parse charge/spin from GJF '{first_input.name}': {e}",
                    err=True,
                )

        if resolved_charge is None and gjf_charge is not None:
            _echo(f"[all] Using total charge from first GJF: {float(gjf_charge):+g}")
            resolved_charge = _round_charge_with_note(float(gjf_charge), prefix="[all]")

        if resolved_charge is None and charge_override is None:
            raise click.ClickException(
                "[all] Total charge could not be resolved. Provide -q/--charge, "
                "--ligand-charge, or a .gjf input with charge metadata."
            )

        if (not user_provided_spin) and (gjf_spin is not None):
            spin = int(gjf_spin)
            _echo(f"[all] Spin multiplicity set from GJF: {spin}")

    if charge_override is not None:
        q_int = int(charge_override)
        if (not skip_extract) and resolved_charge is not None and int(resolved_charge) != q_int:
            raise click.BadParameter(
                f"[all] -q {q_int:+d} does not match extract-derived charge "
                f"{int(resolved_charge):+d}. Either omit -q to auto-derive, "
                f"pass --ligand-charge for non-standard residues, or set -q to "
                f"{int(resolved_charge):+d}."
            )
        override_msg = (
            f"[all] -q/--charge override supplied; using TOTAL system charge {q_int:+d}"
        )
        if resolved_charge is not None:
            override_msg += f" (matches workflow-derived {int(resolved_charge):+d})"
        _echo(override_msg, err=True)
    else:
        q_int = int(resolved_charge) if resolved_charge is not None else 0

    _validate_path = model_outputs[0] if model_outputs else input_paths[0]
    validate_charge_spin_at_path(_validate_path, q_int, spin)

    # Resolve --ref-pdb for topology
    ref_pdb_for_topology: Optional[Path] = None
    if ref_pdb_cli is not None:
        ref_pdb_for_topology = ref_pdb_cli.resolve()
        _echo(f"[all] --ref-pdb provided: {ref_pdb_for_topology}")

    freeze_ref: Optional[Path] = None
    if freeze_links_flag:
        if not skip_extract and model_outputs:
            freeze_ref = model_outputs[0]
        elif ref_pdb_for_topology is not None:
            freeze_ref = ref_pdb_for_topology
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
        coord_type=(
            str(cli_coord_type).lower()
            if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None
            else None
        ),
        precision=precision,
        backend_model=backend_model,
    )

    calc_cfg_shared = _build_calc_cfg(
        q_int,
        spin,
        workers=workers if cli_param_overridden(ctx, "workers") else None,
        workers_per_node=workers_per_node if cli_param_overridden(ctx, "workers_per_node") else None,
        yaml_cfg=yaml_cfg,
        backend=backend if cli_param_overridden(ctx, "backend") else None,
        solvent=solvent if cli_param_overridden(ctx, "solvent") else None,
        solvent_model=solvent_model if cli_param_overridden(ctx, "solvent_model") else None,
    )
    # calc_cfg_shared feeds the run summary (mlip_backend / mlip_model). _build_calc_cfg does
    # not apply --backend-model, so without this the summary records the DEFAULT UMA model even
    # when --backend-model overrides it (the actual calc already honors --backend-model via its
    # own config path; this only corrects the recorded provenance). Use the same canonical helper
    # as path_search for consistency.
    from pdb2reaction.backends import apply_backend_model_to_calc_cfg
    # Unconditional: also pops a raw backend_model token from a --config YAML
    # (the helper no-ops when neither the CLI arg nor the YAML names one).
    apply_backend_model_to_calc_cfg(calc_cfg_shared, backend_model)
    # Same class of fix as --backend-model above: _build_calc_cfg does not apply --precision, so
    # without this the in-process calculators (create_calculator(**calc_cfg_shared) at the TS
    # re-eval / pre-align / endpoint R,P,TS energy sites) run at the default fp32 even when
    # --precision fp64 is requested — mixing an fp32 electronic energy with the fp64 subprocess
    # geometry and thermal correction. Apply via the backend-aware helper (orb/mace key names +
    # aimnet2 fp64 rejection), before the calc-file override may switch to a custom backend.
    from pdb2reaction.backends import apply_effective_precision
    apply_effective_precision(calc_cfg_shared, precision)

    # --calc-file overrides --backend with a user ASE Calculator (custom backend).
    from pdb2reaction.backends import apply_calc_file_to_calc_cfg
    apply_calc_file_to_calc_cfg(calc_cfg_shared, calc_file, calc_factory)

    # Inject backend into freq_overrides so _run_freq_for_state passes it
    _backend_shared = calc_cfg_shared.get("backend", "uma")
    if _backend_shared != "uma":
        freq_overrides["backend"] = _backend_shared
    if calc_cfg_shared.get("calc_file"):
        freq_overrides["calc_file"] = calc_cfg_shared["calc_file"]
        if calc_cfg_shared.get("calc_factory"):
            freq_overrides["calc_factory"] = calc_cfg_shared["calc_factory"]
    # Inject MLIP runtime knobs (workers / workers_per_node / solvent / solvent_model)
    # when CLI-explicit so _run_freq_for_state forwards them to the freq subprocess.
    # Without this, freq silently runs with workers=1 / no solvent even when the parent
    # CLI requested otherwise.
    if cli_param_overridden(ctx, "workers"):
        freq_overrides["workers"] = int(workers)
    if cli_param_overridden(ctx, "workers_per_node"):
        freq_overrides["workers_per_node"] = int(workers_per_node)
    if cli_param_overridden(ctx, "solvent"):
        freq_overrides["solvent"] = str(solvent)
    if cli_param_overridden(ctx, "solvent_model"):
        freq_overrides["solvent_model"] = str(solvent_model)

    if single_tsopt_mode:
        _echo_section("====== [all] TSOPT-only single-structure mode ======")
        tsroot = out_dir / SEGMENTS_DIRNAME / "seg_01"
        ensure_dir(tsroot)

        # In TSOPT-only mode, no MEP search is performed. Use a placeholder for
        # MEP-related fields in downstream summaries.
        mep_mode_kind = "---"

        ts_initial_pdb = model_outputs[0] if not skip_extract else input_paths[0].resolve()

        struct_dir = tsroot / "structures"
        freq_root = _resolve_override_dir(tsroot / "freq", freq_out_dir)
        dft_root = _resolve_override_dir(tsroot / "dft", dft_out_dir)

        _tsopt_ref = ref_pdb_for_topology or (ts_initial_pdb if ts_initial_pdb.suffix.lower() == ".pdb" else None)
        ts_pdb, g_ts = _run_tsopt_on_hei(
            ts_initial_pdb,
            q_int,
            spin,
            calc_cfg_shared,
            args_yaml,
            tsroot,
            freeze_links_flag,
            tsopt_opt_mode_default,
            _tsopt_ref,
            convert_files,
            overrides=tsopt_overrides,
        )

        irc_res = _irc_and_match(
            seg_idx=1,
            seg_dir=tsroot,
            ref_pdb_for_seg=ts_pdb,
            seg_model_pdb=ref_pdb_for_topology or ts_initial_pdb,
            ref_pdb_template=_tsopt_ref,
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

        ensure_dir(struct_dir)
        model_ref = ref_pdb_for_topology or ts_initial_pdb
        _save_single_geom_as_pdb_for_tools(
            g_react_irc, model_ref, struct_dir, "reactant_irc"
        )
        _save_single_geom_as_pdb_for_tools(
            g_prod_irc, model_ref, struct_dir, "product_irc"
        )
        pT = _save_single_geom_as_pdb_for_tools(gT, model_ref, struct_dir, "ts")

        endpoint_opt_dir = tsroot / "endpoint_opt"
        ensure_dir(endpoint_opt_dir)

        # Map IRC left/right Hessians → R/P endpoint (left=forward, right=backward)
        from pdb2reaction.io.hessian_cache import load as _hess_load, store as _hess_store
        _react_hk = "irc_left" if eL >= eR_raw else "irc_right"
        _prod_hk  = "irc_right" if eL >= eR_raw else "irc_left"

        _c = _hess_load(_react_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
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
            _echo(
                f"[post] WARNING: Reactant endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            g_react_opt = g_react_irc

        _c = _hess_load(_prod_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
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
            _echo(
                f"[post] WARNING: Product endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            g_prod_opt = g_prod_irc

        # Clean up endpoint_opt as a temporary working directory
        shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
        _echo_detail("[endpoint-opt] Clean endpoint-opt working dir.")

        pR = _save_single_geom_as_pdb_for_tools(
            g_react_opt, model_ref, struct_dir, "reactant"
        )
        pP = _save_single_geom_as_pdb_for_tools(
            g_prod_opt, model_ref, struct_dir, "product"
        )

        e_react = float(g_react_opt.energy)
        e_prod = float(g_prod_opt.energy)

        diag_payload = _write_segment_energy_diagram(
            tsroot / "energy_diagram_MLIP",
            labels=["R", "TS", "P"],
            energies_au=[e_react, eT, e_prod],
            title_note="(MLIP, TSOPT + IRC)",
        )
        if diag_payload:
            energy_diagrams.append(diag_payload)

        # DO NOT INLINE: (segment path, freq-pre): same closure-capture mechanism: null calculator on retained Geometry before freq subprocess so VRAM reports do not double-count.
        # ── Release GPU memory before freq/thermo/DFT ──
        for _g in [locals().get(n) for n in ("gL", "gR", "gT", "g_react_irc", "g_prod_irc", "g_react_opt", "g_prod_opt")]:
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        thermo_payloads: Dict[str, Dict[str, Any]] = {}

        ref_pdb_for_tsopt_only = ref_pdb_for_topology or (
            ts_initial_pdb if ts_initial_pdb.suffix.lower() == ".pdb" else None
        )

        if do_thermo:
            _echo()
            _echo_detail("[thermo] Single TSOPT: freq on TS/R/P")
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
            from pdb2reaction.io.hessian_cache import clear as _clear_hess_cache
            _clear_hess_cache()  # TS Hessian consumed; R/P need exact computation
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
                    tR.get("sum_EE_and_thermal_free_energy_ha", np.nan)
                )
                GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", np.nan))
                GP = float(
                    tP.get("sum_EE_and_thermal_free_energy_ha", np.nan)
                )
                if not all(np.isfinite([GR, GT, GP])):
                    _echo(
                        "[thermo] NOTE: thermochemistry unavailable (freq --dump "
                        "off, 'thermoanalysis' missing, or freq failed); Gibbs "
                        "diagram skipped.",
                        err=True,
                    )
                elif True:
                    diag_payload = _write_segment_energy_diagram(
                        tsroot / "energy_diagram_G_MLIP",
                        labels=["R", "TS", "P"],
                        energies_au=[GR, GT, GP],
                        title_note="(MLIP + Thermal Correction)",
                        ylabel="ΔG (kcal/mol)",
                    )
                    if diag_payload:
                        energy_diagrams.append(diag_payload)
            except Exception as e:
                _echo(
                    f"[thermo] WARNING: failed to build Gibbs diagram: {e}",
                    err=True,
                )

        if do_dft:
            # DO NOT INLINE: (single-TS path): Geometry retains calculator → torch.nn.Module pinned. Two layers (null calculator + del local name) prevent closure/hook capture from resurrecting model before subprocess fork.
            # ── Aggressively release GPU memory before DFT subprocess ──
            # `del locals()[name]` is a CPython no-op: locals() returns a
            # *copy* of the frame namespace, so deletion never propagated
            # to the real binding and the torch.nn.Module references stayed
            # pinned. Rebind each name to None explicitly instead, so the
            # subsequent gc.collect() + torch.cuda.empty_cache() can actually
            # reclaim the GPU pages before the DFT subprocess fork.
            g_ts = None
            calc = None
            # DO NOT null g_react_opt / g_prod_opt / thermo_payloads / tR / tT /
            # tP here: their calculators were already released above (~L3326),
            # and these objects are consumed downstream — bond-change detection
            # (~L3501), DFT//MLIP Gibbs (~L3466), and the Gibbs log (~L3601).
            # Nulling them raised a swallowed AttributeError that silently
            # dropped the DFT//MLIP Gibbs diagram and forced a bogus
            # "(no covalent changes detected)" summary. They hold no GPU tensors.
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            _echo()
            _echo_detail("[dft] Single TSOPT: DFT on R/TS/P")
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
            dR = dft_payloads.get("R") or {}
            dT = dft_payloads.get("TS") or {}
            dP = dft_payloads.get("P") or {}
            eR_dft = _dft_energy_ha(dR)
            eT_dft = _dft_energy_ha(dT)
            eP_dft = _dft_energy_ha(dP)
            _dft_all_ok = all(e is not None for e in (eR_dft, eT_dft, eP_dft))
            if not _dft_all_ok:
                _failed_states = [s for s, e in zip(["R", "TS", "P"], [eR_dft, eT_dft, eP_dft]) if e is None]
                _echo(f"[dft] WARNING: DFT failed for state(s): {', '.join(_failed_states)}. Skipping DFT diagrams.", err=True)
            if _dft_all_ok:
                try:
                    diag_payload = _write_segment_energy_diagram(
                        tsroot / "energy_diagram_DFT",
                        labels=["R", "TS", "P"],
                        energies_au=[eR_dft, eT_dft, eP_dft],
                        title_note=f"({dft_func_basis_use} // MLIP)",
                    )
                    if diag_payload:
                        energy_diagrams.append(diag_payload)
                except Exception as e:
                    _echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            if do_thermo and _dft_all_ok:
                try:
                    dG_R = float(
                        (thermo_payloads.get("R", {}) or {}).get(
                            "thermal_correction_free_energy_ha", np.nan
                        )
                    )
                    dG_T = float(
                        (thermo_payloads.get("TS", {}) or {}).get(
                            "thermal_correction_free_energy_ha", np.nan
                        )
                    )
                    dG_P = float(
                        (thermo_payloads.get("P", {}) or {}).get(
                            "thermal_correction_free_energy_ha", np.nan
                        )
                    )
                    GR_dftUMA = eR_dft + dG_R
                    GT_dftUMA = eT_dft + dG_T
                    GP_dftUMA = eP_dft + dG_P
                    if not all(np.isfinite([GR_dftUMA, GT_dftUMA, GP_dftUMA])):
                        _echo(
                            "[dft//mlip] NOTE: thermochemistry unavailable; "
                            "DFT//MLIP Gibbs diagram skipped.",
                            err=True,
                        )
                    else:
                        diag_payload = _write_segment_energy_diagram(
                            tsroot / "energy_diagram_G_DFT_plus_MLIP",
                            labels=["R", "TS", "P"],
                            energies_au=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                            title_note=f"({dft_func_basis_use} // MLIP + Thermal Correction)",
                            ylabel="ΔG (kcal/mol)",
                        )
                        if diag_payload:
                            energy_diagrams.append(diag_payload)
                except Exception as e:
                    _echo(
                        f"[dft//mlip] WARNING: failed to build DFT//MLIP Gibbs diagram: {e}",
                        err=True,
                    )

        # Build summary.json for TSOPT-only runs, mirroring path_search/path_opt outputs
        bond_cfg = dict(_path_search.BOND_KW)
        bond_summary = ""
        if g_react_opt is not None and g_prod_opt is not None:
            try:
                changed, bond_summary = _path_search.has_bond_change(
                    g_react_opt, g_prod_opt, bond_cfg
                )
                if not changed:
                    bond_summary = "(no covalent changes detected)"
            except Exception as e:
                _echo(
                    f"[post] WARNING: Failed to detect bond changes for TSOPT-only endpoints: {e}",
                    err=True,
                )
                bond_summary = "(no covalent changes detected)"
        else:
            bond_summary = "(no covalent changes detected)"

        barrier = (eT - e_react) * AU2KCALPERMOL
        delta = (e_prod - e_react) * AU2KCALPERMOL

        n_images = 0
        try:
            irc_trj_path = irc_res.get("irc_trj_path")
            if isinstance(irc_trj_path, Path) and irc_trj_path.exists():
                n_images = len(read_xyz_as_blocks(irc_trj_path))
        except Exception as exc:
            logger.debug("Failed to count IRC trajectory images: %s", exc)
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
        _enrich_summary(
            summary,
            version="",
            pipeline_mode="tsopt-only",
            out_dir=out_dir,
            mlip_backend=calc_cfg_shared.get("backend", "uma"),
            mlip_model=calc_cfg_shared.get("model"),
            charge=q_int,
            spin=spin,
            command=command_str,
            config={
                "refine_path": refine_path,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "dft_status": "failed" if (do_dft and not _dft_all_ok) else ("converged" if (do_dft and _dft_all_ok) else None),
                "dft_func_basis": dft_func_basis_use if do_dft else None,
                "opt_mode": tsopt_opt_mode_default,
                "mep_mode": mep_mode_kind,
            },
            freeze_atoms=_freeze_atoms_for_log(),
        )
        try:
            with open(tsroot / "summary.json", "w") as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            _echo_detail(f"[write] Wrote '{tsroot / 'summary.json'}'.")
            _copy_logged(tsroot / "summary.json", out_dir / "summary.json", label="summary.json")
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
                from pdb2reaction.workflows._all_helpers import build_energy_level_dict
                _structs_seg = {"R": pR, "TS": pT, "P": pP}
                segment_log["uma"] = build_energy_level_dict(
                    labels=["R", "TS", "P"],
                    energies_au=[e_react, eT, e_prod],
                    ref_energy=e_react,
                    au_to_kcal=AU2KCALPERMOL,
                    diagram_path=str(tsroot / "energy_diagram_MLIP.png"),
                    structures=_structs_seg,
                )
                _gkey = "sum_EE_and_thermal_free_energy_ha"
                if do_thermo and thermo_payloads and all(
                    _gkey in (thermo_payloads.get(_s) or {}) for _s in ("R", "TS", "P")
                ):
                    GR = float(thermo_payloads["R"][_gkey])
                    GT = float(thermo_payloads["TS"][_gkey])
                    GP = float(thermo_payloads["P"][_gkey])
                    segment_log["gibbs_uma"] = build_energy_level_dict(
                        labels=["R", "TS", "P"],
                        energies_au=[GR, GT, GP],
                        ref_energy=GR,
                        au_to_kcal=AU2KCALPERMOL,
                        diagram_path=str(tsroot / "energy_diagram_G_MLIP.png"),
                        structures=_structs_seg,
                    )
                if do_dft:
                    if _dft_all_ok:
                        segment_log["dft"] = build_energy_level_dict(
                            labels=["R", "TS", "P"],
                            energies_au=[eR_dft, eT_dft, eP_dft],
                            ref_energy=eR_dft,
                            au_to_kcal=AU2KCALPERMOL,
                            diagram_path=str(tsroot / "energy_diagram_DFT.png"),
                            structures=_structs_seg,
                        )
                    else:
                        segment_log["dft"] = {
                            "status": "failed",
                            "failed_states": [s for s, e in zip(["R", "TS", "P"], [eR_dft, eT_dft, eP_dft]) if e is None],
                        }
                    if do_thermo and _dft_all_ok:
                        segment_log["gibbs_dft_uma"] = build_energy_level_dict(
                            labels=["R", "TS", "P"],
                            energies_au=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                            ref_energy=GR_dftUMA,
                            au_to_kcal=AU2KCALPERMOL,
                            diagram_path=str(tsroot / "energy_diagram_G_DFT_plus_MLIP.png"),
                            structures=_structs_seg,
                        )

                _mlip_backend = calc_cfg_shared.get("backend", "uma")
                _mlip_model = calc_cfg_shared.get("model")
                summary_payload = {
                    "root_out_dir": str(out_dir),
                    "path_dir": str(tsroot),
                    "path_module_dir": tsroot.name,
                    "pipeline_mode": "tsopt-only",
                    "refine_path": refine_path,
                    "tsopt": do_tsopt,
                    "thermo": do_thermo,
                    "dft": do_dft,
                    "opt_mode": tsopt_opt_mode_default,
                    "mep_mode": mep_mode_kind,
                    "mlip_backend": _mlip_backend,
                    "mlip_model": _mlip_model,
                    "command": command_str,
                    "charge": q_int,
                    "spin": spin,
                    "freeze_atoms": _freeze_atoms_for_log(),
                    "mep": {"n_images": n_images, "n_segments": 1},
                    "segments": summary.get("segments", []),
                    "energy_diagrams": summary.get("energy_diagrams", []),
                    "post_segments": [segment_log],
                    "key_files": {},
                }
                # Copy R/TS/P to seg_01/ BEFORE writing summary.log (so tree includes seg_01/)
                try:
                    _input_suffix = first_input.suffix.lower() if first_input else ".xyz"
                    _tsonly_structs = {"R": pR, "TS": pT, "P": pP}
                    _seg_out = _copy_structures_to_seg_dir(
                        _tsonly_structs, out_dir, 1, _input_suffix,
                    )
                    _echo(f"[all] Wrote R/TS/P for segment 01 → {_seg_out}", narrative=True)
                except Exception as e:
                    _echo(
                        f"[all] WARNING: Failed to copy structures in TSOPT-only mode: {e}",
                        err=True,
                    )

                # Refresh summary.json with post_segments and key_output_files
                summary["post_segments"] = _json_safe([segment_log])
                # Rebuild key_output_files now that seg_01/ exists
                try:
                    _kf: Dict[str, Any] = {}
                    for _n, _d in [("summary.log", "Human-readable results summary"),
                                   ("irc_plot_all.png", "Aggregated IRC plot")]:
                        if (out_dir / _n).exists():
                            _kf[_n] = _d
                    for _child in sorted(out_dir.iterdir()) if out_dir.exists() else []:
                        if _child.is_dir() and _child.name.startswith("seg_"):
                            _kf[_child.name] = {"files": sorted(f.name for f in _child.iterdir() if f.is_file())}
                    if _kf:
                        summary["key_output_files"] = _kf
                except Exception:
                    pass
                try:
                    with open(tsroot / "summary.json", "w") as f:
                        json.dump(summary, f, indent=2, ensure_ascii=False)
                    _copy_logged(tsroot / "summary.json", out_dir / "summary.json", label="summary.json", echo=False)
                except Exception:
                    pass

                write_summary_log(tsroot / "summary.log", summary_payload)
                _copy_logged(tsroot / "summary.log", out_dir / "summary.log", label="summary.log", echo=False)
            except Exception as e:
                _echo(f"[write] WARNING: Failed to write summary.log in TSOPT-only mode: {e}", err=True)
        except Exception as e:
            _echo(
                f"[write] WARNING: Failed to write summary.json for TSOPT-only run: {e}",
                err=True,
            )

        try:
            for stem in (
                "energy_diagram_MLIP",
                "energy_diagram_G_MLIP",
                "energy_diagram_DFT",
                "energy_diagram_G_DFT_plus_MLIP",
            ):
                src = tsroot / f"{stem}.png"
                if src.exists():
                    dst = out_dir / f"{stem}_all.png"
                    _copy_logged(src, dst, label=src.name)
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to mirror *_all diagrams in TSOPT-only mode: {e}",
                err=True,
            )

        try:
            if isinstance(irc_plot_path, Path) and irc_plot_path.exists():
                dst = out_dir / "irc_plot_all.png"
                _copy_logged(irc_plot_path, dst, label="irc_plot_all.png")
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to mirror IRC plot in TSOPT-only mode: {e}",
                err=True,
            )

        _echo_section(
            "====== [all] TSOPT-only pipeline successfully finished ======"
        )
        _emit_final_summary(out_dir, time_start)
        return

    # Stage 1b: optional staged scan (single-structure)
    models_for_path: List[Path]
    model_ref_pdbs: Optional[List[Path]] = None
    if is_single and has_scan:
        _echo_section("====== [all] Stage 1b — Staged scan on input ======")
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
            scan_input_pdb = Path(model_outputs[0]).resolve()
            full_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _convert_scan_lists_to_model_indices(
                scan_lists_raw, full_input_pdb, scan_input_pdb
            )
            _echo_detail(
                "[all] Remapped --scan-lists indices from the full PDB to the active site model ordering."
            )

        # Respect the user's explicit --scan-one-based choice symmetrically across
        # both branches. Default to 1-based when unspecified. For extracted models the
        # indices are already 1-based and we decrement to 0-based on explicit False;
        # for skip_extract literals (user-supplied verbatim), an explicit False also
        # decrements so the downstream scan.py sees indices consistent with the
        # --zero-based interpretation we forward below.
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
        scan_endopt_use = False if scan_endopt_override is None else bool(
            scan_endopt_override
        )
        scan_opt_mode_use = opt_mode_norm

        scan_args: List[str] = [
            "-i",
            str(scan_input_pdb),
            "-q",
            str(int(q_int)),
            "-m",
            str(int(spin)),
            "--out-dir",
            str(scan_dir),
            "--freeze-links" if (freeze_links_flag and freeze_ref is not None) else "--no-freeze-links",
            "--convert-files" if convert_files else "--no-convert-files",
            "--preopt" if scan_preopt_use else "--no-preopt",
            "--endopt" if scan_endopt_use else "--no-endopt",
            "--opt-mode",
            str(scan_opt_mode_use),
        ]
        _scan_ref = ref_pdb_for_topology or (scan_input_pdb if scan_input_pdb.suffix.lower() == ".pdb" else None)
        if _scan_ref is not None:
            scan_args.extend(["--ref-pdb", str(_scan_ref)])

        if dump_override_requested:
            scan_args.append("--dump" if dump else "--no-dump")

        if scan_one_based is not None:
            scan_args.append("--one-based" if scan_one_based else "--zero-based")

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        if args_yaml is not None:
            scan_args.extend(["--config", str(args_yaml)])
        if not _forward_calc_file_argv(scan_args, calc_cfg_shared) and cli_param_overridden(ctx, "backend"):
            scan_args.extend(["--backend", str(backend)])
        if cli_param_overridden(ctx, "solvent"):
            scan_args.extend(["--solvent", str(solvent)])
        if cli_param_overridden(ctx, "solvent_model"):
            scan_args.extend(["--solvent-model", str(solvent_model)])
        if scan_stage_literals:
            scan_args.append("--scan-lists")
            scan_args.extend(scan_stage_literals)

        _echo()
        _echo_detail(
            f"[all] dispatch scan: input={scan_input_pdb.name}, "
            f"stages={len(scan_stage_literals)}, preopt={'yes' if scan_preopt_use else 'no'}, "
            f"endopt={'yes' if scan_endopt_use else 'no'}, out={scan_dir}"
        )
        _echo("[all] pdb2reaction scan " + " ".join(scan_args))

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
        _echo_detail("[all] Collected scan stage active site model files:")
        for p in stage_results:
            _echo_detail(f"  - {p}")

        initial_path_for_path = scan_input_pdb
        initial_ref_pdb_for_path = ref_pdb_for_topology or (scan_input_pdb if scan_input_pdb.suffix.lower() == ".pdb" else None)
        if scan_preopt_use:
            preopt_xyz = (scan_dir / "preopt" / "result.xyz").resolve()
            preopt_pdb = (scan_dir / "preopt" / "result.pdb").resolve()
            if preopt_xyz.exists():
                initial_path_for_path = preopt_xyz
                if preopt_pdb.exists():
                    initial_ref_pdb_for_path = preopt_pdb
                _echo_detail(f"[all] Using scan preopt XYZ as initial path endpoint: {initial_path_for_path}")
        models_for_path = [initial_path_for_path] + stage_results

        if initial_ref_pdb_for_path is not None:
            candidate_pdbs: List[Path] = [initial_ref_pdb_for_path]
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
                model_ref_pdbs = candidate_pdbs
            else:
                _echo(
                    "[all] WARNING: active site model PDB snapshots for staged scan were not found; "
                    "full-system merge will use input paths instead.",
                    err=True,
                )
    else:
        if skip_extract:
            # Use post-altloc-fix paths (extract_inputs = final_inputs) so the
            # fix_altloc output is preserved through the skip_extract path.
            # `inputs_for_extract` is the pre-altloc-fix snapshot and would
            # silently discard the altLoc patch.
            models_for_path = [p.resolve() for p in extract_inputs]
        else:
            models_for_path = list(model_outputs)

    # --- Global pre-alignment for coordinate continuity across segments ---
    if not refine_path and len(models_for_path) >= 2:
        _fa = _freeze_atoms_for_log()
        if _fa:
            try:
                _echo("[all] Pre-aligning all input structures to first frame...")
                _align_dir = path_dir / "pre_align"
                ensure_dir(_align_dir)
                _geoms = [geom_loader(str(p), coord_type=DEFAULT_COORD_TYPE) for p in models_for_path]
                for _g in _geoms:
                    _g.freeze_atoms = np.array(_fa, dtype=int)
                _align_calc = create_calculator(**calc_cfg_shared)
                align_and_refine_sequence_inplace(
                    _geoms, shared_calc=_align_calc,
                    out_dir=_align_dir / "refine", verbose=True,
                )
                del _align_calc
                _new_models: List[Path] = []
                for _i, (_g, _orig) in enumerate(zip(_geoms, models_for_path)):
                    _xyz = _align_dir / f"{_i:03d}.xyz"
                    _xyz.write_text(_g.as_xyz() + "\n")
                    if _orig.suffix.lower() == ".pdb":
                        _pdb = _align_dir / f"{_i:03d}.pdb"
                        convert_xyz_to_pdb(_xyz, _orig, _pdb)
                        _new_models.append(_pdb)
                    else:
                        _new_models.append(_xyz)
                models_for_path = _new_models
                _echo("[all] Pre-alignment completed.")
            except Exception as e:
                _echo(
                    f"[all] WARNING: Pre-alignment failed: {e}. "
                    "Continuing with original files.",
                    err=True,
                )

    # Determine availability of full-system templates for downstream merge/copies
    def _is_pdb(path: Path) -> bool:
        return path.suffix.lower() == ".pdb"

    gave_ref_pdb = False

    mep_mode_kind = mep_mode.strip().lower()

    if skip_extract:
        _echo(
            "[all] NOTE: skipping --ref-full-pdb (no --center; inputs already represent full structures)."
        )
    elif is_single and has_scan:
        if _is_pdb(input_paths[0]):
            gave_ref_pdb = True
        else:
            _echo(
                "[all] NOTE: skipping --ref-full-pdb (single+scan: original input is not a PDB)."
            )
    else:
        if all(_is_pdb(p) for p in input_paths):
            gave_ref_pdb = True
        else:
            _echo(
                "[all] NOTE: skipping --ref-full-pdb (one or more original inputs are not PDB)."
            )

    # Stage 2: MEP search
    if not refine_path:
        _echo_section(
            f"====== [all] Stage 2/{stage_total} — Pairwise MEP search via path-opt (no recursive path_search) ======"
        )

        if len(models_for_path) < 2:
            raise click.ClickException("[all] Need at least two structures for path-opt MEP concatenation.")

        combined_blocks: List[str] = []
        path_opt_segments: List[Dict[str, Any]] = []
        for idx, (pL, pR) in enumerate(zip(models_for_path, models_for_path[1:]), start=1):
            # NOTE: internal MEP-engine scratch (3-digit, under _work/); user-facing segment width is 2-digit (segments/seg_NN/).
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
                "--mep-mode",
                mep_mode_kind,
                "--dmf-backend",
                str(dmf_backend).lower(),
                "--max-nodes",
                str(int(max_nodes)),
                "--opt-mode",
                str(opt_mode_norm),
                "--out-dir",
                str(seg_dir),
            ]
            if cli_param_overridden(ctx, "max_cycles"):
                po_args.extend(["--max-cycles", str(int(max_cycles))])
            _append_toggle_arg(po_args, "--freeze-links", bool(freeze_links_flag and freeze_ref is not None))
            if cli_param_overridden(ctx, "climb"):
                _append_toggle_arg(po_args, "--climb", bool(climb))
            _append_toggle_arg(po_args, "--dump", bool(dump))
            _append_toggle_arg(po_args, "--convert-files", bool(convert_files))
            _append_toggle_arg(po_args, "--preopt", bool(preopt))
            ref_pdb_for_seg: Optional[Path] = None
            if model_ref_pdbs and len(model_ref_pdbs) >= idx:
                ref_pdb_for_seg = model_ref_pdbs[idx - 1]
            elif ref_pdb_for_topology is not None:
                ref_pdb_for_seg = ref_pdb_for_topology
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
                po_args.extend(["--config", str(args_yaml)])
            if not _forward_calc_file_argv(po_args, calc_cfg_shared) and cli_param_overridden(ctx, "backend"):
                po_args.extend(["--backend", str(backend)])
            if cli_param_overridden(ctx, "solvent"):
                po_args.extend(["--solvent", str(solvent)])
            if cli_param_overridden(ctx, "solvent_model"):
                po_args.extend(["--solvent-model", str(solvent_model)])

            _echo()
            _echo_detail(
                f"[all] dispatch path-opt segment {idx:02d}/{len(models_for_path) - 1}: "
                f"mode={mep_mode_kind}, preopt={'yes' if preopt else 'no'}, out={seg_dir}"
            )
            _echo("[all] pdb2reaction path-opt " + " ".join(po_args))

            _run_cli_main(
                "path-opt",
                _path_opt.cli,
                po_args,
                on_nonzero="raise",
                prefix=f"all seg {idx:02d}",
            )

            seg_trj = seg_dir / "final_geometries_trj.xyz"
            if not seg_trj.exists():
                raise click.ClickException(
                    f"[all] path-opt segment {idx} did not produce final_geometries_trj.xyz"
                )

            try:
                mirror_dir = path_dir / f"{seg_tag}_mep"
                mirror_trj = mirror_dir / "final_geometries_trj.xyz"

                ensure_dir(mirror_dir)
                if seg_trj.resolve() != mirror_trj.resolve():
                    shutil.copy2(seg_trj, mirror_trj)
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to mirror path-opt trajectory for segment {idx:02d}: {e}",
                    err=True,
                )

            try:
                seg_mep_trj = path_dir / f"mep_seg_{idx:02d}_trj.xyz"
                shutil.copy2(seg_trj, seg_mep_trj)
                if models_for_path[0].suffix.lower() == ".pdb":
                    _path_search._convert_to_pdb_logged(
                        seg_mep_trj,
                        ref_pdb_path=models_for_path[0],
                        out_path=path_dir / f"mep_seg_{idx:02d}.pdb",
                    )
            except Exception as e:
                _echo(
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
                    _echo(
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
                    except Exception as exc:
                        logger.debug("Failed to parse energy from trajectory block: %s", exc)
                        E = np.nan
                energies_seg.append(E)

            first_last = None
            try:
                first_last = xyz_blocks_first_last(raw_blocks, path=seg_trj)
            except Exception as e:
                _echo(
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

        final_trj = path_dir / "mep_trj.xyz"
        try:
            final_trj.write_text("".join(combined_blocks), encoding="utf-8")
            _echo(f"[all] Wrote concatenated MEP trajectory: {final_trj}", narrative=True)
        except Exception as e:
            raise click.ClickException(f"[all] Failed to write concatenated MEP: {e}")

        try:
            run_trj2fig(final_trj, [path_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            close_matplotlib_figures()
            _echo_detail(f"[plot] Saved energy plot → '{path_dir / 'mep_plot.png'}'")
        except Exception as e:
            _echo(f"[plot] WARNING: Failed to plot concatenated MEP: {e}", err=True)

        try:
            if models_for_path[0].suffix.lower() == ".pdb":
                mep_pdb = _path_search._convert_to_pdb_logged(
                    final_trj, ref_pdb_path=models_for_path[0], out_path=path_dir / "mep.pdb"
                )
                if mep_pdb and mep_pdb.exists():
                    dst = out_dir / mep_pdb.name
                    shutil.copy2(mep_pdb, dst)
                    _echo_detail(f"[all] Copied concatenated MEP PDB → {dst}")
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to convert/copy concatenated MEP to PDB: {e}",
                err=True,
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
            _echo(f"[diagram] WARNING: Failed to build GSM diagram for path-opt branch: {e}", err=True)

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
                _echo(
                    f"[all] WARNING: Failed to detect bond changes for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                bond_summary = "(no covalent changes detected)"

            segments_summary.append(
                {
                    "index": seg_idx,
                    "tag": info.get("tag", f"seg_{seg_idx:02d}"),
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
        _enrich_summary(
            summary,
            version="",
            pipeline_mode="path-opt",
            out_dir=out_dir,
            mlip_backend=calc_cfg_shared.get("backend", "uma"),
            mlip_model=calc_cfg_shared.get("model"),
            charge=q_int,
            spin=spin,
            command=command_str,
            config={
                "refine_path": refine_path,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": tsopt_opt_mode_default,
                "mep_mode": mep_mode_kind,
            },
            freeze_atoms=_freeze_atoms_for_log(),
        )
        try:
            with open(path_dir / "summary.json", "w") as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            _echo_detail(f"[write] Wrote '{path_dir / 'summary.json'}'.")
        except Exception as e:
            _echo(f"[write] WARNING: Failed to write summary.json for path-opt branch: {e}", err=True)

        try:
            # MEP deliverables are MOVED to the pipeline root so they are not
            # duplicated under _work/path_opt (path_dir).  The path-opt branch
            # authors the energy diagram lower-cased (energy_diagram_mep.png);
            # glob both casings and land it at root under the canonical name.
            for cand in ("energy_diagram_MEP.png", "energy_diagram_mep.png"):
                src = path_dir / cand
                if src.exists():
                    _move_logged(src, out_dir / "energy_diagram_MEP.png", label=cand)
                    break

            for name in (
                "mep.pdb",
                "mep_w_ref.pdb",
            ):
                src = path_dir / name
                if src.exists():
                    _move_logged(src, out_dir / name, label=name)

            for stem in ("mep", "mep_w_ref"):
                for ext in ("_trj.xyz", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        _move_logged(src, out_dir / src.name, label=src.name)

            # summary.json / summary.log stay COPIES: the path_dir copy is
            # re-read (segments) and re-authored later, and the root copy is
            # consumed by the final-summary banner.
            for name in ("summary.json", "summary.log"):
                src = path_dir / name
                if src.exists():
                    _copy_logged(src, out_dir / name, label=name)
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to relocate path-opt summary files: {e}",
                err=True,
            )
        try:
            diag_for_log: Dict[str, Any] = {}
            for diag in summary.get("energy_diagrams", []) or []:
                if isinstance(diag, dict) and str(diag.get("name", "")).lower().endswith("mep"):
                    diag_for_log = diag
                    break
            # MEP products were moved up to out_dir above; reference the root.
            mep_info = {
                "n_images": summary.get("n_images"),
                "n_segments": summary.get("n_segments"),
                "traj_pdb": str(out_dir / "mep.pdb") if (out_dir / "mep.pdb").exists() else None,
                "mep_plot": str(out_dir / "energy_diagram_MEP.png") if (out_dir / "energy_diagram_MEP.png").exists() else None,
                "diagram": diag_for_log,
            }
            _mlip_backend = calc_cfg_shared.get("backend", "uma")
            _mlip_model = calc_cfg_shared.get("model")
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
                "mlip_backend": _mlip_backend,
                "mlip_model": _mlip_model,
                "command": command_str,
                "charge": q_int,
                "spin": spin,
                "freeze_atoms": _freeze_atoms_for_log(),
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {},
            }
            write_summary_log(path_dir / "summary.log", summary_payload)
            _copy_logged(path_dir / "summary.log", out_dir / "summary.log", label="summary.log")
        except Exception as e:
            _echo(
                f"[write] WARNING: Failed to write summary.log for path-opt branch: {e}",
                err=True,
            )
    if refine_path:
        # --- recursive GSM path_search branch ---
        _echo_section(
            f"====== [all] Stage 2/{stage_total} — MEP search on input structures (recursive GSM) ======"
        )

        ps_args: List[str] = []

        for p in models_for_path:
            ps_args.extend(["-i", str(p)])

        ps_args.extend(["-q", str(q_int)])
        ps_args.extend(["-m", str(int(spin))])

        # Only pass --freeze-links when a PDB reference was available for
        # link-parent detection.  For pure XYZ inputs without --ref-pdb,
        # freeze_ref is None and freeze-links should not be activated.
        _append_toggle_arg(ps_args, "--freeze-links", bool(freeze_links_flag and freeze_ref is not None))
        ps_args.extend(["--mep-mode", mep_mode_kind])
        ps_args.extend(["--max-nodes", str(int(max_nodes))])
        if cli_param_overridden(ctx, "max_cycles"):
            ps_args.extend(["--max-cycles", str(int(max_cycles))])
        if cli_param_overridden(ctx, "climb"):
            _append_toggle_arg(ps_args, "--climb", bool(climb))
        ps_args.extend(["--opt-mode", str(opt_mode_norm)])
        _append_toggle_arg(ps_args, "--dump", bool(dump))
        if thresh is not None:
            ps_args.extend(["--thresh", str(thresh)])
        ps_args.extend(["--out-dir", str(path_dir)])
        _append_toggle_arg(ps_args, "--preopt", bool(preopt))
        _append_toggle_arg(ps_args, "--convert-files", bool(convert_files))
        if args_yaml is not None:
            ps_args.extend(["--config", str(args_yaml)])
        if not _forward_calc_file_argv(ps_args, calc_cfg_shared) and cli_param_overridden(ctx, "backend"):
            ps_args.extend(["--backend", str(backend)])
        if cli_param_overridden(ctx, "solvent"):
            ps_args.extend(["--solvent", str(solvent)])
        if cli_param_overridden(ctx, "solvent_model"):
            ps_args.extend(["--solvent-model", str(solvent_model)])

        if gave_ref_pdb:
            for p in (input_paths if not (is_single and has_scan) else (input_paths[:1] * len(models_for_path))):
                ps_args.extend(["--ref-full-pdb", str(p)])
        # Pass --ref-pdb (active site model PDB snapshots) independently of --ref-full-pdb
        # so that path_search can convert XYZ outputs to PDB even without merge.
        if model_ref_pdbs:
            for p in model_ref_pdbs:
                ps_args.extend(["--ref-pdb", str(p)])
        elif ref_pdb_for_topology is not None:
            # Use CLI --ref-pdb for all inputs when no active site model PDB snapshots are available
            for _ in models_for_path:
                ps_args.extend(["--ref-pdb", str(ref_pdb_for_topology)])

        _echo()
        _echo_detail(
            f"[all] dispatch path-search: inputs={len(models_for_path)}, "
            f"mode={mep_mode_kind}, preopt={'yes' if preopt else 'no'}, "
            f"ref-full={'yes' if gave_ref_pdb else 'no'}, out={path_dir}"
        )
        _echo("[all] pdb2reaction path-search " + " ".join(ps_args))

        _run_cli_main(
            "path_search",
            _path_search.cli,
            ps_args,
            on_nonzero="raise",
            prefix="all",
        )

        try:
            # MEP deliverables are MOVED to the pipeline root so they are not
            # duplicated under _work/path_search (path_dir).  path_search emits
            # the canonical energy_diagram_MEP.png; glob both casings defensively
            # and land it at root under the canonical name.
            for cand in ("energy_diagram_MEP.png", "energy_diagram_mep.png"):
                src = path_dir / cand
                if src.exists():
                    _move_logged(src, out_dir / "energy_diagram_MEP.png", label=cand)
                    break

            for name in (
                "mep.pdb",
                "mep_w_ref.pdb",
            ):
                src = path_dir / name
                if src.exists():
                    _move_logged(src, out_dir / name, label=name)

            for stem in ("mep", "mep_w_ref"):
                for ext in ("_trj.xyz", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        _move_logged(src, out_dir / src.name, label=src.name)

            # summary.json / summary.log stay COPIES: the path_dir copy is
            # re-read (segments) and re-authored later, and the root copy is
            # consumed by the final-summary banner.
            for name in ("summary.json", "summary.log"):
                src = path_dir / name
                if src.exists():
                    _copy_logged(src, out_dir / name, label=name)
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to relocate path_search summary files: {e}",
                err=True,
            )

    # Stage 3: merge to full systems (performed by path_search when enabled)
    _echo_section(f"====== [all] Stage 3/{stage_total} — Merge into full-system templates ======")
    if refine_path and gave_ref_pdb:
        _echo_detail(
            "[all] Merging was carried out by path_search using the original inputs as templates."
        )
        _echo_detail(f"[all] Final products can be found under: {out_dir}")
        _echo_detail(
            "  - mep_w_ref.pdb            (full-system merged trajectory)"
        )
        _echo_detail(f"[all] Raw per-segment merged trajectories stay under: {path_dir}")
        _echo_detail(
            "  - mep_w_ref_seg_XX.pdb     (per-segment merged trajectories for covalent-change segments)"
        )
    elif refine_path:
        _echo_detail(
            "[all] --ref-full-pdb was not provided; full-system merged trajectories are not produced."
        )
        _echo_detail(f"[all] Pocket-only outputs are under: {out_dir}")
    else:
        _echo_detail(
            "[all] path-opt mode produces active site model-level outputs only; full-system merge is not performed."
        )
        _echo_detail(f"[all] Aggregated products are under: {out_dir}")
    _echo_detail("  - summary.json             (segment barriers, ΔE, labels)")
    _echo_detail(
        "  - energy_diagram_MEP.png / energy_diagram.* (MEP energy plot)"
    )
    _echo_section("====== [all] Pipeline (core path) successfully finished ======")

    summary_path = path_dir / "summary.json"
    summary_loaded = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    summary: Dict[str, Any] = summary_loaded if isinstance(summary_loaded, dict) else {}
    segments = _read_summary(summary_path)
    if not energy_diagrams:
        existing_diagrams = summary.get("energy_diagrams", [])
        if isinstance(existing_diagrams, list):
            energy_diagrams.extend(existing_diagrams)

    def _write_pipeline_summary_log(post_segment_logs: Sequence[Dict[str, Any]]) -> None:
        # Payload assembly extracted to pdb2reaction.workflows._all_helpers;
        # the I/O wrapper here keeps the original closure capture + error
        # routing semantics so callers do not change.
        from pdb2reaction.workflows._all_helpers import build_pipeline_summary_payload
        try:
            summary_payload = build_pipeline_summary_payload(
                out_dir=out_dir,
                path_dir=path_dir,
                summary=summary,
                refine_path=refine_path,
                do_tsopt=do_tsopt,
                do_thermo=do_thermo,
                do_dft=do_dft,
                dft_func_basis_use=dft_func_basis_use,
                opt_mode=opt_mode,
                mep_mode_kind=mep_mode_kind,
                mlip_backend=calc_cfg_shared.get("backend", "uma"),
                mlip_model=calc_cfg_shared.get("model"),
                command_str=command_str,
                q_int=q_int,
                spin=spin,
                freeze_atoms=_freeze_atoms_for_log(),
                post_segment_logs=post_segment_logs,
            )
            write_summary_log(path_dir / "summary.log", summary_payload)
            try:
                shutil.copy2(path_dir / "summary.log", out_dir / "summary.log")
                _echo_detail(f"[all] Copied summary.log → {out_dir / 'summary.log'}")
            except (OSError, shutil.Error) as exc:
                logger.debug("Failed to copy summary.log to %s: %s", out_dir, exc)
        except (OSError, KeyError, ValueError, TypeError) as e:
            _echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

    if not (do_tsopt or do_thermo or do_dft):
        if energy_diagrams:
            summary["energy_diagrams"] = list(energy_diagrams)
        _write_pipeline_summary_log([])
        _emit_final_summary(out_dir, time_start)
        return

    # Stage 4: post-processing per reactive segment
    _echo_section(
        f"====== [all] Stage 4/{stage_total} — Post-processing per reactive segment ======"
    )

    if not segments:
        _echo("[post] No segments found in summary; nothing to do.", narrative=True)
        _emit_final_summary(out_dir, time_start)
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
        _echo("[post] No bond-change segments. Skipping TS/thermo/DFT.", narrative=True)
        _emit_final_summary(out_dir, time_start)
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
        _echo_section(f"--- [post] seg_{seg_idx:02d} ({seg_tag}) ---")

        seg_root = path_dir  # MEP-engine scratch root (hei_seg_/mep_seg_/endpoints live here, under _work/)
        seg_dir = out_dir / SEGMENTS_DIRNAME / f"seg_{seg_idx:02d}"  # per-segment deliverables
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
        hei_model_path = _find_with_suffixes(hei_base, [".xyz", ".pdb", ".gjf"])
        if hei_model_path is None:
            _echo(
                f"[post] WARNING: HEI active site model file not found for segment {seg_idx:02d} (searched .pdb/.xyz/.gjf); skipping TSOPT.",
                err=True,
            )
            continue
        ref_pdb_for_seg: Optional[Path] = None
        if hei_model_path.suffix.lower() == ".pdb":
            ref_pdb_for_seg = hei_model_path
        else:
            candidate_ref = hei_base.with_suffix(".pdb")
            if candidate_ref.exists():
                ref_pdb_for_seg = candidate_ref
            elif ref_pdb_for_topology is not None:
                ref_pdb_for_seg = ref_pdb_for_topology

        struct_dir = seg_dir / "structures"
        ensure_dir(struct_dir)
        state_structs: Dict[str, Path] = {}
        uma_ref_energies: Dict[str, float] = {}

        if do_tsopt:
            ts_pdb, g_ts = _run_tsopt_on_hei(
                hei_model_path,
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
                mep_dir=path_dir,
                ref_pdb_for_seg=ts_pdb,
                seg_model_pdb=hei_model_path,
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

            ref_struct_template = ref_pdb_for_seg or hei_model_path
            _save_single_geom_as_pdb_for_tools(
                gL, ref_struct_template, struct_dir, "reactant_irc"
            )
            pT = _save_single_geom_as_pdb_for_tools(
                gT, ref_struct_template, struct_dir, "ts"
            )
            _save_single_geom_as_pdb_for_tools(
                gR, ref_struct_template, struct_dir, "product_irc"
            )

            endpoint_opt_dir = seg_dir / "endpoint_opt"
            ensure_dir(endpoint_opt_dir)

            # Map IRC left/right Hessians → R/P endpoint
            # When reverse_irc is True, _irc_and_match swapped left/right to match GSM endpoints,
            # so "irc_left" (=forward) now corresponds to gR and "irc_right" (=backward) to gL.
            from pdb2reaction.io.hessian_cache import load as _hess_load, store as _hess_store
            _left_hk  = "irc_right" if reverse_irc else "irc_left"
            _right_hk = "irc_left"  if reverse_irc else "irc_right"

            _c = _hess_load(_left_hk)
            if _c:
                _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
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
                _echo(
                    f"[post] WARNING: Reactant endpoint optimization failed for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                g_react_opt = gL

            _c = _hess_load(_right_hk)
            if _c:
                _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
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
                _echo(
                    f"[post] WARNING: Product endpoint optimization failed for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                g_prod_opt = gR

            shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
            _echo_detail("[endpoint-opt] Clean endpoint-opt working dir.")

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
                seg_dir / "energy_diagram_MLIP",
                labels=["R", f"TS{seg_idx}", "P"],
                energies_au=[eR, eT, eP],
                title_note="(MLIP, TSOPT + IRC)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

            tsopt_seg_energies.append((eR, eT, eP))

            segment_log["uma"] = {
                "labels": ["R", f"TS{seg_idx}", "P"],
                "energies_au": [eR, eT, eP],
                "energies_kcal": [0.0, (eT - eR) * AU2KCALPERMOL, (eP - eR) * AU2KCALPERMOL],
                "diagram": str(seg_dir / "energy_diagram_MLIP.png"),
                "structures": state_structs,
                "barrier_kcal": (eT - eR) * AU2KCALPERMOL,
                "delta_kcal": (eP - eR) * AU2KCALPERMOL,
            }

            try:
                _input_suffix = first_input.suffix.lower() if first_input else ".xyz"
                _seg_out = _copy_structures_to_seg_dir(
                    state_structs, out_dir, seg_idx, _input_suffix,
                )
                _echo(f"[all] Wrote R/TS/P for segment {seg_idx:02d} → {_seg_out}", narrative=True)
            except Exception as e:
                _echo(
                    f"[all] WARNING: Failed to copy structures for segment {seg_idx:02d}: {e}",
                    err=True,
                )

        elif do_thermo or do_dft:
            seg_model_path = _find_with_suffixes(
                seg_root / f"mep_seg_{seg_idx:02d}", [".xyz", ".pdb"]
            )

            # Decide reference PDB (if any) for freeze-atoms detection / PDB conversion
            freeze_ref: Optional[Path] = ref_pdb_for_seg
            if freeze_ref is None and seg_model_path is not None and seg_model_path.suffix.lower() == ".pdb":
                freeze_ref = seg_model_path
            elif freeze_ref is None and hei_model_path.suffix.lower() == ".pdb":
                freeze_ref = hei_model_path

            freeze_atoms: List[int] = _get_freeze_atoms(freeze_ref, freeze_links_flag)

            try:
                endpoints = _load_segment_endpoints(seg_root, str(seg_tag), freeze_atoms)
                if endpoints is None:
                    _echo(
                        f"[post] WARNING: final_geometries_trj.xyz not found for segment {seg_idx:02d}; cannot run thermo/DFT without --tsopt. Skipping segment.",
                        err=True,
                    )
                    continue
                gL, gR = endpoints
            except Exception as e:
                _echo(
                    f"[post] WARNING: failed to load segment endpoints from final_geometries_trj.xyz for segment {seg_idx:02d}: {e}. Skipping segment.",
                    err=True,
                )
                continue

            try:
                g_ts = geom_loader(
                    hei_model_path,
                    coord_type=DEFAULT_COORD_TYPE,
                    freeze_atoms=freeze_atoms,
                )
                if freeze_atoms:
                    fa = np.array(freeze_atoms, dtype=int)
                    g_ts.freeze_atoms = fa
            except Exception as e:
                _echo(
                    f"[post] WARNING: failed to load HEI geometry for segment {seg_idx:02d}: {e}. Skipping segment.",
                    err=True,
                )
                continue

            calc_args = dict(calc_cfg_shared)
            calc = create_calculator(**calc_args)
            gL.set_calculator(calc)
            gR.set_calculator(calc)
            g_ts.set_calculator(calc)

            ref_for_structs = ref_pdb_for_seg or (seg_model_path if seg_model_path is not None else hei_model_path)
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

        # ── Release GPU memory before freq/thermo/DFT ──
        # Geometry objects from IRC/endpoint-opt hold calculator references;
        # freq/DFT run as CLI subprocesses and don't need them.
        for _g in (
            locals().get("gL"), locals().get("gR"), locals().get("gT"),
            locals().get("g_ts"), locals().get("g_react_opt"), locals().get("g_prod_opt"),
        ):
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        freq_seg_root = _resolve_override_dir(seg_dir / "freq", freq_out_dir)
        dft_seg_root = _resolve_override_dir(seg_dir / "dft", dft_out_dir)

        if (do_thermo or do_dft) and not state_structs:
            _echo(
                f"[post] WARNING: No segment structures prepared for segment {seg_idx:02d}; skipping thermo/DFT.",
                err=True,
            )
            continue

        p_react = state_structs.get("R")
        p_ts = state_structs.get("TS")
        p_prod = state_structs.get("P")

        if do_thermo:
            if not (p_react and p_ts and p_prod):
                _echo(
                    f"[thermo] WARNING: Missing R/TS/P structures for segment {seg_idx:02d}; skipping thermo.",
                    err=True,
                )
            else:
                _echo()
                _echo_detail(
                    f"[thermo] Segment {seg_idx:02d}: freq on TS/R/P"
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
                from pdb2reaction.io.hessian_cache import clear as _clear_hess_cache
                _clear_hess_cache()  # TS Hessian consumed; R/P need exact computation
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
                        tR.get("sum_EE_and_thermal_free_energy_ha", np.nan)
                    )
                    GT = float(
                        tT.get("sum_EE_and_thermal_free_energy_ha", np.nan)
                    )
                    GP = float(
                        tP.get("sum_EE_and_thermal_free_energy_ha", np.nan)
                    )
                    gibbs_vals = [GR, GT, GP]
                    if all(np.isfinite(gibbs_vals)):
                        diag_payload = _write_segment_energy_diagram(
                            seg_dir / "energy_diagram_G_MLIP",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_au=gibbs_vals,
                            title_note="(MLIP + Thermal Correction)",
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
                            "diagram": str(seg_dir / "energy_diagram_G_MLIP.png"),
                            "structures": state_structs,
                            "barrier_kcal": (GT - GR) * AU2KCALPERMOL,
                            "delta_kcal": (GP - GR) * AU2KCALPERMOL,
                        }
                    else:
                        _echo(
                            "[thermo] NOTE: Gibbs energies non-finite; diagram skipped."
                        )
                except Exception as e:
                    _echo(
                        f"[thermo] WARNING: failed to build Gibbs diagram: {e}",
                        err=True,
                    )

        if do_dft:
            # DO NOT INLINE: (segment path, DFT-pre): same Geometry-retains-calculator mechanism as single-TS path; covers the segment-pipeline branch. Null calculator + del prevents subprocess fork inheriting GPU-pinned model.
            # ── Aggressively release GPU memory before DFT subprocess ──
            for _obj in (
                locals().get("gL"), locals().get("gR"), locals().get("gT"),
                locals().get("g_ts"), locals().get("g_react_opt"), locals().get("g_prod_opt"),
                locals().get("calc"),
            ):
                if _obj is not None:
                    if hasattr(_obj, "calculator"):
                        _obj.calculator = None
                    if hasattr(_obj, "results"):
                        _obj.results = {}
            del _obj
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            if not (p_react and p_ts and p_prod):
                _echo(
                    f"[dft] WARNING: Missing R/TS/P structures for segment {seg_idx:02d}; skipping DFT.",
                    err=True,
                )
            else:
                _echo()
                _echo_detail(f"[dft] Segment {seg_idx:02d}: DFT on R/TS/P")
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
                dR = dft_payloads.get("R") or {}
                dT = dft_payloads.get("TS") or {}
                dP = dft_payloads.get("P") or {}
                eR_dft = _dft_energy_ha(dR)
                eT_dft = _dft_energy_ha(dT)
                eP_dft = _dft_energy_ha(dP)
                _dft_all_ok = all(e is not None for e in (eR_dft, eT_dft, eP_dft))
                if not _dft_all_ok:
                    _failed_states = [s for s, e in zip(["R", "TS", "P"], [eR_dft, eT_dft, eP_dft]) if e is None]
                    _echo(f"[dft] WARNING: DFT failed for state(s): {', '.join(_failed_states)}. Skipping DFT diagrams.", err=True)
                    segment_log["dft"] = {
                        "status": "failed",
                        "failed_states": _failed_states,
                    }
                if _dft_all_ok:
                    try:
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
                    except Exception as e:
                        _echo(
                            f"[dft] WARNING: failed to build DFT diagram: {e}",
                            err=True,
                        )

                if do_thermo and _dft_all_ok:
                    try:
                        dG_R = float(
                            (thermo_payloads.get("R", {}) or {}).get(
                                "thermal_correction_free_energy_ha", np.nan
                            )
                        )
                        dG_T = float(
                            (thermo_payloads.get("TS", {}) or {}).get(
                                "thermal_correction_free_energy_ha", np.nan
                            )
                        )
                        dG_P = float(
                            (thermo_payloads.get("P", {}) or {}).get(
                                "thermal_correction_free_energy_ha", np.nan
                            )
                        )
                        GR_dftUMA = eR_dft + dG_R
                        GT_dftUMA = eT_dft + dG_T
                        GP_dftUMA = eP_dft + dG_P
                        if all(
                            np.isfinite([GR_dftUMA, GT_dftUMA, GP_dftUMA])
                        ):
                            diag_payload = _write_segment_energy_diagram(
                                seg_dir / "energy_diagram_G_DFT_plus_MLIP",
                                labels=["R", f"TS{seg_idx}", "P"],
                                energies_au=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                                title_note=f"({dft_func_basis_use} // MLIP + Thermal Correction)",
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
                                "diagram": str(seg_dir / "energy_diagram_G_DFT_plus_MLIP.png"),
                                "structures": state_structs,
                                "barrier_kcal": (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                                "delta_kcal": (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                            }
                        else:
                            _echo(
                                "[dft//mlip] WARNING: DFT//MLIP Gibbs energies non-finite; diagram skipped.",
                                err=True,
                            )
                    except Exception as e:
                        _echo(
                            f"[dft//mlip] WARNING: failed to build DFT//MLIP Gibbs diagram: {e}",
                            err=True,
                        )

    if tsopt_seg_energies:
        tsopt_all_energies = [e for triple in tsopt_seg_energies for e in triple]
        tsopt_all_labels = _build_global_segment_labels(len(tsopt_seg_energies))
        if tsopt_all_labels and len(tsopt_all_labels) == len(tsopt_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_MLIP_all",
                labels=tsopt_all_labels,
                energies_au=tsopt_all_energies,
                title_note="(MLIP, TSOPT + IRC; all segments)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_thermo and g_uma_seg_energies:
        g_uma_all_energies = [e for triple in g_uma_seg_energies for e in triple]
        g_uma_all_labels = _build_global_segment_labels(len(g_uma_seg_energies))
        if g_uma_all_labels and len(g_uma_all_labels) == len(g_uma_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_G_MLIP_all",
                labels=g_uma_all_labels,
                energies_au=g_uma_all_energies,
                title_note="(MLIP + Thermal Correction; all segments)",
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
                out_dir / "energy_diagram_G_DFT_plus_MLIP_all",
                labels=g_dftuma_all_labels,
                energies_au=g_dftuma_all_energies,
                title_note=f"({dft_func_basis_use} // MLIP + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    # Aggregated IRC plot over all reactive segments (single trj + trj2fig)
    if irc_trj_for_all:
        _merge_irc_trajectories_to_single_plot(
            irc_trj_for_all, out_dir / "irc_plot_all.png"
        )

    # Refresh summary.json with final energy diagram metadata (including aggregated diagrams)
    try:
        summary["energy_diagrams"] = list(energy_diagrams)
        _enrich_summary(
            summary,
            version="",
            pipeline_mode="path-search" if refine_path else "path-opt",
            out_dir=out_dir,
            mlip_backend=calc_cfg_shared.get("backend", "uma"),
            mlip_model=calc_cfg_shared.get("model"),
            charge=q_int,
            spin=spin,
            command=command_str,
            post_segments=post_segment_logs,
            config={
                "refine_path": refine_path,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": tsopt_opt_mode_default,
                "mep_mode": mep_mode_kind,
            },
            freeze_atoms=_freeze_atoms_for_log(),
        )
        with open(path_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        _echo_detail(f"[write] Updated '{path_dir / 'summary.json'}' with energy diagrams.")
        try:
            dst_summary = out_dir / "summary.json"
            shutil.copy2(path_dir / "summary.json", dst_summary)
            _echo_detail(f"[all] Copied summary.json → {dst_summary}")
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to mirror summary.json to {out_dir}: {e}",
                err=True,
            )
        _write_pipeline_summary_log(post_segment_logs)
    except Exception as e:
        _echo(
            f"[write] WARNING: Failed to refresh summary.json with energy diagram metadata: {e}",
            err=True,
        )

    _emit_final_summary(out_dir, time_start)


_configure_all_help_visibility(cli)
