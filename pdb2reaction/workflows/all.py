"""End-to-end extraction, reaction-path, TS, IRC, and analysis workflow."""

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
from pdb2reaction.core.output import emit
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
from pdb2reaction.workflows.align_freeze import (
    align_and_refine_sequence_inplace,
    alignment_failed_pair_indices,
)
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    OUT_DIR_ALL,
    SEGMENTS_DIRNAME,
    THRESH_CHOICES,
    WORK_DIRNAME,
    UMA_CALC_KW as _UMA_CALC_KW,
    fresh_dmf_config,
)
DEFAULT_COORD_TYPE = GEOM_KW_DEFAULT["coord_type"]
from pdb2reaction.io.trj2fig import run_trj2fig
from pdb2reaction.io.summary import (
    emit_method_citations,
    method_references,
    write_summary_log,
)
from pdb2reaction.io.structure_formats import (
    coordinate_template_for,
    register_coordinate_template,
    register_output_template_and_write_cif,
    unregister_coordinate_template,
)
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
    fill_charge_spin_from_gjf_inputs,
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
    cli_param_overridden,
    calculator_provenance,
    verbose_level,
)
from pdb2reaction.core.result_commit import (
    apply_current_run_id,
    commit_exact,
    commit_json,
)
from pdb2reaction.workflows._run_session import (
    CalculatorLease,
    InvocationManifest,
    RunSession,
    claim_public_output as _claim_public_output,
    claim_path_deliverables as _claim_path_deliverables,
    current_key_output_files as _current_key_output_files,
    current_output_paths as _current_output_paths,
    declare_public_output as _declare_public_output,
    declare_path_deliverables as _declare_path_deliverables,
    public_output_key as _public_output_key,
    refresh_current_public_outputs as _refresh_current_public_outputs,
)
from pdb2reaction.cli.common_options import add_coord_type_option, add_precision_option, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
# Advanced-help visibility/callback: one product-local implementation lives in
# cli.help_pages; the `all` command routes through it (a presentation-layer
# dependency) instead of keeping a private copy.
from pdb2reaction.cli.help_pages import _hide_advanced_options, _show_advanced_subcommand_help

logger = logging.getLogger(__name__)


def _validate_postprocessing_dependencies(
    *, do_tsopt: bool, do_thermo: bool, do_dft: bool
) -> None:
    if (do_thermo or do_dft) and not do_tsopt:
        raise click.UsageError(
            "`all --thermo` and `all --dft` require `--tsopt`; an unoptimized "
            "MEP highest-energy image is not a validated transition state."
        )

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
        kwargs.setdefault("narrative", False)
        emit(*args, **kwargs)
        self._started = True

    def section(self, message: str, **kwargs) -> None:
        # Section banners are the narrative spine of the pipeline log: keep
        # them visible at default verbosity (caller may pass narrative=False
        # for a detail-only sub-banner). The leading blank inherits the same
        # tag so the spacing around a shown banner survives the gate.
        narrative = kwargs.setdefault("narrative", True)
        if self._started:
            emit(narrative=narrative)
        emit(message, **kwargs)
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


def _charge_override_message(
    label: str,
    charge: int,
    workflow_charge: Optional[int],
) -> str:
    """Describe an explicit total-charge choice without misstating its relation."""
    message = (
        f"[all] WARNING: {label} supplied; using TOTAL system charge {charge:+d}"
    )
    if workflow_charge is not None:
        relation = "matches" if charge == workflow_charge else "overrides"
        message += f" ({relation} workflow-derived {workflow_charge:+d})"
    return message


def _echo_section(message: str, **kwargs) -> None:
    """Echo a section header (narrative) with a leading blank line unless first."""
    _echo_state.section(message, **kwargs)


def _emit_final_summary(
    out_dir: Path | None,
    time_start: float,
    manifest: Optional[InvocationManifest] = None,
    citation_payload: Optional[Dict[str, Any]] = None,
) -> None:
    """Print a visual `====== Pipeline summary ======` block + Elapsed line.

    Reads ``summary.json`` if present and lifts the most-asked-for numbers
    (status, highest local barrier, reactive-segment count, output dir) so
    the user sees them at the bottom of the log without scrolling back
    through `[diagram] Wrote ...` / `[time] Elapsed Time for X:` clutter.
    Falls back to just the Elapsed line when summary.json is absent
    (dry-run, early failure, TSOPT-only without aggregation).
    """
    summary: Dict[str, Any] = {}
    if out_dir is not None:
        summary_path = (Path(out_dir) / "summary.json").resolve(strict=False)
        current_summary = manifest is None or any(
            path == summary_path for path in manifest.paths("output.public.")
        )
        if current_summary and summary_path.exists():
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
                    f"Highest local barrier: {float(barrier):.2f} kcal/mol (segment {seg_idx}, method {method})",
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
    if citation_payload:
        emit_method_citations(citation_payload)
    _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start), narrative=True)


from pdb2reaction.workflows import scan as _scan_cli
from pdb2reaction.domain.add_elem_info import assign_elements as _assign_elem_info
from pdb2reaction.io.pdb_fix import has_altloc as _has_altloc, fix_altloc_file as _fix_altloc
from pdb2reaction.workflows import irc as _irc_cli




def _copy_logged(src: Path, dst: Path, *, label: Optional[str] = None, echo: bool = True) -> bool:
    """Copy files with consistent warning messages; return success."""
    try:
        if Path(src).name == "summary.json" and Path(dst).name == "summary.json":
            commit_exact(Path(dst), Path(src).read_bytes())
        else:
            shutil.copy2(src, dst)
        template = coordinate_template_for(src)
        if template is not None:
            register_coordinate_template(dst, template)
        if echo:
            shown = label or src.name
            _echo_detail(f"[all] Copied {shown} → {dst}")
        return True
    except Exception as e:
        shown = label or src
        _echo(f"[all] WARNING: Failed to copy {shown} to {dst}: {e}", err=True)
        return False


def _is_pipeline_public_destination(root: Path, path: Path) -> bool:
    """Return whether a lexical destination belongs to the public layout."""

    try:
        _public_output_key(root, path)
    except ValueError:
        return False
    return True


def _write_summary_json(path: Path, summary: dict) -> Path:
    """Publish an aggregate summary atomically with the active run identity."""

    identified = apply_current_run_id(_json_safe(summary))
    summary.clear()
    summary.update(identified)
    return commit_json(Path(path), summary)


def _persist_run_manifest(manifest: InvocationManifest, out_dir: Path) -> Path:
    """Persist current-run ownership under the private ``_work`` tree."""

    return manifest.write_internal(
        Path(out_dir) / WORK_DIRNAME / "_run_manifest.json"
    )


def _publish_manifest_summary(
    path: Path,
    summary: dict,
    *,
    manifest: InvocationManifest,
    key: str,
    out_dir: Path,
) -> Path:
    """Atomically publish and claim one aggregate summary generation."""

    destination = Path(path).resolve(strict=False)
    if key not in manifest.expected:
        manifest.declare(key, [destination])
    summary["run_id"] = manifest.run_id
    published = _write_summary_json(destination, summary)
    manifest.claim_one(key)
    _persist_run_manifest(manifest, out_dir)
    return published


def _move_logged(src: Path, dst: Path, *, label: Optional[str] = None, echo: bool = True) -> bool:
    """Move files with consistent warning messages; return success.

    Used to promote MEP deliverables from the ``_work`` scratch tree up to the
    pipeline root without leaving a duplicate behind. ``shutil.move`` replaces
    an existing destination, so re-runs land cleanly.
    """
    try:
        dst = Path(dst)
        template = coordinate_template_for(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        # shutil.move treats an existing destination *directory* as a target to
        # move into; clear any existing destination *file* first so the result
        # is always the intended path (idempotent across re-runs).
        if dst.is_file():
            dst.unlink()
        shutil.move(str(src), str(dst))
        if template is not None:
            register_coordinate_template(dst, template)
            unregister_coordinate_template(src)
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
) -> Optional[int]:
    """Run a Click command with a temporary argv and consistent error handling.

    Returns the child exit code (0 on success, nonzero on a warned failure, 1 on
    a swallowed exception, ``None`` when unknown). Callers that need explicit
    stage results — e.g. FREQ/DFT thermochemistry — must gate on this code rather
    than infer success from an output file's existence.
    """
    saved = list(sys.argv)
    label = prefix or cmd_name
    exit_code: Optional[int] = 0
    # In-proc subcommand dispatch — flag the child's banner / device echo to
    # stay silent so a 4-stage `all` pipeline doesn't reprint the same lines
    # `pdb2reaction ver. X` / `[calc] Resolved device: cuda` once per stage.
    from pdb2reaction.core.utils import is_child_mode, set_child_mode
    prior_child_mode = is_child_mode()
    set_child_mode(True)
    try:
        sys.argv = ["pdb2reaction", cmd_name] + list(args)
        _echo("")
        cli_obj.main(args=list(args), standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            exit_code = code if isinstance(code, int) else 1
            if on_nonzero == "raise":
                raise click.ClickException(f"[{label}] {cmd_name} exit code {code}.")
            _echo(f"[{label}] WARNING: {cmd_name} exited with code {code}", err=True)
    except Exception as e:
        exit_code = 1
        if on_exception == "raise":
            raise click.ClickException(f"[{label}] {cmd_name} failed: {e}")
        _echo(f"[{label}] WARNING: {cmd_name} failed: {e}", err=True)
    finally:
        sys.argv = saved
        set_child_mode(prior_child_mode)
        # Break calculator cycles before reclaiming the allocator cache between stages.
        # Release GPU memory between pipeline stages to prevent OOM.
        # gc.collect() breaks cyclic refs inside torch.nn.Module.
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        _echo("")
    return exit_code


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


def _append_explicit_child_runtime_argv(
    child_args: List[str],
    *,
    workers: Optional[int] = None,
    workers_per_node: Optional[int] = None,
    dmf_backend: Optional[str] = None,
    backend: Optional[str] = None,
    solvent: Optional[str] = None,
    solvent_model: Optional[str] = None,
    max_nodes: Optional[int] = None,
    preopt: Optional[bool] = None,
) -> None:
    """Forward only parent-explicit runtime values to a child command.

    Omitted values stay absent so the child can resolve its own YAML/default
    precedence.  Explicit values, including ``1``, must be forwarded because
    they may intentionally override a different value in the shared YAML.
    """
    if dmf_backend is not None:
        child_args.extend(["--dmf-backend", str(dmf_backend).lower()])
    if workers is not None:
        child_args.extend(["--workers", str(int(workers))])
    if workers_per_node is not None:
        child_args.extend(["--workers-per-node", str(int(workers_per_node))])
    if backend is not None:
        child_args.extend(["--backend", str(backend).lower()])
    if solvent is not None:
        child_args.extend(["--solvent", str(solvent).lower()])
    if solvent_model is not None:
        child_args.extend(["--solvent-model", str(solvent_model).lower()])
    if max_nodes is not None:
        child_args.extend(["--max-nodes", str(int(max_nodes))])
    _append_toggle_arg(child_args, "--preopt", preopt)


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
    # Apply YAML first (overrides defaults)
    if yaml_cfg:
        apply_yaml_overrides(
            yaml_cfg,
            [
                (cfg, (("calc",),)),
            ],
        )
    # CLI explicit values override YAML
    cfg["charge"] = int(charge)
    cfg["spin"] = int(spin)
    if workers is not None:
        cfg["workers"] = int(workers)
    if workers_per_node is not None:
        cfg["workers_per_node"] = int(workers_per_node)
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
    *,
    one_based: bool = True,
) -> List[List[Tuple[int, int, float]]]:
    """Parse ``--scan-lists`` literals, interpreting them per ``one_based``.

    ``one_based`` is how to READ the user's literal (it follows ``--scan-one-based``); the
    returned stages are always 1-based so the downstream index mapping stays base-agnostic.
    """
    stages: List[List[Tuple[int, int, float]]] = []
    for idx_stage, literal in enumerate(scan_lists_raw, start=1):
        tuples, _ = parse_scan_list_triples(
            literal,
            one_based=one_based,
            atom_meta=atom_meta,
            option_name=f"--scan-lists #{idx_stage}",
            return_one_based=True,
        )
        if not tuples:
            raise click.BadParameter(
                f"--scan-lists #{idx_stage} must contain at least one (i,j,target) triple."
            )
        if any(len(entry) != 3 for entry in tuples):
            raise click.BadParameter(
                "pdb2reaction all accepts only (i,j,target) scan triples; "
                "use standalone scan for bidirectional (i,j,start,end) stages."
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
    *,
    one_based: bool = True,
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
    stages = _parse_scan_lists_literals(
        scan_lists_raw, atom_meta=full_atom_meta, one_based=one_based
    )

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
            "Increase extraction coverage (e.g., --radius/--radius-het2het, --selected-resn, or set --no-exclude-backbone), "
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
    tr_projection: Optional[str] = None,
    precision: Optional[str] = None,
    backend_model: Optional[str] = None,
    session: Optional[RunSession] = None,
) -> Optional[Path]:
    """
    Write ``freeze_atoms`` and (optionally) ``coord_type`` /
    ``tr_projection`` / ``precision`` into a
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
    if (
        not freeze_atoms
        and coord_type is None
        and tr_projection is None
        and precision is None
        and backend_model is None
    ):
        return args_yaml
    if session is None:
        raise ValueError("A RunSession is required when generating child YAML.")

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
    if tr_projection is not None:
        geom_cfg["tr_projection"] = tr_projection
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
    session.resources.own_path(tmp_dir)
    out_path = tmp_dir / "args_freeze_atoms.yaml"
    with out_path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)

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


def _read_path_opt_segment_converged(seg_dir: Path) -> Optional[bool]:
    """Read a path-opt segment child's reported MEP convergence (tri-state).

    Reads the additive ``stage_outcomes`` leaf ``converged`` bit from the child's
    ``result.json`` — sourced from the optimizer for both GSM and DMF (IPOPT
    bit), independent of the byte-compatible legacy ``converged`` field.
    Returns ``None`` when no readable signal exists (fail-closed: a missing or
    unreadable child result never promotes the segment to converged).
    """
    try:
        rj = seg_dir / "result.json"
        if not rj.exists():
            return None
        data = json.loads(rj.read_text(encoding="utf-8")) or {}
        for leaf in data.get("stage_outcomes") or []:
            if isinstance(leaf, dict) and isinstance(leaf.get("converged"), bool):
                return bool(leaf["converged"])
        return None
    except Exception as exc:
        logger.debug("Failed to read path-opt segment convergence %s: %s", seg_dir, exc)
        return None


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
    manifest: Optional[InvocationManifest] = None,
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

    def declare_destination(path: Path) -> None:
        if manifest is not None:
            _declare_public_output(
                manifest,
                out_dir,
                path,
            )

    def claim_destination(path: Path) -> None:
        if manifest is not None:
            _claim_public_output(
                manifest,
                out_dir,
                path,
            )

    def copy_destination(src: Path, dst: Path) -> None:
        declare_destination(dst)
        shutil.copy2(src, dst)
        claim_destination(dst)

    for key, src_xyz in state_structs.items():
        src = Path(src_xyz)
        if not src.exists():
            continue
        dst_name = name_map.get(key, key.lower())

        if input_suffix in {".pdb", ".cif", ".mmcif"}:
            src_pdb = src.with_suffix(".pdb")
            if src_pdb.exists():
                dst_pdb = seg_dir / f"{dst_name}.pdb"
                declare_destination(dst_pdb)
                if not _copy_logged(src_pdb, dst_pdb, echo=False):
                    continue
                claim_destination(dst_pdb)
                template = coordinate_template_for(src_pdb)
                if template is not None:
                    dst_cif = dst_pdb.with_suffix(".cif")
                    declare_destination(dst_cif)
                    register_output_template_and_write_cif(dst_pdb, template)
                    claim_destination(dst_cif)
            else:
                copy_destination(src, seg_dir / f"{dst_name}.xyz")
        elif input_suffix == ".gjf":
            if (
                prepared_input is not None
                and getattr(prepared_input, "gjf_template", None) is not None
            ):
                dst_gjf = seg_dir / f"{dst_name}.gjf"
                declare_destination(dst_gjf)
                try:
                    convert_xyz_to_gjf(
                        src,
                        prepared_input.gjf_template,
                        dst_gjf,
                    )
                    claim_destination(dst_gjf)
                except Exception:
                    copy_destination(src, seg_dir / f"{dst_name}.xyz")
            else:
                copy_destination(src, seg_dir / f"{dst_name}.xyz")
        else:
            copy_destination(src, seg_dir / f"{dst_name}.xyz")

    return seg_dir


def _is_reactive_segment(item: Any) -> bool:
    """Return whether a segment legitimately requires TS post-processing."""
    if not isinstance(item, dict):
        return False
    kind = item.get("kind", "seg")
    if kind == "tsopt":
        return True
    if kind != "seg":
        return False
    # Legacy/directly constructed segment records predate bond-change
    # serialization and remain reactive.  Only an explicit no-change result
    # suppresses post-processing.
    if "bond_changes" not in item:
        return True
    changes = str(item.get("bond_changes", "")).strip()
    return bool(changes and changes != "(no covalent changes detected)")


def _derive_pipeline_status(
    summary: dict,
    *,
    post_segments: Optional[list],
    config: Optional[dict],
) -> Tuple[str, List[str]]:
    """Return the aggregate pipeline status and machine-readable reasons.

    A usable MEP is the baseline for ``success``.  When the final
    post-processing records are available, every requested optional stage must
    also have produced its expected result, and a thermochemistry run must
    confirm exactly one TS imaginary mode.  Optional-stage failures therefore
    remain visible even when an earlier MLIP energy diagram exists.
    """

    segments = summary.get("segments") or []
    has_diagrams = bool(summary.get("energy_diagrams"))
    reasons: List[str] = []

    if not segments and not has_diagrams:
        return "failed", ["no usable path segments or energy diagrams were produced"]
    if not has_diagrams:
        reasons.append("no usable energy diagram was produced")

    cfg = config or {}
    requested = any(bool(cfg.get(name)) for name in ("tsopt", "thermo", "dft"))

    # ``None`` means this is an intermediate summary written before the final
    # per-segment records exist.  An explicit list, including an empty one,
    # means post-processing is complete and can be validated.
    if post_segments is not None and requested:
        logs = [item for item in post_segments if isinstance(item, dict)]
        reactive_ids = {
            item.get("index")
            for item in segments
            if _is_reactive_segment(item) and item.get("index") is not None
        }
        if not logs and reactive_ids:
            reasons.append("requested post-processing produced no segment records")

        for ordinal, item in enumerate(logs, start=1):
            index = item.get("index", ordinal)
            prefix = f"segment {index}"

            if cfg.get("tsopt"):
                if not isinstance(item.get("mlip"), dict):
                    reasons.append(f"{prefix}: TSOPT/IRC refined MLIP energies are missing")
                if not item.get("irc_traj"):
                    reasons.append(f"{prefix}: IRC trajectory is missing")

            if cfg.get("thermo"):
                if not isinstance(item.get("gibbs_mlip"), dict):
                    reasons.append(f"{prefix}: MLIP thermochemistry result is missing")
                ts_imag = item.get("ts_imag")
                if not isinstance(ts_imag, dict) or ts_imag.get("n_imag") is None:
                    reasons.append(f"{prefix}: TS imaginary-mode validation is missing")
                elif int(ts_imag["n_imag"]) != 1:
                    reasons.append(
                        f"{prefix}: TS imaginary-mode validation found "
                        f"n_imag={int(ts_imag['n_imag'])}, expected 1"
                    )

            if cfg.get("dft"):
                dft = item.get("dft")
                if not isinstance(dft, dict):
                    reasons.append(f"{prefix}: DFT result is missing")
                elif dft.get("status") == "failed":
                    failed_states = dft.get("failed_states") or []
                    detail = f" ({', '.join(map(str, failed_states))})" if failed_states else ""
                    reasons.append(f"{prefix}: DFT failed{detail}")

                if cfg.get("thermo") and not isinstance(item.get("gibbs_dft_mlip"), dict):
                    reasons.append(f"{prefix}: DFT//MLIP thermochemistry result is missing")

    # The TS-only branch records the DFT aggregate before its post-segment
    # payload is assembled.  Preserve that explicit failure signal in both the
    # intermediate and final summary.
    if cfg.get("dft") and cfg.get("dft_status") == "failed":
        reasons.append("DFT failed for one or more TS-only states")

    # De-duplicate while preserving diagnostic order.
    reasons = list(dict.fromkeys(reasons))
    return ("partial" if reasons else "success"), reasons


_STATUS_SEVERITY = {"success": 0, "partial": 1, "failed": 2}


def _pipeline_aggregate_truth(
    summary: dict,
    *,
    post_segments: Optional[list],
    config: Optional[dict],
    legacy_status: str,
    legacy_reasons: Optional[Sequence[str]] = None,
):
    """Compose the ``all``-pipeline aggregate from per-segment leaves.

    One required :class:`LeafOutcome` is built per reactive segment.  A path
    segment is usable only when its MEP and every post-processing convergence
    signal are explicitly ``True``.  A direct TSOPT segment has no MEP stage, so
    it is gated by its IRC and, when present, both endpoint optimizations.  A
    dict-present / trajectory-present but nonconverged leaf never counts toward
    fail-closed completeness — a never_stop / max-cycle IRC therefore
    cannot yield ``scientific_status == "success"``.

    The convergence-gated aggregate is then composed with the legacy completeness
    axis (``legacy_status`` from :func:`_derive_pipeline_status`, which already
    covers DFT / thermo / n_imag): ``scientific_status`` is the MORE severe of
    the two so the new field carries at least as much information as the legacy
    ``status``.
    The legacy ``status`` string itself is untouched (byte-compatible).
    """

    from pdb2reaction.workflows._outcomes import (
        AggregateTruth,
        aggregate_workflow_truth,
        make_leaf,
    )

    segments = summary.get("segments") or []
    reactive = [s for s in segments if _is_reactive_segment(s)]
    cfg = config or {}
    tsopt_requested = bool(cfg.get("tsopt"))
    legacy_reasons = list(legacy_reasons or [])

    post_by_idx: Dict[Any, dict] = {}
    if post_segments is not None:
        for ps in post_segments:
            if isinstance(ps, dict) and ps.get("index") is not None:
                post_by_idx[ps.get("index")] = ps

    def _and3(a: Optional[bool], b: Optional[bool]) -> Optional[bool]:
        if a is False or b is False:
            return False
        if a is None or b is None:
            return None
        return True

    leaves: List[Any] = []
    expected: List[str] = []
    for s in reactive:
        idx = s.get("index")
        if idx is None:
            continue
        seg_id = f"segment_{idx}"
        expected.append(seg_id)
        post = post_by_idx.get(idx)
        reason = ""
        artifacts: List[str] = []
        # The segment's own reported convergence, threaded from path_search's
        # SegmentReport / the path-opt child leaf.  A missing field is None
        # (fail-closed), never a silent True.
        _seg_conv = s.get("converged")
        seg_converged: Optional[bool] = _seg_conv if isinstance(_seg_conv, bool) else None
        # ``kind=tsopt`` is the direct-TS branch: no MEP child runs and therefore
        # no MEP convergence field exists.  Only that explicit kind bypasses the
        # MEP gate; unknown/future kinds remain fail-closed like path segments.
        mep_converged: Optional[bool] = (
            True if s.get("kind") == "tsopt" else seg_converged
        )
        if post is not None:
            # Post-processing ran: compose its explicit IRC / endpoint records
            # with the MEP engine's own convergence.  Successful downstream
            # work must never promote a nonconverged/unknown path segment.
            converged: Optional[bool] = mep_converged
            if converged is not True:
                reason = (
                    "mep_not_converged"
                    if converged is False
                    else "mep_convergence_unknown"
                )
            irc = post.get("irc")
            if isinstance(irc, dict):
                _u = irc.get("usable")
                _irc_conv = True if _u is True else (False if _u is False else None)
                converged = _and3(converged, _irc_conv)
                if _irc_conv is not True and not reason:
                    reason = f"irc:{irc.get('reason') or 'not_usable'}"
                _traj = irc.get("traj")
                if _traj:
                    artifacts.append(str(_traj))
            elif tsopt_requested:
                # IRC requested but no directional record: fail closed
                # rather than trust the trajectory file's existence.
                converged = _and3(converged, None)
                if not reason:
                    reason = "irc_missing"
            if tsopt_requested and s.get("kind") != "tsopt":
                endpoint_assignment = post.get("endpoint_assignment")
                connectivity = (
                    endpoint_assignment.get("connectivity_validated")
                    if isinstance(endpoint_assignment, dict)
                    else None
                )
                connectivity_truth = (
                    connectivity if isinstance(connectivity, bool) else None
                )
                converged = _and3(converged, connectivity_truth)
                if connectivity_truth is not True and not reason:
                    reason = "irc_endpoint_connectivity_unvalidated"
            eo = post.get("endpoint_opt")
            if isinstance(eo, dict):
                for _k in ("reactant_converged", "product_converged"):
                    if _k in eo:
                        _v = eo.get(_k)
                        converged = _and3(converged, _v if isinstance(_v, bool) else None)
                        if not (isinstance(_v, bool) and _v) and not reason:
                            reason = f"endpoint_opt:{_k}"
        elif tsopt_requested:
            # tsopt was requested but this segment's IRC/endpoint post-processing
            # has not run yet (the intermediate MEP summary is written before
            # post-processing). Fail closed rather than promote a reactive leaf on
            # the MEP trajectory's existence alone: a required post stage is
            # missing, so the segment is not yet usable.
            converged = None
            reason = "post_missing"
        else:
            # Path-only final summary (no tsopt): the segment's own reported
            # convergence is the whole truth. A missing/unknown field fails
            # closed (None) — never default to True.
            converged = mep_converged
            if converged is not True and not reason:
                reason = "not_converged" if converged is False else "convergence_unknown"
        leaves.append(
            make_leaf(
                "all",
                seg_id,
                required=True,
                executed=True,
                converged=converged,
                reason=reason,
                artifacts=artifacts,
            )
        )

    if leaves:
        agg = aggregate_workflow_truth(leaves, expected)
        agg_sci = agg.scientific_status
        agg_exec = agg.execution_status
        agg_reasons = list(agg.status_reasons)
        observed = list(agg.observed_item_ids)
    else:
        # No reactive-segment leaves to gate on (degenerate/endpoint-only
        # summary): mirror the legacy completeness axis rather than manufacture a
        # spurious failure.
        agg_sci = legacy_status
        agg_exec = "failed" if legacy_status == "failed" else "completed"
        agg_reasons = []
        observed = list(expected)

    # Compose with the legacy completeness axis: keep the MORE severe verdict.
    if _STATUS_SEVERITY.get(legacy_status, 0) >= _STATUS_SEVERITY.get(agg_sci, 0):
        scientific = legacy_status
    else:
        scientific = agg_sci
    execution = "failed" if (legacy_status == "failed" or agg_exec == "failed") else "completed"
    reasons = legacy_reasons + [r for r in agg_reasons if r not in legacy_reasons]

    return AggregateTruth(
        execution_status=execution,
        scientific_status=scientific,
        status_reasons=tuple(reasons),
        expected_item_ids=tuple(expected),
        observed_item_ids=tuple(observed),
    )


def _apply_pipeline_truth(
    summary: dict,
    *,
    post_segments: Optional[list],
    config: Optional[dict],
    legacy_status: str,
    legacy_reasons: Optional[Sequence[str]] = None,
) -> None:
    """Write the outcome axes onto ``summary`` in place.

    Never touches the legacy overloaded ``status`` field; only adds
    ``execution_status`` / ``scientific_status`` / expected+observed IDs and the
    distinct ``scientific_status_reasons`` key.
    """

    truth = _pipeline_aggregate_truth(
        summary,
        post_segments=post_segments,
        config=config,
        legacy_status=legacy_status,
        legacy_reasons=legacy_reasons,
    )
    summary["execution_status"] = truth.execution_status
    summary["scientific_status"] = truth.scientific_status
    summary["expected_item_ids"] = list(truth.expected_item_ids)
    summary["observed_item_ids"] = list(truth.observed_item_ids)
    if truth.scientific_status != "success" and truth.status_reasons:
        summary["scientific_status_reasons"] = list(truth.status_reasons)
    else:
        summary.pop("scientific_status_reasons", None)


def _enrich_summary(
    summary: dict,
    *,
    version: str,
    pipeline_mode: str,
    mlip_backend: str,
    mlip_model: Optional[str],
    mlip_precision: Optional[str] = None,
    charge: int,
    spin: int,
    command: str = "",
    post_segments: Optional[list] = None,
    config: Optional[dict] = None,
    freeze_atoms: Optional[str] = None,
    out_dir: Optional[Path] = None,
    manifest: Optional[InvocationManifest] = None,
) -> dict:
    """Add machine-readable metadata to summary dict for AI agent consumption.

    The resulting dict is written as summary.json and is the machine-readable
    pipeline output consumed by MCP tools and other clients. Formatted tables
    and the filesystem tree remain specific to summary.log.
    """
    # Import through the package fallback: _version.py is generated only by a
    # build/install and is intentionally absent from a fresh source checkout.
    from pdb2reaction import __version__
    from pdb2reaction.core.utils import RESULT_JSON_SCHEMA_VERSION

    segments = summary.get("segments", [])
    reactive = [s for s in segments if _is_reactive_segment(s)]
    n_reactive = len(reactive)
    ts_only = pipeline_mode == "tsopt-only"

    status, status_reasons = _derive_pipeline_status(
        summary,
        post_segments=post_segments,
        config=config,
    )

    # Legacy ``rate_limiting_step`` key: highest independently referenced local
    # barrier among reactive segments. This is not a kinetic RLS assignment.
    # Select the highest-level method that covers every reactive segment; PNG
    # creation and diagram list order must not determine the scientific method.
    best_method = None
    rls = None
    if reactive:
        post_by_idx = {
            ps.get("index"): ps
            for ps in (post_segments or [])
            if isinstance(ps, dict)
        }
        method_key = None
        for candidate_method, candidate_key in (
            ("DFT//MLIP_Gibbs", "gibbs_dft_mlip"),
            ("DFT", "dft"),
            ("MLIP_Gibbs", "gibbs_mlip"),
            ("MLIP", "mlip"),
        ):
            covered = True
            for segment in reactive:
                block = (post_by_idx.get(segment.get("index")) or {}).get(candidate_key)
                barrier = block.get("barrier_kcal") if isinstance(block, dict) else None
                try:
                    covered = barrier is not None and bool(np.isfinite(float(barrier)))
                except (TypeError, ValueError):
                    covered = False
                if not covered:
                    break
            if covered:
                best_method = candidate_method
                method_key = candidate_key
                break
        if best_method is None:
            best_method = "MLIP" if ts_only else "MEP"

        # Prefer the refined TSOPT+IRC barrier (post_segments' mlip/gibbs_mlip) matching
        # best_method; segments[*].barrier_kcal is only the un-refined MEP band and must
        # NOT be reported under an MLIP/Gibbs label (the two can differ by several kcal/mol).
        # An agent reading this field gets the TS/IRC-refined ΔE‡/ΔG‡, with the raw MEP
        # value kept alongside as `mep_barrier_kcal` for transparency.
        max_barrier = -1e9
        for s in reactive:
            idx = s.get("index")
            refined = (post_by_idx.get(idx) or {}).get(method_key) if method_key else None
            if isinstance(refined, dict) and refined.get("barrier_kcal") is not None:
                b = refined.get("barrier_kcal") or 0
                cur_method = best_method
            else:
                # In TS-only mode segments[*] is itself the refined
                # TS-minus-assigned-endpoint result; no MEP was run.  Ordinary
                # path segments still fall back to the raw MEP band.
                b = s.get("barrier_kcal", 0) or 0
                cur_method = (
                    "MLIP"
                    if ts_only and s.get("kind") == "tsopt"
                    else "MEP"
                )
            if b > max_barrier:
                max_barrier = b
                rls = {
                    "segment": idx,
                    "barrier_kcal": round(b, 2),
                    "method": cur_method,
                }
                if not (ts_only and s.get("kind") == "tsopt"):
                    rls["mep_barrier_kcal"] = round(
                        s.get("barrier_kcal", 0) or 0, 2
                    )

    # Overall reaction energy from the best all-segment diagram
    overall_rxn_e = None
    overall_rxn_method = None
    diagrams_by_name = {
        str(diag.get("name", "")): diag
        for diag in summary.get("energy_diagrams", [])
        if isinstance(diag, dict)
    }
    method_rank = {
        "MEP": 0,
        "MLIP": 1,
        "MLIP_Gibbs": 2,
        "DFT": 3,
        "DFT//MLIP_Gibbs": 4,
    }
    max_overall_rank = method_rank.get(best_method or "MEP", 0)
    for diagram_name, method in (
        ("energy_diagram_G_DFT_plus_MLIP_all", "DFT//MLIP_Gibbs"),
        ("energy_diagram_G_DFT_plus_MLIP", "DFT//MLIP_Gibbs"),
        ("energy_diagram_DFT_all", "DFT"),
        ("energy_diagram_DFT", "DFT"),
        ("energy_diagram_G_MLIP_all", "MLIP_Gibbs"),
        ("energy_diagram_G_MLIP", "MLIP_Gibbs"),
        ("energy_diagram_MLIP_all", "MLIP"),
        ("energy_diagram_MLIP", "MLIP"),
        ("energy_diagram_MEP", "MEP"),
        ("MEP", "MEP"),
    ):
        if method_rank[method] > max_overall_rank:
            continue
        diag = diagrams_by_name.get(diagram_name)
        energies = diag.get("energies_kcal", []) if diag else []
        if len(energies) >= 2:
            try:
                first, last = float(energies[0]), float(energies[-1])
            except (TypeError, ValueError):
                continue
            if np.isfinite(first) and np.isfinite(last):
                overall_rxn_e = round(last - first, 2)
                overall_rxn_method = method
                break

    # Core metadata
    summary["pdb2reaction_version"] = __version__
    summary["schema_version"] = RESULT_JSON_SCHEMA_VERSION
    summary["pipeline_mode"] = pipeline_mode
    summary["status"] = status
    if status_reasons:
        summary["status_reasons"] = status_reasons
    else:
        summary.pop("status_reasons", None)
    # Keep the legacy overloaded ``status`` intact and
    # expose the execution/scientific split plus expected/observed segment IDs so
    # a forward-compatible consumer can tell "the pipeline ran" from "the science
    # is complete and usable". ``scientific_status`` is computed from explicit
    # per-segment LeafOutcomes (IRC directional + endpoint-opt convergence)
    # composed with the legacy completeness axis, so a nonconverged IRC/endpoint
    # leaf whose trajectory still exists cannot make the pipeline a success.
    _apply_pipeline_truth(
        summary,
        post_segments=post_segments,
        config=config,
        legacy_status=status,
        legacy_reasons=status_reasons,
    )
    summary["mlip_backend"] = mlip_backend
    summary["mlip_model"] = mlip_model
    summary["mlip_precision"] = mlip_precision
    summary["charge"] = charge
    summary["spin"] = spin
    if manifest is not None:
        summary["run_id"] = manifest.run_id
    summary["n_segments_reactive"] = n_reactive
    if rls:
        summary["rate_limiting_step"] = rls
    else:
        summary.pop("rate_limiting_step", None)
    if overall_rxn_e is not None:
        summary["overall_reaction_energy_kcal"] = overall_rxn_e
        summary["overall_reaction_energy_method"] = overall_rxn_method
    else:
        summary.pop("overall_reaction_energy_kcal", None)
        summary.pop("overall_reaction_energy_method", None)
    if command:
        summary["command"] = command

    # Pipeline config (mirrors summary.log header)
    if config:
        summary["config"] = config
    citation_config = config or {}
    summary["references"] = method_references(
        {
            "pipeline_mode": pipeline_mode,
            "opt_mode": citation_config.get("opt_mode"),
            "opt_mode_post": citation_config.get("opt_mode_post"),
            "path_opt_mode": citation_config.get("path_opt_mode"),
            "post_opt_mode": citation_config.get("post_opt_mode"),
            "ts_opt_mode": citation_config.get("ts_opt_mode"),
            "endpoint_opt_mode": citation_config.get("endpoint_opt_mode"),
            "mep_mode": citation_config.get("mep_mode"),
            "dmf_correlated": citation_config.get("dmf_correlated"),
            "post_segments": post_segments or [],
        }
    )
    if freeze_atoms:
        summary["freeze_atoms"] = freeze_atoms

    # Per-segment post-processing results use backend-neutral machine keys;
    # concrete provenance lives in top-level mlip_backend/mlip_model fields.
    if post_segments:
        summary["post_segments"] = _json_safe(post_segments)

    # Key output file paths for AI agent consumption
    if "out_dir" in summary:
        # Real pipeline root. Fallback to the legacy module_dir.parent for any
        # caller that does not pass out_dir explicitly.
        root = Path(out_dir) if out_dir is not None else Path(summary["out_dir"]).parent
        current_paths = _current_output_paths(manifest, root) if manifest is not None else []
        key_files = _current_key_output_files(manifest, root) if manifest is not None else {}
        if current_paths:
            summary["current_output_paths"] = current_paths
        else:
            summary.pop("current_output_paths", None)
        if key_files:
            summary["key_output_files"] = key_files
        else:
            summary.pop("key_output_files", None)

    # Environment info
    try:
        from pdb2reaction.core.utils import _collect_environment_info
        summary.setdefault("environment", _collect_environment_info())
    except Exception:
        pass

    identified = apply_current_run_id(summary)
    summary.clear()
    summary.update(identified)
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


def _required_xyz_block_energies(
    blocks: Sequence[Sequence[str]],
    *,
    path: Path,
    context: str,
) -> List[float]:
    """Return one finite Hartree energy per XYZ block or fail the segment."""
    from pdb2reaction.io.xyz_trajectory import parse_xyz_energy_comment

    energies: List[float] = []
    for frame_index, block in enumerate(blocks, start=1):
        comment = block[1] if len(block) >= 2 else ""
        energy, provenance = parse_xyz_energy_comment(comment)
        if energy is None or not np.isfinite(energy):
            raise click.ClickException(
                f"[all] {context} trajectory {path} has no unambiguous finite "
                f"energy in frame {frame_index} ({provenance}); write "
                "E=<value> with an optional unit."
            )
        energies.append(float(energy))
    return energies


def _ensure_hei_path_tangent(
    mep_trj: Path,
    hei_path: Path,
    mode_path: Path,
) -> Optional[Path]:
    """Write the central MEP tangent at the image matching ``hei_path``.

    The normalized Cartesian difference is unit-independent.  This fallback
    also makes older path-search artifacts usable; new path-search runs emit
    the same reference-mode file directly when selecting the HEI.
    """
    if not mep_trj.exists() or not hei_path.exists():
        return None
    try:
        from ase.io import read as ase_read

        images = list(ase_read(str(mep_trj), index=":"))
        hei = ase_read(str(hei_path), index=0)
        if len(images) < 2:
            return None
        hei_numbers = np.asarray(hei.numbers)
        if any(
            len(image) != len(hei)
            or not np.array_equal(np.asarray(image.numbers), hei_numbers)
            for image in images
        ):
            return None
        hei_pos = np.asarray(hei.positions, dtype=float)
        distances = [
            float(np.linalg.norm(np.asarray(image.positions) - hei_pos))
            for image in images
        ]
        index = int(np.argmin(distances))
        energies = None
        raw_blocks = read_xyz_as_blocks(mep_trj)
        if len(raw_blocks) == len(images):
            parsed_energies = []
            for block in raw_blocks:
                from pdb2reaction.io.xyz_trajectory import parse_xyz_energy_comment

                comment = block[1] if len(block) >= 2 else ""
                energy, _provenance = parse_xyz_energy_comment(comment)
                parsed_energies.append(energy)
            if all(
                energy is not None and np.isfinite(energy)
                for energy in parsed_energies
            ):
                energies = parsed_energies
        tangent = _path_search._normalized_path_tangent(
            [np.asarray(image.positions, dtype=float) for image in images],
            index,
            energies=energies,
        )
        if tangent is None:
            return None
        mode_path.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(mode_path, tangent, fmt="%.17e")
        _echo(
            f"[tsopt] Derived HEI path-tangent reference mode → {mode_path}"
        )
        return mode_path
    except Exception as exc:
        _echo(
            f"[tsopt] WARNING: Could not derive HEI path tangent: {exc}",
            err=True,
        )
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
    image_written = False
    image_error: Optional[str] = None
    try:
        fig.write_image(str(png), scale=2)
        image_written = True
        _echo(f"[diagram] Wrote energy diagram → {png.name}")
    except Exception as e:
        image_error = str(e)
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
    payload["image"] = str(png) if image_written else None
    payload["image_written"] = image_written
    if image_error:
        payload["image_error"] = image_error
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
    calc_identity_cfg: Optional[Dict[str, Any]] = None,
    reject_uphill: Optional[bool] = None,
) -> Tuple[Any, Path, Optional[bool]]:
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
        ``(optimized_geometry, final_xyz_path, is_converged)``.  ``is_converged``
        is the fail-closed tri-state convergence bit of the final optimizer:
        an endpoint whose optimization did not explicitly converge is retained as
        a geometry/artifact but must not promote its segment to a usable success.
    """
    from pdb2reaction.workflows._outcomes import optimizer_converged_bit
    geom.set_calculator(getattr(geom, "calculator", None))
    mode = (opt_mode_default or "hess").lower()
    if mode in ("grad", "lbfgs", "dimer"):
        run_sequence = ("lbfgs",)
    elif mode in ("hess", "rfo", "rsprfo"):
        run_sequence = ("rfo",)
    else:
        run_sequence = ("rfo",)

    final_xyz: Optional[Path] = None
    _endpoint_conv: Optional[bool] = None
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
        # RFO-only endpoint re-optimization uphill-rejection toggle (min-scoped).
        # ``None`` leaves RFO_KW's own default in place (byte-identical behavior);
        # an explicit bool is threaded from the ``all`` command's
        # --reject-uphill/--no-reject-uphill flag.
        if sopt_kind == "rfo" and reject_uphill is not None:
            cfg["reject_uphill"] = bool(reject_uphill)

        # Seed cached IRC endpoint Hessian for RFO when available, but only on
        # a full evaluation-identity match (run/system/evaluator/active
        # space/potential).  Without a resolved evaluator config we cannot
        # prove compatibility, so a fresh Hessian is computed instead.
        if sopt_kind == "rfo":
            from pdb2reaction.io.hessian_cache import (
                load_matching as _hess_load_matching,
                identity_from_context as _hess_identity,
            )
            _cached = None
            if calc_identity_cfg is not None:
                _cached = _hess_load_matching(
                    "irc_endpoint",
                    _hess_identity(geom, calc_identity_cfg, role="irc_endpoint"),
                )
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
            _endpoint_conv = optimizer_converged_bit(opt)
        except (OptimizationError, ZeroStepLength) as e:
            _echo(
                f"[endpoint-opt] WARNING: optimization for '{tag}' terminated early ({e}); using last geometry.",
                err=True,
            )
            _endpoint_conv = False

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
    return g_final, final_xyz, _endpoint_conv


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
    _append_cli_arg(args, "--symmetry-number", overrides.get("symmetry_number"))
    _append_toggle_arg(args, "--dump", bool(dump_use))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        args.extend(["--hessian-calc-mode", str(hess_mode)])

    # Pass backend so freq uses the same MLIP (or the custom calc-file).
    _freq_custom = _forward_calc_file_argv(args, overrides)

    # Forward MLIP runtime knobs when the parent CLI explicitly set them; values are
    # injected into overrides by the cli() body before this helper is called.
    _append_explicit_child_runtime_argv(
        args,
        workers=overrides.get("workers"),
        workers_per_node=overrides.get("workers_per_node"),
        backend=None if _freq_custom else overrides.get("backend"),
        solvent=overrides.get("solvent"),
        solvent_model=overrides.get("solvent_model"),
    )

    if args_yaml is not None:
        args.extend(["--config", str(args_yaml)])
    _freq_rc = _run_cli_main("freq", _freq_cli.cli, args, on_nonzero="warn", prefix="freq")
    y = fdir / "thermoanalysis.yaml"
    # a nonzero freq exit means the thermochemistry is NOT usable, even if a
    # thermoanalysis.yaml (from a prior run or a partial write) exists with finite
    # fields. Never infer FREQ success from the filename or a finite number.
    if _freq_rc not in (None, 0):
        _echo(
            f"[freq] WARNING: freq exited with code {_freq_rc}; thermochemistry is "
            "unusable and will not enter any Gibbs diagram.",
            err=True,
        )
        return {}
    if not bool(dump_use):
        return {}
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
    """Extract DFT energy in hartree, or None if DFT failed or the value is not finite.

    Non-finite is reported as None at this single chokepoint so that every consumer's
    ``is not None`` check is sufficient; a NaN/inf must never reach a diagram, summary.json
    or logged state energies.
    """
    if not _dft_succeeded(result):
        return None
    try:
        value = float((result.get("energy") or {}).get("hartree"))
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _thermo_value_ha(payload: Any, key: str) -> Optional[float]:
    """Return one finite thermochemistry value in hartree."""

    if not isinstance(payload, dict):
        return None
    try:
        value = float(payload.get(key))
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


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
    # A fresh interpreter isolates the incompatible libcusolver versions linked
    # by torch (UMA/ORB) and gpu4pyscf.
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


def _validate_tsopt_result_payload(payload: Dict[str, Any]) -> None:
    """Reject a TS result that is not a verified first-order saddle."""
    status = str(payload.get("status") or "unknown")
    n_imag = payload.get("n_imaginary_modes")
    if status != "converged" or n_imag != 1:
        raise click.ClickException(
            "[tsopt] TS optimization did not produce a validated first-order "
            f"saddle (status={status!r}, n_imag={n_imag!r}); IRC was not started."
        )


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
    runtime_overrides: Optional[Dict[str, Any]] = None,
    manifest: Optional[InvocationManifest] = None,
    artifact_prefix: str = "tsopt",
    public_root: Optional[Path] = None,
) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI model structure; return (final_geom_path, ts_geom).

    Prefer the XYZ output to preserve coordinate precision between workflow steps, while still writing
    PDB/CIF/GJF companions when requested by the original input type.
    """
    overrides = overrides or {}
    manifest = manifest or InvocationManifest()
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

        reference_mode = overrides.get("reference_mode")
        if reference_mode is not None:
            ts_args.extend(["--ref-mode", str(reference_mode)])

        _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
        _append_toggle_arg(ts_args, "--dump", overrides.get("dump"))
        _append_cli_arg(ts_args, "--thresh", overrides.get("thresh"))
        _append_toggle_arg(ts_args, "--flatten", overrides.get("flatten"))

        hess_mode = overrides.get("hessian_calc_mode")
        if hess_mode:
            ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

        # Pass backend from calc_cfg so tsopt uses the same MLIP (or custom calc-file).
        _ts_custom = _forward_calc_file_argv(ts_args, calc_cfg)

        # Forward MLIP runtime knobs when CLI-explicit (calc_cfg only contains keys
        # supplied via --workers / --workers-per-node / --solvent / --solvent-model;
        # YAML-only values stay in --config args_yaml below). Without these explicit
        # CLI argv entries the subprocess silently falls back to defaults.
        _runtime = runtime_overrides
        if _runtime is None:
            # Backward-compatible internal-call fallback. The composite CLI
            # always supplies a source-aware mapping, including an empty one.
            _runtime = {
                "workers": calc_cfg.get("workers") if calc_cfg.get("workers") != 1 else None,
                "workers_per_node": (
                    calc_cfg.get("workers_per_node")
                    if calc_cfg.get("workers_per_node") != 1 else None
                ),
                "backend": calc_cfg.get("backend") if calc_cfg.get("backend") != "uma" else None,
                "solvent": calc_cfg.get("solvent") if calc_cfg.get("solvent") != "none" else None,
                "solvent_model": (
                    calc_cfg.get("solvent_model")
                    if calc_cfg.get("solvent_model") != "alpb" else None
                ),
            }
        _append_explicit_child_runtime_argv(
            ts_args,
            workers=_runtime.get("workers"),
            workers_per_node=_runtime.get("workers_per_node"),
            backend=None if _ts_custom else _runtime.get("backend"),
            solvent=_runtime.get("solvent"),
            solvent_model=_runtime.get("solvent_model"),
        )

        if args_yaml is not None:
            ts_args.extend(["--config", str(args_yaml)])
        if ref_pdb is not None:
            ts_args.extend(["--ref-pdb", str(ref_pdb)])
        ts_args.append("--out-json")

        _echo()
        _echo_detail(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
        result_key = f"{artifact_prefix}.result"
        geometry_key = f"{artifact_prefix}.geometry"
        result_destination = ts_dir / "result.json"
        geometry_destinations = [
            ts_dir / "final_geometry.xyz",
            ts_dir / "final_geometry.pdb",
            ts_dir / "final_geometry.gjf",
        ]
        manifest.declare(result_key, [result_destination])
        manifest.declare(
            geometry_key,
            geometry_destinations,
        )
        if public_root is not None:
            for destination in [result_destination, *geometry_destinations]:
                if _is_pipeline_public_destination(public_root, destination):
                    _declare_public_output(
                        manifest,
                        public_root,
                        destination,
                    )
        _run_cli_main("tsopt", _tsopt.cli, ts_args, on_nonzero="raise", prefix="tsopt")

        result_path = manifest.claim_one(result_key)
        if public_root is not None and _is_pipeline_public_destination(
            public_root, result_path
        ):
            _claim_public_output(
                manifest,
                public_root,
                result_path,
            )
        try:
            tsopt_result = json.loads(result_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise click.ClickException(
                f"[tsopt] Could not read TS validation result '{result_path}': {exc}"
            ) from exc
        _validate_tsopt_result_payload(tsopt_result)

        ts_pdb = ts_dir / "final_geometry.pdb"
        ts_gjf = ts_dir / "final_geometry.gjf"
        ts_geom_path = manifest.claim_one(geometry_key)

        if ts_geom_path.suffix.lower() == ".xyz":
            try:
                convert_xyz_like_outputs(
                    ts_geom_path,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=ts_pdb if needs_pdb else None,
                    out_gjf_path=ts_gjf if needs_gjf else None,
                )
            except Exception as e:
                _echo(f"[tsopt] WARNING: Failed to convert TS geometry: {e}", err=True)

        if public_root is not None:
            for destination in geometry_destinations:
                if _is_pipeline_public_destination(public_root, destination):
                    _claim_public_output(
                        manifest,
                        public_root,
                        destination,
                    )

        g_ts = geom_loader(
            ts_geom_path,
            coord_type=DEFAULT_COORD_TYPE,
            freeze_atoms=_freeze_atoms_for_log(),
        )
        g_ts._tsopt_result = tsopt_result

        return ts_geom_path, g_ts
    finally:
        prepared_input.cleanup()


def _orient_irc_endpoints(
    g_left: Any,
    g_right: Any,
    *,
    endpoint_trajectory: Optional[Path],
    freeze_atoms: Sequence[int],
    seg_tag: Optional[str],
) -> Tuple[Any, Any, str, str, bool, Dict[str, Any]]:
    """Orient IRC endpoints against one already-claimed MEP trajectory."""

    left_tag = "forward"
    right_tag = "backward"
    reverse_irc = False
    assignment: Dict[str, Any] = {
        "method": "raw",
        "reversed": False,
        "connectivity_validated": None if seg_tag is None else False,
    }
    if seg_tag is None:
        _echo_detail("[irc] TSOPT-only mode: Use raw irc orientation.")
        return (
            g_left,
            g_right,
            left_tag,
            right_tag,
            reverse_irc,
            assignment,
        )
    if endpoint_trajectory is None:
        _echo(
            f"[irc] WARNING: current MEP endpoints were not claimed for segment tag '{seg_tag}'; "
            "using raw IRC orientation.",
            err=True,
        )
        assignment["reason"] = "mep_endpoint_trajectory_missing"
        return (
            g_left,
            g_right,
            left_tag,
            right_tag,
            reverse_irc,
            assignment,
        )

    try:
        elems, c_first, c_last = read_xyz_first_last(endpoint_trajectory)
        gL_end = _geom_from_angstrom(elems, c_first, freeze_atoms)
        gR_end = _geom_from_angstrom(elems, c_last, freeze_atoms)
        bond_cfg = dict(_path_search.BOND_KW)

        def _matches(x, y) -> bool:
            try:
                changed, _ = _path_search.has_bond_change(x, y, bond_cfg)
                return not changed
            except Exception as exc:
                logger.debug(
                    "Bond-change check failed during IRC endpoint matching: %s", exc
                )
                return False

        L_L = _matches(g_left, gL_end)
        L_R = _matches(g_left, gR_end)
        R_L = _matches(g_right, gL_end)
        R_R = _matches(g_right, gR_end)
        direct_match = bool(L_L and R_R)
        swapped_match = bool(L_R and R_L)
        assignment["match_matrix"] = {
            "left_to_mep_left": bool(L_L),
            "left_to_mep_right": bool(L_R),
            "right_to_mep_left": bool(R_L),
            "right_to_mep_right": bool(R_R),
        }

        if direct_match ^ swapped_match:
            assignment["method"] = "bond_topology"
            assignment["connectivity_validated"] = True
            if swapped_match:
                g_left, g_right = g_right, g_left
                left_tag, right_tag = right_tag, left_tag
                reverse_irc = True
        else:
            assignment["method"] = (
                "rmsd_topology_tie"
                if direct_match and swapped_match
                else "rmsd_topology_unmatched"
            )
            assignment["connectivity_validated"] = bool(
                direct_match and swapped_match
            )
            try:
                d_LL = _path_search._rmsd_between(g_left, gL_end)
                d_LR = _path_search._rmsd_between(g_left, gR_end)
                d_RL = _path_search._rmsd_between(g_right, gL_end)
                d_RR = _path_search._rmsd_between(g_right, gR_end)
                direct_score = float(d_LL + d_RR)
                swapped_score = float(d_LR + d_RL)
                assignment["rmsd_direct"] = direct_score
                assignment["rmsd_swapped"] = swapped_score
                if swapped_score < direct_score:
                    g_left, g_right = g_right, g_left
                    left_tag, right_tag = right_tag, left_tag
                    reverse_irc = True
            except Exception as exc:
                assignment["method"] = "unresolved"
                assignment["connectivity_validated"] = False
                assignment["reason"] = f"rmsd_failed:{exc}"
                _echo(
                    f"[irc] WARNING: segment endpoint mapping via RMSD failed: {exc}",
                    err=True,
                )

        assignment["reversed"] = bool(reverse_irc)
        if not assignment["connectivity_validated"]:
            assignment.setdefault(
                "reason",
                "neither endpoint pairing matched the MEP bond topology",
            )
    except Exception as exc:
        assignment["method"] = "unresolved"
        assignment["connectivity_validated"] = False
        assignment["reason"] = f"endpoint_mapping_failed:{exc}"
        _echo(f"[irc] WARNING: segment endpoint mapping failed: {exc}", err=True)
    return (
        g_left,
        g_right,
        left_tag,
        right_tag,
        reverse_irc,
        assignment,
    )


def _read_irc_outcome(irc_dir: Path) -> Dict[str, Any]:
    """Read the IRC child's ``result.json`` into a fail-closed usability record.

    The IRC leaf is *usable* only when the child reports ``scientific_status ==
    "success"`` — i.e. every requested direction explicitly converged. A
    missing / unreadable result, or any nonconverged requested direction, yields
    ``usable=False`` while the endpoint trajectory remains a reportable artifact.
    """

    outcome: Dict[str, Any] = {
        "usable": False,
        "reason": "irc_result_missing",
        "scientific_status": None,
        "forward_converged": None,
        "backward_converged": None,
        "n_frames_forward": None,
        "n_frames_backward": None,
        "traj": None,
    }
    result_path = irc_dir / "result.json"
    if not result_path.exists():
        return outcome
    try:
        data = json.loads(result_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        outcome["reason"] = "irc_result_unreadable"
        return outcome
    if not isinstance(data, dict):
        outcome["reason"] = "irc_result_unreadable"
        return outcome

    sci = data.get("scientific_status")
    outcome["scientific_status"] = sci
    outcome["forward_converged"] = data.get("forward_converged")
    outcome["backward_converged"] = data.get("backward_converged")
    outcome["n_frames_forward"] = data.get("n_frames_forward")
    outcome["n_frames_backward"] = data.get("n_frames_backward")
    _files = data.get("files") if isinstance(data.get("files"), dict) else {}
    outcome["traj"] = _files.get("finished_irc")

    if sci == "success":
        outcome["usable"] = True
        outcome["reason"] = "ok"
    elif isinstance(sci, str):
        outcome["usable"] = False
        reasons = data.get("scientific_status_reasons")
        outcome["reason"] = (
            ";".join(str(r) for r in reasons)
            if isinstance(reasons, list) and reasons
            else f"irc_{sci}"
        )
    else:
        # No explicit status field: fail closed rather than trust file existence.
        outcome["usable"] = False
        outcome["reason"] = "irc_status_unknown"
    return outcome


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
    irc_step_size: Optional[float] = None,
    irc_never_stop: Optional[bool] = None,
    seg_tag: Optional[str] = None,
    endpoint_trajectory: Optional[Path] = None,
    session: Optional[RunSession] = None,
    manifest: Optional[InvocationManifest] = None,
    artifact_prefix: str = "irc",
    public_root: Optional[Path] = None,
    runtime_overrides: Optional[Dict[str, Any]] = None,
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
    manifest = manifest or InvocationManifest()
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
    _irc_custom = _forward_calc_file_argv(irc_args, calc_cfg)

    # Forward MLIP runtime knobs from calc_cfg when CLI-explicit (only set in the
    # shared calc_cfg when user passed them; otherwise downstream irc CLI would
    # fall back to its own defaults regardless of parent --workers / --solvent etc.).
    _runtime = runtime_overrides
    if _runtime is None:
        _runtime = {
            "workers": calc_cfg.get("workers") if calc_cfg.get("workers") != 1 else None,
            "workers_per_node": (
                calc_cfg.get("workers_per_node")
                if calc_cfg.get("workers_per_node") != 1 else None
            ),
            "backend": calc_cfg.get("backend") if calc_cfg.get("backend") != "uma" else None,
            "solvent": calc_cfg.get("solvent") if calc_cfg.get("solvent") != "none" else None,
            "solvent_model": (
                calc_cfg.get("solvent_model")
                if calc_cfg.get("solvent_model") != "alpb" else None
            ),
        }
    _append_explicit_child_runtime_argv(
        irc_args,
        workers=_runtime.get("workers"),
        workers_per_node=_runtime.get("workers_per_node"),
        backend=None if _irc_custom else _runtime.get("backend"),
        solvent=_runtime.get("solvent"),
        solvent_model=_runtime.get("solvent_model"),
    )

    if args_yaml is not None:
        irc_args.extend(["--config", str(args_yaml)])
    if irc_step_size is not None:
        irc_args.extend(["--step-size", str(float(irc_step_size))])
    if irc_never_stop is not None:
        _append_toggle_arg(irc_args, "--never-stop", bool(irc_never_stop))
    # request the child's machine-readable result.json so the aggregate can
    # gate on reported per-direction IRC convergence instead of trajectory-file
    # existence. A never_stop / max-cycle direction still writes its trajectory,
    # but the child reports it as nonconverged and we must not promote it.
    _append_toggle_arg(irc_args, "--out-json", True)
    _echo()
    _echo_detail(f"[irc] Running EulerPC IRC → out={irc_dir}")
    trajectory_key = f"{artifact_prefix}.trajectory"
    plot_key = f"{artifact_prefix}.plot"
    pdb_key = f"{artifact_prefix}.pdb"
    trajectory_destination = irc_dir / "finished_irc_trj.xyz"
    plot_destination = irc_dir / "irc_plot.png"
    pdb_destination = irc_dir / "finished_irc.pdb"
    manifest.declare(trajectory_key, [trajectory_destination])
    manifest.declare(plot_key, [plot_destination])
    manifest.declare(pdb_key, [pdb_destination])
    if public_root is not None:
        for destination in (
            trajectory_destination,
            plot_destination,
            pdb_destination,
        ):
            if _is_pipeline_public_destination(public_root, destination):
                _declare_public_output(
                    manifest,
                    public_root,
                    destination,
                )
    _run_cli_main("irc", _irc_cli.cli, irc_args, on_nonzero="raise", prefix="irc")

    # read the child's per-direction convergence. The IRC leaf is
    # usable only when EVERY requested direction explicitly converged; a
    # trajectory can exist for a nonconverged (never_stop / max-cycle) direction,
    # so promotion must gate on this outcome, not on file existence.
    irc_outcome = _read_irc_outcome(irc_dir)

    finished_pdb = irc_dir / "finished_irc.pdb"
    finished_trj = manifest.claim_one(trajectory_key)
    if public_root is not None and _is_pipeline_public_destination(
        public_root, finished_trj
    ):
        _claim_public_output(
            manifest,
            public_root,
            finished_trj,
        )
    irc_plot = irc_dir / "irc_plot.png"

    # Ensure we have a PDB for visualization if possible
    try:
        if manifest.claim_optional(pdb_key) is None:
            ref_for_conv: Optional[Path] = None
            if seg_model_pdb.suffix.lower() == ".pdb":
                ref_for_conv = seg_model_pdb
            elif ref_pdb_for_seg.suffix.lower() == ".pdb":
                ref_for_conv = ref_pdb_for_seg
            if ref_for_conv is not None:
                _path_search._convert_to_pdb_logged(finished_trj, ref_pdb_path=ref_for_conv, out_path=finished_pdb)
    except Exception as e:
        _echo(f"[irc] WARNING: failed to convert finished_irc_trj.xyz to PDB: {e}", err=True)

    manifest.claim_optional(pdb_key)
    if public_root is not None and _is_pipeline_public_destination(
        public_root, pdb_destination
    ):
        _claim_public_output(
            manifest,
            public_root,
            pdb_destination,
        )
    elems, c_first, c_last = read_xyz_first_last(finished_trj)
    g_left = _geom_from_angstrom(elems, c_first, freeze_atoms)
    g_right = _geom_from_angstrom(elems, c_last, freeze_atoms)
    shared_calc = create_calculator(**calc_cfg)
    lease = CalculatorLease(shared_calc)
    if session is not None:
        session.resources.add(lease.release)
    try:
        lease.attach(g_ts)
        lease.attach(g_left)
        lease.attach(g_right)
        (
            g_left,
            g_right,
            left_tag,
            right_tag,
            reverse_irc,
            endpoint_assignment,
        ) = _orient_irc_endpoints(
            g_left,
            g_right,
            endpoint_trajectory=endpoint_trajectory,
            freeze_atoms=freeze_atoms,
            seg_tag=seg_tag,
        )

        try:
            run_trj2fig(
                finished_trj,
                [irc_plot],
                unit="kcal",
                reference="init",
                reverse_x=False,
            )
            close_matplotlib_figures()
        except Exception as e:
            _echo(f"[irc] WARNING: failed to plot finished IRC trajectory: {e}", err=True)
        current_plot = manifest.claim_optional(plot_key)
        if public_root is not None and _is_pipeline_public_destination(
            public_root, plot_destination
        ):
            _claim_public_output(
                manifest,
                public_root,
                plot_destination,
            )

        return {
            "left_min_geom": g_left,
            "right_min_geom": g_right,
            "ts_geom": g_ts,
            "left_tag": left_tag,
            "right_tag": right_tag,
            "freeze_atoms": freeze_atoms,
            "irc_plot_path": current_plot,
            "irc_trj_path": finished_trj,
            "reverse_irc": reverse_irc,
            "endpoint_assignment": endpoint_assignment,
            "calculator_lease": lease,
            "irc_outcome": irc_outcome,
        }
    except BaseException:
        lease.release()
        raise



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


@click.command(
    help=(
        "Run active site model extraction → optional staged scan → MEP search → full-structure merge in one run.\n"
        "If exactly one input is provided: (a) with --scan-lists, run staged scan on the active site model (or full structure "
        "when extraction is skipped) and use stage results as inputs for path-opt (path_search with --refine-path); "
        "(b) with --tsopt and no --scan-lists, run TSOPT-only mode."
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
    callback=_show_advanced_subcommand_help,
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
        "Two or more **full structures** (PDB/mmCIF/XYZ/GJF) in reaction order (reactant [intermediates ...] product), "
        "or a single **full structure** (with --scan-lists or with --tsopt). "
        "Extraction (-c/--center) accepts PDB/mmCIF. mmCIF is processed through an internal PDB "
        "bridge and emitted again as mmCIF; oversized PDBs use the same safe bridge. "
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
        "a PDB/mmCIF path, a residue-ID list like '123,124' or 'A:123,B:456' "
        "(insertion codes OK: '123A' / 'A:123A'), "
        "a residue-name list like 'GPP,SAM', or chain-qualified 'A:SAM' / 'A:SAM:123'. "
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
        "Force-include residues using the same selectors as -c/--center: IDs ('123', "
        "'A:123A'), names ('SAM'), or chain-qualified names ('A:SAM', 'A:SAM:123'); "
        "comma/space separated."
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
        "(extraction skipped) the same mapping is applied to the full input PDB/mmCIF to "
        "derive the total system charge. A bare number sets the total directly. "
        "PDB/mmCIF inputs only. Without extraction, -q/--charge explicitly sets "
        "the total; with extraction, it asserts the derived total."
    ),
)
@click.option(
    "-q",
    "--charge",
    "charge_override",
    type=int,
    default=None,
    help=(
        "Total system charge. With -c/--center, this is an assertion and must "
        "match the extractor-derived charge; omit it to auto-derive. Without "
        "extraction it explicitly sets/overrides the total and emits a warning."
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
    help="Freeze parent atoms of cap hydrogens (PDB/mmCIF input or XYZ/GJF with --ref-pdb).",
)
@click.option(
    "--tr-projection",
    type=click.Choice(["constrained", "legacy-active"], case_sensitive=False),
    default=GEOM_KW_DEFAULT["tr_projection"],
    show_default=True,
    help=(
        "Rigid translation/rotation treatment forwarded to TSopt, IRC, freq, "
        "and flatten PHVA. The default respects frozen anchors; 'legacy-active' "
        "is deprecated and must not be used for pass/HOSP transition-state certification."
    ),
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
    default=_path_search.GS_KW["max_nodes"],
    show_default=True,
    help=("Movable internal images per GSM/DMF segment; the complete segment has "
          "max_nodes+2 images including endpoints."),
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
        "When --thermo is enabled, freq always retains thermoanalysis.yaml because "
        "the composite workflow consumes that file; --no-dump does not suppress it."
    ),
)
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/CIF/GJF companions based on the input format.",
)
@click.option(
    "--refine-path",
    "refine_path",
    type=click.BOOL,
    default=False,
    show_default=True,
    help=(
        "Run a single-pass path-opt GSM between each adjacent pair and concatenate the "
        "segments (default; no path_search). Use --refine-path to run recursive "
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
        "With --thermo, also generate a DFT//MLIP Gibbs diagram."
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
    "--reject-uphill/--no-reject-uphill",
    "reject_uphill",
    default=True,
    show_default=True,
    help=(
        "Reject uphill RFO trials during post-IRC endpoint re-optimization only "
        "and final-check the retained endpoint at the emergency floor. Does not "
        "affect TS optimization or path search."
    ),
)
@click.option(
    "--irc-step-size",
    type=float,
    default=None,
    help=(
        "Override IRC --step-size (Bohr). If an IRC stops after only a few "
        "frames, retry with a smaller value such as 0.05."
    ),
)
@click.option(
    "--irc-never-stop/--no-irc-never-stop",
    default=None,
    help=(
        "Override IRC energy-stop handling. When enabled, transient energy "
        "rises/plateaus do not stop tracing; physical convergence and the "
        "IRC max-cycle limit still apply. Default off."
    ),
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
    "--freq-symmetry-number",
    type=click.IntRange(min=1),
    default=None,
    help=(
        "Use one rotational symmetry number for every R/TS/P frequency job. "
        "When omitted, each child follows its YAML/default setting."
    ),
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
        "Scan targets: inline Python literal. "
        "Multiple inline literals define sequential stages, e.g. "
        "'[(12,45,1.35)]' '[(10,55,2.20),(23,34,1.80)]'. "
        "Indices refer to the original full input ordering (1-based); atom strings may use "
        "CHAIN:RESNAME:RESSEQ[ICODE]:ATOM. When extraction is used, selections are auto-mapped "
        "to the active site model after extraction."
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
        "Reference PDB/mmCIF for topology when -i provides XYZ inputs. "
        "Enables topology-aware PDB/mmCIF output conversion in TSOPT-only, scan, and path_search pipelines."
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
    tr_projection: str,
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
    reject_uphill: bool,
    irc_step_size: Optional[float],
    irc_never_stop: Optional[bool],
    freq_out_dir: Optional[Path],
    freq_max_write: Optional[int],
    freq_amplitude_ang: Optional[float],
    freq_n_frames: Optional[int],
    freq_sort: Optional[str],
    freq_temperature: Optional[float],
    freq_pressure: Optional[float],
    freq_symmetry_number: Optional[int],
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
    calc_factory: Optional[str],
) -> None:
    """
    The **all** command composes `extract` → (optional `scan` on model or full input) → MEP search
    (single-pass `path-opt` by default, or recursive `path_search` with ``--refine-path``) and hides
    ref-template bookkeeping.
    With single input:
      - with --scan-lists: run staged scan and use stage results as inputs for path-opt (or path_search),
      - with --tsopt and no --scan-lists: run TSOPT-only mode (no MEP search).

    By default, a single-pass ``path-opt`` GSM is run between each adjacent pair of inputs and the segments
    are concatenated into the final MEP.  With ``--refine-path``, the recursive ``path_search`` workflow
    is used instead for automatic multistep discovery.
    """
    from pdb2reaction.core.utils import reject_option_like_extra_args

    # These negative spellings have long been documented for the legacy
    # value-style booleans.  Apply them explicitly before rejecting every
    # other option-like residual token.
    _negative_bool_aliases = {
        "--no-include-h2o",
        "--no-exclude-backbone",
        "--no-add-linkh",
        "--no-freeze-links",
        "--no-climb",
        "--no-dump",
        "--no-convert-files",
        "--no-refine-path",
        "--no-preopt",
        "--no-tsopt",
        "--no-thermo",
        "--no-dft",
    }
    from pdb2reaction.core.utils import current_cli_args

    _argv = current_cli_args(ctx)
    reject_option_like_extra_args(
        ctx.args,
        allowed_options=_negative_bool_aliases,
        allowed_values=collect_option_values(
            _argv, ("-i", "--input", "-s", "--scan-lists")
        ),
        consumed_values=[*input_paths, *scan_lists_raw],
    )
    _negative_bool_args = frozenset(str(value) for value in ctx.args)
    if "--no-include-h2o" in _negative_bool_args:
        include_h2o = False
    if "--no-exclude-backbone" in _negative_bool_args:
        exclude_backbone = False
    if "--no-add-linkh" in _negative_bool_args:
        add_linkh = False
    if "--no-freeze-links" in _negative_bool_args:
        freeze_links_flag = False
    if "--no-climb" in _negative_bool_args:
        climb = False
    if "--no-dump" in _negative_bool_args:
        dump = False
    if "--no-convert-files" in _negative_bool_args:
        convert_files = False
    if "--no-refine-path" in _negative_bool_args:
        refine_path = False
    if "--no-preopt" in _negative_bool_args:
        preopt = False
    if "--no-tsopt" in _negative_bool_args:
        do_tsopt = False
    if "--no-thermo" in _negative_bool_args:
        do_thermo = False
    if "--no-dft" in _negative_bool_args:
        do_dft = False

    global _FREEZE_ATOMS_GLOBAL, _FREEZE_ATOMS_YAML
    from pdb2reaction.core.utils import (
        is_child_mode,
        is_convert_file_enabled,
        pipeline_mode_enabled,
        set_child_mode,
        set_pipeline_mode,
    )

    # Register one invocation owner before the first process-global mutation.
    session = RunSession()
    ctx.call_on_close(session.close)
    session.own_run_id_environment()
    ctx.meta["pdb2reaction_run_session"] = session
    prior_pipeline_mode = pipeline_mode_enabled()
    prior_child_mode = is_child_mode()
    prior_convert_files = is_convert_file_enabled()
    prior_echo_started = bool(_echo_state._started)
    prior_freeze_global = (
        None if _FREEZE_ATOMS_GLOBAL is None else list(_FREEZE_ATOMS_GLOBAL)
    )
    prior_freeze_yaml = (
        None if _FREEZE_ATOMS_YAML is None else list(_FREEZE_ATOMS_YAML)
    )
    prior_sigint = signal.getsignal(signal.SIGINT)

    def _restore_invocation_state() -> None:
        global _FREEZE_ATOMS_GLOBAL, _FREEZE_ATOMS_YAML
        signal.signal(signal.SIGINT, prior_sigint)
        _FREEZE_ATOMS_GLOBAL = (
            None if prior_freeze_global is None else list(prior_freeze_global)
        )
        _FREEZE_ATOMS_YAML = (
            None if prior_freeze_yaml is None else list(prior_freeze_yaml)
        )
        _echo_state._started = prior_echo_started
        set_convert_file_enabled(prior_convert_files)
        set_child_mode(prior_child_mode)
        set_pipeline_mode(prior_pipeline_mode)

    session.resources.add(_restore_invocation_state)

    # Engage pipeline-scoped default-verbosity suppression for this invocation.
    set_pipeline_mode(True)
    argv_all = _argv

    def _sigint_handler(signum, frame):
        _echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    signal.signal(signal.SIGINT, _sigint_handler)

    _echo_state.reset()
    _FREEZE_ATOMS_GLOBAL = None
    _FREEZE_ATOMS_YAML = None
    set_convert_file_enabled(convert_files)
    command_str = "pdb2reaction " + " ".join(argv_all)
    time_start = time.perf_counter()
    energy_diagrams: List[Dict[str, Any]] = []

    dump_override_requested = False
    flatten_override_requested = False
    try:
        dump_source = ctx.get_parameter_source("dump")
        dump_override_requested = dump_source not in (None, ParameterSource.DEFAULT)
        flatten_source = ctx.get_parameter_source("flatten")
        flatten_override_requested = flatten_source not in (None, ParameterSource.DEFAULT)
    except Exception as exc:
        logger.debug("Failed to check dump/flatten parameter source: %s", exc)
        dump_override_requested = False
        flatten_override_requested = False

    args_yaml = config_yaml
    merged_yaml_cfg = load_yaml_dict(config_yaml) if config_yaml is not None else {}
    _dmf_yaml_cfg = (
        merged_yaml_cfg.get("dmf", {})
        if isinstance(merged_yaml_cfg, dict)
        else {}
    )
    dmf_correlated_effective = bool(
        fresh_dmf_config(_dmf_yaml_cfg).get("correlated", False)
    )

    # `calc.charge` / `calc.spin` are workflow inputs below explicit CLI
    # values and --ligand-charge, but above GJF metadata and defaults.
    calc_yaml_cfg = merged_yaml_cfg.get("calc") if isinstance(merged_yaml_cfg, dict) else None
    charge_override_label = "-q/--charge"
    if (
        charge_override is None
        and ligand_charge is None
        and isinstance(calc_yaml_cfg, dict)
        and calc_yaml_cfg.get("charge") is not None
    ):
        try:
            charge_override = int(calc_yaml_cfg["charge"])
        except (TypeError, ValueError) as exc:
            raise click.BadParameter(
                f"calc.charge must be an integer, got {calc_yaml_cfg['charge']!r}."
            ) from exc
        charge_override_label = "YAML calc.charge"
    charge_override_remove_hint = (
        "remove calc.charge from the YAML"
        if charge_override_label == "YAML calc.charge"
        else "omit -q/--charge"
    )

    spin_cli_explicit = cli_param_overridden(ctx, "spin")
    # Post-IRC endpoint re-optimization uphill-rejection toggle. ``None`` unless
    # the flag was explicitly passed, so the default path keeps RFO_KW's own
    # reject_uphill (unchanged behavior); an explicit --reject-uphill/--no-...
    # is threaded into _optimize_endpoint_geom (endpoint re-opt) only.
    _reject_uphill_eff = (
        bool(reject_uphill) if cli_param_overridden(ctx, "reject_uphill") else None
    )
    spin_configured = False
    if (
        not spin_cli_explicit
        and isinstance(calc_yaml_cfg, dict)
        and calc_yaml_cfg.get("spin") is not None
    ):
        try:
            spin = int(calc_yaml_cfg["spin"])
        except (TypeError, ValueError) as exc:
            raise click.BadParameter(
                f"calc.spin must be an integer multiplicity, got {calc_yaml_cfg['spin']!r}."
            ) from exc
        if spin < 1:
            raise click.BadParameter(
                f"calc.spin must be an integer multiplicity >= 1, got {spin}."
            )
        spin_configured = True

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

    scan_vals = collect_single_option_values(
        argv_all,
        ("-s", "--scan-lists"),
        "--scan-lists",
    )
    if scan_vals:
        scan_lists_raw = tuple(scan_vals)

    is_single = len(input_paths) == 1
    has_scan = bool(scan_lists_raw)
    if has_scan and not is_single:
        raise click.BadParameter(
            "--scan-lists requires exactly one input structure; provide one "
            "reactant for a staged scan or omit --scan-lists for an endpoint path."
        )
    single_tsopt_mode = is_single and (not has_scan) and do_tsopt

    if (len(input_paths) < 2) and not (is_single and (has_scan or do_tsopt)):
        raise click.BadParameter(
            "Provide at least two structures with -i/--input in reaction order, "
            "or use a single structure with --scan-lists, or a single structure with --tsopt."
        )

    # Normalize mmCIF and PDBs beyond fixed-column limits once, before any
    # extract/dry-run preflight. The prepared objects remain alive for this
    # composite workflow, while every child stage sees an ordinary PDB path.
    user_input_paths = tuple(Path(path) for path in input_paths)
    prepared_all_inputs = []
    for path in user_input_paths:
        prepared = prepare_input_structure(path)
        prepared_all_inputs.append(prepared)
        session.resources.own_cleanup(prepared)
    input_paths = tuple(prepared.source_path for prepared in prepared_all_inputs)
    for user_path, prepared in zip(user_input_paths, prepared_all_inputs):
        if prepared.source_path != user_path:
            _echo(
                f"[all] Structure bridge: {user_path} → internal PDB "
                f"({prepared.structure_template.reason if prepared.structure_template else 'normalized'})",
                narrative=True,
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

    _validate_postprocessing_dependencies(
        do_tsopt=do_tsopt,
        do_thermo=do_thermo,
        do_dft=do_dft,
    )

    tsopt_opt_mode_default: Optional[str] = None
    if opt_mode_post_set and opt_mode_post_norm is not None:
        tsopt_opt_mode_default = opt_mode_post_norm
    elif opt_mode_set:
        tsopt_opt_mode_default = opt_mode_norm
    else:
        tsopt_opt_mode_default = "hess"

    citation_post_segments: List[Dict[str, Any]] = []

    def _all_method_citation_payload() -> Dict[str, Any]:
        return {
            "pipeline_mode": all_mode,
            "path_opt_mode": opt_mode_norm,
            "post_opt_mode": tsopt_opt_mode_default,
            "ts_opt_mode": tsopt_opt_mode_default,
            "endpoint_opt_mode": tsopt_opt_mode_default,
            "mep_mode": mep_mode,
            "dmf_correlated": dmf_correlated_effective,
            "post_segments": citation_post_segments,
        }

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
    if flatten_override_requested:
        tsopt_overrides["flatten"] = bool(flatten)

    freq_overrides: Dict[str, Any] = {}
    # backend will be injected after calc_cfg_shared is built (see below)
    from pdb2reaction.workflows.freq import _validated_thermo_condition

    if freq_max_write is not None:
        freq_overrides["max_write"] = int(freq_max_write)
    if freq_amplitude_ang is not None:
        freq_overrides["amplitude_ang"] = float(freq_amplitude_ang)
    if freq_n_frames is not None:
        freq_overrides["n_frames"] = int(freq_n_frames)
    if freq_sort is not None:
        freq_overrides["sort"] = freq_sort.lower()
    if freq_temperature is not None:
        freq_overrides["temperature"] = _validated_thermo_condition(
            freq_temperature, name="temperature"
        )
    if freq_pressure is not None:
        freq_overrides["pressure"] = _validated_thermo_condition(
            freq_pressure, name="pressure_atm"
        )
    if freq_symmetry_number is not None:
        freq_overrides["symmetry_number"] = int(freq_symmetry_number)
    # all.py reads thermochemistry EXCLUSIVELY from freq's thermoanalysis.yaml
    # (its in-memory return is discarded via _run_cli_main). Always let freq
    # write that yaml so `--no-dump` cannot silently zero out tR/tT/tP and
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
                "inputs": [str(p) for p in user_input_paths],
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
                "tr_projection": str(tr_projection),
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
        emit(
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
            session.resources.own_path(_dry_tmp)
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
                        f"[all] {charge_override_label} {int(charge_override):+d} does not match extract-derived "
                        f"charge {int(_q_extracted):+d}. Either {charge_override_remove_hint} to auto-derive, "
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
        extract_note = (
            "extract pre-check ran" if not _skip_extract_dry
            else "extraction was not requested"
        )
        _echo(
            f"[all] Dry-run mode: no search/post-processing was executed ({extract_note}).",
            narrative=True,
        )
        planned_stages: List[str] = []
        if not _skip_extract_dry:
            planned_stages.append("extract")
        if has_scan:
            planned_stages.append("scan")
        if not single_tsopt_mode:
            planned_stages.append("path_search" if refine_path else "path_opt")
        if do_tsopt:
            planned_stages.extend(["tsopt", "irc"])
        if do_thermo:
            planned_stages.append("freq")
        if do_dft:
            planned_stages.append("dft")
        _echo(
            "[all] Planned stages: " + (
                " -> ".join(planned_stages) if planned_stages else "validation only"
            ) + ".",
            narrative=True,
        )
        _emit_final_summary(out_dir, time_start, session.manifest)
        return

    yaml_cfg = load_yaml_dict(args_yaml)
    _set_yaml_freeze_atoms(yaml_cfg)

    skip_extract = center_spec is None or str(center_spec).strip() == ""
    first_input = input_paths[0].resolve() if input_paths else None

    out_dir = out_dir.resolve()
    # Public ownership is established at each producer immediately before its
    # exact destination is written.  No root/segments traversal participates
    # in provenance admission.
    manifest = session.manifest

    def _declare_public(path: Path) -> str:
        return _declare_public_output(
            manifest,
            out_dir,
            path,
        )

    def _claim_public(path: Path) -> Optional[Path]:
        return _claim_public_output(
            manifest,
            out_dir,
            path,
        )

    def _copy_public_logged(
        src: Path,
        dst: Path,
        *,
        label: Optional[str] = None,
        echo: bool = True,
    ) -> bool:
        _declare_public(dst)
        copied = _copy_logged(src, dst, label=label, echo=echo)
        if copied:
            _claim_public(dst)
        return copied

    def _move_public_logged(
        src: Path,
        dst: Path,
        *,
        label: Optional[str] = None,
        echo: bool = True,
    ) -> bool:
        _declare_public(dst)
        moved = _move_logged(src, dst, label=label, echo=echo)
        if moved:
            _claim_public(dst)
        return moved

    def _write_public_energy_diagram(
        prefix: Path,
        labels: List[str],
        energies_au: List[float],
        title_note: str,
        ylabel: str = "ΔE (kcal/mol)",
    ) -> Optional[Dict[str, Any]]:
        destination = prefix.with_suffix(".png")
        _declare_public(destination)
        payload = _write_segment_energy_diagram(
            prefix,
            labels=labels,
            energies_au=energies_au,
            title_note=title_note,
            ylabel=ylabel,
        )
        _claim_public(destination)
        return payload

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
        for model_index, model_output in enumerate(model_outputs, start=1):
            manifest.declare(
                f"extract.model.{model_index:02d}",
                [model_output],
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

        for model_index in range(1, len(model_outputs) + 1):
            manifest.claim_one(f"extract.model.{model_index:02d}")
        _persist_run_manifest(manifest, out_dir)

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

        user_provided_spin = spin_cli_explicit or spin_configured

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
                except click.BadParameter:
                    # A rejected input (e.g. multi-MODEL) must abort, not fall through to the
                    # numeric-total fallback below.
                    raise
                except Exception as e:
                    _echo(
                        f"[all] NOTE: failed to derive total charge from full complex: {e}; "
                        "a numeric --ligand-charge will be treated as the total charge; "
                        "a residue mapping still requires successful metadata resolution.",
                        err=True,
                    )
            else:
                _echo(
                    "[all] NOTE: --ligand-charge derivation requires PDB/mmCIF residue metadata; skipping full-complex derivation.",
                    err=True,
                )

        if resolved_charge is None:
            if ligand_charge_numeric is not None:
                _echo(
                    f"[all] Using --ligand-charge as TOTAL system charge: {ligand_charge_numeric:+g}"
                )
                resolved_charge = _round_charge_with_note(ligand_charge_numeric, prefix="[all]")

        if any(
            prepared.gjf_template is not None
            for prepared in prepared_all_inputs
        ) and (resolved_charge is None or not user_provided_spin):
            parsed_charge, parsed_spin = fill_charge_spin_from_gjf_inputs(
                resolved_charge,
                spin if user_provided_spin else None,
                [prepared.gjf_template for prepared in prepared_all_inputs],
            )
            if resolved_charge is None:
                gjf_charge = parsed_charge
            if not user_provided_spin:
                gjf_spin = parsed_spin
            _echo(
                "[all] Parsed consistent GJF metadata: "
                f"charge={gjf_charge if gjf_charge is not None else resolved_charge:+d}, "
                f"spin={gjf_spin if gjf_spin is not None else spin}"
            )

        if resolved_charge is None and gjf_charge is not None:
            _echo(f"[all] Using total charge from GJF headers: {float(gjf_charge):+g}")
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
                f"[all] {charge_override_label} {q_int:+d} does not match extract-derived charge "
                f"{int(resolved_charge):+d}. Either {charge_override_remove_hint} to auto-derive, "
                f"pass --ligand-charge for non-standard residues, or set -q to "
                f"{int(resolved_charge):+d}."
            )
        override_msg = _charge_override_message(
            charge_override_label,
            q_int,
            None if resolved_charge is None else int(resolved_charge),
        )
        _echo(override_msg, err=True, narrative=True)
    else:
        q_int = int(resolved_charge) if resolved_charge is not None else 0

    _validate_path = model_outputs[0] if model_outputs else input_paths[0]
    validate_charge_spin_at_path(_validate_path, q_int, spin)

    # Resolve --ref-pdb for topology
    ref_pdb_for_topology: Optional[Path] = None
    if ref_pdb_cli is not None:
        prepared_ref = prepare_input_structure(ref_pdb_cli.resolve())
        prepared_all_inputs.append(prepared_ref)
        session.resources.own_cleanup(prepared_ref)
        ref_pdb_for_topology = prepared_ref.source_path
        _echo(f"[all] --ref-pdb provided: {ref_pdb_cli.resolve()}")

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
        tr_projection=(
            str(tr_projection).lower()
            if cli_param_overridden(ctx, "tr_projection")
            else None
        ),
        precision=precision,
        backend_model=backend_model,
        session=session,
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
    # Resolve model and precision in the shared calculator mapping so both
    # in-process evaluations and recorded provenance use the same settings.
    from pdb2reaction.backends import apply_backend_model_to_calc_cfg
    apply_backend_model_to_calc_cfg(calc_cfg_shared, backend_model)
    from pdb2reaction.backends import apply_effective_precision
    apply_effective_precision(calc_cfg_shared, precision)

    # --calc-file overrides --backend with a user ASE Calculator (custom backend).
    from pdb2reaction.backends import apply_calc_file_to_calc_cfg
    apply_calc_file_to_calc_cfg(calc_cfg_shared, calc_file, calc_factory)
    _shared_provenance = calculator_provenance(calc_cfg_shared)
    _mlip_backend_shared = str(_shared_provenance["mlip_backend"])
    _mlip_model_shared = _shared_provenance["mlip_model"]
    _mlip_precision_shared = _shared_provenance["mlip_precision"]

    # Preserve the source of parent CLI overrides. Child commands reread the
    # same YAML, so even values equal to child defaults must remain explicit.
    child_runtime_overrides: Dict[str, Any] = {}
    for _name, _value in (
        ("workers", workers),
        ("workers_per_node", workers_per_node),
        ("backend", backend),
        ("solvent", solvent),
        ("solvent_model", solvent_model),
    ):
        if cli_param_overridden(ctx, _name):
            child_runtime_overrides[_name] = _value
    freq_overrides.update(child_runtime_overrides)
    if calc_cfg_shared.get("calc_file"):
        freq_overrides["calc_file"] = calc_cfg_shared["calc_file"]
        if calc_cfg_shared.get("calc_factory"):
            freq_overrides["calc_factory"] = calc_cfg_shared["calc_factory"]
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
            runtime_overrides=child_runtime_overrides,
            manifest=manifest,
            artifact_prefix="post.01.tsopt",
            public_root=out_dir,
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
            irc_step_size=irc_step_size,
            irc_never_stop=irc_never_stop,
            seg_tag=None,
            session=session,
            manifest=manifest,
            artifact_prefix="post.01.irc",
            public_root=out_dir,
            runtime_overrides=child_runtime_overrides,
        )
        _persist_run_manifest(manifest, out_dir)
        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        calculator_lease = irc_res["calculator_lease"]
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

        # Record the endpoint-assignment provenance. In TS-only mode the
        # higher-energy IRC endpoint is presented as the reactant (with a
        # deterministic left-side tie rule). This is an energy-order presentation
        # convention, NOT a chemically established reaction direction; the R/P
        # labels, files, barrier and delta are unchanged, but a consumer must not
        # read chemical direction into them. See docs/all.md / docs/ja/all.md.
        endpoint_assignment = {
            "policy": "higher_energy_endpoint_as_reactant",
            "chemical_direction_known": False,
            "left_role": "reactant" if eL >= eR_raw else "product",
            "right_role": "product" if eL >= eR_raw else "reactant",
            "left_energy_hartree": float(eL),
            "right_energy_hartree": float(eR_raw),
            "tie_rule": "on equal energy (eL == eR) the left endpoint is the reactant",
        }

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
        from pdb2reaction.io.hessian_cache import (
            discard as _hess_discard,
            load as _hess_load,
            store as _hess_store,
        )
        _react_hk = "irc_left" if eL >= eR_raw else "irc_right"
        _prod_hk  = "irc_right" if eL >= eR_raw else "irc_left"

        _hess_discard("irc_endpoint")
        _c = _hess_load(_react_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
        _react_opt_conv: Optional[bool] = None
        try:
            g_react_opt, _, _react_opt_conv = _optimize_endpoint_geom(
                g_react_irc,
                tsopt_opt_mode_default,
                endpoint_opt_dir,
                "reactant",
                dump=dump,
                thresh=thresh_post,
                calc_identity_cfg=calc_cfg_shared,
                reject_uphill=_reject_uphill_eff,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Reactant endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            g_react_opt = g_react_irc
            _react_opt_conv = None

        _hess_discard("irc_endpoint")
        _c = _hess_load(_prod_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
        _prod_opt_conv: Optional[bool] = None
        try:
            g_prod_opt, _, _prod_opt_conv = _optimize_endpoint_geom(
                g_prod_irc,
                tsopt_opt_mode_default,
                endpoint_opt_dir,
                "product",
                dump=dump,
                thresh=thresh_post,
                calc_identity_cfg=calc_cfg_shared,
                reject_uphill=_reject_uphill_eff,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Product endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            g_prod_opt = g_prod_irc
            _prod_opt_conv = None

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

        diag_payload = _write_public_energy_diagram(
            tsroot / "energy_diagram_MLIP",
            labels=["R", "TS", "P"],
            energies_au=[e_react, eT, e_prod],
            title_note="(MLIP, TSOPT + IRC)",
        )
        if diag_payload:
            energy_diagrams.append(diag_payload)

        # Clear the calculator retained by Geometry before the frequency
        # subprocess so VRAM accounting does not include both stages.
        # ── Release GPU memory before freq/thermo/DFT ──
        calculator_lease.release()
        for _g in [locals().get(n) for n in ("gL", "gR", "gT", "g_react_irc", "g_prod_irc", "g_react_opt", "g_prod_opt")]:
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        dG_R = dG_T = dG_P = None
        GR_dft_mlip = GT_dft_mlip = GP_dft_mlip = None

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
            GR = _thermo_value_ha(tR, "sum_EE_and_thermal_free_energy_ha")
            GT = _thermo_value_ha(tT, "sum_EE_and_thermal_free_energy_ha")
            GP = _thermo_value_ha(tP, "sum_EE_and_thermal_free_energy_ha")
            dG_R = _thermo_value_ha(tR, "thermal_correction_free_energy_ha")
            dG_T = _thermo_value_ha(tT, "thermal_correction_free_energy_ha")
            dG_P = _thermo_value_ha(tP, "thermal_correction_free_energy_ha")
            try:
                if any(value is None for value in (GR, GT, GP)):
                    _echo(
                        "[thermo] NOTE: one or more R/TS/P free energies are "
                        "unavailable; Gibbs diagram skipped.",
                        err=True,
                    )
                else:
                    diag_payload = _write_public_energy_diagram(
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
            # Release calculator-owning references before the DFT subprocess.
            # Keep endpoint geometries and thermochemistry payloads: they are
            # still consumed by bond analysis, diagrams, and summary output.
            g_ts = None
            calc = None
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
                    diag_payload = _write_public_energy_diagram(
                        tsroot / "energy_diagram_DFT",
                        labels=["R", "TS", "P"],
                        energies_au=[eR_dft, eT_dft, eP_dft],
                        title_note=f"({dft_func_basis_use} // MLIP)",
                    )
                    if diag_payload:
                        energy_diagrams.append(diag_payload)
                except Exception as e:
                    _echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            if (
                do_thermo
                and _dft_all_ok
                and all(value is not None for value in (dG_R, dG_T, dG_P))
            ):
                try:
                    GR_dft_mlip = eR_dft + dG_R
                    GT_dft_mlip = eT_dft + dG_T
                    GP_dft_mlip = eP_dft + dG_P
                    if not all(np.isfinite([GR_dft_mlip, GT_dft_mlip, GP_dft_mlip])):
                        _echo(
                            "[dft//mlip] NOTE: thermochemistry unavailable; "
                            "DFT//MLIP Gibbs diagram skipped.",
                            err=True,
                        )
                    else:
                        diag_payload = _write_public_energy_diagram(
                            tsroot / "energy_diagram_G_DFT_plus_MLIP",
                            labels=["R", "TS", "P"],
                            energies_au=[GR_dft_mlip, GT_dft_mlip, GP_dft_mlip],
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

        n_images: Optional[int] = None
        try:
            irc_trj_path = irc_res.get("irc_trj_path")
            if isinstance(irc_trj_path, Path) and irc_trj_path.exists():
                n_images = len(read_xyz_as_blocks(irc_trj_path, strict=True))
        except Exception as exc:
            _echo(
                f"[post] WARNING: could not count complete IRC trajectory frames: {exc}",
                err=True,
            )

        summary = {
            "out_dir": str(tsroot),
            "n_images": n_images,
            "n_segments": 1,
            # Presentation-convention provenance for the R/P assignment.
            "endpoint_assignment": endpoint_assignment,
            "segments": [
                {
                    "index": 1,
                    "tag": "seg_01",
                    "kind": "tsopt",
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": _path_search._bond_changes_block(bond_summary),
                    "endpoint_assignment": endpoint_assignment,
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
            mlip_backend=_mlip_backend_shared,
            mlip_model=_mlip_model_shared,
            mlip_precision=_mlip_precision_shared,
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
                "path_opt_mode": opt_mode_norm,
                "post_opt_mode": tsopt_opt_mode_default,
                "ts_opt_mode": tsopt_opt_mode_default,
                "endpoint_opt_mode": tsopt_opt_mode_default,
                "mep_mode": mep_mode_kind,
                "dmf_correlated": dmf_correlated_effective,
            },
            freeze_atoms=_freeze_atoms_for_log(),
            manifest=manifest,
        )
        try:
            _publish_manifest_summary(
                tsroot / "summary.json",
                summary,
                manifest=manifest,
                key="ts.summary.01",
                out_dir=out_dir,
            )
            _echo_detail(f"[write] Wrote '{tsroot / 'summary.json'}'.")
            _copy_public_logged(tsroot / "summary.json", out_dir / "summary.json", label="summary.json")
            try:
                ts_freq_info = (
                    _read_imaginary_frequency(freq_root / "TS") if do_thermo else None
                )
                segment_log = {
                    "index": 1,
                    "tag": "seg_01",
                    "kind": "tsopt",
                    "bond_changes": summary["segments"][0].get("bond_changes"),
                    "post_dir": str(tsroot),
                    "irc_plot": str(irc_plot_path) if isinstance(irc_plot_path, Path) else None,
                    "irc_traj": str(irc_trj_path) if isinstance(irc_trj_path, Path) else None,
                    # IRC directional convergence and endpoint-opt
                    # convergence so the aggregate gates on convergence, not on
                    # trajectory-file existence.
                    "irc": irc_res.get("irc_outcome"),
                    "endpoint_opt": {
                        "reactant_converged": _react_opt_conv,
                        "product_converged": _prod_opt_conv,
                    },
                    # Presentation-convention provenance (energy-order R/P).
                    "endpoint_assignment": endpoint_assignment,
                }
                if ts_freq_info is not None:
                    segment_log["ts_imag"] = ts_freq_info
                    if ts_freq_info.get("nu_imag_max_cm") is not None:
                        segment_log["ts_imag_freq_cm"] = ts_freq_info["nu_imag_max_cm"]
                from pdb2reaction.workflows._all_helpers import (
                    build_energy_level_dict,
                    build_thermo_symmetry_provenance,
                )
                _thermo_symmetry = build_thermo_symmetry_provenance(
                    thermo_payloads
                )
                if _thermo_symmetry:
                    segment_log["thermo_symmetry"] = _thermo_symmetry
                _structs_seg = {"R": pR, "TS": pT, "P": pP}
                segment_log["mlip"] = build_energy_level_dict(
                    labels=["R", "TS", "P"],
                    energies_au=[e_react, eT, e_prod],
                    ref_energy=e_react,
                    au_to_kcal=AU2KCALPERMOL,
                    diagram_path=str(tsroot / "energy_diagram_MLIP.png"),
                    structures=_structs_seg,
                )
                if do_thermo and all(
                    value is not None for value in (GR, GT, GP)
                ):
                    segment_log["gibbs_mlip"] = build_energy_level_dict(
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
                    if (
                        do_thermo
                        and _dft_all_ok
                        and all(
                            value is not None
                            for value in (
                                GR_dft_mlip,
                                GT_dft_mlip,
                                GP_dft_mlip,
                            )
                        )
                    ):
                        segment_log["gibbs_dft_mlip"] = build_energy_level_dict(
                            labels=["R", "TS", "P"],
                            energies_au=[GR_dft_mlip, GT_dft_mlip, GP_dft_mlip],
                            ref_energy=GR_dft_mlip,
                            au_to_kcal=AU2KCALPERMOL,
                            diagram_path=str(tsroot / "energy_diagram_G_DFT_plus_MLIP.png"),
                            structures=_structs_seg,
                        )

                summary_payload = {
                    "root_out_dir": str(out_dir),
                    "path_dir": str(tsroot),
                    "path_module_dir": tsroot.name,
                    "pipeline_mode": "tsopt-only",
                    "n_images": n_images,
                    "n_segments": 1,
                    "refine_path": refine_path,
                    "tsopt": do_tsopt,
                    "thermo": do_thermo,
                    "dft": do_dft,
                    "opt_mode": tsopt_opt_mode_default,
                    "opt_mode_post": (
                        opt_mode_post.lower() if opt_mode_post else None
                    ),
                    "post_opt_mode": tsopt_opt_mode_default,
                    "ts_opt_mode": tsopt_opt_mode_default,
                    "endpoint_opt_mode": tsopt_opt_mode_default,
                    "mep_mode": mep_mode_kind,
                    "dmf_correlated": dmf_correlated_effective,
                    "mlip_backend": _mlip_backend_shared,
                    "mlip_model": _mlip_model_shared,
                    "mlip_precision": _mlip_precision_shared,
                    "command": command_str,
                    "charge": q_int,
                    "spin": spin,
                    "freeze_atoms": _freeze_atoms_for_log(),
                    "mep": {"n_images": n_images, "n_segments": 1},
                    "segments": summary.get("segments", []),
                    "energy_diagrams": summary.get("energy_diagrams", []),
                    "post_segments": [segment_log],
                    "key_files": {},
                    "current_output_paths": [
                        path.relative_to(out_dir).as_posix()
                        for path in _refresh_current_public_outputs(
                            manifest, out_dir
                        )
                        if path.is_relative_to(out_dir)
                    ],
                }
                # Copy R/TS/P to seg_01/ BEFORE writing summary.log (so tree includes seg_01/)
                try:
                    _input_suffix = first_input.suffix.lower() if first_input else ".xyz"
                    _tsonly_structs = {"R": pR, "TS": pT, "P": pP}
                    _seg_out = _copy_structures_to_seg_dir(
                        _tsonly_structs,
                        out_dir,
                        1,
                        _input_suffix,
                        manifest=manifest,
                    )
                    _echo(f"[all] Wrote R/TS/P for segment 01 → {_seg_out}", narrative=True)
                except Exception as e:
                    _echo(
                        f"[all] WARNING: Failed to copy structures in TSOPT-only mode: {e}",
                        err=True,
                    )

                # Re-enrich after the real TS/IRC post-processing payload exists.
                # The first pass cannot select thermochemistry/DFT provenance or
                # derive the corresponding local barrier/reaction energy.
                _enrich_summary(
                    summary,
                    version="",
                    pipeline_mode="tsopt-only",
                    out_dir=out_dir,
                    mlip_backend=_mlip_backend_shared,
                    mlip_model=_mlip_model_shared,
                    mlip_precision=_mlip_precision_shared,
                    charge=q_int,
                    spin=spin,
                    command=command_str,
                    post_segments=[segment_log],
                    config={
                        "refine_path": refine_path,
                        "tsopt": do_tsopt,
                        "thermo": do_thermo,
                        "dft": do_dft,
                        "dft_status": (
                            "failed"
                            if do_dft and not _dft_all_ok
                            else (
                                "converged"
                                if do_dft and _dft_all_ok
                                else None
                            )
                        ),
                        "dft_func_basis": (
                            dft_func_basis_use if do_dft else None
                        ),
                        "opt_mode": tsopt_opt_mode_default,
                        "path_opt_mode": opt_mode_norm,
                        "post_opt_mode": tsopt_opt_mode_default,
                        "ts_opt_mode": tsopt_opt_mode_default,
                        "endpoint_opt_mode": tsopt_opt_mode_default,
                        "mep_mode": mep_mode_kind,
                        "dmf_correlated": dmf_correlated_effective,
                    },
                    freeze_atoms=_freeze_atoms_for_log(),
                    manifest=manifest,
                )
                try:
                    _publish_manifest_summary(
                        tsroot / "summary.json",
                        summary,
                        manifest=manifest,
                        key="ts.summary.01",
                        out_dir=out_dir,
                    )
                    _copy_public_logged(tsroot / "summary.json", out_dir / "summary.json", label="summary.json", echo=False)
                except Exception as exc:
                    _echo(
                        f"[write] WARNING: Failed to refresh summary.json in TSOPT-only mode: {exc}",
                        err=True,
                    )

                write_summary_log(tsroot / "summary.log", summary_payload)
                _copy_public_logged(tsroot / "summary.log", out_dir / "summary.log", label="summary.log", echo=False)
            except Exception as e:
                _echo(f"[write] WARNING: Failed to write summary.log in TSOPT-only mode: {e}", err=True)
        except Exception as e:
            _echo(
                f"[write] WARNING: Failed to write summary.json for TSOPT-only run: {e}",
                err=True,
            )

        try:
            current_public_set = set(
                _refresh_current_public_outputs(manifest, out_dir)
            )
            for stem in (
                "energy_diagram_MLIP",
                "energy_diagram_G_MLIP",
                "energy_diagram_DFT",
                "energy_diagram_G_DFT_plus_MLIP",
            ):
                src = tsroot / f"{stem}.png"
                if src.resolve(strict=False) in current_public_set:
                    dst = out_dir / f"{stem}_all.png"
                    _copy_public_logged(src, dst, label=src.name)
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to mirror *_all diagrams in TSOPT-only mode: {e}",
                err=True,
            )

        try:
            if isinstance(irc_plot_path, Path):
                dst = out_dir / "irc_plot_all.png"
                _copy_public_logged(irc_plot_path, dst, label="irc_plot_all.png")
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to mirror IRC plot in TSOPT-only mode: {e}",
                err=True,
            )

        summary["key_output_files"] = _current_key_output_files(
            manifest, out_dir
        )
        summary["current_output_paths"] = _current_output_paths(
            manifest, out_dir
        )
        _publish_manifest_summary(
            tsroot / "summary.json",
            summary,
            manifest=manifest,
            key="ts.summary.01",
            out_dir=out_dir,
        )
        _copy_public_logged(
            tsroot / "summary.json",
            out_dir / "summary.json",
            label="summary.json",
            echo=False,
        )
        summary_payload["current_output_paths"] = [
            path.relative_to(out_dir).as_posix()
            for path in _refresh_current_public_outputs(manifest, out_dir)
            if path.is_relative_to(out_dir)
        ]
        write_summary_log(tsroot / "summary.log", summary_payload)
        _copy_public_logged(
            tsroot / "summary.log",
            out_dir / "summary.log",
            label="summary.log",
            echo=False,
        )
        _refresh_current_public_outputs(manifest, out_dir)
        _persist_run_manifest(manifest, out_dir)

        _echo_section(
            "====== [all] TSOPT-only pipeline successfully finished ======"
        )
        citation_post_segments = [segment_log]
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    # Stage 1b: optional staged scan (single-structure)
    models_for_path: List[Path]
    model_ref_pdbs: Optional[List[Path]] = None
    if is_single and has_scan:
        _echo_section("====== [all] Stage 1b — Staged scan on input ======")
        ensure_dir(scan_dir)

        # --scan-one-based decides how the user's literal is READ, so it has to be known
        # before parsing. Default to 1-based when unspecified.
        scan_one_based_effective = True if scan_one_based is None else bool(
            scan_one_based
        )

        if skip_extract:
            scan_input_pdb = Path(input_paths[0]).resolve()
            scan_atom_meta = None
            if scan_input_pdb.suffix.lower() == ".pdb":
                scan_atom_meta = load_pdb_atom_metadata(scan_input_pdb)
            converted_scan_stages = _parse_scan_lists_literals(
                scan_lists_raw,
                atom_meta=scan_atom_meta,
                one_based=scan_one_based_effective,
            )
        else:
            scan_input_pdb = Path(model_outputs[0]).resolve()
            full_input_pdb = Path(input_paths[0]).resolve()
            converted_scan_stages = _convert_scan_lists_to_model_indices(
                scan_lists_raw,
                full_input_pdb,
                scan_input_pdb,
                one_based=scan_one_based_effective,
            )
            _echo_detail(
                "[all] Remapped --scan-lists indices from the full PDB to the active site model ordering."
            )

        # Both parse paths return 1-based indices regardless of how the literal was read, so
        # on an explicit 0-based request we decrement once here and forward --zero-based to
        # scan.py below. The two conversions are inverse, not additive.
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
            "--out-json",
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
        _append_explicit_child_runtime_argv(
            scan_args,
            workers=workers if cli_param_overridden(ctx, "workers") else None,
            workers_per_node=(
                workers_per_node
                if cli_param_overridden(ctx, "workers_per_node")
                else None
            ),
        )
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

        for stage_idx in range(1, len(scan_stage_literals) + 1):
            stage_root = scan_dir / f"stage_{stage_idx:02d}" / "result"
            manifest.declare(
                f"scan.stage.{stage_idx:02d}",
                [stage_root.with_suffix(suffix) for suffix in (".xyz", ".pdb", ".gjf")],
            )
            manifest.declare(
                f"scan.stage_ref.{stage_idx:02d}",
                [stage_root.with_suffix(".pdb")],
            )
        if scan_preopt_use:
            preopt_root = scan_dir / "preopt" / "result"
            manifest.declare(
                "scan.preopt",
                [preopt_root.with_suffix(suffix) for suffix in (".xyz", ".pdb", ".gjf")],
            )
            manifest.declare("scan.preopt_ref", [preopt_root.with_suffix(".pdb")])
        manifest.declare("scan.result", [scan_dir / "result.json"])

        _run_cli_main("scan", _scan_cli.cli, scan_args, on_nonzero="raise", prefix="all")

        scan_result_path = manifest.claim_one("scan.result")
        try:
            scan_result = json.loads(
                scan_result_path.read_text(encoding="utf-8")
            )
        except (OSError, ValueError, TypeError) as exc:
            raise click.ClickException(
                f"[all] Could not read scan outcome from {scan_result_path}: {exc}"
            ) from exc
        if scan_result.get("scientific_status") != "success":
            reasons = scan_result.get("scientific_status_reasons") or []
            detail = "; ".join(str(reason) for reason in reasons) or "unknown reason"
            raise click.ClickException(
                f"[all] Staged scan did not produce a scientifically usable path: {detail}"
            )
        scan_preopt_usable = scan_result.get("preopt_converged") is True

        stage_results = [
            manifest.claim_one(f"scan.stage.{stage_idx:02d}")
            for stage_idx in range(1, len(scan_stage_literals) + 1)
        ]
        stage_ref_results = {
            stage_idx: manifest.claim_optional(f"scan.stage_ref.{stage_idx:02d}")
            for stage_idx in range(1, len(scan_stage_literals) + 1)
        }
        preopt_result = manifest.claim_optional("scan.preopt") if scan_preopt_use else None
        preopt_ref_result = (
            manifest.claim_optional("scan.preopt_ref") if scan_preopt_use else None
        )
        _persist_run_manifest(manifest, out_dir)
        _echo_detail("[all] Collected scan stage active site model files:")
        for p in stage_results:
            _echo_detail(f"  - {p}")

        initial_path_for_path = scan_input_pdb
        initial_ref_pdb_for_path = ref_pdb_for_topology or (scan_input_pdb if scan_input_pdb.suffix.lower() == ".pdb" else None)
        if preopt_result is not None and scan_preopt_usable:
            initial_path_for_path = preopt_result
            if preopt_ref_result is not None:
                initial_ref_pdb_for_path = preopt_ref_result
            _echo_detail(f"[all] Using current scan preopt result as initial path endpoint: {initial_path_for_path}")
        models_for_path = [initial_path_for_path] + stage_results

        if initial_ref_pdb_for_path is not None:
            candidate_pdbs: List[Path] = [initial_ref_pdb_for_path]
            missing_pdb = False
            for stage_idx, stage_path in enumerate(stage_results, start=1):
                stage_ref = stage_ref_results.get(stage_idx)
                if stage_ref is not None:
                    candidate_pdbs.append(stage_ref)
                elif stage_path.suffix.lower() == ".pdb":
                    candidate_pdbs.append(stage_path)
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
                alignment_results = align_and_refine_sequence_inplace(
                    _geoms, shared_calc=_align_calc,
                    out_dir=_align_dir / "refine", verbose=True,
                )
                failed_pairs = alignment_failed_pair_indices(alignment_results)
                if failed_pairs:
                    raise click.ClickException(
                        "Input alignment did not converge for pair(s): "
                        + ", ".join(str(index) for index in failed_pairs)
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
                raise click.ClickException(
                    f"Path pre-alignment failed: {e}"
                ) from e

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

        _declare_path_deliverables(manifest, path_dir)
        if "path.summary" not in manifest.expected:
            manifest.declare("path.summary", [path_dir / "summary.json"])
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
                "--opt-mode",
                str(opt_mode_norm),
                "--out-dir",
                str(seg_dir),
                # Pipeline-owned machine contract: the aggregate reads this
                # child's real MEP convergence from result.json.
                "--out-json",
            ]
            if cli_param_overridden(ctx, "max_cycles"):
                po_args.extend(["--max-cycles", str(int(max_cycles))])
            _append_toggle_arg(po_args, "--freeze-links", bool(freeze_links_flag and freeze_ref is not None))
            if cli_param_overridden(ctx, "climb"):
                _append_toggle_arg(po_args, "--climb", bool(climb))
            _append_toggle_arg(po_args, "--dump", bool(dump))
            _append_toggle_arg(po_args, "--convert-files", bool(convert_files))
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
            _append_explicit_child_runtime_argv(
                po_args,
                dmf_backend=(
                    dmf_backend if cli_param_overridden(ctx, "dmf_backend") else None
                ),
                workers=workers if cli_param_overridden(ctx, "workers") else None,
                workers_per_node=(
                    workers_per_node
                    if cli_param_overridden(ctx, "workers_per_node")
                    else None
                ),
                max_nodes=(
                    max_nodes if cli_param_overridden(ctx, "max_nodes") else None
                ),
                preopt=preopt if cli_param_overridden(ctx, "preopt") else None,
            )
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

            seg_trj_expected = seg_dir / "final_geometries_trj.xyz"
            child_hei_base = seg_dir / "hei"
            manifest.declare(
                f"path.segment.{idx:02d}.trajectory",
                [seg_trj_expected],
            )
            manifest.declare(
                f"path.segment.{idx:02d}.endpoint_trajectory",
                [seg_trj_expected],
            )
            manifest.declare(
                f"path.segment.{idx:02d}.hei_child",
                [child_hei_base.with_suffix(suffix) for suffix in (".xyz", ".pdb", ".gjf")],
            )
            manifest.declare(
                f"path.segment.{idx:02d}.hei_child_ref",
                [child_hei_base.with_suffix(".pdb")],
            )
            manifest.declare(
                f"path.segment.{idx:02d}.hei_child_gjf",
                [child_hei_base.with_suffix(".gjf")],
            )
            manifest.declare(
                f"path.segment.{idx:02d}.result",
                [seg_dir / "result.json"],
            )

            _run_cli_main(
                "path-opt",
                _path_opt.cli,
                po_args,
                on_nonzero="raise",
                prefix=f"all seg {idx:02d}",
            )

            seg_trj = manifest.claim_one(f"path.segment.{idx:02d}.trajectory")
            child_hei = manifest.claim_optional(f"path.segment.{idx:02d}.hei_child")
            child_hei_ref = manifest.claim_optional(
                f"path.segment.{idx:02d}.hei_child_ref"
            )
            child_hei_gjf = manifest.claim_optional(
                f"path.segment.{idx:02d}.hei_child_gjf"
            )
            manifest.claim_one(f"path.segment.{idx:02d}.result")

            try:
                mirror_dir = path_dir / f"{seg_tag}_mep"
                mirror_trj = mirror_dir / "final_geometries_trj.xyz"

                ensure_dir(mirror_dir)
                if seg_trj.resolve() != mirror_trj.resolve():
                    shutil.copy2(seg_trj, mirror_trj)
                manifest.claim_one(f"path.segment.{idx:02d}.endpoint_trajectory")
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to mirror path-opt trajectory for segment {idx:02d}: {e}",
                    err=True,
                )

            try:
                seg_mep_trj = path_dir / f"mep_seg_{idx:02d}_trj.xyz"
                manifest.declare(f"path.mep.{idx:02d}.trajectory", [seg_mep_trj])
                shutil.copy2(seg_trj, seg_mep_trj)
                manifest.claim_one(f"path.mep.{idx:02d}.trajectory")
                if models_for_path[0].suffix.lower() == ".pdb":
                    seg_mep_pdb = path_dir / f"mep_seg_{idx:02d}.pdb"
                    manifest.declare(f"path.mep.{idx:02d}.pdb", [seg_mep_pdb])
                    _path_search._convert_to_pdb_logged(
                        seg_mep_trj,
                        ref_pdb_path=models_for_path[0],
                        out_path=seg_mep_pdb,
                    )
                    manifest.claim_one(f"path.mep.{idx:02d}.pdb")
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to emit per-segment trajectory copies for segment {idx:02d}: {e}",
                    err=True,
                )

            if child_hei is not None:
                try:
                    hei_destination = path_dir / f"hei_seg_{idx:02d}{child_hei.suffix}"
                    manifest.declare(f"path.hei.{idx:02d}", [hei_destination])
                    shutil.copy2(child_hei, hei_destination)
                    manifest.claim_one(f"path.hei.{idx:02d}")
                    if child_hei_ref is not None and child_hei_ref != child_hei:
                        hei_pdb_destination = path_dir / f"hei_seg_{idx:02d}.pdb"
                        manifest.declare(
                            f"path.hei_ref.{idx:02d}", [hei_pdb_destination]
                        )
                        shutil.copy2(child_hei_ref, hei_pdb_destination)
                        manifest.claim_one(f"path.hei_ref.{idx:02d}")
                    if child_hei_gjf is not None and child_hei_gjf != child_hei:
                        hei_gjf_destination = path_dir / f"hei_seg_{idx:02d}.gjf"
                        manifest.declare(
                            f"path.hei_gjf.{idx:02d}", [hei_gjf_destination]
                        )
                        shutil.copy2(child_hei_gjf, hei_gjf_destination)
                        manifest.claim_one(f"path.hei_gjf.{idx:02d}")
                except Exception as e:
                    _echo(
                        f"[all] WARNING: failed to prepare HEI artifacts for segment {idx:02d}: {e}",
                        err=True,
                    )

            _persist_run_manifest(manifest, out_dir)

            raw_blocks = read_xyz_as_blocks(seg_trj, strict=True)
            blocks = ["\n".join(b) + "\n" for b in raw_blocks]
            if not blocks:
                raise click.ClickException(
                    f"[all] No frames read from path-opt segment {idx} trajectory: {seg_trj}"
                )
            if idx > 1:
                blocks = blocks[1:]
            combined_blocks.extend(blocks)

            energies_seg = _required_xyz_block_energies(
                raw_blocks,
                path=seg_trj,
                context=f"path-opt segment {idx}",
            )

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
                    # Per-segment convergence from the path-opt child, so
                    # the no-tsopt aggregate fails closed on a nonconverged MEP
                    # segment instead of assuming it converged.
                    "converged": _read_path_opt_segment_converged(seg_dir),
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
                current_mep_pdb = manifest.claim_optional(
                    "path.deliverable.mep.pdb"
                )
                if mep_pdb and current_mep_pdb is not None:
                    dst = out_dir / current_mep_pdb.name
                    _copy_public_logged(current_mep_pdb, dst, echo=False)
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
                energies_chain.append(max(Es))
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
            Es = [float(x) for x in info.get("energies", [])]
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
                    # Use the path-opt child's explicit convergence so
                    # the no-tsopt path-opt aggregate gates on real per-segment
                    # convergence rather than assuming the segment converged.
                    "converged": info.get("converged"),
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
            mlip_backend=_mlip_backend_shared,
            mlip_model=_mlip_model_shared,
            mlip_precision=_mlip_precision_shared,
            charge=q_int,
            spin=spin,
            command=command_str,
            config={
                "refine_path": refine_path,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": tsopt_opt_mode_default,
                "path_opt_mode": opt_mode_norm,
                "post_opt_mode": tsopt_opt_mode_default,
                "ts_opt_mode": tsopt_opt_mode_default,
                "endpoint_opt_mode": tsopt_opt_mode_default,
                "mep_mode": mep_mode_kind,
                "dmf_correlated": dmf_correlated_effective,
            },
            freeze_atoms=_freeze_atoms_for_log(),
            manifest=manifest,
        )
        _publish_manifest_summary(
            path_dir / "summary.json",
            summary,
            manifest=manifest,
            key="path.summary",
            out_dir=out_dir,
        )
        _echo_detail(f"[write] Wrote '{path_dir / 'summary.json'}'.")

        try:
            # MEP deliverables are MOVED to the pipeline root so they are not
            # duplicated under _work/path_opt (path_dir).  The path-opt branch
            # authors the energy diagram lower-cased (energy_diagram_mep.png);
            # glob both casings and land it at root under the canonical name.
            path_deliverables = _claim_path_deliverables(manifest)
            diagram = path_deliverables.get("diagram")
            if diagram is not None:
                _move_public_logged(
                    diagram,
                    out_dir / "energy_diagram_MEP.png",
                    label=diagram.name,
                )
            for name, src in path_deliverables.items():
                if name == "diagram":
                    continue
                _move_public_logged(src, out_dir / src.name, label=src.name)

            # summary.json / summary.log stay COPIES: the path_dir copy is
            # re-read (segments) and re-authored later, and the root copy is
            # consumed by the final-summary banner.
            _copy_public_logged(
                manifest.path("path.summary"),
                out_dir / "summary.json",
                label="summary.json",
            )
            _refresh_current_public_outputs(manifest, out_dir)
            _persist_run_manifest(manifest, out_dir)
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
            current_public = _refresh_current_public_outputs(
                manifest, out_dir
            )
            current_public_set = set(current_public)
            # MEP products were moved up to out_dir above; reference the root.
            mep_info = {
                "n_images": summary.get("n_images"),
                "n_segments": summary.get("n_segments"),
                "traj_pdb": str(out_dir / "mep.pdb") if (out_dir / "mep.pdb").resolve(strict=False) in current_public_set else None,
                "traj_cif": str(out_dir / "mep.cif") if (out_dir / "mep.cif").resolve(strict=False) in current_public_set else None,
                "mep_plot": str(out_dir / "energy_diagram_MEP.png") if (out_dir / "energy_diagram_MEP.png").resolve(strict=False) in current_public_set else None,
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
                "opt_mode_post": (
                    opt_mode_post.lower() if opt_mode_post else None
                ),
                "path_opt_mode": opt_mode_norm,
                "post_opt_mode": tsopt_opt_mode_default,
                "ts_opt_mode": tsopt_opt_mode_default,
                "endpoint_opt_mode": tsopt_opt_mode_default,
                "mep_mode": mep_mode_kind,
                "dmf_correlated": dmf_correlated_effective,
                "mlip_backend": _mlip_backend_shared,
                "mlip_model": _mlip_model_shared,
                "mlip_precision": _mlip_precision_shared,
                "command": command_str,
                "charge": q_int,
                "spin": spin,
                "freeze_atoms": _freeze_atoms_for_log(),
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {},
                "current_output_paths": [
                    path.relative_to(out_dir).as_posix()
                    for path in current_public
                    if path.is_relative_to(out_dir)
                ],
            }
            write_summary_log(path_dir / "summary.log", summary_payload)
            _copy_public_logged(path_dir / "summary.log", out_dir / "summary.log", label="summary.log")
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
        if cli_param_overridden(ctx, "max_cycles"):
            ps_args.extend(["--max-cycles", str(int(max_cycles))])
        if cli_param_overridden(ctx, "climb"):
            _append_toggle_arg(ps_args, "--climb", bool(climb))
        ps_args.extend(["--opt-mode", str(opt_mode_norm)])
        _append_toggle_arg(ps_args, "--dump", bool(dump))
        if thresh is not None:
            ps_args.extend(["--thresh", str(thresh)])
        ps_args.extend(["--out-dir", str(path_dir)])
        _append_toggle_arg(ps_args, "--convert-files", bool(convert_files))
        if args_yaml is not None:
            ps_args.extend(["--config", str(args_yaml)])
        if not _forward_calc_file_argv(ps_args, calc_cfg_shared) and cli_param_overridden(ctx, "backend"):
            ps_args.extend(["--backend", str(backend)])
        _append_explicit_child_runtime_argv(
            ps_args,
            dmf_backend=(
                dmf_backend if cli_param_overridden(ctx, "dmf_backend") else None
            ),
            workers=workers if cli_param_overridden(ctx, "workers") else None,
            workers_per_node=(
                workers_per_node
                if cli_param_overridden(ctx, "workers_per_node")
                else None
            ),
            max_nodes=(
                max_nodes if cli_param_overridden(ctx, "max_nodes") else None
            ),
            preopt=preopt if cli_param_overridden(ctx, "preopt") else None,
        )
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

        path_stage_candidates: List[Path] = [path_dir / "summary.json"]
        for pattern in (
            "hei_seg_*.*",
            "mep_seg_*_trj.xyz",
            "seg_*_mep/final_geometries_trj.xyz",
            "seg_*_refine_mep/final_geometries_trj.xyz",
        ):
            path_stage_candidates.extend(path_dir.glob(pattern))
        for name in (
            "energy_diagram_MEP.png",
            "energy_diagram_mep.png",
            "mep.pdb",
            "mep.cif",
            "mep_w_ref.pdb",
            "mep_w_ref.cif",
            "mep_trj.xyz",
            "mep.xyz",
            "mep_w_ref_trj.xyz",
            "mep_w_ref.xyz",
        ):
            path_stage_candidates.append(path_dir / name)
        path_stage_snapshot = InvocationManifest.snapshot(path_stage_candidates)
        if "path.summary" not in manifest.expected:
            manifest.declare(
                "path.summary",
                [path_dir / "summary.json"],
                snapshot=path_stage_snapshot,
            )
        _declare_path_deliverables(
            manifest,
            path_dir,
            snapshot=path_stage_snapshot,
        )

        _run_cli_main(
            "path_search",
            _path_search.cli,
            ps_args,
            on_nonzero="raise",
            prefix="all",
        )

        claimed_summary = manifest.claim_one("path.summary")
        try:
            path_summary_payload = json.loads(
                claimed_summary.read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as exc:
            raise click.ClickException(
                f"[all] Could not read current path-search summary '{claimed_summary}': {exc}"
            ) from exc
        if not isinstance(path_summary_payload, dict):
            raise click.ClickException(
                f"[all] Current path-search summary is not a JSON object: {claimed_summary}"
            )
        path_segments = path_summary_payload.get("segments") or []
        if not isinstance(path_segments, list):
            raise click.ClickException(
                f"[all] Current path-search summary has a non-list segments field: {claimed_summary}"
            )
        for segment in path_segments:
            if not isinstance(segment, dict):
                continue
            try:
                current_idx = int(segment.get("index", 0) or 0)
            except (TypeError, ValueError):
                continue
            if current_idx <= 0:
                continue
            current_tag = str(segment.get("tag") or f"seg_{current_idx:03d}")
            base_tag = _path_search._segment_base_id(current_tag)
            hei_base = path_dir / f"hei_seg_{current_idx:02d}"
            manifest.declare(
                f"path.hei.{current_idx:02d}",
                [hei_base.with_suffix(suffix) for suffix in (".xyz", ".pdb", ".gjf")],
                snapshot=path_stage_snapshot,
            )
            manifest.declare(
                f"path.hei_ref.{current_idx:02d}",
                [hei_base.with_suffix(".pdb")],
                snapshot=path_stage_snapshot,
            )
            manifest.declare(
                f"path.hei_gjf.{current_idx:02d}",
                [hei_base.with_suffix(".gjf")],
                snapshot=path_stage_snapshot,
            )
            manifest.declare(
                f"path.mep.{current_idx:02d}.trajectory",
                [path_dir / f"mep_seg_{current_idx:02d}_trj.xyz"],
                snapshot=path_stage_snapshot,
            )
            manifest.declare(
                f"path.segment.{current_idx:02d}.endpoint_trajectory",
                [
                    path_dir / f"{base_tag}_refine_mep" / "final_geometries_trj.xyz",
                    path_dir / f"{base_tag}_mep" / "final_geometries_trj.xyz",
                ],
                snapshot=path_stage_snapshot,
            )
            manifest.claim_optional(f"path.hei.{current_idx:02d}")
            manifest.claim_optional(f"path.hei_ref.{current_idx:02d}")
            manifest.claim_optional(f"path.hei_gjf.{current_idx:02d}")
            manifest.claim_optional(f"path.mep.{current_idx:02d}.trajectory")
            manifest.claim_optional(
                f"path.segment.{current_idx:02d}.endpoint_trajectory"
            )
        path_deliverables = _claim_path_deliverables(manifest)
        _persist_run_manifest(manifest, out_dir)

        try:
            # MEP deliverables are MOVED to the pipeline root so they are not
            # duplicated under _work/path_search (path_dir).  path_search emits
            # the canonical energy_diagram_MEP.png; glob both casings defensively
            # and land it at root under the canonical name.
            diagram = path_deliverables.get("diagram")
            if diagram is not None:
                _move_public_logged(
                    diagram,
                    out_dir / "energy_diagram_MEP.png",
                    label=diagram.name,
                )
            for name, src in path_deliverables.items():
                if name == "diagram":
                    continue
                _move_public_logged(src, out_dir / src.name, label=src.name)

            # summary.json / summary.log stay COPIES: the path_dir copy is
            # re-read (segments) and re-authored later, and the root copy is
            # consumed by the final-summary banner.
            _copy_public_logged(
                claimed_summary,
                out_dir / "summary.json",
                label="summary.json",
            )
            _refresh_current_public_outputs(manifest, out_dir)
            _persist_run_manifest(manifest, out_dir)
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
            "  - mep_w_ref.pdb/.cif       (full-system merged trajectory)"
        )
        _echo_detail(f"[all] Raw per-segment merged trajectories stay under: {path_dir}")
        _echo_detail(
            "  - mep_w_ref_seg_XX.pdb/.cif (per-segment merged trajectories for covalent-change segments)"
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

    summary_path = manifest.path("path.summary")
    try:
        summary_loaded = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise click.ClickException(
            f"[all] Could not read current path summary '{summary_path}': {exc}"
        ) from exc
    if not isinstance(summary_loaded, dict):
        raise click.ClickException(
            f"[all] Current path summary is not a JSON object: {summary_path}"
        )
    summary: Dict[str, Any] = summary_loaded
    _publish_manifest_summary(
        summary_path,
        summary,
        manifest=manifest,
        key="path.summary",
        out_dir=out_dir,
    )
    _copy_public_logged(summary_path, out_dir / "summary.json", label="summary.json", echo=False)
    _refresh_current_public_outputs(manifest, out_dir)
    _persist_run_manifest(manifest, out_dir)
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
        nonlocal citation_post_segments
        citation_post_segments = list(post_segment_logs)
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
            opt_mode_post=opt_mode_post,
            path_opt_mode=opt_mode_norm,
            post_opt_mode=tsopt_opt_mode_default,
            ts_opt_mode=tsopt_opt_mode_default,
            endpoint_opt_mode=tsopt_opt_mode_default,
            mep_mode_kind=mep_mode_kind,
            dmf_correlated=dmf_correlated_effective,
            mlip_backend=_mlip_backend_shared,
            mlip_model=_mlip_model_shared,
            mlip_precision=_mlip_precision_shared,
            command_str=command_str,
            q_int=q_int,
            spin=spin,
            freeze_atoms=_freeze_atoms_for_log(),
            post_segment_logs=post_segment_logs,
        )
        current_public = _refresh_current_public_outputs(manifest, out_dir)
        current_public_set = set(current_public)
        summary_payload["current_output_paths"] = [
            path.relative_to(out_dir).as_posix()
            for path in current_public
            if path.is_relative_to(out_dir)
        ]
        mep_payload = summary_payload.get("mep") or {}
        for field, name in (
            ("traj_pdb", "mep.pdb"),
            ("mep_plot", "energy_diagram_MEP.png"),
        ):
            candidate = (out_dir / name).resolve(strict=False)
            mep_payload[field] = (
                str(out_dir / name) if candidate in current_public_set else None
            )
        write_summary_log(path_dir / "summary.log", summary_payload)
        copied = _copy_public_logged(
            path_dir / "summary.log",
            out_dir / "summary.log",
            label="summary.log",
            echo=False,
        )
        if not copied:
            raise click.ClickException(
                f"Failed to publish summary.log to {out_dir / 'summary.log'}."
            )
        _echo_detail(f"[all] Copied summary.log → {out_dir / 'summary.log'}")

    def _finalize_current_summary() -> None:
        summary["key_output_files"] = _current_key_output_files(
            manifest, out_dir
        )
        summary["current_output_paths"] = _current_output_paths(
            manifest, out_dir
        )
        _publish_manifest_summary(
            summary_path,
            summary,
            manifest=manifest,
            key="path.summary",
            out_dir=out_dir,
        )
        _copy_public_logged(
            summary_path,
            out_dir / "summary.json",
            label="summary.json",
            echo=False,
        )
        _refresh_current_public_outputs(manifest, out_dir)
        _persist_run_manifest(manifest, out_dir)

    if not (do_tsopt or do_thermo or do_dft):
        if energy_diagrams:
            summary["energy_diagrams"] = list(energy_diagrams)
        _write_pipeline_summary_log([])
        _finalize_current_summary()
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    # Stage 4: post-processing per reactive segment
    _echo_section(
        f"====== [all] Stage 4/{stage_total} — Post-processing per reactive segment ======"
    )

    if not segments:
        _echo("[post] No segments found in summary; nothing to do.", narrative=True)
        _write_pipeline_summary_log([])
        _finalize_current_summary()
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    reactive = [s for s in segments if _is_reactive_segment(s)]
    if not reactive:
        _echo("[post] No bond-change segments. Skipping TS/thermo/DFT.", narrative=True)
        _write_pipeline_summary_log([])
        _finalize_current_summary()
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    # Per-category per-segment energies
    tsopt_seg_energies: List[Tuple[float, float, float]] = []
    g_mlip_seg_energies: List[Tuple[float, float, float]] = []
    dft_seg_energies: List[Tuple[float, float, float]] = []
    g_dft_mlip_seg_energies: List[Tuple[float, float, float]] = []
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

        hei_key = f"path.hei.{seg_idx:02d}"
        hei_model_path = (
            manifest.produced[hei_key][0] if hei_key in manifest.produced else None
        )
        if hei_model_path is None:
            _echo(
                f"[post] WARNING: current invocation did not claim an HEI model for segment {seg_idx:02d}; skipping TSOPT.",
                err=True,
            )
            continue
        ref_pdb_for_seg: Optional[Path] = None
        if hei_model_path.suffix.lower() == ".pdb":
            ref_pdb_for_seg = hei_model_path
        else:
            hei_ref_key = f"path.hei_ref.{seg_idx:02d}"
            if hei_ref_key in manifest.produced:
                ref_pdb_for_seg = manifest.produced[hei_ref_key][0]
            elif ref_pdb_for_topology is not None:
                ref_pdb_for_seg = ref_pdb_for_topology

        struct_dir = seg_dir / "structures"
        ensure_dir(struct_dir)
        state_structs: Dict[str, Path] = {}

        if do_tsopt:
            _seg_tsopt_overrides = dict(tsopt_overrides)
            _hei_mode_path = seg_root / f"hei_mode_seg_{seg_idx:02d}.txt"
            mep_key = f"path.mep.{seg_idx:02d}.trajectory"
            current_mep = (
                manifest.produced[mep_key][0] if mep_key in manifest.produced else None
            )
            _hei_mode_path = (
                _ensure_hei_path_tangent(
                    current_mep,
                    hei_model_path,
                    _hei_mode_path,
                )
                if current_mep is not None
                else None
            )
            if _hei_mode_path is not None:
                _seg_tsopt_overrides["reference_mode"] = _hei_mode_path
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
                overrides=_seg_tsopt_overrides,
                runtime_overrides=child_runtime_overrides,
                manifest=manifest,
                artifact_prefix=f"post.{seg_idx:02d}.tsopt",
                public_root=out_dir,
            )

            endpoint_key = f"path.segment.{seg_idx:02d}.endpoint_trajectory"
            current_endpoint_trajectory = (
                manifest.produced[endpoint_key][0]
                if endpoint_key in manifest.produced
                else None
            )

            irc_res = _irc_and_match(
                seg_idx=seg_idx,
                seg_dir=seg_dir,
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
                irc_step_size=irc_step_size,
                irc_never_stop=irc_never_stop,
                seg_tag=str(seg_tag),
                endpoint_trajectory=current_endpoint_trajectory,
                session=session,
                manifest=manifest,
                artifact_prefix=f"post.{seg_idx:02d}.irc",
                public_root=out_dir,
                runtime_overrides=child_runtime_overrides,
            )
            _persist_run_manifest(manifest, out_dir)

            gL = irc_res["left_min_geom"]
            gR = irc_res["right_min_geom"]
            gT = irc_res["ts_geom"]
            calculator_lease = irc_res["calculator_lease"]
            irc_plot_path = irc_res.get("irc_plot_path")
            irc_trj_path = irc_res.get("irc_trj_path")
            reverse_irc = bool(irc_res.get("reverse_irc", False))

            if isinstance(irc_plot_path, Path) and irc_plot_path.exists():
                segment_log["irc_plot"] = str(irc_plot_path)
            if isinstance(irc_trj_path, Path) and irc_trj_path.exists():
                segment_log["irc_traj"] = str(irc_trj_path)

            # IRC directional convergence for the aggregate gate.
            segment_log["irc"] = irc_res.get("irc_outcome")
            segment_log["endpoint_assignment"] = irc_res.get(
                "endpoint_assignment"
            )

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
            from pdb2reaction.io.hessian_cache import (
                discard as _hess_discard,
                load as _hess_load,
                store as _hess_store,
            )
            _left_hk  = "irc_right" if reverse_irc else "irc_left"
            _right_hk = "irc_left"  if reverse_irc else "irc_right"

            _hess_discard("irc_endpoint")
            _c = _hess_load(_left_hk)
            if _c:
                _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
            _react_opt_conv: Optional[bool] = None
            try:
                g_react_opt, _, _react_opt_conv = _optimize_endpoint_geom(
                    gL,
                    tsopt_opt_mode_default,
                    endpoint_opt_dir,
                    f"seg_{seg_idx:02d}_reactant",
                    dump=dump,
                    thresh=thresh_post,
                    calc_identity_cfg=calc_cfg_shared,
                    reject_uphill=_reject_uphill_eff,
                )
            except Exception as e:
                _echo(
                    f"[post] WARNING: Reactant endpoint optimization failed for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                g_react_opt = gL
                _react_opt_conv = None

            _hess_discard("irc_endpoint")
            _c = _hess_load(_right_hk)
            if _c:
                _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
            _prod_opt_conv: Optional[bool] = None
            try:
                g_prod_opt, _, _prod_opt_conv = _optimize_endpoint_geom(
                    gR,
                    tsopt_opt_mode_default,
                    endpoint_opt_dir,
                    f"seg_{seg_idx:02d}_product",
                    dump=dump,
                    thresh=thresh_post,
                    calc_identity_cfg=calc_cfg_shared,
                    reject_uphill=_reject_uphill_eff,
                )
            except Exception as e:
                _echo(
                    f"[post] WARNING: Product endpoint optimization failed for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                g_prod_opt = gR
                _prod_opt_conv = None

            # record endpoint-opt convergence so a nonconverged endpoint
            # (whose geometry is still used for the diagram) does not silently
            # promote its segment to a usable success.
            segment_log["endpoint_opt"] = {
                "reactant_converged": _react_opt_conv,
                "product_converged": _prod_opt_conv,
            }

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
            diag_payload = _write_public_energy_diagram(
                seg_dir / "energy_diagram_MLIP",
                labels=["R", f"TS{seg_idx}", "P"],
                energies_au=[eR, eT, eP],
                title_note="(MLIP, TSOPT + IRC)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

            tsopt_seg_energies.append((eR, eT, eP))

            segment_log["mlip"] = {
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
                    state_structs,
                    out_dir,
                    seg_idx,
                    _input_suffix,
                    manifest=manifest,
                )
                _echo(f"[all] Wrote R/TS/P for segment {seg_idx:02d} → {_seg_out}", narrative=True)
            except Exception as e:
                _echo(
                    f"[all] WARNING: Failed to copy structures for segment {seg_idx:02d}: {e}",
                    err=True,
                )

        elif do_thermo or do_dft:
            seg_model_key = f"path.mep.{seg_idx:02d}.pdb"
            seg_model_path = (
                manifest.produced[seg_model_key][0]
                if seg_model_key in manifest.produced
                else None
            )

            # Decide reference topology (if any) for freeze-atoms detection and conversion.
            freeze_ref: Optional[Path] = ref_pdb_for_seg
            if freeze_ref is None and seg_model_path is not None and seg_model_path.suffix.lower() == ".pdb":
                freeze_ref = seg_model_path
            elif freeze_ref is None and hei_model_path.suffix.lower() == ".pdb":
                freeze_ref = hei_model_path

            freeze_atoms: List[int] = _get_freeze_atoms(freeze_ref, freeze_links_flag)

            try:
                endpoint_key = f"path.segment.{seg_idx:02d}.endpoint_trajectory"
                endpoint_path = (
                    manifest.produced[endpoint_key][0]
                    if endpoint_key in manifest.produced
                    else None
                )
                if endpoint_path is None:
                    _echo(
                        f"[post] WARNING: current endpoint trajectory was not claimed for segment {seg_idx:02d}; cannot run thermo/DFT without --tsopt. Skipping segment.",
                        err=True,
                    )
                    continue
                elems, c_first, c_last = read_xyz_first_last(endpoint_path)
                endpoints = (
                    _geom_from_angstrom(elems, c_first, freeze_atoms),
                    _geom_from_angstrom(elems, c_last, freeze_atoms),
                )
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
        if do_tsopt:
            calculator_lease.release()
        for _g in (
            locals().get("gL"), locals().get("gR"), locals().get("gT"),
            locals().get("g_ts"), locals().get("g_react_opt"), locals().get("g_prod_opt"),
        ):
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        calc = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        dG_R = dG_T = dG_P = None
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
                from pdb2reaction.workflows._all_helpers import (
                    build_thermo_symmetry_provenance,
                )
                GR = _thermo_value_ha(tR, "sum_EE_and_thermal_free_energy_ha")
                GT = _thermo_value_ha(tT, "sum_EE_and_thermal_free_energy_ha")
                GP = _thermo_value_ha(tP, "sum_EE_and_thermal_free_energy_ha")
                dG_R = _thermo_value_ha(tR, "thermal_correction_free_energy_ha")
                dG_T = _thermo_value_ha(tT, "thermal_correction_free_energy_ha")
                dG_P = _thermo_value_ha(tP, "thermal_correction_free_energy_ha")
                _thermo_symmetry = build_thermo_symmetry_provenance(
                    thermo_payloads
                )
                if _thermo_symmetry:
                    segment_log["thermo_symmetry"] = _thermo_symmetry
                ts_freq_info = _read_imaginary_frequency(freq_seg_root / "TS")
                if ts_freq_info is not None:
                    segment_log["ts_imag"] = ts_freq_info
                    if ts_freq_info.get("nu_imag_max_cm") is not None:
                        segment_log["ts_imag_freq_cm"] = ts_freq_info["nu_imag_max_cm"]
                try:
                    if all(value is not None for value in (GR, GT, GP)):
                        gibbs_vals = [GR, GT, GP]
                        diag_payload = _write_public_energy_diagram(
                            seg_dir / "energy_diagram_G_MLIP",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_au=gibbs_vals,
                            title_note="(MLIP + Thermal Correction)",
                            ylabel="ΔG (kcal/mol)",
                        )
                        if diag_payload:
                            energy_diagrams.append(diag_payload)
                        g_mlip_seg_energies.append((GR, GT, GP))
                        segment_log["gibbs_mlip"] = {
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
                            "[thermo] NOTE: one or more R/TS/P free energies "
                            "are unavailable; Gibbs diagram skipped."
                        )
                except Exception as e:
                    _echo(
                        f"[thermo] WARNING: failed to build Gibbs diagram: {e}",
                        err=True,
                    )

        if do_dft:
            # Clear the calculator retained by Geometry before forking the DFT
            # subprocess so it does not inherit the GPU-pinned model.
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
                        diag_payload = _write_public_energy_diagram(
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

                if (
                    do_thermo
                    and _dft_all_ok
                    and all(value is not None for value in (dG_R, dG_T, dG_P))
                ):
                    try:
                        GR_dft_mlip = eR_dft + dG_R
                        GT_dft_mlip = eT_dft + dG_T
                        GP_dft_mlip = eP_dft + dG_P
                        if all(
                            np.isfinite([GR_dft_mlip, GT_dft_mlip, GP_dft_mlip])
                        ):
                            diag_payload = _write_public_energy_diagram(
                                seg_dir / "energy_diagram_G_DFT_plus_MLIP",
                                labels=["R", f"TS{seg_idx}", "P"],
                                energies_au=[GR_dft_mlip, GT_dft_mlip, GP_dft_mlip],
                                title_note=f"({dft_func_basis_use} // MLIP + Thermal Correction)",
                                ylabel="ΔG (kcal/mol)",
                            )
                            if diag_payload:
                                energy_diagrams.append(diag_payload)
                            g_dft_mlip_seg_energies.append(
                                (GR_dft_mlip, GT_dft_mlip, GP_dft_mlip)
                            )
                            segment_log["gibbs_dft_mlip"] = {
                                "labels": ["R", f"TS{seg_idx}", "P"],
                                "energies_au": [GR_dft_mlip, GT_dft_mlip, GP_dft_mlip],
                                "energies_kcal": [
                                    (GR_dft_mlip - GR_dft_mlip) * AU2KCALPERMOL,
                                    (GT_dft_mlip - GR_dft_mlip) * AU2KCALPERMOL,
                                    (GP_dft_mlip - GR_dft_mlip) * AU2KCALPERMOL,
                                ],
                                "diagram": str(seg_dir / "energy_diagram_G_DFT_plus_MLIP.png"),
                                "structures": state_structs,
                                "barrier_kcal": (GT_dft_mlip - GR_dft_mlip) * AU2KCALPERMOL,
                                "delta_kcal": (GP_dft_mlip - GR_dft_mlip) * AU2KCALPERMOL,
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

    if len(tsopt_seg_energies) == len(reactive):
        tsopt_all_energies = [e for triple in tsopt_seg_energies for e in triple]
        tsopt_all_labels = _build_global_segment_labels(len(tsopt_seg_energies))
        if tsopt_all_labels and len(tsopt_all_labels) == len(tsopt_all_energies):
            diag_payload = _write_public_energy_diagram(
                out_dir / "energy_diagram_MLIP_all",
                labels=tsopt_all_labels,
                energies_au=tsopt_all_energies,
                title_note="(MLIP, TSOPT + IRC; all segments)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_thermo and len(g_mlip_seg_energies) == len(reactive):
        g_mlip_all_energies = [e for triple in g_mlip_seg_energies for e in triple]
        g_mlip_all_labels = _build_global_segment_labels(len(g_mlip_seg_energies))
        if g_mlip_all_labels and len(g_mlip_all_labels) == len(g_mlip_all_energies):
            diag_payload = _write_public_energy_diagram(
                out_dir / "energy_diagram_G_MLIP_all",
                labels=g_mlip_all_labels,
                energies_au=g_mlip_all_energies,
                title_note="(MLIP + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_dft and len(dft_seg_energies) == len(reactive):
        dft_all_energies = [e for triple in dft_seg_energies for e in triple]
        dft_all_labels = _build_global_segment_labels(len(dft_seg_energies))
        if dft_all_labels and len(dft_all_labels) == len(dft_all_energies):
            diag_payload = _write_public_energy_diagram(
                out_dir / "energy_diagram_DFT_all",
                labels=dft_all_labels,
                energies_au=dft_all_energies,
                title_note=f"({dft_func_basis_use}; all segments)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_dft and do_thermo and len(g_dft_mlip_seg_energies) == len(reactive):
        g_dft_mlip_all_energies = [e for triple in g_dft_mlip_seg_energies for e in triple]
        g_dft_mlip_all_labels = _build_global_segment_labels(len(g_dft_mlip_seg_energies))
        if g_dft_mlip_all_labels and len(g_dft_mlip_all_labels) == len(g_dft_mlip_all_energies):
            diag_payload = _write_public_energy_diagram(
                out_dir / "energy_diagram_G_DFT_plus_MLIP_all",
                labels=g_dft_mlip_all_labels,
                energies_au=g_dft_mlip_all_energies,
                title_note=f"({dft_func_basis_use} // MLIP + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    # Aggregated IRC plot over all reactive segments (single trj + trj2fig)
    if irc_trj_for_all:
        _declare_public(out_dir / "irc_plot_all.png")
        _merge_irc_trajectories_to_single_plot(
            irc_trj_for_all, out_dir / "irc_plot_all.png"
        )
        _claim_public(out_dir / "irc_plot_all.png")

    # Refresh summary.json with final energy diagram metadata (including aggregated diagrams)
    summary["energy_diagrams"] = list(energy_diagrams)
    _enrich_summary(
        summary,
        version="",
        pipeline_mode="path-search" if refine_path else "path-opt",
        out_dir=out_dir,
        mlip_backend=_mlip_backend_shared,
        mlip_model=_mlip_model_shared,
        mlip_precision=_mlip_precision_shared,
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
            "path_opt_mode": opt_mode_norm,
            "post_opt_mode": tsopt_opt_mode_default,
            "ts_opt_mode": tsopt_opt_mode_default,
            "endpoint_opt_mode": tsopt_opt_mode_default,
            "mep_mode": mep_mode_kind,
            "dmf_correlated": dmf_correlated_effective,
        },
        freeze_atoms=_freeze_atoms_for_log(),
        manifest=manifest,
    )
    _publish_manifest_summary(
        path_dir / "summary.json",
        summary,
        manifest=manifest,
        key="path.summary",
        out_dir=out_dir,
    )
    _echo_detail(f"[write] Updated '{path_dir / 'summary.json'}' with energy diagrams.")
    dst_summary = out_dir / "summary.json"
    if not _copy_public_logged(
        path_dir / "summary.json",
        dst_summary,
        label="summary.json",
        echo=False,
    ):
        raise click.ClickException(
            f"Failed to publish summary.json to {dst_summary}."
        )
    _echo_detail(f"[all] Copied summary.json → {dst_summary}")
    _write_pipeline_summary_log(post_segment_logs)
    # summary.log becomes current-run owned only after the writer returns;
    # republish JSON once more so key_output_files includes that final root
    # artifact as well as every earlier producer-owned output.
    _finalize_current_summary()

    _emit_final_summary(
        out_dir,
        time_start,
        manifest,
        citation_payload=_all_method_citation_payload(),
    )


_hide_advanced_options(cli, _ALL_PRIMARY_HELP_OPTIONS)
