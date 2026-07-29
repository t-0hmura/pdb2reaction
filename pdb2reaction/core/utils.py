"""
Utilities for configuration, plotting, coordinates, Gaussian templates, and link-freezing.

Categories:
    - Generic helpers: pretty_block, format_elapsed, deep_update, load_yaml_dict, etc.
    - Gaussian (.gjf): parse_gjf_template, prepare_input_structure, convert_xyz_to_gjf
    - Plotting: build_energy_diagram (Plotly-based energy diagrams)
    - Coordinate conversion: convert_xyz_to_pdb (XYZ to PDB with reference topology)
    - Link-freezing: detect_freeze_links, resolve_freeze_atoms (for PDB link hydrogens)

Dependencies: PyYAML, ASE, Plotly, Click
"""

from __future__ import annotations

import ast
import functools
import logging
import math
import os
import re
import sys
import tempfile
import time
from collections import Counter
from collections.abc import Iterable as _Iterable, Mapping, Sequence as _Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from numbers import Real, Integral
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, List, Tuple, Callable, Iterator

import click
from click.core import ParameterSource
import numpy as np
import yaml
from ase.data import chemical_symbols
from ase.io import read
import plotly.graph_objs as go

from pdb2reaction.domain.add_elem_info import guess_element
from pdb2reaction.core.defaults import RFO_KW
from pdb2reaction.core.output import _TAG_AWARE_MARKER, emit
from pdb2reaction.io.structure_formats import (
    CIF_SUFFIXES,
    CoordinateTemplate,
    cleanup_normalized_structure,
    coordinate_template_for,
    is_cif_path,
    normalize_structure_to_pdb,
    pdb_requires_normalization,
    register_coordinate_template,
    render_mmcif_frames,
    render_pdb_coordinate_frames,
    unregister_coordinate_template,
)
from pdb2reaction.core.result_commit import commit_payloads
from pdb2reaction.io.charge import compute_charge_summary, _format_echo_message
from pysisyphus.constants import ANG2BOHR
from pysisyphus.helpers import geom_loader

logger = logging.getLogger(__name__)

# CLI verbosity state (set by the top-level --verbose callback in cli/app.py).
# `pretty_block` and any other config-echo helpers consult `is_verbose()` so
# that running without `-v` keeps stdout focused on milestones / errors /
# user-set parameters, and `-v` (any positive count) restores the full
# config dump for debugging.
#
# Two console verbosity levels drive `_patch_click_echo` (the single output
# chokepoint):
#   default (no -v)     -> NARRATIVE only: stage banners, per-stage one-line
#                          status, and the final summary. Detail lines (the bulk
#                          of the per-stage chatter) are suppressed.
#   -v / --verbose      -> full detail (every echo), plus pretty_block config
#                          dumps and DEBUG stdlib logging.
# A message is marked NARRATIVE by passing `narrative=True` to `click.echo`
# (the patched echo pops the flag before delegating); the `emit` helper below
# is the thin wrapper that does this for narrative-tagged lines.
_VERBOSE_LEVEL: int = 0

# Console gating is OFF until the real CLI entry point turns it on (the
# `-v` group callback calls `set_console_gating(True)`). While OFF the
# patched echo only shortens paths and dedupes blank lines — every message
# prints, so library/test callers that monkeypatch echo via `set_base_dir`
# but never run the CLI keep their full output. Once ON, `-v` applies
# globally; default-verbosity narrative-suppression applies only in the `all`
# pipeline (see `_PIPELINE_MODE` / child mode below).
_GATE_ACTIVE: bool = False

# Pipeline-mode state: True while the `all` command runs. Default-verbosity
# suppression of DETAIL lines is scoped to the pipeline (the `all` parent +
# its in-proc child stages) so that a standalone leaf/report command — whose
# stdout IS its deliverable (bond-summary, energy-diagram, …) — keeps full
# output at default verbosity. `-v` still applies everywhere.
_PIPELINE_MODE: bool = False

# Child-invocation state: True when this process is dispatching one of its
# subcommands through `_run_cli_main` (in-proc), so the child's group-callback
# banner / `[calc] Resolved device:` echo can be skipped to avoid repeating
# the same line 4-8x per pipeline run.
_CHILD_MODE: bool = False


def set_verbose_level(level: int) -> None:
    """Record the CLI --verbose level (0=silent .. 3=full) for gating."""
    global _VERBOSE_LEVEL
    _VERBOSE_LEVEL = max(0, min(3, int(level)))


def verbose_level() -> int:
    """Current console verbosity: 0=silent, 1=milestones, 2=default(+detail), 3=full."""
    return _VERBOSE_LEVEL


def is_verbose() -> bool:
    """True iff verbosity reaches the detail tier (level >= 2). Back-compat shim
    for call sites that gate optimizer/SCF detail (cycle tables, PySCF logger)."""
    return _VERBOSE_LEVEL >= 2


def set_console_gating(value: bool) -> None:
    """Enable/disable the narrative/detail console verbosity gate.

    Called once from the real CLI group callbacks so that only an actual
    `pdb2reaction` invocation engages default-suppression; direct library or
    test use of ``click.echo`` is left untouched.
    """
    global _GATE_ACTIVE
    _GATE_ACTIVE = bool(value)


def is_console_gating() -> bool:
    """True iff the narrative/detail console gate is currently engaged."""
    return _GATE_ACTIVE


def set_pipeline_mode(value: bool) -> None:
    """Mark that the `all` pipeline is running (scopes default-suppression)."""
    global _PIPELINE_MODE
    _PIPELINE_MODE = bool(value)


def pipeline_mode_enabled() -> bool:
    """Return the raw parent-pipeline flag without child-mode aggregation."""

    return _PIPELINE_MODE


def is_pipeline_mode() -> bool:
    """True iff the `all` pipeline (parent or child stage) is running."""
    return _PIPELINE_MODE or _CHILD_MODE


def set_child_mode(value: bool) -> None:
    """Toggle child-invocation mode for in-proc subcommand dispatch."""
    global _CHILD_MODE
    _CHILD_MODE = bool(value)


def _echo_info(msg: str, *args: Any, level: int = 2) -> None:
    # Console output mapped to the unified `-v` gate by tier:
    #   level 1 = milestone (atom counts, net/total charge, output path)
    #   level 2 = detail [default] (options, per-residue/per-resname, context)
    #   level 3 = debug (terminal-cap charge corrections, etc.)
    # The click.echo gate only fires inside the `all` pipeline (a standalone
    # leaf keeps full stdout). extract is the exception: its deliverable is the
    # pocket PDB, not stdout, so it should still tier by -v standalone. Apply
    # the same per-record gate here for the standalone case (silent at -v 0).
    if is_console_gating() and not is_pipeline_mode():
        _lvl = verbose_level()
        if _lvl <= 0:
            return
        required = 1 if level <= 1 else (2 if level == 2 else 3)
        if _lvl < required:
            return
    rendered = _format_echo_message(msg, *args)
    if is_console_gating():
        emit(rendered, narrative=(level <= 1), detail=(level == 2))
    else:
        # Programmatic extract_api callers do not install the CLI echo shim,
        # so passing private verbosity tags to vanilla Click would raise.
        click.echo(rendered)


def log_charge_summary(prefix: str,
                       summary: Dict[str, Any]):
    """
    Emit concise charge summary logs.
    """
    total = summary["total_charge"]
    protein = summary["protein_charge"]
    ligand = summary.get("ligand_total_charge", 0.0)
    ion_list: List[Tuple[str, float]] = summary.get("ion_charges", [])
    ion_total = summary.get("ion_total_charge", sum(q for _, q in ion_list))
    unk_map: Dict[str, float] = summary.get("unknown_residue_charges", {}) or {}

    if unk_map:
        items = ", ".join(f"{res}: {q:g}" for res, q in sorted(unk_map.items()))
        _echo_info("%s Per-resname ligand charges: %s", prefix, items)
    else:
        _echo_info("%s Per-resname ligand charges: (none)", prefix)

    _echo_info("%s Net protein charge: %+g", prefix, protein, level=1)
    _echo_info("%s Net ligand charge: %+g", prefix, ligand, level=1)
    if ion_list:
        _echo_info("%s Ion charges (each):", prefix)
        for tag, q in ion_list:
            _echo_info("  %s  ->  %+g", tag, q)
        _echo_info("%s Net ion charge: %+g", prefix, ion_total)
    else:
        _echo_info("%s Ion charges: (none)", prefix)
    _echo_info("%s Total active site model charge: %+g", prefix, total, level=1)
    for _rn, _rs, _dq in summary.get("terminal_corrections", []):
        _lbl = "C-terminal carboxylate" if _dq < 0 else "N-terminal ammonium"
        _echo_info("%s   %s %s %s: %+d", prefix, _lbl, _rn, _rs, _dq, level=3)


def is_child_mode() -> bool:
    """True iff we are dispatching one of our own subcommands in-process."""
    return _CHILD_MODE


def echo_run_summary(items: Dict[str, Any]) -> None:
    """Echo a compact `[key] value` run summary, then a blank line.

    Skipped in child mode (`all` already printed its own summary) and when
    items is empty. Each subcommand calls this at its entry point so a
    default-verbosity run still surfaces the input file, backend, opt mode,
    and output dir without dumping the full per-stage config block (which
    only fires under `-v`).
    """
    if is_child_mode() or not items:
        return
    for key, value in items.items():
        if value is None or value == "":
            continue
        emit(f"[{key}] {value}", narrative=True)
    emit("", narrative=False)


# YAML helpers (shared representers)


class YamlLiteralStr(str):
    """String marker to force literal block style when dumping YAML."""


class YamlFlowList(list):
    """List marker to force flow style when dumping YAML."""


def echo_resolved_device() -> None:
    """Print the resolved torch device (cuda/cpu) as a cosmetic CLI breadcrumb.

    Best-effort: silently no-ops on torch-import or CUDA-probe failure
    (the echo is informational only; a missing torch is a separate
    error path that fires elsewhere when the calculator is built).
    Replaces a try/except block that was duplicated across 10 workflow files.
    """
    # Skip in child mode: the parent `pdb2reaction all` already printed this
    # once; reprinting at every stage entry adds 4-8 identical lines per run.
    if is_child_mode():
        return
    try:
        import torch as _torch
        _resolved_dev = "cuda" if _torch.cuda.is_available() else "cpu"
    except (ImportError, AttributeError, RuntimeError):
        return
    import click as _click
    # Device breadcrumb is a level-3 (debug) line: untagged so it only shows at -v 3.
    _click.echo(f"[calc] Resolved device: {_resolved_dev}")


def optimizer_cycle_count(optimizer: Any) -> Optional[int]:
    """Return the number of optimizer cycles spent, if the optimizer exposes it."""
    cur_cycle = getattr(optimizer, "cur_cycle", None)
    if cur_cycle is None:
        return None
    try:
        return max(int(cur_cycle) + 1, 0)
    except (TypeError, ValueError):
        return None


def optimizer_terminal_status(optimizer: Any) -> str:
    """Map a pysisyphus optimizer (or a product-local runner) terminal state to
    the public status vocabulary.

    Returns ``"stalled"`` for an energy-plateau outcome,
    ``"converged"`` for a genuine stationary point, and ``"not_converged"``
    otherwise.  ``stalled`` takes precedence so a plateau is never reported as
    converged; legacy callers that only read ``is_converged`` still see
    ``False`` for a stall.
    """
    status = getattr(optimizer, "termination_status", None)
    if status in ("stalled", "converged", "not_converged"):
        return status
    if getattr(optimizer, "is_stalled", False):
        return "stalled"
    return "converged" if getattr(optimizer, "is_converged", False) else "not_converged"


def emit_optimizer_terminal_status(
    label: str,
    *,
    converged: Optional[bool],
    cycles: Optional[int],
    max_cycles: Optional[int],
    stalled: bool = False,
    stop_reason: Optional[str] = None,
) -> None:
    """Emit a consistent optimizer terminal status at detail verbosity.

    ``stalled`` renders the energy-plateau outcome and takes
    precedence over the convergence/max-cycle branches so a stalled run is
    never printed as ``Converged!``.
    """
    prefix = f"[{label}]"
    if stalled:
        if stop_reason:
            emit(
                f"{prefix} Stalled (energy plateau; not converged): {stop_reason}",
                narrative=False,
                detail=True,
            )
        else:
            emit(
                f"{prefix} Stalled (energy plateau; not converged).",
                narrative=False,
                detail=True,
            )
    elif converged is True:
        emit(f"{prefix} Converged!", narrative=False, detail=True)
    elif cycles is not None and max_cycles is not None and cycles >= max_cycles:
        emit(
            f"{prefix} Reached max cycles ({cycles}/{max_cycles}).",
            narrative=False,
            detail=True,
        )
    elif converged is False:
        if cycles is None:
            emit(f"{prefix} Stopped without convergence.", narrative=False, detail=True)
        else:
            emit(
                f"{prefix} Stopped without convergence (cycles={cycles}).",
                narrative=False,
                detail=True,
            )
    elif cycles is not None:
        emit(f"{prefix} Finished (cycles={cycles}).", narrative=False, detail=True)

    if cycles is not None:
        emit(f"{prefix} Total cycles: {cycles}", narrative=False, detail=True)


def register_yaml_representers() -> None:
    """Register shared YAML representers (literal strings and flow lists)."""
    yaml.add_representer(
        YamlLiteralStr,
        lambda dumper, data: dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|"),
        Dumper=yaml.SafeDumper
    )
    yaml.SafeDumper.add_representer(
        YamlFlowList,
        lambda dumper, data: dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)
    )


register_yaml_representers()



def pretty_block(title: str, content: Dict[str, Any], *, force: bool = False) -> str:
    """Return a YAML block with an underlined title.

    Returns an empty string below verbosity level 3 so that the default
    CLI output stays focused on milestones and user-set parameters; the
    full config dump is restored under `-v 3` for debugging.

    ``force=True`` bypasses that gate. Use it for output the user asked for
    explicitly (``--show-config``, ``--print-parsed``): a flag whose whole purpose is
    to print something must not render nothing at the default verbosity.
    """
    if not force and verbose_level() < 3:
        return ""
    if not content:
        return ""  # suppress empty blocks entirely
    if _base_dir is not None:
        content = _shorten_paths(content)
    body = yaml.safe_dump(_to_yaml_safe(content), sort_keys=False, allow_unicode=True).strip()
    return f"\n{title}\n" + "-" * len(title) + "\n" + body


# Module-level base directory for relative path display.
_base_dir: Path | None = None
_original_click_echo = None

# Raw-stdout verbosity tap for bundled optimizers whose per-cycle tables and
# summaries are written via sys.stdout (bypassing the click.echo gate). At -v 1
# only convergence verdicts are kept. At -v 2 optimizer tables are kept but
# high-volume DLC/trust-radius chatter and IPOPT metric rows that grep as
# "error" are hidden. At -v 3 raw optimizer stdout is passed through unchanged.
_PYSIS_L1_ALLOW = re.compile(
    r"^(?:Converged!|Final summary:"
    r"|max\(forces,\s*\w+\):|rms\(forces,\s*\w+\):|energy:\s"
    r"|Path with \d+ moving images\.|Number of cycles exceeded!"
    r"|Operator indicated convergence!|Insignificant coordinate change"
    r"|Energy plateau detected|Wrote final geometr)"
)
_PYSIS_L2_DENY = re.compile(
    r"^(?:\d+\s+(?:[-\d.]|nan)"          # compact table data row (cycle 0 col = nan*)
    r"|cycle\s+\S.*energy"               # compact table header
    r"|[-=]{5,}\s*$"                      # separator rule
    r"|If not specified otherwise, all quantities"
    r"|Spent\s+[\d.]+\s+s\s+preparing"
    r"|Convergence thresholds|'Superscript"
    r"|max\(\|force\|\)\s*<=|rms\(force\)\s*<="
    r"|max\(\|step\|\)\s*<=|rms\(step\)\s*<="
    r"|Rebuilt internal coordinates|Interfragment distances increased"
    r"|Dumped latest coordinates|String=|Overall NLP error\.+:)"
)
_PYSIS_V2_DENY = re.compile(
    r"^(?:Rebuilt internal coordinates|Interfragment distances increased"
    r"|Dumped latest coordinates"
    r"|Overall NLP error\.+:"
    r"|Unexpected energy increase"
    r"|Current trust radius:|Decreasing trust radius\.|Increasing trust radius\."
    r"|Keeping current trust radius|Updated trust radius:)"
)


def _pysis_stdout_visible(stripped: str) -> bool:
    """Whether a raw-stdout optimizer line is visible at the current level."""
    if _VERBOSE_LEVEL <= 0:
        return False                       # -v 0: silent
    if _VERBOSE_LEVEL >= 3:
        return True                        # -v 3: full raw optimizer output
    if _VERBOSE_LEVEL >= 2:
        return not _PYSIS_V2_DENY.match(stripped)
    if _PYSIS_L1_ALLOW.match(stripped):    # -v 1: keep the convergence verdict
        return True
    return not _PYSIS_L2_DENY.match(stripped)  # drop only the table noise


def _patch_click_echo() -> None:
    """Monkey-patch click.echo and sys.stdout to shorten absolute paths
    and suppress consecutive blank lines in output."""
    import click as _click
    import sys as _sys
    global _original_click_echo
    if _original_click_echo is not None:
        return  # already patched
    _original_click_echo = _click.echo

    _last_was_blank = [False]
    _last_visible_line = [""]
    _raw_path_echo_depth = [0]

    def _is_hessian_status_line(line: str) -> bool:
        return (
            line.startswith("[hessian]")
            or line.startswith("[HessianTiming]")
            or line.startswith("[HessianVRAM]")
        )

    def _wants_blank_before(first_visible: str, prev_visible: str) -> bool:
        if not first_visible:
            return False
        starts_hessian_status = _is_hessian_status_line(first_visible)
        return (
            first_visible.startswith("======")
            or first_visible.startswith("[time]")
            or first_visible.startswith("[stage]")
            or (
                starts_hessian_status
                and not _is_hessian_status_line(prev_visible)
            )
            or (
                first_visible.startswith("[Imaginary modes]")
                and not prev_visible.startswith("[Imaginary modes]")
            )
            or first_visible.startswith("Convergence thresholds (non mass-weighted gradient)")
            or first_visible.startswith("IRC steps exceeded.")
            or first_visible.startswith("Transition vector is mode")
            or first_visible.startswith("Wrote final geometry")
            or (
                _is_hessian_status_line(prev_visible)
                and not _is_hessian_status_line(first_visible)
            )
        )

    def _wants_blank_after(first_visible: str) -> bool:
        return first_visible.startswith("======") or first_visible.startswith("[stage]")

    def _patched_echo(message=None, **kwargs):
        # Pull the verbosity tags out before delegating (click.echo rejects
        # unknown kwargs). A line's required level: narrative -> 1 (milestone),
        # detail -> 2 (cycle tables / timing / deliverable paths), untagged -> 3
        # (config dumps / per-file path bullets / device breadcrumbs). Untagged
        # lines fall through to level 3, never to "dropped", so `-v 3`
        # reproduces every original line.
        narrative = bool(kwargs.pop("narrative", False))
        detail = bool(kwargs.pop("detail", False))
        # `force=True` bypasses the verbosity gate entirely: it tags an
        # explicitly requested machine-readable deliverable (e.g. `--json`
        # output) that must always reach stdout (except at -v 0).
        # Path shortening below still applies.
        force = bool(kwargs.pop("force", False))
        raw_path = bool(kwargs.pop("raw_path", False))
        is_err = bool(kwargs.get("err", False))
        is_blank = (message is None or (isinstance(message, str) and message.strip() == ""))
        # Verbosity gate (single chokepoint for the whole CLI), engaged only
        # once the real CLI entry point calls set_console_gating(True):
        #   - level 0 (-v 0): completely silent (even errors/forced output).
        #   - errors/warnings (err=True) and force=True: print at level >= 1.
        #   - otherwise, inside the `all` pipeline a line prints iff the current
        #     level >= its required level (narrative 1 / detail 2 / untagged 3).
        #     Blank lines pass through (spacing). A standalone leaf/report
        #     command keeps full output (its stdout is the deliverable).
        if _GATE_ACTIVE:
            if _VERBOSE_LEVEL <= 0:
                return  # -v 0: completely silent
            if (
                not is_err
                and not force
                and not is_blank
                and is_pipeline_mode()
            ):
                required = 1 if narrative else (2 if detail else 3)
                if _VERBOSE_LEVEL < required:
                    return
        if not raw_path and message is not None and _base_dir is not None and isinstance(message, str):
            bd = str(_base_dir)
            if bd in message:
                message = message.replace(bd + "/", "").replace(bd, ".")
        if isinstance(message, str):
            first_visible = next((line.strip() for line in message.splitlines() if line.strip()), "")
            prev_visible = _last_visible_line[0]
            if not _last_was_blank[0] and not message.startswith("\n") and _wants_blank_before(first_visible, prev_visible):
                message = "\n" + message
            if first_visible and _wants_blank_after(first_visible) and not message.endswith("\n"):
                message += "\n"
        # Suppress consecutive blank lines
        if isinstance(message, str) and _last_was_blank[0] and message.startswith("\n"):
            message = message.lstrip("\n")
        if is_blank and _last_was_blank[0]:
            return
        ends_with_nl = isinstance(message, str) and message.endswith("\n")
        _last_was_blank[0] = is_blank or ends_with_nl
        if raw_path:
            _raw_path_echo_depth[0] += 1
        try:
            _original_click_echo(message, **kwargs)
            if isinstance(message, str):
                for line in reversed(message.splitlines()):
                    if line.strip():
                        _last_visible_line[0] = line.strip()
                        break
        finally:
            if raw_path:
                _raw_path_echo_depth[0] -= 1

    setattr(_patched_echo, _TAG_AWARE_MARKER, True)
    _click.echo = _patched_echo

    # Wrap sys.stdout to suppress consecutive blank lines and shorten paths
    _real_stdout = _sys.stdout

    class _FilteredStdout:
        def __init__(self, stream):
            self._stream = stream
            self._last_was_blank = False
            self._current_line_has_content = False
            self._has_visible_output = False
            self._suppress_next_nl = False

        @staticmethod
        def _starts_raw_section(stripped: str) -> bool:
            return (
                stripped.startswith("Spent ")
                or stripped.startswith("Path with ")
                or stripped.startswith("Final summary:")
                or stripped.startswith("Convergence thresholds")
                or stripped.startswith("IRC steps exceeded.")
                or stripped.startswith("Transition vector is mode")
                or stripped.startswith("Wrote final geometry")
                or stripped.startswith("======")
            )

        def write(self, s):
            original_len = len(s) if isinstance(s, str) else None
            # Verbosity tap for raw optimizer stdout (pysisyphus prints the
            # per-cycle table + Converged! via sys.stdout, bypassing the
            # click.echo gate). Only active under the CLI gate; library/test
            # callers keep full output. A suppressed content line also swallows
            # its trailing newline so no blank gap is left behind.
            if _GATE_ACTIVE and isinstance(s, str):
                stripped = s.strip()
                if not stripped:
                    # blank/whitespace: swallow entirely at -v 0, or when it
                    # trails a line we just suppressed (avoid a blank gap).
                    if _VERBOSE_LEVEL <= 0 or self._suppress_next_nl:
                        self._suppress_next_nl = False
                        return len(s)
                elif not _pysis_stdout_visible(stripped):
                    self._suppress_next_nl = True
                    return len(s)
                else:
                    self._suppress_next_nl = False
                    prev_visible = _last_visible_line[0]
                    visible_before = self._has_visible_output or bool(prev_visible)
                    needs_blank = (
                        self._starts_raw_section(stripped)
                        or _wants_blank_before(stripped, prev_visible)
                    )
                    if (
                        visible_before
                        and needs_blank
                        and not self._last_was_blank
                        and not _last_was_blank[0]
                    ):
                        self._stream.write("\n")
                        self._last_was_blank = True
                        _last_was_blank[0] = True
                        self._current_line_has_content = False
            if _raw_path_echo_depth[0] <= 0 and _base_dir is not None and isinstance(s, str):
                bd = str(_base_dir)
                if bd in s:
                    s = s.replace(bd + "/", "").replace(bd, ".")
            if isinstance(s, str):
                if self._last_was_blank and s.startswith("\n"):
                    s = s.lstrip("\n")
                    if not s:
                        return original_len
                s = re.sub(r"\n{3,}", "\n\n", s)
            # Suppress consecutive blank lines
            if s == "\n":
                if self._current_line_has_content:
                    self._current_line_has_content = False
                    self._last_was_blank = False
                    self._has_visible_output = True
                    _last_was_blank[0] = False
                    return self._stream.write(s)
                if self._last_was_blank:
                    _last_was_blank[0] = True
                    return len(s)
                self._last_was_blank = True
            elif s.strip() == "":
                pass  # whitespace-only, don't reset
            else:
                self._current_line_has_content = True
                self._last_was_blank = False
                self._has_visible_output = True
                if s.endswith("\n\n"):
                    self._current_line_has_content = False
                    self._last_was_blank = True
                elif s.endswith("\n"):
                    self._current_line_has_content = False
            _last_was_blank[0] = self._last_was_blank
            if isinstance(s, str):
                for line in reversed(s.splitlines()):
                    if line.strip():
                        _last_visible_line[0] = line.strip()
                        break
            return self._stream.write(s)

        def flush(self):
            self._stream.flush()

        def __getattr__(self, name):
            return getattr(self._stream, name)

    _sys.stdout = _FilteredStdout(_real_stdout)


def set_base_dir(path: Path | str | None) -> None:
    """Set the base directory for relative path display.

    Also monkey-patches ``click.echo`` so that any absolute path under
    *base_dir* is automatically shortened to a relative path in all
    CLI output.
    """
    global _base_dir
    _base_dir = Path(path).resolve() if path else None
    _patch_click_echo()


def rel_display(path: Path | str) -> str:
    """Return a display string for *path*, relative to the base dir when possible."""
    p = Path(path)
    if _base_dir is not None:
        try:
            return str(p.resolve().relative_to(_base_dir))
        except ValueError:
            pass
    return str(p)


def _shorten_paths(content: Dict[str, Any]) -> Dict[str, Any]:
    """Replace absolute path strings with relative paths in a config dict."""
    out: Dict[str, Any] = {}
    for k, v in content.items():
        if isinstance(v, str) and v.startswith("/") and ("/" in v[1:]):
            out[k] = rel_display(v)
        else:
            out[k] = v
    return out


def _to_yaml_safe(value: Any) -> Any:
    """Recursively convert NumPy values/containers into YAML-safe builtins."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_to_yaml_safe(v) for v in value.tolist()]
    if isinstance(value, Mapping):
        out: Dict[Any, Any] = {}
        for k, v in value.items():
            nk = _to_yaml_safe(k)
            if isinstance(nk, (list, tuple, set, dict)):
                nk = str(nk)
            out[nk] = _to_yaml_safe(v)
        return out
    if isinstance(value, tuple):
        return [_to_yaml_safe(v) for v in value]
    if isinstance(value, YamlFlowList):
        return YamlFlowList([_to_yaml_safe(v) for v in value])
    if isinstance(value, list):
        return [_to_yaml_safe(v) for v in value]
    if isinstance(value, set):
        return [_to_yaml_safe(v) for v in sorted(value, key=lambda x: str(x))]
    return value


def strip_inherited_keys(
    child_cfg: Dict[str, Any],
    base_cfg: Mapping[str, Any],
    *,
    mode: str = "present",
) -> Dict[str, Any]:
    """Return child_cfg without inherited keys (for concise logs)."""
    trimmed: Dict[str, Any] = {}
    if mode not in {"present", "same"}:
        raise ValueError(f"Unknown strip_inherited_keys mode: {mode}")
    for key, value in child_cfg.items():
        if key in base_cfg:
            if mode == "present":
                continue
            if base_cfg.get(key) == value:
                continue
        trimmed[key] = value
    return trimmed


def format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize geometry configuration for CLI echo output.
    """
    g = dict(geom_cfg)
    freeze_atoms = g.get("freeze_atoms")
    if freeze_atoms is None:
        return g

    if isinstance(freeze_atoms, str):
        return g

    try:
        items = list(freeze_atoms)
    except TypeError:
        return g

    # Display as 1-based (internal is 0-based)
    items_1based = [i + 1 for i in items]
    joined = ",".join(map(str, items_1based))
    g["freeze_atoms"] = f"[{joined}]" if items_1based else "[]"
    return g


def format_elapsed(prefix: str, start_time: float, end_time: Optional[float] = None) -> str:
    """Return a formatted elapsed-time string with the provided ``prefix`` label."""
    elapsed = max(0.0, (end_time if end_time is not None else time.perf_counter()) - start_time)
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{prefix}: {int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"


def xyz_string_with_energy(geom: Any, energy: Optional[float] = None) -> str:
    """Return an XYZ string, optionally overwriting the comment line with an energy value."""
    s = geom.as_xyz()
    lines = s.splitlines()
    if energy is not None and len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{energy:.12f}"
        s = "\n".join(lines)
    if not s.endswith("\n"):
        s += "\n"
    return s


def distance_A_from_coords(coords_bohr: "np.ndarray", i: int, j: int) -> float:
    """Return interatomic distance in Å given coords in Bohr."""
    diff = coords_bohr[i] - coords_bohr[j]
    return float(np.linalg.norm(diff) / ANG2BOHR)


def distance_tag(value_A: float, *, digits: int = 2, pad: int = 3) -> str:
    """Format a distance in Å as a zero-padded integer tag (default: ×10^2)."""
    scale = 10 ** digits
    return f"{int(round(value_A * scale)):0{pad}d}"


def as_list(raw: Any) -> List[Any]:
    """Return ``raw`` as a list, or [] when not iterable/None."""
    if raw is None:
        return []
    try:
        return list(raw)
    except Exception as exc:
        logger.debug("as_list: cannot convert %r to list: %s", type(raw).__name__, exc)
        return []


def ensure_dir(path: Path) -> None:
    """Create a directory (parents ok); noop if it already exists."""
    path.mkdir(parents=True, exist_ok=True)


def collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
    """
    Collect a variadic option's values in raw occurrence order.

    Supports grouped/repeated spellings plus Click's ``--long=value`` and
    attached-short ``-ivalue`` forms.  This is used only by the compatibility
    commands whose public grammar historically allowed ``-i A B C``.
    """
    vals: List[str] = []
    names_set = set(names)
    long_names = tuple(name for name in names_set if name.startswith("--"))
    short_names = tuple(
        name for name in names_set if name.startswith("-") and not name.startswith("--")
    )
    i = 0
    while i < len(argv):
        tok = argv[i]
        matched = tok in names_set
        inline: Optional[str] = None
        if not matched:
            for name in long_names:
                if tok.startswith(name + "="):
                    matched = True
                    inline = tok.split("=", 1)[1]
                    break
        if not matched:
            for name in short_names:
                if tok.startswith(name) and tok != name:
                    matched = True
                    inline = tok[len(name):]
                    break
        if not matched:
            i += 1
            continue
        if inline is not None:
            vals.append(inline)
        i += 1
        while i < len(argv) and not argv[i].startswith("-"):
            vals.append(argv[i])
            i += 1
    return vals


def current_cli_args(ctx: Optional[click.Context] = None) -> List[str]:
    """Return normalized arguments for the current top-level CLI invocation."""
    if ctx is None:
        ctx = click.get_current_context(silent=True)
    if ctx is not None:
        recorded = ctx.meta.get("pdb2reaction.cli.raw_args")
        if recorded is not None:
            return [str(value) for value in recorded]
    return list(sys.argv[1:])


def reject_option_like_extra_args(
    extra_args: Sequence[str],
    *,
    allowed_options: Sequence[str] = (),
    allowed_values: Sequence[str] = (),
    consumed_values: Sequence[Any] = (),
) -> None:
    """Reject residual tokens not claimed by a legacy variadic option.

    A few public commands accept the historical ``-i A B`` spelling by
    enabling Click's ``allow_extra_args`` mode.  That mode must not turn an
    option typo or unrelated bare token into a successful no-op, so only
    values recovered from a declared variadic option are allowed through.
    """
    allowed = frozenset(str(value) for value in allowed_options)
    remaining = Counter(str(value) for value in allowed_values)
    for raw in consumed_values:
        value = str(raw)
        if remaining[value] > 0:
            remaining[value] -= 1
    for raw in extra_args:
        value = str(raw)
        if remaining[value] > 0:
            remaining[value] -= 1
            continue
        if value.startswith("-") and value not in allowed:
            raise click.UsageError(f"No such option: {value}")
        if value not in allowed:
            raise click.UsageError(f"Unexpected extra argument: {value}")


def collect_single_option_values(
    argv: Sequence[str],
    names: Sequence[str],
    label: str,
) -> List[str]:
    """Collect values following a flag that must appear at most once."""
    vals: List[str] = []
    seen = 0
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            seen += 1
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    if seen > 1:
        raise click.BadParameter(
            f"Use a single {label} followed by multiple values; repeated flags are not accepted."
        )
    return vals


def geom_from_xyz_string(
    xyz_text: str,
    *,
    coord_type: str,
    freeze_atoms: Optional[Sequence[int]] = None,
) -> Any:
    """Load a pysisyphus Geometry from an XYZ text string (tempfile-backed)."""
    s = xyz_text if xyz_text.endswith("\n") else (xyz_text + "\n")
    freeze_atoms = list(freeze_atoms) if freeze_atoms is not None else []
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()

        g = geom_loader(
            Path(tmp.name),
            coord_type=coord_type,
            freeze_atoms=freeze_atoms,
        )
        try:
            g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        except Exception:
            click.echo(
                "[geom] WARNING: Failed to attach freeze_atoms to geometry.",
                err=True,
            )
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def snapshot_geometry(geom: Any, *, coord_type_default: str) -> Any:
    """Create an independent pysisyphus Geometry snapshot from the given Geometry."""
    s = geom.as_xyz()
    return geom_from_xyz_string(
        s,
        coord_type=getattr(geom, "coord_type", coord_type_default),
        freeze_atoms=getattr(geom, "freeze_atoms", []),
    )


def make_snapshot_geometry(coord_type_default: str) -> Callable[[Any], Any]:
    """Return a snapshot helper bound to a default coord_type (scan helpers)."""
    return functools.partial(snapshot_geometry, coord_type_default=coord_type_default)


def normalize_freeze_atoms(raw: Any) -> List[int]:
    """Normalize freeze_atoms values (string/list/iterable) into a list of integers."""
    if raw is None:
        return []
    if isinstance(raw, str):
        tokens = re.findall(r"-?\d+", raw)
        return [int(tok) for tok in tokens]
    try:
        items = list(raw)
    except TypeError as exc:
        logger.debug("normalize_freeze_atoms: %r is not a sequence of indices: %s", raw, exc)
        return []
    out: List[int] = []
    for item in items:
        try:
            out.append(int(item))
        except (TypeError, ValueError) as exc:
            # An unparsable entry must not degrade to an empty list: that let
            # `geom.freeze_atoms: [1, 2, three]` run with NOTHING frozen, silently turning a
            # hard-freeze / PHVA result into an unconstrained one.
            raise ValueError(
                f"freeze_atoms: cannot interpret {item!r} as an atom index"
            ) from exc
    return out


def merge_freeze_atom_indices(
    geom_cfg: Dict[str, Any],
    *indices: _Iterable[int],
) -> List[int]:
    """Merge one or more iterables of indices into ``geom_cfg['freeze_atoms']``.

    Existing entries are preserved, duplicates removed, and the result sorted.
    The updated list is returned.
    """
    merged: set[int] = set()
    base = geom_cfg.get("freeze_atoms", None)
    merged.update(normalize_freeze_atoms(base))
    for seq in indices:
        merged.update(normalize_freeze_atoms(seq))
    result = sorted(merged)
    geom_cfg["freeze_atoms"] = result
    return result


def merge_freeze_atom_groups(*groups: Sequence[int]) -> List[int]:
    """Merge multiple freeze_atoms groups into a sorted list of ints."""
    merged: set[int] = set()
    for group in groups:
        merged.update(normalize_freeze_atoms(group))
    return sorted(merged)


def yaml_freeze_to_internal(indices: List[int]) -> List[int]:
    """Convert 1-based YAML freeze_atoms to 0-based internal indices."""
    return sorted(i - 1 for i in indices if i > 0)


def _parse_freeze_atoms(arg: Optional[str]) -> List[int]:
    """Parse comma-separated 1-based indices (e.g., "1,3,5") into 0-based ints."""
    if arg is None:
        return []

    items = [chunk.strip() for chunk in str(arg).split(",")]
    indices: List[int] = []
    for idx, chunk in enumerate(items, start=1):
        if not chunk:
            continue
        try:
            value = int(chunk)
        except ValueError as exc:
            raise click.BadParameter(
                f"Invalid integer in --freeze-atoms entry #{idx}: '{chunk}'"
            ) from exc
        if value <= 0:
            raise click.BadParameter(
                f"--freeze-atoms expects 1-based positive indices; got {value}"
            )
        indices.append(value - 1)
    return sorted(set(indices))


def build_sopt_kwargs(
    kind: str,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    max_step_bohr: float,
    relax_max_cycles: Optional[int] = None,
    relax_override_requested: bool = False,
    out_dir: Optional[Path] = None,
    prefix: Optional[str] = None,
) -> Dict[str, Any]:
    """Build LBFGS/RFO kwargs from one already-resolved scan config.

    ``relax_max_cycles`` and ``relax_override_requested`` retain the previous
    callable signature for embedded users. Scan workflows resolve that value
    into ``opt_cfg`` once and omit the compatibility arguments.
    """
    if out_dir is None or prefix is None:
        raise TypeError("out_dir and prefix are required")
    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    if kind == "lbfgs":
        args = {**lbfgs_cfg, **common}
        args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
    else:
        args = {**rfo_cfg, **common}
        tr = float(rfo_cfg.get("trust_radius", RFO_KW["trust_radius"]))
        args["trust_radius"] = min(tr, max_step_bohr)
        args["trust_max"] = min(float(rfo_cfg.get("trust_max", RFO_KW["trust_max"])), max_step_bohr)
    if relax_override_requested:
        if relax_max_cycles is None:
            raise TypeError(
                "relax_max_cycles is required when relax_override_requested is true"
            )
        args["max_cycles"] = int(relax_max_cycles)
    return args


def make_sopt_optimizer(
    geom: Any,
    kind: str,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    max_step_bohr: float,
    relax_max_cycles: Optional[int] = None,
    relax_override_requested: bool = False,
    out_dir: Optional[Path] = None,
    prefix: Optional[str] = None,
):
    """Construct an LBFGS/RFO optimizer based on shared settings."""
    args = build_sopt_kwargs(
        kind,
        lbfgs_cfg,
        rfo_cfg,
        opt_cfg,
        max_step_bohr,
        relax_max_cycles,
        relax_override_requested,
        out_dir,
        prefix,
    )
    from pysisyphus.optimizers.LBFGS import LBFGS
    from pysisyphus.optimizers.RFOptimizer import RFOptimizer

    if kind == "lbfgs":
        return LBFGS(geom, **args)
    return RFOptimizer(geom, **args)


def normalize_choice(
    value: str,
    *,
    param: str,
    alias_groups: _Sequence[tuple[_Sequence[str], str]],
    allowed_hint: str,
) -> str:
    """Normalize *value* using alias groups and raise ``click.BadParameter`` on failure."""
    key = (value or "").strip().lower()
    for aliases, canonical in alias_groups:
        if any(key == alias.lower() for alias in aliases):
            return canonical

    hint = allowed_hint.strip()
    detail = f" Allowed: {hint}." if hint else ""
    raise click.BadParameter(f"Unknown value for {param} '{value}'.{detail}")


def cli_param_overridden(ctx: click.Context, name: str) -> bool:
    """Return True when a CLI parameter value was explicitly provided."""
    try:
        source = ctx.get_parameter_source(name)
    except Exception as exc:
        # On a failed source query the param name is unknown to the context, so
        # treat it as NOT explicitly provided: falling through to --config
        # YAML / defaults is safer than letting a Click DEFAULT value outrank an
        # explicit YAML entry.
        logger.debug("cli_param_overridden: failed to query source for %r: %s", name, exc)
        return False
    return source not in (None, ParameterSource.DEFAULT)


def deep_update(dst: Dict[str, Any], src: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Recursively update mapping *dst* with *src*, returning *dst*.
    """
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _get_mapping_section(cfg: Mapping[str, Any], path: _Sequence[str]) -> Optional[Dict[str, Any]]:
    cur: Any = cfg
    for key in path:
        if not isinstance(cur, Mapping):
            return None
        cur = cur.get(key)
        if cur is None:
            return None
    return cur if isinstance(cur, dict) else None


def apply_yaml_overrides(
    yaml_cfg: Mapping[str, Any],
    overrides: _Sequence[Tuple[Dict[str, Any], _Sequence[_Sequence[str]]]],
) -> None:
    """Apply YAML overrides to multiple target dictionaries.

    Parameters
    ----------
    yaml_cfg : Mapping[str, Any]
        Parsed YAML configuration (root-level mapping).
    overrides : Sequence[Tuple[Dict[str, Any], Sequence[Sequence[str]]]]
        Each entry consists of the target dictionary to update followed by one or
        more candidate key paths. The first existing path is used. For example::

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (lbfgs_cfg, (("opt", "lbfgs"), ("lbfgs",))),
                ],
            )

        This mirrors the previous ``deep_update(..., yaml_cfg.get(...))`` pattern
        while centralizing the shared logic.
    """
    for target, paths in overrides:
        for path in paths:
            norm_path = tuple(path)
            section = _get_mapping_section(yaml_cfg, norm_path)
            if section is not None:
                deep_update(target, section)
                break


def load_yaml_dict(path: Optional[Path]) -> Dict[str, Any]:
    """
    Load a YAML file whose root must be a mapping. Return an empty dict if *path* is None.
    """
    if not path:
        return {}

    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")

    return data


def build_scan_configs(
    yaml_cfg: Mapping[str, Any],
    *,
    geom_kw: Dict[str, Any],
    calc_kw: Dict[str, Any],
    opt_kw: Dict[str, Any],
    lbfgs_kw: Dict[str, Any],
    rfo_kw: Dict[str, Any],
    bias_kw: Dict[str, Any],
    extra_overrides: Sequence[Tuple[Dict[str, Any], _Sequence[_Sequence[str]]]] = (),
    charge: Optional[int] = None,
    spin: Optional[int] = None,
    workers: int = 1,
    workers_per_node: int = 1,
    out_dir: str = ".",
    thresh: Optional[str] = None,
    bias_k: Optional[float] = None,
    relax_max_cycles: Optional[int] = None,
    relax_max_cycles_overridden: bool = False,
    set_charge_spin: bool = True,
    workers_overridden: bool = True,
    workers_per_node_overridden: bool = True,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Build common scan configs with uniform precedence ``defaults < --config (YAML) < CLI``.

    The YAML configuration is applied first (over the built-in defaults), then
    explicit CLI-derived values override it, matching the precedence used by the
    other subcommands (e.g. ``opt``).  ``workers`` / ``workers_per_node`` carry a
    non-``None`` CLI default, so they are only treated as a CLI override when the
    caller signals it via ``workers_overridden`` / ``workers_per_node_overridden``
    (typically ``cli_param_overridden(ctx, "workers")``). ``thresh`` / ``bias_k``
    default to ``None`` and are self-gating. ``relax_max_cycles`` is applied only
    when its source flag says the user supplied it explicitly.  The returned
    dictionaries are therefore the sole effective configuration consumed by
    optimizer construction and config echoing.
    """
    geom_cfg = dict(geom_kw)
    calc_cfg = dict(calc_kw)
    opt_cfg = dict(opt_kw)
    lbfgs_cfg = dict(lbfgs_kw)
    rfo_cfg = dict(rfo_kw)
    bias_cfg = dict(bias_kw)

    # 1. YAML (--config) over built-in defaults.
    apply_yaml_overrides(
        yaml_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (opt_cfg, (("opt",),)),
            (lbfgs_cfg, (("lbfgs",),)),
            (rfo_cfg, (("rfo",),)),
            (bias_cfg, (("bias",),)),
            *list(extra_overrides),
        ],
    )

    # 2. Explicit CLI-derived values over YAML (uniform with other subcommands).
    if set_charge_spin:
        if charge is not None:
            calc_cfg["charge"] = int(charge)
        if spin is not None:
            calc_cfg["spin"] = int(spin)
    if workers_overridden:
        calc_cfg["workers"] = int(workers)
    if workers_per_node_overridden:
        calc_cfg["workers_per_node"] = int(workers_per_node)
    # out_dir / dump are run-scoped (not YAML-tunable) and always set.
    opt_cfg["out_dir"] = out_dir
    opt_cfg["dump"] = False
    if thresh is not None:
        opt_cfg["thresh"] = str(thresh)
    if relax_max_cycles_overridden and relax_max_cycles is not None:
        opt_cfg["max_cycles"] = int(relax_max_cycles)
    if bias_k is not None:
        bias_cfg["k"] = float(bias_k)

    return geom_cfg, calc_cfg, opt_cfg, lbfgs_cfg, rfo_cfg, bias_cfg



# Plotly: Energy diagram builder
def build_energy_diagram(
    energies: Sequence[float],
    labels: Sequence[str],
    ylabel: str = "ΔE",
    baseline: bool = False,
    showgrid: bool = False,
) -> go.Figure:
    """
    Plot an energy diagram using Plotly.

    Parameters
    ----------
    energies : Sequence[float]
        Energies for each state (same unit). Values are plotted without conversion.
    labels : Sequence[str]
        Labels corresponding to each state (for example, ["R", "TS1", "IM1", "TS2", "P"]).
        Must be the same length as ``energies``.
    ylabel : str, optional
        Y-axis label (for example, "ΔE" or "ΔG"). Defaults to ``"ΔE"``.
    baseline : bool, optional
        If ``True``, draw a dotted baseline at the energy of the first state across the plot.
    showgrid : bool, optional
        If ``True``, show grid lines on both axes. Defaults to ``False``.

    Returns
    -------
    plotly.graph_objs.Figure
        Figure containing the energy diagram.

    Notes
    -----
    - Each state is rendered as a thick horizontal segment (width ``HLINE_WIDTH``).
    - Adjacent states are connected by dotted diagonal segments from the right end of
      the left state to the left end of the right state.
    - Segment length automatically shrinks with additional states so that gaps remain
      between neighbors.
    - X-axis ticks are centered on each state and labeled using ``labels``.
    """
    if len(energies) == 0:
        raise ValueError("`energies` must contain at least one value.")
    if len(energies) != len(labels):
        raise ValueError("`energies` and `labels` must have the same length.")

    n = len(energies)
    energies = [float(e) for e in energies]

    AXIS_WIDTH = 3
    FONT_SIZE = 18
    AXIS_TITLE_SIZE = 20
    HLINE_WIDTH = 6           # Width of the horizontal state segments
    CONNECTOR_WIDTH = 2       # Width of the dotted connectors
    LINE_COLOR = "#1C1C1C"
    GRID_COLOR = "lightgrey"

    # Geometry along the X axis (centers and segment lengths)
    # Place segment centers at 0.5, 1.5, 2.5, ... (equally spaced)
    centers = [i + 0.5 for i in range(n)]

    # Shorten the segment as n grows (min 0.35, max 0.85)
    # Examples: n=5 -> 0.7, n=10 -> 0.5, n>=20 -> 0.35
    seg_width = min(0.85, max(0.35, 0.90 - 0.04 * n))
    half = seg_width / 2.0

    lefts = [c - half for c in centers]
    rights = [c + half for c in centers]

    fig = go.Figure()

    # Baseline (dotted line at the first energy level)
    if baseline:
        fig.add_trace(
            go.Scatter(
                x=[lefts[0], rights[-1]],
                y=[energies[0], energies[0]],
                mode="lines",
                line=dict(color=GRID_COLOR, dash="dot", width=2),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Horizontal segments for each state
    for i, (e, lab) in enumerate(zip(energies, labels)):
        fig.add_trace(
            go.Scatter(
                x=[lefts[i], rights[i]],
                y=[e, e],
                mode="lines",
                line=dict(color=LINE_COLOR, width=HLINE_WIDTH),
                hovertemplate=f"{lab}: %{{y:.6f}}<extra></extra>",
                showlegend=False,
            )
        )

    # Dotted diagonals between adjacent states (right end -> left end)
    for i in range(n - 1):
        fig.add_trace(
            go.Scatter(
                x=[rights[i], lefts[i + 1]],
                y=[energies[i], energies[i + 1]],
                mode="lines",
                line=dict(color=LINE_COLOR, width=CONNECTOR_WIDTH, dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Add a small margin beyond the first/last segments on X
    xpad = max(0.08, 0.15 * (1.0 - seg_width))
    x_min = lefts[0] - xpad
    x_max = rights[-1] + xpad

    # Add vertical padding above and below
    y_min = min(energies)
    y_max = max(energies)
    span = max(1e-6, y_max - y_min)  # Avoid zero span even if all values match
    ypad_low = 0.10 * span
    ypad_high = 0.20 * span
    y_range = [y_min - ypad_low, y_max + ypad_high]

    xaxis_config = dict(
        range=[x_min, x_max],
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        tickmode="array",
        tickvals=centers,
        ticktext=list(labels),
        title=dict(text="", font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    yaxis_config = dict(
        range=y_range,
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        title=dict(text=ylabel, font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    fig.update_layout(
        xaxis=xaxis_config,
        yaxis=yaxis_config,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=80, r=40, t=40, b=80),
    )

    return fig


def convert_xyz_to_pdb(xyz_path: Path, ref_pdb_path: Path, out_pdb_path: Path) -> None:
    """Overlay coordinates from *xyz_path* onto the topology of *ref_pdb_path* and write to *out_pdb_path*.

    The reference PDB is used as a text template: only the coordinate columns
    (31–54) of ATOM/HETATM records are replaced with coordinates from the XYZ
    frames.  All other PDB metadata (atom names, residue info, element columns,
    chain IDs, B-factors, etc.) are preserved verbatim, avoiding element
    misidentification bugs in external PDB parsers (e.g., ASE reading ``ZN``
    atom names as nitrogen).

    Notes
    -----
        - *xyz_path* may contain one or many frames. For multi-frame trajectories,
          MODEL/ENDMDL blocks are written for each frame.
        - The complete trajectory is validated before the destination is atomically replaced.
    """
    from ase.io import read as ase_read

    traj = ase_read(str(xyz_path), index=":", format="xyz")
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")
    symbols = [frame.get_chemical_symbols() for frame in traj]
    frames = [np.asarray(frame.get_positions(), dtype=float) for frame in traj]
    pdb_payload = render_pdb_coordinate_frames(ref_pdb_path, symbols, frames).encode(
        "utf-8"
    )

    # A CIF/oversized-PDB input is represented by a registered internal PDB.
    # Render its companion before publishing either path, so a validation or
    # serialization failure leaves every existing artifact and registry entry
    # untouched.
    template = coordinate_template_for(ref_pdb_path)
    cif_path = Path(out_pdb_path).with_suffix(".cif")
    cif_payload: Optional[bytes] = None
    if template is not None:
        expected = tuple(record.element.title() for record in template.records)
        for frame_index, frame_symbols in enumerate(symbols, start=1):
            actual = tuple(str(symbol).title() for symbol in frame_symbols)
            if actual != expected:
                raise ValueError(
                    f"Ordered elements differ from the retained coordinate template "
                    f"in XYZ frame {frame_index}."
                )
        cif_payload = render_mmcif_frames(frames, template).encode("utf-8")

    output_payloads = {Path(out_pdb_path): pdb_payload}
    if cif_payload is not None:
        output_payloads[cif_path] = cif_payload
    commit_payloads(Path(out_pdb_path), output_payloads)
    if template is not None and cif_payload is not None:
        register_coordinate_template(out_pdb_path, template)
    else:
        unregister_coordinate_template(out_pdb_path)


# Gaussian input (.gjf) helpers

_FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
_GJF_COORD_RE = re.compile(
    rf"^(?P<prefix>.*?)(?P<sep0>\s*)(?P<x>{_FLOAT_PATTERN})(?P<sep1>\s+)"
    rf"(?P<y>{_FLOAT_PATTERN})(?P<sep2>\s+)(?P<z>{_FLOAT_PATTERN})(?P<suffix>\s*)$"
)


def _count_decimals(token: str) -> int:
    if "." not in token:
        return 6
    mantissa = token.split(".", 1)[1]
    mantissa = mantissa.split("e", 1)[0].split("E", 1)[0]
    mantissa = mantissa.split("d", 1)[0].split("D", 1)[0]
    return sum(ch.isdigit() for ch in mantissa)


def _format_like(template: str, value: float) -> str:
    stripped = template.strip()
    width = len(template)
    if not stripped:
        formatted = f"{value:.10f}"
    elif any(ch in stripped for ch in "eEdD"):
        mantissa = stripped.replace("d", "e").replace("D", "E")
        mantissa_only, _, _ = mantissa.partition("E")
        decimals = _count_decimals(mantissa_only)
        formatted = f"{value:.{decimals}E}"
        if "d" in template:
            formatted = formatted.replace("E", "d")
        elif "D" in template:
            formatted = formatted.replace("E", "D")
        elif "e" in template:
            formatted = formatted.replace("E", "e")
    else:
        decimals = _count_decimals(stripped)
        formatted = f"{value:.{decimals}f}"
    if len(formatted) < width:
        formatted = formatted.rjust(width)
    return formatted


def _token_to_symbol(text: str) -> str:
    stripped = text.strip()
    if not stripped:
        raise ValueError("Empty atom token in Gaussian coordinate line.")
    token = stripped.split()[0]
    m = re.match(r"([A-Za-z]{1,2})", token)
    if m:
        return m.group(1).title()
    if token.isdigit():
        z = int(token)
        if 0 <= z < len(chemical_symbols):
            return chemical_symbols[z]
    raise ValueError(f"Could not determine element symbol from '{token}'.")


def _convert_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


@dataclass
class GjfCoordinateLine:
    prefix: str
    sep0: str
    sep1: str
    sep2: str
    suffix: str
    x_template: str
    y_template: str
    z_template: str
    symbol: str
    x: float
    y: float
    z: float

    def render(self, coords: Tuple[float, float, float]) -> str:
        x_str = _format_like(self.x_template, coords[0])
        y_str = _format_like(self.y_template, coords[1])
        z_str = _format_like(self.z_template, coords[2])
        return f"{self.prefix}{self.sep0}{x_str}{self.sep1}{y_str}{self.sep2}{z_str}{self.suffix}"


@dataclass
class GjfTemplate:
    path: Path
    prefix_lines: List[str]
    suffix_lines: List[str]
    coord_lines: List[GjfCoordinateLine]
    charge: int
    spin: int

    @property
    def natoms(self) -> int:
        return len(self.coord_lines)

    def as_xyz_string(self) -> str:
        lines = [str(self.natoms), f"converted from {self.path.name}"]
        for atom in self.coord_lines:
            lines.append(f"{atom.symbol}  {atom.x:.10f}  {atom.y:.10f}  {atom.z:.10f}")
        return "\n".join(lines) + "\n"


def write_xyz_trj_with_energy(images: Sequence[Any], energies: Sequence[float], path: Path) -> None:
    """Write an XYZ `_trj.xyz` with the energy on line 2 of each block."""
    blocks: List[str] = []
    E = np.array(energies, dtype=float)
    for geom, e in zip(images, E):
        if hasattr(geom, "as_xyz"):
            blocks.append(xyz_string_with_energy(geom, energy=float(e)))
            continue
        # ASE Atoms fallback
        symbols = geom.get_chemical_symbols()
        coords = geom.get_positions()
        lines = [str(len(symbols)), f"{float(e):.12f}"]
        lines.extend(
            f"{sym} {x:.15f} {y:.15f} {z:.15f}"
            for sym, (x, y, z) in zip(symbols, coords)
        )
        blocks.append("\n".join(lines) + "\n")
    with open(path, "w") as f:
        f.write("".join(blocks))


def set_freeze_atoms_or_warn(
    geom: Any,
    freeze_atoms: Sequence[int],
    *,
    context: str,
) -> None:
    """Attach freeze_atoms to a geometry; warn once on failure."""
    if freeze_atoms is None:
        return
    try:
        if isinstance(freeze_atoms, np.ndarray):
            if freeze_atoms.size == 0:
                return
        elif len(freeze_atoms) == 0:
            return
    except TypeError:
        if not freeze_atoms:
            return
    try:
        geom.freeze_atoms = np.array(sorted({int(i) for i in freeze_atoms}), dtype=int)
    except Exception:
        click.echo(f"[{context}] WARNING: Failed to attach freeze_atoms to geometry.", err=True)


def read_xyz_energies(path: Path | str) -> List[float]:
    """Extract Hartree energies with the strict shared XYZ parser."""
    from pdb2reaction.io.xyz_trajectory import read_xyz_trajectory

    parsed = read_xyz_trajectory(path, require_energies=True)
    return [float(value) for value in parsed["energies_ha"]]


def parse_xyz_block(
    block: Sequence[str],
    *,
    path: Path,
    frame_idx: int,
) -> Tuple[List[str], np.ndarray]:
    if not block:
        raise click.ClickException(f"[xyz] Empty XYZ frame in {path}")
    try:
        nat = int(block[0].strip().split()[0])
    except Exception:
        raise click.ClickException(
            f"[xyz] Malformed XYZ/TRJ header in frame {frame_idx} of {path}"
        )
    if len(block) < 2 + nat:
        raise click.ClickException(
            f"[xyz] Incomplete XYZ frame {frame_idx} in {path} (expected {nat} atoms)."
        )
    elems: List[str] = []
    coords: List[List[float]] = []
    for k in range(nat):
        parts = block[2 + k].split()
        if len(parts) < 4:
            raise click.ClickException(
                f"[xyz] Malformed atom line in frame {frame_idx} of {path}"
            )
        elems.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords, dtype=float)


def xyz_blocks_first_last(
    blocks: Sequence[Sequence[str]],
    *,
    path: Path,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    if not blocks:
        raise click.ClickException(f"[xyz] No frames found in {path}")
    first_elems, first_coords = parse_xyz_block(blocks[0], path=path, frame_idx=1)
    last_elems, last_coords = parse_xyz_block(blocks[-1], path=path, frame_idx=len(blocks))
    if first_elems != last_elems:
        raise click.ClickException(f"[xyz] Element list changed across frames in {path}")
    return first_elems, first_coords, last_coords


def read_xyz_first_last(trj_path: Path) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """
    Lightweight XYZ trajectory reader: return (elements, first_coords[Å], last_coords[Å]).
    Assumes standard multi-frame XYZ: natoms line, comment line, natoms atom lines.
    """
    blocks = read_xyz_as_blocks(trj_path, strict=True)
    return xyz_blocks_first_last(blocks, path=trj_path)


def read_xyz_as_blocks(trj_path: Path, *, strict: bool = False) -> List[List[str]]:
    """
    Read a multi-frame XYZ/TRJ file and return a list of frames, each as a list of lines.

    When *strict* is True, malformed headers or truncated frames raise a ClickException.
    """
    try:
        lines = trj_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as e:
        raise click.ClickException(f"[xyz] Failed to read XYZ/TRJ: {trj_path} ({e})")

    blocks: List[List[str]] = []
    i = 0
    n = len(lines)
    while i < n:
        while i < n and not lines[i].strip():
            i += 1
        if i >= n:
            break
        header = lines[i].strip()
        try:
            nat = int(header.split()[0])
        except Exception:
            if strict:
                raise click.ClickException(f"[xyz] Malformed header at line {i+1} in {trj_path}")
            break
        block = lines[i : i + 2 + nat]
        if len(block) < 2 + nat:
            if strict:
                raise click.ClickException(
                    f"[xyz] Incomplete frame at line {i+1} in {trj_path} (expected {nat} atoms)."
                )
            break
        blocks.append(block)
        i += 2 + nat
    return blocks


@dataclass
class PreparedInputStructure:
    source_path: Path
    geom_path: Path
    gjf_template: Optional[GjfTemplate] = None
    original_path: Optional[Path] = None
    structure_template: Optional[CoordinateTemplate] = None
    _tmp_geom_path: Optional[Path] = None
    _normalized_structures: List[Tuple[Path, Path]] = field(default_factory=list)

    @property
    def is_gjf(self) -> bool:
        return self.gjf_template is not None

    @property
    def is_cif(self) -> bool:
        return bool(
            self.original_path is not None
            and self.original_path.suffix.lower() in CIF_SUFFIXES
        )

    @property
    def display_path(self) -> Path:
        return self.original_path or self.source_path

    def cleanup(self) -> None:
        if self._tmp_geom_path and self._tmp_geom_path.exists():
            try:
                self._tmp_geom_path.unlink()
            except FileNotFoundError:
                pass
        for internal_path, tmp_dir in self._normalized_structures:
            cleanup_normalized_structure(internal_path, tmp_dir)
        self._normalized_structures.clear()

    def __enter__(self) -> "PreparedInputStructure":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.cleanup()

    def __del__(self) -> None:
        # Best-effort safety net for callers that miss explicit/context-manager
        # cleanup; interpreter shutdown does not guarantee finalizer timing.
        try:
            self.cleanup()
        except Exception:
            pass


def _parse_coord_line(line: str) -> GjfCoordinateLine:
    match = _GJF_COORD_RE.match(line)
    if not match:
        raise ValueError(f"Could not parse Gaussian coordinate line: '{line}'.")
    prefix = match.group("prefix")
    sep0 = match.group("sep0")
    sep1 = match.group("sep1")
    sep2 = match.group("sep2")
    suffix = match.group("suffix")
    x_token = match.group("x")
    y_token = match.group("y")
    z_token = match.group("z")
    symbol = _token_to_symbol(f"{prefix}{sep0}")
    return GjfCoordinateLine(
        prefix=prefix,
        sep0=sep0,
        sep1=sep1,
        sep2=sep2,
        suffix=suffix,
        x_template=x_token,
        y_template=y_token,
        z_template=z_token,
        symbol=symbol,
        x=_convert_float(x_token),
        y=_convert_float(y_token),
        z=_convert_float(z_token),
    )


def parse_gjf_template(path: Path) -> GjfTemplate:
    lines = path.read_text().splitlines()
    section = 0
    charge_line_idx = None
    charge = None
    spin = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            if section < 2:
                section += 1
            continue
        if section < 2:
            continue
        tokens = stripped.split()
        if len(tokens) < 2:
            continue
        try:
            charge = int(tokens[0])
            spin = int(tokens[1])
        except ValueError:
            continue
        charge_line_idx = idx
        break
    if charge_line_idx is None or charge is None or spin is None:
        raise ValueError(f"Failed to locate charge/multiplicity line in '{path}'.")

    coord_start = charge_line_idx + 1
    while coord_start < len(lines) and not lines[coord_start].strip():
        coord_start += 1

    prefix_lines = lines[:coord_start]
    coord_lines_raw: List[str] = []
    idx = coord_start
    while idx < len(lines):
        line = lines[idx]
        if not line.strip():
            break
        coord_lines_raw.append(line)
        idx += 1
    if not coord_lines_raw:
        raise ValueError(f"No coordinates found in '{path}'.")
    suffix_lines = lines[idx:]

    coord_lines = [_parse_coord_line(line) for line in coord_lines_raw]
    return GjfTemplate(
        path=path,
        prefix_lines=prefix_lines,
        suffix_lines=suffix_lines,
        coord_lines=coord_lines,
        charge=charge,
        spin=spin,
    )


def prepare_input_structure(path: Path) -> PreparedInputStructure:
    path = Path(path)
    suffix = path.suffix.lower()
    registered_template = coordinate_template_for(path)
    if registered_template is not None:
        return PreparedInputStructure(
            source_path=path,
            geom_path=path,
            original_path=registered_template.source_path,
            structure_template=registered_template,
        )
    if is_cif_path(path) or (suffix == ".pdb" and pdb_requires_normalization(path)):
        internal, structure_template, tmp_dir = normalize_structure_to_pdb(path)
        return PreparedInputStructure(
            source_path=internal,
            geom_path=internal,
            original_path=path.resolve(),
            structure_template=structure_template,
            _normalized_structures=[(internal, tmp_dir)],
        )
    if suffix != ".gjf":
        return PreparedInputStructure(source_path=path, geom_path=path, original_path=path)
    template = parse_gjf_template(path)
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(template.as_xyz_string())
        tmp.flush()
    finally:
        tmp.close()
    tmp_path = Path(tmp.name)
    return PreparedInputStructure(
        source_path=path,
        geom_path=tmp_path,
        gjf_template=template,
        original_path=path,
        _tmp_geom_path=tmp_path,
    )


def _read_atom_count(path: Path) -> int:
    try:
        atoms = read(path, index=0)
    except Exception as e:
        raise click.ClickException(f"Failed to read '{path}' for --ref-pdb validation: {e}")
    return len(atoms)


def _validate_ref_pdb_atom_count(geom_path: Path, ref_pdb_path: Path) -> None:
    geom_count = _read_atom_count(geom_path)
    ref_count = _read_atom_count(ref_pdb_path)
    if geom_count != ref_count:
        raise click.ClickException(
            f"--ref-pdb atom count ({ref_count}) does not match input ({geom_count})."
        )


def apply_ref_pdb_override(
    prepared_input: PreparedInputStructure,
    ref_pdb: Optional[Path],
) -> Optional[Path]:
    """Use a reference PDB/mmCIF topology while keeping input coordinates."""
    if ref_pdb is None:
        return None
    ref_pdb = Path(ref_pdb).resolve()
    if ref_pdb.suffix.lower() not in ({".pdb"} | set(CIF_SUFFIXES)):
        raise click.BadParameter("--ref-pdb must be a .pdb, .cif, or .mmcif file.")
    prepared_ref = prepare_input_structure(ref_pdb)
    try:
        _validate_ref_pdb_atom_count(prepared_input.geom_path, prepared_ref.geom_path)
    except BaseException:
        prepared_ref.cleanup()
        raise
    prepared_input.source_path = prepared_ref.source_path
    prepared_input.structure_template = prepared_ref.structure_template
    if prepared_ref.structure_template is not None:
        prepared_input.original_path = ref_pdb
        prepared_input._normalized_structures.extend(prepared_ref._normalized_structures)
        prepared_ref._normalized_structures.clear()
    return prepared_input.source_path


def fill_charge_spin_from_gjf(
    charge: Optional[int],
    spin: Optional[int],
    template: Optional[GjfTemplate],
) -> Tuple[Optional[int], Optional[int]]:
    if template is not None:
        if charge is None:
            charge = template.charge
        if spin is None:
            spin = template.spin
    return charge, spin


def fill_charge_spin_from_gjf_inputs(
    charge: Optional[int],
    spin: Optional[int],
    templates: Sequence[Optional[GjfTemplate]],
) -> Tuple[Optional[int], Optional[int]]:
    """Fill unresolved fields only when all GJF headers agree.

    Explicitly supplied fields remain authoritative.  Conflicts in an
    unresolved field are rejected with the contributing files and values
    rather than silently inheriting the first image's electronic surface.
    """
    present = [template for template in templates if template is not None]

    def _resolve_field(current: Optional[int], field: str) -> Optional[int]:
        if current is not None or not present:
            return current
        by_value: Dict[int, List[str]] = {}
        for template in present:
            value = int(getattr(template, field))
            by_value.setdefault(value, []).append(str(template.path))
        if len(by_value) > 1:
            details = "; ".join(
                f"{value}: {', '.join(paths)}"
                for value, paths in sorted(by_value.items())
            )
            label = "charge" if field == "charge" else "multiplicity"
            raise click.ClickException(
                f"Conflicting GJF {label} headers across reaction structures "
                f"({details}). Supply an explicit {'-q/--charge' if field == 'charge' else '-m/--multiplicity'} "
                "only if one value is intended for every image."
            )
        return next(iter(by_value))

    return _resolve_field(charge, "charge"), _resolve_field(spin, "spin")


def resolve_configured_charge_spin(
    yaml_cfg: Mapping[str, Any],
    *,
    charge: Optional[int],
    spin: Optional[int],
    ligand_charge: Optional[float | str | Dict[str, float]] = None,
) -> Tuple[Optional[int], Optional[int]]:
    """Apply ``calc.charge`` / ``calc.spin`` below explicit CLI inputs.

    The public YAML schema has always exposed these two calculator fields, but
    the workflow-level charge resolver used to overwrite them before it ever
    saw their configured values.  Keep the precedence explicit and shared:

    ``CLI -q/-m or --ligand-charge > YAML calc.charge/spin > GJF > defaults``.

    ``--ligand-charge`` is a user request to derive the total charge from the
    selected residues, so a configured ``calc.charge`` must not bypass it.
    """
    calc_cfg = yaml_cfg.get("calc") if isinstance(yaml_cfg, Mapping) else None
    if not isinstance(calc_cfg, Mapping):
        return charge, spin

    if charge is None and ligand_charge is None and "charge" in calc_cfg:
        raw_charge = calc_cfg.get("charge")
        if raw_charge is not None:
            try:
                charge = int(raw_charge)
            except (TypeError, ValueError) as exc:
                raise click.BadParameter(
                    f"calc.charge must be an integer, got {raw_charge!r}."
                ) from exc

    if spin is None and "spin" in calc_cfg:
        raw_spin = calc_cfg.get("spin")
        if raw_spin is not None:
            try:
                spin = int(raw_spin)
            except (TypeError, ValueError) as exc:
                raise click.BadParameter(
                    f"calc.spin must be an integer multiplicity, got {raw_spin!r}."
                ) from exc
            if spin < 1:
                raise click.BadParameter(
                    f"calc.spin must be an integer multiplicity >= 1, got {spin}."
                )

    return charge, spin


def _round_charge_with_note(q: float, prefix: str) -> int:
    if not math.isfinite(q):
        raise click.BadParameter(f"{prefix} Computed total charge is non-finite: {q!r}")
    q_int = int(round(q))
    if not math.isclose(q_int, q):
        click.echo(
            f"{prefix} NOTE: extractor total charge = {q:+g} → rounded to integer {q_int:+d}."
        )
    return q_int


def _derive_charge_from_ligand_charge(
    prepared: PreparedInputStructure,
    ligand_charge: Optional[float | str | Dict[str, float]],
    *,
    prefix: str,
) -> Optional[int]:
    if ligand_charge is None:
        return None
    try:
        from Bio import PDB

        parser = PDB.PDBParser(QUIET=True)
        complex_struct = parser.get_structure("complex", str(prepared.source_path))
        models = list(complex_struct.get_models())
        if len(models) > 1:
            # Refuse rather than derive from the first model: the geometry loader used for the
            # actual calculation does not honour MODEL/ENDMDL at all — it concatenates every
            # ATOM/HETATM record into one system — so a "first model" charge would describe a
            # structure the run never computes. `extract` may warn and continue only because it
            # detaches the other models first.
            raise click.BadParameter(
                f"Input '{prepared.source_path}' contains {len(models)} MODELs; "
                "--ligand-charge total-charge derivation requires a single-model PDB. "
                "Split the intended model into its own PDB before running the calculation.",
                param_hint="-i / --input",
            )
        selected_ids = {res.get_full_id() for res in complex_struct.get_residues()}
        from pdb2reaction.io.charge import infer_present_terminal_cap_ids

        keep_ncap_ids, keep_ccap_ids = infer_present_terminal_cap_ids(
            complex_struct,
            selected_ids,
        )
        summary = compute_charge_summary(
            complex_struct,
            selected_ids,
            set(),
            ligand_charge,
            keep_ncap_ids=keep_ncap_ids,
            keep_ccap_ids=keep_ccap_ids,
        )
        log_charge_summary(prefix, summary)
        q_total = float(summary.get("total_charge", 0.0))
        click.echo(
            f"{prefix} Charge summary from full complex (--ligand-charge without extraction):"
        )
        click.echo(
            f"  Protein: {summary.get('protein_charge', 0.0):+g},  "
            f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
            f"Ions: {summary.get('ion_total_charge', 0.0):+g},  "
            f"Total: {q_total:+g}"
        )
        return _round_charge_with_note(q_total, prefix)
    except click.BadParameter:
        # A rejected input is a decision, not a derivation failure: it must not degrade to the
        # "could not derive, carry on" path below.
        raise
    except Exception as e:
        click.echo(
            f"{prefix} NOTE: failed to derive total charge from --ligand-charge: {e}",
            err=True,
        )
        return None


def resolve_charge_spin(
    prepared_inputs: PreparedInputStructure | Sequence[PreparedInputStructure],
    charge: Optional[int],
    spin: Optional[int],
    *,
    spin_default: int = 1,
    charge_default: int = 0,
    ligand_charge: Optional[float | str | Dict[str, float]] = None,
    prefix: str = "[charge]",
    cleanup_on_error: Optional[Callable[[], None]] = None,
) -> Tuple[int, int]:
    """Resolve charge/spin from CLI args, ligand metadata, and GJF templates.

    Accepts either a single PreparedInputStructure or a sequence of them.
    """
    # Normalize to sequence
    if isinstance(prepared_inputs, PreparedInputStructure):
        inputs = [prepared_inputs]
        cleanup_on_error = cleanup_on_error or prepared_inputs.cleanup
    else:
        inputs = list(prepared_inputs)

    resolved_charge = charge
    resolved_spin = spin

    if ligand_charge is not None and resolved_charge is None:
        for prepared in inputs:
            if prepared.source_path.suffix.lower() in {".xyz", ".gjf"}:
                if cleanup_on_error:
                    cleanup_on_error()
                raise click.ClickException(
                    "--ligand-charge is only supported for PDB/mmCIF inputs; it cannot be used with .xyz or .gjf files without --ref-pdb."
                )
        resolved_charge = _derive_charge_from_ligand_charge(
            inputs[0], ligand_charge, prefix=prefix
        )

    resolved_charge, resolved_spin = fill_charge_spin_from_gjf_inputs(
        resolved_charge,
        resolved_spin,
        [prepared.gjf_template for prepared in inputs],
    )

    if resolved_charge is None:
        if any(not p.is_gjf for p in inputs):
            if cleanup_on_error:
                cleanup_on_error()
            raise click.ClickException(
                "-q/--charge is required unless the input is a .gjf template with charge metadata."
            )
        if cleanup_on_error:
            cleanup_on_error()
        raise click.ClickException(
            "Charge metadata was not found in the .gjf input(s). Provide -q/--charge "
            "or add charge/spin metadata to the .gjf template."
        )

    if resolved_spin is None:
        resolved_spin = spin_default
    return int(resolved_charge), int(resolved_spin)

@contextmanager
def prepared_cli_input(
    input_path: Path,
    *,
    ref_pdb: Optional[Path],
    charge: Optional[int],
    spin: Optional[int],
    ligand_charge: Optional[float | str | Dict[str, float]] = None,
    prefix: str = "[charge]",
) -> Iterator[Tuple[PreparedInputStructure, int, int]]:
    """Context-managed input preparation with charge/spin resolution."""
    with prepare_input_structure(input_path) as prepared:
        apply_ref_pdb_override(prepared, ref_pdb)
        charge_res, spin_res = resolve_charge_spin(
            prepared,
            charge,
            spin,
            ligand_charge=ligand_charge,
            prefix=prefix,
        )
        yield prepared, charge_res, spin_res


_CONVERT_FILES_ENABLED: bool = True


def set_convert_file_enabled(enabled: bool) -> None:
    """Globally enable or disable XYZ/TRJ conversions to PDB/CIF/GJF outputs.

    The toggle mirrors the ``--convert-files {True|False}`` CLI flag used
    by every workflow. When disabled, format-aware conversions are skipped even
    if reference templates are available.
    """

    global _CONVERT_FILES_ENABLED
    _CONVERT_FILES_ENABLED = bool(enabled)


def is_convert_file_enabled() -> bool:
    """Return the raw process-wide output-conversion toggle."""

    return _CONVERT_FILES_ENABLED


def convert_xyz_to_gjf(xyz_path: Path, template: GjfTemplate, out_path: Path) -> None:
    """Render single- or multi-frame XYZ/TRJ coordinates into a Gaussian template.

    A single frame produces a normal template-based Gaussian input. Multi-frame
    trajectories retain the historical one-header, blank-separated coordinate
    archive layout; that layout is not by itself an executable Gaussian
    QST/Link1 input and callers must split/select a frame before submission.
    """
    traj = read(xyz_path, index=":", format="xyz")
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")
    if len(traj) > 1:
        logger.warning(
            "Writing multi-frame GJF coordinate archive %s from %s. The "
            "one-header, blank-separated layout is not directly executable "
            "as a Gaussian QST/Link1 job; select one frame before submission.",
            out_path,
            xyz_path,
        )
    new_lines: List[str] = list(template.prefix_lines)
    for frame_idx, atoms in enumerate(traj):
        if len(atoms) != template.natoms:
            raise ValueError(
                f"Atom count mismatch for '{xyz_path}': xyz has {len(atoms)} atoms, "
                f"but template has {template.natoms}."
            )
        coords = atoms.get_positions()
        if frame_idx > 0:
            new_lines.append("")  # Historical coordinate-archive separator.
        for idx, coord_line in enumerate(template.coord_lines):
            new_lines.append(coord_line.render(tuple(map(float, coords[idx]))))
    new_lines.extend(template.suffix_lines)
    text = "\n".join(new_lines)
    if not text.endswith("\n"):
        text += "\n"
    out_path.write_text(text)


def convert_xyz_like_outputs(
    xyz_path: Path,
    prepared_input: PreparedInputStructure,
    *,
    ref_pdb_path: Optional[Path],
    out_pdb_path: Optional[Path] = None,
    out_gjf_path: Optional[Path] = None,
    context: str = "outputs",
    on_error: str = "raise",
) -> bool:
    """Convert XYZ/TRJ to requested PDB/CIF/GJF companion outputs.

    Parameters
    ----------
    xyz_path:
        Source XYZ-like file (single or multi-frame).
    prepared_input:
        Prepared input structure returned by :func:`prepare_input_structure`.
    ref_pdb_path:
        Reference PDB topology (required for PDB conversions).
    out_pdb_path / out_gjf_path:
        Targets for the converted files. Conversions are skipped when the
        corresponding output path is ``None`` or when the input type does not
        request that format.
    Returns True when at least one conversion was attempted and succeeded; False otherwise.
    """

    if not _CONVERT_FILES_ENABLED:
        return False

    needs_pdb = out_pdb_path is not None and ref_pdb_path is not None
    needs_gjf = (
        xyz_path.suffix.lower() == ".xyz"
        and prepared_input.is_gjf
        and prepared_input.gjf_template is not None
        and out_gjf_path is not None
    )

    if not (needs_pdb or needs_gjf):
        return False

    try:
        if needs_pdb:
            convert_xyz_to_pdb(xyz_path, ref_pdb_path, out_pdb_path)
        if needs_gjf:
            convert_xyz_to_gjf(xyz_path, prepared_input.gjf_template, out_gjf_path)
    except Exception as e:
        if on_error == "warn":
            click.echo(f"[convert] WARNING: Failed to convert {context}: {e}", err=True)
            return False
        raise click.ClickException(f"[convert] Failed to convert {context}: {e}") from e
    return True


def _convert_to_pdb_logged(
    src_path: Path, ref_pdb_path: Optional[Path], out_path: Optional[Path] = None
) -> Optional[Path]:
    """Convert an XYZ/TRJ to PDB when conversion is enabled; return path or None."""
    try:
        if ref_pdb_path is None or not _CONVERT_FILES_ENABLED:
            return None
        src_path = Path(src_path)
        if (not src_path.exists()) or src_path.suffix.lower() != ".xyz":
            return None
        out_path = out_path if out_path is not None else src_path.with_suffix(".pdb")
        convert_xyz_to_pdb(src_path, ref_pdb_path, out_path)
        if out_path.exists():
            click.echo(f"[convert] Wrote '{out_path}'.")
            return out_path
        return None
    except Exception as e:
        click.echo(
            f"[convert] WARNING: Failed to convert '{src_path}' to PDB: {e}",
            err=True,
        )
        return None


def parse_pdb_coords(pdb_path):
    """Parse ATOM/HETATM records from *pdb_path* and separate link hydrogen (HL) atoms.

    Returns:
        A tuple (others, lkhs) where:
            - others: list of tuples (index, x, y, z, line) for all atoms except the
              'HL' atom of residue 'LKH'. ``index`` is the 0-based position in the
              atom sequence as loaded from the *first* MODEL (or the full file if no
              MODEL records are present).
            - lkhs: list of tuples (x, y, z, line) for atoms where residue name is
              'LKH' and atom name is 'HL' in the same MODEL selection.

    Notes
    -----
        - Coordinates are read from standard PDB columns:
          X: columns 31–38, Y: 39–46, Z: 47–54 (1-based indexing).
        - If multiple MODEL blocks are present, only the first model is considered.
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    others = []
    lkhs = []
    model_seen = False
    in_first_model = True
    atom_index = 0
    for line in lines:
        if line.startswith("MODEL"):
            if not model_seen:
                model_seen = True
                in_first_model = True
            else:
                in_first_model = False
            continue
        if line.startswith("ENDMDL"):
            if model_seen and in_first_model:
                break
            continue
        if model_seen and not in_first_model:
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        current_index = atom_index
        atom_index += 1

        name = line[12:16].strip()
        resname = line[17:20].strip()
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        if resname == "LKH" and name == "HL":
            lkhs.append((x, y, z, line))
        else:
            others.append((current_index, x, y, z, line))
    return others, lkhs


def nearest_index(point, pool):
    """Find the nearest point in *pool* to *point* using Euclidean distance.

    Args:
        point: Tuple (x, y, z) representing the query coordinate.
        pool: Iterable of tuples (index, x, y, z, line) to search.

    Returns:
        A tuple (index, distance) where:
            - index is the 0-based index of the nearest entry in *pool* (or -1 if *pool* is empty).
            - distance is the Euclidean distance to that entry (``inf`` if *pool* is empty).
    """
    x, y, z = point
    best_i = -1
    best_d2 = float("inf")
    for atom_index, a, b, c, _ in pool:
        d2 = (a - x) ** 2 + (b - y) ** 2 + (c - z) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_i = atom_index
    return best_i, math.sqrt(best_d2)


def load_pdb_atom_metadata(pdb_path: Path) -> List[Dict[str, Any]]:
    """Return per-atom PDB metadata in file order, restoring original CIF IDs."""

    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            serial_txt = line[6:11].strip()
            resseq_txt = line[22:26].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            icode = line[26:27].strip()
            element_txt = line[76:78].strip()
            is_hetatm = line.startswith("HETATM")

            try:
                serial = int(serial_txt) if serial_txt else None
            except ValueError:
                serial = None
            try:
                resseq = int(resseq_txt) if resseq_txt else None
            except ValueError:
                resseq = None

            if not element_txt:
                inferred = guess_element(atom_name, res_name, is_hetatm)
                element_txt = inferred or ""

            atoms.append(
                {
                    "serial": serial,
                    "name": atom_name,
                    "resname": res_name,
                    "resseq": resseq,
                    "chain": chain_id,
                    "icode": icode,
                    "element": element_txt,
                }
            )

    template = coordinate_template_for(pdb_path)
    if template is not None:
        if len(atoms) != template.natoms:
            raise ValueError(
                f"PDB metadata atom count ({len(atoms)}) does not match retained "
                f"template ({template.natoms}) for {pdb_path}."
            )
        for meta, record in zip(atoms, template.records):
            meta.update(
                {
                    "name": record.atom_name,
                    "resname": record.resname,
                    "resseq": int(record.resseq) if re.fullmatch(r"[-+]?\d+", record.resseq) else record.resseq,
                    "chain": record.chain_id,
                    "icode": record.icode,
                    "element": record.element,
                }
            )

    return atoms


def _split_atom_spec_tokens(spec: str) -> List[str]:
    # Accept whitespace, comma, slash, backtick, or backslash separators.
    tokens = [t for t in re.split(r"[\s/:`,\\]+", spec.strip().replace(' ',',')) if t]
    return tokens


def resolve_atom_spec_index(spec: str, atom_meta: Sequence[Dict[str, Any]]) -> int:
    """
    Resolve an atom selector string into a 0-based atom index using PDB metadata.

    Three fields select by residue name, residue number, and atom name; their
    order remains flexible for compatibility. A fourth chain-ID field uses the
    unambiguous positional form ``CHAIN:RESNAME:RESSEQ[ICODE]:ATOM``.
    """
    tokens = _split_atom_spec_tokens(spec)
    if len(tokens) not in {3, 4}:
        raise ValueError(
            f"Atom spec '{spec}' must have 3 fields (resname, resseq, atomname) "
            "or 4 fields including chain ID."
        )

    tokens_upper = [t.upper() for t in tokens]
    canonical_four = spec.count(":") == 3 and all(
        part.strip() for part in spec.split(":")
    )
    canonical_parts = [part.strip() for part in spec.split(":")] if canonical_four else []
    matches: List[int] = []
    for idx, meta in enumerate(atom_meta):
        resname = (meta.get("resname") or "").strip().upper()
        resseq = meta.get("resseq")
        atom = (meta.get("name") or "").strip().upper()
        chain_text = (meta.get("chain") or "").strip()
        chain = chain_text.upper()
        if resseq is None:
            continue
        resseq_text = str(resseq)
        if canonical_four:
            chain_token = canonical_parts[0]
            resname_token = canonical_parts[1].upper()
            resseq_token = canonical_parts[2].upper()
            atom_token = canonical_parts[3].upper()
            numbered = re.fullmatch(r"(?P<number>[-+]?\d+)(?P<icode>[A-Z]?)", resseq_token)
            if numbered is not None:
                try:
                    same_resseq = int(numbered.group("number")) == int(resseq_text)
                except ValueError:
                    same_resseq = numbered.group("number") == resseq_text.upper()
                requested_icode = numbered.group("icode")
                if requested_icode:
                    same_resseq = same_resseq and requested_icode == str(meta.get("icode") or "").upper()
            else:
                same_resseq = resseq_token == resseq_text.upper()
            is_match = (
                chain_token == chain_text
                and resname_token == resname
                and same_resseq
                and atom_token == atom
            )
        else:
            normalized_tokens = [
                str(int(token)) if re.fullmatch(r"[-+]?\d+", token) else token
                for token in tokens_upper
            ]
            expected = [resname, resseq_text, atom]
            if len(tokens) == 4:
                expected.append(chain)
            is_match = Counter(normalized_tokens) == Counter(expected)
        if is_match:
            matches.append(idx)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Atom spec '{spec}' matches {len(matches)} atoms; add chain ID as "
            "CHAIN:RESNAME:RESSEQ[ICODE]:ATOM or use an explicit atom index."
        )

    if len(tokens) == 4 and not canonical_four:
        raise ValueError(
            f"Atom spec '{spec}' did not match any atom. Use the positional "
            "CHAIN:RESNAME:RESSEQ[ICODE]:ATOM form for chain-qualified selectors."
        )
    raise ValueError(f"Atom spec '{spec}' did not match any atom.")


def values_from_bounds(low: float, high: float, h: float) -> "np.ndarray":
    """Return evenly spaced values from low→high with step cap h (inclusive)."""
    if h <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    delta = abs(high - low)
    if delta < 1e-12:
        return np.array([low], dtype=float)
    N = int(math.ceil(delta / h))
    return np.linspace(low, high, N + 1, dtype=float)


def atom_label_from_meta(atom_meta: Sequence[Dict[str, Any]], index: int) -> str:
    if index < 0 or index >= len(atom_meta):
        return f"idx{index}"
    meta = atom_meta[index]
    resname = (meta.get("resname") or "?").strip() or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    icode = (meta.get("icode") or "").strip()
    if icode:
        resseq_txt += icode
    atom = (meta.get("name") or "?").strip() or "?"
    chain = (meta.get("chain") or "").strip()
    prefix = f"{chain}:" if chain else ""
    return f"{prefix}{resname}:{resseq_txt}:{atom}"


def axis_label_csv(
    axis_name: str,
    i_idx: int,
    j_idx: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
    pair_raw: Optional[Tuple[Any, Any, float, float]] = None,
) -> str:
    if pair_raw and (isinstance(pair_raw[0], str) or isinstance(pair_raw[1], str)) and atom_meta:
        i_label = atom_label_from_meta(atom_meta, i_idx)
        j_label = atom_label_from_meta(atom_meta, j_idx)
        return f"{axis_name}_{i_label}_{j_label}_A"
    i_disp = i_idx + 1 if one_based else i_idx
    j_disp = j_idx + 1 if one_based else j_idx
    return f"{axis_name}_{i_disp}_{j_disp}_A"


def axis_label_html(label: str) -> str:
    parts = label.split("_")
    if len(parts) >= 4 and parts[-1] == "A":
        axis = parts[0]
        i_disp = parts[1]
        j_disp = parts[2]
        return f"{axis} ({i_disp},{j_disp}) (Å)"
    return label


def parse_scan_list_quads_checked(
    raw: str,
    *,
    expected_len: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str,
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]]]:
    parsed, raw_pairs = parse_scan_list_quads(
        raw,
        expected_len=expected_len,
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=option_name,
    )
    for i, j, low, high in parsed:
        if low <= 0.0 or high <= 0.0:
            raise click.BadParameter(f"Distances must be positive: {(i, j, low, high)}")
    return parsed, raw_pairs


def parse_scan_list_triples(
    raw: str,
    *,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str,
    return_one_based: bool = False,
) -> Tuple[List[Tuple[int, int, float]], List[Tuple[Any, Any, float]]]:
    """Parse --scan-lists entries into indices (0-based by default).

    Accepts both 3-tuples ``(i, j, target)`` and 4-tuples
    ``(i, j, start, end)`` for bidirectional scans.  4-tuples are
    expanded into two 3-tuple stages (initial→start, then initial→end)
    by the caller in scan.py.

    The returned *parsed* list contains tuples of length 3 **or** 4:
    ``(i, j, target)`` or ``(i, j, start, end)``.
    """
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not isinstance(obj, (list, tuple)):
        raise click.BadParameter(
            f"{option_name} must be a list/tuple of (i,j,target) or (i,j,start,end)."
        )
    if len(obj) == 0:
        raise click.BadParameter(f"{option_name} must contain at least one atom pair.")

    parsed: list = []
    seen_pairs: set[tuple[int, int]] = set()
    for entry_idx, t in enumerate(obj, start=1):
        is_3 = (
            isinstance(t, (list, tuple))
            and len(t) == 3
            and isinstance(t[2], Real)
        )
        is_4 = (
            isinstance(t, (list, tuple))
            and len(t) == 4
            and isinstance(t[2], Real)
            and isinstance(t[3], Real)
        )
        if not (is_3 or is_4):
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i,j,target) or (i,j,start,end): got {t}"
            )

        i = resolve_scan_index(
            t[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            t[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if i == j:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} selects the same atom twice."
            )
        pair_key = tuple(sorted((i, j)))
        if pair_key in seen_pairs:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} repeats atom pair "
                f"{pair_key}; each simultaneous scan axis must be unique."
            )
        seen_pairs.add(pair_key)
        if return_one_based:
            i += 1
            j += 1
        if is_4:
            parsed.append((i, j, float(t[2]), float(t[3])))
        else:
            parsed.append((i, j, float(t[2])))

    return parsed, list(obj)


def parse_dist_freeze_list(
    raw: str,
    *,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str = "--dist-freeze",
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse ``--dist-freeze`` entries: ``(i,j)`` or ``(i,j,target_A)``.

    Uses the same :func:`resolve_scan_index` as ``--scan-lists``, so string
    atom specs (e.g. ``'A:SER:123:OG'``) are supported when PDB metadata is
    available.
    """
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not isinstance(obj, (list, tuple)):
        raise click.BadParameter(f"{option_name} must be a list/tuple of (i,j) or (i,j,target).")

    # Single tuple → wrap in list
    if obj and not isinstance(obj[0], (list, tuple)):
        obj = [obj]

    parsed: List[Tuple[int, int, Optional[float]]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for entry_idx, t in enumerate(obj, start=1):
        if not (isinstance(t, (list, tuple)) and len(t) in (2, 3)):
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i,j) or (i,j,target): got {t}"
            )
        i = resolve_scan_index(
            t[0], one_based=one_based, atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            t[1], one_based=one_based, atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if i == j:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} selects the same atom twice."
            )
        pair_key = tuple(sorted((i, j)))
        if pair_key in seen_pairs:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} repeats atom pair {pair_key}."
            )
        seen_pairs.add(pair_key)
        target: Optional[float] = None
        if len(t) == 3:
            if not isinstance(t[2], Real):
                raise click.BadParameter(
                    f"Target distance must be numeric in {option_name} entry {entry_idx}: {t}"
                )
            target = float(t[2])
            if target <= 0.0:
                raise click.BadParameter(
                    f"Target distance must be > 0 in {option_name} entry {entry_idx}: {t}"
                )
        parsed.append((i, j, target))
    return parsed


def parse_dist_freeze_spec(
    spec_path: Path,
    *,
    one_based_default: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str = "--dist-freeze",
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse a YAML/JSON dist-freeze spec file.

    Expected format::

        constraints:       # or "pairs" / "stages"
          - [1, 5, 1.4]   # (i, j, target_A) — target optional
          - [2, 6]         # freeze at current distance
        one_based: true    # optional, defaults to CLI value
    """
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    key, raw_list = _first_spec_field(spec_cfg, ("constraints", "pairs", "stages"))
    if key is None:
        raise click.BadParameter(
            f"{option_name} spec must define 'constraints', 'pairs', or 'stages'."
        )
    if not isinstance(raw_list, (list, tuple)) or len(raw_list) == 0:
        raise click.BadParameter(
            f"{option_name} field '{key}' must be a non-empty list."
        )

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name,
    )
    return parse_dist_freeze_list(
        repr(list(raw_list)),
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=f"{option_name} {key}",
    )


def unbiased_energy_hartree(geom, base_calc) -> float:
    """Evaluate the underlying MLIP energy (Hartree) without harmonic bias."""
    import numpy as np

    # geom.coords3d is always Cartesian (Bohr); geom.coords returns the
    # internal-coordinate vector for redund/dlc/tric coord types, which would
    # feed get_energy() a non-Cartesian array -> NaN/garbage energy (scan2d/3d).
    coords_bohr = np.asarray(geom.coords3d)
    elems = getattr(geom, "atoms", None)
    if elems is None:
        return float("nan")
    try:
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception as exc:
        click.echo(
            f"[energy] WARNING: energy evaluation failed: {exc}",
            err=True,
        )
        return float("nan")


def _canonical_mlip_precision(backend: str, value: Any) -> Optional[str]:
    """Return the backend-neutral public precision token (``fp32``/``fp64``)."""
    if backend == "custom":
        return None
    token = "" if value is None else str(value).strip().lower()
    if backend == "aimnet2":
        return "fp32"
    if token in {"fp64", "float64", "double", "highest"}:
        return "fp64"
    if token in {"fp32", "float32", "float32-high", "float32-highest", "single"}:
        return "fp32"
    return token or None


def calculator_provenance(calc_cfg: Mapping[str, Any]) -> Dict[str, Any]:
    """Return backend-neutral MLIP provenance for machine-readable outputs."""
    from pdb2reaction.core.defaults import CALC_KW_DEFAULT

    backend = str(calc_cfg.get("backend") or CALC_KW_DEFAULT["backend"]).lower()
    if backend == "custom":
        calc_file = calc_cfg.get("calc_file")
        factory = calc_cfg.get("calc_factory") or "get_calculator"
        model = f"{Path(calc_file).name}:{factory}" if calc_file else str(factory)
        precision = None
    else:
        model = calc_cfg.get("model")
        if model is None:
            from pdb2reaction.core.defaults import apply_backend_defaults

            resolved = dict(CALC_KW_DEFAULT)
            resolved.update(calc_cfg)
            apply_backend_defaults(resolved)
            model = resolved.get("model")
        raw_precision = (
            calc_cfg.get("default_dtype")
            if backend == "mace"
            else calc_cfg.get("precision")
        )
        if raw_precision is None or str(raw_precision).strip().lower() == "auto":
            raw_precision = {
                "uma": "fp32",
                "orb": "fp64",
                "mace": "fp64",
                "aimnet2": "fp32",
            }.get(backend)
        precision = _canonical_mlip_precision(backend, raw_precision)

    return {
        "mlip_backend": backend,
        "mlip_model": None if model is None else str(model),
        "mlip_precision": precision,
    }


def close_matplotlib_figures() -> None:
    """Best-effort cleanup for matplotlib figures to avoid open-figure warnings."""
    try:
        import matplotlib.pyplot as plt  # type: ignore
        plt.close("all")
    except Exception:
        pass


def resolve_scan_index(
    value: Any,
    *,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    context: str,
) -> int:
    """Resolve an index or atom-spec string for scan lists with consistent errors."""
    if isinstance(value, Integral):
        idx_val = int(value)
        if one_based:
            idx_val -= 1
        if idx_val < 0:
            raise click.BadParameter(
                f"Negative atom index after base conversion in {context}: {idx_val} (0-based expected)."
            )
        # Out-of-range upper bound: when atom metadata is available we know
        # the atom count, so reject an index past the end here (clean error,
        # consistent for bracketed and unbracketed specs) instead of letting
        # --dry-run pass and only crashing in the real run.
        if atom_meta and idx_val >= len(atom_meta):
            raise click.BadParameter(
                f"{context}: atom index {idx_val} (0-based) is out of range "
                f"for a {len(atom_meta)}-atom structure."
            )
        return idx_val
    if isinstance(value, str):
        if not atom_meta:
            raise click.BadParameter(
                f"{context} uses a string atom spec, but no PDB metadata is available."
            )
        try:
            return resolve_atom_spec_index(value, atom_meta)
        except ValueError as exc:
            raise click.BadParameter(f"{context} {exc}")
    raise click.BadParameter(f"{context} must be an int index or atom spec string.")


def parse_scan_list_quads(
    raw: str,
    *,
    expected_len: int,
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
    option_name: str,
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]]]:
    """Parse --scan-lists quadruples into 0-based indices."""
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not (isinstance(obj, (list, tuple)) and len(obj) == expected_len):
        quads = ",".join([f"(i{n},j{n},low{n},high{n})" for n in range(1, expected_len + 1)])
        raise click.BadParameter(
            f"{option_name} must contain exactly {expected_len} quadruples: [{quads}]"
        )

    parsed: List[Tuple[int, int, float, float]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for entry_idx, q in enumerate(obj, start=1):
        if not (
            isinstance(q, (list, tuple))
            and len(q) == 4
            and isinstance(q[2], Real)
            and isinstance(q[3], Real)
        ):
            arity_hint = ""
            if isinstance(q, (list, tuple)):
                if len(q) == 3:
                    arity_hint = (
                        " — 3-tuple is the `scan` command's format (i, j, target); "
                        "scan2d/scan3d need the 4-tuple form (i, j, low, high)"
                    )
                elif len(q) == 5:
                    arity_hint = (
                        " — 5-tuple is not accepted; step count per axis is set "
                        "via --max-step-size, not inside the tuple"
                    )
                elif len(q) != 4:
                    arity_hint = f" — expected 4-tuple, got {len(q)}-tuple"
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i, j, low, high): "
                f"got {q}{arity_hint}"
            )

        i = resolve_scan_index(
            q[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            q[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if i == j:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} selects the same atom twice."
            )
        pair_key = tuple(sorted((i, j)))
        if pair_key in seen_pairs:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} repeats atom pair "
                f"{pair_key}; each scan axis must be unique."
            )
        seen_pairs.add(pair_key)
        parsed.append((i, j, float(q[2]), float(q[3])))

    return parsed, list(obj)


def _load_scan_spec_root(
    spec_path: Path,
    *,
    option_name: str = "--scan-lists",
) -> Mapping[str, Any]:
    """Load a scan spec (YAML/JSON) and ensure mapping root."""
    try:
        with open(spec_path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except Exception as exc:
        raise click.BadParameter(
            f"Failed to parse {option_name} file '{spec_path}': {exc}"
        )

    if data is None:
        raise click.BadParameter(f"{option_name} file '{spec_path}' is empty.")
    if not isinstance(data, Mapping):
        raise click.BadParameter(
            f"{option_name} file '{spec_path}' must have a mapping at the YAML/JSON root."
        )
    # A misspelled root key used to be dropped in silence, so the whole spec file became a
    # no-op and the run continued with no scan/freeze constraints at all.
    unknown = sorted(set(data) - {"stages", "pairs", "constraints", "one_based"})
    if unknown:
        raise click.BadParameter(
            f"{option_name} file '{spec_path}' has unrecognized root key(s): "
            f"{', '.join(unknown)}. Expected any of: stages, pairs, constraints, one_based."
        )
    return data


def _spec_one_based(
    value: Any,
    *,
    default: bool,
    option_name: str = "--scan-lists",
) -> bool:
    """Resolve one_based value from spec with CLI fallback."""
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        key = value.strip().lower()
        if key in {"1", "true", "yes", "y", "on"}:
            return True
        if key in {"0", "false", "no", "n", "off"}:
            return False
    raise click.BadParameter(
        f"{option_name} field 'one_based' must be a boolean (true/false)."
    )


def _first_spec_field(
    spec_cfg: Mapping[str, Any],
    candidates: _Sequence[str],
) -> Tuple[Optional[str], Any]:
    for key in candidates:
        if key in spec_cfg:
            return key, spec_cfg[key]
    return None, None


def is_scan_spec_file(value: str) -> bool:
    """Return True if *value* looks like an existing YAML/JSON scan spec file."""
    p = Path(value)
    return p.is_file() and p.suffix.lower() in {".yaml", ".yml", ".json"}


def parse_scan_spec_stages(
    spec_path: Path,
    *,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--scan-lists",
    return_bidirectional_markers: bool = False,
) -> Any:
    """Parse staged 1D scan spec into 0-based executable stages.

    A stage containing only ``(i, j, target)`` entries stays simultaneous.
    The legacy ``(i, j, start, end)`` form expands exactly like the inline
    CLI form: snapshot before the first leg, restore before the second.
    """
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    stages_key, stages_raw = _first_spec_field(spec_cfg, ("stages",))
    if stages_key is None:
        raise click.BadParameter(f"{option_name} must define 'stages'.")
    if not isinstance(stages_raw, (list, tuple)) or len(stages_raw) == 0:
        raise click.BadParameter(f"{option_name} field '{stages_key}' must be a non-empty list.")

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name
    )
    stages: List[List[Tuple[int, int, float]]] = []
    reset_before: set[int] = set()
    snapshot_before: set[int] = set()
    for stage_idx, stage_raw in enumerate(stages_raw, start=1):
        if not isinstance(stage_raw, (list, tuple)):
            raise click.BadParameter(
                f"{option_name} {stages_key}[{stage_idx}] must be a list of (i,j,target) entries."
            )
        parsed, _ = parse_scan_list_triples(
            repr(list(stage_raw)),
            one_based=one_based,
            atom_meta=atom_meta,
            option_name=f"{option_name} {stages_key}[{stage_idx}]",
        )
        if not parsed:
            raise click.BadParameter(
                f"{option_name} {stages_key}[{stage_idx}] must contain at least one (i,j,target) triple."
            )
        for entry in parsed:
            if any(float(distance) <= 0.0 for distance in entry[2:]):
                raise click.BadParameter(
                    f"Non-positive target distance in {option_name} "
                    f"{stages_key}[{stage_idx}]: {entry}."
                )
        if any(len(entry) == 4 for entry in parsed):
            for entry in parsed:
                if len(entry) == 4:
                    i, j, start, end = entry
                    first_leg = len(stages)
                    stages.append([(i, j, start)])
                    snapshot_before.add(first_leg)
                    reset_before.add(first_leg + 1)
                    stages.append([(i, j, end)])
                else:
                    stages.append([entry])
        else:
            stages.append(parsed)
    if return_bidirectional_markers:
        return (
            stages,
            one_based,
            frozenset(snapshot_before),
            frozenset(reset_before),
        )
    return stages, one_based


def parse_scan_spec_quads(
    spec_path: Path,
    *,
    expected_len: int,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--scan-lists",
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]], bool]:
    """Parse 2D/3D scan spec into 0-based quad tuples."""
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    pairs_key, pairs_raw = _first_spec_field(spec_cfg, ("pairs",))
    if pairs_key is None:
        raise click.BadParameter(f"{option_name} must define 'pairs'.")
    if not isinstance(pairs_raw, (list, tuple)):
        raise click.BadParameter(f"{option_name} field '{pairs_key}' must be a list.")

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name
    )
    parsed, raw_pairs = parse_scan_list_quads_checked(
        repr(list(pairs_raw)),
        expected_len=expected_len,
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=f"{option_name} {pairs_key}",
    )
    return parsed, raw_pairs, one_based


def format_pdb_atom_metadata_header() -> str:
    """Column legend for :func:`format_pdb_atom_metadata`, aligned to match values."""

    return f"{'id':>5} {'atom':<4} {'res':<4} {'resid':>6} {'chain':<8} {'el':<2}"


def format_pdb_atom_metadata(atom_meta: Sequence[Dict[str, Any]], index: int) -> str:
    """Format metadata for atom *index* as aligned text: serial name resname resseq element."""

    fallback_serial = index + 1
    if index < 0 or index >= len(atom_meta):
        return f"{fallback_serial:>5} {'?':<4} {'?':<4} {'?':>6} {'?':<8} {'?':<2}"

    meta = atom_meta[index]
    serial = meta.get("serial") or fallback_serial
    name = meta.get("name") or "?"
    resname = meta.get("resname") or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    chain = (meta.get("chain") or "?").strip() or "?"
    element = (meta.get("element") or "?").strip() or "?"

    return f"{serial:>5} {name:<4} {resname:<4} {resseq_txt:>6} {chain:<8} {element:<2}"


def detect_freeze_links(pdb_path):
    """Identify link-parent atom indices for 'LKH'/'HL' link hydrogens.

    For each 'HL' atom in residue 'LKH', find the nearest atom among all other
    ATOM/HETATM records and return the indices of those nearest neighbors in the
    same atom ordering used by geometry loading (first MODEL if present).

    Args:
        pdb_path: Path to the input PDB file.

    Returns:
        List of 0-based indices into the full atom sequence (including any link H atoms)
        corresponding to the nearest neighbors (link parents). Returns an empty list if
        no LKH/HL atoms are present or if link hydrogens exist without any other atoms.
    """
    others, lkhs = parse_pdb_coords(pdb_path)

    if not lkhs or not others:
        # --freeze-links defaults on. Returning [] freezes nothing, and the run
        # then reads as "cap parents frozen" when no cap hydrogen was found at
        # all -- usually a structure that was never passed through `extract`.
        if not lkhs:
            click.echo(
                f"[freeze-links] WARNING: no LKH/HL cap hydrogen found in "
                f"'{pdb_path}'; no cap parent is frozen despite --freeze-links.",
                err=True,
            )
        return []

    indices = []
    for (x, y, z, _line) in lkhs:
        idx, dist = nearest_index((x, y, z), others)
        if idx >= 0:
            indices.append(idx)
    return indices


def detect_freeze_links_logged(pdb_path: Path) -> List[int]:
    """Return link-parent indices and raise a user-facing error on failure."""
    try:
        return list(detect_freeze_links(pdb_path))
    except Exception as e:  # pragma: no cover - defensive logging helper
        raise click.ClickException(
            f"[freeze-links] Failed to detect link parents for '{pdb_path.name}': {e}"
        ) from e


def merge_detected_freeze_links(
    geom_cfg: Dict[str, Any],
    pdb_path: Path,
    *,
    prefix: str = "[freeze-links]",
) -> List[int]:
    """Detect link-parent atoms and merge them into ``geom_cfg['freeze_atoms']``."""
    detected = detect_freeze_links_logged(pdb_path)
    merged = merge_freeze_atom_indices(geom_cfg, detected)
    if merged:
        click.echo(f"{prefix} Freeze atoms (1-based): {','.join(str(i+1) for i in merged)}")
    return merged


def resolve_freeze_atoms(
    geom_cfg: Dict[str, Any],
    source_path: Optional[Path],
    freeze_links: bool,
    *,
    prefix: str = "[freeze-links]",
    on_error: str = "raise",
) -> List[int]:
    """Normalize freeze_atoms and optionally merge detected link-parent atoms."""
    merge_freeze_atom_indices(geom_cfg)
    if not freeze_links or source_path is None or source_path.suffix.lower() != ".pdb":
        return list(geom_cfg.get("freeze_atoms", []))
    try:
        return merge_detected_freeze_links(geom_cfg, source_path, prefix=prefix)
    except Exception as e:
        if on_error == "warn":
            click.echo(f"{prefix} WARNING: Could not detect link parents: {e}", err=True)
            return list(geom_cfg.get("freeze_atoms", []))
        raise


def load_prepared_geometries(
    prepared_inputs: Sequence["PreparedInputStructure"],
    *,
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
    prefix: str = "[freeze-links]",
) -> List[Any]:
    """Load multiple PreparedInputStructure geometries and apply freeze atom logic."""
    geoms: List[Any] = []

    for prepared in prepared_inputs:
        src_path = prepared.source_path
        geom_path = prepared.geom_path
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        freeze = resolve_freeze_atoms(
            cfg,
            src_path,
            auto_freeze_links,
            prefix=f"{prefix} {src_path.name}:",
        )

        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)

    return geoms




def _collect_environment_info() -> dict:
    """Collect compute environment info (CPU, RAM, GPU, VRAM, resolved device)."""
    import platform
    env: dict = {}
    try:
        import torch
        cuda_ok = torch.cuda.is_available()
        env["device"] = "cuda" if cuda_ok else "cpu"
        if cuda_ok:
            try:
                env["gpu_name"] = torch.cuda.get_device_name(0)
                props = torch.cuda.get_device_properties(0)
                vram = getattr(props, "total_memory", None) or getattr(props, "total_mem", None)
                if vram:
                    env["gpu_vram_gb"] = round(vram / 1e9, 1)
            except Exception:
                pass
            env["cuda_version"] = getattr(torch.version, "cuda", None) or "unknown"
    except Exception:
        env["device"] = "cpu"
    try:
        import os
        # CPU model
        cpu_info = platform.processor()
        if not cpu_info or cpu_info == "x86_64":
            try:
                with open("/proc/cpuinfo") as f:
                    for line in f:
                        if "model name" in line:
                            cpu_info = line.split(":")[1].strip()
                            break
            except Exception:
                pass
        env["cpu"] = cpu_info or "unknown"
        env["n_cpus"] = os.cpu_count()
        # RAM
        try:
            import psutil
            env["ram_gb"] = round(psutil.virtual_memory().total / 1e9, 1)
        except ImportError:
            pass
    except Exception:
        pass
    return env


# Schema version for result/summary JSON. Version 2.0 removes the UMA-specific
# all-workflow energy keys in favor of backend-neutral MLIP keys.
RESULT_JSON_SCHEMA_VERSION = "2.0"

# Union of public command-specific values for the ``status`` field. Individual
# commands intentionally use narrower enums; see docs/json-output.md.
RESULT_JSON_STATUS_VALUES = (
    "completed",
    "converged",
    "error",
    "failed",
    "not_converged",
    "ok",
    "partial",
    "stalled",
    "success",
    "unknown",
)


def write_result_json(
    out_dir: Path,
    data: dict,
    *,
    command: str,
    elapsed_seconds: Optional[float] = None,
    filename: str = "result.json",
    also_write_summary_json: bool = True,
) -> Path:
    """Write a machine-readable result.json for a subcommand.

    The ``data`` dict is augmented with common envelope fields
    (``command``, ``pdb2reaction_version``, ``schema_version``,
    ``status``, ``elapsed_seconds``, ``files``, ``environment``) and
    serialized as indented JSON.

    When ``also_write_summary_json`` is True (default) the same payload
    is mirrored to ``summary.json`` so downstream consumers can converge
    on a single filename across every subcommand (``all`` and
    ``path-search`` already use ``summary.json``).

    Returns the authoritative primary path.  Serialization, staging, and
    publication failures raise ``ResultCommitError``; they never report a
    successful path or silently suppress a failed compatibility mirror.
    """
    from pdb2reaction.core.result_commit import apply_current_run_id, commit_json

    try:
        from pdb2reaction._version import __version__
    except ImportError:
        __version__ = "unknown"

    data = dict(data)
    data.setdefault("command", command)
    data.setdefault("pdb2reaction_version", __version__)
    data.setdefault("schema_version", RESULT_JSON_SCHEMA_VERSION)
    data.setdefault("status", "unknown")
    # Preserve MLIP provenance as two independent, backend-neutral fields.
    # Existing leaf workflows commonly populate backend/model; normalize them
    # here so every result.json has the same machine-readable contract.
    if data.get("backend") is not None:
        data.setdefault("mlip_backend", data.get("backend"))
    if "model" in data:
        data.setdefault("mlip_model", data.get("model"))
    if elapsed_seconds is not None:
        data["elapsed_seconds"] = round(elapsed_seconds, 3)
    data.setdefault("environment", _collect_environment_info())

    # Convert non-serializable objects for json.dump
    def _to_json(obj):
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, dict):
            return {k: _to_json(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_to_json(i) for i in obj]
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        # numpy scalar / array
        try:
            import numpy as _np
            if isinstance(obj, _np.generic):
                return _to_json(obj.item())
            if isinstance(obj, _np.ndarray):
                return _to_json(obj.tolist())
        except ImportError:
            pass
        # torch tensor
        try:
            import torch as _th
            if isinstance(obj, _th.Tensor):
                return _to_json(obj.detach().cpu().tolist())
        except ImportError:
            pass
        return obj

    data = _to_json(data)

    data = apply_current_run_id(data)
    dest = Path(out_dir) / filename
    mirrors: tuple[Path, ...] = ()
    if also_write_summary_json and Path(filename).name != "summary.json":
        mirrors = (Path(out_dir) / "summary.json",)
    return commit_json(dest, data, mirrors=mirrors)


_ALLOW_CHARGE_MULT_MISMATCH = False


def set_allow_charge_mult_mismatch(value: bool = True) -> None:
    """Process-global toggle for ``--allow-charge-mult-mismatch`` (set by the CLI eager callback)."""
    global _ALLOW_CHARGE_MULT_MISMATCH
    _ALLOW_CHARGE_MULT_MISMATCH = bool(value)


def validate_charge_spin(elements, charge, multiplicity):
    """Raise ValueError if sum_Z(elements) - charge has the wrong parity for multiplicity,
    unless ``--allow-charge-mult-mismatch`` was set (then log a warning and skip)."""
    from pysisyphus.elem_data import ATOMIC_NUMBERS

    multiplicity = int(multiplicity)
    if multiplicity < 1:
        raise ValueError(
            f"Spin multiplicity must be an integer >= 1, got {multiplicity}."
        )
    sum_z = sum(ATOMIC_NUMBERS[str(e).lower()] for e in elements)
    total = sum_z - int(charge)
    unpaired = multiplicity - 1
    if total < unpaired or (total - unpaired) % 2:
        if _ALLOW_CHARGE_MULT_MISMATCH:
            import logging
            logging.getLogger(__name__).warning(
                "Cluster electron-parity check SKIPPED (--allow-charge-mult-mismatch): "
                "sum_Z=%d, charge=%d, total_electrons=%d, multiplicity=%d -- proceeding; "
                "make sure this charge/multiplicity is intentional.",
                sum_z, charge, total, multiplicity,
            )
            return
        raise ValueError(
            f"Cluster electron count inconsistent: sum_Z={sum_z}, charge={charge}, "
            f"total_electrons={total}, multiplicity={multiplicity}. Adjust charge (-q / -l) or "
            f"multiplicity (-m) so the electron count matches the spin state (e.g. -m 2 for an odd "
            f"electron count). Common cause: a modified/covalently-linked residue whose cluster cut "
            f"left an unpaired electron, or an open-shell system. If this charge/multiplicity is "
            f"intentional, pass --allow-charge-mult-mismatch to skip this check."
        )


def validate_charge_spin_at_path(path, charge, multiplicity):
    """Read elements at the given geometry path and run validate_charge_spin.

    This wrapper is the canonical CLI entry point — bare ValueError from
    the underlying parity check (or from ase.io.read on a malformed
    geometry) is re-raised as click.UsageError so the standard CLI
    exception handler (cli_utils.run_cli) renders a single-line error
    instead of a full Python traceback.
    """
    from ase.io import read as _ase_read

    try:
        atoms = _ase_read(path)
    except FileNotFoundError as e:
        raise click.UsageError(f"Geometry file not found: {path} ({e})")
    except Exception as e:
        raise click.UsageError(f"Failed to read geometry {path}: {e}")
    try:
        validate_charge_spin(atoms.get_chemical_symbols(), charge, multiplicity)
    except ValueError as e:
        raise click.UsageError(str(e))


def validate_charge_spin_for_prepared(prepared_or_list, charge, multiplicity):
    """Read elements from PreparedInputStructure(s) and run validate_charge_spin per geometry.

    Same wrapping as validate_charge_spin_at_path — see that function's
    docstring for the rationale.
    """
    items = (
        list(prepared_or_list)
        if isinstance(prepared_or_list, (list, tuple))
        else [prepared_or_list]
    )
    for prepared in items:
        validate_charge_spin_at_path(prepared.geom_path, charge, multiplicity)


# Compatibility re-export: the bounded-peak Hessian symmetrizer lives in the
# bundled-engine layer (``pysisyphus.normal_modes``) so the pure
# normal-mode kernel there stays free of any upward ``pdb2reaction`` import. It is
# re-exported here so existing callers of
# ``pdb2reaction.core.utils.symmetrize_inplace`` (backends/uma, workflows/tsopt,
# tests) keep resolving to the SAME function object.
from pysisyphus.normal_modes import symmetrize_inplace  # noqa: F401,E402


# ---------------------------------------------------------------------------
# XYZTrajectoryWriter — per-frame streaming XYZ writer (tail-able)
# ---------------------------------------------------------------------------


class XYZTrajectoryWriter:
    """Per-frame streaming XYZ writer (open / write+flush / close).

    Usage::

        with XYZTrajectoryWriter(path, mode="w") as w:
            for geom in geoms:
                w.write_raw(_coords3d_to_xyz_string(geom))
    """

    __slots__ = ("path", "_mode", "_handle", "frames_written")

    def __init__(self, path, *, mode: str = "w") -> None:
        if mode not in ("w", "a"):
            raise ValueError(
                f"XYZTrajectoryWriter mode must be 'w' or 'a', got {mode!r}"
            )
        self.path = Path(path)
        self._mode = mode
        self._handle = None
        self.frames_written = 0

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open(self._mode, encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._handle is not None:
            try:
                self._handle.flush()
            finally:
                self._handle.close()
            self._handle = None

    def write_block(self, frame_xyz: str, energy=None) -> None:
        if energy is not None:
            lines = frame_xyz.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{float(energy):.12f}"
            frame_xyz = "\n".join(lines)
        self.write_raw(frame_xyz)

    def write_raw(self, block: str, energy=None) -> None:
        if self._handle is None:
            raise RuntimeError(
                "XYZTrajectoryWriter must be used as a context manager."
            )
        if energy is not None:
            self.write_block(block, energy=energy)
            return
        if not block.endswith("\n"):
            block = block + "\n"
        self._handle.write(block)
        self._handle.flush()
        self.frames_written += 1
