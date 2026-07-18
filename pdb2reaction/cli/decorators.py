"""CLI utilities for standardized exception handling and shared boilerplate."""

from __future__ import annotations

import gc
import sys
import textwrap
import time
import traceback
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple, Type

import click

from pdb2reaction.core.utils import deep_update, load_yaml_dict


def resolve_yaml_sources(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], bool]:
    """Resolve internal YAML-source compatibility inputs.

    Public commands expose one ``--config`` layer. ``override_yaml`` and
    ``args_yaml_legacy`` remain callable compatibility inputs for embedded
    callers; they are not additional public CLI layers.
    """
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use a single YAML source option."
        )
    if args_yaml_legacy is not None:
        return config_yaml, args_yaml_legacy, True
    return config_yaml, override_yaml, False


def load_merged_yaml_cfg(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Load and merge YAML config and override files.

    Returns ``(merged, config_dict, override_dict)`` so that callers can
    use the individual layers for staged ``apply_yaml_overrides`` and the
    merged dict for ``show_config`` display without re-reading the files.
    """
    config_dict = load_yaml_dict(config_yaml)
    override_dict = load_yaml_dict(override_yaml)
    # Keep the returned provenance layers immutable from the perspective of
    # callers. ``deep_update`` assigns previously unseen nested mappings by
    # reference, so both inputs must be copied before forming the effective
    # tree.
    merged: Dict[str, Any] = deepcopy(config_dict)
    deep_update(merged, deepcopy(override_dict))
    return merged, config_dict, override_dict


def _write_error_json_with_chain(exc, klass, label):
    """Build the structured error payload (used by _write_error_json).

    Returns a dict containing the legacy keys (``status`` / ``error`` /
    ``error_type``) plus ``error_class_chain`` / ``error_module`` /
    ``error_label`` so MCP clients can pattern-match the class hierarchy.
    """
    return {
        "status": "error",
        "error": str(exc),
        "error_type": klass.__name__,
        "error_class_chain": [
            c.__name__ for c in klass.__mro__ if c is not object
        ],
        "error_module": klass.__module__,
        "error_label": label,
    }


def _write_error_json(
    out_dir: Path, command: str, exc: Exception, label: str,
    time_start: Optional[float] = None,
) -> None:
    """Write result.json + summary.json with a structured error envelope.

    Adds ``error_class_chain`` / ``error_module`` / ``error_label`` so
    MCP clients can pattern-match the class hierarchy without parsing
    text. Legacy ``status`` / ``error`` / ``error_type`` keys are
    preserved for backward compatibility.
    """
    try:
        from pdb2reaction.core.utils import write_result_json
        elapsed = time.perf_counter() - time_start if time_start else None
        write_result_json(
            out_dir,
            _write_error_json_with_chain(exc, type(exc), label),
            command=command,
            elapsed_seconds=elapsed,
        )
    except OSError:
        pass  # Best-effort; don't mask the original error


def render_cli_exception(
    e,
    *,
    label: str,
    out_dir=None,
    command=None,
    time_start=None,
):
    """Shared terminal renderer for a subcommand's top-level exception.

    User-input errors (Click UsageError/BadParameter/ClickException, or a
    malformed --config/override YAML) print a clean one-line ``Error:`` and
    exit with the conventional code (2). Everything else falls through to
    the full traceback (exit 1) so genuine internal bugs stay visible — no
    exception masking. Always calls sys.exit (never returns).
    """
    if out_dir is not None and command:
        _write_error_json(Path(out_dir).resolve(), command, e, "UnhandledError", time_start)
    if isinstance(e, click.ClickException):
        e.show()
        click.echo(
            f"Try 'pdb2reaction {command} -h' for help." if command else "Try 'pdb2reaction -h' for help.",
            err=True,
        )
        sys.exit(e.exit_code)
    try:
        import yaml as _yaml
        _is_yaml_err = isinstance(e, _yaml.YAMLError)
    except ImportError:
        _is_yaml_err = False
    if _is_yaml_err:
        click.echo(
            f"Error: invalid YAML in a configuration/override file: {e}",
            err=True,
        )
        click.echo(
            f"Try 'pdb2reaction {command} -h' for help." if command else "Try 'pdb2reaction -h' for help.",
            err=True,
        )
        sys.exit(2)
    # Emit the full traceback once via the human-readable terminal echo so
    # genuine internal bugs stay visible without duplicate output.
    tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
    click.echo(
        f"Unhandled error during {label}:\n" + textwrap.indent(tb, "  "),
        err=True,
    )
    sys.exit(1)


def run_cli(
    fn: Callable[[], None],
    *,
    label: str,
    zero_step_exc: Optional[Type[BaseException]] = None,
    zero_step_msg: Optional[str] = None,
    opt_exc: Optional[Type[BaseException]] = None,
    opt_msg: Optional[str] = None,
    out_dir: Optional[Path] = None,
    command: Optional[str] = None,
    time_start: Optional[float] = None,
) -> None:
    """Standard CLI exception handling with consistent messaging."""
    try:
        fn()
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        if out_dir and command:
            _write_error_json(out_dir, command, e, label, time_start)
        if zero_step_exc is not None and isinstance(e, zero_step_exc):
            click.echo(
                zero_step_msg
                or "ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).",
                err=True,
            )
            sys.exit(2)
        if opt_exc is not None and isinstance(e, opt_exc):
            msg = opt_msg or "ERROR: Optimization failed - {exc}"
            click.echo(msg.format(exc=e), err=True)
            sys.exit(3)
        # User-input errors raised in the command body (after Click's own
        # option parsing) should render as a clean one-line message, not a
        # raw Python traceback. Click's UsageError/BadParameter/ClickException
        # carry a user-facing message and conventional exit code (2);
        # malformed --config YAML is likewise a user error. Generic
        # exceptions still fall through to the full traceback so real
        # internal bugs remain visible (no exception masking).
        if isinstance(e, click.exceptions.ClickException):
            e.show()
            sys.exit(e.exit_code)
        try:
            import yaml as _yaml
            _is_yaml_err = isinstance(e, _yaml.YAMLError)
        except ImportError:
            _is_yaml_err = False
        if _is_yaml_err:
            click.echo(
                f"Error: invalid YAML in a configuration/override file: {e}",
                err=True,
            )
            sys.exit(2)
        # Emit the full traceback once via the human-readable terminal echo
        # so genuine internal bugs stay visible without duplicate output.
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo(
            f"Unhandled error during {label}:\n" + textwrap.indent(tb, "  "),
            err=True,
        )
        sys.exit(1)
    finally:
        # DO NOT INLINE: per-subcommand wrapper; runs when subcommand invoked standalone (e.g. `pdb2reaction tsopt ...`). Both this and the all.py finally block are required: different invocation paths.
        # Release GPU memory (model + Hessian) after CLI command finishes
        # so that subsequent pipeline stages (e.g. tsopt → irc) don't OOM.
        # gc.collect() breaks cyclic refs inside torch.nn.Module.
        gc.collect()
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except ImportError:
            pass
