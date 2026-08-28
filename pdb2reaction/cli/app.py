from __future__ import annotations

import logging
import shlex
import sys
import warnings
from pathlib import Path

import click

from pdb2reaction.cli.help_pages import (
    _configure_subcommand_help_visibility,
    _ensure_help_advanced_option,
)
from pdb2reaction.cli.bool_compat import _parse_bool_literal, normalize_bool_argv
from pdb2reaction.cli.default_group import DefaultGroup
from pdb2reaction import __version__

_PDB2REACTION_BANNER = r"""
 ____  ____  ____  ____  ____                 _   _
|  _ \|  _ \| __ ||___ \|  _ \ ___  __ _  ___| |_|_| ___  _ __
| |_| | | | |  _ \  __| | |_| / _ \/ _` |/ __| __| |/ _ \| '_ \
|  __/| |_| | |_| |/ __/|  _ <  __/ |_| | |__| |_| | |_| | | | |
|_|   |____/|____/|_____|_| \_\___|\__,_|\___|\__|_|\___/|_| |_|
""".strip("\n")

_CONSOLE_SCRIPT_NAMES = {"pdb2reaction", "p2r"}


def _command_argv(argv: list[str]) -> list[str]:
    if not argv:
        return []
    argv0_name = Path(argv[0]).name
    if argv0_name == "__main__.py":
        return [sys.executable, "-m", "pdb2reaction", *argv[1:]]
    if argv0_name in _CONSOLE_SCRIPT_NAMES:
        return [argv0_name, *argv[1:]]
    return list(argv)


def _quoted_argv(argv: list[str]) -> str:
    """Return a shell-safe representation of the executed argv."""
    return shlex.join([str(arg) for arg in _command_argv(argv)])


def _has_help_or_version_request(argv: list[str]) -> bool:
    return any(arg in {"-h", "--help", "--version", "--help-advanced"} for arg in argv[1:])


def _toggle_enabled(args: list[str], flag: str) -> bool:
    """Resolve a `--flag` / `--no-flag` pair, accepting the legacy value style."""
    enabled = False
    i = 0
    while i < len(args):
        name, separator, value = args[i].partition("=")
        name = name.lower()
        if name == f"--no-{flag}":
            enabled = False
        elif name == f"--{flag}":
            parsed = _parse_bool_literal(value) if separator else None
            if parsed is None and not separator and i + 1 < len(args):
                parsed = _parse_bool_literal(args[i + 1])
                if parsed is not None:
                    i += 1
            enabled = True if parsed is None else parsed
        i += 1
    return enabled


def _requests_stdout_json(argv: list[str]) -> bool:
    """Return whether bond-summary promises JSON-only stdout."""
    args = argv[1:]
    if "bond-summary" not in args:
        return False
    return _toggle_enabled(args, "json")


def _requests_dry_run(argv: list[str]) -> bool:
    """Return whether the run only reports what it would do."""
    return _toggle_enabled(argv[1:], "dry-run")


def _emit_start_header(ctx: click.Context) -> None:
    from pdb2reaction.core.utils import emit, is_child_mode, verbose_level

    if (
        is_child_mode()
        or _has_help_or_version_request(sys.argv)
        or _requests_stdout_json(sys.argv)
    ):
        return

    # A dry run reports a plan, so keep the version line but drop the artwork.
    if verbose_level() >= 2 and not _requests_dry_run(sys.argv):
        emit(f"{_PDB2REACTION_BANNER}\n\npdb2reaction ver. {__version__}\n", narrative=True)
    else:
        emit(f"pdb2reaction ver. {__version__}\n", narrative=True)

    subcommand = (
        getattr(ctx, "invoked_subcommand", None)
        or getattr(ctx, "info_name", None)
        or (ctx.command.name if ctx.command is not None else None)
        or "all"
    )
    if verbose_level() >= 2:
        emit(f"[command] {_quoted_argv(sys.argv)}", narrative=True, raw_path=True)
        if subcommand != "all":
            emit(f"[mode] {subcommand}", narrative=True)
    emit("", narrative=True)


_LAZY_SUBCOMMANDS: dict[str, tuple[str, str, str]] = {
    "all": ("pdb2reaction.workflows.all", "cli", "End-to-end reaction workflow with MEP, scan, or TS-only entry routes."),
    "scan": ("pdb2reaction.workflows.scan", "cli", "Run staged 1D scan with harmonic restraints."),
    "opt": ("pdb2reaction.workflows.opt", "cli", "Optimize one structure."),
    "path-opt": ("pdb2reaction.workflows.path_opt", "cli", "Optimize a reaction path segment."),
    "path-search": ("pdb2reaction.workflows.path_search", "cli", "Search reaction pathways recursively."),
    "tsopt": ("pdb2reaction.workflows.tsopt", "cli", "Optimize a transition-state candidate."),
    "freq": ("pdb2reaction.workflows.freq", "cli", "Run vibrational analysis and thermochemistry."),
    "irc": ("pdb2reaction.workflows.irc", "cli", "Run IRC integration from a TS geometry."),
    "trj2fig": ("pdb2reaction.io.trj2fig", "cli", "Plot energy profile from trajectory."),
    "add-elem-info": ("pdb2reaction.domain.add_elem_info", "cli", "Repair/add PDB element columns."),
    "dft": ("pdb2reaction.workflows.dft", "cli", "Run single-point DFT."),
    "sp": ("pdb2reaction.workflows.sp", "cli", "Run single-point MLIP energy + forces."),
    "scan2d": ("pdb2reaction.workflows.scan2d", "cli", "Run 2D distance scan."),
    "scan3d": ("pdb2reaction.workflows.scan3d", "cli", "Run 3D distance scan."),
    "extract": ("pdb2reaction.workflows.extract", "cli", "Extract an active site model."),
    "fix-altloc": ("pdb2reaction.io.pdb_fix", "cli", "Resolve PDB alternate locations."),
    "energy-diagram": ("pdb2reaction.io.energy_diagram", "cli", "Draw energy diagrams from values."),
    "bond-summary": ("pdb2reaction.domain.bond_summary", "cli", "Detect bond changes between structures."),
}

# ``all`` retains value-style Click booleans for its established interface;
# toggle-style options on other commands are discovered at runtime.
_COMMAND_BOOL_VALUE_OPTIONS: dict[str, frozenset[str]] = {
    "all": frozenset(
        {
            "--add-linkh",
            "--climb",
            "--convert-files",
            "--dft",
            "--dump",
            "--exclude-backbone",
            "--freeze-links",
            "--include-h2o",
            "--preopt",
            "--refine-path",
            "--scan-endopt",
            "--scan-one-based",
            "--scan-preopt",
            "--thermo",
            "--tsopt",
            "--write-ref-merge",
        }
    ),
}

# Manual toggle-option hints.  ``DefaultGroup._resolve_bool_options()``
# auto-detects toggle options from Click's ``is_bool_flag`` attribute,
# but entries here ensure correct normalization *before* the lazy
# subcommand is imported (needed for early argv rewriting).
_COMMAND_BOOL_TOGGLE_OPTIONS: dict[str, frozenset[str]] = {
    "add-elem-info": frozenset(
        {
            "--overwrite",
        }
    ),
    "all": frozenset(
        {
            "--dry-run",
            "--flatten",
            "--show-config",
        }
    ),
    "bond-summary": frozenset(
        {
            "--json",
            "--one-based",
        }
    ),
    "dft": frozenset(
        {
            "--dry-run",
            "--lowmem",
            "--out-json",
            "--show-config",
        }
    ),
    "energy-diagram": frozenset(
        {
            "--out-json",
        }
    ),
    "extract": frozenset(
        {
            "--add-linkh",
            "--exclude-backbone",
            "--include-h2o",
            "--out-json",
        }
    ),
    "fix-altloc": frozenset(
        {
            "--force",
            "--inplace",
            "--overwrite",
            "--recursive",
        }
    ),
    "freq": frozenset(
        {
            "--convert-files",
            "--dry-run",
            "--dump",
            "--freeze-links",
            "--out-json",
            "--show-config",
        }
    ),
    "irc": frozenset(
        {
            "--backward",
            "--convert-files",
            "--dry-run",
            "--forward",
            "--freeze-links",
            "--never-stop",
            "--out-json",
            "--show-config",
        }
    ),
    "opt": frozenset(
        {
            "--convert-files",
            "--dry-run",
            "--dump",
            "--flatten",
            "--freeze-links",
            "--one-based",
            "--out-json",
            "--show-config",
        }
    ),
    "path-opt": frozenset(
        {
            "--climb",
            "--convert-files",
            "--dry-run",
            "--dump",
            "--fix-ends",
            "--freeze-links",
            "--out-json",
            "--preopt",
            "--show-config",
        }
    ),
    "path-search": frozenset(
        {
            "--align",
            "--climb",
            "--convert-files",
            "--dry-run",
            "--dump",
            "--freeze-links",
            "--preopt",
            "--show-config",
            "--write-ref-merge",
        }
    ),
    "scan": frozenset(
        {
            "--convert-files",
            "--dry-run",
            "--dump",
            "--endopt",
            "--freeze-links",
            "--one-based",
            "--out-json",
            "--preopt",
            "--print-parsed",
        }
    ),
    "scan2d": frozenset(
        {
            "--convert-files",
            "--dry-run",
            "--dump",
            "--freeze-links",
            "--one-based",
            "--out-json",
            "--preopt",
            "--print-parsed",
        }
    ),
    "scan3d": frozenset(
        {
            "--convert-files",
            "--dry-run",
            "--dump",
            "--freeze-links",
            "--one-based",
            "--out-json",
            "--preopt",
            "--print-parsed",
        }
    ),
    "trj2fig": frozenset(
        {
            "--out-json",
            "--reverse-x",
        }
    ),
    "tsopt": frozenset(
        {
            "--convert-files",
            "--dry-run",
            "--dump",
            "--flatten",
            "--freeze-links",
            "--out-json",
            "--show-config",
        }
    ),
}

_COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES: dict[str, dict[str, str]] = {
    "scan": {
        "--one-based": "--zero-based",
    },
    "scan2d": {
        "--one-based": "--zero-based",
    },
    "scan3d": {
        "--one-based": "--zero-based",
    },
    "opt": {
        "--one-based": "--zero-based",
    },
}

_SUBCOMMAND_PRIMARY_HELP_OPTIONS: dict[str, frozenset[str]] = {
    "scan": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "-s", "--scan-lists",
            "--config",
            "-o", "--out-dir",
            "--bias-k",
            "--max-step-size",
            "--thresh",
            "--print-every",
            "--help-advanced",
        }
    ),
    "scan2d": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "-s", "--scan-lists",
            "--config",
            "-o", "--out-dir",
            "--bias-k",
            "--max-step-size",
            "--thresh",
            "--print-every",
            "--help-advanced",
        }
    ),
    "scan3d": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "-s", "--scan-lists",
            "--csv",
            "--config",
            "-o", "--out-dir",
            "--bias-k",
            "--max-step-size",
            "--thresh",
            "--print-every",
            "--help-advanced",
        }
    ),
    "opt": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--opt-mode",
            "--thresh",
            "--dump", "--no-dump",
            "--dist-freeze",
            "--bias-k",
            "--config",
            "-o", "--out-dir",
            "--max-cycles",
            "--print-every",
            "--help-advanced",
        }
    ),
    "path-opt": frozenset(
        {
            "-i",
            "--input",
            "--mep-mode",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--config",
            "--max-nodes",
            "--fix-ends",
            "-o", "--out-dir",
            "--max-cycles-gsm",
            "--max-cycles-dmf",
            "--help-advanced",
        }
    ),
    "path-search": frozenset(
        {
            "-i",
            "--input",
            "--mep-mode",
            "--refine-mode",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--config",
            "--max-nodes",
            "-o", "--out-dir",
            "--max-cycles-gsm",
            "--max-cycles-dmf",
            "--help-advanced",
        }
    ),
    "tsopt": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--config",
            "--max-cycles",
            "--opt-mode",
            "-o", "--out-dir",
            "--hessian-calc-mode",
            "--thresh",
            "--dump", "--no-dump",
            "--print-every",
            "--help-advanced",
        }
    ),
    "freq": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--temperature",
            "--pressure",
            "--config",
            "-o", "--out-dir",
            "--hessian-calc-mode",
            "--help-advanced",
        }
    ),
    "irc": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--max-cycles",
            "--step-size",
            "--forward",
            "--backward",
            "--config",
            "-o", "--out-dir",
            "--hessian-calc-mode",
            "--help-advanced",
        }
    ),
    "dft": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--func-basis",
            "--engine",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "sp": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--hess",
            "--hessian-calc-mode",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "trj2fig": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--out",
            "--unit",
            "-r",
            "--reference",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--reverse-x", "--no-reverse-x",
            "--help-advanced",
        }
    ),
    "add-elem-info": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--out",
            "--help-advanced",
        }
    ),
    "energy-diagram": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--output",
            "--help-advanced",
        }
    ),
    "extract": frozenset(
        {
            "-i",
            "--input",
            "-c",
            "--center",
            "-o",
            "--output",
            "-r",
            "--radius",
            "--selected-resn",
            "-l",
            "--ligand-charge",
            "--help-advanced",
        }
    ),
    "fix-altloc": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--out",
            "--recursive", "--no-recursive",
            "--inplace", "--no-inplace",
            "--help-advanced",
        }
    ),
    "bond-summary": frozenset(
        {
            "-i",
            "--input",
            "--device",
            "--bond-factor",
            "--help-advanced",
        }
    ),
}

_PARSER_WRAPPER_SUBCOMMANDS: frozenset[str] = frozenset()

_PARSER_WRAPPER_BOOL_OPTION_PROVIDERS: dict[str, object] = {}

_COMMAND_GROUPS: dict[str, tuple[str, ...]] = {
    "Pipelines": ("all",),
    "Pipeline stages": (
        "opt", "sp", "tsopt", "freq", "irc", "dft",
        "scan", "scan2d", "scan3d", "path-opt", "path-search",
    ),
    "Inputs & topology": ("extract", "fix-altloc", "add-elem-info"),
    "Analysis": ("energy-diagram", "bond-summary", "trj2fig"),
}

_COMMAND_BOOL_SINGLE_FLAG_OPTIONS: dict[str, frozenset[str]] = {}

_DEFAULT_GROUP_KWARGS = {
    "command_bool_value_options": _COMMAND_BOOL_VALUE_OPTIONS,
    "command_bool_toggle_options": _COMMAND_BOOL_TOGGLE_OPTIONS,
    "command_bool_toggle_negative_aliases": _COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES,
    "command_bool_single_flag_options": _COMMAND_BOOL_SINGLE_FLAG_OPTIONS,
    "parser_wrapper_subcommands": _PARSER_WRAPPER_SUBCOMMANDS,
    "parser_wrapper_bool_option_providers": _PARSER_WRAPPER_BOOL_OPTION_PROVIDERS,
    "subcommand_primary_help_options": _SUBCOMMAND_PRIMARY_HELP_OPTIONS,
    "normalize_bool_argv": normalize_bool_argv,
    "ensure_help_advanced_option": _ensure_help_advanced_option,
    "configure_subcommand_help_visibility": _configure_subcommand_help_visibility,
    "command_groups": _COMMAND_GROUPS,
}


def _verbose_callback(ctx: click.Context, param: click.Parameter, value: int) -> int:
    # `-v/--verbose LEVEL` (0-3) is the single, unified verbosity control:
    #   0 silent · 1 milestones · 2 (default) +cycle tables/timing/VRAM/paths
    #   · 3 everything (config dumps, per-file paths, DEBUG logging).
    # Injected into every subcommand (see `_ensure_verbose_option`), so this
    # eager callback is where each invocation sets its level, enables the
    # console gate, and emits the start header once the level is known.
    from pdb2reaction.core.utils import set_verbose_level, set_console_gating
    set_verbose_level(value)
    set_console_gating(True)
    if value >= 3:
        from pdb2reaction.core.logging import setup_logging
        setup_logging(value)
    _emit_start_header(ctx)
    return value


# Subcommands that run an iterative optimizer, so their -v 2 genuinely adds
# per-cycle optimizer tables. Other subcommands are single-shot compute or
# lightweight IO/analysis steps. A backend may still emit its own Hessian VRAM
# detail for any command that evaluates a Hessian, including freq and sp.
_OPTIMIZER_SUBCOMMANDS = frozenset({
    "all", "opt", "tsopt", "irc",
    "scan", "scan2d", "scan3d",
    "path-opt", "path-search",
})


def _verbose_help(name: str | None) -> str:
    """`-v/--verbose` help text, specialised per subcommand so the level-2
    description never claims optimizer cycle tables for non-optimizer commands
    (extract, sp, freq, dft, file IO, analysis)."""
    if name in _OPTIMIZER_SUBCOMMANDS:
        level2 = "optimizer cycle tables, per-stage timing, VRAM, deliverable paths"
    else:
        level2 = "detailed step logging and deliverable paths"
    return (
        "Console verbosity 0-3 (default 2). 0=silent; 1=milestones only; "
        f"2=+{level2}; "
        "3=everything (full config blocks, per-file paths, DEBUG logging)."
    )


def _ensure_verbose_option(
    command: click.Command, cmd_name: str | None = None
) -> click.Command:
    """Attach the unified `-v/--verbose LEVEL` control to a subcommand when absent.

    Verbosity is a per-subcommand option (not a root-group option) so it can be
    written in the natural position, e.g. ``pdb2reaction opt -v 0 ...``.  The
    eager callback sets the level, enables the console gate, and emits the start
    header. ``cmd_name`` is the registry key (``command.name`` is unreliable:
    lazy subcommands share the underlying ``cli`` function name) and selects the
    per-subcommand help text.
    """
    if any(
        isinstance(param, click.Option)
        and ("--verbose" in param.opts or "-v" in param.opts)
        for param in command.params
    ):
        return command

    option = click.Option(
        ["-v", "--verbose"],
        type=click.IntRange(0, 3),
        default=2,
        metavar="LEVEL",
        is_eager=True,
        expose_value=False,
        callback=_verbose_callback,
        help=_verbose_help(cmd_name),
    )
    command.params.insert(0, option)
    return command


@click.group(
    cls=DefaultGroup,
    default="all",
    lazy_subcommands=_LAZY_SUBCOMMANDS,
    **_DEFAULT_GROUP_KWARGS,
    ensure_verbose_option=_ensure_verbose_option,
    help="pdb2reaction: Run workflow steps via subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=__version__, prog_name="pdb2reaction")
@click.pass_context
def cli(ctx: click.Context) -> None:
    # `-v/--verbose` is a per-subcommand option (injected by
    # `_ensure_verbose_option`); the start header is emitted by that option's
    # eager callback once the level is known, not from this group body.
    if not ctx.resilient_parsing:
        from pdb2reaction.core.utils import set_base_dir
        set_base_dir(Path.cwd())


# Pysisyphus log suppression is handled by DefaultGroup._silence_pysisyphus_loggers()
# which runs after lazy subcommand import (when pysisyphus __init__ handlers are created).

# Filter noisy UMA/pydmf warnings that clutter CLI output
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"var\(\): degrees of freedom is <= 0\. Correction should be strictly less than the reduction factor.*",
    module=r"fairchem\.core\.models\.uma\.escn_moe"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"t_eval update skipped due to insufficient candidates",
    module=r"dmf"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"Setting global torch default dtype to torch\.float32\.",
)
# Suppress fairchem dataset_list deprecation warning (baked into UMA checkpoint
# config; cannot be fixed caller-side until fairchem removes dataset_list from
# checkpoints). The message is emitted from `escn_md.resolve_dataset_mapping`
# via the bare `logging.warning(...)` call, which routes through the root
# logger (visible as `WARNING:root:...`), not a `fairchem.*` named logger.
# Filter the root logger but require the verbatim opening clause so unrelated
# user log messages that happen to contain "dataset_list" survive.
_DATASET_LIST_MSG = "If 'dataset_list' is provided in the config"
logging.getLogger().addFilter(
    lambda record: _DATASET_LIST_MSG not in record.getMessage()
)
