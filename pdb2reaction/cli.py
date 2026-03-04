# pdb2reaction/cli.py

from __future__ import annotations

import logging
import warnings

import click

from .advanced_help import (
    _configure_subcommand_help_visibility,
    _ensure_help_advanced_option,
)
from .bool_compat import normalize_bool_argv
from .default_group import DefaultGroup
from pdb2reaction import __version__

_LAZY_SUBCOMMANDS: dict[str, tuple[str, str, str]] = {
    "all": (".all", "cli", "End-to-end workflow (extract -> MEP -> TS -> IRC -> freq -> DFT)."),
    "scan": (".scan", "cli", "Run staged 1D scan with harmonic restraints."),
    "opt": (".opt", "cli", "Optimize one structure."),
    "path-opt": (".path_opt", "cli", "Optimize a reaction path segment."),
    "path-search": (".path_search", "cli", "Search reaction pathways recursively."),
    "tsopt": (".tsopt", "cli", "Optimize a transition-state candidate."),
    "freq": (".freq", "cli", "Run vibrational analysis and thermochemistry."),
    "irc": (".irc", "cli", "Run IRC integration from a TS geometry."),
    "trj2fig": (".trj2fig", "cli", "Plot energy profile from trajectory."),
    "add-elem-info": (".add_elem_info", "cli", "Repair/add PDB element columns."),
    "dft": (".dft", "cli", "Run single-point DFT."),
    "scan2d": (".scan2d", "cli", "Run 2D distance scan."),
    "scan3d": (".scan3d", "cli", "Run 3D distance scan."),
    "extract": (".extract", "cli", "Extract a binding pocket."),
    "fix-altloc": (".fix_altloc", "cli", "Resolve PDB alternate locations."),
    "energy-diagram": (".energy_diagram", "cli", "Draw energy diagrams from values."),
}

# Only the ``all`` subcommand is listed here because it uses Click's
# ``type=click.BOOL`` (value-style) booleans that cannot be auto-detected
# from ``is_bool_flag``.  For all other subcommands the ``DefaultGroup``
# in ``default_group.py`` inspects the Click command's parameters at
# runtime and auto-discovers ``is_bool_flag`` / ``BoolParamType`` options,
# so they do not need to be repeated in these manual registries.
# All subcommands now use native Click options. The parser-wrapper
# infrastructure below is kept empty for forward compatibility.
_COMMAND_BOOL_VALUE_OPTIONS: dict[str, frozenset[str]] = {
    "all": frozenset(
        {
            "--include-H2O",
            "--include-h2o",
            "--exclude-backbone",
            "--add-linkH",
            "--verbose",
            "--freeze-links",
            "--climb",
            "--dump",
            "--convert-files",
            "--refine-path",
            "--preopt",
            "--tsopt",
            "--thermo",
            "--dft",
            "--scan-one-based",
            "--scan-preopt",
            "--scan-endopt",
        }
    ),
}

# Manual toggle-option hints.  ``DefaultGroup._resolve_bool_options()``
# auto-detects toggle options from Click's ``is_bool_flag`` attribute,
# but entries here ensure correct normalization *before* the lazy
# subcommand is imported (needed for early argv rewriting).
_COMMAND_BOOL_TOGGLE_OPTIONS: dict[str, frozenset[str]] = {
    "all": frozenset({"--flatten", "--show-config", "--dry-run", "--resume"}),
    "trj2fig": frozenset({"--reverse-x"}),
    "add-elem-info": frozenset({"--overwrite"}),
    "scan": frozenset(
        {
            "--one-based",
            "--freeze-links",
            "--dump",
            "--convert-files",
            "--preopt",
            "--endopt",
        }
    ),
    "scan2d": frozenset({"--one-based", "--freeze-links", "--dump", "--convert-files", "--preopt"}),
    "scan3d": frozenset({"--one-based", "--freeze-links", "--dump", "--convert-files", "--preopt"}),
    "opt": frozenset(
        {
            "--one-based",
            "--freeze-links",
            "--convert-files",
            "--dump",
            "--show-config",
            "--dry-run",
        }
    ),
    "dft": frozenset({"--convert-files", "--show-config", "--dry-run"}),
    "tsopt": frozenset(
        {
            "--freeze-links",
            "--convert-files",
            "--flatten",
            "--dump",
            "--show-config",
            "--dry-run",
        }
    ),
    "path-opt": frozenset(
        {
            "--freeze-links",
            "--climb",
            "--dump",
            "--convert-files",
            "--preopt",
            "--fix-ends",
            "--show-config",
            "--dry-run",
        }
    ),
    "path-search": frozenset(
        {
            "--freeze-links",
            "--climb",
            "--dump",
            "--convert-files",
            "--preopt",
            "--align",
            "--show-config",
            "--dry-run",
        }
    ),
    "freq": frozenset(
        {
            "--freeze-links",
            "--convert-files",
            "--dump",
            "--show-config",
            "--dry-run",
        }
    ),
    "irc": frozenset(
        {
            "--forward",
            "--backward",
            "--freeze-links",
            "--convert-files",
            "--show-config",
            "--dry-run",
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
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--spec",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "scan2d": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--spec",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "scan3d": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--spec",
            "--csv",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "opt": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--opt-mode",
            "--config",
            "--dry-run",
            "--out-dir",
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
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--config",
            "--dry-run",
            "--max-nodes",
            "--out-dir",
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
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--config",
            "--dry-run",
            "--max-nodes",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "tsopt": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--config",
            "--dry-run",
            "--max-cycles",
            "--opt-mode",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "freq": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--temperature",
            "--pressure",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "irc": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--max-cycles",
            "--step-size",
            "--forward",
            "--backward",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "dft": frozenset(
        {
            "-i",
            "--input",
            "-q",
            "--charge",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "--func-basis",
            "--engine",
            "--config",
            "--dry-run",
            "--out-dir",
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
}

_PARSER_WRAPPER_SUBCOMMANDS: frozenset[str] = frozenset()

_PARSER_WRAPPER_BOOL_OPTION_PROVIDERS: dict[str, object] = {}

_DEFAULT_GROUP_KWARGS = {
    "command_bool_value_options": _COMMAND_BOOL_VALUE_OPTIONS,
    "command_bool_toggle_options": _COMMAND_BOOL_TOGGLE_OPTIONS,
    "command_bool_toggle_negative_aliases": _COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES,
    "parser_wrapper_subcommands": _PARSER_WRAPPER_SUBCOMMANDS,
    "parser_wrapper_bool_option_providers": _PARSER_WRAPPER_BOOL_OPTION_PROVIDERS,
    "subcommand_primary_help_options": _SUBCOMMAND_PRIMARY_HELP_OPTIONS,
    "normalize_bool_argv": normalize_bool_argv,
    "ensure_help_advanced_option": _ensure_help_advanced_option,
    "configure_subcommand_help_visibility": _configure_subcommand_help_visibility,
}


@click.group(
    cls=DefaultGroup,
    default="all",
    lazy_subcommands=_LAZY_SUBCOMMANDS,
    **_DEFAULT_GROUP_KWARGS,
    help="pdb2reaction: Run workflow steps via subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=__version__, prog_name="pdb2reaction")
@click.pass_context
def cli(ctx: click.Context) -> None:
    if not ctx.resilient_parsing:
        click.echo(f"pdb2reaction ver. {__version__}\n")


# Silence pysisyphus logger without muting application/global logging.
_pysisyphus_logger = logging.getLogger("pysisyphus")
_pysisyphus_logger.setLevel(logging.CRITICAL)
_pysisyphus_logger.propagate = False

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
