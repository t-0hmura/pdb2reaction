# pdb2reaction/cli.py

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
    "init": (".init", "cli", "Generate starter YAML templates."),
    "scan2d": (".scan2d", "cli", "Run 2D distance scan."),
    "scan3d": (".scan3d", "cli", "Run 3D distance scan."),
    "fix-altloc": (".fix_altloc", "cli", "Resolve PDB alternate locations."),
    "energy-diagram": (".energy_diagram", "cli", "Draw energy diagrams from values."),
}

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
            "--flatten-imag-mode",
            "--scan-one-based",
            "--scan-preopt",
            "--scan-endopt",
        }
    ),
}

_COMMAND_BOOL_TOGGLE_OPTIONS: dict[str, frozenset[str]] = {
    "extract": frozenset(
        {
            "--include-H2O",
            "--include-h2o",
            "--exclude-backbone",
            "--add-linkH",
            "--verbose",
        }
    ),
    "trj2fig": frozenset({"--reverse-x"}),
    "add-elem-info": frozenset({"--overwrite"}),
    "fix-altloc": frozenset({"--recursive", "--inplace", "--overwrite", "--force"}),
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
            "--flatten-imag-mode",
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

_PARSER_WRAPPER_SUBCOMMANDS = frozenset({"extract", "fix-altloc"})

_DEFAULT_GROUP_KWARGS = {
    "command_bool_value_options": _COMMAND_BOOL_VALUE_OPTIONS,
    "command_bool_toggle_options": _COMMAND_BOOL_TOGGLE_OPTIONS,
    "command_bool_toggle_negative_aliases": _COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES,
    "parser_wrapper_subcommands": _PARSER_WRAPPER_SUBCOMMANDS,
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
def cli() -> None:
    pass


@click.command(
    name="extract",
    help="Extract a binding pocket.",
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "help_option_names": [],
    },
)
@click.pass_context
def extract_cmd(ctx: click.Context) -> None:
    from . import extract as _extract_mod
    args = _extract_mod.parse_args(list(ctx.args))
    _extract_mod.extract(args)


cli.add_command(extract_cmd, name="extract")

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
