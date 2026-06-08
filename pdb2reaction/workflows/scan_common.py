from __future__ import annotations

import logging
from pathlib import Path
from typing import Callable, Dict, Any

import click

from pdb2reaction.core.defaults import THRESH_CHOICES
from pdb2reaction.core.utils import deep_update, load_yaml_dict

logger = logging.getLogger(__name__)


def add_scan_common_options(
    *,
    workers_default: int,
    workers_per_node_default: int,
    out_dir_default: str,
    baseline_help: str,
    dump_help: str,
    max_step_help: str = "Maximum step size per scanned distance [Å].",
    thresh_default: str | None = "baker",
    max_step_size_default: float = 0.20,
    bias_k_default: float | None = None,
    relax_max_cycles_default: int = 10000,
    opt_mode_default: str = "grad",
    freeze_links_default: bool = True,
    dump_default: bool = False,
    convert_files_default: bool = True,
    preopt_default: bool = False,
    one_based_default: bool = True,
    include_baseline: bool = True,
    include_zmin_zmax: bool = True,
    args_yaml_sections: str = "geom, calc, opt, lbfgs, rfo, bias",
) -> Callable[[Callable], Callable]:
    """Attach the shared scan2d/scan3d CLI options to a Click command."""
    thresh_note = f" Defaults to '{thresh_default}'." if thresh_default is not None else ""
    options = [
        click.option(
            "-q",
            "--charge",
            type=int,
            required=False,
            help=(
                "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
                "(PDB inputs or XYZ/GJF with --ref-pdb)."
            ),
        ),
        click.option(
            "--workers",
            type=int,
            default=workers_default,
            show_default=True,
            help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: the analytical Hessian raises a RuntimeError when workers>1; run with --workers 1 for Hessian-based modes.",
        ),
        click.option(
            "--workers-per-node",
            "workers_per_node",
            type=int,
            default=workers_per_node_default,
            show_default=True,
            help="Workers per node when using a parallel MLIP predictor (workers>1).",
        ),
        click.option(
            "-l",
            "--ligand-charge",
            type=str,
            default=None,
            show_default=False,
            help=(
                "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
                "when -q is omitted (requires PDB input or --ref-pdb)."
            ),
        ),
        click.option(
            "-m",
            "--multiplicity",
            "spin",
            type=int,
            default=None,
            show_default=False,
            help="Spin multiplicity (2S+1).",
        ),
        click.option(
            "--one-based/--zero-based",
            "one_based",
            default=one_based_default,
            show_default=True,
            help="Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.",
        ),
        click.option(
            "--max-step-size",
            type=float,
            default=max_step_size_default,
            show_default=True,
            help=max_step_help,
        ),
        click.option(
            "--bias-k",
            type=float,
            default=bias_k_default,
            show_default=False,
            help=(
                "Harmonic well strength k [eV/Å^2]. "
                "Defaults to YAML bias.k (BIAS_KW['k']=300 in defaults.py) when omitted; "
                "explicit CLI value overrides YAML."
            ),
        ),
        click.option(
            "--relax-max-cycles",
            type=int,
            default=relax_max_cycles_default,
            show_default=True,
            help=(
                "Maximum optimizer cycles per grid relaxation. When explicitly provided, "
                "used unless YAML sets opt.max_cycles."
            ),
        ),
        click.option(
            "--opt-mode",
            type=click.Choice(["grad", "hess"], case_sensitive=False),
            default=opt_mode_default,
            show_default=True,
            help="Relaxation mode: grad (=LBFGS) or hess (=RFO).",
        ),
        click.option(
            "--freeze-links/--no-freeze-links",
            "freeze_links",
            default=freeze_links_default,
            show_default=True,
            help="Freeze parent atoms of link hydrogens (PDB input or XYZ/GJF with --ref-pdb).",
        ),
        click.option(
            "--freeze-atoms",
            "freeze_atoms_text",
            type=str,
            default=None,
            show_default=False,
            help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
        ),
        click.option(
            "--dump/--no-dump",
            "dump",
            default=dump_default,
            show_default=True,
            help=dump_help,
        ),
        click.option(
            "--convert-files/--no-convert-files",
            "convert_files",
            default=convert_files_default,
            show_default=True,
            help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
        ),
        click.option(
            "--ref-pdb",
            type=click.Path(path_type=Path, exists=True, dir_okay=False),
            default=None,
            help="Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
        ),
        click.option(
            "-o",
            "--out-dir",
            type=str,
            default=out_dir_default,
            show_default=True,
            help="Base output directory.",
        ),
        click.option(
            "--thresh",
            type=click.Choice(THRESH_CHOICES, case_sensitive=False),
            default=thresh_default,
            show_default=False,
            help=(
                "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
                f"{thresh_note}"
            ),
        ),
        click.option(
            "--config",
            "config_yaml",
            type=click.Path(path_type=Path, exists=True, dir_okay=False),
            default=None,
            show_default=False,
            help="Base YAML configuration file applied before explicit CLI options.",
        ),
        click.option(
            "--preopt/--no-preopt",
            "preopt",
            default=preopt_default,
            show_default=True,
            help="Pre-optimize the initial structure without bias before the scan.",
        ),
        click.option(
            "-b", "--backend",
            type=click.Choice(["uma", "orb", "mace", "aimnet2"]),
            default="uma",
            show_default=True,
            help="MLIP backend.",
        ),
        click.option(
            "--solvent",
            default="none",
            show_default=True,
            help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.",
        ),
        click.option(
            "--solvent-model",
            "solvent_model",
            default="alpb",
            type=click.Choice(["alpb", "cpcmx"]),
            show_default=True,
            help="xTB solvent model.",
        ),
    ]
    if include_baseline:
        options.append(
            click.option(
                "--baseline",
                type=click.Choice(["min", "first"]),
                default="min",
                show_default=True,
                help=baseline_help,
            )
        )
    if include_zmin_zmax:
        options.extend(
            [
                click.option(
                    "--zmin",
                    type=float,
                    default=None,
                    show_default=False,
                    help="Lower bound of color scale for plots (kcal/mol).",
                ),
                click.option(
                    "--zmax",
                    type=float,
                    default=None,
                    show_default=False,
                    help="Upper bound of color scale for plots (kcal/mol).",
                ),
            ]
        )

    def decorator(func):
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator



def resolve_yaml_sources(
    config_yaml: Path | None,
    override_yaml: Path | None,
    args_yaml_legacy: Path | None,
) -> tuple[Path | None, Path | None, bool]:
    """Resolve YAML layers and legacy alias usage for scan-family commands."""
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use a single YAML source option."
        )
    if args_yaml_legacy is not None:
        return config_yaml, args_yaml_legacy, True
    return config_yaml, override_yaml, False


def load_merged_yaml_cfg(
    config_yaml: Path | None,
    override_yaml: Path | None,
) -> Dict[str, Any]:
    """Load and deep-merge scan YAML layers as config < override."""
    merged: Dict[str, Any] = {}
    deep_update(merged, load_yaml_dict(config_yaml))
    deep_update(merged, load_yaml_dict(override_yaml))
    return merged
