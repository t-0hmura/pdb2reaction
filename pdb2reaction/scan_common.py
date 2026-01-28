# pdb2reaction/scan_common.py

from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Tuple, Any

import click


def add_scan_common_options(
    *,
    workers_default: int,
    workers_per_node_default: int,
    out_dir_default: str,
    baseline_help: str,
    dump_help: str,
    max_step_help: str = "Maximum step size in either distance [Å].",
    thresh_default: str | None = "baker",
    max_step_size_default: float = 0.20,
    bias_k_default: float = 100.0,
    relax_max_cycles_default: int = 10000,
    opt_mode_default: str = "light",
    freeze_links_default: bool = True,
    dump_default: bool = False,
    convert_files_default: bool = True,
    preopt_default: bool = True,
    one_based_default: bool = True,
    include_baseline: bool = True,
    include_zmin_zmax: bool = True,
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
            help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
        ),
        click.option(
            "--workers-per-node",
            "workers_per_node",
            type=int,
            default=workers_per_node_default,
            show_default=True,
            help="Workers per node when using a parallel UMA predictor (workers>1).",
        ),
        click.option(
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
            default=1,
            show_default=True,
            help="Spin multiplicity (2S+1) for the ML region.",
        ),
        click.option(
            "--one-based",
            "one_based",
            type=click.BOOL,
            default=one_based_default,
            show_default=True,
            help="Interpret (i,j) indices in --scan-list as 1-based (default) or 0-based.",
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
            show_default=True,
            help="Harmonic well strength k [eV/Å^2].",
        ),
        click.option(
            "--relax-max-cycles",
            type=int,
            default=relax_max_cycles_default,
            show_default=True,
            help=(
                "Maximum optimizer cycles per grid relaxation. When explicitly provided, "
                "overrides opt.max_cycles from YAML."
            ),
        ),
        click.option(
            "--opt-mode",
            type=click.Choice(["light", "heavy"], case_sensitive=False),
            default=opt_mode_default,
            show_default=True,
            help="Relaxation mode: light (=LBFGS) or heavy (=RFO).",
        ),
        click.option(
            "--freeze-links",
            type=click.BOOL,
            default=freeze_links_default,
            show_default=True,
            help="If input is PDB, freeze parent atoms of link hydrogens.",
        ),
        click.option(
            "--dump",
            type=click.BOOL,
            default=dump_default,
            show_default=True,
            help=dump_help,
        ),
        click.option(
            "--convert-files",
            "convert_files",
            type=click.BOOL,
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
            "--out-dir",
            type=str,
            default=out_dir_default,
            show_default=True,
            help="Base output directory.",
        ),
        click.option(
            "--thresh",
            type=str,
            default=thresh_default,
            show_default=False,
            help=(
                "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
                f"{thresh_note}"
            ),
        ),
        click.option(
            "--args-yaml",
            type=click.Path(path_type=Path, exists=True, dir_okay=False),
            default=None,
            help="YAML file with extra args (sections: geom, calc, opt, lbfgs, rfo, bias).",
        ),
        click.option(
            "--preopt",
            type=click.BOOL,
            default=preopt_default,
            show_default=True,
            help="Pre-optimize the initial structure without bias before the scan.",
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


def build_scan_defaults(
    *,
    geom_kw_default: Dict[str, Any],
    calc_kw_default: Dict[str, Any],
    opt_base_kw: Dict[str, Any],
    lbfgs_kw: Dict[str, Any],
    rfo_kw: Dict[str, Any],
    out_dir: str,
    thresh: str = "baker",
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Return (geom_kw, calc_kw, opt_kw, lbfgs_kw, rfo_kw) with shared scan defaults applied."""
    geom_kw = dict(geom_kw_default)
    calc_kw = dict(calc_kw_default)
    opt_cfg = dict(opt_base_kw)
    opt_cfg.update(
        {
            "out_dir": out_dir,
            "dump": False,
            "max_cycles": 10000,
            "thresh": thresh,
        }
    )
    lbfgs_cfg = dict(lbfgs_kw)
    lbfgs_cfg.update({"out_dir": out_dir})
    rfo_cfg = dict(rfo_kw)
    rfo_cfg.update({"out_dir": out_dir})
    return geom_kw, calc_kw, opt_cfg, lbfgs_cfg, rfo_cfg
