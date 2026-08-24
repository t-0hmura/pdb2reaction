from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Any, Callable, Sequence

import click

from pdb2reaction.core.defaults import THRESH_CHOICES
from pdb2reaction.core.utils import (
    is_scan_spec_file,
    parse_scan_list_quads_checked,
    parse_scan_list_triples,
    parse_scan_spec_quads,
    parse_scan_spec_stages,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class StagedScanRequest:
    """One fully parsed 1D scan request shared by planning and execution."""

    stages: tuple[tuple[tuple[int, int, float], ...], ...]
    one_based: bool
    source: str
    raw_values: tuple[str, ...]
    bidirectional_reset_before: frozenset[int] = frozenset()
    bidirectional_snapshot_before: frozenset[int] = frozenset()


@dataclass(frozen=True)
class GridScanRequest:
    """One fully parsed 2D/3D scan request shared by planning and execution."""

    pairs: tuple[tuple[int, int, float, float], ...]
    raw_pairs: tuple[tuple[Any, Any, float, float], ...]
    one_based: bool
    source: str
    raw_value: str


def collect_staged_scan_values(
    option_values: Sequence[str],
    extra_args: Sequence[str],
    *,
    option_name: str = "--scan-lists",
) -> tuple[str, ...]:
    """Recover the legacy one-flag/many-values form without reading ``sys.argv``.

    Click consumes the first value after a ``--scan-lists`` option and leaves
    subsequent positional stage literals in ``Context.args``. Both the legacy
    grouped form and Click's natural repeated-option form describe stages.
    """

    values = tuple(str(value) for value in option_values)
    extras = tuple(str(value) for value in extra_args)
    unexpected = next((value for value in extras if value.startswith("-")), None)
    if unexpected is not None:
        raise click.BadParameter(f"Unexpected option or argument: {unexpected}")
    if extras and len(values) > 1:
        raise click.BadParameter(
            f"Do not mix repeated {option_name} options with the grouped "
            "one-option/many-values form."
        )
    combined = values + extras
    if not combined:
        raise click.BadParameter(f"{option_name} is required.")
    return combined


def parse_staged_scan_request(
    raw_values: Sequence[str],
    *,
    one_based: bool,
    atom_meta: Sequence[dict[str, Any]],
    option_name: str = "--scan-lists",
) -> StagedScanRequest:
    """Parse one complete staged request before dry-run or execution branches."""

    values = tuple(str(value) for value in raw_values)
    if not values:
        raise click.BadParameter(f"{option_name} is required.")

    scan_one_based = bool(one_based)
    source = option_name
    reset_before: set[int] = set()
    snapshot_before: set[int] = set()
    stages: list[list[tuple[int, int, float]]]
    if len(values) == 1 and is_scan_spec_file(values[0]):
        spec_path = Path(values[0])
        (
            stages,
            scan_one_based,
            snapshot_before_spec,
            reset_before_spec,
        ) = parse_scan_spec_stages(
            spec_path,
            one_based_default=one_based,
            atom_meta=atom_meta,
            option_name=option_name,
            return_bidirectional_markers=True,
        )
        snapshot_before.update(snapshot_before_spec)
        reset_before.update(reset_before_spec)
        source = f"{option_name} ({spec_path})"
    else:
        stages = []
        for value_index, raw in enumerate(values, start=1):
            parsed, _ = parse_scan_list_triples(
                raw,
                one_based=scan_one_based,
                atom_meta=atom_meta,
                option_name=f"{option_name} #{value_index}",
            )
            for entry in parsed:
                if any(float(distance) <= 0.0 for distance in entry[2:]):
                    raise click.BadParameter(
                        f"Non-positive target length in {option_name} #{value_index}: {entry}."
                    )
            if any(len(entry) == 4 for entry in parsed):
                for entry in parsed:
                    if len(entry) == 4:
                        i, j, start, end = entry
                        stage_index = len(stages)
                        stages.append([(i, j, start)])
                        snapshot_before.add(stage_index)
                        reset_before.add(stage_index + 1)
                        stages.append([(i, j, end)])
                    else:
                        stages.append([entry])
            else:
                stages.append(parsed)

    return StagedScanRequest(
        stages=tuple(tuple(stage) for stage in stages),
        one_based=scan_one_based,
        source=source,
        raw_values=values,
        bidirectional_reset_before=frozenset(reset_before),
        bidirectional_snapshot_before=frozenset(snapshot_before),
    )


def parse_grid_scan_request(
    raw_value: str | None,
    *,
    dimensions: int,
    one_based: bool,
    atom_meta: Sequence[dict[str, Any]],
    option_name: str = "--scan-lists",
) -> GridScanRequest:
    """Parse one literal or YAML/JSON grid request exactly once."""

    if raw_value is None:
        raise click.BadParameter(f"{option_name} is required.")
    raw = str(raw_value)
    scan_one_based = bool(one_based)
    source = option_name
    if is_scan_spec_file(raw):
        spec_path = Path(raw)
        parsed, raw_pairs, scan_one_based = parse_scan_spec_quads(
            spec_path,
            expected_len=dimensions,
            one_based_default=one_based,
            atom_meta=atom_meta,
            option_name=option_name,
        )
        source = f"{option_name} ({spec_path})"
    else:
        parsed, raw_pairs = parse_scan_list_quads_checked(
            raw,
            expected_len=dimensions,
            one_based=scan_one_based,
            atom_meta=atom_meta,
            option_name=option_name,
        )
    return GridScanRequest(
        pairs=tuple(parsed),
        raw_pairs=tuple(raw_pairs),
        one_based=scan_one_based,
        source=source,
        raw_value=raw,
    )


def add_scan_common_options(
    *,
    workers_default: int,
    workers_per_node_default: int,
    out_dir_default: str,
    baseline_help: str,
    dump_help: str,
    max_step_help: str = "Maximum step size per scanned distance [Å].",
    thresh_default: str | None = "baker",
    # Display-only: a command that resolves `--thresh` downstream must keep its
    # declared default None so `cli_param_overridden` still sees an omission,
    # while help/docs/the Colab option list still name the effective value.
    thresh_shown: str | None = None,
    max_step_size_default: float = 0.20,
    bias_k_default: float | None = None,
    # Display-only fallback when the command resolves --bias-k downstream.
    bias_k_shown: str | None = "300.0",
    relax_max_cycles_default: int | None = None,
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
    thresh_effective = thresh_shown if thresh_shown is not None else thresh_default
    options = [
        click.option(
            "-q",
            "--charge",
            type=int,
            required=False,
            help=(
                "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
                "(PDB/mmCIF inputs or XYZ/GJF with --ref-pdb)."
            ),
        ),
        click.option(
            "--workers",
            type=int,
            default=workers_default,
            show_default=True,
            help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
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
                "when -q is omitted (requires PDB/mmCIF input or --ref-pdb)."
            ),
        ),
        click.option(
            "-m",
            "--multiplicity",
            "spin",
            type=int,
            default=None,
            show_default="1",
            help="Spin multiplicity (2S+1).",
        ),
        click.option(
            "--one-based/--zero-based",
            "one_based",
            default=one_based_default,
            show_default=True,
            help="Interpret (i,j) indices in --scan-lists as 1-based or 0-based.",
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
            show_default=(bias_k_default if bias_k_default is not None else bias_k_shown),
            help=(
                "Harmonic well strength k [eV/Å^2]. "
                "YAML bias.k applies when this option is omitted; explicit CLI wins."
            ),
        ),
        click.option(
            "--relax-max-cycles",
            type=click.IntRange(min=1),
            default=relax_max_cycles_default,
            show_default="100000",
            help=(
                "Maximum optimizer cycles per grid relaxation. An explicitly "
                "provided value overrides YAML opt.max_cycles."
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
            help="Freeze parent atoms of cap hydrogens (PDB/mmCIF input or XYZ/GJF with --ref-pdb).",
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
            help="Convert XYZ/TRJ outputs into PDB/CIF/GJF companions based on the input topology/template.",
        ),
        click.option(
            "--ref-pdb",
            type=click.Path(path_type=Path, exists=True, dir_okay=False),
            default=None,
            help="Reference PDB/mmCIF topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
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
            # Per-command value, so the rendered default matches the command.
            show_default=(thresh_effective if thresh_effective is not None else False),
            help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
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
            help="Experimental, computationally expensive xTB solvent delta correction. Examples: water, methanol, acetonitrile, dmso, thf, toluene. 'none' disables it.",
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
                    show_default="the surface minimum",
                    help="Lower bound of color scale for plots (kcal/mol).",
                ),
                click.option(
                    "--zmax",
                    type=float,
                    default=None,
                    show_default="the surface maximum",
                    help="Upper bound of color scale for plots (kcal/mol).",
                ),
            ]
        )

    def decorator(func):
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator
