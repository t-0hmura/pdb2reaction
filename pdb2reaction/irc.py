# pdb2reaction/irc.py

"""
IRC calculations using the EulerPC predictor-corrector integrator with UMA.

Example:
    pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/

For detailed documentation, see: docs/irc.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, List

import os
import shutil
import sys
import textwrap

import click
from click.core import ParameterSource
import yaml
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC
from pdb2reaction.uma_pysis import uma_pysis
from pdb2reaction.defaults import CALC_KW_DEFAULT, GEOM_KW_DEFAULT, UMA_CALC_KW, IRC_KW
from pdb2reaction.utils import (
    load_yaml_dict,
    deep_update,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    resolve_freeze_atoms,
    prepared_cli_input,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
)


def _resolve_yaml_sources(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], bool]:
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use a single YAML source option."
        )
    if args_yaml_legacy is not None:
        return config_yaml, args_yaml_legacy, True
    return config_yaml, override_yaml, False


def _load_merged_yaml_cfg(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    deep_update(merged, load_yaml_dict(config_yaml))
    deep_update(merged, load_yaml_dict(override_yaml))
    return merged


def _echo_convert_trj_if_exists(
    trj_path: Path,
    prepared_input: "PreparedInputStructure",
    *,
    out_pdb: Optional[Path] = None,
) -> None:
    if trj_path.exists():
        ref_pdb = prepared_input.source_path if prepared_input.source_path.suffix.lower() == ".pdb" else None
        if convert_xyz_like_outputs(
            trj_path,
            prepared_input,
            ref_pdb_path=ref_pdb,
            out_pdb_path=out_pdb,
            context=f"'{trj_path.name}' outputs",
        ):
            targets = [p for p in (out_pdb,) if p is not None and p.exists()]
            if targets:
                written = ", ".join(f"'{p.name}'" for p in targets)
                click.echo(f"[convert] Wrote {written}.")


def _first_existing_artifact(out_dir: Path, patterns: Sequence[str]) -> Optional[Path]:
    """Resolve the first existing artifact for a list of relative patterns."""
    for pattern in patterns:
        if any(ch in pattern for ch in "*?[]"):
            for candidate in sorted(out_dir.glob(pattern)):
                if candidate.is_file():
                    return candidate.resolve()
            continue
        candidate = out_dir / pattern
        if candidate.is_file():
            return candidate.resolve()
    return None


def _link_or_copy_file(src: Path, dst: Path) -> bool:
    """Create a symlink when possible; fall back to copy."""
    try:
        if dst.exists() or dst.is_symlink():
            if dst.is_dir():
                return False
            dst.unlink()
        rel = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel)
        return True
    except Exception:
        try:
            shutil.copy2(src, dst)
            return True
        except Exception:
            return False


def _write_output_summary_md(out_dir: Path) -> None:
    """summary.md and key_* outputs are disabled."""
    return None

# --------------------------
# CLI
# --------------------------

@click.command(
    help="Run an IRC calculation with EulerPC. Only the documented CLI options are accepted; all other settings come from YAML.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, _trj.xyz, etc.).",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=False,
    help=(
        "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
        "(PDB inputs or XYZ/GJF with --ref-pdb)."
    ),
)
@click.option(
    "--workers",
    type=int,
    default=UMA_CALC_KW["workers"],
    show_default=True,
    help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=UMA_CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel UMA predictor (workers>1).",
)
@click.option(
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
        "when -q is omitted (requires PDB input or --ref-pdb)."
    ),
)
@click.option("-m", "--multiplicity", "spin", type=int, default=None, show_default=False, help="Spin multiplicity (2S+1) for the ML region.")
@click.option(
    "--max-cycles",
    type=int,
    default=None,
    help=(
        "Maximum number of IRC steps; used unless YAML sets irc.max_cycles. "
        "Defaults to 125 when not provided."
    ),
)
@click.option(
    "--step-size",
    type=float,
    default=None,
    help=(
        "Step length in mass-weighted coordinates; used unless YAML sets irc.step_length. "
        "Defaults to 0.10."
    ),
)
@click.option(
    "--root",
    type=int,
    default=None,
    help=(
        "Imaginary mode index used for the initial displacement; used unless YAML sets irc.root. "
        "Defaults to 0."
    ),
)
@click.option(
    "--forward/--no-forward",
    "forward",
    default=None,
    help=(
        "Run the forward IRC; used unless YAML sets irc.forward. "
        "Defaults to True."
    ),
)
@click.option(
    "--backward/--no-backward",
    "backward",
    default=None,
    help=(
        "Run the backward IRC; used unless YAML sets irc.backward. "
        "Defaults to True."
    ),
)
@click.option(
    "--freeze-links/--no-freeze-links",
    "freeze_links_flag",
    default=True,
    show_default=True,
    help="Freeze parent atoms of link hydrogens (PDB only).",
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option(
    "--out-dir",
    type=str,
    default="./result_irc/",
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="How UMA builds the Hessian (Analytical or FiniteDifference); used unless YAML sets calc.hessian_calc_mode. Defaults to 'FiniteDifference'.",
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
    help="Validate options and print the execution plan without running IRC.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    max_cycles: Optional[int],
    step_size: Optional[float],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    freeze_links_flag: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    hessian_calc_mode: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
) -> None:
    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception:
            return False

    config_yaml, override_yaml, used_legacy_yaml = _resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg = _load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )

    set_convert_file_enabled(convert_files)
    with prepared_cli_input(
        input_path,
        ref_pdb=ref_pdb,
        charge=charge,
        spin=spin,
        ligand_charge=ligand_charge,
        prefix="[irc]",
    ) as (prepared_input, resolved_charge, resolved_spin):
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path
        try:
            time_start = time.perf_counter()

            # --------------------------
            # 1) Assemble configuration: defaults < config < CLI(explicit) < override
            # --------------------------
            config_layer_cfg = load_yaml_dict(config_yaml)
            override_layer_cfg = load_yaml_dict(override_yaml)

            geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
            calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
            irc_cfg: Dict[str, Any] = dict(IRC_KW)

            apply_yaml_overrides(
                config_layer_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",),)),
                    (irc_cfg, (("irc",),)),
                ],
            )

            if _is_param_explicit("workers"):
                calc_cfg["workers"] = int(workers)
            if _is_param_explicit("workers_per_node"):
                calc_cfg["workers_per_node"] = int(workers_per_node)
            if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
                calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
            if _is_param_explicit("max_cycles") and max_cycles is not None:
                irc_cfg["max_cycles"] = int(max_cycles)
            if _is_param_explicit("step_size") and step_size is not None:
                irc_cfg["step_length"] = float(step_size)
            if _is_param_explicit("root") and root is not None:
                irc_cfg["root"] = int(root)
            if _is_param_explicit("forward") and forward is not None:
                irc_cfg["forward"] = bool(forward)
            if _is_param_explicit("backward") and backward is not None:
                irc_cfg["backward"] = bool(backward)
            if _is_param_explicit("out_dir"):
                irc_cfg["out_dir"] = str(out_dir)

            charge_value = calc_cfg.get("charge", resolved_charge)
            if charge_value is None:
                charge_value = resolved_charge
            calc_cfg["charge"] = int(charge_value)
            if _is_param_explicit("charge"):
                calc_cfg["charge"] = int(resolved_charge)

            spin_value = calc_cfg.get("spin", resolved_spin)
            if spin_value is None:
                spin_value = resolved_spin
            calc_cfg["spin"] = int(spin_value)
            if _is_param_explicit("spin"):
                calc_cfg["spin"] = int(resolved_spin)

            apply_yaml_overrides(
                override_layer_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",),)),
                    (irc_cfg, (("irc",),)),
                ],
            )

            # Normalize freeze_atoms and optionally add link-parent indices for PDB inputs
            resolve_freeze_atoms(geom_cfg, source_path, freeze_links_flag)

            # EulerPC currently only supports Cartesian coordinates
            geom_cfg["coord_type"] = "cart"

            # Ensure the calculator receives the freeze list used by geometry
            # (so FD Hessian can skip frozen DOF, etc.)
            calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
            calc_cfg["return_partial_hessian"] = False

            out_dir_path = Path(irc_cfg["out_dir"]).resolve()
            if show_config:
                click.echo(
                    pretty_block(
                        "yaml_layers",
                        {
                            "config": None if config_yaml is None else str(config_yaml),
                            "override_yaml": None if override_yaml is None else str(override_yaml),
                            "merged_keys": sorted(merged_yaml_cfg.keys()),
                        },
                    )
                )
            if dry_run:
                click.echo(
                    pretty_block(
                        "dry_run_plan",
                        {
                            "input_geometry": str(geom_input_path),
                            "output_dir": str(out_dir_path),
                            "freeze_links": bool(freeze_links_flag),
                            "convert_files": bool(convert_files),
                            "will_run_irc": True,
                            "will_write_trajectories": True,
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. IRC execution was skipped.")
                return

            out_dir_path.mkdir(parents=True, exist_ok=True)

            # Pretty-print configuration (expand freeze_atoms for readability)
            click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
            click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
            click.echo(pretty_block("irc",  {**irc_cfg, "out_dir": str(out_dir_path)}))

            # --------------------------
            # 2) Load geometry and configure UMA calculator
            # --------------------------
            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)

            geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

            calc = uma_pysis(**calc_cfg)
            geometry.set_calculator(calc)

            # --------------------------
            # 3) Construct and run EulerPC
            # --------------------------
            # EulerPC.__init__ forwards **kwargs directly to IRC.__init__
            eulerpc = EulerPC(geometry, **irc_cfg)

            click.echo("\n====== IRC (EulerPC) started ======\n")
            eulerpc.run()
            click.echo("\n====== IRC (EulerPC) finished ======\n")

            # --------------------------
            # 4) Convert trajectories to PDB based on input type
            # --------------------------
            suffix_prefix = irc_cfg.get("prefix", "")
            _echo_convert_trj_if_exists(
                out_dir_path / f"{suffix_prefix}{'finished_irc_trj.xyz'}",
                prepared_input,
                out_pdb=out_dir_path / f"{suffix_prefix}{'finished_irc.pdb'}" if prepared_input.source_path.suffix.lower() == ".pdb" else None,
            )
            _echo_convert_trj_if_exists(
                out_dir_path / f"{suffix_prefix}{'forward_irc_trj.xyz'}",
                prepared_input,
                out_pdb=out_dir_path / f"{suffix_prefix}{'forward_irc.pdb'}" if prepared_input.source_path.suffix.lower() == ".pdb" else None,
            )
            _echo_convert_trj_if_exists(
                out_dir_path / f"{suffix_prefix}{'backward_irc_trj.xyz'}",
                prepared_input,
                out_pdb=out_dir_path / f"{suffix_prefix}{'backward_irc.pdb'}" if prepared_input.source_path.suffix.lower() == ".pdb" else None,
            )

            # summary.md and key_* outputs are disabled.
            click.echo(format_elapsed("[time] Elapsed Time for IRC", time_start))

        except KeyboardInterrupt:
            click.echo("Interrupted by user.")
            sys.exit(130)
        except Exception as e:
            tb = textwrap.indent("".join(__import__("traceback").format_exception(type(e), e, e.__traceback__)), "  ")
            click.echo("Unhandled exception during IRC:\n" + tb)
            sys.exit(1)
