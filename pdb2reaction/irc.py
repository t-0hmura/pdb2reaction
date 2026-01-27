# pdb2reaction/irc.py

"""
irc — IRC calculations with the EulerPC algorithm
====================================================================

Usage (CLI)
-----------
    pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q <charge>] [--ligand-charge <number|'RES:Q,...'>] \
        [--workers <int>] [--workers-per-node <int>] [-m <multiplicity>] \
        [--max-cycles <int>] [--step-size <float>] [--root <int>] \
        [--forward {True|False}] [--backward {True|False}] \
        [--freeze-links {True|False}] [--convert-files {True|False}] [--ref-pdb <file>] \
        [--out-dir <dir>] [--hessian-calc-mode {Analytical|FiniteDifference}] \
        [--args-yaml <file>]

Examples
--------
    # Forward-only with finite-difference Hessian and custom step size
    pdb2reaction irc -i ts.xyz -q -1 -m 2 --forward True --backward False \
        --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

    # Use a PDB input so trajectories are also exported as PDB
    pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/

Description
-----------
- Purpose: Run Intrinsic Reaction Coordinate (IRC) calculations using the EulerPC predictor–corrector integrator.
- Inputs: Any structure readable by `pysisyphus.helpers.geom_loader` (.pdb, .xyz, .trj, ...).
  If the input is `.pdb`, trajectory files written by the run are additionally converted to PDB (when conversion is enabled).
- For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling
  format-aware PDB/GJF output conversion.
- Configuration model: Only the CLI options listed above are accepted. All other parameters
  (geometry options, UMA calculator configuration, and detailed EulerPC/IRC settings) must be provided via YAML.
  Final configuration precedence: built-in defaults → CLI → YAML.
- Charge/spin defaults: `-q/--charge` and `-m/--multiplicity` inherit values from `.gjf` templates when provided. For non-`.gjf`
  inputs, omitting `-q/--charge` is allowed only when ``--ligand-charge`` is set: the full complex is treated as an
  enzyme–substrate system and its total charge is derived with ``extract.py``’s residue-aware logic. Otherwise the CLI aborts;
  multiplicity still defaults to 1 when unspecified, and an explicit `-q` overrides any derived charge.

CLI options
-----------
  - `-i/--input PATH` (required): Structure file (.pdb/.xyz/.trj/…).
  - `-q/--charge INT`: Total charge; sets `calc.charge`. Required for non-`.gjf` inputs; `.gjf` templates
    supply defaults when available.
  - `--workers`, `--workers-per-node`: UMA predictor parallelism (workers > 1 disables analytic Hessians).
  - `-m/--multiplicity INT` (default 1): Spin multiplicity (2S+1); sets `calc.spin` and defaults to the template multiplicity or `1`.
  - `--max-cycles INT`: Max number of IRC steps; sets `irc.max_cycles`.
  - `--step-size FLOAT`: Step length in mass-weighted coordinates; sets `irc.step_length`.
  - `--root INT`: Imaginary mode index for the initial displacement; sets `irc.root`.
  - `--forward BOOL`: Run the forward IRC (explicit `True`/`False`); sets `irc.forward`.
  - `--backward BOOL`: Run the backward IRC (explicit `True`/`False`); sets `irc.backward`.
  - `--freeze-links BOOL` (default `True`): Freeze parent atoms of link hydrogens when the input is PDB.
  - `--convert-files {True|False}` (default `True`): Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.
  - `--ref-pdb PATH`: Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).
  - `--out-dir STR` (default `./result_irc/`): Output directory; sets `irc.out_dir`.
  - `--hessian-calc-mode {Analytical,FiniteDifference}`: How UMA builds the Hessian; sets `calc.hessian_calc_mode`.
  - `--args-yaml PATH`: YAML file with sections `geom`, `calc`, and `irc`.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_irc/)
  ├─ irc_data.h5                    # Optional HDF5 dump (if enabled by EulerPC)
  ├─ <prefix>finished_irc.trj       # Full IRC trajectory (TRJ)
  ├─ <prefix>forward_irc.trj        # Forward-only segment (TRJ)
  ├─ <prefix>backward_irc.trj       # Backward-only segment (TRJ)
  ├─ <prefix>finished_irc.pdb       # PDB conversions (written when the input was .pdb and conversion is enabled)
  ├─ <prefix>forward_irc.pdb        # PDB conversions (written when the input was .pdb and conversion is enabled)
  └─ <prefix>backward_irc.pdb       # PDB conversions (written when the input was .pdb and conversion is enabled)

All files honor ``irc.prefix`` when it is set; the directory is created if missing.

Notes
-----
- YAML overrides CLI for overlapping keys such as `calc.charge`, `calc.spin`, IRC step control, and output paths.
- UMA options are passed directly to `pdb2reaction.uma_pysis.uma_pysis`. With `device: "auto"`,
  the calculator selects GPU/CPU automatically. If `hessian_calc_mode: "FiniteDifference"`,
  `geom.freeze_atoms` can be used to skip frozen DOF in FD Hessian construction.
  The geometry's freeze list is also forwarded to the calculator as `calc.freeze_atoms`.
- `--step-size` is in mass-weighted coordinates; `--root` selects the imaginary-frequency index used
  for the initial displacement.
- Output conversion steps can be disabled via `--convert-files False`.
- Standard output includes progress and timing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import sys
import textwrap

import click
import yaml
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC
from pdb2reaction.uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW
from pdb2reaction.utils import (
    convert_xyz_like_outputs,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    detect_freeze_links,
    merge_freeze_atom_indices,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
)


# --------------------------
# Default configuration
# --------------------------

CALC_KW_DEFAULT: Dict[str, Any] = dict(_UMA_CALC_KW)

IRC_KW_DEFAULT: Dict[str, Any] = {
    # Arguments for IRC.__init__ (forwarded to EulerPC via **kwargs)
    "step_length": 0.10,         # float, default step length in mass-weighted coordinates (overridden by CLI)
    "max_cycles": 125,           # int, maximum IRC steps (overridden by CLI)
    "downhill": False,           # bool, follow downhill potential (debug option)
    "forward": True,             # bool, integrate forward branch (CLI override --forward)
    "backward": True,            # bool, integrate backward branch (CLI override --backward)
    "root": 0,                   # int, imaginary mode index for initial displacement (CLI override --root)
    "hessian_init": "calc",      # str, initial Hessian source ("calc" = calculator-provided TS Hessian)
    "displ": "energy",           # str, displacement metric (energy|length)
    "displ_energy": 1.0e-3,      # float, energy step (Hartree) when displ == "energy"
    "displ_length": 0.10,        # float, length step in mass-weighted coordinates when displ == "length"
    "rms_grad_thresh": 1.0e-3,   # float, RMS gradient threshold for convergence (Hartree/bohr)
    "hard_rms_grad_thresh": None,# Optional[float], stricter RMS gradient cutoff
    "energy_thresh": 1.0e-6,     # float, energy-change threshold for convergence (Hartree)
    "imag_below": 0.0,           # float, treat imaginary frequency below this as zero
    "force_inflection": True,    # bool, stop when force inflection detected
    "check_bonds": False,        # bool, enable bond-change detection during IRC
    "out_dir": "./result_irc/",  # str, output directory
    "prefix": "",                # str, file name prefix

    # EulerPC-specific options
    "hessian_update": "bofill",  # str, Hessian update algorithm
    "hessian_recalc": None,      # Optional[int], force Hessian recalculation every N steps
    "max_pred_steps": 500,       # int, predictor steps per segment
    "loose_cycles": 3,           # int, cycles using looser thresholds
    "corr_func": "mbs",          # str, correction function selection
}


def _echo_convert_trj_if_exists(
    trj_path: Path,
    prepared_input: "PreparedInputStructure",
    *,
    out_pdb: Optional[Path] = None,
) -> None:
    if trj_path.exists():
        try:
            ref_pdb = prepared_input.source_path if prepared_input.source_path.suffix.lower() == ".pdb" else None
            convert_xyz_like_outputs(
                trj_path,
                prepared_input,
                ref_pdb_path=ref_pdb,
                out_pdb_path=out_pdb,
            )
            targets = [p for p in (out_pdb,) if p is not None and p.exists()]
            if targets:
                written = ", ".join(f"'{p.name}'" for p in targets)
                click.echo(f"[convert] Wrote {written}.")
        except Exception as e:
            click.echo(
                f"[convert] WARNING: Failed to convert '{trj_path.name}' outputs: {e}",
                err=True,
            )


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
    help="Input structure file (.pdb, .xyz, .trj, etc.).",
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
    default=CALC_KW_DEFAULT["workers"],
    show_default=True,
    help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=CALC_KW_DEFAULT["workers_per_node"],
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
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1) for the ML region.")
@click.option("--max-cycles", type=int, default=None, help="Maximum number of IRC steps; overrides irc.max_cycles from YAML.")
@click.option("--step-size", type=float, default=None, help="Step length in mass-weighted coordinates; overrides irc.step_length from YAML.")
@click.option("--root", type=int, default=None, help="Imaginary mode index used for the initial displacement; overrides irc.root from YAML.")
@click.option("--forward", type=bool, default=None, help="Run the forward IRC; overrides irc.forward from YAML. Specify True/False explicitly.")
@click.option("--backward", type=bool, default=None, help="Run the backward IRC; overrides irc.backward from YAML. Specify True/False explicitly.")
@click.option(
    "--freeze-links",
    "freeze_links_flag",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Freeze parent atoms of link hydrogens when the input is PDB.",
)
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
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
@click.option("--out-dir", type=str, default="./result_irc/", show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="How UMA builds the Hessian (Analytical or FiniteDifference); overrides calc.hessian_calc_mode from YAML. Defaults to 'FiniteDifference'.",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file providing extra parameters (sections: geom, calc, irc).",
)
def cli(
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
    args_yaml: Optional[Path],
) -> None:
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    apply_ref_pdb_override(prepared_input, ref_pdb)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input,
        charge,
        spin,
        ligand_charge=ligand_charge,
        prefix="[irc]",
    )
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble configuration: defaults -> CLI overrides -> YAML overrides
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg:  Dict[str, Any] = dict(IRC_KW_DEFAULT)

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)
        calc_cfg["workers"] = int(workers)
        calc_cfg["workers_per_node"] = int(workers_per_node)

        if hessian_calc_mode is not None:
            calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

        if max_cycles is not None:
            irc_cfg["max_cycles"] = int(max_cycles)
        if step_size is not None:
            irc_cfg["step_length"] = float(step_size)
        if root is not None:
            irc_cfg["root"] = int(root)
        if forward is not None:
            irc_cfg["forward"] = bool(forward)
        if backward is not None:
            irc_cfg["backward"] = bool(backward)
        if out_dir:
            irc_cfg["out_dir"] = str(out_dir)

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (irc_cfg, (("irc",),)),
            ],
        )

        # Normalize any existing freeze list and optionally augment with link parents
        merged_freeze = merge_freeze_atom_indices(geom_cfg)
        if freeze_links_flag and source_path.suffix.lower() == ".pdb":
            try:
                detected = detect_freeze_links(source_path)
            except Exception as e:
                click.echo(
                    f"[freeze-links] WARNING: Could not detect link parents: {e}",
                    err=True,
                )
                detected = []
            merged_freeze = merge_freeze_atom_indices(geom_cfg, detected)
            if merged_freeze:
                click.echo(
                    f"[freeze-links] Freeze atoms (0-based): {','.join(map(str, merged_freeze))}"
                )

        # EulerPC currently only supports Cartesian coordinates
        geom_cfg["coord_type"] = "cart"

        # Ensure the calculator receives the freeze list used by geometry
        # (so FD Hessian can skip frozen DOF, etc.)
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
        calc_cfg["return_partial_hessian"] = False

        out_dir_path = Path(irc_cfg["out_dir"]).resolve()
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Pretty-print configuration (expand freeze_atoms for readability)
        click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
        click.echo(pretty_block("calc", format_freeze_atoms_for_echo(calc_cfg)))
        click.echo(pretty_block("irc",  {**irc_cfg, "out_dir": str(out_dir_path)}))

        # --------------------------
        # 2) Load geometry and configure UMA calculator
        # --------------------------
        coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)

        geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

        calc_builder_or_instance = uma_pysis(**calc_cfg)
        if callable(calc_builder_or_instance):
            geometry.set_calculator(calc_builder_or_instance())
        else:
            geometry.set_calculator(calc_builder_or_instance)

        # --------------------------
        # 3) Construct and run EulerPC
        # --------------------------
        # EulerPC.__init__ forwards **kwargs directly to IRC.__init__
        eulerpc = EulerPC(geometry, **irc_cfg)

        click.echo("\n=== IRC (EulerPC) started ===\n")
        eulerpc.run()
        click.echo("\n=== IRC (EulerPC) finished ===\n")

        # --------------------------
        # 4) Convert trajectories to PDB based on input type
        # --------------------------
        suffix_prefix = irc_cfg.get("prefix", "")
        _echo_convert_trj_if_exists(
            out_dir_path / f"{suffix_prefix}{'finished_irc.trj'}",
            prepared_input,
            out_pdb=out_dir_path / f"{suffix_prefix}{'finished_irc.pdb'}" if prepared_input.source_path.suffix.lower() == ".pdb" else None,
        )
        _echo_convert_trj_if_exists(
            out_dir_path / f"{suffix_prefix}{'forward_irc.trj'}",
            prepared_input,
            out_pdb=out_dir_path / f"{suffix_prefix}{'forward_irc.pdb'}" if prepared_input.source_path.suffix.lower() == ".pdb" else None,
        )
        _echo_convert_trj_if_exists(
            out_dir_path / f"{suffix_prefix}{'backward_irc.trj'}",
            prepared_input,
            out_pdb=out_dir_path / f"{suffix_prefix}{'backward_irc.pdb'}" if prepared_input.source_path.suffix.lower() == ".pdb" else None,
        )

        click.echo(format_elapsed("[time] Elapsed Time for IRC", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = textwrap.indent("".join(__import__("traceback").format_exception(type(e), e, e.__traceback__)), "  ")
        click.echo("Unhandled exception during IRC:\n" + tb, err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()
