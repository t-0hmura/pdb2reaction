# pdb2reaction/irc.py

"""
irc — Concise CLI for IRC calculations with the EulerPC integrator
====================================================================

Usage (CLI)
-----
    # Minimal and full examples shown below.
    # Minimal (explicitly set charge/spin to avoid unintended conditions)
    pdb2reaction irc -i a.pdb -q 0 -s 1

    # Full example
    pdb2reaction irc \
      -i a.pdb -q 0 -s 1 \
      --max-cycles 20 \
      --step-size 0.5 \
      --root 0 \
      --forward True \
      --backward False \
      --freeze-links True \
      --out-dir "./result_irc/" \
      --hessian-calc-mode Analytical \
      --args-yaml args.yaml

Examples::
    # Forward-only with finite-difference Hessian and custom step size
    pdb2reaction irc -i ts.xyz -q -1 -s 2 --forward True --backward False \
      --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

    # Use a PDB input so trajectories are also exported as PDB
    pdb2reaction irc -i ts.pdb -q 0 -s 1 --max-cycles 50 --out-dir ./result_irc/

Description
-----
- Purpose: Run Intrinsic Reaction Coordinate (IRC) calculations using the EulerPC predictor–corrector integrator.
- Inputs: Any structure readable by `pysisyphus.helpers.geom_loader` (.pdb, .xyz, .trj, ...).
  If the input is `.pdb`, the generated trajectories are additionally converted to PDB.
- Configuration model: Only the CLI options listed above are accepted. All other parameters
  (geometry options, UMA calculator configuration, and detailed EulerPC/IRC settings) must be provided via YAML.
  Final configuration precedence: built-in defaults → YAML → CLI.
- Charge/spin defaults: `-q/--charge` and `-s/--spin` inherit values from `.gjf` templates when provided and otherwise fall back
  to `0`/`1`. Set them explicitly to avoid running under unintended conditions.
- CLI options:
  - `-i/--input PATH` (required): Structure file (.pdb/.xyz/.trj/…).
  - `-q/--charge INT`: Total charge; overrides `calc.charge` from YAML and defaults to a `.gjf` template value or `0` when omitted.
  - `-s/--spin INT` (default 1): Spin multiplicity (2S+1); overrides `calc.spin` and defaults to the template multiplicity or `1`.
  - `--max-cycles INT`: Max number of IRC steps; overrides `irc.max_cycles`.
  - `--step-size FLOAT`: Step length in mass-weighted coordinates; overrides `irc.step_length`.
  - `--root INT`: Imaginary mode index for the initial displacement; overrides `irc.root`.
  - `--forward BOOL`: Run the forward IRC (explicit `True`/`False`); overrides `irc.forward`.
  - `--backward BOOL`: Run the backward IRC (explicit `True`/`False`); overrides `irc.backward`.
  - `--freeze-links BOOL` (default `True`): Freeze parent atoms of link hydrogens when the input is PDB.
  - `--out-dir STR` (default `./result_irc/`): Output directory; overrides `irc.out_dir`.
  - `--hessian-calc-mode {Analytical,FiniteDifference}`: How UMA builds the Hessian; overrides `calc.hessian_calc_mode`.
  - `--args-yaml PATH`: YAML file with sections `geom`, `calc`, and `irc`.

Outputs (& Directory Layout)
-----
- All files are written under `--out-dir` (default: `./result_irc/`). The directory is created if missing.
- If `irc.prefix` is set, it is prepended to filenames.

    <out_dir>/
      ├─ irc_data.h5                    # HDF5 dump written every `irc.dump_every` steps
      ├─ <prefix>finished_irc.trj       # Full IRC trajectory (XYZ/TRJ)
      ├─ <prefix>forward_irc.trj        # Forward path segment
      ├─ <prefix>backward_irc.trj       # Backward path segment
      ├─ <prefix>finished_irc.pdb       # PDB conversion (only if input was .pdb)
      ├─ <prefix>forward_irc.pdb        # PDB conversion (only if input was .pdb)
      └─ <prefix>backward_irc.pdb       # PDB conversion (only if input was .pdb)

Notes:
-----
- CLI overrides YAML for:
  `charge`→`calc.charge`, `spin`→`calc.spin`, `step-size`→`irc.step_length`, `max-cycles`→`irc.max_cycles`,
  `root`→`irc.root`, `forward`→`irc.forward`, `backward`→`irc.backward`, `out-dir`→`irc.out_dir`,
  `hessian-calc-mode`→`calc.hessian_calc_mode`.
- UMA options are passed directly to `pdb2reaction.uma_pysis.uma_pysis`. With `device: "auto"`,
  the calculator selects GPU/CPU automatically. If `hessian_calc_mode: "FiniteDifference"`,
  `geom.freeze_atoms` can be used to skip frozen DOF in FD Hessian construction.
- `--step-size` is in mass-weighted coordinates; `--root` selects the imaginary-frequency index used
  for the initial displacement.
- Standard output includes progress and timing. Exit codes: `0` on success, `130` on `KeyboardInterrupt`,
  `1` on unhandled exceptions.
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
    convert_xyz_to_pdb,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    detect_freeze_links,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    charge_option,
    spin_option,
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


def _echo_convert_trj_to_pdb_if_exists(trj_path: Path, ref_pdb: Path, out_path: Path) -> None:
    if trj_path.exists():
        try:
            convert_xyz_to_pdb(trj_path, ref_pdb, out_path)
            click.echo(f"[convert] Wrote '{out_path}'.")
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert '{trj_path.name}' to PDB: {e}", err=True)


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
@charge_option(
    "Total charge; overrides calc.charge from YAML. Defaults to the .gjf template value when present, otherwise 0."
)
@spin_option(
    "Spin multiplicity (2S+1); overrides calc.spin from YAML. Defaults to the .gjf template value when present, otherwise 1."
)
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
@click.option("--out-dir", type=str, default="./result_irc/", show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None,
    help="How UMA builds the Hessian (Analytical or FiniteDifference); overrides calc.hessian_calc_mode from YAML.",
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
    spin: Optional[int],
    max_cycles: Optional[int],
    step_size: Optional[float],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    freeze_links_flag: bool,
    out_dir: str,
    hessian_calc_mode: Optional[str],
    args_yaml: Optional[Path],
) -> None:
    prepared_input = prepare_input_structure(input_path)
    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble configuration: defaults -> YAML overrides -> CLI overrides
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg:  Dict[str, Any] = dict(IRC_KW_DEFAULT)

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (irc_cfg, (("irc",),)),
            ],
        )

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)

        if hessian_calc_mode is not None:
            # pass through exactly as chosen; uma_pysis normalizes internally
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

        # Normalize any existing freeze list and optionally augment with link parents
        merged_freeze = merge_freeze_atom_indices(geom_cfg)
        if freeze_links_flag and input_path.suffix.lower() == ".pdb":
            try:
                detected = detect_freeze_links(input_path)
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

        # Ensure the calculator receives the freeze list used by geometry
        #      (so FD Hessian can skip frozen DOF, etc.)
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))

        out_dir_path = Path(irc_cfg["out_dir"]).resolve()
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Pretty-print configuration (expand freeze_atoms for readability)
        click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
        click.echo(pretty_block("calc", calc_cfg))
        click.echo(pretty_block("irc",  {**irc_cfg, "out_dir": str(out_dir_path)}))

        # --------------------------
        # 2) Load geometry and configure UMA calculator
        # --------------------------
        coord_type = geom_cfg.get("coord_type", "cart")
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)

        geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

        calc_builder_or_instance = uma_pysis(**calc_cfg)
        try:
            # Case 1: builder function returned
            geometry.set_calculator(calc_builder_or_instance())
        except TypeError:
            # Case 2: calculator instance returned directly
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
        # 4) Convert trajectories to PDB when the input was PDB
        # --------------------------
        if input_path.suffix.lower() == ".pdb":
            ref_pdb = input_path.resolve()

            # Whole IRC trajectory
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.trj'}",
                ref_pdb,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.pdb'}",
            )
            # Forward/backward trajectories
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.trj'}",
                ref_pdb,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.pdb'}",
            )
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.trj'}",
                ref_pdb,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.pdb'}",
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


# Script entry point
if __name__ == "__main__":
    cli()
