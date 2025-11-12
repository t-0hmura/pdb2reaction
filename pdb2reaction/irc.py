# pdb2reaction/irc.py
# -*- coding: utf-8 -*-
"""Command-line interface for IRC calculations using the EulerPC integrator.

Example
-------
::

  pdb2reaction irc \
    -i a.pdb -q 0 -s 1 \
    --max-cycles 20 \
    --step-size 0.5 \
    --root 0 \
    --forward True \
    --backward False \
    --out-dir "./result_irc/" \
    --args-yaml args.yaml

Overview
--------
* The CLI accepts only the options shown above. All other parameters (geometry,
  UMA configuration, and detailed EulerPC/IRC settings) must be provided via
  the YAML file.
* When the input is a ``.pdb`` file, the generated trajectories (``finished_irc.trj``),
  ``forward_irc.trj``, and ``backward_irc.trj`` are converted to PDB format.

Optional YAML layout
--------------------
::

  geom:
    coord_type: cart
    freeze_atoms: []

  calc:
    charge: 0
    spin: 1
    model: "uma-s-1p1"
    task_name: "omol"
    device: "auto"
    max_neigh: null
    radius: null
    r_edges: false
    out_hess_torch: false

  irc:
    # Base IRC settings
    downhill: false
    forward: true
    backward: true
    hessian_init: "calc"
    displ: "energy"
    displ_energy: 1.0e-3
    displ_length: 0.1
    rms_grad_thresh: 1.0e-3
    hard_rms_grad_thresh: null
    energy_thresh: 1.0e-6
    imag_below: 0.0
    force_inflection: true
    check_bonds: false
    prefix: ""
    dump_fn: "irc_data.h5"
    dump_every: 5

    # EulerPC-specific settings
    hessian_update: "bofill"
    hessian_recalc: null
    max_pred_steps: 500
    loose_cycles: 3
    corr_func: "mbs"

Notes
-----
* CLI arguments override values loaded from YAML (``charge``/``spin``,
  ``step_length`` via ``--step-size``, ``max_cycles`` via ``--max-cycles``,
  ``root``, ``forward``/``backward``, and ``out_dir``).
* UMA options are passed directly to :func:`pdb2reaction.uma_pysis.uma_pysis`.
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
from pdb2reaction.uma_pysis import uma_pysis
from pdb2reaction.utils import (
    convert_xyz_to_pdb,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
)


# --------------------------
# Default configuration
# --------------------------

GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",
    "freeze_atoms": [],  # 0-based indices
}

CALC_KW_DEFAULT: Dict[str, Any] = {
    "charge": 0,
    "spin": 1,
    "model": "uma-s-1p1",
    "task_name": "omol",
    "device": "auto",         # "cuda" | "cpu" | "auto"
    "max_neigh": None,
    "radius": None,           # Å
    "r_edges": False,
    "out_hess_torch": True,   # IRC supports Hessian input on the GPU
}

IRC_KW_DEFAULT: Dict[str, Any] = {
    # Arguments for IRC.__init__ (forwarded to EulerPC via **kwargs)
    "step_length": 0.10,         # Overridden by CLI --step-size
    "max_cycles": 125,           # Overridden by CLI --max-cycles
    "downhill": False,
    "forward": True,             # Overridden by CLI --forward
    "backward": True,            # Overridden by CLI --backward
    "root": 0,                   # Overridden by CLI --root
    "hessian_init": "calc",      # Default assumption: start from TS Hessian
    "displ": "energy",
    "displ_energy": 1.0e-3,
    "displ_length": 0.10,
    "rms_grad_thresh": 1.0e-3,
    "hard_rms_grad_thresh": None,
    "energy_thresh": 1.0e-6,
    "imag_below": 0.0,
    "force_inflection": True,
    "check_bonds": False,
    "out_dir": "./result_irc/",
    "prefix": "",
    "dump_fn": "irc_data.h5",
    "dump_every": 5,

    # EulerPC-specific options
    "hessian_update": "bofill",
    "hessian_recalc": None,
    "max_pred_steps": 500,
    "loose_cycles": 3,
    "corr_func": "mbs",
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
@click.option("-q", "--charge", type=int, required=True, help="Total charge; overrides calc.charge from YAML.")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1); overrides calc.spin from YAML.")
@click.option("--max-cycles", type=int, default=None, help="Maximum number of IRC steps; overrides irc.max_cycles from YAML.")
@click.option("--step-size", type=float, default=None, help="Step length in mass-weighted coordinates; overrides irc.step_length from YAML.")
@click.option("--root", type=int, default=None, help="Imaginary mode index used for the initial displacement; overrides irc.root from YAML.")
@click.option("--forward", type=bool, default=None, help="Run the forward IRC; overrides irc.forward from YAML. Specify True/False explicitly.")
@click.option("--backward", type=bool, default=None, help="Run the backward IRC; overrides irc.backward from YAML. Specify True/False explicitly.")
@click.option("--out-dir", type=str, default="./result_irc/", show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file providing extra parameters (sections: geom, calc, irc).",
)
def cli(
    input_path: Path,
    charge: int,
    spin: int,
    max_cycles: Optional[int],
    step_size: Optional[float],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
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

        geometry = geom_loader(input_path, coord_type=coord_type, **coord_kwargs)

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

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for IRC: {hh:02d}:{mm:02d}:{ss:06.3f}")

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = textwrap.indent("".join(__import__("traceback").format_exception(type(e), e, e.__traceback__)), "  ")
        click.echo("Unhandled exception during IRC:\n" + tb, err=True)
        sys.exit(1)


# Script entry point
if __name__ == "__main__":
    cli()
