# pdb2reaction/commands/opt.py
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Literal

import click
import numpy as np
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from .uma_pysis import uma_pysis
from .utils import convert_xyz_to_pdb, freeze_links

@click.command(help="Optimization")
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Input structure (pdb/xyz)"
)
def cli(input_path: Path) -> None:
    print(f"opt {input_path.resolve()}")

# -----------------------------------------------
# Defaults variables (overridable via YAML/CLI)
# -----------------------------------------------
CALC_KW: Dict[str, Any] = dict(
    charge=0,
    spin=1,                # multiplicity
    model="uma-s-1p1",
    task_name="omol",
    device="auto",         # "cuda" or "cpu" or "auto"
    # out_hess_torch is set depending on optimizer
)

BASE_OPT_KW: Dict[str, Any] = dict(
    thresh="gau",
    max_cycles=10000,
    print_every=1,
    dump=True,
    dump_restart=False,
    align=False,
    out_dir="./dump/opt/",
)

LBFGS_KW: Dict[str, Any] = dict(
    keep_last=7,
    beta=1.0,
    max_step=0.20,
    double_damp=True,
    gamma_mult=False,
    line_search=False,
    mu_reg=None,
)

RFO_KW: Dict[str, Any] = dict(
    line_search=True,
    gdiis=True,
    gediis=False,
    gdiis_thresh=2.5e-3,
    gediis_thresh=1.0e-2,
    max_micro_cycles=25,
    adapt_step_func=False,
    hessian_update="bfgs",
    hessian_init="calc",
    hessian_recalc=None,
    hessian_recalc_adapt=None,
)

# -----------------------------------------------
# Defaults variables (overridable via YAML/CLI)
# -----------------------------------------------