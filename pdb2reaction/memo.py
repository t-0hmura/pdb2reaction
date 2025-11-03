# pdb2reaction/commands/opt.py
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Literal

import click
import numpy as np
import torch

from .uma_pysis import uma_pysis
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

# -----------------------
# Helpers
# -----------------------
def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    """Load a YAML file into a dictionary (empty dict if path is None)."""
    if path is None:
        return {}
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise click.ClickException(
            f"PyYAML is required. Please install with `pip install pyyaml`: {e}"
        )
    try:
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        if not isinstance(data, dict):
            raise TypeError("Root of YAML must be a mapping (dict).")
        return data
    except Exception as e:
        raise click.ClickException(f"Failed to read YAML: {path}\n  -> {e}")

def _parse_bool(v: Optional[object], default: bool = False) -> bool:
    """Parse a truthy/falsey value robustly (supports true/false/1/0/yes/no...)."""
    if v is None:
        return default
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("1", "true", "t", "yes", "y", "on"):
        return True
    if s in ("0", "false", "f", "no", "n", "off"):
        return False
    raise click.ClickException(f"Cannot parse boolean value: {v!r}")

def _numpy_int_array(x: Optional[Sequence[int]]) -> Optional[np.ndarray]:
    if x is None:
        return None
    arr = np.asarray(list(x), dtype=int)
    if arr.size == 0:
        return None
    return arr

# -----------------------
# Defaults (overridable via YAML/CLI)
# -----------------------
DEFAULT_CALC_KW: Dict[str, Any] = dict(
    charge=0,
    spin=1,                # multiplicity
    model="uma-s-1p1",
    task_name="omol",
    device="auto",         # "cuda" or "cpu" or "auto"
    # out_hess_torch is set depending on optimizer
)

DEFAULT_BASE_OPT_KW: Dict[str, Any] = dict(
    thresh="gau_loose",
    max_cycles=500,
    print_every=1,
    dump=True,
    dump_restart=False,
    align=False,
    out_dir="./dump/opt/",
)

DEFAULT_LBFGS_KW: Dict[str, Any] = dict(
    keep_last=7,
    beta=1.0,
    max_step=0.20,
    double_damp=True,
    gamma_mult=False,
    line_search=False,
    mu_reg=None,
)

DEFAULT_RFO_KW: Dict[str, Any] = dict(
    line_search=True,
    gdiis=True,
    gediis=False,          # keep both keys available if your RFOptimizer supports either
    gdiis_thresh=2.5e-3,   # rms(step)
    gediis_thresh=1.0e-2,  # rms(force)
    max_micro_cycles=25,
    adapt_step_func=False,
    hessian_update="bfgs",
    hessian_init="calc",
    hessian_recalc=None,
    hessian_recalc_adapt=None,
)

def _build_calculator(calc_kw: Dict[str, Any], optimizer_kind: Literal["lbfgs", "rfo"]) -> uma_pysis:
    """Build UMA calculator; keep Hessian on GPU (out_hess_torch) for RFO when supported."""
    out_hess_torch = (optimizer_kind == "rfo")
    kw = dict(calc_kw)
    # Accept both 'mult' and 'spin' for multiplicity; prefer 'spin'
    if "mult" in kw and "spin" not in kw:
        kw["spin"] = kw.pop("mult")
    # Some wrappers may not accept out_hess_torch; fallback cleanly
    try:
        return uma_pysis(**kw, out_hess_torch=out_hess_torch)
    except TypeError:
        return uma_pysis(**kw)

def _build_geometry(input_file: Path, calc: uma_pysis, freeze_atoms: Optional[Sequence[int]] = None):
    """Load geometry (cartesian) and attach the calculator; optionally freeze atoms."""
    geom = geom_loader(str(input_file), coord_type="cart")
    if freeze_atoms:
        geom.freeze_atoms = _numpy_int_array(freeze_atoms)
    geom.set_calculator(calc)
    return geom

def _run_lbfgs(geom, base_kw: Dict[str, Any], lbfgs_kw: Dict[str, Any]) -> None:
    kw = dict(base_kw)
    kw.update(lbfgs_kw)
    LBFGS(geometry=geom, **kw).run()

def _run_rfo(geom, base_kw: Dict[str, Any], rfo_kw: Dict[str, Any]) -> None:
    kw = dict(base_kw)
    kw.update(rfo_kw)
    RFOptimizer(geometry=geom, **kw).run()

def _finalize_out_dir(base_kw: Dict[str, Any], subdir: str) -> Path:
    """Ensure output directory exists; preserve user YAML value if provided."""
    out_dir = Path(base_kw.get("out_dir", f"./dump/{subdir}/"))
    _ensure_dir(out_dir)
    base_kw["out_dir"] = str(out_dir)
    return out_dir

def _choose_optimizer(opt_mode: str) -> Literal["lbfgs", "rfo"]:
    m = opt_mode.strip().lower()
    if m in ("light", "lbfgs"):
        return "lbfgs"
    if m in ("heavy", "rfo"):
        return "rfo"
    raise click.ClickException(
        f"--opt-mode must be 'light(lbfgs)' or 'heavy(rfo)', got: {opt_mode!r}"
    )

# -----------------------
# Click command
# -----------------------
@click.command(help="Single-structure optimization (LBFGS or RFO). All UMA/optimizer kwargs can be overridden via YAML.")
@click.option(
    "-i", "--input", "input_file",
    type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True,
    help="Input structure file (.pdb, .xyz, ...)",
)
@click.option("-q", "--charge", type=int, default=None, help="Total charge (priority: YAML > CLI > defaults)")
@click.option("-s", "--mult", "mult", type=int, default=None, help="Multiplicity (spin), same priority as above")
# Two forms so both '--freeze-links' and '--freeze_links true' work
@click.option("--freeze-links/--no-freeze-links", default=None, help="Freeze link atoms. Prefer specifying indices via YAML (freeze.freeze_atoms).")
@click.option("--freeze_links", type=click.BOOL, default=None, help="Same as --freeze-links (true/false form).")
@click.option("--max-cycles", "max_cycles", type=int, default=None, help="Override optimizer.base.max_cycles")
@click.option("--max_cycles", "max_cycles_us", type=int, default=None, help="Alias of --max-cycles")
@click.option("--opt-mode", "opt_mode", type=str, default="heavy", show_default=True, help="light (=LBFGS) or heavy (=RFO)")
@click.option("--opt_mode", "opt_mode_us", type=str, default=None, help="Alias of --opt-mode")
@click.option(
    "--args-yaml", "args_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None,
    help="YAML with extra args: uma.*, optimizer.base/lbfgs/rfo.*, freeze.freeze_atoms, ...",
)
@click.option(
    "--args_yaml", "args_yaml_us",
    type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None,
    help="Alias of --args-yaml",
)
def cli(
    input_file: Path,
    charge: Optional[int],
    mult: Optional[int],
    freeze_links: Optional[bool],
    freeze_links_bool: Optional[bool],
    max_cycles: Optional[int],
    max_cycles_us: Optional[int],
    opt_mode: Optional[str],
    opt_mode_us: Optional[str],
    args_yaml: Optional[Path],
    args_yaml_us: Optional[Path],
) -> None:
    # Resolve aliases
    args_yaml = args_yaml or args_yaml_us
    max_cycles = max_cycles if max_cycles is not None else max_cycles_us
    if opt_mode_us:
        opt_mode = opt_mode_us

    # Load YAML (if provided)
    y = _load_yaml(args_yaml)
    y_uma = dict(y.get("uma", {}))
    y_opt = dict(y.get("optimizer", {}))
    y_base = dict(y_opt.get("base", {}))
    y_lbfgs = dict(y_opt.get("lbfgs", {}))
    y_rfo = dict(y_opt.get("rfo", {}))
    y_freeze = dict(y.get("freeze", {}))

    # Final booleans
    if freeze_links is None:
        freeze_links = freeze_links_bool  # may still be None
    freeze_links = _parse_bool(freeze_links, default=bool(y_freeze.get("freeze_links", False)))

    # Calculator kwargs (priority: YAML > CLI > defaults)
    calc_kw = dict(DEFAULT_CALC_KW)
    calc_kw.update(y_uma)
    if charge is not None:
        calc_kw["charge"] = charge
    if mult is not None:
        calc_kw["spin"] = mult

    # Base optimizer kwargs
    base_kw = dict(DEFAULT_BASE_OPT_KW)
    base_kw.update(y_base)
    if max_cycles is not None:
        base_kw["max_cycles"] = int(max_cycles)

    # Choose optimizer
    optimizer_kind: Literal["lbfgs", "rfo"] = _choose_optimizer(opt_mode or "heavy")

    # Method-specific kwargs
    if optimizer_kind == "lbfgs":
        method_kw = dict(DEFAULT_LBFGS_KW)
        method_kw.update(y_lbfgs)
    else:
        method_kw = dict(DEFAULT_RFO_KW)
        method_kw.update(y_rfo)

    # Output directory (preserve YAML if provided)
    out_dir = _finalize_out_dir(base_kw, subdir="opt")

    # Build calculator & geometry
    torch.set_default_dtype(torch.float64)
    calc = _build_calculator(calc_kw, optimizer_kind)
    freeze_atoms: Optional[Sequence[int]] = y_freeze.get("freeze_atoms")
    if freeze_links and not freeze_atoms:
        click.echo("[WARN] --freeze-links=True but freeze.freeze_atoms is not set; no atoms will be frozen.", err=True)
    geom = _build_geometry(input_file, calc, freeze_atoms)

    # Banner
    device = "cuda" if torch.cuda.is_available() else "cpu"
    try:
        n_atoms = len(geom.atoms)
    except Exception:
        n_atoms = getattr(geom, "n_atoms", -1)
    click.echo(f"Device: {device}")
    click.echo(f"Atoms:  {n_atoms}")
    try:
        click.echo(f"Start energy (au): {geom.energy:.8f}")
    except Exception:
        pass

    # Run optimization
    if optimizer_kind == "lbfgs":
        click.echo("Running LBFGS ...")
        _run_lbfgs(geom, base_kw, method_kw)
    else:
        click.echo("Running RFOptimizer (RFO) ...")
        _run_rfo(geom, base_kw, method_kw)

    click.echo(f"\nDone. Outputs are under: {out_dir.resolve()}")

以上のスクリプトを参考にして、

pdb2reaction opt -i a.pdb -charge 0 -mult 1 --freeze_links true --max_cycles 10000 --opt_mode light(lbfgs) or heavy(rfo) --args_yaml additional_args.yaml
このコマンドのように動くようにしてください。その他の、umaやoptimizerの変数はadditional_args.yamlで入力できるようにしてください。