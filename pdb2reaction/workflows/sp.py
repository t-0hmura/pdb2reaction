"""
Single-point MLIP energy / forces (and optional Hessian) calculation.

Example:
    pdb2reaction sp -i input.pdb -q 0 -m 1
    pdb2reaction sp -i input.pdb -q 0 -m 1 --hess

For detailed documentation, see: docs/sp.md
"""
# DOMAIN_PURE

from __future__ import annotations

import gc
import logging
import sys
import textwrap
import time
import traceback
from pathlib import Path
from typing import Optional

import click
import numpy as np
import torch
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.constants import ANG2BOHR, AU2EV

from pdb2reaction.backends import create_calculator, apply_calc_file_to_calc_cfg
from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    OUT_DIR_SP,
    UMA_CALC_KW,
    apply_backend_defaults,
)
from pdb2reaction.core.utils import (
    apply_yaml_overrides,
    cli_param_overridden,
    format_elapsed,
    prepare_input_structure,
    resolve_charge_spin,
    resolve_freeze_atoms,
    set_convert_file_enabled,
    validate_charge_spin_for_prepared,
)
from pdb2reaction.cli.common_options import (
    add_ml_charge_spin_options,
    add_precision_option,
    add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_print_every_option,
)
from pdb2reaction.cli.decorators import (
    _write_error_json,
    load_merged_yaml_cfg,
    resolve_yaml_sources,
)
from pdb2reaction.workflows.freq import _calc_full_hessian_torch, _torch_device

logger = logging.getLogger(__name__)

EV2AU = 1.0 / AU2EV
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)


@click.command(
    name="sp",
)
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Input structure (PDB / XYZ / GJF).",
)
@click.option(
    "--workers", type=int, default=UMA_CALC_KW["workers"], show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: analytical Hessian raises a RuntimeError when workers>1; pass --hessian-calc-mode FiniteDifference explicitly.",
)
@click.option(
    "--workers-per-node", "workers_per_node",
    type=int, default=UMA_CALC_KW["workers_per_node"], show_default=True,
    help="Workers per node when using a parallel MLIP predictor (workers>1).",
)
@click.option(
    "-o", "--out-dir", type=str, default=OUT_DIR_SP,
    show_default=True, help="Output directory.",
)
@click.option(
    "--hess/--no-hess", "do_hess", default=False, show_default=True,
    help="Also compute the full Hessian and save to hessian.npy.",
)
@click.option(
    "--hessian-calc-mode", "hessian_calc_mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None, show_default=False,
    help="Hessian backend when --hess is set. Analytical only works for UMA; other backends fall back to FiniteDifference.",
)
@click.option(
    "--convert-files/--no-convert-files", "convert_files",
    default=True, show_default=True,
    help="Auto-convert output XYZ-like files into companion PDB files written alongside them when the input had PDB metadata.",
)
@click.option(
    "--config", "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None, help="YAML config file with sections (calc:, geom:, …).",
)
@click.option(
    "--show-config/--no-show-config", "show_config",
    default=False, help="Print resolved configuration and continue execution.",
)
@click.option(
    "--dry-run/--no-dry-run", "dry_run",
    default=False, help="Validate options and print the plan without running the calculation.",
)
@click.option(
    "--out-json/--no-out-json", "out_json",
    default=False, show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "-b", "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"]),
    default="uma", show_default=True, help="MLIP backend.",
)
@click.option(
    "--solvent", default="none", show_default=True,
    help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.",
)
@click.option(
    "--solvent-model", "solvent_model",
    default="alpb", type=click.Choice(["alpb", "cpcmx"]),
    show_default=True, help="xTB solvent model.",
)
@add_ml_charge_spin_options()
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@add_print_every_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    workers: int,
    workers_per_node: int,
    out_dir: str,
    do_hess: bool,
    hessian_calc_mode: Optional[str],
    convert_files: bool,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    backend: str,
    solvent: str,
    solvent_model: str,
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: str,
    print_every: Optional[int],
) -> None:
    """Compute a single-point MLIP energy + forces (and optionally Hessian)."""
    set_convert_file_enabled(convert_files)

    config_yaml, _override, _legacy = resolve_yaml_sources(
        config_yaml=config_yaml, override_yaml=None, args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, _override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml, override_yaml=None,
    )

    prepared_inputs = [prepare_input_structure(input_path)]
    out_dir_path: Optional[Path] = None
    time_start: Optional[float] = time.perf_counter()

    try:
        geom_cfg: dict = dict(GEOM_KW_DEFAULT)
        calc_cfg: dict = dict(UMA_CALC_KW)
        sp_cfg: dict = {"out_dir": out_dir, "hess": False, "hessian_calc_mode": None}

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (sp_cfg, (("sp",),)),
            ],
        )

        # CLI overrides
        if cli_param_overridden(ctx, "out_dir"):
            sp_cfg["out_dir"] = out_dir
        if cli_param_overridden(ctx, "do_hess"):
            sp_cfg["hess"] = bool(do_hess)
        if cli_param_overridden(ctx, "hessian_calc_mode") and hessian_calc_mode is not None:
            sp_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
        if cli_param_overridden(ctx, "workers"):
            calc_cfg["workers"] = int(workers)
        if cli_param_overridden(ctx, "workers_per_node"):
            calc_cfg["workers_per_node"] = int(workers_per_node)
        if cli_param_overridden(ctx, "backend"):
            calc_cfg["backend"] = str(backend).lower()
        if cli_param_overridden(ctx, "solvent"):
            calc_cfg["solvent"] = str(solvent)
        if cli_param_overridden(ctx, "solvent_model"):
            calc_cfg["solvent_model"] = str(solvent_model)
        if cli_param_overridden(ctx, "precision") and precision is not None:
            calc_cfg["precision"] = str(precision)
        if cli_param_overridden(ctx, "backend_model") and backend_model is not None:
            calc_cfg["model"] = str(backend_model)
        # --calc-file overrides --backend with a user ASE Calculator (custom backend).
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        if cli_param_overridden(ctx, "print_every") and print_every is not None:
            calc_cfg["print_every"] = int(print_every)

        apply_backend_defaults(calc_cfg)

        resolved_charge, resolved_spin = resolve_charge_spin(
            prepared_inputs, charge=charge, spin=spin,
            ligand_charge=ligand_charge, prefix="[sp]",
        )
        validate_charge_spin_for_prepared(prepared_inputs, resolved_charge, resolved_spin)
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        out_dir_path = Path(sp_cfg["out_dir"]).resolve()

        if show_config:
            click.echo(yaml.safe_dump(
                {"calc": calc_cfg, "geom": geom_cfg, "sp": sp_cfg},
                sort_keys=False, allow_unicode=True,
            ).rstrip())
            # Matches the rest of the family (opt/tsopt/freq/...): print and
            # continue. Use --dry-run if you want to stop after the print.

        if dry_run:
            click.echo(f"[sp] dry-run: would compute SP on {input_path} → {out_dir_path}")
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)
        resolve_freeze_atoms(geom_cfg, prepared_inputs[0].geom_path, freeze_links=False)

        coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geom = geom_loader(prepared_inputs[0].geom_path, coord_type=coord_type, **coord_kwargs)

        # Default-verbosity entry summary (skipped in child mode).
        from pdb2reaction.core.utils import echo_run_summary
        _model = calc_cfg.get("model")
        echo_run_summary({
            "input": str(input_path),
            "backend": (
                f"{calc_cfg.get('backend', '?')} ({_model}, {calc_cfg.get('precision', 'fp32')})"
                if _model else calc_cfg.get("backend", "?")
            ),
            "charge/spin": f"{calc_cfg.get('charge')}/{calc_cfg.get('spin')}",
            "out": str(out_dir_path),
        })
        calc = create_calculator(**calc_cfg)
        geom.set_calculator(calc)

        # Energy + forces
        t0 = time.perf_counter()
        energy_au = float(geom.energy)
        gradient_au = np.asarray(geom.gradient, dtype=float)  # Hartree / Bohr
        forces_au = -gradient_au.reshape(-1, 3)
        elapsed_ef = time.perf_counter() - t0
        click.echo(f"[sp] energy = {energy_au:.10f} a.u.  |force|_max = {np.max(np.abs(forces_au)):.4e} a.u./bohr  ({elapsed_ef:.2f} s)")

        forces_path = out_dir_path / "forces.npy"
        np.save(forces_path, forces_au)

        # Optional Hessian
        hessian_path: Optional[Path] = None
        if sp_cfg["hess"]:
            mode = sp_cfg.get("hessian_calc_mode") or ("Analytical" if calc_cfg["backend"] == "uma" else "FiniteDifference")
            click.echo(f"[sp] computing full Hessian (mode={mode}) …")
            t0 = time.perf_counter()
            if mode.lower() == "analytical" and calc_cfg["backend"] == "uma":
                device = _torch_device()
                H_tensor = _calc_full_hessian_torch(geom, calc_cfg, device)
                H_np = H_tensor.detach().cpu().numpy()
                # C-NEED: sp dumps Hessian to npy; release GPU copy now.
                del H_tensor
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            else:
                # FiniteDifference via pysisyphus geom.hessian path
                H_np = np.asarray(geom.hessian, dtype=float)
            elapsed_h = time.perf_counter() - t0
            hessian_path = out_dir_path / "hessian.npy"
            np.save(hessian_path, H_np)
            click.echo(f"[sp] Hessian {H_np.shape} written to {hessian_path}  ({elapsed_h:.2f} s)")

        # Summary
        elapsed_total = format_elapsed("[sp] Elapsed Time", time_start)
        summary = {
            "stage": "sp",
            "status": "ok",
            "input": str(prepared_inputs[0].source_path),
            "backend": calc_cfg["backend"],
            "charge": calc_cfg["charge"],
            "spin": calc_cfg["spin"],
            "energy_au": energy_au,
            "forces_path": str(forces_path),
            "hessian_path": str(hessian_path) if hessian_path else None,
            "elapsed": elapsed_total,
        }
        if out_json:
            from pdb2reaction.core.utils import write_result_json
            write_result_json(
                out_dir_path, summary,
                command="sp",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        click.echo(f"[sp] {elapsed_total}.")

    except click.ClickException:
        # Click usage errors (missing flag, bad value) carry their own clean
        # render. Let Click handle the exit so the user sees a one-line
        # "Error: ..." instead of a full Python traceback.
        raise
    except Exception as exc:
        _write_error_json(
            out_dir_path or Path(sp_cfg.get("out_dir", OUT_DIR_SP)).resolve(),
            "sp", exc, "UnhandledError", time_start,
        )
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        click.echo(
            "Unhandled error during single-point calculation:\n"
            + textwrap.indent(tb, "  "),
            err=True,
        )
        sys.exit(1)
    finally:
        try:
            del geom  # type: ignore[name-defined]
        except NameError:
            pass
        try:
            del calc  # type: ignore[name-defined]
        except NameError:
            pass
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


