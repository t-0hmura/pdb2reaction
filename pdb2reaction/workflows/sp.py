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

from pdb2reaction.backends import (
    BackendError,
    apply_calc_file_to_calc_cfg,
    create_calculator,
    normalize_hessian_calc_mode,
)
from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    OUT_DIR_SP,
    UMA_CALC_KW,
    apply_backend_defaults,
)
from pdb2reaction.core.utils import (
    apply_yaml_overrides,
    cli_param_overridden,
    emit_dry_run_complete,
    format_elapsed,
    prepare_input_structure,
    resolve_charge_spin,
    resolve_freeze_atoms,
    validate_charge_spin_for_prepared,
    yaml_freeze_to_internal,
)
from pdb2reaction.cli.common_options import (
    add_ml_charge_spin_options,
    add_precision_option,
    add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from pdb2reaction.cli.decorators import (
    _write_error_json,
    load_merged_yaml_cfg,
    resolve_yaml_sources,
)
logger = logging.getLogger(__name__)


@click.command(
    name="sp",
)
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Input structure (PDB / mmCIF / XYZ / GJF).",
)
@click.option(
    "--workers", type=int, default=UMA_CALC_KW["workers"], show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
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
    help="Also compute a Hessian and save it to hessian.npy (active block when YAML geom.freeze_atoms is non-empty).",
)
@click.option(
    "--hessian-calc-mode", "hessian_calc_mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None, show_default="FiniteDifference",
    help=(
        "Hessian backend when --hess is set. FiniteDifference is the default "
        "for every backend; UMA, ORB, MACE, and AIMNet2 also support an "
        "explicit Analytical request."
    ),
)
@click.option(
    "--config", "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None, help="YAML config file with sections (calc:, geom:, …).",
)
@click.option(
    "--show-config/--no-show-config", "show_config",
    default=False, show_default=True, help="Print resolved configuration and continue execution.",
)
@click.option(
    "--dry-run/--no-dry-run", "dry_run",
    default=False, show_default=True, help="Validate options and print the plan without running the calculation.",
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
    help="Experimental, computationally expensive xTB solvent delta correction. Examples: water, methanol, acetonitrile, dmso, thf, toluene. 'none' disables it.",
)
@click.option(
    "--solvent-model", "solvent_model",
    default="alpb", type=click.Choice(["alpb", "cpcmx"]),
    show_default=True, help="xTB solvent model.",
)
@add_ml_charge_spin_options(allow_ref_pdb=False)
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
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
    calc_factory: Optional[str],
) -> None:
    """Compute a single-point MLIP energy + forces (and optionally Hessian)."""
    config_yaml, _override, _legacy = resolve_yaml_sources(
        config_yaml=config_yaml, override_yaml=None, args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, _override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml, override_yaml=None,
    )
    from pdb2reaction.core.utils import resolve_configured_charge_spin
    charge, spin = resolve_configured_charge_spin(
        merged_yaml_cfg, charge=charge, spin=spin, ligand_charge=ligand_charge,
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
        if sp_cfg["hess"]:
            try:
                calc_cfg["hessian_calc_mode"] = normalize_hessian_calc_mode(
                    sp_cfg.get("hessian_calc_mode")
                    or calc_cfg.get("hessian_calc_mode", "FiniteDifference")
                )
            except BackendError as exc:
                raise click.ClickException(str(exc)) from exc
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
        from pdb2reaction.backends import apply_backend_model_to_calc_cfg
        # Use the canonical helper like the other subcommands (was an inline
        # calc_cfg["model"]=... that never popped a --config YAML backend_model
        # token). Unconditional: the helper no-ops when neither the CLI arg nor
        # the YAML names a model.
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        # --calc-file overrides --backend with a user ASE Calculator (custom backend).
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        from pdb2reaction.backends import apply_effective_precision
        apply_effective_precision(
            calc_cfg,
            precision if cli_param_overridden(ctx, "precision") else None,
        )
        apply_backend_defaults(calc_cfg)

        # YAML uses human-facing 1-based atom indices; pysisyphus and the
        # calculator adapters use 0-based indices internally.
        if geom_cfg.get("freeze_atoms"):
            geom_cfg["freeze_atoms"] = yaml_freeze_to_internal(geom_cfg["freeze_atoms"])

        resolved_charge, resolved_spin = resolve_charge_spin(
            prepared_inputs, charge=charge, spin=spin,
            ligand_charge=ligand_charge, prefix="[sp]",
        )
        validate_charge_spin_for_prepared(prepared_inputs, resolved_charge, resolved_spin)
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        # Resolve this before --show-config/--dry-run and make the calculator
        # contract explicit. Geometry.set_calculator() also propagates its
        # freeze set at runtime; carrying it here keeps reported and runtime
        # configuration identical and avoids relying on that later side effect.
        resolve_freeze_atoms(
            geom_cfg, prepared_inputs[0].geom_path, freeze_links=False
        )
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
        if sp_cfg["hess"]:
            calc_cfg["return_partial_hessian"] = True

        out_dir_path = Path(sp_cfg["out_dir"]).resolve()

        if (
            sp_cfg["hess"]
            and str(calc_cfg.get("hessian_calc_mode", "")).lower() == "analytical"
            and calc_cfg.get("backend") == "uma"
            and int(calc_cfg.get("workers", 1)) > 1
        ):
            raise click.ClickException(
                "Analytical Hessian cannot be combined with UMA workers>1: "
                "the parallel predictor exposes no autograd model. Use "
                "workers=1 or select hessian_calc_mode='FiniteDifference'."
            )

        if show_config:
            click.echo(yaml.safe_dump(
                {"calc": calc_cfg, "geom": geom_cfg, "sp": sp_cfg},
                sort_keys=False, allow_unicode=True,
            ).rstrip())
            # Matches the rest of the family (opt/tsopt/freq/...): print and
            # continue. Use --dry-run if you want to stop after the print.

        if dry_run:
            click.echo(f"[sp] dry-run: would compute SP on {input_path} → {out_dir_path}")
            emit_dry_run_complete()
            return

        owned_outputs = tuple(
            out_dir_path / name
            for name in ("forces.npy", "hessian.npy", "result.json", "summary.json")
        )
        for protected in (
            input_path,
            prepared_inputs[0].original_path,
            prepared_inputs[0].source_path,
            prepared_inputs[0].geom_path,
            config_yaml,
        ):
            if protected is None:
                continue
            protected_path = Path(protected)
            if any(
                protected_path.resolve() == output.resolve()
                or (
                    protected_path.exists()
                    and output.exists()
                    and protected_path.samefile(output)
                )
                for output in owned_outputs
            ):
                raise click.UsageError(
                    f"Input {protected_path} collides with a reserved SP output "
                    f"path under {out_dir_path}."
                )

        out_dir_path.mkdir(parents=True, exist_ok=True)
        # These public files are owned by this invocation. Remove the previous
        # generation before calculator evaluation so a failed rerun cannot
        # leave old numerical results beside the current error envelope.
        for name in ("forces.npy", "hessian.npy", "result.json", "summary.json"):
            (out_dir_path / name).unlink(missing_ok=True)

        coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geom = geom_loader(prepared_inputs[0].geom_path, coord_type=coord_type, **coord_kwargs)

        # Default-verbosity entry summary (skipped in child mode).
        from pdb2reaction.core.utils import calculator_run_label, echo_run_summary
        echo_run_summary({
            "input": str(input_path),
            "backend": calculator_run_label(calc_cfg),
            "charge/spin": f"{calc_cfg.get('charge')}/{calc_cfg.get('spin')}",
            "out": str(out_dir_path),
        })
        hessian_mode: Optional[str] = None
        if sp_cfg["hess"]:
            hessian_mode = str(calc_cfg["hessian_calc_mode"])
            # Every built-in pysis calculator implements the same get_hessian
            # contract. Propagate the resolved request before constructing the
            # calculator; otherwise ORB/MACE/AIMNet2 silently retain their FD
            # default even when the CLI explicitly requested Analytical.
            calc_cfg["out_hess_torch"] = True

        calc = create_calculator(**calc_cfg)
        geom.set_calculator(calc)

        # Energy + forces
        t0 = time.perf_counter()
        ef_results = geom.get_energy_and_cart_forces_at(geom.cart_coords)
        geom.set_results(ef_results)
        energy_au = float(ef_results["energy"])
        forces_au = np.asarray(ef_results["forces"], dtype=float).reshape(-1, 3)
        if not np.isfinite(energy_au) or not np.all(np.isfinite(forces_au)):
            raise ValueError("Single-point energy and forces must be finite.")
        elapsed_ef = time.perf_counter() - t0
        click.echo(f"[sp] energy = {energy_au:.10f} a.u.  |force|_max = {np.max(np.abs(forces_au)):.4e} a.u./bohr  ({elapsed_ef:.2f} s)")

        forces_path = out_dir_path / "forces.npy"
        np.save(forces_path, forces_au)

        # Optional Hessian
        hessian_path: Optional[Path] = None
        if sp_cfg["hess"]:
            assert hessian_mode is not None
            click.echo(f"[sp] computing Hessian (mode={hessian_mode}) …")
            t0 = time.perf_counter()
            hess_results = calc.get_hessian(geom.atoms, geom.cart_coords)
            H_raw = hess_results["hessian"]
            if isinstance(H_raw, torch.Tensor):
                H_np = H_raw.detach().cpu().numpy()
                del H_raw
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            else:
                H_np = np.asarray(H_raw, dtype=float)
            n_full_dof = 3 * len(geom.atoms)
            n_active_dof = n_full_dof - 3 * len(calc_cfg.get("freeze_atoms", []))
            valid_shapes = {(n_full_dof, n_full_dof), (n_active_dof, n_active_dof)}
            if H_np.shape not in valid_shapes or not np.all(np.isfinite(H_np)):
                expected = " or ".join(str(shape) for shape in sorted(valid_shapes))
                raise ValueError(
                    "Single-point Hessian must be finite with full or active "
                    f"Cartesian shape {expected}; got {H_np.shape}."
                )
            try:
                geom.set_results(hess_results)
            except Exception:
                logger.debug("Failed to cache SP Hessian on Geometry", exc_info=True)
            elapsed_h = time.perf_counter() - t0
            hessian_path = out_dir_path / "hessian.npy"
            np.save(hessian_path, H_np)
            click.echo(f"[sp] Hessian {H_np.shape} written to {hessian_path}  ({elapsed_h:.2f} s)")

        # Summary
        elapsed_total = format_elapsed("[time] Elapsed Time for SP", time_start)
        _custom_calculator = (
            f"{Path(str(calc_cfg['calc_file'])).name}:"
            f"{calc_cfg.get('calc_factory', 'get_calculator')}"
            if calc_cfg["backend"] == "custom" and calc_cfg.get("calc_file")
            else None
        )
        from pdb2reaction.core.utils import calculator_provenance

        summary = {
            "stage": "sp",
            "status": "ok",
            "input": str(prepared_inputs[0].display_path),
            "backend": calc_cfg["backend"],
            "model": _custom_calculator or calc_cfg.get("model"),
            "custom_calculator": _custom_calculator,
            **calculator_provenance(calc_cfg),
            "charge": calc_cfg["charge"],
            "spin": calc_cfg["spin"],
            "n_atoms": len(geom.atoms),
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

        click.echo(format_elapsed("[time] Elapsed Time for SP", time_start))

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
        for prepared in prepared_inputs:
            prepared.cleanup()
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
