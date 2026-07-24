"""
IRC calculations using the EulerPC predictor-corrector integrator.

Example:
    pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/

For detailed documentation, see: docs/irc.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import logging
import sys

import gc

import click
from pdb2reaction.core.output import emit
import numpy as np
import time
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import GEOM_KW_DEFAULT, UMA_CALC_KW, IRC_KW, OUT_DIR_IRC, apply_backend_defaults
from pdb2reaction.core.utils import (
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    resolve_freeze_atoms,
    prepared_cli_input,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    cli_param_overridden,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
)
from pdb2reaction.cli.common_options import (
    add_ml_charge_spin_options,
    add_precision_option, add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_irc_pos_def_option,
)
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, render_cli_exception

logger = logging.getLogger(__name__)


def _directional_endpoint_energy_fields(all_energies: Any, ts_energy: Any) -> Dict[str, Any]:
    """Build standalone IRC energy fields without assigning chemical R/P identity."""
    first = float(all_energies[0]) if len(all_energies) > 0 else None
    last = float(all_energies[-1]) if len(all_energies) > 0 else None
    ts = float(ts_energy) if ts_energy is not None else None
    return {
        "energy_first_hartree": first,
        "energy_ts_hartree": ts,
        "energy_last_hartree": last,
        "endpoint_energy_orientation": "finished_first_to_finished_last",
        # Compatibility aliases: directional only, not chemical assignments.
        "energy_reactant_hartree": first,
        "energy_product_hartree": last,
    }


def _consume_mw_hessian_to_cartesian_active(
    mw_hessian: Any,
    mass_sqrt_active: Any,
) -> np.ndarray:
    """Mass-unweight a terminal active Hessian in place and return it on CPU.

    The input is consumed: callers must first detach it from the EulerPC owner
    and must not reuse it. This avoids a second dense device allocation at the
    IRC-to-cache stage boundary.
    """
    if isinstance(mw_hessian, torch.Tensor):
        ms_t = torch.as_tensor(
            mass_sqrt_active,
            dtype=mw_hessian.dtype,
            device=mw_hessian.device,
        )
        with torch.no_grad():
            mw_hessian.mul_(ms_t.unsqueeze(1))
            mw_hessian.mul_(ms_t.unsqueeze(0))
        return mw_hessian.detach().cpu().numpy()

    result = np.asarray(mw_hessian)
    masses = np.asarray(mass_sqrt_active, dtype=result.dtype)
    result *= masses[:, None]
    result *= masses[None, :]
    return result


def _irc_output_path(eulerpc: EulerPC, filename: str) -> Path:
    """Resolve an engine-authored IRC filename, including normalized prefix."""
    return Path(eulerpc.get_path_for_fn(filename))


_IRC_GENERATION_FILENAMES = tuple(
    f"{stem}{suffix}"
    for stem in (
        "finished_irc",
        "forward_irc",
        "backward_irc",
        "finished_first",
        "finished_last",
        "forward_first",
        "forward_last",
        "backward_first",
        "backward_last",
    )
    for suffix in (
        ("_trj.xyz", ".pdb", ".cif")
        if stem.endswith("_irc")
        else (".xyz", ".pdb", ".cif")
    )
)


def _prepare_irc_output_dir(
    path: Path,
    *,
    prefix: str = "",
    protected_inputs: tuple[Optional[Path], ...] = (),
) -> Path:
    """Invalidate command-owned IRC artifacts before a real generation."""
    resolved = Path(path).resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    normalized_prefix = f"{prefix}_" if prefix else ""
    owned = [
        *(resolved / f"{normalized_prefix}{name}" for name in _IRC_GENERATION_FILENAMES),
        resolved / "result.json",
        resolved / "summary.json",
    ]
    reserved = {candidate.resolve() for candidate in owned}
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise click.UsageError(
                f"Input {protected} collides with a reserved IRC output path "
                f"under {resolved}."
            )
    for candidate in owned:
        candidate.unlink(missing_ok=True)
    return resolved


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





@click.command(
    help="Run an IRC calculation with EulerPC. Only the documented CLI options are accepted; all other settings come from YAML.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .cif, .mmcif, .xyz, .gjf, _trj.xyz, etc.).",
)
@click.option(
    "--workers",
    type=int,
    default=UMA_CALC_KW["workers"],
    show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=UMA_CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel MLIP predictor (workers>1).",
)
@click.option(
    "--max-cycles",
    type=int,
    default=None,
    help=(
        "Maximum number of IRC steps; an explicit value overrides YAML irc.max_cycles. "
        "Defaults to 125 when not provided."
    ),
)
@click.option(
    "--step-size",
    type=float,
    default=None,
    help=(
        "Step length in Bohr (unweighted Cartesian coordinates); an explicit value overrides YAML irc.step_length. "
        "Default: 0.10 Bohr."
    ),
)
@click.option(
    "--never-stop/--no-never-stop",
    "never_stop",
    default=None,
    help=(
        "Ignore transient energy increases/plateaus and keep tracing through "
        "small shoulders. Gradient/integrator convergence and --max-cycles "
        "still stop the run. An explicit toggle overrides YAML irc.never_stop; default off."
    ),
)
@click.option(
    "--root",
    type=int,
    default=None,
    help=(
        "Imaginary mode index used for the initial displacement; an explicit value overrides YAML irc.root. "
        "Defaults to 0."
    ),
)
@click.option(
    "--forward/--no-forward",
    "forward",
    default=None,
    help=(
        "Run the forward IRC; an explicit toggle overrides YAML irc.forward. "
        "Defaults to True."
    ),
)
@click.option(
    "--backward/--no-backward",
    "backward",
    default=None,
    help=(
        "Run the backward IRC; an explicit toggle overrides YAML irc.backward. "
        "Defaults to True."
    ),
)
@click.option(
    "--freeze-links/--no-freeze-links",
    "freeze_links_flag",
    default=True,
    show_default=True,
    help="Freeze parent atoms of cap hydrogens (PDB/mmCIF input or XYZ/GJF with --ref-pdb).",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--tr-projection",
    type=click.Choice(["constrained", "legacy-active"], case_sensitive=False),
    default=None,
    help=(
        "Rigid-mode treatment for a frozen/partial Hessian. 'constrained' "
        "removes only full-system rigid motions compatible with the anchors "
        "(default); 'legacy-active' is deprecated comparison-only behavior and "
        "must not be used for pass/HOSP transition-state certification."
    ),
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/CIF/GJF companions based on the input format.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB/mmCIF topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option(
    "-o", "--out-dir",
    type=str,
    default=OUT_DIR_IRC,
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
    default=None,
    help="How the ML backend builds the Hessian (Analytical or FiniteDifference); an explicit value overrides YAML calc.hessian_calc_mode. Defaults to 'FiniteDifference'.",
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
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              show_default=True, help="xTB solvent model.")
@add_ml_charge_spin_options()
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@add_irc_pos_def_option()
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
    never_stop: Optional[bool],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    freeze_links_flag: bool,
    freeze_atoms_text: Optional[str],
    tr_projection: Optional[str],
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
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
    irc_pos_def: Optional[bool],
) -> None:
    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )
    from pdb2reaction.core.utils import resolve_configured_charge_spin
    charge, spin = resolve_configured_charge_spin(
        merged_yaml_cfg, charge=charge, spin=spin, ligand_charge=ligand_charge,
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
        validate_charge_spin_for_prepared(prepared_input, resolved_charge, resolved_spin)
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path
        calc = eulerpc = geometry = None
        time_start = time.perf_counter()
        out_dir_path = Path(out_dir).resolve()
        try:
            geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
            calc_cfg: Dict[str, Any] = dict(UMA_CALC_KW)
            irc_cfg: Dict[str, Any] = dict(IRC_KW)

            apply_yaml_overrides(
                config_layer_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",),)),
                    (irc_cfg, (("irc",),)),
                ],
            )

            if cli_param_overridden(ctx, "workers"):
                calc_cfg["workers"] = int(workers)
            if cli_param_overridden(ctx, "workers_per_node"):
                calc_cfg["workers_per_node"] = int(workers_per_node)
            if cli_param_overridden(ctx, "backend"):
                calc_cfg["backend"] = backend
            if cli_param_overridden(ctx, "solvent"):
                calc_cfg["solvent"] = solvent
            if cli_param_overridden(ctx, "solvent_model"):
                calc_cfg["solvent_model"] = solvent_model
            from pdb2reaction.backends import apply_effective_precision
            apply_effective_precision(calc_cfg, precision)
            from pdb2reaction.backends import apply_backend_model_to_calc_cfg
            # Unconditional: also pops a raw backend_model token from a --config YAML
            # (the helper no-ops when neither the CLI arg nor the YAML names one).
            apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
            from pdb2reaction.backends import apply_calc_file_to_calc_cfg
            apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
            apply_backend_defaults(calc_cfg)
            if cli_param_overridden(ctx, "hessian_calc_mode") and hessian_calc_mode is not None:
                calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
            if cli_param_overridden(ctx, "max_cycles") and max_cycles is not None:
                irc_cfg["max_cycles"] = int(max_cycles)
            if cli_param_overridden(ctx, "step_size") and step_size is not None:
                irc_cfg["step_length"] = float(step_size)
            if cli_param_overridden(ctx, "never_stop") and never_stop is not None:
                irc_cfg["never_stop"] = bool(never_stop)
            if cli_param_overridden(ctx, "root") and root is not None:
                irc_cfg["root"] = int(root)
            if cli_param_overridden(ctx, "forward") and forward is not None:
                irc_cfg["forward"] = bool(forward)
            if cli_param_overridden(ctx, "backward") and backward is not None:
                irc_cfg["backward"] = bool(backward)
            if cli_param_overridden(ctx, "out_dir"):
                irc_cfg["out_dir"] = str(out_dir)
            if cli_param_overridden(ctx, "tr_projection") and tr_projection is not None:
                geom_cfg["tr_projection"] = str(tr_projection).lower()
            # CLI knobs → irc_cfg. require_pos_def_hessian gates PSD-Hessian convergence.
            if cli_param_overridden(ctx, "irc_pos_def") and irc_pos_def is not None:
                irc_cfg["require_pos_def_hessian"] = bool(irc_pos_def)

            # resolved_charge / resolved_spin already incorporate CLI
            # -q/-m, gjf metadata, and --ligand-charge derivation. Direct
            # assign; an earlier .get("charge", resolved) silently kept
            # the UMA_CALC_KW default 0 when -q was not passed.
            calc_cfg["charge"] = int(resolved_charge)
            calc_cfg["spin"] = int(resolved_spin)

            apply_yaml_overrides(
                override_layer_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",),)),
                    (irc_cfg, (("irc",),)),
                ],
            )
            from pysisyphus.tr_projection import normalize_tr_projection_mode
            geom_cfg["tr_projection"] = normalize_tr_projection_mode(
                geom_cfg.get("tr_projection")
            )

            # Convert 1-based YAML freeze_atoms to 0-based internal
            if geom_cfg.get("freeze_atoms"):
                geom_cfg["freeze_atoms"] = yaml_freeze_to_internal(geom_cfg["freeze_atoms"])
            # Merge CLI --freeze-atoms (already 0-based)
            try:
                freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            if freeze_atoms_cli:
                merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)
            # Normalize freeze_atoms and optionally add link-parent indices for PDB inputs
            resolve_freeze_atoms(geom_cfg, source_path, freeze_links_flag)

            # EulerPC currently only supports Cartesian coordinates
            geom_cfg["coord_type"] = "cart"

            # Ensure the calculator receives the freeze list used by geometry
            # (so FD Hessian can skip frozen DOF, etc.)
            calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
            calc_cfg["return_partial_hessian"] = True

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
                            "tr_projection": geom_cfg["tr_projection"],
                            "will_run_irc": True,
                            "will_write_trajectories": True,
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. IRC execution was skipped.")
                return

            out_dir_path = _prepare_irc_output_dir(
                out_dir_path,
                prefix=str(irc_cfg.get("prefix") or ""),
                protected_inputs=(
                    prepared_input.source_path,
                    geom_input_path,
                    config_yaml,
                    override_yaml,
                ),
            )

            # Default-verbosity entry summary (skipped in child mode).
            from pdb2reaction.core.utils import echo_run_summary
            _model = calc_cfg.get("model")
            echo_run_summary({
                "input": str(input_path),
                "backend": (
                    f"{calc_cfg.get('backend', '?')} ({_model}, {calc_cfg.get('precision', 'fp32')})"
                    if _model else calc_cfg.get("backend", "?")
                ),
                "out": str(out_dir_path),
            })

            # Pretty-print configuration (expand freeze_atoms for readability)
            click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
            click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
            click.echo(pretty_block("irc",  {**irc_cfg, "out_dir": str(out_dir_path)}))

            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)

            geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

            calc = create_calculator(**calc_cfg)
            geometry.set_calculator(calc)
            echo_resolved_device()

            # Seed cached TS Hessian if available (from tsopt in ``all`` workflow)
            from pdb2reaction.io.hessian_cache import (
                discard as _hess_discard,
                load_matching as _hess_load_matching,
                store as _hess_store,
                identity_from_context as _hess_identity,
            )
            # Reuse the tsopt TS Hessian only on a full evaluation-identity
            # match; the all-workflow may pass the TS through a three-decimal
            # PDB, so the coordinate field keeps the wider bohr tolerance.
            cached = _hess_load_matching(
                "ts",
                _hess_identity(geometry, calc_cfg, role="ts"),
                atol=1.1e-3,
            )
            if cached is not None:
                emit("[irc] Reusing cached TS Hessian from tsopt.", narrative=True)
                _dev = calc_cfg.get("device", "auto")
                if _dev == "auto":
                    _dev = "cuda" if torch.cuda.is_available() else "cpu"
                active_dofs = cached.get("active_dofs")
                h_raw = cached["hessian"]
                if isinstance(h_raw, torch.Tensor):
                    h_init = h_raw.to(device=torch.device(_dev))
                else:
                    h_init = torch.as_tensor(h_raw, dtype=torch.float64, device=torch.device(_dev))
                if active_dofs is not None:
                    geometry.within_partial_hessian = {
                        "active_n_dof": len(active_dofs),
                        "full_n_dof": geometry.cart_coords.size,
                        "active_dofs": active_dofs,
                        "active_atoms": sorted(set(d // 3 for d in active_dofs)),
                    }
                geometry.cart_hessian = h_init
                click.echo(f"[irc] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
                del h_init

            eulerpc = EulerPC(geometry, **irc_cfg)
            from pysisyphus.tr_projection import active_tr_basis
            _basis, _rigid_info = active_tr_basis(
                torch.as_tensor(geometry.coords3d, dtype=torch.float64),
                torch.as_tensor(geometry.masses, dtype=torch.float64),
                eulerpc._act_atoms,
                mode=geometry.tr_projection,
            )
            del _basis
            eulerpc.rigid_projection_info = _rigid_info
            click.echo(
                "[irc] Rigid projection: "
                f"treatment={_rigid_info.treatment}, "
                f"rank={_rigid_info.effective_rank}, "
                f"full_rigid_rank={_rigid_info.full_rigid_rank}."
            )

            emit("\n====== IRC (EulerPC) ======\n", narrative=True)
            # Clear per-direction values before running: a failed or one-sided
            # IRC must never expose an endpoint from an earlier segment.
            _hess_discard("irc_left")
            _hess_discard("irc_right")
            eulerpc.run()

            # A one- or two-step branch usually indicates that the Euler step
            # skipped across a shallow feature or immediately triggered an
            # energy stop. Give an actionable remedy in the normal CLI output.
            _quick_directions = []
            for _direction in ("forward", "backward"):
                if not getattr(eulerpc, _direction, False):
                    continue
                _n_frames = len(getattr(eulerpc, f"{_direction}_energies", []))
                if 0 < _n_frames <= 3:
                    _quick_directions.append(_direction)
            if _quick_directions:
                click.echo(
                    "[irc] IRC stopped after only a few frames in "
                    + ", ".join(_quick_directions)
                    + ". Retry with a smaller maximum step, for example "
                    "--step-size 0.05. If a small uphill/flat section is "
                    "intentional, also consider --never-stop; it is opt-in.",
                    err=True,
                )

            # Cache IRC endpoint Hessians (Bofill-updated mw → Cartesian)
            def _unmw_and_store(mw_H, key, endpoint_cart_coords, direction):
                act = eulerpc._act_dofs
                m_sqrt = geometry.masses_rep ** 0.5
                ms_act = m_sqrt[act]
                H_cart_act_np = _consume_mw_hessian_to_cartesian_active(
                    mw_H, ms_act
                )
                _hess_store(
                    key,
                    H_cart_act_np,
                    active_dofs=list(act),
                    meta={
                        "cart_coords": endpoint_cart_coords,
                        "irc_direction": direction,
                    },
                    identity=_hess_identity(
                        geometry,
                        calc_cfg,
                        role=key,
                        cart_coords=endpoint_cart_coords,
                    ),
                )

            # M42: cache an endpoint Hessian ONLY for a requested direction that
            # explicitly CONVERGED. A nonconverged (e.g. max-cycle) direction may
            # still carry a Bofill-updated Hessian, but promoting it would let a
            # nonconverged endpoint seed downstream RFO as if it were a real
            # minimum. Keys were discarded above; keep them discarded when the
            # direction did not converge so no stale/never-converged Hessian is
            # reused.
            from pdb2reaction.workflows._outcomes import (
                irc_hessian_cache_eligible as _irc_hess_eligible,
            )
            _fwd_conv = _irc_hess_eligible(eulerpc, "forward_is_converged")
            _bwd_conv = _irc_hess_eligible(eulerpc, "backward_is_converged")
            if (
                eulerpc.forward
                and _fwd_conv
                and getattr(eulerpc, "forward_mw_hessian", None) is not None
            ):
                _forward_mw_hessian = eulerpc.forward_mw_hessian
                eulerpc.forward_mw_hessian = None
                _forward_endpoint = (
                    np.asarray(eulerpc.forward_mw_coords[0], dtype=float)
                    / np.asarray(eulerpc.m_sqrt, dtype=float)
                )
                try:
                    _unmw_and_store(
                        _forward_mw_hessian,
                        "irc_left",
                        _forward_endpoint,
                        "forward",
                    )
                finally:
                    del _forward_mw_hessian
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                click.echo("[irc] Cached forward endpoint Hessian as 'irc_left'.")
            else:
                eulerpc.forward_mw_hessian = None
                _hess_discard("irc_left")
            if (
                eulerpc.backward
                and _bwd_conv
                and getattr(eulerpc, "mw_hessian", None) is not None
            ):
                _backward_mw_hessian = eulerpc.mw_hessian
                eulerpc.mw_hessian = None
                _backward_endpoint = (
                    np.asarray(eulerpc.backward_mw_coords[-1], dtype=float)
                    / np.asarray(eulerpc.m_sqrt, dtype=float)
                )
                try:
                    _unmw_and_store(
                        _backward_mw_hessian,
                        "irc_right",
                        _backward_endpoint,
                        "backward",
                    )
                finally:
                    del _backward_mw_hessian
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                click.echo("[irc] Cached backward endpoint Hessian as 'irc_right'.")
            else:
                eulerpc.mw_hessian = None
                _hess_discard("irc_right")
            eulerpc.mw_hessian = None
            eulerpc.forward_mw_hessian = None
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            for _stem in ("finished", "forward", "backward"):
                _echo_convert_trj_if_exists(
                    _irc_output_path(eulerpc, f"{_stem}_irc_trj.xyz"),
                    prepared_input,
                    out_pdb=(
                        _irc_output_path(eulerpc, f"{_stem}_irc.pdb")
                        if prepared_input.source_path.suffix.lower() == ".pdb"
                        else None
                    ),
                )

            # summary.md and key_* outputs are disabled.
            emit(format_elapsed("[time] Elapsed Time for IRC", time_start), narrative=True)

            # result.json (if --out-json)
            if out_json:
                from pdb2reaction.core.utils import calculator_provenance, write_result_json
                _all_e = eulerpc.all_energies
                _n_fwd = len(getattr(eulerpc, "forward_energies", [])) if hasattr(eulerpc, "forward_energies") else 0
                _n_bwd = len(getattr(eulerpc, "backward_energies", [])) if hasattr(eulerpc, "backward_energies") else 0
                _ts_e = eulerpc.ts_energy if hasattr(eulerpc, "ts_energy") else None
                _irc_files = {}
                for _fn in ("finished_irc_trj.xyz", "forward_irc_trj.xyz", "backward_irc_trj.xyz"):
                    _fp = _irc_output_path(eulerpc, _fn)
                    if _fp.exists():
                        _irc_files[_fn.replace("_trj.xyz", "")] = _fp.name
                for _fn in (
                    "finished_irc.pdb", "forward_irc.pdb", "backward_irc.pdb",
                    "finished_irc.cif", "forward_irc.cif", "backward_irc.cif",
                ):
                    _fp = _irc_output_path(eulerpc, _fn)
                    if _fp.exists():
                        _irc_files[_fn.replace(".", "_")] = _fp.name
                result_data = {
                    "status": "completed",
                    "n_frames_forward": _n_fwd,
                    "n_frames_backward": _n_bwd,
                    "n_frames_total": len(_all_e),
                    "forward_converged": getattr(eulerpc, 'forward_is_converged', None),
                    "backward_converged": getattr(eulerpc, 'backward_is_converged', None),
                    "forward_energy_increased": getattr(eulerpc, 'forward_energy_increased', None),
                    "backward_energy_increased": getattr(eulerpc, 'backward_energy_increased', None),
                    "backend": calc_cfg.get("backend", backend),
                    "charge": calc_cfg["charge"],
                    "spin": calc_cfg["spin"],
                    "model": calc_cfg.get("model"),
                    **calculator_provenance(calc_cfg),
                    "never_stop": bool(irc_cfg.get("never_stop", False)),
                    "never_stop_energy_bypasses": int(
                        getattr(eulerpc, "never_stop_energy_increase_bypasses", 0)
                        + getattr(eulerpc, "never_stop_energy_convergence_bypasses", 0)
                    ),
                    "rigid_projection": {
                        **getattr(eulerpc, "rigid_projection_info", _rigid_info).as_dict(),
                        "hessian_space": (
                            "active" if len(eulerpc._act_atoms) < len(geometry.atoms) else "full"
                        ),
                        "hessian_shape": list(eulerpc.init_hessian_shape),
                        "hessian_source": "cache" if cached is not None else "fresh",
                        "hessian_representation": "cartesian-unweighted-unprojected",
                    },
                    "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                    "solvent": calc_cfg.get("solvent", "none"),
                    "step_length": irc_cfg.get("step_length"),
                    "max_cycles": irc_cfg.get("max_cycles"),
                    "input_file": str(input_path),
                    "files": _irc_files,
                }
                result_data.update(_directional_endpoint_energy_fields(_all_e, _ts_e))

                # M42: one truthful LeafOutcome per requested IRC direction. A
                # requested direction is usable only when it explicitly
                # converged; a disabled direction is optional (not a failure).
                # Legacy ``status`` stays "completed" (the IRC process ran).
                from pdb2reaction.workflows._outcomes import (
                    aggregate_workflow_truth as _agg_truth,
                    attach_outcomes as _attach,
                    irc_direction_leaves as _irc_dir_leaves,
                )
                _dir_leaves, _dir_expected = _irc_dir_leaves(
                    (
                        (
                            "forward",
                            bool(getattr(eulerpc, "forward", False)),
                            getattr(eulerpc, "forward_is_converged", None),
                            _n_fwd,
                            [_irc_files["forward_irc"]] if "forward_irc" in _irc_files else [],
                        ),
                        (
                            "backward",
                            bool(getattr(eulerpc, "backward", False)),
                            getattr(eulerpc, "backward_is_converged", None),
                            _n_bwd,
                            [_irc_files["backward_irc"]] if "backward_irc" in _irc_files else [],
                        ),
                    )
                )
                _irc_truth = _agg_truth(_dir_leaves, _dir_expected)
                _attach(result_data, truth=_irc_truth, stage_outcomes=_dir_leaves)

                # Bond changes between IRC endpoints
                try:
                    from pdb2reaction.domain.bond_changes import compare_structures
                    _irc_first = _irc_output_path(eulerpc, "finished_first.xyz")
                    _irc_last = _irc_output_path(eulerpc, "finished_last.xyz")
                    if _irc_first.exists() and _irc_last.exists():
                        _g1 = geom_loader(str(_irc_first))
                        _g2 = geom_loader(str(_irc_last))
                        _bc = compare_structures(_g1, _g2, device="cpu")
                        _elems = [a.capitalize() for a in _g1.atoms]
                        result_data["bond_changes"] = {
                            "formed": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.formed_covalent)],
                            "broken": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.broken_covalent)],
                        }
                        result_data["bond_changes_direction"] = "finished_first_to_finished_last"
                except Exception:
                    pass

                write_result_json(
                    out_dir_path, result_data,
                    command="irc",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

        except KeyboardInterrupt:
            click.echo("Interrupted by user.", err=True)
            sys.exit(130)
        except Exception as e:
            render_cli_exception(e, label="IRC", out_dir=out_dir_path, command="irc", time_start=time_start)
        finally:
            # DO NOT INLINE: IRC eulerpc retains a Hessian-sized tensor; explicit
            # `= None` followed by `del` breaks local-frame variable bindings so
            # torch.nn.Module cyclic refs are released before gc.collect() reaps
            # them. Assigning None alone leaves stale frame entries that retain
            # torch hooks, which keeps the Hessian resident across stages.
            calc = eulerpc = geometry = None
            del calc, eulerpc, geometry
            gc.collect()  # break cyclic refs inside torch.nn.Module
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
