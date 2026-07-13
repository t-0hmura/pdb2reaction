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
    help="Input structure file (.pdb, .xyz, _trj.xyz, etc.).",
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
        "Maximum number of IRC steps; used unless YAML sets irc.max_cycles. "
        "Defaults to 125 when not provided."
    ),
)
@click.option(
    "--step-size",
    type=float,
    default=None,
    help=(
        "Step length in Bohr (unweighted Cartesian coordinates); used unless YAML sets irc.step_length. "
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
        "still stop the run. Used unless YAML sets irc.never_stop; default off."
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
    help="Freeze parent atoms of cap hydrogens (PDB input or XYZ/GJF with --ref-pdb).",
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
    help="How the ML backend builds the Hessian (Analytical or FiniteDifference); used unless YAML sets calc.hessian_calc_mode. Defaults to 'FiniteDifference'.",
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
    calc_factory: str,
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
        try:
            time_start = time.perf_counter()

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
                            "will_run_irc": True,
                            "will_write_trajectories": True,
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. IRC execution was skipped.")
                return

            out_dir_path.mkdir(parents=True, exist_ok=True)

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
                load as _hess_load,
                matches_cart_coords as _hess_matches_coords,
                store as _hess_store,
            )
            cached = _hess_load("ts")
            if cached is not None and not _hess_matches_coords(
                cached,
                geometry.cart_coords,
                # The all-workflow may pass the TS through a three-decimal PDB.
                atol=1.1e-3,
            ):
                click.echo(
                    "[irc] Cached TS Hessian does not match the IRC start "
                    "geometry; calculating a fresh Hessian.",
                    err=True,
                )
                cached = None
            if cached is not None:
                click.echo("[irc] Reusing cached TS Hessian from tsopt.")
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

            click.echo("\n====== IRC (EulerPC) ======\n", narrative=True)
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
                if isinstance(mw_H, torch.Tensor):
                    ms_t = torch.as_tensor(ms_act, dtype=mw_H.dtype, device=mw_H.device)
                    H_cart_act = ms_t.unsqueeze(1) * mw_H * ms_t.unsqueeze(0)
                    H_cart_act_np = H_cart_act.detach().cpu().numpy()
                    # C-NEED: disk cache contract; free GPU copy after npy dump.
                    del H_cart_act
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                else:
                    H_cart_act_np = np.diag(ms_act) @ mw_H @ np.diag(ms_act)
                _hess_store(
                    key,
                    H_cart_act_np,
                    active_dofs=list(act),
                    meta={
                        "cart_coords": endpoint_cart_coords,
                        "irc_direction": direction,
                    },
                )

            if eulerpc.forward and getattr(eulerpc, "forward_mw_hessian", None) is not None:
                _forward_endpoint = (
                    np.asarray(eulerpc.forward_mw_coords[0], dtype=float)
                    / np.asarray(eulerpc.m_sqrt, dtype=float)
                )
                _unmw_and_store(
                    eulerpc.forward_mw_hessian,
                    "irc_left",
                    _forward_endpoint,
                    "forward",
                )
                click.echo("[irc] Cached forward endpoint Hessian as 'irc_left'.")
            if eulerpc.backward and getattr(eulerpc, "mw_hessian", None) is not None:
                _backward_endpoint = (
                    np.asarray(eulerpc.backward_mw_coords[-1], dtype=float)
                    / np.asarray(eulerpc.m_sqrt, dtype=float)
                )
                _unmw_and_store(
                    eulerpc.mw_hessian,
                    "irc_right",
                    _backward_endpoint,
                    "backward",
                )
                click.echo("[irc] Cached backward endpoint Hessian as 'irc_right'.")

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
            click.echo(format_elapsed("[time] Elapsed Time for IRC", time_start), narrative=True)

            # result.json (if --out-json)
            if out_json:
                from pdb2reaction.core.utils import write_result_json
                _all_e = eulerpc.all_energies
                _n_fwd = len(getattr(eulerpc, "forward_energies", [])) if hasattr(eulerpc, "forward_energies") else 0
                _n_bwd = len(getattr(eulerpc, "backward_energies", [])) if hasattr(eulerpc, "backward_energies") else 0
                _ts_e = float(eulerpc.ts_energy) if hasattr(eulerpc, "ts_energy") else None
                # Forward endpoint = first element of all_energies (reactant side)
                # Backward endpoint = last element of all_energies (product side)
                _e_reactant = float(_all_e[0]) if len(_all_e) > 0 else None
                _e_product = float(_all_e[-1]) if len(_all_e) > 0 else None
                _irc_files = {}
                for _fn in ("finished_irc_trj.xyz", "forward_irc_trj.xyz", "backward_irc_trj.xyz"):
                    _fp = out_dir_path / f"{suffix_prefix}{_fn}"
                    if _fp.exists():
                        _irc_files[_fn.replace("_trj.xyz", "")] = _fp.name
                for _fn in ("finished_irc.pdb", "forward_irc.pdb", "backward_irc.pdb"):
                    _fp = out_dir_path / f"{suffix_prefix}{_fn}"
                    if _fp.exists():
                        _irc_files[_fn.replace(".pdb", "_pdb")] = _fp.name
                result_data = {
                    "status": "completed",
                    "n_frames_forward": _n_fwd,
                    "n_frames_backward": _n_bwd,
                    "n_frames_total": len(_all_e),
                    "energy_reactant_hartree": _e_reactant,
                    "energy_ts_hartree": _ts_e,
                    "energy_product_hartree": _e_product,
                    "forward_converged": getattr(eulerpc, 'forward_is_converged', None),
                    "backward_converged": getattr(eulerpc, 'backward_is_converged', None),
                    "forward_energy_increased": getattr(eulerpc, 'forward_energy_increased', None),
                    "backward_energy_increased": getattr(eulerpc, 'backward_energy_increased', None),
                    "backend": calc_cfg.get("backend", backend),
                    "charge": calc_cfg["charge"],
                    "spin": calc_cfg["spin"],
                    "model": calc_cfg.get("model"),
                    "never_stop": bool(irc_cfg.get("never_stop", False)),
                    "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                    "solvent": calc_cfg.get("solvent", "none"),
                    "step_length": irc_cfg.get("step_length"),
                    "max_cycles": irc_cfg.get("max_cycles"),
                    "input_file": str(input_path),
                    "files": _irc_files,
                }

                # Bond changes between IRC endpoints
                try:
                    from pdb2reaction.domain.bond_changes import compare_structures
                    _irc_first = out_dir_path / "finished_first.xyz"
                    _irc_last = out_dir_path / "finished_last.xyz"
                    if _irc_first.exists() and _irc_last.exists():
                        _g1 = geom_loader(str(_irc_first))
                        _g2 = geom_loader(str(_irc_last))
                        _bc = compare_structures(_g1, _g2, device="cpu")
                        _elems = [a.capitalize() for a in _g1.atoms]
                        result_data["bond_changes"] = {
                            "formed": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.formed_covalent)],
                            "broken": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.broken_covalent)],
                        }
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
