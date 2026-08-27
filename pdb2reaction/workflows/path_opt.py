"""
Pairwise MEP optimization via GSM or DMF with an MLIP calculator.

Example:
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1

For detailed documentation, see: docs/path-opt.md
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import gc
import inspect
import logging
import math
import sys
import time

import click
from pdb2reaction.core.output import emit
import numpy as np
import torch

from pysisyphus.constants import ANG2BOHR
from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from pdb2reaction.backends import create_calculator, create_ase_calculator
from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    THRESH_CHOICES,
    UMA_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    OPT_MODE_ALIASES,
    DMF_KW,
    fresh_dmf_config,
    GS_KW,
    STOPT_KW,
    OUT_DIR_PATH_OPT,
    apply_backend_defaults,
)
from pdb2reaction.core.utils import (
    YamlFlowList,
    apply_yaml_overrides,
    deep_update,
    pretty_block,
    format_geom_for_echo,
    emit_dry_run_complete,
    format_elapsed,
    prepare_input_structure,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    convert_xyz_to_pdb,
    PreparedInputStructure,
    apply_ref_pdb_override,
    resolve_charge_spin,
    validate_charge_spin_for_prepared,
    load_prepared_geometries,
    write_xyz_trj_with_energy,
    normalize_choice,
    cli_param_overridden,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
    lossless_int,
    optional_positive_int,
)
from pdb2reaction.workflows.align_freeze import (
    align_and_refine_sequence_inplace,
    alignment_failed_pair_indices,
)
from pdb2reaction.workflows._path_yaml_helpers import apply_single_opt_yaml_layer
from pdb2reaction.cli.common_options import add_coord_type_option, add_precision_option, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, _write_error_json, render_cli_exception

logger = logging.getLogger(__name__)


def _select_hei_index(energies: Sequence[float]) -> int:
    """Return the global highest-energy-image index."""
    E = np.array(energies, dtype=float)
    if E.size == 0:
        raise ValueError("Cannot select an HEI from an empty energy profile.")
    if not np.all(np.isfinite(E)):
        raise ValueError("Cannot select an HEI from non-finite energies.")
    return int(np.argmax(E))


@dataclass
class DMFMepResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int
    # Explicit convergence of the DMF (IPOPT) solve. ``None`` means the
    # engine exposed no readable convergence signal (fail-closed: not "converged"
    # and never promoted to success by artifact existence).
    is_converged: Optional[bool] = None
    reason: str = ""
    ipopt_status: Optional[int] = None


def _is_cuda_oom(exc: BaseException) -> bool:
    """True if `exc` looks like a CUDA out-of-memory (torch.cuda.OutOfMemoryError or a
    RuntimeError carrying 'out of memory'), so the DMF gpu backend can advise --dmf-backend cpu."""
    if type(exc).__name__ == "OutOfMemoryError":
        return True
    msg = str(exc).lower()
    return "out of memory" in msg or "cuda oom" in msg


def _create_dmf_ase_calculator(calc_cfg: Dict[str, Any]):
    """Create the DMF evaluator from the complete resolved calculator map.

    The central ASE factory owns backend-specific key filtering. Copying here
    keeps the caller's canonical configuration immutable while preserving every
    accepted ORB, MACE, UMA, or custom-calculator field.
    """

    return create_ase_calculator(**dict(calc_cfg))


def _main_dmf_ipopt_options(
    ipopt_options: Dict[str, Any],
    out_dir_path: Path,
    max_cycles: Any,
) -> Dict[str, Any]:
    """Resolve IPOPT options for the main DMF solve."""
    resolved = dict(ipopt_options)
    resolved["output_file"] = str(out_dir_path / "dmf_ipopt.out")
    if max_cycles is not None:
        max_iter = lossless_int(max_cycles, "dmf.max_cycles")
        if max_iter < 1:
            raise click.BadParameter("dmf.max_cycles must be a positive integer.")
        resolved["max_iter"] = max_iter
    return resolved


DMF_TOL_PRESETS = ("tight", "middle", "loose")


def resolve_dmf_solve_tol(
    dmf_cfg: Mapping[str, Any], prefix: str = "[path-opt]"
) -> Any:
    """Resolve the ``tol`` argument of the DMF solve from the DMF configuration.

    ``dmf.tol`` accepts the pydmf presets ``tight`` / ``middle`` / ``loose``
    (IPOPT ``dual_inf_tol`` 0.04 / 0.10 / 0.20) or a positive float.  When it is
    unset, an explicitly pinned ``dmf.ipopt_options.dual_inf_tol`` is honoured by
    returning ``None``: pydmf's ``solve`` applies its own ``tol`` after the
    caller's IPOPT options, so a preset passed here would silently replace that
    value.  With neither set, the historical ``tight`` default applies.
    """
    raw = dmf_cfg.get("tol")
    if raw is None:
        pinned = (dmf_cfg.get("ipopt_options") or {}).get("dual_inf_tol")
        return None if pinned is not None else "tight"
    if isinstance(raw, str):
        text = raw.strip().lower()
        if text in DMF_TOL_PRESETS:
            return text
    if isinstance(raw, bool):
        raise click.ClickException(
            f"{prefix} Invalid DMF tolerance '{raw}': expected "
            f"{'|'.join(DMF_TOL_PRESETS)} or a positive float (IPOPT "
            "dual_inf_tol)."
        )
    try:
        value = float(raw)
    except (TypeError, ValueError):
        raise click.ClickException(
            f"{prefix} Invalid DMF tolerance '{raw}': expected "
            f"{'|'.join(DMF_TOL_PRESETS)} or a positive float (IPOPT "
            "dual_inf_tol). Gaussian convergence presets apply to --thresh / "
            "--thresh-gsm, not to the DMF optimizer."
        ) from None
    if not math.isfinite(value) or value <= 0.0:
        raise click.ClickException(
            f"{prefix} Invalid DMF tolerance '{raw}': the IPOPT "
            "dual_inf_tol must be a finite positive number."
        )
    return value


def _combine_path_opt_outcomes(preopt_outcomes, mep_outcome):
    """Combine requested endpoint preoptimizations with the MEP outcome."""
    from pdb2reaction.workflows._outcomes import aggregate_workflow_truth

    outcomes = [*preopt_outcomes, mep_outcome]
    expected = [outcome.item_id for outcome in outcomes]
    return aggregate_workflow_truth(outcomes, expected), outcomes


def _release_dmf_interpolation_cache(mxflx_fbenm: Any) -> None:
    """Release an optional PyDMF device cache at the interpolation boundary."""
    release_cache = getattr(mxflx_fbenm, "release_device_cache", None)
    if callable(release_cache):
        release_cache(empty_cache=False)


def _torch_dmf_runtime_kwargs(
    dmf_backend: str,
    dmf_options: Mapping[str, Any],
    fbenm_options: Mapping[str, Any],
    cfbenm_options: Mapping[str, Any],
    *,
    supports_keep_history: bool = True,
) -> Dict[str, Any]:
    """Resolve Torch-only path settings shared by both DMF stages."""
    if dmf_backend != "gpu":
        return {}

    resolved: Dict[str, Any] = {}
    if supports_keep_history:
        resolved["keep_history"] = bool(dmf_options.get("keep_history", False))
    for name in ("device", "dtype"):
        value = dmf_options.get(name)
        value = fbenm_options.get(name, value)
        value = cfbenm_options.get(name, value)
        if value is not None:
            resolved[name] = value
    if "device" not in resolved:
        # The public gpu backend must run on CUDA.  Without an explicit expert
        # device, dmf.torch would pick its own default and silently compute on
        # the CPU when CUDA is unavailable.
        if not torch.cuda.is_available():
            raise click.ClickException(
                "dmf.backend 'gpu' requires CUDA, which is unavailable. Select "
                "dmf.backend 'cpu' (or --dmf-backend cpu), or set an explicit "
                "dmf.dmf_options.device."
            )
        resolved["device"] = "cuda"
    return resolved


def _validate_dmf_solvent_compatibility(calc_cfg: Mapping[str, Any]) -> None:
    """Reject the unsupported DMF/implicit-solvent PES mismatch."""

    solvent = str(calc_cfg.get("solvent", "none") or "none").strip().lower()
    if solvent not in ("", "none"):
        raise click.ClickException(
            f"--mep-mode dmf is not compatible with --solvent '{solvent}': the DMF path "
            "optimizer runs on the gas-phase ASE PES (no implicit-solvent wrapper), so the "
            "path would be optimized without solvent while the rest of the pipeline uses the "
            "solvent PES. Use --mep-mode gsm (its eval uses the solvent-corrected pysisyphus "
            "calculator), or drop --solvent."
        )


def _run_dmf_mep(
    geoms: Sequence[Any],
    calc_cfg: Dict[str, Any],
    out_dir_path: Path,
    prepared_inputs: Sequence[PreparedInputStructure],
    max_nodes: int,
    fix_atoms: Sequence[int],
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> DMFMepResult:
    """Run Direct Max Flux (DMF) MEP optimization between two endpoints.

    Uses pydmf with harmonic constraints for frozen atoms. The backend is selected by
    ``dmf_cfg["backend"]``: ``"gpu"`` (default) imports the PyTorch backend ``dmf.torch``
    (runs on CUDA), ``"cpu"`` imports the NumPy backend ``dmf``.

    References:
    [1] S.-i. Koda and  S. Saito, Locating Transition States by Variational Reaction Path Optimization with an Energy-Derivative-Free Objective Function, JCTC, 20, 2798–2811 (2024). [doi: 10.1021/acs.jctc.3c01246]
    [2] S.-i. Koda and  S. Saito, Flat-bottom Elastic Network Model for Generating Improved Plausible Reaction Paths, JCTC, 20, 7176−7187 (2024). [doi: 10.1021/acs.jctc.4c00792]
    [3] S.-i. Koda and  S. Saito, Correlated Flat-bottom Elastic Network Model for Improved Bond Rearrangement in Reaction Paths, JCTC, 21, 3513−3522 (2025). [doi: 10.1021/acs.jctc.4c01549]
    """

    # DMF optimizes on the ASE-path PES, which has no implicit-solvent wrapper
    # (solvent is a pysisyphus-path-only xTB correction). Under --solvent the DMF
    # path would be optimized gas-phase while the rest of the pipeline uses the
    # solvent PES — an opt-PES ≠ score-PES mismatch. Refuse and direct the user to
    # GSM, whose per-image eval goes through the solvent-aware pysisyphus calculator.
    _validate_dmf_solvent_compatibility(calc_cfg)

    dmf_backend = str((dmf_cfg or {}).get("backend", "gpu")).strip().lower()
    if dmf_backend not in {"cpu", "gpu"}:
        raise click.ClickException(
            "dmf.backend must be either 'cpu' or 'gpu'; "
            f"received {dmf_backend!r}."
        )
    try:
        from ase.io import read as ase_read
        from ase.io import write as ase_write
        from ase.calculators.mixing import SumCalculator
        if dmf_backend == "cpu":
            from dmf import DirectMaxFlux, interpolate_fbenm
        else:
            from dmf.torch import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise click.ClickException(
            "DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2 "
            "(`pip install 'pydmf[torch]'` for the default GPU backend; the `cpu` "
            "backend needs only `pip install pydmf`). "
            f"Import error: {e}"
        ) from e

    from pdb2reaction.workflows.restraints import HarmonicFixAtoms

    def _geom_to_ase(g: Any):
        from io import StringIO

        return ase_read(StringIO(g.as_xyz()), format="xyz")

    fix_atoms = list(sorted(set(map(int, fix_atoms))))

    ref_images = [_geom_to_ase(g) for g in geoms]
    primary_prepared = prepared_inputs[0] if prepared_inputs else None
    ref_pdb = (
        primary_prepared.source_path.resolve()
        if primary_prepared and primary_prepared.source_path.suffix.lower() == ".pdb"
        else None
    )
    needs_pdb = ref_pdb is not None
    needs_gjf = bool(primary_prepared and primary_prepared.is_gjf)

    charge = int(calc_cfg.get("charge", 0))
    spin = int(calc_cfg.get("spin", 1))
    for img in ref_images:
        img.info["charge"] = charge
        img.info["spin"] = spin

    ase_calc = _create_dmf_ase_calculator(calc_cfg)

    dmf_cfg = fresh_dmf_config(dmf_cfg)
    fbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("fbenm_options", {}))
    cfbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("cfbenm_options", {}))
    dmf_opts: Dict[str, Any] = dict(dmf_cfg.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
    supports_keep_history = (
        dmf_backend == "gpu"
        and "keep_history" in inspect.signature(DirectMaxFlux.__init__).parameters
    )
    torch_dmf_kwargs = _torch_dmf_runtime_kwargs(
        dmf_backend, dmf_opts, fbenm_opts, cfbenm_opts,
        supports_keep_history=supports_keep_history,
    )
    if supports_keep_history:
        # Product workflows consume only the final path.  Keeping every
        # iteration would retain force arrays and Atoms copies unnecessarily.
        dmf_opts.setdefault("keep_history", False)
    k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))

    # Default-mode IPOPT options: print_level=0 silences the per-iteration
    # IPOPT banner + MUMPS license header that would otherwise print 16+
    # times per pipeline (one per align+scan refinement). Under `-v` we
    # let dmf use its own default (print_level=5 = full table). User-
    # supplied `dmf_cfg["ipopt_options"]` always wins.
    from pdb2reaction.core.utils import is_verbose
    ipopt_opts: Dict[str, Any] = dict(dmf_cfg.get("ipopt_options", {}))
    if "print_level" not in ipopt_opts and not is_verbose():
        ipopt_opts["print_level"] = 0

    # Run FB-ENM interpolation (pydmf CPU version)
    mxflx_fbenm = interpolate_fbenm(
        ref_images,
        nmove=max(1, int(max_nodes)),
        fbenm_only_endpoints=bool(dmf_cfg.get("fbenm_only_endpoints", False)),
        correlated=bool(dmf_cfg.get("correlated", False)),
        sequential=bool(dmf_cfg.get("sequential", False)),
        output_file=str(out_dir_path / "dmf_fbenm_ipopt.out"),
        fbenm_options=fbenm_opts,
        cfbenm_options=cfbenm_opts,
        dmf_options=dmf_opts,
        ipopt_options=ipopt_opts,
    )

    initial_trj = out_dir_path / "dmf_initial_trj.xyz"
    ase_write(initial_trj, mxflx_fbenm.images, format="xyz")
    if primary_prepared is not None and needs_pdb:
        convert_xyz_like_outputs(
            initial_trj,
            primary_prepared,
            ref_pdb_path=ref_pdb,
            out_pdb_path=initial_trj.with_suffix(".pdb") if needs_pdb else None,
        )
    coefs = mxflx_fbenm.coefs.copy()

    # The interpolation stage owns a separate, large FB-ENM/CFB-ENM cache.
    # Release it only at this phase boundary: keeping it through the accurate
    # PES stage multiplies peak VRAM, while releasing it inside force sweeps
    # would repeat the CPU-to-GPU constant uploads.
    _release_dmf_interpolation_cache(mxflx_fbenm)
    del mxflx_fbenm
    gc.collect()
    if dmf_backend == "gpu" and torch.cuda.is_available():
        torch.cuda.empty_cache()

    # Create DirectMaxFlux object (pydmf CPU version)
    mxflx = DirectMaxFlux(
        ref_images,
        coefs=coefs,
        nmove=max(1, int(max_nodes)),
        update_teval=update_teval,
        remove_rotation_and_translation=bool(
            dmf_opts.get("remove_rotation_and_translation", False)
        ),
        mass_weighted=bool(dmf_opts.get("mass_weighted", False)),
        parallel=bool(dmf_opts.get("parallel", False)),
        eps_vel=float(dmf_opts.get("eps_vel", 0.01)),
        eps_rot=float(dmf_opts.get("eps_rot", 0.01)),
        beta=float(dmf_opts.get("beta", 10.0)),
        **torch_dmf_kwargs,
    )

    # Assign calculators to images
    # For frozen atoms, use HarmonicFixAtoms combined with UMA via SumCalculator
    for image in mxflx.images:
        if "charge" not in image.info:
            image.info["charge"] = charge
        if "spin" not in image.info:
            image.info["spin"] = spin

        if fix_atoms:
            # Create harmonic constraint calculator for frozen atoms
            ref_positions = image.get_positions()[fix_atoms]
            harmonic_calc = HarmonicFixAtoms(
                indices=fix_atoms,
                ref_positions=ref_positions,
                k_fix=k_fix,
            )
            # Combine the configured MLIP calculator with harmonic constraints.
            image.calc = SumCalculator([ase_calc, harmonic_calc])
        else:
            image.calc = ase_calc

    max_cycles = dmf_cfg.get("max_cycles") if isinstance(dmf_cfg, dict) else None
    main_ipopt_opts = _main_dmf_ipopt_options(
        ipopt_opts,
        out_dir_path,
        max_cycles,
    )
    mxflx.add_ipopt_options(main_ipopt_opts)
    # IPOPT returns (x_opt, info); info["status"] == 0 is Solve_Succeeded and
    # 1 is Solved_To_Acceptable_Level. Anything else (max-iter, infeasible) is
    # NOT convergence. A missing/unreadable status fails closed to unknown.
    _solve_ret = mxflx.solve(tol=resolve_dmf_solve_tol(dmf_cfg))
    _info = (
        _solve_ret[1]
        if isinstance(_solve_ret, (tuple, list)) and len(_solve_ret) >= 2
        else None
    )
    from pdb2reaction.workflows._outcomes import ipopt_status_to_converged
    _dmf_status: Optional[int] = None
    if isinstance(_info, dict) and _info.get("status") is not None:
        try:
            _dmf_status = int(_info["status"])
        except (TypeError, ValueError):
            _dmf_status = None
    _dmf_conv, _dmf_reason = ipopt_status_to_converged(_dmf_status)

    # Free DMF calculator before creating eval calculator to avoid GPU OOM
    for image in mxflx.images:
        image.calc = None
    del ase_calc
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    calc_eval_kw = dict(calc_cfg)
    if fix_atoms:
        calc_eval_kw.setdefault("freeze_atoms", fix_atoms)

    calc_eval = create_calculator(**calc_eval_kw)

    energies = []
    for image in mxflx.images:
        elems = image.get_chemical_symbols()
        coords_bohr = np.asarray(image.get_positions(), dtype=float).reshape(-1, 3) * ANG2BOHR
        energies.append(float(calc_eval.get_energy(elems, coords_bohr)["energy"]))
    hei_idx = _select_hei_index(energies)

    final_trj = out_dir_path / "final_geometries_trj.xyz"
    write_xyz_trj_with_energy(mxflx.images, energies, final_trj)
    if primary_prepared is not None and needs_pdb:
        convert_xyz_like_outputs(
            final_trj,
            primary_prepared,
            ref_pdb_path=ref_pdb,
            out_pdb_path=final_trj.with_suffix(".pdb") if needs_pdb else None,
        )
    if primary_prepared is not None and (needs_pdb or needs_gjf):
        hei_tmp = out_dir_path / "hei.xyz"
        write_xyz_trj_with_energy([mxflx.images[hei_idx]], [energies[hei_idx]], hei_tmp)
        convert_xyz_like_outputs(
            hei_tmp,
            primary_prepared,
            ref_pdb_path=ref_pdb,
            out_pdb_path=out_dir_path / "hei.pdb" if needs_pdb else None,
            out_gjf_path=out_dir_path / "hei.gjf" if needs_gjf else None,
        )

    # Free eval calculator before returning (caller may create a new one)
    result = DMFMepResult(
        images=list(mxflx.images),
        energies=list(energies),
        hei_idx=int(hei_idx),
        is_converged=_dmf_conv,
        reason=_dmf_reason,
        ipopt_status=_dmf_status,
    )
    del calc_eval, mxflx
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return result


def _optimize_single(
    g,
    shared_calc,
    opt_kind: str,
    single_opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    prepared_input: Optional[PreparedInputStructure],
    ref_pdb: Optional[Path],
):
    """
    Single-structure optimization (LBFGS or RFO) shared by path_opt and path_search.

    Returns ``(geometry, converged)`` where ``converged`` is the optimizer's
    fail-closed tri-state convergence bit (``optimizer_converged_bit``): a
    non-raising ``run()`` is not convergence, so the caller must gate on
    this bit rather than assume a completed optimization converged.
    """
    g.set_calculator(shared_calc)

    seg_dir = out_dir / f"{tag}_{opt_kind}_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(single_opt_cfg)
    args["out_dir"] = str(seg_dir)

    if opt_kind == "lbfgs":
        opt = LBFGS(g, **args)
    else:
        opt = RFOptimizer(g, **args)

    emit(f"\n====== [{tag}] Single-structure {opt_kind.upper()} ======\n", narrative=True)
    opt.run()
    # Capture the optimizer's explicit convergence bit before any
    # geometry post-processing. A normal (non-raising) run() is NOT convergence;
    # thread this to the caller so a nonconverged single-structure optimization
    # cannot silently become a usable (e.g. kink-segment) reactive leaf.
    from pdb2reaction.workflows._outcomes import optimizer_converged_bit
    converged = optimizer_converged_bit(opt)

    try:
        final_xyz = Path(opt.final_fn)
        if prepared_input is not None:
            ref_pdb_path = ref_pdb
            if ref_pdb_path is None and prepared_input.source_path.suffix.lower() == ".pdb":
                ref_pdb_path = prepared_input.source_path.resolve()
            needs_pdb = ref_pdb_path is not None
            needs_gjf = prepared_input.is_gjf
            if needs_pdb or needs_gjf:
                converted = convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb_path,
                    out_pdb_path=final_xyz.with_suffix(".pdb") if needs_pdb else None,
                    out_gjf_path=final_xyz.with_suffix(".gjf") if needs_gjf else None,
                )
                # Fallback: if the source wasn't PDB (e.g. .xyz from a scan),
                # convert_xyz_like_outputs skips PDB output.  Use ref_pdb directly.
                if not converted and needs_pdb and not final_xyz.with_suffix(".pdb").exists():
                    try:
                        convert_xyz_to_pdb(final_xyz, ref_pdb_path, final_xyz.with_suffix(".pdb"))
                    except Exception as e:
                        click.echo(
                            f"[{tag}] WARNING: PDB fallback conversion failed: {e}",
                            err=True,
                        )
        elif ref_pdb is not None:
            # No prepared_input but ref_pdb is available — still try PDB conversion
            try:
                out_pdb = final_xyz.with_suffix(".pdb")
                if not out_pdb.exists():
                    convert_xyz_to_pdb(final_xyz, ref_pdb, out_pdb)
            except Exception as e:
                click.echo(
                    f"[{tag}] WARNING: PDB conversion failed: {e}",
                    err=True,
                )
        g_final = geom_loader(
            final_xyz,
            coord_type=g.coord_type,
            freeze_atoms=getattr(g, "freeze_atoms", []),
        )
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception as e:
            click.echo(
                f"[path-opt] WARNING: Failed to propagate freeze_atoms to final geometry: {e}",
                err=True,
            )
        g_final.set_calculator(shared_calc)
        return g_final, converged
    except Exception as exc:
        click.echo(
            f"[path-opt] WARNING: Failed to load optimized geometry; returning input: {exc}",
            err=True,
        )
        # Could not even load the optimized geometry: return the input geometry
        # and fail closed on convergence (None) rather than claim the discarded
        # optimization converged.
        return g, None



@click.command(
    help="MEP optimization via GSM or DMF.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    nargs=2,
    required=True,
    help="Two endpoint structures (reactant and product); accepts PDB, mmCIF, XYZ, or GJF.",
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
)
@click.option(
    "--dmf-backend",
    type=click.Choice(["cpu", "gpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="DMF compute backend (--mep-mode dmf only): gpu (dmf.torch / CUDA) or cpu (dmf / NumPy). "
    "On a GPU out-of-memory error, retry with cpu.",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=False,
    help="Total charge. Required unless a .gjf template provides charge metadata or --ligand-charge is supplied for PDB/mmCIF inputs.",
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
    "-l",
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
        "when -q is omitted (requires PDB/mmCIF input or --ref-pdb)."
    ),
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default="1",
    help="Spin multiplicity (2S+1).",
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
    "--max-nodes",
    type=int,
    default=GS_KW["max_nodes"],
    show_default=True,
    help=("Number of movable internal images for GSM or DMF; the complete path "
          "has max_nodes+2 images including the two endpoints."),
)
@click.option(
    "--gsm-param",
    type=click.Choice(["equi", "energy"], case_sensitive=False),
    default=None,
    show_default="equi",
    help=(
        "GSM node parameterization after string growth. The energy scheme "
        "concentrates nodes in high-energy regions and may be tried when an "
        "equidistant path skips the reaction-coordinate region near the HEI."
    ),
)
@click.option(
    "--max-cycles-gsm",
    type=click.IntRange(min=1),
    default=None,
    show_default="300",
    help="Maximum GSM string-optimizer cycles for the MEP stage.",
)
@click.option(
    "--max-cycles-dmf",
    type=click.IntRange(min=1),
    default=None,
    show_default="300",
    help=(
        "Maximum IPOPT iterations for the DMF MEP stage. This is a solver "
        "iteration count, not a string-optimizer cycle count."
    ),
)
@click.option(
    "--climb/--no-climb",
    default=True,
    show_default=True,
    help="Enable the GSM climbing-image search after path growth (accepted but unused by DMF).",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="grad",
    show_default=True,
    help="Single-structure optimizer for endpoint preoptimization: grad (=LBFGS) or hess (=RFO).",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write GSM/single-optimizer trajectory and restart data (accepted but unused by DMF).",
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
    "out_dir",
    type=str,
    default=OUT_DIR_PATH_OPT,
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau",
    help=(
        "Convergence preset for endpoint preoptimization only "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option(
    "--thresh-gsm",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau_loose",
    help=(
        "Convergence preset for the GSM string optimizer "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option(
    "--thresh-dmf",
    type=str,
    default=None,
    show_default="tight",
    help=(
        "IPOPT dual-infeasibility tolerance for the DMF path optimizer: "
        "tight (0.04) | middle (0.10) | loose (0.20) or a positive float. "
        "This is not a Gaussian preset."
    ),
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
    help="Validate options and print the execution plan without running path optimization.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "--preopt/--no-preopt",
    default=True,
    show_default=True,
    help="Preoptimize each endpoint via the selected single-structure optimizer (LBFGS/RFO) before alignment and path optimization.",
)
@click.option(
    "--preopt-max-cycles",
    type=click.IntRange(min=1),
    default=None,
    show_default="100000",
    help="Maximum cycles for each endpoint preoptimization pass.",
)
@click.option(
    "--fix-ends/--no-fix-ends",
    default=True,
    show_default=True,
    help="Fix input endpoints during GSM path optimization (accepted but unused by DMF).",
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Experimental, computationally expensive xTB solvent delta correction. Examples: water, methanol, acetonitrile, dmso, thf, toluene. 'none' disables it.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              show_default=True, help="xTB solvent model.")
@add_coord_type_option(choices=("cart", "dlc"))
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    mep_mode: str,
    dmf_backend: str,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links_flag: bool,
    freeze_atoms_text: Optional[str],
    max_nodes: int,
    gsm_param: Optional[str],
    max_cycles_gsm: Optional[int],
    max_cycles_dmf: Optional[int],
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    thresh_gsm: Optional[str],
    thresh_dmf: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    preopt: bool,
    preopt_max_cycles: Optional[int],
    fix_ends: bool,
    backend: str,
    solvent: str,
    solvent_model: str,
    cli_coord_type: Optional[str],
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
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

    input_paths = tuple(Path(p) for p in input_paths)
    set_convert_file_enabled(convert_files)
    prepared_inputs = []
    try:
        for path in input_paths:
            prepared_inputs.append(prepare_input_structure(path))
    except BaseException:
        for prepared in reversed(prepared_inputs):
            prepared.cleanup()
        raise
    # Initialize early so that early-failure except handlers do not raise
    # an UnboundLocalError when referencing out_dir_path / time_start, which
    # would mask the real exception. Reassigned to the resolved value once
    # the run config is parsed.
    out_dir_path: Optional[Path] = None
    time_start: Optional[float] = time.perf_counter()
    try:
        if ref_pdb is not None:
            for prepared in prepared_inputs:
                apply_ref_pdb_override(prepared, ref_pdb)

        geom_cfg = dict(GEOM_KW_DEFAULT)
        calc_cfg = dict(UMA_CALC_KW)
        dmf_cfg = fresh_dmf_config()
        gs_cfg = dict(GS_KW)
        stopt_cfg = dict(STOPT_KW)
        stopt_cfg["out_dir"] = out_dir
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg = dict(RFO_KW)

        def _apply_single_opt_yaml_layer(layer_cfg: Dict[str, Any]) -> None:
            apply_single_opt_yaml_layer(
                layer_cfg,
                lbfgs_cfg=lbfgs_cfg,
                rfo_cfg=rfo_cfg,
                opt_base_kw=OPT_BASE_KW,
                deep_update=deep_update,
                apply_yaml_overrides=apply_yaml_overrides,
            )

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",),)),
            ],
        )
        _apply_single_opt_yaml_layer(config_layer_cfg)

        # Resolve charge/spin from templates/ligand charge and then apply precedence.
        resolved_charge, resolved_spin = resolve_charge_spin(
            prepared_inputs,
            charge=charge,
            spin=spin,
            ligand_charge=ligand_charge,
            prefix="[path-opt]",
        )
        validate_charge_spin_for_prepared(prepared_inputs, resolved_charge, resolved_spin)
        # resolved_charge / resolved_spin already incorporate CLI -q/-m,
        # gjf metadata, and --ligand-charge derivation. Direct assign;
        # an earlier .get("charge", resolved) silently kept the UMA_CALC_KW
        # default 0 when -q was not passed.
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        if cli_param_overridden(ctx, "workers"):
            calc_cfg["workers"] = int(workers)
        if cli_param_overridden(ctx, "workers_per_node"):
            calc_cfg["workers_per_node"] = int(workers_per_node)
        if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
            geom_cfg["coord_type"] = str(cli_coord_type).lower()
        if cli_param_overridden(ctx, "max_nodes"):
            gs_cfg["max_nodes"] = int(max_nodes)
        if cli_param_overridden(ctx, "gsm_param") and gsm_param is not None:
            gs_cfg["param"] = str(gsm_param).lower()
        # The GSM cycle budget also bounds the fully-grown string; DMF's budget
        # is a separate IPOPT iteration count.
        if cli_param_overridden(ctx, "max_cycles_gsm") and max_cycles_gsm is not None:
            stopt_cfg["max_cycles"] = int(max_cycles_gsm)
            stopt_cfg["stop_in_when_full"] = int(max_cycles_gsm)
        if cli_param_overridden(ctx, "max_cycles_dmf") and max_cycles_dmf is not None:
            dmf_cfg["max_cycles"] = int(max_cycles_dmf)
        if cli_param_overridden(ctx, "dmf_backend"):
            dmf_cfg["backend"] = str(dmf_backend).lower()
        if cli_param_overridden(ctx, "climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if cli_param_overridden(ctx, "fix_ends"):
            if str(mep_mode).lower() == "dmf" and not fix_ends:
                raise click.BadParameter(
                    "--no-fix-ends is supported only for GSM paths.",
                    param_hint="--no-fix-ends",
                )
            gs_cfg["fix_first"] = bool(fix_ends)
            gs_cfg["fix_last"] = bool(fix_ends)
        if cli_param_overridden(ctx, "dump"):
            stopt_cfg["dump"] = bool(dump)
            lbfgs_cfg["dump"] = bool(dump)
            rfo_cfg["dump"] = bool(dump)
        if cli_param_overridden(ctx, "out_dir"):
            stopt_cfg["out_dir"] = out_dir
            lbfgs_cfg["out_dir"] = out_dir
            rfo_cfg["out_dir"] = out_dir
        if cli_param_overridden(ctx, "thresh") and thresh is not None:
            lbfgs_cfg["thresh"] = str(thresh)
            rfo_cfg["thresh"] = str(thresh)
        if cli_param_overridden(ctx, "thresh_gsm") and thresh_gsm is not None:
            stopt_cfg["thresh"] = str(thresh_gsm)
        if cli_param_overridden(ctx, "thresh_dmf") and thresh_dmf is not None:
            dmf_cfg["tol"] = str(thresh_dmf)
        if cli_param_overridden(ctx, "preopt_max_cycles") and preopt_max_cycles is not None:
            lbfgs_cfg["max_cycles"] = int(preopt_max_cycles)
            rfo_cfg["max_cycles"] = int(preopt_max_cycles)

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",),)),
            ],
        )
        _apply_single_opt_yaml_layer(override_layer_cfg)

        mep_mode_kind = mep_mode.strip().lower()

        # A dormant YAML DMF section does not affect GSM. An explicit CLI
        # tolerance is still validated as user input, regardless of MEP mode.
        if mep_mode_kind == "dmf" or cli_param_overridden(ctx, "thresh_dmf"):
            resolve_dmf_solve_tol(dmf_cfg)

        # Convert 1-based YAML freeze_atoms to 0-based internal
        if geom_cfg.get("freeze_atoms"):
            geom_cfg["freeze_atoms"] = yaml_freeze_to_internal(geom_cfg["freeze_atoms"])
        # Merge CLI --freeze-atoms (already 0-based)
        try:
            freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
        except click.BadParameter:
            raise
        if freeze_atoms_cli:
            merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)

        # Use external Kabsch alignment; keep internal align disabled.
        stopt_cfg["align"] = False

        opt_kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=OPT_MODE_ALIASES,
            allowed_hint="grad|hess",
        )
        if opt_kind == "lbfgs":
            single_opt_kind = "lbfgs"
            single_opt_cfg = lbfgs_cfg
        else:
            single_opt_kind = "rfo"
            single_opt_cfg = rfo_cfg

        single_opt_cfg = dict(single_opt_cfg)
        preopt_max_cycles_effective = single_opt_cfg.get(
            "max_cycles", preopt_max_cycles
        )
        if mep_mode_kind == "dmf":
            dmf_cycles = optional_positive_int(dmf_cfg.get("max_cycles"), "dmf.max_cycles")
            dmf_cfg["max_cycles"] = dmf_cycles
            try:
                k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))
            except (TypeError, ValueError, OverflowError) as exc:
                raise click.BadParameter(
                    "dmf.k_fix must be finite and non-negative."
                ) from exc
            if not np.isfinite(k_fix) or k_fix < 0.0:
                raise click.BadParameter(
                    "dmf.k_fix must be finite and non-negative."
                )
            dmf_cfg["k_fix"] = k_fix
        else:
            string_cycles = optional_positive_int(
                stopt_cfg.get("max_cycles"), "stopt.max_cycles"
            )
            stopt_cfg["max_cycles"] = string_cycles
            stopt_cfg["stop_in_when_full"] = string_cycles
        if preopt:
            preopt_max_cycles_effective = optional_positive_int(
                single_opt_cfg.get("max_cycles"), "preopt max_cycles"
            )

        # Apply backend/solvent CLI overrides early (before display)
        if cli_param_overridden(ctx, "backend"):
            calc_cfg["backend"] = backend
        if cli_param_overridden(ctx, "solvent"):
            calc_cfg["solvent"] = solvent
        if cli_param_overridden(ctx, "solvent_model"):
            calc_cfg["solvent_model"] = solvent_model
        from pdb2reaction.backends import apply_backend_model_to_calc_cfg
        # Unconditional: also pops a raw backend_model token from a --config YAML
        # (the helper no-ops when neither the CLI arg nor the YAML names one).
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        from pdb2reaction.backends import apply_calc_file_to_calc_cfg
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        from pdb2reaction.backends import apply_effective_precision
        apply_effective_precision(calc_cfg, precision)
        apply_backend_defaults(calc_cfg)

        # For display: resolved configuration
        out_dir_path = Path(stopt_cfg["out_dir"]).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_geom_for_echo(calc_cfg)
        echo_gs = dict(gs_cfg)
        echo_stopt = dict(stopt_cfg)
        echo_stopt["out_dir"] = str(out_dir_path)

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs", echo_gs))
        click.echo(pretty_block("stopt", echo_stopt))
        if mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg))
        echo_opt = dict(single_opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)
        echo_opt["out_dir_per_tag"] = f"{out_dir_path}/<tag>_{single_opt_kind}_opt"
        click.echo(pretty_block("opt." + single_opt_kind, echo_opt))
        click.echo(
            pretty_block(
                "run_flags",
                {
                    "preopt": bool(preopt),
                    "preopt_max_cycles": preopt_max_cycles_effective,
                    "mep_mode": mep_mode_kind,
                },
            )
        )

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override_yaml": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                force=True)
            )
        if dry_run:
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_endpoints": [str(p) for p in input_paths],
                        "output_dir": str(out_dir_path),
                        "mep_mode": mep_mode_kind,
                        "gsm_param": str(gs_cfg.get("param", GS_KW["param"])),
                        "opt_mode": ("grad" if opt_kind == "lbfgs" else "hess"),
                        "preopt": bool(preopt),
                        "preopt_max_cycles": preopt_max_cycles_effective,
                        "freeze_links": bool(freeze_links_flag),
                        "convert_files": bool(convert_files),
                        "will_run_path_opt": True,
                        "will_write_summary": True,
                    },
                )
            )
            emit_dry_run_complete()
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)
        for name in (
            "hei.pdb",
            "hei.cif",
            "hei.gjf",
            "final_geometries.pdb",
            "final_geometries.cif",
        ):
            (out_dir_path / name).unlink(missing_ok=True)
        if not out_json:
            for name in ("result.json", "summary.json"):
                (out_dir_path / name).unlink(missing_ok=True)

        source_paths = [prep.source_path for prep in prepared_inputs]

        geoms = load_prepared_geometries(
            prepared_inputs,
            coord_type=geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"]),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )
        if geoms:
            freeze_effective: Dict[str, List[int]] = {}
            for prepared, g in zip(prepared_inputs, geoms):
                try:
                    freeze_list = list(getattr(g, "freeze_atoms", []))
                except Exception:
                    freeze_list = list(geom_cfg.get("freeze_atoms", []))
                freeze_effective[prepared.source_path.name] = YamlFlowList([i + 1 for i in freeze_list])
            click.echo(pretty_block("freeze_atoms (effective, 1-based)", freeze_effective))

        # Keep UMA freeze_atoms aligned with the resolved geometry freeze list.
        if geoms:
            freeze_union = sorted(
                {int(i) for g in geoms for i in getattr(g, "freeze_atoms", [])}
            )
            calc_cfg["freeze_atoms"] = freeze_union

        # (backend/solvent CLI overrides already applied above, before display)

        # Shared calculator (reuse the same instance for all images)
        shared_calc = create_calculator(**calc_cfg)
        echo_resolved_device()

        # Optional endpoint pre-optimization (LBFGS/RFO) before alignment/GSM
        preopt_outcomes = []
        if preopt:
            from pdb2reaction.workflows._outcomes import make_leaf as _mk_leaf

            emit("\n====== Preoptimizing endpoints via single-structure optimizer ======\n", narrative=True)
            ref_pdb_for_preopt: Optional[Path] = None
            for p in source_paths:
                if p.suffix.lower() == ".pdb":
                    ref_pdb_for_preopt = p.resolve()
                    break
            if ref_pdb_for_preopt is None and ref_pdb is not None:
                ref_pdb_for_preopt = Path(ref_pdb).resolve()

            preopt_out_dir = out_dir_path
            _seg_prefix = ""
            if (
                out_dir_path.name.startswith("seg_")
                and out_dir_path.parent.name == "path_opt"
            ):
                preopt_out_dir = out_dir_path.parent
                preopt_out_dir.mkdir(parents=True, exist_ok=True)
                # Include segment name in tag to avoid overwriting across segments
                _seg_prefix = out_dir_path.name.split("_mep")[0] + "_"

            new_geoms = []
            for i, g in enumerate(geoms):
                tag = f"{_seg_prefix}init{i:02d}"
                try:
                    g_opt, preopt_converged = _optimize_single(
                        g,
                        shared_calc,
                        single_opt_kind,
                        single_opt_cfg,
                        preopt_out_dir,
                        tag=tag,
                        prepared_input=prepared_inputs[i] if i < len(prepared_inputs) else None,
                        ref_pdb=ref_pdb_for_preopt,
                    )
                    new_geoms.append(g_opt)
                    preopt_outcomes.append(
                        _mk_leaf(
                            "path-opt",
                            f"preopt_endpoint_{i}",
                            executed=True,
                            converged=preopt_converged,
                        )
                    )
                except Exception as e:
                    click.echo(
                        f"[preopt] WARNING: Failed to preoptimize endpoint {i}: {e}",
                        err=True,
                    )
                    new_geoms.append(g)
                    preopt_outcomes.append(
                        _mk_leaf(
                            "path-opt",
                            f"preopt_endpoint_{i}",
                            executed=True,
                            converged=False,
                            reason="optimization_error",
                        )
                    )
            geoms = new_geoms
        else:
            click.echo(
                "[preopt] Skipping endpoint preoptimization as requested by --no-preopt."
            )

        # External Kabsch alignment (if freeze_atoms exist, use only them)
        align_thresh = str(single_opt_cfg.get("thresh", "gau"))
        try:
            click.echo(
                "\n====== Aligning all inputs to the first structure "
                "(freeze-guided scan + relaxation) ======\n"
            )
            alignment_results = align_and_refine_sequence_inplace(
                geoms,
                thresh=align_thresh,
                shared_calc=shared_calc,
                out_dir=out_dir_path / "align_refine",
                verbose=True,
            )
            failed_pairs = alignment_failed_pair_indices(alignment_results)
            if failed_pairs:
                raise click.ClickException(
                    "Input alignment did not converge for pair(s): "
                    + ", ".join(str(index) for index in failed_pairs)
                )
            click.echo("[align] Completed input alignment.")
        except Exception as e:
            raise click.ClickException(f"Input alignment failed: {e}") from e

        fix_atoms: List[int] = []
        try:
            fix_atoms = sorted(
                {int(i) for g in geoms for i in getattr(g, "freeze_atoms", [])}
            )
        except Exception:
            logger.warning("Failed to collect freeze_atoms from geometries", exc_info=True)

        if mep_mode_kind == "dmf":
            try:
                dmf_res = _run_dmf_mep(
                    geoms,
                    calc_cfg,
                    out_dir_path,
                    prepared_inputs,
                    int(gs_cfg["max_nodes"]),
                    fix_atoms,
                    dmf_cfg=dmf_cfg,
                )
            except Exception as e:
                if str(dmf_cfg.get("backend", "gpu")).lower() != "cpu" and _is_cuda_oom(e):
                    click.echo(
                        "[dmf] GPU out of memory. Retry with `--dmf-backend cpu` "
                        "(NumPy backend; slower but not limited by GPU memory).",
                        err=True,
                    )
                else:
                    click.echo(f"[dmf] ERROR: DMF optimization failed: {e}", err=True)
                sys.exit(3)

            try:
                hei_idx = int(dmf_res.hei_idx)
                hei_xyz = out_dir_path / "hei.xyz"
                write_xyz_trj_with_energy(
                    [dmf_res.images[hei_idx]], [dmf_res.energies[hei_idx]], hei_xyz
                )
                click.echo(f"[write] Wrote '{hei_xyz}'.")
                main_prepared = prepared_inputs[0] if prepared_inputs else None
                if main_prepared is not None:
                    ref_pdb = (
                        main_prepared.source_path.resolve()
                        if main_prepared.source_path.suffix.lower() == ".pdb"
                        else None
                    )
                    needs_pdb = ref_pdb is not None
                    needs_gjf = main_prepared.is_gjf
                    if needs_pdb or needs_gjf:
                        try:
                            convert_xyz_like_outputs(
                                hei_xyz,
                                main_prepared,
                                ref_pdb_path=ref_pdb,
                                out_pdb_path=out_dir_path / "hei.pdb" if needs_pdb else None,
                                out_gjf_path=out_dir_path / "hei.gjf" if needs_gjf else None,
                            )
                            click.echo("[convert] Wrote 'hei' outputs.")
                        except Exception as e:
                            click.echo(
                                f"[convert] WARNING: Failed to convert HEI to requested formats: {e}",
                                err=True,
                            )
            except Exception as e:
                click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
                sys.exit(5)

            # result.json (if --out-json) — DMF path
            if out_json:
                from pdb2reaction.core.utils import calculator_provenance, write_result_json
                from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
                _dmf_energies = list(dmf_res.energies)
                _dmf_hei = int(dmf_res.hei_idx)
                _dmf_hei_E = float(_dmf_energies[_dmf_hei])
                _dmf_e0 = float(_dmf_energies[0])
                _dmf_eN = float(_dmf_energies[-1])
                _barrier = (_dmf_hei_E - _dmf_e0) * _AU2KCAL
                _delta = (_dmf_eN - _dmf_e0) * _AU2KCAL
                _dmf_converged = getattr(dmf_res, 'is_converged', None)
                _dmf_reason = getattr(dmf_res, 'reason', "") or ""
                result_data: Dict[str, Any] = {
                    # Legacy byte compatibility: DMFMepResult exposed no
                    # readable is_converged, so both the legacy `status` and the
                    # legacy `converged` fields always read "completed" / null.
                    # The new IPOPT convergence truth is carried ONLY by the
                    # additive scientific_status / stage_outcomes (attached below);
                    # both legacy values are pinned to their base so no downstream
                    # consumer of `status`/`converged` observes a non-additive flip
                    # on a genuinely converged run.
                    "status": "completed",
                    "converged": None,
                    "mep_mode": "dmf",
                    "backend": calc_cfg.get("backend", backend),
                    "charge": calc_cfg["charge"],
                    "spin": calc_cfg["spin"],
                    "model": calc_cfg.get("model"),
                    **calculator_provenance(calc_cfg),
                    "solvent": calc_cfg.get("solvent", "none"),
                    "preopt": bool(preopt),
                    "reactant_energy_hartree": float(_dmf_e0),
                    "product_energy_hartree": float(_dmf_eN),
                    "image_energies_hartree": [float(e) for e in _dmf_energies],
                    "n_images": len(_dmf_energies),
                    "hei_index": _dmf_hei,
                    "hei_energy_hartree": _dmf_hei_E,
                    "barrier_kcal": round(_barrier, 6),
                    "delta_kcal": round(_delta, 6),
                    "files": {
                        "final_geometries_trj_xyz": "final_geometries_trj.xyz",
                        "hei_xyz": "hei.xyz",
                    },
                }
                for ext in (".pdb", ".cif", ".gjf"):
                    f = out_dir_path / f"hei{ext}"
                    if f.exists():
                        result_data["files"][f"hei_{ext[1:]}"] = f.name
                for ext in (".pdb", ".cif"):
                    f = out_dir_path / f"final_geometries{ext}"
                    if f.exists():
                        result_data["files"][f"final_geometries_{ext[1:]}"] = f.name
                # Additive outcome fields: the DMF path is a required leaf that
                # is usable only when the IPOPT solve explicitly converged. A
                # nonconverged solve keeps its trajectory artifact but is not
                # promoted to a usable path.
                from pdb2reaction.workflows._outcomes import (
                    attach_outcomes as _attach,
                    make_leaf as _mk_leaf,
                )
                _dmf_leaf = _mk_leaf(
                    "path-opt",
                    "dmf_mep",
                    executed=True,
                    converged=_dmf_converged,
                    artifacts=["final_geometries_trj.xyz"],
                    reason=_dmf_reason,
                )
                _dmf_truth, _dmf_outcomes = _combine_path_opt_outcomes(
                    preopt_outcomes,
                    _dmf_leaf,
                )
                _attach(
                    result_data,
                    truth=_dmf_truth,
                    stage_outcomes=_dmf_outcomes,
                )
                write_result_json(
                    out_dir_path, result_data,
                    command="path-opt",
                    elapsed_seconds=time.perf_counter() - time_start,
                )
            emit(
                format_elapsed("[time] Elapsed Time for Path Opt", time_start),
                narrative=True,
            )
            return

        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # Used when GrowingString generates new nodes
            return shared_calc

        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        opt_args = dict(stopt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"},
        )

        emit("\n====== Growing String optimization ======\n", narrative=True)
        optimizer.run()

        final_trj = out_dir_path / "final_geometries_trj.xyz"
        try:
            try:
                energies = np.array(gs.energy, dtype=float)
                blocks = []
                for _idx, (geom, E) in enumerate(zip(gs.images, energies)):
                    s = geom.as_xyz()
                    lines = s.splitlines()
                    if len(lines) >= 2 and lines[0].strip().isdigit():
                        lines[1] = f"{E:.12f}"
                    s_mod = "\n".join(lines)
                    if not s_mod.endswith("\n"):
                        s_mod += "\n"
                    blocks.append(s_mod)
                annotated = "".join(blocks)
                with open(final_trj, "w") as f:
                    f.write(annotated)
                click.echo(f"[write] Wrote '{final_trj}' with energy.")
            except Exception:
                with open(final_trj, "w") as f:
                    f.write(gs.as_xyz())
                click.echo(f"[write] Wrote '{final_trj}'.")

            main_prepared = prepared_inputs[0]
            needs_pdb = main_prepared.source_path.suffix.lower() == ".pdb"
            needs_gjf = main_prepared.is_gjf
            ref_pdb = main_prepared.source_path.resolve() if needs_pdb else None
            if needs_pdb or needs_gjf:
                try:
                    convert_xyz_like_outputs(
                        final_trj,
                        main_prepared,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=out_dir_path / "final_geometries.pdb" if needs_pdb else None,
                        out_gjf_path=out_dir_path / "final_geometries.gjf" if needs_gjf else None,
                    )
                    click.echo("[convert] Wrote 'final_geometries' outputs.")
                except Exception as e:
                    click.echo(
                        f"[convert] WARNING: Failed to convert MEP path trajectory: {e}",
                        err=True,
                    )

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        try:
            energies = np.array(gs.energy, dtype=float)
            hei_idx = _select_hei_index(energies)

            hei_geom = gs.images[int(hei_idx)]
            hei_E = float(energies[int(hei_idx)])

            hei_xyz = out_dir_path / "hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            main_prepared = prepared_inputs[0]
            needs_pdb = main_prepared.source_path.suffix.lower() == ".pdb"
            needs_gjf = main_prepared.is_gjf
            ref_pdb = main_prepared.source_path.resolve() if needs_pdb else None
            if needs_pdb or needs_gjf:
                try:
                    convert_xyz_like_outputs(
                        hei_xyz,
                        main_prepared,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=out_dir_path / "hei.pdb" if needs_pdb else None,
                        out_gjf_path=out_dir_path / "hei.gjf" if needs_gjf else None,
                    )
                    click.echo("[convert] Wrote 'hei' outputs.")
                except Exception as e:
                    click.echo(
                        f"[convert] WARNING: Failed to convert HEI structure: {e}",
                        err=True,
                    )
            else:
                click.echo("[convert] Skipped HEI conversion (no PDB/CIF/GJF template).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

        # result.json (if --out-json) — GSM path
        if out_json:
            from pdb2reaction.core.utils import calculator_provenance, write_result_json
            from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
            _gsm_energies = list(map(float, energies))
            _gsm_hei = int(hei_idx)
            _gsm_hei_E = float(_gsm_energies[_gsm_hei])
            _gsm_e0 = float(_gsm_energies[0])
            _gsm_eN = float(_gsm_energies[-1])
            _barrier = (_gsm_hei_E - _gsm_e0) * _AU2KCAL
            _delta = (_gsm_eN - _gsm_e0) * _AU2KCAL
            result_data_gsm: Dict[str, Any] = {
                "mep_mode": "gsm",
                "backend": calc_cfg.get("backend", backend),
                "charge": calc_cfg["charge"],
                "spin": calc_cfg["spin"],
                "model": calc_cfg.get("model"),
                **calculator_provenance(calc_cfg),
                "solvent": calc_cfg.get("solvent", "none"),
                "preopt": bool(preopt),
                "reactant_energy_hartree": float(_gsm_e0),
                "product_energy_hartree": float(_gsm_eN),
                "image_energies_hartree": [float(e) for e in _gsm_energies],
                "n_images": len(_gsm_energies),
                "hei_index": _gsm_hei,
                "hei_energy_hartree": _gsm_hei_E,
                "barrier_kcal": round(_barrier, 6),
                "delta_kcal": round(_delta, 6),
                "files": {
                    "final_geometries_trj_xyz": "final_geometries_trj.xyz",
                    "hei_xyz": "hei.xyz",
                },
            }
            for ext in (".pdb", ".cif", ".gjf"):
                f = out_dir_path / f"hei{ext}"
                if f.exists():
                    result_data_gsm["files"][f"hei_{ext[1:]}"] = f.name
            for ext in (".pdb", ".cif"):
                f = out_dir_path / f"final_geometries{ext}"
                if f.exists():
                    result_data_gsm["files"][f"final_geometries_{ext[1:]}"] = f.name
            # Additive outcome fields: the GSM path is usable only when the
            # optimizer explicitly converged.
            from pdb2reaction.workflows._outcomes import (
                attach_outcomes as _attach,
                make_leaf as _mk_leaf,
                optimizer_converged_bit as _optimizer_converged_bit,
            )
            _converged = _optimizer_converged_bit(optimizer)
            result_data_gsm["status"] = (
                "converged"
                if _converged is True
                else ("not_converged" if _converged is False else "completed")
            )
            result_data_gsm["converged"] = _converged
            _gsm_leaf = _mk_leaf(
                "path-opt",
                "gsm_mep",
                executed=True,
                converged=_converged,
                artifacts=["final_geometries_trj.xyz"],
            )
            _gsm_truth, _gsm_outcomes = _combine_path_opt_outcomes(
                preopt_outcomes,
                _gsm_leaf,
            )
            _attach(
                result_data_gsm,
                truth=_gsm_truth,
                stage_outcomes=_gsm_outcomes,
            )
            write_result_json(
                out_dir_path, result_data_gsm,
                command="path-opt",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        emit(
            format_elapsed("[time] Elapsed Time for Path Opt", time_start),
            narrative=True,
        )

    except OptimizationError as e:
        _write_error_json(out_dir_path, "path-opt", e, "OptimizationError", time_start)
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except click.ClickException:
        raise
    except Exception as e:
        render_cli_exception(e, label="path optimization", out_dir=out_dir_path, command="path-opt", time_start=time_start)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM
        shared_calc = gs = geoms = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
