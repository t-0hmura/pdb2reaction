"""
Pairwise MEP optimization via GSM or DMF with UMA calculator.

Example:
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1

For detailed documentation, see: docs/path-opt.md
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import gc
import logging
import sys
import time

import click
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
    GS_KW,
    STOPT_KW,
    OUT_DIR_PATH_OPT,
    DEFAULT_UMA_MODEL,
    apply_backend_defaults,
)
from pdb2reaction.core.utils import (
    YamlFlowList,
    apply_yaml_overrides,
    deep_update,
    pretty_block,
    format_geom_for_echo,
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
)
from pdb2reaction.workflows.align_freeze import align_and_refine_sequence_inplace
from pdb2reaction.workflows._path_yaml_helpers import apply_single_opt_yaml_layer
from pdb2reaction.cli.common_options import add_coord_type_option, add_precision_option, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, _write_error_json, render_cli_exception

logger = logging.getLogger(__name__)


def _select_hei_index(energies: Sequence[float]) -> int:
    """Pick an HEI index preferring internal local maxima."""

    E = np.array(energies, dtype=float)
    nE = int(len(E))
    hei_idx = None
    if nE >= 3:
        candidates = [i for i in range(1, nE - 1) if (E[i] > E[i - 1] and E[i] > E[i + 1])]
        if candidates:
            hei_idx = int(max(candidates, key=lambda i: E[i]))
        else:
            hei_idx = 1 + int(np.argmax(E[1:-1]))
    if hei_idx is None:
        hei_idx = int(np.argmax(E))
    return hei_idx


@dataclass
class DMFMepResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int


def _is_cuda_oom(exc: BaseException) -> bool:
    """True if `exc` looks like a CUDA out-of-memory (torch.cuda.OutOfMemoryError or a
    RuntimeError carrying 'out of memory'), so the DMF gpu backend can advise --dmf-backend cpu."""
    if type(exc).__name__ == "OutOfMemoryError":
        return True
    msg = str(exc).lower()
    return "out of memory" in msg or "cuda oom" in msg


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
    _solvent = str(calc_cfg.get("solvent", "none") or "none").strip().lower()
    if _solvent not in ("", "none"):
        raise click.ClickException(
            f"--mep-mode dmf is not compatible with --solvent '{_solvent}': the DMF path "
            "optimizer runs on the gas-phase ASE PES (no implicit-solvent wrapper), so the "
            "path would be optimized without solvent while the rest of the pipeline uses the "
            "solvent PES. Use --mep-mode gsm (its eval uses the solvent-corrected pysisyphus "
            "calculator), or drop --solvent."
        )

    dmf_backend = str((dmf_cfg or {}).get("backend", "gpu")).strip().lower()
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

    calc_uma = create_ase_calculator(
        backend=calc_cfg.get("backend", "uma"),
        model=str(calc_cfg.get("model", DEFAULT_UMA_MODEL)),
        device=str(calc_cfg.get("device", "auto")),
        task_name=str(calc_cfg.get("task_name", "omol")),
        workers=int(calc_cfg.get("workers", 1)),
        workers_per_node=int(calc_cfg.get("workers_per_node", 1)),
        # Match the DMF optimizer PES to the pysisyphus HEI-ranking PES (calc_eval
        # below carries the full calc_cfg incl. precision). `solvent` has no ASE
        # equivalent, so DMF cannot run solvent-corrected — that case is rejected
        # up front (see the --solvent guard at the top of this function), hence
        # only precision needs matching here.
        precision=str(calc_cfg.get("precision", "fp32")),
    )

    dmf_cfg = deep_update(dict(DMF_KW), dmf_cfg)
    fbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("fbenm_options", {}))
    cfbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("cfbenm_options", {}))
    dmf_opts: Dict[str, Any] = dict(dmf_cfg.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
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
            # Combine UMA calculator with harmonic constraints
            image.calc = SumCalculator([calc_uma, harmonic_calc])
        else:
            image.calc = calc_uma

    mxflx.add_ipopt_options({"output_file": str(out_dir_path / "dmf_ipopt.out")})
    # Same default-mode IPOPT quiet as the FB-ENM interpolation above:
    # the banner + MUMPS license header would otherwise reprint per solve.
    if not is_verbose() and "print_level" not in ipopt_opts:
        mxflx.add_ipopt_options({"print_level": 0})
    max_cycles = dmf_cfg.get("max_cycles") if isinstance(dmf_cfg, dict) else None
    if max_cycles is not None:
        try:
            max_iter = int(max_cycles)
            if max_iter > 0:
                mxflx.add_ipopt_options({"max_iter": max_iter})
        except Exception:
            pass
    mxflx.solve(tol="tight")

    # Free DMF calculator before creating eval calculator to avoid GPU OOM
    for image in mxflx.images:
        image.calc = None
    calc_uma = None
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
    result = DMFMepResult(images=list(mxflx.images), energies=list(energies), hei_idx=int(hei_idx))
    calc_eval = mxflx = None
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

    click.echo(f"\n====== [{tag}] Single-structure {opt_kind.upper()} ======\n", narrative=True)
    opt.run()

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
        return g_final
    except Exception as exc:
        click.echo(
            f"[path-opt] WARNING: Failed to load optimized geometry; returning input: {exc}",
            err=True,
        )
        return g



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
    help="Two endpoint structures (reactant and product); accepts .pdb or .xyz.",
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
    help="Total charge. Required unless a .gjf template provides charge metadata or --ligand-charge is supplied for PDB inputs.",
)
@click.option(
    "--workers",
    type=int,
    default=UMA_CALC_KW["workers"],
    show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: when workers>1 the analytical Hessian is unavailable and auto-downgrades to finite differences with a warning; run with --workers 1 to force analytical.",
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
        "when -q is omitted (requires PDB input or --ref-pdb)."
    ),
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1).",
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
    "--max-nodes",
    type=int,
    default=20,
    show_default=True,
    help="Maximum number of internal nodes (string has up to max_nodes+2 images including endpoints).",
)
@click.option(
    "--max-cycles",
    type=int,
    default=300,
    show_default=True,
    help="Maximum string optimizer cycles (GSM/DMF path optimization).",
)
@click.option(
    "--climb/--no-climb",
    default=True,
    show_default=True,
    help="Search for a transition state (climbing image) after path growth.",
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
    help="Write optimizer trajectory and restarts during the run.",
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
    show_default=False,
    help=(
        "Convergence preset for endpoint preoptimization only "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
    ),
)
@click.option(
    "--thresh-stopt",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Convergence preset for the string optimizer (stopt) "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau_loose' when not provided."
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
    type=int,
    default=10000,
    show_default=True,
    help="Maximum cycles for each endpoint preoptimization pass (LBFGS or RFO; only used when --preopt is enabled).",
)
@click.option(
    "--fix-ends/--no-fix-ends",
    default=True,
    show_default=True,
    help="Fix structures of input endpoints during path optimization (GSM/DMF).",
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.")
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
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    thresh_stopt: str,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    preopt: bool,
    preopt_max_cycles: int,
    fix_ends: bool,
    backend: str,
    solvent: str,
    solvent_model: str,
    cli_coord_type: Optional[str],
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: str,
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

    input_paths = tuple(Path(p) for p in input_paths)
    set_convert_file_enabled(convert_files)
    prepared_inputs = [prepare_input_structure(p) for p in input_paths]
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
        dmf_cfg = dict(DMF_KW)
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
        if cli_param_overridden(ctx, "max_cycles"):
            stopt_cfg["max_cycles"] = int(max_cycles)
            stopt_cfg["stop_in_when_full"] = int(max_cycles)
            dmf_cfg["max_cycles"] = int(max_cycles)
        if cli_param_overridden(ctx, "dmf_backend"):
            dmf_cfg["backend"] = str(dmf_backend).lower()
        if cli_param_overridden(ctx, "climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if cli_param_overridden(ctx, "fix_ends"):
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
        if cli_param_overridden(ctx, "thresh_stopt") and thresh_stopt is not None:
            stopt_cfg["thresh"] = str(thresh_stopt)
        if cli_param_overridden(ctx, "preopt_max_cycles"):
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

        # Use external Kabsch alignment; keep internal align disabled.
        stopt_cfg["align"] = False
        stopt_cfg["stop_in_when_full"] = int(
            stopt_cfg.get("max_cycles", STOPT_KW["max_cycles"])
        )

        opt_kind = normalize_choice(
            opt_mode,
            param="--opt-mode",
            alias_groups=OPT_MODE_ALIASES,
            allowed_hint="grad|hess",
        )
        mep_mode_kind = mep_mode.strip().lower()
        if opt_kind == "lbfgs":
            single_opt_kind = "lbfgs"
            single_opt_cfg = lbfgs_cfg
        else:
            single_opt_kind = "rfo"
            single_opt_cfg = rfo_cfg

        single_opt_cfg = dict(single_opt_cfg)
        preopt_max_cycles_effective = int(
            single_opt_cfg.get("max_cycles", preopt_max_cycles)
        )

        # Apply backend/solvent CLI overrides early (before display)
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
                    "preopt_max_cycles": int(preopt_max_cycles_effective),
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
                )
            )
        if dry_run:
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_endpoints": [str(p) for p in input_paths],
                        "output_dir": str(out_dir_path),
                        "mep_mode": mep_mode_kind,
                        "opt_mode": ("grad" if opt_kind == "lbfgs" else "hess"),
                        "preopt": bool(preopt),
                        "preopt_max_cycles": int(preopt_max_cycles_effective),
                        "freeze_links": bool(freeze_links_flag),
                        "convert_files": bool(convert_files),
                        "will_run_path_opt": True,
                        "will_write_summary": True,
                    },
                )
            )
            click.echo("[dry-run] Validation complete. Path optimization execution was skipped.")
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)

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
        if preopt:
            click.echo("\n====== Preoptimizing endpoints via single-structure optimizer ======\n", narrative=True)
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
                    g_opt = _optimize_single(
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
                except Exception as e:
                    click.echo(
                        f"[preopt] WARNING: Failed to preoptimize endpoint {i}: {e}",
                        err=True,
                    )
                    new_geoms.append(g)
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
            _ = align_and_refine_sequence_inplace(
                geoms,
                thresh=align_thresh,
                shared_calc=shared_calc,
                out_dir=out_dir_path / "align_refine",
                verbose=True,
            )
            click.echo("[align] Completed input alignment.")
        except Exception as e:
            click.echo(f"[align] WARNING: alignment skipped: {e}", err=True)

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
                    max_nodes,
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

            click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start), narrative=True)

            # result.json (if --out-json) — DMF path
            if out_json:
                from pdb2reaction.core.utils import write_result_json
                from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
                _dmf_energies = list(dmf_res.energies)
                _dmf_hei = int(dmf_res.hei_idx)
                _dmf_hei_E = float(_dmf_energies[_dmf_hei])
                _dmf_e0 = float(_dmf_energies[0])
                _dmf_eN = float(_dmf_energies[-1])
                _barrier = (_dmf_hei_E - _dmf_e0) * _AU2KCAL
                _delta = (_dmf_eN - _dmf_e0) * _AU2KCAL
                _dmf_converged = getattr(dmf_res, 'is_converged', None)
                result_data: Dict[str, Any] = {
                    "status": "converged" if _dmf_converged else ("not_converged" if _dmf_converged is False else "completed"),
                    "converged": _dmf_converged,
                    "mep_mode": "dmf",
                    "backend": calc_cfg.get("backend", backend),
                    "charge": calc_cfg["charge"],
                    "spin": calc_cfg["spin"],
                    "model": calc_cfg.get("model"),
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
                for ext in (".pdb", ".gjf"):
                    f = out_dir_path / f"hei{ext}"
                    if f.exists():
                        result_data["files"][f"hei_{ext[1:]}"] = f.name
                write_result_json(
                    out_dir_path, result_data,
                    command="path-opt",
                    elapsed_seconds=time.perf_counter() - time_start,
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

        click.echo("\n====== Growing String optimization ======\n", narrative=True)
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
                click.echo("[convert] Skipped HEI conversion (no PDB/GJF template).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

        click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start), narrative=True)

        # result.json (if --out-json) — GSM path
        if out_json:
            from pdb2reaction.core.utils import write_result_json
            from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
            _gsm_energies = list(map(float, energies))
            _gsm_hei = int(hei_idx)
            _gsm_hei_E = float(_gsm_energies[_gsm_hei])
            _gsm_e0 = float(_gsm_energies[0])
            _gsm_eN = float(_gsm_energies[-1])
            _barrier = (_gsm_hei_E - _gsm_e0) * _AU2KCAL
            _delta = (_gsm_eN - _gsm_e0) * _AU2KCAL
            _converged = getattr(optimizer, 'is_converged', None) if 'optimizer' in dir() else None
            result_data_gsm: Dict[str, Any] = {
                "status": "converged" if _converged else ("not_converged" if _converged is False else "completed"),
                "converged": _converged,
                "mep_mode": "gsm",
                "backend": calc_cfg.get("backend", backend),
                "charge": calc_cfg["charge"],
                "spin": calc_cfg["spin"],
                "model": calc_cfg.get("model"),
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
            for ext in (".pdb", ".gjf"):
                f = out_dir_path / f"hei{ext}"
                if f.exists():
                    result_data_gsm["files"][f"hei_{ext[1:]}"] = f.name
            write_result_json(
                out_dir_path, result_data_gsm,
                command="path-opt",
                elapsed_seconds=time.perf_counter() - time_start,
            )

    except OptimizationError as e:
        _write_error_json(out_dir_path, "path-opt", e, "OptimizationError", time_start)
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
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
