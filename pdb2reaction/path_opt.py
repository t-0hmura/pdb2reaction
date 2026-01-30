# pdb2reaction/path_opt.py

"""
Pairwise MEP optimization via GSM or DMF with UMA calculator.

Example:
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1

For detailed documentation, see: docs/path_opt.md
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import sys
import traceback
import textwrap
import time

import click
import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from .uma_pysis import uma_ase, uma_pysis
from .defaults import (
    GEOM_KW_DEFAULT,
    UMA_CALC_KW,
    LBFGS_KW,
    RFO_KW,
    DMF_KW,
    GS_KW,
    STOPT_KW,
)
from .utils import (
    load_yaml_dict,
    apply_yaml_overrides,
    deep_update,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    prepare_input_structure,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    PreparedInputStructure,
    apply_ref_pdb_override,
    resolve_charge_spin_multi,
    load_prepared_geometries,
    write_xyz_trj_with_energy,
)
from .align_freeze_atoms import align_and_refine_sequence_inplace


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

    Uses pydmf (CPU version) with harmonic constraints for frozen atoms.

    References:
    [1] S.-i. Koda and  S. Saito, Locating Transition States by Variational Reaction Path Optimization with an Energy-Derivative-Free Objective Function, JCTC, 20, 2798–2811 (2024). [doi: 10.1021/acs.jctc.3c01246]
    [2] S.-i. Koda and  S. Saito, Flat-bottom Elastic Network Model for Generating Improved Plausible Reaction Paths, JCTC, 20, 7176−7187 (2024). [doi: 10.1021/acs.jctc.4c00792]
    [3] S.-i. Koda and  S. Saito, Correlated Flat-bottom Elastic Network Model for Improved Bond Rearrangement in Reaction Paths, JCTC, 21, 3513−3522 (2025). [doi: 10.1021/acs.jctc.4c01549]
    """

    try:
        from ase.io import read as ase_read
        from ase.io import write as ase_write
        from ase.calculators.mixing import SumCalculator
        from dmf import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise RuntimeError(
            "DMF mode requires ase, fairchem (via uma_pysis), cyipopt, and pydmf "
            "to be installed."
        ) from e

    from .harmonic_constraints import HarmonicFixAtoms

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

    calc_uma = uma_ase(
        model=str(calc_cfg.get("model", "uma-s-1p1")),
        device=str(calc_cfg.get("device", "auto")),
        task_name=str(calc_cfg.get("task_name", "omol")),
        workers=int(calc_cfg.get("workers", 1)),
        workers_per_node=int(calc_cfg.get("workers_per_node", 1)),
    )

    dmf_cfg = deep_update(dict(DMF_KW), dmf_cfg)
    fbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("fbenm_options", {}))
    cfbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("cfbenm_options", {}))
    dmf_opts: Dict[str, Any] = dict(dmf_cfg.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
    k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))

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
    )

    initial_trj = out_dir_path / "dmf_initial.trj"
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
    mxflx.solve(tol="tight")

    calc_eval_kw = dict(calc_cfg)
    if fix_atoms:
        calc_eval_kw.setdefault("freeze_atoms", fix_atoms)

    calc_eval = uma_pysis(**calc_eval_kw)

    energies = []
    for image in mxflx.images:
        elems = image.get_chemical_symbols()
        coords_bohr = np.asarray(image.get_positions(), dtype=float).reshape(-1, 3) * ANG2BOHR
        energies.append(float(calc_eval.get_energy(elems, coords_bohr)["energy"]))
    hei_idx = _select_hei_index(energies)

    final_trj = out_dir_path / "final_geometries.trj"
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

    return DMFMepResult(images=list(mxflx.images), energies=list(energies), hei_idx=int(hei_idx))


def _optimize_single(
    g,
    shared_calc,
    sopt_kind: str,
    sopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    prepared_input: Optional[PreparedInputStructure],
    ref_pdb: Optional[Path],
):
    """
    Single-structure optimization (LBFGS or RFO) shared by path_opt and path_search.
    """
    g.set_calculator(shared_calc)

    seg_dir = out_dir / f"{tag}_{sopt_kind}_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(sopt_cfg)
    args["out_dir"] = str(seg_dir)

    if sopt_kind == "lbfgs":
        opt = LBFGS(g, **args)
    else:
        opt = RFOptimizer(g, **args)

    click.echo(f"\n====== [{tag}] Single-structure {sopt_kind.upper()} started ======\n")
    opt.run()
    click.echo(f"\n====== [{tag}] Single-structure {sopt_kind.upper()} finished ======\n")

    try:
        final_xyz = Path(opt.final_fn)
        if prepared_input is not None:
            ref_pdb_path = ref_pdb
            if ref_pdb_path is None and prepared_input.source_path.suffix.lower() == ".pdb":
                ref_pdb_path = prepared_input.source_path.resolve()
            needs_pdb = ref_pdb_path is not None
            needs_gjf = prepared_input.is_gjf
            if needs_pdb or needs_gjf:
                convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb_path,
                    out_pdb_path=final_xyz.with_suffix(".pdb") if needs_pdb else None,
                    out_gjf_path=final_xyz.with_suffix(".gjf") if needs_gjf else None,
                )
        g_final = geom_loader(final_xyz, coord_type=g.coord_type)
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception as e:
            click.echo(
                f"[path-opt] WARNING: Failed to propagate freeze_atoms to final geometry: {e}",
                err=True,
            )
        g_final.set_calculator(shared_calc)
        return g_final
    except Exception:
        return g


# -----------------------------------------------
# CLI
# -----------------------------------------------

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
    help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=UMA_CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel UMA predictor (workers>1).",
)
@click.option(
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive charge "
        "when -q is omitted (PDB inputs only; otherwise it may fall back to 0)."
    ),
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=1,
    show_default=True,
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--freeze-links",
    "freeze_links_flag",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If a PDB is provided, freeze the parent atoms of link hydrogens.",
)
@click.option(
    "--max-nodes",
    type=int,
    default=10,
    show_default=True,
    help="Number of internal nodes (string has max_nodes+2 images including endpoints).",
)
@click.option(
    "--max-cycles",
    type=int,
    default=300,
    show_default=True,
    help="Maximum optimization cycles.",
)
@click.option(
    "--climb",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Search for a transition state (climbing image) after path growth.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Single-structure optimizer for endpoint preoptimization: light (=LBFGS) or heavy (=RFO).",
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Dump optimizer trajectory/restarts during the run.",
)
@click.option(
    "--convert-files",
    "convert_files",
    type=click.BOOL,
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
    "--out-dir",
    "out_dir",
    type=str,
    default="./result_path_opt/",
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Convergence preset for the string optimizer, pre-alignment refinement, "
        "and endpoint preoptimization (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
    ),
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt, dmf, sopt.lbfgs, sopt.rfo).",
)
@click.option(
    "--preopt",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="If True, preoptimize each endpoint via the selected single-structure optimizer (LBFGS/RFO) before alignment and GSM.",
)
@click.option(
    "--preopt-max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum cycles for endpoint preoptimization (applies to the chosen optimizer; only used when --preopt True).",
)
@click.option(
    "--fix-ends",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Fix structures of input endpoints during GSM.",
)
def cli(
    input_paths: Sequence[Path],
    mep_mode: str,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    preopt_max_cycles: int,
    fix_ends: bool,
) -> None:
    input_paths = tuple(Path(p) for p in input_paths)
    set_convert_file_enabled(convert_files)
    prepared_inputs = [prepare_input_structure(p) for p in input_paths]
    try:
        time_start = time.perf_counter()
        if ref_pdb is not None:
            for prepared in prepared_inputs:
                apply_ref_pdb_override(prepared, ref_pdb)

        # --------------------------
        # 1) Assemble final config (defaults ← CLI ← YAML)
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW_DEFAULT)
        calc_cfg = dict(UMA_CALC_KW)
        dmf_cfg = dict(DMF_KW)
        gs_cfg = dict(GS_KW)
        opt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg = dict(RFO_KW)

        # Resolve charge/spin (defaults ← CLI/GJF/ligand-charge)
        resolved_charge, resolved_spin = resolve_charge_spin_multi(
            prepared_inputs,
            charge=charge,
            spin=spin,
            ligand_charge=ligand_charge,
            prefix="[path-opt]",
        )
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)
        calc_cfg["workers"] = int(workers)
        calc_cfg["workers_per_node"] = int(workers_per_node)

        gs_cfg["max_nodes"] = int(max_nodes)
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)
        gs_cfg["fix_first"] = bool(fix_ends)
        gs_cfg["fix_last"] = bool(fix_ends)

        # Lanczos tangent estimation follows the CLI --climb flag by default but
        # can still be overridden via YAML (`gs.climb_lanczos`).
        opt_cfg["dump"] = bool(dump)
        opt_cfg["out_dir"] = out_dir
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["thresh"] = str(thresh)
            rfo_cfg["thresh"] = str(thresh)

        # Use external Kabsch alignment; keep internal align disabled.
        opt_cfg["align"] = False

        lbfgs_cfg["dump"] = bool(dump)
        rfo_cfg["dump"] = bool(dump)
        lbfgs_cfg["out_dir"] = out_dir
        rfo_cfg["out_dir"] = out_dir

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("sopt", "lbfgs"), ("opt", "lbfgs"), ("lbfgs",))),
                (rfo_cfg, (("sopt", "rfo"), ("opt", "rfo"), ("rfo",))),
            ],
        )

        opt_kind = opt_mode.strip().lower()
        mep_mode_kind = mep_mode.strip().lower()
        if opt_kind == "light":
            sopt_kind = "lbfgs"
            sopt_cfg = lbfgs_cfg
        elif opt_kind == "heavy":
            sopt_kind = "rfo"
            sopt_cfg = rfo_cfg
        else:
            raise click.BadParameter(f"Unknown --opt-mode '{opt_mode}'.")

        sopt_cfg = dict(sopt_cfg)
        sopt_cfg["max_cycles"] = int(preopt_max_cycles)

        # For display: resolved configuration
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_geom_for_echo(calc_cfg)
        echo_gs = dict(gs_cfg)
        echo_opt = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs", echo_gs))
        click.echo(pretty_block("opt", echo_opt))
        if mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg))
        echo_sopt = dict(sopt_cfg)
        echo_sopt["out_dir"] = str(out_dir_path)
        echo_sopt["out_dir_per_tag"] = f"{out_dir_path}/<tag>_{sopt_kind}_opt"
        click.echo(pretty_block("sopt." + sopt_kind, echo_sopt))
        click.echo(
            pretty_block(
                "run_flags",
                {
                    "preopt": bool(preopt),
                    "preopt_max_cycles": int(preopt_max_cycles),
                    "mep_mode": mep_mode_kind,
                },
            )
        )

        # --------------------------
        # 2) Prepare structures (load two endpoints and apply freezing)
        # --------------------------
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
                freeze_effective[prepared.source_path.name] = freeze_list
            click.echo(pretty_block("freeze_atoms (effective)", freeze_effective))

        # Shared UMA calculator (reuse the same instance for all images)
        shared_calc = uma_pysis(**calc_cfg)

        # Optional endpoint pre-optimization (LBFGS/RFO) before alignment/GSM
        if preopt:
            click.echo("\n====== Preoptimizing endpoints via single-structure optimizer started ======\n")
            ref_pdb_for_preopt: Optional[Path] = None
            for p in source_paths:
                if p.suffix.lower() == ".pdb":
                    ref_pdb_for_preopt = p.resolve()
                    break
            if ref_pdb_for_preopt is None and ref_pdb is not None:
                ref_pdb_for_preopt = Path(ref_pdb).resolve()

            preopt_out_dir = out_dir_path
            if (
                out_dir_path.name.startswith("seg_")
                and out_dir_path.parent.name == "path_opt"
            ):
                preopt_out_dir = out_dir_path.parent
                preopt_out_dir.mkdir(parents=True, exist_ok=True)

            new_geoms = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                try:
                    g_opt = _optimize_single(
                        g,
                        shared_calc,
                        sopt_kind,
                        sopt_cfg,
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
                "[preopt] Skipping endpoint preoptimization (use --preopt True to enable)."
            )

        # External Kabsch alignment (if freeze_atoms exist, use only them)
        align_thresh = str(opt_cfg.get("thresh", "gau"))
        try:
            click.echo(
                "\n====== Aligning all inputs to the first structure "
                "(freeze-guided scan + relaxation) started ==="
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
            pass

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

            click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start))
            return

        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # Used when GrowingString generates new nodes
            return shared_calc

        # --------------------------
        # 3) Build path object and optimizer
        # --------------------------
        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        opt_args = dict(opt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"},
        )

        # --------------------------
        # 4) Run optimization
        # --------------------------
        click.echo("\n====== Growing String optimization started ======\n")
        optimizer.run()
        click.echo("\n====== Growing String optimization finished ======\n")

        # --------------------------
        # 5) Write final path (final_geometries.trj)
        # --------------------------
        final_trj = out_dir_path / "final_geometries.trj"
        try:
            try:
                energies = np.array(gs.energy, dtype=float)
                blocks = []
                for idx, (geom, E) in enumerate(zip(gs.images, energies)):
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

        # --------------------------
        # 6) Identify and write HEI (hei.xyz[/pdb/.gjf])
        # --------------------------
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

        click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start))

    except OptimizationError as e:
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo(
            "Unhandled error during path optimization:\n"
            + textwrap.indent(tb, "  "),
            err=True,
        )
        sys.exit(1)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
