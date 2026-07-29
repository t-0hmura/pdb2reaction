# pdb2reaction/scan.py

"""
Staged bond-length scan with harmonic restraints and full relaxation.

Example:
    pdb2reaction scan -i input.pdb -q 0 --scan-lists '[(12,45,1.35)]' --preopt --endopt

For detailed documentation, see: docs/scan.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import logging
import math
import sys
import textwrap

import click
from pdb2reaction.core.output import emit
import numpy as np
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR

from pdb2reaction.core.defaults import (
    GEOM_KW_DEFAULT,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    BIAS_KW,
    BOND_KW,
    OPT_MODE_ALIASES,
    UMA_CALC_KW,
    OUT_DIR_SCAN,
    apply_backend_defaults,
)
from pdb2reaction.backends import create_calculator
from pdb2reaction.workflows.restraints import HarmonicBiasCalculator
from pdb2reaction.workflows._outcomes import optimizer_converged_bit
from pdb2reaction.core.utils import (
    optimizer_terminal_status,
    build_sopt_kwargs,
    make_sopt_optimizer,
    build_scan_configs,
    cli_param_overridden,
    pretty_block,
    strip_inherited_keys,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    make_snapshot_geometry,
    resolve_freeze_atoms,
    xyz_string_with_energy,
    unbiased_energy_hartree,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
)
from pdb2reaction.cli.decorators import (
    load_merged_yaml_cfg,
    resolve_yaml_sources,
    run_cli,
)
from pdb2reaction.cli.common_options import (
    add_coord_type_option,
    add_print_every_option,
    add_precision_option, add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from pdb2reaction.domain.bond_changes import has_bond_change
from pdb2reaction.workflows.scan_common import (
    add_scan_common_options,
    collect_staged_scan_values,
    parse_staged_scan_request,
)

logger = logging.getLogger(__name__)


# All defaults imported from defaults.py


def _ensure_stage_dir(base: Path, k: int) -> Path:
    d = base / f"stage_{k:02d}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _echo_scan_summary(stages: List[Dict[str, Any]]) -> None:
    """Print a readable end-of-run summary."""
    emit("\n====== Scan summary ======\n", narrative=True)
    for idx_s, s in enumerate(stages):
        idx = int(s.get("index", 0))
        pairs_1b = list(s.get("pairs_1based", []))
        r0 = list(s.get("initial_distances_A", []))
        rT = list(s.get("target_distances_A", []))
        dA = list(s.get("per_pair_step_A", []))
        N = int(s.get("num_steps", 0))
        bchg = s.get("bond_change", {}) or {}
        changed = bool(bchg.get("changed"))
        summary_txt = (bchg.get("summary") or "").strip()

        # Build the targets triplet string inline
        triples = [f"({i}, {j}, {f'{t:.3f}'.rstrip('0').rstrip('.')})" for (i, j), t in zip(pairs_1b, rT)]
        targets_str = "[" + ", ".join(triples) + "]"

        emit(f"[stage {idx}] Targets (i,j,target Å): {targets_str}", narrative=True)
        emit(f"[stage {idx}] initial distances (Å) = [" + ", ".join(f"'{v:.3f}'" for v in r0) + "]", narrative=True)
        emit(f"[stage {idx}] target distances  (Å) = [" + ", ".join(f"'{v:.3f}'" for v in rT) + "]", narrative=True)
        emit(f"[stage {idx}] per_pair_step     (Å) = [" + ", ".join(f"'{v:.3f}'" for v in dA) + "]", narrative=True)
        emit(f"[stage {idx}] steps N = {N}", narrative=True)
        emit(f"[stage {idx}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}", narrative=True)
        if changed and summary_txt:
            click.echo(textwrap.indent(summary_txt, prefix="  "))
        if not changed:
            click.echo("  (no covalent changes detected)")
        if idx_s != len(stages) - 1:
            click.echo("")  # blank line between stages


def _pair_distances(coords_ang: np.ndarray, pairs: Iterable[Tuple[int, int]]) -> List[float]:
    """
    coords_ang: (N,3) in Å; returns a list of distances (Å) for the given pairs.
    """
    dists: List[float] = []
    for i, j in pairs:
        v = coords_ang[i] - coords_ang[j]
        d = float(np.linalg.norm(v))
        dists.append(d)
    return dists


def _schedule_for_stage(
    coords_ang: np.ndarray,
    tuples: List[Tuple[int, int, float]],
    max_step_size_ang: float,
) -> Tuple[int, List[float], List[float], List[float]]:
    """
    Given current *Å* coords and stage tuples, compute:
      N: number of steps
      r0: initial distances per tuple (Å)
      rT: target distances per tuple (Å)
      step_widths: δ_k per tuple (Å, signed)
    """
    pairs = [(i, j) for (i, j, _) in tuples]
    r0 = _pair_distances(coords_ang, pairs)
    rT = [t for (_, _, t) in tuples]
    deltas = [RT - R0 for (R0, RT) in zip(r0, rT)]
    d_max = max((abs(d) for d in deltas), default=0.0)
    if d_max <= 0.0:
        return 0, r0, rT, [0.0] * len(tuples)
    if max_step_size_ang <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    N = int(math.ceil(d_max / max_step_size_ang))
    step_widths = [d / N for d in deltas]
    return N, r0, rT, step_widths


_COORD_TYPE_DEFAULT = GEOM_KW_DEFAULT["coord_type"]
_snapshot_geometry = make_snapshot_geometry(_COORD_TYPE_DEFAULT)


@click.command(
    help="Bond-length driven scan with staged harmonic restraints and relaxation.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .cif, .mmcif, .xyz, _trj.xyz, ...).",
)
@click.option(
    "-s", "--scan-lists",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help="Scan targets: inline Python literal (e.g. '[(1,5,1.4)]') or a YAML/JSON spec file path. "
         "Atom strings accept positional CHAIN:RESNAME:RESSEQ[ICODE]:ATOM. "
         "Multiple inline literals define sequential stages.",
)
@add_scan_common_options(
    workers_default=UMA_CALC_KW["workers"],
    workers_per_node_default=UMA_CALC_KW["workers_per_node"],
    out_dir_default=OUT_DIR_SCAN,
    baseline_help="(unused)",
    dump_help="Write per-step optimizer trajectory files. scan_trj.xyz is always written; "
              "PDB/CIF companions additionally require conversion topology.",
    max_step_help="Maximum change in any scanned bond length per step [Å].",
    thresh_default=None,
    include_baseline=False,
    include_zmin_zmax=False,
    args_yaml_sections="geom, calc, opt, lbfgs, rfo, bias, bond",
)
@click.option(
    "--print-parsed/--no-print-parsed",
    "print_parsed",
    default=False,
    show_default=True,
    help="Print parsed scan targets after resolving --scan-lists.",
)
@click.option(
    "--endopt/--no-endopt",
    "endopt",
    default=False,
    show_default=True,
    help="After each stage, run an additional unbiased optimization of the stage result.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help=(
        "Resolve and validate options (input, charge/spin parity, "
        "--scan-lists parse) and print the planned scan, then exit "
        "without running any optimization."
    ),
)
@add_coord_type_option()
@add_print_every_option()
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
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    scan_lists_raw: Sequence[str],
    one_based: bool,
    max_step_size: float,
    bias_k: Optional[float],
    relax_max_cycles: int,
    opt_mode: str,
    freeze_links: bool,
    freeze_atoms_text: Optional[str],
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    preopt: bool,
    print_parsed: bool,
    endopt: bool,
    out_json: bool,
    dry_run: bool,
    backend: str,
    solvent: str,
    solvent_model: str,
    cli_coord_type: Optional[str],
    print_every: Optional[int],
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    set_convert_file_enabled(convert_files)
    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, _, _ = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )
    from pdb2reaction.core.utils import resolve_configured_charge_spin
    charge, spin = resolve_configured_charge_spin(
        merged_yaml_cfg, charge=charge, spin=spin, ligand_charge=ligand_charge,
    )

    cycles_overridden = cli_param_overridden(ctx, "relax_max_cycles")
    thresh_overridden = cli_param_overridden(ctx, "thresh")

    with prepared_cli_input(
        input_path,
        ref_pdb=ref_pdb,
        charge=charge,
        spin=spin,
        ligand_charge=ligand_charge,
        prefix="[scan]",
    ) as (prepared_input, resolved_charge, resolved_spin):
        validate_charge_spin_for_prepared(prepared_input, resolved_charge, resolved_spin)
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path
        needs_pdb = source_path.suffix.lower() == ".pdb"
        needs_gjf = prepared_input.is_gjf
        ref_pdb = source_path.resolve() if needs_pdb else None
        time_start = time.perf_counter()
        out_dir_path = Path(out_dir).resolve()
        pdb_atom_meta: List[Dict[str, Any]] = []
        if source_path.suffix.lower() == ".pdb":
            pdb_atom_meta = load_pdb_atom_metadata(source_path)
        cli_scan_values = collect_staged_scan_values(scan_lists_raw, ctx.args)
        scan_request = parse_staged_scan_request(
            cli_scan_values,
            one_based=one_based,
            atom_meta=pdb_atom_meta,
        )
        if dry_run:
            click.echo("[scan] --dry-run: input, charge/spin parity, and --scan-lists parse OK.")
            click.echo(f"[scan] input geometry  : {geom_input_path}")
            click.echo(f"[scan] resolved charge : {int(resolved_charge):+d}")
            click.echo(f"[scan] resolved spin   : {int(resolved_spin)} (multiplicity)")
            click.echo(f"[scan] out_dir         : {out_dir_path}")
            click.echo(
                f"[scan] --scan-lists    : {scan_request.raw_values} "
                f"→ {len(scan_request.stages)} stage(s)"
            )
            click.echo(f"[scan] preopt={bool(preopt)}  endopt={bool(endopt)}  freeze_links={bool(freeze_links)}")
            click.echo("[scan] No scan / preopt was executed.")
            return

        def _run() -> None:

            yaml_cfg = merged_yaml_cfg
            geom_cfg = dict(GEOM_KW_DEFAULT)
            calc_cfg = dict(UMA_CALC_KW)
            opt_cfg  = dict(OPT_BASE_KW)
            lbfgs_cfg = dict(LBFGS_KW)
            rfo_cfg   = dict(RFO_KW)
            bias_cfg  = dict(BIAS_KW)
            bond_cfg  = dict(BOND_KW)

            geom_cfg, calc_cfg, opt_cfg, lbfgs_cfg, rfo_cfg, bias_cfg = build_scan_configs(
                yaml_cfg,
                geom_kw=geom_cfg,
                calc_kw=calc_cfg,
                opt_kw=opt_cfg,
                lbfgs_kw=lbfgs_cfg,
                rfo_kw=rfo_cfg,
                bias_kw=bias_cfg,
                extra_overrides=((bond_cfg, (("bond",),)),),
                charge=resolved_charge,
                spin=resolved_spin,
                workers=workers,
                workers_per_node=workers_per_node,
                out_dir=out_dir,
                thresh=thresh if thresh_overridden else None,
                bias_k=bias_k,
                relax_max_cycles=relax_max_cycles,
                relax_max_cycles_overridden=cycles_overridden,
                workers_overridden=cli_param_overridden(ctx, "workers"),
                workers_per_node_overridden=cli_param_overridden(ctx, "workers_per_node"),
            )

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
            if cli_param_overridden(ctx, "dump"):
                opt_cfg["dump"] = bool(dump)
            if cli_param_overridden(ctx, "print_every") and print_every is not None:
                opt_cfg["print_every"] = int(print_every)
            if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
                geom_cfg["coord_type"] = str(cli_coord_type).lower()
            apply_backend_defaults(calc_cfg)

            kind = normalize_choice(
                opt_mode,
                param="--opt-mode",
                alias_groups=OPT_MODE_ALIASES,
                allowed_hint="grad|hess",
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
            # Resolve freeze list before logging so printed config matches runtime.
            freeze = resolve_freeze_atoms(geom_cfg, source_path, freeze_links)
            calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
            calc_cfg["return_partial_hessian"] = True

            # Present final config
            out_dir_path = Path(opt_cfg["out_dir"]).resolve()
            echo_geom = format_geom_for_echo(geom_cfg)
            echo_calc = format_geom_for_echo(calc_cfg)
            echo_opt  = dict(opt_cfg)
            echo_opt["out_dir"] = str(out_dir_path)
            echo_bias = dict(bias_cfg)
            echo_bond = dict(bond_cfg)
            click.echo(pretty_block("geom", echo_geom))
            click.echo(pretty_block("calc", echo_calc))
            click.echo(pretty_block("opt",  echo_opt))
            max_step_bohr_for_log = float(max_step_size) * ANG2BOHR
            echo_sopt = build_sopt_kwargs(
                kind,
                lbfgs_cfg,
                rfo_cfg,
                opt_cfg,
                max_step_bohr_for_log,
                out_dir=out_dir_path,
                prefix=str(opt_cfg.get("prefix", "")),
            )
            echo_sopt = strip_inherited_keys(echo_sopt, opt_cfg)
            click.echo(pretty_block("lbfgs" if kind == "lbfgs" else "rfo", echo_sopt))
            click.echo(pretty_block("bias", echo_bias))
            click.echo(pretty_block("bond", echo_bond))

            stages = [list(stage) for stage in scan_request.stages]
            scan_one_based = scan_request.one_based
            scan_source = scan_request.source
            _bidir_reset_before = set(scan_request.bidirectional_reset_before)
            _bidir_snapshot_before = set(scan_request.bidirectional_snapshot_before)
            K = len(stages)
            emit(f"[scan] Received {K} stage(s).", narrative=True)
            if print_parsed:
                click.echo(
                    pretty_block(
                        "scan-parsed",
                        {
                            "source": scan_source,
                            "one_based": bool(scan_one_based),
                            "stages_0based": stages,
                        },
                    force=True)
                )

            if pdb_atom_meta:
                emit("[scan] PDB atom details for scanned pairs:", detail=True)
                legend = format_pdb_atom_metadata_header()
                emit(f"        legend: {legend}", detail=True)
                for stage_idx, tuples in enumerate(stages, start=1):
                    emit(f"  Stage {stage_idx}:", detail=True)
                    for pair_idx, (i, j, _) in enumerate(tuples, start=1):
                        emit(
                            f"    pair {pair_idx} i: {format_pdb_atom_metadata(pdb_atom_meta, i)}",
                            detail=True,
                        )
                        emit(
                            f"           j: {format_pdb_atom_metadata(pdb_atom_meta, j)}",
                            detail=True,
                        )

            # Prepare end-of-run summary collector
            stages_summary: List[Dict[str, Any]] = []

            out_dir_path.mkdir(parents=True, exist_ok=True)

            # Load
            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            geom = geom_loader(geom_input_path, coord_type=coord_type, freeze_atoms=freeze)

            max_step_bohr = float(max_step_size) * ANG2BOHR  # shared cap for LBFGS step / RFO trust radii

            # Attach freeze indices to Geometry for optimizer awareness
            if freeze:
                try:
                    geom.freeze_atoms = np.array(freeze, dtype=int)
                except Exception as e:
                    click.echo(
                        f"[scan] WARNING: Failed to attach freeze_atoms to geometry: {e}",
                        err=True,
                    )

            # Build the configured MLIP calculator.
            base_calc = create_calculator(**calc_cfg)
            echo_resolved_device()

            if preopt:
                pre_dir = out_dir_path / "preopt"
                pre_dir.mkdir(parents=True, exist_ok=True)
                geom.set_calculator(base_calc)
                click.echo(f"[preopt] Unbiased relaxation ({kind}) ...")
                optimizer0 = make_sopt_optimizer(
                    geom,
                    kind,
                    lbfgs_cfg,
                    rfo_cfg,
                    opt_cfg,
                    max_step_bohr,
                    out_dir=pre_dir,
                    prefix="preopt",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.")
                except OptimizationError as e:
                    click.echo(f"[preopt] OptimizationError — {e}")

                # Write preopt result
                pre_xyz = pre_dir / "result.xyz"
                with open(pre_xyz, "w") as f:
                    f.write(xyz_string_with_energy(geom))
                click.echo(f"[write] Wrote '{pre_xyz}'.")
                if convert_xyz_like_outputs(
                    pre_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=pre_dir / "result.pdb" if needs_pdb else None,
                    out_gjf_path=pre_dir / "result.gjf" if needs_gjf else None,
                    context="preopt result",
                ):
                    if needs_pdb or needs_gjf:
                        written = []
                        if needs_pdb:
                            written.append("'result.pdb'")
                        if needs_gjf:
                            written.append("'result.gjf'")
                        click.echo(f"[convert] Wrote {', '.join(written)}.")

            # Wrap with bias calculator for the scan
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))
            geom.set_calculator(biased)


            all_trj_blocks: List[str] = []
            _bidir_saved_geom = None
            _bidir_pass1_trj: List[str] = []

            # Iterate stages
            for k, tuples in enumerate(stages, start=1):
                stage_idx_0 = k - 1  # 0-based
                # Bidirectional support: snapshot before pass 1
                if stage_idx_0 in _bidir_snapshot_before:
                    _bidir_saved_geom = _snapshot_geometry(geom)
                    _bidir_pass1_trj = []
                # Bidirectional support: restore geometry before pass 2
                if stage_idx_0 in _bidir_reset_before and _bidir_saved_geom is not None:
                    emit("[bidir] Restoring initial geometry for reverse-direction pass.", narrative=True)
                    geom.coords = _bidir_saved_geom.coords.copy()

                stage_dir = _ensure_stage_dir(out_dir_path, k)
                emit(f"[stage] Stage {k}/{K}", narrative=True)
                tuples_1b = [(i+1, j+1, t) for (i, j, t) in tuples]
                emit(f"Targets (i,j,target Å, 1-based): {tuples_1b}", narrative=True)

                # Snapshot beginning geometry of this stage for bond-change comparison
                start_geom_for_stage = _snapshot_geometry(geom)

                # Current coordinates (Bohr) and schedule computed in Å
                R_bohr = np.array(geom.coords3d, dtype=float)      # (N,3) Bohr
                R_ang  = R_bohr * BOHR2ANG                         # (N,3) Å
                Nsteps, r0, rT, step_widths = _schedule_for_stage(R_ang, tuples, float(max_step_size))
                emit(f"[stage {k}] initial distances (Å) = {['{:.3f}'.format(x) for x in r0]}", narrative=True)
                emit(f"[stage {k}] target distances  (Å) = {['{:.3f}'.format(x) for x in rT]}", narrative=True)
                emit(f"[stage {k}] steps N = {Nsteps}", narrative=True)

                # Record per-stage summary
                srec: Dict[str, Any] = {
                    "index": int(k),
                    "pairs_1based": [(int(i)+1, int(j)+1) for (i, j, _) in tuples],
                    "initial_distances_A": [float(f"{x:.3f}") for x in r0],
                    "target_distances_A": [float(f"{x:.3f}") for x in rT],
                    "per_pair_step_A": [float(f"{x:.3f}") for x in step_widths],
                    "num_steps": int(Nsteps),
                    "converged": None,  # updated after relaxation / endopt
                    # additive terminal status of the last optimizer
                    # (stalled/converged/not_converged) plus its stop_reason.
                    # A stalled leaf keeps ``converged`` False, so aggregation
                    # already treats it as a non-successful required leaf.
                    "optimizer_status": None,
                    "stop_reason": None,
                    # per-step convergence (not just the last optimizer) so a
                    # middle step that fails is not hidden by a converged final step.
                    "step_converged": [],
                    "endopt_converged": None,
                    "endopt_requested": bool(endopt),
                    "bond_change": {"changed": None, "summary": ""},
                }
                stages_summary.append(srec)

                trj_blocks: List[str] = []
                stage_energies: List[Optional[float]] = []
                stage_trj_path = stage_dir / "scan_trj.xyz"
                stage_trj_path.write_text("")  # truncate

                pairs = [(i, j) for (i, j, _) in tuples]

                if Nsteps == 0:
                    # No stepping; optionally perform end-of-stage unbiased optimization
                    if endopt:
                        geom.set_calculator(base_calc)
                        emit(f"[stage {k}] endopt (unbiased) ...", narrative=True)
                        end_optimizer = None
                        try:
                            end_optimizer = make_sopt_optimizer(
                                geom,
                                kind,
                                lbfgs_cfg,
                                rfo_cfg,
                                opt_cfg,
                                max_step_bohr,
                                out_dir=stage_dir,
                                prefix="endopt",
                            )
                            end_optimizer.run()
                        except ZeroStepLength:
                            emit(f"[stage {k}] endopt ZeroStepLength — continuing.", narrative=True)
                        except OptimizationError as e:
                            emit(f"[stage {k}] endopt OptimizationError — {e}", narrative=True)
                        srec["converged"] = getattr(end_optimizer, 'is_converged', None) if end_optimizer is not None else None
                        srec["endopt_converged"] = srec["converged"]
                        if end_optimizer is not None:
                            srec["optimizer_status"] = optimizer_terminal_status(end_optimizer)
                            srec["stop_reason"] = getattr(end_optimizer, "stop_reason", "") or None

                    # No scan steps: empty energy trajectory
                    srec["energies_hartree"] = []

                    # Bond changes: start vs final (possibly endopt)
                    try:
                        changed, summary = has_bond_change(start_geom_for_stage, geom, bond_cfg)
                        emit(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}", narrative=True)
                        if changed and summary and summary.strip():
                            click.echo(textwrap.indent(summary.strip(), prefix="  "))
                        if not changed:
                            click.echo("  (no covalent changes detected)")
                        srec["bond_change"]["changed"] = bool(changed)
                        srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                    except Exception as e:
                        click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                    # Write current (possibly endopted) geometry as the stage result
                    final_energy_h = unbiased_energy_hartree(geom, base_calc)
                    final_xyz = stage_dir / "result.xyz"
                    with open(final_xyz, "w") as f:
                        f.write(
                            xyz_string_with_energy(
                                geom, energy=final_energy_h,
                            )
                        )
                    click.echo(f"[write] Wrote '{final_xyz}'.")
                    srec["final_energy_hartree"] = (
                        float(final_energy_h)
                        if np.isfinite(final_energy_h)
                        else None
                    )
                    if convert_xyz_like_outputs(
                        final_xyz,
                        prepared_input,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=stage_dir / "result.pdb" if needs_pdb else None,
                        out_gjf_path=stage_dir / "result.gjf" if needs_gjf else None,
                        context="stage result",
                    ):
                        if needs_pdb or needs_gjf:
                            written = []
                            if needs_pdb:
                                written.append("'result.pdb'")
                            if needs_gjf:
                                written.append("'result.gjf'")
                            click.echo(f"[convert] Wrote {', '.join(written)}.")
                    continue

                # Run N step(s) with bias
                for s in range(1, Nsteps + 1):
                    # Compute per-pair step target (Å) for this step
                    step_targets = [r0_i + s * dw for (r0_i, dw) in zip(r0, step_widths)]

                    # Update bias well targets (still in Å; wrapper converts internally)
                    biased.set_pairs([(i, j, t) for ((i, j), t) in zip(pairs, step_targets)])
                    # Flushing Geometry caches by re-attaching the calculator
                    geom.set_calculator(biased)

                    # Build optimizer and relax (with bias)
                    prefix = f"scan_s{s:04d}"
                    optimizer = make_sopt_optimizer(
                        geom,
                        kind,
                        lbfgs_cfg,
                        rfo_cfg,
                        opt_cfg,
                        max_step_bohr,
                        out_dir=stage_dir,
                        prefix=prefix,
                    )
                    emit(f"\n[stage {k}] step {s}/{Nsteps}: relaxation ({kind}) ...", narrative=True)
                    try:
                        optimizer.run()
                        srec["step_converged"].append(optimizer_converged_bit(optimizer))
                    except ZeroStepLength:
                        emit(f"[stage {k}] step {s}: ZeroStepLength — continuing to next step.", narrative=True)
                        srec["step_converged"].append(optimizer_converged_bit(optimizer))
                    except OptimizationError as e:
                        emit(f"[stage {k}] step {s}: OptimizationError — {e}", narrative=True)
                        srec["step_converged"].append(False)

                    # Record the UNBIASED (bare PES) energy: geom.energy here
                    # carries the harmonic-bias penalty (the biased calculator is
                    # attached during scan steps), which would inflate the barrier
                    # non-uniformly. scan2d/scan3d already use unbiased_energy_hartree.
                    E_unbiased_h = unbiased_energy_hartree(geom, base_calc)
                    trj_blocks.append(xyz_string_with_energy(geom, energy=E_unbiased_h))
                    stage_energies.append(E_unbiased_h)
                    with open(stage_trj_path, "a") as _tf:
                        _tf.write(trj_blocks[-1])

                # Optional end-of-stage UNBIASED optimization
                if endopt:
                    geom.set_calculator(base_calc)
                    emit(f"[stage {k}] endopt (unbiased) ...", narrative=True)
                    end_optimizer = None
                    try:
                        end_optimizer = make_sopt_optimizer(
                            geom,
                            kind,
                            lbfgs_cfg,
                            rfo_cfg,
                            opt_cfg,
                            max_step_bohr,
                            out_dir=stage_dir,
                            prefix="endopt",
                        )
                        end_optimizer.run()
                    except ZeroStepLength:
                        emit(f"[stage {k}] endopt ZeroStepLength — continuing.", narrative=True)
                    except OptimizationError as e:
                        emit(f"[stage {k}] endopt OptimizationError — {e}", narrative=True)

                # Record convergence of the last optimizer (endopt if used, else last scan step)
                _last_opt = end_optimizer if (endopt and end_optimizer is not None) else optimizer
                srec["converged"] = getattr(_last_opt, 'is_converged', None)
                if _last_opt is not None:
                    srec["optimizer_status"] = optimizer_terminal_status(_last_opt)
                    srec["stop_reason"] = getattr(_last_opt, "stop_reason", "") or None
                if endopt and end_optimizer is not None:
                    srec["endopt_converged"] = getattr(end_optimizer, 'is_converged', None)

                # Store per-step energies in stage record
                srec["energies_hartree"] = stage_energies

                # Bond changes: start vs final (possibly endopt)
                try:
                    changed, summary = has_bond_change(start_geom_for_stage, geom, bond_cfg)
                    emit(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}", narrative=True)
                    if changed and summary and summary.strip():
                        click.echo(textwrap.indent(summary.strip(), prefix="  "))
                    if not changed:
                        click.echo("  (no covalent changes detected)")
                    srec["bond_change"]["changed"] = bool(changed)
                    srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                except Exception as e:
                    click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                # Stage outputs
                if trj_blocks:
                    trj_path = stage_dir / "scan_trj.xyz"
                    click.echo(f"[write] Wrote '{trj_path}'.")
                    # Bidirectional trajectory assembly:
                    # pass 1 (initial→start) is saved; pass 2 (initial→end)
                    # triggers assembly: reversed(pass1) + pass2 → start→initial→end
                    if stage_idx_0 in _bidir_snapshot_before:
                        _bidir_pass1_trj = list(trj_blocks)
                    elif stage_idx_0 in _bidir_reset_before:
                        all_trj_blocks.extend(reversed(_bidir_pass1_trj))
                        all_trj_blocks.extend(trj_blocks)
                        _bidir_pass1_trj = []
                    else:
                        all_trj_blocks.extend(trj_blocks)
                    if convert_xyz_like_outputs(
                        trj_path,
                        prepared_input,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=stage_dir / "scan.pdb" if needs_pdb else None,
                        context="stage trajectory",
                    ):
                        if needs_pdb:
                            click.echo("[convert] Wrote 'scan.pdb'.")

                final_energy_h = unbiased_energy_hartree(geom, base_calc)
                final_xyz = stage_dir / "result.xyz"
                with open(final_xyz, "w") as f:
                    f.write(
                        xyz_string_with_energy(geom, energy=final_energy_h)
                    )
                click.echo(f"[write] Wrote '{final_xyz}'.")
                srec["final_energy_hartree"] = (
                    float(final_energy_h)
                    if np.isfinite(final_energy_h)
                    else None
                )
                if convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=stage_dir / "result.pdb" if needs_pdb else None,
                    out_gjf_path=stage_dir / "result.gjf" if needs_gjf else None,
                    context="stage result",
                ):
                    if needs_pdb or needs_gjf:
                        written = []
                        if needs_pdb:
                            written.append("'result.pdb'")
                        if needs_gjf:
                            written.append("'result.gjf'")
                        click.echo(f"[convert] Wrote {', '.join(written)}.")

            # 4b) Write combined scan_trj.xyz + scan.pdb to out_dir
            if all_trj_blocks:
                combined_trj = out_dir_path / "scan_trj.xyz"
                with open(combined_trj, "w") as f:
                    f.write("".join(all_trj_blocks))
                click.echo(f"[write] Wrote '{combined_trj}'.")
                if convert_xyz_like_outputs(
                    combined_trj,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=out_dir_path / "scan.pdb" if needs_pdb else None,
                    context="combined scan trajectory",
                ):
                    if needs_pdb:
                        click.echo(f"[convert] Wrote '{out_dir_path / 'scan.pdb'}'.")

            _echo_scan_summary(stages_summary)

            emit("\n====== Scan finished ======\n", narrative=True)
            emit(format_elapsed("[time] Elapsed Time for Scan", time_start), narrative=True)

            # result.json (if --out-json)
            if out_json:
                from pdb2reaction.core.utils import calculator_provenance, write_result_json
                from pdb2reaction.workflows._outcomes import (
                    aggregate_workflow_truth,
                    attach_outcomes,
                    combine_step_convergence as _combine_step_convergence,
                    make_leaf,
                )

                json_stages = []
                _stage_leaves = []
                for srec in stages_summary:
                    stage_entry: Dict[str, Any] = {
                        "index": srec["index"],
                        "n_steps": srec["num_steps"],
                        "converged": srec.get("converged"),
                        "bond_changes": srec.get("bond_change", {}),
                        "pairs_1based": srec.get("pairs_1based"),
                        "initial_distances_angstrom": srec.get("initial_distances_A"),
                        "target_distances_angstrom": srec.get("target_distances_A"),
                    }
                    # Final energy: use value captured directly from geometry object
                    stage_entry["final_energy_hartree"] = srec.get("final_energy_hartree")
                    # Per-step energy trajectory
                    stage_entry["energies_hartree"] = srec.get("energies_hartree", [])
                    # Surface the optimizer terminal status and stall
                    # reason so a stalled scan stage is not silently dropped from
                    # result.json (additive; absent for stages that never set it).
                    if srec.get("optimizer_status"):
                        stage_entry["optimizer_status"] = srec["optimizer_status"]
                    if srec.get("stop_reason"):
                        stage_entry["stop_reason"] = srec["stop_reason"]
                    json_stages.append(stage_entry)

                    # the stage leaf is usable only when EVERY step converged
                    # and (when requested) the endopt converged. A converged final
                    # step must not hide an earlier failed step.
                    _steps = list(srec.get("step_converged") or [])
                    _stage_conv = _combine_step_convergence(_steps)
                    if srec.get("endopt_requested"):
                        _eo = srec.get("endopt_converged")
                        _stage_conv = _combine_step_convergence(
                            _steps + [_eo if isinstance(_eo, bool) else None]
                        )
                    _fe = srec.get("final_energy_hartree")
                    _energy_valid = (
                        isinstance(_fe, (int, float, np.integer, np.floating))
                        and np.isfinite(float(_fe))
                    )
                    _stage_leaves.append(
                        make_leaf(
                            "scan",
                            f"stage_{srec['index']}",
                            executed=True,
                            converged=_stage_conv,
                            energy_valid=_energy_valid,
                        )
                    )
                _truth = aggregate_workflow_truth(
                    _stage_leaves,
                    [f"stage_{srec['index']}" for srec in stages_summary],
                )
                result_data: Dict[str, Any] = {
                    "status": "completed",
                    "charge": resolved_charge,
                    "spin": resolved_spin,
                    "backend": calc_cfg.get("backend", backend),
                    "model": calc_cfg.get("model"),
                    **calculator_provenance(calc_cfg),
                    "solvent": calc_cfg.get("solvent", "none"),
                    "preopt": bool(preopt),
                    "max_step_size_angstrom": float(max_step_size),
                    "n_stages": len(stages_summary),
                    "stages": json_stages,
                    "files": {
                        "scan_trj_xyz": "scan_trj.xyz",
                    },
                }
                for ext in (".pdb", ".cif"):
                    f = out_dir_path / f"scan{ext}"
                    if f.exists():
                        result_data["files"][f"scan_{ext[1:]}"] = f.name
                # Additive outcome fields; legacy ``status`` stays "completed".
                attach_outcomes(result_data, truth=_truth, stage_outcomes=_stage_leaves)
                write_result_json(
                    out_dir_path, result_data,
                    command="scan",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

        run_cli(_run, label="scan", out_dir=out_dir_path, command="scan", time_start=time_start)
