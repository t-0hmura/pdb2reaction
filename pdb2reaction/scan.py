# pdb2reaction/scan.py

"""
Staged bond-length scan with harmonic restraints and full relaxation using UMA calculator.

Example:
    pdb2reaction scan -i input.pdb -q 0 --scan-lists '[(12,45,1.35)]' --preopt True --endopt True

For detailed documentation, see: docs/scan.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import math
import sys
import textwrap

import click
import numpy as np
import yaml
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR

from .defaults import (
    GEOM_KW_DEFAULT,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    BIAS_KW,
    BOND_KW,
    OPT_MODE_ALIASES,
    UMA_CALC_KW,
    OUT_DIR_SCAN,
)
from .uma_pysis import uma_pysis
from .opt import HarmonicBiasCalculator
from .utils import (
    build_sopt_kwargs,
    make_sopt_optimizer,
    load_yaml_dict,
    build_scan_configs,
    collect_single_option_values,
    cli_param_overridden,
    pretty_block,
    strip_inherited_keys,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    parse_scan_list_triples,
    prepared_cli_input,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    load_pdb_atom_metadata,
    format_pdb_atom_metadata,
    format_pdb_atom_metadata_header,
    make_snapshot_geometry,
    resolve_freeze_atoms,
    xyz_string_with_energy,
)
from .cli_utils import run_cli
from .bond_changes import has_bond_change
from .scan_common import add_scan_common_options


# All defaults imported from defaults.py


def _ensure_stage_dir(base: Path, k: int) -> Path:
    d = base / f"stage_{k:02d}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _echo_scan_summary(stages: List[Dict[str, Any]]) -> None:
    """Print a readable end-of-run summary."""
    click.echo("\nSummary")
    click.echo("------------------")
    for s in stages:
        idx = int(s.get("index", 0))
        pairs_1b = list(s.get("pairs_1based", []))
        r0 = list(s.get("initial_distances_A", []))
        rT = list(s.get("target_distances_A", []))
        dA = list(s.get("per_pair_step_A", []))
        N = int(s.get("num_steps", 0))
        bchg = s.get("bond_change", {}) or {}
        changed = bool(bchg.get("changed"))
        summary_txt = (bchg.get("summary") or "").strip()

        # Inline _targets_triplet_str and _fmt_target_value
        triples = [f"({i}, {j}, {f'{t:.3f}'.rstrip('0').rstrip('.')})" for (i, j), t in zip(pairs_1b, rT)]
        targets_str = "[" + ", ".join(triples) + "]"

        click.echo(f"[stage {idx}] Targets (i,j,target Å): {targets_str}")
        click.echo(f"[stage {idx}] initial distances (Å) = [" + ", ".join(f"'{v:.3f}'" for v in r0) + "]")
        click.echo(f"[stage {idx}] target distances  (Å) = [" + ", ".join(f"'{v:.3f}'" for v in rT) + "]")
        click.echo(f"[stage {idx}] per_pair_step     (Å) = [" + ", ".join(f"'{v:.3f}'" for v in dA) + "]")
        click.echo(f"[stage {idx}] steps N = {N}")
        click.echo(f"[stage {idx}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
        if changed and summary_txt:
            click.echo(textwrap.indent(summary_txt, prefix="  "))
        if not changed:
            click.echo("  (no covalent changes detected)")
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
    context_settings={"help_option_names": ["-h", "--help"], "allow_extra_args": True},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, ...).",
)
@click.option(
    "--scan-lists",
    "--scan-list",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help="Python-like list of (i,j,target) per stage. Pass a single --scan-list(s) followed by "
         "multiple literals to run sequential stages, e.g. --scan-lists '[(0,1,1.50)]' '[(5,7,1.20)]'.",
)
@add_scan_common_options(
    workers_default=UMA_CALC_KW["workers"],
    workers_per_node_default=UMA_CALC_KW["workers_per_node"],
    out_dir_default="./result_scan/",
    baseline_help="(unused)",
    dump_help="Write stage trajectory as scan.trj (and scan.pdb for PDB input).",
    max_step_help="Maximum change in any scanned bond length per step [Å].",
    thresh_default=None,
    include_baseline=False,
    include_zmin_zmax=False,
)
@click.option("--endopt", type=click.BOOL, default=True, show_default=True,
              help="After each stage, run an additional unbiased optimization of the stage result.")
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
    dump: bool,
    convert_files: bool,
    ref_pdb: Optional[Path],
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    endopt: bool,
) -> None:
    set_convert_file_enabled(convert_files)

    cycles_overridden = cli_param_overridden(ctx, "relax_max_cycles")

    with prepared_cli_input(
        input_path,
        ref_pdb=ref_pdb,
        charge=charge,
        spin=spin,
        ligand_charge=ligand_charge,
        prefix="[scan]",
    ) as (prepared_input, resolved_charge, resolved_spin):
        geom_input_path = prepared_input.geom_path
        source_path = prepared_input.source_path
        needs_pdb = source_path.suffix.lower() == ".pdb"
        needs_gjf = prepared_input.is_gjf
        ref_pdb = source_path.resolve() if needs_pdb else None
        def _run() -> None:
            time_start = time.perf_counter()

            # ------------------------------------------------------------------
            # 1) Assemble configuration (defaults ← CLI ← YAML) - create fresh copies for merging
            # ------------------------------------------------------------------
            yaml_cfg = load_yaml_dict(args_yaml)

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
                thresh=thresh,
                bias_k=bias_k,
            )

            kind = normalize_choice(
                opt_mode,
                param="--opt-mode",
                alias_groups=OPT_MODE_ALIASES,
                allowed_hint="light|heavy",
            )

            # Resolve freeze list before logging so printed config matches runtime.
            freeze = resolve_freeze_atoms(geom_cfg, source_path, freeze_links)

            # Present final config
            out_dir_path = Path(opt_cfg["out_dir"]).resolve()
            echo_geom = format_geom_for_echo(geom_cfg)
            echo_calc = format_geom_for_echo(calc_cfg)
            echo_opt  = dict(opt_cfg)
            if cycles_overridden:
                echo_opt["max_cycles"] = int(relax_max_cycles)
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
                relax_max_cycles,
                cycles_overridden,
                out_dir_path,
                str(opt_cfg.get("prefix", "")),
            )
            echo_sopt = strip_inherited_keys(echo_sopt, opt_cfg)
            click.echo(pretty_block("lbfgs" if kind == "lbfgs" else "rfo", echo_sopt))
            click.echo(pretty_block("bias", echo_bias))
            click.echo(pretty_block("bond", echo_bond))

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            # ------------------------------------------------------------------
            # 2) Parse scan lists
            # ------------------------------------------------------------------
            scan_lists_raw = collect_single_option_values(
                sys.argv[1:], ("--scan-lists", "--scan-list"), "--scan-list/--scan-lists"
            )
            if not scan_lists_raw:
                raise click.BadParameter("--scan-list(s) must be provided at least once.")
            stages: List[List[Tuple[int, int, float]]] = []
            for idx, raw in enumerate(scan_lists_raw, start=1):
                parsed, _ = parse_scan_list_triples(
                    raw,
                    one_based=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name=f"--scan-lists #{idx}",
                )
                for i, j, r in parsed:
                    if r <= 0.0:
                        raise click.BadParameter(
                            f"Non-positive target length in --scan-lists #{idx}: {(i, j, r)}."
                        )
                stages.append(parsed)
            K = len(stages)
            click.echo(f"[scan] Received {K} stage(s).")

            if pdb_atom_meta:
                click.echo("[scan] PDB atom details for scanned pairs:")
                legend = format_pdb_atom_metadata_header()
                click.echo(f"        legend: {legend}")
                for stage_idx, tuples in enumerate(stages, start=1):
                    click.echo(f"  Stage {stage_idx}:")
                    for pair_idx, (i, j, _) in enumerate(tuples, start=1):
                        click.echo(
                            f"    pair {pair_idx} i: {format_pdb_atom_metadata(pdb_atom_meta, i)}"
                        )
                        click.echo(
                            f"           j: {format_pdb_atom_metadata(pdb_atom_meta, j)}"
                        )

            # Prepare end-of-run summary collector
            stages_summary: List[Dict[str, Any]] = []

            # ------------------------------------------------------------------
            # 3) Load geometry (Cartesian) and set calculator (UMA → harmonic-bias wrapper)
            # ------------------------------------------------------------------
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

            # Build UMA calculator (only uma_pysis is supported)
            base_calc = uma_pysis(**calc_cfg)

            # ------------------------------------------------------------------
            # Optional preoptimization WITHOUT bias
            # ------------------------------------------------------------------
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
                    relax_max_cycles,
                    cycles_overridden,
                    pre_dir,
                    "preopt_",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo(f"[preopt] ZeroStepLength — continuing.", err=True)
                except OptimizationError as e:
                    click.echo(f"[preopt] OptimizationError — {e}", err=True)

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

            # ------------------------------------------------------------------
            # 4) Stage-by-stage scan
            # ------------------------------------------------------------------

            # Iterate stages
            for k, tuples in enumerate(stages, start=1):
                stage_dir = _ensure_stage_dir(out_dir_path, k)
                click.echo(f"\n--- Stage {k}/{K} ---")
                click.echo(f"Targets (i,j,target Å): {tuples}")

                # Snapshot beginning geometry of this stage for bond-change comparison
                start_geom_for_stage = _snapshot_geometry(geom)

                # Current coordinates (Bohr) and schedule computed in Å
                R_bohr = np.array(geom.coords3d, dtype=float)      # (N,3) Bohr
                R_ang  = R_bohr * BOHR2ANG                         # (N,3) Å
                Nsteps, r0, rT, step_widths = _schedule_for_stage(R_ang, tuples, float(max_step_size))
                click.echo(f"[stage {k}] initial distances (Å) = {['{:.3f}'.format(x) for x in r0]}")
                click.echo(f"[stage {k}] target distances  (Å) = {['{:.3f}'.format(x) for x in rT]}")
                click.echo(f"[stage {k}] steps N = {Nsteps}")

                # Record per-stage summary
                srec: Dict[str, Any] = {
                    "index": int(k),
                    "pairs_1based": [(int(i)+1, int(j)+1) for (i, j, _) in tuples],
                    "initial_distances_A": [float(f"{x:.3f}") for x in r0],
                    "target_distances_A": [float(f"{x:.3f}") for x in rT],
                    "per_pair_step_A": [float(f"{x:.3f}") for x in step_widths],
                    "num_steps": int(Nsteps),
                    "bond_change": {"changed": None, "summary": ""},
                }
                stages_summary.append(srec)

                trj_blocks: List[str] = [] if dump else None

                pairs = [(i, j) for (i, j, _) in tuples]

                if Nsteps == 0:
                    # No stepping; optionally perform end-of-stage unbiased optimization
                    if endopt:
                        geom.set_calculator(base_calc)
                        click.echo(f"[stage {k}] endopt (unbiased) ...")
                        try:
                            end_optimizer = make_sopt_optimizer(
                                geom,
                                kind,
                                lbfgs_cfg,
                                rfo_cfg,
                                opt_cfg,
                                max_step_bohr,
                                relax_max_cycles,
                                cycles_overridden,
                                stage_dir,
                                "endopt_",
                            )
                            end_optimizer.run()
                        except ZeroStepLength:
                            click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                        except OptimizationError as e:
                            click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)

                    # Bond changes: start vs final (possibly endopt)
                    try:
                        changed, summary = has_bond_change(start_geom_for_stage, geom, bond_cfg)
                        click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                        if changed and summary and summary.strip():
                            click.echo(textwrap.indent(summary.strip(), prefix="  "))
                        if not changed:
                            click.echo("  (no covalent changes detected)")
                        srec["bond_change"]["changed"] = bool(changed)
                        srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                    except Exception as e:
                        click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                    # Write current (possibly endopted) geometry as the stage result
                    final_xyz = stage_dir / "result.xyz"
                    with open(final_xyz, "w") as f:
                        f.write(xyz_string_with_energy(geom))
                    click.echo(f"[write] Wrote '{final_xyz}'.")
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
                    prefix = f"scan_s{s:04d}_"
                    optimizer = make_sopt_optimizer(
                        geom,
                        kind,
                        lbfgs_cfg,
                        rfo_cfg,
                        opt_cfg,
                        max_step_bohr,
                        relax_max_cycles,
                        cycles_overridden,
                        stage_dir,
                        prefix,
                    )
                    click.echo(f"[stage {k}] step {s}/{Nsteps}: relaxation ({kind}) ...")
                    try:
                        optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] step {s}: ZeroStepLength — continuing to next step.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] step {s}: OptimizationError — {e}", err=True)

                    # Record trajectory block only when requested (biased result)
                    if dump and trj_blocks is not None:
                        trj_blocks.append(xyz_string_with_energy(geom))

                # Optional end-of-stage UNBIASED optimization
                if endopt:
                    geom.set_calculator(base_calc)
                    click.echo(f"[stage {k}] endopt (unbiased) ...")
                    try:
                        end_optimizer = make_sopt_optimizer(
                            geom,
                            kind,
                            lbfgs_cfg,
                            rfo_cfg,
                            opt_cfg,
                            max_step_bohr,
                            relax_max_cycles,
                            cycles_overridden,
                            stage_dir,
                            "endopt_",
                        )
                        end_optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)

                # Bond changes: start vs final (possibly endopt)
                try:
                    changed, summary = has_bond_change(start_geom_for_stage, geom, bond_cfg)
                    click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                    if changed and summary and summary.strip():
                        click.echo(textwrap.indent(summary.strip(), prefix="  "))
                    if not changed:
                        click.echo("  (no covalent changes detected)")
                    srec["bond_change"]["changed"] = bool(changed)
                    srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                except Exception as e:
                    click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                # Stage outputs
                if dump and trj_blocks:
                    trj_path = stage_dir / "scan.trj"
                    with open(trj_path, "w") as f:
                        f.write("".join(trj_blocks))
                    click.echo(f"[write] Wrote '{trj_path}'.")
                    if convert_xyz_like_outputs(
                        trj_path,
                        prepared_input,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=stage_dir / "scan.pdb" if needs_pdb else None,
                        out_gjf_path=stage_dir / "scan.gjf" if needs_gjf else None,
                        context="stage trajectory",
                    ):
                        if needs_pdb or needs_gjf:
                            written = []
                            if needs_pdb:
                                written.append("'scan.pdb'")
                            if needs_gjf:
                                written.append("'scan.gjf'")
                            click.echo(f"[convert] Wrote {', '.join(written)}.")

                final_xyz = stage_dir / "result.xyz"
                with open(final_xyz, "w") as f:
                    f.write(xyz_string_with_energy(geom))
                click.echo(f"[write] Wrote '{final_xyz}'.")
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

            # ------------------------------------------------------------------
            # 5) Final summary echo (human‑friendly)
            # ------------------------------------------------------------------
            _echo_scan_summary(stages_summary)
            # ------------------------------------------------------------------

            click.echo("\n=== Scan finished ===\n")

            click.echo(format_elapsed("[time] Elapsed Time for Scan", time_start))

        run_cli(_run, label="scan")
