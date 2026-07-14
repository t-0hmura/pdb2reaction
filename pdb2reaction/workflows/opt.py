"""
Single-structure geometry optimization using LBFGS or RFO with an MLIP calculator.

Example:
    pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess

For detailed documentation, see: docs/opt.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import logging
import sys

import click
import numpy as np
import torch
import time
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, BOHR2ANG, AU2EV
from pysisyphus.tr_projection import normalize_tr_projection_mode
from pdb2reaction.workflows.restraints import HarmonicBiasCalculator

from pdb2reaction.core.defaults import (
    BIAS_KW,
    GEOM_KW_DEFAULT,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    OPT_MODE_ALIASES,
    THRESH_CHOICES,
    UMA_CALC_KW,
    OUT_DIR_OPT,
    apply_backend_defaults,
)
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.utils import (
    resolve_freeze_atoms,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    cli_param_overridden,
    is_scan_spec_file,
    parse_dist_freeze_list,
    parse_dist_freeze_spec,
    load_pdb_atom_metadata,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
    emit_optimizer_terminal_status,
    optimizer_cycle_count,
)
from pdb2reaction.cli.common_options import (
    add_ml_charge_spin_options,
    add_coord_type_option,
    add_print_every_option,
    add_precision_option, add_backend_model_option,
    add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from pdb2reaction.cli.decorators import run_cli, resolve_yaml_sources, load_merged_yaml_cfg
from pdb2reaction.workflows.freq import (
    _torch_device,
    _safe_masses_amu,
    _calc_full_hessian_torch,
    _frequencies_cm_and_modes,
    _calc_energy,
)

logger = logging.getLogger(__name__)

EV2AU = 1.0 / AU2EV  # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)
OPT_FLATTEN_NEG_FREQ_THRESH_CM = 5.0
OPT_FLATTEN_AMP_ANG = 0.10
OPT_FLATTEN_MAX_ITER = 50


def _parse_dist_freeze_args(
    raw_args: Sequence[str],
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse all ``--dist-freeze`` arguments (inline literal or spec file).

    Accepts the same format as ``--scan-lists``: inline Python literal
    (e.g. ``'[(1,5,1.4)]'``) or a YAML/JSON spec file path.  String atom
    specs (e.g. ``'A:SER123:OG'``) are supported when *atom_meta* is
    available.  Target distance is optional — omit to freeze at the current
    distance.
    """
    all_pairs: List[Tuple[int, int, Optional[float]]] = []
    for raw in raw_args:
        if is_scan_spec_file(raw):
            all_pairs.extend(parse_dist_freeze_spec(
                Path(raw),
                one_based_default=one_based,
                atom_meta=atom_meta,
            ))
        else:
            all_pairs.extend(parse_dist_freeze_list(
                raw,
                one_based=one_based,
                atom_meta=atom_meta,
            ))
    return all_pairs


def _resolve_dist_freeze_targets(
    geometry,
    tuples: List[Tuple[int, int, Optional[float]]],
) -> List[Tuple[int, int, float]]:
    coords_bohr = np.array(geometry.coords3d, dtype=float).reshape(-1, 3)
    coords_ang = coords_bohr * BOHR2ANG
    n = coords_ang.shape[0]
    resolved: List[Tuple[int, int, float]] = []
    for (i, j, target) in tuples:
        if not (0 <= i < n and 0 <= j < n):
            raise click.BadParameter(
                f"--dist-freeze indices {(i, j)} are out of bounds for the loaded geometry (N={n})."
            )
        if target is None:
            vec = coords_ang[i] - coords_ang[j]
            dist = float(np.linalg.norm(vec))
        else:
            dist = float(target)
        resolved.append((i, j, dist))
    return resolved


def _convert_outputs(
    prepared_input: "PreparedInputStructure",
    out_dir: Path,
    dump: bool,
    get_trj_fn,
    final_xyz_path: Path,
) -> None:
    """Convert outputs to PDB/CIF/GJF companions requested by the input type."""
    needs_pdb = prepared_input.source_path.suffix.lower() == ".pdb"
    needs_gjf = prepared_input.is_gjf
    if not (needs_pdb or needs_gjf):
        return

    ref_pdb = prepared_input.source_path.resolve() if needs_pdb else None

    # final_geometry.xyz → final_geometry.{pdb|gjf}
    if convert_xyz_like_outputs(
        final_xyz_path,
        prepared_input,
        ref_pdb_path=ref_pdb,
        out_pdb_path=out_dir / "final_geometry.pdb" if needs_pdb else None,
        out_gjf_path=out_dir / "final_geometry.gjf" if needs_gjf else None,
        context="final geometry",
    ):
        click.echo("[convert] Wrote 'final_geometry' outputs.")

    # optimization_trj.xyz → optimization.pdb (if dump)
    if dump and needs_pdb:
        trj_path = get_trj_fn("optimization_trj.xyz")
        if trj_path.exists():
            if convert_xyz_like_outputs(
                trj_path,
                prepared_input,
                ref_pdb_path=ref_pdb,
                out_pdb_path=out_dir / "optimization.pdb" if needs_pdb else None,
                context="optimization trajectory",
            ):
                click.echo("[convert] Wrote 'optimization' outputs.")
        else:
            click.echo("[convert] WARNING: 'optimization_trj.xyz' not found; skipping conversion.", err=True)


def _flatten_all_imag_modes_for_geom(
    geom,
    masses_amu: np.ndarray,
    calc_kwargs: dict,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    flatten_amp_ang: float,
) -> bool:
    """
    Flatten all imaginary modes for a geometry in a single pass.
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) == 0:
        return False

    order = np.argsort(freqs_cm[neg_idx_all])  # most negative first
    targets = [int(x) for x in neg_idx_all[order]]
    mass_scale = np.sqrt(12.011 / masses_amu)[:, None]
    amp_bohr = float(flatten_amp_ang) / BOHR2ANG
    E_ref = _calc_energy(geom, calc_kwargs)

    m3 = np.repeat(masses_amu, 3).reshape(-1, 3)
    for idx in targets:
        v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
        v_cart = v_mw / np.sqrt(m3)
        v_cart /= np.linalg.norm(v_cart)

        disp = amp_bohr * mass_scale * v_cart
        ref = geom.cart_coords.reshape(-1, 3)

        plus = ref + disp
        minus = ref - disp

        geom.coords = plus.reshape(-1)
        E_plus = _calc_energy(geom, calc_kwargs)

        geom.coords = minus.reshape(-1)
        E_minus = _calc_energy(geom, calc_kwargs)

        use_plus = E_plus <= E_minus
        geom.coords = (plus if use_plus else minus).reshape(-1)
        E_keep = E_plus if use_plus else E_minus
        delta_e = E_keep - E_ref
        click.echo(
            f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
            f"E_disp={E_keep:.8f} Ha ΔE={delta_e:+.8f} Ha"
        )

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True



@click.command(
    help="Single-structure geometry optimization using LBFGS or RFO.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .cif, .mmcif, .xyz, .gjf, _trj.xyz, ...).",
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
    "--dist-freeze",
    "dist_freeze_raw",
    type=str,
    multiple=True,
    default=(),
    show_default=False,
    help="Distance restraints: inline Python literal (e.g. '[(1,5,1.4)]') or a YAML/JSON spec file path. "
         "Same format as --scan-lists: (i,j,target_A) triples. "
         "Target may be omitted to freeze at the current distance: (i,j).",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret --dist-freeze / --scan-lists indices as 1-based (default) or 0-based.",
)
@click.option(
    "--bias-k",
    type=float,
    default=BIAS_KW["k"],
    show_default=True,
    help="Harmonic restraint strength k [eV/Å^2] for --dist-freeze.",
)
@click.option(
    "--freeze-links/--no-freeze-links",
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
    default=GEOM_KW_DEFAULT["tr_projection"],
    show_default=True,
    help=(
        "Rigid translation/rotation treatment used by --flatten PHVA. "
        "'constrained' respects frozen anchors; 'legacy-active' treats the "
        "active fragment as isolated."
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
    "--max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum number of optimization cycles.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess", "lbfgs", "rfo"], case_sensitive=False),
    default="grad",
    show_default=True,
    help="Optimization mode: grad (lbfgs) or hess (rfo). Aliases lbfgs/rfo are accepted.",
)
@click.option(
    "--flatten/--no-flatten",
    "flatten",
    default=False,
    show_default=True,
    help="Enable/disable imaginary-mode flatten loop after optimization.",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write optimization trajectory to 'optimization_trj.xyz'.",
)
@click.option(
    "-o", "--out-dir",
    type=str,
    default=OUT_DIR_OPT,
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default=False,
    help=(
        "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
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
    help="Validate options and print the execution plan without running optimization.",
)
@click.option("-b", "--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              show_default=True, help="MLIP backend.")
@click.option("--solvent", default="none", show_default=True,
              help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              show_default=True, help="xTB solvent model.")
@add_coord_type_option()
@add_ml_charge_spin_options()
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
    dist_freeze_raw: Sequence[str],
    one_based: bool,
    bias_k: float,
    freeze_links: bool,
    freeze_atoms_text: Optional[str],
    tr_projection: str,
    convert_files: bool,
    ref_pdb: Optional[Path],
    max_cycles: int,
    opt_mode: str,
    flatten: bool,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    out_json: bool,
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
    time_start = time.perf_counter()

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

    def _run() -> None:
        with prepared_cli_input(
            input_path,
            ref_pdb=ref_pdb,
            charge=charge,
            spin=spin,
            ligand_charge=ligand_charge,
            prefix="[opt]",
        ) as (prepared_input, resolved_charge, resolved_spin):
            validate_charge_spin_for_prepared(prepared_input, resolved_charge, resolved_spin)
            geom_input_path = prepared_input.geom_path
            source_path = prepared_input.source_path

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            dist_freeze = _parse_dist_freeze_args(
                dist_freeze_raw, one_based=bool(one_based), atom_meta=pdb_atom_meta,
            )

            geom_cfg = dict(GEOM_KW_DEFAULT)
            calc_cfg = dict(UMA_CALC_KW)
            opt_cfg = dict(OPT_BASE_KW)
            lbfgs_cfg = dict(LBFGS_KW)
            rfo_cfg = dict(RFO_KW)

            apply_yaml_overrides(
                config_layer_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",),)),
                    (opt_cfg, (("opt",),)),
                    (lbfgs_cfg, (("lbfgs",),)),
                    (rfo_cfg, (("rfo",),)),
                ],
            )

            # CLI explicit overrides
            calc_cfg["charge"] = resolved_charge
            calc_cfg["spin"] = resolved_spin
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
            if cli_param_overridden(ctx, "max_cycles"):
                opt_cfg["max_cycles"] = int(max_cycles)
            if cli_param_overridden(ctx, "dump"):
                opt_cfg["dump"] = bool(dump)
            if cli_param_overridden(ctx, "out_dir"):
                opt_cfg["out_dir"] = out_dir
            if cli_param_overridden(ctx, "thresh") and thresh is not None:
                opt_cfg["thresh"] = str(thresh)
            if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
                geom_cfg["coord_type"] = str(cli_coord_type).lower()
            if cli_param_overridden(ctx, "tr_projection"):
                geom_cfg["tr_projection"] = str(tr_projection).lower()
            if cli_param_overridden(ctx, "print_every") and print_every is not None:
                opt_cfg["print_every"] = int(print_every)

            apply_yaml_overrides(
                override_layer_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",),)),
                    (opt_cfg, (("opt",),)),
                    (lbfgs_cfg, (("lbfgs",),)),
                    (rfo_cfg, (("rfo",),)),
                ],
            )
            try:
                geom_cfg["tr_projection"] = normalize_tr_projection_mode(
                    geom_cfg.get("tr_projection")
                )
            except ValueError as exc:
                raise click.ClickException(str(exc)) from exc

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
            resolve_freeze_atoms(geom_cfg, source_path, freeze_links)
            calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
            calc_cfg["return_partial_hessian"] = True

            # Normalize and select optimizer kind
            kind = normalize_choice(
                opt_mode,
                param="--opt-mode",
                alias_groups=OPT_MODE_ALIASES,
                allowed_hint="grad|hess|lbfgs|rfo",
            )
            main_kind = kind
            flatten_kind = kind

            # Pretty-print the resolved configuration
            out_dir_path = Path(opt_cfg["out_dir"]).resolve()

            # Default-verbosity entry summary (skipped in child mode).
            from pdb2reaction.core.utils import echo_run_summary
            _model = calc_cfg.get("model")
            echo_run_summary({
                "input": str(input_path),
                "backend": (
                    f"{calc_cfg.get('backend', '?')} ({_model}, {calc_cfg.get('precision', 'fp32')})"
                    if _model else calc_cfg.get("backend", "?")
                ),
                "opt": f"{kind}, max_cycles={opt_cfg.get('max_cycles', '?')}",
                "out": str(out_dir_path),
            })

            click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
            click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
            click.echo(pretty_block("opt", {**opt_cfg, "out_dir": str(out_dir_path)}))
            echo_sopt = dict(lbfgs_cfg if kind == "lbfgs" else rfo_cfg)
            echo_sopt = strip_inherited_keys(echo_sopt, opt_cfg)
            click.echo(pretty_block(kind, echo_sopt))
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
            if dist_freeze:
                display_pairs = []
                for (i, j, target) in dist_freeze:
                    label = (f"{target:.4f}" if target is not None else "<current>")
                    display_pairs.append((int(i) + 1, int(j) + 1, label))
                click.echo(
                    pretty_block(
                        "dist_freeze (input)",
                        {
                            "k (eV/Å^2)": float(bias_k),
                            "pairs_1based": display_pairs,
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
                            "optimizer_mode": str(kind),
                            "flatten_optimizer_mode": str(flatten_kind) if bool(flatten) else None,
                            "flatten": bool(flatten),
                            "convert_files": bool(convert_files),
                            "freeze_links": bool(freeze_links),
                            "tr_projection": geom_cfg["tr_projection"],
                            "will_run_optimization": True,
                            "will_write_final_geometry": True,
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. Optimization execution was skipped.")
                return

            out_dir_path.mkdir(parents=True, exist_ok=True)

            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            # Pass all geometry kwargs except coord_type as coord_kwargs
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)
            geometry = geom_loader(
                geom_input_path,
                coord_type=coord_type,
                **coord_kwargs,
            )

            # Attach the configured MLIP calculator.
            base_calc = create_calculator(**calc_cfg)
            geometry.set_calculator(base_calc)
            echo_resolved_device()

            resolved_dist_freeze: List[Tuple[int, int, float]] = []
            if dist_freeze:
                try:
                    resolved_dist_freeze = _resolve_dist_freeze_targets(geometry, dist_freeze)
                except click.BadParameter as e:
                    click.echo(f"ERROR: {e}", err=True)
                    sys.exit(1)
                click.echo(
                    pretty_block(
                        "dist_freeze (active)",
                        {
                            "k (eV/Å^2)": float(bias_k),
                            "pairs_1based": [
                                (int(i) + 1, int(j) + 1, float(f"{t:.4f}"))
                                for (i, j, t) in resolved_dist_freeze
                            ],
                        },
                    )
                )
                bias_calc = HarmonicBiasCalculator(base_calc, k=float(bias_k))
                bias_calc.set_pairs(resolved_dist_freeze)
                geometry.set_calculator(bias_calc)

            opt_calc = geometry.calculator
            common_kwargs = dict(opt_cfg)
            # Ensure paths (strings) are OK; Optimizer expects str, not Path
            common_kwargs["out_dir"] = str(out_dir_path)

            def _build_optimizer(run_kind: str):
                if run_kind == "lbfgs":
                    lbfgs_args = {**lbfgs_cfg, **common_kwargs}
                    return LBFGS(geometry, **lbfgs_args)
                if run_kind == "rfo":
                    rfo_args = {**rfo_cfg, **common_kwargs}
                    return RFOptimizer(geometry, **rfo_args)
                raise click.BadParameter(f"Unknown optimizer kind '{run_kind}'.")

            # Seed cached IRC endpoint Hessian for RFO when available
            if main_kind == "rfo":
                from pdb2reaction.io.hessian_cache import (
                    load as _hess_load,
                    matches_cart_coords as _hess_matches_coords,
                )
                _cached = _hess_load("irc_endpoint")
                if _cached is not None and not _hess_matches_coords(
                    _cached, geometry.cart_coords
                ):
                    click.echo(
                        "[opt] Cached IRC Hessian does not match the input "
                        "geometry; calculating a fresh Hessian.",
                        err=True,
                    )
                    _cached = None
                if _cached is not None:
                    click.echo("[opt] Reusing IRC endpoint Hessian for RFO seeding.")
                    _active_dofs = _cached.get("active_dofs")
                    _h_raw = _cached["hessian"]
                    if isinstance(_h_raw, torch.Tensor):
                        _h_init = _h_raw.clone()
                    else:
                        _h_init = torch.as_tensor(_h_raw, dtype=torch.float64)
                    if _active_dofs is not None:
                        geometry.within_partial_hessian = {
                            "active_n_dof": len(_active_dofs),
                            "full_n_dof": geometry.cart_coords.size,
                            "active_dofs": _active_dofs,
                            "active_atoms": sorted(set(d // 3 for d in _active_dofs)),
                        }
                    geometry.cart_hessian = _h_init
                    click.echo(f"[opt] Initial Hessian seeded (shape={_h_init.shape[0]}x{_h_init.shape[1]}).")
                    del _h_init

            main_label = "LBFGS" if main_kind == "lbfgs" else "RFO"
            optimizer = _build_optimizer(main_kind)
            click.echo(f"\n====== Optimization ({main_label}) ======\n", narrative=True)
            optimizer.run()
            emit_optimizer_terminal_status(
                "opt",
                converged=getattr(optimizer, "is_converged", None),
                cycles=optimizer_cycle_count(optimizer),
                max_cycles=int(opt_cfg.get("max_cycles", 0)) or None,
            )
            last_optimizer = optimizer
            rigid_projection_info: Dict[str, Any] = {}

            if flatten:
                click.echo("\n====== Optimization (Flatten loop) ======\n", narrative=True)

                geometry.set_calculator(None)
                calc_kwargs_for_flatten = dict(calc_cfg)
                calc_kwargs_for_flatten["out_hess_torch"] = True
                device = _torch_device(calc_cfg.get("device", "auto"))
                freeze_idx = list(geom_cfg.get("freeze_atoms", [])) if len(geom_cfg.get("freeze_atoms", [])) > 0 else None
                masses_amu = _safe_masses_amu(geometry.atomic_numbers)

                def _attach_opt_calc() -> None:
                    geometry.set_calculator(opt_calc)

                def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                    H = _calc_full_hessian_torch(geometry, calc_kwargs_for_flatten, device)
                    freqs_local, modes_local = _frequencies_cm_and_modes(
                        H,
                        geometry.atomic_numbers,
                        geometry.cart_coords.reshape(-1, 3),
                        device,
                        freeze_idx=freeze_idx,
                        tr_projection=geom_cfg["tr_projection"],
                        projection_info=rigid_projection_info,
                    )
                    rigid_projection_info.update({
                        "hessian_space": (
                            "full" if H.shape[0] == 3 * len(geometry.atomic_numbers)
                            else "active"
                        ),
                        "raw_hessian_shape": list(H.shape),
                        "source": "opt_flatten",
                    })
                    del H
                    return freqs_local, modes_local

                freqs_cm, modes = _calc_freqs_and_modes()
                neg_mask = freqs_cm < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm if x < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)]
                click.echo(f"[Imaginary modes] n={n_imag}  ({ims})", narrative=True)

                for it in range(OPT_FLATTEN_MAX_ITER):
                    if n_imag == 0:
                        break
                    click.echo(f"[flatten] iteration {it + 1}/{OPT_FLATTEN_MAX_ITER}")
                    did_flatten = _flatten_all_imag_modes_for_geom(
                        geometry,
                        masses_amu,
                        calc_kwargs_for_flatten,
                        freqs_cm,
                        modes,
                        OPT_FLATTEN_NEG_FREQ_THRESH_CM,
                        OPT_FLATTEN_AMP_ANG,
                    )
                    if not did_flatten:
                        click.echo("[flatten] No eligible imaginary modes to flatten; stopping.")
                        break

                    _attach_opt_calc()
                    optimizer = _build_optimizer(flatten_kind)
                    restart_label = "LBFGS" if flatten_kind == "lbfgs" else "RFO"
                    click.echo(f"\n====== Optimization ({restart_label}, flatten retry) ======\n", narrative=True)
                    optimizer.run()
                    emit_optimizer_terminal_status(
                        "opt",
                        converged=getattr(optimizer, "is_converged", None),
                        cycles=optimizer_cycle_count(optimizer),
                        max_cycles=int(opt_cfg.get("max_cycles", 0)) or None,
                    )
                    last_optimizer = optimizer

                    geometry.set_calculator(None)
                    freqs_cm, modes = _calc_freqs_and_modes()
                    neg_mask = freqs_cm < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)
                    n_imag = int(np.sum(neg_mask))
                    ims = [float(x) for x in freqs_cm if x < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)]
                    click.echo(f"[Imaginary modes] n={n_imag}  ({ims})", narrative=True)

                if n_imag > 0:
                    click.echo(
                        f"[flatten] WARNING: Remaining imaginary modes after {OPT_FLATTEN_MAX_ITER} iterations: {n_imag}",
                        err=True,
                    )
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                click.echo(pretty_block("rigid_projection", rigid_projection_info))

            # Final geometry location (Optimizer sets final_fn during run)
            final_xyz_path = (
                last_optimizer.final_fn
                if isinstance(last_optimizer.final_fn, Path)
                else Path(last_optimizer.final_fn)
            )
            _convert_outputs(
                prepared_input=prepared_input,
                out_dir=out_dir_path,
                dump=bool(opt_cfg["dump"]),
                get_trj_fn=last_optimizer.get_path_for_fn,
                final_xyz_path=final_xyz_path,
            )

            click.echo(format_elapsed("[time] Elapsed Time for Opt", time_start), narrative=True)

            if out_json:
                from pdb2reaction.core.utils import write_result_json
                result_data = {
                    "status": "converged" if last_optimizer.is_converged else "not_converged",
                    "energy_hartree": float(geometry.energy) if geometry.energy is not None else None,
                    "n_opt_cycles": last_optimizer.cur_cycle if hasattr(last_optimizer, "cur_cycle") else None,
                    "opt_mode": opt_cfg.get("opt_mode", opt_mode),
                    "backend": backend,
                    "charge": calc_cfg["charge"],
                    "spin": calc_cfg["spin"],
                    "model": calc_cfg.get("model"),
                    "n_atoms": len(geometry.atoms),
                    "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                    "solvent": calc_cfg.get("solvent", "none"),
                    "thresh": opt_cfg.get("thresh", "gau"),
                    "max_cycles": opt_cfg.get("max_cycles"),
                    "input_file": str(prepared_input.display_path),
                    "files": {
                        "final_geometry_xyz": str(final_xyz_path.name),
                    },
                }
                if rigid_projection_info:
                    result_data["rigid_projection"] = dict(rigid_projection_info)
                if hasattr(last_optimizer, 'max_forces') and last_optimizer.max_forces:
                    result_data["final_max_force"] = float(last_optimizer.max_forces[-1])
                    result_data["final_rms_force"] = float(last_optimizer.rms_forces[-1])
                # Convergence thresholds (numeric values for the named preset)
                if hasattr(last_optimizer, 'convergence') and last_optimizer.convergence:
                    result_data["convergence_thresholds"] = {k: float(v) for k, v in last_optimizer.convergence.items()}
                # Final step convergence values
                if hasattr(last_optimizer, 'max_steps') and last_optimizer.max_steps:
                    result_data["final_max_step"] = float(last_optimizer.max_steps[-1])
                    result_data["final_rms_step"] = float(last_optimizer.rms_steps[-1])
                # Add PDB/CIF/GJF companions if generated.
                for ext in (".pdb", ".cif", ".gjf"):
                    f = out_dir_path / f"final_geometry{ext}"
                    if f.exists():
                        result_data["files"][f"final_geometry_{ext[1:]}"] = f.name
                # Add trajectory files if they exist
                for name in ("optimization_trj.xyz", "optimization.pdb", "optimization.cif"):
                    _tf = out_dir_path / name
                    if _tf.exists():
                        key = name.replace(".", "_").replace("-", "_")
                        result_data["files"][key] = name
                write_result_json(
                    out_dir_path, result_data,
                    command="opt",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

    out_dir_path = Path(out_dir).resolve()
    run_cli(
        _run,
        label="optimization",
        zero_step_exc=ZeroStepLength,
        zero_step_msg="ERROR: Step length fell below the minimum allowed (ZeroStepLength).",
        opt_exc=OptimizationError,
        opt_msg="ERROR: Optimization failed - {exc}",
        out_dir=out_dir_path,
        command="opt",
        time_start=time_start,
    )
