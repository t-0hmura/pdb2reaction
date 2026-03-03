# pdb2reaction/opt.py

"""
Single-structure geometry optimization using LBFGS or RFO with UMA calculator.

Example:
    pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess

For detailed documentation, see: docs/opt.md
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import ast
import inspect
import os
import shutil
import sys

import click
import numpy as np
import torch
import yaml
import time
from click.core import ParameterSource
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, BOHR2ANG, AU2EV

from .defaults import (
    GEOM_KW_DEFAULT,
    CALC_KW_DEFAULT,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    OPT_MODE_ALIASES,
    UMA_CALC_KW,
)
from .backends import create_calculator
from .utils import (
    resolve_freeze_atoms,
    deep_update,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
)
from .cli_utils import run_cli, resolve_yaml_sources, load_merged_yaml_cfg, link_or_copy_file
from .freq import (
    _torch_device,
    _safe_masses_amu,
    _calc_full_hessian_torch,
    _frequencies_cm_and_modes,
    _calc_energy,
)

EV2AU = 1.0 / AU2EV  # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)
OPT_FLATTEN_NEG_FREQ_THRESH_CM = 5.0
OPT_FLATTEN_AMP_ANG = 0.10
OPT_FLATTEN_MAX_ITER = 50

# Note: All defaults imported from defaults.py - no local copies needed



_link_or_copy_file = link_or_copy_file  # backward compat alias


class HarmonicBiasCalculator:
    """Wrap a base UMA calculator with harmonic distance restraints."""

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])

    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        coords = np.array(coords_bohr, dtype=float).reshape(-1, 3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        k = self.k_au_bohr2
        for (i, j, target_ang) in self._pairs:
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]
            rij = float(np.linalg.norm(rij_vec))
            if rij < 1e-14:
                continue
            target_bohr = float(target_ang) * ANG2BOHR
            diff_bohr = rij - target_bohr
            E_bias += 0.5 * k * diff_bohr * diff_bohr
            u = rij_vec / max(rij, 1e-14)
            Fi = -k * diff_bohr * u
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_forces(elem, coords_bohr)
        E0 = float(base["energy"])
        F0 = np.asarray(base["forces"], dtype=float).reshape(-1)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias, "forces": F0 + Fbias}

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        E0 = float(self.base.get_energy(elem, coords_bohr)["energy"])
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias}

    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

    def __getattr__(self, name: str):
        return getattr(self.base, name)


def _parse_dist_freeze(
    args: Sequence[str],
    one_based: bool,
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse --dist-freeze arguments into 0-based pairs with optional targets."""
    parsed: List[Tuple[int, int, Optional[float]]] = []
    for idx, raw in enumerate(args, start=1):
        try:
            obj = ast.literal_eval(raw)
        except Exception as e:
            raise click.BadParameter(f"Invalid literal for --dist-freeze #{idx}: {e}")
        if isinstance(obj, (list, tuple)) and not obj:
            iterable = []
        elif isinstance(obj, (list, tuple)) and isinstance(obj[0], (list, tuple)):
            iterable = obj
        else:
            iterable = [obj]
        for entry in iterable:
            if not (isinstance(entry, (list, tuple)) and len(entry) in (2, 3)):
                raise click.BadParameter(
                    f"--dist-freeze #{idx} entries must be (i,j) or (i,j,target_A): {entry}"
                )
            if not (
                isinstance(entry[0], (int, np.integer))
                and isinstance(entry[1], (int, np.integer))
            ):
                raise click.BadParameter(f"Atom indices in --dist-freeze #{idx} must be integers: {entry}")
            i = int(entry[0])
            j = int(entry[1])
            target = None
            if len(entry) == 3:
                if not isinstance(entry[2], (int, float, np.floating)):
                    raise click.BadParameter(f"Target distance must be numeric in --dist-freeze #{idx}: {entry}")
                target = float(entry[2])
                if target <= 0.0:
                    raise click.BadParameter(f"Target distance must be > 0 in --dist-freeze #{idx}: {entry}")
            if one_based:
                i -= 1
                j -= 1
            if i < 0 or j < 0:
                raise click.BadParameter(f"--dist-freeze #{idx} produced negative index after conversion: {entry}")
            parsed.append((i, j, target))
    return parsed


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
    """Convert outputs (final geometry and trajectory) to PDB/GJF when requested by the input type."""
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
    uma_kwargs: dict,
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
    E_ref = _calc_energy(geom, uma_kwargs)

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
        E_plus = _calc_energy(geom, uma_kwargs)

        geom.coords = minus.reshape(-1)
        E_minus = _calc_energy(geom, uma_kwargs)

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


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Single-structure geometry optimization using LBFGS or RFO.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, _trj.xyz, ...).",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=False,
    help=(
        "Total charge. Required for non-.gjf inputs unless --ligand-charge is provided "
        "(PDB inputs or XYZ/GJF with --ref-pdb)."
    ),
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
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--dist-freeze",
    "dist_freeze_raw",
    type=str,
    multiple=True,
    default=(),
    show_default=False,
    help="Python-like list(s) of (i,j,target_A) to restrain distances (target optional).",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret --dist-freeze indices as 1-based (default) or 0-based.",
)
@click.option(
    "--bias-k",
    type=float,
    default=10.0,
    show_default=True,
    help="Harmonic restraint strength k [eV/Å^2] for --dist-freeze.",
)
@click.option(
    "--freeze-links/--no-freeze-links",
    default=True,
    show_default=True,
    help="Freeze parent atoms of link hydrogens (PDB only).",
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
    "--out-dir",
    type=str,
    default="./result_opt/",
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--thresh",
    type=str,
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
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running optimization.",
)
@click.option("--backend", type=click.Choice(["uma", "orb", "mace", "aimnet2"]), default="uma",
              help="MLIP backend.")
@click.option("--solvent", default="none",
              help="Implicit solvent name for xTB correction (e.g. 'water'). 'none' to disable.")
@click.option("--solvent-model", "solvent_model", default="alpb", type=click.Choice(["alpb", "cpcmx"]),
              help="xTB solvent model.")
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
    backend: str,
    solvent: str,
    solvent_model: str,
) -> None:
    time_start = time.perf_counter()

    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception:
            return False

    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
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
            geom_input_path = prepared_input.geom_path
            source_path = prepared_input.source_path

            dist_freeze = _parse_dist_freeze(dist_freeze_raw, one_based=bool(one_based))

            # --------------------------
            # 1) Assemble configuration (defaults < config < CLI(explicit) < override)
            # --------------------------
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
            if _is_param_explicit("workers"):
                calc_cfg["workers"] = int(workers)
            if _is_param_explicit("workers_per_node"):
                calc_cfg["workers_per_node"] = int(workers_per_node)
            if _is_param_explicit("backend"):
                calc_cfg["backend"] = backend
            if _is_param_explicit("solvent"):
                calc_cfg["solvent"] = solvent
            if _is_param_explicit("solvent_model"):
                calc_cfg["solvent_model"] = solvent_model
            if _is_param_explicit("max_cycles"):
                opt_cfg["max_cycles"] = int(max_cycles)
            if _is_param_explicit("dump"):
                opt_cfg["dump"] = bool(dump)
            if _is_param_explicit("out_dir"):
                opt_cfg["out_dir"] = out_dir
            if _is_param_explicit("thresh") and thresh is not None:
                opt_cfg["thresh"] = str(thresh)

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

            # Normalize freeze_atoms and optionally add link-parent indices for PDB inputs
            resolve_freeze_atoms(geom_cfg, source_path, freeze_links)
            calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))

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
                            "will_run_optimization": True,
                            "will_write_final_geometry": True,
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. Optimization execution was skipped.")
                return

            # --------------------------
            # 2) Prepare geometry
            # --------------------------
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

            # Attach UMA calculator
            base_calc = create_calculator(**calc_cfg)
            geometry.set_calculator(base_calc)

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

            # --------------------------
            # 3) Build optimizer
            # --------------------------
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

            # --------------------------
            # 4) Run optimization
            # --------------------------
            main_label = "LBFGS" if main_kind == "lbfgs" else "RFO"
            optimizer = _build_optimizer(main_kind)
            click.echo(f"\n====== Optimization ({main_label}) started ======\n")
            optimizer.run()
            click.echo(f"\n====== Optimization ({main_label}) finished ======\n")
            last_optimizer = optimizer

            # --------------------------
            # 5) Flatten loop (all imaginary modes)
            # --------------------------
            if flatten:
                click.echo("\n====== Optimization (Flatten loop) started ======\n")

                geometry.set_calculator(None)
                uma_kwargs_for_flatten = dict(calc_cfg)
                uma_kwargs_for_flatten["out_hess_torch"] = True
                device = _torch_device(calc_cfg.get("device", "auto"))
                freeze_idx = list(geom_cfg.get("freeze_atoms", [])) if len(geom_cfg.get("freeze_atoms", [])) > 0 else None
                masses_amu = _safe_masses_amu(geometry.atomic_numbers)

                def _attach_opt_calc() -> None:
                    geometry.set_calculator(opt_calc)

                def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                    H = _calc_full_hessian_torch(geometry, uma_kwargs_for_flatten, device)
                    freqs_local, modes_local = _frequencies_cm_and_modes(
                        H,
                        geometry.atomic_numbers,
                        geometry.cart_coords.reshape(-1, 3),
                        device,
                        freeze_idx=freeze_idx,
                    )
                    del H
                    return freqs_local, modes_local

                freqs_cm, modes = _calc_freqs_and_modes()
                neg_mask = freqs_cm < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm if x < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)]
                click.echo(f"[Imaginary modes] n={n_imag}  ({ims})")

                for it in range(OPT_FLATTEN_MAX_ITER):
                    if n_imag == 0:
                        break
                    click.echo(f"[flatten] iteration {it + 1}/{OPT_FLATTEN_MAX_ITER}")
                    did_flatten = _flatten_all_imag_modes_for_geom(
                        geometry,
                        masses_amu,
                        uma_kwargs_for_flatten,
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
                    click.echo(f"\n====== Optimization ({restart_label}) restarted ======\n")
                    optimizer.run()
                    click.echo(f"\n====== Optimization ({restart_label}) finished ======\n")
                    last_optimizer = optimizer

                    geometry.set_calculator(None)
                    freqs_cm, modes = _calc_freqs_and_modes()
                    neg_mask = freqs_cm < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)
                    n_imag = int(np.sum(neg_mask))
                    ims = [float(x) for x in freqs_cm if x < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)]
                    click.echo(f"[Imaginary modes] n={n_imag}  ({ims})")

                if n_imag > 0:
                    click.echo(
                        f"[flatten] WARNING: Remaining imaginary modes after {OPT_FLATTEN_MAX_ITER} iterations: {n_imag}",
                        err=True,
                    )
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                click.echo("\n====== Optimization (Flatten loop) finished ======\n")

            # --------------------------
            # 6) Post-processing: PDB conversions (if input is PDB)
            # --------------------------
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

            # summary.md and key_* outputs are disabled.
            click.echo(format_elapsed("[time] Elapsed Time for Opt", time_start))

    run_cli(
        _run,
        label="optimization",
        zero_step_exc=ZeroStepLength,
        zero_step_msg="ERROR: Step length fell below the minimum allowed (ZeroStepLength).",
        opt_exc=OptimizationError,
        opt_msg="ERROR: Optimization failed - {e}",
    )
