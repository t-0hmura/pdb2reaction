# pdb2reaction/opt.py

"""
Single-structure geometry optimization using LBFGS or RFO with UMA calculator.

Example:
    pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode light

For detailed documentation, see: docs/opt.md
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import ast
import inspect
import sys

import click
import numpy as np
import yaml
import time

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
from .uma_pysis import uma_pysis
from .utils import (
    resolve_freeze_atoms,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    normalize_choice,
    prepared_cli_input,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
)
from .cli_utils import run_cli

EV2AU = 1.0 / AU2EV  # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)

# Note: All defaults imported from defaults.py - no local copies needed


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

    # optimization.trj → optimization.pdb (if dump)
    if dump and needs_pdb:
        trj_path = get_trj_fn("optimization.trj")
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
            click.echo("[convert] WARNING: 'optimization.trj' not found; skipping conversion.", err=True)


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
    help="Input structure file (.pdb, .xyz, .trj, ...).",
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
    default=1,
    show_default=True,
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
    "--one-based",
    "one_based",
    type=click.BOOL,
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
    "--freeze-links",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Freeze the parent atoms of link hydrogens (PDB only).",
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
    "--max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum number of optimization cycles.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Optimization mode: 'light' (=LBFGS) or 'heavy' (=RFO).",
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Write optimization trajectory to 'optimization.trj'.",
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
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra arguments (sections: geom, calc, opt, lbfgs, rfo).",
)
def cli(
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
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
) -> None:
    time_start = time.perf_counter()
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
            # 1) Assemble configuration (create fresh copies for merging)
            # --------------------------
            yaml_cfg = load_yaml_dict(args_yaml)
            geom_cfg = dict(GEOM_KW_DEFAULT)
            calc_cfg = dict(UMA_CALC_KW)
            opt_cfg = dict(OPT_BASE_KW)
            lbfgs_cfg = dict(LBFGS_KW)
            rfo_cfg = dict(RFO_KW)

            # CLI overrides (defaults ← CLI)
            calc_cfg["charge"] = resolved_charge
            calc_cfg["spin"] = resolved_spin
            calc_cfg["workers"] = int(workers)
            calc_cfg["workers_per_node"] = int(workers_per_node)
            opt_cfg["max_cycles"] = int(max_cycles)
            opt_cfg["dump"] = bool(dump)
            opt_cfg["out_dir"] = out_dir
            if thresh is not None:
                opt_cfg["thresh"] = str(thresh)

            # YAML has highest precedence (defaults ← CLI ← YAML)
            apply_yaml_overrides(
                yaml_cfg,
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

            # Normalize and select optimizer kind
            kind = normalize_choice(
                opt_mode,
                param="--opt-mode",
                alias_groups=OPT_MODE_ALIASES,
                allowed_hint="light|heavy",
            )

            # Pretty-print the resolved configuration
            out_dir_path = Path(opt_cfg["out_dir"]).resolve()
            click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
            click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
            click.echo(pretty_block("opt", {**opt_cfg, "out_dir": str(out_dir_path)}))
            echo_sopt = dict(lbfgs_cfg if kind == "lbfgs" else rfo_cfg)
            echo_sopt.update(opt_cfg)
            echo_sopt["out_dir"] = str(out_dir_path)
            click.echo(pretty_block(kind, echo_sopt))
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
            base_calc = uma_pysis(**calc_cfg)
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
            common_kwargs = dict(opt_cfg)
            # Ensure paths (strings) are OK; Optimizer expects str, not Path
            common_kwargs["out_dir"] = str(out_dir_path)

            if kind == "lbfgs":
                lbfgs_args = {**lbfgs_cfg, **common_kwargs}
                optimizer = LBFGS(geometry, **lbfgs_args)
            else:
                rfo_args = {**rfo_cfg, **common_kwargs}
                optimizer = RFOptimizer(geometry, **rfo_args)

            # --------------------------
            # 4) Run optimization
            # --------------------------
            click.echo("\n=== Optimization started ===\n")
            optimizer.run()
            click.echo("\n=== Optimization finished ===\n")

            # --------------------------
            # 5) Post-processing: PDB conversions (if input is PDB)
            # --------------------------
            # Final geometry location (Optimizer sets final_fn during run)
            final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
            _convert_outputs(
                prepared_input=prepared_input,
                out_dir=out_dir_path,
                dump=bool(opt_cfg["dump"]),
                get_trj_fn=optimizer.get_path_for_fn,
                final_xyz_path=final_xyz_path,
            )

            click.echo(format_elapsed("[time] Elapsed Time for Opt", time_start))

    run_cli(
        _run,
        label="optimization",
        zero_step_exc=ZeroStepLength,
        zero_step_msg="ERROR: Step length fell below the minimum allowed (ZeroStepLength).",
        opt_exc=OptimizationError,
        opt_msg="ERROR: Optimization failed - {e}",
    )
