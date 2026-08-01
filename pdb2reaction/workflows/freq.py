# pdb2reaction/workflows/freq.py

"""
Vibrational frequency analysis with PHVA support and thermochemistry output.

Example:
    pdb2reaction freq -i a.pdb -q 0 -m 1

For detailed documentation, see: docs/freq.md
"""

from __future__ import annotations

import gc
import logging
import sys
import textwrap
from pathlib import Path
from typing import Any, Optional

import click
from pdb2reaction.core.output import emit
import numpy as np
import torch
import yaml
import time

# ---------------- pysisyphus / pdb2reaction imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AMU2AU

# Compatibility re-exports: the pure normal-mode kernel lives in the
# lower bundled-engine module ``pysisyphus.normal_modes`` (a sibling of
# ``pysisyphus.tr_projection``), which does not import ``pdb2reaction``.
# It is re-exported here so existing callers of
# ``pdb2reaction.workflows.freq`` (opt/tsopt/tests) keep working unchanged and
# resolve to the SAME function objects as the lower implementation.
from pysisyphus.normal_modes import (  # noqa: F401
    _safe_masses_amu,
    _mw_projected_hessian,
    _mass_weighted_hessian,
    _frequencies_cm_and_modes,
    _mw_mode_to_cart,
)

# local helpers from pdb2reaction
from pdb2reaction.backends import create_calculator
from pdb2reaction.core.defaults import GEOM_KW_DEFAULT, FREQ_CALC_KW, FREQ_KW, THERMO_KW, apply_backend_defaults
from pdb2reaction.core.utils import (
    apply_yaml_overrides,
    convert_xyz_like_outputs,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin,
    validate_charge_spin_for_prepared,
    set_convert_file_enabled,
    resolve_freeze_atoms,
    cli_param_overridden,
    yaml_freeze_to_internal,
    _parse_freeze_atoms,
    merge_freeze_atom_indices,
    echo_resolved_device,
)
from pdb2reaction.cli.common_options import add_ml_charge_spin_options, add_precision_option, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option, add_coord_type_option
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, render_cli_exception

logger = logging.getLogger(__name__)


def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)


def _calc_full_hessian_torch(
    geom,
    calc_kwargs: dict,
    device: torch.device,
    *,
    refresh_geom_meta: bool = False,
    calculator=None,
    cache_geometry: bool = True,
) -> torch.Tensor:
    """
    Return Hessian as torch.Tensor in Hartree/Bohr^2 (3N or 3N_act square).

    Reuses Geometry-level Hessian cache when available so repeated requests on the
    same coordinates do not trigger extra calculator evaluations.

    Parameters
    ----------
    refresh_geom_meta : bool, default False
        Accepted for compatibility; currently a no-op here. When
        ``cache_geometry`` is true, ``geom.set_results(results)``
        refreshes ``geom.within_partial_hessian`` via ``Geometry.set_results``;
        standalone frequency analysis deliberately leaves that dense cache off.
    calculator : optional
        An already-resolved evaluator for the active PES. When supplied, the
        Geometry cache is intentionally bypassed because a coordinate-only
        cache cannot prove that it represents this evaluator's restraints.
    cache_geometry : bool, default True
        Store the result on ``geom`` and return a protective clone. Standalone
        frequency analysis sets this to false because it consumes the Hessian
        once and must not keep a second dense GPU copy alive.
    """
    def _to_torch(H_obj: Any, clone: bool = True) -> torch.Tensor:
        if isinstance(H_obj, torch.Tensor):
            t = H_obj.to(device=device)
            return t.clone() if clone else t
        return torch.as_tensor(H_obj, device=device)

    cached = getattr(geom, "_hessian", None) if calculator is None else None
    if cached is not None:
        return _to_torch(cached, clone=True)

    calc = calculator
    if calc is None:
        kw = dict(calc_kwargs or {})
        kw["out_hess_torch"] = True
        calc = create_calculator(**kw)
        echo_resolved_device()
    results = calc.get_hessian(geom.atoms, geom.cart_coords)

    # Keep Geometry cache in sync so optimizers/freq analysis can share one Hessian
    # evaluation on unchanged coordinates.
    if cache_geometry:
        try:
            geom.set_results(results)
        except Exception:
            import logging
            logging.getLogger(__name__).warning(
                "Failed to set Hessian results on Geometry cache", exc_info=True
            )

    # Clone so downstream in-place mass-weighting / TR projection does not
    # poison the Hessian cached on the Geometry object.
    return _to_torch(results["hessian"], clone=cache_geometry)


def _calc_energy(geom, calc_kwargs: dict, calc=None) -> float:
    """
    Compute electronic energy (Hartree) from the configured calculator.
    """
    if calc is None:
        calc = create_calculator(**calc_kwargs)
    geom.set_calculator(calc)
    try:
        E = float(geom.energy)
    finally:
        geom.set_calculator(None)
    if not np.isfinite(E):
        raise ValueError("Electronic energy must be finite for thermochemistry.")
    return E


def _fmt_ha(x: float) -> str:
    return f"{float(x): .6f} Ha"


def _fmt_cal(x: float) -> str:
    return f"{float(x): .2f} cal/mol"


def _fmt_calK(x: float) -> str:
    return f"{float(x): .2f} cal/(mol*K)"


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            amplitude_ang: float = 0.8,
                            n_frames: int = 20,
                            comment: str = "mode",
                            ref_pdb: Optional[Path] = None,
                            write_pdb: bool = True,
                            prepared_input: Optional["PreparedInputStructure"] = None,
                            out_pdb: Optional[Path] = None) -> None:
    """
    Write a single mode animation as _trj.xyz (XYZ-like) and optionally .pdb.

    If `ref_pdb` is provided and is a .pdb file, the .pdb is generated by
    converting the _trj.xyz using the input PDB as the template.
    Set `write_pdb=False` to skip PDB generation.
    """
    ref_ang = geom.cart_coords.reshape(-1, 3) * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # _trj.xyz (concatenated XYZ-like trajectory)
    try:
        from pysisyphus.xyzloader import make_trj_str  # type: ignore
        amp_ang = amplitude_ang
        steps = np.sin(2.0 * np.pi * np.arange(n_frames) / n_frames)[:, None, None] * (amp_ang * mode[None, :, :])
        traj_ang = ref_ang[None, :, :] + steps  # (T,N,3) in Å
        comments = [f"{comment}  frame={i+1}/{n_frames}" for i in range(n_frames)]
        trj_str = make_trj_str(geom.atoms, traj_ang, comments=comments)
        out_trj.write_text(trj_str, encoding="utf-8")
    except Exception:
        with out_trj.open("w", encoding="utf-8") as f:
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                coords = ref_ang + phase * amplitude_ang * mode
                f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
                for sym, (x, y, z) in zip(geom.atoms, coords):
                    f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")

    needs_pdb = write_pdb and out_pdb is not None

    if not needs_pdb:
        return

    ref_for_conv = ref_pdb if (ref_pdb and ref_pdb.suffix.lower() == ".pdb") else None
    try:
        convert_xyz_like_outputs(
            out_trj,
            prepared_input,  # type: ignore[arg-type]
            ref_pdb_path=ref_for_conv,
            out_pdb_path=out_pdb if needs_pdb else None,
        )
    except Exception as e:
        click.echo(
            f"[convert] WARNING: Failed to convert mode trajectory '{out_trj.name}' to PDB: {e}",
            err=True,
        )


# Geometry defaults (local copy for CLI)
GEOM_KW = dict(GEOM_KW_DEFAULT)

# Calc defaults for freq (from defaults.py)
CALC_KW = FREQ_CALC_KW


def _validated_symmetry_number(value: object) -> int:
    """Return an external rotational symmetry number accepted by thermoanalysis."""
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise click.UsageError(
            "thermo.symmetry_number must be an integer greater than or equal to 1."
        )
    return int(value)


def _validated_thermo_condition(value: object, *, name: str) -> float:
    """Return a finite, strictly positive thermochemistry state variable."""
    if isinstance(value, bool):
        raise click.UsageError(
            f"thermo.{name} must be a finite number greater than zero."
        )
    try:
        resolved = float(value)
    except (TypeError, ValueError) as exc:
        raise click.UsageError(
            f"thermo.{name} must be a finite number greater than zero."
        ) from exc
    if not np.isfinite(resolved) or resolved <= 0.0:
        raise click.UsageError(
            f"thermo.{name} must be a finite number greater than zero."
        )
    return resolved


def _prepare_thermo_output_paths(
    out_dir: Path,
    *,
    protected_inputs: tuple[Optional[Path], ...] = (),
) -> tuple[Path, Path]:
    """Invalidate a prior thermochemistry generation before real work starts."""
    thermo_yaml = Path(out_dir) / "thermoanalysis.yaml"
    thermo_yaml_tmp = Path(out_dir) / "thermoanalysis.yaml.tmp"
    reserved = {thermo_yaml.resolve(), thermo_yaml_tmp.resolve()}
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise click.UsageError(
                f"Configuration input {protected} collides with a reserved "
                f"frequency output path under {out_dir}."
            )
    thermo_yaml.unlink(missing_ok=True)
    thermo_yaml_tmp.unlink(missing_ok=True)
    return thermo_yaml, thermo_yaml_tmp


def _prepare_frequency_output_paths(
    out_dir: Path,
    *,
    protected_inputs: tuple[Optional[Path], ...] = (),
) -> tuple[Path, Path]:
    """Invalidate every public artifact owned by one real frequency run."""
    out_dir = Path(out_dir)
    owned = [
        out_dir / "frequencies_cm-1.txt",
        out_dir / "result.json",
        out_dir / "summary.json",
        *out_dir.glob("mode_*cm-1_trj.xyz"),
        *out_dir.glob("mode_*cm-1.pdb"),
        *out_dir.glob("mode_*cm-1.cif"),
    ]
    reserved = {path.resolve() for path in owned}
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise click.UsageError(
                f"Input {protected} collides with a reserved frequency output "
                f"path under {out_dir}."
            )
    thermo_paths = _prepare_thermo_output_paths(
        out_dir,
        protected_inputs=protected_inputs,
    )
    for path in owned:
        path.unlink(missing_ok=True)
    return thermo_paths



@click.command(
    help="Vibrational frequency analysis and mode writer (+ default thermochemistry summary).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
# NOTE: the workers/workers_per_node help text below is intentionally explicit:
# freq DOES support workers>1 (FD Hessian path). Do not shorten the parenthetical
# when trimming docs — the verbose form is what disambiguates the supported case.
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
    default=CALC_KW["workers"],
    show_default=True,
    help="MLIP predictor workers; >1 spawns a parallel predictor. NOTE: with UMA, workers>1 plus an explicit Analytical Hessian request is an error; use workers=1 or FiniteDifference.",
)
@click.option(
    "--workers-per-node",
    "workers_per_node",
    type=int,
    default=CALC_KW["workers_per_node"],
    show_default=True,
    help="Workers per node when using a parallel MLIP predictor (workers>1).",
)
@click.option(
    "--freeze-links/--no-freeze-links",
    "freeze_links",
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
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB companions when a PDB template is available.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB/mmCIF topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option("--max-write", type=int, default=FREQ_KW["max_write"], show_default=True,
              help="How many modes to export (after sorting per --sort).")
@click.option("--amplitude-ang", type=float, default=FREQ_KW["amplitude_ang"], show_default=True,
              help="Animation amplitude (Å) used for both _trj.xyz and .pdb.")
@click.option("--n-frames", type=int, default=FREQ_KW["n_frames"], show_default=True,
              help="Number of frames per mode animation.")
@click.option("--sort", type=click.Choice(["value", "abs"]), default="value", show_default=True,
              help="Sort modes by 'value' (cm^-1) or by absolute value.")
@click.option("-o", "--out-dir", type=str, default=FREQ_KW["out_dir"], show_default=True, help="Output directory.")
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Base YAML configuration file applied before explicit CLI options.",
)
# Thermochemistry options
@click.option("--temperature", type=float, default=THERMO_KW["temperature"], show_default=True,
              help="Temperature (K) for thermochemistry summary.")
@click.option("--pressure", "pressure_atm",
              type=float, default=THERMO_KW["pressure_atm"], show_default=True,
              help="Pressure (atm) for thermochemistry summary.")
@click.option(
    "--symmetry-number",
    type=click.IntRange(min=1),
    default=THERMO_KW["symmetry_number"],
    show_default=True,
    help="External rotational symmetry number used in the thermochemistry partition function.",
)
@click.option(
    "--dump/--no-dump",
    "dump",
    default=THERMO_KW["dump"],
    show_default=True,
    help="When True, write 'thermoanalysis.yaml' under out-dir.",
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
    help="Validate options and print the execution plan without running frequency analysis.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
# Hessian calculation mode
@click.option("--hessian-calc-mode",
              type=click.Choice(["FiniteDifference", "Analytical"], case_sensitive=False),
              default=None,
              help="How the ML backend computes Hessian. Defaults to 'FiniteDifference' (can also be set via YAML).")
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
@add_coord_type_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_node: int,
    spin: Optional[int],
    freeze_links: bool,
    freeze_atoms_text: Optional[str],
    convert_files: bool,
    ref_pdb: Optional[Path],
    max_write: int,
    amplitude_ang: float,
    n_frames: int,
    sort: str,
    out_dir: str,
    config_yaml: Optional[Path],
    # thermo
    temperature: float,
    pressure_atm: float,
    symmetry_number: int,
    dump: bool,
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    # hessian
    hessian_calc_mode: Optional[str],
    # backend
    backend: str,
    solvent: str,
    solvent_model: str,
    precision: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
    cli_coord_type: Optional[str],
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

    time_start = time.perf_counter()
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    # Click closes the context for normal returns, dry-run, UsageError/SystemExit,
    # and execution exceptions. Register ownership immediately so every path
    # after successful preparation has the same cleanup boundary.
    ctx.call_on_close(prepared_input.cleanup)
    apply_ref_pdb_override(prepared_input, ref_pdb)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin(
        prepared_input,
        charge,
        spin,
        ligand_charge=ligand_charge,
        prefix="[freq]",
    )
    validate_charge_spin_for_prepared(prepared_input, charge, spin)

    geom_cfg = dict(GEOM_KW)
    calc_cfg = dict(CALC_KW)
    freq_cfg = dict(FREQ_KW)
    thermo_cfg = dict(THERMO_KW)

    apply_yaml_overrides(
        config_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",),)),
        ],
    )
    _config_thermo = config_layer_cfg.get("thermo")
    symmetry_number_source = (
        "config"
        if isinstance(_config_thermo, dict)
        and "symmetry_number" in _config_thermo
        else "default"
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
    # Precision: the `--precision` CLI flag wins, else the config's calc.precision
    # (see apply_effective_precision — the `all` pipeline propagates precision via
    # the config, invoking children with --config and no --precision, so orb child
    # stages silently fell back to float32-high/TF32 before this).
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
    if cli_param_overridden(ctx, "max_write"):
        freq_cfg["max_write"] = int(max_write)
    if cli_param_overridden(ctx, "amplitude_ang"):
        freq_cfg["amplitude_ang"] = float(amplitude_ang)
    if cli_param_overridden(ctx, "n_frames"):
        freq_cfg["n_frames"] = int(n_frames)
    if cli_param_overridden(ctx, "sort"):
        freq_cfg["sort"] = str(sort)
    if cli_param_overridden(ctx, "out_dir"):
        freq_cfg["out_dir"] = str(out_dir)
    if cli_param_overridden(ctx, "cli_coord_type") and cli_coord_type is not None:
        geom_cfg["coord_type"] = str(cli_coord_type).lower()
    if cli_param_overridden(ctx, "temperature"):
        thermo_cfg["temperature"] = float(temperature)
    if cli_param_overridden(ctx, "pressure_atm"):
        thermo_cfg["pressure_atm"] = float(pressure_atm)
    if cli_param_overridden(ctx, "symmetry_number"):
        thermo_cfg["symmetry_number"] = int(symmetry_number)
        symmetry_number_source = "cli"
    if cli_param_overridden(ctx, "dump"):
        thermo_cfg["dump"] = bool(dump)

    # `charge` / `spin` were already resolved (CLI -q/-m, gjf metadata,
    # --ligand-charge derivation) at the top of cli via resolve_charge_spin
    # and validated via validate_charge_spin_for_prepared. Assign directly;
    # an earlier .get("charge", charge) idiom silently returned the
    # CALC_KW default 0 when -q was not passed.
    calc_cfg["charge"] = int(charge)
    calc_cfg["spin"] = int(spin)

    apply_yaml_overrides(
        override_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",),)),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",),)),
        ],
    )
    _override_thermo = override_layer_cfg.get("thermo")
    if (
        isinstance(_override_thermo, dict)
        and "symmetry_number" in _override_thermo
    ):
        symmetry_number_source = "override"
    thermo_cfg["symmetry_number"] = _validated_symmetry_number(
        thermo_cfg.get("symmetry_number")
    )
    thermo_cfg["temperature"] = _validated_thermo_condition(
        thermo_cfg.get("temperature"), name="temperature"
    )
    thermo_cfg["pressure_atm"] = _validated_thermo_condition(
        thermo_cfg.get("pressure_atm"), name="pressure_atm"
    )
    max_write_value = freq_cfg.get("max_write")
    n_frames_value = freq_cfg.get("n_frames")
    amplitude_value = freq_cfg.get("amplitude_ang")
    if (
        isinstance(max_write_value, bool)
        or not isinstance(max_write_value, (int, np.integer))
        or int(max_write_value) < 0
    ):
        raise click.BadParameter(
            f"freq.max_write must be a non-negative integer, got {max_write_value!r}."
        )
    if (
        isinstance(n_frames_value, bool)
        or not isinstance(n_frames_value, (int, np.integer))
        or int(n_frames_value) < 1
    ):
        raise click.BadParameter(
            f"freq.n_frames must be a positive integer, got {n_frames_value!r}."
        )
    try:
        amplitude_value = float(amplitude_value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise click.BadParameter(
            f"freq.amplitude_ang must be finite, got {amplitude_value!r}."
        ) from exc
    if not np.isfinite(amplitude_value):
        raise click.BadParameter(
            f"freq.amplitude_ang must be finite, got {amplitude_value!r}."
        )
    freq_cfg["max_write"] = int(max_write_value)
    freq_cfg["n_frames"] = int(n_frames_value)
    freq_cfg["amplitude_ang"] = amplitude_value
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
    resolve_freeze_atoms(geom_cfg, source_path, freeze_links, on_error="warn")

    # Ensure calc config reflects the geometry freeze list used in the run.
    calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
    calc_cfg["return_partial_hessian"] = True

    out_dir_path = Path(freq_cfg.get("out_dir", out_dir)).resolve()

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
                    "input_geometry": str(geom_input_path),
                    "output_dir": str(out_dir_path),
                    "freeze_links": bool(freeze_links),
                    "convert_files": bool(convert_files),
                    "will_run_frequency_analysis": True,
                    "will_write_modes": True,
                    "will_dump_thermo_yaml": bool(thermo_cfg.get("dump", False)),
                    "tr_projection": geom_cfg["tr_projection"],
                },
            )
        )
        click.echo("[dry-run] Validation complete. Frequency execution was skipped.")
        emit(
            format_elapsed("[time] Elapsed Time for Freq", time_start),
            narrative=True,
        )
        return

    out_dir_path.mkdir(parents=True, exist_ok=True)
    # A real invocation owns this exact output generation.  Invalidate both
    # published and staged thermochemistry before any calculator or Hessian
    # operation can fail, including under --no-dump.
    _thermo_yaml, _thermo_yaml_tmp = _prepare_frequency_output_paths(
        out_dir_path,
        protected_inputs=(
            input_path,
            geom_input_path,
            source_path,
            ref_pdb,
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

    # Pretty-print config summary
    click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
    click.echo(pretty_block("calc", format_geom_for_echo(calc_cfg)))
    click.echo(pretty_block("freq", {**freq_cfg, "out_dir": str(out_dir_path)}))
    thermo_block = {
        "temperature": thermo_cfg["temperature"],
        "pressure_atm": thermo_cfg["pressure_atm"],
        "symmetry_number": thermo_cfg["symmetry_number"],
        "dump": thermo_cfg["dump"],
    }
    click.echo(pretty_block("thermo", thermo_block))

    coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
    coord_kwargs = dict(geom_cfg)
    coord_kwargs.pop("coord_type", None)
    geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

    # Masses (AU tensor for TR projection & MW->Cart conversion)
    masses_amu = _safe_masses_amu(geometry.atomic_numbers)
    masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=torch.float32)
    device = _torch_device(calc_cfg.get("device", "auto"))
    masses_au_t = masses_au_t.to(device=device)

    try:
        freeze_list = list(calc_cfg.get("freeze_atoms", []))
        from pdb2reaction.io.hessian_cache import (
            load_matching as _hess_load_matching,
            identity_from_context as _hess_identity,
        )
        _ts_calc_cfg = dict(calc_cfg)
        _ts_calc_cfg.setdefault("freeze_atoms", freeze_list)
        # A cached TS Hessian is reused only when the full evaluation identity
        # (run/system/evaluator/active space/potential) matches; the
        # all-workflow may pass the TS through a three-decimal PDB, so the
        # coordinate field keeps the wider bohr tolerance.
        _cached_ts = _hess_load_matching(
            "ts",
            _hess_identity(geometry, _ts_calc_cfg, role="ts"),
            atol=1.1e-3,
        )
        if _cached_ts is not None:
            emit("[freq] Reusing cached TS Hessian.", narrative=True)
            H = _cached_ts["hessian"]
            if isinstance(H, torch.Tensor):
                H = H.to(device=device)
            else:
                H = torch.as_tensor(H, device=device)
        else:
            H = _calc_full_hessian_torch(
                geometry,
                calc_cfg,
                device,
                cache_geometry=False,
            )
        coords_bohr = geometry.cart_coords.reshape(-1, 3)

        # PHVA: use the freeze list to carve out the active subspace and apply TR projection there.
        _n_atoms = len(geometry.atomic_numbers)
        _n_frozen = len(set(int(i) for i in freeze_list))
        _n_active = max(_n_atoms - _n_frozen, 0)
        emit(
            f"[freq] Hessian ready: shape={tuple(H.shape)}, active_atoms={_n_active}/{_n_atoms}, "
            f"frozen_atoms={_n_frozen}, active_dof={3 * _n_active}",
            detail=True,
        )
        _rigid_projection = {}
        freqs_cm, modes_mw = _frequencies_cm_and_modes(
            H,
            geometry.atomic_numbers,
            coords_bohr,
            device,
            freeze_idx=freeze_list if len(freeze_list) > 0 else None,
            tr_projection=geom_cfg["tr_projection"],
            projection_info=_rigid_projection,
        )
        _rigid_projection.update(
            {
                "hessian_space": (
                    "full" if H.shape[0] == 3 * len(geometry.atomic_numbers)
                    else "active"
                ),
                "hessian_shape": list(H.shape),
                "hessian_source": "cache" if _cached_ts is not None else "fresh",
                "hessian_representation": "cartesian-unweighted-unprojected",
            }
        )
        click.echo(
            "[freq] Rigid projection: "
            f"treatment={_rigid_projection['treatment']}, "
            f"rank={_rigid_projection['effective_rank']}, "
            f"full_rigid_rank={_rigid_projection['full_rigid_rank']}."
        )

        # Normal-mode data no longer needs the Cartesian Hessian. Geometry
        # also owns the evaluator result, so clear that cache before the
        # thermochemistry energy evaluation creates another model.
        geometry.clear()
        del H
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        if freq_cfg["sort"] == "abs":
            order = np.argsort(np.abs(freqs_cm))
        else:
            order = np.argsort(freqs_cm)

        n_write = int(min(freq_cfg["max_write"], len(order)))
        _imag = [f for f in freqs_cm if f < 0]
        _preview_n = min(20, len(order))
        _freq_preview = ", ".join(f"{float(freqs_cm[j]):+.1f}" for j in order[:_preview_n])
        _suffix = ", ..." if len(order) > _preview_n else ""
        emit(
            f"[freq] {len(freqs_cm)} modes ({len(_imag)} imaginary); "
            f"first {_preview_n} by {freq_cfg['sort']}: [{_freq_preview}{_suffix}] cm⁻¹; "
            f"full list: {out_dir_path / 'frequencies_cm-1.txt'}",
            narrative=True,
        )
        if len(order) > _preview_n:
            _freq_str = ", ".join(f"{float(freqs_cm[j]):+.1f}" for j in order)
            click.echo(f"[freq:all] [{_freq_str}] cm⁻¹")
        emit(f"[INFO] Writing first {n_write} modes ({freq_cfg['sort']} ascending).", detail=True)

        # Reference PDB (only when input is PDB)
        ref_pdb = source_path if source_path.suffix.lower() == ".pdb" else None
        write_pdb = ref_pdb is not None

        # Write modes and retain their authoritative relative paths for
        # result.json/Results consumers.
        _mode_output_files: list[str] = []
        for k, idx in enumerate(order[:n_write], start=1):
            freq = float(freqs_cm[idx])
            mode_cart_3N = _mw_mode_to_cart(modes_mw[idx], masses_au_t)
            out_trj = out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1_trj.xyz"
            out_pdb = (
                out_dir_path / f"mode_{k:04d}_{freq:+.2f}cm-1.pdb"
                if write_pdb
                else None
            )
            _write_mode_trj_and_pdb(
                geometry,
                mode_cart_3N,
                out_trj,
                amplitude_ang=freq_cfg["amplitude_ang"],
                n_frames=freq_cfg["n_frames"],
                comment=f"mode {k}  {freq:+.2f} cm-1",
                ref_pdb=ref_pdb,
                write_pdb=write_pdb,
                prepared_input=prepared_input,
                out_pdb=out_pdb,
            )
            _mode_output_files.append(out_trj.name)
            if out_pdb is not None and out_pdb.is_file():
                _mode_output_files.append(out_pdb.name)
                out_cif = out_pdb.with_suffix(".cif")
                if out_cif.is_file():
                    _mode_output_files.append(out_cif.name)
        (out_dir_path / "frequencies_cm-1.txt").write_text(
            "\n".join(f"{i+1:4d}  {float(freqs_cm[j]):+12.4f}" for i, j in enumerate(order)),
            encoding="utf-8"
        )

        del modes_mw
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        _thermo_data = None  # populated below if thermoanalysis succeeds
        try:
            from thermoanalysis.QCData import QCData
            from thermoanalysis.thermo import thermochemistry
            from thermoanalysis.constants import J2AU, NA, J2CAL
            from thermoanalysis.config import WORKFLOW_THERMO_POLICY

            qc_data = {
                "coords3d": geometry.cart_coords.reshape(-1, 3) * BOHR2ANG,  # Å
                "wavenumbers": freqs_cm,                                 # cm^-1
                "scf_energy": _calc_energy(geometry, calc_cfg),          # Hartree
                "masses": masses_amu,
                "mult": int(calc_cfg["spin"]),
            }
            qc = QCData(qc_data, point_group="c1", mult=int(calc_cfg["spin"]))
            qc.symmetry_number = int(thermo_cfg["symmetry_number"])

            T = float(thermo_cfg["temperature"])
            p_atm = float(thermo_cfg["pressure_atm"])
            symmetry_number = int(thermo_cfg["symmetry_number"])
            p_pa = p_atm * 101325.0  # Pa

            # The standalone-freq policy is library-default QRRHO with no
            # imaginary inversion and NO positive-frequency floor. Pass every value
            # explicitly and serialize the effective policy below; this reproduces
            # the historical bare thermochemistry(qc, T, pressure=p) numbers.
            _thermo_policy = WORKFLOW_THERMO_POLICY
            tr = thermochemistry(
                qc, T, pressure=p_pa, **_thermo_policy.thermochemistry_kwargs()
            )

            au2CalMol = (1.0 / J2AU) * NA * J2CAL
            to_cal_per_mol = lambda x: float(x) * au2CalMol
            J_per_Kmol_to_cal_per_Kmol = lambda j: float(j) * J2CAL

            n_imag = int(np.sum(freqs_cm < 0.0))

            EE = float(tr.U_el)
            ZPE = float(tr.ZPE)
            dE_therm = float(tr.U_therm)
            dH_therm = float(tr.H - tr.U_el)
            dG_therm = float(tr.dG)

            sum_EE_ZPE = EE + ZPE
            sum_EE_thermal_E = float(tr.U_tot)
            sum_EE_thermal_H = float(tr.H)
            sum_EE_thermal_G = float(tr.G)

            E_thermal_cal = to_cal_per_mol(tr.U_therm)
            Cv_cal_per_Kmol = J_per_Kmol_to_cal_per_Kmol(tr.c_tot)
            S_cal_per_Kmol  = to_cal_per_mol(tr.S_tot)

            emit("\n====== Thermochemistry summary ======\n", narrative=True)
            click.echo(f"Structure               = {input_path}")
            click.echo(f"Temperature (K)         = {T:.2f}")
            click.echo(f"Pressure    (atm)       = {p_atm:.4f}")
            click.echo(
                f"Rotational symmetry no. = {symmetry_number:d} "
                f"({symmetry_number_source})"
            )
            if freeze_list:
                emit("[NOTE] Thermochemistry uses active DOF (PHVA) due to frozen atoms.", narrative=True)
            click.echo(f"Number of Imaginary Freq = {n_imag:d}")
            click.echo("")

            click.echo(f"Electronic Energy (E)                  = {_fmt_ha(EE)}")
            click.echo(f"Zero-point Energy Correction           = {_fmt_ha(ZPE)}")
            click.echo(f"Thermal Correction to Energy           = {_fmt_ha(dE_therm)}")
            click.echo(f"Thermal Correction to Enthalpy         = {_fmt_ha(dH_therm)}")
            click.echo(f"Gibbs Free Energy Correction (G_corr)  = {_fmt_ha(dG_therm)}")
            click.echo(f"EE + Zero-point Energy                 = {_fmt_ha(sum_EE_ZPE)}")
            click.echo(f"EE + Thermal Energy Correction         = {_fmt_ha(sum_EE_thermal_E)}")
            click.echo(f"EE + Thermal Enthalpy Correction       = {_fmt_ha(sum_EE_thermal_H)}")
            click.echo(f"Gibbs Free Energy (G = E + G_corr)      = {_fmt_ha(sum_EE_thermal_G)}")
            click.echo("")
            click.echo(f"E (Thermal)                            = {_fmt_cal(E_thermal_cal)}")
            click.echo(f"Heat Capacity (Cv)                     = {_fmt_calK(Cv_cal_per_Kmol)}")
            click.echo(f"Entropy (S)                            = {_fmt_calK(S_cal_per_Kmol)}")

            if bool(thermo_cfg["dump"]):
                payload = {
                    "structure": str(input_path),
                    "temperature_K": T,
                    "pressure_atm": p_atm,
                    "symmetry_number": symmetry_number,
                    "symmetry_number_source": symmetry_number_source,
                    "num_imag_freq": n_imag,
                    "n_freeze_atoms": len(freeze_list),
                    "thermo_policy": _thermo_policy.as_dict(),
                    "rigid_projection": _rigid_projection,
                    "electronic_energy_ha": EE,
                    "zpe_correction_ha": ZPE,
                    "thermal_correction_energy_ha": dE_therm,
                    "thermal_correction_enthalpy_ha": dH_therm,
                    "thermal_correction_free_energy_ha": dG_therm,
                    "sum_EE_and_ZPE_ha": sum_EE_ZPE,
                    "sum_EE_and_thermal_energy_ha": sum_EE_thermal_E,
                    "sum_EE_and_thermal_enthalpy_ha": sum_EE_thermal_H,
                    "sum_EE_and_thermal_free_energy_ha": sum_EE_thermal_G,
                    "E_thermal_cal_per_mol": E_thermal_cal,
                    "Cv_cal_per_mol_K": Cv_cal_per_Kmol,
                    "S_cal_per_mol_K": S_cal_per_Kmol,
                }
                try:
                    with _thermo_yaml_tmp.open("w", encoding="utf-8") as f:
                        yaml.safe_dump(
                            payload, f, sort_keys=False, allow_unicode=True
                        )
                    _thermo_yaml_tmp.replace(_thermo_yaml)
                finally:
                    _thermo_yaml_tmp.unlink(missing_ok=True)
                emit(
                    f"[dump] Wrote thermoanalysis summary → {_thermo_yaml}",
                    detail=True,
                )

            _thermo_data = {
                "thermo_policy": _thermo_policy.as_dict(),
                "symmetry_number": symmetry_number,
                "symmetry_number_source": symmetry_number_source,
                "electronic_energy_ha": EE,
                "zpe_correction_ha": ZPE,
                "thermal_correction_energy_ha": dE_therm,
                "thermal_correction_enthalpy_ha": dH_therm,
                "thermal_correction_free_energy_ha": dG_therm,
                "sum_EE_and_ZPE_ha": sum_EE_ZPE,
                "sum_EE_and_thermal_energy_ha": sum_EE_thermal_E,
                "sum_EE_and_thermal_enthalpy_ha": sum_EE_thermal_H,
                "sum_EE_and_thermal_free_energy_ha": sum_EE_thermal_G,
                "E_thermal_cal_per_mol": E_thermal_cal,
                "Cv_cal_per_mol_K": Cv_cal_per_Kmol,
                "S_cal_per_mol_K": S_cal_per_Kmol,
            }

        except ImportError as e:
            raise click.ClickException(
                "Thermochemistry failed because 'thermoanalysis' is unavailable."
            ) from e
        except Exception as e:
            raise click.ClickException(f"Thermochemistry failed: {e}") from e

        # summary.md and key_* outputs are disabled.
        emit(f"[DONE] Wrote modes and list → {out_dir_path}", detail=True)

        # result.json (if --out-json)
        if out_json:
            from pdb2reaction.core.utils import calculator_provenance, write_result_json
            _all_freqs = [float(f) for f in freqs_cm]
            _imag_freqs = [f for f in _all_freqs if f < 0.0]
            result_data = {
                "status": "completed",
                "n_modes": len(_all_freqs),
                "n_imaginary": len(_imag_freqs),
                "frequencies_cm": _all_freqs,
                "imaginary_frequencies_cm": _imag_freqs,
                "thermochemistry": _thermo_data,
                "rigid_projection": _rigid_projection,
                "backend": calc_cfg.get("backend", backend),
                "charge": calc_cfg["charge"],
                "spin": calc_cfg["spin"],
                "model": calc_cfg.get("model"),
                **calculator_provenance(calc_cfg),
                "n_atoms": len(geometry.atomic_numbers),
                "n_freeze_atoms": len(freeze_list) if 'freeze_list' in dir() else 0,
                "solvent": calc_cfg.get("solvent", "none"),
                "temperature_K": thermo_cfg.get("temperature", 298.15) if 'thermo_cfg' in dir() else 298.15,
                "pressure_atm": thermo_cfg.get("pressure_atm", 1.0) if 'thermo_cfg' in dir() else 1.0,
                "input_file": str(input_path),
                "files": {
                    "frequencies_txt": "frequencies_cm-1.txt",
                    "mode_files": _mode_output_files,
                },
            }
            if _thermo_data is not None and bool(thermo_cfg.get("dump", False)):
                result_data["files"]["thermoanalysis_yaml"] = "thermoanalysis.yaml"
            write_result_json(
                out_dir_path, result_data,
                command="freq",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        emit(
            format_elapsed("[time] Elapsed Time for Freq", time_start),
            narrative=True,
        )

    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        render_cli_exception(e, label="frequency analysis", out_dir=out_dir_path, command="freq", time_start=time_start)
    finally:
        # Release GPU memory so subsequent pipeline stages don't OOM:
        # rebind to None first (decrements heavy torch.nn.Module refs),
        # then del to drop the local names; gc.collect breaks cycles.
        calc = geometry = H = H_analysis = modes = modes_mw = None
        del calc, geometry, H, H_analysis, modes, modes_mw
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
