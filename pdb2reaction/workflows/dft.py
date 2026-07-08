"""
Single-point DFT calculation with GPU acceleration (GPU4PySCF or CPU PySCF).

Example:
    pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

For detailed documentation, see: docs/dft.md
"""
# DOMAIN_PURE

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, List, Union

import sys
import traceback
import textwrap
import time
from functools import reduce

import click
import yaml
import numpy as np

from pysisyphus.helpers import geom_loader
from pysisyphus.constants import AU2KCALPERMOL

from pdb2reaction.core.utils import (
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    prepared_cli_input,
    set_convert_file_enabled,
    YamlFlowList,
    cli_param_overridden,
)
from pdb2reaction.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, _write_error_json
from pdb2reaction.core.defaults import GEOM_KW_DEFAULT, OUT_DIR_DFT

logger = logging.getLogger(__name__)


# Defaults (override via CLI / YAML)

DFT_DEFAULT_FUNC = "wb97m-v"
# CHEMISTRY-RULE:5 def2 family default basis (= auto-ECP for heavy atoms via def2).
DFT_DEFAULT_BASIS = "def2-tzvpd"

DFT_KW: Dict[str, Any] = {
    "conv_tol": 1e-9,          # SCF convergence tolerance (Eh)
    "max_cycle": 100,          # Maximum number of SCF iterations
    "grid_level": 3,           # Numerical integration grid level (PySCF grids.level)
    "verbose": 0,              # PySCF verbosity (0..9)
    "out_dir": OUT_DIR_DFT,    # Output directory
    "func": DFT_DEFAULT_FUNC,  # XC functional (can be overridden via YAML)
    "basis": DFT_DEFAULT_BASIS,# Basis set (can be overridden via YAML)
    # CHEMISTRY-RULE:4 gpu4pyscf rks_lowmem (closed-shell GPU only) triple-guard.
    "lowmem": True,            # Use gpu4pyscf rks_lowmem (closed-shell GPU only); auto-fallback otherwise
}



def _parse_func_basis(s: str) -> Tuple[str, str]:
    """
    Parse 'FUNC/BASIS' into (xc, basis).
    Mixed case is accepted (PySCF is case-insensitive for common names).
    """
    if not s or "/" not in s:
        raise click.BadParameter("Expected 'FUNC/BASIS' (e.g., 'wb97m-v/def2-tzvpd').")
    func, basis = s.split("/", 1)
    func = func.strip()
    basis = basis.strip()
    if not func or not basis:
        raise click.BadParameter("Functional or basis is empty. Example: --func-basis 'wb97m-v/6-31g**'")
    return func, basis


def _geometry_to_pyscf_atoms_string(geometry) -> Tuple[str, Sequence[Tuple[str, Tuple[float, float, float]]]]:
    """
    Convert a pysisyphus Geometry to (xyz_string, PySCF atom list).
    The atom list is [(symbol, (x, y, z)), ...] in Angstrom.
    """
    s = geometry.as_xyz()
    lines = s.splitlines()
    atoms: list[Tuple[str, Tuple[float, float, float]]] = []
    for ln in lines[2:]:
        parts = ln.split()
        if len(parts) >= 4:
            sym = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((sym, (x, y, z)))
    return s, atoms


def _AU2KCALPERMOL(Eh: float) -> float:
    return float(Eh * AU2KCALPERMOL)


def _configure_scf_object(mf, dft_cfg: Dict[str, Any], xc: str, *, use_density_fit: bool = True):
    """Apply common SCF settings (XC, DF, tolerances, grids) to an SCF object.

    `use_density_fit=False` is required for `gpu4pyscf.dft.rks_lowmem.RKS`,
    which intentionally does not implement `density_fit`.
    """
    mf.xc = xc
    mf.max_cycle = int(dft_cfg["max_cycle"])
    mf.conv_tol = float(dft_cfg["conv_tol"])
    mf.grids.level = int(dft_cfg["grid_level"])
    mf.chkfile = None
    if use_density_fit:
        mf = mf.density_fit()

    return mf


def _build_cpu_surrogate_for_analysis(mol, mf, xc: str, spin2s: int):
    """Construct a CPU `pyscf.dft` SCF object carrying converged orbitals from `mf`.

    Used when the GPU SCF object's `to_cpu()` is not implemented
    (e.g. `gpu4pyscf.dft.rks_lowmem.RKS`). Population analysis routines in
    `pyscf.scf.hf` / `pyscf.scf.uhf` operate on numpy arrays, so MOs are
    transferred from CuPy if needed.
    """
    from pyscf import dft as pdft

    def _to_np(x):
        try:
            import cupy as _cp
            if isinstance(x, _cp.ndarray):
                return _cp.asnumpy(x)
        except Exception:
            pass
        return np.asarray(x)

    surrogate = pdft.RKS(mol) if spin2s == 0 else pdft.UKS(mol)
    surrogate.xc = xc
    surrogate.chkfile = None
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    mo_energy = getattr(mf, "mo_energy", None)
    if spin2s == 0:
        surrogate.mo_coeff = _to_np(mo_coeff)
        surrogate.mo_occ = _to_np(mo_occ)
        if mo_energy is not None:
            surrogate.mo_energy = _to_np(mo_energy)
    else:
        surrogate.mo_coeff = (_to_np(mo_coeff[0]), _to_np(mo_coeff[1]))
        surrogate.mo_occ = (_to_np(mo_occ[0]), _to_np(mo_occ[1]))
        if mo_energy is not None:
            surrogate.mo_energy = (_to_np(mo_energy[0]), _to_np(mo_energy[1]))
    surrogate.converged = bool(getattr(mf, "converged", False))
    surrogate.e_tot = float(getattr(mf, "e_tot", float("nan")))
    return surrogate


def _format_row_for_echo(row: List[Union[int, str, float, None]]) -> str:
    """Format a row like: [0, H, 0.0, 0.0, 0.0]."""
    return "[" + ", ".join(_format_row_value_for_echo(v) for v in row) + "]"


def _format_row_value_for_echo(x: Union[int, str, float, None]) -> str:
    if x is None:
        return "null"
    if isinstance(x, float):
        return f"{x:.10g}"
    return str(x)


# This function is based on https://pyscf.org/_modules/pyscf/lo/iao.html
def fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=None):
    """
    IAO-basis Mulliken spin population analysis.

    Args:
        mol : Mole or Cell object
        dm  : AO density matrix; for UKS/UHF, a (2, nao, nao) array
        iaos: 2D array of IAO orbitals (orthogonal or non-orthogonal)
        verbose: PySCF logger level

    Returns:
        (spin_pop_ao, Ms_by_atom)
        spin_pop_ao : Mulliken spin population on each IAO
        Ms_by_atom  : Mulliken spin density per atom (sum over IAOs on atom)
    """
    import numpy
    from pyscf.lib import logger as pyscf_logger
    from pyscf.lo.iao import reference_mol
    from pyscf.scf import uhf as scf_uhf

    if verbose is None:
        verbose = pyscf_logger.DEBUG

    pmol = reference_mol(mol)
    # Overlap in the large basis
    if getattr(mol, 'pbc_intor', None):  # cell?
        ovlpS = mol.pbc_intor('int1e_ovlp')
    else:
        ovlpS = mol.intor_symmetric('int1e_ovlp')

    # Transform DM in big basis to IAO basis
    # |IAO> = |big> C
    # DM_IAO = C^{-1} DM (C^{-1})^T = S_IAO^{-1} C^T S DM S C S_IAO^{-1}
    cs = numpy.dot(iaos.T.conj(), ovlpS)
    s_iao = numpy.dot(cs, iaos)
    iao_inv = numpy.linalg.solve(s_iao, cs)

    # Restricted case: spin density is identically zero
    if isinstance(dm, numpy.ndarray) and dm.ndim == 2:
        spin_pop_ao = numpy.zeros(s_iao.shape[0], dtype=float)
        Ms = numpy.zeros(pmol.natm, dtype=float)
        return spin_pop_ao, Ms

    # Unrestricted: transform alpha/beta DM to IAO basis
    dm_a = reduce(numpy.dot, (iao_inv, dm[0], iao_inv.conj().T))
    dm_b = reduce(numpy.dot, (iao_inv, dm[1], iao_inv.conj().T))

    return scf_uhf.mulliken_spin_pop(pmol, [dm_a, dm_b], s_iao, verbose)


def _get_occupied_orbitals(mf) -> np.ndarray:
    """Return occupied MO coefficients (AO→MO) for RKS or UKS/UHF."""
    mo = mf.mo_coeff
    mo_occ = mf.mo_occ
    if isinstance(mo, np.ndarray) and mo.ndim == 2:
        occ_idx = np.asarray(mo_occ) > 0
        return mo[:, occ_idx]
    else:
        occ_idx = np.asarray(mo_occ[0]) > 0
        return mo[0][:, occ_idx]


def _safe_array(label: str, what: str, func):
    try:
        vals = func()
        return np.asarray(vals, dtype=float).tolist()
    except Exception as e:
        click.echo(f"{label} WARNING: Failed to compute {what}: {e}", err=True)
        return None


def _compute_atomic_charges(mol, mf) -> Dict[str, Optional[List[float]]]:
    """
    Compute atomic charges by three schemes:
      - 'mulliken' : scf.hf.mulliken_pop (failure → None)
      - 'lowdin'   : scf.hf.mulliken_pop_meta_lowdin_ao (failure → None)
      - 'iao'      : lo.iao.fast_iao_mullikan_pop (failure → None)
    """
    from pyscf.scf import hf as scf_hf
    from pyscf.lo import iao as lo_iao

    dm = mf.make_rdm1()
    S = mf.get_ovlp()
    # Total density (for charges)
    dm_tot = dm[0] + dm[1] if (isinstance(dm, np.ndarray) and dm.ndim == 3) else dm

    mull_q: Optional[List[float]] = _safe_array(
        "[Mulliken]",
        "charges",
        lambda: scf_hf.mulliken_pop(mol, dm_tot, s=S, verbose=0)[1],
    )
    low_q: Optional[List[float]] = _safe_array(
        "[Löwdin]",
        "charges",
        lambda: scf_hf.mulliken_pop_meta_lowdin_ao(mol, dm_tot, verbose=0, s=S)[1],
    )
    iao_q: Optional[List[float]] = _safe_array(
        "[IAO]",
        "charges",
        lambda: lo_iao.fast_iao_mullikan_pop(
            mol,
            dm,
            lo_iao.iao(mol, _get_occupied_orbitals(mf), minao="minao"),
            verbose=0,
        )[1],
    )

    return {
        "mulliken": mull_q,
        "lowdin": low_q,
        "iao": iao_q,
    }


def _compute_atomic_spin_densities(mol, mf) -> Dict[str, Optional[List[float]]]:
    """
    Compute atomic spin densities (Ms) by three schemes:
      - 'mulliken' : scf.uhf.mulliken_spin_pop (RKS → zeros; failure → None)
      - 'lowdin'   : scf.uhf.mulliken_spin_pop_meta_lowdin_ao (RKS → zeros; failure → None)
      - 'iao'      : fast_iao_mullikan_spin_pop (RKS → zeros; failure → None)
    """
    from pyscf.scf import uhf as scf_uhf
    from pyscf.lo import iao as lo_iao

    dm = mf.make_rdm1()
    S = mf.get_ovlp()
    nat = mol.natm

    # RKS (restricted) → spin densities are zero
    if not (isinstance(dm, np.ndarray) and dm.ndim == 3):
        zeros = [0.0] * nat
        return {"mulliken": zeros, "lowdin": zeros, "iao": zeros}

    mull: Optional[List[float]] = _safe_array(
        "[Spin Mulliken]",
        "spin densities",
        lambda: scf_uhf.mulliken_spin_pop(mol, dm, s=S, verbose=0)[1],
    )
    low: Optional[List[float]] = _safe_array(
        "[Spin Löwdin]",
        "spin densities",
        lambda: scf_uhf.mulliken_spin_pop_meta_lowdin_ao(mol, dm, verbose=0, s=S)[1],
    )
    iao_ms: Optional[List[float]] = _safe_array(
        "[Spin IAO]",
        "spin densities",
        lambda: fast_iao_mullikan_spin_pop(
            mol,
            dm,
            lo_iao.iao(mol, _get_occupied_orbitals(mf), minao="minao"),
            verbose=0,
        )[1],
    )

    return {"mulliken": mull, "lowdin": low, "iao": iao_ms}



@click.command(
    help="Single-point DFT using GPU4PySCF (CPU PySCF backend).",
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
@click.option("-m", "--multiplicity", "spin", type=int, default=None, show_default=False, help="Spin multiplicity (2S+1; inherits from .gjf when available; otherwise defaults to 1).")
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Accepted for interface consistency; dft does not emit PDB/GJF outputs.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates).",
)
@click.option(
    "--func-basis",
    "func_basis",
    type=str,
    default=f"{DFT_DEFAULT_FUNC}/{DFT_DEFAULT_BASIS}",
    show_default=True,
    help="Exchange–correlation functional and basis set as 'FUNC/BASIS' (e.g., 'wb97m-v/6-31g**', 'wb97m-v/def2-tzvpd').",
)
@click.option("--max-cycle", type=int, default=DFT_KW["max_cycle"], show_default=True, help="Maximum SCF iterations.")
@click.option("--conv-tol", type=float, default=DFT_KW["conv_tol"], show_default=True, help="SCF convergence tolerance (Eh).")
@click.option("--grid-level", type=int, default=DFT_KW["grid_level"], show_default=True, help="Numerical integration grid level (PySCF grids.level).")
@click.option("-o", "--out-dir", type=str, default=DFT_KW["out_dir"], show_default=True, help="Output directory.")
@click.option(
    "--engine",
    type=click.Choice(["gpu", "cpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="SCF backend: gpu (GPU4PySCF, raises error if unavailable) or cpu (PySCF).",
)
@click.option(
    "--lowmem/--no-lowmem",
    "lowmem",
    default=DFT_KW["lowmem"],
    show_default=True,
    help="Use gpu4pyscf rks_lowmem.RKS for closed-shell GPU runs "
         "(memory-efficient direct JK; mutually exclusive with density fitting). "
         "Open-shell, CPU, or pre-rks_lowmem GPU4PySCF installs auto-fall back "
         "to standard RKS/UKS with density fitting.",
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
    help="Validate options and print the execution plan without running DFT.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    convert_files: bool,
    ref_pdb: Optional[Path],
    func_basis: str,
    max_cycle: int,
    conv_tol: float,
    grid_level: int,
    out_dir: str,
    engine: str,
    lowmem: bool,
    config_yaml: Optional[Path],
    show_config: bool,
    out_json: bool,
    dry_run: bool,
) -> None:
    config_yaml, override_yaml, _ = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, _, _ = load_merged_yaml_cfg(
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
        prefix="[dft]",
    ) as (prepared_input, resolved_charge, resolved_spin):
        geom_input_path = prepared_input.geom_path
        out_dir_path = Path(out_dir).resolve()
        try:
            time_start = time.perf_counter()
            geom_cfg = dict(GEOM_KW_DEFAULT)
            dft_cfg = dict(DFT_KW)

            apply_yaml_overrides(
                merged_yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (dft_cfg, (("dft",),)),
                ],
            )

            if cli_param_overridden(ctx, "conv_tol"):
                dft_cfg["conv_tol"] = float(conv_tol)
            if cli_param_overridden(ctx, "max_cycle"):
                dft_cfg["max_cycle"] = int(max_cycle)
            if cli_param_overridden(ctx, "grid_level"):
                dft_cfg["grid_level"] = int(grid_level)
            if cli_param_overridden(ctx, "out_dir"):
                dft_cfg["out_dir"] = out_dir
            if cli_param_overridden(ctx, "lowmem"):
                dft_cfg["lowmem"] = bool(lowmem)

            # Combined "FUNC/BASIS" only overrides the split-form dft.func /
            # dft.basis when it was actually supplied (CLI --func-basis or
            # config dft.func_basis); otherwise preserve YAML split-form values
            # (or the DFT_KW defaults, which always seed func/basis).
            if cli_param_overridden(ctx, "func_basis"):
                func_basis_value = func_basis
            else:
                func_basis_value = dft_cfg.get("func_basis")
            if func_basis_value:
                cfg_func, cfg_basis = _parse_func_basis(func_basis_value)
                dft_cfg["func"] = cfg_func
                dft_cfg["basis"] = cfg_basis

            xc = str(dft_cfg.get("func", "")).strip()
            basis = str(dft_cfg.get("basis", "")).strip()
            if not xc or not basis:
                raise click.BadParameter("Functional and basis must be non-empty (set via --func-basis or YAML dft.func/basis)")
            multiplicity = int(resolved_spin)
            if multiplicity < 1:
                raise click.BadParameter("Multiplicity (spin) must be >= 1.")
            spin2s = multiplicity - 1  # PySCF expects 2S

            # Echo resolved config
            engine_name = str(dft_cfg.get("engine", engine if engine else "gpu")).strip().lower()
            if cli_param_overridden(ctx, "engine"):
                engine_name = (engine or "gpu").strip().lower()
            out_dir_path = Path(dft_cfg["out_dir"]).resolve()
            lowmem_requested = bool(dft_cfg.get("lowmem", True))
            echo_cfg = {
                "charge": int(resolved_charge),
                "multiplicity": multiplicity,
                "spin (PySCF expects 2S)": spin2s,
                "xc": xc,
                "basis": basis,
                "conv_tol": dft_cfg["conv_tol"],
                "max_cycle": dft_cfg["max_cycle"],
                "grid_level": dft_cfg["grid_level"],
                "out_dir": str(out_dir_path),
                "engine": engine_name,
                "lowmem": lowmem_requested,
            }
            click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
            click.echo(pretty_block("dft", echo_cfg))
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
                            "engine": engine_name,
                            "convert_files": bool(convert_files),
                            "will_run_scf": True,
                            "will_write_result_yaml": True,
                            "will_run_population_analysis": True,
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. DFT execution was skipped.")
                return

            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)
            geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)
            xyz_s, atoms_list = _geometry_to_pyscf_atoms_string(geometry)

            out_dir_path.mkdir(parents=True, exist_ok=True)
            # Write a provenance snapshot of the input geometry
            input_xyz = out_dir_path / "input_geometry.xyz"
            input_xyz.write_text(xyz_s if xyz_s.endswith("\n") else (xyz_s + "\n"))
            click.echo(f"[write] Wrote '{input_xyz}'.")

            try:
                from pyscf import gto
            except (ModuleNotFoundError, ImportError) as e:
                click.echo(f"ERROR: PySCF import failed: {e}", err=True)
                sys.exit(2)

            from pdb2reaction.core.utils import is_verbose
            mol = gto.Mole()
            # PySCF verbose level: 0 silent / 2 warnings / 3 info / 4 prints the
            # full `[INPUT]` per-atom coordinate dump + SCF iteration banner.
            # Default runs stay quiet (DFT_KW seeds a low value); `-v` raises the
            # level to >=4 to restore the dump for debugging (~600 lines saved
            # per pipeline).
            mol.verbose = int(dft_cfg.pop("verbose", 0))
            if is_verbose():
                mol.verbose = max(mol.verbose, 4)
            # DO NOT INLINE: omitting this auto-attach gives silent wrong DFT energies for Z>=21; PySCF basis= alone does NOT imply ecp= even for def2 family.
            # def2 family includes Stuttgart ECPs for heavy elements (Z>=21);
            # must set mol.ecp explicitly or PySCF uses all-electron treatment.
            _ecp = dft_cfg.get("ecp", None)
            if _ecp is None and basis.lower().startswith("def2"):
                _ecp = basis
            _build_kw: Dict[str, Any] = dict(
                atom=atoms_list,
                unit="Angstrom",
                charge=int(resolved_charge),
                spin=int(spin2s),
                basis=basis,
            )
            if _ecp:
                _build_kw["ecp"] = _ecp
                click.echo(f"[dft] Using ECP: {_ecp}")
            mol.build(**_build_kw)

            engine = engine_name
            using_gpu = False
            using_lowmem = False
            engine_label = "pyscf(cpu)"
            make_ks = (lambda mod: mod.RKS(mol) if spin2s == 0 else mod.UKS(mol))


            # --- Detect Blackwell GPU and emit warning ---
            is_blackwell_gpu = False
            try:
                import cupy as cp
                dev_id = cp.cuda.runtime.getDevice()
                props = cp.cuda.runtime.getDeviceProperties(dev_id)
                name = props["name"]
                if isinstance(name, bytes):
                    name = name.decode()
                if ("rtx 50" in name.lower()) or ("nvidia b" in name.lower()):
                    is_blackwell_gpu = True
            except Exception:
                is_blackwell_gpu = False

            if is_blackwell_gpu:
                click.echo("[gpu] WARNING: Detected a Blackwell GPU; GPU4PySCF may be unsupported.", err=True)


            if engine == "gpu":
                try:
                    import cupy as _cp
                    _dev_id = _cp.cuda.runtime.getDevice()
                    _props = _cp.cuda.runtime.getDeviceProperties(_dev_id)
                    _dev_name = _props["name"]
                    if isinstance(_dev_name, bytes):
                        _dev_name = _dev_name.decode()
                    click.echo(f"[dft] Using GPU device {_dev_id}: {_dev_name}")
                except Exception:
                    pass

                try:
                    from gpu4pyscf import dft as gdf

                    # DO NOT INLINE: triple-guard (feature flag + closed-shell + ImportError) keeps the CLI working across the wide gpu4pyscf install matrix; rks_lowmem only exists in recent gpu4pyscf and is RKS-only.
                    rks_lowmem_mod = None
                    if lowmem_requested and spin2s == 0:
                        try:
                            from gpu4pyscf.dft import rks_lowmem as rks_lowmem_mod  # type: ignore
                        except ImportError:
                            click.echo(
                                "[lowmem] WARNING: gpu4pyscf.dft.rks_lowmem is not available "
                                "in this gpu4pyscf install; falling back to standard RKS.",
                                err=True,
                            )
                            rks_lowmem_mod = None

                    if rks_lowmem_mod is not None:
                        mf = rks_lowmem_mod.RKS(mol, xc=xc)
                        using_gpu = True
                        using_lowmem = True
                        engine_label = "gpu4pyscf(rks_lowmem)"
                        # density_fit is intentionally NotImplemented in rks_lowmem.
                        mf = _configure_scf_object(mf, dft_cfg, xc, use_density_fit=False)
                    else:
                        if lowmem_requested and spin2s != 0:
                            click.echo(
                                "[lowmem] NOTE: rks_lowmem only supports closed-shell; "
                                "open-shell run uses standard UKS.",
                                err=True,
                            )
                        mf = make_ks(gdf)
                        using_gpu = True
                        engine_label = "gpu4pyscf"
                        mf = _configure_scf_object(mf, dft_cfg, xc)
                    e_tot = mf.kernel()

                except Exception as e:
                    raise click.ClickException(
                        f"[gpu] GPU backend failed: {e}. "
                        "Use --engine cpu to explicitly run on CPU."
                    )

            if engine == "cpu":
                from pyscf import lib as _pyscf_lib
                click.echo(f"[dft] PySCF is using {_pyscf_lib.num_threads()} threads on CPU.")

                from pyscf import dft as pdft
                mf = make_ks(pdft)
                mf = _configure_scf_object(mf, dft_cfg, xc)
                e_tot = mf.kernel()



            converged = bool(getattr(mf, "converged", False))
            if e_tot is None:
                e_tot = float(getattr(mf, "e_tot", np.nan))

            e_h = float(e_tot)
            e_kcal = _AU2KCALPERMOL(e_h)

            if using_lowmem:
                # rks_lowmem.RKS.to_cpu() raises NotImplementedError; build a CPU
                # surrogate carrying the converged orbitals for population analysis.
                mf_for_analysis = _build_cpu_surrogate_for_analysis(mol, mf, xc, spin2s)
            else:
                try:
                    mf_for_analysis = mf.to_cpu()  # GPU → CPU (no-op on CPU backend)
                except Exception:
                    mf_for_analysis = mf

            charges = _compute_atomic_charges(mol, mf_for_analysis)
            spins   = _compute_atomic_spin_densities(mol, mf_for_analysis)

            def _round_list(xs, tol=1e-10):
                return [0.0 if (x == x) and abs(x) < tol else float(x) for x in xs]  # keep NaN as-is

            # Round tiny numbers (None-safe)
            for dct in (charges, spins):
                for key in ("mulliken", "lowdin", "iao"):
                    dct[key] = None if dct[key] is None else _round_list(dct[key])

            # Build per-atom tables
            charges_table: List[List[Any]] = []
            spins_table:   List[List[Any]] = []
            for i in range(mol.natm):
                elem = mol.atom_symbol(i)
                q_mull = None if charges["mulliken"] is None else charges["mulliken"][i]
                q_low  = None if charges["lowdin"]   is None else charges["lowdin"][i]
                q_iao  = None if charges["iao"]      is None else charges["iao"][i]
                charges_table.append([i, elem, q_mull, q_low, q_iao])

                s_mull = None if spins["mulliken"] is None else spins["mulliken"][i]
                s_low  = None if spins["lowdin"]   is None else spins["lowdin"][i]
                s_iao  = None if spins["iao"]      is None else spins["iao"][i]
                spins_table.append([i, elem, s_mull, s_low, s_iao])

            # ---- Echo charges/spins to stdout in flow style lines ----
            click.echo("")
            click.echo("charges [index, element, mulliken, lowdin, iao]:")
            for row in charges_table:
                click.echo(f"- {_format_row_for_echo(row)}")

            # Skip the per-atom spin_densities echo for closed-shell systems —
            # spin=0 forces every row to all-zero and the resulting N zero rows
            # contribute nothing but noise. The data still lands in result.yaml.
            if spin2s != 0:
                click.echo("")
                click.echo("spin_densities [index, element, mulliken, lowdin, iao]:")
                for row in spins_table:
                    click.echo(f"- {_format_row_for_echo(row)}")

            charges_rows_flow = [YamlFlowList(r) for r in charges_table]
            spins_rows_flow = [YamlFlowList(r) for r in spins_table]

            result_yaml = {
                "input": dict(echo_cfg),  # configuration snapshot
                "energy": {
                    "hartree": e_h,
                    "kcal_per_mol": e_kcal,
                    "converged": converged,
                    "engine": engine_label,
                    "used_gpu": bool(using_gpu),
                    "used_lowmem": bool(using_lowmem),
                },
                # Table-style outputs (flow lists)
                "charges [index, element, mulliken, lowdin, iao]": charges_rows_flow,
                "spin_densities [index, element, mulliken, lowdin, iao]": spins_rows_flow,
            }
            (out_dir_path / "result.yaml").write_text(
                yaml.safe_dump(result_yaml, sort_keys=False, allow_unicode=True)
            )
            click.echo(f"[write] Wrote '{out_dir_path / 'result.yaml'}'.", detail=True)
            # summary.md and key_* outputs are disabled.
            click.echo("")
            click.echo(f"E_total (Hartree): {e_h:.12f}")
            click.echo(f"E_total (kcal/mol): {e_kcal:.6f}")

            # Exit codes: 0 if converged, 3 otherwise
            if not converged:
                click.echo("WARNING: SCF did not converge to the requested tolerance.", err=True)
                sys.exit(3)

            click.echo(format_elapsed("[time] Elapsed Time for DFT", time_start), narrative=True)

            # result.json (if --out-json)
            if out_json:
                from pdb2reaction.core.utils import write_result_json
                result_data: Dict[str, Any] = {
                    "converged": converged,
                    "charge": resolved_charge,
                    "spin": multiplicity,
                    "n_atoms": mol.natm,
                    "grid_level": dft_cfg.get("grid_level"),
                    "conv_tol": dft_cfg.get("conv_tol"),
                    "max_cycle": dft_cfg.get("max_cycle"),
                    "input_file": str(input_path),
                    "energy_hartree": e_h,
                    "energy_kcal_per_mol": e_kcal,
                    "xc_functional": xc,
                    "basis_set": basis,
                    "engine": engine_label,
                    "used_gpu": bool(using_gpu),
                    "used_lowmem": bool(using_lowmem),
                    "charges": {k: v for k, v in charges.items()},
                    "spin_densities": {k: v for k, v in spins.items()},
                    "files": {
                        "result_yaml": "result.yaml",
                        "input_geometry_xyz": "input_geometry.xyz",
                    },
                }
                write_result_json(
                    out_dir_path, result_data,
                    command="dft",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

        except KeyboardInterrupt:
            click.echo("Interrupted by user.", err=True)
            sys.exit(130)
        except click.ClickException:
            raise
        except Exception as e:
            _write_error_json(out_dir_path, "dft", e, "UnhandledError", time_start)
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            click.echo("Unhandled error during DFT single-point:\n" + textwrap.indent(tb, "  "), err=True)
            sys.exit(1)
