# pdb2reaction/dft.py

"""
dft — Single-point DFT calculation
====================================================================

Usage (CLI)
-----------
    pdb2reaction dft -i INPUT.{pdb|xyz|gjf|...} [-q <charge>] [--ligand-charge <number|'RES:Q,...'>] [-m <multiplicity>] \
        [--func-basis 'FUNC/BASIS'] [--max-cycle <int>] [--conv-tol <hartree>] \
        [--grid-level <int>] [--out-dir <dir>] [--engine {gpu|cpu|auto}] \
        [--convert-files {True|False}] [--args-yaml <file>]

Examples
--------
    # Default GPU-first policy with an explicit functional/basis pair
    pdb2reaction dft -i input.pdb -q 0 -m 1 --func-basis 'wb97m-v/6-31g**'

    # Tight SCF controls with a larger basis and CPU-only fallback
    pdb2reaction dft -i input.pdb -q 0 -m 2 --func-basis 'wb97m-v/def2-tzvpd' \
        --max-cycle 150 --conv-tol 1e-9 --engine cpu

Description
-----------
- Single-point DFT engine with optional GPU acceleration (GPU4PySCF) and a CPU PySCF backend.
  The backend policy is controlled by --engine:
  * gpu  (default): try GPU4PySCF first; on import/runtime errors, automatically fall back to CPU PySCF.
                    Blackwell GPUs emit a warning on detection.
  * cpu           : use CPU PySCF only.
  * auto          : try GPU4PySCF first and fall back to CPU PySCF if unavailable (same behavior as "gpu").
- RKS/UKS is selected automatically from the spin multiplicity (2S+1).
- Inputs: any structure format supported by pysisyphus.helpers.geom_loader (.pdb, .xyz, .trj, …).
  The geometry is written back unchanged as input_geometry.xyz.
- Functional/basis specified as 'FUNC/BASIS' via --func-basis (e.g., 'wb97m-v/6-31g**', 'wb97m-v/def2-svp', 'wb97m-v/def2-tzvpd').
  Names are case-insensitive in PySCF.
- Density fitting (DF) is enabled via PySCF's density_fit(); the auxiliary basis is left to
  PySCF's default selection.
- SCF controls: --conv-tol (Eh), --max-cycle, --grid-level (mapped to PySCF grids.level).
  Verbosity defaults to 0 and can be overridden via YAML (dft.verbose). Output directory
  selection is handled separately via --out-dir.
- VV10 / other nonlocal corrections are **not** configured explicitly; backends run with their
  defaults for the chosen functional.
- Charge/spin are resolved by internal helpers; .gjf templates may supply charge/spin when available.
  Provide explicit -q/--charge and -m/--multiplicity values whenever possible to enforce the intended state
  (multiplicity > 1 selects UKS). When ``-q`` is omitted but ``--ligand-charge`` is given for a parseable
  complex structure, the total system charge is inferred using ``extract.py``’s residue-aware logic; an
  explicit ``-q`` always takes precedence.
- **Atomic properties:** from the final density, **atomic charges** and **atomic spin densities** are reported by three schemes:
    * Mulliken (charges: scf.hf.mulliken_pop; spins: scf.uhf.mulliken_spin_pop for UKS; RKS → zeros; failure → null)
    * meta‑Löwdin (charges: scf.hf.mulliken_pop_meta_lowdin_ao; spins: scf.uhf.mulliken_spin_pop_meta_lowdin_ao for UKS; RKS → zeros; not available or failure → null)
    * IAO (charges: lo.iao.fast_iao_mullikan_pop; spins: fast_iao_mullikan_spin_pop implemented here; RKS → zeros; failure → null)
- Energies are reported in Hartree and kcal/mol; SCF convergence status and backend selection are recorded.
  Elapsed time is printed to stdout.
- The per-atom tables are echoed to stdout and saved to YAML in flow-style rows for readability.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_dft/)
  ├─ result.yaml                # Input metadata, SCF energy (Eh/kcal), convergence status, backend used, and per-atom charge/spin tables
  └─ input_geometry.xyz         # Geometry snapshot passed to PySCF (as read; unchanged)

Notes
-----
- Charge/spin resolution: -q/--charge and -m/--multiplicity may inherit values from .gjf templates when present.
  For non-.gjf inputs, omitting -q/--charge is allowed only when ``--ligand-charge`` is supplied; the full complex
  is treated as an enzyme–substrate system and the total charge is derived with the same logic as ``extract.py``.
  Otherwise the CLI aborts. Explicit ``-q`` overrides any derived charge. When a `.gjf` template omits charge/spin,
  defaults are `0` and `1`. Provide explicit values whenever possible to enforce the intended state
  (multiplicity > 1 selects UKS).
- YAML overrides: --args-yaml points to a file with top-level keys "dft" (func, basis, conv_tol,
  max_cycle, grid_level, verbose, out_dir, or combined func_basis) and "geom" (passed to
  pysisyphus.helpers.geom_loader).
- Grids: sets grids.level when supported.
- Units: input coordinates are in Å.
- If any population analysis (Mulliken, meta‑Löwdin, IAO) fails, a WARNING is printed and the corresponding column is null.
"""

from __future__ import annotations

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

from .utils import (
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    format_elapsed,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
)
from .uma_pysis import GEOM_KW_DEFAULT


# -----------------------------------------------
# Defaults (override via CLI / YAML)
# -----------------------------------------------

DFT_DEFAULT_FUNC = "wb97m-v"
DFT_DEFAULT_BASIS = "def2-tzvpd"

DFT_KW: Dict[str, Any] = {
    "conv_tol": 1e-9,          # SCF convergence tolerance (Eh)
    "max_cycle": 100,          # Maximum number of SCF iterations
    "grid_level": 3,           # Numerical integration grid level (PySCF grids.level)
    "verbose": 0,              # PySCF verbosity (0..9)
    "out_dir": "./result_dft/",# Output directory
    "func": DFT_DEFAULT_FUNC,  # XC functional (can be overridden via YAML)
    "basis": DFT_DEFAULT_BASIS,# Basis set (can be overridden via YAML)
}


# -----------------------------------------------
# Utilities
# -----------------------------------------------

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


def _configure_scf_object(mf, dft_cfg: Dict[str, Any], xc: str):
    """Apply common SCF settings (XC, DF, tolerances, grids) to an SCF object."""
    mf.xc = xc
    mf.max_cycle = int(dft_cfg["max_cycle"])
    mf.conv_tol = float(dft_cfg["conv_tol"])
    mf.grids.level = int(dft_cfg["grid_level"])
    mf.chkfile = None
    mf = mf.density_fit()

    return mf


# ---------------- Flow-style YAML helper (only for inner row lists) -----------
class FlowList(list):
    """A list that will be dumped in YAML flow style: [a, b, c]."""
    pass


def _flow_seq_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.SafeDumper.add_representer(FlowList, _flow_seq_representer)


def _format_row_for_echo(row: List[Union[int, str, float, None]]) -> str:
    """Format a row like: [0, H, 0.0, 0.0, 0.0]."""
    def _fmt(x):
        if x is None:
            return "null"
        if isinstance(x, float):
            return f"{x:.10g}"
        return str(x)
    return "[" + ", ".join(_fmt(v) for v in row) + "]"


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

    # return scf_uhf.mulliken_pop(pmol, [dm_a, dm_b], s_iao, verbose)
    # -->
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

    # Mulliken charges
    try:
        _, mull_chg = scf_hf.mulliken_pop(mol, dm_tot, s=S, verbose=0)
        mull_q: Optional[List[float]] = np.asarray(mull_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Mulliken] WARNING: Failed to compute Mulliken charges: {e}", err=True)
        mull_q = None

    # meta-Löwdin charges
    try:
        _, low_chg = scf_hf.mulliken_pop_meta_lowdin_ao(mol, dm_tot, verbose=0, s=S)
        low_q: Optional[List[float]] = np.asarray(low_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Löwdin] WARNING: Failed to compute meta-Löwdin charges: {e}", err=True)
        low_q = None

    # IAO charges
    iao_q: Optional[List[float]]
    try:
        iaos = lo_iao.iao(mol, _get_occupied_orbitals(mf), minao="minao")
        _, iao_chg = lo_iao.fast_iao_mullikan_pop(mol, dm, iaos, verbose=0)
        iao_q = np.asarray(iao_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[IAO] WARNING: Failed to compute IAO charges: {e}", err=True)
        iao_q = None

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

    try:
        _, Ms_mull = scf_uhf.mulliken_spin_pop(mol, dm, s=S, verbose=0)
        mull: Optional[List[float]] = np.asarray(Ms_mull, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin Mulliken] WARNING: Failed to compute Mulliken spin densities: {e}", err=True)
        mull = None

    try:
        _, Ms_low = scf_uhf.mulliken_spin_pop_meta_lowdin_ao(mol, dm, verbose=0, s=S)
        low: Optional[List[float]] = np.asarray(Ms_low, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin Löwdin] WARNING: Failed to compute meta-Löwdin spin densities: {e}", err=True)
        low = None

    iao_ms: Optional[List[float]]
    try:
        iaos = lo_iao.iao(mol, _get_occupied_orbitals(mf), minao="minao")
        _, Ms_iao = fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=0)
        iao_ms = np.asarray(Ms_iao, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin IAO] WARNING: Failed to compute IAO spin densities: {e}", err=True)
        iao_ms = None

    return {"mulliken": mull, "lowdin": low, "iao": iao_ms}


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Single-point DFT using GPU4PySCF (CPU PySCF backend).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, etc.; loaded via pysisyphus.helpers.geom_loader).",
)
@click.option("-q", "--charge", type=int, required=False, help="Charge of the ML region.")
@click.option(
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help="Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) for unknown residues.",
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default="GJF template or 1",
    help="Spin multiplicity (2S+1) for the ML region (inherits from .gjf when available; otherwise defaults to 1).",
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
@click.option("--out-dir", type=str, default=DFT_KW["out_dir"], show_default=True, help="Output directory.")
@click.option(
    "--engine",
    type=click.Choice(["gpu", "cpu", "auto"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="Preferred SCF backend: GPU (GPU4PySCF), CPU, or auto (try GPU then CPU if GPU is unavailable).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help='Optional YAML overrides under key "dft" (func/basis, conv_tol, max_cycle, grid_level, verbose, out_dir).',
)
def cli(
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
    args_yaml: Optional[Path],
) -> None:
    set_convert_file_enabled(convert_files)
    prepared_input = prepare_input_structure(input_path)
    apply_ref_pdb_override(prepared_input, ref_pdb)
    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input,
        charge,
        spin,
        ligand_charge=ligand_charge,
        prefix="[dft]",
    )
    try:
        time_start = time.perf_counter()
        # --------------------------
        # 1) Assemble configuration (defaults ← CLI ← YAML)
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)
        geom_cfg = dict(GEOM_KW_DEFAULT)
        dft_cfg = dict(DFT_KW)

        # CLI overrides
        dft_cfg["conv_tol"] = float(conv_tol)
        dft_cfg["max_cycle"] = int(max_cycle)
        dft_cfg["grid_level"] = int(grid_level)
        dft_cfg["out_dir"] = out_dir
        cli_xc, cli_basis = _parse_func_basis(func_basis)
        dft_cfg["func"] = cli_xc
        dft_cfg["basis"] = cli_basis

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (dft_cfg, (("dft",),)),
            ],
        )

        if "func_basis" in dft_cfg:
            # Allow a combined "FUNC/BASIS" field in YAML for convenience.
            yaml_func, yaml_basis = _parse_func_basis(str(dft_cfg["func_basis"]))
            dft_cfg["func"] = yaml_func
            dft_cfg["basis"] = yaml_basis

        xc = str(dft_cfg.get("func", "")).strip()
        basis = str(dft_cfg.get("basis", "")).strip()
        if not xc or not basis:
            raise click.BadParameter("Functional and basis must be non-empty (set via --func-basis or YAML dft.func/basis)")
        multiplicity = int(spin)
        if multiplicity < 1:
            raise click.BadParameter("Multiplicity (spin) must be >= 1.")
        spin2s = multiplicity - 1  # PySCF expects 2S

        # Echo resolved config
        out_dir_path = Path(dft_cfg["out_dir"]).resolve()
        echo_cfg = {
            "charge": int(charge),
            "multiplicity": multiplicity,
            "spin (PySCF expects 2S)": spin2s,
            "xc": xc,
            "basis": basis,
            "conv_tol": dft_cfg["conv_tol"],
            "max_cycle": dft_cfg["max_cycle"],
            "grid_level": dft_cfg["grid_level"],
            "out_dir": str(out_dir_path),
            "engine": engine,
        }
        click.echo(pretty_block("geom", format_geom_for_echo(geom_cfg)))
        click.echo(pretty_block("dft", echo_cfg))

        # --------------------------
        # 2) Load geometry
        # --------------------------
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

        # --------------------------
        # 3) Build PySCF molecule
        # --------------------------
        try:
            from pyscf import gto
        except Exception as e:
            click.echo(f"ERROR: PySCF import failed: {e}", err=True)
            sys.exit(2)

        mol = gto.Mole()
        mol.verbose = int(dft_cfg.get("verbose", 4))
        mol.build(
            atom=atoms_list,
            unit="Angstrom",
            charge=int(charge),
            spin=int(spin2s),
            basis=basis,
        )

        # --------------------------
        # 4) Activate GPU & build SCF object
        # --------------------------
        engine = (engine or "gpu").strip().lower()
        using_gpu = False
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
            click.echo("[gpu] WARNING: Detected a Blackwell GPU; GPU4PySCF may be unsupported.")
        # --------------------------------------------------


        if engine in ("gpu", "auto"):
            try:
                from gpu4pyscf import dft as gdf
                mf = make_ks(gdf)
                using_gpu = True
                engine_label = "gpu4pyscf"
                mf = _configure_scf_object(mf, dft_cfg, xc)
                e_tot = mf.kernel()

            except Exception as e:
                click.echo(
                    f"[gpu] WARNING: GPU backend unavailable ({e}); falling back to CPU.",
                )
                using_gpu = False
                engine_label = "pyscf(cpu)"
                engine = "cpu"

        if not using_gpu:
            from pyscf import dft as pdft
            mf = make_ks(pdft)
            mf = _configure_scf_object(mf, dft_cfg, xc)
            e_tot = mf.kernel()

        # --------------------------
        # 5) Run SCF
        # --------------------------
        

        converged = bool(getattr(mf, "converged", False))
        if e_tot is None:
            e_tot = float(getattr(mf, "e_tot", np.nan))

        e_h = float(e_tot)
        e_kcal = _AU2KCALPERMOL(e_h)

        # --------------------------
        # 6) Charges (Mulliken / meta-Löwdin-AO / IAO) & Spin densities
        # --------------------------
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
        click.echo("\ncharges [index, element, mulliken, lowdin, iao]:")
        for row in charges_table:
            click.echo(f"- {_format_row_for_echo(row)}")

        click.echo("\nspin_densities [index, element, mulliken, lowdin, iao]:")
        for row in spins_table:
            click.echo(f"- {_format_row_for_echo(row)}")

        # --------------------------
        # 7) Save result.yaml (flow style rows for readability)
        # --------------------------
        charges_rows_flow = [FlowList(r) for r in charges_table]
        spins_rows_flow   = [FlowList(r) for r in spins_table]

        result_yaml = {
            "input": dict(echo_cfg),  # configuration snapshot
            "energy": {
                "hartree": e_h,
                "kcal_per_mol": e_kcal,
                "converged": converged,
                "engine": engine_label,
                "used_gpu": bool(using_gpu),
            },
            # Table-style outputs (flow lists)
            "charges [index, element, mulliken, lowdin, iao]": charges_rows_flow,
            "spin_densities [index, element, mulliken, lowdin, iao]": spins_rows_flow,
        }
        (out_dir_path / "result.yaml").write_text(
            yaml.safe_dump(result_yaml, sort_keys=False, allow_unicode=True)
        )
        click.echo(f"[write] Wrote '{out_dir_path / 'result.yaml'}'.")

        # --------------------------
        # 8) Final print: energies
        # --------------------------
        click.echo(f"\nE_total (Hartree): {e_h:.12f}")
        click.echo(f"E_total (kcal/mol): {e_kcal:.6f}")

        # Exit codes: 0 if converged, 3 otherwise
        if not converged:
            click.echo("WARNING: SCF did not converge to the requested tolerance.", err=True)
            sys.exit(3)

        click.echo(format_elapsed("[time] Elapsed Time for DFT", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except click.ClickException:
        raise
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during DFT single-point:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()
