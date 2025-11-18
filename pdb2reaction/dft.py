# pdb2reaction/dft.py

"""
dft — Single-point DFT (GPU4PySCF with CPU PySCF fallback)
====================================================================

Usage (CLI)
-----
    pdb2reaction dft -i INPUT -q CHARGE [-s SPIN] [--func-basis "FUNC/BASIS"] [--max-cycle N] [--conv-tol Eh] [--grid-level L] [--out-dir OUT_DIR] [--args-yaml YAML]

    # -q/--charge and -s/--spin fall back to .gjf template values when available (otherwise 0/1),
    # but set them explicitly to avoid unphysical states.

Examples:::
    pdb2reaction dft -i input.pdb -q 0 -s 1 --func-basis "wb97m-v/6-31g**"
    pdb2reaction dft -i input.pdb -q 0 -s 2 --func-basis "wb97m-v/def2-tzvpd" --max-cycle 150 --conv-tol 1e-9

Description
-----
- Single-point DFT engine that activates GPU4PySCF when available and otherwise uses CPU PySCF.
- RKS/UKS is selected automatically from the spin multiplicity (2S+1).
- Inputs: any structure format supported by pysisyphus.helpers.geom_loader (.pdb, .xyz, .trj, …). The geometry is written back unchanged as `input_geometry.xyz`.
- Functional/basis specified as "FUNC/BASIS" via --func-basis (e.g., "wb97m-v/6-31g**", "wb97m-v/def2-tzvpd"). Names are case-insensitive in PySCF.
- SCF controls: --conv-tol (Eh), --max-cycle, --grid-level (mapped to PySCF `grids.level`), --out-dir. Verbosity can be overridden via YAML.
- Nonlocal VV10 is enabled automatically when the functional ends with "-v" or contains "vv10".
- **Atomic properties:** from the final density, **atomic charges** and **atomic spin densities** are reported by three schemes:
    * Mulliken (charges: `scf.hf.mulliken_pop`; spins: `scf.uhf.mulliken_spin_pop` for UKS; failure → null)
    * meta‑Löwdin (charges: `scf.hf.mulliken_pop_meta_lowdin_ao`; spins: `scf.uhf.mulliken_spin_pop_meta_lowdin_ao` for UKS; not available or failure → null)
    * IAO (charges: `lo.iao.fast_iao_mullikan_pop`; spins: `fast_iao_mullikan_spin_pop` implemented here; failure → null)
- Energies are reported in Hartree and kcal/mol; SCF convergence metadata and timing are recorded.

Outputs (& Directory Layout)
-----
    OUT_DIR/  (default: ./result_dft/)
    ├── result.yaml
    │     energy:
    │       hartree: <float>
    │       kcal_per_mol: <float>
    │       converged: <bool>
    │       scf_time_sec: <float>
    │       engine: "gpu4pyscf" or "pyscf(cpu)"
    │       used_gpu: <bool>
    │     charges [index, element, mulliken, lowdin, iao]:
    │       - [0, H, <float or null>, <float or null>, <float or null>]
    │       - ...
    │     spin_densities [index, element, mulliken, lowdin, iao]:
    │       - [0, H, <float or null>, <float or null>, <float or null>]
    │       - ...
    └── input_geometry.xyz  # geometry snapshot used for SCF (as read; unchanged)

Notes:
-----
- Charge/spin resolution: `-q/--charge` and `-s/--spin` inherit values from `.gjf` templates when present and otherwise fall back
  to `0`/`1`. Provide explicit values whenever possible to enforce the intended state (multiplicity > 1 selects UKS).
- YAML overrides: --args-yaml points to a file with top-level key "dft" (conv_tol, max_cycle, grid_level, verbose, out_dir).
- Grids and checkpointing: sets `grids.level` when supported; disables SCF checkpoint files when possible.
- Units: input coordinates are in Å. Functional/basis names are PySCF-style and case-insensitive for common sets.
- Exit codes: 0 if SCF converged; 3 if not converged; 2 if PySCF import fails; 1 on unhandled errors; 130 on user interrupt.
- If any population analysis (Mulliken, meta‑Löwdin, IAO) fails, a WARNING is printed and the corresponding column becomes `null`.
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

from .utils import (
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_elapsed,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    maybe_convert_xyz_to_gjf,
    charge_option,
    spin_option,
)


# -----------------------------------------------
# Defaults (override via CLI / YAML)
# -----------------------------------------------

DFT_KW: Dict[str, Any] = {
    "conv_tol": 1e-9,          # SCF convergence tolerance (Eh)
    "max_cycle": 100,          # Maximum number of SCF iterations
    "grid_level": 3,           # Numerical integration grid level (PySCF grids.level)
    "verbose": 4,              # PySCF verbosity (0..9)
    "out_dir": "./result_dft/",# Output directory
}


# -----------------------------------------------
# Utilities
# -----------------------------------------------

HARTREE_TO_KCALMOL = 627.5094740631  # Commonly used Hartree → kcal/mol conversion factor
def _parse_func_basis(s: str) -> Tuple[str, str]:
    """
    Parse "FUNC/BASIS" into (xc, basis).
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
    s = geometry.as_xyz()  # trusted by other tools in this package
    lines = s.splitlines()
    atoms: list[Tuple[str, Tuple[float, float, float]]] = []
    for ln in lines[2:]:
        parts = ln.split()
        if len(parts) >= 4:
            sym = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((sym, (x, y, z)))
    return s, atoms


def _hartree_to_kcalmol(Eh: float) -> float:
    return float(Eh * HARTREE_TO_KCALMOL)


# ---- Aux-basis chooser for DF (JKFIT) ----
def _choose_auxbasis_for_orbital_basis(basis: str) -> Optional[str]:
    """
    Pick a practical JKFIT auxiliary basis for density fitting from the orbital basis name.
    Policy:
      - def2 family:
          * SVP  → def2-SVP-JKFIT
          * TZ*  → def2-TZVPP-JKFIT
          * QZ*  → def2-QZVPP-JKFIT
          * otherwise → def2-universal-jkfit
      - (aug-)cc-pVXZ family: use cc-pVXZ-JKFIT
      - Pople (6-31G**, 6-311G** ...): use "weigend"
      - Unknown → None (delegate to PySCF default)
    """
    if not basis:
        return None
    b = basis.strip().lower()
    # def2 family
    if "def2" in b:
        if "svp" in b and "tz" not in b and "qz" not in b:
            return "def2-svp-jkfit"
        if "tz" in b:
            return "def2-tzvpp-jkfit"
        if "qz" in b:
            return "def2-qzvpp-jkfit"
        return "def2-universal-jkfit"
    # (aug-)cc-pVXZ family
    if "cc-pv" in b:  # matches both cc-pv and aug-cc-pv
        if "dz" in b: return "cc-pvdz-jkfit"
        if "tz" in b: return "cc-pvtz-jkfit"
        if "qz" in b: return "cc-pvqz-jkfit"
        if "5z" in b: return "cc-pv5z-jkfit"
    # Pople family
    if "6-31" in b or "6-311" in b:
        return "weigend"
    return None


# ---------------- Flow-style YAML helper (only for inner row lists) -----------
class FlowList(list):
    """A list that will be dumped in YAML flow style: [a, b, c]."""
    pass

def _flow_seq_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

yaml.SafeDumper.add_representer(FlowList, _flow_seq_representer)


def _format_row_for_echo(row: List[Union[int, str, float, None]]) -> str:
    """Format a row like: [0, H, 0.0, 0.0, 0.0]"""
    def _fmt(x):
        if x is None:
            return "null"
        if isinstance(x, float):
            # keep concise; similar to YAML default float rendering
            return f"{x:.10g}"
        return str(x)
    return "[" + ", ".join(_fmt(v) for v in row) + "]"


def fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=None):
    """
    IAO-basis Mulliken spin population analysis.

    Args:
        mol : Mole or Cell object
        dm  : AO density matrix; for UKS/UHF, a (2, nao, nao) array
        iaos: 2D array of IAO orbitals (orthogonal or non-orthogonal)
        verbose: PySCF logger level (defaults to logger.DEBUG if None)

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


# ---- Small helpers to remove duplication ------------------------------------
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
    help="Single-point DFT using GPU4PySCF (CPU PySCF fallback).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, etc.; loaded via pysisyphus.helpers.geom_loader).",
)
@charge_option()
@spin_option()
@click.option(
    "--func-basis",
    "func_basis",
    type=str,
    default="wb97m-v/6-31g**",
    show_default=True,
    help='Exchange–correlation functional and basis set as "FUNC/BASIS" (e.g., "wb97m-v/6-31g**", "wb97m-v/def2-tzvpd").',
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
    help="Preferred SCF backend: GPU (GPU4PySCF when available), CPU, or auto (try GPU then CPU).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help='Optional YAML overrides under key "dft" (conv_tol, max_cycle, grid_level, verbose, out_dir).',
)

def cli(
    input_path: Path,
    charge: Optional[int],
    spin: Optional[int],
    func_basis: str,
    max_cycle: int,
    conv_tol: float,
    grid_level: int,
    out_dir: str,
    engine: str,
    args_yaml: Optional[Path],
) -> None:
    prepared_input = prepare_input_structure(input_path)
    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)
    try:
        time_start = time.perf_counter()
        # --------------------------
        # 1) Assemble configuration
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)
        dft_cfg = dict(DFT_KW)
        apply_yaml_overrides(yaml_cfg, [(dft_cfg, (("dft",),))])

        # CLI overrides
        dft_cfg["conv_tol"] = float(conv_tol)
        dft_cfg["max_cycle"] = int(max_cycle)
        dft_cfg["grid_level"] = int(grid_level)
        dft_cfg["out_dir"] = out_dir

        xc, basis = _parse_func_basis(func_basis)
        multiplicity = int(spin)
        if multiplicity < 1:
            raise click.BadParameter("Multiplicity (spin) must be >= 1.")
        spin2s = multiplicity - 1  # PySCF expects 2S

        # ---- choose aux-basis for DF from orbital basis ----
        aux_basis_guess = _choose_auxbasis_for_orbital_basis(basis)

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
        click.echo(pretty_block("dft", echo_cfg))

        # --------------------------
        # 2) Load geometry
        # --------------------------
        geometry = geom_loader(geom_input_path, coord_type="cart")
        xyz_s, atoms_list = _geometry_to_pyscf_atoms_string(geometry)

        out_dir_path.mkdir(parents=True, exist_ok=True)
        # Write a provenance snapshot of the input geometry
        input_xyz = out_dir_path / "input_geometry.xyz"
        input_xyz.write_text(xyz_s if xyz_s.endswith("\n") else (xyz_s + "\n"))
        click.echo(f"[write] Wrote '{input_xyz}'.")
        gjf_written = maybe_convert_xyz_to_gjf(input_xyz, prepared_input.gjf_template, out_dir_path / "input_geometry.gjf")
        if gjf_written:
            click.echo(f"[convert] Wrote '{gjf_written}'.")

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
        gpu_exc: Optional[str] = None
        if engine in ("gpu", "auto"):
            try:
                import gpu4pyscf
                gpu4pyscf.activate()  # patch PySCF backends to GPU where supported
                from gpu4pyscf import dft as gdf
                mf = make_ks(gdf)
                using_gpu = True
                engine_label = "gpu4pyscf"
            except Exception as e:
                gpu_exc = str(e)
                if engine == "gpu":
                    click.echo(
                        f"[gpu] WARNING: GPU backend requested but unavailable ({gpu_exc}); falling back to CPU.",
                        err=True,
                    )
        if not using_gpu:
            from pyscf import dft as pdft
            mf = make_ks(pdft)
            if engine == "gpu" and gpu_exc is None:
                click.echo(
                    "[gpu] WARNING: GPU backend requested but unavailable (unknown error); falling back to CPU.",
                    err=True,
                )

        # ---- Enable density fitting (RI/DF) & set aux-basis if available ----
        try:
            mf = mf.density_fit()
            if aux_basis_guess:
                try:
                    mf.with_df.auxbasis = aux_basis_guess
                except Exception as e:
                    click.echo(f"[df] WARNING: Could not set auxbasis='{aux_basis_guess}': {e}", err=True)
        except Exception as e:
            click.echo(f"[df] WARNING: density_fit() failed or not available: {e}", err=True)

        # SCF settings
        mf.xc = xc
        mf.max_cycle = int(dft_cfg["max_cycle"])
        mf.conv_tol = float(dft_cfg["conv_tol"])
        try:
            # grids.level is the standard PySCF knob; supported by GPU4PySCF
            mf.grids.level = int(dft_cfg["grid_level"])
        except Exception as e:
            click.echo(f"[grids] WARNING: Could not set grids.level={dft_cfg['grid_level']}: {e}", err=True)

        # Disable checkpoint file if possible
        try:
            mf.chkfile = None
        except Exception:
            pass

        # Enable VV10 for "-v" functionals
        xcl = xc.lower()
        if xcl.endswith("-v") or "vv10" in xcl:
            # PySCF style for VV10 nonlocal correlation
            mf.nlc = "vv10"

        # --------------------------
        # 5) Run SCF
        # --------------------------
        click.echo("\n=== DFT single-point started (GPU if available) ===\n")
        tic_scf = time.time()
        e_tot = mf.kernel()
        toc_scf = time.time()
        click.echo("\n=== DFT single-point finished ===\n")

        converged = bool(getattr(mf, "converged", False))
        if e_tot is None:
            # Some PySCF versions return None on non-convergence
            e_tot = float(getattr(mf, "e_tot", np.nan))

        e_h = float(e_tot)
        e_kcal = _hartree_to_kcalmol(e_h)

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
            "input": dict(echo_cfg),  # reuse echoed configuration
            "energy": {
                "hartree": e_h,
                "kcal_per_mol": e_kcal,
                "converged": converged,
                "scf_time_sec": round(toc_scf - tic_scf, 3),
                "engine": engine_label,
                "used_gpu": bool(using_gpu),
            },
            # Requested table-style outputs (flow lists)
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

        # Exit codes: 0 if converged, 3 otherwise (for compatibility with prior behavior)
        if not converged:
            click.echo("WARNING: SCF did not converge to the requested tolerance.", err=True)
            sys.exit(3)

        click.echo(format_elapsed("[time] Elapsed Time for DFT", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except click.ClickException:
        # Re-raise click-specific errors
        raise
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during DFT single-point:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()


# Enable `python -m pdb2reaction.dft` execution
if __name__ == "__main__":
    cli()
