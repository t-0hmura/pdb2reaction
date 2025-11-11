# pdb2reaction/dft.py

"""
Single-point DFT using GPU4PySCF (falls back to CPU PySCF).

Features
--------
- Input formats: any format supported by pysisyphus.helpers.geom_loader (.pdb, .xyz, .trj, …).
- Engine: GPU4PySCF + PySCF; RKS/UKS is chosen automatically from the spin multiplicity.
- Functional / basis: pass as "FUNC/BASIS" via --func-basis (e.g., "wb97m-v/6-31g**", "wb97m-v/def2-tzvpd").
- Outputs (in out_dir):
    - result.yaml        : total energy (Hartree & kcal/mol), SCF metadata, and atomic charges (Mulliken / Löwdin / IAO)
    - input_geometry.xyz : the geometry snapshot used for SCF (as read; unchanged)

Examples
--------
pdb2reaction dft -i input.pdb -q 0 -s 1 --func-basis "wb97m-v/6-31g**"
pdb2reaction dft -i input.pdb -q 0 -s 2 --func-basis "wb97m-v/def2-tzvpd" --max-cycle 150 --conv-tol 1e-9
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, List

import sys
import traceback
import textwrap
import time

import click
import yaml
import numpy as np

from pysisyphus.helpers import geom_loader

from .utils import deep_update, load_yaml_dict


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


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


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


def _maybe_enable_vv10(mf, xc: str) -> None:
    """
    Attempt to enable VV10 for "-v" functionals (e.g., wb97m-v).
    If the GPU backend does not support nlc, print a warning and continue.
    """
    xcl = xc.lower()
    if xcl.endswith("-v") or "vv10" in xcl:
        try:
            # PySCF style for VV10 nonlocal correlation
            mf.nlc = "vv10"
            # Use library defaults for parameters; users can override via args-yaml if needed
        except Exception as e:
            click.echo(f"[vv10] WARNING: Could not enable VV10 nonlocal correlation on the GPU backend: {e}", err=True)


def _compute_atomic_charges(mol, mf) -> Dict[str, Optional[List[float]]]:
    """
    Compute atomic charges:
      - 'mulliken': Mulliken charges per atom
      - 'lowdin'  : Löwdin charges per atom
      - 'iao'     : IAO charges per atom (None if IAO construction fails)
    """
    # Total density matrix in the AO basis
    dm = mf.make_rdm1()
    if isinstance(dm, np.ndarray) and dm.ndim == 3:
        # UKS: dm has shape (2, nao, nao); sum α and β
        dm_tot = dm[0] + dm[1]
    else:
        dm_tot = dm

    Z = mol.atom_charges()
    aoslice = mol.aoslice_by_atom()  # (natm, 4): [ish0, ish1, iao0, iao1]
    S = mol.intor_symmetric("int1e_ovlp")

    # ---------- Mulliken ----------
    # AO populations = diag(D * S)
    pop_ao_mull = np.einsum("ij,ji->i", dm_tot, S)
    mull = []
    for _, _, p0, p1 in aoslice:
        pop_atom = float(np.sum(pop_ao_mull[p0:p1]))
        mull.append(float(pop_atom))
    mull_q = (Z - np.array(mull)).astype(float).tolist()

    # ---------- Löwdin ----------
    # AO populations = diag(S^{1/2} * D * S^{1/2})
    w, U = np.linalg.eigh(S)
    w = np.clip(w, 0.0, None)
    S_half = (U * np.sqrt(w)) @ U.T
    pop_ao_low = np.diag(S_half @ dm_tot @ S_half)
    low = []
    for _, _, p0, p1 in aoslice:
        pop_atom = float(np.sum(pop_ao_low[p0:p1]))
        low.append(float(pop_atom))
    low_q = (Z - np.array(low)).astype(float).tolist()

    # ---------- IAO ----------
    iao_q: Optional[List[float]] = None
    try:
        from pyscf.lo import iao, orth
        # For UKS, use the alpha coefficients as the reference to build IAOs
        mo_coeff = mf.mo_coeff
        if isinstance(mo_coeff, (list, tuple)) or (hasattr(mo_coeff, "ndim") and getattr(mo_coeff, "ndim", 0) == 3):
            mo_for_iao = mo_coeff[0]
        else:
            mo_for_iao = mo_coeff
        C_iao = iao.iao(mol, mo_for_iao, minao="minao")
        C_iao = orth.vec_lowdin(C_iao, S)  # orthonormalize IAOs in the AO metric

        # Population in the IAO basis: diag(C^T D C)
        D_iao = C_iao.T @ dm_tot @ C_iao
        pop_iao_orb = np.diag(D_iao)

        # Assign each IAO to an atom by the largest AO-squared weight on that atom
        nat = mol.natm
        ao2atom = np.empty(S.shape[0], dtype=int)
        for a, (_, __, p0, p1) in enumerate(aoslice):
            ao2atom[p0:p1] = a

        iao_atom = np.zeros(nat, dtype=float)
        for j in range(C_iao.shape[1]):
            col = C_iao[:, j]
            w_per_atom = np.zeros(nat, dtype=float)
            # AO-squared weights by atom
            for mu, c in enumerate(col):
                w_per_atom[ao2atom[mu]] += float(c * c)
            amax = int(np.argmax(w_per_atom))
            iao_atom[amax] += float(pop_iao_orb[j])

        iao_q = (Z - iao_atom).astype(float).tolist()
    except Exception as e:
        click.echo(f"[IAO] WARNING: Failed to compute IAO charges: {e}", err=True)
        iao_q = None

    return {
        "mulliken": mull_q,
        "lowdin": low_q,
        "iao": iao_q,
    }


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
@click.option("-q", "--charge", type=int, required=True, help="Total charge.")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1).")
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
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help='Optional YAML overrides under key "dft" (conv_tol, max_cycle, grid_level, verbose, out_dir).',
)

def cli(
    input_path: Path,
    charge: int,
    spin: int,
    func_basis: str,
    max_cycle: int,
    conv_tol: float,
    grid_level: int,
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
    try:
        time_start = time.perf_counter()
        # --------------------------
        # 1) Assemble configuration
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)
        dft_cfg = dict(DFT_KW)
        deep_update(dft_cfg, yaml_cfg.get("dft", {}))

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
        }
        click.echo(_pretty_block("dft", echo_cfg))

        # --------------------------
        # 2) Load geometry
        # --------------------------
        geometry = geom_loader(input_path, coord_type="cart")
        xyz_s, atoms_list = _geometry_to_pyscf_atoms_string(geometry)

        out_dir_path.mkdir(parents=True, exist_ok=True)
        # Write a provenance snapshot of the input geometry
        (out_dir_path / "input_geometry.xyz").write_text(xyz_s if xyz_s.endswith("\n") else (xyz_s + "\n"))
        click.echo(f"[write] Wrote '{out_dir_path / 'input_geometry.xyz'}'.")

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
        using_gpu = False
        engine_label = "pyscf(cpu)"
        try:
            import gpu4pyscf
            gpu4pyscf.activate()  # patch PySCF backends to GPU where supported
            from gpu4pyscf import dft as gdf
            if spin2s == 0:
                mf = gdf.RKS(mol)
            else:
                mf = gdf.UKS(mol)
            using_gpu = True
            engine_label = "gpu4pyscf"
        except Exception:
            from pyscf import dft as pdft
            if spin2s == 0:
                mf = pdft.RKS(mol)
            else:
                mf = pdft.UKS(mol)
            using_gpu = False
            engine_label = "pyscf(cpu)"

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

        # Attempt to enable VV10 for "-v" functionals; ignore if unsupported on GPU
        _maybe_enable_vv10(mf, xc)

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
        # 6) Charges (Mulliken / Löwdin / IAO)
        # --------------------------
        charges = _compute_atomic_charges(mol, mf)

        # --------------------------
        # 7) Save result.yaml
        # --------------------------
        # Per-atom metadata for mapping indices to elements
        atoms_meta = [{"index": i, "element": mol.atom_symbol(i)} for i in range(mol.natm)]

        result_yaml = {
            "energy": {
                "hartree": e_h,
                "kcal_per_mol": e_kcal,
                "converged": converged,
                "scf_time_sec": round(toc_scf - tic_scf, 3),
                "engine": engine_label,
                "used_gpu": bool(using_gpu),
            },
            "charges": {
                "atoms": atoms_meta,
                "mulliken": charges["mulliken"],
                "lowdin": charges["lowdin"],
                "iao": charges["iao"],  # may be None
            },
        }
        (out_dir_path / "result.yaml").write_text(yaml.safe_dump(result_yaml, sort_keys=False, allow_unicode=True))
        click.echo(f"[write] Wrote '{out_dir_path / 'result.yaml'}'.")

        # --------------------------
        # 8) Final print: energies
        # --------------------------
        click.echo(f"E_total (Hartree): {e_h:.12f}")
        click.echo(f"E_total (kcal/mol): {e_kcal:.6f}")

        # Exit codes: 0 if converged, 3 otherwise (for compatibility with prior behavior)
        if not converged:
            click.echo("WARNING: SCF did not converge to the requested tolerance.", err=True)
            sys.exit(3)
        
        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for DFT: {hh:02d}:{mm:02d}:{ss:06.3f}")

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


# Enable `python -m pdb2reaction.dft` execution
if __name__ == "__main__":
    cli()
