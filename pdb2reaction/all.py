# pdb2reaction/all.py

"""
Pipeline driver:  **extract active‑site pockets → run MEP search → merge back to full systems**
+ (optional) TS optimization / pseudo-IRC / thermo / DFT post-processing per reactive segment

Usage (example)
---------------
pdb2reaction all -i a.pdb b.pdb c.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
  --tsopt True --thermo True --dft True

# Also works with residue-id lists as the substrate specification:
#   -c "308,309"

What this does
--------------
1) Runs the binding‑pocket extractor (`extract.extract_api`) on all input PDBs *together*,
   producing per‑structure pocket PDBs under `<out-dir>/pockets/`.
   - All key extractor behaviors can be configured via this CLI (radius, backbone handling, etc.).
2) Reads the **total pocket charge** computed by `extract` and uses it as `-q/--charge`
   for the subsequent MEP search (rounded to nearest integer, with a console note if rounding).
3) Runs the recursive GSM minimum-energy path search (`path_search.cli`) **on the pocket PDBs**.
   - All major path‑search behaviors can be configured via this CLI (spin, max_nodes, optimizer, etc.).
4) **Automatically merges** the pocket MEP back into the **original full PDBs** as reference templates
   (no `--ref-pdb` from the user required), producing full‑system trajectories in the path_search output dir.

5) (Optional; new)
   Per **bond-change** segment:
     - (--tsopt True) Run `ts_opt` on the segment HEI (pocket), then perform a **pseudo-IRC**:
         * compute the imaginary mode at the TS,
         * displace along ±mode (small amplitude) and optimize both to minima (LBFGS),
         * decide forward/backward by **bond-state match** to the segment endpoints from the MEP,
         * rebuild an energy diagram using these (R, TS, P) per reactive segment.
       A small "IRC plot" (TS↔min) is also emitted for each direction.
     - (--thermo True) Run `freq` on (R, TS, P) to get vibrational analysis & **thermochemistry**,
       and build a **Gibbs free-energy diagram** (UMA).
     - (--dft True) Run `dft` single-point on (R, TS, P) and build a **DFT energy diagram**.
       If **both** (--dft True and --thermo True), also build a **DFT//UMA thermal** Gibbs diagram
       (DFT electronic energy + UMA thermal correction to free energy).

Outputs (directory layout)
--------------------------
<out-dir>/
  pockets/
    pocket_<input1_basename>.pdb
    pocket_<input2_basename>.pdb
    ...
  path_search/
    mep.trj
    energy.png
    mep.pdb                     (pocket-only trajectory if pocket inputs were PDB)
    mep_w_ref.pdb               (full-system merged trajectory; references = original inputs)
    mep_w_ref_seg_XX.pdb        (per-segment merged trajectories with covalent changes)
    summary.yaml                (segment barriers, ΔE, etc.)
    ... (segment subfolders)
    tsopt_seg_XX/               (NEW when --tsopt True; per-segment TS/IRC results)
      ts/ ...
      irc/ ...
      freq/ ...                 (when --thermo True)
      dft/  ...                 (when --dft True)
      energy_diagram_tsopt.(html|png)
      energy_diagram_G_UMA.(html|png)                  (when --thermo True)
      energy_diagram_DFT.(html|png)                    (when --dft True)
      energy_diagram_G_DFT_plus_UMA.(html|png)         (when --dft True and --thermo True)

Notes
-----
- Requires Python ≥ 3.10.
- This subcommand intentionally **does not expose** `--ref-pdb`; the original input PDBs are used automatically.
- The extractor runs in multi‑structure union mode to ensure a consistent pocket topology across the ensemble.
- The total charge passed to the path search is taken from the extractor’s **first model** charge summary,
  which matches the extractor’s documented behavior.

CLI parity with underlying tools
--------------------------------
- Extractor options exposed:
  center (PDB path / residue-ID list / residue-name list), radius, radius_het2het, include_H2O, exclude_backbone, add_linkH,
  selected_resn, ligand_charge, verbose.
- Path search options exposed:
  spin, freeze-links, max-nodes, max-cycles, climb, sopt-mode, dump,
  args-yaml, pre-opt, out-dir (the path_search subdir is created inside this).
- (NEW) Post options exposed:
  tsopt (True/False), thermo (True/False), dft (True/False).
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Optional, Tuple, Dict, Any

import sys
import math
import click
import time  # timing
import yaml
import numpy as np
import torch

# Biopython for PDB parsing (post-processing helpers)
from Bio import PDB

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL

# Local imports from the package
from .extract import extract_api
from . import path_search as _path_search
from . import ts_opt as _ts_opt
from . import freq as _freq_cli
from . import dft as _dft_cli
from .uma_pysis import uma_pysis
from .trj2fig import run_trj2fig
from .utils import build_energy_diagram


# -----------------------------
# Helpers
# -----------------------------

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _collect_option_values(argv: Sequence[str], names: Sequence[str]) -> List[str]:
    """
    Robustly collect values following a flag that may appear **once** followed by multiple space-separated values,
    e.g., "-i A B C". This mirrors the behavior implemented in `path_search.cli`.
    """
    vals: List[str] = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    return vals


def _round_charge_with_note(q: float) -> int:
    """
    Cast the extractor's total charge (float) to an integer suitable for the path search.
    If it is not already an integer within 1e-6, round to the nearest integer with a console note.
    """
    q_rounded = int(round(float(q)))
    if not math.isfinite(q):
        raise click.BadParameter(f"Computed total charge is non-finite: {q!r}")
    if abs(float(q) - q_rounded) > 1e-6:
        click.echo(f"[all] NOTE: extractor total charge = {q:g} → rounded to integer {q_rounded} for the path search.")
    return q_rounded


# ---------- Post-processing helpers (minimal, reuse internals) ----------

def _read_summary(summary_yaml: Path) -> List[Dict[str, Any]]:
    """Read path_search/summary.yaml and return segments list (empty if not found)."""
    try:
        if not summary_yaml.exists():
            return []
        data = yaml.safe_load(summary_yaml.read_text(encoding="utf-8")) or {}
        segs = data.get("segments", []) or []
        if not isinstance(segs, list):
            return []
        return segs
    except Exception:
        return []


def _pdb_models_to_coords_and_elems(pdb_path: Path) -> Tuple[List[np.ndarray], List[str]]:
    """
    Return ([coords_model1, coords_model2, ...] in Å), [elements] from a multi-model PDB.
    """
    parser = PDB.PDBParser(QUIET=True)
    st = parser.get_structure("seg", str(pdb_path))
    models = list(st.get_models())
    if not models:
        raise click.ClickException(f"[post] No MODEL found in PDB: {pdb_path}")
    # atom order taken from first model
    atoms0 = [a for a in models[0].get_atoms()]
    elems: List[str] = []
    for a in atoms0:
        el = (a.element or "").strip()
        if not el:
            # fall back: derive from atom name
            nm = a.get_name().strip()
            el = "".join([c for c in nm if c.isalpha()])[:2].title() or "C"
        elems.append(el)
    coords_list: List[np.ndarray] = []
    for m in models:
        atoms = [a for a in m.get_atoms()]
        if len(atoms) != len(atoms0):
            raise click.ClickException(f"[post] Atom count mismatch across models in {pdb_path}")
        coords = np.array([a.get_coord() for a in atoms], dtype=float)
        coords_list.append(coords)
    return coords_list, elems


def _geom_from_angstrom(elems: Sequence[str],
                        coords_ang: np.ndarray,
                        freeze_atoms: Sequence[int]) -> Any:
    """Create a Geometry from Å coordinates using _path_search._new_geom_from_coords (expects Bohr)."""
    coords_bohr = np.asarray(coords_ang, dtype=float) / BOHR2ANG
    return _path_search._new_geom_from_coords(elems, coords_bohr, coord_type="cart", freeze_atoms=freeze_atoms)


def _load_segment_end_geoms(seg_pdb: Path, freeze_atoms: Sequence[int]) -> Tuple[Any, Any]:
    """Load first/last model as Geometries from a per-segment pocket PDB."""
    coords_list, elems = _pdb_models_to_coords_and_elems(seg_pdb)
    gL = _geom_from_angstrom(elems, coords_list[0], freeze_atoms)
    gR = _geom_from_angstrom(elems, coords_list[-1], freeze_atoms)
    return gL, gR


def _compute_imag_mode_direction(ts_geom: Any,
                                 uma_kwargs: Dict[str, Any],
                                 freeze_atoms: Sequence[int]) -> np.ndarray:
    """
    Compute imaginary mode direction (N×3, unit vector in Cartesian space) at TS geometry.
    Uses ts_opt internal helpers to minimize new code.
    """
    # full analytic Hessian (torch tensor)
    H_t = _ts_opt._calc_full_hessian_torch(ts_geom, uma_kwargs=uma_kwargs,
                                           device=torch.device(uma_kwargs.get("device", "cuda" if torch.cuda.is_available() else "cpu")))
    coords_bohr_t = torch.as_tensor(ts_geom.coords.reshape(-1, 3), dtype=H_t.dtype, device=H_t.device)
    # masses in a.u.
    from ase.data import atomic_masses
    masses_amu = np.array([atomic_masses[z] for z in ts_geom.atomic_numbers])
    masses_au_t = torch.as_tensor(masses_amu * _ts_opt.AMU2AU, dtype=H_t.dtype, device=H_t.device)
    mode = _ts_opt._mode_direction_by_root(H_t, coords_bohr_t, masses_au_t,
                                           root=0,
                                           freeze_idx=list(freeze_atoms) if len(freeze_atoms) > 0 else None)
    # ensure unit length
    norm = float(np.linalg.norm(mode.reshape(-1)))
    if norm <= 0:
        raise click.ClickException("[post] Imaginary mode direction has zero norm.")
    return (mode / norm)


def _displaced_geometry_along_mode(geom: Any,
                                   mode_xyz: np.ndarray,
                                   amplitude_ang: float,
                                   freeze_atoms: Sequence[int]) -> Any:
    """
    Displace geometry along mode by ± amplitude (Å). Returns new Geometry.
    """
    coords_bohr = np.asarray(geom.coords3d, dtype=float)  # Bohr
    disp_bohr = (amplitude_ang / BOHR2ANG) * np.asarray(mode_xyz, dtype=float)  # (N,3)
    new_coords_bohr = coords_bohr + disp_bohr
    return _path_search._new_geom_from_coords(geom.atoms, new_coords_bohr, coord_type=geom.coord_type, freeze_atoms=freeze_atoms)


def _save_single_geom_as_pdb_for_tools(g: Any, ref_pdb: Path, out_dir: Path, name: str) -> Path:
    """
    Write a single-geometry XYZ/TRJ with energy and convert to PDB using the pocket ref (for downstream CLI tools).
    Returns PDB path.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_trj = out_dir / f"{name}.trj"
    _path_search._write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)
    pdb_out = out_dir / f"{name}.pdb"
    _path_search._maybe_convert_to_pdb(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
    return pdb_out


def _run_tsopt_on_hei(hei_pdb: Path,
                      charge: int,
                      spin: int,
                      args_yaml: Optional[Path],
                      out_dir: Path,
                      freeze_links: bool) -> Tuple[Path, Any]:
    """
    Run ts_opt CLI on a HEI pocket PDB; return (final_ts_pdb_path, ts_geom)
    """
    ts_dir = out_dir / "ts"
    _ensure_dir(ts_dir)
    ts_args: List[str] = [
        "-i", str(hei_pdb),
        "-q", str(int(charge)),
        "-s", str(int(spin)),
        "--freeze-links", "True" if freeze_links else "False",
        "--max-cycles", "10000",
        "--opt-mode", "light",
        "--dump", "False",
        "--out-dir", str(ts_dir),
    ]
    if args_yaml is not None:
        ts_args.extend(["--args-yaml", str(args_yaml)])

    click.echo(f"[tsopt] Running ts_opt on HEI → out={ts_dir}")
    _saved = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "ts_opt"] + ts_args
        _ts_opt.cli.main(args=ts_args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            raise click.ClickException(f"[tsopt] ts_opt exit code {code}.")
    finally:
        sys.argv = _saved

    # Prefer PDB (ts_opt converts when input is PDB)
    ts_pdb = ts_dir / "final_geometry.pdb"
    if not ts_pdb.exists():
        # fallback: use final .xyz and convert
        xyz_path = ts_dir / "final_geometry.xyz"
        if not xyz_path.exists():
            raise click.ClickException("[tsopt] TS outputs not found.")
        _path_search._maybe_convert_to_pdb(xyz_path, hei_pdb, ts_dir / "final_geometry.pdb")
    ts_pdb = ts_dir / "final_geometry.pdb"
    g_ts = geom_loader(ts_pdb, coord_type="cart")

    # Ensure calculator to have energy on g_ts
    calc = uma_pysis(charge=int(charge), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
    g_ts.set_calculator(calc)
    _ = float(g_ts.energy)

    return ts_pdb, g_ts


def _pseudo_irc_and_match(seg_idx: int,
                          seg_dir: Path,
                          ref_pdb_for_seg: Path,
                          seg_pocket_pdb: Path,
                          g_ts: Any,
                          q_int: int,
                          spin: int,
                          freeze_links_flag: bool) -> Dict[str, Any]:
    """
    From a TS pocket geometry, perform pseudo-IRC:
      - compute imag. mode
      - displace ± (0.25 Å) and optimize both to minima (LBFGS)
      - map each min to left/right segment endpoint by bond-change check
    Returns dict with paths/energies/geoms for {left, ts, right}, and small IRC plots.
    """
    # Freeze parents of link-H if requested
    freeze_atoms = []
    if freeze_links_flag and seg_pocket_pdb.suffix.lower() == ".pdb":
        try:
            freeze_atoms = list(_path_search._freeze_links_for_pdb(seg_pocket_pdb))
        except Exception:
            freeze_atoms = []

    # Mode direction
    uma_kwargs = dict(charge=int(q_int), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
    mode_xyz = _compute_imag_mode_direction(g_ts, uma_kwargs=uma_kwargs, freeze_atoms=freeze_atoms)

    # Displace ± and optimize
    irc_dir = seg_dir / "irc"
    _ensure_dir(irc_dir)
    amp = 0.25  # Å; small stable displacement
    g_plus0 = _displaced_geometry_along_mode(g_ts,  mode_xyz, +amp, freeze_atoms)
    g_minus0 = _displaced_geometry_along_mode(g_ts, mode_xyz, -amp, freeze_atoms)

    # Shared UMA calc
    shared_calc = uma_pysis(charge=int(q_int), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
    # LBFGS settings (reuse defaults)
    sopt_cfg = dict(_path_search.LBFGS_KW)
    sopt_cfg["dump"] = True
    sopt_cfg["out_dir"] = str(irc_dir)

    # Optimize
    g_plus  = _path_search._optimize_single(g_plus0, shared_calc, "lbfgs", sopt_cfg, irc_dir, tag=f"seg_{seg_idx:02d}_irc_plus",  ref_pdb_path=seg_pocket_pdb)
    g_minus = _path_search._optimize_single(g_minus0, shared_calc, "lbfgs", sopt_cfg, irc_dir, tag=f"seg_{seg_idx:02d}_irc_minus", ref_pdb_path=seg_pocket_pdb)

    # IRC mini plots (TS→min)
    try:
        trj_plus  = irc_dir / f"seg_{seg_idx:02d}_irc_plus_opt/optimization.trj"
        trj_minus = irc_dir / f"seg_{seg_idx:02d}_irc_minus_opt/optimization.trj"
        if trj_plus.exists():
            run_trj2fig(trj_plus, [irc_dir / f"irc_plus_plot.png"], unit="kcal", reference="init", reverse_x=False)
        if trj_minus.exists():
            run_trj2fig(trj_minus, [irc_dir / f"irc_minus_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to plot IRC mini plots: {e}", err=True)

    # Load segment endpoints (pocket-only)
    seg_pocket_path = seg_dir.parent / f"mep_seg_{seg_idx:02d}.pdb"
    if not seg_pocket_path.exists():
        # fallback to TRJ if PDB missing
        raise click.ClickException(f"[post] segment pocket PDB not found: {seg_pocket_path}")
    gL_end, gR_end = _load_segment_end_geoms(seg_pocket_path, freeze_atoms)

    # Decide mapping by bond-change
    bond_cfg = dict(_path_search.BOND_KW)
    def _matches(x, y) -> bool:
        try:
            chg, _ = _path_search._has_bond_change(x, y, bond_cfg)
            return (not chg)
        except Exception:
            # fallback: small RMSD threshold
            return (_path_search._rmsd_between(x, y, align=True) < 1e-3)

    # Try to assign (g_plus, g_minus) to (left,right)
    candidates = [("plus", g_plus), ("minus", g_minus)]
    mapping: Dict[str, Any] = {"left": None, "right": None}
    # First pass: exact match on bond changes
    for tag, g in candidates:
        if _matches(g, gL_end) and not _matches(g, gR_end):
            mapping["left"] = (tag, g)
        elif _matches(g, gR_end) and not _matches(g, gL_end):
            mapping["right"] = (tag, g)
    # Second pass: fill missing by RMSD
    for side, g_end in (("left", gL_end), ("right", gR_end)):
        if mapping[side] is None:
            # pick closer one that's not already used
            remain = [(t, gg) for (t, gg) in candidates if mapping.get("left", (None, None))[0] != t and mapping.get("right", (None, None))[0] != t]
            if not remain:
                remain = candidates
            best = min(remain, key=lambda p: _path_search._rmsd_between(p[1], g_end, align=True))
            mapping[side] = best

    # Energies (ensure calculator)
    for _, g in candidates:
        _path_search._ensure_calc_on_geom(g, shared_calc)
        _ = float(g.energy)
    _path_search._ensure_calc_on_geom(g_ts, shared_calc); _ = float(g_ts.energy)

    # Dump tiny TS↔min trj for each direction
    try:
        for side in ("left", "right"):
            tag, gmin = mapping[side]
            trj = irc_dir / f"irc_{side}.trj"
            _path_search._write_xyz_trj_with_energy([g_ts, gmin], [float(g_ts.energy), float(gmin.energy)], trj)
            run_trj2fig(trj, [irc_dir / f"irc_{side}_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception:
        pass

    return {
        "left_min_geom": mapping["left"][1],
        "right_min_geom": mapping["right"][1],
        "ts_geom": g_ts,
        "left_tag": mapping["left"][0],
        "right_tag": mapping["right"][0],
        "freeze_atoms": freeze_atoms,
    }


def _write_segment_energy_diagram(prefix: Path,
                                  labels: List[str],
                                  energies_eh: List[float],
                                  title_note: str) -> None:
    """Write energy diagram (HTML + PNG) using utils.build_energy_diagram."""
    if not energies_eh:
        return
    e0 = energies_eh[0]
    energies_kcal = [(e - e0) * AU2KCALPERMOL for e in energies_eh]
    fig = build_energy_diagram(
        energies=energies_kcal,
        labels=labels,
        ylabel="ΔE (kcal/mol)",
        baseline=True,
        showgrid=False,
        title=f"{prefix.name} {title_note}",
    )
    html = prefix.with_suffix(".html")
    png = prefix.with_suffix(".png")
    fig.write_html(str(html))
    try:
        fig.write_image(str(png), scale=2)
    except Exception:
        pass
    click.echo(f"[diagram] Wrote energy diagram → {html.name} / {png.name}")


def _run_freq_for_state(pdb_path: Path,
                        q_int: int,
                        spin: int,
                        out_dir: Path,
                        args_yaml: Optional[Path],
                        freeze_links: bool) -> Dict[str, Any]:
    """Run freq CLI; return parsed thermo dict (may be empty)."""
    fdir = out_dir
    _ensure_dir(fdir)
    args = [
        "-i", str(pdb_path),
        "-q", str(int(q_int)),
        "-s", str(int(spin)),
        "--freeze-links", "True" if freeze_links else "False",
        "--max-write", "20",
        "--out-dir", str(fdir),
        "--dump", "True",
    ]
    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _saved = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "freq"] + args
        _freq_cli.cli.main(args=args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            click.echo(f"[freq] WARNING: freq exited with code {code}", err=True)
    finally:
        sys.argv = _saved
    # parse thermoanalysis.yaml if any
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _run_dft_for_state(pdb_path: Path,
                       q_int: int,
                       spin: int,
                       out_dir: Path,
                       args_yaml: Optional[Path],
                       func_basis: str = "wb97x-v/def2-tzvp") -> Dict[str, Any]:
    """Run dft CLI; return parsed result.yaml dict (may be empty)."""
    ddir = out_dir
    _ensure_dir(ddir)
    args = [
        "-i", str(pdb_path),
        "-q", str(int(q_int)),
        "-s", str(int(spin)),
        "--func-basis", str(func_basis),
        "--out-dir", str(ddir),
        "--dump", "False",
    ]
    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _saved = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "dft"] + args
        _dft_cli.cli.main(args=args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            click.echo(f"[dft] WARNING: dft exited with code {code}", err=True)
    finally:
        sys.argv = _saved
    y = out_dir / "result.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


# -----------------------------
# CLI
# -----------------------------

@click.command(
    help="Run pocket extraction → MEP search → merge to full PDBs in one shot.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
# ===== Inputs =====
@click.option(
    "-i", "--input", "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True, required=True,
    help=("Two or more **full** PDBs in reaction order (reactant [intermediates ...] product). "
          "You may pass a single '-i' followed by multiple space-separated files (e.g., '-i A.pdb B.pdb C.pdb').")
)
@click.option(
    "-c", "--center", "center_spec",
    type=str, required=True,
    help=("Substrate specification for the extractor: "
          "a PDB path, a residue-ID list like '123,124' or 'A:123,B:456' "
          "(insertion codes OK: '123A' / 'A:123A'), "
          "or a residue-name list like 'GPP,MMT'.")
)
@click.option(
    "--out-dir", "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("./result_all/"), show_default=True,
    help="Top-level output directory for the pipeline."
)
# ===== Extractor knobs (subset of extract.parse_args) =====
@click.option("-r", "--radius", type=float, default=2.6, show_default=True,
              help="Inclusion cutoff (Å) around substrate atoms.")
@click.option("--radius_het2het", type=float, default=0.0, show_default=True,
              help="Independent hetero–hetero cutoff (Å) for non‑C/H pairs.")
@click.option("--include_H2O", "--include-h2o", "include_h2o", type=click.BOOL, default=True, show_default=True,
              help="Include waters (HOH/WAT/TIP3/SOL) in the pocket.")
@click.option("--exclude_backbone", "--exclude-backbone", "exclude_backbone", type=click.BOOL, default=True, show_default=True,
              help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).")
@click.option("--add_linkH", "--add-linkH", "add_linkh", type=click.BOOL, default=True, show_default=True,
              help="Add link hydrogens for severed bonds (carbon-only) in pockets.")
@click.option("--selected_resn", "--selected-resn", type=str, default="", show_default=True,
              help="Force-include residues (comma/space separated; chain/insertion codes allowed).")
@click.option("--ligand_charge", "--ligand-charge", type=str, default=None,
              help=("Either a total charge (number) to distribute across unknown residues "
                    "or a mapping like 'GPP:-3,MMT:-1'."))
@click.option("--verbose", type=click.BOOL, default=True, show_default=True, help="Enable INFO-level logging inside extractor.")
# ===== Path search knobs (subset of path_search.cli) =====
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1).")
@click.option("--freeze-links", "freeze_links_flag", type=click.BOOL, default=True, show_default=True,
              help="For pocket PDB input, freeze parent atoms of link hydrogens.")
@click.option("--max-nodes", type=int, default=10, show_default=True,
              help="Max internal nodes for **segment** GSM (String has max_nodes+2 images including endpoints).")
@click.option("--max-cycles", type=int, default=100, show_default=True, help="Maximum GSM optimization cycles.")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Enable transition-state climbing after growth for the **first** segment in each pair.")
@click.option("--sopt-mode", type=click.Choice(["lbfgs", "rfo", "light", "heavy"], case_sensitive=False),
              default="lbfgs", show_default=True,
              help="Single-structure optimizer kind for HEI±1 and kink nodes.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump GSM / single-structure trajectories during the run.")
@click.option("--args-yaml", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="YAML with extra args for path_search (sections: geom, calc, gs, opt, sopt, bond, search).")
@click.option("--pre-opt", "--pre_opt", "pre_opt", type=click.BOOL, default=True, show_default=True,
              help="If False, skip initial single-structure optimizations of the pocket inputs.")
# ===== NEW: Post-processing toggles =====
@click.option("--tsopt", "do_tsopt", type=click.BOOL, default=False, show_default=True,
              help="TS optimization + pseudo-IRC per reactive segment, and rebuild segment-level energy diagram.")
@click.option("--thermo", "do_thermo", type=click.BOOL, default=False, show_default=True,
              help="Run freq on (R,TS,P) per reactive segment and build Gibbs free-energy diagram (UMA).")
@click.option("--dft", "do_dft", type=click.BOOL, default=False, show_default=True,
              help="Run DFT single-point on (R,TS,P) per reactive segment and build DFT energy diagram. "
                   "If also --thermo True, build DFT//UMA thermal Gibbs diagram.")
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    center_spec: str,
    out_dir: Path,
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    ligand_charge: Optional[str],
    verbose: bool,
    spin: int,
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    sopt_mode: str,
    dump: bool,
    args_yaml: Optional[Path],
    pre_opt: bool,
    do_tsopt: bool,
    do_thermo: bool,
    do_dft: bool,
) -> None:
    """
    The **all** command composes `extract` → `path_search` and hides ref-template bookkeeping.
    It also accepts the sloppy `-i A B C` style like `path_search` does.
    """
    time_start = time.perf_counter()

    # --- Robustly accept a single "-i" followed by multiple paths (like path_search.cli) ---
    argv_all = sys.argv[1:]
    i_vals = _collect_option_values(argv_all, ("-i", "--input"))
    if i_vals:
        i_parsed: List[Path] = []
        for tok in i_vals:
            p = Path(tok)
            if (not p.exists()) or p.is_dir():
                raise click.BadParameter(
                    f"Input path '{tok}' not found or is a directory. "
                    f"When using '-i', list only existing file paths (multiple paths may follow a single '-i')."
                )
            i_parsed.append(p)
        input_paths = tuple(i_parsed)

    # --------------------------
    # Validate input count
    # --------------------------
    if len(input_paths) < 2:
        raise click.BadParameter("Provide at least two PDBs with -i/--input in reaction order.")

    # --------------------------
    # Prepare directories
    # --------------------------
    out_dir = out_dir.resolve()
    pockets_dir = out_dir / "pockets"
    path_dir = out_dir / "path_search"
    _ensure_dir(out_dir)
    _ensure_dir(pockets_dir)
    _ensure_dir(path_dir)

    click.echo("\n=== [all] Stage 1/3 — Active-site pocket extraction (multi-structure union) ===\n")

    # Build per-structure pocket output file list
    pocket_outputs: List[Path] = []
    for p in input_paths:
        pocket_outputs.append((pockets_dir / f"pocket_{p.stem}.pdb").resolve())

    # Run extractor via its public API (multi-structure union mode)
    try:
        ex_res = extract_api(
            complex_pdb=[str(p) for p in input_paths],
            center=center_spec,
            output=[str(p) for p in pocket_outputs],
            radius=float(radius),
            radius_het2het=float(radius_het2het),
            include_H2O=bool(include_h2o),
            exclude_backbone=bool(exclude_backbone),
            add_linkH=bool(add_linkh),
            selected_resn=selected_resn or "",
            ligand_charge=ligand_charge,
            verbose=bool(verbose),
        )
    except Exception as e:
        raise click.ClickException(f"[all] Extractor failed: {e}")

    # Report extractor outputs and charge breakdown
    click.echo("[all] Pocket files:")
    for op in pocket_outputs:
        click.echo(f"  - {op}")

    try:
        cs = ex_res.get("charge_summary", {})
        q_total = float(cs.get("total_charge", 0.0))
        q_prot = float(cs.get("protein_charge", 0.0))
        q_lig = float(cs.get("ligand_total_charge", 0.0))
        q_ion = float(cs.get("ion_total_charge", 0.0))
        click.echo("\n[all] Charge summary from extractor (model #1):")
        click.echo(f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}")
        q_int = _round_charge_with_note(q_total)
    except Exception as e:
        raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")

    # --------------------------
    # Stage 2: Path search on pockets (auto-supplying ref templates = original full PDBs)
    # --------------------------
    click.echo("\n=== [all] Stage 2/3 — MEP search on pocket structures (recursive GSM) ===\n")

    # Build path_search CLI args using *repeated* options (robust for Click)
    ps_args: List[str] = []

    # Inputs: repeat "-i" per pocket to satisfy Click even without their argv aggregator
    for p in pocket_outputs:
        ps_args.extend(["-i", str(p)])

    # Charge & spin
    ps_args.extend(["-q", str(q_int)])
    ps_args.extend(["-s", str(int(spin))])

    # Freeze-links, nodes, cycles, climb, optimizer, dump, out-dir, pre-opt, args-yaml
    ps_args.extend(["--freeze-links", "True" if freeze_links_flag else "False"])
    ps_args.extend(["--max-nodes", str(int(max_nodes))])
    ps_args.extend(["--max-cycles", str(int(max_cycles))])
    ps_args.extend(["--climb", "True" if climb else "False"])
    ps_args.extend(["--sopt-mode", str(sopt_mode)])
    ps_args.extend(["--dump", "True" if dump else "False"])
    ps_args.extend(["--out-dir", str(path_dir)])
    ps_args.extend(["--pre-opt", "True" if pre_opt else "False"])
    if args_yaml is not None:
        ps_args.extend(["--args-yaml", str(args_yaml)])

    # Auto-provide ref templates (original full PDBs) for full-system merge. Repeat "--ref-pdb" per file.
    for p in input_paths:
        ps_args.extend(["--ref-pdb", str(p)])

    click.echo("[all] Invoking path_search with arguments:")
    click.echo("  " + " ".join(ps_args))

    # CRITICAL FIX:
    _saved_argv = list(sys.argv)
    try:
        sys.argv = ["pdb2reaction", "path_search"] + ps_args
        _path_search.cli.main(args=ps_args, standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            raise click.ClickException(f"[all] path_search terminated with exit code {code}.")
    except Exception as e:
        raise click.ClickException(f"[all] path_search failed: {e}")
    finally:
        sys.argv = _saved_argv

    # --------------------------
    # Stage 3: Merge (performed by path_search when --ref-pdb was supplied)
    # --------------------------
    click.echo("\n=== [all] Stage 3/3 — Merge into full-system templates ===\n")
    click.echo("[all] Merging was carried out by path_search using the original inputs as templates.")
    click.echo(f"[all] Final products can be found under: {path_dir}")
    click.echo("  - mep_w_ref.pdb            (full-system merged trajectory)")
    click.echo("  - mep_w_ref_seg_XX.pdb     (per-segment merged trajectories for covalent-change segments)")
    click.echo("  - summary.yaml             (segment barriers, ΔE, labels)")
    click.echo("  - energy.png / energy_diagram.*")
    click.echo("\n=== [all] Pipeline finished successfully ===\n")

    # --------------------------
    # Optional Stage 4: TSOPT / THERMO / DFT (per reactive segment)
    # --------------------------
    if not (do_tsopt or do_thermo or do_dft):
        # Elapsed time
        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[all] Total Elapsed: {hh:02d}:{mm:02d}:{ss:06.3f}")
        return

    click.echo("\n=== [all] Stage 4 — Post-processing per reactive segment ===\n")

    # Load segment summary
    summary_yaml = path_dir / "summary.yaml"
    segments = _read_summary(summary_yaml)
    if not segments:
        click.echo("[post] No segments found in summary; nothing to do.")
        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[all] Total Elapsed: {hh:02d}:{mm:02d}:{ss:06.3f}")
        return

    # Iterate only bond-change segments (kind='seg' and bond_changes not empty and not '(no covalent...)')
    reactive = [s for s in segments if (s.get("kind", "seg") == "seg" and str(s.get("bond_changes", "")).strip() and str(s.get("bond_changes", "")).strip() != "(no covalent changes detected)")]
    if not reactive:
        click.echo("[post] No bond-change segments. Skipping TS/thermo/DFT.")
        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[all] Total Elapsed: {hh:02d}:{mm:02d}:{ss:06.3f}")
        return

    # For each reactive segment
    for s in reactive:
        seg_idx = int(s.get("index", 0) or 0)
        seg_tag = s.get("tag", f"seg_{seg_idx:02d}")
        click.echo(f"\n--- [post] Segment {seg_idx:02d} ({seg_tag}) ---")

        seg_root = path_dir  # base
        seg_dir = seg_root / f"tsopt_seg_{seg_idx:02d}"
        _ensure_dir(seg_dir)

        # HEI pocket file prepared by path_search (only for bond-change segments)
        hei_pocket_pdb = seg_root / f"hei_seg_{seg_idx:02d}.pdb"
        if not hei_pocket_pdb.exists():
            click.echo(f"[post] WARNING: HEI pocket PDB not found for segment {seg_idx:02d}; skipping TSOPT.", err=True)
            continue

        # 4.1 TS optimization (optional; still needed to drive IRC & diagrams)
        if do_tsopt:
            ts_pdb, g_ts = _run_tsopt_on_hei(hei_pocket_pdb, q_int, spin, args_yaml, seg_dir, freeze_links_flag)
        else:
            # If TSOPT off: use the GSM HEI (pocket) as TS geometry
            ts_pdb = hei_pocket_pdb
            g_ts = geom_loader(ts_pdb, coord_type="cart")
            calc = uma_pysis(charge=int(q_int), spin=int(spin), model="uma-s-1p1", task_name="omol", device="auto")
            g_ts.set_calculator(calc); _ = float(g_ts.energy)

        # 4.2 Pseudo-IRC & mapping to (left,right)
        irc_res = _pseudo_irc_and_match(seg_idx=seg_idx,
                                        seg_dir=seg_dir,
                                        ref_pdb_for_seg=ts_pdb,
                                        seg_pocket_pdb=hei_pocket_pdb,
                                        g_ts=g_ts,
                                        q_int=q_int,
                                        spin=spin,
                                        freeze_links_flag=freeze_links_flag)

        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        # Save standardized PDBs for tools
        struct_dir = seg_dir / "structures"
        _ensure_dir(struct_dir)
        pL = _save_single_geom_as_pdb_for_tools(gL, hei_pocket_pdb, struct_dir, "reactant_like")
        pT = _save_single_geom_as_pdb_for_tools(gT, hei_pocket_pdb, struct_dir, "ts")
        pR = _save_single_geom_as_pdb_for_tools(gR, hei_pocket_pdb, struct_dir, "product_like")

        # 4.3 Segment-level energy diagram from UMA (R,TS,P)
        eR = float(gL.energy)
        eT = float(gT.energy)
        eP = float(gR.energy)
        _write_segment_energy_diagram(seg_dir / "energy_diagram_tsopt",
                                      labels=["R", f"TS{seg_idx}", "P"],
                                      energies_eh=[eR, eT, eP],
                                      title_note="(UMA, TSOPT/IRC)")

        # 4.4 Thermochemistry (UMA freq) and Gibbs diagram
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        if do_thermo:
            click.echo(f"[thermo] Segment {seg_idx:02d}: freq on R/TS/P")
            tR = _run_freq_for_state(pL, q_int, spin, seg_dir / "freq" / "R", args_yaml, freeze_links_flag)
            tT = _run_freq_for_state(pT, q_int, spin, seg_dir / "freq" / "TS", args_yaml, freeze_links_flag)
            tP = _run_freq_for_state(pR, q_int, spin, seg_dir / "freq" / "P", args_yaml, freeze_links_flag)
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            try:
                GR = float(tR.get("sum_EE_and_thermal_free_energy_ha", eR))
                GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", eT))
                GP = float(tP.get("sum_EE_and_thermal_free_energy_ha", eP))
                _write_segment_energy_diagram(seg_dir / "energy_diagram_G_UMA",
                                              labels=["R", f"TS{seg_idx}", "P"],
                                              energies_eh=[GR, GT, GP],
                                              title_note="(Gibbs, UMA)")
            except Exception as e:
                click.echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)

        # 4.5 DFT single-point and (optionally) DFT//UMA Gibbs
        if do_dft:
            click.echo(f"[dft] Segment {seg_idx:02d}: DFT on R/TS/P")
            dR = _run_dft_for_state(pL, q_int, spin, seg_dir / "dft" / "R", args_yaml, func_basis="wb97x-v/def2-tzvp")
            dT = _run_dft_for_state(pT, q_int, spin, seg_dir / "dft" / "TS", args_yaml, func_basis="wb97x-v/def2-tzvp")
            dP = _run_dft_for_state(pR, q_int, spin, seg_dir / "dft" / "P", args_yaml, func_basis="wb97x-v/def2-tzvp")
            try:
                eR_dft = float(((dR or {}).get("energy", {}) or {}).get("hartree", np.nan))
                eT_dft = float(((dT or {}).get("energy", {}) or {}).get("hartree", np.nan))
                eP_dft = float(((dP or {}).get("energy", {}) or {}).get("hartree", np.nan))
                if all(map(np.isfinite, [eR_dft, eT_dft, eP_dft])):
                    _write_segment_energy_diagram(seg_dir / "energy_diagram_DFT",
                                                  labels=["R", f"TS{seg_idx}", "P"],
                                                  energies_eh=[eR_dft, eT_dft, eP_dft],
                                                  title_note="(DFT wb97x-v/def2-tzvp)")
                else:
                    click.echo("[dft] WARNING: some DFT energies missing; diagram skipped.", err=True)
            except Exception as e:
                click.echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            # DFT//UMA thermal Gibbs (E_DFT + ΔG_therm(UMA))
            if do_thermo:
                try:
                    dG_R = float((thermo_payloads.get("R", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_T = float((thermo_payloads.get("TS", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_P = float((thermo_payloads.get("P", {}) or {}).get("thermal_correction_free_energy_ha", 0.0))
                    eR_dft = float(((dR or {}).get("energy", {}) or {}).get("hartree", eR))
                    eT_dft = float(((dT or {}).get("energy", {}) or {}).get("hartree", eT))
                    eP_dft = float(((dP or {}).get("energy", {}) or {}).get("hartree", eP))
                    GR_dftUMA = eR_dft + dG_R
                    GT_dftUMA = eT_dft + dG_T
                    GP_dftUMA = eP_dft + dG_P
                    _write_segment_energy_diagram(seg_dir / "energy_diagram_G_DFT_plus_UMA",
                                                  labels=["R", f"TS{seg_idx}", "P"],
                                                  energies_eh=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                                                  title_note="(Gibbs, DFT//UMA)")
                except Exception as e:
                    click.echo(f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}", err=True)

    # --------------------------
    # Elapsed time
    # --------------------------
    elapsed = time.perf_counter() - time_start
    hh = int(elapsed // 3600)
    mm = int((elapsed % 3600) // 60)
    ss = elapsed - (hh * 3600 + mm * 60)
    click.echo(f"[all] Total Elapsed: {hh:02d}:{mm:02d}:{ss:06.3f}")


if __name__ == "__main__":
    cli()
