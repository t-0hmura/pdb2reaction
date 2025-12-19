# pdb2reaction/path_opt.py

"""
path_opt — Minimum-energy path (MEP) optimization via GSM or DMF with a UMA calculator
=====================================================================================

Usage (CLI)
-----------
    pdb2reaction path-opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} -q CHARGE -m MULT \
                        [--mep-mode {gsm|dmf}] [--freeze-links BOOL] [--max-nodes N] [--max-cycles N] \
                        [--climb BOOL] [--dump BOOL] [--thresh PRESET] \
                        [--convert-files/--no-convert-files] \
                        [--out-dir DIR] [--args-yaml FILE]

Examples
--------
    # Minimal: two endpoints, neutral singlet
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1

    # Typical full run with YAML overrides, dumps, and a convergence preset
    pdb2reaction path-opt -i reac.pdb prod.pdb -q 0 -m 1 \
      --freeze-links True --max-nodes 10 --max-cycles 100 \
      --climb True --dump True --out-dir ./result_path_opt/ \
      --thresh gau_tight --args-yaml ./args.yaml

Description
-----------
- Optimizes a minimum-energy path between two endpoints using GSM (pysisyphus `GrowingString` + `StringOptimizer`) or DMF, with UMA as the calculator (via `uma_pysis`).
- Inputs: two structures (.pdb or .xyz). If a PDB is provided and `--freeze-links=True` (default), parent atoms of link hydrogens are added to `freeze_atoms` (0-based indices).
- Configuration via YAML with sections `geom`, `calc`, `gs`, `opt`, and single-structure optimizer sections such as `sopt.lbfgs` / `sopt.rfo` (also accepting `opt.lbfgs` / `opt.rfo` and top-level `lbfgs` / `rfo`). Precedence: YAML > CLI > built-in defaults.
- Optional endpoint pre-optimization: with `--preopt=True` (default False), each endpoint is relaxed individually via single-structure LBFGS ("light", default) or RFO ("heavy") before alignment and MEP search. The iteration limit for this pre-optimization is controlled independently by `--preopt-max-cycles` (default: 10000) for whichever optimizer is selected.
- Path generator: `--mep-mode` accepts GSM or DMF, with GSM enabled by default to match the CLI default.
- Alignment: before optimization, all inputs after the first are rigidly Kabsch-aligned to the first structure using an external routine with a short relaxation. `StringOptimizer.align` is disabled. If either endpoint specifies `freeze_atoms`, the RMSD fit uses only those atoms and the resulting rigid transform is applied to all atoms.
- With `--climb=True` (default), a climbing-image step refines the highest-energy image. Lanczos-based tangent estimation (`gs.climb_lanczos`) is enabled by default and follows the `--climb` flag; YAML can still override it.
- `--thresh` sets the convergence preset used by the string optimizer, the optional endpoint pre-optimization, and the pre-alignment refinement (e.g., `gau_loose|gau|gau_tight|gau_vtight|baker|never`).
- `--fix-ends=True` fixes both endpoint geometries during GSM (`fix_first=True`, `fix_last=True`).
- After optimization, the highest-energy image (HEI) is identified as the highest-energy internal local maximum (preferring internal nodes). If none exist, the maximum among internal nodes is used; if there are no internal nodes, the global maximum is used. The selected HEI is exported.

Outputs (& Directory Layout)
----------------------------
out_dir/ (default: ./result_path_opt/)
  ├─ final_geometries.trj        # XYZ trajectory of the optimized path; comment line carries per-image energy when available
  ├─ final_geometries.pdb        # Converted from .trj when the *first* endpoint is a PDB and conversion is enabled
  ├─ hei.xyz                     # Highest-energy image with energy on the comment line
  ├─ hei.pdb                     # HEI converted to PDB when the *first* endpoint is a PDB and conversion is enabled
  ├─ hei.gjf                     # HEI written using a detected .gjf template when available and conversion is enabled
  ├─ align_refine/               # Files from external alignment/refinement
  └─ <optimizer dumps / restarts>  # Emitted when dumping is enabled (e.g., via `--dump` and/or YAML `dump_restart` settings)

Notes
-----
- Charge/spin: `-q/--charge` (recommended for non-`.gjf` inputs; otherwise defaults to 0) and `-m/--multiplicity`
  are reconciled with any `.gjf` template values: explicit CLI options win. When ``-q`` is omitted but ``--ligand-charge`` is set,
  the full complex is treated as an enzyme–substrate system and the total charge is inferred using ``extract.py``’s residue-aware
  logic. If neither `-q` nor `--ligand-charge` is supplied, the charge falls back to 0; set it explicitly to avoid unphysical
  conditions.
- Coordinates are Cartesian; `freeze_atoms` use 0-based indices. With `--freeze-links=True` and PDB inputs, link-hydrogen parents are added automatically.
- `--max-nodes` sets the number of internal nodes; the string has (max_nodes + 2) images including endpoints.
- `--max-cycles` limits optimization; after full growth, the same bound applies to additional refinement.
- `--preopt-max-cycles` limits only the optional endpoint single-structure preoptimization (LBFGS or RFO, selected via `--opt-mode`) and does not affect `--max-cycles`.
- Format-aware XYZ/TRJ → PDB/GJF conversions respect the global `--convert-files/--no-convert-files` toggle (default: enabled).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import sys
import traceback
import textwrap
import time

import click
import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from .uma_pysis import uma_pysis, GEOM_KW_DEFAULT, CALC_KW as _UMA_CALC_KW

from .utils import (
    convert_xyz_to_pdb,
    detect_freeze_links,
    detect_freeze_links_safe,
    load_yaml_dict,
    apply_yaml_overrides,
    deep_update,
    pretty_block,
    format_geom_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    fill_charge_spin_from_gjf,
    _derive_charge_from_ligand_charge,
    maybe_convert_xyz_to_gjf,
    set_convert_file_enabled,
    convert_xyz_like_outputs,
    PreparedInputStructure,
)
from .opt import (
    LBFGS_KW as _LBFGS_KW,
    RFO_KW as _RFO_KW,
)
from .align_freeze_atoms import align_and_refine_sequence_inplace


# -----------------------------------------------
# Defaults (overridden by YAML/CLI)
# -----------------------------------------------

# Geometry (input handling)
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)

CALC_KW: Dict[str, Any] = dict(_UMA_CALC_KW)

# DMF (Direct Max Flux + (C)FB-ENM)
DMF_KW: Dict[str, Any] = {
    # Top-level interpolate_fbenm options
    "correlated": True,             # Add CFB_ENM for correlated paths
    "sequential": True,             # Enable staged barrier construction during optimization
    "fbenm_only_endpoints": False,  # If False, use all images (not just endpoints) for ENM references

    # FB_ENM_Bonds options (fbenm_options)
    "fbenm_options": {
        "delta_scale": 0.2,         # Scale for the distance penalty width
        "bond_scale": 1.25,         # Bond test: d_ij < bond_scale * (r_cov_i + r_cov_j)
        "fix_planes": True,         # Add plane constraints to preserve planarity
        "two_hop_mode": "sparse",   # Two-hop neighbor construction ("sparse" | "dense")
    },

    # CFB_ENM options (cfbenm_options)
    "cfbenm_options": {
        "bond_scale": 1.25,         # neighbor cutoff multiplier for CFB-ENM graph construction
        "corr0_scale": 1.10,        # d_corr0 ~ corr0_scale * d_bond
        "corr1_scale": 1.50,        # scale for first correlation shell distance
        "corr2_scale": 1.60,        # scale for second correlation shell distance
        "eps": 0.05,                # sqrt(pp^2 + eps^2) term's epsilon
        "pivotal": True,            # enable pivotal constraints in the CFB-ENM
        "single": True,             # enforce single-path constraint in correlation graph
        "remove_fourmembered": True,# prune four-membered rings in the correlation network
        "two_hop_mode": "sparse",   # Two-hop construction on the CFB_ENM side
    },

    # DirectMaxFlux core options (forwarded via dmf_options)
    "dmf_options": {
        "remove_rotation_and_translation": False,  # Do not explicitly constrain rigid-body motions
        "mass_weighted": False,                    # Whether to use mass-weighted velocity norms
        "parallel": False,                         # allow parallel execution inside DMF core
        "eps_vel": 0.01,                           # stabilization epsilon for velocity norms
        "eps_rot": 0.01,                           # stabilization epsilon for rotational terms
        "beta": 10.0,                              # "Beta" for the geometric action
        "update_teval": False,                     # Control node relocation from interpolate_fbenm
    },

    # Strength of fix_atoms harmonic restraints
    "k_fix": 100.0,                                # [eV/Å^2]
}

# GrowingString (path representation)
GS_KW: Dict[str, Any] = {
    "fix_first": True,           # keep the first image fixed during optimization
    "fix_last": True,            # keep the last image fixed during optimization
    "max_nodes": 10,            # int, internal nodes; total images = max_nodes + 2 including endpoints
    "perp_thresh": 5e-3,        # float, frontier growth criterion (RMS/NORM of perpendicular force)
    "reparam_check": "rms",     # str, "rms" | "norm"; convergence check metric after reparam
    "reparam_every": 1,         # int, reparametrize every N steps while growing
    "reparam_every_full": 1,    # int, reparametrize every N steps after fully grown
    "param": "equi",            # str, "equi" (even spacing) | "energy" (weight by energy)
    "max_micro_cycles": 10,     # int, micro-optimization cycles per macro iteration
    "reset_dlc": True,          # bool, reset DLC coordinates when appropriate
    "climb": True,              # bool, enable climbing image
    "climb_rms": 5e-4,          # float, RMS force threshold to start climbing image
    "climb_lanczos": True,      # bool, use Lanczos to estimate the HEI tangent (enabled by default; tied to --climb)
    "climb_lanczos_rms": 5e-4,  # float, RMS force threshold for Lanczos tangent
    "climb_fixed": False,       # bool, fix the HEI image index instead of adapting it
    "scheduler": None,          # Optional[str], execution scheduler; None = serial (shared calculator)
}

# StringOptimizer (optimization control)
STOPT_KW: Dict[str, Any] = {
    "type": "string",           # str, tag for bookkeeping / output labelling
    "stop_in_when_full": 300,   # int, allow N extra cycles after the string is fully grown
    "align": False,             # bool, keep internal align disabled; use external Kabsch alignment instead
    "scale_step": "global",     # str, "global" | "per_image" scaling policy
    "max_cycles": 300,          # int, maximum macro cycles for the optimizer
    "dump": False,              # bool, write optimizer trajectory to disk
    "dump_restart": False,      # bool | int, write restart YAML every N cycles (False disables)
    "reparam_thresh": 0.0,      # float, convergence threshold for reparametrization
    "coord_diff_thresh": 0.0,   # float, tolerance for coordinate difference before pruning
    "out_dir": "./result_path_opt/",  # str, output directory for optimizer artifacts
    "print_every": 10,          # int, status print frequency (cycles)
}


def _load_two_endpoints(
    inputs: Sequence[PreparedInputStructure],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> Sequence:
    """
    Load the two endpoint structures and set `freeze_atoms` as needed.
    """
    geoms = []
    for prepared in inputs:
        geom_path = prepared.geom_path
        src_path = prepared.source_path
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        if auto_freeze_links and src_path.suffix.lower() == ".pdb":
            detected = detect_freeze_links_safe(src_path)
            freeze = merge_freeze_atom_indices(cfg, detected)
            if detected and freeze:
                click.echo(
                    f"[freeze-links] {src_path.name}: Freeze atoms (0-based): "
                    f"{','.join(map(str, freeze))}"
                )
        else:
            freeze = merge_freeze_atom_indices(cfg)
        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


def _maybe_convert_to_pdb(
    src_path: Path, ref_pdb_path: Optional[Path], out_path: Optional[Path] = None
) -> Optional[Path]:
    """
    Convert an XYZ/TRJ path to PDB using a reference topology when available.

    Returns the written path on success; otherwise None. Conversion is skipped when
    the source does not exist, lacks an XYZ/TRJ suffix, or no reference PDB is
    provided.
    """
    try:
        if ref_pdb_path is None:
            return None
        src_path = Path(src_path)
        if (not src_path.exists()) or src_path.suffix.lower() not in (".xyz", ".trj"):
            return None

        out_path = out_path if out_path is not None else src_path.with_suffix(".pdb")
        convert_xyz_to_pdb(src_path, ref_pdb_path, out_path)
        click.echo(f"[convert] Wrote '{out_path}'.")
        return out_path
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert '{src_path}' to PDB: {e}", err=True)
        return None


def _select_hei_index(energies: Sequence[float]) -> int:
    """Pick an HEI index preferring internal local maxima."""

    E = np.array(energies, dtype=float)
    nE = int(len(E))
    hei_idx = None
    if nE >= 3:
        candidates = [i for i in range(1, nE - 1) if (E[i] > E[i - 1] and E[i] > E[i + 1])]
        if candidates:
            hei_idx = int(max(candidates, key=lambda i: E[i]))
        else:
            hei_idx = 1 + int(np.argmax(E[1:-1]))
    if hei_idx is None:
        hei_idx = int(np.argmax(E))
    return hei_idx


def _write_ase_trj_with_energy(images: Sequence[Any], energies: Sequence[float], path: Path) -> None:
    """Write an XYZ `.trj` from ASE Atoms with the energy on line 2."""

    blocks = []
    for atoms, e in zip(images, np.array(energies, dtype=float)):
        symbols = atoms.get_chemical_symbols()
        coords = atoms.get_positions()
        lines = [str(len(symbols)), f"{e:.12f}"]
        lines.extend(
            f"{sym} {x:.15f} {y:.15f} {z:.15f}" for sym, (x, y, z) in zip(symbols, coords)
        )
        blocks.append("\n".join(lines) + "\n")

    with open(path, "w") as f:
        f.write("".join(blocks))


@dataclass
class DMFMepResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int


def _run_dmf_mep(
    geoms: Sequence[Any],
    calc_cfg: Dict[str, Any],
    out_dir_path: Path,
    prepared_inputs: Sequence[PreparedInputStructure],
    max_nodes: int,
    fix_atoms: Sequence[int],
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> DMFMepResult:
    """Run Direct Max Flux (DMF) MEP optimization between two endpoints.

    References:
    [1] S.-i. Koda and  S. Saito, Locating Transition States by Variational Reaction Path Optimization with an Energy-Derivative-Free Objective Function, JCTC, 20, 2798–2811 (2024). [doi: 10.1021/acs.jctc.3c01246]
    [2] S.-i. Koda and  S. Saito, Flat-bottom Elastic Network Model for Generating Improved Plausible Reaction Paths, JCTC, 20, 7176−7187 (2024). [doi: 10.1021/acs.jctc.4c00792]
    [3] S.-i. Koda and  S. Saito, Correlated Flat-bottom Elastic Network Model for Improved Bond Rearrangement in Reaction Paths, JCTC, 21, 3513−3522 (2025). [doi: 10.1021/acs.jctc.4c01549]
    """

    try:
        import torch
        from ase.io import read as ase_read
        from ase.io import write as ase_write
        from fairchem.core import pretrained_mlip, FAIRChemCalculator
        from torch_dmf import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise RuntimeError(
            "DMF mode requires torch, ase, fairchem, cyiopt, and torch_dmf to be installed."
        ) from e

    def _geom_to_ase(g: Any):
        from io import StringIO

        return ase_read(StringIO(g.as_xyz()), format="xyz")

    device = str(calc_cfg.get("device", "auto"))
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"

    fix_atoms = list(sorted(set(map(int, fix_atoms))))

    ref_images = [_geom_to_ase(g) for g in geoms]
    primary_prepared = prepared_inputs[0] if prepared_inputs else None
    ref_pdb = (
        primary_prepared.source_path.resolve()
        if primary_prepared and primary_prepared.source_path.suffix.lower() == ".pdb"
        else None
    )
    needs_pdb = ref_pdb is not None
    needs_gjf = bool(primary_prepared and primary_prepared.is_gjf)

    charge = int(calc_cfg.get("charge", 0))
    spin = int(calc_cfg.get("spin", 1))
    for img in ref_images:
        img.info["charge"] = charge
        img.info["spin"] = spin

    predictor = pretrained_mlip.get_predict_unit(
        calc_cfg.get("model", "uma-s-1p1"),
        device=device,
    )

    calc_uma = FAIRChemCalculator(
        predictor,
        task_name=str(calc_cfg.get("task_name", "omol")),
    )

    dmf_cfg = deep_update(dict(DMF_KW), dmf_cfg)
    fbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("fbenm_options", {}))
    cfbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("cfbenm_options", {}))
    dmf_opts: Dict[str, Any] = dict(dmf_cfg.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
    k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))

    mxflx_fbenm = interpolate_fbenm(
        ref_images,
        nmove=max(1, int(max_nodes)),
        fbenm_only_endpoints=bool(dmf_cfg.get("fbenm_only_endpoints", False)),
        correlated=bool(dmf_cfg.get("correlated", False)),
        sequential=bool(dmf_cfg.get("sequential", False)),
        output_file=str(out_dir_path / "dmf_fbenm_ipopt.out"),
        device=device,
        dtype="float64",
        fix_atoms=fix_atoms,
        fbenm_options=fbenm_opts,
        cfbenm_options=cfbenm_opts,
        dmf_options=dmf_opts,
        k_fix=k_fix,
    )

    initial_trj = out_dir_path / "dmf_initial.trj"
    ase_write(initial_trj, mxflx_fbenm.images, format="xyz")
    if primary_prepared is not None and needs_pdb:
        convert_xyz_like_outputs(
            initial_trj,
            primary_prepared,
            ref_pdb_path=ref_pdb,
            out_pdb_path=initial_trj.with_suffix(".pdb") if needs_pdb else None,
        )
    coefs = mxflx_fbenm.coefs.copy()

    mxflx = DirectMaxFlux(
        ref_images,
        coefs=coefs,
        nmove=max(1, int(max_nodes)),
        update_teval=update_teval,
        device=device,
        dtype="float64",
        fix_atoms=fix_atoms,
        remove_rotation_and_translation=bool(
            dmf_opts.get("remove_rotation_and_translation", False)
        ),
        mass_weighted=bool(dmf_opts.get("mass_weighted", False)),
        parallel=bool(dmf_opts.get("parallel", False)),
        eps_vel=float(dmf_opts.get("eps_vel", 0.01)),
        eps_rot=float(dmf_opts.get("eps_rot", 0.01)),
        beta=float(dmf_opts.get("beta", 10.0)),
        k_fix=k_fix,
    )

    for image in mxflx.images:
        if "charge" not in image.info:
            image.info["charge"] = charge
        if "spin" not in image.info:
            image.info["spin"] = spin
        image.calc = calc_uma

    mxflx.add_ipopt_options({"output_file": str(out_dir_path / "dmf_ipopt.out")})
    mxflx.solve(tol="tight")

    calc_eval_kw = dict(calc_cfg)
    if fix_atoms:
        calc_eval_kw.setdefault("freeze_atoms", fix_atoms)

    calc_eval = uma_pysis(**calc_eval_kw)

    energies = []
    for image in mxflx.images:
        elems = image.get_chemical_symbols()
        coords_bohr = np.asarray(image.get_positions(), dtype=float).reshape(-1, 3) * ANG2BOHR
        energies.append(float(calc_eval.get_energy(elems, coords_bohr)["energy"]))
    hei_idx = _select_hei_index(energies)

    final_trj = out_dir_path / "final_geometries.trj"
    _write_ase_trj_with_energy(mxflx.images, energies, final_trj)
    if primary_prepared is not None and needs_pdb:
        convert_xyz_like_outputs(
            final_trj,
            primary_prepared,
            ref_pdb_path=ref_pdb,
            out_pdb_path=final_trj.with_suffix(".pdb") if needs_pdb else None,
        )
    if primary_prepared is not None and (needs_pdb or needs_gjf):
        hei_tmp = out_dir_path / "hei.xyz"
        _write_ase_trj_with_energy([mxflx.images[hei_idx]], [energies[hei_idx]], hei_tmp)
        convert_xyz_like_outputs(
            hei_tmp,
            primary_prepared,
            ref_pdb_path=ref_pdb,
            out_pdb_path=out_dir_path / "hei.pdb" if needs_pdb else None,
            out_gjf_path=out_dir_path / "hei.gjf" if needs_gjf else None,
        )

    return DMFMepResult(images=list(mxflx.images), energies=list(energies), hei_idx=int(hei_idx))


def _optimize_single(
    g,
    shared_calc,
    sopt_kind: str,
    sopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    prepared_input: Optional[PreparedInputStructure],
):
    """
    Single-structure optimization (LBFGS or RFO) shared by path_opt and path_search.
    """
    g.set_calculator(shared_calc)

    seg_dir = out_dir / f"{tag}_{sopt_kind}_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(sopt_cfg)
    args["out_dir"] = str(seg_dir)

    if sopt_kind == "lbfgs":
        opt = LBFGS(g, **args)
    else:
        opt = RFOptimizer(g, **args)

    click.echo(f"\n=== [{tag}] Single-structure {sopt_kind.upper()} started ===\n")
    opt.run()
    click.echo(f"\n=== [{tag}] Single-structure {sopt_kind.upper()} finished ===\n")

    try:
        final_xyz = Path(opt.final_fn)
        if prepared_input is not None:
            ref_pdb = (
                prepared_input.source_path.resolve()
                if prepared_input.source_path.suffix.lower() == ".pdb"
                else None
            )
            needs_pdb = ref_pdb is not None
            needs_gjf = prepared_input.is_gjf
            if needs_pdb or needs_gjf:
                convert_xyz_like_outputs(
                    final_xyz,
                    prepared_input,
                    ref_pdb_path=ref_pdb,
                    out_pdb_path=final_xyz.with_suffix(".pdb") if needs_pdb else None,
                    out_gjf_path=final_xyz.with_suffix(".gjf") if needs_gjf else None,
                )
        g_final = geom_loader(final_xyz, coord_type=g.coord_type)
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        g_final.set_calculator(shared_calc)
        return g_final
    except Exception:
        return g


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="MEP optimization via the Growing String method.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    nargs=2,
    required=True,
    help="Two endpoint structures (reactant and product); accepts .pdb or .xyz.",
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
)
@click.option("-q", "--charge", type=int, required=False, help="Charge of the ML region.")
@click.option(
    "--workers",
    type=int,
    default=CALC_KW["workers"],
    show_default=True,
    help="UMA predictor workers; >1 spawns a parallel predictor (disables analytic Hessian).",
)
@click.option(
    "--workers-per-nodes",
    "workers_per_nodes",
    type=int,
    default=CALC_KW["workers_per_nodes"],
    show_default=True,
    help="Workers per node when using a parallel UMA predictor (workers>1).",
)
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
    default=1,
    show_default=True,
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--freeze-links",
    "freeze_links_flag",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="If a PDB is provided, freeze the parent atoms of link hydrogens.",
)
@click.option(
    "--max-nodes",
    type=int,
    default=10,
    show_default=True,
    help="Number of internal nodes (string has max_nodes+2 images including endpoints).",
)
@click.option(
    "--max-cycles",
    type=int,
    default=300,
    show_default=True,
    help="Maximum optimization cycles.",
)
@click.option(
    "--climb",
    type=click.BOOL,
    default=True,
    show_default=True,
    help="Search for a transition state (climbing image) after path growth.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Single-structure optimizer for endpoint preoptimization: light (=LBFGS) or heavy (=RFO).",
)
@click.option(
    "--dump",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Dump optimizer trajectory/restarts during the run.",
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB/GJF companions based on the input format.",
)
@click.option(
    "--out-dir",
    "out_dir",
    type=str,
    default="./result_path_opt/",
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help="Convergence preset for the string optimizer, pre-alignment refinement, and endpoint preoptimization (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt, sopt.lbfgs, sopt.rfo).",
)
@click.option(
    "--preopt",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="If True, preoptimize each endpoint via the selected single-structure optimizer (LBFGS/RFO) before alignment and GSM.",
)
@click.option(
    "--preopt-max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum cycles for endpoint preoptimization (applies to the chosen optimizer; only used when --preopt True).",
)
@click.option(
    "--fix-ends",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Fix structures of input endpoints during GSM.",
)
def cli(
    input_paths: Sequence[Path],
    mep_mode: str,
    charge: Optional[int],
    ligand_charge: Optional[str],
    workers: int,
    workers_per_nodes: int,
    spin: Optional[int],
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    dump: bool,
    convert_files: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    preopt: bool,
    preopt_max_cycles: int,
    fix_ends: bool,
) -> None:
    input_paths = tuple(Path(p) for p in input_paths)
    set_convert_file_enabled(convert_files)
    prepared_inputs = [prepare_input_structure(p) for p in input_paths]
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble final config (defaults ← YAML ← CLI)
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        dmf_cfg = dict(DMF_KW)
        gs_cfg = dict(GS_KW)
        opt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(_LBFGS_KW)
        rfo_cfg = dict(_RFO_KW)

        # CLI overrides (defaults ← CLI)
        resolved_charge = charge
        resolved_spin = spin
        for prepared in prepared_inputs:
            resolved_charge, resolved_spin = fill_charge_spin_from_gjf(
                resolved_charge, resolved_spin, prepared.gjf_template
            )
        if resolved_charge is None and ligand_charge is not None:
            resolved_charge = _derive_charge_from_ligand_charge(
                prepared_inputs[0], ligand_charge, prefix="[path-opt]"
            )
        if resolved_charge is None:
            resolved_charge = 0
        if resolved_spin is None:
            resolved_spin = 1
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)
        calc_cfg["workers"] = int(workers)
        calc_cfg["workers_per_nodes"] = int(workers_per_nodes)

        gs_cfg["max_nodes"] = int(max_nodes)
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)
        gs_cfg["fix_first"] = bool(fix_ends)
        gs_cfg["fix_last"] = bool(fix_ends)

        # Lanczos tangent estimation follows the CLI --climb flag by default but
        # can still be overridden via YAML (`gs.climb_lanczos`).
        opt_cfg["dump"] = bool(dump)
        opt_cfg["out_dir"] = out_dir
        if thresh is not None:
            opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["thresh"] = str(thresh)
            rfo_cfg["thresh"] = str(thresh)

        # Use external Kabsch alignment; keep internal align disabled.
        opt_cfg["align"] = False

        lbfgs_cfg["dump"] = bool(dump)
        rfo_cfg["dump"] = bool(dump)
        lbfgs_cfg["out_dir"] = out_dir
        rfo_cfg["out_dir"] = out_dir

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (dmf_cfg, (("dmf",),)),
                (gs_cfg, (("gs",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("sopt", "lbfgs"), ("opt", "lbfgs"), ("lbfgs",))),
                (rfo_cfg, (("sopt", "rfo"), ("opt", "rfo"), ("rfo",))),
            ],
        )

        opt_kind = opt_mode.strip().lower()
        mep_mode_kind = mep_mode.strip().lower()
        if opt_kind == "light":
            sopt_kind = "lbfgs"
            sopt_cfg = lbfgs_cfg
        elif opt_kind == "heavy":
            sopt_kind = "rfo"
            sopt_cfg = rfo_cfg
        else:
            raise click.BadParameter(f"Unknown --opt-mode '{opt_mode}'.")

        sopt_cfg = dict(sopt_cfg)
        sopt_cfg["max_cycles"] = int(preopt_max_cycles)

        # For display: resolved configuration
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = format_freeze_atoms_for_echo(calc_cfg)
        echo_gs = dict(gs_cfg)
        echo_opt = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs", echo_gs))
        click.echo(pretty_block("opt", echo_opt))
        if mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg))
        click.echo(pretty_block("sopt." + sopt_kind, sopt_cfg))
        click.echo(
            pretty_block(
                "run_flags",
                {
                    "preopt": bool(preopt),
                    "preopt_max_cycles": int(preopt_max_cycles),
                    "mep_mode": mep_mode_kind,
                },
            )
        )

        # --------------------------
        # 2) Prepare structures (load two endpoints and apply freezing)
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        source_paths = [prep.source_path for prep in prepared_inputs]

        geoms = _load_two_endpoints(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"]),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )

        # Shared UMA calculator (reuse the same instance for all images)
        shared_calc = uma_pysis(**calc_cfg)

        # Optional endpoint pre-optimization (LBFGS/RFO) before alignment/GSM
        if preopt:
            click.echo(
                "\n=== Preoptimizing endpoints via single-structure optimizer ===\n"
            )
            ref_pdb_for_preopt: Optional[Path] = None
            for p in source_paths:
                if p.suffix.lower() == ".pdb":
                    ref_pdb_for_preopt = p.resolve()
                    break

            preopt_out_dir = out_dir_path
            if (
                out_dir_path.name.startswith("seg_")
                and out_dir_path.parent.name == "path_opt"
            ):
                preopt_out_dir = out_dir_path.parent
                preopt_out_dir.mkdir(parents=True, exist_ok=True)

            new_geoms = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                try:
                    g_opt = _optimize_single(
                        g,
                        shared_calc,
                        sopt_kind,
                        sopt_cfg,
                        preopt_out_dir,
                        tag=tag,
                        prepared_input=prepared_inputs[i] if i < len(prepared_inputs) else None,
                    )
                    new_geoms.append(g_opt)
                except Exception as e:
                    click.echo(
                        f"[preopt] WARNING: Failed to preoptimize endpoint {i}: {e}",
                        err=True,
                    )
                    new_geoms.append(g)
            geoms = new_geoms
        else:
            click.echo(
                "[preopt] Skipping endpoint preoptimization (use --preopt True to enable)."
            )

        # External Kabsch alignment (if freeze_atoms exist, use only them)
        align_thresh = str(opt_cfg.get("thresh", "gau"))
        try:
            click.echo(
                "\n=== Aligning all inputs to the first structure "
                "(freeze-guided scan + relaxation) ===\n"
            )
            _ = align_and_refine_sequence_inplace(
                geoms,
                thresh=align_thresh,
                shared_calc=shared_calc,
                out_dir=out_dir_path / "align_refine",
                verbose=True,
            )
            click.echo("[align] Completed input alignment.")
        except Exception as e:
            click.echo(f"[align] WARNING: alignment skipped: {e}", err=True)

        fix_atoms: List[int] = []
        try:
            fix_atoms = sorted(
                {int(i) for g in geoms for i in getattr(g, "freeze_atoms", [])}
            )
        except Exception:
            pass

        if mep_mode_kind == "dmf":
            try:
                dmf_res = _run_dmf_mep(
                    geoms,
                    calc_cfg,
                    out_dir_path,
                    prepared_inputs,
                    max_nodes,
                    fix_atoms,
                    dmf_cfg=dmf_cfg,
                )
            except Exception as e:
                click.echo(f"[dmf] ERROR: DMF optimization failed: {e}", err=True)
                sys.exit(3)

            try:
                hei_idx = int(dmf_res.hei_idx)
                hei_xyz = out_dir_path / "hei.xyz"
                _write_ase_trj_with_energy(
                    [dmf_res.images[hei_idx]], [dmf_res.energies[hei_idx]], hei_xyz
                )
                click.echo(f"[write] Wrote '{hei_xyz}'.")
                main_prepared = prepared_inputs[0] if prepared_inputs else None
                if main_prepared is not None:
                    ref_pdb = (
                        main_prepared.source_path.resolve()
                        if main_prepared.source_path.suffix.lower() == ".pdb"
                        else None
                    )
                    needs_pdb = ref_pdb is not None
                    needs_gjf = main_prepared.is_gjf
                    if needs_pdb or needs_gjf:
                        try:
                            convert_xyz_like_outputs(
                                hei_xyz,
                                main_prepared,
                                ref_pdb_path=ref_pdb,
                                out_pdb_path=out_dir_path / "hei.pdb" if needs_pdb else None,
                                out_gjf_path=out_dir_path / "hei.gjf" if needs_gjf else None,
                            )
                            click.echo("[convert] Wrote 'hei' outputs.")
                        except Exception as e:
                            click.echo(
                                f"[convert] WARNING: Failed to convert HEI to requested formats: {e}",
                                err=True,
                            )
            except Exception as e:
                click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
                sys.exit(5)

            click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start))
            return

        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # Used when GrowingString generates new nodes
            return shared_calc

        # --------------------------
        # 3) Build path object and optimizer
        # --------------------------
        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        opt_args = dict(opt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"},
        )

        # --------------------------
        # 4) Run optimization
        # --------------------------
        click.echo("\n=== Growing String optimization started ===\n")
        optimizer.run()
        click.echo("\n=== Growing String optimization finished ===\n")

        # --------------------------
        # 5) Write final path (final_geometries.trj)
        # --------------------------
        final_trj = out_dir_path / "final_geometries.trj"
        try:
            try:
                energies = np.array(gs.energy, dtype=float)
                blocks = []
                for idx, (geom, E) in enumerate(zip(gs.images, energies)):
                    s = geom.as_xyz()
                    lines = s.splitlines()
                    if len(lines) >= 2 and lines[0].strip().isdigit():
                        lines[1] = f"{E:.12f}"
                    s_mod = "\n".join(lines)
                    if not s_mod.endswith("\n"):
                        s_mod += "\n"
                    blocks.append(s_mod)
                annotated = "".join(blocks)
                with open(final_trj, "w") as f:
                    f.write(annotated)
                click.echo(f"[write] Wrote '{final_trj}' with energy.")
            except Exception:
                with open(final_trj, "w") as f:
                    f.write(gs.as_xyz())
                click.echo(f"[write] Wrote '{final_trj}'.")

            main_prepared = prepared_inputs[0]
            needs_pdb = main_prepared.source_path.suffix.lower() == ".pdb"
            needs_gjf = main_prepared.is_gjf
            ref_pdb = main_prepared.source_path.resolve() if needs_pdb else None
            if needs_pdb or needs_gjf:
                try:
                    convert_xyz_like_outputs(
                        final_trj,
                        main_prepared,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=out_dir_path / "final_geometries.pdb" if needs_pdb else None,
                        out_gjf_path=out_dir_path / "final_geometries.gjf" if needs_gjf else None,
                    )
                    click.echo("[convert] Wrote 'final_geometries' outputs.")
                except Exception as e:
                    click.echo(
                        f"[convert] WARNING: Failed to convert MEP path trajectory: {e}",
                        err=True,
                    )

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        # --------------------------
        # 6) Identify and write HEI (hei.xyz[/pdb/.gjf])
        # --------------------------
        try:
            energies = np.array(gs.energy, dtype=float)
            hei_idx = _select_hei_index(energies)

            hei_geom = gs.images[int(hei_idx)]
            hei_E = float(energies[int(hei_idx)])

            hei_xyz = out_dir_path / "hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            main_prepared = prepared_inputs[0]
            needs_pdb = main_prepared.source_path.suffix.lower() == ".pdb"
            needs_gjf = main_prepared.is_gjf
            ref_pdb = main_prepared.source_path.resolve() if needs_pdb else None
            if needs_pdb or needs_gjf:
                try:
                    convert_xyz_like_outputs(
                        hei_xyz,
                        main_prepared,
                        ref_pdb_path=ref_pdb,
                        out_pdb_path=out_dir_path / "hei.pdb" if needs_pdb else None,
                        out_gjf_path=out_dir_path / "hei.gjf" if needs_gjf else None,
                    )
                    click.echo("[convert] Wrote 'hei' outputs.")
                except Exception as e:
                    click.echo(
                        f"[convert] WARNING: Failed to convert HEI structure: {e}",
                        err=True,
                    )
            else:
                click.echo("[convert] Skipped HEI conversion (no PDB/GJF template).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

        click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start))

    except OptimizationError as e:
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo(
            "Unhandled error during path optimization:\n"
            + textwrap.indent(tb, "  "),
            err=True,
        )
        sys.exit(1)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()


def freeze_links_helper(pdb_path: Path):
    """
    Expose detect_freeze_links for external callers/tests.
    """
    return detect_freeze_links(pdb_path)
