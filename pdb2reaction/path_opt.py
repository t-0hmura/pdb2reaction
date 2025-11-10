# pdb2reaction/path_opt.py

"""
MEP optimization via the Growing String method.

- Endpoints: two structures (e.g., reactant/product) given by -i/--input.
- Calculator: UMA (through the uma_pysis wrapper).
- Path optimizer: GrowingString + StringOptimizer (pysisyphus).
- Configuration: YAML with sections {geom, calc, gs, opt}. Precedence: CLI > YAML > defaults.
- If inputs are PDB and --freeze-links=True, freeze the parent atoms of link hydrogens.
- **Default behavior:** Before optimization, rigidly align the second endpoint to the first
  using an external Kabsch alignment. Do **not** use StringOptimizer’s internal `align` option
  (known to be fragile).
- **If either endpoint specifies `freeze_atoms`, the alignment is determined using *only* those
  atoms (RMSD over `freeze_atoms`), and the resulting rigid transform is applied to *all* atoms.**
- Dump the highest-energy image (HEI) as gsm_hei.xyz and gsm_hei.pdb (if a PDB reference is available).

Examples
--------
pdb2reaction path_opt -i reac.pdb prod.pdb -q 0 -s 1 \
  --freeze-links True --max-nodes 10 --max-cycles 100 \
  --dump True --dump_dir ./result_path_opt/ --args-yaml ./args.yaml
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import sys
import traceback
import textwrap

import click
import numpy as np
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError

from .uma_pysis import uma_pysis
from .utils import convert_xyz_to_pdb, freeze_links


# -----------------------------------------------
# Defaults (overridden by YAML/CLI)
# -----------------------------------------------

# Geometry (input handling)
GEOM_KW: Dict[str, Any] = {
    "coord_type": "cart",   # GrowingString works best with Cartesian coordinates
    "freeze_atoms": [],     # 0-based atom indices
}

# UMA calculator settings
CALC_KW: Dict[str, Any] = {
    "charge": 0,
    "spin": 1,                  # multiplicity (= 2S + 1)
    "model": "uma-s-1p1",
    "task_name": "omol",
    "device": "auto",
    "max_neigh": None,
    "radius": None,
    "r_edges": False,
}

# GrowingString (path representation)
GS_KW: Dict[str, Any] = {
    "max_nodes": 30,            # Including endpoints, the string has max_nodes + 2 images
    "perp_thresh": 5e-3,        # Frontier growth criterion (RMS/NORM of perpendicular force)
    "reparam_check": "rms",     # "rms" | "norm"; convergence check after reparam (RMS of structural change)
    "reparam_every": 1,         # Reparametrize every N steps
    "reparam_every_full": 1,    # After the path is fully grown, reparametrize every N steps
    "param": "equi",            # "equi" (even spacing) | "energy" (denser near the peak via weighting)
    "max_micro_cycles": 10,
    "reset_dlc": True,          # Reset DLC coordinates when appropriate
    "climb": True,              # Enable climbing image
    "climb_rms": 5e-4,          # Threshold (RMS of force) to start CI
    "climb_lanczos": True,      # Use Lanczos to estimate the HEI tangent
    "climb_lanczos_rms": 5e-4,
    "climb_fixed": False,       # Fix the HEI image index
    "scheduler": None,          # Serial execution (assumes a shared calculator instance)
}

# StringOptimizer (optimization control)
OPT_KW: Dict[str, Any] = {
    "type": "string",           # Tag for bookkeeping
    "stop_in_when_full": 100,  # After fully grown, stop after N additional cycles
    "align": False,             # Keep internal align disabled; use external Kabsch alignment instead
    "scale_step": "global",     # "global" | "per_image"
    "max_cycles": 100,
    "dump": False,
    "dump_restart": False,
    "reparam_thresh": 1e-3,     # Convergence after reparam (RMS of step)
    "coord_diff_thresh": 0.0,   # Near-duplicate image check (0 disables)
    "out_dir": "./result_path_opt/",
    "print_every": 1,
}


# -----------------------------------------------
# Utilities
# -----------------------------------------------

def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively update dict *dst* with *src*, returning *dst*."""
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    """Load a YAML mapping from *path*. Return an empty dict if *path* is None."""
    if not path:
        return {}
    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")
    return data


def _pretty_block(title: str, content: Dict[str, Any]) -> str:
    """Format a titled YAML block for echoing resolved configuration."""
    body = yaml.safe_dump(content, sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _format_geom_for_echo(geom_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Pretty-print helper for geom config (joins freeze_atoms as CSV)."""
    g = dict(geom_cfg)
    fa = g.get("freeze_atoms")
    if isinstance(fa, (list, tuple, np.ndarray)):
        g["freeze_atoms"] = ",".join(map(str, fa)) if len(fa) else ""
    return g


def _freeze_links_for_pdb(pdb_path: Path) -> Sequence[int]:
    """Detect the parent atoms of link hydrogens in a PDB and return 0‑based indices."""
    try:
        return freeze_links(pdb_path)
    except Exception as e:
        click.echo(f"[freeze-links] WARNING: Could not detect link parents for '{pdb_path.name}': {e}", err=True)
        return []


def _load_two_endpoints(
    paths: Sequence[Path],
    coord_type: str,
    base_freeze: Sequence[int],
    auto_freeze_links: bool,
) -> Sequence:
    """Load the two endpoint structures and set `freeze_atoms` as needed."""
    geoms = []
    for p in paths:
        g = geom_loader(p, coord_type=coord_type)
        # Attach freeze_atoms to the Geometry object
        freeze = list(base_freeze)
        if auto_freeze_links and p.suffix.lower() == ".pdb":
            detected = _freeze_links_for_pdb(p)
            if detected:
                freeze = sorted(set(freeze).union(detected))
                click.echo(f"[freeze-links] {p.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# ---------- Kabsch alignment (external, robust, minimal) ----------

def _kabsch_R_t(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute optimal rotation R and translation t that aligns Q onto P
    (minimizing ||P - (Q R + t)||_F) via Kabsch.

    Parameters
    ----------
    P, Q : (N, 3) arrays

    Returns
    -------
    R : (3, 3) rotation matrix
    t : (3,) translation vector
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("Kabsch expects P, Q with shape (N, 3).")
    mu_P = P.mean(axis=0)
    mu_Q = Q.mean(axis=0)
    Pc = P - mu_P
    Qc = Q - mu_Q
    # NOTE: covariance order is Pc.T @ Qc (maps Q -> P). Do not swap!
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # Ensure right-handed rotation
    if np.linalg.det(R) < 0.0:
        Vt[-1, :] *= -1.0
        R = Vt.T @ U.T
    t = mu_P - mu_Q @ R
    return R, t


def _align_second_to_first_kabsch(geom_ref, geom_to_align) -> Tuple[float, float, int]:
    """
    Rigidly align geom_to_align onto geom_ref using Kabsch (no scaling).

    - If either endpoint provides `freeze_atoms`, determine the transform using the
      union of those atoms only, then apply the same rigid transform to *all* atoms.
      Otherwise, use all atoms to determine the transform.

    Returns
    -------
    rmsd_before, rmsd_after, n_used
    """
    P = np.array(geom_ref.coords3d, dtype=float)    # (N, 3)
    Q = np.array(geom_to_align.coords3d, dtype=float)
    if P.shape != Q.shape:
        raise ValueError(f"Different atom counts for endpoints: {P.shape[0]} vs {Q.shape[0]}")

    N = P.shape[0]
    # Indices to use for alignment
    fa0 = getattr(geom_ref, "freeze_atoms", np.array([], dtype=int))
    fa1 = getattr(geom_to_align, "freeze_atoms", np.array([], dtype=int))
    freeze_union = sorted(set(map(int, fa0)) | set(map(int, fa1)))

    if len(freeze_union) > 0:
        # Use ONLY the frozen atoms to determine the rigid transform
        use_mask = np.zeros(N, dtype=bool)
        # Clip any out-of-range indices defensively
        valid_idx = [i for i in freeze_union if 0 <= i < N]
        use_mask[valid_idx] = True
    else:
        # No frozen atoms specified -> use all atoms
        use_mask = np.ones(N, dtype=bool)

    P_sel = P[use_mask]
    Q_sel = Q[use_mask]
    n_used = int(P_sel.shape[0])

    # Pre-align RMSD (on selected set)
    def _rmsd(A: np.ndarray, B: np.ndarray) -> float:
        return float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1)))) if len(A) else float("nan")

    rmsd_before = _rmsd(P_sel, Q_sel)

    R, t = _kabsch_R_t(P_sel, Q_sel)
    Q_aligned = (Q @ R) + t  # apply to all atoms

    # Write back (expects flattened 1D array of length 3N)
    geom_to_align.set_coords(Q_aligned.reshape(-1))

    rmsd_after = _rmsd(P[use_mask], Q_aligned[use_mask])

    return rmsd_before, rmsd_after, n_used


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
@click.option("-q", "--charge", type=int, required=True, help="Total charge.")
@click.option("-s", "--spin", type=int, default=1, show_default=True, help="Spin multiplicity (2S+1).")
@click.option("--freeze-links", "freeze_links_flag", type=click.BOOL, default=True, show_default=True,
              help="If a PDB is provided, freeze the parent atoms of link hydrogens.")
@click.option("--max-nodes", type=int, default=30, show_default=True,
              help="Number of internal nodes (string has max_nodes+2 images including endpoints).")
@click.option("--max-cycles", type=int, default=1000, show_default=True, help="Maximum optimization cycles.")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Search for a transition state (climbing image) after path growth.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump optimizer trajectory/restarts during the run.")
@click.option("--out-dir", "out_dir", type=str, default="./result_path_opt/", show_default=True,
              help="Output directory.")
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML with extra args (sections: geom, calc, gs, opt).",
)
def cli(
    input_paths: Sequence[Path],
    charge: int,
    spin: int,
    freeze_links_flag: bool,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    dump: bool,
    out_dir: str,
    args_yaml: Optional[Path],
) -> None:
    try:
        # --------------------------
        # 1) Assemble final config (defaults ← YAML ← CLI)
        # --------------------------
        yaml_cfg = _load_yaml(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg   = dict(GS_KW)
        opt_cfg  = dict(OPT_KW)

        _deep_update(geom_cfg, yaml_cfg.get("geom", {}))
        _deep_update(calc_cfg, yaml_cfg.get("calc", {}))
        _deep_update(gs_cfg,   yaml_cfg.get("gs",   {}))
        _deep_update(opt_cfg,  yaml_cfg.get("opt",  {}))

        # CLI overrides
        calc_cfg["charge"] = int(charge)
        calc_cfg["spin"]   = int(spin)

        gs_cfg["max_nodes"] = int(max_nodes)
        opt_cfg["max_cycles"] = int(max_cycles)
        opt_cfg["stop_in_when_full"] = int(max_cycles)
        gs_cfg["climb"] = bool(climb)
        gs_cfg["climb_lanczos"] = bool(climb)

        opt_cfg["dump"]       = bool(dump)
        opt_cfg["out_dir"]    = out_dir  # Pass --out-dir to the optimizer via "out_dir"

        # Important: do not use internal alignment; use external Kabsch alignment instead
        opt_cfg["align"] = False

        # For display: resolved configuration
        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        echo_geom = _format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_opt  = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        click.echo(_pretty_block("geom", echo_geom))
        click.echo(_pretty_block("calc", echo_calc))
        click.echo(_pretty_block("gs",   echo_gs))
        click.echo(_pretty_block("opt",  echo_opt))

        # --------------------------
        # 2) Prepare structures (load two endpoints and apply freezing)
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Input paths
        p0, p1 = Path(input_paths[0]), Path(input_paths[1])

        # Load endpoints (if PDB, merge in link-parent freezing)
        geoms = _load_two_endpoints(
            paths=[p0, p1],
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
            auto_freeze_links=bool(freeze_links_flag),
        )

        # By default, apply external Kabsch alignment (if freeze_atoms exist, use only them)
        try:
            rmsd_before, rmsd_after, n_used = _align_second_to_first_kabsch(geoms[0], geoms[1])
            click.echo(f"[align] Applied Kabsch alignment: 2nd endpoint -> 1st endpoint "
                       f"(used {n_used} atoms; RMSD {rmsd_before:.4f} Å -> {rmsd_after:.4f} Å).")
        except Exception as e:
            click.echo(f"[align] WARNING: Kabsch alignment skipped: {e}", err=True)

        # Shared UMA calculator (reuse the same instance for all images)
        shared_calc = uma_pysis(**calc_cfg)
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

        # StringOptimizer expects 'out_dir' under the key "out_dir"
        opt_args = dict(opt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"}  # 'type' is just a tag
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

            if input_paths[0].suffix.lower() == ".pdb":
                ref_pdb = input_paths[0].resolve()

                try:
                    out_pdb = out_dir_path / "final_geometries.pdb"
                    convert_xyz_to_pdb(final_trj, ref_pdb, out_pdb)
                    click.echo(f"[convert] Wrote '{out_pdb}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert MEP path trajectory to PDB: {e}", err=True)

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        try:
            energies = np.array(gs.energy, dtype=float)
            # --- HEI identification logic ---
            # Choose the internal local maximum (exclude endpoints) with the highest energy,
            # i.e., nodes whose immediate neighbors have lower energy.
            # Fallback 1: if none exist, pick the maximum among internal nodes (exclude endpoints).
            # Fallback 2: if internal nodes are unavailable, pick the global maximum.
            nE = int(len(energies))
            hei_idx = None
            if nE >= 3:
                # Strict internal local maxima (both neighbors lower)
                candidates = [i for i in range(1, nE - 1)
                              if energies[i] > energies[i - 1] and energies[i] > energies[i + 1]]
                if candidates:
                    cand_es = energies[candidates]
                    rel = int(np.argmax(cand_es))
                    hei_idx = int(candidates[rel])
                else:
                    # Fallback 1: maximum over internal nodes (exclude endpoints)
                    if nE > 2:
                        rel = int(np.argmax(energies[1:-1]))
                        hei_idx = 1 + rel
            if hei_idx is None:
                # Fallback 2: global maximum
                hei_idx = int(np.argmax(energies))

            hei_geom = gs.images[hei_idx]
            hei_E = float(energies[hei_idx])

            hei_xyz = out_dir_path / "gsm_hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            ref_pdb = None
            if input_paths[0].suffix.lower() == ".pdb":
                ref_pdb = input_paths[0].resolve()
            if ref_pdb is not None:
                hei_pdb = out_dir_path / "gsm_hei.pdb"
                convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
                click.echo(f"[convert] Wrote '{hei_pdb}'.")
            else:
                click.echo("[convert] Skipped 'gsm_hei.pdb' (no PDB reference among inputs).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

    except OptimizationError as e:
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during path optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


def freeze_links_helper(pdb_path: Path):
    """Expose freeze_links for external callers/tests."""
    return freeze_links(pdb_path)


if __name__ == "__main__":
    cli()
