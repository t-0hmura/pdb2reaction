# pdb2reaction/path_opt.py

"""
path_opt — Minimum-energy path (MEP) optimization via the Growing String method (pysisyphus) with a UMA calculator
====================================================================

Usage
-----
    pdb2reaction path_opt -i REACTANT.{pdb|xyz} PRODUCT.{pdb|xyz} \
        -q <charge> [-s <multiplicity>] [--freeze-links {True|False}] \
        [--max-nodes <int>] [--max-cycles <int>] [--climb {True|False}] \
        [--dump {True|False}] [--out-dir <dir>] [--args-yaml <file>]

Examples::
    # Minimal: two endpoints, neutral singlet
    pdb2reaction path_opt -i reac.pdb prod.pdb -q 0 -s 1

    # Typical full run with YAML overrides and dumps
    pdb2reaction path_opt -i reac.pdb prod.pdb -q 0 -s 1 \
      --freeze-links True --max-nodes 10 --max-cycles 100 \
      --dump True --out-dir ./result_path_opt/ --args-yaml ./args.yaml


Description
-----
- Optimize a minimum-energy path between two endpoints using pysisyphus `GrowingString` and `StringOptimizer`, with UMA as the calculator (via `uma_pysis`).
- Inputs: two structures (.pdb or .xyz). If a PDB is provided and `--freeze-links=True` (default), parent atoms of link hydrogens are added to `freeze_atoms` (0-based indices).
- Configuration via YAML with sections {geom, calc, gs, opt}. Precedence: CLI > YAML > built-in defaults.
- Default alignment: before optimization, the second endpoint is rigidly Kabsch-aligned to the first using an external routine. `StringOptimizer.align` is disabled (known to be fragile). If either endpoint specifies `freeze_atoms`, the RMSD fit uses only those atoms and the resulting rigid transform is applied to all atoms.
- With `--climb=True` (default), a climbing-image step refines the highest-energy image; Lanczos-based tangent estimation is enabled by default (configurable via YAML).
- After optimization, the highest-energy image (HEI) is identified and exported.

Outputs (& Directory Layout)
-----
out_dir/ (default: ./result_path_opt/)
  ├─ final_geometries.trj        # XYZ trajectory of the optimized path; comment line holds per-image energy
  ├─ final_geometries.pdb        # Converted from .trj when a PDB reference is available
  ├─ gsm_hei.xyz                 # Highest-energy image with energy on the comment line
  ├─ gsm_hei.pdb                 # HEI converted to PDB when a PDB reference is available
  ├─ align_refine/               # Files from external alignment/refinement
  └─ <optimizer dumps / restarts>  # Emitted when --dump=True (from pysisyphus)

Notes:
-----
- Always set correct total charge (`-q/--charge`) and spin multiplicity (`-s/--spin`) to avoid unphysical conditions.
- Coordinates are Cartesian; `freeze_atoms` use 0-based indices. With `--freeze-links=True` and PDB inputs, link-hydrogen parents are added automatically.
- `--max-nodes` sets the number of internal nodes; the string has (max_nodes + 2) images including endpoints.
- `--max-cycles` limits optimization; after full growth, the same bound applies to additional refinement.
- Exit codes: 0 (success); 3 (optimization failed); 4 (final trajectory write error); 5 (HEI dump error); 130 (keyboard interrupt); 1 (unhandled error).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import sys
import traceback
import textwrap

import click
import numpy as np
import yaml
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError

from .uma_pysis import uma_pysis
from .utils import (
    convert_xyz_to_pdb,
    freeze_links as detect_freeze_links,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_geom_for_echo,
    merge_freeze_atom_indices,
)
from .align_freeze_atoms import align_and_refine_sequence_inplace


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
    "stop_in_when_full": 100,   # After fully grown, stop after N additional cycles
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

def _freeze_links_for_pdb(pdb_path: Path) -> Sequence[int]:
    """Detect the parent atoms of link hydrogens in a PDB and return 0‑based indices."""
    try:
        return detect_freeze_links(pdb_path)
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
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        if auto_freeze_links and p.suffix.lower() == ".pdb":
            detected = _freeze_links_for_pdb(p)
            freeze = merge_freeze_atom_indices(cfg, detected)
            if detected and freeze:
                click.echo(f"[freeze-links] {p.name}: Freeze atoms (0-based): {','.join(map(str, freeze))}")
        else:
            freeze = merge_freeze_atom_indices(cfg)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


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
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble final config (defaults ← YAML ← CLI)
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg   = dict(GS_KW)
        opt_cfg  = dict(OPT_KW)

        # Prefer centralized override helper (imports from utils)
        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (gs_cfg,   (("gs",),)),
                (opt_cfg,  (("opt",),)),
            ],
        )

        # CLI overrides (highest precedence)
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
        echo_geom = format_geom_for_echo(geom_cfg)
        echo_calc = dict(calc_cfg)
        echo_gs   = dict(gs_cfg)
        echo_opt  = dict(opt_cfg)
        echo_opt["out_dir"] = str(out_dir_path)

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs",   echo_gs))
        click.echo(pretty_block("opt",  echo_opt))

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

        # Shared UMA calculator (reuse the same instance for all images)
        shared_calc = uma_pysis(**calc_cfg)

        # By default, apply external Kabsch alignment (if freeze_atoms exist, use only them)
        try:
            click.echo("\n=== Aligning all inputs to the first structure (freeze-guided scan + relaxation) ===\n")
            _ = align_and_refine_sequence_inplace(
                geoms,
                shared_calc=shared_calc,
                out_dir=out_dir_path / "align_refine",
                verbose=True,
            )
            click.echo("[align] Completed input alignment.")
        except Exception as e:
            click.echo(f"[align] WARNING: alignment skipped: {e}", err=True)

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

        elapsed = time.perf_counter() - time_start
        hh = int(elapsed // 3600)
        mm = int((elapsed % 3600) // 60)
        ss = elapsed - (hh * 3600 + mm * 60)
        click.echo(f"[time] Elapsed Time for Path Opt: {hh:02d}:{mm:02d}:{ss:06.3f}")

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
