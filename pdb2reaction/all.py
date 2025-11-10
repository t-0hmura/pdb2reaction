# pdb2reaction/all.py
# -*- coding: utf-8 -*-
"""
Pipeline driver:  **extract active‑site pockets → run MEP search → merge back to full systems**

Usage (example)
---------------
pdb2reaction all -i a.pdb b.pdb c.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \

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
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Optional

import sys
import math
import click
import time  # timing

# Local imports from the package
from .extract import extract_api
from . import path_search as _path_search


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
@click.option("--max-cycles", type=int, default=1000, show_default=True, help="Maximum GSM optimization cycles.")
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
    # path_search.cli re-parses sys.argv internally to support the style "-i A B C".
    # When calling it programmatically, we temporarily replace sys.argv so its internal
    # collector sees *our* args (pockets + ref templates), not the outer "all" CLI args.
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
    # Elapsed time
    # --------------------------
    elapsed = time.perf_counter() - time_start
    hh = int(elapsed // 3600)
    mm = int((elapsed % 3600) // 60)
    ss = elapsed - (hh * 3600 + mm * 60)
    click.echo(f"[all] Total Elapsed: {hh:02d}:{mm:02d}:{ss:06.3f}")


if __name__ == "__main__":
    cli()
