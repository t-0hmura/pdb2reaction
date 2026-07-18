"""
bond-summary — Detect and summarize bond changes between molecular structures.

Usage:
    pdb2reaction bond-summary -i A.xyz B.xyz [C.xyz ...]

Compares consecutive pairs of structures (A→B, B→C, …) and reports
covalent bonds formed and broken. Supports XYZ, PDB, mmCIF, and GJF formats.

For detailed documentation, see: docs/bond-summary.md
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import click
from pdb2reaction.core.output import emit

from pdb2reaction.domain.bond_changes import compare_structures, summarize_changes


def _load_geom(path: str):
    """Load a geometry from XYZ, PDB, mmCIF, or GJF file."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix in {".pdb", ".cif", ".mmcif", ".gjf"}:
        from pdb2reaction.core.utils import prepare_input_structure
        with prepare_input_structure(p) as prepared:
            if suffix == ".gjf":
                # The existing p2r preparation path is the canonical ordinary
                # Gaussian-input parser and owns its temporary XYZ bridge.
                from pysisyphus.helpers import geom_from_xyz_file

                return geom_from_xyz_file(str(prepared.geom_path))
            from pysisyphus.io.pdb import geom_from_pdb

            return geom_from_pdb(str(prepared.geom_path))
    else:
        from pysisyphus.helpers import geom_from_xyz_file
        return geom_from_xyz_file(str(p))


@click.command("bond-summary", help="Detect bond changes between consecutive structures.")
@click.option(
    "-i", "--input", "inputs", multiple=True,
    help="Input structure files (XYZ/PDB/mmCIF/GJF). Repeat -i for each file.",
)
@click.argument("extra_inputs", nargs=-1, type=click.Path(exists=True))
@click.option(
    "--device", default="cpu", show_default=True,
    help="Compute device for distance calculations.",
)
@click.option(
    "--bond-factor", default=1.20, show_default=True,
    help="Scaling factor for covalent radii sum.",
)
@click.option(
    "--one-based/--zero-based", default=True, show_default=True,
    help="Use 1-based atom indices in output.",
)
@click.option(
    "--json/--no-json",
    "out_json",
    default=False,
    show_default=True,
    help="Emit JSON to stdout (no file written).",
)
def cli(inputs: tuple, extra_inputs: tuple, device: str, bond_factor: float, one_based: bool, out_json: bool) -> None:
    """Detect and summarize bond changes between consecutive structure pairs.

    \b
    Usage:
      pdb2reaction bond-summary -i A.xyz -i B.xyz -i C.xyz
      pdb2reaction bond-summary A.xyz B.xyz C.xyz
      pdb2reaction bond-summary -i A.xyz B.xyz C.xyz
    """
    files: List[str] = list(inputs) + list(extra_inputs)
    if len(files) < 2:
        raise click.BadParameter("At least two input files are required.", param_hint="'-i'")

    geoms = []
    for f in files:
        p = Path(f)
        if not p.exists():
            raise click.FileError(f, hint="File not found.")
        geoms.append((p.name, _load_geom(f)))

    comparisons_json: List[dict] = []
    # Track per-pair failures so the exit code and JSON envelope status reflect
    # any dropped comparison instead of silently reporting "ok" / exit 0.
    n_pairs = 0
    n_failed = 0

    for idx in range(len(geoms) - 1):
        name1, g1 = geoms[idx]
        name2, g2 = geoms[idx + 1]
        n_pairs += 1

        if not out_json:
            click.echo(f"{'=' * 60}")
            click.echo(f"  {name1}  →  {name2}")
            click.echo(f"{'=' * 60}")

        try:
            result = compare_structures(g1, g2, device=device, bond_factor=bond_factor)
            if out_json:
                comparisons_json.append({
                    "structure_a": name1,
                    "structure_b": name2,
                    "bonds_formed": len(result.formed_covalent),
                    "bonds_broken": len(result.broken_covalent),
                })
            else:
                summary = summarize_changes(g2, result, one_based=one_based)
                click.echo(summary)
        except AssertionError as e:
            n_failed += 1
            click.echo(f"  ERROR: {e}", err=True)
        except Exception as e:
            n_failed += 1
            click.echo(f"  ERROR: {e}", err=True)

        if not out_json:
            click.echo()

    # Envelope status: "ok" when every pair compared cleanly; "partial" when at
    # least one pair succeeded but others failed; "failed" when no pair compared.
    if n_failed == 0:
        status = "ok"
    elif n_failed < n_pairs:
        status = "partial"
    else:
        status = "failed"

    if out_json:
        import json as _json
        # `--json` is an explicit machine-readable deliverable: it must always
        # reach stdout, even at default verbosity (which otherwise suppresses
        # detail console output via the core.utils click.echo chokepoint).
        # `force=True` bypasses the verbosity gate; see
        # pdb2reaction.core.utils._patch_click_echo.
        emit(
            _json.dumps(
                {"status": status, "comparisons": comparisons_json},
                indent=2, ensure_ascii=False,
            ),
            force=True,
        )

    # Surface dropped comparisons via a non-zero exit code so callers (shell /
    # MCP) cannot mistake a partial/failed run for success.
    if n_failed > 0:
        raise SystemExit(1)
