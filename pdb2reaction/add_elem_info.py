# pdb2reaction/add_elem_info.py

"""
add-elem-info — repair PDB element symbols (columns 77–78) with Biopython
===========================================================================

Usage (CLI)
-----------
    pdb2reaction add-elem-info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite {True|False}]

Examples
--------
    # Populate element fields and write to "<input>_add_elem.pdb"
    pdb2reaction add-elem-info -i 1abc.pdb

    # Write to a specific output file
    pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

    # Overwrite the input file in-place
    pdb2reaction add-elem-info -i 1abc.pdb --overwrite True

Output behavior
---------------------------
- If `-o/--out` is omitted and `--overwrite` is not `True`, write to `<input>_add_elem.pdb` (replace a
  trailing `.pdb` with `_add_elem.pdb`; otherwise append `_add_elem.pdb`).
- If `--overwrite True` is set without `-o/--out`, overwrite the input file in-place.
- When `-o/--out` is supplied, write to that path and ignore `--overwrite`.

Workflow
--------
- Parse the input with `Bio.PDB.PDBParser`, sharing residue definitions with `extract.py`
  (`AMINO_ACIDS`, `WATER_RES`, `ION`).
- Infer elements per atom using the atom name, residue name, and HETATM flag:
  - Ion residues from the `ION` dict: use residue-derived elements (polyatomic ions handled per
    atom; D* → H).
  - Proteins/nucleic acids/water: special handling for H/D, Se, and first-letter mapping for
    C/N/O/P/S; carbon sidechain labels default to C.
  - Other ligands: use atom-name prefixes (C*/P*, excluding CL) and fall back to element-symbol
    normalization (recognizing halogens and D → H).
- Write the structure through `PDBIO` and print a summary: totals for processed/assigned atoms,
  per-element counts, and up to 50 unresolved atoms.

Outputs
-------
- PDB with element columns 77–78 populated/corrected at the path determined above.
- Console report with totals, per-element counts, and truncated unresolved-atom list.

Notes
-----
- Only element columns are modified; coordinates, occupancies, B-factors, charges, altlocs,
  insertion codes, and record ordering stay untouched.
- Supports ATOM and HETATM records across all models/chains/residues.
- Depends on Biopython (`Bio.PDB`) and Click; deuterium labels map to hydrogen; selenium (`SE*`) and
  halogens are recognized automatically.
"""

from __future__ import annotations

import argparse
import collections
import os
import re
import sys
from pathlib import Path
from typing import Optional

import click
from Bio.PDB import PDBParser, PDBIO

# Reuse residue/ion dictionaries from extract.py to keep definitions in sync
from .extract import AMINO_ACIDS, ION, WATER_RES

# -----------------------------
# Element symbols (IUPAC, 1–118)
# -----------------------------
ELEMENTS: set[str] = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db",
    "Sg","Bh","Hs","Mt","Ds","Rg","Cn","Fl","Lv","Ts","Og"
}

# Common residue classes
PROTEIN_RES = set(AMINO_ACIDS.keys())
NUCLEIC_RES = {
    "DA","DT","DG","DC","DI",
    "A","U","G","C","I",
}

# -----------------------------
# Helper: normalize strings to element symbols
# -----------------------------
_re_letters = re.compile(r"[A-Za-z]+")

def _normalize_symbol(s: str) -> Optional[str]:
    """Remove non-letters; prefer a 2-letter match, then 1-letter, against known elements.
    Returns the correctly cased symbol if matched.
    Treat deuterium 'D' as hydrogen 'H' (PDB often uses D interchangeably with H).
    """
    if not s:
        return None
    m = _re_letters.findall(s)
    if not m:
        return None
    letters = "".join(m)
    if len(letters) >= 2:
        cand2 = (letters[:2][0].upper() + letters[:2][1].lower())
        if cand2 in ELEMENTS:
            return cand2
    cand1 = letters[0].upper()
    if cand1 in ELEMENTS:
        return cand1
    if letters[0].upper() == "D":
        return "H"
    return None

def _symbol_from_resname(resname: str) -> Optional[str]:
    """Extract an element symbol from an ion residue name (e.g., CA, FE2, Cl-, YB2)."""
    res = resname.strip()
    return _normalize_symbol(res)

def _default_out_pdb_path(in_pdb: str) -> str:
    """Default output path when -o/--out is omitted:
    replace trailing '.pdb' (case-insensitive) with '_add_elem.pdb';
    if no trailing '.pdb', append '_add_elem.pdb'.
    """
    p = Path(in_pdb)
    name = p.name
    if name.lower().endswith(".pdb"):
        name = name[:-4] + "_add_elem.pdb"
    else:
        name = name + "_add_elem.pdb"
    return str(p.with_name(name))

# -----------------------------
# Element inference (use residue to disambiguate)
# -----------------------------
def guess_element(atom_name: str, resname: str, is_het: bool) -> Optional[str]:
    """
    Infer the element from atom name + residue name.
    Priority:
      1) Ion residues: prefer the residue name (NH4 / H3O+ handled per-atom as H/N/O)
      2) Polymers (protein/nucleic acid) and water: follow convention (H/C/N/O/S/P/Se)
      3) Other ligands: use atom-name prefix; then normalization
      4) Unresolved → None
    """
    name_u = atom_name.strip().upper()
    res_u = resname.strip().upper()

    # 1) Ion residues — strongly prefer the residue-derived element
    if res_u in {k.upper() for k in ION.keys()}:
        # Polyatomic ions (NH4, H3O+, …): decide per atom name (treat D* as H)
        if name_u.startswith(("H", "D")):
            return "H"
        if name_u.startswith("N"):
            return "N"
        if name_u.startswith("O"):
            return "O"
        # Monatomic metals/halogens: from residue name
        sym = _symbol_from_resname(res_u)
        if sym:
            return sym

    # 2) Polymers (protein/nucleic) / water
    is_protein = res_u in PROTEIN_RES
    is_nucl = res_u in NUCLEIC_RES
    is_water = res_u in WATER_RES
    if is_protein or is_nucl or is_water:
        # Water: only O and H (treat D* as H)
        if is_water:
            if name_u.startswith(("H", "D")):
                return "H"
            return "O"

        # Hydrogen (including D*)
        if name_u.startswith(("H", "D")):
            return "H"

        # Selenium (e.g., selenomethionine/selenocysteine)
        if name_u.startswith("SE"):
            return "Se"

        # P, N, O, S map directly by first letter
        if name_u.startswith("P"):
            return "P"
        if name_u.startswith("N"):
            return "N"
        if name_u.startswith("O"):
            return "O"
        if name_u.startswith("S"):
            return "S"

        # Carbon for Cα/sidechain labels (CA, CB, CG, CD, CE, CZ, CH*, etc.)
        if name_u.startswith("C"):
            return "C"

        # Fallback to normalization (rare halogens or atypical labels)
        sym = _normalize_symbol(name_u)
        if sym:
            return sym

    # 3) Non-polymers (ligands / cofactors)
    if name_u.startswith("C") and not name_u.startswith("CL"):
        return "C"
    if name_u.startswith("P"):
        return "P"

    sym = _normalize_symbol(name_u)
    if sym:
        return sym

    # 4) Unresolved
    return None

def _get_atom_serial(atom) -> Optional[int]:
    """Safely obtain the serial number from a Biopython Atom, handling version differences."""
    sn = getattr(atom, "serial_number", None)
    if sn is None and hasattr(atom, "get_serial_number"):
        try:
            sn = atom.get_serial_number()
        except Exception:
            sn = None
    return sn

# -----------------------------
# Main processing
# -----------------------------
def assign_elements(in_pdb: str, out_pdb: Optional[str], overwrite: bool = False) -> None:
    # If an explicit output path is provided, never overwrite in-place even when --overwrite is
    # passed. This keeps -o/-\-out as the higher-priority choice.
    effective_overwrite = overwrite and out_pdb is None

    parser = PDBParser(QUIET=True)
    structure_id = os.path.splitext(os.path.basename(in_pdb))[0]
    structure = parser.get_structure(structure_id, in_pdb)

    total = 0
    assigned_or_updated = 0   # atoms whose inferred element differs from the parsed value
    unknown = []              # could not infer (left unchanged)

    by_element = collections.Counter()

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0].strip()  # '' (standard); 'W' = water; 'H_' = HETATM
                is_het = (hetflag != "")
                resname = residue.get_resname()
                for atom in residue:
                    total += 1
                    name = atom.get_name()

                    sym = guess_element(name, resname, is_het)
                    if sym is None:
                        serial = _get_atom_serial(atom)
                        unknown.append((model.id, chain.id, residue.id, resname, name, serial))
                        continue

                    prev = getattr(atom, "element", None)
                    atom.element = sym
                    by_element[sym] += 1
                    if prev != sym:
                        assigned_or_updated += 1

    io = PDBIO()
    io.set_structure(structure)
    if effective_overwrite:
        out_path = in_pdb
    else:
        out_path = out_pdb if out_pdb else _default_out_pdb_path(in_pdb)
    io.save(out_path)

    # Summary
    print(f"[OK] Wrote: {out_path}")
    print(f"  total atoms                 : {total}")
    print(f"  assigned/updated            : {assigned_or_updated}")
    if by_element:
        top = ", ".join(f"{k}:{v}" for k, v in by_element.most_common())
        print(f"  assignment breakdown        : {top}")
    if unknown:
        print(f"[WARN] Could not confidently assign {len(unknown)} atoms; left unchanged.")
        for (mid, chid, resid, resn, aname, serial) in unknown[:50]:
            if isinstance(resid, tuple):
                resseq = resid[1]
                icode = resid[2].strip()
            else:
                resseq, icode = "?", ""
            s_str = f" serial {serial}" if serial is not None else ""
            print(f"    model {mid} chain {chid} {resn} {resseq}{icode} : {aname}{s_str}")
    if len(unknown) > 50:
        print("    ... (truncated) ...")


def _parse_bool(value: str) -> bool:
    value_lower = value.strip().lower()
    if value_lower in {"true", "1", "yes", "y", "t"}:
        return True
    if value_lower in {"false", "0", "no", "n", "f"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value!r}. Use True/False.")


def main():
    ap = argparse.ArgumentParser(
        description="Add/repair element columns (77–78) in a PDB using Biopython."
    )
    ap.add_argument(
        "-i",
        "--input",
        dest="in_pdb",
        required=True,
        help="Input PDB filepath",
    )
    ap.add_argument(
        "-o",
        "--out",
        help='output PDB filepath (default: replace ".pdb" with "_add_elem.pdb"; when provided, --overwrite is ignored)',
    )
    ap.add_argument(
        "--overwrite",
        type=_parse_bool,
        default=False,
        help="Overwrite the input file in-place when -o/--out is omitted. Use True/False.",
    )
    args = ap.parse_args()

    if not os.path.isfile(args.in_pdb):
        print(f"[ERR] Input not found: {args.in_pdb}", file=sys.stderr)
        sys.exit(1)

    try:
        assign_elements(args.in_pdb, args.out, overwrite=args.overwrite)
    except Exception as e:
        print(f"[ERR] Failed: {e}", file=sys.stderr)
        sys.exit(2)

# -----------------------------
# Click subcommand (pdb2reaction add-elem-info)
# -----------------------------
@click.command(
    help="Add/repair element columns (77–78) in a PDB using Biopython.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "in_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input PDB filepath",
)
@click.option(
    "-o", "--out",
    "out_pdb",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help='Output PDB filepath (default: replace ".pdb" with "_add_elem.pdb"; when provided, --overwrite is ignored)',
)
@click.option(
    "--overwrite",
    type=click.BOOL,
    default=False,
    show_default=True,
    help="Overwrite the input file in-place when -o/--out is omitted.",
)
def cli(in_pdb: Path, out_pdb: Optional[Path], overwrite: bool) -> None:
    """Click wrapper to run via the `pdb2reaction add-elem-info` subcommand."""
    try:
        assign_elements(str(in_pdb), (str(out_pdb) if out_pdb else None), overwrite=overwrite)
    except SystemExit as e:
        raise e
    except Exception as e:
        click.echo(f"[ERR] Failed: {e}", err=True)
        sys.exit(2)
