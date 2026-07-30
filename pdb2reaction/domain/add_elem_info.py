# pdb2reaction/domain/add_elem_info.py

"""
Repair PDB element symbols (columns 77-78) without rewriting other records.

Example:
    pdb2reaction add-elem-info -i 1abc.pdb

For detailed documentation, see: docs/add-elem-info.md
"""

from __future__ import annotations

import collections
import re
import sys
import time
from pathlib import Path
from typing import Optional

import click

# Reuse the canonical residue/ion/water tables (leaf data module) to keep
# element inference and charge inference reading one source of truth.
from pdb2reaction.domain.residue_data import AMINO_ACIDS, ION, WATER_RES

# Element symbols (IUPAC, 1–118)
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

# Helper: normalize strings to element symbols
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


def _symbol_from_aligned_atom_name(atom_name: str) -> Optional[str]:
    """Infer an element from the four-column PDB atom-name alignment."""
    if len(atom_name) < 4:
        return None
    raw = atom_name[:4]
    if raw[0].isspace():
        return _normalize_symbol(raw.lstrip()[:1])
    if raw[0].isdigit():
        return _normalize_symbol(raw.lstrip("0123456789")[:1])
    return _normalize_symbol(raw[:2])


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

# Element inference (use residue to disambiguate)
def guess_element(atom_name: str, resname: str, _is_het: bool = False) -> Optional[str]:
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

    if res_u in {k.upper() for k in ION.keys()}:
        # Genuinely polyatomic ions (NH4, H3O+) contain more than one element,
        # so decide per atom name (treat D* as H). Monatomic metal/halogen ions
        # fall through to residue-name resolution below, so an ion whose symbol
        # starts with H/N/O (Na, Ni, Hg, Nd, Hf, He, …) is not mislabelled.
        if res_u in {"NH4", "H3O+", "H3O"}:
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

    aligned = _symbol_from_aligned_atom_name(atom_name)
    if aligned is not None:
        return aligned

    # Unaligned programmatic inputs retain the historical prefix fallback.
    if name_u.startswith(("H", "D")):
        return "H"
    if name_u.startswith("C") and not name_u.startswith("CL"):
        return "C"
    if name_u.startswith("N"):
        return "N"
    if name_u.startswith("O"):
        return "O"
    if name_u.startswith("P"):
        return "P"

    sym = _normalize_symbol(name_u)
    if sym:
        return sym

    return None


def _replace_element_field(line: str, symbol: str) -> str:
    if line.endswith("\r\n"):
        content, ending = line[:-2], "\r\n"
    elif line.endswith(("\n", "\r")):
        content, ending = line[:-1], line[-1:]
    else:
        content, ending = line, ""
    content = content.ljust(78)
    return content[:76] + f"{symbol:>2}" + content[78:] + ending


def assign_elements(in_pdb: str, out_pdb: Optional[str], overwrite: bool = False) -> None:
    # If an explicit output path is provided, never overwrite in-place even when --overwrite is
    # passed. This keeps -o/-\-out as the higher-priority choice.
    effective_overwrite = overwrite and out_pdb is None

    total = 0
    assigned_or_updated = 0
    unknown = []
    by_element = collections.Counter()
    with open(in_pdb, "r", encoding="utf-8", errors="surrogateescape", newline="") as handle:
        lines = handle.readlines()

    model_id: object = 0
    rewritten = []
    for line in lines:
        if line.startswith("MODEL"):
            token = line[10:14].strip()
            model_id = int(token) if token.isdigit() else token or model_id
        if not line.startswith(("ATOM  ", "HETATM")):
            rewritten.append(line)
            continue

        total += 1
        atom_name = line[12:16]
        resname = line[17:20]
        symbol = guess_element(atom_name, resname, line.startswith("HETATM"))
        serial_text = line[6:11].strip()
        serial = int(serial_text) if serial_text.isdigit() else None
        if symbol is None:
            unknown.append(
                (
                    model_id,
                    line[21:22].strip(),
                    resname.strip(),
                    line[22:26].strip(),
                    line[26:27].strip(),
                    atom_name.strip(),
                    serial,
                )
            )
            rewritten.append(line)
            continue

        previous = line[76:78].strip() if len(line.rstrip("\r\n")) >= 78 else ""
        by_element[symbol] += 1
        if previous != symbol:
            assigned_or_updated += 1
        rewritten.append(_replace_element_field(line, symbol))

    if effective_overwrite:
        out_path = in_pdb
    else:
        out_path = out_pdb if out_pdb else _default_out_pdb_path(in_pdb)
    with open(out_path, "w", encoding="utf-8", errors="surrogateescape", newline="") as handle:
        handle.writelines(rewritten)

    # Summary
    click.echo(f"[OK] Wrote: {out_path}")
    click.echo(f"  total atoms                 : {total}")
    click.echo(f"  assigned/updated            : {assigned_or_updated}")
    if by_element:
        top = ", ".join(f"{k}:{v}" for k, v in by_element.most_common())
        click.echo(f"  assignment breakdown        : {top}")
    if unknown:
        click.echo(f"[WARN] Could not confidently assign {len(unknown)} atoms; left unchanged.")
        for mid, chid, resn, resseq, icode, aname, serial in unknown[:50]:
            s_str = f" serial {serial}" if serial is not None else ""
            click.echo(f"    model {mid} chain {chid} {resn} {resseq}{icode} : {aname}{s_str}")
    if len(unknown) > 50:
        click.echo("    ... (truncated) ...")


# Click subcommand (pdb2reaction add-elem-info)
@click.command(
    help="Add/repair element columns (77–78) in a PDB.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "in_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input PDB filepath.",
)
@click.option(
    "-o", "--out",
    "out_pdb",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help='Output PDB filepath (default: replace ".pdb" with "_add_elem.pdb"; when provided, --overwrite is ignored).',
)
@click.option(
    "--overwrite/--no-overwrite",
    default=False,
    show_default=True,
    help="Overwrite the input file in-place when -o/--out is omitted.",
)
def cli(in_pdb: Path, out_pdb: Optional[Path], overwrite: bool) -> None:
    """Click wrapper to run via the `pdb2reaction add-elem-info` subcommand."""
    time_start = time.perf_counter()
    try:
        assign_elements(str(in_pdb), (str(out_pdb) if out_pdb else None), overwrite=overwrite)
    except SystemExit as e:
        raise e
    except Exception as e:
        click.echo(f"[ERR] Failed: {e}", err=True)
        sys.exit(2)
    from pdb2reaction.core.output import emit
    from pdb2reaction.core.utils import format_elapsed

    emit(
        format_elapsed("[time] Elapsed Time for Add Element Info", time_start),
        narrative=True,
    )
