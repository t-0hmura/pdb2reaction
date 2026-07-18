#!/usr/bin/env python3
"""Verify exact chain/residue selection across the internal PDB bridge."""

from __future__ import annotations

import sys
from pathlib import Path


pdb_path = Path(sys.argv[1])
cif_path = Path(sys.argv[2])
pdb_atoms = [
    line for line in pdb_path.read_text(encoding="utf-8").splitlines()
    if line.startswith(("ATOM  ", "HETATM"))
]
if len(pdb_atoms) != 2:
    raise SystemExit(f"expected exactly two selected atoms, found {len(pdb_atoms)}")
if {line[12:16].strip() for line in pdb_atoms} != {"C1", "O1"}:
    raise SystemExit("selected PDB does not contain exactly PRE:C1/O1")
if any(line[17:20].strip() != "PRE" for line in pdb_atoms):
    raise SystemExit("selected PDB contains a residue other than PRE")

cif_lines = cif_path.read_text(encoding="utf-8").splitlines()
cif_rows = [line.split() for line in cif_lines if line.startswith(("ATOM ", "HETATM "))]
if len(cif_rows) != 2:
    raise SystemExit(f"expected exactly two selected CIF rows, found {len(cif_rows)}")

# Resolve each field from the file's own `_atom_site.` header rather than
# hard-coding a column index: mmCIF loops declare their own column order, so a
# fixed index silently checks the wrong field the moment the writer adds or
# reorders one.
header = [line.strip() for line in cif_lines if line.strip().startswith("_atom_site.")]
def column(name: str) -> int:
    try:
        return next(i for i, tag in enumerate(header) if tag == f"_atom_site.{name}")
    except StopIteration:
        raise SystemExit(f"CIF declares no _atom_site.{name}")

expected = {"auth_seq_id": "10001", "auth_comp_id": "PRE", "auth_asym_id": "LONG_CHAIN"}
for row in cif_rows:
    for name, want in expected.items():
        got = row[column(name)]
        if got != want:
            raise SystemExit(
                f"CIF auth identity not preserved: _atom_site.{name} is {got!r}, expected {want!r}"
            )
