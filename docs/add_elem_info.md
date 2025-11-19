# `add-elem-info` subcommand

## Overview
`add-elem-info` populates the element-symbol columns (77–78) of ATOM/HETATM
records using the same heuristics as the extractor module. The command relies on
Biopython to parse/write the structure, preserves coordinates and metadata, and
only rewrites element fields unless `--overwrite` is supplied.

## Usage
```bash
pdb2reaction add-elem-info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite]
```

### Examples
```bash
# Repair an input PDB in-place
pdb2reaction add-elem-info -i 1abc.pdb

# Write to a new file and force re-assignment for all atoms
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb --overwrite
```

## Workflow
1. Parse the input file with `Bio.PDB.PDBParser`, mirroring the residue
   definitions used elsewhere (`extract.AMINO_ACIDS`, `WATER_RES`, `ION`).
2. For each atom, guess the element by combining the atom name, residue name,
   and whether the record is HETATM:
   - Ion residues: prefer the residue name for monatomic ions; polyatomic ions
     (e.g., NH4, H3O+) are split into per-atom assignments.
   - Proteins/nucleic acids/water: apply special handling for H/D, Se, and
     first-letter mapping for C/N/O/P/S; carbon sidechain labels default to C.
   - Other ligands: use atom-name prefixes and fall back to element-symbol
     normalization (recognising halogens, deuterium → hydrogen, etc.).
3. Preserve existing element fields unless `--overwrite` is set; in overwrite
   mode, all atoms are re-inferred even if the original fields were populated.
4. Write the structure through `PDBIO`, touching only columns 77–78. When no
   output path is given, the input file is overwritten.
5. Print a summary reporting how many atoms were newly assigned, kept, or
   overwritten plus per-element totals and a truncated list of unresolved atoms.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output path; overwrites the input file when omitted. | _None_ |
| `--overwrite` | Re-infer elements even if the original fields are present. | `False` |

## Outputs
- A PDB with element symbols populated (`--out` destination or the input file in
  place).
- Console report with totals for processed/assigned/kept/overwritten atoms,
  per-element counts, and up to 50 unresolved atoms.

## Notes
- Only columns 77–78 are modified; coordinates, occupancies, altlocs, insertion
  codes, and record ordering stay untouched.
- ATOM and HETATM records across all models/chains/residues are supported.
- Deuterium labels map to hydrogen; selenium (`SE*`) and halogens are recognized
  automatically.
- Requires Biopython and Click at runtime.
