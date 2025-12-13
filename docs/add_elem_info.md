# `add-elem-info` subcommand

## Overview
`add-elem-info` populates (or repairs) the element-symbol columns (77–78) of ATOM/HETATM
records in a PDB file. The command relies on Biopython to parse/write the structure,
preserves coordinates and metadata, and only rewrites element fields unless `--overwrite`
is supplied.

### Output behavior (important)
- If `-o/--out` is **omitted** and `--overwrite` is **not** set, the output is written to
  `<input>_add_elem.pdb` (i.e., it replaces a trailing `.pdb` with `_add_elem.pdb`;
  if the input does not end with `.pdb`, `_add_elem.pdb` is appended).
- If `--overwrite` **is set** **and** `-o/--out` is **omitted**, the **input file is overwritten
  in-place**. When `-o/--out` is supplied, `--overwrite` is ignored.

## Usage
```bash
pdb2reaction add-elem-info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite]
```

## Examples
```bash
# Populate element fields and write to "<input>_add_elem.pdb"
pdb2reaction add-elem-info -i 1abc.pdb

# Write to a specific output file (preserve existing symbols by default)
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

# Overwrite the input file in-place and re-infer all element symbols
pdb2reaction add-elem-info -i 1abc.pdb --overwrite
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
4. Write the structure through `PDBIO`:
   - default output: `<input>_add_elem.pdb` (when `-o/--out` is omitted and `--overwrite` is not set)
   - `-o/--out`: write to the specified path; `--overwrite` is ignored when this is provided
   - `--overwrite` (without `-o/--out`): overwrite the input path in-place
5. Print a summary reporting how many atoms were newly assigned, kept, or
   overwritten plus per-element totals and a truncated list of unresolved atoms.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output path. When set, `--overwrite` is ignored. | _None_ → `<input>_add_elem.pdb` |
| `--overwrite` | Re-infer elements even if the original fields are present **and** overwrite the input file in-place when `-o/--out` is omitted. | `False` |

## Outputs
- A PDB file with element symbols populated/corrected:
  - `<input>_add_elem.pdb` by default (when `-o/--out` is omitted and `--overwrite` is not set)
  - `OUTPUT.pdb` if `-o/--out` is provided (regardless of `--overwrite`)
  - `INPUT.pdb` overwritten in-place if `--overwrite` is set without `-o/--out`
- Console report with totals for processed/assigned/kept/overwritten atoms,
  per-element counts, and up to 50 unresolved atoms.

## Notes
- Existing element fields are detected by scanning the original file’s ATOM/HETATM lines
  (serials 7–11, elements 77–78) to reflect the true presence/absence and avoid parser side effects.
- Only columns 77–78 are modified; coordinates, occupancies, B-factors, charges, altlocs,
  insertion codes, and record ordering stay untouched.
- ATOM and HETATM records across all models/chains/residues are supported.
- Deuterium labels map to hydrogen; selenium (`SE*`) and halogens are recognized automatically.
- Requires Biopython (`Bio.PDB`) and Click at runtime.
