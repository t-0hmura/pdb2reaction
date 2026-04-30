# `add-elem-info`

## Overview
`add-elem-info` repairs the element-symbol columns (77–78) of ATOM/HETATM
records in a PDB file.

## Usage
```bash
pdb2reaction add-elem-info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite/--no-overwrite]
```

## Examples
```bash
# Populate element fields and write to "<input>_add_elem.pdb"
pdb2reaction add-elem-info -i 1abc.pdb

# Write to a specific output file
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

# Overwrite the input file in-place
pdb2reaction add-elem-info -i 1abc.pdb --overwrite
```

## Workflow
1. Parse the input file with `Bio.PDB.PDBParser`, mirroring the residue
    definitions used in `extract.py` (`AMINO_ACIDS`, `WATER_RES`, `ION`).
2. For each atom, guess the element by combining the atom name, residue name,
    and whether the record is HETATM:
 - Monatomic ion residues in the `ION` dict: use the corresponding element.
 - Proteins/nucleic acids/water: apply special handling for H/D, Se, and
  first-letter mapping for C/N/O/P/S; carbon side-chain labels default to C.
 - Other ligands: use atom-name prefixes and fall back to element-symbol
  normalization (recognizing halogens, deuterium → hydrogen, etc.).
3. Write the structure through `PDBIO`:
 - default output: `<input>_add_elem.pdb` (when `-o/--out` is omitted and `--overwrite` is not `True`)
 - `-o/--out`: write to the specified path; `--overwrite` is ignored when this is provided
 - `--overwrite` (without `-o/--out`): overwrite the input path in-place
4. Print a summary reporting how many atoms were assigned/reassigned, plus
    per-element totals and a truncated list of unresolved atoms.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output path. When set, `--overwrite` is ignored. | _None_ → `<input>_add_elem.pdb` |
| `--overwrite/--no-overwrite` | Overwrite the input file in-place when `-o/--out` is omitted. | `False` |

## Outputs
- A PDB file with element symbols populated/corrected:
 - `<input>_add_elem.pdb` by default (when `-o/--out` is omitted and `--overwrite` is not `True`)
 - `OUTPUT.pdb` if `-o/--out` is provided (regardless of `--overwrite`)
 - `INPUT.pdb` overwritten in-place if `--overwrite` is set without `-o/--out`
- Console report with totals for processed/assigned atoms,
  per-element counts, and up to 50 unresolved atoms.

## Integration with `all` workflow

When running `pdb2reaction all`, `add-elem-info` is **automatically invoked as a
preflight step** on every PDB input before any other processing (followed by
`fix-altloc`, then active site model extraction). This ensures element symbols
are populated before bond detection or charge assignment runs. You therefore do
**not** need to run `add-elem-info` manually before `pdb2reaction all`; it is
only required when you use individual subcommands such as `extract` or `opt`
directly on a PDB that lacks element columns.

## Notes
- Only columns 77–78 are modified; coordinates, occupancies, B-factors, charges, altlocs,
  insertion codes, and record ordering stay untouched.
- ATOM and HETATM records across all models/chains/residues are supported.
- Deuterium labels map to hydrogen; selenium (`SE*`) and halogens are recognized automatically.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [extract](extract.md) -- Active site model extraction after element-column repair
- [all](all.md) -- End-to-end workflow entrypoint
