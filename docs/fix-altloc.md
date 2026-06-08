# `fix-altloc`

Remove alternate-location (altLoc) indicators from PDB files by selecting the best conformer for each atom based on occupancy and dropping duplicates. For each atom, the highest-occupancy conformer is kept (ties broken by file order) and column 17 is blanked. Run it before any cluster extraction or geometry stage (`extract`, `opt`, `tsopt`, ...) that cannot consume multi-conformer input.

## Examples

```bash
# Command form
pdb2reaction fix-altloc -i INPUT.pdb [-o OUTPUT.pdb] [OPTIONS]

# Process a single file (output: INPUT_clean.pdb)
pdb2reaction fix-altloc -i 1abc.pdb

# Specify output file
pdb2reaction fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb

# Process a directory recursively
pdb2reaction fix-altloc -i ./structures -o ./cleaned --recursive

# Overwrite input files in-place (creates .bak backups)
pdb2reaction fix-altloc -i ./structures --inplace --recursive
```

## Workflow

1. Check if the input file contains any non-blank altLoc characters (column 17).
 - If no altLoc is found and `--force` is not set, skip the file.
2. For each ATOM/HETATM record, build an identity key ignoring the altLoc field:
 - record name, atom name, residue name, chain ID, residue sequence, insertion code, segID
3. Among atoms with the same identity key, select the one with:
 - Highest occupancy (columns 55–60)
 - If tied, the earliest appearance in the file
4. Write output with:
 - Only the selected atoms retained
 - altLoc column (17) blanked to a single space
 - ANISOU records filtered to match retained atoms

### Handled records

- `ATOM` / `HETATM`: altLoc selection and blanking
- `ANISOU`: kept only if the corresponding ATOM/HETATM line (same serial) is kept

### Handling different atom counts between altLoc states

When different altLoc states contain different atoms (e.g., altLoc A has atoms
N, CA, CB, CG while altLoc B has N, CA, CB, CD), `fix-altloc` handles this correctly:

- **Duplicate atoms** (same residue + atom name in multiple altLocs, e.g., N, CA, CB):
  The best one is selected based on occupancy (highest first, then earliest in file).
- **Unique atoms** (only present in one altLoc, e.g., CG in A, CD in B):
  ALL unique atoms are preserved in the output.

This ensures the output structure contains all atoms from all altLoc states,
with only true duplicates resolved to a single conformer.

**Example:**
```
Input:
 ATOM 1 N AALA A 1... 0.50 # altLoc A
 ATOM 2 CA AALA A 1... 0.50 # altLoc A
 ATOM 3 CG AALA A 1... 0.50 # altLoc A only
 ATOM 4 N BALA A 1... 0.40 # altLoc B
 ATOM 5 CA BALA A 1... 0.40 # altLoc B
 ATOM 6 CD BALA A 1... 0.40 # altLoc B only

Output:
 ATOM 1 N ALA A 1... 0.50 # from A (higher occ)
 ATOM 2 CA ALA A 1... 0.50 # from A (higher occ)
 ATOM 3 CG ALA A 1... 0.50 # kept (A only)
 ATOM 6 CD ALA A 1... 0.40 # kept (B only)
```

## Outputs

- A PDB file with alternate locations removed:
 - File input: `<input>_clean.pdb` by default (when `-o/--out` is omitted)
 - Directory input: `<input>_clean/` directory by default (mirrors subpaths)
 - `OUTPUT.pdb` if `-o/--out` is provided
 - Original file overwritten if `--inplace` is set (backup saved as `<input>.pdb.bak`)

## Python API

For programmatic use, the module exports:
```python
from pdb2reaction.io.pdb_fix import has_altloc, fix_altloc_file

# Check if a file has altLoc
if has_altloc(Path("input.pdb")):
 # Fix altLoc
 was_processed = fix_altloc_file("input.pdb", "output.pdb", overwrite=True)
```

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file or directory. | Required |
| `-o, --out PATH` | Output file (if input is a file) or directory (if input is a directory). | File input: `<input>_clean.pdb`; directory input: `<input>_clean/` |
| `--recursive/--no-recursive` | Process `*.pdb` files recursively when input is a directory. | `False` |
| `--inplace/--no-inplace` | Overwrite input file(s) in-place (creates `.bak` backup). | `False` |
| `--overwrite/--no-overwrite` | Allow overwriting existing output files. | `False` |
| `--force/--no-force` | Process files even if no altLoc is detected. | `False` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes

- By default, if a file contains **no altLoc characters** (all column 17 positions are blank), the file is **skipped** and no output is written. Use `--force` to process files regardless of altLoc presence.
- Atom serial numbers are **NOT renumbered** (gaps may remain after duplicate removal).
- `CONECT` and other connectivity/annotation records are **NOT updated**.
- Only column 17 (altLoc) is modified; coordinates, occupancies, B-factors, charges,
  insertion codes, and record ordering stay untouched (except for duplicate removal).
- MODEL/ENDMDL blocks are processed independently.

## See Also

- [Common Error Recipes](recipes-common-errors.md)
- [Troubleshooting](troubleshooting.md)
- [all](all.md) -- end-to-end workflow that auto-invokes `add-elem-info` then `fix-altloc` as preflight
