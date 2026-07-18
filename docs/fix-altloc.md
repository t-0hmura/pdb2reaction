# `fix-altloc`

Remove alternate-location (altLoc) indicators from PDB files by selecting one
coherent non-blank altLoc label per residue. The label with the highest mean
occupancy across that residue's labeled atoms is selected; ties are broken by
the label's first appearance. Blank/shared atoms are retained, coordinates from
other labels are dropped, and column 17 is blanked on surviving records.
Geometry workflows apply the same residue-coherent rule automatically through
the common input bridge. Use this command when a cleaned PDB file itself is a
deliverable, for directory processing, or to inspect the choice.

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
2. Group labeled ATOM/HETATM records by residue (residue name, chain ID,
   residue sequence, insertion code, and segID).
3. Select one non-blank label per residue using:
 - Highest mean parsed occupancy across that label's atoms (columns 55–60)
 - A label with no parsed occupancies ranks below every label with a parsed mean
 - Equal scores, including the all-missing case, are resolved by first appearance
4. Keep blank/shared atoms plus atoms from the selected label. Resolve any
   remaining duplicate atom identities by occupancy and file order.
5. Write output with:
 - Only blank/shared atoms and the selected residue conformer retained
 - altLoc column (17) blanked to a single space
 - ANISOU records filtered to match retained atoms

### Handled records

- `ATOM` / `HETATM`: altLoc selection and blanking
- `ANISOU`: kept only if the corresponding ATOM/HETATM line (same serial) is kept

### Handling different atom counts between altLoc states

When altLoc states contain different atoms (for example A has N, CA, CB, CG
while B has N, CA, CB, CD), only atoms belonging to the selected residue label
are retained. An atom unique to an unselected label is dropped. This avoids the
old per-atom behavior that produced an A/B hybrid corresponding to no deposited
conformer.

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
- Surviving coordinate records retain their coordinates, occupancies, B-factors,
  charges, insertion codes, and relative order; records from unselected altLoc
  labels are removed and column 17 is blanked.
- MODEL/ENDMDL blocks are processed independently.
- The residue-level occupancy rule is still a heuristic. If the active-site
  conformer must be selected by chemical contacts or a deposited ensemble
  interpretation, select it explicitly in a structure editor and inspect it.

## See Also

- [Common Error Recipes](recipes-common-errors.md)
- [Troubleshooting](troubleshooting.md)
- [all](all.md) -- end-to-end workflow using the common altloc/large-structure bridge
