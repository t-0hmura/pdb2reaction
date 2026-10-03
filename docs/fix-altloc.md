# `fix-altloc` (resolve PDB alternate locations)

## Overview

`fix-altloc` **removes alternate locations (altLoc)** from PDB files. For each residue it keeps one altLoc label, the one with the highest mean occupancy, so each residue is one conformer that was actually deposited. Commands that read a PDB apply the same rule on their own when they read it. Use `fix-altloc` when you need the cleaned file itself.

### What it is for

* **A clean PDB file to keep**: one conformer per residue, for other programs or for your records.
* **Many files at once**: every `.pdb` in a directory, optionally with its subdirectories.
* **Inspecting the choice**: open the cleaned file next to the input in a structure viewer to see which conformer each residue keeps, before you run a calculation.

---

## Examples

### 1. One file

Clean one file and write `1abc_clean.pdb`.

```bash
pdb2reaction fix-altloc -i 1abc.pdb
```

The console prints `[fix-altloc] Fixed altLoc → 1abc_clean.pdb`.

### 2. Choose the output file

```bash
pdb2reaction fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb
```

### 3. A directory, recursively

Clean every `.pdb` under `./structures` and write the results to `./cleaned` with the same subdirectories.

```bash
pdb2reaction fix-altloc -i ./structures -o ./cleaned --recursive
```

### 4. In place, with a backup

```bash
pdb2reaction fix-altloc -i ./structures --inplace --recursive
```

---

## How it works

1. **Detecting altLoc**:
`fix-altloc` looks for non-blank altLoc characters (column 17) in each file.
2. **Grouping by residue**:
Labeled ATOM and HETATM records are grouped by residue: chain ID, residue number, insertion code, and segID. The residue name is not part of the key.
3. **Choosing one label per residue**:
The label whose atoms have the highest mean occupancy (columns 55–60) is chosen; a tie goes to the label that appears first.
4. **Writing**:
Blank (shared) atoms and the atoms of the chosen label are kept, and column 17 is blanked.

ANISOU records are kept only for the atoms that remain (same serial number); every other record is written unchanged.

### Different atom counts between altLoc states

When the altLoc states contain different atoms, only the atoms of the chosen label remain, and an atom found only in the other label is dropped.

```text
Input:
 ATOM 1 N ALYS A 1... 0.50 # altLoc A
 ATOM 2 CA ALYS A 1... 0.50 # altLoc A
 ATOM 3 CB ALYS A 1... 0.50 # altLoc A
 ATOM 4 CG ALYS A 1... 0.50 # altLoc A
 ATOM 5 N BLYS A 1... 0.40 # altLoc B
 ATOM 6 CA BLYS A 1... 0.40 # altLoc B
 ATOM 7 CB BLYS A 1... 0.40 # altLoc B
 ATOM 8 CG BLYS A 1... 0.40 # altLoc B
 ATOM 9 CD BLYS A 1... 0.40 # altLoc B only

Output:
 ATOM 1 N LYS A 1... 0.50 # from A (higher occupancy)
 ATOM 2 CA LYS A 1... 0.50 # from A
 ATOM 3 CB LYS A 1... 0.50 # from A
 ATOM 4 CG LYS A 1... 0.50 # from A
 (CD, in altLoc B only, is dropped)
```

---

## Output files

* **File input**: `<input>_clean.pdb` by default, or the path given with `-o`.
* **Directory input**: `<input>_clean/` by default, or the directory given with `-o`, with the same relative paths as the input. The console prints `[fix-altloc] Processed N file(s) → …` and, for files without altLoc, `Skipped N file(s)`.
* **`--inplace`**: the input files are overwritten, and each original is saved as `<name>.pdb.bak`.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input PDB file or directory |
| `-o, --output` | path | `None` | Output file (file input) or directory (directory input); without it, `<input>_clean.pdb` or `<input>_clean/` |
| `--recursive/--no-recursive` | flag | `False` | For a directory, also process `.pdb` files in subdirectories |
| `--inplace/--no-inplace` | flag | `False` | Overwrite the input files, keeping `.bak` backups |
| `--overwrite/--no-overwrite` | flag | `False` | Allow overwriting existing output files; without it, an existing output stops the run with `Output exists: <path> (use --overwrite to overwrite)` |
| `--force/--no-force` | flag | `False` | Process files even when no altLoc is found |

See the [generated CLI reference](reference/commands/fix_altloc.md) for every option.

---

## Notes

* **Files without altLoc**: a file whose column 17 is blank everywhere is skipped and nothing is written. `--force` processes it anyway.
* **`--inplace` and `-o`** cannot be used together. An existing `.bak` file is not replaced, so it keeps the file from before the first in-place run.
* **Serial numbers** are not renumbered, so gaps can remain where atoms were removed. `CONECT` and other connectivity or annotation records are not updated.
* **Kept records** keep their coordinates, occupancies, B-factors, charges, insertion codes, and order.
* **MODEL/ENDMDL blocks** are processed one by one.
* **The occupancy rule is a heuristic**: when the active-site conformer must be chosen by chemical contacts or by how the deposited ensemble is interpreted, choose it yourself in a structure editor and check it.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [extract](extract.md) — active-site model extraction, which applies the same altLoc rule when it reads a PDB
* [add-elem-info](add-elem-info.md) — repair the element columns of a PDB
* [Before you run: the input structures](getting-started.md#before-you-run-the-input-structures) — add hydrogens before a calculation
* [all](all.md) — the full workflow
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
