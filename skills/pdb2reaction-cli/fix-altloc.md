# `pdb2reaction fix-altloc`

## Purpose

Blank the PDB altLoc column (col 17) without shifting any other field and
select one coherent non-blank altLoc label per residue. The default rule is
**highest mean occupancy across that residue's labelled atoms, then earliest
label appearance** — there is no hard-coded preference for A. Blank/shared
atoms are retained; atoms unique to an unselected conformer are dropped rather
than merged into a hybrid residue. Geometry workflows apply the same
residue-coherent selection through the common input bridge. Call this utility
when a cleaned PDB itself is needed, for directory processing, or to inspect
the selected conformer before calculation.

## Synopsis

```bash
pdb2reaction fix-altloc -i in.pdb [-o out.pdb] [--inplace] [--overwrite] [--recursive]
```

`-i` can be a single PDB or a directory; `-o` matches accordingly
(file → file, directory → directory).

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB file or directory |
| `-o, --out` | path | derived | Output file (file input) or directory (directory input). Without `--inplace`, omission creates `<stem>_clean.pdb` or `<directory>_clean/`. |
| `--inplace / --no-inplace` | flag | `--no-inplace` | Overwrite input file(s) in place; writes a `.bak` next to each file |
| `--overwrite / --no-overwrite` | flag | `--no-overwrite` | Overwrite existing files in non-inplace output mode; `--inplace` already replaces each input after creating its backup |
| `--recursive / --no-recursive` | flag | `--no-recursive` | Recurse into sub-directories when `-i` is a directory |
| `--force / --no-force` | flag | `--no-force` | Process files even if no altLoc detected (default skips them) |

## Examples

```bash
# Single PDB
pdb2reaction fix-altloc -i raw.pdb -o cleaned.pdb

# Whole directory
pdb2reaction fix-altloc -i raw_pdbs/ -o cleaned_pdbs/

# In place (overwrite)
pdb2reaction fix-altloc -i raw.pdb --inplace
```

## Output

- Single-file mode: one PDB at `-o` (or in place).
- Directory mode: one PDB per input file, mirrored under `-o/`.
- altLoc column (col 17) is blanked for all surviving atoms.

## Caveats

- Selection is residue-level, but it is still an occupancy heuristic. If the
  active-site conformer must be chosen by ligand contacts or mechanism rather
  than occupancy, select it with a structure editor and validate the result.
- Use `add-elem-info` only when element fields are missing/wrong and
  `fix-altloc` only when alternate locations are present; neither cleanup is a
  mandatory rewrite of every RCSB file.

## See also

- `add-elem-info.md` — usually run together.
- `extract.md` — automatically resolves one coherent altloc per residue; use
  this utility only when a separate cleaned PDB deliverable is needed.
- `../pdb2reaction-structure-io/pdb.md` — col 17 (altLoc) reference.
