# `pdb2reaction fix-altloc`

## Purpose

Blank the PDB altLoc column (col 17) without shifting any other field,
and keep **one** altLoc per atom. The default rule is **highest
occupancy first, then earliest appearance** — there is no letter-based
selection. `pdb2reaction all` preflight-runs this only when altLoc is
detected. Call it explicitly when invoking `extract` / `opt` / `tsopt`
directly on a raw RCSB PDB.

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
| `-o, --out` | path | derived | Output file (if input is file) or directory (if directory). Omit to overwrite in place via `--inplace`. |
| `--inplace` | flag | False | Overwrite the input file in place |
| `--overwrite` | flag | False | Overwrite an existing output file |
| `--recursive` | flag | False | Recurse into sub-directories when `-i` is a directory |

The selection rule is **fixed**: highest occupancy, then earliest
appearance. There is no `--keep <letter>` flag.

## Examples

```bash
# Single PDB
pdb2reaction fix-altloc -i raw.pdb -o cleaned.pdb

# Whole directory
pdb2reaction fix-altloc -i raw_pdbs/ -o cleaned_pdbs/

# In place (overwrite)
pdb2reaction fix-altloc -i raw.pdb --inplace --overwrite
```

## Output

- Single-file mode: one PDB at `-o` (or in place).
- Directory mode: one PDB per input file, mirrored under `-o/`.
- altLoc column (col 17) is blanked for all surviving atoms.

## Caveats

- The selection rule (highest occupancy, then earliest) is the only
  rule. If the active site has a high-occupancy alt-loc that is **not**
  the chemistry you want, edit the PDB by hand or run a different tool.
- `add-elem-info` and `fix-altloc` are typically run together on raw
  RCSB downloads; order does not matter.

## See also

- `add-elem-info.md` — usually run together.
- `extract.md` — expects altloc-resolved PDB input.
- `../pdb2reaction-structure-io/pdb.md` — col 17 (altLoc) reference.