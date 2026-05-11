# `pdb2reaction add-elem-info`

## Purpose

Repair the element column (PDB cols 77–78) when it is blank or
inconsistent with the atom name (a common state for tleap-emitted PDBs;
without it, `extract`'s element-aware truncation logic fails).

`pdb2reaction all` preflight-runs this only when the element field is
missing. Call it explicitly when invoking `extract` / `opt` / `tsopt`
directly on a raw RCSB PDB.

## Synopsis

```bash
pdb2reaction add-elem-info -i in.pdb -o out.pdb
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB |
| `-o, --out` | path | `<input>_add_elem.pdb` (auto) | Output PDB with element column populated |
| `--overwrite / --no-overwrite` | flag | `--no-overwrite` | Overwrite input file in-place (only honored when `-o/--out` omitted) |

## Examples

```bash
pdb2reaction add-elem-info -i raw.pdb -o cleaned.pdb
```

## Algorithm

The element is inferred from atom name + residue name with the
following priority (`add_elem_info.guess_element`):

1. **Ion residues** (residue name is in the internal `ION` dict):
   the residue name is the element source. Polyatomic ions
   (`NH4`, `H3O+`, …) dispatch per atom-name prefix (`H`/`D` → H,
   `N` → N, `O` → O); monatomic metals/halogens use the residue.
2. **Polymers and water** (protein, nucleic acid, water): use the
   PDB convention element subset (`H`/`C`/`N`/`O`/`S`/`P`/`Se`).
3. **Other ligands**: 1–2 alphabetic prefix of the atom name with a
   case-sensitivity table for two-letter symbols (`Mg`, `Mn`, `Fe`,
   `Zn`, `Ca`, …); 4-character names starting with a digit (e.g.
   `1HG2`) use the second character.
4. Unresolved → reported in the diagnostic summary (truncated at 50
   entries); the existing element field is left unchanged.

## Caveats

- Element-column values in the output are recomputed from atom / residue
  names (existing values are replaced). The input file is not modified
  unless `--overwrite` is passed with no `-o/--out`. Use a diff to confirm
  the changes are sensible.
- Atom names that don't follow the standard convention (e.g.
  exotic ligand names) may be misclassified; verify by spot-check.

## See also

- `fix-altloc.md` — typically run together as a "PDB cleanup" step.
- `extract.md` — depends on a well-formed element column.
- [`pdb2reaction-structure-io/pdb.md`](../pdb2reaction-structure-io/pdb.md) — PDB column layout reference.
