# `add-elem-info`

Repair the element-symbol columns (77–78) of ATOM/HETATM records in a PDB file. Each element is inferred from the fixed-column atom name and residue context. The `all` preflight repairs blank element fields only; it preserves nonblank values without chemically validating them. Run this utility explicitly for wrong nonblank symbols or before a standalone command that receives missing fields.

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
1. Read the raw PDB records and classify atoms with the residue definitions
    used in `extract.py` (`AMINO_ACIDS`, `WATER_RES`, `ION`).
2. For each atom, guess the element by combining the atom name, residue name,
    and whether the record is HETATM:
 - Monatomic ion residues in the `ION` dict: use the corresponding element.
 - Proteins/nucleic acids/water: apply special handling for H/D, Se, and
  first-letter mapping for C/N/O/P/S; carbon side-chain labels default to C.
 - Other ligands: use atom-name prefixes and fall back to element-symbol
  normalization (recognizing halogens, deuterium → hydrogen, etc.).
3. Replace only columns 77–78 on ATOM/HETATM records and write all other
   columns and records unchanged (see [Outputs](#outputs) for path precedence).
4. Print a summary reporting how many atoms were assigned/reassigned, plus
    per-element totals and a truncated list of unresolved atoms.

## Outputs
- A PDB file with element symbols populated/corrected:
 - `<input>_add_elem.pdb` by default (when `-o/--out` is omitted and `--overwrite` is not `True`)
 - `OUTPUT.pdb` if `-o/--out` is provided (regardless of `--overwrite`)
 - `INPUT.pdb` overwritten in-place if `--overwrite` is set without `-o/--out`
- Console report with totals for processed/assigned atoms,
  per-element counts, and up to 50 unresolved atoms.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output path. When set, `--overwrite` is ignored. | _None_ → `<input>_add_elem.pdb` |
| `--overwrite/--no-overwrite` | Overwrite the input file in-place when `-o/--out` is omitted. | `False` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes
- Every input line is preserved byte-for-byte except columns 77–78 of
  ATOM/HETATM records selected for repair. HEADER, REMARK, CONECT, ANISOU,
  and the legacy charge column (79–80) are retained.
- ATOM and HETATM records across all models/chains/residues are supported.
- Deuterium labels map to hydrogen; selenium (`SE*`) and halogens are recognized automatically.
- Existing nonblank element symbols are preserved by inference; the writer may still normalize record formatting. See [all](all.md) for the blank-field-only preflight boundary.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [extract](extract.md) -- Active site model extraction after element-column repair
- [all](all.md) -- End-to-end workflow entrypoint
