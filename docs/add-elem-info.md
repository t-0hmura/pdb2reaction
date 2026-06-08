# `add-elem-info`

Repair the element-symbol columns (77–78) of ATOM/HETATM records in a PDB file. Per-atom elements are inferred from atom name + residue context, using `Bio.PDB.PDBParser` re-parse + per-atom element inference; only columns 77–78 are rewritten.

## When to use

- Use when a PDB file has missing or wrong element columns (77–78) and downstream subcommands (`extract`, `opt`, `tsopt`, ...) reject it.
- `all` auto-invokes `add-elem-info` as a preflight, so manual use is only needed before standalone subcommands.

## Quick examples

```bash
# Populate element fields and write to "<input>_add_elem.pdb"
pdb2reaction add-elem-info -i 1abc.pdb

# Write to a specific output file
pdb2reaction add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb

# Overwrite the input file in-place
pdb2reaction add-elem-info -i 1abc.pdb --overwrite
```

## Inputs

Command form:

```bash
pdb2reaction add-elem-info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite/--no-overwrite]
```

| Input | Required | Notes |
| --- | --- | --- |
| `-i, --input PATH` | yes | Input PDB file. |
| `-o, --out PATH` | no | Output path. When set, `--overwrite` is ignored; defaults to `<input>_add_elem.pdb`. |

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
3. Write the structure through `PDBIO` to the chosen output path (see [Outputs](#outputs) for the default / `-o` / `--overwrite` precedence).
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
- Only columns 77–78 are modified; coordinates, occupancies, B-factors, charges, altlocs,
  insertion codes, and record ordering stay untouched.
- ATOM and HETATM records across all models/chains/residues are supported.
- Deuterium labels map to hydrogen; selenium (`SE*`) and halogens are recognized automatically.
- Re-running on a PDB that already carries valid element symbols is a no-op (atoms pass through unchanged). See [all](all.md) for how the `all` preflight invokes `add-elem-info` automatically only when element columns are missing.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [extract](extract.md) -- Active site model extraction after element-column repair
- [all](all.md) -- End-to-end workflow entrypoint
