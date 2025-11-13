# `add_elem_info` subcommand

## Purpose
Repairs or populates element-symbol columns (77–78) in PDB files using Biopython heuristics aligned with the extractor’s residue definitions.

## Usage
```bash
pdb2reaction add_elem_info -i INPUT.pdb [-o OUTPUT.pdb] [--overwrite]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output PDB file. If omitted, the input is overwritten. | _None_ |
| `--overwrite` | Re-infer all element fields even when present. | `False` |

## YAML configuration (`--args-yaml`)
Not supported.

## Outputs
- A PDB file with element symbols filled (`--out` path or the original file when omitted).
- Console summary including counts of assigned/kept/overwritten elements and up to 50 unresolved atoms.

## Notes
- Recognises standard protein, nucleic acid, water, ion, and common ligand nomenclature; deuterium (`D`) is mapped to hydrogen.
- Depends on Biopython’s `PDBParser`/`PDBIO`; only ATOM/HETATM records are modified.
- Keeps coordinates and other columns unchanged; only columns 77–78 are rewritten when assignments exist.
