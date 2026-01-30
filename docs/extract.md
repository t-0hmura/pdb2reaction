# `extract` subcommand

## Overview
Automatically extract binding pockets (active sites) from protein–substrate complexes. The tool applies chemically-aware residue selection (distance cutoffs plus heuristics for disulfides, PRO adjacency, etc.), truncates side chains/backbone segments, optionally appends link hydrogens, and can process single structures or ensembles.

## Usage
```bash
pdb2reaction extract -i COMPLEX.pdb [COMPLEX2.pdb ...]
                     -c SUBSTRATE_SPEC
                     [-o POCKET.pdb [POCKET2.pdb ...]]
                     [--radius Å] [--radius-het2het Å]
                     [--include-H2O {True\|False}]
                     [--exclude-backbone {True\|False}]
                     [--add-linkH {True\|False}]
                     [--selected-resn LIST]
                     [--ligand-charge MAP_OR_NUMBER]
                     [--verbose {True\|False}]
```

### Examples
```bash
# Minimal (ID-based substrate) with explicit total ligand charge
pdb2reaction extract -i complex.pdb -c '123' -o pocket.pdb --ligand-charge -3

# Substrate provided as a PDB; per-resname charge mapping (others remain 0)
pdb2reaction extract -i complex.pdb -c substrate.pdb -o pocket.pdb --ligand-charge 'GPP:-3,SAM:1'

# Name-based substrate selection including all matches (WARNING is logged)
pdb2reaction extract -i complex.pdb -c 'GPP,SAM' -o pocket.pdb --ligand-charge 'GPP:-3,SAM:1'

# Multi-structure to single multi-MODEL output with hetero-hetero proximity enabled
pdb2reaction extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket_multi.pdb --radius-het2het 2.6 --ligand-charge 'GPP:-3,SAM:1'

# Multi-structure to multiple outputs with hetero-hetero proximity enabled
pdb2reaction extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket1.pdb pocket2.pdb --radius-het2het 2.6 --ligand-charge 'GPP:-3,SAM:1'
```

## Workflow
### Residue inclusion
- Always include the substrate residues from `-c/--center`.
- **Standard cutoff (`--radius`, default 2.6 Å):**
  - When `--exclude-backbone False`, any atom within the cutoff qualifies a residue.
  - When `--exclude-backbone True`, amino-acid residues must contact the substrate with a **non-backbone** atom (not N/H*/CA/HA*/C/O). Non-amino acids use any atom.
- **Independent hetero–hetero cutoff (`--radius-het2het`):** adds residues when a substrate hetero atom (non C/H) lies within the specified Å of a protein hetero atom. With backbone exclusion enabled the protein atom must be non-backbone.
- **Water handling:** HOH/WAT/H2O/DOD/TIP/TIP3/SOL are included by default (`--include-H2O True`).
- **Forced inclusion:** `--selected-resn` accepts IDs with chains/insertion codes (e.g., `A:123A`).
- **Neighbor safeguards:**
  - When backbone exclusion is off and a residue contacts the substrate with a backbone atom, auto-include the peptide-adjacent N/C neighbors (C–N ≤ 1.9 Å). Termini keep caps (N/H* or C/O/OXT).
  - Disulfide bonds (SG–SG ≤ 2.5 Å) bring both cysteines.
  - Non-terminal PRO residues always pull in the N-side amino acid; CA is preserved even if backbone atoms are removed, and when `--exclude-backbone True`, the neighbor’s C/O/OXT remain to maintain the peptide bond.

### Truncation/capping
- Isolated residues retain only side-chain atoms; amino-acid backbone atoms (N, CA, C, O, OXT plus N/CA hydrogens) are removed except for PRO/HYP safeguards.
- Continuous peptide stretches keep internal backbone atoms; only terminal caps (N/H* or C/O/OXT) are removed. TER awareness prevents capping across chain breaks.
- With `--exclude-backbone True`, main-chain atoms on all **non-substrate** amino acids are stripped (subject to PRO/HYP safeguards and PRO neighbor retention).
- Non-amino-acid residues never lose atoms named like backbone (N/CA/HA/H/H1/H2/H3).

### Link hydrogens (`--add-linkH True`)
- Adds carbon-only link hydrogens at 1.09 Å along severed bond vectors (CB–CA, CA–N, CA–C; PRO/HYP use CA–C only).
- Inserted after a `TER` as contiguous `HETATM` records named `HL` in residue `LKH` (chain `L`). Serial numbers continue from the main block.
- In multi-structure mode the same bonds are capped across all models; coordinates remain model-specific.

### Charge summary (`--ligand-charge`)
- Amino acids and common ions draw charges from internal dictionaries; waters are zero.
- Unknown residues default to 0 unless `--ligand-charge` supplies either a total charge (distributed across unknown substrate residues, or all unknowns when no unknown substrate) or a per-resname mapping like `GPP:-3,SAM:1`.
- Summaries (protein/ligand/ion/total) are logged for the first input when verbose mode is enabled.

### Multi-structure ensembles
- Accepts multiple input PDBs (identical atom ordering is validated at the head/tail of each file). Each structure is processed independently and the **union** of selected residues is applied to every model so that outputs remain consistent.
- Output policy:
  - No `-o`, multiple inputs → per-file `pocket_<original_basename>.pdb`.
  - One `-o` path → single multi-MODEL PDB.
  - N outputs where N == number of inputs → N individual PDBs.
- Diagnostics echo raw vs. kept atom counts per model along with residue IDs.

### Substrate specification (`-c/--center`)
- PDB path: the coordinates must match the first input exactly (tolerance 1e-3 Å); residue IDs propagate to other structures.
- Residue IDs: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'` (insertion codes supported).
- Residue names: comma-separated list (case insensitive). If multiple residues share a name, **all** matches are included and a warning is logged.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | One or more protein–ligand PDB files (identical atom ordering required). | Required |
| `-c, --center SPEC` | Substrate specification (PDB path, residue IDs, or residue names). | Required |
| `-o, --output PATH...` | Pocket PDB output(s). One path ⇒ multi-MODEL, N paths ⇒ per input. | Auto (`pocket.pdb` or `pocket_<input>.pdb`) |
| `-r, --radius FLOAT` | Atom–atom distance cutoff (Å) for inclusion. | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å, non C/H). | `0.0` (internally 0.001 Å when zero) |
| `--include-H2O {True\|False}` | Include HOH/WAT/H2O/DOD/TIP/TIP3/SOL waters. | `True` |
| `--exclude-backbone {True\|False}` | Remove backbone atoms on non-substrate amino acids (PRO/HYP safeguards). | `True` |
| `--add-linkH {True\|False}` | Add carbon-only link hydrogens at 1.09 Å along severed bonds. | `True` |
| `--selected-resn TEXT` | Force-include residues (IDs with optional chains/insertion codes). | `""` |
| `--ligand-charge TEXT` | Total charge or per-resname mapping (e.g., `GPP:-3,SAM:1`). | _None_ |
| `-v, --verbose` | Emit INFO-level logging (`true`) or keep warnings only (`false`). | `true` |

## Outputs
```
<output>.pdb  # Pocket PDB(s) with optional link hydrogens after a TER record
               # Single input → pocket.pdb by default
               # Multiple inputs without -o → pocket_<original_basename>.pdb per structure
               # One -o path with multiple inputs → single multi-MODEL PDB
               # Output directories are not created automatically; ensure they exist
```
- Charge summary (protein/ligand/ion/total) is logged for model #1 when verbose mode is enabled.
- Programmatic use (`extract_api`) returns `{"outputs": [...], "counts": [...], "charge_summary": {...}}`.

## Notes
- `--radius` defaults to 2.6 Å; `0` is nudged to 0.001 Å to avoid empty selections. `--radius-het2het` is off by default (also nudged to 0.001 Å when zero is provided).
- Waters can be excluded with `--include-H2O False`.
- Backbone trimming plus capping respect chain breaks and PRO/HYP safeguards as outlined above; non-amino residues never lose backbone-like atom names.
- Link hydrogens are inserted only on carbon cuts and reuse identical bonding patterns across models in ensemble mode.
- INFO logs summarize residue selection, truncation counts, and charge breakdowns.
