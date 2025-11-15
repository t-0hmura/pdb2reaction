# `extract` subcommand

## Purpose
Extracts binding pockets from protein–substrate complexes by selecting residues within distance cutoffs, optionally pruning backbone atoms, and adding link hydrogens.

## Usage
```bash
pdb2reaction extract -i COMPLEX.pdb [COMPLEX2.pdb ...]
                     -c SUBSTRATE_SPEC
                     [-o POCKET.pdb [POCKET2.pdb ...]]
                     [--radius Å] [--radius-het2het Å]
                     [--include-H2O true|false]
                     [--exclude-backbone true|false]
                     [--add-linkH true|false]
                     [--selected-resn LIST]
                     [--ligand-charge MAP_OR_NUMBER]
                     [--verbose true|false]
```

`--verbose true|false` mirrors the Click wrapper (`pdb2reaction cli extract`) and now behaves identically when the script is invoked directly via argparse.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | One or more full protein–ligand PDB files (identical atom ordering required). | Required |
| `-c, --center SPEC` | Substrate specification: PDB path, residue-ID list (e.g. `123,124` or `A:123,B:456`), or residue-name list (e.g. `GPP,MMT`). | Required |
| `-o, --output PATH...` | Pocket PDB output(s). Provide one path for a single multi-MODEL file or N paths matching the number of inputs for per-structure outputs. | Auto-generated (`pocket.pdb` for single input / `pocket_<input>.pdb` per input) |
| `-r, --radius FLOAT` | Atom–atom distance cutoff (Å) for inclusion. | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å) for non-C/H pairs. | `0.0` (internally treated as 0.001 Å) |
| `--include-H2O BOOL` | Include water residues (HOH/WAT/TIP3/SOL). | `true` |
| `--exclude-backbone BOOL` | Remove backbone atoms (N/H*/CA/HA*/C/O) on non-substrate amino acids (PRO/HYP safeguards). | `true` |
| `--add-linkH BOOL` | Add carbon-only link hydrogens at 1.09 Å along severed bonds. | `true` |
| `--selected-resn TEXT` | Residue IDs to force-include (comma/space separated; chain/insertion codes supported). | `""` |
| `--ligand-charge TEXT` | Either a total charge (number) or mapping like `GPP:-3,MMT:-1` to distribute across unknown residues. | `None` |
| `-v, --verbose` | Enable INFO-level logging (set `true` to emit INFO via `click.echo`, `false` keeps WARNING-only output). | `false` |

## Outputs
- Pocket PDB(s) containing the extracted residues, with optional link hydrogens appended after a `TER` record.
  - Single input: defaults to `pocket.pdb`.
  - Multiple inputs + no `-o`: defaults to `pocket_<original_basename>.pdb` per structure.
  - Supplying one `-o` path with multiple inputs writes a single multi-MODEL PDB.
- Charge summary (protein/ligand/ion/total) is logged for model #1 when verbose mode is enabled (Click wrapper echoes INFO lines by default; the argparse script now follows the same behaviour via `--verbose true`).

## Notes
- Multi-structure mode (`-i` with multiple files) unions the selected residues across all inputs (after verifying identical atom ordering). That union is applied consistently to every structure, and the outputs can be written either as a single multi-MODEL PDB or one file per input.
- `--radius` defines the standard cutoff around substrate atoms; `--radius-het2het` independently includes residues whose hetero atoms approach substrate hetero atoms (non-C/H) within the specified Å.
- Waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL) are included by default. `--include-H2O false` removes them from the pocket.
- Backbone trimming removes N/H*/CA/HA*/C/O from non-substrate amino acids when `--exclude-backbone true` (default); PRO/HYP safeguards keep the atoms needed to retain the pyrrolidine ring and neighbouring peptide bond.
- Substrate residue lists support insertion codes (e.g. `123A`, `A:123A`). When residue-name mode yields duplicates, all matches are included and a warning is logged.
- Link-hydrogen placement (`--add-linkH true`) adds 1.09 Å C–H vectors for severed bonds and appends them as contiguous `HL` atoms in residue `LKH` after a `TER`.
