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

_Subcommand `extract` is invoked via the Click wrapper in `pdb2reaction cli`, which forwards to the argparse-based implementation shown above._

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | One or more full protein–ligand PDB files (identical atom ordering required). | Required |
| `-c, --center SPEC` | Substrate specification: PDB path, residue-ID list (e.g. `123,124` or `A:123,B:456`), or residue-name list (e.g. `GPP,MMT`). | Required |
| `-o, --output PATH...` | Pocket PDB output(s). Provide one path for a multi-MODEL file or N paths matching the number of inputs. | Auto-generated (`pocket.pdb` / `pocket_<input>.pdb`) |
| `-r, --radius FLOAT` | Atom–atom distance cutoff (Å) for inclusion. | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å) for non-C/H pairs. | `0.0` (internally treated as 0.001 Å) |
| `--include-H2O BOOL` | Include water residues (HOH/WAT/TIP3/SOL). | `true` |
| `--exclude-backbone BOOL` | Remove backbone atoms (N/H*/CA/HA*/C/O) on non-substrate amino acids (PRO/HYP safeguards). | `true` |
| `--add-linkH BOOL` | Add carbon-only link hydrogens at 1.09 Å along severed bonds. | `true` |
| `--selected-resn TEXT` | Residue IDs to force-include (comma/space separated; chain/insertion codes supported). | `""` |
| `--ligand-charge TEXT` | Either a total charge (number) or mapping like `GPP:-3,MMT:-1` to distribute across unknown residues. | `None` |
| `-v, --verbose` | Enable INFO-level logging. | `true` |

## Outputs
- Pocket PDB(s) containing the extracted residues, with optional link hydrogens appended after a `TER` record.
- Charge summary (protein/ligand/ion/total) printed to stdout when extraction succeeds.

## Notes
- Multi-structure mode (`-i` with multiple files) produces a pocket per input; atom ordering must match across inputs.
- Substrate residue lists support insertion codes (e.g. `123A`, `A:123A`). When residue-name mode yields duplicates, all matches are included and a warning is logged.
- Water selection recognises HOH/WAT/TIP3/SOL residues; backbone trimming skips PRO/HYP nitrogen atoms.
