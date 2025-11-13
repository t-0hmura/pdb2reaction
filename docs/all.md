# `all` subcommand

## Purpose
Runs the full pipeline: pocket extraction across multiple full PDBs, minimum-energy-path search on the pockets, automatic merge back into the originals, and optional per-segment post-processing (TS optimisation, thermochemistry, DFT).

## Usage
```bash
pdb2reaction all -i R.pdb [I.pdb ...] P.pdb -c SUBSTRATE_SPEC
                 [--ligand-charge MAP_OR_NUMBER]
                 [--spin 2S+1]
                 [--freeze-links/--no-freeze-links]
                 [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
                 [--sopt-mode lbfgs|rfo|light|heavy]
                 [--dump/--no-dump] [--out-dir DIR]
                 [--pre-opt/--no-pre-opt]
                 [--args-yaml FILE]
                 [--tsopt/--no-tsopt] [--thermo/--no-thermo] [--dft/--no-dft]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full PDBs in reaction order. A single `-i` may be followed by multiple files. | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names) passed to the extractor. | Required |
| `--out-dir PATH` | Top-level output directory. | `./result_all/` |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). | `0.0` |
| `--include-H2O / --no-include-H2O` | Include waters. | `--include-H2O` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate residues. | `--exclude-backbone` |
| `--add-linkH / --no-add-linkH` | Add carbon-only link hydrogens. | `--add-linkH` |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--ligand-charge TEXT` | Total charge or mapping for unknown residues (recommended). | `None` |
| `--verbose / --no-verbose` | Extractor logging. | `--verbose` |
| `-s, --spin INT` | Spin multiplicity forwarded to path search and post-processing. | `1` |
| `--freeze-links / --no-freeze-links` | Freeze link parents in pocket PDBs. | `--freeze-links` |
| `--max-nodes INT` | GSM internal nodes per segment. | `10` |
| `--max-cycles INT` | GSM max cycles. | `100` |
| `--climb / --no-climb` | Enable climbing image for the first segment per pair. | `--climb` |
| `--sopt-mode TEXT` | Single-structure optimiser for HEI±1/kink nodes. | `lbfgs` |
| `--dump / --no-dump` | Dump GSM and single-structure trajectories. | `--no-dump` |
| `--args-yaml FILE` | YAML forwarded to `path_search` (see below). | _None_ |
| `--pre-opt / --no-pre-opt` | Pre-optimise pocket endpoints before GSM. | `--pre-opt` |
| `--tsopt / --no-tsopt` | Run TS optimisation and pseudo-IRC per reactive segment. | `--no-tsopt` |
| `--thermo / --no-thermo` | Run vibrational analysis (freq) on R/TS/P and build Gibbs diagram. | `--no-thermo` |
| `--dft / --no-dft` | Run single-point DFT on R/TS/P and build DFT energy diagram. | `--no-dft` |

## YAML configuration (`--args-yaml`)
The YAML file is forwarded unchanged to `path_search`. See [`path_search`](path_search.md#yaml-configuration-args-yaml) for accepted sections (`geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search`).

## Outputs
- `<out-dir>/pockets/`: Pocket PDBs for each input.
- `<out-dir>/path_search/`: GSM results (trajectory, merged full-system PDBs, energy diagrams, `summary.yaml`, segment folders).
- Optional `<out-dir>/path_search/tsopt_seg_XX/` subtrees with TS optimisation, pseudo-IRC, frequency, and DFT results depending on toggles.
- Console logs covering pocket charge summary, resolved configuration blocks, and per-stage timing.

## Notes
- The total pocket charge from the extractor (first model) is rounded to the nearest integer and used as the GSM charge.
- Reference PDB templates for merging are taken automatically from the original inputs; the explicit `--ref-pdb` option of `path_search` is intentionally hidden in this wrapper.
- When both `--thermo` and `--dft` are enabled, the post-processing stage also produces a DFT//UMA Gibbs diagram (DFT energy + UMA thermal correction).
- Always provide `--ligand-charge` when formal charges are not inferable to ensure the correct total charge propagates through the pipeline.
