# `trj2fig` subcommand

## Purpose
Extracts energies from XYZ trajectory comments, computes either absolute energies or ΔE relative to a reference frame, and exports figures and/or CSV tables.

## Usage
```bash
pdb2reaction trj2fig -i TRAJ.xyz [-o OUT ...] [--unit kcal|hartree]
                     [-r init|None|INDEX] [--reverse-x]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | XYZ trajectory file whose comment line (2nd line) contains the energy. | Required |
| `-o, --out PATH` | Output filenames. Repeatable; supports .png/.html/.svg/.pdf/.csv. | _None_ → `energy.png` |
| `extra_outs` | Positional filenames after options; combined with `-o`. | _None_ |
| `--unit CHOICE` | Energy unit for output values (`kcal` or `hartree`). | `kcal` |
| `-r, --reference TEXT` | Reference specification: `init`, `None`, or integer frame index. | `init` |
| `--reverse-x` | Reverse x-axis so the last frame appears on the left (also flips `init`). | `False` |

## YAML configuration (`--args-yaml`)
Not supported.

## Outputs
- Figure files (`.png`, `.html`, `.svg`, `.pdf`) showing energy or ΔE versus frame.
- Optional `.csv` with columns `frame`, `energy_hartree`, and either `delta_kcal`, `energy_kcal`, `delta_hartree`, or `energy_hartree` depending on mode.
- Console warning/error messages when energies cannot be parsed or unsupported extensions are requested.

## Notes
- Energies are parsed as the first floating-point number on each comment line; malformed comments raise an error.
- The CSV/figure data respect `--unit` and whether a reference is used.
- `--reverse-x` swaps the interpretation of `-r init` to use the last frame as the reference.
