# `trj2fig` subcommand

## Overview
`trj2fig` converts an XYZ trajectory into polished energy profiles. By default it
reads the Hartree energies encoded in each frame’s comment line, converts them
to kcal/mol or Hartree, optionally references all values to a chosen frame, and
exports the resulting series as static/interactive figures and CSV tables. The
reference can be the first frame (`init`), the last frame when `--reverse-x` is
used, or any explicit index. When you supply `-q/--charge` and/or
`-m/--multiplicity`, all energies are recomputed for every frame with the
`uma_pysis` calculator using the provided charge/spin instead of the comment
lines.

## Usage
```bash
pdb2reaction trj2fig -i TRAJECTORY.xyz [-o OUTPUTS...] [-q CHARGE] [-m MULT] [options]
```

### Examples
```bash
# Default PNG, ΔE relative to the first frame
pdb2reaction trj2fig -i traj.xyz

# CSV + SVG with ΔE relative to frame 5, reported in Hartree
pdb2reaction trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

# Multiple figure formats with the x-axis reversed (reference becomes last frame)
pdb2reaction trj2fig -i traj.xyz --reverse-x -o energy.png energy.html energy.pdf

# Recompute all frame energies with UMA before plotting
pdb2reaction trj2fig -i traj.xyz -q 0 -m 1 -o energy.png
```

## Workflow
1. Parse the XYZ trajectory. By default, read the first floating-point number
   found in every frame comment (scientific notation is not supported). If
   `-q/-m` is present, recompute Hartree energies for each frame with
   `uma_pysis` using those charge/spin values instead of the comment.
   If no energies are found or produced, the run aborts.
2. Normalize the reference specification:
   - `init` → frame `0` (or the last frame when `--reverse-x` is active).
   - `None`/`none`/`null` → absolute energies (no referencing).
   - Integer literal → the corresponding 0-based frame index.
3. Convert energies to either kcal/mol (default) or Hartree and, when a
   reference is active, subtract the reference value to produce ΔE.
4. Build the Plotly figure (strong ticks, spline interpolation, markers, no
   title) and export it to every requested extension.
5. Optionally emit a CSV table with columns `frame`, `energy_hartree`, and the
   appropriate ΔE/E column in the requested unit.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | XYZ trajectory whose second line stores energies. | Required |
| `-o, --out PATH` | Repeatable output filenames; supports `.png/.jpg/.jpeg/.html/.svg/.pdf/.csv`. | `energy.png` |
| _extra arguments_ | Positional filenames listed after options; merged with the `-o` list. | _None_ |
| `--unit {kcal,hartree}` | Target unit for the plotted/exported values. | `kcal` |
| `-r, --reference TEXT` | Reference specification (`init`, `None`, or 0-based integer). | `init` |
| `-q, --charge INT` | Total charge; triggers energy recomputation with `uma_pysis` when provided. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); triggers energy recomputation with `uma_pysis` when provided. | _None_ |
| `--reverse-x` | Reverse the x-axis so the last frame appears on the left (and `init` becomes the last frame). | `False` |

## Outputs
```
<output>.[png|jpg|jpeg|html|svg|pdf]  # Plotly export for every requested extension (defaults to energy.png)
<output>.csv                          # Optional energy table when CSV is requested
```
- When no `-o` or positional outputs are provided, a single `energy.png` is written
  to the current directory. CSV exports include `frame`, `energy_hartree`, and either
  a ΔE column (`delta_kcal`/`delta_hartree`) or absolute column (`energy_kcal`/`energy_hartree`
  when no reference is applied).
- Console diagnostics describing parsing failures or unsupported extensions.

## Notes
- Energies are taken from the first decimal number in each comment; malformed
  comments raise an error.
- Unsupported extensions abort the run; `.png` uses Plotly’s PNG export with
  `scale=2` for sharper output.
- `--reverse-x` flips both the axis direction and the behavior of `-r init` so
  the visualized path reads backwards.
- The legacy `--output-peak` option was removed.
