# `trj2fig`

Plot energies from an XYZ trajectory's comment lines (or recompute them with an MLIP backend) and export the series as static / interactive figures and CSV tables. Use it to visualize energies along an XYZ trajectory produced by `opt`, `scan`, `path-opt`, `path-search`, or `irc`. By default it reads the per-frame energies (in hartree) from each comment line, converts them to kcal/mol or hartree, and optionally references all values to a chosen frame to produce ΔE. Supplying `-q/--charge` and/or `-m/--multiplicity` instead recomputes every frame's energy with the MLIP backend (default UMA) using the given charge/spin.

## Examples

Command form:

```bash
pdb2reaction trj2fig -i TRAJECTORY.xyz [-o OUTPUTS...] [-q CHARGE] [-m MULT] [options]
```

```bash
# Default PNG, ΔE relative to the first frame
pdb2reaction trj2fig -i traj.xyz

# CSV + SVG with ΔE relative to frame 5, reported in hartree
pdb2reaction trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

# Multiple figure formats with the x-axis reversed (reference becomes last frame)
pdb2reaction trj2fig -i traj.xyz --reverse-x -o energy.png energy.html energy.pdf

# Recompute all frame energies with UMA before plotting
pdb2reaction trj2fig -i traj.xyz -q 0 -m 1 -o energy.png
```

## Workflow

1. Parse the XYZ trajectory. By default, read the energy from every frame
    comment: an explicit `E=`/`Energy:` token takes precedence, otherwise the lone numeric token is used (scientific notation such as `1.5e-3` is supported; with multiple bare numbers the LAST is taken and a warning is emitted). If
    `-q/-m` is present, recompute energies (in hartree) for each frame with
    the MLIP backend using those charge/spin values instead of the comment.
    If no energies are found or produced, the run aborts.
2. Normalize the reference specification:
 - `init` → frame `0` (or the last frame when `--reverse-x` is active).
 - `None`/`none`/`null` → absolute energies (no referencing).
 - Integer literal → the corresponding 0-based frame index.
3. Convert energies to either kcal/mol (default) or hartree and, when a
    reference is active, subtract the reference value to produce ΔE.
4. Build the Plotly figure (strong ticks, spline interpolation, markers, no
    title) and export it to every requested extension.
5. Optionally emit a CSV table with columns `frame`, `energy_hartree`, and the
    appropriate ΔE/E column in the requested unit.

## Outputs

```
<output>.[png|jpg|jpeg|html|svg|pdf] # Plotly export for every requested extension (defaults to energy.png)
<output>.csv                         # Optional energy table when CSV is requested
result.json                          # Machine-readable summary when --out-json is set
```
- CSV exports include `frame`, `energy_hartree`, and either a ΔE column
  (`delta_kcal`/`delta_hartree`) or an absolute column (`energy_kcal`/`energy_hartree`)
  when no reference is applied.
- Console diagnostics describing parsing failures or unsupported extensions.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | XYZ trajectory whose second line stores energies. | Required |
| `-o, --out PATH` | Repeatable output filenames; supports `.png/.jpg/.jpeg/.html/.svg/.pdf/.csv`. | `energy.png` |
| _extra arguments_ | Positional filenames listed after options; merged with the `-o` list. | _None_ |
| `--unit {kcal,hartree}` | Target unit for the plotted/exported values. | `kcal` |
| `-r, --reference TEXT` | Reference specification (`init`, `None`, or 0-based integer). | `init` |
| `-q, --charge INT` | Total charge; triggers energy recomputation when provided. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); triggers energy recomputation when provided. | _None_ |
| `--reverse-x/--no-reverse-x` | Reverse the x-axis so the last frame appears on the left (and `init` becomes the last frame). | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend for energy recomputation. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` next to the output(s). See [JSON Output Schema](json-output.md) for the schema. | `False` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes
- A comment with no parseable energy raises an error (see Workflow step 1 for the `E=`/`Energy:` vs. bare-number precedence).
- Unsupported extensions abort the run; `.png` uses Plotly’s PNG export with
  `scale=2` for sharper output.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [path-search](path-search.md) -- MEP trajectories for profiling
- [irc](irc.md) -- IRC trajectories for profiling
- [all](all.md) -- End-to-end workflow
