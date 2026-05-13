# `energy-diagram`

## Overview

Draw a state energy diagram directly from numeric values (no structure file, no quantum/thermo calculation).

### At a glance
- **Use when:** Numeric state energies are already known (e.g. from `summary.json`) and only the formatted diagram is needed.
- **Method:** Matplotlib renderer — no structure file is read, no QM / thermo / MLIP call is performed.
- **Outputs:** One image file (`.png` / `.jpg` / `.jpeg` / `.svg` / `.pdf`); optional `result.json` with `--out-json`. Default output is `energy_diagram.png`.
- **Defaults:** `-o energy_diagram.png`, `--label-x` auto-numbered (`S1`, `S2`, ...), `--label-y "ΔE (kcal/mol)"`, `--out-json False`.
- **Next step:** Pair with [`trj2fig`](trj2fig.md) when an energy trajectory is also on hand, or with [`all`](all.md) outputs (`summary.json`) for end-to-end pipelines.

`pdb2reaction energy-diagram` only visualizes numbers you provide. It does not read structure files and does not run `--thermo` / `--dft`-style calculations.

## Usage
```bash
pdb2reaction energy-diagram -i VALUES... [-o OUTPUT] [--label-x...] [--label-y...]
```

## Examples
```bash
# Multiple numeric arguments
pdb2reaction energy-diagram -i 0 12.5 4.3 -o energy.png

# List string
pdb2reaction energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png

# X/Y labels
pdb2reaction energy-diagram -i 0 12.5 4.3 --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## Workflow
1. Collect values from `-i/--input` (supports repeated flags, multiple values after one flag, and list-like strings).
2. Parse all input values as floats and fail early if fewer than two values are provided.
3. Parse optional `--label-x` values. If omitted, labels are auto-generated as `S1`, `S2`,...
4. Validate label count (`--label-x`) against value count, then render the diagram.
5. Save the image to `-o/--output` and print the saved path.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input TEXT...` | Numeric values (multiple args or list-like string). | Required |
| `-o, --output PATH` | Output image path (`.png/.jpg/.jpeg/.svg/.pdf`). | `energy_diagram.png` |
| `--label-x TEXT...` | X-axis state labels. Count must match input value count. | `S1, S2,...` |
| `--label-y TEXT` | Y-axis label. | `ΔE (kcal/mol)` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` next to the output image. See [JSON Output Schema](json-output.md) for the schema. | `False` |

## Outputs
```
OUTPUT.(png|jpg|jpeg|svg|pdf)
result.json   # Optional sidecar with input energies and labels (when --out-json is set)
```
- If `-o/--output` is omitted, `energy_diagram.png` is written to the current directory.
- When output has no extension, `.png` is appended automatically.
- Parent directories are created automatically when needed.

## Notes
- Input order is used directly as plotting order.
- At least two numeric values are required.
- This command does not read structure files and does not compute energies.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [trj2fig](trj2fig.md) -- Plot profile from trajectory energies
- [all](all.md) -- End-to-end workflow with built-in energy diagram output
