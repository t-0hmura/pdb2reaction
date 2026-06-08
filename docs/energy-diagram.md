# `energy-diagram`

Draw a state energy diagram directly from numeric values you provide — it does not read structure files and runs no quantum / thermo / MLIP calculation (`--thermo` / `--dft`). Use it when numeric state energies are already known (e.g. from `summary.json`) and only the formatted diagram is needed; it produces one image file and an optional machine-readable sidecar.

## Examples

```bash
# Command form
pdb2reaction energy-diagram -i VALUES... [-o OUTPUT] [--label-x...] [--label-y...]
```

```bash
# List string (recommended for ad-hoc plots)
pdb2reaction energy-diagram -i "[0, 12.5, 4.3]" -o energy.png
```

```bash
# Repeated flag form (use one `-i` per value)
pdb2reaction energy-diagram -i 0 -i 12.5 -i 4.3 -o energy.png
```

```bash
# X/Y labels — same applies to --label-x / --label-y
pdb2reaction energy-diagram -i "[0, 12.5, 4.3]" \
    --label-x "['R','TS','P']" --label-y "ΔE (kcal/mol)" -o energy.png
```

## Workflow
1. Collect values from `-i/--input` (supports repeated flags and a single list-like string).
2. Parse all input values as floats and fail early if fewer than two values are provided.
3. Parse optional `--label-x` values. If omitted, labels are auto-generated as `S1`, `S2`,...
4. Validate label count (`--label-x`) against value count, then render the diagram.
5. Save the image to `-o/--output` and print the saved path.

## Outputs
```
OUTPUT.(png|jpg|jpeg|svg|pdf)
result.json   # Optional sidecar with status, n_points, files plus standard envelope fields (command/version/schema/environment); no per-point energies or labels (when --out-json is set)
```
- If `-o/--output` is omitted, `energy_diagram.png` is written to the current directory.
- When output has no extension, `.png` is appended automatically.
- Parent directories are created automatically when needed.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input TEXT...` | Numeric values (multiple args or list-like string). | Required |
| `-o, --output PATH` | Output image path (`.png/.jpg/.jpeg/.svg/.pdf`). | `energy_diagram.png` |
| `--label-x TEXT...` | X-axis state labels. Count must match input value count. | `S1, S2,...` |
| `--label-y TEXT` | Y-axis label. | `ΔE (kcal/mol)` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` next to the output image. See [JSON Output Schema](json-output.md) for the schema. | `False` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes
- Input order is used directly as plotting order.
- At least two numeric values are required.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [trj2fig](trj2fig.md) -- Plot profile from trajectory energies
- [all](all.md) -- End-to-end workflow with built-in energy diagram output
