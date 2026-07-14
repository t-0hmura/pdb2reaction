# `pdb2reaction energy-diagram`

## Purpose

Plot an energy diagram from numeric inputs. Use this when you have
energy values from one or more `pdb2reaction` runs (e.g. R from one
run, TS / IM from another) and want a single composite figure.

## Synopsis

```bash
pdb2reaction energy-diagram -i "[0, 12.5, 4.3]" \
    [-o energy_diagram.png]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | numeric sequence | required | Energy values. Accepts a Python-list literal (`-i "[0, 12.5, 4.3]"`) or repeated `-i` calls (`-i 0 -i 12.5 -i 4.3`). A bare space-separated list after one flag (`-i 0 12.5 4.3`) is rejected |
| `-o, --output` | path | `energy_diagram.png` | Output image path |
| `--label-x` | sequence | `S1, S2, ...` | Per-state labels. Repeat the flag (`--label-x R --label-x TS`) or pass one list literal; the count must equal the energy count. |
| `--label-y` | str | `ΔE (kcal/mol)` | Y-axis label |
| `--out-json / --no-out-json` | flag | `--no-out-json` | Write `result.json` next to the image |

Without `--label-x`, points are plotted in input order with `S1, S2, …` labels.

## Examples

### Five points (two-step mechanism)

```bash
pdb2reaction energy-diagram -i 0.0 -i 21.5 -i -0.7 -i 2.2 -i -18.2 -o diagram.png
```

### Bracketed list literal

```bash
pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7, 2.2, -18.2]" -o diagram.png
```

### Explicit state labels and JSON

```bash
pdb2reaction energy-diagram -i "[0.0, 21.5, -18.2]" \
    --label-x R --label-x TS --label-x P \
    --label-y 'ΔG (kcal/mol)' --out-json -o diagram.png
```

With `--out-json`, `result.json` is written in the image directory with
`status`, `n_points`, and `files` (the resolved image path).

## Caveats

- Energies are taken **as-is**; the command does no unit conversion.
  Make sure all values share a unit (kcal/mol is typical) before
  calling.
- Do not write `-i 0 12.5 4.3` or `--label-x R TS P`; each option consumes
  one token. Repeat the option, or quote one list literal.
- For a profile along a continuous trajectory (XYZ frames with energies
  in the comment line), use `trj2fig.md`.

## See also

- `trj2fig.md` — plot from a trajectory file.
- `../pdb2reaction-workflows-output/SKILL.md` — extracting per-segment
  energies from `summary.json` to feed this command.
