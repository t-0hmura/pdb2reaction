# `pdb2reaction energy-diagram`

## Purpose

Plot an energy diagram from numeric inputs. Use this when you have
energy values from one or more `pdb2reaction` runs (e.g. R from one
run, TS / IM from another) and want a single composite figure.

## Synopsis

```bash
pdb2reaction energy-diagram -i 0 12.5 4.3 [-i 18.2 -2.0] \
    [-o energy_diagram.png]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | numeric sequence | required | Energy values. Accepts space-separated numbers (`-i 0 12.5 4.3`), a Python-list literal (`-i "[0,12.5,4.3]"`), or repeated `-i` calls |
| `-o, --output` | path | `energy_diagram.png` | Output image path |
| `--label-x` | sequence | `S1, S2, ...` | Per-state labels on the x-axis |
| `--label-y` | str | `Energy (kcal/mol)` | Y-axis label |

Without `--label-x`, points are plotted in input order with `S1, S2, …` labels.

## Examples

### Five points (two-step mechanism)

```bash
pdb2reaction energy-diagram -i 0.0 21.5 -0.7 2.2 -18.2 -o diagram.png
```

### Bracketed list literal

```bash
pdb2reaction energy-diagram -i "[0.0, 21.5, -0.7, 2.2, -18.2]" -o diagram.png
```

## Caveats

- Energies are taken **as-is**; the command does no unit conversion.
  Make sure all values share a unit (kcal/mol is typical) before
  calling.
- For a profile along a continuous trajectory (XYZ frames with energies
  in the comment line), use `trj2fig.md`.
- Per-state labels are passed via `--label-x` and the y-axis caption via `--label-y`.

## See also

- `trj2fig.md` — plot from a trajectory file.
- `../pdb2reaction-workflows-output/SKILL.md` — extracting per-segment
  energies from `summary.json` to feed this command.