# `pdb2reaction energy-diagram`

## Purpose

Build an ad-hoc energy diagram from a list of state names + energy
values. Use this when you've collected energies from multiple
`pdb2reaction` runs (e.g. R from one run, TS / IM from another) and
want a single composite figure.

## Synopsis

```bash
pdb2reaction energy-diagram \
    --states 'R:0.0' 'TS1:21.5' 'IM:-0.7' 'TS2:2.2' 'P:-18.2' \
    [-o diagram.png] [--html]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `--states` | repeated `'NAME:VALUE'` | required | Each state name + energy value (kcal/mol relative to first) |
| `-o, --output` | path | `diagram.png` | Output figure |
| `--html` | flag | off | Plotly HTML instead of static PNG |
| `--units` | str | `kcal/mol` | y-axis units label |

## Examples

### Two-step energy diagram

```bash
pdb2reaction energy-diagram \
    --states 'R:0.0' 'TS1:21.5' 'IM:-0.7' 'TS2:2.2' 'P:-18.2' \
    -o diagram.png
```

## Caveats

- Energies are taken **as-is**; the command does no conversion. Make
  sure your values are in consistent units (kcal/mol typical) before
  calling.
- For a profile *along* a trajectory (continuous coordinate), use
  `trj2fig.md` instead.

## See also

- `trj2fig.md` — plot continuous trajectories.
- `pdb2reaction-workflows-output/SKILL.md` — extracting per-segment
  energies from `summary.json`.
