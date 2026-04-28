# `pdb2reaction trj2fig`

## Purpose

Plot an energy profile from an XYZ trajectory. Reads ASE-style `energy=...`
metadata from the comment line of each frame and produces a static PNG
or HTML plot. Useful for quickly visualizing IRC, MEP, or scan output.

## Synopsis

```bash
pdb2reaction trj2fig -i trajectory.xyz [-o energy.png|energy.html|energy.csv]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | XYZ trajectory with energy in comment line |
| `-o, --out` | path | `energy.png` | Output figure path; suffix selects format (`.png` / `.svg` / `.pdf` / `.jpg` / `.html` / `.csv`) |
| `--unit` | choice | `kcal` | `kcal` or `hartree` |
| `-r, --reference` | int / `init` / `None` | `init` | Reference frame for ΔE: `init` (initial frame; last frame if `--reverse-x`), `None` (absolute E), or a 1-based integer index |
| `-q, --charge` / `-m, --multiplicity` / `-b, --backend` / `--solvent` / `--solvent-model` | — | — | Recompute energies via MLIP if the input XYZ has no energies in its comment lines |
| `--reverse-x/--no-reverse-x` | flag | `--no-reverse-x` | Flip the x-axis |

## Examples

### Static PNG

```bash
pdb2reaction trj2fig -i finished_irc_trj.xyz -o irc_profile.png
```

### Interactive HTML (suffix selects format)

```bash
pdb2reaction trj2fig -i mep.xyz -o mep.html
```

## Caveats

- The XYZ comment line must encode the energy. The reader pulls the
  first numeric token (any decimal / scientific / negative form), so
  bare floats like `-12345.67` and ASE-style `... energy=-1234.56`
  both work. If the comment line has no numeric token at all, the
  reader raises (no silent flat plot).
- For a labeled energy diagram (R / TS / IM / P), use `energy-diagram.md`
  instead.

## See also

- `energy-diagram.md` — composed energy diagrams from explicit values.
- `irc.md`, `path-search.md`, `scan.md` — produce trajectories that
  feed `trj2fig`.
