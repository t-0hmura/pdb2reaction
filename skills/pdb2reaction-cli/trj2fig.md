# `pdb2reaction trj2fig`

## Purpose

Plot an energy profile from an XYZ trajectory. Reads an explicit `E=`/`Energy:`
token or a bare numeric energy from each frame's comment line and produces a static PNG
or HTML plot. Useful for quickly visualizing IRC, MEP, or scan output.

## Synopsis

```bash
pdb2reaction trj2fig -i trajectory.xyz \
    [-o energy.png] [-o energy.html] [-o energy.csv] [--out-json]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | XYZ trajectory with energy in comment line |
| `-o, --out` | path(s) | `energy.png` | Output path; repeat `-o` to write several formats (`.png` / `.svg` / `.pdf` / `.jpg` / `.jpeg` / `.html` / `.csv`) |
| `--unit` | choice | `kcal` | `kcal` or `hartree` |
| `-r, --reference` | int / `init` / `None` | `init` | Reference frame for ΔE: `init` (initial frame; last frame if `--reverse-x`), `None` (absolute E), or a 0-based integer index |
| `-q, --charge` / `-m, --multiplicity` / `-b, --backend` / `--solvent` / `--solvent-model` | — | — | Recompute every frame's energy via the MLIP backend when `-q/--charge` and/or `-m/--multiplicity` are supplied (instead of reading the comment line) |
| `--reverse-x/--no-reverse-x` | flag | `--no-reverse-x` | Flip the x-axis |
| `--out-json / --no-out-json` | flag | `--no-out-json` | Write `result.json` in the first output's directory |

## Examples

### Static PNG

```bash
pdb2reaction trj2fig -i finished_irc_trj.xyz -o irc_profile.png
```

### Interactive HTML (suffix selects format)

```bash
pdb2reaction trj2fig -i mep.xyz -o mep.html
```

### Several formats plus JSON

```bash
pdb2reaction trj2fig -i mep.xyz \
    -o mep.png -o mep.html -o mep.csv --out-json
```

`result.json` contains `status`, `n_frames`, the minimum/maximum Hartree
energies, `energy_source`, calculator provenance, and a `files` mapping for
every output. With comment-line energies, `energy_source` is
`trajectory_comment` and backend/charge/multiplicity/solvent fields are null;
the default `-b uma` is not falsely recorded as a calculator that ran. With
`-q` or `-m`, the source is `mlip_recomputed` and those fields record the
effective recomputation request.

## Caveats

- The XYZ comment line must encode the energy. An explicit keyed token
  (`E=` / `Energy:`, any decimal / scientific / negative form) takes
  precedence; otherwise exactly one bare numeric token is accepted. Several
  bare numbers are rejected as ambiguous. So
  bare floats like `-12345.67` and ASE-style `... energy=-1234.56`
  both work. If the comment line has no numeric token at all, the
  reader raises an error (rather than emitting a silent flat plot).
- Supplying either `-q` or `-m` switches from comment-line energies to an
  MLIP recomputation of every frame. Confirm the omitted charge or
  multiplicity default is appropriate; for a charged/open-shell system,
  normally provide both.
- For a labeled energy diagram (R / TS / IM / P), use `energy-diagram.md`
  instead.

## See also

- `energy-diagram.md` — composed energy diagrams from explicit values.
- `irc.md`, `path-search.md`, `scan.md` — produce trajectories that
  feed `trj2fig`.
