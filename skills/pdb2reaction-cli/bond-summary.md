# `pdb2reaction bond-summary`

## Purpose

Detect bond changes between consecutive structures. Reports formed /
broken bonds with the same algorithm that `path-search` uses for
recursive segmentation and that `irc` reports under `bond_changes`.

Use it as a sanity check on R vs P, or to understand how a recursive
`path-search` decided where to split a multi-step mechanism.

## Synopsis

```bash
pdb2reaction bond-summary -i a.pdb -i b.pdb [-i c.pdb ...] \
    [--bond-factor 1.2] [--device cpu] [--one-based|--zero-based] [--json]
# Positional file args also accepted:
pdb2reaction bond-summary A.xyz B.xyz [C.xyz ...]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required (≥2) | Repeat `-i` for each structure (XYZ / PDB / GJF). Atom ordering must be identical across all inputs. |
| `--device` | str | `cpu` | Compute device for distance calculations |
| `--bond-factor` | float | `1.2` | Covalent-radius multiplier for bond cutoff |
| `--one-based / --zero-based` | flag | `--one-based` | Atom-index numbering convention in the report |
| `--json / --no-json` | flag | `--no-json` | Emit machine-readable JSON to stdout in place of the text report |

(Internal `margin_fraction` of 0.05 further shrinks the threshold; see
`bond_changes.py`.)

## Examples

```bash
# R vs P
pdb2reaction bond-summary -i 1.R.pdb -i 3.P.pdb

# Multi-frame check across an MEP
pdb2reaction bond-summary -i frame_01.xyz -i frame_05.xyz -i frame_10.xyz
```

## Output

Text on stdout, e.g.:

```
============================================================
  1.R.pdb  →  2.P.pdb
============================================================
Bond formed (2):
  - C320-C321 : 3.170 Å --> 1.680 Å
  - O186-H321 : 2.230 Å --> 0.980 Å
Bond broken (2):
  - S320-C320 : 1.800 Å --> 3.430 Å
  - C321-H321 : 1.100 Å --> 2.320 Å
```

## Caveats

- Atom ordering must match across all inputs. If it doesn't, run
  `pdb2reaction extract` first to canonicalize.
- The default `1.2` × covalent-radius cutoff (with internal margin
  fraction `0.05`) is geometry-only — it does not classify covalent
  vs ionic vs hydrogen-bonded; metal–ligand interactions may hover
  near the cutoff. Raise `--bond-factor` to 1.5 / 1.6 for permissive
  detection.
- Bond-change blocks are also embedded in `summary.json`'s per-segment
  output as a list of `{"Bond formed (k)": [...], "Bond broken (k)": [...]}`
  dicts (one per consecutive pair). See `../pdb2reaction-workflows-output/SKILL.md`.

## See also

- `irc.md` — same algorithm applied to IRC endpoints.
- `path-search.md` — uses bond changes for recursive segmentation.
- `../pdb2reaction-workflows-output/SKILL.md` — interpreting
  `summary.json["segments"][i]["bond_changes"]`.
- Defaults: bond-summary uses its own CLI defaults (`device=cpu`, `bond_factor=1.20`). `BOND_KW` (`device=auto`) only governs the internal bond checks in scan / path-search / all.