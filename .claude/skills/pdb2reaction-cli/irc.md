# `pdb2reaction irc`

## Purpose

Intrinsic Reaction Coordinate (IRC) integration from a TS geometry.
Uses **EulerPC** in mass-weighted Cartesians (the only supported
integrator; not exposed as a CLI flag). Forward and backward branches
are run, then each endpoint is LBFGS-optimized to the nearest minimum.
Output: a stitched IRC trajectory plus the optimized R and P geometries.

## Synopsis

```bash
pdb2reaction irc -i ts.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--max-cycles 125] [--step-size 0.1] \
    [-b uma|orb|mace|aimnet2] [-o ./result_irc/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Optimized TS geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--max-cycles` | int | 125 | Max IRC steps per branch (forward + backward) |
| `--step-size` | float | (live default) | Step in Bohr; check `IRC_KW.step_size` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_irc/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default IRC from a tsopt'd geometry

```bash
pdb2reaction irc -i result_tsopt/final_geometry.xyz -q 0 -m 1 -b uma --out-json -o result_irc
```

`--out-json` enables the `result.json` examples below; omit it for trajectory files only.

### Tighter step / longer integration for shallow surfaces

```bash
pdb2reaction irc -i ts.xyz -q -1 -m 1 \
    --max-cycles 250 --step-size 0.05 \
    -b uma -o result_irc_long
```

## Output

```
result_irc/
├── forward_irc_trj.xyz             # IRC forward trajectory
├── backward_irc_trj.xyz            # IRC backward trajectory
├── finished_irc_trj.xyz            # stitched: backward (reversed) + TS + forward
├── finished_first.xyz              # backward endpoint (raw IRC)
├── finished_last.xyz               # forward endpoint (raw IRC)
├── forward_irc.pdb                 # PDB companions (only when input is .pdb)
├── backward_irc.pdb
├── finished_irc.pdb
└── result.json                     # only when --out-json is passed
```

`result.json` keys:

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_reactant_hartree"], d["energy_ts_hartree"], d["energy_product_hartree"])
print(d["bond_changes"])           # {"formed": [...], "broken": [...]}
print(d["status"])                  # "completed" / "diverged" / ...
```

## R/TS/P canonical geometries (downstream — `pdb2reaction all`)

Standalone `irc` writes only the raw IRC endpoints (`finished_first.xyz`
/ `finished_last.xyz`). The post-IRC LBFGS-optimized minima and the
canonical `reactant.{xyz,pdb}` / `product.{xyz,pdb}` filenames live
under `<out_dir>/seg_NN/` and are produced by the `pdb2reaction all`
pipeline (which calls IRC internally and then optimizes the endpoints
under `path_search/post_seg_NN/`); see
`pdb2reaction-workflows-output/SKILL.md`.

## Bond-change check

`bond_changes` records which bonds are different between R and P
according to a 1.20× covalent-radius cutoff (with internal margin
0.05). This is the same algorithm used by `bond-summary` and
`path-search` segmentation.

```python
import json
bc = json.load(open("result_irc/result.json"))["bond_changes"]
for b in bc["formed"]: print("FORMED ", b)
for b in bc["broken"]: print("BROKEN ", b)
```

## Caveats

- IRC starts from a **single imaginary mode** TS. If `tsopt` produced
  multiple imaginary modes, IRC may follow the wrong one — re-tsopt
  first.
- `--max-cycles 125` is enough for most clusters. If forward / backward
  hits the cap, the surface is probably very shallow; try a smaller
  `--step-size`.
- The bond-change detector is geometry-based (covalent-radius cutoff),
  not physics-based. Metal–ligand bonds may flicker on the borderline.

## See also

- `tsopt.md` — produces the IRC starting geometry.
- `freq.md`, `dft.md` — downstream.
- `bond-summary.md` — same bond-change algorithm, standalone.
- `pdb2reaction-workflows-output/SKILL.md` — R/TS/P path conventions.
- Defaults: `import pdb2reaction.defaults as d; print(d.IRC_KW)`
