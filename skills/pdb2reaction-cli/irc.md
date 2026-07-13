# `pdb2reaction irc`

## Purpose

Intrinsic Reaction Coordinate (IRC) integration from a TS geometry.
Uses **EulerPC** (the only supported integrator; not exposed as a CLI
flag). Forward and backward branches are run from there; output
is a stitched IRC trajectory and raw endpoint geometries. (For LBFGS-refined
R / P, use `pdb2reaction all`.)

## Synopsis

```bash
pdb2reaction irc -i ts.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--max-cycles 125] [--step-size 0.1] [--never-stop] \
    [-b uma|orb|mace|aimnet2] [-o ./result_irc/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Optimized TS geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `--max-cycles` | int | 125 | Max IRC steps per branch (forward + backward) |
| `--step-size` | float | `0.10` | Step in Bohr (unweighted Cartesian); maps to `IRC_KW["step_length"]` |
| `--never-stop / --no-never-stop` | bool | `False` | Ignore energy-rise/plateau stops so a small shoulder can be crossed; gradient/integrator convergence and `max_cycles` still apply |
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

If IRC stops after only a few frames, **reduce `--step-size` first**
(usually from `0.10` to `0.05`). If the trajectory visibly contains a small
uphill/flat shoulder, retry with `--never-stop`; this is opt-in and the
trajectory/endpoints must be inspected because it may pass the nearest basin.

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/forward_irc_trj.xyz`, `backward_irc_trj.xyz` | always | IRC forward / backward trajectories |
| `<out_dir>/finished_irc_trj.xyz` | always | stitched: backward (reversed) + TS + forward |
| `<out_dir>/finished_first.xyz`, `finished_last.xyz` | always | raw IRC endpoints (backward / forward) |
| `<out_dir>/{forward,backward,finished}_irc.pdb` | input is `.pdb` | PDB companions |
| `<out_dir>/result.json` | `--out-json` | machine-readable result |

`result.json` keys:

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_reactant_hartree"], d["energy_ts_hartree"], d["energy_product_hartree"])
print(d["bond_changes"])           # {"formed": [...], "broken": [...]}
print(d["status"])                  # "completed" (success) / "error" (on failure)
```

For LBFGS-optimized canonical `reactant.{xyz,pdb}` / `product.{xyz,pdb}`
under `<out_dir>/segments/seg_NN/`, run `pdb2reaction all` (it calls `irc`
internally and post-processes the endpoints).

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
- Defaults: `import pdb2reaction.core.defaults as d; print(d.IRC_KW)`
