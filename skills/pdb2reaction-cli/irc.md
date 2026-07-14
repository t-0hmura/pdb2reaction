# `pdb2reaction irc`

## Purpose

Intrinsic Reaction Coordinate (IRC) integration from a TS geometry.
Uses **EulerPC** (the only supported integrator; not exposed as a CLI
flag). Forward and backward branches are run from there; output
is a stitched IRC trajectory and raw endpoint geometries. (For RFO-refined
R / P, use `pdb2reaction all`.)

## Synopsis

```bash
pdb2reaction irc -i ts.{pdb,cif,xyz,gjf} \
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
| `--tr-projection` | str | `constrained` | Rigid-mode treatment for the initial Hessian; `legacy-active` is for isolated-active comparison only |
| `--workers`, `--workers-per-node` | int | `1`, `1` | UMA predictor workers. `workers > 1` plus an explicit `Analytical` Hessian raises `BackendError`; use one worker or finite differences. |
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
| `<out_dir>/forward_irc_trj.xyz`, `backward_irc_trj.xyz` | corresponding direction enabled and produces output | directional IRC trajectories |
| `<out_dir>/finished_irc_trj.xyz` | IRC produces stitched output | stitched path from first endpoint through TS to last endpoint |
| `<out_dir>/finished_first.xyz`, `finished_last.xyz` | stitched trajectory has endpoints | first / last raw stitched-path endpoints; do not infer R/P identity without comparison |
| `<out_dir>/{forward,backward,finished}_irc.pdb` | `--convert-files` and PDB/mmCIF topology/reference available | PDB companions |
| `<out_dir>/{forward,backward,finished}_irc.cif` | `--convert-files` and input/reference required the mmCIF or oversized-PDB bridge | public trajectories with original IDs |
| `<out_dir>/result.json` | `--out-json` | machine-readable result |

`result.json` keys:

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_first_hartree"], d["energy_ts_hartree"], d["energy_last_hartree"])
print(d.get("bond_changes"))        # omitted if endpoint comparison was unavailable/failed
print(d["endpoint_energy_orientation"])
print(d.get("bond_changes_direction"))
print(d["status"])                 # "completed" means the runner returned and wrote output
print(d["forward_converged"], d["backward_converged"])
print(d["forward_energy_increased"], d["backward_energy_increased"])
print(d["rigid_projection"]["treatment"], d["rigid_projection"]["effective_rank"])
```

`completed` is not an IRC convergence verdict. Check the enabled direction's
`*_converged`, `*_energy_increased`, frame count, endpoints, and bond changes.
An energy-stop can produce `completed` with `*_converged == false`.
The older `energy_reactant_hartree` / `energy_product_hartree` keys are retained
as directional first/last aliases only. Standalone IRC has no endpoint
references, so compare `finished_first.xyz` and `finished_last.xyz` with the
intended states before assigning chemical R/P labels.
`endpoint_energy_orientation` and, when present, `bond_changes_direction` are
both `finished_first_to_finished_last`.

For RFO-optimized canonical `reactant.{xyz,pdb}` / `product.{xyz,pdb}`
under `<out_dir>/segments/seg_NN/`, run `pdb2reaction all` (it calls `irc`
internally and post-processes the endpoints). Set `opt-mode-post: grad` in the
`all` YAML only when you deliberately want LBFGS endpoint refinement.

## Bond-change check

`bond_changes` records the directed change from `finished_first.xyz` to
`finished_last.xyz`
according to a 1.20× covalent-radius cutoff (with internal margin
0.05). This is the same algorithm used by `bond-summary` and
`path-search` segmentation. This direction is not necessarily chemical R→P;
after structural assignment, reverse formed/broken labels if the chemical
orientation is last→first.

```python
import json
bc = json.load(open("result_irc/result.json")).get("bond_changes")
if bc is None:
    raise RuntimeError("IRC endpoint comparison was not available")
for b in bc["formed"]: print("FORMED ", b)
for b in bc["broken"]: print("BROKEN ", b)
```

## Caveats

- IRC starts from a **single imaginary mode** TS. If `tsopt` produced
  multiple imaginary modes, IRC may follow the wrong one — re-tsopt
  first.
- If a forward/backward branch reaches `--max-cycles`, inspect its gradient,
  energy trend, and geometry. A smaller `--step-size` can help integration
  stability but may require more cycles; there is no universal cycle count for
  every surface.
- The bond-change detector is geometry-based (covalent-radius cutoff),
  not physics-based. Metal–ligand bonds may flicker on the borderline.
- In the in-process `all` workflow, IRC may reuse TSOPT's cached Hessian only
  when its stored Cartesian coordinate fingerprint matches the IRC start
  (1.1e-3 bohr absolute tolerance to allow three-decimal PDB round-tripping).
  A missing or mismatched fingerprint is rejected and a fresh Hessian is
  calculated; endpoint Hessians are stored separately for endpoint RFO.
- `--tr-projection constrained` removes only full-system rigid motions that
  leave frozen anchors fixed. It does not select a reaction path; see
  `freeze-atoms.md`. An all-frozen structure raises an explicit error.

## See also

- `tsopt.md` — produces the IRC starting geometry.
- `freq.md`, `dft.md` — downstream.
- `bond-summary.md` — same bond-change algorithm, standalone.
- `pdb2reaction-workflows-output/SKILL.md` — R/TS/P path conventions.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.IRC_KW)`
