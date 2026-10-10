# `pdb2reaction irc`

## When to use

`irc` traces the intrinsic reaction coordinate from an optimized TS in both
directions with EulerPC, the only integrator, and writes the path and its two
raw end structures. Optimize those ends with `opt` to see which R and P the TS
connects. For optimized `reactant` / `product` files, run `pdb2reaction all`,
which calls `irc` and optimizes the ends for you.

## Minimal run

```bash
pdb2reaction irc -i result_tsopt/final_geometry.xyz -q 0 -m 1 -b uma --out-json -o result_irc
```

Success: `finished_irc_trj.xyz` shows the reaction, and the optimized ends are
the intended R and P (next section). `--out-json` writes `result.json`.

Smaller step and longer integration for a shallow surface:

```bash
pdb2reaction irc -i ts.xyz -q -1 -m 1 \
    --max-cycles 250 --step-size 0.05 \
    -b uma -o result_irc_long
```

If IRC stops after only a few frames, **reduce `--step-size` first** (usually
from `0.10` to `0.05`). Use `--never-stop` only when you want to trace to the
cycle cap unconditionally. It is opt-in and may pass the nearest basin, so
inspect the trajectory and endpoints.

## Judge success

Even if the IRC does not converge, the result is usable when the endpoints,
optimized with `opt`, reach the intended R and P.

1. Optimize both ends. They are `.xyz` only, so pass the TS PDB (for example
   `result_tsopt/final_geometry.pdb` when tsopt read a PDB) with
   `--ref-pdb` to keep the boundary frozen:

   ```bash
   pdb2reaction opt -i result_irc/finished_first.xyz --ref-pdb ts.pdb -q 0 -m 1 --out-json -o result_opt_first
   pdb2reaction opt -i result_irc/finished_last.xyz --ref-pdb ts.pdb -q 0 -m 1 --out-json -o result_opt_last
   ```

2. Both runs must end `converged`. Compare their bond patterns and coordinates
   with the intended R and P (for a TS from a path search, the MEP ends) to
   decide which is R and which is P. Do not infer chemical identity from
   first / last or from energy alone.

Files: `forward_irc_trj.xyz`, `backward_irc_trj.xyz`, the stitched
`finished_irc_trj.xyz` (first end → TS → last end), and the raw ends
`finished_first.xyz` and `finished_last.xyz`. For PDB/mmCIF input,
`--convert-files` (on by default) adds `.pdb` trajectories, and `.cif` with
the original IDs when needed. With a non-empty YAML `irc.prefix`, EulerPC
inserts one underscore before each filename (`prefix: trial` →
`trial_finished_irc_trj.xyz`); read the names from `files` in `result.json`.

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["forward_integration_stop_reason"], d["backward_integration_stop_reason"])
print(d["forward_integration_converged"], d["backward_integration_converged"])
print(d["energy_first_hartree"], d["energy_ts_hartree"], d["energy_last_hartree"])
print(d["never_stop"], d["never_stop_energy_bypasses"])
```

For `execution_status`, `completed` is not an IRC convergence verdict; it
means the runner returned.
IRC has no independent scientific success verdict, and `scientific_status` is
`success` whenever the integration ran without an error. Each requested
direction records its frame count and `*_integration_stop_reason`.
`*_integration_converged` records whether the RMS-gradient criterion fired, so
`--never-stop` leaves it false. It and `*_downhill_departure_valid` are
diagnostics, not gates for endpoint optimization: finite ends can go to `opt`
whatever the stop reason. Missing or non-finite coordinates and execution
errors must still be reported. `never_stop` records whether the opt-in mode was
on; `never_stop_energy_bypasses` counts the energy stops it bypassed.
Standalone IRC has no endpoint references. `endpoint_energy_orientation` and
`bond_changes_direction` are both `finished_first_to_finished_last`, not
chemical R → P.

`all` writes optimized, oriented `segments/seg_NN/reactant.*`, `ts.*`, and
`product.*`. `all --opt-mode-post grad` switches both its TS optimization (to
Dimer) and its endpoint optimization (to L-BFGS) from the default `hess`.

## Bond changes

`bond_changes` records the change from `finished_first.xyz` to
`finished_last.xyz` with a 1.20× covalent-radius cutoff and a 0.05 margin, the
same algorithm as `bond-summary` and `path-search` segmentation. After you
assign R and P, swap formed and broken if the chemical direction is last →
first. The key is absent when the endpoint comparison was not available.

```python
import json
bc = json.load(open("result_irc/result.json")).get("bond_changes")
if bc is None:
    raise RuntimeError("IRC endpoint comparison was not available")
for b in bc["formed"]: print("FORMED ", b)
for b in bc["broken"]: print("BROKEN ", b)
```

## Pitfalls and recovery

- **More than one imaginary mode at the start.** IRC needs a TS with one
  imaginary mode; otherwise it may follow the wrong one. Re-run `tsopt` first.
- **A branch reaches `--max-cycles`.** Inspect its gradient, energy trend, and
  geometry. A smaller `--step-size` stabilizes the integration but may need
  more cycles; no cycle count fits every surface.
- **Flickering bonds.** The bond-change check is geometric, so metal–ligand
  bonds near the cutoff may appear and disappear.
- **Hessian reuse in `all`.** IRC reuses the TS Hessian only when the geometry
  (within 1.1e-3 bohr, for three-decimal PDB round trips), calculator
  settings, and frozen atoms all match; otherwise it computes a new one.
  Endpoint Hessians are kept separately for endpoint optimization.
- **Frozen atoms.** Only rigid motions that keep every frozen atom in place are
  removed; this does not choose the path. See
  [PHVA treatment](extract.md#phva-treatment)
  and the [JSON record](../../docs/json-output.md#rigid-projection-provenance).

## Next step

- [tsopt.md](tsopt.md): makes the starting TS.
- [opt.md](opt.md), [freq.md](freq.md), [dft.md](dft.md): endpoints and downstream.
- [bond-summary](utilities.md#bond-summary): the same bond-change check on its own.
- [Reading outputs](../pdb2reaction-overview/outputs.md): R/TS/P paths in `all`.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.IRC_KW)`
