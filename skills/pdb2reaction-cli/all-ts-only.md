# `pdb2reaction all` — TS-only mode

## When to use

You already have a **TS candidate** (typically from another QM code, an
older `pdb2reaction` run, or a manual guess) and want to run only the
validation + thermochemistry stages — `tsopt → irc → freq → (dft)` —
without the upstream extract / path-search.

## Synopsis

```bash
pdb2reaction all -i ts_candidate.xyz \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    [--dft --dft-func-basis 'wb97m-v/def2-tzvpd'] \
    -o result_ts_only
```

Or with a PDB that carries residue / charge info:

```bash
pdb2reaction all -i ts_candidate.pdb \
    -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_ts_only
```

## How it differs from the other two modes

`pdb2reaction all` falls into TS-only mode when **all three** hold:

- exactly **one** `-i` input is given,
- **no** `--scan-lists` is provided,
- `--tsopt` (or `--tsopt True`) is passed.

Without `--scan-lists` or `--tsopt` the CLI raises `BadParameter`
(`all.md` covers the orchestrator's input gate).

For finer control, run the underlying subcommands directly:

```bash
pdb2reaction tsopt -i ts.xyz -q ... -m 1 -o result_tsopt -b uma
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -o result_irc -b uma
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -o result_freq -b uma
```

`extract` and `path-search` are skipped entirely; the chain collapses
to `tsopt → irc → freq → (dft)`. The output tree:

| Path | When | Content |
|---|---|---|
| `<out_dir>/summary.json` | always | machine-readable result |
| `<out_dir>/summary.log` | always | human-readable text + dir tree |
| `<out_dir>/segments/seg_01/{reactant,ts,product}.{pdb,xyz}` | always | canonical R/TS/P (TS-only mode; extension follows the `-i` input format) |
| `<out_dir>/segments/seg_01/structures/` | always | `.xyz` working copies of R/TS/P (plus `.pdb` when the input/ref is a PDB) |
| `<out_dir>/segments/seg_01/ts/final_geometry.{xyz,pdb}`, `optimization_trj.xyz` | always (`_trj.xyz` with `--dump`) | tsopt output |
| `<out_dir>/segments/seg_01/irc/{forward,backward,finished}_irc_trj.xyz` | always | IRC trajectories |
| `<out_dir>/segments/seg_01/freq/{R,TS,P}/{frequencies_cm-1.txt, thermoanalysis.yaml}` | always | per-state freq + thermo |
| `<out_dir>/segments/seg_01/dft/{R,TS,P}/result.yaml` | `--dft` | per-state DFT (`result.json` only when `dft` is run standalone with `--out-json`) |

## Output keys

```python
import json
d = json.load(open("result_ts_only/summary.json"))
seg  = d["segments"][0]            # MEP-style block (barrier/delta/bond_changes)
post = d["post_segments"][0]       # post-processing block (ts_imag / uma / dft …)
print(seg["barrier_kcal"], seg["delta_kcal"])
print(seg["bond_changes"])         # what bonds broke / formed along the IRC
print(post["ts_imag"]["n_imag"])   # should be 1 for a true TS
print(post["uma"]["energies_au"])  # [E(R), E(TS), E(P)] in hartree
```

If `post["ts_imag"]["n_imag"] != 1`, the geometry is **not a true first-order
saddle**; see "Distinctive failure modes" below.

## Distinctive failure modes

| Symptom | Likely cause | Fix |
|---|---|---|
| `tsopt.status == "not_converged"` | Initial Hessian misleading or step size too large | `pdb2reaction tsopt -i ts.xyz --opt-mode rsirfo --max-cycles 200` standalone, then re-run downstream stages |
| `post["ts_imag"]["n_imag"] == 0` | Geometry collapsed to a minimum during refinement | TS guess was not a real saddle; re-do `path-search` instead |
| `post["ts_imag"]["n_imag"] >= 2` | Two near-degenerate negative modes | Normal for some metalloenzyme TSs; check whether the second imaginary mode is a translation / rotation residue (often resolved by tightening `freeze_atoms`) |
| `irc.bond_changes == {}` (no bonds change) | TS connects two essentially identical wells (numerical ringing) | Verify the imaginary mode visualization in `freq/`; this is sometimes a non-physical TS |

## When *not* to use TS-only mode

- You do not yet have a TS candidate. Run `path-search` (or the
  full `all` in endpoint-MEP / scan-list mode) instead.
- You have a candidate but suspect the connectivity is wrong (i.e.
  you're not sure whether your "TS" sits between the right reactant
  and product). Use `path-search` to discover the connectivity.

## Caveats

- For an XYZ TS candidate, you must supply `-q` and `-m` explicitly
  (XYZ has no header). Use `--ref-pdb cluster.pdb` if you want
  `-l 'RES:Q'` to work.
- The IRC step here is the **canonical validation** that the TS
  connects the expected R and P. Always read `seg_01/{reactant,product}.pdb`
  to confirm the IRC ended up where you thought.

## See also

- `all.md` — base orientation.
- `tsopt.md`, `irc.md`, `freq.md`, `dft.md` — the underlying
  subcommands (which you can also run standalone if you want
  fine-grained control).
- `pdb2reaction-workflows-output/SKILL.md` — IRC interpretation
  and bond-change conventions.