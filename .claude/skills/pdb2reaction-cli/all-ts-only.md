# `pdb2reaction all` — TS-only mode

## When to use

You already have a **TS candidate** (typically from another QM code, an
older `pdb2reaction` run, or a manual guess) and want to run only the
validation + thermochemistry stages — `tsopt → irc → freq → (dft)` —
without the upstream extract / path-search.

This is the "I trust this geometry, just check it for me" mode.

## Synopsis

```bash
pdb2reaction all -i ts_candidate.xyz \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    [--dft --dft-func-basis 'wb97m-v/def2-svp'] \
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

If a single `-i` is given without either `--scan-lists` or `--tsopt`,
the CLI raises `BadParameter` instead of running (see `all.py`:
"Provide at least two structures... or a single structure with
--scan-lists, or a single structure with --tsopt True"). With the
required `--tsopt` set, the orchestrator skips path-search and starts
the pipeline at `tsopt`.

For finer control, run the underlying subcommands directly:

```bash
pdb2reaction tsopt -i ts.xyz -q ... -m 1 -o result_tsopt -b uma
pdb2reaction irc   -i result_tsopt/final_geometry.xyz -o result_irc -b uma
pdb2reaction freq  -i result_tsopt/final_geometry.xyz -o result_freq -b uma
```

## Pipeline collapses to

```
ts_candidate.{xyz,pdb,gjf}
       │
       ▼
   [tsopt]            (RS-I-RFO default; Dimer alternative)
       │
       ▼
   [irc]              (forward + backward; LBFGS endpoint refinement)
       │
       ▼
   [freq]             (Hessian + thermo)
       │
       ▼
   [dft]              (optional)
```

`extract` and `path-search` are skipped entirely. The output tree
collapses to one segment:

```
result_ts_only/
├── summary.json
├── summary.log
├── seg_01/                # canonical R/TS/P (same layout as multi-input mode)
│   ├── reactant.{pdb,xyz}
│   ├── ts.{pdb,xyz}
│   └── product.{pdb,xyz}
└── tsopt_single/
    ├── ts/            final_geometry.{xyz,pdb}, optimization_trj.xyz (with --dump)
    ├── irc/           forward_irc_trj.xyz, backward_irc_trj.xyz, finished_irc_trj.xyz
    ├── freq/          frequencies_cm-1.txt, thermoanalysis.yaml (with --dump)
    ├── dft/           (when --dft)
    ├── endpoint_opt/  intermediate optimizer scratch
    └── structures/    reactant.pdb, ts.pdb, product.pdb (mirror of seg_01/)
```

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
| `tsopt.n_imaginary_modes == 0` | Geometry collapsed to a minimum during refinement | TS guess was not a real saddle; re-do `path-search` instead |
| `tsopt.n_imaginary_modes == 2+` | Two near-degenerate negative modes | Normal for some metalloenzyme TSs; check whether the second imaginary mode is a translation / rotation residue (often resolved by tightening `freeze_atoms`) |
| `irc.bond_changes == {}` (no bonds change) | TS connects two essentially identical wells (numerical ringing) | Verify the imaginary mode visualization in `freq/`; this is sometimes a non-physical TS |

## When *not* to use TS-only mode

- You do not yet have a TS candidate. Run `path-search` (or the
  full `all` in endpoint-MEP / scan-list mode) instead.
- You have a candidate but suspect the connectivity is wrong (i.e.
  you're not sure whether your "TS" sits between the right reactant
  and product). Use `path-search` to discover the connectivity.

## Caveats

- Passing `--tsopt false` skips TS optimization entirely, which is
  rarely what you want in this mode.
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