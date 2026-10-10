# path-opt and path-search

## When to use

- `path-opt` optimizes one MEP between two endpoints with GSM (default) or
  DMF. It does not split the path. Use it when one endpoint pair is the
  intended candidate unit.
- `path-search` takes two or more endpoints in reaction order, finds bond
  changes along the path, and splits it recursively into candidate steps, each
  with its own HEI. `all --refine-path` runs it; by default `all` runs one
  `path-opt` pass. It does not run `tsopt`, `irc`, `freq`, or `dft`; `all`
  does that under `segments/seg_NN/`.

Both give a candidate path and HEI; `tsopt`, `freq`, and `irc` decide whether
a segment is one elementary step.

## Minimal run

All inputs follow one `-i`.

### path-opt

```bash
pdb2reaction path-opt -i R.xyz P.xyz -q 0 -m 1 -b uma --out-json -o result_path_opt

# DMF engine
pdb2reaction path-opt -i R.pdb P.pdb -l 'GPP:-3' --mep-mode dmf -b mace \
    -o result_path_opt_dmf
```

### path-search

```bash
pdb2reaction path-search -i 1.R.pdb 3.P.pdb -l 'SAM:1,GPP:-3' -b uma \
    -o result_path_search

# With an intermediate you supply
pdb2reaction path-search -i 1.R.pdb 2.IM.pdb 3.P.pdb -l 'SAM:1,GPP:-3' \
    -b uma --max-nodes 30 -o result_path_search_im

# DMF engine
pdb2reaction path-search -i 1.R.pdb 3.P.pdb --mep-mode dmf --refine-mode minima \
    -l 'SAM:1,GPP:-3' -b uma -o result_path_search_dmf
```

## Judge success

### path-opt

The console prints `[write] Wrote '…/hei.xyz'.` when the HEI, the TS
candidate, is written. In `result.json` (`--out-json`), `scientific_status` is
`success` when every requested stage (endpoint pre-optimization, MEP)
converged, otherwise `partial` or `failed`. `optimization_status` is
`converged` or `not_converged` when the engine reports convergence, and
`completed` when it reports nothing; do not read `completed` as converged.
`barrier_kcal` is the HEI energy relative to the first image.

`hei_index` between 1 and `n_images − 2` gives a TS candidate. At 0 or
`n_images − 1` the HEI is an endpoint and not a TS candidate: check the
endpoints, or get a candidate another way
([ts-strategy](../pdb2reaction-overview/ts-strategy.md#when-the-ts-does-not-come-out)).

```text
result_path_opt/
├─ final_geometries_trj.xyz   # final path, energies on the comment lines
├─ final_geometries.pdb       # PDB/mmCIF input; DMF names it final_geometries_trj.pdb
├─ hei.xyz, hei.pdb           # HEI, the TS candidate
├─ dmf_initial_trj.xyz        # DMF only, with dmf_ipopt.out
├─ align_refine/              # endpoint alignment and relaxation
└─ result.json, summary.json  # --out-json
```

Give `hei.pdb` to `tsopt` so that `-l` and `--freeze-links` apply; with
`hei.xyz`, add `--ref-pdb`.

### path-search

`summary.json` and `summary.log` are written once the path is built, fully
or in part; open the
`[2] Segment-level MEP summary` section of `summary.log` or read
`summary.json`. `scientific_status` is `success` when the pre-optimizations
and every path run converged, otherwise `partial` or `failed`. `segments`
lists each segment in path order:

```python
{
  "index": 1, "tag": "seg_01", "kind": "seg", "converged": True,
  "barrier_kcal": 21.5, "delta_kcal": -0.7,
  "bond_changes": [
    {"Bond formed (1)": ["C508-C567 : 3.166 Å --> 1.675 Å"]},
    {"Bond broken (1)": ["S507-C508 : 1.798 Å --> 3.459 Å"]}
  ]
}
```

- `kind` is `seg` (reactive), `kink` (conformation only), or `bridge` (short
  connecting path). A `seg` entry with bond changes, and its
  `hei_seg_NN.xyz`, is a TS candidate.
- A `tag` ending in `_maxdepth` means splitting stopped at the depth limit or
  after repeated kinks; the segment may hold more than one step. Raise
  `--max-depth` (default 10) or give intermediates.
- Only `kink` segments mean that no bond change was found between the ends;
  check the inputs or add intermediates.
- The warning `HEI is at an endpoint` means that in that interval no image lies
  above the higher end, so, as for `path-opt` above, there is no TS candidate:
  `path-search` returns the interval as a raw path with no `segments` entry or
  `hei_seg_NN.xyz`, adds `endpoint_hei` to `scientific_status_reasons`, and
  `all --refine-path` runs no TS optimization on it. It is not a zero barrier
  and does not show that no bond changes. Read the energies on the comment
  lines of that run's `seg_NNN_*/final_geometries_trj.xyz` and its bond changes
  ([`bond-summary`](utilities.md#bond-summary)); a frame at an interior local
  maximum can start `tsopt` like an HEI. Other routes:
  [ts-strategy](../pdb2reaction-overview/ts-strategy.md#when-the-ts-does-not-come-out).

A successful TS optimization gives one imaginary mode along the reaction
coordinate; confirm each HEI with `tsopt` (n_imag = 1) and IRC before reading
it as a step of the mechanism.

```text
result_path_search/
├─ mep_trj.xyz                       # whole stitched MEP
├─ mep_trj.pdb                       # PDB/mmCIF input or --ref-pdb
├─ mep_plot.png, energy_diagram_MEP.png  # the PNG diagram needs Kaleido
├─ summary.json, summary.log
├─ mep_seg_NN_trj.xyz, hei_seg_NN.xyz  # reactive segment NN and its HEI
├─ mep_w_ref*.pdb                    # --write-ref-merge with --ref-full-pdb
└─ seg_NNN_*/                        # working files of each GSM/DMF run
```

NN is the segment `index` (from 01); NNN counts the GSM/DMF runs from 000, so
the two numbers differ. Working-directory names carry the tags `_mep`,
`_maxdepth`, and `_bridge`.

## Pitfalls and recovery

- All inputs need the same atoms, elements, and order; equal atom counts are
  not enough.
- Convergence depends on the endpoints. Endpoint pre-optimization is on by
  default; on `not_converged`, inspect its result and the string before
  choosing a separate `opt`, a different endpoint, or other path settings.
  When it converges but changes the bonding, see
  [Run stage by stage](../pdb2reaction-overview/SKILL.md#run-stage-by-stage-and-judge-each-stage).
- With frozen atoms, `path-opt`, `path-search`, and `all` align each input to
  the previous one and then relax it with L-BFGS (up to 10000 cycles), even
  when the frozen atoms already coincide and whatever `--preopt/--no-preopt`
  says. If that relaxation does not converge, the run stops with
  `Input alignment did not converge for pair(s)`. Only `path-search` can skip
  the step (`--align/--no-align`); use `--no-align` for endpoints cut together
  that share identical frozen coordinates.
- If GSM stalls, DMF (`--mep-mode dmf`) is an alternative. It needs
  `cyipopt` from conda-forge; after a GPU out-of-memory error with the default
  `--dmf-backend gpu`, retry with `--dmf-backend cpu`. If a DMF run makes
  almost no progress, set `BLIS_NUM_THREADS=1` in the job script before Python
  starts: nested BLIS threads in IPOPT can stall it, and a change inside the
  running process does not reach the solver. DMF holds frozen atoms with a
  stiff restraint, so they can move slightly.
- Another engine or a larger `--max-nodes` costs more and does not repair an
  inconsistent endpoint pair or a wrong mechanism.
- `path-search` can return more segments than inputs minus one by proposing
  intermediates you did not give; validate each one. A segment is not
  guaranteed to be one elementary step or to contain exactly one TS.
- `path-search` writes no optimized TS; `all --tsopt` writes it as
  `segments/seg_NN/ts.pdb`.

## Next step

- Optimize `hei.pdb` or `hei_seg_NN.xyz`: [tsopt.md](tsopt.md).
- Relax the endpoints first: [opt.md](opt.md).
- Bond changes between two structures:
  [`bond-summary`](utilities.md#bond-summary).
- The same path inside the full workflow: [all-endpoint-mep.md](all-endpoint-mep.md).
- Flags and defaults: `--help-advanced` and [SKILL.md](SKILL.md#where-flags-and-defaults-live).
