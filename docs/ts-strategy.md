# Transition-state strategy and barrier interpretation

A decision guide for getting a clean transition state (TS) and a correctly-read
reaction barrier out of a `pdb2reaction` (pure-MLIP) enzyme run. It covers six
choices you make in almost every barrier study: precision, how to build a TS
candidate, what to do when the imaginary-mode count is wrong, how to read a
barrier when the scan started from the product, staged vs concerted scans, and
how to keep a mutant-vs-WT comparison controlled.

For per-subcommand flags and examples, keep [`tsopt`](tsopt.md),
[`scan`](scan.md), and [`path-search`](path-search.md) open in parallel.

## 1. Precision: pick `fp64` on datacenter GPUs, `fp32` elsewhere

`--precision` selects the floating-point precision of MLIP inference. It is
backend-agnostic — the CLI routes the value into each backend's own key — and
defaults to `fp32` (the established, screening-speed baseline).

| GPU class | Recommended | Why |
| --- | --- | --- |
| HPC datacenter (H100 / H200 / A100) | `--precision fp64` | Deterministic, low numerical-noise TS optimization and Hessians; the fp64 throughput cost is small on these cards |
| Consumer (RTX 50xx / 40xx) | `--precision fp32` (default) | fp64 is substantially slower here; fp32 is the speed / screening baseline |

How the value is routed per backend:

| Backend | `fp32` | `fp64` |
| --- | --- | --- |
| `uma` | `precision='fp32'` | `precision='fp64'` |
| `orb` | `precision='float32-high'` | `precision='float64'` |
| `mace` | `default_dtype='float32'` | `default_dtype='float64'` |
| `aimnet2` | no-op (already fp32) | rejected — its model inputs are cast to float32 upstream |

```bash
# datacenter H200: production TS + Hessian
pdb2reaction tsopt -i ts_candidate.xyz -q -1 -m 1 -b uma --precision fp64 -o result_tsopt
```

fp64 has a non-trivial effect on TS optimization and Hessians for the
OMol-trained UMA backend, so use it for final / production numbers — not only
for screening. For bit-identical reruns, combine it with `--deterministic`
(see [Reproducibility](reproducibility.md)).

## 2. Two routes to a TS candidate

There are two complementary ways to obtain a TS candidate before you run
`tsopt → freq → irc` to confirm it.

| Route | Subcommand | Use when | What it does |
| --- | --- | --- | --- |
| (a) MEP / path search | [`path-search`](path-search.md) | You have both endpoints (reactant **and** product) and want the TS bracketed automatically | Recursive minimum-energy-path search (GSM / DMF) with bond-change detection; it auto-segments a multi-step path, refines each reactive segment, and returns the highest-energy image per segment (`hei_seg_NN.xyz`) |
| (b) Distance-restrained scan | [`scan`](scan.md) | You have only the reactant, or want to drive a specific reacting distance directly | Staged harmonic distance restraints, `E = ½k(r − target)²`, drive each reacting distance with full relaxation, walking the system up to a TS candidate |

There is no separate `opt --restraint` flag: `opt` is plain unrestrained
minimization, and the distance-restrained build-up route is `scan`. `scan` can
relax the endpoints around the driven path with `--preopt` / `--endopt`.

```bash
# route (a): bracket the TS between two endpoints
pdb2reaction path-search -i reactant.pdb product.pdb -q 0 -m 1 -o result_mep

# route (b): drive the reacting distance from the reactant only
pdb2reaction scan -i reactant.pdb -s '[("Ca RES 10","Cb RES 11",1.6)]' -q 0 -m 1 -o result_scan
```

## 3. Wrong imaginary-mode count after TS optimization

A true first-order saddle has **exactly one** imaginary frequency, and its mode
displaces along the reaction coordinate (detection cutoff
`hessian_dimer.neg_freq_thresh_cm`, default 5 cm⁻¹). If TS-opt instead gives a
spurious second small imaginary mode, or no dominant reaction mode, escalate the
following levers — they are complementary, so you can combine them.

| Lever | Flag | Effect |
| --- | --- | --- |
| Raise precision | `--precision fp64` | A cleaner Hessian removes numerical-noise imaginary modes |
| Internal coordinates | `--coord-type dlc` | Delocalised internal coordinates — slower, but more robust convergence to a clean first-order saddle on torsion-rich systems |
| Flatten small modes | `--flatten` | Displaces along and re-relaxes residual small imaginary modes (available on `tsopt`, `opt`, and `all`) |

`--coord-type` accepts `cart` (the robust default behind the published numbers),
`redund`, `dlc`, and `tric`; on `path-opt` / `path-search` only `cart` and `dlc`
are accepted. Try `--precision fp64` and/or `--coord-type dlc` first, then add
`--flatten` to clean up any residual small modes.

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q -1 -m 1 \
    --precision fp64 --coord-type dlc --flatten -o result_tsopt
```

## 4. A product-start scan reports the *reverse* barrier

If your `scan` (or path) starts from the **product** side, the raw barrier it
reports is the **reverse** barrier — `E(TS) − E(product)`. To quote the forward
barrier, compute it from the reactant:

| You ran | Forward barrier |
| --- | --- |
| A product-start scan | `E(TS) − E(reactant)` — **not** the raw product-start number |

This is a read-time interpretation, not a flag. Always confirm which endpoint
the scan started from before quoting a barrier, especially when the workflow was
seeded from a crystallographic product complex.

## 5. Staged vs concerted scans

`-s/--scan-lists` can be passed more than once. A single `-s` flag defines one
stage; the coordinate tuples inside it move together. Repeating `-s` defines
sequential stages, each written to its own `stage_NN/` directory.

| Mode | Syntax | Use when |
| --- | --- | --- |
| Concerted | one `-s` with several coordinate tuples | The coordinates move together in a single step; you do not need to break the mechanism into stages |
| Staged | `-s` repeated (one per sequential stage) | The mechanism is known up front and you want clean per-step control and per-stage output |

When the mechanism **is** known, the staged form is generally preferred — it
gives per-step barriers and per-stage geometries. When the mechanism is unknown
or multi-step, let [`path-search`](path-search.md) auto-segment the path instead
of guessing the stages yourself. A 4-tuple `(i, j, low, high)` expands into a
bidirectional 2-stage scan.

```bash
# concerted: two coordinates move together in one stage
pdb2reaction scan -i reactant.pdb \
    -s '[("Ca RES 10","Cb RES 11",1.6),("H RES 11","O GLU 20",1.0)]' -o result_concerted

# staged: stage_01 then stage_02
pdb2reaction scan -i reactant.pdb \
    -s '[("Ca RES 10","Cb RES 11",1.6)]' \
    -s '[("H RES 11","O GLU 20",1.0)]' -o result_staged
```

## 6. Controlled comparisons need the same atom set

For a mutant-vs-WT (or mechanism-vs-mechanism) barrier comparison, **every
compared model must use the identical atom set** — the same atom count and the
same residues. Otherwise the comparison is not controlled and the barrier
difference is not interpretable.

`pdb2reaction` is a pure-MLIP cluster tool: there is no ML/MM layer concept, no
layer-detection flag, and no geometric layer split. The same-atom-set rule is
therefore enforced by construction:

- Prepare **one** cluster atom set, then apply the mutation (or change the
  mechanism) **on that same set**, so the atom count and residues stay identical
  across every compared run. Edit the shared cluster in place — do **not**
  re-extract each variant independently, because a different `--radius` or
  residue inclusion silently changes the atom set and breaks the comparison.
- Keep the non-standard ligand charge consistent across all compared runs with
  `-l 'RES:Q'` (e.g. `-l 'GPP:-3,SAM:1'`), so a charge difference never
  confounds the barrier comparison.

```bash
# WT and mutant share one prepared cluster (identical atom count + residues)
pdb2reaction all -i wt_cluster.pdb     -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_wt
pdb2reaction all -i mutant_cluster.pdb -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_mutant
```

## See also

- [`tsopt`](tsopt.md) — TS optimization flags (`--coord-type`, `--flatten`, optimizer modes).
- [`scan`](scan.md) — staged distance-restrained scans.
- [`path-search`](path-search.md) — recursive MEP search and bond-change segmentation.
- [Reproducibility](reproducibility.md) — `--precision` together with `--deterministic`.
- [Quickstart: TS-only](quickstart-tsopt-freq.md) — confirm a TS candidate end-to-end.
