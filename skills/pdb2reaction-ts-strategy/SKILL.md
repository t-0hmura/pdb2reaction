---
name: pdb2reaction-ts-strategy
description: "Decision know-how for pure-MLIP enzyme reaction-barrier runs with pdb2reaction — backend-aware fp32/fp64 choices, the two TS-candidate routes (path-search MEP vs distance-restrained scan), fixing a bad imaginary-mode count, reading a P-start scan barrier as the reverse direction, staged vs concerted scans using one -s occurrence, and controlled cross-variant comparisons without confusing them with within-path atom-order requirements. TRIGGER when choosing precision, building a TS candidate, debugging imaginary modes, interpreting a barrier, choosing staged vs concerted scans, or setting up a mutant-vs-WT or mechanism comparison. SKIP for install, HPC scheduler, output parsing, or structure-format editing questions."
---

# pdb2reaction TS strategy

All flags below verified against `pdb2reaction/cli/common_options.py`, `core/defaults.py`,
`workflows/{scan,tsopt,path_search,restraints}.py`. Verify version-sensitive values with
`pdb2reaction <sub> --help` / `--help-advanced` and `pdb2reaction.core.defaults`.

## 1. Precision: `--precision fp32|fp64`

`--precision` is backend-agnostic; the CLI routes it per backend. Unset = the `auto` sentinel
(`CALC_KW_DEFAULT['precision'] = 'auto'`), resolved **per backend** by
`backends._BACKEND_DEFAULT_PRECISION`: `uma` → fp32, `orb` → fp64, `mace` → fp64, `aimnet2` → fp32.
So ORB and MACE run fp64 Hessians by default; UMA and AIMNet2 default to fp32.

Leave precision unset (`auto`) to preserve the backend default. When precision
matters, compare supported settings on the target model, system, and hardware.
AIMNet2 accepts fp32 only. Every setting still requires independent frequency
and IRC validation; precision is not proof of a saddle.

| Backend | fp32 routes to | fp64 routes to |
|---|---|---|
| `uma` | `precision='fp32'` | `precision='fp64'` |
| `orb` | `precision='float32-high'` | `precision='float64'` |
| `mace` | `default_dtype='float32'` | `default_dtype='float64'` |
| `aimnet2` | no-op (already fp32) | REJECTED (inputs cast to float32 upstream) |

- fp64 can materially change TSopt/Hessian behavior for UMA, but its speed cost
  is hardware/model dependent; validate the chosen production setting.
- `--deterministic` requests same-stack repeatability; it does not increase
  numerical precision, remove physical imaginary modes, or guarantee identical
  output across PyTorch/backend versions, hardware, or custom calculators.

## 2. Two routes to a TS candidate

| Route | Subcommand | When | Mechanics |
|---|---|---|---|
| (a) MEP / path search | `path-search` | Have both endpoints (R and P); want the TS bracketed automatically | Recursive GSM/DMF MEP with bond-change detection; auto-segments a multi-step path, refines each reactive segment, returns `hei_seg_NN.xyz` (highest-energy image per segment) |
| (b) Distance-restrained build-up | `scan` | Have only the reactant (or want to drive a specific reacting distance) | Staged harmonic restraints `E=½k(r−target)²` (scan default `k=300` via `BIAS_KW`; the `10.0` in `HarmonicBiasCalculator` is only an unused constructor fallback) optimize the remaining degrees of freedom while driving each distance; frames are biased geometries, not stationary points on the bare PES |

- There is **no** `opt --restraint` flag, but `opt` **does** support harmonic distance restraints via `--dist-freeze` (with `--bias-k`); `scan` is the route for *driving/walking* a reacting coordinate up to a TS candidate.
- `scan` supports `--preopt` (unbiased optimization of the **initial structure** before the scan) and `--endopt` (unbiased optimization of **each stage's result**, run after that stage); both default off.
- Feed either route's TS candidate into `tsopt → freq → irc` to confirm it (see `pdb2reaction-cli`).

## 3. Wrong imaginary-mode count after TS-opt

Symptom: a second imaginary mode, or no dominant reaction mode (a certified
first-order saddle has exactly **one** imaginary frequency whose mode displaces
along the reaction coordinate).

| Lever | Flag | Effect |
|---|---|---|
| Compare precision where supported | `--precision fp32|fp64` | Numerical behavior is backend/model/system dependent; neither setting removes a genuine second negative-curvature direction |
| Coordinate representation | `--coord-type dlc` | Delocalized internal coordinates can change conditioning relative to Cartesian coordinates; benchmark both on the problematic seed because neither is uniformly faster or more reliable |
| Flatten spurious modes | `--flatten` | Extra-imaginary-mode flattening loop (`grad`: dimer loop; `hess`: post-RS-P-RFO); `--no-flatten` forces `flatten_max_iter=0`. On `tsopt`, `opt`, `all` |

- `--coord-type` choices: `cart` (default) | `redund` | `dlc` | `tric`. On `path-opt` / `path-search` only `cart` / `dlc` are accepted.
- Inspect all mode displacements, the MEP seed, and optimizer stop reason before retrying an appropriate precision/coordinate/flattening setting. AIMNet2 rejects fp64. In every case, independently recompute frequencies and verify IRC connectivity.
- If a path-derived HEI is simply poor, rerun the parent `all` command with
  `--refine-path` so recursive `path-search` resolves the MEP before
  TSOPT. This is deliberately off by default: a bad/noisy path can be split
  into unnecessary segments and multiply MEP, TSOPT, IRC, and frequency cost.
- `--ref-mode` is not a normal standalone TSOPT remedy. It is the
  advanced/internal Cartesian path direction that `all` derives from its MEP
  and passes automatically for mode tracking and bounded saddle recovery.
  Only external-path expert workflows should provide it manually.

## 4. P-start scan → barrier is the REVERSE direction

If a `scan` (or path) STARTS from the Product side, the raw reported barrier is the **reverse**
barrier. This is a read-time interpretation, not a flag.

| You have | Forward barrier |
|---|---|
| P-start run | `E(TS) − E(reactant)` — NOT the raw P-start number (that is `E(TS) − E(product)`, the reverse barrier) |

- Always confirm which endpoint the scan started from before quoting a barrier.

## 5. Staged vs concerted scan

Pass a **single** `-s/--scan-lists` followed by one or more space-separated Python literals: each
literal is one stage, and the tuples inside a literal move together within that stage. Repeating the
flag is rejected (`Use a single --scan-lists followed by multiple values; repeated flags are not
accepted.`).

| Mode | Syntax | Use |
|---|---|---|
| Concerted | one `-s`, one literal holding several coordinate tuples | Tests a coupled-coordinate hypothesis by restraining the coordinates in one stage |
| Staged | one `-s`, several literals separated by spaces (each literal = one sequential stage, `stage_NN/`) | Tests an explicitly ordered sequence, with a separate restrained optimization/output directory per stage |

- Choose staged versus concerted from the mechanistic hypothesis; staging a
  truly concerted event can create an artificial intermediate, while coupling
  unrelated coordinates can hide a stepwise route.
- When the mechanism is unknown and both endpoints are available,
  `path-search` can propose bond-change segments. Those are candidates that
  still require TS/frequency/IRC validation, not an automatic mechanism proof.
- A 4-tuple `(i,j,low,high)` expands to a bidirectional 2-stage scan.

```bash
# concerted (one stage, two coords move together)
pdb2reaction scan -i R.pdb -q 0 -m 1 -s '[("Ca RES 10","Cb RES 11",1.6),("H RES 11","O GLU 20",1.0)]' -o out
# staged (stage_01 then stage_02) — ONE -s, two space-separated literals
pdb2reaction scan -i R.pdb -q 0 -m 1 \
  -s '[("Ca RES 10","Cb RES 11",1.6)]' \
     '[("H RES 11","O GLU 20",1.0)]' -o out
```

## 6. Controlled comparisons: distinguish two atom-set rules

Within **one** R→IM→P path, every structure must have exactly the same atoms,
elements, and ordering; this is an input contract of the MEP machinery. Across
WT and mutant systems, however, a real mutation can change residue identity and
atom count. Their raw total energies must not be subtracted from one another,
but their separately computed activation barriers can be compared:

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`.

Keep the modeling protocol controlled: use the same residue-position selection
and cluster-boundary policy except for the intended mutation, preserve the same
protonation/charge convention and cap placement rules, and use the same
backend/model/precision/thermochemistry settings. Independent automatic
extractions can choose different boundary residues near a radius cutoff, so
compare and harmonize the selected residue positions rather than demanding an
impossible identical atom count after every mutation. For mechanism comparisons
on the **same chemical composition**, keep one common atom set and order.

## See also

- `pdb2reaction-cli/SKILL.md` — per-subcommand flags (`scan` / `tsopt` / `path-search` / `freq` / `irc`).
- `pdb2reaction-workflows-output/SKILL.md` — six canonical workflows, R/TS/P paths, `--deterministic`.
- docs `tsopt.md` (TS routes / imaginary-mode fix / controlled comparison), `scan.md` (staged vs concerted, barrier direction), `reproducibility.md` (precision by GPU class).
