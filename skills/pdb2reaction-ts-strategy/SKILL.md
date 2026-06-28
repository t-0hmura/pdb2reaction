---
name: pdb2reaction-ts-strategy
description: Decision know-how for pure-MLIP enzyme reaction-barrier runs with pdb2reaction — precision (fp32 vs fp64) per GPU class, the two TS-candidate routes (path-search MEP vs distance-restrained `scan`), fixing a bad imaginary-mode count (fp64 / `--coord-type dlc` / `--flatten`), reading a P-start scan barrier as the REVERSE direction, staged (`-s` repeated) vs concerted (one `-s`, many tuples) scans, and the same-atom-set rule for any mutant-vs-WT / mechanism-vs-mechanism comparison. TRIGGER when choosing precision, building a TS candidate, debugging imaginary modes, interpreting a barrier number, choosing staged vs concerted, or setting up a controlled barrier comparison. SKIP for install / HPC scheduler / output-parsing / structure-format-editing questions (use the dedicated skills). pdb2reaction is a PURE-MLIP cluster tool, so the comparison rule is enforced by feeding the SAME prepared atom set to every compared run.
---

# pdb2reaction TS strategy

All flags below verified against `pdb2reaction/cli/common_options.py`, `core/defaults.py`,
`workflows/{scan,tsopt,path_search,restraints}.py`. Verify version-sensitive values with
`pdb2reaction <sub> --help` / `--help-advanced` and `pdb2reaction.core.defaults`.

## 1. Precision: `--precision fp32|fp64`

`--precision` is backend-agnostic; the CLI routes it per backend. Effective default `fp32`
(`UMA_CALC_KW['precision']='fp32'`, option default `None` falls through).

| GPU class | Use | Why |
|---|---|---|
| HPC datacenter (H100 / H200 / A100) | `--precision fp64` | Deterministic, low numerical-noise TSopt / Hessian; fp64 throughput penalty is small on these cards |
| Consumer (RTX 50xx / 40xx) | `--precision fp32` (default) | fp64 is much slower here; fp32 is the speed / screening baseline |

| Backend | fp32 routes to | fp64 routes to |
|---|---|---|
| `uma` | `precision='fp32'` | `precision='fp64'` |
| `orb` | `precision='float32-high'` | `precision='float64'` |
| `mace` | `default_dtype='float32'` | `default_dtype='float64'` |
| `aimnet2` | no-op (already fp32) | REJECTED (inputs cast to float32 upstream) |

- fp64 has non-trivial TSopt / Hessian impact for OMol-trained UMA — use it for final / production numbers, not just screening.
- Combine with `--deterministic` (see `pdb2reaction-workflows-output` / docs `reproducibility.md`) for bit-identical reruns.

## 2. Two routes to a TS candidate

| Route | Subcommand | When | Mechanics |
|---|---|---|---|
| (a) MEP / path search | `path-search` | Have both endpoints (R and P); want the TS bracketed automatically | Recursive GSM/DMF MEP with bond-change detection; auto-segments a multi-step path, refines each reactive segment, returns `hei_seg_NN.xyz` (highest-energy image per segment) |
| (b) Distance-restrained build-up | `scan` | Have only the reactant (or want to drive a specific reacting distance) | Staged harmonic restraints `E=½k(r−target)²` (scan default `k=300` via `BIAS_KW`; the `10.0` in `HarmonicBiasCalculator` is only an unused constructor fallback) drive each reacting distance with full relaxation, walking up to a TS candidate |

- There is **no** `opt --restraint` flag. `opt` is plain unrestrained minimization; the restrained route is `scan`.
- `scan` supports `--preopt` / `--endopt` to relax the endpoints around the driven path.
- Feed either route's TS candidate into `tsopt → freq → irc` to confirm it (see `pdb2reaction-cli`).

## 3. Wrong imaginary-mode count after TS-opt

Symptom: a spurious 2nd small imaginary mode, or no dominant reaction mode (true first-order
saddle = exactly **one** imaginary frequency whose mode displaces along the reaction coordinate).
Detection cutoff `hessian_dimer.neg_freq_thresh_cm` (default 5 cm⁻¹).

| Lever | Flag | Effect |
|---|---|---|
| Raise precision | `--precision fp64` | Cleaner Hessian → removes numerical-noise imaginary modes |
| Internal coordinates | `--coord-type dlc` | Delocalised internal coords — slower but more robust convergence to a clean first-order saddle on torsion-rich systems |
| Flatten spurious modes | `--flatten` | Extra-imaginary-mode flattening loop (`grad`: dimer loop; `hess`: post-RS-I-RFO); `--no-flatten` forces `flatten_max_iter=0`. On `tsopt`, `opt`, `all` |

- `--coord-type` choices: `cart` (default, robust, published numbers) | `redund` | `dlc` | `tric`. On `path-opt` / `path-search` only `cart` / `dlc` are accepted.
- Try `--precision fp64` and/or `--coord-type dlc` first; add `--flatten` for residual small modes. They are complementary.
- Example: a mutant CM TS came out as the dominant Claisen mode −223 cm⁻¹ plus a residual −12.5 cm⁻¹; `--flatten` clears the spurious extra mode to a clean single-imaginary saddle.

## 4. P-start scan → barrier is the REVERSE direction

If a `scan` (or path) STARTS from the Product side, the raw reported barrier is the **reverse**
barrier. This is a read-time interpretation, not a flag.

| You have | Forward barrier |
|---|---|
| P-start run | `E(TS) − E(reactant)` — NOT the raw P-start number (that is `E(TS) − E(product)`, the reverse barrier) |

- Always confirm which endpoint the scan started from before quoting a barrier.

## 5. Staged vs concerted scan

`-s/--scan-lists` is `multiple=True`. One flag = one stage; tuples inside it run together.

| Mode | Syntax | Use |
|---|---|---|
| Concerted | one `-s` with several coordinate tuples | Coordinates move together in one stage; no mechanism breakdown needed |
| Staged | `-s` repeated (each `-s` = one sequential stage, `stage_NN/`) | Mechanism known up front; cleaner per-step control |

- When the mechanism IS known, **staged is generally preferred** (per-step barriers, per-stage output dirs).
- When the mechanism is unknown / multi-step, let `path-search` auto-segment instead of guessing stages.
- A 4-tuple `(i,j,low,high)` expands to a bidirectional 2-stage scan.

```bash
# concerted (one stage, two coords move together)
pdb2reaction scan -i R.pdb -s '[("Ca RES 10","Cb RES 11",1.6),("H RES 11","O GLU 20",1.0)]' -o out
# staged (stage_01 then stage_02)
pdb2reaction scan -i R.pdb \
  -s '[("Ca RES 10","Cb RES 11",1.6)]' \
  -s '[("H RES 11","O GLU 20",1.0)]' -o out
```

## 6. Controlled comparison: SAME atom set across all models

For a mutant-vs-WT (or mechanism-vs-mechanism) barrier comparison, **every compared model must
use the identical atom set** — same atom count and same residues. Otherwise it is not a controlled
experiment and the barrier difference is not interpretable.

| Aspect | pdb2reaction (pure-MLIP) |
|---|---|
| How the rule is enforced | Prepare ONE cluster atom set, then mutate / change mechanism on THAT same set so atom count + residues stay identical across all compared runs |
| Non-standard ligand charge | Keep `-l 'RES:Q'` (or `-l 'GPP:-3,SAM:1'`) consistent across compared runs so charge differences don't confound the comparison |

- Do NOT re-extract each variant independently (different `--radius` / residue inclusion silently changes the atom set → mismatch, uncontrolled comparison). Build the shared cluster once and edit in place.

## See also

- `pdb2reaction-cli/SKILL.md` — per-subcommand flags (`scan` / `tsopt` / `path-search` / `freq` / `irc`).
- `pdb2reaction-workflows-output/SKILL.md` — six canonical workflows, R/TS/P paths, `--deterministic`.
- docs `tsopt.md` (TS routes / imaginary-mode fix / controlled comparison), `scan.md` (staged vs concerted, barrier direction), `reproducibility.md` (precision by GPU class) — the user-facing form of this know-how, folded from the former standalone page.
