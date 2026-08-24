# `pdb2reaction tsopt`

## Purpose

Transition-state optimization. The default algorithm is RS-P-RFO
(full-Hessian; `--opt-mode hess`/`rsprfo`). Hessian-Guided Dimer
(`--opt-mode grad`/`dimer`) is the **alternative** TS optimizer when
RS-P-RFO fails to converge (e.g. unstable full-Hessian eigenstructure or very
large clusters where full-Hessian recomputation is prohibitive). The
dimer still uses an initial Hessian to set the search direction, then
updates the lowest mode via dimer rotation rather than recomputing the
full Hessian each cycle. Use after `path-search` or `scan` to refine a
HEI to a true first-order saddle, or as a standalone validator on an
externally-generated TS guess.

## Synopsis

```bash
pdb2reaction tsopt -i ts_guess.{pdb,cif,mmcif,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--opt-mode grad|hess|dimer|rsirfo|trim|rsprfo] \
    [--max-cycles 100000] \
    [-b uma|orb|mace|aimnet2] [-o ./result_tsopt/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | TS candidate; `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `hess` | `grad`/`dimer` (Hessian-Guided Dimer), `hess`/`rsprfo` (RS-P-RFO), `rsirfo` (RS-I-RFO), or `trim` (TRIM/Helgaker) |
| `--max-cycles` | int | 100000 | Optimization step cap |
| `--hessian-calc-mode` | str | (live default) | `Analytical` or `FiniteDifference` (default: `FiniteDifference`); selects how the initial Hessian is computed |
| `--workers`, `--workers-per-node` | int | `1`, `1` | UMA predictor workers. `workers > 1` with explicit `Analytical` raises `BackendError`; it does not fall back. Other built-in backends ignore these worker kwargs. |
| `--ref-mode` | path | none | Advanced/internal Cartesian 3N reaction direction used for initial-root selection and overlap tracking. `all` supplies it by default; with `all --no-tsopt-from-mep-tan`, TSOPT selects from the initial-structure Hessian modes. Ordinary standalone `tsopt` runs should omit it. Visible only in `--help-advanced`. |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_tsopt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

`tsopt` always forces `reject_uphill=False`, regardless of optimizer mode or
YAML. Uphill trial steps can be part of saddle-point mode following. The
`--reject-uphill/--no-reject-uphill` toggle belongs only to minimum
optimization (`opt`) and post-IRC endpoint refinement (`all`).

## Examples

### Default RS-P-RFO

```bash
pdb2reaction tsopt -i hei.xyz -q 0 -m 1 -b uma --out-json -o result_tsopt
```

### Dimer mode (alternative when RS-P-RFO fails to converge)

```bash
pdb2reaction tsopt -i hei.xyz -q 0 -m 1 \
    --opt-mode dimer -b uma -o result_tsopt_dimer
```

### Alternative RFO family on a difficult seed

```bash
pdb2reaction tsopt -i hei.pdb -l 'SAM:1,GPP:-3' \
    --opt-mode rsirfo --max-cycles 200 -b mace \
    -o result_tsopt_rsirfo
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/final_geometry.xyz` | completed optimizer run | final geometry; it is a validated first-order saddle only when `result.json["status"] == "converged"` and the downstream IRC is correct |
| `<out_dir>/final_geometry.pdb` | `--convert-files` (default on) and PDB/mmCIF topology/reference available | normalized PDB companion used between pipeline stages |
| `<out_dir>/final_geometry.cif` | `--convert-files` and input/reference required the mmCIF or oversized-PDB bridge | public companion with original IDs |
| `<out_dir>/final_geometry.gjf` | `--convert-files` and input is `.gjf` | Gaussian companion |
| `<out_dir>/optimization_trj.xyz`, `optimization.{pdb,cif}` | `--dump`, Hessian-family mode; companions require `--convert-files`, topology, and CIF bridge metadata | full optimization trajectory |
| `<out_dir>/optimization_all_trj.xyz`, `optimization_all.{pdb,cif}` | `--dump`, Dimer-family mode; companions require `--convert-files`, topology, and CIF bridge metadata | full optimization trajectory |
| `<out_dir>/vib/imag_<freq>cm-1_trj.xyz`, `.pdb`, `.cif` | imaginary modes exported; companions require `--convert-files` plus topology and CIF needs a bridge | imaginary-mode displacement |

`result.json` (only when `--out-json` is passed) keys:

```python
import json
d = json.load(open("result_tsopt/result.json"))
print(d["status"])                      # "converged" / "not_converged" / "stalled"; "error" on failure
print(d["energy_hartree"])
print(d["n_imaginary_modes"])           # should be 1 for a real TS
print(d["imaginary_frequencies_cm"])    # list of cm⁻¹
print(d["rigid_projection"]["treatment"], d["rigid_projection"]["effective_rank"])
print(d["rigid_projection"]["source"], d["rigid_projection"]["raw_hessian_shape"])
print(d["files"]["final_geometry_xyz"]) # path under out_dir
```

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `hess` / `rsprfo` (default) | RS-P-RFO with full Hessian | Tested default. Full-Hessian cost and convergence relative to the other modes are system/backend dependent |
| `grad` / `dimer` | Hessian-Guided Dimer | Alternative when RS-P-RFO fails to converge or full-Hessian recomputation is prohibitive on large clusters. |
| `trim` | TRIM / Helgaker | Standalone `--opt-mode` value (not an alias of grad/hess) |
| `rsirfo` | RS-I-RFO | Standalone `--opt-mode` value (not an alias of grad/hess) |


## Validation: imaginary modes

A real TS has exactly one imaginary frequency that corresponds to the
reaction coordinate.

```python
import json
d = json.load(open("result_tsopt/result.json"))
if d["status"] != "converged":
    print("NOT CONVERGED:", d["status"])
elif d["n_imaginary_modes"] == 1:
    print("OK: single imaginary mode at", d["imaginary_frequencies_cm"][0], "cm-1")
elif d["n_imaginary_modes"] == 0:
    print("BAD: collapsed to a minimum during refinement")
elif d["n_imaginary_modes"] is not None and d["n_imaginary_modes"] > 1:
    print("NOT VALIDATED: multiple imaginary modes; inspect vib/imag_*_trj.xyz")
```

For multi-imaginary cases, `optimization_status` remains the numerical optimizer
outcome while `saddle_validation` is `higher_order`. Visualize every
`vib/imag_*_trj.xyz`; `.pdb`/`.cif` companions exist only with conversion
topology. Diagnose constraints, precision, and mode character, then obtain a
better seed or use the documented flatten workflow rather than accepting it as
a first-order saddle. `all` may run only warning-labelled diagnostic IRC when a
valid negative root is available.

## Caveats

- A converged `tsopt` is **not** a complete validation; always follow
  with `irc.md` to confirm the TS connects the expected R and P.
- Do not add `--ref-mode` merely because standalone TSOPT is difficult. It
  is an MEP handoff, not a general convergence switch. Supply it manually
  only when an external path provides a deliberate non-zero Cartesian 3N
  reaction direction.
- The fixed constrained treatment controls rigid-mode removal in Cartesian
  PHVA. It is distinct from `--ref-mode`, which identifies the intended
  reaction direction. See `freeze-atoms.md` for constrained-rank behavior.
- `--max-cycles` defaults to 100000 as a safety upper bound; if a run
  burns through many cycles without converging, inspect the trajectory and
  seed rather than assuming that more cycles will repair it. Add `--flatten`
  for surplus imaginary modes. If the seed came
  from `all`, `--refine-path` can produce a better HEI, but recursive
  segmentation may multiply MEP/TSOPT/IRC/freq cost and is therefore off by
  default; inspect the coarse MEP before enabling it.
- Curvature is precision-sensitive on every backend. Keep the backend's
  effective precision (Orb defaults to fp64), inspect the actual mode, and
  confirm the result with an independent `freq` plus IRC rather than ranking a
  saddle from optimizer convergence alone.
- All four built-in backends implement analytical Hessians. For UMA only, its
  multi-worker predictor cannot expose the autograd model, so combine
  `Analytical` with `--workers 1` or use explicit finite differences.

## See also

- `path-search.md` — produces TS candidates for `tsopt`.
- `irc.md`, `freq.md` — downstream validation.
- [`pdb2reaction-install-backends/uma.md`](../pdb2reaction-install-backends/uma.md) / `mace.md` — backend-specific setup and validation notes.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.RSIRFO_KW, d.DIMER_KW, d.HESSIAN_DIMER_KW)`
