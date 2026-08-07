# `pdb2reaction opt`

## Purpose

Single-structure geometry optimization with L-BFGS or RFO.
Use this to relax a starting geometry toward a local minimum
before feeding it to `path-search` / `path-opt`, or as a post-IRC
endpoint refinement.

## Synopsis

```bash
pdb2reaction opt -i input.pdb [-q 0 -m 1] \
    [--opt-mode grad|hess|lbfgs|rfo] \
    [-b uma|orb|mace|aimnet2] [-o ./result_opt/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | `.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `--opt-mode` | str | `grad` | `grad` (L-BFGS) or `hess` (RFO); aliases `lbfgs` / `rfo` |
| `--max-cycles` | int | `10000` | Stop after N cycles; see `OPT_BASE_KW["max_cycles"]` |
| `--reject-uphill / --no-reject-uphill` | toggle | off | Opt in to rejecting an energy-raising Hessian/RFO trial above `1e-4` Hartree, restoring the lower-energy geometry and shrinking the trust radius. At the emergency floor, run one final convergence check on the retained geometry. Ignored in L-BFGS mode. |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_opt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

With `--thresh baker`, convergence requires ALL of `max(|force|) <= 3e-4`,
`rms(force) <= 2e-4`, `max(|step|) <= 3e-4`, `rms(step) <= 2e-4` and
`|delta E| < 1e-6`. This is a deliberately tightened variant of the published
criterion (Bakken and Helgaker, J. Chem. Phys. 117, 9160 (2002)), which requires
only `max(|force|)` and (`|delta E|` or `max(|step|)`); the looser form accepts
geometries whose remaining RMS force still displaces the structure.

## Examples

### Default L-BFGS

```bash
pdb2reaction opt -i my.pdb -l 'SAM:1' -b uma -o result_opt
```

### RFO alternative when L-BFGS is problematic

```bash
pdb2reaction opt -i my.xyz -q -1 -m 1 --opt-mode rfo -b mace -o result_opt_rfo
```

### Pre-relax endpoints before path-opt

```bash
pdb2reaction opt -i 1.R.pdb -q 0 -m 1 -o /tmp/relax_R
pdb2reaction opt -i 3.P.pdb -q 0 -m 1 -o /tmp/relax_P
pdb2reaction path-opt -i /tmp/relax_R/final_geometry.pdb /tmp/relax_P/final_geometry.pdb \
    -q 0 -m 1 -o result_path_opt
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/final_geometry.xyz` | completed optimizer run | final geometry; inspect `result.json["status"]` before calling it converged |
| `<out_dir>/final_geometry.pdb` | `--convert-files` (default on) and PDB/mmCIF topology/reference available | normalized PDB companion used between pipeline stages |
| `<out_dir>/final_geometry.cif` | `--convert-files` and input/reference required the mmCIF or oversized-PDB bridge | public companion with original chain/residue IDs |
| `<out_dir>/optimization_trj.xyz` | `--dump` | full optimization trajectory |
| `<out_dir>/optimization.{pdb,cif}` | `--dump`, `--convert-files`, and conversion topology; CIF only for bridge inputs | topology-bearing trajectory companions |

`result.json` (only when `--out-json` is passed) keys: `status`
(`converged` / `not_converged`; `error` on failure), `n_opt_cycles`, `energy_hartree`,
`final_max_force`, `final_rms_force`, and the `files` block whose
`final_geometry_xyz` entry points at the final attempted geometry. When
`--flatten` runs, `rigid_projection` records the treatment, effective rank,
and raw Hessian source and shape.

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `grad` / `lbfgs` | L-BFGS | Software default; gradient-history method with no full initial Hessian |
| `hess` / `rfo` | RFO with Hessian updates | Alternative when L-BFGS oscillates or its step history is poorly conditioned; relative cost/convergence is system dependent |

## Caveats

- Not a TS optimizer — for TS use `tsopt.md`.
- Optimizer convergence alone does not prove a minimum. Run `freq`; if a
  chemically meaningful imaginary mode remains, improve the starting
  geometry or retry with `--opt-mode rfo`, then verify again.
- `--config` YAML is the way to override less-common settings (step
  limits, trust radius, etc.); inspect `OPT_BASE_KW` and `LBFGS_KW`
  in `pdb2reaction.core.defaults`.
- The fixed constrained rigid-mode treatment applies only when `--flatten`
  runs; it does not change L-BFGS or RFO steps. See `freeze-atoms.md`.

## See also

- `tsopt.md` — TS analog.
- `freq.md` — verify the optimized minimum (zero imaginary modes).
- Defaults: `import pdb2reaction.core.defaults as d; print(d.OPT_BASE_KW, d.LBFGS_KW, d.RFO_KW)`
