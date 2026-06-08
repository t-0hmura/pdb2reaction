# `pdb2reaction opt`

## Purpose

Single-structure geometry optimization with L-BFGS or RFO.
Use this to relax a starting geometry to its nearest local minimum
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
| `-i, --input` | path | required | `.pdb` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `grad` | `grad` (L-BFGS) or `hess` (RFO); aliases `lbfgs` / `rfo` |
| `--max-cycles` | int | `10000` | Stop after N cycles; see `OPT_BASE_KW["max_cycles"]` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_opt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default L-BFGS

```bash
pdb2reaction opt -i my.pdb -l 'SAM:1' -b uma -o result_opt
```

### RFO for stiffer convergence

```bash
pdb2reaction opt -i my.xyz -q -1 -m 1 --opt-mode rfo -b mace -o result_opt_rfo
```

### Pre-relax endpoints before path-opt

```bash
pdb2reaction opt -i 1.R.pdb -l '...' -o /tmp/relax_R
pdb2reaction opt -i 3.P.pdb -l '...' -o /tmp/relax_P
pdb2reaction path-opt -i /tmp/relax_R/final_geometry.xyz /tmp/relax_P/final_geometry.xyz ...
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/final_geometry.xyz` | always | converged geometry |
| `<out_dir>/final_geometry.pdb` | `--convert-files` (default on) | PDB companion |
| `<out_dir>/optimization_trj.xyz` | `--dump` | full optimization trajectory |

`result.json` (only when `--out-json` is passed) keys: `status`
(`converged` / `not_converged`; `error` on failure), `n_opt_cycles`, `energy_hartree`,
`final_max_force`, `final_rms_force`, and the `files` block whose
`final_geometry_xyz` entry points at the converged geometry.

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `grad` / `lbfgs` | L-BFGS | Default, fast, robust for most well-conditioned minima |
| `hess` / `rfo` | RFO with Hessian updates | Stiffer convergence; useful when L-BFGS oscillates |

## Caveats

- Not a TS optimizer — for TS use `tsopt.md`.
- L-BFGS occasionally walks past a saddle on shallow surfaces; if the
  resulting geometry has imaginary frequencies (run `freq` to check),
  re-run with `--opt-mode rfo`.
- `--config` YAML is the way to override less-common settings (step
  limits, trust radius, etc.); inspect `OPT_BASE_KW` and `LBFGS_KW`
  in `pdb2reaction.core.defaults`.

## See also

- `tsopt.md` — TS analog.
- `freq.md` — verify the optimized minimum (zero imaginary modes).
- Defaults: `import pdb2reaction.core.defaults as d; print(d.OPT_BASE_KW, d.LBFGS_KW, d.RFO_KW)`