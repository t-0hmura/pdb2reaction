# `pdb2reaction opt`

## When to use

`opt` relaxes one structure to a local minimum with L-BFGS (default) or RFO.
Use it to relax R and P before `path-opt` / `path-search`, and to optimize
IRC endpoints.

## Minimal run

```bash
pdb2reaction opt -i my.pdb -l 'SAM:1' -b uma --out-json -o result_opt
```

Success: the console prints `[opt] Converged!`.

RFO when L-BFGS has trouble:

```bash
pdb2reaction opt -i my.xyz -q -1 -m 1 --opt-mode rfo -b mace -o result_opt_rfo
```

Relax the endpoints before `path-opt`:

```bash
pdb2reaction opt -i 1.R.pdb -q 0 -m 1 -o result_opt_R
pdb2reaction opt -i 3.P.pdb -q 0 -m 1 -o result_opt_P
pdb2reaction path-opt -i result_opt_R/final_geometry.pdb result_opt_P/final_geometry.pdb \
    -q 0 -m 1 -o result_path_opt
```

## Judge success

`optimization_status` in `result.json` is `converged`, `not_converged`
(reached `--max-cycles`), or `stalled` (energy plateau with `--stop-plateau`).
Only `converged` counts. The file also records `n_opt_cycles`,
`energy_hartree`, `final_max_force`, `final_rms_force`, and
`files.final_geometry_xyz`.

The default `--thresh gau` matches Gaussian's default. With `--thresh baker`,
convergence requires ALL of `max(|force|) <= 3e-4`, `rms(force) <= 2e-4`,
`max(|step|) <= 3e-4`, `rms(step) <= 2e-4` and `|delta E| < 1e-6`. This is a
deliberately tightened variant of the published criterion.

Convergence gives a stationary point, not necessarily a minimum: run `freq`
and check n_imag = 0.

Files: `final_geometry.xyz`; with `--convert-files` (on by default), also
`.pdb` for PDB/mmCIF input and `.cif` for mmCIF or very large PDB input.
`--dump` adds `optimization_trj.xyz`, and `--out-json` adds `result.json`.

## Choosing --opt-mode

| Mode | Algorithm | When |
|---|---|---|
| `grad` / `lbfgs` (default) | L-BFGS | Gradient history, no initial Hessian |
| `hess` / `rfo` | RFO with Hessian updates | L-BFGS oscillates or its history is poorly conditioned; cost depends on the system |

## Pitfalls and recovery

- **Not a TS optimizer.** For a TS, use [tsopt.md](tsopt.md).
- **Imaginary mode after convergence.** If `freq` shows a chemically
  meaningful imaginary mode, improve the starting geometry or retry with
  `--opt-mode rfo`, then check again.
- **`--reject-uphill`** (off by default) rejects an RFO trial that raises the
  energy by more than 1e-4 Hartree, restores the lower-energy geometry, and
  shrinks the trust radius. At the smallest trust radius it runs one final
  convergence check on the retained geometry. L-BFGS ignores it.
- **Other settings.** Override step limits, trust radius, and the like with
  `--config` YAML; see `OPT_BASE_KW` and `LBFGS_KW`.
- **Frozen atoms.** The removal of rigid motions applies only with `--flatten`
  and does not change L-BFGS or RFO steps. See
  [PHVA treatment](extract.md#phva-treatment).

## Next step

- [freq.md](freq.md): confirm the minimum (n_imag = 0).
- [path.md](path.md): MEP between relaxed endpoints. [tsopt.md](tsopt.md): the TS counterpart.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.OPT_BASE_KW, d.LBFGS_KW, d.RFO_KW)`
