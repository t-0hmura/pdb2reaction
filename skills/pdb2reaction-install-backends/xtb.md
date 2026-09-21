# xTB solvent correction

`--solvent` adds an xTB solvent-minus-vacuum delta to the base MLIP surface:

```text
ΔE = E_xTB(solvent) - E_xTB(vacuum)
E_total = E_base + ΔE
```

Forces and Hessians use the corresponding difference. The correction is
disabled by default and computationally expensive because it runs both solvent
and vacuum xTB calculations. Use it mainly for small-molecule solution
calculations. Roughly 200–300 atoms is a practical upper range, not a hard
limit; benchmark the actual system and hardware before a production run.

For a matched reaction, the solution-phase barrier can be compared with the
enzyme-cluster barrier when the reacting species, charge, multiplicity,
backend/model, and energy references are kept consistent. Continuum models are
also used with enzyme cluster models, but this option represents a named bulk
solvent. Apply it to an enzyme cluster only when that dielectric model is
scientifically justified.

## Install

`pdb2reaction` calls the standalone `xtb` executable, not its Python bindings:

```bash
conda install -c conda-forge xtb
xtb --version
```

The executable must remain on `PATH` in batch jobs. Configure the solvent and
executable under `calc` as listed in `docs/yaml-reference.md`; use
`--help-advanced` for the corresponding CLI options.

If xTB SCC convergence is poor, pass additional xTB options through the
advanced command override, for example
`--solvent-xtb-cmd 'xtb --etemp 1000'`.

`trj2fig` also accepts it, but only uses it when `-q` or `-m` triggers frame
energy recomputation. It is not accepted by `dft`, `extract`, or
structure-only utilities.
