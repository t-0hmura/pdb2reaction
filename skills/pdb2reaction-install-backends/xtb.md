# Experimental xTB solvent correction

`--solvent` adds an xTB solvent-minus-vacuum delta to the base MLIP surface:

```text
ΔE = E_xTB(solvent) - E_xTB(vacuum)
E_total = E_base + ΔE
```

Forces and Hessians use the corresponding difference. The correction is
experimental, disabled by default, and computationally expensive because it
runs both solvent and vacuum xTB calculations.

## Install

`pdb2reaction` calls the standalone `xtb` executable, not its Python bindings:

```bash
conda install -c conda-forge xtb
xtb --version
```

The executable must remain on `PATH` in batch jobs. Configure the solvent and
executable under `calc` as listed in `docs/yaml-reference.md`; use
`--help-advanced` for the corresponding CLI options.

`trj2fig` also accepts it, but only uses it when `-q` or `-m` triggers frame
energy recomputation. It is not accepted by `dft`, `extract`, or
structure-only utilities.
