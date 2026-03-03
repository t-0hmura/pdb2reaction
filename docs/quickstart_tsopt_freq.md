# Quickstart: `pdb2reaction tsopt`

## Goal

Optimize a TS candidate and verify that it is a first-order saddle point.

## Prerequisites

- TS candidate geometry: `ts_guess.pdb`
- Charge and multiplicity (`-q`, `-m`) for the target state

## 1. TS optimization

```bash
pdb2reaction tsopt -i ts_guess.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

`tsopt` performs a final Hessian evaluation and imaginary-frequency check automatically at the end of optimization. Check the terminal output for lines like:

```
[Imaginary modes] n=1  ([-593.1])
```

## What to check

- `result_tsopt/final_geometry.pdb` — optimized TS structure
- `result_tsopt/vib/` — animation files for the imaginary-frequency normal mode (`final_imag_mode_*.xyz`, `.pdb`)
- Terminal output: **n=1** with a sufficiently large imaginary frequency (|ν| ≥ 100 cm⁻¹) indicates a good TS candidate

## 2. (Optional) Separate frequency analysis

A standalone `freq` run is useful when you want full vibrational frequency output or thermochemistry corrections (`--thermo` in the `all` command). If you only need the imaginary-frequency check, the `tsopt` output above is sufficient.

```bash
pdb2reaction freq -i ./result_tsopt/final_geometry.pdb -q 0 -m 1 --out-dir ./result_freq
```

## Tips

- Use `--hessian-calc-mode Analytical` when VRAM is sufficient.
- Check full options with `pdb2reaction tsopt --help-advanced` and `pdb2reaction freq --help-advanced`.

## Next step

- Trace the path with [irc](irc.md), or run the end-to-end route with [all](all.md).
