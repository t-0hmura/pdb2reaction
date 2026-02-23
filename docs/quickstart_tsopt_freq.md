# Quickstart: `pdb2reaction tsopt` -> `pdb2reaction freq`

## Goal

Optimize a TS candidate and validate it by frequency analysis.

## Prerequisites

- TS candidate geometry: `ts_guess.pdb`
- Charge and multiplicity (`-q`, `-m`) for the target state

## 1. TS optimization

```bash
pdb2reaction tsopt -i ts_guess.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## 2. Frequency check on optimized TS

```bash
pdb2reaction freq -i ./result_tsopt/final_geometry.pdb -q 0 -m 1 --out-dir ./result_freq
```

## What to check

- `result_tsopt/final_geometry.pdb`
- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz` and `result_freq/mode_*.pdb`

For a valid first-order saddle, frequencies should contain exactly one imaginary mode (negative cm^-1).

## Tips

- Use `--hessian-calc-mode Analytical` when VRAM is sufficient.
- Check full options with `pdb2reaction tsopt --help-advanced` and `pdb2reaction freq --help-advanced`.

## Next step

- Trace the path with [irc](irc.md), or run the end-to-end route with [all](all.md).
