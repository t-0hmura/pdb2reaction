# Quickstart: `pdb2reaction tsopt`

## Goal

Optimize a TS candidate and verify that it is a first-order saddle point.

## Prerequisites

- TS candidate geometry: `.pdb`
- Charge (`-q/--charge` or `-l/--ligand-charge`) and multiplicity (`-m`) for the target state

## 1. TS optimization

```bash
pdb2reaction tsopt -i ts_guess.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

`tsopt` performs a final Hessian evaluation and imaginary-frequency check automatically at the end of optimization. Check the terminal output for lines like:

```
[Imaginary modes] n=1  ([-593.1])
```

## Expected output

```text
result_tsopt/
├── final_geometry.xyz     # Optimized TS geometry
├── final_geometry.pdb     # PDB format (if input was PDB)
└── vib/
    ├── imag_-593.10cm-1.pdb       # Imaginary mode animation
    └── imag_-593.10cm-1_trj.xyz   # Trajectory format
```

**What to check:**

1. Terminal output: **`n=1`** with |frequency| >= 100 cm⁻¹ indicates a valid first-order saddle point
2. `vib/imag_*.pdb` — open in PyMOL and animate; the mode should correspond to the expected bond-breaking/forming
3. If `n=0` (no imaginary mode): the optimization converged to a minimum, not a TS. Try a different initial guess
4. If `n>1` (multiple imaginary modes): use `--flatten-max-iter 3` to attempt flattening extra modes

## 2. (Optional) Separate frequency analysis

A standalone `freq` run is useful when you want full vibrational frequency output or thermochemistry corrections (zero-point energy (ZPE), Gibbs free energy, etc.; `--thermo` in the `all` command). If you only need the imaginary-frequency check, the `tsopt` output above is sufficient.

```bash
pdb2reaction freq -i ./result_tsopt/final_geometry.pdb -q 0 -m 1 --out-dir ./result_freq
```

## Tips

- Using `--hessian-calc-mode Analytical` is recommended when VRAM is sufficient (default: `FiniteDifference`).
- Check full options with `pdb2reaction tsopt --help-advanced` and `pdb2reaction freq --help-advanced`.

## Next step

- Trace the path with [irc](irc.md), or run the end-to-end route with [all](all.md).
