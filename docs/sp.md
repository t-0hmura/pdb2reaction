# `sp`

`pdb2reaction sp` evaluates the MLIP energy and atomic forces (optionally the full Hessian) at a single geometry. Use it for a quick energy / forces / Hessian sanity check on a structure before running an optimization, for comparing backends head-to-head, or for generating reference values and Hessians outside the optimizer loop.

## Examples

Command form:

```bash
pdb2reaction sp -i FILE -q INT -m INT [-b uma|orb|mace|aimnet2] [--hess] [options]
```

Energy + forces (UMA backend, neutral closed-shell):

```bash
# energy + forces (UMA backend, neutral closed-shell)
pdb2reaction sp -i structure.pdb -q 0 -m 1
```

Also compute the full Hessian (UMA → analytical; other backends → finite-difference):

```bash
# also compute the full Hessian (UMA → analytical; other backends → finite-difference)
pdb2reaction sp -i structure.pdb -q 0 -m 1 --hess
```

## Outputs

`sp` writes its outputs under `result_sp/` by default. The scalar energy and `|force|_max` are printed to stdout; `forces.npy` (and `hessian.npy` with `--hess`) are always written there.

| file | contents | written |
|---|---|---|
| _stdout_ | scalar energy (a.u.) and `|force|_max`, printed as `[sp] energy = …` | always |
| `forces.npy` | `(N, 3)` array of forces in atomic units (Hartree / Bohr) | always |
| `hessian.npy` | `(3N, 3N)` mass-unweighted Hessian (Hartree / Bohr²) | only with `--hess` |
| `result.json` / `summary.json` | machine-readable energy (a.u.), backend, charge/spin, paths to npy outputs, elapsed time | only with `--out-json` |

`sp` writes no human-readable `summary.log`.

### Hessian backend

When `--hess` is set, the backend choice picks the Hessian computation strategy:

- `--backend uma` (default) → `Analytical` Hessian via the UMA torch autograd path
- `--backend orb` / `mace` / `aimnet2` → falls back to `FiniteDifference` (no analytical Hessian implementation upstream)

`--hessian-calc-mode` lets you force `FiniteDifference` even on UMA if you want a sanity check against the analytical implementation. Forcing `Analytical` on a non-UMA backend has no effect — those backends always fall back to `FiniteDifference`.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| flag | default | meaning |
|---|---|---|
| `-i, --input FILE` | — | PDB / XYZ / GJF structure file (required) |
| `-q, --charge INT` | — | total charge of the system (required for non-GJF; GJF inherits the template) |
| `-l, --ligand-charge TEXT` | — | per-residue charge mapping (e.g. `SAM:1,GPP:-3`), used to derive `-q` automatically |
| `-m, --multiplicity INT` | `1` | spin multiplicity, 2S+1 (required for non-GJF; GJF inherits the template) |
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | MLIP backend selection |
| `--hess / --no-hess` | `--no-hess` | also compute and write `hessian.npy` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | auto | force a specific Hessian mode (only applies with `--hess`) |
| `-o, --out-dir PATH` | `./result_sp/` | output directory |
| `--precision [fp32\|fp64]` | `fp32` | numeric precision passed to the backend |
| `--config PATH` | — | YAML config providing `calc.*`, `geom.*` defaults |
| `--out-json / --no-out-json` | `--no-out-json` | also write a machine-readable `result.json` (mirrored to `summary.json`) into the output directory |
| `--show-config / --dry-run` | off | print effective merged config / validate without running |

Run `pdb2reaction sp --help-advanced` for the full list (workers, solvent corrections, etc.).

## Notes

- For single-point DFT (gpu4pyscf / PySCF) benchmarking use [`dft`](dft.md) instead.

## See Also

- [`opt`](opt.md) — optimize the structure
- [`tsopt`](tsopt.md) — refine a TS candidate
- [`freq`](freq.md) — vibrational analysis with thermochemistry
- [`dft`](dft.md) — single-point DFT counterpart (uses PySCF / gpu4pyscf)
