# `sp`

`pdb2reaction sp` evaluates the MLIP energy and atomic forces (optionally a Hessian) at a single geometry. Use it for a quick energy / forces / Hessian sanity check on a structure before running an optimization, for comparing backends head-to-head, or for generating reference values and Hessians outside the optimizer loop.

## Examples

Command form:

```bash
pdb2reaction sp -i FILE [-q INT | -l 'RES:Q,...'] [-m INT] [-b uma|orb|mace|aimnet2] [--hess] [options]
```

Energy + forces (UMA backend, neutral closed-shell):

```bash
# energy + forces (UMA backend, neutral closed-shell)
pdb2reaction sp -i structure.pdb -q 0 -m 1
```

Also compute the full Hessian (finite differences are used by default for every backend):

```bash
# also compute the full Hessian (FiniteDifference by default)
pdb2reaction sp -i structure.pdb -q 0 -m 1 --hess
```

## Outputs

`sp` writes its outputs under `result_sp/` by default. After a successful
calculation, the scalar energy and `|force|_max` are printed to stdout and
`forces.npy` (plus `hessian.npy` with `--hess`) is written there.

| file | contents | written |
|---|---|---|
| _stdout_ | scalar energy (a.u.) and `|force|_max`, printed as `[sp] energy = …` | successful calculation |
| `forces.npy` | `(N, 3)` array of forces in atomic units (Hartree / Bohr) | successful calculation |
| `hessian.npy` | mass-unweighted Hessian (Hartree / Bohr²): `(3N, 3N)` without frozen atoms, or the active block with YAML `geom.freeze_atoms` | only with `--hess` |
| `result.json` / `summary.json` | machine-readable energy (a.u.), backend, charge/spin, paths to npy outputs, elapsed time | only with `--out-json` |

`sp` writes no human-readable `summary.log`.

### Hessian backend

When `--hess` is set, `--hessian-calc-mode` selects the Hessian computation strategy:

- Every backend uses `FiniteDifference` by default.
- Pass `--hessian-calc-mode Analytical` to use a supported backend's autograd path explicitly.

UMA, ORB, MACE, and AIMNet2 all implement analytical Hessians. Use `--hessian-calc-mode Analytical` to request one explicitly, or force `FiniteDifference` for a numerical cross-check. With UMA, `workers > 1` cannot be combined with an explicit analytical Hessian request and raises an error; use `workers = 1` or finite differences.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| flag | default | meaning |
|---|---|---|
| `-i, --input FILE` | — | PDB / mmCIF / XYZ / GJF structure file (required) |
| `-q, --charge INT` | — | total charge; alternatively derive it with `-l` for residue-bearing PDB/mmCIF, while a valid GJF can inherit its header value |
| `-l, --ligand-charge TEXT` | — | per-residue charge mapping (e.g. `SAM:1,GPP:-3`), used to derive `-q` automatically |
| `-m, --multiplicity INT` | `1` | spin multiplicity, 2S+1 (optional; defaults to 1. GJF inherits the template) |
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | MLIP backend selection |
| `--hess / --no-hess` | `--no-hess` | also compute and write `hessian.npy` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | `FiniteDifference` | select the Hessian mode (only applies with `--hess`) |
| `-o, --out-dir PATH` | `./result_sp/` | output directory |
| `--precision [fp32\|fp64]` | backend-dependent | numeric precision passed to the backend |
| `--config PATH` | — | YAML config providing `calc.*`, `geom.*` defaults |
| `--out-json / --no-out-json` | `--no-out-json` | also write a machine-readable `result.json` (mirrored to `summary.json`) into the output directory |
| `--show-config / --dry-run` | off | print effective merged config / validate without running |

Run `pdb2reaction sp --help-advanced` for the full list (workers, solvent corrections, etc.).

## Notes

- `sp` has no freeze CLI flag, but honors the 1-based YAML
  `geom.freeze_atoms` list. Frozen forces are zeroed by the backend geometry
  contract and `--hess` writes the active partial-Hessian block by default.
- For single-point DFT (gpu4pyscf / PySCF) benchmarking use [`dft`](dft.md) instead.

## See Also

- [`opt`](opt.md) — optimize the structure
- [`tsopt`](tsopt.md) — refine a TS candidate
- [`freq`](freq.md) — vibrational analysis with thermochemistry
- [`dft`](dft.md) — single-point DFT counterpart (uses PySCF / gpu4pyscf)
