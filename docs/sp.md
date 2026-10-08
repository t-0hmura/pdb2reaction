# `sp` (single point)

`sp` computes the **energy and atomic forces** of one structure with the selected backend, and with `--hess` also the **Hessian**. It runs no optimization: the geometry stays as given.

---

## What it is for

* **Check before an optimization**: confirm that the charge and multiplicity are accepted and that the backend returns a finite energy and forces.
* **Compare backends**: evaluate the same structure with the [MLIPs](backends.md) (machine-learning interatomic potentials) UMA, ORB, MACE, and AIMNet2, or with DFT (`-b dft`).
* **Reference values**: forces and Hessians as `.npy` files, and the energy in the console or `result.json`, for your own analysis.

---

## Examples

### 1. Energy and forces

Evaluate a neutral singlet with the default backend (UMA).

```bash
pdb2reaction sp -i structure.pdb -q 0 -m 1 --out-json
```

The console prints `[sp] energy = … a.u.  |force|_max = … a.u./bohr`, and `result_sp/` has `forces.npy` and `result.json` with `energy_au`.

### 2. Add the full Hessian

`--hess` also computes the Hessian.

```bash
pdb2reaction sp -i structure.pdb -q 0 -m 1 --hess
```

---

## How it works

1. **Reading the structure**:
PDB, mmCIF, XYZ, and GJF inputs are read. The charge comes from `-q`, `-l` (PDB/mmCIF input), YAML `calc.charge`, or a `.gjf` header. Atoms given with `--freeze-atoms` are frozen.
2. **Energy and forces**:
The backend is called once at the input geometry. `sp` prints the energy and the largest force component and saves the forces to `forces.npy`.
3. **Hessian (with `--hess`)**:
`--hessian-calc-mode FiniteDifference` differentiates the forces numerically; `Analytical` uses the analytical Hessian of UMA, ORB, MACE, AIMNet2, or DFT. With UMA, `Analytical` {ref}`cannot run <workers-analytical-error>` with `--uma-workers` above 1.

---

## Output files

`sp` writes these files to `--out-dir`:

| File | Contents | Written |
| --- | --- | --- |
| `forces.npy` | Forces as an `(N, 3)` array in Hartree/bohr | Always |
| `hessian.npy` | Cartesian Hessian without mass weighting (Hartree/bohr²): `(3N, 3N)`, or `(3M, 3M)` for the M moving atoms in input order when atoms are frozen | With `--hess` |
| `result.json` | Energy (`energy_au`), backend, model, charge, multiplicity, atom count, paths to the `.npy` files, elapsed time | With `--out-json` |
| `summary.json` | Copy of `result.json`; read `result.json` | With `--out-json` |

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input structure (`.pdb`, `.cif`, `.xyz`, `.gjf`, ...) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l`, YAML `calc.charge`, or a `.gjf` input gives it |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1); a `.gjf` input supplies its own |
| `-l, --ligand-charge` | text | `None` | Per-residue formal charges (e.g. `'SAM:1,GPP:-3'`) or one total ligand charge. Needs PDB/mmCIF input |
| `-b, --backend` | text | `uma` | Calculator (`uma`, `orb`, `mace`, `aimnet2`, `dft`); for the `-b dft` settings see [Refine an MLIP TS with DFT](dft-backend.md) |
| `--hess/--no-hess` | flag | `False` | Also compute the Hessian and write `hessian.npy` |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | Hessian method (finite difference / analytical); used with `--hess` |
| `--freeze-atoms` | text | `None` | Atoms to freeze (1-based, comma-separated, e.g. `'1,3,5'`) |
| `-o, --out-dir` | path | `./result_sp/` | Output directory |
| `--out-json/--no-out-json` | flag | `False` | Write `result.json` and `summary.json` |

See the [generated CLI reference](reference/commands/sp.md) for every option.

> **Note:** In YAML (`--config`), `calc` sets the backend and `geom.freeze_atoms` adds frozen atoms (1-based), merged with `--freeze-atoms`.

---

## Notes

* **Energy looks wrong**: re-check the {ref}`charge and multiplicity <charge-spin-problems>`.
* **Frozen atoms** get zero force.
* **Cap hydrogens**: `sp` does not freeze the parent atoms of the cap hydrogens that `extract` adds; list them in `--freeze-atoms` if you want them fixed.
* **Atomic charges**: `sp -b dft` gives the DFT energy and forces only. For Mulliken, meta-Löwdin, and IAO charges, use [`dft`](dft.md).
* **A failed run** prints a one-line `Error: …` or `Unhandled error during single-point calculation:` with a traceback, and exits with a nonzero code; see [Error handling](json-output.md#error-handling).
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [opt](opt.md) — optimize the structure
* [tsopt](tsopt.md) — optimize a transition-state (TS) candidate
* [freq](freq.md) — vibrational analysis and thermochemistry
* [dft](dft.md) — DFT single point with atomic charges
* [MLIP Backends](backends.md) — choosing a backend, and the Hugging Face login that UMA needs
* [Refine an MLIP TS with DFT](dft-backend.md) — `-b dft` settings and GPU memory
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
