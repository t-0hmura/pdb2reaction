# `irc` (intrinsic reaction coordinate)

`irc` traces the intrinsic reaction coordinate (IRC) from an optimized transition state (TS) in both directions with EulerPC (an Euler predictor–corrector integrator). It writes the trajectory of each branch and the two endpoint candidates. Optimizing those endpoints with [`opt`](opt.md) shows which reactant (R) and product (P) the TS connects.

---

## What it is for

* **Checking a TS**: confirm that the TS from [`tsopt`](tsopt.md) connects the intended R and P.
* **Getting R and P**: optimize the endpoints with [`opt`](opt.md) to obtain the R and P structures of this TS.
* **Rerunning the IRC step of `all`**: trace the IRC of an [`all`](all.md) run again on its own, with different settings.

The default backend is **UMA**, Meta's pretrained [machine-learning interatomic potential (MLIP)](backends.md); `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or **DFT**.

---

## Examples

### 1. Both branches

Trace both directions from the TS, and write a summary with `--out-json`.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --out-json --out-dir ./result_irc
```

The endpoint candidates are `finished_first.xyz` and `finished_last.xyz`.

### 2. Forward branch only with a larger step

Trace only the forward branch, with a maximum step of 0.2 bohr.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward --step-size 0.2 --out-dir ./result_irc_forward
```

### 3. Analytical Hessian

Compute the starting Hessian analytically instead of by finite differences.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

### 4. Retry with a smaller step

When a branch stops after only a few frames, retry with a maximum step of 0.05 bohr.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 --out-dir ./result_irc_small_step
```

### 5. Trace to the cycle limit

Add `--never-stop` to ignore the gradient and energy stop criteria and trace each branch until `--max-cycles`.

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 --never-stop \
    --max-cycles 250 --out-dir ./result_irc_continue
```

---

## How it works

1. **Starting direction**: `irc` computes the Hessian at the TS (or reads it with `--read-hess`), removes rigid motions as [`freq`](freq.md#rigid-modes-with-frozen-boundaries) does, and takes the eigenvector `--root` (`0` = lowest eigenvalue) as the reaction mode. If that mode is not imaginary, the run stops with an error.
2. **EulerPC integration**: each branch (forward, then backward) starts from the TS. Every step is an Euler predictor along the mass-weighted steepest-descent direction, with the gradient estimated from a second-order Taylor expansion with the current Hessian (Bofill update). Each predictor step is followed by a modified Bulirsch–Stoer corrector on a DWI (distance-weighted interpolation) surface. A branch stops when the RMS gradient falls below 1 × 10⁻³ hartree/bohr after leaving the TS region, when the energy rises, when the energy changes by 1 × 10⁻⁶ hartree or less in one step, or at `--max-cycles`.
3. **Writing the path**: `irc` writes each branch, the whole path through the TS, and the two end structures of that path. For PDB/mmCIF input, the trajectories are also converted to PDB.

---

## Judging the IRC

Even if the IRC does not converge, the result is usable when the endpoints, optimized with `opt`, reach the intended R and P.

| What to check | Where to look |
| --- | --- |
| The start was a TS | The console line `Transition vector is mode 0 with wavenumber … cm⁻¹.` shows a negative wavenumber |
| How each branch stopped | `forward_integration_converged` / `backward_integration_converged` in `result.json`: `true` when the RMS gradient fell below the threshold, `false` for an energy stop or the cycle limit |
| Bonds that change along the path | `bond_changes` in `result.json` (`formed` and `broken`, from `finished_first` to `finished_last`) |
| Which end is R and which is P | Optimize `finished_first.xyz` and `finished_last.xyz` with [`opt`](opt.md) and compare them with the intended R and P. The order first / last does not decide it |

`result.json` records the outcome as `scientific_status`: `success` (exit code 0) means that the integration ran without an error, however each branch stopped. `irc` does not judge the endpoints; whether they are the intended R and P is for you to check.

If the endpoints are not the intended R and P, see {ref}`When the TS search fails <ts-search-fails>`.

---

## Output files

When the run finishes, `--out-dir` contains:

```text
result_irc/
├─ finished_irc_trj.xyz    # Whole IRC path: forward end → TS → backward end
├─ finished_irc.pdb        # Same path as PDB (PDB/mmCIF input)
├─ finished_first.xyz      # First frame: forward end (the TS with --no-forward)
├─ finished_last.xyz       # Last frame: backward end (the TS with --no-backward)
├─ {forward,backward}_{first,last}.xyz  # Ends of each branch
├─ forward_irc_trj.xyz     # Forward branch (when it runs)
├─ forward_irc.pdb         # Same branch as PDB (PDB/mmCIF input)
├─ backward_irc_trj.xyz    # Backward branch (when it runs)
├─ backward_irc.pdb        # Same branch as PDB (PDB/mmCIF input)
└─ result.json             # Summary (--out-json)
```

{ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

* **Endpoint candidates**: `finished_first.xyz` and `finished_last.xyz` are the structures to optimize with [`opt`](opt.md). Of the ends of each branch, `forward_first.xyz` and `backward_last.xyz` are the far ends, and `forward_last.xyz` and `backward_first.xyz` are the first steps from the TS.
* **Path**: open `finished_irc_trj.xyz` or `finished_irc.pdb` in PyMOL or VMD to watch the reaction.
* **Summary**: with `--out-json`, [`result.json`](json-output.md) records the number of frames of each branch, how each branch stopped, `bond_changes`, and the energies of the two ends and the TS.
* **File prefix**: with YAML `irc.prefix: trial`, the names of the trajectory and structure files start with `trial_`; `files` in `result.json` also records the prefixed names.
* **Console**: the step table of each branch and the elapsed time; with {ref}`-v 3 <verbosity-levels>`, also the `geom`, `calc`, and `irc` settings actually used.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | TS structure (`.pdb`, `.cif`, `.mmcif`, `.xyz`, `.gjf`). For a trajectory, extract one frame to `.xyz` first (see {ref}`Extract one frame from a trajectory <trajectory-one-frame>`) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l` is given or the input is `.gjf` |
| `-l, --ligand-charge` | text | `None` | Total ligand charge (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `-b, --backend` | text | `uma` | Backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--max-cycles` | integer | `125` | Maximum number of IRC steps per branch |
| `--step-size` | float | `0.10` | Maximum step length in bohr (unweighted Cartesian coordinates) |
| `--never-stop/--no-never-stop` | flag | `False` | Ignore the gradient and energy stop criteria and trace until `--max-cycles` |
| `--forward/--no-forward` | flag | `True` | Run the forward branch |
| `--backward/--no-backward` | flag | `True` | Run the backward branch |
| `--root` | integer | `0` | Hessian eigenvector used as the reaction mode, counted from 0 in ascending order of eigenvalue |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | How the starting Hessian is computed |
| `--read-hess` | path | `None` | Start from the Hessian in a `.npy` file (for example from `freq` or `tsopt --dump-hess`) instead of computing it |
| `--freeze-links/--no-freeze-links` | flag | `True` | Freeze the parent atoms of cap hydrogens (PDB/mmCIF input or `--ref-pdb`) |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |
| `-o, --out-dir` | path | `./result_irc/` | Output directory |

See the [generated CLI reference](reference/commands/irc.md) for every option.

> **Note:** in YAML (`--config`), every key of the `irc` block, including the stop thresholds, is listed under [`irc`](yaml-reference.md#irc-section) in the YAML Reference.

---

## Notes

* **A branch that stops at once**: when a branch ends after three frames or fewer, the console warns `[irc] IRC stopped after only a few frames in …`. Try example 4 first; a large step can make EulerPC unstable.
* **`--never-stop` is off by default**: numerical failures and interruptions still stop the run. Inspect the trajectory and optimize the endpoints, and raise `--max-cycles` only when the extra path is useful.
* **`--root` counts from 0**: a successful TS optimization gives one imaginary mode along the reaction coordinate, so for a TS with n_imag = 1 keep `--root 0` (the only negative eigenvalue). Use `1`, `2`, … only when you know that spurious modes with lower (more negative) eigenvalues come before the reaction mode.
* **Fixed settings**: `irc` always uses Cartesian coordinates (`geom.coord_type: cart`) and the Hessian of the movable atoms only (`calc.return_partial_hessian: true`), whatever the YAML says.
* **The `--read-hess` file** is the same `.npy` file as in [`freq`](freq.md). It needs `irc.hessian_init: calc` (the default); when the file is used, `rigid_projection.hessian_source` in `result.json` is `"file"`.
* **Analytical Hessian and `--uma-workers`**: with UMA, `--hessian-calc-mode Analytical` cannot run with `--uma-workers` above 1 and stops with an error. Use `--uma-workers 1` for an analytical Hessian. Its speed and memory use depend on the backend, the model, and the system size, so test it on your system first.
* **Optimizing the endpoints**: `finished_first.xyz` and `finished_last.xyz` are written only as `.xyz`, so pass the TS PDB with `--ref-pdb` to keep the parent atoms of the cap hydrogens {ref}`frozen <freeze-atoms-and-restraints>`:

  ```bash
  pdb2reaction opt -i result_irc/finished_first.xyz --ref-pdb ts.pdb -q 0 -m 1 --out-dir ./result_opt_first
  pdb2reaction opt -i result_irc/finished_last.xyz --ref-pdb ts.pdb -q 0 -m 1 --out-dir ./result_opt_last
  ```

* **Frozen atoms**: `--freeze-atoms` (1-based) freezes atoms in addition to `--freeze-links`. `result.json` records the removed rigid motions and the starting Hessian under `rigid_projection`.
* **Large systems**: `--hess-device cpu` keeps the starting Hessian and the IRC Hessian operations on the CPU, to stay within GPU memory.
* **At least one branch**: `--no-forward` together with `--no-backward` stops with an error.

---

## See also

* [tsopt](tsopt.md) — optimize the TS before running IRC
* [opt](opt.md) — optimize the IRC endpoints to R and P
* [freq](freq.md) — full vibrational analysis and thermochemistry
* [all](all.md) — the full workflow, which runs IRC after `tsopt` and optimizes the endpoints
* [Troubleshooting](troubleshooting.md) — when a run fails
* [YAML Reference](yaml-reference.md) — every `irc` setting
* [Glossary](glossary.md) — IRC and other terms
* {ref}`Exit codes <exit-codes>` — what each exit status means
