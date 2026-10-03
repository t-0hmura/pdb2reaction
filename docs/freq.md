# `freq` (vibrational analysis and thermochemistry)

## Overview

`freq` computes **harmonic vibrational frequencies** and **thermochemical corrections** (ZPE, enthalpy, Gibbs free energy) for a structure.

### What it is for

* **Checking a stationary point**: count the imaginary frequencies (n_imag) to confirm a minimum (n_imag = 0) or a transition state (TS, n_imag = 1).
* **Thermochemistry**: free energies and other thermodynamic quantities from the QRRHO (quasi-rigid-rotor harmonic oscillator) model.
* **Seeing the modes**: atomic-displacement animations of the imaginary mode or any other mode.

The default backend is **UMA**, Meta's pretrained [machine-learning interatomic potential (MLIP)](backends.md); `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or [**DFT**](dft-backend.md).

---

## Examples

### 1. Minimal run (explicit charge and multiplicity)

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

The console summary prints n_imag as `Number of Imaginary Freq = N`.

### 2. Detailed thermochemistry file

Add `--dump` to also write the detailed thermochemistry file `thermoanalysis.yaml`.

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --dump --out-dir ./result_freq_dump
```

### 3. Analytical Hessian

Use this to avoid the step-size error of finite differences.

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 \
  --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

---

## How it works

1. **Reading the structure and freezing the boundary (PHVA)**:
With `--freeze-links` (on by default), `freq` finds the cap hydrogens that `extract` adds (atom `HL` in residue `LKH`) and freezes their parent atoms. When any atoms are frozen, the vibrational analysis runs on the movable atoms only (PHVA, partial Hessian vibrational analysis).
2. **Hessian mode**:
`--hessian-calc-mode` selects `FiniteDifference` (finite differences, the default) or `Analytical`.
3. **Thermochemistry (QRRHO)**:
From the frequencies, the QRRHO model, which corrects the entropy of low frequencies (rotor cutoff 100 cm⁻¹), gives the Gibbs free-energy correction `G_corr`. G is E + `G_corr`, where E is the electronic energy, and always includes the translational and rotational terms of the whole structure besides the vibrational terms of the positive frequencies. The console summary prints `G_corr` and G in Hartree on the lines `Gibbs Free Energy Correction (G_corr)` and `Gibbs Free Energy (G = E + G_corr)`; `--dump` also writes them to `thermoanalysis.yaml`, and `--out-json` to `thermochemistry` in `result.json`, where G is `sum_EE_and_thermal_free_energy_ha`. The point group and rotational symmetry number are detected from the structure; YAML `thermo.symmetry_number` overrides the number.
4. **Writing the modes**:
Up to `--max-write` mode animations are written, starting from the imaginary or lowest modes.

### Rigid modes with frozen boundaries

Without frozen atoms, `freq` removes the six rigid motions (three translations and three rotations), which are not vibrations, before reporting frequencies. With frozen atoms, it removes only the rigid motions that keep every frozen atom in place. With three or more frozen atoms that do not lie on one line, the normal case for a cluster model, nothing is removed and every vibrational mode of the movable atoms is kept. With one frozen atom, three motions are removed (rotations about that atom); with two, one is removed (rotation about the axis through them).

`irc`, the TS frequency check and the Dimer direction in `tsopt`, and `--flatten` (removing extra imaginary modes) in `opt` and `tsopt` treat rigid motions the same way. With `--out-json`, `result.json` records the number of removed motions and the Hessian used under `rigid_projection`; see [JSON Output Reference](json-output.md#rigid-projection-provenance).

---

## Reading the frequencies

How `frequencies_cm-1.txt` and the JSON record treat each case:

| Item | Value / behavior | Meaning |
| --- | --- | --- |
| **Imaginary modes** | Negative values (ν < 0 cm⁻¹) | An imaginary mode is listed as a negative frequency. |
| **Imaginary threshold** | ν < −5.00 cm⁻¹ | Such a mode counts as imaginary (n_imag; the JSON field is `n_imaginary`). The cutoff is YAML `freq.zero_cutoff_cm` (default `5.0`). |
| **Tiny negative modes** | −5.00 ≤ ν < 0 cm⁻¹ | Small negative modes from numerical noise do not count toward n_imag. `n_negative_modes` counts every negative frequency, including these. |
| **Thermochemistry** | No inversion, no floor | Imaginary modes are not flipped and small positive modes are not raised. QRRHO uses only the positive modes, so imaginary modes are left out of ZPE and G. |

---

## Output files

`freq` writes these files to `--out-dir`:

```text
result_freq/
├─ frequencies_cm-1.txt          # All frequencies (cm⁻¹)
├─ mode_0001_-385.20cm-1_trj.xyz # Animation of each mode (XYZ)
├─ mode_0001_-385.20cm-1.pdb     # Same animation as PDB (PDB or mmCIF input)
├─ mode_0001_-385.20cm-1.cif     # Same animation as mmCIF (mmCIF or very large PDB input)
├─ thermoanalysis.yaml           # Detailed thermochemistry (with --dump)
└─ result.json                   # Summary (--out-json)
```

* **Minimum or TS?** The console summary prints n_imag as `Number of Imaginary Freq = N`; `freq` does not judge it, so `scientific_status` in `result.json` is `success` whatever n_imag is. A successful TS optimization gives one imaginary mode along the reaction coordinate: **exactly one** clear negative value at the top of `frequencies_cm-1.txt`, and every later value positive or within the tolerance. A TS then goes to [`irc`](irc.md). If a structure meant to be a minimum has imaginary modes, optimize it again with [`opt`](opt.md) `--flatten`; if a TS has none or several, see {ref}`When the TS search fails <ts-search-fails>`.
* **Watching the motion**: open `mode_*_trj.xyz` or `.pdb` in PyMOL, VMD, OVITO, or another viewer to animate the vibration.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input structure (`.pdb`, `.cif`, `.xyz`, ...) |
| `-q, --charge` | integer | `None` | Total charge. Required unless `-l` is given or the input is `.gjf` |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `-l, --ligand-charge` | text | `None` | Per-residue formal charges (e.g. `'SAM:1,GPP:-3'`) |
| `--ref-pdb` | path | `None` | Reference PDB/mmCIF topology for `.xyz` / `.gjf` input; the coordinates still come from `-i` |
| `-o, --out-dir` | path | `./result_freq/` | Output directory |
| `-b, --backend` | text | `uma` | Calculator backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | How the Hessian is computed (finite differences / analytical) |
| `--read-hess` | path | `None` | Read the Hessian from a `.npy` file (for example one saved by `freq` or `tsopt --dump-hess`) instead of computing it |
| `--dump-hess` | path | `None` | Save the Hessian as a `.npy` file for `--read-hess` in `freq`, `tsopt`, or `irc` |
| `--freeze-links/--no-freeze-links` | flag | `True` | Freeze the parent atoms of cap hydrogens at the cluster boundary |
| `--freeze-atoms` | text | `None` | Atoms to freeze (1-based, comma-separated, e.g. `'1,3,5'`) |
| `--max-write` | integer | `10` | Maximum number of mode animations to write |
| `--sort` | `value` / `abs` | `value` | Order of the modes (by value / by absolute value) |
| `--temperature` | float | `298.15` | Temperature for thermochemistry (K) |
| `--pressure` | float | `1.0` | Pressure for thermochemistry (atm) |
| `--dump/--no-dump` | flag | `False` | Write the detailed thermochemistry file `thermoanalysis.yaml` |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` |

See the [generated CLI reference](reference/commands/freq.md) for every option.

> **Note:** In YAML (`--config`), the [`freq`](yaml-reference.md#freq-section) section sets the imaginary threshold `zero_cutoff_cm` and the number and amplitude of the written modes, and the [`thermo`](yaml-reference.md#thermo) section sets temperature and pressure.

---

## Notes

* **`freq` or `tsopt`?** `tsopt` already checks the imaginary frequencies. Run `freq` on its own when you need detailed thermochemistry or mode animations.
* **At least one atom must move.** If every atom is frozen there is no vibration to analyze, and `freq` stops with an error.
* **Analytical Hessian and `--uma-workers`**: with UMA, `--hessian-calc-mode Analytical` cannot run with `--uma-workers` above 1 and stops with an error. Use `--uma-workers 1` for an analytical Hessian. It uses more GPU memory, so test it on your system first.
* **`all --thermo` keeps the thermochemistry file**: `all` reads its thermochemistry from `thermoanalysis.yaml`, so with `--thermo` it writes this file even under `--no-dump`.
* **The `--read-hess` / `--dump-hess` file** is one NumPy array (`numpy.save`): the Cartesian Hessian in Hartree/bohr², not mass-weighted, with atoms in input order. It covers all atoms (3N × 3N) or, when atoms are frozen, only the movable ones. `--read-hess` checks only the size, finiteness, and symmetry of the matrix, so pass a Hessian computed for the same geometry, charge, multiplicity, and calculator.

---

## See also

* [opt](opt.md) — geometry optimization to a minimum
* [tsopt](tsopt.md) — TS optimization
* [irc](irc.md) — IRC from a TS
* [dft](dft.md) — DFT single-point energies
* [all](all.md) — the full workflow: extraction, path search, TS optimization, and vibrational analysis
* [YAML Reference](yaml-reference.md) — configuration file format
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
* {ref}`Exit codes <exit-codes>` — what each exit status means
