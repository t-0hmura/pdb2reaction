# `pdb2reaction freq`

## When to use

`freq` builds the Hessian, gives harmonic frequencies and a displacement file
per mode, and computes QRRHO thermochemistry (298.15 K and 1 atm by default).
Use it to check a minimum (n_imag = 0) or a TS (n_imag = 1), or to get free
energies. With frozen atoms it runs a partial Hessian vibrational analysis
(PHVA) on its own.

## Minimal run

```bash
pdb2reaction freq -i ts.xyz -q 0 -m 1 -b uma --out-json -o result_freq
```

Success: the console prints `Number of Imaginary Freq = N` with the N you
expect, and the Gibbs free energy. Without `--out-json` you still get the text
and trajectory files, but no `result.json`.

Higher temperature for activation enthalpy:

```bash
pdb2reaction freq -i ts.pdb -l 'SAM:1' \
    --temperature 310.15 --pressure 1.0 \
    -b uma -o result_freq_310K
```

Pass the Hessian on to IRC:

```bash
pdb2reaction freq -i ts.xyz -q 0 -m 1 --dump-hess result_freq/hessian.npy -o result_freq
pdb2reaction irc -i ts.xyz -q 0 -m 1 --read-hess result_freq/hessian.npy -o result_irc
```

`--dump-hess` writes to that exact path, relative to the current directory,
not under `--out-dir`. The file is one `numpy.save` array: the Cartesian
Hessian in Hartree/bohr², not mass-weighted, atoms in input order, 3N×3N or
only the movable atoms when some are frozen. `--read-hess` checks only size,
symmetry, and finiteness, so pass a Hessian computed for the same geometry,
charge, multiplicity, and calculator.

## Judge success

n_imag is printed as `Number of Imaginary Freq = N` and stored as
`n_imaginary`. A minimum has n_imag = 0. A successful TS optimization gives one
imaginary mode along the reaction coordinate: one clear negative value at the
top of `frequencies_cm-1.txt`. The default imaginary criterion is
ν < −5.00 cm⁻¹; `freq.zero_cutoff_cm` sets the cutoff magnitude. `freq` does not
judge n_imag, so `scientific_status` is `success` whatever n_imag is.

Files:

- `frequencies_cm-1.txt`: all modes in cm⁻¹.
- `mode_NNNN_<freq>cm-1_trj.xyz`: mode animations at the top level, not under
  `vib/`, up to `--max-write` of them. `--convert-files`, on by default, adds
  `.pdb` for PDB/mmCIF input and `.cif` for mmCIF or very large PDB input.
- `thermoanalysis.yaml` with `--dump`, and `result.json` with `--out-json`.

```python
import json
d = json.load(open("result_freq/result.json"))
print(d["n_imaginary"], d["n_negative_modes"])
print(d["frequencies_cm"][:5])                # first five frequencies (cm-1)
t = d["thermochemistry"]
print(t["electronic_energy_ha"])              # E (Hartree)
print(t["zpe_correction_ha"])                 # ZPE correction (Hartree)
print(t["thermal_correction_free_energy_ha"]) # G_corr (Hartree)
print(t["sum_EE_and_thermal_free_energy_ha"]) # G (Hartree)
print(t["S_cal_per_mol_K"])                   # entropy (cal/mol·K)
```

## Thermochemistry

QRRHO uses a 100 cm⁻¹ rotor cutoff: vibrations below it are interpolated
toward the free-rotor entropy, and higher ones use the harmonic oscillator.
The free energy is E + G_corr = G, where E is the electronic energy. `THERMO_KW`
(`pdb2reaction.core.defaults`) exposes `temperature`, `pressure_atm`, an
optional `symmetry_number` override, and `dump`. The point group and
rotational symmetry are detected from each structure, and the `1/sigma`
correction is always included. The rotor cutoff is not a `THERMO_KW` setting.

## PHVA

When atoms are frozen, `freq` builds and diagonalizes only the block of the
movable atoms. This shrinks the dense matrix; the actual time and memory
saving depends on the number of movable atoms and the backend. Only rigid
motions that keep every frozen atom in place are removed, so a normal cluster
boundary with three or more frozen atoms off one line loses no mode. `tsopt`
and `irc` use the same treatment; the [JSON record](../../docs/json-output.md#rigid-projection-provenance)
shows what was removed.

`pdb2reaction` does not read PDB B-factors as a freeze list. The frozen set is
the union of `--freeze-atoms` with 1-based indices, YAML `geom.freeze_atoms`,
and `--freeze-links`, which freezes the parents of the `LKH` / `HL` cap
hydrogens written by `extract`. See
[Freeze atoms at the cluster boundary](extract.md#freeze-atoms-at-the-cluster-boundary).

## Pitfalls and recovery

- **Imaginary counts.** `all` can still follow a valid negative mode from a
  converged higher-order candidate as a diagnostic IRC; that is not a
  first-order TS. Imaginary modes in R or P do not block thermochemistry.
- **All atoms frozen.** No vibration is left, and `freq` raises an error.
- **Signed modes.** `freq` retains every signed physical mode; imaginary modes
  are not flipped, and positive modes between 0 and 5 cm⁻¹ stay in
  thermochemistry. Raw negative counts are diagnostic and do not add a pipeline
  failure gate. Inspect all of `frequencies_cm`, `n_negative_modes`, and the
  displacements.
- **Small imaginary frequency.** It may be numerical noise or a real shallow
  mode. Inspect its displacement and repeat the Hessian at a suitable
  precision. The 100 cm⁻¹ QRRHO cutoff does not turn an imaginary mode real or
  validate a stationary point.
- **Hessian mode.** `--hessian-calc-mode FiniteDifference` usually lowers peak
  memory. It can be faster or slower than `Analytical` depending on backend,
  model, system, precision, and hardware; both build a dense Hessian of the
  movable atoms.
- **Analytical Hessian with UMA workers.** UMA, ORB, MACE, and AIMNet2 support
  analytical Hessians with one calculator. With UMA, `--uma-workers` above 1
  plus `Analytical` raises `BackendError` instead of changing the method; use
  `--uma-workers 1` or finite differences.
- **Charge and spin.** `-q` / `-m` select the electronic state of the
  Hessian. A wrong state invalidates the frequencies, ZPE, and thermochemistry.

## Next step

- [tsopt.md](tsopt.md), [irc.md](irc.md): the usual stages before and after.
- [UMA notes](../pdb2reaction-install/backends.md#uma): `--hessian-calc-mode`.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.FREQ_KW, d.THERMO_KW, d.FREQ_CALC_KW)`
