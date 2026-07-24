# `pdb2reaction freq`

## Purpose

Vibrational analysis: build the Hessian, diagonalize for normal-mode
frequencies, write per-mode geometry displacements, and compute
QRRHO thermochemistry. Default temperature 298.15 K, 1 atm.
Partial-Hessian variant (PHVA) activates automatically when
`freeze_atoms` is non-empty.

## Synopsis

```bash
pdb2reaction freq -i geom.{pdb,cif,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--temperature 298.15] [--pressure 1.0] [--symmetry-number 1] \
    [-b uma|orb|mace|aimnet2] [-o ./result_freq/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--temperature` | float | 298.15 | K, for thermochemistry |
| `--pressure` | float | 1.0 | atm, for thermochemistry |
| `--symmetry-number` | int ≥ 1 | 1 | External rotational symmetry number. Point-group symmetry is not inferred. |
| `--hessian-calc-mode` | str | `FiniteDifference` | `Analytical` / `FiniteDifference`; check `UMA_CALC_KW` |
| `--tr-projection` | str | `constrained` | PHVA rigid-mode treatment. `legacy-active` is deprecated comparison-only behavior; never use it for pass/HOSP transition-state certification. |
| `--workers`, `--workers-per-node` | int | `1`, `1` | UMA predictor workers. An explicit `Analytical` request with `workers > 1` raises `BackendError`; use one worker or finite differences. Other built-in backends ignore these worker kwargs. |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_freq/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default (298.15 K, 1 atm)

```bash
pdb2reaction freq -i ts.xyz -q 0 -m 1 -b uma --out-json -o result_freq
```

`--out-json` enables the `result.json` shown below; without it you still get the text / trajectory outputs (just no `result.json`).

### Higher temperature for activation enthalpy

```bash
pdb2reaction freq -i ts.pdb -l 'SAM:1' \
    --temperature 310.15 --pressure 1.0 \
    -b uma -o result_freq_310K
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/frequencies_cm-1.txt` | successful Hessian/frequency evaluation | all modes, sorted (cm⁻¹) |
| `<out_dir>/thermoanalysis.yaml` | `--dump` | ZPE, thermal corrections, S, H, G |
| `<out_dir>/mode_NNNN_<freq>cm-1_trj.xyz` | up to `--max-write` selected modes | displacement trajectory (top-level, NOT under `vib/`); none are written when `--max-write 0` |
| `<out_dir>/mode_NNNN_<freq>cm-1.pdb` | `--convert-files` and PDB/mmCIF topology/reference available | PDB companion |
| `<out_dir>/mode_NNNN_<freq>cm-1.cif` | `--convert-files` and input/reference required the mmCIF or oversized-PDB bridge | public companion with original IDs |

`result.json` keys:

```python
import json
d = json.load(open("result_freq/result.json"))
print(d["n_imaginary"])                       # count of every frequency < 0 cm-1
print(d["rigid_projection"]["treatment"], d["rigid_projection"]["effective_rank"])
print(d["frequencies_cm"][:5])                # first five frequencies (cm-1)
t = d["thermochemistry"]
print(t["electronic_energy_ha"])              # EE (Hartree)
print(t["zpe_correction_ha"])                 # ZPE correction (Hartree)
print(t["thermal_correction_free_energy_ha"]) # dG_therm (Hartree)
print(t["sum_EE_and_thermal_free_energy_ha"]) # EE + dG_therm (Hartree)
print(t["S_cal_per_mol_K"])                   # entropy (cal/mol·K)
print(t["symmetry_number"], t["symmetry_number_source"])
```

`thermoanalysis.yaml` contains the same `rigid_projection` provenance when it is written.

## QRRHO thermochemistry

Default thermochemistry uses the QRRHO (Grimme) treatment with a
100 cm⁻¹ rotor cutoff:

- low-frequency vibrations (< 100 cm⁻¹) are interpolated toward the
  free-rotor entropy formula,
- high-frequency vibrations use the standard harmonic-oscillator
  partition function.

`THERMO_KW` (`pdb2reaction.core.defaults`) exposes `temperature`,
`pressure_atm`, `symmetry_number`, and `dump`. Supply the external rotational
symmetry number explicitly when it is not 1; the workflow does not infer a
point group. The QRRHO rotor cutoff is internal to `thermoanalysis` and is not
a `THERMO_KW` knob.

## Partial-Hessian Vibrational Analysis (PHVA)

When the input has frozen atoms, `freq` automatically computes the
**partial Hessian**: only the mobile-atom block is built and
diagonalized; frozen atoms are projected out. This reduces the active dense
block, with the actual time/memory benefit set by the number of mobile degrees
of freedom and backend implementation.

The default `--tr-projection constrained` removes only full-system rigid
motions that leave every frozen anchor fixed. A normal multi-anchor cluster
boundary therefore usually has effective rank 0. `legacy-active` is an
isolated-active comparison treatment, not a bitwise replay guarantee for
near-linear or degenerate structures. It is deprecated and must not be used
for pass/HOSP transition-state certification. See `freeze-atoms.md`.

`pdb2reaction` does **not** read PDB B-factors as a freeze list. The
freeze set is assembled as a union of CLI `--freeze-atoms`
(1-based atom indices), CLI `--freeze-links/--no-freeze-links` (auto
freeze of `LKH/HL` cap-H parents written by `extract`), and YAML
`geom.freeze_atoms`. See `freeze-atoms.md`.

## Caveats

- A minimum should have **0 imaginary frequencies**, a TS should have
  **exactly 1**.
- An all-frozen structure has no active vibrational DOF and raises an error.
- `freq` reports every strictly negative frequency in `n_imaginary`.
  `tsopt` deliberately applies its configured small-negative threshold
  (5 cm⁻¹ by default), so the two counters can differ for near-zero numerical
  modes. Inspect `frequencies_cm`, not just the integer.
- A small-magnitude imaginary frequency may be numerical or a real shallow
  mode; inspect its displacement and repeat the Hessian at suitable precision.
  The 100 cm⁻¹ QRRHO cutoff regularizes low-frequency thermochemistry but
  does **not** turn an imaginary mode into a real one or validate a stationary point.
- `--hessian-calc-mode FiniteDifference` usually lowers peak model/autograd
  memory. It can be faster or slower than `Analytical` depending on backend,
  model, system, precision, and hardware; both paths still materialize a dense
  active-space Hessian.
- UMA's multi-worker predictor has no analytical autograd model. Therefore
  `--workers > 1 --hessian-calc-mode Analytical` is rejected rather than
  silently changing the requested numerical method. UMA, ORB, MACE, and
  AIMNet2 all support explicit analytical Hessians in their supported
  single-calculator configurations.
- Charge/spin select the MLIP electronic state from which the Hessian is built.
  Verify `-q`/`-m`; an incorrect state invalidates the frequencies and derived
  ZPE/thermochemistry.

## See also

- `tsopt.md`, `irc.md` — usual upstream stages.
- [`pdb2reaction-install-backends/uma.md`](../pdb2reaction-install-backends/uma.md) — `--hessian-calc-mode` knob.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.FREQ_KW, d.THERMO_KW, d.FREQ_CALC_KW)`
