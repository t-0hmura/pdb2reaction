# `pdb2reaction freq`

## Purpose

Vibrational analysis: build the Hessian, diagonalize for normal-mode
frequencies, write per-mode geometry displacements, and compute
QRRHO thermochemistry. Default temperature 298.15 K, 1 atm.
Partial-Hessian variant (PHVA) activates automatically when
`freeze_atoms` is non-empty.

## Synopsis

```bash
pdb2reaction freq -i geom.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--temperature 298.15] [--pressure 1.0] \
    [-b uma|orb|mace|aimnet2] [-o ./result_freq/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--temperature` | float | 298.15 | K, for thermochemistry |
| `--pressure` | float | 1.0 | atm, for thermochemistry |
| `--hessian-calc-mode` | str | `FiniteDifference` | `Analytical` / `FiniteDifference`; check `UMA_CALC_KW` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_freq/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default (298.15 K, 1 atm)

```bash
pdb2reaction freq -i ts.xyz -q 0 -m 1 -b uma --out-json -o result_freq
```

`--out-json` enables the `result.json` shown below; omit it for `frequencies_cm-1.txt` only.

### Higher temperature for activation enthalpy

```bash
pdb2reaction freq -i ts.xyz -l 'SAM:1' \
    --temperature 310.15 --pressure 1.0 \
    -b uma -o result_freq_310K
```

## Output

| Path | When | Content |
|---|---|---|
| `<out_dir>/result.json` | `--out-json` | machine-readable result |
| `<out_dir>/frequencies_cm-1.txt` | always | all modes, sorted (cm⁻¹) |
| `<out_dir>/thermoanalysis.yaml` | `--dump` | ZPE, thermal corrections, S, H, G |
| `<out_dir>/mode_NNNN_<freq>cm-1_trj.xyz` | always | per-mode displacement trajectory (top-level, NOT under `vib/`) |
| `<out_dir>/mode_NNNN_<freq>cm-1.pdb` | input has PDB topology | PDB companion |

`result.json` keys:

```python
import json
d = json.load(open("result_freq/result.json"))
print(d["n_imaginary"])                       # 0 for minimum, 1 for TS
print(d["frequencies_cm"][:5])                # first five frequencies (cm-1)
t = d["thermochemistry"]
print(t["electronic_energy_ha"])              # EE (Hartree)
print(t["zpe_correction_ha"])                 # ZPE correction (Hartree)
print(t["thermal_correction_free_energy_ha"]) # dG_therm (Hartree)
print(t["sum_EE_and_thermal_free_energy_ha"]) # EE + dG_therm (Hartree)
print(t["S_cal_per_mol_K"])                   # entropy (cal/mol·K)
```

## QRRHO thermochemistry

Default thermochemistry uses the QRRHO (Grimme) treatment with a
100 cm⁻¹ rotor cutoff:

- low-frequency vibrations (< 100 cm⁻¹) are interpolated toward the
  free-rotor entropy formula,
- high-frequency vibrations use the standard harmonic-oscillator
  partition function.

Inspect the active QRRHO knob via `pdb2reaction.core.defaults.THERMO_KW`.

## Partial-Hessian Vibrational Analysis (PHVA)

When the input has frozen atoms, `freq` automatically computes the
**partial Hessian**: only the mobile-atom block is built and
diagonalized; frozen atoms are projected out. This is much cheaper for
large clusters.

`pdb2reaction` does **not** read PDB B-factors as a freeze list. The
freeze set is assembled (in priority order) from CLI `--freeze-atoms`
(1-based atom indices), CLI `--freeze-links/--no-freeze-links` (auto
freeze of `LKH/HL` link-H parents written by `extract`), and YAML
`geom.freeze_atoms`. See `freeze-atoms.md`.

## Caveats

- A minimum should have **0 imaginary frequencies**, a TS should have
  **exactly 1**.
- Imaginary frequencies < ~50 cm⁻¹ are often numerical noise, not real
  modes. The QRRHO cutoff (100 cm⁻¹) is one safeguard.
- `--hessian-calc-mode FiniteDifference` is more memory-friendly but
  slower than `Analytical`.
- Thermochemistry depends on charge / spin — make sure `-q`/`-m` are
  correct or ZPE will be off.

## See also

- `tsopt.md`, `irc.md` — usual upstream stages.
- [`pdb2reaction-install-backends/uma.md`](../pdb2reaction-install-backends/uma.md) — `--hessian-calc-mode` knob.
- Defaults: `import pdb2reaction.core.defaults as d; print(d.FREQ_KW, d.THERMO_KW, d.FREQ_CALC_KW)`
