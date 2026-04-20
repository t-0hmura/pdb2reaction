# `freq`

## Overview

> **Summary:** Compute vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) with an MLIP backend (UMA by default; `-b/--backend` also supports ORB, MACE, AIMNet2). When VRAM permits, `--hessian-calc-mode Analytical` speeds Hessian evaluation. Imaginary frequencies appear as negative values.

### At a glance
- **Use when:** You need full vibrational analysis (e.g., confirm a stationary point is a true minimum with no imaginary frequencies) and/or compute thermochemistry corrections from an MLIP backend (UMA by default). Note: `tsopt` already includes an imaginary-frequency check, so a separate `freq` run is mainly for thermochemistry or detailed vibrational mode inspection.
- **Frozen atoms:** Supported via PHVA (Partial Hessian Vibrational Analysis).
- **Outputs:** `frequencies_cm-1.txt`, per-mode `_trj.xyz` animations (and optional `.pdb`), plus `thermoanalysis.yaml` when enabled/available.
- **TS check:** A properly converged first-order saddle point (TS) is expected to have **exactly one** imaginary frequency (negative cm⁻¹ value). Note: the 5 cm⁻¹ internal detection threshold and the ~100 cm⁻¹ TS-quality gate answer different questions and must not be conflated — see {ref}`imaginary-mode-thresholds` for the canonical definition.
- **Performance:** For Hessian evaluation modes, see {ref}`hessian-evaluation`.

`pdb2reaction freq` performs vibrational analysis with an MLIP backend (UMA by default), honoring frozen atoms via PHVA. It exports normal-mode animations as `_trj.xyz` (and `.pdb` when a PDB template is available and conversion is enabled), and prints a Gaussian-style thermochemistry summary when the optional `thermoanalysis` package is installed.


## Minimal example

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

## Output checklist

- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz`
- `result_freq/mode_*.pdb` (for PDB inputs with conversion enabled)

## Common examples

1. Limit exported modes for quick inspection.

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

2. Run PHVA with link freezing and dump thermo payload.

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --freeze-links --dump --out-dir ./result_freq_phva
```

3. Use analytical Hessian mode on VRAM-rich nodes.

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## Usage
```bash
pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] \
 [--freeze-links/--no-freeze-links] \
 [--max-write N] [--amplitude-ang Å] [--n-frames N] [--sort value|abs] \
 [--out-dir DIR] [--config FILE] [--show-config] [--dry-run] \
 [--temperature K] [--pressure FLOAT] [--dump/--no-dump] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### Examples
```bash
# Minimal run with explicit charge and spin
pdb2reaction freq -i a.pdb -q 0 -m 1

# PHVA with YAML overrides and a custom output directory
pdb2reaction freq -i a.xyz -q -1 --config ./freq.yaml --out-dir ./result_freq/
```

## Workflow
- **Geometry loading & freeze handling**: structures are read via
  `pysisyphus.helpers.geom_loader`. For PDB inputs, `--freeze-links` detects link
  hydrogens and freezes their parent atoms, then merges the resulting indices with
  `geom.freeze_atoms`; the merged list is echoed and propagated to the MLIP backend and PHVA.
- **MLIP backend**: `--hessian-calc-mode` selects analytical or finite-difference Hessians.
  The MLIP backend may return a partial (active) Hessian block whenever atoms are frozen.
  For Hessian evaluation modes, see {ref}`hessian-evaluation`.
- **PHVA & TR projection**: with frozen atoms, eigenanalysis occurs inside the active
  subspace with translation/rotation modes projected there. Both 3N×3N and active-block
  Hessians are accepted, and frequencies are reported in cm⁻¹ (negatives = imaginary).
- **Mode export**: `--max-write` limits how many modes are animated. Modes are sorted by
  value (or absolute value with `--sort abs`). The sinusoidal animation amplitude
  (`--amplitude-ang`) and frame count (`--n-frames`) match the YAML defaults. `_trj.xyz`
  animations are produced for every input; `.pdb` animations are written only when a PDB
  template exists **and** `--convert-files` remains enabled (ASE conversion is used as a
  fallback).
- **Thermochemistry**: if `thermoanalysis` is installed, a QRRHO-like summary (EE, ZPE, E/H/G
  corrections, heat capacities, entropies) is printed using PHVA frequencies. CLI pressure in
  atm is converted internally to Pa. When `--dump`, a `thermoanalysis.yaml` snapshot is
  also written.
- **Performance & exit behavior**: the implementation minimizes GPU memory usage by keeping
  a single Hessian resident.
  Keyboard interrupts exit with code 130; other failures print a traceback and exit with code 1.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. When omitted, charge can be inferred from `--ligand-charge`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge` supplies it |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | PDB-only. Freeze parents of link hydrogens and merge with `geom.freeze_atoms`. See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-write INT` | Number of modes to export. | `10` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `-o, --out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). On the CLI this flag is `--pressure`; the matching YAML key under `thermo:` is `pressure_atm` (explicit unit suffix). Both are in atm and get converted to Pa internally. | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. Standalone `freq` defaults to `False`; when invoked as part of `pdb2reaction all --thermo` the wrapper flips this to `True` unless you pass `--no-dump`. | `False` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB companions when a PDB template is available (GJF is not written). | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running frequency analysis. | `False` |

## Outputs
```
out_dir/ (default:./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz # Per-mode animations
├─ mode_XXXX_±freqcm-1.pdb # Only when a PDB template exists and conversion is enabled
├─ frequencies_cm-1.txt # Full frequency list using the selected sort order
└─ thermoanalysis.yaml # Present when `thermoanalysis` is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 130 | Keyboard interrupt |
| 1 | Unexpected error |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Imaginary frequencies are reported as negative values in cm⁻¹. `freq` prints how many were detected
  and dumps details when `--dump`.
- `--hessian-calc-mode` follows the standard precedence (defaults < config < explicit CLI); an explicit CLI `--hessian-calc-mode` value takes precedence over `calc.hessian_calc_mode` in the config YAML.

Provide mappings with merge order **defaults < config < explicit CLI**.
Shared sections reuse [YAML Reference](yaml-reference.md).
An additional `thermo` section is supported for thermochemistry controls.

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 task_name: omol # UMA task name
 device: auto # MLIP device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true # allow partial Hessians
 backend: uma # MLIP backend: uma, orb, mace, aimnet2
 solvent: none # implicit solvent name (e.g. water) or none
 solvent_model: alpb # xTB solvent model: alpb or cpcmx
freq:
 amplitude_ang: 0.8 # displacement amplitude for modes (Å)
 n_frames: 20 # number of frames per mode
 max_write: 10 # maximum number of modes to write
 sort: value # sort order: value vs abs
thermo:
 temperature: 298.15 # thermochemistry temperature (K)
 pressure_atm: 1.0 # thermochemistry pressure (atm)
 dump: false # write thermoanalysis.yaml when true
```

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing

- [tsopt](tsopt.md) — Optimize TS candidates (includes imaginary-frequency check; follow with IRC for endpoint validation)
- [irc](irc.md) — IRC from TS (freq is often run on IRC endpoints for thermochemistry)
- [dft](dft.md) — Single-point DFT for higher-level energy refinement
- [all](all.md) — End-to-end workflow with `--thermo`
- [YAML Reference](yaml-reference.md) — Full `freq` and `thermo` configuration options
- [Glossary](glossary.md) — Definitions of ZPE, Gibbs Energy, Enthalpy, Entropy
