# `freq`

## Overview

> **Summary:** Compute vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) with UMA. When VRAM permits, `--hessian-calc-mode Analytical` speeds Hessian evaluation. Imaginary frequencies appear as negative values.

### At a glance
- **Use when:** You want to validate a minimum/TS candidate and/or compute thermo corrections from UMA.
- **Frozen atoms:** Supported via PHVA (partial Hessian vibrational analysis).
- **Outputs:** `frequencies_cm-1.txt`, per-mode `_trj.xyz` animations (and optional `.pdb`), plus `thermoanalysis.yaml` when enabled/available.
- **TS check:** A properly converged TS is expected to have **exactly one** imaginary frequency (negative cm⁻¹).
- **Performance:** If you have ample VRAM, `--hessian-calc-mode Analytical` is usually recommended.

`pdb2reaction freq` performs vibrational analysis with the UMA calculator, honoring frozen atoms via PHVA. It exports normal-mode animations as `_trj.xyz` (and `.pdb` when a PDB template is available and conversion is enabled), and prints a Gaussian-style thermochemistry summary when the optional `thermoanalysis` package is installed.


## Minimal example

```bash
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

## Output checklist

- `result_freq/summary.md`
- `result_freq/key_frequencies.txt`
- `result_freq/key_mode_1_trj.xyz`
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
pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [--freeze-links/--no-freeze-links] \
 [--max-write N] [--amplitude-ang Å] [--n-frames N] \
 [--show-config] [--dry-run] \
 [--temperature K] [--pressure atm] [--dump/--no-dump] \
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
 `geom.freeze_atoms`; the merged list is echoed and propagated to UMA and PHVA.
- **UMA calculator**: `--hessian-calc-mode` selects analytical or finite-difference Hessians.
 UMA may return a partial (active) Hessian block whenever atoms are frozen.
 When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.
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
 a single Hessian resident, preferring upper-triangular eigendecompositions (`UPLO="U"`).
 Keyboard interrupts exit with code 130; other failures print a traceback and exit with code 1.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. When omitted, charge can be inferred from `--ligand-charge`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge` supplies it |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | PDB-only. Freeze parents of link hydrogens and merge with `geom.freeze_atoms`. | `True` |
| `--max-write INT` | Number of modes to export. | `10` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `--out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. | `False` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB companions when a PDB template is available (GJF is not written). | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running frequency analysis. | `False` |

## Outputs
```
out_dir/ (default:./result_freq/)
├─ summary.md # Quick index of key outputs
├─ key_frequencies.txt # Shortcut to frequencies_cm-1.txt
├─ key_mode_1_trj.xyz # Shortcut to a representative mode trajectory
├─ key_mode_1.pdb # Shortcut to representative mode PDB (when available)
├─ key_thermo.yaml # Shortcut to thermoanalysis.yaml (when available)
├─ mode_XXXX_±freqcm-1_trj.xyz # Per-mode animations
├─ mode_XXXX_±freqcm-1.pdb # Only when a PDB template exists and conversion is enabled
├─ frequencies_cm-1.txt # Full frequency list using the selected sort order
└─ thermoanalysis.yaml # Present when `thermoanalysis` is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Imaginary modes are reported as negative frequencies. `freq` prints how many were detected
 and dumps details when `--dump`.
- `--hessian-calc-mode` follows the standard precedence (defaults < config < explicit CLI < override); if YAML
 specifies `calc.hessian_calc_mode`, it overrides the CLI value.

Provide mappings with merge order **defaults < config < explicit CLI < override**.
Shared sections reuse [YAML Reference](yaml_reference.md).
An additional `thermo` section is supported for thermochemistry controls.

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # UMA model tag
 task_name: omol # UMA task name
 device: auto # UMA device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true # allow partial Hessians
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

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [tsopt](tsopt.md) — Optimize TS candidates (validate with freq/IRC; expected: one imaginary frequency)
- [irc](irc.md) — IRC from TS (often paired with freq on endpoints)
- [dft](dft.md) — Single-point DFT for higher-level energy refinement
- [all](all.md) — End-to-end workflow with `--thermo`
- [YAML Reference](yaml_reference.md) — Full `freq` and `thermo` configuration options
- [Glossary](glossary.md) — Definitions of ZPE, Gibbs Energy, Enthalpy, Entropy
