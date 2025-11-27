# `freq` subcommand

## Overview
`pdb2reaction freq` performs vibrational analysis with the UMA calculator, honoring any
frozen atoms via partial Hessian vibrational analysis (PHVA). It exports mass-weighted
normal modes as `.trj`/`.pdb` animations, prints a Gaussian-style thermochemistry summary
when the optional `thermoanalysis` package is installed, and can emit a YAML summary when
`--dump True`. Configuration values are read from YAML (`geom`, `calc`, `freq`) and then
overridden by CLI switches, so the same template can drive both standalone runs and
workflows launched by other subcommands.

## Usage
```bash
pdb2reaction freq -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-m 2S+1] \
                  [--freeze-links {True|False}] \
                  [--max-write N] [--amplitude-ang Å] [--n-frames N] \
                  [--sort value|abs] [--out-dir DIR] [--args-yaml FILE] \
                  [--temperature K] [--pressure atm] [--dump {True|False}] \
                  [--hessian-calc-mode Analytical|FiniteDifference] \
                  [--convert-files/--no-convert-files]
```

### Examples
```bash
# Minimal run with explicit charge and spin
pdb2reaction freq -i a.pdb -q 0 -m 1

# PHVA with YAML overrides and a custom output directory
pdb2reaction freq -i a.xyz -q -1 --args-yaml ./args.yaml --out-dir ./result_freq/
```

## Workflow
- **Geometry loading & freeze handling**: structures are read via
  `pysisyphus.helpers.geom_loader`. For PDB inputs, `--freeze-links True` detects link
  hydrogens and freezes their parent atoms, then merges the resulting indices with
  `geom.freeze_atoms`; the merged list is echoed and propagated to UMA and PHVA.
- **UMA calculator**: `--hessian-calc-mode` selects analytical or finite-difference Hessians.
  UMA may return a partial (active) Hessian block whenever atoms are frozen.
- **PHVA & TR projection**: with frozen atoms, eigenanalysis occurs inside the active
  subspace with translation/rotation modes projected there. Both 3N×3N and active-block
  Hessians are accepted, and frequencies are reported in cm⁻¹ (negatives = imaginary).
- **Mode export**: `--max-write` limits how many modes are animated. Modes are sorted by
  value (or absolute value with `--sort abs`). The sinusoidal animation amplitude
  (`--amplitude-ang`) and frame count (`--n-frames`) match the YAML defaults. `.trj`
  animations are produced for every input; `.pdb` animations mirror them when a PDB template
  is available (ASE conversion is used as a fallback).
- **Thermochemistry**: if `thermoanalysis` is installed, a QRRHO-like summary (EE, ZPE, E/H/G
  corrections, heat capacities, entropies) is printed using PHVA frequencies. CLI pressure in
  atm is converted internally to Pa. When `--dump True`, a `thermoanalysis.yaml` snapshot is
  also written.
- **Performance & exit behavior**: the implementation minimizes GPU memory usage by keeping
  a single Hessian resident, preferring upper-triangular eigendecompositions (`UPLO="U"`).
  Keyboard interrupts exit with code 130; other failures print a traceback and exit with code 1.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Required unless the input is a `.gjf` template with charge metadata. | Required when not in template |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | PDB-only. Freeze parents of link hydrogens and merge with `geom.freeze_atoms`. | `True` |
| `--max-write INT` | Number of modes to export. | `20` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `--out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump BOOL` | Explicit `True`/`False`. Write `thermoanalysis.yaml`. | `False` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ (use YAML/default) |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `--convert-files` |
| `--args-yaml FILE` | YAML overrides (sections: `geom`, `calc`, `freq`). | _None_ |

## Outputs
- `<out-dir>/mode_XXXX_±freqcm-1.trj` and `.pdb` animations for each written mode (mirroring follows `--convert-files`).
- `<out-dir>/frequencies_cm-1.txt` listing every computed frequency according to the sort
  order.
- `<out-dir>/thermoanalysis.yaml` when both `thermoanalysis` is importable and `--dump True`.
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Notes
- Imaginary modes are reported as negative frequencies. `freq` prints how many were detected
  and dumps details when `--dump True`.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after YAML merging.
- Charge/spin inherit `.gjf` metadata when available; otherwise `-q/--charge` is required
  and multiplicity defaults to `1`. Override them explicitly to ensure the intended state.

## YAML configuration (`--args-yaml`)
Provide a mapping; CLI values override YAML. Shared sections reuse
[`opt`](opt.md#yaml-configuration-args-yaml).

```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
  model: uma-s-1p1
  task_name: omol
  device: auto
  max_neigh: null
  radius: null
  r_edges: false
  out_hess_torch: true
  freeze_atoms: null
  hessian_calc_mode: Analytical
  return_partial_hessian: true
freq:
  amplitude_ang: 0.8
  n_frames: 20
  max_write: 20
  sort: value
```
