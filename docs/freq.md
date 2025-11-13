# `freq` subcommand

## Purpose
Computes vibrational frequencies using UMA, performs partial Hessian vibrational analysis (PHVA) when atoms are frozen, exports animated modes, and prints a thermochemistry summary.

## Usage
```bash
pdb2reaction freq -i INPUT -q CHARGE [--spin 2S+1]
                  [--freeze-links/--no-freeze-links]
                  [--max-write N] [--amplitude-ang Å] [--n-frames N]
                  [--sort value|abs] [--out-dir DIR]
                  [--temperature K] [--pressure atm] [--dump/--no-dump]
                  [--hessian-calc-mode Analytical|FiniteDifference]
                  [--args-yaml FILE]
```

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. | Required |
| `-s, --spin INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-links / --no-freeze-links` | For PDB inputs, freeze link-hydrogen parents (merged with `geom.freeze_atoms`). | `--freeze-links` |
| `--max-write INT` | Number of modes to export. | `20` |
| `--amplitude-ang FLOAT` | Animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `--out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump / --no-dump` | Write `thermoanalysis.yaml`. | `--no-dump` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ |
| `--args-yaml FILE` | YAML overrides (see below). | _None_ |

## YAML configuration (`--args-yaml`)
Accepts a mapping; CLI overrides YAML. Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

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

### Shared sections
- `geom`, `calc`: same keys as [`opt`](opt.md#yaml-configuration-args-yaml). `--freeze-links` augments `geom.freeze_atoms`, which are also forwarded to UMA for Hessian calculations.

### Section `freq`
Controls mode export.

- `amplitude_ang` (`0.8`): Animation amplitude (Å).
- `n_frames` (`20`): Frames per vibrational animation.
- `max_write` (`20`): Number of modes to export.
- `sort` (`"value"`): Ordering (`"value"` ascending or `"abs"` for absolute value).

_The thermochemistry parameters (`temperature`, `pressure_atm`, `dump`) are currently CLI-only and not read from YAML._

## Outputs
- `<out-dir>/mode_XXXX_±freqcm-1.(trj|pdb)` animations for exported modes.
- `<out-dir>/frequencies_cm-1.txt` sorted list of frequencies.
- Optional `<out-dir>/thermoanalysis.yaml` when `--dump` is enabled.
- Console blocks summarising resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Notes
- Imaginary modes are reported as negative frequencies; PHVA restricts the Hessian to active degrees of freedom when atoms are frozen.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode`; Analytical mode may require more GPU memory than finite differences.
- Mode animations use sinusoidal displacements; PDB animations employ MODEL/ENDMDL records for multi-model files.
