# YAML Schema

## Top-level Keys

| Key |
|---|
| `extract` |
| `calc` |
| `path_search` |
| `scan` |
| `tsopt` |
| `freq` |
| `dft` |

## Starter Template

```yaml
# Starter config for `pdb2reaction all`

extract:
  radius: 2.6
  radius_het2het: 0.0

calc:
  workers: 1
  workers_per_node: 1

path_search:
  mep_mode: gsm
  max_nodes: 20
  max_cycles: 300

scan:
  max_step_size: 0.2
  bias_k: 300.0
  relax_max_cycles: 10000

tsopt:
  max_cycles: 10000

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value

thermo:
  temperature: 298.15
  pressure_atm: 1.0

dft:
  func_basis: wb97m-v/def2-tzvpd
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
```

## Scalar Defaults

| Key | Type | Default |
|---|---|---|
| `extract.radius` | `float` | `2.6` |
| `extract.radius_het2het` | `float` | `0.0` |
| `calc.workers` | `int` | `1` |
| `calc.workers_per_node` | `int` | `1` |
| `path_search.mep_mode` | `str` | `'gsm'` |
| `path_search.max_nodes` | `int` | `20` |
| `path_search.max_cycles` | `int` | `300` |
| `scan.max_step_size` | `float` | `0.2` |
| `scan.bias_k` | `float` | `300.0` |
| `scan.relax_max_cycles` | `int` | `10000` |
| `tsopt.max_cycles` | `int` | `10000` |
| `freq.max_write` | `int` | `10` |
| `freq.amplitude_ang` | `float` | `0.8` |
| `freq.n_frames` | `int` | `20` |
| `freq.sort` | `str` | `'value'` |
| `thermo.temperature` | `float` | `298.15` |
| `thermo.pressure_atm` | `float` | `1.0` |
| `dft.func_basis` | `str` | `'wb97m-v/def2-tzvpd'` |
| `dft.max_cycle` | `int` | `100` |
| `dft.conv_tol` | `float` | `1e-09` |
| `dft.grid_level` | `int` | `3` |

## Scan Spec Shapes

Accepted by `scan`, `scan2d`, and `scan3d` with `-s/--scan-lists`.

```yaml
# scan (1D staged)
one_based: false
stages:
  - - [1, 2, 1.65]
  - - [2, 3, 2.30]

# scan2d / scan3d
one_based: false
pairs:
  - [1, 2, 1.40, 2.20]
  - [2, 3, 1.20, 2.00]
  - [3, 4, 1.00, 1.80]  # required only for scan3d
```
