# YAML Schema

## Top-level Keys

| Key |
|---|
| `geom` |
| `calc` |
| `opt` |
| `lbfgs` |
| `rfo` |
| `gs` |
| `stopt` |
| `search` |
| `dmf` |
| `hessian_dimer` |
| `rsirfo` |
| `irc` |
| `freq` |
| `thermo` |
| `bias` |
| `bond` |
| `dft` |

## Starter Template

```yaml
# Starter config for `pdb2reaction all`

calc:
  workers: 1
  workers_per_node: 1

gs:
  max_nodes: 20

stopt:
  max_cycles: 300

search:
  max_depth: 10

bias:
  k: 300.0

hessian_dimer:
  thresh: baker

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value

thermo:
  temperature: 298.15
  pressure_atm: 1.0

dft:
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
```

## Scalar Defaults

| Key | Type | Default |
|---|---|---|
| `calc.workers` | `int` | `1` |
| `calc.workers_per_node` | `int` | `1` |
| `gs.max_nodes` | `int` | `20` |
| `stopt.max_cycles` | `int` | `300` |
| `search.max_depth` | `int` | `10` |
| `bias.k` | `float` | `300.0` |
| `hessian_dimer.thresh` | `str` | `'baker'` |
| `freq.max_write` | `int` | `10` |
| `freq.amplitude_ang` | `float` | `0.8` |
| `freq.n_frames` | `int` | `20` |
| `freq.sort` | `str` | `'value'` |
| `thermo.temperature` | `float` | `298.15` |
| `thermo.pressure_atm` | `float` | `1.0` |
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
