# `freq` subcommand

## Overview
`pdb2reaction freq` performs vibrational analysis with the UMA calculator, honoring any
frozen atoms via partial Hessian vibrational analysis (PHVA). It exports mass-weighted
normal modes as `.trj`/`.pdb` animations, prints a Gaussian-style thermochemistry summary
when the optional `thermoanalysis` package is installed, and can emit a YAML summary when
`--dump True`. Configuration starts from defaults, applies CLI switches, and finally
applies YAML overrides (`geom`, `calc`, `freq`) with highest precedence, so the same
template can drive both standalone runs and workflows launched by other subcommands.

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
  animations are produced for every input; `.pdb` animations are written only when a PDB
  template exists **and** `--convert-files` remains enabled (ASE conversion is used as a
  fallback).
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
| `-q, --charge INT` | Total charge. When omitted, charge can be inferred from `--ligand-charge`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge` supplies it |
| `--ligand-charge TEXT` | Total charge or per-resname mapping used when `-q` is omitted. Triggers extract-style charge derivation on the full complex. | `None` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links BOOL` | PDB-only. Freeze parents of link hydrogens and merge with `geom.freeze_atoms`. | `True` |
| `--max-write INT` | Number of modes to export. | `10` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `--out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump BOOL` | Explicit `True`/`False`. Write `thermoanalysis.yaml`. | `False` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _None_ (uses YAML/default of `FiniteDifference`) |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB companions when a PDB template is available (GJF is not written). | `--convert-files` |
| `--args-yaml FILE` | YAML overrides (sections: `geom`, `calc`, `freq`). | _None_ |

## Outputs
```
out_dir/ (default: ./result_freq/)
├─ mode_XXXX_±freqcm-1.trj  # Per-mode animations
├─ mode_XXXX_±freqcm-1.pdb  # Only when a PDB template exists and conversion is enabled
├─ frequencies_cm-1.txt     # Full frequency list using the selected sort order
└─ thermoanalysis.yaml      # Present when `thermoanalysis` is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Notes
- Imaginary modes are reported as negative frequencies. `freq` prints how many were detected
  and dumps details when `--dump True`.
- `--hessian-calc-mode` overrides `calc.hessian_calc_mode` after YAML merging.
- Charge/spin inherit `.gjf` metadata when available. If `-q` is omitted but
  `--ligand-charge` is provided, the input is treated as an enzyme–substrate
  complex and `extract.py`’s charge summary computes the total charge; explicit
  `-q` still overrides. Otherwise charge defaults to `0` and multiplicity to `1`.
  Override them explicitly to ensure the intended state.

## YAML configuration (`--args-yaml`)
Provide a mapping; YAML values override both defaults and CLI switches (highest
precedence). Shared sections reuse [`opt`](opt.md#yaml-configuration-args-yaml).

```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  task_name: omol            # UMA task name
  device: auto               # UMA device selection
  max_neigh: null            # maximum neighbors for graph construction
  radius: null               # cutoff radius for neighbor search
  r_edges: false             # store radial edges
  out_hess_torch: true       # request torch-form Hessian
  freeze_atoms: null         # calculator-level frozen atoms
  hessian_calc_mode: FiniteDifference   # Hessian mode selection
  return_partial_hessian: true          # allow partial Hessians
freq:
  amplitude_ang: 0.8         # displacement amplitude for modes (Å)
  n_frames: 20               # number of frames per mode
  max_write: 10              # maximum number of modes to write
  sort: value                # sort order: value vs abs
```
