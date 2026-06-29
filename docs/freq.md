# `freq`

Compute vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) using an MLIP backend (UMA by default; `-b/--backend` also supports ORB, MACE, AIMNet2). Use it when full vibrational analysis is required — for example, to confirm that a stationary point is a true minimum with no imaginary frequencies, or that a TS has exactly one — or when thermochemistry corrections (ZPE, Gibbs energy) are needed. When VRAM permits, `--hessian-calc-mode Analytical` speeds up Hessian evaluation. Imaginary frequencies appear as negative values.

## Examples

Minimal run with explicit charge and spin:

```bash
# Minimal run with explicit charge and spin
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

PHVA with link freezing and dump thermo payload:

```bash
# PHVA with link freezing and dump thermo payload
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --freeze-links --dump --out-dir ./result_freq_phva
```

Analytical Hessian mode on VRAM-rich nodes:

```bash
# Analytical Hessian mode on VRAM-rich nodes
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## Workflow

- **Geometry loading & freeze handling**: structures are read via
  `pysisyphus.helpers.geom_loader`. For PDB inputs, `--freeze-links` detects cap
  hydrogens and freezes their parent atoms, then merges the resulting indices with
  `geom.freeze_atoms`; the merged list is echoed and propagated to the MLIP backend and PHVA.
- **MLIP backend**: `--hessian-calc-mode` selects analytical or finite-difference Hessians.
  The MLIP backend may return a partial (active) Hessian block whenever atoms are frozen.
  For Hessian evaluation modes, see {ref}`hessian-evaluation`.
- **PHVA & TR projection**: with frozen atoms, eigenanalysis occurs inside the active
  subspace, where translation/rotation modes are projected out. Both 3N×3N and active-block
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
- **Performance**: the implementation minimizes GPU memory usage by keeping a single Hessian resident.

## Outputs

```text
out_dir/ (default:./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz # Per-mode animations
├─ mode_XXXX_±freqcm-1.pdb # Only when a PDB template exists and conversion is enabled
├─ frequencies_cm-1.txt # Full frequency list using the selected sort order
└─ thermoanalysis.yaml # Present when `thermoanalysis` is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Exit codes

See {ref}`exit-codes` in CLI Conventions.

## CLI options

The tables below cover the options that need explanation; the full flag list is in the generated [command reference](reference/commands/index.md) — do not hand-duplicate it here.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader` (`.pdb` / `.xyz` / `.trj` / ...). | Required |
| `-q, --charge INT` | Total charge. When omitted, charge can be inferred from `--ligand-charge/-l`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge/-l` supplies it |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | PDB input (or XYZ/GJF with `--ref-pdb`). Freeze parents of cap hydrogens and merge with `geom.freeze_atoms`. See [extract](extract.md) for cap-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-write INT` | Number of modes to export. | `10` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `-o, --out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). On the CLI this flag is `--pressure`; the matching YAML key under `thermo:` is `pressure_atm` (explicit unit suffix). Both are in atm and get converted to Pa internally. | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. Standalone `freq` defaults to `False` (use `--no-dump` to keep it off explicitly); when invoked as part of `pdb2reaction all --thermo` the freq stage uses `dump=True` by default — pass `--dump False` to `all` to disable it. | `False` |
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

## YAML configuration

The `geom`, `calc`, `freq`, and `thermo` sections are unchanged from the canonical definitions in [YAML Reference](yaml-reference.md): see [`geom`](yaml-reference.md#geom), [`calc`](yaml-reference.md#calc), [`freq`](yaml-reference.md#freq-section), and [`thermo`](yaml-reference.md#thermo). `freq` automatically sets `calc.return_partial_hessian = true` (PHVA) by default; YAML can override.

The only `freq`-specific default that differs from the canonical block is the output directory:

```yaml
freq:
 out_dir: ./result_freq/ # freq default
```

## Notes

- `tsopt` already includes an imaginary-frequency check, so a separate `freq` run is mainly for thermochemistry or detailed mode inspection.
- A properly converged first-order saddle point (TS) is expected to have **exactly one** imaginary frequency (detection cutoff `hessian_dimer.neg_freq_thresh_cm`, default 5 cm⁻¹).
- Imaginary frequencies are reported as negative values in cm⁻¹. `freq` prints how many were detected
  and dumps details when `--dump`.
- `--hessian-calc-mode` follows the standard precedence (defaults < config < explicit CLI); an explicit CLI `--hessian-calc-mode` value takes precedence over `calc.hessian_calc_mode` in the config YAML.

## See Also

- [tsopt](tsopt.md) — Optimize TS candidates (includes imaginary-frequency check; follow with IRC for endpoint validation)
- [irc](irc.md) — IRC from TS (freq is often run on IRC endpoints for thermochemistry)
- [dft](dft.md) — Single-point DFT for higher-level energy refinement
- [all](all.md) — End-to-end workflow with `--thermo`
- [YAML Reference](yaml-reference.md) — Full `freq` and `thermo` configuration options
- [Glossary](glossary.md) — Definitions of ZPE, Gibbs Energy, Enthalpy, Entropy
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed fixes for common failure modes
