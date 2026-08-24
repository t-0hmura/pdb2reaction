# `freq`

Compute vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) using an MLIP backend (UMA by default; `-b/--backend` also supports ORB, MACE, AIMNet2). Use it when full vibrational analysis is required — for example, to confirm that a stationary point is a true minimum with no imaginary frequencies, or that a TS has exactly one — or when thermochemistry corrections (ZPE, Gibbs energy) are needed. Finite differences are the default. `--hessian-calc-mode Analytical` avoids displacement error but may be faster or slower and usually needs more accelerator memory; benchmark and validate the selected backend/model on the target system. Imaginary frequencies appear as negative values.

## Examples

Minimal run with explicit charge and spin:

```bash
# Minimal run with explicit charge and spin
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --out-dir ./result_freq
```

PHVA with cap-hydrogen parent freezing and dump thermo payload:

```bash
# PHVA with cap-hydrogen parent freezing and dump thermo payload
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 --freeze-links --dump --out-dir ./result_freq_phva
```

Explicit analytical Hessian mode (after validating memory and runtime):

```bash
# Explicit analytical Hessian mode
pdb2reaction freq -i ts_or_min.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## Workflow

- **Geometry loading & freeze handling**: structures pass through the common
  structure bridge before `pysisyphus.helpers.geom_loader`. For PDB/mmCIF topology inputs, `--freeze-links` detects cap
  hydrogens and freezes their parent atoms, then merges the resulting indices with
  `geom.freeze_atoms`; the merged list is echoed and propagated to the MLIP backend and PHVA.
- **MLIP backend**: `--hessian-calc-mode` selects analytical or finite-difference Hessians.
  The MLIP backend may return a partial (active) Hessian block whenever atoms are frozen.
  For Hessian evaluation modes, see {ref}`hessian-evaluation`.
- **PHVA & rigid modes**: with frozen atoms, eigenanalysis occurs inside the active
  subspace. The default `constrained` treatment removes only full-system rigid motions that
  leave every frozen anchor fixed; a normal multi-anchor cluster boundary therefore usually
  has effective rank 0. Both 3N×3N and active-block Hessians are accepted. See
  [Frozen Atoms](freeze-atoms.md#rigid-modes-with-frozen-boundaries).
- **Mode export**: `--max-write` limits how many modes are animated. Modes are sorted by
  value (or absolute value with `--sort abs`). The sinusoidal animation amplitude
  (`--amplitude-ang`) and frame count (`--n-frames`) match the YAML defaults. `_trj.xyz`
  animations are produced for every input; topology inputs also receive `.pdb` animations
  when `--convert-files` remains enabled, and mmCIF/oversized-PDB bridge inputs additionally
  receive `.cif` animations with the original identifiers.
- **Thermochemistry**: if `thermoanalysis` is installed, a QRRHO-like summary (E, ZPE, E/H/G
  corrections, heat capacities, entropies) is printed using PHVA frequencies. CLI pressure in
  atm is converted internally to Pa. When `--dump`, a `thermoanalysis.yaml` snapshot is
  also written. The console reports the structure energy in Hartree as
  `E + G_corr = G` (electronic energy + Gibbs free-energy correction = Gibbs
  free energy). The molecular point group and external rotational symmetry number are
  detected from each analyzed structure and the resulting `1/sigma` correction is always
  included. An expert can override the detected number with
  `thermo.symmetry_number` in YAML.
- **Frequency-treatment policy**: `freq` applies the **standalone-freq policy** — QRRHO with a
  100 cm⁻¹ rotor cutoff, unit frequency/ZPE scaling, **no** imaginary-frequency inversion, and
  **no** positive-frequency floor. This is deliberately different from the internal
  `Geometry.get_thermoanalysis` policy used by some bundled-engine paths, which additionally
  inverts small imaginaries (from −15 cm⁻¹) and floors positive frequencies below 25 cm⁻¹.
  Neither is a universal scientific default; each is tied to its entry point. The effective
  policy (`kind`, `rotor_cutoff_cm`, `frequency_scale`, `zpe_scale`, `invert_imag_from_cm`,
  `positive_frequency_floor_cm`) is serialized under `thermo_policy` in `thermoanalysis.yaml`
  and in `result.json`.
- **Performance**: the implementation minimizes GPU memory usage by keeping a single Hessian resident.

## Outputs

```text
out_dir/ (default:./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz # Per-mode animations
├─ mode_XXXX_±freqcm-1.pdb # PDB/mmCIF topology exists and conversion is enabled
├─ mode_XXXX_±freqcm-1.cif # mmCIF/oversized-PDB bridge input
├─ frequencies_cm-1.txt # Full frequency list using the selected sort order
└─ thermoanalysis.yaml # Present when `thermoanalysis` is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Exit codes

See {ref}`exit-codes` in CLI Conventions.

## CLI options

The tables below cover the options that need explanation; the full flag list is in the generated [command reference](reference/commands/index.md).

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by the input bridge (`.pdb` / `.cif` / `.mmcif` / `.xyz` / `.trj` / ...). | Required |
| `-q, --charge INT` | Total charge. Explicit `-q` overrides YAML `calc.charge` and `--ligand-charge/-l`; when omitted, YAML, residue derivation, or `.gjf` metadata may supply it. | Required unless YAML/template/derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB/mmCIF residue metadata. Used when `-q` is omitted (PDB/mmCIF inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | UMA predictor parallelism. `workers > 1` cannot be combined with an explicit analytical Hessian request; use `workers = 1` or finite differences. See {ref}`workers-analytical-error`. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Explicit `-m` overrides YAML `calc.spin`; otherwise YAML, `.gjf`, or `1` is used. | YAML/`.gjf`/`1` |
| `--freeze-links/--no-freeze-links` | PDB/mmCIF input (or XYZ/GJF with `--ref-pdb`). Freeze parents of cap hydrogens and merge with `geom.freeze_atoms`. See [extract](extract.md) for cap-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-write INT` | Number of modes to export. | `10` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm⁻¹) or `abs`. | `value` |
| `-o, --out-dir TEXT` | Output directory. | `./result_freq/` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). On the CLI this flag is `--pressure`; the matching YAML key under `thermo:` is `pressure_atm` (explicit unit suffix). Both are in atm and get converted to Pa internally. | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. Standalone `freq` defaults to off. `pdb2reaction all --thermo` always retains this internal file because the composite workflow consumes it; `all --no-dump` still controls optional scan/MEP/TS trajectories but does not suppress the thermochemistry channel. | `False` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/CIF companions when a PDB/mmCIF topology is available (GJF is not written). | `True` |
| `--ref-pdb FILE` | Reference PDB or mmCIF topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running frequency analysis. | `False` |

## YAML configuration

The `geom`, `calc`, `freq`, and `thermo` sections follow the canonical definitions in [YAML Reference](yaml-reference.md): see [`geom`](yaml-reference.md#geom), [`calc`](yaml-reference.md#calc), [`freq`](yaml-reference.md#freq-section), and [`thermo`](yaml-reference.md#thermo). `freq` forces `calc.return_partial_hessian = true` (PHVA) regardless of YAML.

The only `freq`-specific default that differs from the canonical block is the output directory:

```yaml
freq:
 zero_cutoff_cm: 5.0 # remove |frequency| <= 5.0 cm^-1
 out_dir: ./result_freq/ # freq default
```

## Notes

- `tsopt` already includes an imaginary-frequency check, so a separate `freq` run is mainly for thermochemistry or detailed mode inspection.
- A properly converged first-order saddle point (TS) is expected to have **exactly one** imaginary frequency. `freq.zero_cutoff_cm` removes `|frequency| <= cutoff` consistently before standalone and TS-side classification.
- Imaginary frequencies are reported as negative values in cm⁻¹. `freq` prints how many were detected
  and dumps details when `--dump`.
- An all-frozen structure has no active vibrational DOF and raises an explicit error.
- `--hessian-calc-mode` follows the standard precedence (defaults < config < explicit CLI); an explicit CLI `--hessian-calc-mode` value takes precedence over `calc.hessian_calc_mode` in the config YAML.

## See Also

- [tsopt](tsopt.md) — Optimize TS candidates (includes imaginary-frequency check; follow with IRC for endpoint validation)
- [irc](irc.md) — IRC from TS (freq is often run on IRC endpoints for thermochemistry)
- [dft](dft.md) — Single-point DFT for higher-level energy evaluation
- [all](all.md) — End-to-end workflow with `--thermo`
- [YAML Reference](yaml-reference.md) — Full `freq` and `thermo` configuration options
- [Glossary](glossary.md) — Definitions of ZPE, Gibbs Energy, Enthalpy, Entropy
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed fixes for common failure modes
