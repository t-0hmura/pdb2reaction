# `opt`

## Overview

> **Summary:** Optimizes a single structure toward a local minimum using L-BFGS (`--opt-mode grad`, default) or RFO (`--opt-mode hess`). Optional imaginary-frequency flattening can be enabled with `--flatten`.

### At a glance
- **Use when:** Relaxing a single structure (PDB/XYZ/GJF/`_trj.xyz`) to a local minimum, optionally with distance restraints or imaginary-mode flattening.
- **Method:** pysisyphus L-BFGS (`--opt-mode grad`, default) or RFOptimizer (`--opt-mode hess`) driven by an MLIP backend (UMA by default; ORB/MACE/AIMNet2 via `-b`); link-hydrogen parents on PDB inputs are auto-frozen with `--freeze-links`.
- **Outputs:** `final_geometry.xyz` (always), plus `final_geometry.pdb` (PDB inputs) and `final_geometry.gjf` (Gaussian templates) when `--convert-files` is enabled; `optimization_trj.xyz`/`optimization.pdb` when `--dump` is on; optional restart YAML.
- **Defaults:** Backend `uma`, `--opt-mode grad`, `--thresh gau`, `--max-cycles 10000`, `--bias-k 300`, `--freeze-links True`, `--convert-files True`, `--flatten False`, `--dump False`, `--out-dir ./result_opt/`.
- **Next step:** Confirm the minimum with [`freq`](freq.md), or proceed to [`tsopt`](tsopt.md) for saddle points and [`path-search`](path-search.md) / [`all`](all.md) for full reaction pathways.

`pdb2reaction opt` optimizes a single structure to a local minimum using L-BFGS (`--opt-mode grad`, default) or RFO (`--opt-mode hess`). For PDB inputs, link-hydrogen parents are automatically frozen.

The command uses pysisyphus LBFGS (`lbfgs`) or RFOptimizer (`rfo`) while an MLIP backend (UMA by default; ORB, MACE, and AIMNet2 also available via `-b/--backend`) provides energies, gradients, and Hessians. Input structures can be `.pdb`, `.xyz`, `_trj.xyz`, or any format supported by `geom_loader`. Settings follow precedence: **defaults < config < explicit CLI**.

When the starting structure is a PDB or Gaussian template, the command also writes `.pdb` (PDB inputs) and `.gjf` (Gaussian templates) companions, controlled by `--convert-files/--no-convert-files` (enabled by default). PDB-specific conveniences include:
- With `--freeze-links` (default `True`), parent atoms of link hydrogens are detected and merged into `geom.freeze_atoms` (1-based indices).
- Output conversion produces `final_geometry.pdb` (and `optimization.pdb` when dumping trajectories) using the input PDB as the topology reference.
For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates,
enabling format-aware PDB/GJF output conversion.

A Gaussian `.gjf` template seeds the charge/spin defaults and enables automatic export of the optimized structure as `.gjf` when conversion is enabled.

## Minimal example

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --out-dir ./result_opt
```

## Output checklist

- `result_opt/final_geometry.xyz`
- `result_opt/final_geometry.pdb` (for PDB input with conversion enabled)
- `result_opt/optimization_trj.xyz` (when `--dump` is enabled)

## Common examples

1. Optimize with a tighter threshold and keep trajectory dumps.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --thresh gau_tight --dump \
 --out-dir ./result_opt_tight
```

2. Add a harmonic distance restraint.

```bash
# Restrain atoms 1 and 5 to 2.0 Å. The example uses --bias-k 20.0 (a loose restraint
# suitable for a light guide near the target distance); the default `bias.k` is 300
# eV·Å⁻² and is better when you want the restraint to dominate during optimization.
pdb2reaction opt -i input.pdb -q 0 -m 1 \
 --dist-freeze '[(1,5,2.0)]' --bias-k 20.0 --out-dir ./result_opt_rest

# 2-tuple form: restrain atoms 1 and 5 to their current distance
pdb2reaction opt -i input.pdb -q 0 -m 1 \
 --dist-freeze '[(1,5)]' --out-dir ./result_opt_freeze
```

3. Switch explicitly to RFO mode.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess \
 --out-dir ./result_opt_hess
```

4. Run LBFGS mode and flatten imaginary modes after optimization.

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode grad --flatten \
 --out-dir ./result_opt_grad_flat
```

## Usage
```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--opt-mode grad|hess|lbfgs|rfo] [--flatten/--no-flatten] [--freeze-links/--no-freeze-links] \
 [--dist-freeze '[(i,j,target_Å),...]'] [--one-based|--zero-based] \
 [--bias-k K_eV_per_Å²] [--dump/--no-dump] [-o/--out-dir DIR] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

## Workflow
- **Optimizers**: `--opt-mode grad` (alias: `lbfgs`, default) → L-BFGS; `--opt-mode hess` (alias: `rfo`) → RFOptimizer.
  > **Naming note:** The CLI accepts `grad|lbfgs` and `hess|rfo`. In YAML, use `lbfgs` or `rfo` directly.
- **Flatten loop**: `--flatten` enables post-optimization flattening of imaginary vibrational modes. In `opt`, all detected imaginary modes are flattened each iteration until none remain or the internal loop cap is reached.
- **Restraints**: `--dist-freeze` consumes Python-literal tuples `(i, j, target_Å)` where `target_Å` is the target distance in Å; omitting the third element restrains the starting distance. `--bias-k` sets a global harmonic strength (eV·Å⁻²). Indices default to 1-based but can be flipped to 0-based with `--zero-based`.
- **Charge/spin resolution**: Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>` for details).
- **Freeze atoms**: When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- **Dumping & conversion**: `--dump` mirrors `opt.dump=True` and writes `optimization_trj.xyz`; when conversion is enabled, trajectories are mirrored to `.pdb` for PDB inputs. `opt.dump_restart` can emit restart YAML snapshots.
- **Exit codes**: See {ref}`exit-codes` in CLI Conventions.

## CLI options

> **Note:** Default values shown are used when the option is not specified.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Net charge. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Falls back to `.gjf` template or `1`. | Template/`1` |
| `--dist-freeze TEXT` | Repeatable string parsed as Python literal describing `(i,j,target_Å)` tuples for harmonic restraints. | _None_ |
| `--one-based/--zero-based` | Interpret `--dist-freeze` indices as 1-based (default) or 0-based. | `True` |
| `--bias-k FLOAT` | Harmonic bias strength applied to every `--dist-freeze` tuple (eV·Å⁻²). | `300` |
| `--freeze-links/--no-freeze-links` | Toggle link-hydrogen parent freezing (PDB inputs only). See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-cycles INT` | Hard limit on optimization iterations (`opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Optimizer preset: `grad` (`lbfgs`) or `hess` (`rfo`). Aliases `lbfgs`/`rfo` are accepted. For the full subcommand-dependent table (`opt` uses L-BFGS/RFO; `tsopt` uses Dimer/RS-I-RFO), see {ref}`opt-mode-semantics`. | `grad` |
| `--flatten/--no-flatten` | Enable/disable the post-optimization imaginary-mode flattening loop. | `False` |
| `--dump/--no-dump` | Emit trajectory dumps (`optimization_trj.xyz`). | `False` |
| `--convert-files/--no-convert-files` | Enable or disable XYZ/TRJ → PDB companions for PDB inputs and XYZ → GJF companions for Gaussian templates. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `-o, --out-dir TEXT` | Output directory for all files. | `./result_opt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--config FILE` | Base YAML configuration file. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layer information before execution. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running optimization. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |

## Outputs
```
out_dir/
├─ final_geometry.xyz # Always written
├─ final_geometry.pdb # Only when the input was a PDB and conversion is enabled
├─ final_geometry.gjf # When a Gaussian template was detected and conversion is enabled
├─ optimization_trj.xyz # Only if dumping is enabled
├─ optimization.pdb # PDB conversion of the trajectory (PDB inputs, conversion enabled)
└─ restart*.yml # Optional restarts when opt.dump_restart is set
```
The console prints the resolved `geom`, `calc`, `opt`, `lbfgs`/`rfo` blocks plus cycle-by-cycle progress and total runtime.

See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

### `geom`
- `coord_type` (`"cart"`): Cartesian vs. `"dlc"` delocalized internal coordinates.
- `freeze_atoms` (`[]`): 1-based frozen indices; automatically merged with CLI link detection.

### `calc`
- MLIP backend configuration (`model`, `task_name`, device selection, neighbor radii, Hessian format, etc.).
- `charge`/`spin` mirror the CLI options; defaults come from `.gjf` when present.

### `opt`
Shared optimizer controls used by both LBFGS and RFO:
- `thresh` presets (Gaussian-like or Baker rule). Presets translate to the force/step thresholds documented in `pdb2reaction/opt.py`.
- `max_cycles`, `print_every` (`100`), `min_step_norm` (`1e-8`), `assert_min_step`, convergence toggles (`rms_force`, etc.), RMSD-based `converge_to_geom_rms_thresh`, `overachieve_factor`, `check_eigval_structure`, `line_search`.
- Energy plateau fallback (`energy_plateau`, `energy_plateau_thresh`, `energy_plateau_window`) — declares convergence when the recent energy range flattens (useful when MLIP force noise prevents force-based convergence; see note below).
- Dumping/bookkeeping fields (`dump`, `dump_restart`, `prefix`, `out_dir`).

### `lbfgs`
Extends `opt` with L-BFGS specifics: `keep_last`, `beta`, `gamma_mult`, `max_step`, `control_step`, `double_damp`, and the optional regularization parameters `mu_reg`/`max_mu_reg_adaptions`.

### `rfo`
Extends `opt` with RFOptimizer fields: trust-region sizing (`trust_radius`, `trust_min`, `trust_max`, `trust_update`), `max_energy_incr`, Hessian management (`hessian_update`, `hessian_init`, `hessian_recalc`, `hessian_recalc_adapt`, `small_eigval_thresh`), micro-iteration controls (`alpha0`, `max_micro_cycles`, `rfo_overlaps`), DIIS helpers (`gdiis`, `gediis`, thresholds, `gdiis_test_direction`), and `adapt_step_func`.

### Opt-specific defaults

The `opt` subcommand sets `opt.out_dir` (and `lbfgs.out_dir` / `rfo.out_dir`) to `./result_opt/`.

For the full YAML schema of `geom`, `calc`, `opt`, `lbfgs`, and `rfo`, see [YAML Reference](yaml-reference.md).

```{note}
**Energy plateau fallback (new in v0.3.5).** When `energy_plateau: true`, the optimizer
is declared converged if the energy range (max − min) over the last
`energy_plateau_window` steps falls below `energy_plateau_thresh`
(default `1×10⁻⁴ au ≈ 0.06 kcal/mol` over 50 steps). This prevents wasted cycles when
the MLIP force noise floor (~4×10⁻⁴ au) exceeds the force-based convergence threshold
(e.g. `baker` max_force = 3×10⁻⁴ au). The fallback is skipped for chain-of-states
optimizers, which store per-image energy arrays.
```

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [tsopt](tsopt.md) — Optimize transition states (saddle points) instead of minima
- [freq](freq.md) — Vibrational analysis to confirm optimization reached a minimum
- [extract](extract.md) — Generate active site model (binding pocket) PDBs before optimization
- [all](all.md) — End-to-end workflow that pre-optimizes endpoints
- [YAML Reference](yaml-reference.md) — Full `opt`, `lbfgs`, `rfo` configuration options
- [Glossary](glossary.md) — Definitions of L-BFGS, RFO
