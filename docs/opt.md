# `opt`

Optimizes a single structure toward a local minimum using L-BFGS (`--opt-mode grad`, default) or RFO (`--opt-mode hess`). The command uses pysisyphus L-BFGS (`lbfgs`) or RFOptimizer (`rfo`) while an MLIP backend (UMA by default; ORB, MACE, and AIMNet2 also available via `-b/--backend`) provides energies, gradients, and Hessians. Input structures can be `.pdb`, `.xyz`, `_trj.xyz`, or any format supported by `geom_loader`. Optional imaginary-frequency flattening can be enabled with `--flatten`.

## When to use

- Relax a single structure (PDB/XYZ/GJF/`_trj.xyz`) to a local minimum, optionally with distance restraints or imaginary-mode flattening.
- Pick `--opt-mode grad` (alias `lbfgs`, default) for L-BFGS minimisation; pick `--opt-mode hess` (alias `rfo`) for RFOptimizer.

## Quick examples

```bash
pdb2reaction opt -i input.pdb -q 0 -m 1 --out-dir ./result_opt
```

```bash
# Tighter threshold and keep trajectory dumps
pdb2reaction opt -i input.pdb -q 0 -m 1 --thresh gau_tight --dump \
 --out-dir ./result_opt_tight
```

```bash
# Add a harmonic distance restraint. The example uses --bias-k 20.0 (a loose restraint
# suitable for a light guide near the target distance); the default `bias.k` is 300
# eV·Å⁻² and is better when you want the restraint to dominate during optimization.
pdb2reaction opt -i input.pdb -q 0 -m 1 \
 --dist-freeze '[(1,5,2.0)]' --bias-k 20.0 --out-dir ./result_opt_rest
# 2-tuple form restrains atoms 1 and 5 to their current distance: --dist-freeze '[(1,5)]'
```

```bash
# Switch explicitly to RFO mode
pdb2reaction opt -i input.pdb -q 0 -m 1 --opt-mode hess \
 --out-dir ./result_opt_hess
```

## Inputs

Command form:

```bash
pdb2reaction opt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--opt-mode grad|hess|lbfgs|rfo] [--flatten/--no-flatten] [--freeze-links/--no-freeze-links] \
 [--dist-freeze '[(i,j,target_Å),...]'] [--one-based|--zero-based] \
 [--bias-k K_eV_per_Å²] [--dump/--no-dump] [-o/--out-dir DIR] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

| Input | Required | Notes |
| --- | --- | --- |
| `-i, --input` | yes | Input structure accepted by `geom_loader` (`.pdb`, `.xyz`, `_trj.xyz`, `.gjf`). |
| `-q, --charge` | unless template/derivation applies | Net charge; supplied by a `.gjf` template or derived from `-l/--ligand-charge` when omitted. |
| `-m, --multiplicity` | no | Spin multiplicity (2S+1); falls back to `.gjf` template or `1`. |
| `--ref-pdb` | for XYZ/GJF | Reference PDB topology when the input is XYZ/GJF, enabling format-aware PDB/GJF output conversion. |

## Workflow

- **Optimizer naming**: the CLI accepts `grad|lbfgs` and `hess|rfo`; in the YAML `opt_mode` key, use `lbfgs` or `rfo` directly. See {ref}`opt-mode-semantics` for the per-subcommand token→algorithm mapping.
- **Flatten loop**: `--flatten` enables post-optimization flattening of imaginary vibrational modes. In `opt`, all detected imaginary modes are flattened each iteration until none remain or the internal loop cap is reached.
- **Restraints**: `--dist-freeze` consumes Python-literal tuples `(i, j, target_Å)` where `target_Å` is the target distance in Å; omitting the third element restrains the starting distance. `--bias-k` sets a global harmonic strength (eV·Å⁻²). Indices default to 1-based but can be flipped to 0-based with `--zero-based`.
- **Charge/spin resolution**: Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>` for details).
- **Freeze atoms**: When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- **Dumping & conversion**: `--dump` mirrors `opt.dump=True` and writes `optimization_trj.xyz`; when conversion is enabled, trajectories are mirrored to `.pdb` for PDB inputs. `opt.dump_restart` can emit restart YAML snapshots.
- **Exit codes**: See {ref}`exit-codes` in CLI Conventions.

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

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation, and the exhaustive list is not hand-duplicated here.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Net charge. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Falls back to `.gjf` template or `1`. | Template/`1` |
| `--dist-freeze TEXT` | Repeatable string parsed as Python literal describing `(i,j,target_Å)` tuples for harmonic restraints. | _None_ |
| `--one-based/--zero-based` | Interpret `--dist-freeze` indices as 1-based (default) or 0-based. | `True` |
| `--bias-k FLOAT` | Harmonic bias strength applied to every `--dist-freeze` tuple (eV·Å⁻²). | `300` |
| `--freeze-links/--no-freeze-links` | Toggle link-hydrogen parent freezing (PDB inputs only). See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-cycles INT` | Hard limit on optimization iterations (`opt.max_cycles`). | `10000` |
| `--opt-mode TEXT` | Optimizer preset: `grad` (`lbfgs`) or `hess` (`rfo`). Aliases `lbfgs`/`rfo` are accepted. On `opt`, `grad` = L-BFGS minimisation; on `tsopt`, `grad` = Hessian-Guided Dimer TS search. For the full subcommand-dependent table, see {ref}`opt-mode-semantics`. | `grad` |
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

## YAML configuration

Shared sections reuse [YAML Reference](yaml-reference.md); adjust only the values you need to change. `geom`, `calc`, `opt`, and the optimizer-specific `lbfgs`/`rfo` blocks use the canonical keys and defaults — see [`geom`](yaml-reference.md#geom), [`calc`](yaml-reference.md#calc), [`opt`](yaml-reference.md#opt), [`lbfgs`](yaml-reference.md#lbfgs), [`rfo`](yaml-reference.md#rfo). A minimal representative configuration:

```yaml
geom:
  coord_type: cart        # or `dlc` for delocalized internal coordinates
  freeze_atoms: []        # 1-based frozen indices; merged with CLI link detection
calc:
  charge: 0               # mirrors the CLI option; defaults from `.gjf` when present
  spin: 1
opt:
  thresh: gau
  max_cycles: 10000
  out_dir: ./result_opt/  # opt-specific default
```

The only `opt`-specific default is `out_dir: ./result_opt/` (also applied to `lbfgs.out_dir`/`rfo.out_dir`). The CLI accepts `grad|lbfgs` and `hess|rfo`; in YAML, use `lbfgs` or `rfo` directly.

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## Notes

```{note}
**Energy plateau fallback.** When `energy_plateau: true`, the optimizer
is declared converged if the energy range (max − min) over the last
`energy_plateau_window` steps falls below `energy_plateau_thresh`
(default `1×10⁻⁴ au ≈ 0.06 kcal/mol` over 50 steps). This prevents wasted cycles when
the MLIP force noise floor (~4×10⁻⁴ au) exceeds the force-based convergence threshold
(e.g. `baker` max_force = 3×10⁻⁴ au). The fallback is skipped for chain-of-states
optimizers, which store per-image energy arrays.
```

- **Naming note:** The CLI accepts `grad|lbfgs` and `hess|rfo`. In YAML, use `lbfgs` or `rfo` directly.

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [tsopt](tsopt.md) — Optimize transition states (saddle points) instead of minima
- [freq](freq.md) — Vibrational analysis to confirm optimization reached a minimum
- [extract](extract.md) — Generate active site model (binding pocket) PDBs before optimization
- [all](all.md) — End-to-end workflow that pre-optimizes endpoints
- [YAML Reference](yaml-reference.md) — Full `opt`, `lbfgs`, `rfo` configuration options
- [Glossary](glossary.md) — Definitions of L-BFGS, RFO
