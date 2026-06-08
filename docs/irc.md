# `irc`

Runs EulerPC (Euler Predictor-Corrector)-based intrinsic reaction coordinate (IRC) integration from an optimized transition state (validated by `tsopt`) toward reactants and products, confirming endpoint connectivity (R ↔ TS ↔ P). By default both forward and backward branches are computed; use `--no-backward` (or `--no-forward`) to follow only one direction. Setting `--hessian-calc-mode Analytical` is strongly recommended when VRAM permits. For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB output conversion. A typical workflow is `tsopt` → `irc`.

## Examples

Command synopsis:

```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] [-m 2S+1] \
 [--max-cycles N] [--step-size Δs] [--root k] \
 [--forward/--no-forward] [--backward/--no-backward] \
 [--freeze-links/--no-freeze-links] \
 [--out-dir DIR] [--config FILE] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--show-config] [--dry-run]
```

Basic run with both branches:

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

Forward-only branch, finite-difference Hessian, larger step size:

```bash
# Forward-only branch, finite-difference Hessian, larger step size
pdb2reaction irc -i ts.pdb -q 0 -m 1 --no-backward \
 --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/
```

Increase step size and use analytical Hessians:

```bash
# Increase step size and use analytical Hessians
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.20 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

## Workflow

1. **Input preparation** – Any format supported by `geom_loader` is accepted. When a reference PDB is available (input is `.pdb` or `--ref-pdb` is supplied), EulerPC trajectories are converted to PDB using that topology, and `--freeze-links` augments `geom.freeze_atoms` by freezing parents of link hydrogens for PDB inputs.
2. **EulerPC integration** – The EulerPC predictor-corrector integrator traces the IRC path from the transition state. Forward and/or backward branches are run according to `--forward`/`--backward` flags. Each step uses an Euler predictor along the mass-weighted steepest-descent direction (with the gradient approximated via a second-order Taylor expansion using the current Hessian), followed by a modified-Bulirsch–Stoer corrector on a distance-weighted-interpolation surface.
3. **Trajectory output** – Finished, forward, and backward IRC trajectories are written as XYZ files. When a reference PDB is available, PDB companions are also generated (`--convert-files`).

## Outputs

```text
out_dir/ (default:./result_irc/)
├─ <prefix>finished_irc_trj.xyz   # Complete IRC trajectory
├─ <prefix>finished_irc.pdb       # PDB companion (when ref PDB available + conversion enabled)
├─ <prefix>forward_irc_trj.xyz    # Present when the forward branch runs
├─ <prefix>forward_irc.pdb        # Forward-branch PDB companion (same gating)
├─ <prefix>backward_irc_trj.xyz   # Present when the backward branch runs
└─ <prefix>backward_irc.pdb       # Backward-branch PDB companion (same gating)
```
- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus wall-clock timing.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge; used unless YAML sets `calc.charge`. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Explicit `-q` still overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); used unless YAML sets `calc.spin`. | `.gjf` template value or `1` |
| `--max-cycles INT` | Maximum IRC steps; used unless YAML sets `irc.max_cycles`. | `125` |
| `--step-size FLOAT` | Step length in unweighted Cartesian coordinates (Bohr); used unless YAML sets `irc.step_length`. | `0.10` |
| `--root INT` | **0-based** index into the projected Hessian's eigenvalues sorted in **ascending order** (most-negative first), used to pick the mode for the initial IRC displacement. For a validated TS with exactly one imaginary mode, leave `--root 0` (the sole negative eigenvalue). Use `--root 1`, `--root 2`, … only if you know the active imaginary mode is ranked above more-negative spurious modes. Used unless YAML sets `irc.root`. | `0` |
| `--forward/--no-forward` | Run forward branch (`irc.forward`), used unless YAML sets `irc.forward`. | `True` |
| `--backward/--no-backward` | Run backward branch (`irc.backward`), used unless YAML sets `irc.backward`. | `True` |
| `--freeze-links/--no-freeze-links` | For PDB inputs, freeze link-H parents (merged with `geom.freeze_atoms`). See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `-o, --out-dir TEXT` | Output directory (`irc.out_dir`), used unless YAML sets `irc.out_dir`. | `./result_irc/` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB companions when a reference PDB is available. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`calc.hessian_calc_mode`), used unless YAML sets `calc.hessian_calc_mode`. | `FiniteDifference` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. | `False` |

## YAML configuration

See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

The `geom`, `calc`, and `irc` sections are unchanged from the canonical definitions in [YAML Reference](yaml-reference.md): see [`geom`](yaml-reference.md#geom), [`calc`](yaml-reference.md#calc), and [`irc`](yaml-reference.md#irc-section). `--freeze-links` augments `geom.freeze_atoms` for PDB inputs, and `--hessian-calc-mode` plus CLI charge/spin values supplement the merged `calc` block.

**`irc`-specific hard overrides** (applied after YAML/CLI merging, regardless of YAML values):

```yaml
geom:
 coord_type: cart # forced to cart for irc (YAML value ignored)
calc:
 return_partial_hessian: true # forced true for irc (partial Hessian with active-DOF processing)
```

## Exit codes

See {ref}`exit-codes` in CLI Conventions.

## Notes

- The MLIP backend (UMA by default) is reused throughout the IRC; aggressive `step_length` values can destabilize EulerPC.
- When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed fixes for common failure modes
- [tsopt](tsopt.md) — Optimize the TS before running IRC
- [freq](freq.md) — Full vibrational analysis and thermochemistry
- [opt](opt.md) — Optimize IRC endpoints to true minima
- [all](all.md) — End-to-end workflow that runs IRC after tsopt
- [YAML Reference](yaml-reference.md) — Full `irc` configuration options
- [Glossary](glossary.md) — Definition of IRC (Intrinsic Reaction Coordinate)
