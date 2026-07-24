# `irc`

Runs EulerPC (Euler Predictor-Corrector)-based intrinsic reaction coordinate (IRC) integration from an optimized transition state (validated by `tsopt`) toward reactants and products, confirming endpoint connectivity (R ↔ TS ↔ P). By default both forward and backward branches are computed; use `--no-backward` (or `--no-forward`) to follow only one direction. Hessians default to finite differences; analytical autograd is an explicit alternative whose speed and memory cost depend on the backend, model, system size, precision, and GPU, so validate it on the target setup instead of treating it as a blanket recommendation. For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB/mmCIF topology while keeping the XYZ coordinates, enabling format-aware companion output. A typical workflow is `tsopt` → `irc`.

## Examples

Command synopsis:

```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--workers N] [--workers-per-node N] [-m 2S+1] \
 [--max-cycles N] [--step-size Δs] [--never-stop/--no-never-stop] [--root k] \
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

Request an analytical Hessian explicitly:

```bash
# Keep workers at 1 when requesting an analytical UMA Hessian
pdb2reaction irc -i ts.pdb -q 0 -m 1 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

If a branch stops after only one or two steps, first reduce the maximum
EulerPC step. `0.05` Bohr is a useful retry value:

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 \
 --out-dir ./result_irc_small_step
```

For a path that contains a small numerical uphill or flat shoulder, add
`--never-stop`. This opt-in mode ignores energy-increase and energy-change
stops, but still obeys gradient/integrator convergence and `--max-cycles`:

```bash
pdb2reaction irc -i ts.pdb -q 0 -m 1 --step-size 0.05 --never-stop \
 --max-cycles 250 --out-dir ./result_irc_continue
```

## Workflow

1. **Input preparation** – The common bridge accepts PDB, mmCIF, and formats supported by `geom_loader`. When a PDB/mmCIF reference topology is available, EulerPC trajectories are converted to PDB; bridge inputs also receive CIF with original IDs. `--freeze-links` augments `geom.freeze_atoms` by freezing parents of cap hydrogens for topology inputs.
2. **EulerPC integration** – The EulerPC predictor-corrector integrator traces the IRC path from the transition state. Forward and/or backward branches are run according to `--forward`/`--backward` flags. The default `constrained` rigid-mode treatment removes only full-system translations/rotations compatible with the frozen anchors before selecting the initial mode. Each step then uses an Euler predictor along the mass-weighted steepest-descent direction (with the gradient approximated via a second-order Taylor expansion using the current Hessian), followed by a modified-Bulirsch–Stoer corrector on a distance-weighted-interpolation surface.
3. **Trajectory output** – Finished, forward, and backward IRC trajectories are written as XYZ files. With a reference topology and `--convert-files`, PDB companions are generated; mmCIF/oversized-PDB bridge inputs also receive CIF companions.

## Outputs

```text
out_dir/ (default:./result_irc/)
├─ <prefix>finished_irc_trj.xyz   # Complete IRC trajectory
├─ <prefix>finished_irc.pdb       # PDB companion (when ref PDB available + conversion enabled)
├─ <prefix>finished_irc.cif       # Bridge-input companion with original IDs
├─ <prefix>forward_irc_trj.xyz    # Present when the forward branch runs
├─ <prefix>forward_irc.pdb        # Forward-branch PDB companion (same gating)
├─ <prefix>forward_irc.cif        # Bridge-input companion
├─ <prefix>backward_irc_trj.xyz   # Present when the backward branch runs
├─ <prefix>backward_irc.pdb       # Backward-branch PDB companion (same gating)
└─ <prefix>backward_irc.cif       # Bridge-input companion
```

When `irc.prefix` is non-empty, EulerPC inserts one underscore before the
filename; for example, `prefix: trial` produces
`trial_finished_irc_trj.xyz`. `result.json.files` records the normalized names.

Low-level periodic HDF5 checkpointing is available only through YAML:
`irc.dump_every` defaults to `null` and a positive value writes
`<prefix>irc_data.h5`. The file is overwritten with the current direction's
coordinates, energies, and gradients; it is not a final combined IRC artifact,
contains no Hessian, and is omitted from `result.json.files`.

- Console summaries of resolved `geom`, `calc`, and `irc` configurations plus wall-clock timing.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Transition-state structure accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge. Explicit `-q` overrides YAML `calc.charge` and `--ligand-charge/-l`; when omitted, YAML, residue derivation, or `.gjf` metadata may supply it. | Required unless YAML/template/derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB/mmCIF residue metadata. Used when `-q` is omitted (PDB/mmCIF inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | UMA predictor parallelism. `workers > 1` cannot be combined with an explicit analytical Hessian request; use `workers = 1` or finite differences. See {ref}`workers-analytical-error`. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). Explicit `-m` overrides YAML `calc.spin`; otherwise YAML, `.gjf`, or `1` is used. | YAML/`.gjf`/`1` |
| `--max-cycles INT` | Maximum IRC steps. An explicit value overrides YAML `irc.max_cycles`. | `125` |
| `--step-size FLOAT` | Step length in unweighted Cartesian coordinates (Bohr). An explicit value overrides YAML `irc.step_length`. | `0.10` |
| `--never-stop/--no-never-stop` | Ignore transient energy increases/plateaus so EulerPC can cross a small shoulder. It does not disable gradient/integrator convergence or the `max_cycles` hard cap. | `False` |
| `--root INT` | **0-based** index into the projected Hessian's eigenvalues sorted in **ascending order** (most-negative first), used to pick the mode for the initial IRC displacement. For a validated TS with exactly one imaginary mode, leave `--root 0` (the sole negative eigenvalue). Use `--root 1`, `--root 2`, … only if you know the active imaginary mode is ranked above more-negative spurious modes. An explicit value overrides YAML `irc.root`. | `0` |
| `--forward/--no-forward` | Run forward branch (`irc.forward`); an explicit toggle overrides YAML. | `True` |
| `--backward/--no-backward` | Run backward branch (`irc.backward`); an explicit toggle overrides YAML. | `True` |
| `--irc-pos-def/--no-irc-pos-def` | Opt in to requiring a positive-definite projected Hessian before accepting IRC convergence. Enable this guard when a shoulder could otherwise look converged. | `False` |
| `--freeze-links/--no-freeze-links` | For PDB/mmCIF topology inputs, freeze cap-H parents (merged with `geom.freeze_atoms`). See [extract](extract.md) for cap-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--tr-projection [constrained\|legacy-active]` | Rigid-mode treatment for the initial frozen/partial Hessian. `legacy-active` is deprecated comparison-only behavior and must not be used for pass/HOSP transition-state certification. | `constrained` |
| `-o, --out-dir TEXT` | Output directory (`irc.out_dir`); an explicit value overrides YAML. | `./result_irc/` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/CIF companions when a reference PDB/mmCIF topology is available. | `True` |
| `--ref-pdb FILE` | Reference PDB or mmCIF topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`calc.hessian_calc_mode`); an explicit value overrides YAML. | `FiniteDifference` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. | `False` |

## YAML configuration

See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

The `geom`, `calc`, and `irc` sections are unchanged from the canonical definitions in [YAML Reference](yaml-reference.md): see [`geom`](yaml-reference.md#geom), [`calc`](yaml-reference.md#calc), and [`irc`](yaml-reference.md#irc-section). `--freeze-links` augments `geom.freeze_atoms` for PDB/mmCIF topology inputs, and `--hessian-calc-mode` plus CLI charge/spin values supplement the merged `calc` block.

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

- The MLIP backend (UMA by default) is reused throughout the IRC; aggressive `step_length` values can destabilize EulerPC. A branch that stops almost immediately should be retried with a smaller `--step-size` (for example `0.05`) before changing other controls.
- `--never-stop` is intentionally off by default. It is useful for a small numerical hill/shoulder, but it can also trace farther than the chemically intended basin; inspect the trajectory and endpoint connectivity.
- When `--freeze-links` is active, cap-hydrogen parent atoms are automatically frozen (see {ref}`Cap hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- `result.json["rigid_projection"]` records the treatment, effective rank, and initial Hessian source and shape. See [Frozen Atoms](freeze-atoms.md#rigid-modes-with-frozen-boundaries).

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed fixes for common failure modes
- [tsopt](tsopt.md) — Optimize the TS before running IRC
- [freq](freq.md) — Full vibrational analysis and thermochemistry
- [opt](opt.md) — Optimize IRC endpoints to true minima
- [all](all.md) — End-to-end workflow that runs IRC after tsopt
- [YAML Reference](yaml-reference.md) — Full `irc` configuration options
- [Glossary](glossary.md) — Definition of IRC (Intrinsic Reaction Coordinate)
