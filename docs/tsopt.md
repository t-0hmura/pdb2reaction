# `tsopt`

`pdb2reaction tsopt` optimizes a transition-state *candidate* into a first-order saddle point. The default optimizer is **RS-I-RFO** (Restricted-Step Image Rational Function Optimization, `--opt-mode hess`); when RS-I-RFO struggles (or when full-Hessian recomputation is prohibitive), switch to the **Hessian-Guided Dimer** (`--opt-mode grad`). After convergence, `tsopt` performs a final Hessian calculation and imaginary-frequency check automatically — a validated TS should show **exactly one** imaginary frequency. A separate [`freq`](freq.md) run is only needed for full vibrational analysis or thermochemistry; the imaginary-frequency check is already included in `tsopt`. Always confirm endpoint connectivity with [`irc`](irc.md).

> **Naming note:** the CLI accepts `grad` / `dimer` (Dimer) and `hess` / `rsirfo` (RS-I-RFO, default). In YAML, use the top-level `hessian_dimer:` (Dimer) or `rsirfo:` (RS-I-RFO) blocks directly.

Use it to refine a TS guess — the HEI from `path-opt` / `path-search`, or a user-supplied structure — into an optimized first-order saddle point with a built-in imaginary-frequency check. Pick `--opt-mode hess` (RS-I-RFO, default) in most systems, since it uses a full Hessian and is more reliable. Switch to `--opt-mode grad` (Hessian-Guided Dimer) when RS-I-RFO fails to converge or full-Hessian recomputation is prohibitive. Enable `--flatten` (disabled by default) when the candidate has multiple imaginary frequencies and you need surplus-mode cleanup. For XYZ / GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping the XYZ coordinates, enabling format-aware PDB / GJF output conversion. If you need a TS guess first, run [`path-opt`](path-opt.md) (two structures) or [`path-search`](path-search.md) (two or more structures) and then optimize the HEI with `tsopt` → `irc`.

## Examples

Default RS-I-RFO optimization of a PDB candidate:

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

Dimer mode with analytical Hessian (VRAM permitting):

```bash
# Dimer mode with analytical Hessian (VRAM permitting)
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
    --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

RS-I-RFO mode driven by YAML overrides:

```bash
# RS-I-RFO mode driven by YAML overrides
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
    --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
```

RS-I-RFO mode with surplus-imaginary-mode flattening enabled:

```bash
# RS-I-RFO mode with surplus-imaginary-mode flattening enabled
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
    --opt-mode hess --flatten --out-dir ./result_tsopt_flatten
```

Add `--dump` to keep the full optimization trajectory for inspection.

## Workflow

- **Charge / spin** are resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>`).
- **Geometry loading + freeze-links** — structures are read through `pysisyphus.helpers.geom_loader`. When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- **MLIP Hessians** (default UMA) — `--hessian-calc-mode` toggles analytical vs finite-difference; both honour active (PHVA) subspaces. The MLIP backend may return only the active block when frozen atoms are present. See {ref}`hessian-evaluation` for the full Hessian-evaluation matrix.
- **Dimer mode** — the Hessian-Guided Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian (active subspace, TR-projected); the lowest eigenpair uses `torch.lobpcg` when `root == 0`, falling back to `torch.linalg.eigh`. With `--flatten`, the active Hessian is updated via Bofill (SR1/MS ↔ PSB blend; toggle via `hessian_dimer.flatten_loop_bofill`) using displacements Δx and gradient differences Δg. Each loop estimates imaginary modes, flattens once, refreshes the dimer direction, runs a Dimer + L-BFGS micro-segment, and (optionally) performs a Bofill update. Once only one imaginary mode remains, a final exact Hessian is computed for frequency analysis. If `root != 0`, that root seeds only the initial dimer direction; subsequent refreshes follow the most negative mode (`root = 0`).
- **RS-I-RFO mode** — runs the RS-I-RFO optimizer with optional Hessian reference files, R+S splitting safeguards, and micro-cycle controls defined in the `rsirfo` YAML section. With `--flatten`, when more than one imaginary mode remains after convergence the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten-iteration cap is reached.
- **Mode export + conversion** — all detected imaginary modes are written to `vib/imag_*_trj.xyz` and mirrored to `.pdb` when the input was PDB and conversion is enabled. The optimization trajectory and final geometry are also converted to PDB via the input template when `--dump`; Gaussian templates receive a `.gjf` companion for the final geometry only.

## Outputs

Validate a run by opening `final_geometry.*` (the optimized saddle point) and the `vib/imag_*` modes (expect exactly one for a valid TS).

```text
out_dir/   (default: ./result_tsopt/)
├─ final_geometry.xyz              # Always written
├─ final_geometry.pdb              # When the input was PDB (conversion enabled)
├─ final_geometry.gjf              # When the input was Gaussian (conversion enabled)
├─ optimization_all_trj.xyz        # Dimer-mode dump (--dump)
├─ optimization_all.pdb            # Dimer-mode PDB companion (--dump, PDB input)
├─ optimization_trj.xyz            # RSIRFO-mode trajectory (--dump)
├─ optimization.pdb                # RSIRFO-mode PDB companion (--dump)
├─ vib/
│  ├─ imag_±XXXX.Xcm-1_trj.xyz
│  └─ imag_±XXXX.Xcm-1.pdb
└─ .dimer_mode.dat                 # Dimer-mode orientation seed
```

Exit codes: see {ref}`exit-codes` in CLI Conventions.

## CLI options

Command form:

```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l 'RES:Q,...'] [-m 2S+1] \
    [-b uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [--opt-mode grad|hess|dimer|rsirfo|trim|rsprfo] [--flatten / --no-flatten] \
    [--freeze-links / --no-freeze-links] [--max-cycles N] [--thresh PRESET] \
    [--hessian-calc-mode Analytical|FiniteDifference] \
    [--convert-files / --no-convert-files] [--ref-pdb FILE]
```

`pdb2reaction tsopt --help` shows core options; `pdb2reaction tsopt --help-advanced` shows the full option list. For full input-file requirements (hydrogens, element columns, atom-order parity, charge specification), see [CLI Conventions](cli-conventions.md).

The tables below cover the options that need explanation. The full flag list is in the generated [command reference](reference/commands/index.md); do not hand-duplicate it here.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader` (`.pdb` / `.xyz` / `.gjf` / `.trj`). | Required |
| `-q, --charge INT` | Net charge. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB inputs or XYZ / GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template / derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g. `-1`) for the total ligand charge, or a per-residue mapping (e.g. `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted. | _None_ |
| `--workers INT`, `--workers-per-node INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; see {ref}`workers-fd-downgrade`). | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links / --no-freeze-links` | PDB-only. Freeze parents of link hydrogens (merged into `geom.freeze_atoms`). | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g. `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | TS optimizer preset (Choice: `grad` / `hess` / `dimer` / `rsirfo` / `trim` / `rsprfo`). `grad` and `dimer` → Hessian-Guided Dimer; `hess` and `rsirfo` → RS-I-RFO (default); `trim` → TRIM (Helgaker, non-microiter); `rsprfo` → RS-P-RFO (Banerjee, non-microiter). On `opt`, the same `grad` token picks L-BFGS minimization instead — see {ref}`opt-mode-semantics`. | `hess` |
| `--dump / --no-dump` | Dump trajectories. | `False` |
| `-o, --out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--flatten / --no-flatten` | Enable the surplus-imaginary-mode flattening loop (`False` forces `flatten_max_iter = 0`). After TS optimization converges, iteratively flattens surplus negative-eigenvalue modes until only one imaginary frequency remains (or the iteration cap is reached). Applies to both Dimer (dimer loop) and RS-I-RFO (post-convergence). | `False` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--convert-files / --no-convert-files` | Toggle XYZ / TRJ → PDB / GJF companions for PDB or Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology when the input is XYZ / GJF (keeps XYZ coordinates). | _None_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config / --no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--out-json / --no-out-json` | Write a machine-readable `result.json` to `out_dir`. Schema: [JSON Output Schema](json-output.md). | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run / --no-dry-run` | Validate inputs / config and print the execution plan without running TS optimization. | `False` |

(flatten-precedence-caveat)=
### `--flatten` precedence caveat

```{note}
**`--flatten` is disabled by default (precedence caveat).** Although `defaults.py` defines `flatten_max_iter: 50`, the CLI initializer seeds `flatten_max_iter = 0`. Effective resolution:

- CLI `--flatten` **not** passed → `flatten_max_iter = 0` unless **explicitly set in YAML** via `hessian_dimer.flatten_max_iter` (the flatten counter is read only from the `hessian_dimer` block for both Dimer and RS-I-RFO paths). The `defaults.py` value of 50 is ignored.
- CLI `--flatten` passed → the YAML / `defaults.py` value applies (default `flatten_max_iter = 50`); you can still override via YAML.

If your TS candidate has multiple imaginary frequencies, add `--flatten` to enable the surplus-mode cleanup loop.
```

## YAML configuration

Shared sections reuse [YAML Reference](yaml-reference.md); adjust only the values you need to change. `geom` and `calc` are unchanged from canonical — see [`geom`](yaml-reference.md#geom) and [`calc`](yaml-reference.md#calc). The `opt` block uses the same keys as [`opt`](yaml-reference.md#opt) with these `tsopt`-specific defaults:

```yaml
opt:
  thresh: baker                # tsopt default (vs. `gau` for `opt`)
  out_dir: ./result_tsopt/     # tsopt default (vs. `./result_opt/` for `opt`)
```

```{note}
**Energy-plateau fallback.** The RS-I-RFO optimizer honours the shared `energy_plateau` setting: if the energy range (max − min) over the last `energy_plateau_window` (default 50) steps falls below `energy_plateau_thresh` (default `1×10⁻⁴ au ≈ 0.06 kcal/mol`), the run is declared converged. This matters for large TS systems where the MLIP force noise floor (~4×10⁻⁴ au) exceeds the `baker` `max_force` threshold (3×10⁻⁴ au), making the force criterion unreachable even though the energy has plainly flattened. Set `energy_plateau: false` to disable.
```

### Dimer mode (`--opt-mode grad`)

Used with `--opt-mode grad` (Hessian-Guided Dimer + L-BFGS translation). The full `hessian_dimer` block (including the inner `dimer:` and its nested `lbfgs:`) is documented in [`hessian_dimer`](yaml-reference.md#hessian_dimer); the inner `lbfgs:` inherits from [`lbfgs`](yaml-reference.md#lbfgs), with this `tsopt`-specific override:

```yaml
hessian_dimer:
  dimer:
    lbfgs:
      out_dir: ./result_tsopt/   # tsopt override (defaults.py value is ./result_opt/)
```

### RS-I-RFO mode (`--opt-mode hess`, default)

Used with `--opt-mode hess` (RS-I-RFO, the default). The full `rsirfo` block is documented in [`rsirfo`](yaml-reference.md#rsirfo) (which inherits trust-region and Hessian-update keys from [`rfo`](yaml-reference.md#rfo)). `tsopt`-specific overrides:

```yaml
rsirfo:
  trust_max: 0.10              # maximum trust radius (bohr)
  out_dir: ./result_tsopt/     # tsopt override (defaults.py value is ./result_opt/)
  hessian_recalc: 500          # rebuild exact Hessian every N macro steps
```

```{tip}
Set `rsirfo.track_mode_by_overlap: true` if the TS mode switches root during optimization (e.g. when multiple imaginary frequencies are present). If TS convergence is slow or the TS mode is lost, lowering `hessian_recalc` (e.g. to 50–200) helps — more frequent exact Hessian recalculations improve robustness at the cost of additional Hessian evaluations.
```

## Notes

- Imaginary-frequency **detection** threshold defaults to 5.0 cm⁻¹ (configurable via `hessian_dimer.neg_freq_thresh_cm`); frequencies with magnitudes below this threshold are not counted as imaginary. The selected `root` controls which vibrational mode is followed during optimization; it is set via YAML (`rsirfo.root` or `hessian_dimer.root`; default `0`) — `tsopt` has no `--root` CLI flag (unlike [`irc`](irc.md)).
- Use `--opt-mode` to choose the algorithm directly (`rsirfo` by default) rather than editing YAML mode mappings.
- Dimer mode applies translation / rotation projection (PHVA when frozen atoms are present) before the initial Hessian diagonalization, matching the `freq` implementation. RS-I-RFO mode operates directly on the active-DOF Cartesian Hessian without TR projection (frozen atoms remove the rigid-body symmetry).
- See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

## See Also

[Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [path-search](path-search.md) · [irc](irc.md) · [freq](freq.md) · [all](all.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
