# `tsopt`

## Overview

> **Summary:** Optimize a transition-state *candidate* using RS‑I‑RFO (Restricted-Step Image Rational Function Optimization) (`--opt-mode hess`, default) or, as an alternative when RS‑I‑RFO struggles, Hessian-Guided Dimer (`--opt-mode grad`). `tsopt` performs a final Hessian calculation and imaginary-frequency check automatically; a validated TS (first-order saddle point) should show **exactly one** imaginary frequency. Always confirm endpoint connectivity with `irc`.

### At a glance
- **Use when:** You have a TS guess (HEI from `path-opt`/`path-search`, or your own structure) and need to refine it into an optimized first-order saddle point with a built-in imaginary-frequency check.
- **Method:** `--opt-mode hess` (`rsirfo`) = RS‑I‑RFO with full Hessian (default, more reliable for most systems); `--opt-mode grad` (`dimer`) = Hessian-Guided Dimer (alternative when RS‑I‑RFO fails to converge or full-Hessian recomputation is prohibitive). `--flatten` (default disabled) controls surplus-imaginary-mode cleanup.
- **Outputs:** `final_geometry.{xyz,pdb,gjf}`, `vib/imag_*_trj.xyz` (and `.pdb` for PDB inputs); optimization trajectories with `--dump`.
- **Defaults:** `--opt-mode hess` (RS-I-RFO), `--thresh baker`, `--hessian-calc-mode FiniteDifference`, `--max-cycles 10000`, `--flatten` disabled, backend `uma`.
- **Next step:** The result is still a *candidate* until [irc](irc.md) confirms endpoint connectivity. A separate [freq](freq.md) run is only needed for full vibrational analysis or thermochemistry.

### Choosing `--opt-mode`
- Use **`--opt-mode hess` (RS‑I‑RFO)** when you want the default, conservative optimizer and you can afford Hessian work.
- Use **`--opt-mode grad` (Dimer)** when you want a lighter-weight search, or when you plan to iterate quickly from several TS guesses.

> **Naming note:** The CLI accepts `grad|dimer` (= Dimer) and `hess|rsirfo` (= RS-I-RFO, default). In YAML, use the top-level `hessian_dimer:` (Dimer) or `rsirfo:` (RS-I-RFO) blocks directly.

For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB/GJF output conversion. If you need a TS guess first, run [path-opt](path-opt.md) (two structures) or [path-search](path-search.md) (two or more structures) and then optimize the HEI with `tsopt` (which includes an imaginary-frequency check) → `irc`.

## Minimal example

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## Output checklist

- `result_tsopt/final_geometry.pdb` (or `final_geometry.xyz`)
- `result_tsopt/vib/imag_*_trj.xyz`
- `result_tsopt/vib/imag_*.pdb` (for PDB inputs)

## Common examples

1. Use dimer mode with analytical Hessian when VRAM is sufficient.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

2. Keep the optimization trajectory for inspection.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. Run rsirfo mode with YAML overrides.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
```

4. Run rsirfo mode with flattening enabled.

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
 --opt-mode hess --flatten --out-dir ./result_tsopt_flatten
```

## Usage
```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m 2S+1] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [--opt-mode grad|hess|dimer|rsirfo] [--flatten/--no-flatten] \
 [--freeze-links/--no-freeze-links] [--max-cycles N] [--thresh PRESET] \
 [--hessian-calc-mode Analytical|FiniteDifference] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### Examples
```bash
# Recommended baseline: specify charge/multiplicity and use rsirfo
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode hess --out-dir ./result_tsopt/

# Dimer mode with YAML overrides, finite-difference Hessian, and freeze-links handling
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --freeze-links \
 --opt-mode grad --max-cycles 10000 --no-dump \
 --out-dir ./result_tsopt/ --config ./args.yaml \
 --hessian-calc-mode FiniteDifference

# RS-I-RFO mode driven entirely by YAML
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --opt-mode hess \
 --config ./args.yaml --out-dir ./result_tsopt/
```

## Workflow
- **Charge/spin resolution**: Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>` for details).
- **Geometry loading & freeze-links**: structures are read via
  `pysisyphus.helpers.geom_loader`. When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- **MLIP Hessians (default: UMA)**: `--hessian-calc-mode` toggles between analytical and finite-difference
  evaluations; both honor active (PHVA) subspaces. The MLIP backend may return only the active block when
  frozen atoms are present.
  For Hessian evaluation modes, see {ref}`hessian-evaluation`.
- **Dimer mode details**:
 - The Hessian Guided Dimer stage periodically refreshes the dimer direction by evaluating an exact
  Hessian (active subspace, TR-projected) and prefers `torch.lobpcg` for the lowest
  eigenpair when `root == 0` (falling back to `torch.linalg.eigh`).
 - When enabled (`--flatten`), the flatten loop updates the stored active Hessian via
  Bofill (SR1/MS ↔ PSB blend; toggle via `hessian_dimer.flatten_loop_bofill`) using
  displacements Δx and gradient differences Δg. Each loop estimates imaginary modes, flattens
  once, refreshes the dimer direction, runs a dimer+LBFGS micro-segment, and (optionally)
  performs a Bofill update. Once only one imaginary mode remains, a final exact Hessian is
  computed for frequency analysis.
 - If `root != 0`, that root seeds only the initial dimer direction; subsequent refreshes
  follow the most negative mode (`root = 0`).
- **RS-I-RFO mode**: runs the RS-I-RFO optimizer with optional Hessian reference files,
  R+S splitting safeguards, and micro-cycle controls defined in the `rsirfo` YAML section.
  When `--flatten` is enabled and more than one imaginary mode remains after
  convergence, the workflow flattens extra modes and reruns RS-I-RFO until only one
  imaginary mode remains or the flatten iteration cap is reached.
- **Mode export & conversion**: all detected imaginary modes are written to `vib/imag_*_trj.xyz`
  and mirrored to `.pdb` when the input was PDB and conversion is enabled. The optimization
  trajectory and final geometry are also converted to PDB via the input template when `--dump`;
  Gaussian templates receive a `.gjf` companion for the final geometry only.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Net charge. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB inputs or XYZ/GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers INT` | MLIP predictor parallelism (workers > 1 disables analytic Hessians). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1` |
| `--workers-per-node INT` | Workers per node, forwarded to the parallel predictor. | `1` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--freeze-links/--no-freeze-links` | PDB-only. Freeze parents of link hydrogens (merged into `geom.freeze_atoms`). See [extract](extract.md) for link-hydrogen details. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | Optimizer preset: `grad` (`dimer`) or `hess` (`rsirfo`). Aliases `dimer`/`rsirfo` are accepted. For the full subcommand-dependent table (`opt` uses L-BFGS/RFO; `tsopt` uses Dimer/RS-I-RFO), see {ref}`opt-mode-semantics`. | `hess` |
| `--dump/--no-dump` | Dump trajectories. | `False` |
| `-o, --out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--flatten/--no-flatten` | Enable the surplus-imaginary-mode flattening loop (`False` forces `flatten_max_iter=0`). After TS optimization converges, iteratively flattens surplus negative-eigenvalue modes of the Hessian matrix until only one imaginary frequency remains (or the iteration cap is reached). Applies to both dimer (dimer loop) and RS-I-RFO (post-convergence). | `False` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB or Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--dry-run/--no-dry-run` | Validate inputs/config and print the execution plan without running TS optimization. | `False` |

(flatten-precedence-caveat)=
### `--flatten` precedence caveat

```{note}
**`--flatten` is disabled by default (precedence caveat).** Although `defaults.py` defines `flatten_max_iter: 50` (and the YAML table below shows `flatten_max_iter: 50`), the CLI initializer forces `flatten_max_iter = 0` unless `--flatten` is **explicitly** passed on the command line. In other words, the effective value is:

- CLI `--flatten` **not** passed → `flatten_max_iter = 0` (surplus-mode cleanup disabled). The YAML value of 50 is **ignored**.
- CLI `--flatten` passed → the YAML / `defaults.py` value applies (default `flatten_max_iter = 50`); you can override via YAML `hessian_dimer.flatten_max_iter` or `rsirfo.flatten_max_iter`.

If your TS candidate has multiple imaginary frequencies, add `--flatten` to enable the surplus-mode cleanup loop.
```

## Outputs
```
out_dir/ (default:./result_tsopt/)
├─ final_geometry.xyz # Always written
├─ final_geometry.pdb # When the input was PDB (conversion enabled)
├─ final_geometry.gjf # When the input was Gaussian (conversion enabled)
├─ optimization_all_trj.xyz # Dimer-mode dump when --dump is True
├─ optimization_all.pdb # Dimer-mode companion for PDB inputs (conversion enabled, --dump)
├─ optimization_trj.xyz # RSIRFO-mode trajectory when --dump is True
├─ optimization.pdb # RSIRFO-mode PDB companion when conversion is enabled and --dump is True
├─ vib/
│ ├─ imag_±XXXX.Xcm-1_trj.xyz
│ └─ imag_±XXXX.Xcm-1.pdb
└─.dimer_mode.dat # Dimer-mode orientation seed
```

## Exit codes

See {ref}`exit-codes` in CLI Conventions.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.
- Imaginary-frequency **detection** threshold defaults to 5.0 cm⁻¹ (configurable via
  `hessian_dimer.neg_freq_thresh_cm`); frequencies with magnitudes below this threshold are not counted as imaginary. The selected `root` controls which vibrational mode is followed during optimization. See {ref}`imaginary-mode-thresholds` in the glossary for the 5 cm⁻¹ detection threshold vs 100 cm⁻¹ quality gate.
- Use `--opt-mode` to choose the algorithm workflow directly (`rsirfo` by default), instead of
  manually editing YAML mode mappings.
- PHVA translation/rotation projection follows the same implementation as `freq`, while reducing
  memory usage and preserving correct active-space eigenvectors.
- See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

Shared sections reuse
[YAML Reference](yaml-reference.md). Adjust only the values you need to change.

```{note}
**Reference duplication.** The YAML keys for `geom`, `calc`, `opt`, `hessian_dimer`, and `rsirfo` listed below mirror the canonical definitions in [YAML Reference](yaml-reference.md). When the two pages disagree, the canonical [YAML Reference](yaml-reference.md) entries (and `pdb2reaction/defaults.py`) take precedence; the inline appendix on this page is reproduced only for `tsopt`-specific defaults (e.g. `out_dir: ./result_tsopt/`, the `--flatten` interaction documented above). Note that `flatten_max_iter` is forced to `0` by the CLI initializer unless `--flatten` is passed, regardless of the value shown in the inline YAML.
```

### Shared configuration (common to both modes)

`geom` and `calc` keys are unchanged from the canonical definitions; see [`geom`](yaml-reference.md#geom) and [`calc`](yaml-reference.md#calc) in the YAML Reference.

The `opt` block uses the same keys as [`opt`](yaml-reference.md#opt), with these `tsopt`-specific defaults:

```yaml
opt:
 thresh: baker # tsopt default (vs. `gau` for `opt`)
 out_dir: ./result_tsopt/ # tsopt default (vs. `./result_opt/` for `opt`)
```

```{note}
**Energy plateau fallback (new in v0.3.5).** The RS-I-RFO optimizer honours the
shared `energy_plateau` setting: if the energy range (max − min) over the last
`energy_plateau_window` (default 50) steps falls below `energy_plateau_thresh`
(default `1×10⁻⁴ au ≈ 0.06 kcal/mol`), the run is declared converged. This is
especially useful for large TS systems where the MLIP force noise floor
(~4×10⁻⁴ au) exceeds the `baker` max_force threshold (3×10⁻⁴ au), making the
force criterion unreachable even though the energy has plainly flattened. Set
`energy_plateau: false` to disable.
```

### Dimer mode (`--opt-mode grad`)

Used with `--opt-mode grad` (Hessian Guided Dimer + LBFGS translation).

The full `hessian_dimer` block (including the inner `dimer:` and its nested `lbfgs:`) is documented in [`hessian_dimer`](yaml-reference.md#hessian_dimer). The inner `lbfgs:` inherits from the [`lbfgs`](yaml-reference.md#lbfgs) section, with this `tsopt`-specific override:

```yaml
hessian_dimer:
 dimer:
   lbfgs:
     out_dir: ./result_tsopt/ # tsopt override (defaults.py value is ./result_opt/)
```

Recall that `flatten_max_iter` is forced to `0` by the CLI initializer unless `--flatten` is passed (see {ref}`flatten-precedence-caveat` above).

### RS-I-RFO mode (`--opt-mode hess`, default)

Used with `--opt-mode hess` (RS-I-RFO, the default).

The full `rsirfo` block is documented in [`rsirfo`](yaml-reference.md#rsirfo) (which also inherits trust-region and Hessian-update keys from [`rfo`](yaml-reference.md#rfo)). The `tsopt`-specific overrides are:

```yaml
rsirfo:
 trust_max: 0.10 # maximum trust radius (bohr); reduced from 0.20 in v0.3.5 to damp near-saddle oscillations on large systems
 out_dir: ./result_tsopt/ # tsopt override (defaults.py value is ./result_opt/)
 hessian_recalc: 500 # rebuild exact Hessian every N macro steps; lower to 50-200 if TS convergence is slow (see tip below)
```

```{tip}
Set `rsirfo.track_mode_by_overlap: true` if the TS mode switches root during optimization (e.g., when multiple imaginary frequencies are present).
```

```{tip}
If TS convergence is slow or the TS mode is lost during optimization, try lowering `hessian_recalc` (e.g., to 50--200) in the `rsirfo` section. More frequent exact Hessian recalculations improve robustness at the cost of additional Hessian evaluations.
```

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing

- [path-search](path-search.md) — MEP search that identifies TS candidates (HEI)
- [irc](irc.md) — Trace the reaction path from an optimized TS
- [freq](freq.md) — Full vibrational analysis and thermochemistry (imaginary-frequency check is already included in `tsopt`)
- [all](all.md) — End-to-end workflow that chains extraction → MEP → tsopt → IRC (→ optional freq/DFT)
- [YAML Reference](yaml-reference.md) — Full `hessian_dimer` (Hessian Guided Dimer) and `rsirfo` configuration options
- [Glossary](glossary.md) — Definitions of TS, Dimer, RS-I-RFO, Hessian
