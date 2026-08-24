# `tsopt`

`pdb2reaction tsopt` refines a transition-state (TS) *candidate* into an optimized first-order saddle point, with a built-in imaginary-frequency check. The candidate can be the highest-energy image (HEI) from `path-opt` / `path-search`, or a user-supplied structure.

Pick the optimizer with `--opt-mode`. The default `hess` mode is **RS-P-RFO**
(Restricted-Step Partitioned Rational Function Optimization, Banerjee) and uses
a full Hessian. The `grad` mode is the **Hessian-Guided Dimer** alternative
when full-Hessian recomputation is prohibitive or a separate optimization route
is needed. Compare the methods on the actual TS seed and validate the result
with one imaginary mode and IRC connectivity. Enable `--flatten` (disabled by
default) when the candidate has multiple imaginary frequencies and surplus-mode
cleanup is needed.

`tsopt` always sets `reject_uphill: false` for its RFO-family and Dimer
optimizers, including after YAML overrides. A saddle search must be able to
raise the physical energy along the reaction mode. The
`--reject-uphill/--no-reject-uphill` toggle belongs only to minimum
optimization (`opt` and post-IRC endpoint re-optimization in `all`).

At optimizer termination, `tsopt` retains the final geometry. Terminal exact PHVA runs only after numerical convergence; a non-converged or stalled run stops without PHVA. A PHVA failure is recorded without discarding the structure or fabricating frequencies. Numerical optimizer status is independent of saddle order: first-order certification still requires **exactly one** imaginary frequency, the intended displacement, and correct [`irc`](irc.md) connectivity. A separate [`freq`](freq.md) run is needed only for full vibrational analysis or thermochemistry.

### Terminal outcomes and fatal errors

| Condition | `tsopt` artifacts | Composite `all` behavior |
| --- | --- | --- |
| Convergence criteria unmet, explicit cycle limit reached, or opt-in energy plateau | Retain the final geometry and trajectory; skip terminal PHVA | Register the TS result and stop before IRC |
| Terminal PHVA fails | Retain the geometry and set `hessian_status: failed` with the error; do not invent frequencies | Stop before IRC after artifact registration |
| Invalid input/geometry or an unrecoverable optimizer exception such as `ZeroStepLength` / `OptimizationError` | Follow the structured error-envelope path; only files already written are retained on a best-effort basis | Abort the stage rather than relabelling it as ordinary non-convergence |


If you need a TS guess first, run [`path-opt`](path-opt.md) (two structures) or [`path-search`](path-search.md) (two or more structures), then optimize the HEI with `tsopt` → `irc`. For XYZ / GJF inputs, `--ref-pdb` supplies a reference PDB/mmCIF topology while keeping the XYZ coordinates, enabling format-aware PDB / CIF / GJF companion output.

`--ref-mode` is an advanced/internal handoff, not a normal requirement for
standalone `tsopt`. It accepts one or more atom-order-matched Cartesian 3N
candidate directions from `.npz`, `.npy`, or whitespace text. The `all`
workflow supplies CPU/file-cached MEP tangent candidates automatically for
Hessian TS optimizers; legacy paths without readable energies use normalized
secants. Dimer does not consume `--ref-mode`. With
`all --no-tsopt-from-mep-tan`, cache creation/use is disabled and TSOPT selects
its initial root from the initial-structure Hessian modes.

The reference direction guides Hessian-root identity and overlap tracking; it
is **not** an initial Hessian replacement and does not make a failed TS search
successful. Terminal exact PHVA remains authoritative for saddle order.
`n_imag = 0` is `no_imaginary`; `n_imag > 1` is `higher_order`. Neither state
rewrites a numerically converged optimizer result as numerical non-convergence.
A higher-order result may be used only for warning-labelled diagnostic IRC by
`all`, never as first-order TS certification.

`--flatten` is a separate, explicit cleanup for surplus imaginary modes. It can
remove extra negative directions but cannot create a missing reaction mode.

> **Naming note:** the CLI accepts `grad` / `dimer` (Dimer), `hess` / `rsprfo` (RS-P-RFO, default), and `rsirfo` (RS-I-RFO) / `trim` (TRIM). In YAML, use the top-level `hessian_dimer:` (Dimer) block, or the `rsirfo:` block (shared by RS-P-RFO, RS-I-RFO, and TRIM), directly.

## Two routes to a TS candidate

`tsopt` refines a candidate you already have; there are two complementary ways to *build* that candidate first. Pick the route that matches the information you have.

| Route | Subcommand | Use when | What it does |
| --- | --- | --- | --- |
| (a) MEP / path search | [`path-search`](path-search.md) | You have both endpoints (reactant **and** product) and want the TS bracketed automatically | Recursive minimum-energy-path search (GSM / DMF) with bond-change detection; it auto-segments a multi-step path, refines each reactive segment, and returns the highest-energy image per segment (`hei_seg_NN.xyz`) |
| (b) Distance-restrained scan | [`scan`](scan.md) | You have only the reactant, or want to drive a specific reacting distance directly | Harmonic distance restraints, `E = ½k(r − target)²`, drive each reacting distance with full relaxation, advancing the system toward a TS candidate |

There is no `opt --restraint` flag: `opt` restrains distances with `--dist-freeze` (harmonic, `--bias-k`) rather than driving them, and the distance-driven build-up route is `scan` (which can relax the endpoints around the driven path with `--preopt` / `--endopt`). Feed the candidate from either route into `tsopt → freq → irc` to optimize and validate it.

## Examples

Default RS-P-RFO optimization of a PDB candidate:

```bash
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

Dimer mode with analytical Hessian (VRAM permitting):

```bash
# Dimer mode with analytical Hessian (VRAM permitting)
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
    --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

RS-P-RFO mode driven by YAML overrides:

```bash
# RS-P-RFO mode driven by YAML overrides
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
    --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
```

RS-P-RFO mode with surplus-imaginary-mode flattening enabled:

```bash
# RS-P-RFO mode with surplus-imaginary-mode flattening enabled
pdb2reaction tsopt -i ts_cand.pdb -q 0 -m 1 \
    --opt-mode hess --flatten --out-dir ./result_tsopt_flatten
```

Add `--dump` to keep the full optimization trajectory for inspection.

## Workflow

- **Charge / spin** are resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>`).
- **Geometry loading + freeze-links** — structures are read through `pysisyphus.helpers.geom_loader`. When `--freeze-links` is active, cap-hydrogen parent atoms are automatically frozen (see {ref}`Cap hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- **MLIP Hessians** (default UMA) — `--hessian-calc-mode` toggles analytical vs finite-difference; both honor active (PHVA) subspaces. The MLIP backend may return only the active block when frozen atoms are present. See {ref}`hessian-evaluation` for the full Hessian-evaluation matrix.
- **Dimer mode** — the Hessian-Guided Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian in the active subspace. It uses the fixed constrained rigid-mode treatment, which removes only full-system rigid motions compatible with frozen anchors. Every stored, rotated, and trial orientation has frozen Cartesian components set to zero, and every off-center force evaluation retains the central image's frozen coordinates exactly. The lowest eigenpair uses `torch.lobpcg` when `root == 0`, falling back to `torch.linalg.eigh`. With `--flatten`, the active Hessian is updated via a Bofill update (an SR1/MS ↔ PSB blend; toggle via `hessian_dimer.flatten_loop_bofill`) using displacements Δx and gradient differences Δg. Each flatten loop:
  - estimates imaginary modes, flattens once, and refreshes the dimer direction;
  - runs a Dimer + L-BFGS micro-segment;
  - optionally performs a Bofill update.

  At termination, one exact PHVA is produced after numerical convergence; a non-converged or stalled run retains the final geometry and skips PHVA. If `root != 0`, that root seeds only the initial dimer direction; subsequent refreshes follow the most negative mode (`root = 0`).
- **RS-I-RFO mode** — runs the RS-I-RFO optimizer with optional Hessian reference files, R+S splitting safeguards, and micro-cycle controls defined in the `rsirfo` YAML section. With `--flatten`, when more than one imaginary mode remains after convergence the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten-iteration cap is reached.
- **Mode export + conversion** — imaginary modes smaller than the configured magnitude threshold (5 cm⁻¹ by default) are ignored; the remaining modes are written to `vib/imag_*_trj.xyz`. With conversion enabled, PDB inputs receive `.pdb` companions and mmCIF/oversized-PDB bridge inputs receive `.pdb` plus `.cif`; Gaussian templates receive a `.gjf` companion for the final geometry only.

## Outputs

Validate a run from `result.json`, the final geometry in `final_geometry.*`, and the `vib/imag_*` modes (expect exactly one for a valid TS).

```text
out_dir/   (default: ./result_tsopt/)
├─ final_geometry.xyz              # Always written
├─ final_geometry.pdb              # PDB/mmCIF topology input (conversion enabled)
├─ final_geometry.cif              # mmCIF/oversized-PDB bridge input
├─ final_geometry.gjf              # When the input was Gaussian (conversion enabled)
├─ optimization_all_trj.xyz        # Dimer-mode dump (--dump)
├─ optimization_all.pdb            # Dimer-mode PDB companion (--dump, topology input)
├─ optimization_all.cif            # Bridge-input companion with original IDs
├─ optimization_trj.xyz            # RS-P-RFO/RS-I-RFO/TRIM trajectory (--dump)
├─ optimization.pdb                # RFO-mode PDB companion (--dump)
├─ optimization.cif                # Bridge-input companion with original IDs
├─ vib/
│  ├─ imag_±XXXX.Xcm-1_trj.xyz
│  ├─ imag_±XXXX.Xcm-1.pdb
│  └─ imag_±XXXX.Xcm-1.cif         # Bridge input
└─ .dimer_mode.dat                 # Dimer-mode orientation seed
```

Exit codes: see {ref}`exit-codes` in CLI Conventions.

## CLI options

Command form:

```bash
pdb2reaction tsopt -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l 'RES:Q,...'] [-m 2S+1] \
    [-b uma|orb|mace|aimnet2] [--opt-mode grad|hess|dimer|rsirfo|trim|rsprfo] [--flatten / --no-flatten] \
    [--freeze-links / --no-freeze-links] [--max-cycles N] [--thresh PRESET] \
    [--hessian-calc-mode Analytical|FiniteDifference] \
    [--convert-files / --no-convert-files] [--ref-pdb FILE]
```

`pdb2reaction tsopt --help` shows core options; `pdb2reaction tsopt --help-advanced` shows the full option list. For full input-file requirements (hydrogens, element columns, atom-order parity, charge specification), see [CLI Conventions](cli-conventions.md).

The tables below cover the options that need explanation. The full flag list is in the generated [command reference](reference/commands/index.md).

| Option | Description | Default |
| --- | --- | --- |
| **Input & charge** | | |
| `-i, --input PATH` | One geometry (`.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf`). Extract a desired trajectory frame to `.xyz` first. | Required |
| `-q, --charge INT` | Net charge. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB/mmCIF inputs or XYZ / GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template / derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g. `-1`) for the total ligand charge, or a per-residue mapping (e.g. `GPP:-3,SAM:1`) that derives the total from PDB/mmCIF residue metadata. Used when `-q` is omitted. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--ref-pdb FILE` | Reference PDB/mmCIF topology when the input is XYZ / GJF (keeps XYZ coordinates). | _None_ |
| **Backend & compute** | | |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--workers INT`, `--workers-per-node INT` | UMA predictor parallelism. `workers > 1` cannot be combined with an explicit analytical Hessian request; use `workers = 1` or finite differences. See {ref}`workers-analytical-error`. | `1`, `1` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| **Active-region freezing** | | |
| `--freeze-links / --no-freeze-links` | PDB/mmCIF input (or XYZ/GJF with `--ref-pdb`). Freeze parents of cap hydrogens (merged into `geom.freeze_atoms`). | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g. `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| **TS optimizer & mode** | | |
| `--opt-mode TEXT` | TS optimizer preset (Choice: `grad` / `hess` / `dimer` / `rsirfo` / `trim` / `rsprfo`). `grad` and `dimer` → Hessian-Guided Dimer; `hess` and `rsprfo` → RS-P-RFO (Banerjee, default, non-microiter); `rsirfo` → RS-I-RFO; `trim` → TRIM (Helgaker, non-microiter). On `opt`, the same `grad` token picks L-BFGS minimization instead — see {ref}`opt-mode-semantics`. | `hess` |
| `--ref-mode PATH` | Advanced/internal Cartesian reference candidate(s) from `.npz`, `.npy`, or whitespace text (one 3N vector or a 2-D candidate table). Guides Hessian-root identity/overlap; does not replace the Hessian and is unsupported by Dimer. `all` supplies it from the MEP for Hessian TS optimizers. | _None_ |
| `--flatten / --no-flatten` | Enable surplus-imaginary-mode flattening for Dimer and the RS-P-RFO / RS-I-RFO / TRIM Hessian family. `--ref-mode` identifies which negative mode must be retained but does not enable flattening by itself. | `False` |
| `--coord-type TEXT` | Optimization coordinate system (`cart` / `redund` / `dlc` / `tric`). `cart` is the default. `dlc` changes the conditioning, but neither representation is uniformly faster or more reliable; compare them on the problematic seed. Hessian-based `tsopt` modes support all four, while `path-opt` / `path-search` accept only `cart` / `dlc`. | `cart` |
| `--precision [fp32\|fp64]` | MLIP backend precision, routed to the backend-native kwarg (UMA `precision` / ORB `precision` / MACE `default_dtype`; `aimnet2`: `fp32` no-op, `fp64` rejected). Compare supported settings on the target system; see [Reproducibility](reproducibility.md#choosing-precision-by-backend-and-purpose). | per backend (uma `fp32`; orb, mace `fp64`) |
| **Thresholds & cycles** | | |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `100000` |
| **Output & config** | | |
| `-o, --out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--convert-files / --no-convert-files` | Toggle XYZ/TRJ → PDB/CIF/GJF companions according to the input topology/template. | `True` |
| `--dump / --no-dump` | Dump trajectories. | `False` |
| `--out-json / --no-out-json` | Write a machine-readable `result.json` to `out_dir`. Schema: [JSON Output Schema](json-output.md). | `False` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config / --no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--dry-run / --no-dry-run` | Validate inputs / config and print the execution plan without running TS optimization. | `False` |

(flatten-precedence-caveat)=
### `--flatten` precedence caveat

```{note}
**`--flatten` is disabled by default (precedence caveat).** Although `defaults.py` defines `flatten_max_iter: 50`, the CLI initializer seeds `flatten_max_iter = 0`. Effective resolution:

- CLI `--flatten` **not** passed → `flatten_max_iter = 0` unless **explicitly set in YAML** via `hessian_dimer.flatten_max_iter` (the flatten counter is read only from the `hessian_dimer` block for both Dimer and RS-I-RFO paths). The `defaults.py` value of 50 is ignored.
- CLI `--flatten` passed → the YAML / `defaults.py` value applies (default `flatten_max_iter = 50`); you can still override via YAML.
- CLI `--no-flatten` passed → `flatten_max_iter = 0`, overriding YAML.

If your TS candidate has multiple imaginary frequencies, add `--flatten` to enable the surplus-mode cleanup loop.

If TS optimization still fails from a path HEI, there are two distinct retries:

1. Add `--flatten` when the candidate retains surplus imaginary modes.
2. In the `all` workflow, use `--refine-path` to run recursive
   `path-search` and obtain a better-resolved HEI before TSOPT.

The second retry is intentionally not the default: a poor or noisy path can be
split into unnecessary elementary segments, multiplying MEP, TSOPT, IRC, and
frequency work. Inspect the unrefined MEP first and enable recursive refinement
only when the coarse HEI is the likely cause.
```

### Wrong imaginary-mode count after optimization

A true first-order saddle has **exactly one** imaginary frequency, and its mode displaces along the reaction coordinate. If `tsopt` instead reports a spurious second small imaginary mode, or no dominant reaction mode, escalate the following levers — they are complementary, so you can combine them:

| Lever | Flag | Effect |
| --- | --- | --- |
| Compare precision | `--precision fp32|fp64` | Numerical behavior is backend/model/system dependent; AIMNet2 rejects fp64, and neither setting removes genuine negative curvature. |
| Internal coordinates | `--coord-type dlc` | Changes the optimization conditioning. Benchmark `cart` and `dlc` on the problematic seed because neither is uniformly faster or more reliable. |
| Flatten small modes | `--flatten` | Runs an extra-imaginary-mode flattening loop (`grad`: dimer loop; `hess`: post-RS-P-RFO step); `--no-flatten` forces `flatten_max_iter = 0`. |

Inspect the mode displacements and optimizer stop reason, then retry an appropriate supported precision/coordinate setting and use `--flatten` only for surplus modes. For example:

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q -1 -m 1 \
    --precision fp64 --coord-type dlc --flatten -o result_tsopt
```

See also [Common Error Recipes → Convergence and post-processing failures](recipes-common-errors.md).

### Reading a barrier scanned from the product side

If the `scan` (or path) that produced this TS candidate started from the **product**, the raw barrier it reports is the **reverse** barrier, `E(TS) − E(product)`. The forward barrier you usually want is computed from the reactant:

| You ran | Forward barrier |
| --- | --- |
| A product-start scan | `E(TS) − E(reactant)` — **not** the raw product-start number |

This is a read-time interpretation, not a flag. Always confirm which endpoint the scan started from before quoting a barrier, especially when the workflow was seeded from a crystallographic product complex. See also [`scan` → Scan direction and barrier sign](scan.md#scan-direction-and-barrier-sign).

### Controlled mutant-vs-WT comparison

Do not confuse the MEP input contract with a cross-variant comparison. Within
each R→IM→P path, all structures must contain the same atoms in the same order.
A real WT→mutant substitution may change residue identity and atom count, so
raw total energies from WT and mutant must not be subtracted directly. Compare
the activation energy or free energy computed within each system instead:

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`.

- Match the selected residue **positions** and cluster-boundary/cap policy as
  closely as chemically meaningful, with the intended mutation as the only
  designed composition change. Compare independent radius-based extractions,
  because a boundary residue can enter one cluster but not the other.
- Keep protonation, charge-assignment convention, backend/model, precision,
  constraints, and thermochemistry settings controlled. The verified total
  charge itself may legitimately differ if the mutation changes protonation or
  formal charge; do not force equal `-q` values to make the comparison look symmetric.
- For two mechanisms of the **same composition**, use one common atom set and
  ordering across the compared paths.

```bash
# Each path is internally atom-consistent; the two clusters use matched boundaries.
pdb2reaction all -i wt_cluster.pdb     -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_wt
pdb2reaction all -i mutant_cluster.pdb -l 'GPP:-3,SAM:1' --tsopt --thermo -o result_mutant
```

## YAML configuration

Shared sections reuse [YAML Reference](yaml-reference.md); adjust only the values you need to change. `geom` and `calc` are unchanged from canonical — see [`geom`](yaml-reference.md#geom) and [`calc`](yaml-reference.md#calc). The `opt` block uses the same keys as [`opt`](yaml-reference.md#opt) with these `tsopt`-specific defaults:

```yaml
opt:
  thresh: baker                # tsopt default (vs. `gau` for `opt`)
  out_dir: ./result_tsopt/     # tsopt default (vs. `./result_opt/` for `opt`)
```

```{note}
**Energy-plateau stop (opt-in, default off).** Hessian-family TS optimizers
(RS-P-RFO, RS-I-RFO, TRIM, and Dimer) honor the shared `energy_plateau` setting, which
`--stop-plateau` turns on. An energy range below `--stop-plateau-thresh`
(default `1×10⁻⁴ au` over the last 50 steps) stops the search as `stalled` and
skips terminal PHVA, as does reaching `max_cycles` without convergence. This can save cycles when a
backend/model/system-specific force floor prevents the selected force threshold
from being reached. It is off unless you ask for it, because a TS search that
stops on a flat energy typically still carries extra imaginary modes.
```

### Dimer mode (`--opt-mode grad`)

Used with `--opt-mode grad` (Hessian-Guided Dimer + L-BFGS translation). The full `hessian_dimer` block, including its sibling `dimer:` and `lbfgs:` sections, is documented in [`hessian_dimer`](yaml-reference.md#hessian_dimer). The option names under `hessian_dimer.lbfgs` follow the [`lbfgs`](yaml-reference.md#lbfgs) schema, and `tsopt` reads their values from this sibling section:

```yaml
hessian_dimer:
  lbfgs:
    out_dir: ./result_tsopt/   # tsopt override (defaults.py value is ./result_opt/)
```

### RS-P-RFO / RS-I-RFO mode (`--opt-mode hess`, default → RS-P-RFO)

Used with `--opt-mode hess` (RS-P-RFO, the default; `rsirfo` selects RS-I-RFO and `trim` selects TRIM — all three share this block). The full `rsirfo` block is documented in [`rsirfo`](yaml-reference.md#rsirfo) (which inherits trust-region and Hessian-update keys from [`rfo`](yaml-reference.md#rfo)). `tsopt`-specific overrides:

```yaml
rsirfo:
  trust_max: 0.10              # maximum trust radius (bohr)
  out_dir: ./result_tsopt/     # tsopt override (defaults.py value is ./result_opt/)
  hessian_recalc: 500          # rebuild exact Hessian every N macro steps
  saddle_recovery_max_cycles: 0      # automatic n_imag=0 recovery is disabled
```

```{tip}
Set `rsirfo.track_mode_by_overlap: true` if the TS mode switches root during optimization (e.g. when multiple imaginary frequencies are present). If TS convergence is slow, lowering `hessian_recalc` (e.g. to 50–200) gives more frequent exact Hessian updates at the cost of additional evaluations. Automatic curvature recovery can be enabled through YAML by setting `saddle_recovery_max_cycles` above zero, but it is not part of the default search.
```

## Notes

- Imaginary frequencies smaller than the configured threshold (5 cm⁻¹ by default) are ignored consistently by final TS validation, mode-file output, and flattening.
- Hessian-family optimizers follow exactly one root for a first-order TS. Set it as a one-item YAML list (for example, `rsirfo.roots: [0]`); empty or multi-root lists are rejected. Dimer uses the separate singular `hessian_dimer.root` key (default `0`). `tsopt` has no `--root` CLI flag, unlike [`irc`](irc.md).
- Use `--opt-mode` to choose the algorithm directly (`rsprfo` by default) rather than editing YAML mode mappings.
- Dimer orientation, rotation forces, flattening, and final exact PHVA validation use the same constrained projector as `freq`. The Dimer rebuilds this basis whenever its central image changes. It never subtracts translations of the active fragment unless they are actual rigid null directions compatible with every frozen anchor. Hessian RFO optimization itself operates on the active-DOF Cartesian Hessian without this projection. See [Frozen Atoms](freeze-atoms.md#rigid-modes-with-frozen-boundaries).
- See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

## See Also

[Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [path-search](path-search.md) · [irc](irc.md) · [freq](freq.md) · [all](all.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
