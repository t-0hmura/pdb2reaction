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

After convergence, `tsopt` performs a final Hessian calculation and imaginary-frequency check automatically — a validated TS should show **exactly one** imaginary frequency. A separate [`freq`](freq.md) run is only needed for full vibrational analysis or thermochemistry. Always confirm endpoint connectivity with [`irc`](irc.md).

If you need a TS guess first, run [`path-opt`](path-opt.md) (two structures) or [`path-search`](path-search.md) (two or more structures), then optimize the HEI with `tsopt` → `irc`. For XYZ / GJF inputs, `--ref-pdb` supplies a reference PDB/mmCIF topology while keeping the XYZ coordinates, enabling format-aware PDB / CIF / GJF companion output.

`--ref-mode` is an advanced/internal handoff, not a normal requirement for
standalone `tsopt`. The `all` workflow supplies it automatically from the MEP
as the standard energy-upwinding Cartesian tangent at each HEI; legacy paths
without readable energies use a normalized secant bisector. An expert
standalone run may pass
`--ref-mode PATH` only when it deliberately has a non-zero Cartesian 3N reaction direction
from an external path calculation. The mode selects the initial Hessian root,
using the softest curvature among roots whose overlap is within 90% of the
maximum when a discrete tangent spans several near-tied modes. It is then
transported by overlap across the complete positive/negative Hessian spectrum
as it rotates. Before that target has ever acquired negative curvature,
`n_imag=0` recovery retains the complete path vector instead of collapsing it
to one positive-curvature normal mode; subsequent recovery uses the transported
mode. Neither an initial raw negative root nor a transient negative
quasi-Newton root arms mode-loss rollback until exact PHVA or an explicit
recovery crossing confirms the physical target. The reference therefore
supplies both the recovery and exact-identity direction until PHVA first
confirms the requested negative curvature; only after that crossing does
continuous eigenmode transport become authoritative. If exact validation
finds `n_imag=0`, a structure with no path information cannot uniquely identify
which nearby saddle is intended. Bounded recovery checks an exact physical
Hessian every 50 steps for up to 200 steps, so it can resume as soon as the
physical mode crosses negative even when the quasi-Newton model lags behind.
Exact validation matches that authoritative mode
against all physical normal modes; a different imaginary mode cannot substitute
for the tracked mode if it has become positive. Persistent higher-order candidates stop after three exact
terminal checks. When `--flatten` is explicitly enabled, the subsequent cleanup
removes non-path negative modes within the configured iteration budget. Because
individual eigenvectors can exchange identity inside a multi-negative subspace,
each higher-order exact check re-anchors the mode to retain against the immutable
MEP tangent before flattening; after the first
physical crossing, order-zero and first-order checks use overlap transport so
genuine mode rotation is preserved. A retry keeps the standard three-check exact window so the retained
saddle direction can be re-established after the orthogonal displacement.
For every cleanup retry, if the energy-selected sign does not already give
the requested first-order saddle, the opposite sign is also optimized and the
branch that best retains the requested mode while approaching first order is
kept; displacement energy alone does not select the branch.
Each retry inherits that latest PHVA-identified path mode rather than
selecting an arbitrary negative root after the flatten displacement. A separate
snapshot keeps the exact PHVA Cartesian mode distinct
while the raw Hessian eigenvector continues to update for per-cycle root
tracking. If the first path-guided
Hessian run instead stops at order
zero/one without a verified saddle, `tsopt` makes bounded same-optimizer
restarts from ±0.10 Å and then ±0.20 Å along that mode. If a kinked discrete
tangent finds only unrelated modes, the same bounded shells are tried along the
soft, path-correlated Hessian root selected at the original HEI. The source and
outcome of at most eight attempts are recorded in `result.json`; if all fail,
the original result is retained.

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

  Once only one imaginary mode remains, a final exact Hessian is computed for frequency analysis. If `root != 0`, that root seeds only the initial dimer direction; subsequent refreshes follow the most negative mode (`root = 0`).
- **RS-I-RFO mode** — runs the RS-I-RFO optimizer with optional Hessian reference files, R+S splitting safeguards, and micro-cycle controls defined in the `rsirfo` YAML section. With `--flatten`, when more than one imaginary mode remains after convergence the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten-iteration cap is reached.
- **Mode export + conversion** — all detected imaginary modes are written to `vib/imag_*_trj.xyz`. With conversion enabled, PDB inputs receive `.pdb` companions and mmCIF/oversized-PDB bridge inputs receive `.pdb` plus `.cif`; Gaussian templates receive a `.gjf` companion for the final geometry only.

## Outputs

Validate a run by opening `final_geometry.*` (the optimized saddle point) and the `vib/imag_*` modes (expect exactly one for a valid TS).

```text
out_dir/   (default: ./result_tsopt/)
├─ final_geometry.xyz              # Always written
├─ final_geometry.pdb              # PDB/mmCIF topology input (conversion enabled)
├─ final_geometry.cif              # mmCIF/oversized-PDB bridge input
├─ final_geometry.gjf              # When the input was Gaussian (conversion enabled)
├─ optimization_all_trj.xyz        # Dimer-mode dump (--dump)
├─ optimization_all.pdb            # Dimer-mode PDB companion (--dump, topology input)
├─ optimization_all.cif            # Bridge-input companion with original IDs
├─ optimization_trj.xyz            # RSIRFO-mode trajectory (--dump)
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
    [-b uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [--opt-mode grad|hess|dimer|rsirfo|trim|rsprfo] [--flatten / --no-flatten] \
    [--freeze-links / --no-freeze-links] [--max-cycles N] [--thresh PRESET] \
    [--hessian-calc-mode Analytical|FiniteDifference] \
    [--convert-files / --no-convert-files] [--ref-pdb FILE]
```

`pdb2reaction tsopt --help` shows core options; `pdb2reaction tsopt --help-advanced` shows the full option list. For full input-file requirements (hydrogens, element columns, atom-order parity, charge specification), see [CLI Conventions](cli-conventions.md).

The tables below cover the options that need explanation. The full flag list is in the generated [command reference](reference/commands/index.md).

| Option | Description | Default |
| --- | --- | --- |
| **Input & charge** | | |
| `-i, --input PATH` | Structure file accepted by the input bridge (`.pdb` / `.cif` / `.mmcif` / `.xyz` / `.gjf` / `.trj`). | Required |
| `-q, --charge INT` | Net charge. Required unless a `.gjf` template or `--ligand-charge/-l` (PDB/mmCIF inputs or XYZ / GJF with `--ref-pdb`) supplies it. Overrides `--ligand-charge/-l` when both are set. | Required unless template / derivation applies |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g. `-1`) for the total ligand charge, or a per-residue mapping (e.g. `GPP:-3,SAM:1`) that derives the total from PDB/mmCIF residue metadata. Used when `-q` is omitted. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `.gjf` template value or `1` |
| `--ref-pdb FILE` | Reference PDB topology when the input is XYZ / GJF (keeps XYZ coordinates). | _None_ |
| **Backend & compute** | | |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--workers INT`, `--workers-per-node INT` | UMA predictor parallelism. `workers > 1` cannot be combined with an explicit analytical Hessian request; use `workers = 1` or finite differences. See {ref}`workers-analytical-error`. | `1`, `1` |
| `--hessian-calc-mode CHOICE` | MLIP Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| **Active-region freezing** | | |
| `--freeze-links / --no-freeze-links` | PDB input (or XYZ/GJF with `--ref-pdb`). Freeze parents of cap hydrogens (merged into `geom.freeze_atoms`). | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g. `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| **TS optimizer & mode** | | |
| `--opt-mode TEXT` | TS optimizer preset (Choice: `grad` / `hess` / `dimer` / `rsirfo` / `trim` / `rsprfo`). `grad` and `dimer` → Hessian-Guided Dimer; `hess` and `rsprfo` → RS-P-RFO (Banerjee, default, non-microiter); `rsirfo` → RS-I-RFO; `trim` → TRIM (Helgaker, non-microiter). On `opt`, the same `grad` token picks L-BFGS minimization instead — see {ref}`opt-mode-semantics`. | `hess` |
| `--ref-mode PATH` | Advanced/internal MEP handoff containing a Cartesian 3N direction as whitespace text or `.npy`. `all` supplies it automatically; ordinary standalone runs omit it. Expert use covers external-path root selection, overlap tracking, and `n_imag=0` recovery. | _None_ |
| `--flatten / --no-flatten` | Enable general surplus-imaginary-mode flattening. After TS optimization, iteratively flattens surplus negative-eigenvalue modes until only one imaginary frequency remains (or the iteration cap is reached). Applies to both Dimer (dimer loop) and RS-P-RFO / RS-I-RFO (post-convergence). `--ref-mode` identifies which negative mode must be retained but does not enable flattening by itself. | `False` |
| `--coord-type TEXT` | Optimization coordinate system (`cart` / `redund` / `dlc` / `tric`). `cart` is the default. `dlc` changes the conditioning, but neither representation is uniformly faster or more reliable; compare them on the problematic seed. Hessian-based `tsopt` modes support all four, while `path-opt` / `path-search` accept only `cart` / `dlc`. | `cart` |
| `--precision [fp32\|fp64]` | MLIP backend precision, routed to the backend-native kwarg (UMA `precision` / ORB `precision` / MACE `default_dtype`; `aimnet2`: `fp32` no-op, `fp64` rejected). Compare supported settings on the target system; see [Reproducibility](reproducibility.md#choosing-precision-by-backend-and-purpose). | per backend (uma `fp32`; orb, mace `fp64`) |
| **Thresholds & cycles** | | |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--max-cycles INT` | Macro-cycle cap forwarded to `opt.max_cycles`. | `10000` |
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

A true first-order saddle has **exactly one** imaginary frequency, and its mode displaces along the reaction coordinate (detection cutoff `hessian_dimer.neg_freq_thresh_cm`, default 5 cm⁻¹). If `tsopt` instead reports a spurious second small imaginary mode, or no dominant reaction mode, escalate the following levers — they are complementary, so you can combine them:

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
**Energy-plateau fallback.** Hessian-family TS optimizers (RS-P-RFO,
RS-I-RFO, and TRIM) honor the shared `energy_plateau` setting. An energy
range below `energy_plateau_thresh` (default `1×10⁻⁴ au` over the last 50
steps) triggers exact-Hessian/PHVA terminal validation; it does not bypass the
required first-order saddle and physical convergence checks. This can avoid
wasted cycles when a backend/model/system-specific force floor prevents the
selected force threshold from being reached. Set `energy_plateau: false` to
disable the trigger.
```

### Dimer mode (`--opt-mode grad`)

Used with `--opt-mode grad` (Hessian-Guided Dimer + L-BFGS translation). The full `hessian_dimer` block (including the inner `dimer:` and its nested `lbfgs:`) is documented in [`hessian_dimer`](yaml-reference.md#hessian_dimer); the inner `lbfgs:` inherits from [`lbfgs`](yaml-reference.md#lbfgs), with this `tsopt`-specific override:

```yaml
hessian_dimer:
  dimer:
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
  saddle_recovery_check_interval: 50 # exact PHVA cadence during n_imag=0 recovery
  saddle_recovery_max_cycles: 200    # bounded n_imag=0 recovery cap
```

```{tip}
Set `rsirfo.track_mode_by_overlap: true` if the TS mode switches root during optimization (e.g. when multiple imaginary frequencies are present). If TS convergence is slow or the TS mode is lost, lowering `hessian_recalc` (e.g. to 50–200) helps — more frequent exact Hessian recalculations improve TS-mode tracking and convergence at the cost of additional Hessian evaluations.
```

## Notes

- Imaginary-frequency **detection** threshold defaults to 5.0 cm⁻¹ (configurable via `hessian_dimer.neg_freq_thresh_cm`). Frequencies with magnitudes below this threshold are not counted as imaginary.
- The selected `root` controls which vibrational mode is followed during optimization. It is set via YAML (`rsirfo.root` or `hessian_dimer.root`; default `0`); `tsopt` has no `--root` CLI flag, unlike [`irc`](irc.md).
- Use `--opt-mode` to choose the algorithm directly (`rsprfo` by default) rather than editing YAML mode mappings.
- Dimer orientation, rotation forces, flattening, and final exact PHVA validation use the same constrained projector as `freq`. The Dimer rebuilds this basis whenever its central image changes. It never subtracts translations of the active fragment unless they are actual rigid null directions compatible with every frozen anchor. Hessian RFO optimization itself operates on the active-DOF Cartesian Hessian without this projection. See [Frozen Atoms](freeze-atoms.md#rigid-modes-with-frozen-boundaries).
- See {ref}`CLI Conventions: Configuration precedence <configuration-precedence>` for the full resolution order.

## See Also

[Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [path-search](path-search.md) · [irc](irc.md) · [freq](freq.md) · [all](all.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
