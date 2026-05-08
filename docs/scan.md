# `scan`

## Overview

> **Summary:** Drive a reaction coordinate by scanning bond distances with harmonic restraints. Use `--scan-lists/-s` to define targets as either a YAML/JSON spec file path (recommended) or inline Python literals.

### At a glance
- **Use when:** A single structure needs specific distances driven to explore a plausible path (often before `path-search`/`path-opt`). Input is one structure + `-s/--scan-lists scan.yaml` (recommended), or one or more `--scan-lists/-s` inline literals (each literal = one stage). YAML/JSON file paths avoid shell-quoting pitfalls and version better; inline literals are fine for simple single-stage scans.
- **Method:** MLIP backend (UMA by default; selectable via `-b/--backend`) with harmonic restraints `E = Σ ½ k (|ri − rj| − target)²` and LBFGS (`--opt-mode grad`) or RFOptimizer (`--opt-mode hess`) per step.
- **Outputs:** Per-stage `result.xyz` (+ optional `.pdb`/`.gjf`), and concatenated scan trajectories. See "Notes" for the `--dump` caveat.
- **Defaults:** `--opt-mode grad` (LBFGS), `--no-preopt`, `--no-endopt`, `--max-step-size 0.20 Å`, `--bias-k 300 eV·Å⁻²`, `--thresh gau`, `--out-dir ./result_scan/`.
- **Next step:** Feed the staged endpoints (`stage_XX/result.pdb`) to `path-search`/`path-opt` for MEP refinement, or use `pdb2reaction all -s ...` to chain scan → MEP → TSOPT/IRC/freq/DFT in one command.

`pdb2reaction scan` performs a staged, bond-length–driven scan using an MLIP backend (UMA by default) and harmonic restraints. At each step, the temporary targets are updated, restraint wells are applied, and the structure is relaxed with LBFGS (`--opt-mode grad`) or RFOptimizer (`--opt-mode hess`).


For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB/GJF output conversion.

## Minimal example

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml -o ./result_scan
```

## Output checklist

- `result_scan/stage_01/result.pdb` (or `result.xyz`)
- `result_scan/stage_02/result.pdb` (or `result.xyz`)
- `result_scan/stage_*/scan_trj.xyz` (always written; `scan.pdb` companion when input is PDB and conversion is enabled)

## Common examples

1. Run from a YAML spec.

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml
```

2. Use literal input.

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 -s '[("TYR,285,CA","SAM,309,C10",1.35)]'
```

3. Dump trajectories for stage-by-stage inspection.

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml --dump -o ./result_scan_dump
```

> **Note:** Add `--print-parsed` when you want to verify parsed stage targets from `--scan-lists/-s`.

## Usage
```bash
pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan.yaml | '[(i,j,targetÅ),...]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### Examples
```bash
# Recommended: YAML/JSON spec file
cat > scan.yaml << 'YAML'
one_based: true
stages:
 - [["TYR,285,CA", "SAM,309,C10", 1.35]]
 - [["TYR,285,CA", "SAM,309,C10", 2.20], ["TYR,285,CB", "SAM,309,C11", 1.80]]
YAML
pdb2reaction scan -i input.pdb -q 0 -s scan.yaml

# Alternative: inline Python literal
pdb2reaction scan -i input.pdb -q 0 -s '[("TYR,285,CA","SAM,309,C10",1.35)]'

# Two stages, LBFGS relaxations, and trajectory dumping
pdb2reaction scan -i input.pdb -q 0 -s \
 '[("TYR,285,CA","SAM,309,C10",1.35)]' \
 '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
 --max-step-size 0.20 --dump -o ./result_scan/ --opt-mode grad \
 --preopt --endopt

# Supply multiple stage literals after a single -s/--scan-lists
pdb2reaction scan -i input.pdb -q 0 -s \
 '[("TYR,285,CA","SAM,309,C10",1.35)]' \
 '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]'
```

## Scan-list spec

For the YAML/JSON file format, inline Python literal syntax, atom selectors, and quoting rules,
see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`.

### Multiple stages

Pass multiple literals after a single `--scan-lists/-s` flag. Each literal becomes one stage:

```bash
# Stage 1: drive one bond to 1.35 Å
# Stage 2: drive two bonds simultaneously
-s \
 '[("TYR,285,CA","SAM,309,C10",1.35)]' \
 '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]'
```

Stages run sequentially; each starts from the previous stage's relaxed result.

### Bidirectional scan (4-tuple)

Instead of a 3-tuple `(i, j, target)`, you can pass a **4-tuple** `(i, j, start, end)` to scan in both directions from the current geometry. The CLI automatically expands each 4-tuple into two stages:

1. **Pass 1:** Drive `i`--`j` from the current distance toward `start`.
2. **Pass 2:** Restore the initial geometry and drive `i`--`j` toward `end`.

The concatenated trajectory is assembled as `start → initial → end`, giving a continuous path through the starting structure.

```bash
# Bidirectional scan: drive bond 12--45 from current geometry
# toward 1.35 Å (pass 1) and toward 2.50 Å (pass 2)
pdb2reaction scan -i input.pdb -q 0 -s '[(12, 45, 1.35, 2.50)]'
```

This is equivalent to two manual stages with a geometry reset between them, but avoids the need to script it yourself. Mixed 3-tuples and 4-tuples are accepted in the same literal.

```{note}
**Stage counter with 4-tuples.** A 4-tuple expands into **two** stages in the output tree: the `start` pass is written under `stage_NN/` and the `end` pass under `stage_NN+1/`. So if you pass a single 4-tuple as your first literal, you will see `stage_01/` and `stage_02/`, not one combined `stage_01/`. When mixing 3-tuples and 4-tuples, the counter advances by `+1` per 3-tuple and `+2` per 4-tuple.
```

## Workflow
1. Load the structure through `geom_loader`. Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>` for details).
2. Optionally run an unbiased preoptimization (`--preopt`) before any
    biasing so the starting point is relaxed.
3. Parse stage targets from `--scan-lists/-s` (YAML/JSON file or inline literal), then normalize the
    `(i, j)` indices (1-based by default). When the input is a PDB, each entry
    may be either an integer index or an atom selector string like `'TYR,285,CA'`;
    selector fields can be separated by spaces, commas, slashes, backticks, or
    backslashes and may be in any order (fallback assumes resname, resseq, atom).
    Compute the per-bond displacement
    `Δ = target − current` and split it into `N = ceil(max(|Δ|) / h)` steps using
    `h = --max-step-size`. Every bond receives its own `δ = Δ / N` increment.
4. March through all steps, updating the temporary targets, applying the
    harmonic wells `E = Σ ½ k (|ri − rj| − target)²`, and minimizing with the MLIP backend.
    Optimizer cycles are capped by `--relax-max-cycles` unless YAML specifies `opt.max_cycles`.
5. After the last step of each stage, optionally run an unbiased relaxation
    (`--endopt`) before reporting covalent bond changes and writing the
    `result.*` files.
6. Repeat for every stage. Concatenated scan trajectories (`scan_trj.xyz` and
    `scan.pdb`) are always written. Per-step optimizer trajectory files are not
    currently emitted from the CLI; set `opt.dump: true` in YAML if you need
    them.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge (CLI > template). When omitted, charge can be inferred from `--ligand-charge/-l`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge/-l` supplies it |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only; `workers_per_node` forwarded to the parallel predictor). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1. Inherits the `.gjf` template value when available; defaults to `1` when omitted. | `.gjf` template value or `1` |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (recommended) or inline Python literal with `(i,j,targetÅ)` triples or `(i,j,start,end)` 4-tuples for bidirectional scans. Each inline literal is one stage; supply multiple literals after a single flag. `i`/`j` can be integer indices or PDB atom selectors like `'TYR,285,CA'`. | Required |
| `--one-based/--zero-based` | Interpret atom indices as 1- or 0-based. These are mutually exclusive toggle aliases for the same flag (`--one-based` sets it to `True`, `--zero-based` sets it to `False`). | `True` |
| `--print-parsed/--no-print-parsed` | Print parsed stage tuples after `--scan-lists/-s` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per step (Å). Controls the number of integration steps. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². | `300` |
| `--relax-max-cycles INT` | Cap on optimizer cycles during preopt, each biased step, and end-of-stage cleanups. Used unless YAML sets `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `grad` → LBFGS, `hess` → RFOptimizer. See {ref}`opt-mode-semantics` for how the same token maps to different optimizers under `tsopt`. | `grad` |
| `--freeze-links/--no-freeze-links` | When the input is PDB, freeze the parents of link hydrogens. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--dump/--no-dump` | Reserved flag; the CLI does not currently forward this to the optimizer (`opt_cfg["dump"]` is always set to `False` for `scan`). To enable per-step optimizer trajectory output, set `opt.dump: true` in YAML. `scan_trj.xyz`/`scan.pdb` are always written regardless. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs (trajectory conversion only writes PDB). | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `-o, --out-dir TEXT` | Output directory root. | `./result_scan/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. **Scope-dependent default:** `False` standalone; flipped to `True` when invoked via `pdb2reaction all` (see [`all` → Scan Options](all.md#scan-options-single-input-runs)). | `False` |
| `--endopt/--no-endopt` | Run an unbiased optimization after each stage. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical keys to those documented in
  [YAML Reference](yaml-reference.md). `opt.dump` can be set in YAML for optimizer dumps;
  use `--dump` to control scan-stage trajectories.
- `--relax-max-cycles` applies only when explicitly provided **and** YAML does not set `opt.max_cycles` (default `10000`).

### Section `bias`
- `k` (`300`): Harmonic strength in eV·Å⁻².

(section-bond)=
### Section `bond`
MLIP-based bond-change detection shared with `path-search`:
- `device` (`"auto"`): MLIP device for bond analysis.
- `bond_factor` (`1.20`): Covalent-radius scaling for cutoff.
- `margin_fraction` (`0.05`): Fractional tolerance for comparisons.
- `delta_fraction` (`0.05`): Minimum relative change to flag formation/breaking.

## Outputs
```
out_dir/ (default:./result_scan/)
├─ preopt/ # Present when --preopt is True
│ ├─ result.xyz
│ ├─ result.pdb # PDB companion for PDB inputs when conversion is enabled
│ └─ result.gjf # When a Gaussian template exists and conversion is enabled
├─ stage_XX/ # One folder per stage
│ ├─ result.xyz
│ ├─ result.pdb # PDB mirror of the final structure (conversion enabled)
│ ├─ result.gjf # Gaussian mirror when templates exist and conversion is enabled
│ ├─ scan_trj.xyz # Always written (concatenated biased trajectory)
│ └─ scan.pdb # Always written for PDB inputs when conversion is enabled (no scan.gjf is produced)
├─ scan_trj.xyz # Combined trajectory across all stages
└─ scan.pdb # Combined PDB trajectory (when conversion is enabled)
```
- Console summaries of the resolved `geom`, `calc`, `opt`, `bias`, `bond`, and optimizer blocks plus per-stage bond-change reports.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Provide multiple literals after a single `--scan-lists/-s` flag.
  Tuples must have positive targets. Atom indices are normalized to 0-based internally for computation. For
  PDB inputs, `i`/`j` can be selector strings with flexible delimiters
  (space/comma/slash/backtick/backslash) and unordered tokens.
- When `--freeze-links` is active, link-hydrogen parent atoms are automatically frozen (see {ref}`Link hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).
- Stage results (`result.xyz` plus optional PDB/GJF companions) are always
  written. Concatenated scan trajectories (`scan_trj.xyz` and `scan.pdb` for
  PDB inputs with conversion enabled) are also always written. Per-step
  optimizer trajectory files are not currently emitted from the CLI — the
  `--dump` flag is reserved; use `opt.dump: true` in YAML to enable the
  optimizer-level dump.


```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 task_name: omol # UMA task name
 device: auto # MLIP device selection
 max_neigh: null # maximum neighbors for graph construction
 radius: null # cutoff radius for neighbor search
 r_edges: false # store radial edges
 out_hess_torch: true # request torch-form Hessian
 freeze_atoms: null # calculator-level frozen atoms
 hessian_calc_mode: FiniteDifference # Hessian mode selection
 return_partial_hessian: true  # partial Hessian over active DOFs
opt:
 thresh: gau # convergence preset (Gaussian/Baker-style)
 max_cycles: 10000 # optimizer cycle cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum norm for step acceptance
 assert_min_step: true # stop if steps fall below threshold
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # geom RMS threshold when converging to ref
 overachieve_factor: 0.0 # factor to tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 energy_plateau: true # enable plateau-based early stop
 energy_plateau_thresh: 1.0e-04 # plateau detection threshold
 energy_plateau_window: 50 # plateau detection window (cycles)
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir: ./result_scan/ # output directory
lbfgs:
 thresh: gau # LBFGS convergence preset
 max_cycles: 10000 # iteration limit
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 energy_plateau: true # enable plateau-based early stop
 energy_plateau_thresh: 1.0e-04 # plateau detection threshold
 energy_plateau_window: 50 # plateau detection window (cycles)
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir: ./result_scan/ # output directory
 keep_last: 7 # history size for LBFGS buffers
 beta: 1.0 # initial damping beta
 gamma_mult: false # multiplicative gamma update toggle
 max_step: 0.3 # maximum step length
 control_step: true # control step length adaptively
 double_damp: true # double damping safeguard
 mu_reg: null # regularization strength
 max_mu_reg_adaptions: 10 # cap on mu adaptations
rfo:
 thresh: gau # RFOptimizer convergence preset
 max_cycles: 10000 # iteration cap
 print_every: 100 # logging stride
 min_step_norm: 1.0e-08 # minimum accepted step norm
 assert_min_step: true # assert when steps stagnate
 rms_force: null # explicit RMS force target
 rms_force_only: false # rely only on RMS force convergence
 max_force_only: false # rely only on max force convergence
 force_only: false # skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when targeting geometry
 overachieve_factor: 0.0 # tighten thresholds
 check_eigval_structure: false # validate Hessian eigenstructure
 energy_plateau: true # enable plateau-based early stop
 energy_plateau_thresh: 1.0e-04 # plateau detection threshold
 energy_plateau_window: 50 # plateau detection window (cycles)
 line_search: true # enable line search
 dump: false # dump trajectory/restart data
 dump_restart: false # dump restart checkpoints
 prefix: "" # filename prefix
 out_dir: ./result_scan/ # output directory
 trust_radius: 0.10 # trust-region radius
 trust_update: true # enable trust-region updates
 trust_min: 0.0001 # minimum trust radius
 trust_max: 0.10 # maximum trust radius
 max_energy_incr: null # allowed energy increase per step
 hessian_update: bfgs # Hessian update scheme
 hessian_init: calc # Hessian initialization source
 hessian_recalc: 500 # rebuild Hessian every N steps
 hessian_recalc_adapt: null # adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # eigenvalue threshold for stability
 alpha0: 1.0 # initial micro step
 max_micro_cycles: 50 # micro-iteration limit
 rfo_overlaps: false # enable RFO overlaps
 gediis: false # enable GEDIIS
 gdiis: true # enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # test descent direction before DIIS
 adapt_step_func: true # adaptive step scaling toggle
bias:
 k: 300 # harmonic bias strength (eV·Å⁻²)
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # covalent-radius scaling
 margin_fraction: 0.05 # tolerance margin for comparisons
 delta_fraction: 0.05 # minimum relative change to flag bonds
```

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing

- [all](all.md) — End-to-end workflow with `--scan-lists/-s` for single-structure inputs
- [scan2d](scan2d.md) — Two-distance grid scan (d₁, d₂) with the same MLIP backend and YAML controls
- [scan3d](scan3d.md) — Three-distance grid scan (d₁, d₂, d₃) with isosurface output
- [path-search](path-search.md) — MEP search using scan endpoints as intermediates
- [extract](extract.md) — Generate active site model (binding pocket) PDBs before scanning
- [YAML Reference](yaml-reference.md) — Full `bias` and `bond` configuration options
- [Glossary](glossary.md) — Definitions of MEP, Segment
