# `scan`

Drive a reaction coordinate by scanning bond distances with harmonic restraints. Use `pdb2reaction scan` to drive specific distances in a single structure and explore a plausible path (often before `path-search`/`path-opt`). It performs a staged, bond-length–driven scan using an MLIP backend (UMA by default) and harmonic restraints. At each step, the temporary targets are updated, restraint wells are applied, and the structure is relaxed with L-BFGS (`--opt-mode grad`) or RFOptimizer (`--opt-mode hess`). For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB/GJF output conversion.

## Examples

```bash
# Minimal: run from a YAML spec
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml -o ./result_scan
```

```bash
# Inline Python literal
pdb2reaction scan -i input.pdb -q 0 -m 1 -s '[("TYR,285,CA","SAM,309,C10",1.35)]'
```

```bash
# Dump trajectories for stage-by-stage inspection
pdb2reaction scan -i input.pdb -q 0 -m 1 -s scan.yaml --dump -o ./result_scan_dump
```

Command form:

```bash
pdb2reaction scan -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan.yaml | '[(i,j,targetÅ),...]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

> **Note:** Add `--print-parsed` when you want to verify parsed stage targets from `--scan-lists/-s`.

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
    `scan.pdb`) are always written. Pass `--dump` to additionally emit per-step
    optimizer trajectory files (`opt.dump` from YAML is run-scoped and ignored).

## Outputs

```text
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

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge (CLI > template). When omitted, charge can be inferred from `--ligand-charge/-l`; explicit `-q` overrides any derived value. | Required unless a `.gjf` template or `--ligand-charge/-l` supplies it |
| `-l, --ligand-charge TEXT` | Either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB residue charges. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only; `workers_per_node` forwarded to the parallel predictor). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1. Inherits the `.gjf` template value when available; defaults to `1` when omitted. | `.gjf` template value or `1` |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (recommended) or inline Python literal with `(i,j,targetÅ)` triples or `(i,j,start,end)` 4-tuples for bidirectional scans. Each inline literal is one stage; supply multiple literals after a single flag. `i`/`j` can be integer indices or PDB atom selectors like `'TYR,285,CA'`. | Required |
| `--one-based/--zero-based` | Interpret atom indices as 1- or 0-based. These are mutually exclusive toggle aliases for the same flag (`--one-based` sets it to `True`, `--zero-based` sets it to `False`). | `True` |
| `--print-parsed/--no-print-parsed` | Print parsed stage tuples after `--scan-lists/-s` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per step (Å). Controls the number of integration steps. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². | `300` |
| `--relax-max-cycles INT` | Cap on optimizer cycles during preopt, each biased step, and end-of-stage cleanups. Used unless YAML sets `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `grad` → L-BFGS, `hess` → RFOptimizer. See {ref}`opt-mode-semantics` for how the same token maps to different optimizers under `tsopt`. | `grad` |
| `--freeze-links/--no-freeze-links` | When the input is PDB, freeze the parents of cap hydrogens. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--dump/--no-dump` | Forward to the per-step optimizer (`opt_cfg["dump"]`), emitting per-step optimizer trajectory files. `scan_trj.xyz`/`scan.pdb` are always written regardless. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs (trajectory conversion only writes PDB). | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `-o, --out-dir TEXT` | Output directory root. | `./result_scan/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. **Scope-dependent default:** `False` standalone; flipped to `True` when invoked via `pdb2reaction all` (see [`all` → TSOPT / freq / DFT / scan overrides](all.md#tsopt--freq--dft--scan-overrides)). | `False` |
| `--endopt/--no-endopt` | Run an unbiased optimization after each stage. | `False` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical keys to those documented in
  [YAML Reference](yaml-reference.md). Per-step optimizer trajectories are controlled
  by `--dump` (CLI) only — `opt.dump` and `opt.out_dir` from YAML are run-scoped and
  overwritten (not YAML-tunable); the scan-stage trajectories `scan_trj.xyz`/`scan.pdb`
  are always written regardless.
- `--relax-max-cycles` applies only when explicitly provided **and** YAML does not set `opt.max_cycles` (default `10000`).

### Section `bias`
- `k` (`300`): Harmonic strength in eV·Å⁻².

(section-bond)=
### Section `bond`
MLIP-based bond-change detection shared with `path-search`. Full keys and
defaults (`device`, `bond_factor`, `margin_fraction`, `delta_fraction`): see
[YAML Reference](yaml-reference.md#bond).

## YAML configuration

```yaml
geom:
 coord_type: cart        # cartesian vs dlc internals
calc:
 model: uma-s-1p1        # uma-s-1p1 | uma-m-1p1
 task_name: omol         # UMA task name
opt:
 thresh: gau             # convergence preset
 max_cycles: 10000       # optimizer cycle cap
 # out_dir is run-scoped: set via -o/--out-dir, not YAML (a YAML value here is ignored)
lbfgs:
 max_step: 0.3           # maximum step length (grad mode)
rfo:
 trust_radius: 0.10      # trust-region radius (hess mode)
bias:
 k: 300                  # harmonic bias strength (eV·Å⁻²)
bond:
 bond_factor: 1.2        # covalent-radius scaling
 margin_fraction: 0.05   # tolerance margin
 delta_fraction: 0.05    # minimum relative change to flag bonds
```

More YAML options for `opt`/`lbfgs`/`rfo`/`bias`/`bond` and their defaults are in [YAML Reference](yaml-reference.md).

## Scan-list spec

For the YAML/JSON file format, inline Python literal syntax, atom selectors, and quoting rules,
see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`.

### Staged vs concerted scans

The number of `(i, j, target)` tuples inside one literal and the number of literals together decide whether the coordinates are driven *together* (concerted) or *in sequence* (staged):

| Mode | Syntax | Use when |
| --- | --- | --- |
| Concerted | one `-s` with several coordinate tuples | The coordinates move together in a single step; you do not need to break the mechanism into stages |
| Staged | `-s` repeated (one literal per sequential stage) | The mechanism is known up front and you want clean per-step control and per-stage output |

When the mechanism **is** known, the staged form is generally preferred — it gives per-step barriers and per-stage geometries. When the mechanism is unknown or multi-step, let [`path-search`](path-search.md) auto-segment the path instead of guessing the stages yourself. (A 4-tuple `(i, j, low, high)` expands into a bidirectional 2-stage scan; see [Bidirectional scan](#bidirectional-scan-4-tuple).)

```bash
# Concerted: two coordinates move together in one stage
pdb2reaction scan -i reactant.pdb \
    -s '[("Ca RES 10","Cb RES 11",1.6),("H RES 11","O GLU 20",1.0)]' -o result_concerted
```

Pass multiple literals after a single `--scan-lists/-s` flag for a staged scan. Each literal becomes one stage:

```bash
# Stage 1: drive one bond to 1.35 Å
# Stage 2: drive two bonds simultaneously
-s \
 '[("TYR,285,CA","SAM,309,C10",1.35)]' \
 '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]'
```

Stages run sequentially; each starts from the previous stage's relaxed result.

(scan-direction-and-barrier-sign)=
### Scan direction and barrier sign

If a `scan` (or path) **starts from the product** side, the raw barrier it reports is the **reverse** barrier, `E(TS) − E(product)`. To quote the forward barrier, compute it from the reactant:

| You ran | Forward barrier |
| --- | --- |
| A product-start scan | `E(TS) − E(reactant)` — **not** the raw product-start number |

This is something to interpret when reading results, not a CLI flag. Always confirm which endpoint the scan started from before quoting a barrier, especially when the workflow was seeded from a crystallographic product complex.

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

This is equivalent to two manual stages with a geometry reset between them. Mixed 3-tuples and 4-tuples are accepted in the same literal.

```{note}
**Stage counter with 4-tuples.** A 4-tuple expands into **two** stages in the output tree: the `start` pass is written under `stage_NN/` and the `end` pass under `stage_NN+1/`. So if you pass a single 4-tuple as your first literal, you will see `stage_01/` and `stage_02/`, not one combined `stage_01/`. When mixing 3-tuples and 4-tuples, the counter advances by `+1` per 3-tuple and `+2` per 4-tuple.
```

## Notes

- The scan input is one structure plus `-s/--scan-lists scan.yaml` (recommended) or one or more `--scan-lists/-s` inline literals (each literal = one stage). YAML/JSON file paths avoid shell-quoting pitfalls and version better; inline literals are fine for simple single-stage scans.
- Provide multiple literals after a single `--scan-lists/-s` flag. Tuples must have positive targets. Atom indices are normalized to 0-based internally for computation. For PDB inputs, `i`/`j` can be integer indices or selector strings (see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`).
- When `--freeze-links` is active, cap-hydrogen parent atoms are automatically frozen (see {ref}`Cap hydrogen and frozen atoms <link-hydrogen-and-frozen-atoms>`).

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed fixes for common failure modes
- [all](all.md) — End-to-end workflow with `--scan-lists/-s` for single-structure inputs
- [scan2d](scan2d.md) — Two-distance grid scan (d₁, d₂) with the same MLIP backend and YAML controls
- [scan3d](scan3d.md) — Three-distance grid scan (d₁, d₂, d₃) with isosurface output
- [path-search](path-search.md) — MEP search using scan endpoints as intermediates
- [extract](extract.md) — Generate active site model (binding pocket) PDBs before scanning
- [YAML Reference](yaml-reference.md) — Full `bias` and `bond` configuration options
- [Glossary](glossary.md) — Definitions of MEP, Segment
