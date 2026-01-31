# `all`

## Overview

`pdb2reaction all` runs the entire workflow end-to-end:

pocket extraction → (optional) staged UMA scan → recursive MEP search (`path_search`, GSM/DMF) → merge back into the full system → (optional) TS optimization + IRC (`tsopt`) → (optional) vibrational analysis / thermochemistry (`freq`) → (optional) single-point DFT (`dft`).

It supports three common modes:

- **Multi-structure workflow** — Provide ≥2 structures (PDB/GJF/XYZ) in reaction order plus a substrate definition. `all` extracts pockets, runs GSM/DMF MEP search, merges the optimized path back into the full-system template(s), and optionally runs TSOPT/freq/DFT per reactive segment.
- **Single-structure + staged scan** — Provide one structure plus one or more `--scan-lists`. The scan generates an ordered set of intermediates that become MEP endpoints.
  - One `--scan-lists` literal runs a single scan stage.
  - Multiple stages are passed as multiple values after a single `--scan-lists` flag (the flag itself cannot be repeated).
- **TSOPT-only pocket TS optimization** — Provide a single input structure, omit `--scan-lists`, and set `--tsopt True`. `all` extracts the pocket (if `-c/--center` is given) and runs TS optimization + IRC, with optional freq/DFT, on that single system.

PDB/GJF companion files are generated when templates are available, controlled by `--convert-files {True\|False}` (enabled by default).


## Usage
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [options]
```

### Examples
```bash
# Multi-structure ensemble with explicit ligand charges and post-processing
pdb2reaction all -i reactant.pdb product.pdb -c 'GPP,MMT' \
    --ligand-charge 'GPP:-3,MMT:-1' --mult 1 --freeze-links True \
    --max-nodes 10 --max-cycles 100 --climb True --opt-mode light \
    --out-dir result_all_${date} --tsopt True --thermo True --dft True

# Single-structure staged scan followed by GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c '308,309' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
    --opt-mode heavy --tsopt True --thermo True --dft True

# TSOPT-only workflow (no path search)
pdb2reaction all -i reactant.pdb -c 'GPP,MMT' \
    --ligand-charge 'GPP:-3,MMT:-1' --tsopt True --thermo True --dft True
```

## Workflow
1. **Active-site pocket extraction** (if `-c/--center` is provided)
   - Substrates may be specified via PDB paths, residue IDs (`123,124` or `A:123,B:456`), or residue names (`GPP,MMT`).
   - Optional toggles forward to the extractor: `--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`, and `--verbose`.
   - Per-input pocket PDBs are saved under `<out-dir>/pockets/`. When multiple structures are supplied, their pockets are unioned per residue selection.
   - The **first pocket’s total charge** is propagated to scan/MEP/TSOPT.

2. **Optional staged scan (single-input only)**
   - Each `--scan-lists` argument is a Python-like list of `(i,j,target_Å)` tuples describing a UMA scan stage. Atom indices refer to the original input ordering (1-based) and are remapped to the pocket ordering. For PDB inputs, `i`/`j` can be integer indices or selector strings like `'TYR,285,CA'`; selectors accept spaces/commas/slashes/backticks/backslashes (` ` `,` `/` `` ` `` `\`) as delimiters and allow unordered tokens (fallback assumes resname, resseq, atom).
   - A single literal runs a one-stage scan; multiple literals run **sequentially** so stage 2 begins from stage 1's result, and so on. Supply multiple literals after a single flag (repeated flags are not accepted).
   - Scan inherits charge/spin, `--freeze-links`, the UMA optimizer preset (`--opt-mode`), `--args-yaml`, and `--preopt`. The `--dump` flag is forwarded to scan only when explicitly set on this command; otherwise scan uses its own default (`False`). Overrides such as `--scan-out-dir`, `--scan-one-based`, `--scan-max-step-size`, `--scan-bias-k`, `--scan-relax-max-cycles`, `--scan-preopt`, and `--scan-endopt` apply per run.
   - Stage endpoints (`stage_XX/result.pdb`) become the ordered intermediates that feed the subsequent MEP step.

3. **MEP search on pockets (recursive GSM/DMF)**
   - Executes `path_search` by default using the extracted pockets (or the original entire structures if extraction is skipped). Relevant options: `--mult`, `--freeze-links`, `--max-nodes`, `--max-cycles`, `--climb`, `--opt-mode`, `--dump`, `--preopt`, `--args-yaml`, and `--out-dir`.
   - Use `--refine-path False` to switch to a single-pass `path-opt` GSM/DMF chain without the recursive refiner.
   - For multi-input PDB runs, the full-system templates are automatically passed to `path_search` for reference merging. Single-structure scan runs reuse the original full PDB template for every stage.

4. **Merge pockets back to the full systems**
   - When reference PDB templates exist, merged `mep_w_ref*.pdb` and per-segment `mep_w_ref_seg_XX.pdb` files are emitted under `<out-dir>/path_search/`.

5. **Optional per-segment post-processing**
   - `--tsopt True`: run TS optimization on each HEI pocket, follow with EulerPC IRC, and emit segment energy diagrams.
   - `--thermo True`: call `freq` on (R, TS, P) to obtain vibrational/thermochemistry data and a UMA Gibbs diagram.
   - `--dft True`: launch single-point DFT on (R, TS, P) and build a DFT diagram. When combined with `--thermo True`, a DFT//UMA Gibbs diagram (DFT energies + UMA thermal correction) is also produced.
   - Shared overrides include `--opt-mode`, `--opt-mode-post` (overrides TSOPT/post-IRC optimization mode), `--flatten-imag-mode`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, and `--dft-engine` (GPU-first by default).
   - When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.

6. **TSOPT-only mode** (single input, `--tsopt True`, no `--scan-lists`)
   - Skips the MEP/merge stages. Runs `tsopt` on the pocket (or full input if extraction is skipped), performs EulerPC IRC, identifies the higher-energy endpoint as reactant (R), and generates the same set of energy diagrams plus optional freq/DFT outputs.

### Charge and spin precedence

**Charge resolution (highest to lowest priority):**

| Priority | Source | When Used |
|----------|--------|-----------|
| 1 | `-q/--charge` | Explicit CLI override |
| 2 | Pocket extraction | When `-c` is provided (sums amino acids, ions, `--ligand-charge`) |
| 3 | `--ligand-charge` (numeric) | Fallback when extraction fails or is skipped |
| 4 | `.gjf` template | Embedded charge/spin metadata |
| 5 | Default | 0 |

**Spin resolution:** `--mult` (CLI) → `.gjf` template → default (1)

> **Tip:** Always provide `--ligand-charge` for non-standard substrates to ensure correct charge propagation.

### Input expectations
- Extraction enabled (`-c/--center`): inputs must be **PDB** files so residues can be located.
- Extraction skipped: inputs may be **PDB/XYZ/GJF**.
- Multi-structure runs require ≥2 structures.

## CLI Options

> **Note:** Default values shown are used when the option is not specified.

### Input/Output Options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full structures in reaction order (single input allowed only with `--scan-lists` or `--tsopt True`). | Required |
| `--out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--convert-files {True\|False}` | Global toggle for XYZ/TRJ → PDB/GJF companions when templates are available. | `True` |
| `--dump BOOLEAN` | Dump MEP (GSM/DMF) trajectories. | `False` |
| `--args-yaml FILE` | YAML forwarded unchanged to all subcommands. | _None_ |

### Charge/Spin Options

| Option | Description | Default |
| --- | --- | --- |
| `--ligand-charge TEXT` | Total charge or residue-specific mapping for unknown residues (recommended). | _None_ |
| `-q, --charge INT` | Force the total system charge (overrides `--ligand-charge`). | _None_ |
| `-m, --mult INT` | Spin multiplicity forwarded to all downstream steps. | `1` |

### Extraction Options

| Option | Description | Default |
| --- | --- | --- |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). | Required for extraction |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). | `0.0` |
| `--include-H2O BOOLEAN` | Include waters (HOH/WAT/TIP3/SOL). | `True` |
| `--exclude-backbone BOOLEAN` | Remove backbone atoms on non-substrate amino acids. | `True` |
| `--add-linkH BOOLEAN` | Add link hydrogens for severed bonds. | `True` |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--freeze-links BOOLEAN` | Freeze link parents in pocket PDBs. | `True` |
| `--verbose BOOLEAN` | Enable INFO-level extractor logging. | `True` |

### MEP Search Options

| Option | Description | Default |
| --- | --- | --- |
| `--max-nodes INT` | MEP internal nodes per segment. | `10` |
| `--max-cycles INT` | MEP maximum optimization cycles. | `300` |
| `--climb BOOLEAN` | Enable TS climbing for the first segment. | `True` |
| `--opt-mode [light\|heavy]` | Optimizer preset (light → LBFGS/Dimer, heavy → RFO/RSIRFO). | `light` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`). | `gau` |
| `--preopt BOOLEAN` | Pre-optimize pocket endpoints before MEP search. | `True` |

### UMA Calculator Options

| Option | Description | Default |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians). | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | Shared UMA Hessian engine. | `FiniteDifference` |

### Post-Processing Options

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt BOOLEAN` | Run TS optimization + IRC per reactive segment. | `False` |
| `--thermo BOOLEAN` | Run vibrational analysis (`freq`) on R/TS/P. | `False` |
| `--dft BOOLEAN` | Run single-point DFT on R/TS/P. | `False` |
| `--opt-mode-post [light\|heavy]` | Optimizer preset for TSOPT and post-IRC optimization. | _None_ |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations. | `baker` |
| `--flatten-imag-mode {True\|False}` | Enable extra-imaginary-mode flattening in `tsopt`. | `False` |

### TSOPT Overrides

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles`. | `10000` |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory. | _None_ |

### Freq Overrides

| Option | Description | Default |
| --- | --- | --- |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Maximum modes to write. | `10` |
| `--freq-amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--freq-n-frames INT` | Frames per mode animation. | `20` |
| `--freq-sort [value\|abs]` | Mode sorting behavior. | `value` |
| `--freq-temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--freq-pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |

### DFT Overrides

| Option | Description | Default |
| --- | --- | --- |
| `--dft-engine [gpu\|cpu\|auto]` | Preferred backend (`auto` tries GPU then CPU). | `gpu` |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Functional/basis pair. | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | Maximum SCF iterations. | `100` |
| `--dft-conv-tol FLOAT` | SCF convergence tolerance. | `1e-9` |
| `--dft-grid-level INT` | PySCF grid level. | `3` |

### Scan Options (Single-Input Runs)

| Option | Description | Default |
| --- | --- | --- |
| `--scan-lists TEXT...` | Staged scans: `(i,j,target_Å)` tuples. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | _None_ |
| `--scan-one-based BOOLEAN` | Force scan indexing to 1-based or 0-based. | `True` |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | `0.20` |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV/Å²). | `100` |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | `10000` |
| `--scan-preopt BOOLEAN` | Override the scan preoptimization toggle. | `True` |
| `--scan-endopt BOOLEAN` | Override the scan end-of-stage optimization toggle. | `True` |

## Outputs
```text
out_dir/ (default: ./result_all/)
├─ summary.log               # formatted summary for quick inspection
├─ summary.yaml              # YAML version summary
├─ pockets/                  # Per-input pocket PDBs when extraction runs
├─ scan/                     # Staged pocket scan results (present when --scan-lists is provided)
├─ path_search/              # MEP results (GSM/DMF): trajectories, merged PDBs, diagrams, summary.yaml, per-segment folders
├─ path_search/post_seg_XX/  # Post-processing outputs (TS optimization, IRC, freq, DFT, diagrams)
└─ tsopt_single/             # TSOPT-only outputs with IRC endpoints and optional freq/DFT directories
```
- Console logs summarizing pocket charge resolution, YAML contents, scan stages, MEP progress (GSM/DMF), and per-stage timing.

### Reading `summary.log`
The log is organized into numbered sections:
- **[1] Global MEP overview** – image/segment counts, MEP trajectory plot paths, and the aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (UMA path)** – per-segment barriers (`ΔE‡`), reaction energies (`ΔE`), and bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** – per-segment TS imaginary frequency checks, IRC outputs, and UMA/thermo/DFT energy tables.
- **[4] Energy diagrams (overview)** – diagram tables for MEP/UMA/Gibbs/DFT series plus an optional cross-method summary table.
- **[5] Output directory structure** – a compact tree of generated files with inline annotations.

### Reading `summary.yaml`
The YAML is a compact, machine-readable summary. Common top-level keys include:
- `out_dir`, `n_images`, `n_segments` – run metadata and total counts.
- `segments` – list of per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, and `bond_changes`.
- `energy_diagrams` (optional) – diagram payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, and `image` paths.

`summary.yaml` intentionally omits the formatted tables and filesystem tree that appear in `summary.log`.

## Notes
- Always provide `--ligand-charge` (numeric or per-residue mapping) when formal charges cannot be inferred so the correct total charge propagates to scan/MEP/TSOPT/DFT.
- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-full-pdb` option of `path_search` is intentionally hidden in this wrapper.
- Convergence presets: `--thresh` defaults to `gau`; `--thresh-post` defaults to `baker`.
- Extraction radii: passing `0` to `--radius` or `--radius-het2het` is internally clamped to `0.001 Å` by the extractor.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to the MEP/tsopt/freq/DFT stages; single-structure runs still require either `--scan-lists` or `--tsopt True`.
- `--args-yaml` lets you coordinate all calculators from a single configuration file. YAML values override CLI flags.

## YAML configuration (`--args-yaml`)

The same YAML file is forwarded unchanged to **every** invoked subcommand. Each tool reads the sections described in its own documentation:

| Subcommand | YAML Sections |
|------------|---------------|
| [`path_search`](path_search.md) | `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

> **Note:** YAML contents take precedence over CLI values when both are provided.

**Minimal example:**
```yaml
calc:
  model: uma-s-1p1
  hessian_calc_mode: Analytical  # recommended when VRAM permits
gs:
  max_nodes: 12
  climb: true
dft:
  grid_level: 6
```

For a complete reference of all YAML options, see **[YAML Configuration Reference](yaml-reference.md)**.

---

## See Also

- [Getting Started](getting-started.md) — Installation and first run tutorial
- [Concepts & Workflow](concepts.md) — Mental model of pockets, segments, and stages
- [extract](extract.md) — Standalone pocket extraction (called internally by `all`)
- [path-search](path_search.md) — Standalone MEP search (called internally by `all`)
- [tsopt](tsopt.md) — Standalone TS optimization
- [freq](freq.md) — Standalone vibrational analysis
- [dft](dft.md) — Standalone DFT calculations
- [Troubleshooting](troubleshooting.md) — Common errors and fixes
- [YAML Reference](yaml-reference.md) — Complete YAML configuration options
- [Glossary](glossary.md) — Definitions of MEP, TS, IRC, GSM, DMF