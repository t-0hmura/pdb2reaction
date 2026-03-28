# `all`

## Overview

`pdb2reaction all` runs the entire workflow end-to-end:

pocket extraction → (optional) staged MLIP scan → recursive MEP search (`path-search`, GSM/DMF) → merge back into the full system → (optional) TS optimization + IRC (`tsopt`) → (optional) vibrational analysis / thermochemistry (`freq`) → (optional) single-point DFT (`dft`). The default MLIP backend is UMA; select an alternative with `-b/--backend`.

```{important}
`--tsopt` produces **TS candidates**. `all` automatically runs IRC for endpoint validation (`tsopt` itself includes an imaginary-frequency check), but always inspect the results (imaginary-frequency count + endpoint connectivity) before mechanistic interpretation.
```

It supports three common modes:

- **Multi-structure workflow** — Provide ≥2 structures (PDB/GJF/XYZ) in reaction order plus a substrate definition. `all` extracts pockets, runs GSM/DMF MEP search, merges the optimized path back into the full-system template(s), and optionally runs TSOPT+IRC/freq/DFT per reactive segment.
- **Single-structure + staged scan** — Provide one structure plus one or more `--scan-lists`. The scan generates an ordered set of intermediates that become MEP endpoints.
 - One `--scan-lists` literal runs a single scan stage.
 - Multiple stages are passed by repeating `--scan-lists`.
- **TSOPT-only pocket TS optimization** — Provide a single input structure, omit `--scan-lists`, and set `--tsopt`. `all` extracts the pocket (if `-c/--center` is given) and runs TS optimization + IRC, with optional freq/DFT, on that single system.

## Minimal example

```bash
pdb2reaction all -i R.pdb -i P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" --out-dir ./result_all
```

## Output checklist

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb` (or `result_all/path_search/seg_*/`)

## Common examples

1. Run full post-processing in one command.

```bash
pdb2reaction all -i R.pdb -i P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_all
```

2. Single-structure staged scan route.

```bash
pdb2reaction all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" --scan-lists "[(10,55,2.20)]" \
 --multiplicity 1 --out-dir ./result_scan_all
```

PDB/GJF companion files are generated when templates are available, controlled by `--convert-files` (enabled by default).


## Usage
```bash
pdb2reaction all -i INPUT1 -i [INPUT2 ...] -c SUBSTRATE [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
```

For help output, `pdb2reaction all --help` shows core options and `pdb2reaction all --help-advanced` shows the full option list.

### Examples
```bash
# Multi-structure ensemble with explicit ligand charges and post-processing
pdb2reaction all -i reactant.pdb -i product.pdb -c 'GPP,SAM' \
 -l 'GPP:-3,SAM:1' --multiplicity 1 --freeze-links \
 --max-nodes 10 --max-cycles 100 --climb --opt-mode grad \
 --out-dir ./result_all --tsopt --thermo --dft

# Single-structure staged scan followed by GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c '308,309' \
 --scan-lists '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
 --opt-mode hess --tsopt --thermo --dft

# TSOPT-only workflow (no path search)
pdb2reaction all -i reactant.pdb -c 'GPP,SAM' \
 -l 'GPP:-3,SAM:1' --tsopt --thermo --dft
```

## Workflow

0. **Preflight checks** (automatic)
 - `all` automatically runs `add-elem-info` (fills missing element symbols in PDB columns 77–78) and `fix-altloc` (resolves alternate conformations) on every PDB input before any other processing. When using individual subcommands (e.g., `extract`, `opt`), you must run these manually if needed.

1. **Active-site pocket extraction** (if `-c/--center` is provided)
 - Substrates may be specified via PDB paths, residue IDs (`123,124` or `A:123,B:456`), or residue names (`GPP,SAM`).
 - Optional toggles forward to the extractor: `--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`, and `--verbose`.
 - Per-input pocket PDBs are saved under `<out-dir>/pockets/`. When multiple structures are supplied, their pockets are unioned per residue selection.
 - The **first pocket’s net charge** is propagated to scan/MEP/TSOPT.

2. **Optional staged scan (single-input only)**
 - Each `--scan-lists` argument is a Python-like list of `(i,j,target_Å)` tuples describing an MLIP scan stage. Atom indices refer to the original input ordering (1-based) and are remapped to the pocket ordering. For PDB inputs, `i`/`j` can be integer indices or selector strings like `'TYR,285,CA'`; selectors accept spaces/commas/slashes/backticks/backslashes (` ` `,` `/` `` ` `` `\`) as delimiters and allow unordered tokens (fallback assumes resname, resseq, atom).
 - A single literal runs a one-stage scan; multiple literals run **sequentially** so stage 2 begins from stage 1's result, and so on. Supply multiple literals by repeating `--scan-lists`.
 - Stage endpoints (`stage_XX/result.pdb`) become the ordered intermediates that feed the subsequent MEP step.

3. **MEP search on pockets (recursive GSM/DMF)**
 - Use `--refine-path False` to switch to a single-pass `path-opt` GSM/DMF chain without the recursive refiner.
 - For multi-input PDB runs, the full-system templates are automatically passed to `path-search` for reference merging. Single-structure scan runs reuse the original full PDB template for every stage.

4. **Merge pockets back to the full systems**
 - When reference PDB templates exist, merged `mep_w_ref*.pdb` and per-segment `mep_w_ref_seg_XX.pdb` files are emitted under `<out-dir>/path_search/`.

5. **Optional per-segment post-processing** (only for reactive segments — segments with bond changes; bridge segments are skipped)
 - `--tsopt`: run TS optimization on each HEI pocket, follow with EulerPC-based IRC, then re-optimize IRC endpoints with `--thresh-post` (default `baker`). The endpoint optimization working directory is automatically deleted after completion.
 - `--thermo`: call `freq` on (R, TS, P) to obtain vibrational/thermochemistry data and an MLIP Gibbs diagram.
 - `--dft`: launch single-point DFT on (R, TS, P) and build a DFT diagram. When combined with `--thermo`, a DFT//MLIP Gibbs diagram (DFT energies + MLIP thermal correction) is also produced.
  - Shared overrides include `--opt-mode`, `--opt-mode-post` (overrides TSOPT/post-IRC optimization mode), `--flatten/--no-flatten`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, and `--dft-engine` (GPU-first by default).
 - For Hessian evaluation modes, see [MLIP Calculator](uma_pysis.md#hessian-evaluation).

6. **TSOPT-only mode** (single input, `--tsopt`, no `--scan-lists`)
 - Skips the MEP/merge stages. Runs `tsopt` on the pocket (or full input if extraction is skipped), performs EulerPC IRC, identifies the higher-energy endpoint as reactant (R), and generates the same set of energy diagrams plus optional freq/DFT outputs.

### Charge and spin precedence

Charge is resolved via the standard priority chain (see [CLI Conventions: Charge specification](cli_conventions.md#charge-specification) for details).

**Spin resolution:** `--multiplicity` (CLI) → `.gjf` template → default (1)

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
| `-i, --input PATH...` | Two or more full structures in reaction order (single input allowed only with `--scan-lists` or `--tsopt`). | Required |
| `--ref-pdb FILE` | Reference PDB for topology when `-i` provides XYZ inputs. | _None_ |
| `-o, --out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--convert-files/--no-convert-files` | Global toggle for XYZ/TRJ → PDB/GJF companions when templates are available. | `True` |
| `--dump/--no-dump` | Dump MEP (GSM/DMF) trajectories. Always forwarded to `path-search`/`path-opt`; forwarded to `scan`/`tsopt` only when explicitly set here. `freq` defaults to dump=True unless you pass `--no-dump`. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration before execution. | `False` |
| `--dry-run/--no-dry-run` | Validate and print plan without running stages. | `False` |
| `--resume/--no-resume` | Resume a previous run from `--out-dir`. Completed stages whose output files already exist are skipped. | `False` |

### Charge/Spin Options

| Option | Description | Default |
| --- | --- | --- |
| `-l, --ligand-charge TEXT` | Net charge or per-resname mapping used when `-q` is omitted (recommended). Triggers extract-style charge derivation on the full complex (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `-q, --charge INT` | Force the net system charge (overrides `--ligand-charge`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity forwarded to all downstream steps. | `1` |

### Extraction Options

| Option | Description | Default |
| --- | --- | --- |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). | Required for extraction |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). | `0.0` |
| `--include-H2O/--no-include-H2O` | Include waters (HOH/WAT/TIP3/SOL). | `True` |
| `--exclude-backbone/--no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `False` |
| `--add-linkH/--no-add-linkH` | Add link hydrogens for severed bonds. | `True` |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--freeze-links/--no-freeze-links` | Freeze link parents in pocket PDBs. | `True` |
| `--verbose/--no-verbose` | Enable INFO-level extractor logging. | `True` |

### MEP Search Options

| Option | Description | Default |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP search algorithm: GSM (Growing String Method) or DMF (Direct Max Flux). | `gsm` |
| `--max-nodes INT` | MEP internal nodes per segment. | `20` |
| `--max-cycles INT` | MEP maximum optimization cycles. | `300` |
| `--climb/--no-climb` | Enable TS climbing for the first segment. | `True` |
| `--opt-mode [grad\|hess]` | Workflow preset (`grad` → LBFGS/Dimer, `hess` → RFO/RSIRFO). For direct commands, prefer `opt --opt-mode grad|hess` and `tsopt --opt-mode grad|hess`. | `grad` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--preopt/--no-preopt` | Pre-optimize pocket endpoints before MEP search. | `True` |
| `--refine-path/--no-refine-path` | If True, run recursive `path-search`; if False, chain `path-opt` segments without recursive refinement. | `True` |

### MLIP Calculator Options

| Option | Description | Default |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only). | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | Shared MLIP Hessian engine. | `FiniteDifference` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |

### Post-Processing Options

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | Run TS optimization + IRC per reactive segment. | `False` |
| `--thermo/--no-thermo` | Run vibrational analysis (`freq`) on R/TS/P. | `False` |
| `--dft/--no-dft` | Run single-point DFT on R/TS/P. | `False` |
| `--opt-mode-post [grad\|hess]` | Optimizer preset override for TSOPT and post-IRC optimization (`grad` → Dimer/LBFGS, `hess` → RSIRFO/RFO). | `hess` |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--flatten/--no-flatten` | Enable surplus-imaginary-mode flattening in `tsopt`. | `False` |

TSOPT optimizer selection order: `--opt-mode-post` (if set) → `--opt-mode` (only when explicitly provided) → TSOPT default (`hess` → `rsirfo`).

Example: `--opt-mode grad --opt-mode-post hess` uses LBFGS for path optimization and RS-I-RFO for TS refinement.

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
| `-s, --scan-lists TEXT...` | Staged scans: `(i,j,target_Å)` tuples. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | _None_ |
| `--scan-one-based/--no-scan-one-based` | Force scan indexing to 1-based or 0-based. | _None_ |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | `0.20` |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV·Å⁻²). | `300` |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | `10000` |
| `--scan-preopt/--no-scan-preopt` | Override the scan preoptimization toggle. | _None_ |
| `--scan-endopt/--no-scan-endopt` | Override the scan end-of-stage optimization toggle. | _None_ |

## Outputs
```text
out_dir/ (default:./result_all/)
├─ summary.log # formatted summary for quick inspection
├─ summary.yaml # YAML version summary
├─ pockets/ # Per-input pocket PDBs when extraction runs
├─ scan/ # Staged pocket scan results (present when --scan-lists is provided)
├─ path_search/ # MEP results (GSM/DMF): trajectories, merged PDBs, diagrams, summary.yaml, per-segment folders
├─ path_search/post_seg_XX/ # Post-processing outputs (TS optimization, IRC, freq, DFT, diagrams)
└─ tsopt_single/ # TSOPT-only outputs with IRC endpoints and optional freq/DFT directories
```
- Console logs summarizing pocket charge resolution, YAML contents, scan stages, MEP progress (GSM/DMF), and per-stage timing.

### Energy diagram naming convention

Energy diagram files are named by method and scope:

| File name | Generated when | Content |
|---|---|---|
| `energy_diagram_MEP.png` | path-opt/path-search completes | All-segment MEP barriers (raw GSM/DMF values) |
| `energy_diagram_UMA.png` | per-segment tsopt+IRC completes | R→TS→P (MLIP energy) |
| `energy_diagram_G_UMA.png` | per-segment thermo completes | R→TS→P (MLIP Gibbs free energy) |
| `energy_diagram_DFT.png` | per-segment DFT completes | R→TS→P (DFT energy) |
| `energy_diagram_G_DFT_plus_UMA.png` | per-segment DFT+thermo completes | R→TS→P (DFT energy + MLIP thermal correction) |
| `energy_diagram_UMA_all.png` | all segments aggregated | All segments combined (MLIP) |
| `energy_diagram_G_UMA_all.png` | all segments + thermo | All segments combined (MLIP Gibbs) |
| `energy_diagram_DFT_all.png` | all segments + DFT | All segments combined (DFT) |
| `energy_diagram_G_DFT_plus_UMA_all.png` | all segments + DFT + thermo | All segments combined (DFT//MLIP Gibbs) |

### Reading `summary.log`
The log is organized into numbered sections:
- **[1] Global MEP overview** – image/segment counts, MEP trajectory plot paths, and the aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (MLIP path)** – per-segment barriers (`ΔE‡`), reaction energies (`ΔE`), and bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** – per-segment TS imaginary frequency checks, IRC outputs, and MLIP/thermo/DFT energy tables.
- **[4] Energy diagrams (overview)** – diagram tables for MEP/MLIP/Gibbs/DFT series plus an optional cross-method summary table.
- **[5] Output directory structure** – a compact tree of generated files with inline annotations.

### Reading `summary.yaml`
The YAML is a compact, machine-readable summary. Common top-level keys include:
- `out_dir`, `n_images`, `n_segments` – run metadata and total counts.
- `segments` – list of per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, and `bond_changes`.
- `energy_diagrams` (optional) – diagram payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, and `image` paths.

`summary.yaml` intentionally omits the formatted tables and filesystem tree that appear in `summary.log`.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Always provide `--ligand-charge` (numeric or per-residue mapping) when formal charges cannot be inferred so the correct net charge propagates to scan/MEP/TSOPT/DFT.
- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-full-pdb` option of `path-search` is intentionally hidden in this wrapper.
- Convergence presets: `--thresh` defaults to `gau`; `--thresh-post` defaults to `baker`.
- Extraction radii: passing `0` to `--radius` or `--radius-het2het` is internally clamped to `0.001 Å` by the extractor.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to the MEP/tsopt/freq/DFT stages; single-structure runs still require either `--scan-lists` or `--tsopt`.
- **`--resume`**: Re-run the same command with `--resume` to skip stages whose output files already exist. Each stage is guarded by sentinel-file checks (e.g. `summary.yaml` for MEP, `final_geometry.*` + `finished_irc_trj.xyz` for TSOPT/IRC, `R/`+`TS/`+`P/` directories for freq/DFT). When extraction is skipped on resume, provide `-q/--charge` or `--ligand-charge` explicitly so the charge can be resolved without re-running the extractor.


`all` supports layered YAML:

- `--config FILE`: base settings.

`defaults < config < CLI < override-yaml`

The effective YAML is forwarded to **every** invoked subcommand. Each tool reads the sections described in its own documentation:

| Subcommand | YAML Sections |
|------------|---------------|
| [`path-search`](path_search.md) | `geom`, `calc`, `gs`, `stopt`, `opt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

> **Note:** Applied after CLI values.

**Minimal example:**
```yaml
calc:
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 hessian_calc_mode: Analytical # recommended when VRAM permits
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

For a complete reference of all YAML options, see **[YAML Configuration Reference](yaml_reference.md)**.

---

## See Also

- [Installation](installation.md) — Setup and dependency installation
- [Getting Started](getting_started.md) — First run and workflow overview
- [Concepts & Workflow](concepts.md) — Mental model of pockets, segments, and stages
- [extract](extract.md) — Standalone pocket extraction (called internally by `all`)
- [path-search](path_search.md) — Standalone MEP search (called internally by `all`)
- [tsopt](tsopt.md) — Standalone TS optimization
- [freq](freq.md) — Standalone vibrational analysis
- [dft](dft.md) — Standalone DFT calculations
- [Common Error Recipes](recipes_common_errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Common errors and fixes
- [YAML Reference](yaml_reference.md) — Complete YAML configuration options
- [Glossary](glossary.md) — Definitions of MEP, TS, IRC, GSM, DMF
