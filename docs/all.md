# `all`

## Overview

`pdb2reaction all` runs the entire workflow end-to-end:

Active site model extraction → (optional) staged scan → MEP search (recursive `path-search` by default) → (optional) TS optimization + IRC (`tsopt`) → (optional) vibrational analysis / thermochemistry (`freq`) → (optional) single-point DFT (`dft`). Use `--refine-path False` to fall back to single-pass `path-opt` (GSM/DMF). The default MLIP backend is UMA; select an alternative with `-b/--backend`.

```{important}
The `all` workflow **without `--tsopt`** produces **TS candidates** (Highest-Energy Images from MEP search). Adding `--tsopt` refines these into optimized TS structures validated by imaginary-frequency check, followed by IRC for endpoint validation. Always inspect the results (imaginary-frequency count + endpoint connectivity) before mechanistic interpretation.
```

### At a glance
- **Use when:** End-to-end pipeline from PDB(s) (extraction → MEP → TS optimization → IRC → thermo → DFT).
- **Method:** Three modes — multi-structure MEP, single-structure + staged scan, or TSOPT-only — selected by the inputs and flags you provide.
- **Outputs:** `summary.log`, `summary.json`, and `path_search/mep.pdb` (or `path_opt/` when `--refine-path False`); per-segment `seg_XX/` and post-processing `path_search/post_seg_XX/` when `--tsopt`/`--thermo`/`--dft` are enabled.
- **Defaults:** Backend `uma`, `--mep-mode gsm`, `--opt-mode grad`, `--refine-path True`, `--preopt True`, `--thresh gau`, `--thresh-post baker`; `--tsopt`/`--thermo`/`--dft` are off.
- **Next step:** Without `--tsopt`, results are TS *candidates* (HEIs); add `--tsopt` (imaginary-frequency check) and IRC for validation, then optionally `--thermo` and `--dft`.

## Workflow at a glance

Most workflows follow this flow:

```text
Full system(s) (PDB/XYZ/GJF)
 │
 ├─ (optional) active site model extraction [`extract`](extract.md) ← requires PDB when you use --center/-c
 │ ↓
 │ Active site model/cluster model(s) (PDB)
 │ │
 │ ├─ (optional) staged scan [`scan`](scan.md) ← single-structure workflows
 │ │ ↓
 │ │ Ordered intermediates
 │ │ ↓
 │ └─ MEP search [`path-search`](path-search.md) or [`path-opt`](path-opt.md)
 │ ↓
 │ MEP trajectory (mep_trj.xyz) + energy diagrams
 │ ↓
 └─ (optional) TS optimization + IRC [`tsopt`](tsopt.md) → [`irc`](irc.md)
 └─ (optional) thermo [`freq`](freq.md)
 └─ (optional) single-point DFT [`dft`](dft.md)
```

Each stage is available as an individual subcommand. The `pdb2reaction all` command runs many stages end-to-end.

It supports three modes:

- **Multi-structure workflow** — Provide ≥2 structures (PDB/GJF/XYZ) in reaction order plus a substrate definition. `all` extracts active site models, runs GSM/DMF MEP search, merges the optimized path back into the full-system template(s), and optionally runs TSOPT+IRC/freq/DFT per reactive segment.
- **Single-structure + staged scan** — Provide one structure plus one or more `--scan-lists/-s`. The (staged) scan generates an ordered set of intermediates that become MEP endpoints.

  - One `--scan-lists/-s` literal runs a single scan stage.
  - Multiple stages are passed as multiple arguments to a single `--scan-lists/-s` (e.g. `-s '[(…)]' '[(…)]'`).

- **TSOPT-only active site model TS optimization** — Provide a single input structure, omit `--scan-lists/-s`, and set `--tsopt`. `all` extracts the active site model (if `-c/--center` is given) and runs TS optimization + IRC, with optional freq/DFT, on that single system.

```{tip}
For large active site models, the single-structure scan workflow (`--scan-lists/-s`) tends to produce more reliable reaction barriers than the multi-structure MEP workflow. When multiple full PDB structures are provided, structural differences in regions unrelated to the reaction coordinate can accumulate, leading to overestimated barriers. The scan workflow avoids this by starting from a single structure and driving only the relevant coordinates, minimizing irrelevant structural noise. This effect becomes more pronounced as the model size increases.
```

> **Working examples:** The [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) directory contains complete `all` workflow scripts for GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)), covering both multi-structure MEP and scan-based pipelines.

## Minimal example

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 --out-dir ./result_all
```

## Output checklist

- `result_all/summary.log`
- `result_all/summary.json`
- `result_all/path_search/mep.pdb` (or `result_all/path_opt/` when `--refine-path False` is used)

## Common examples

1. Run full post-processing in one command.

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_mep
```

2. Single-structure staged scan route.

```bash
pdb2reaction all -i 1.R.pdb -c "SAM,GPP,MG" -l "SAM:1,GPP:-3" \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --tsopt --thermo --out-dir ./result_scan
```

PDB/GJF companion files are generated when templates are available, controlled by `--convert-files` (enabled by default).


## Usage
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
```

For help output, `pdb2reaction all --help` shows core options and `pdb2reaction all --help-advanced` shows the full option list.

### Examples
```bash
# Multi-structure MEP with explicit ligand charges and post-processing
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' \
 -l 'SAM:1,GPP:-3' --multiplicity 1 --freeze-links \
 --max-nodes 10 --max-cycles 100 --climb --opt-mode grad \
 --out-dir ./result_mep --tsopt --thermo --dft

# Single-structure staged scan followed by GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --opt-mode hess --tsopt --thermo --dft

# TSOPT-only workflow (no path search)
pdb2reaction all -i TS_candidate.pdb -c 'SAM,GPP,MG' \
 -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```

## Workflow

0. **Preflight checks** (automatic)
 - `all` automatically runs `add-elem-info` (fills missing element symbols in PDB columns 77–78) and `fix-altloc` (resolves alternate conformations) on every PDB input before any other processing. When using individual subcommands (e.g., `extract`, `opt`), you must run these manually if needed.

1. **Active site model (binding pocket) extraction** (if `-c/--center` is provided)
 - Substrates may be specified via PDB paths, residue IDs (`123,124` or `A:123,B:456`), or residue names (`GPP,SAM`).
 - Optional toggles forward to the extractor: `--radius`, `--radius-het2het`, `--include-h2o`, `--exclude-backbone`, `--add-linkh`, `--selected-resn`, and `--verbose`.
 - Per-input active site model PDBs are saved under `<out-dir>/models/`. When multiple structures are supplied, their active site models are unioned per residue selection.
 - The **first active site model’s net charge** is propagated to scan/MEP/TSOPT.

2. **Optional staged scan (single-input only)**
 - Each `--scan-lists/-s` argument is a Python-like list of `(i,j,target_Å)` tuples describing an MLIP scan stage. Atom indices refer to the original input ordering (1-based) and are remapped to the active site model ordering. For PDB inputs, `i`/`j` can be integer indices or selector strings like `'TYR,285,CA'`; selectors accept spaces/commas/slashes/backticks/backslashes (` ` `,` `/` `` ` `` `\`) as delimiters and allow unordered tokens (fallback assumes resname, resseq, atom).
 - A single literal runs a one-stage scan; multiple literals run **sequentially** so stage 2 begins from stage 1's result, and so on. Supply multiple literals as arguments to a single `--scan-lists/-s` (e.g. `-s '[(…)]' '[(…)]'`).
 - Stage endpoints (`stage_XX/result.pdb`) become the ordered intermediates that feed the subsequent MEP step.

3. **MEP search on active site models (recursive `path-search`)**
 - By default, runs recursive `path-search`, which automatically detects multistep reactions and builds a detailed MEP for each elementary step (outputs go to `<out-dir>/path_search/`).  Complex multistep mechanisms may require manual trial-and-error to obtain a satisfactory pathway.
 - Use `--refine-path False` to fall back to single-pass `path-opt` GSM/DMF on each adjacent pair (outputs go to `<out-dir>/path_opt/`).
 - For multi-input PDB runs, the full-system templates are automatically passed to `path-search` for reference merging. Single-structure scan runs reuse the original full PDB template for every stage.

4. **Merge active site models back to the full systems** (default with `--refine-path`)
 - When `--refine-path` is True (default) and reference PDB templates exist, merged `mep_w_ref*.pdb` and per-segment `mep_w_ref_seg_XX.pdb` files are emitted under `<out-dir>/path_search/`. With `--refine-path False` (`path-opt` mode), full-system merge is not performed.

5. **Optional per-segment post-processing** (only for reactive segments — segments with bond changes; bridge segments are skipped)
 - `--tsopt`: run TS optimization on each HEI active site model, follow with EulerPC-based IRC, then re-optimize IRC endpoints with `--thresh-post` (default `baker`). The endpoint optimization working directory is automatically deleted after completion.
 - `--thermo`: call `freq` on (R, TS, P) to obtain vibrational/thermochemistry data and an MLIP Gibbs diagram.
 - `--dft`: launch single-point DFT on (R, TS, P) and build a DFT diagram. When combined with `--thermo`, a DFT//MLIP Gibbs diagram (DFT energies + MLIP thermal correction) is also produced.
  - Shared overrides include `--opt-mode`, `--opt-mode-post` (overrides TSOPT/post-IRC optimization mode), `--flatten/--no-flatten`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, and `--dft-engine` (GPU-first by default).
 - For Hessian evaluation modes, see {ref}`hessian-evaluation`.

6. **TSOPT-only mode** (single input, `--tsopt`, no `--scan-lists/-s`)
 - Skips the MEP/merge stages. Runs `tsopt` on the active site model (or full input if extraction is skipped), performs EulerPC IRC, identifies the higher-energy endpoint as reactant (R), and generates the same set of energy diagrams plus optional freq/DFT outputs.

### Charge and spin precedence

Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>` for details). In the `all` command, charge derivation from active site model extraction (when `-c` is specified) acts as an additional priority layer.

**Spin resolution:** `--multiplicity` (CLI) → `.gjf` template → default (1)

> **Tip:** Always provide `--ligand-charge/-l` for non-standard substrates to ensure correct charge propagation.

### Input expectations
- Extraction enabled (`-c/--center`): inputs must be **PDB** files so residues can be located.
- Extraction skipped: inputs may be **PDB/XYZ/GJF**.
- Multi-structure runs require ≥2 structures.

## CLI Options

> **Note:** Default values shown are used when the option is not specified.

### Input/Output Options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full structures in reaction order (single input allowed only with `--scan-lists/-s` or `--tsopt`). | Required |
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
| `-q, --charge INT` | Force the net system charge (overrides `--ligand-charge/-l`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity forwarded to all downstream steps. | `1` |

### Extraction Options

| Option | Description | Default |
| --- | --- | --- |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). | Required for extraction |
| `-r, --radius FLOAT` | Active site model inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). Passing `0` is internally nudged to `0.001 Å` to avoid empty selections (same behavior as standalone `extract`). | `0.0` |
| `--include-h2o/--no-include-h2o` | Include waters (HOH/WAT/TIP3/SOL). | `True` |
| `--exclude-backbone/--no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `False` |
| `--add-linkh/--no-add-linkh` | Add link hydrogens for severed bonds. | `True` |
| `--selected-resn TEXT` | Residues to force include. **Despite the name, this flag accepts residue IDs (colon-separated integers with optional chains/insertion codes, e.g. `A:123A`), not 3-letter residue names.** Use `-c/--center 'GPP,SAM'` for residue-name-based selection. | `""` |
| `--modified-residue TEXT` | Comma-separated residue names (with optional charge) to treat as amino acids for backbone truncation and charge assignment (e.g., `HD1,HD2,HD3` or `HD1:0,SEP:-2`). | `""` |
| `--freeze-links/--no-freeze-links` | Freeze link parents in active site model PDBs. | `True` |
| `--verbose/--no-verbose` | Enable INFO-level extractor logging. | `True` |

### MEP Search Options

| Option | Description | Default |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP search algorithm: GSM (Growing String Method) or DMF (Direct Max Flux). | `gsm` |
| `--max-nodes INT` | MEP internal nodes per segment. **GSM:** total images = `max_nodes + 2` (endpoints fixed). **DMF:** number of *movable* images along the chain (no implicit endpoint expansion). | `20` |
| `--max-cycles INT` | MEP maximum optimization cycles. | `300` |
| `--climb/--no-climb` | Enable TS climbing for the first segment. | `True` |
| `--opt-mode [grad\|hess]` | Workflow preset (`grad` → LBFGS/Dimer, `hess` → RFO/RSIRFO). For direct commands, prefer `opt --opt-mode grad|hess` and `tsopt --opt-mode grad|hess`. The token-to-algorithm mapping depends on the scope; see {ref}`opt-mode-semantics` for the per-subcommand table and note that `all`'s pre-opt default (`grad`) is not the same as `tsopt`'s default (`hess`). | `grad` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--preopt/--no-preopt` | Pre-optimize active site model endpoints before MEP search. **Note:** `all` overrides the child-subcommand default here. Standalone `path-search`, `path-opt`, `scan`, `scan2d`, and `scan3d` default `--preopt` to `False`. | `True` |
| `--refine-path BOOL` | If `True` (default), run recursive `path-search`; if `False`, chain `path-opt` segments without recursive refinement. Pass `--refine-path False` to disable (no `--no-refine-path` auto-form: `click.BOOL` parameter, not a flag). | `True` |

### MLIP Calculator Options

| Option | Description | Default |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
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

```{warning}
The `--dft` single-point calculations (powered by PySCF/GPU4PySCF) are very expensive for models exceeding ~300 atoms. For such systems, HPC clusters with high-end GPUs (e.g. A100, H200) are typically required.
```

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
| `--dft-engine [gpu\|cpu]` | DFT backend: gpu (GPU4PySCF) or cpu (PySCF). In the `all` wrapper this is named `--dft-engine` (prefix-disambiguated); on the standalone `dft` subcommand the same option is named `--engine`. | `gpu` |
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
├─ summary.log                  # Text summary
├─ summary.json                 # JSON results
├─ models/                      # Extracted active site model PDBs when extraction runs
│  └─ model_<input_stem>.pdb   #   One per -i input
├─ scan/                        # Staged scan results (present when --scan-lists is provided)
├─ seg_XX/                      # Refined TS and optimized IRC endpoints of segment XX
│  ├─ reactant.{pdb,xyz,gjf}   #   Output format matches input format
│  ├─ ts.{pdb,xyz,gjf}
│  └─ product.{pdb,xyz,gjf}
├─ path_search/                 # MEP results (recursive path-search, default); path_opt/ when --refine-path False
│  └─ post_seg_XX/              # Post-processing: TSOPT, IRC, freq, DFT per segment
│     ├─ structures/            # Optimized R/TS/P structures (IRC endpoints)
│     ├─ irc/                   # IRC trajectories and plots
│     ├─ ts/                    # TS optimization output and vibrational analysis
│     ├─ freq/{R,TS,P}/         # Per-state frequencies_cm-1.txt + thermoanalysis.yaml
│     └─ dft/                   # DFT single-point results (when --dft is enabled)
└─ tsopt_single/                # TSOPT-only outputs with IRC endpoints
   ├─ structures/               # Working copy of R/TS/P (canonical version lives at top-level seg_NN/)
   ├─ ts/                       # TS optimization output and vibrational analysis
   ├─ irc/                      # IRC trajectories and plots
   ├─ freq/{R,TS,P}/            # Per-state frequencies_cm-1.txt + thermoanalysis.yaml
   └─ dft/                      # DFT single-point results (when --dft is enabled)
```

```{note}
**`seg_XX/` vs `path_search/post_seg_XX/`.** The two per-segment trees serve different purposes and are **not** duplicates:

- **`seg_XX/`** (top level of `out_dir`) is the *segment-level merged result* aggregating the final reactant/TS/product structures for reactive segment `XX`. It is written after the full post-processing pipeline converges and holds the canonical `reactant.{pdb,xyz,gjf}`, `ts.{pdb,xyz,gjf}`, `product.{pdb,xyz,gjf}` you should cite when reporting mechanisms.
- **`path_search/post_seg_XX/`** is the *per-segment post-processing workspace*. It stores the intermediate products of each stage (`ts/`, `irc/`, `structures/`, `freq/`, `dft/`) and is the right place to inspect when debugging a single stage (e.g. checking `ts/vib/imag_*_trj.xyz` or `irc/*_trj.xyz`). Reactive segments get populated; bridge segments (no bond change) are skipped.

When `--refine-path False` is passed, the workspace moves under `path_opt/post_seg_XX/` instead.
```

- Console logs summarizing active site model charge resolution, YAML contents, scan stages, MEP progress (GSM/DMF), and per-stage timing.

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

### Reading `summary.json`
The JSON summary contains structured data. Common top-level keys include:
- `out_dir`, `n_images`, `n_segments` – run metadata and total counts.
- `segments` – list of per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, and `bond_changes`.
- `energy_diagrams` (optional) – diagram payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, and `image` paths.

`summary.json` intentionally omits the formatted tables and filesystem tree that appear in `summary.log`.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Always provide `--ligand-charge/-l` (numeric or per-residue mapping) when formal charges cannot be inferred so the correct net charge propagates to scan/MEP/TSOPT/DFT.
- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-full-pdb` option of `path-search` is hidden in this wrapper.
- Convergence presets: `--thresh` defaults to `gau`; `--thresh-post` defaults to `baker`.
- Extraction radii: passing `0` to `--radius` or `--radius-het2het` is internally clamped to `0.001 Å` by the extractor.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to the MEP/tsopt/freq/DFT stages; single-structure runs still require either `--scan-lists/-s` or `--tsopt`.
- **`--resume`**: Re-run the same command with `--resume` to skip stages whose output files already exist. Each stage is guarded by sentinel-file checks (e.g. `summary.json` for MEP, `final_geometry.*` + `finished_irc_trj.xyz` for TSOPT/IRC, `R/`+`TS/`+`P/` directories for freq/DFT). When extraction is skipped on resume, provide `-q/--charge` or `--ligand-charge/-l` explicitly so the charge can be resolved without re-running the extractor. **Sentinel-corruption caveat:** `--resume` only checks for the *presence* of the sentinel files, not their integrity. If a stage was killed mid-write (SIGKILL, OOM, cluster preemption) and the sentinel file was already on disk but is truncated or corrupted, `--resume` will still consider the stage complete. Delete the stage's output directory (e.g. `path_search/post_seg_XX/ts/`) before resuming if you suspect a partially written result.


`all` supports layered YAML:

- `--config FILE`: base settings.

`defaults < config < CLI`

The effective YAML is forwarded to **every** invoked subcommand. Each tool reads the sections described in its own documentation:

| Subcommand | YAML Sections |
|------------|---------------|
| [`path-search`](path-search.md) | `geom`, `calc`, `gs`, `stopt`, `opt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |
| [`irc`](irc.md) | `geom`, `calc`, `irc` |

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

For a complete reference of all YAML options, see **[YAML Configuration Reference](yaml-reference.md)**.

## See Also

- [Installation](installation.md) — Setup and dependency installation
- [Getting Started](getting-started.md) — First run, workflow overview, and key concepts
- [extract](extract.md) — Standalone active site model extraction (called internally by `all`)
- [scan](scan.md) — Standalone staged distance scan
- [path-opt](path-opt.md) — Single-pass MEP optimization (GSM/DMF)
- [path-search](path-search.md) — Recursive MEP search (called internally by `all`)
- [tsopt](tsopt.md) — Standalone TS optimization
- [irc](irc.md) — Standalone IRC calculation
- [freq](freq.md) — Standalone vibrational analysis
- [dft](dft.md) — Standalone DFT calculations
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Common errors and fixes
- [YAML Reference](yaml-reference.md) — Complete YAML configuration options
- [Glossary](glossary.md) — Definitions of MEP, TS, IRC, GSM, DMF
