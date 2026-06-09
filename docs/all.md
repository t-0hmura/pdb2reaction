# `all`

`pdb2reaction all` runs the entire workflow end-to-end so you can go from structures to a validated mechanism in one command, instead of chaining `extract` → `scan` / `path-search` → `tsopt` → `irc` / `freq` / `dft` by hand. Starting from one or more PDB inputs, it extracts an active-site cluster model, runs an optional staged scan, performs an MEP search (recursive `path-search` by default; `--refine-path False` falls back to single-pass `path-opt`), and optionally chains TS optimization, IRC, vibrational analysis, and single-point DFT. The default MLIP backend is UMA; choose an alternative with `-b/--backend`.

`all` runs in one of three modes, chosen by what you pass:

- **Multi-structure MEP** (`[mode] all (mep)`) — give ≥ 2 structures in reaction order plus a substrate definition. `all` extracts active-site models, runs GSM / DMF MEP search, merges the optimized path back into the full-system template(s), and optionally runs TSOPT + IRC / freq / DFT per reactive segment.
- **Single-structure staged scan** (`[mode] all (scan-lists)`) — give one structure plus one or more `--scan-lists/-s` literals, each defining a scan stage; the staged scan produces the ordered intermediates that drive the MEP step.
- **TSOPT-only** — give a single input and set `--tsopt` (no `--scan-lists`). `all` skips the MEP / merge stages, runs `tsopt` + EulerPC IRC on the active-site model (or the full input if extraction is skipped), and identifies the higher-energy endpoint as the reactant.

```{important}
Without `--tsopt`, the workflow produces **TS candidates** (highest-energy images from MEP search). Adding `--tsopt` refines them into optimized TS structures validated by an imaginary-frequency check, followed by IRC for endpoint validation. Always inspect the results (imaginary-frequency count and IRC endpoint connectivity) before mechanistic interpretation.
```

## Examples

Working examples for GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)) covering both multi-structure MEP and scan-based pipelines: [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples).

Command form:

```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [-b uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
```

Multi-structure MEP with TS + thermo + DFT:

```bash
# Multi-structure MEP with TS + thermo + DFT
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --out-dir ./result_mep
```

Single-structure staged scan (two stages):

```bash
# Single-structure staged scan (two stages)
pdb2reaction all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
    --tsopt --thermo --out-dir ./result_scan
```

TSOPT-only validation of a single TS candidate:

```bash
# TSOPT-only validation of a single TS candidate
pdb2reaction all -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

PDB / GJF companion files are generated automatically when reference templates are available; control with `--convert-files` (on by default).

`pdb2reaction all --help` shows core options; `pdb2reaction all --help-advanced` shows the full option list.

## Workflow

```text
Full system(s) (PDB / XYZ / GJF)
  ├─ (optional) active-site extraction `extract` — requires PDB when `-c` is used
  │   └─ active-site cluster model(s)
  │       ├─ (optional) staged scan `scan` — single-structure workflows
  │       │   └─ ordered intermediates
  │       └─ MEP search `path-search` (recursive) or `path-opt`
  │           └─ MEP trajectory `mep_trj.xyz` + energy diagrams
  └─ (optional) TS optimisation + IRC `tsopt` → `irc`
      ├─ (optional) thermochemistry `freq`
      └─ (optional) single-point DFT `dft`
```

`all` runs the following stages in order. Stages 0 and 1 are automatic preprocessing; the rest fire based on the flags you pass.

0. **Preflight** (automatic) — `add-elem-info` fills missing element symbols and runs only on PDB inputs that lack them (an empty element field in cols 77–78), and `fix-altloc` resolves alternate conformations and runs only on PDB inputs that contain alternate conformations (altLoc). When you invoke individual subcommands (e.g. `extract`, `opt`) you must run these manually if needed.
1. **Active-site model (binding-pocket) extraction** (when `-c/--center` is set) — runs `extract` to build the active-site cluster from the substrate selection (see Inputs for the `-c/--center` syntax). Forwarded extractor toggles: `--radius`, `--radius-het2het`, `--include-h2o`, `--exclude-backbone`, `--add-linkh`, `--selected-resn`, `--verbose`. Per-input PDBs are saved under `<out-dir>/_work/models/`; when multiple structures are supplied, the active-site models are unioned per residue selection. The **first active-site model's net charge** is propagated to scan / MEP / TSOPT.
2. **Optional staged scan** (single-input only) — each `--scan-lists/-s` literal is a list of `(i, j, target_Å)` tuples. Atom indices use the original input ordering (1-based) and are remapped to the active-site model ordering. PDB selector strings like `'TYR,285,CA'` are also accepted (space / comma / slash / backtick / backslash delimiters; token order is flexible). Stages run sequentially (stage 2 starts from stage 1's result), and the stage endpoints become the ordered intermediates that feed the MEP step.
3. **MEP search** — by default runs recursive `path-search`, which automatically detects multistep reactions and builds a detailed MEP per elementary step. Complex multistep mechanisms may need manual trial-and-error to converge a satisfactory pathway. `--refine-path False` falls back to single-pass `path-opt` GSM / DMF on each adjacent pair. The raw engine output is written under `<out-dir>/_work/path_search/` (or `_work/path_opt/`); the merged products (`mep.pdb`, `mep_trj.xyz`, `energy_diagram_MEP.png`) are promoted to the top level. For multi-input runs, full-system PDB templates are forwarded automatically for reference merging.
4. **Merge to full systems** (default with `--refine-path`) — when reference templates exist, the merged `mep_w_ref.pdb` is promoted to `<out-dir>/`, and per-segment `mep_w_ref_seg_NN.pdb` files remain under `<out-dir>/_work/path_search/`. `--refine-path False` skips the full-system merge.
5. **Per-segment post-processing** (reactive segments only — bridge segments without bond changes are skipped):
   - `--tsopt` — TS optimization on each HEI active-site model, followed by EulerPC IRC, then IRC-endpoint re-optimization with `--thresh-post` (default `baker`). The endpoint optimization working directory is deleted automatically after completion.
   - `--thermo` — `freq` on (R, TS, P) for vibrational + thermochemistry data and an MLIP Gibbs diagram.
   - `--dft` — single-point DFT on (R, TS, P) and a DFT diagram. With `--thermo`, a DFT//MLIP Gibbs diagram (DFT energies + MLIP thermal correction) is also produced.
   - Shared overrides: `--opt-mode`, `--opt-mode-post`, `--flatten`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, `--dft-engine` (GPU-first by default). For Hessian evaluation modes see {ref}`hessian-evaluation`.
6. **TSOPT-only mode** (single input + `--tsopt`, no `--scan-lists`) — skips MEP / merge; runs `tsopt` + EulerPC IRC and generates the same energy diagrams plus optional freq / DFT outputs.

## Outputs

The tree has three top-level zones: **deliverables at the root**, **per-segment deliverables under `segments/seg_NN/`**, and **pipeline scratch under `_work/`** (safe to `rm -rf` once you have the results you need).

```text
out_dir/   (default: ./result_all/)
├─ summary.log                 # Text summary (authored at the root)
├─ summary.json                # JSON results
├─ mep.pdb                     # Merged MEP path (promoted from the engine)
├─ mep_w_ref.pdb               # MEP merged into the full-system template (multi-input runs)
├─ mep_trj.xyz                 # Full MEP trajectory
├─ energy_diagram_MEP.png      # All-segment MEP barriers
├─ energy_diagram_*.png        # Aggregated post-processing diagrams (UMA / Gibbs / DFT, with --tsopt etc.)
├─ segments/                   # Per-reactive-segment deliverables (bridge segments are skipped)
│  └─ seg_NN/                  # 2-digit index, e.g. seg_01, seg_02
│     ├─ reactant.{pdb,xyz,gjf}   # Canonical R/TS/P (output format matches input format)
│     ├─ ts.{pdb,xyz,gjf}
│     ├─ product.{pdb,xyz,gjf}
│     ├─ ts/                   # TS optimisation output + vibrational analysis (--tsopt)
│     ├─ irc/                  # IRC trajectories + plots (--tsopt)
│     ├─ freq/{R,TS,P}/        # frequencies_cm-1.txt + thermoanalysis.yaml (--thermo)
│     └─ dft/                  # DFT single-point results (--dft)
└─ _work/                      # Pipeline scratch (safe to delete)
   ├─ models/                  # Extracted active-site model PDBs (model_<input_stem>.pdb, when extraction runs)
   ├─ scan/                    # Staged scan results (with --scan-lists)
   ├─ add_elem_info/           # Preflight element-symbol fills
   ├─ fix_altloc/              # Preflight altLoc resolution
   └─ path_search/             # Raw MEP-engine output (path_opt/ when --refine-path False)
```

In **TSOPT-only mode** (single input + `--tsopt`, no `--scan-lists`) there is no MEP stage: the optimized R/TS/P plus `ts/`, `irc/`, `freq/`, and `dft/` land directly under `segments/seg_01/`, and `_work/path_search/` is absent.

```{note}
**The canonical structures are `segments/seg_NN/reactant.*`, `ts.*`, `product.*`** — cite these when reporting mechanisms. The `ts/`, `irc/`, `freq/`, and `dft/` subdirectories inside the same `seg_NN/` hold the per-stage working files (e.g. `ts/vib/imag_*_trj.xyz`, `irc/*_trj.xyz`) for debugging a single stage. The raw MEP-search engine output under `_work/path_search/` is scratch — the products you need (`mep.pdb`, `mep_trj.xyz`, `energy_diagram_MEP.png`) are already promoted to the root.
```

At `-v 2` the console summarises active-site charge resolution, YAML contents, scan stages, MEP progress (GSM / DMF), and per-stage timing; see {ref}`verbosity-levels`.

### Plot file naming

Energy-diagram filenames encode method and scope:

| File | Generated when | Content |
|---|---|---|
| `energy_diagram_MEP.png` | `path-opt` / `path-search` completes | All-segment MEP barriers (raw GSM / DMF values) |
| `energy_diagram_UMA.png` | per-segment `tsopt` + IRC completes | R → TS → P (MLIP energy) |
| `energy_diagram_G_UMA.png` | per-segment thermo completes | R → TS → P (MLIP Gibbs) |
| `energy_diagram_DFT.png` | per-segment DFT completes | R → TS → P (DFT energy) |
| `energy_diagram_G_DFT_plus_UMA.png` | per-segment DFT + thermo | R → TS → P (DFT energy + MLIP thermal correction) |
| `energy_diagram_UMA_all.png` / `_G_UMA_all.png` / `_DFT_all.png` / `_G_DFT_plus_UMA_all.png` | all segments aggregated (variants for MLIP / Gibbs / DFT / DFT//MLIP Gibbs) | Combined across all segments |
| `irc_plot.png` (per `segments/seg_NN/irc/`) | per-segment IRC completes | IRC profile (MLIP energy along the trajectory) |
| `irc_plot_all.png` | all segments aggregated | IRC profiles concatenated across segments |

### Reading `summary.log`

The log is organised into numbered sections:

- **[1] Global MEP overview** — image / segment counts, MEP trajectory plot paths, aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (UMA path)** — per-segment barriers (ΔE‡), reaction energies (ΔE), bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** — TS imaginary-frequency checks, IRC outputs, MLIP / thermo / DFT energy tables.
- **[4] Energy diagrams (overview)** — diagram tables for MEP / MLIP / Gibbs / DFT plus an optional cross-method summary.
- **[5] Output directory structure** — a compact tree of generated files with inline annotations.

### Reading `summary.json`

Top-level keys: `out_dir`, `n_images`, `n_segments` (run metadata and counts); `segments` (per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes`); `energy_diagrams` (optional payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` paths). `summary.json` intentionally omits the formatted tables and filesystem tree from `summary.log`.

## CLI options

Defaults shown are used when the option is not specified. The full flag list is in the generated [command reference](reference/commands/index.md); the tables below cover the options that need explanation.

Input expectations:

- Extraction enabled (`-c/--center`): inputs must be **PDB** so residues can be located.
- Extraction skipped: inputs may be **PDB / XYZ / GJF**.
- Multi-structure runs require ≥ 2 structures. For full input-file requirements (hydrogens, element columns, atom-order parity), see [CLI Conventions](cli-conventions.md).

Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>`). In `all`, the charge derivation from active-site model extraction (when `-c` is set) acts as an additional priority layer. Spin resolution: `--multiplicity` CLI → `.gjf` template → default `1`. Always provide `--ligand-charge/-l` for non-standard substrates so the correct net charge propagates to scan / MEP / TSOPT / DFT.

### Input / output

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full structures in reaction order (single input allowed only with `--scan-lists/-s` or `--tsopt`). | Required |
| `--ref-pdb FILE` | Reference PDB for topology when `-i` provides XYZ inputs. | _None_ |
| `-o, --out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--convert-files / --no-convert-files` | Global toggle for XYZ / TRJ → PDB / GJF companions when templates are available. | `True` |
| `--dump / --no-dump` | Dump MEP (GSM / DMF) trajectories. Always forwarded to `path-search` / `path-opt`; forwarded to `scan` / `tsopt` only when explicitly set. `freq` defaults to `dump=True` unless you pass `--no-dump`. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config / --no-show-config` | Print resolved configuration before execution. | `False` |
| `--dry-run / --no-dry-run` | Validate and print plan without running stages. | `False` |

### Charge / spin

| Option | Description | Default |
| --- | --- | --- |
| `-l, --ligand-charge TEXT` | Net charge or per-resname mapping used when `-q` is omitted (recommended). Triggers extract-style charge derivation on the full complex (PDB inputs, or XYZ / GJF with `--ref-pdb`). | _None_ |
| `-q, --charge INT` | Force the net system charge (overrides `--ligand-charge/-l`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity forwarded to all downstream steps. | `1` |

### Extraction

| Option | Description | Default |
| --- | --- | --- |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). | Required for extraction |
| `-r, --radius FLOAT` | Active-site model inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). `0` is internally nudged to `0.001 Å` to avoid empty selections (same as standalone `extract`). | `0.0` |
| `--include-h2o / --no-include-h2o` | Include waters (HOH / WAT / TIP3 / SOL). | `True` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `False` |
| `--add-linkh / --no-add-linkh` | Add link hydrogens for severed bonds. | `True` |
| `--selected-resn TEXT` | Residues to force include. **Despite the name, accepts residue IDs (colon-separated integers with optional chains / insertion codes, e.g. `A:123A`), not 3-letter residue names.** Use `-c/--center 'GPP,SAM'` for residue-name selection. | `""` |
| `--modified-residue TEXT` | Comma-separated residue names (with optional charge) to treat as amino acids for backbone truncation and charge assignment (e.g. `HD1,HD2,HD3` or `HD1:0,SEP:-2`). | `""` |
| `--freeze-links / --no-freeze-links` | Freeze link parents in active-site model PDBs. | `True` |

### MEP search

| Option | Description | Default |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP algorithm: GSM (Growing String Method) or DMF (Direct Max Flux). | `gsm` |
| `--max-nodes INT` | MEP internal nodes per segment. **GSM**: total images = `max_nodes + 2` (endpoints fixed). **DMF**: number of *movable* images along the chain (no implicit endpoint expansion). | `20` |
| `--max-cycles INT` | MEP maximum optimization cycles. | `300` |
| `--climb / --no-climb` | Enable climbing image for standard GSM segments (bridge segments always disable climbing). | `True` |
| `--opt-mode [grad\|hess]` | Workflow preset (`grad` → L-BFGS / Dimer, `hess` → RFO / RS-I-RFO). Token-to-algorithm mapping depends on scope — see {ref}`opt-mode-semantics` for the per-subcommand table; note that `all`'s pre-opt default (`grad`) differs from `tsopt`'s default (`hess`). | `grad` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--preopt / --no-preopt` | Pre-optimize active-site model endpoints before MEP search. Standalone `scan` / `scan2d` / `scan3d` default `--preopt` to `False`. | `True` |
| `--refine-path / --no-refine-path` | On (default) → recursive `path-search`; off → chain `path-opt` segments without recursive refinement. | `True` |

### MLIP calculator

| Option | Description | Default |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only). See {ref}`workers-fd-downgrade`. | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | Shared MLIP Hessian engine. | `FiniteDifference` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |

### Post-processing

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt / --no-tsopt` | Run TS optimization + IRC per reactive segment. | `False` |
| `--thermo / --no-thermo` | Run vibrational analysis (`freq`) on R / TS / P. | `False` |
| `--dft / --no-dft` | Run single-point DFT on R / TS / P. | `False` |
| `--opt-mode-post [grad\|hess]` | Optimizer preset for TSOPT + post-IRC (`grad` → Dimer / L-BFGS, `hess` → RS-I-RFO / RFO). | `hess` |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations. | `baker` |
| `--flatten / --no-flatten` | Enable surplus-imaginary-mode flattening in `tsopt`. | `False` |

```{warning}
`--dft` single-point calculations (PySCF / GPU4PySCF) are very expensive for models above ~300 atoms — HPC clusters with high-end GPUs (e.g. A100, H200) are typically required.
```

TSOPT optimizer selection order: `--opt-mode-post` (if set) → `--opt-mode` (only when explicitly provided) → TSOPT default (`hess` → `rsirfo`). Example: `--opt-mode grad --opt-mode-post hess` uses L-BFGS for path optimization and RS-I-RFO for TS refinement.

### TSOPT / freq / DFT / scan overrides

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles`. | `10000` |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory. | _None_ |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Maximum modes to write. | `10` |
| `--freq-amplitude-ang FLOAT` | Mode animation amplitude (Å). | `0.8` |
| `--freq-n-frames INT` | Frames per mode animation. | `20` |
| `--freq-sort [value\|abs]` | Mode sorting behavior. | `value` |
| `--freq-temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--freq-pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dft-engine [gpu\|cpu]` | DFT backend (GPU4PySCF or PySCF). In `all` the option is named `--dft-engine`; the standalone `dft` subcommand uses `--engine`. | `gpu` |
| `--dft-out-dir PATH` | DFT outputs base directory override. | _None_ |
| `--dft-func-basis TEXT` | Functional / basis pair. | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | Maximum SCF iterations. | `100` |
| `--dft-conv-tol FLOAT` | SCF convergence tolerance. | `1e-9` |
| `--dft-grid-level INT` | PySCF grid level. | `3` |
| `-s, --scan-lists TEXT...` | Staged scans: `(i, j, target_Å)` tuples (single-input runs). | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | _None_ |
| `--scan-one-based / --no-scan-one-based` | Force scan indexing. | _None_ |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | `0.20` |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV · Å⁻²). | `300` |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | `10000` |
| `--scan-preopt / --no-scan-preopt` | Override the scan preoptimization toggle. | _None_ |
| `--scan-endopt / --no-scan-endopt` | Override the scan end-of-stage optimization toggle. | _None_ |

## YAML configuration

`all` supports layered YAML — base settings via `--config FILE`, with the precedence `defaults < config < CLI`. The effective YAML is forwarded to every invoked subcommand, and each subcommand reads its own sections:

| Subcommand | YAML sections |
|---|---|
| [`path-search`](path-search.md) | `geom`, `calc`, `gs`, `stopt`, `opt`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |
| [`irc`](irc.md) | `geom`, `calc`, `irc` |

```yaml
# Minimal example
calc:
  model: uma-s-1p1            # uma-s-1p1 | uma-m-1p1
  hessian_calc_mode: Analytical   # recommended when VRAM permits
gs:
  max_nodes: 12
  climb: true
dft:
  grid_level: 6
```

Full schema: [YAML Reference](yaml-reference.md).

## Notes

```{tip}
For large active-site models, the single-structure scan workflow tends to produce more reliable reaction barriers than the multi-structure MEP workflow. When multiple full PDB structures are provided, structural differences in regions unrelated to the reaction coordinate can accumulate and overestimate the barrier. The scan workflow avoids this by driving only the relevant coordinates from a single starting structure. The effect becomes more pronounced as the model size grows.
```

- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-full-pdb` option of `path-search` is hidden in this wrapper.
- Extraction radii: passing `0` to `--radius` or `--radius-het2het` is internally clamped to `0.001 Å` by the extractor.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to MEP / `tsopt` / `freq` / `dft`; single-structure runs still require either `--scan-lists/-s` or `--tsopt`.

## See Also

[Installation](installation.md) · [Getting Started](getting-started.md) · [extract](extract.md) · [scan](scan.md) · [path-opt](path-opt.md) · [path-search](path-search.md) · [tsopt](tsopt.md) · [irc](irc.md) · [freq](freq.md) · [dft](dft.md) · [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
