# `all`

`pdb2reaction all` runs the workflow end-to-end from PDB/mmCIF/XYZ/GJF.
mmCIF and oversized PDB inputs are converted to a safely reindexed internal
PDB; coordinate outputs also include CIF with the original identifiers.

`all` runs in one of three modes, chosen by what you pass:

- **Multi-structure MEP** (`[mode] all (mep)`) — give ≥ 2 structures in reaction order. With `-c`, `all` first extracts active-site models; without it, the full supplied structures are used. It then runs GSM / DMF MEP search and optionally runs TSOPT + IRC / freq / DFT per reactive segment. With `--refine-path` and full-system templates, it also merges the optimized path back into those templates.
- **Single-structure staged scan** (`[mode] all (scan-lists)`) — give one structure plus one or more `--scan-lists/-s` literals, each defining a scan stage; the staged scan produces the ordered intermediates that drive the MEP step.
- **TSOPT-only** — give a single input and set `--tsopt` (no `--scan-lists`). `all` skips the MEP / merge stages, runs `tsopt` + EulerPC IRC on the active-site model (or the full input if extraction is skipped), and identifies the higher-energy endpoint as the reactant.

```{note}
The TSOPT-only reactant/product labels follow an **energy-order presentation convention**, not a chemically established reaction direction: the higher-energy IRC endpoint is presented as the reactant (on an exact energy tie the left endpoint is the reactant). The R/P labels, `reactant_irc`/`product_irc` filenames, barrier and delta are computed under this convention. The machine-readable summary records it explicitly under `endpoint_assignment` (`policy = "higher_energy_endpoint_as_reactant"`, `chemical_direction_known = false`); read that field and do not infer chemical direction from the labels alone.
```

```{important}
Without `--tsopt`, the workflow produces **TS candidates** (highest-energy images from MEP search). Adding `--tsopt` refines them into optimized TS structures validated by an imaginary-frequency check. IRC starts only when the optimizer reports a converged first-order saddle (`n_imag = 1`); otherwise `all` stops before downstream post-processing. Always inspect the results (imaginary-frequency count and IRC endpoint connectivity) before mechanistic interpretation.
```

## Examples

Working examples for GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)) covering both multi-structure MEP and scan-based pipelines: [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples).

Command form:

```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] [-c SUBSTRATE] [-b uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] [options]
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

PDB / CIF / GJF companion files are generated automatically when reference templates are available; CIF is emitted for mmCIF/oversized-PDB bridge inputs. Control with `--convert-files` (on by default).

`pdb2reaction all --help` shows core options; `pdb2reaction all --help-advanced` shows the full option list.

## Workflow

```text
Full system(s) (PDB / mmCIF / XYZ / GJF)
  ├─ (optional) active-site extraction `extract` — requires PDB/mmCIF when `-c` is used
  │   └─ active-site cluster model(s)
  │       ├─ (optional) staged scan `scan` — single-structure workflows
  │       │   └─ ordered intermediates
  │       └─ MEP search `path-opt` (default) or `path-search` (recursive, `--refine-path`)
  │           └─ MEP trajectory `mep_trj.xyz` + energy diagrams
  └─ (optional) TS optimization + IRC `tsopt` → `irc`
      ├─ (optional) thermochemistry `freq`
      └─ (optional) single-point DFT `dft`
```

`all` runs the following stages in order. Stages 0 and 1 are automatic preprocessing; the rest fire based on the flags you pass.

0. **Structure bridge and preflight** (automatic) — mmCIF, oversized/nonstandard PDB, and PDB altloc input are converted once to a safely reindexed internal PDB; altloc is selected coherently per residue. For an ordinary PDB with blank element columns, `all` runs `add-elem-info`. Standalone `fix-altloc` is only needed when you want a cleaned PDB deliverable; standalone commands use the same bridge. Missing element data must still be repaired for an ordinary PDB or supplied as mmCIF `_atom_site.type_symbol`.
1. **Active-site model extraction** (when `-c/--center` is set) — accepts PDB/mmCIF paths, IDs/names, `CHAIN:RESNAME`, and `CHAIN:RESNAME:RESSEQ`. Per-input internal PDBs are saved under `<out-dir>/_work/models/`; bridge inputs also produce CIF companions.
2. **Optional staged scan** (single-input only) — each `--scan-lists/-s` literal is a list of `(i, j, target_Å)` tuples. Atom indices use the original input ordering, 1-based by default (pass `--no-scan-one-based` to interpret them as 0-based), and are remapped to the active-site model ordering. Three-field selectors like `'TYR,285,CA'` are order-flexible; use positional `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` for repeated names or numbering. Stages run sequentially (stage 2 starts from stage 1's result), and the stage endpoints become the ordered intermediates that feed the MEP step.
3. **MEP search** — by default runs single-pass `path-opt`; `--refine-path` switches to recursive `path-search`. Recursive refinement can improve a poor HEI but can also split a noisy/bad path into unnecessary segments and increase cost, so it is off by default. Segmentation is only a candidate mechanism until TS/frequency/IRC validation. Raw engine output stays under `_work`; `mep.pdb`, bridge-input `mep.cif`, `mep_trj.xyz`, and the diagram are promoted to the top level.
4. **Merge to full systems** (with `--refine-path`) — writes `mep_w_ref.pdb` and, for mmCIF/oversized-PDB templates, `mep_w_ref.cif`; per-segment equivalents remain under `_work/path_search/`.
5. **Per-segment post-processing** (reactive segments only — bridge segments without bond changes are skipped):
   - `--tsopt` — TS optimization on each HEI active-site model. A machine-readable exact-Hessian result must report `status=converged` and `n_imag=1`; otherwise the workflow stops before IRC. A validated TS is followed by EulerPC IRC and IRC-endpoint re-optimization with `--thresh-post` (default `baker`). For Hessian TS optimization, the MEP's energy-upwinding Cartesian tangent is passed by default as the reaction reference mode (a normalized secant bisector is used only for legacy trajectories without readable energies). With `--no-tsopt-from-mep-tan`, TSOPT instead computes the initial-structure Hessian and selects the initial mode from its vibrational modes. The endpoint optimization working directory is deleted automatically after completion. Endpoint RFO uphill rejection is disabled by default; pass `--reject-uphill` to enable it for endpoint re-optimization only.
   - `--thermo` — `freq` on (R, TS, P) for vibrational + thermochemistry data and an MLIP Gibbs diagram.
   - `--dft` — single-point DFT on (R, TS, P) and a DFT diagram. With `--thermo`, a DFT//MLIP Gibbs diagram (DFT energies + MLIP thermal correction) is also produced.
   - Shared overrides: `--opt-mode`, `--opt-mode-post`, `--flatten`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, `--dft-engine` (GPU-first by default). Frozen-boundary PHVA always uses the constrained rigid-mode treatment; it is unrelated to the MEP-derived `--ref-mode`. For Hessian evaluation modes see {ref}`hessian-evaluation`.
6. **TSOPT-only mode** (single input + `--tsopt`, no `--scan-lists`) — skips MEP / merge; runs `tsopt` + EulerPC IRC and generates the same energy diagrams plus optional freq / DFT outputs.

## Outputs

The tree has three top-level zones: **deliverables at the root**, **per-segment deliverables under `segments/seg_NN/`**, and **pipeline scratch under `_work/`** (safe to `rm -rf` once you have the results you need).

```text
out_dir/   (default: ./result_all/)
├─ summary.log                 # Text summary (authored at the root)
├─ summary.json                # JSON results
├─ mep.pdb                     # Merged MEP path (promoted from the engine)
├─ mep.cif                     # Bridge inputs only; original identifiers restored
├─ mep_w_ref.pdb               # --refine-path + reference template only
├─ mep_w_ref.cif               # --refine-path + bridge reference template only
├─ mep_trj.xyz                 # Full MEP trajectory
├─ energy_diagram_MEP.png      # All-segment MEP barriers
├─ energy_diagram_*.png        # Aggregated post-processing diagrams (MLIP / Gibbs / DFT, with --tsopt etc.)
├─ segments/                   # Per-reactive-segment deliverables (bridge segments are skipped)
│  └─ seg_NN/                  # 2-digit index, e.g. seg_01, seg_02
│     ├─ reactant.{pdb,cif,xyz,gjf} # Canonical R/TS/P; bridge inputs write PDB + CIF
│     ├─ ts.{pdb,cif,xyz,gjf}
│     ├─ product.{pdb,cif,xyz,gjf}
│     ├─ ts/                   # TS optimization output + vibrational analysis (--tsopt)
│     ├─ irc/                  # IRC trajectories + plots (--tsopt)
│     ├─ freq/{R,TS,P}/        # frequencies_cm-1.txt + thermoanalysis.yaml (--thermo)
│     └─ dft/                  # DFT single-point results (--dft)
└─ _work/                      # Pipeline scratch (safe to delete)
   ├─ models/                  # Extracted active-site model PDBs (model_<input_stem>.pdb, when extraction runs)
   ├─ scan/                    # Staged scan results (with --scan-lists)
   ├─ add_elem_info/           # Preflight element-symbol fills
   └─ path_opt/                # Raw MEP-engine output (path_search/ with --refine-path)
```

In **TSOPT-only mode** (single input + `--tsopt`, no `--scan-lists`) there is no MEP stage: the optimized R/TS/P plus `ts/`, `irc/`, `freq/`, and `dft/` land directly under `segments/seg_01/`, and the MEP work directory (`_work/path_opt/`) is absent.

```{note}
**The canonical structures are `segments/seg_NN/reactant.*`, `ts.*`, `product.*`** — cite these when reporting mechanisms. The `ts/`, `irc/`, `freq/`, and `dft/` subdirectories inside the same `seg_NN/` hold the per-stage working files (e.g. `ts/vib/imag_*_trj.xyz`, `irc/*_trj.xyz`) for debugging a single stage. The raw MEP-search engine output under `_work/path_opt/` is scratch — the products you need (`mep.pdb`, bridge-input `mep.cif`, `mep_trj.xyz`, `energy_diagram_MEP.png`) are already promoted to the root.
```

At `-v 2` the console summarizes active-site charge resolution, YAML contents, scan stages, MEP progress (GSM / DMF), and per-stage timing; see {ref}`verbosity-levels`.

### Plot file naming

Energy-diagram filenames encode method and scope:

| File | Generated when | Content |
|---|---|---|
| `energy_diagram_MEP.png` | `path-opt` / `path-search` completes | All-segment MEP barriers (raw GSM / DMF values) |
| `energy_diagram_MLIP.png` | per-segment `tsopt` + IRC completes | R → TS → P (MLIP energy) |
| `energy_diagram_G_MLIP.png` | per-segment thermo completes | R → TS → P (MLIP Gibbs) |
| `energy_diagram_DFT.png` | per-segment DFT completes | R → TS → P (DFT energy) |
| `energy_diagram_G_DFT_plus_MLIP.png` | per-segment DFT + thermo | R → TS → P (DFT energy + MLIP thermal correction) |
| `energy_diagram_MLIP_all.png` / `_G_MLIP_all.png` / `_DFT_all.png` / `_G_DFT_plus_MLIP_all.png` | all segments aggregated (variants for MLIP / Gibbs / DFT / DFT//MLIP Gibbs) | Combined across all segments |
| `irc_plot.png` (per `segments/seg_NN/irc/`) | per-segment IRC completes | IRC profile (MLIP energy along the trajectory) |
| `irc_plot_all.png` | all segments aggregated | IRC profiles concatenated across segments |

### Reading `summary.log`

The log is organized into numbered sections:

- **[1] Global MEP overview** — image / segment counts, MEP trajectory plot paths, aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (MLIP path)** — per-segment barriers (ΔE‡), reaction energies (ΔE), bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** — TS imaginary-frequency checks, IRC outputs, MLIP / thermo / DFT energy tables.
- **[4] Energy diagrams (overview)** — diagram tables for MEP / MLIP / Gibbs / DFT plus an optional cross-method summary.
- **[5] Output directory structure** — a compact tree of generated files with inline annotations.

### Reading `summary.json`

Top-level keys: `out_dir`, `n_images`, `n_segments` (run metadata and counts); `segments` (per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes`); `energy_diagrams` (optional payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` paths). `summary.json` intentionally omits the formatted tables and filesystem tree from `summary.log`.

## CLI options

Defaults shown are used when the option is not specified. The full flag list is in the generated [command reference](reference/commands/index.md); the tables below cover the options that need explanation.

Input expectations:

- Extraction enabled (`-c/--center`): inputs must be **PDB / mmCIF** so residues can be located.
- Extraction skipped: inputs may be **PDB / mmCIF / XYZ / GJF**.
- Multi-structure runs require ≥ 2 structures. For full input-file requirements (hydrogens, element columns, atom-order parity), see [CLI Conventions](cli-conventions.md).

Charge is resolved via the standard priority chain (see {ref}`CLI Conventions: Charge specification <charge-specification>`). In `all`, the charge derivation from active-site model extraction (when `-c` is set) acts as an additional priority layer. Spin resolution: `--multiplicity` CLI → `.gjf` template → default `1`. Always provide `--ligand-charge/-l` for non-standard substrates so the correct net charge propagates to scan / MEP / TSOPT / DFT.

### Input / output

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full structures in reaction order (single input allowed only with `--scan-lists/-s` or `--tsopt`). | Required |
| `--ref-pdb FILE` | Reference PDB/mmCIF topology when `-i` provides XYZ/GJF coordinates. | _None_ |
| `-o, --out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--convert-files / --no-convert-files` | Global toggle for XYZ / TRJ → PDB / CIF / GJF companions when templates are available. CIF requires a bridged mmCIF/oversized-PDB topology. | `True` |
| `--dump / --no-dump` | Dump MEP (GSM / DMF) trajectories. An explicit parent toggle is forwarded to `path-search` / `path-opt` and `scan` / `tsopt`; when omitted, each child resolves its YAML/default. With `--thermo`, freq always retains the internal `thermoanalysis.yaml` channel even when `--no-dump` is supplied. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config / --no-show-config` | Print the effective YAML plus parent settings before execution. Child-governing parent options are `null` when omitted, because each child then resolves its own YAML/default. | `False` |
| `--dry-run / --no-dry-run` | Validate options and print the plan. With `--center`, run extraction in a temporary directory to validate derived charge and electron parity. No scan/MEP/TSOPT/freq/DFT stage runs, and no persistent output is produced. | `False` |

### Charge / spin

| Option | Description | Default |
| --- | --- | --- |
| `-l, --ligand-charge TEXT` | Net charge or per-resname mapping used when `-q` is omitted (PDB/mmCIF metadata or `--ref-pdb`). | _None_ |
| `-q, --charge INT` | Explicit net system charge with highest priority. If it differs from the extracted/workflow-derived value, warn and use `-q`; omit it for automatic derivation. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity forwarded to all downstream steps. | `1` |

### Extraction

| Option | Description | Default |
| --- | --- | --- |
| `-c, --center TEXT` | PDB/mmCIF path, IDs/names, `CHAIN:RESNAME`, or `CHAIN:RESNAME:RESSEQ`. | Required for extraction |
| `-r, --radius FLOAT` | Active-site model inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). `0` is internally nudged to `0.001 Å` to avoid empty selections (same as standalone `extract`). | `0.0` |
| `--include-h2o / --no-include-h2o` | Include waters (HOH / WAT / TIP3 / SOL). | `True` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `False` |
| `--add-linkh / --no-add-linkh` | Add cap hydrogens for severed bonds. | `True` |
| `--selected-resn TEXT` | Force-include using the same ID/name/chain-qualified selector forms as `--center`. | `""` |
| `--modified-residue TEXT` | Comma-separated residue names (with optional charge) to treat as amino acids for backbone truncation and charge assignment (e.g. `HD1,HD2,HD3` or `HD1:0,SEP:-2`). | `""` |
| `--freeze-links / --no-freeze-links` | Freeze cap parents in active-site model PDBs. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices frozen in every stage; indices refer to the extracted model. Merged with `--freeze-links` and YAML `geom.freeze_atoms`. | _None_ |

### MEP search

```{note}
Do not set `--max-cycles` on `pdb2reaction all`; let each stage use its own
default. Set `--max-cycles` only when running a single-stage subcommand
directly, such as `opt`, `tsopt`, or `path-opt`.
```

| Option | Description | Default |
| --- | --- | --- |
| `--mep-mode [gsm\|dmf]` | MEP algorithm: GSM (Growing String Method) or DMF (Direct Max Flux). | `gsm` |
| `--max-nodes INT` | Movable internal images per GSM/DMF segment. Both engines retain two endpoints, so total images = `max_nodes + 2`. | `20` |
| `--max-cycles INT` | MEP maximum optimization cycles. | `300` |
| `--climb / --no-climb` | Enable climbing image for standard GSM segments (bridge segments always disable climbing). | `True` |
| `--opt-mode [grad\|hess]` | Workflow preset (`grad` → L-BFGS / Dimer, `hess` → RFO / RS-P-RFO). Token-to-algorithm mapping depends on scope — see {ref}`opt-mode-semantics` for the per-subcommand table; note that `all`'s pre-opt default (`grad`) differs from `tsopt`'s default (`hess`). | `grad` |
| `--thresh TEXT` | Convergence preset for single-structure optimizations and scan relaxations (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--thresh-gsm TEXT` | Convergence preset for the GSM string optimizer of the MEP stage (same presets as `--thresh`). | `gau_loose` |
| `--thresh-dmf TEXT` | IPOPT dual-infeasibility tolerance of the DMF MEP stage: `tight` (0.04), `middle` (0.10), `loose` (0.20), or a positive float. Not a Gaussian preset. | `tight` |
| `--preopt / --no-preopt` | Pre-optimize active-site model endpoints before MEP search. Standalone `scan` / `scan2d` / `scan3d` default `--preopt` to `False`. | `True` |
| `--refine-path / --no-refine-path` | Enable recursive `path-search` with automatic bond-change segmentation / use the default single-pass `path-opt` per adjacent pair. | disabled |

### MLIP calculator

| Option | Description | Default |
| --- | --- | --- |
| `--workers`, `--workers-per-node` | UMA predictor parallelism. `workers > 1` cannot be combined with an explicit analytical Hessian request; use `workers = 1` or finite differences. See {ref}`workers-analytical-error`. | `1`, `1` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | Shared MLIP Hessian engine. | `FiniteDifference` |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |

### Post-processing

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt / --no-tsopt` | Run TS optimization + IRC per reactive segment. | `False` |
| `--tsopt-from-mep-tan / --no-tsopt-from-mep-tan` | Select the initial TS root from the HEI MEP tangent; when off, select it from the initial-structure Hessian modes. | `True` |
| `--thermo / --no-thermo` | Run vibrational analysis (`freq`) on R / TS / P; requires `--tsopt`. | `False` |
| `--dft / --no-dft` | Run single-point DFT on R / TS / P; requires `--tsopt`. | `False` |
| `--opt-mode-post [grad\|hess]` | Optimizer preset for TSOPT + post-IRC (`grad` → Dimer / L-BFGS, `hess` → RS-P-RFO / RFO). | `hess` |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations. | `baker` |
| `--flatten / --no-flatten` | Enable surplus-imaginary-mode flattening in `tsopt`. | `False` |
| `--reject-uphill / --no-reject-uphill` | Opt in to rejecting energy-raising RFO steps during post-IRC **endpoint re-optimization only**, using a `1e-4` Hartree tolerance (roll back to the lower-energy geometry and shrink the trust radius); TS optimization forces rejection off, and path search is unaffected. At the emergency floor, the retained endpoint receives a final normal convergence check. | `False` |
| `--irc-step-size FLOAT` | Override the IRC maximum EulerPC step (Bohr). If IRC stops after only a few frames, retry with a smaller value such as `0.05`. | IRC default `0.10` |
| `--irc-never-stop / --no-irc-never-stop` | Ignore IRC gradient and energy stop conditions and trace each branch until its max-cycle cap. Numerical/integration failures and external interruption still stop propagation. | `False` |

```{warning}
`--dft` cost and memory depend on basis-function count, elements, functional,
grid, engine, and hardware; atom count alone is not a reliable cutoff. Pilot a
representative state and monitor peak memory before selecting resources.
```

TSOPT optimizer selection order: `--opt-mode-post` (if set) → `--opt-mode` (only when explicitly provided) → TSOPT default (`hess` → `rsprfo`). Example: `--opt-mode grad --opt-mode-post hess` uses L-BFGS for path optimization and RS-P-RFO for TS refinement.

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
| `--scan-one-based / --no-scan-one-based` | How to read the `--scan-lists` atom indices: `True` = 1-based, `False` = 0-based. | _None_ (1-based) |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | `0.20` |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV · Å⁻²). | `300` |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | `10000` |
| `--scan-preopt / --no-scan-preopt` | Override the scan preoptimization toggle. | _None_ |
| `--scan-endopt / --no-scan-endopt` | Override the scan end-of-stage optimization toggle. | _None_ |

## YAML configuration

`all` supports layered YAML — base settings via `--config FILE`, with the precedence `defaults < config < CLI`. The effective YAML is forwarded to every invoked subcommand, and each subcommand reads its own sections:

| Subcommand | YAML sections |
|---|---|
| [`path-opt`](path-opt.md) | `geom`, `calc`, `gs`, `dmf`, `stopt`, `opt`, `lbfgs`, `rfo` |
| [`path-search`](path-search.md) | `geom`, `calc`, `gs`, `dmf`, `stopt`, `opt`, `lbfgs`, `rfo`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond` |
| [`tsopt`](tsopt.md) | `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |
| [`irc`](irc.md) | `geom`, `calc`, `irc` |

```yaml
# Minimal example
calc:
  model: uma-s-1p2            # uma-s-1p2 | uma-m-1p1
  hessian_calc_mode: FiniteDifference   # default; benchmark Analytical before opting in
gs:
  max_nodes: 12
  climb: true
dft:
  grid_level: 6
```

Full schema: [YAML Reference](yaml-reference.md).

## Notes

Structural differences outside the reaction coordinate can affect barriers
obtained from independently prepared full-system structures. Inspect the
structures and validate the selected path workflow for the modeled system.

- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-full-pdb` option of `path-search` is hidden in this wrapper.
- Extraction radii: passing `0` to `--radius` or `--radius-het2het` is internally clamped to `0.001 Å` by the extractor.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to MEP / `tsopt` / `freq` / `dft`; single-structure runs still require either `--scan-lists/-s` or `--tsopt`.

## See Also

[Installation](installation.md) · [Getting Started](getting-started.md) · [extract](extract.md) · [scan](scan.md) · [path-opt](path-opt.md) · [path-search](path-search.md) · [tsopt](tsopt.md) · [irc](irc.md) · [freq](freq.md) · [dft](dft.md) · [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
