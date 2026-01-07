# `all` subcommand

## Overview
`pdb2reaction all` is the umbrella command that orchestrates **every pipeline stage**: pocket extraction → optional staged UMA scan → recursive MEP search (`path_search`, GSM/DMF) → full-system merging → optional TS optimization + IRC (`tsopt`) → optional vibrational analysis (`freq`) → optional single-point DFT (`dft`). The command accepts multi-structure ensembles, converts single-structure scans into ordered intermediates, and can fall back to a TSOPT-only pocket workflow. All downstream tools share a single CLI surface so you can coordinate long reaction campaigns from one invocation. Format-aware XYZ/TRJ → PDB/GJF conversions across every stage are controlled by the shared `--convert-files {True|False}` flag (enabled by default).

Key modes:
- **End-to-end ensemble** – Supply ≥2 PDBs/GJFs/XYZ files in reaction order plus a substrate definition; the command extracts pockets, runs GSM/DMF MEP search, merges to the parent PDB(s), and optionally runs TSOPT/freq/DFT per reactive segment.
- **Single-structure + staged scan** – Provide one PDB plus one or more `--scan-lists`; UMA scans on the extracted pocket generate intermediates that become MEP endpoints.
- **TSOPT-only pocket refinement** – Provide one input structure, omit `--scan-lists`, and enable `--tsopt True`; the command extracts the pocket (if `-c/--center` is given) and only runs TS optimization + IRC (with optional freq/DFT) on that single system.

## Usage
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [options]
```

### Examples
```bash
# Multi-structure ensemble with explicit ligand charges and post-processing
pdb2reaction all -i reactant.pdb product.pdb -c "GPP,MMT" \
    --ligand-charge "GPP:-3,MMT:-1" --mult 1 --freeze-links True \
    --max-nodes 10 --max-cycles 100 --climb True --opt-mode light \
    --out-dir result_all_${date} --tsopt True --thermo True --dft True

# Single-structure staged scan followed by GSM/DMF + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c "308,309" \
    --scan-lists "[('TYR,285,CA','MMT,309,C10',2.20),('TYR,285,CB','MMT,309,C11',1.80)]" \
    --opt-mode heavy --tsopt True --thermo True --dft True

# TSOPT-only workflow (no path search)
pdb2reaction all -i reactant.pdb -c "GPP,MMT" \
    --ligand-charge "GPP:-3,MMT:-1" --tsopt True --thermo True --dft True
```

## Workflow
1. **Active-site pocket extraction** (if `-c/--center` is provided)
   - Substrates may be specified via PDB paths, residue IDs (`123,124` or `A:123,B:456`), or residue names (`GPP,MMT`).
   - Optional toggles forward to the extractor: `--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected_resn`, and `--verbose`.
   - Per-input pocket PDBs are saved under `<out-dir>/pockets/`. When multiple structures are supplied, their pockets are unioned per residue selection.
   - The **first pocket’s total charge** is rounded to the nearest integer and propagated to scan/MEP/TSOPT (a console note appears when rounding occurs).

2. **Optional staged scan (single-input only)**
   - Each `--scan-lists` argument is a Python-like list of `(i,j,target_Å)` tuples describing a UMA scan stage. Atom indices refer to the original input PDB (1-based) and are remapped to the pocket ordering. For PDB inputs, `i`/`j` can be integer indices or selector strings like `"TYR,285,CA"`; selectors accept spaces/commas/slashes/backticks/backslashes as delimiters and allow unordered tokens (fallback assumes resname, resseq, atom).
   - When `--scan-lists` is repeated, stages run **sequentially**: stage 1 relaxes the starting structure, stage 2 begins from stage 1's result, and so on.
   - Scan inherits charge/spin, `--freeze-links`, the UMA optimizer preset (`--opt-mode`), `--dump`, `--args-yaml`, and `--preopt`. Overrides such as `--scan-out-dir`, `--scan-one-based`, `--scan-max-step-size`, `--scan-bias-k`, `--scan-relax-max-cycles`, `--scan-preopt`, and `--scan-endopt` apply per run.
   - Stage endpoints (`stage_XX/result.pdb`) become the ordered intermediates that feed the subsequent MEP step.

3. **MEP search on pockets (recursive GSM/DMF)**
   - Executes `path_search` by default using the extracted pockets (or the original structures if extraction is skipped). Relevant options: `--mult`, `--freeze-links`, `--max-nodes`, `--max-cycles`, `--climb`, `--opt-mode`, `--dump`, `--preopt`, `--args-yaml`, and `--out-dir`.
   - Use `--refine-path False` to switch to a single-pass `path-opt` GSM/DMF chain without the recursive refiner.
   - For multi-input PDB runs, the full-system templates are automatically passed to `path_search` for reference merging. Single-structure scan runs reuse the original full PDB template for every stage.

4. **Merge pockets back to the full systems**
   - When reference PDB templates exist, merged `mep_w_ref*.pdb` and per-segment `mep_w_ref_seg_XX.pdb` files are emitted under `<out-dir>/path_search/`.

5. **Optional per-segment post-processing**
   - `--tsopt True`: run TS optimization on each HEI pocket, follow with EulerPC IRC, and emit segment energy diagrams.
   - `--thermo True`: call `freq` on (R, TS, P) to obtain vibrational/thermochemistry data and a UMA Gibbs diagram.
   - `--dft True`: launch single-point DFT on (R, TS, P) and build a DFT diagram. When combined with `--thermo True`, a DFT//UMA Gibbs diagram (DFT energies + UMA thermal correction) is also produced.
   - Shared overrides include `--opt-mode`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, and `--dft-engine` (GPU-first by default).

6. **TSOPT-only mode** (single input, `--tsopt True`, no `--scan-lists`)
   - Skips the MEP/merge stages. Runs `tsopt` on the pocket (or full input if extraction is skipped), performs EulerPC IRC, identifies the higher-energy endpoint as reactant (R), and generates the same set of energy diagrams plus optional freq/DFT outputs.

### Charge and spin precedence
- With extraction: pocket charge = rounded extractor charge; spin comes from `--mult` (default 1).
- Without extraction: explicit `-q/--charge` wins. If omitted but `--ligand-charge` is provided, the **full complex is treated as an enzyme–substrate adduct** and `extract.py`’s charge summary logic derives the total charge from the supplied substrate charge(s); otherwise the first `.gjf` supplies the charge or the default is 0. Spin precedence becomes explicit `--mult`, else `.gjf`, else 1.

### Input expectations
- Extraction enabled (`-c/--center`): inputs must be **PDB** files so residues can be located.
- Extraction skipped: inputs may be **PDB/XYZ/GJF**; no staged scan is available unless the input is PDB.
- Multi-structure runs require ≥2 structures unless TSOPT-only mode is triggered as described above.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full structures in reaction order (single input allowed only with `--scan-lists` or `--tsopt True`). | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs like `123,124` / `A:123,B:456`, or residue names like `GPP,MMT`). | Required for extraction |
| `--out-dir PATH` | Top-level output directory. | `./result_all/` |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero–hetero cutoff (Å). | `0.0` |
| `--include-H2O BOOLEAN` | Include waters (set `False` to drop HOH/WAT/TIP3/SOL). | `True` |
| `--exclude-backbone BOOLEAN` | Remove backbone atoms on non-substrate amino acids. | `True` |
| `--add-linkH BOOLEAN` | Add link hydrogens for severed bonds (carbon-only). | `True` |
| `--selected_resn TEXT` | Residues to force include (comma/space separated; chain/insertion codes allowed). | `""` |
| `--ligand-charge TEXT` | Total charge or residue-specific mapping for unknown residues (recommended). When `-q` is omitted, triggers extract-style charge derivation on the full complex. | `None` |
| `-q, --charge INT` | Force the total system charge, overriding extractor rounding / `.gjf` metadata / `--ligand-charge` (logs a warning). | _None_ |
| `--workers`, `--workers-per-node` | UMA predictor parallelism (workers > 1 disables analytic Hessians; `workers_per_node` forwarded to the parallel predictor). | `1`, `1` |
| `--verbose BOOLEAN` | Enable INFO-level extractor logging. | `True` |
| `-m, --mult INT` | Spin multiplicity forwarded to all downstream steps. | `1` |
| `--freeze-links BOOLEAN` | Freeze link parents in pocket PDBs (reused by scan/tsopt/freq). | `True` |
| `--max-nodes INT` | MEP internal nodes per segment (GSM string images or DMF images). | `10` |
| `--max-cycles INT` | MEP maximum optimization cycles (GSM/DMF). | `300` |
| `--climb BOOLEAN` | Enable TS climbing for the first segment in each pair. | `True` |
| `--opt-mode [light\|heavy]` | Optimizer preset shared across scan, tsopt, and path_search (light → LBFGS/Dimer, heavy → RFO/RSIRFO). | `light` |
| `--dump BOOLEAN` | Dump MEP (GSM/DMF) and single-structure trajectories (propagates to scan/tsopt/freq). | `False` |
| `--convert-files {True|False}` | Global toggle for XYZ/TRJ → PDB/GJF companions when templates are available. | `True` |
| `--args-yaml FILE` | YAML forwarded unchanged to `path_search`, `scan`, `tsopt`, `freq`, and `dft`. | _None_ |
| `--preopt BOOLEAN` | Pre-optimize pocket endpoints before MEP search (also the default for scan preopt). | `True` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | Shared UMA Hessian engine forwarded to tsopt and freq. | _None_ (uses YAML/default of `FiniteDifference`) |
| `--tsopt BOOLEAN` | Run TS optimization + IRC per reactive segment, or enable TSOPT-only mode (single input). | `False` |
| `--thermo BOOLEAN` | Run vibrational analysis (`freq`) on R/TS/P and build UMA Gibbs diagram. | `False` |
| `--dft BOOLEAN` | Run single-point DFT on R/TS/P and build DFT energy + optional DFT//UMA diagrams. | `False` |
| `--dft-engine [gpu\|cpu\|auto]` | Preferred backend for the DFT stage (`auto` tries GPU then CPU). | `gpu` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles` for each refinement. | _None_ |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory (resolved against `<out-dir>` when relative). | _None_ |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Override `freq --max-write`. | _None_ |
| `--freq-amplitude-ang FLOAT` | Override `freq --amplitude-ang` (Å). | _None_ |
| `--freq-n-frames INT` | Override `freq --n-frames`. | _None_ |
| `--freq-sort [value\|abs]` | Override freq mode sorting behavior. | _None_ |
| `--freq-temperature FLOAT` | Override freq thermochemistry temperature (K). | _None_ |
| `--freq-pressure FLOAT` | Override freq thermochemistry pressure (atm). | _None_ |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Override `dft --func-basis`. | _None_ |
| `--dft-max-cycle INT` | Override `dft --max-cycle`. | _None_ |
| `--dft-conv-tol FLOAT` | Override `dft --conv-tol`. | _None_ |
| `--dft-grid-level INT` | Override `dft --grid-level`. | _None_ |
| `--scan-lists TEXT...` | One or more Python-like lists describing staged scans on the extracted pocket (single-input runs only). Each element is `(i,j,target_Å)`; `i`/`j` can be integer indices or PDB atom selectors like `"TYR,285,CA"` and are remapped internally. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory (`<out-dir>/scan`). | _None_ |
| `--scan-one-based BOOLEAN` | Force scan indexing to 1-based (`True`) or 0-based (`False`); `None` keeps the scan default (1-based). | _None_ |
| `--scan-max-step-size FLOAT` | Override scan `--max-step-size` (Å). | _None_ |
| `--scan-bias-k FLOAT` | Override the harmonic bias strength `k` (eV/Å²). | _None_ |
| `--scan-relax-max-cycles INT` | Override scan relaxation max cycles per step. | _None_ |
| `--scan-preopt BOOLEAN` | Override the scan preoptimization toggle (otherwise `--preopt` propagates). | _None_ |
| `--scan-endopt BOOLEAN` | Override the scan end-of-stage optimization toggle. | _None_ |

## Outputs
```
out_dir/ (default: ./result_all/)
├─ summary.log               # formatted summary for quick inspection
├─ summary.yaml              # YAML version summary
├─ pockets/                  # Per-input pocket PDBs when extraction runs
├─ scan/                     # Staged pocket scan results (present when --scan-lists is provided)
├─ path_search/              # MEP results (GSM/DMF): trajectories, merged PDBs, diagrams, summary.yaml, per-segment folders
├─ path_search/tsopt_seg_XX/ # Post-processing outputs (TS optimization, IRC, freq, DFT, diagrams)
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
- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-pdb` option of `path_search` is intentionally hidden in this wrapper.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to the MEP/tsopt/freq/DFT stages; single-structure runs still require either `--scan-lists` or `--tsopt True`.
- `--args-yaml` lets you coordinate all calculators from a single configuration file. YAML values override CLI flags.

## YAML configuration (`--args-yaml`)
The same YAML file is forwarded unchanged to **every** invoked subcommand. Each tool reads the sections described in its own documentation:

- [`path_search`](path_search.md#yaml-configuration-args-yaml): `geom`, `calc`, `gs`, `opt`, `sopt`, `bond`, `search`.
- [`scan`](scan.md#yaml-configuration-args-yaml): `geom`, `calc`, `opt`, `lbfgs`, `rfo`, `bias`, `bond`.
- [`tsopt`](tsopt.md#yaml-configuration-args-yaml): `geom`, `calc`, `opt`, `hessian_dimer`, `rsirfo`.
- [`freq`](freq.md#yaml-configuration-args-yaml): `geom`, `calc`, `freq`, `thermo`.
- [`dft`](dft.md#yaml-configuration-args-yaml): `dft`.

Include whichever sections you need at the YAML root; overlapping names such as `geom`, `calc`, or `opt` are shared across modules, while module-specific blocks (for example `freq` or `dft`) apply only where supported. YAML contents take precedence over CLI values when both are provided.

Example snippet combining shared and module-specific sections:

```yaml
geom:
  coord_type: cart                     # coordinate type: cartesian vs dlc internals
calc:
  model: uma-s-1p1                     # UMA model tag
  hessian_calc_mode: FiniteDifference  # Hessian mode selection
gs:
  max_nodes: 12                        # maximum string nodes
freq:
  max_write: 8                         # maximum modes written
dft:
  grid_level: 6                        # PySCF grid level
```

```yaml
geom:
  coord_type: cart           # coordinate type: cartesian vs dlc internals
  freeze_atoms: []           # 0-based frozen atoms merged with CLI/link detection
calc:
  charge: 0                  # total charge (CLI/template override)
  spin: 1                    # spin multiplicity 2S+1
  model: uma-s-1p1           # UMA model tag
  task_name: omol            # UMA task name
  device: auto               # UMA device selection
  max_neigh: null            # maximum neighbors for graph construction
  radius: null               # cutoff radius for neighbor search
  r_edges: false             # store radial edges
  out_hess_torch: true       # request torch-form Hessian
  freeze_atoms: null         # calculator-level frozen atoms
  hessian_calc_mode: FiniteDifference   # Hessian mode selection
  return_partial_hessian: false         # full Hessian (avoids shape mismatches)
gs:
  max_nodes: 10              # maximum string nodes
  perp_thresh: 0.005         # perpendicular displacement threshold
  reparam_check: rms         # reparametrization check metric
  reparam_every: 1           # reparametrization stride
  reparam_every_full: 1      # full reparametrization stride
  param: equi                # parametrization scheme
  max_micro_cycles: 10       # micro-iteration limit
  reset_dlc: true            # rebuild delocalized coordinates each step
  climb: true                # enable climbing image
  climb_rms: 0.0005          # climbing RMS threshold
  climb_lanczos: true        # Lanczos refinement for climbing
  climb_lanczos_rms: 0.0005  # Lanczos RMS threshold
  climb_fixed: false         # keep climbing image fixed
  scheduler: null            # optional scheduler backend
opt:
  type: string               # optimizer type label
  stop_in_when_full: 100     # early stop threshold when string is full
  align: false               # alignment toggle (kept off)
  scale_step: global         # step scaling mode
  max_cycles: 100            # maximum optimization cycles
  dump: false                # dump trajectory/restart data
  dump_restart: false        # dump restart checkpoints
  reparam_thresh: 0.001      # reparametrization threshold
  coord_diff_thresh: 0.0     # coordinate difference threshold
  out_dir: ./result_path_search/   # output directory
  print_every: 10            # logging stride
sopt:
  lbfgs:
    thresh: gau                # LBFGS convergence preset
    max_cycles: 10000          # iteration limit
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_path_search/   # output directory
    keep_last: 7               # history size for LBFGS buffers
    beta: 1.0                  # initial damping beta
    gamma_mult: false          # multiplicative gamma update toggle
    max_step: 0.3              # maximum step length
    control_step: true         # control step length adaptively
    double_damp: true          # double damping safeguard
    mu_reg: null               # regularization strength
    max_mu_reg_adaptions: 10   # cap on mu adaptations
  rfo:
    thresh: gau                # RFOptimizer convergence preset
    max_cycles: 10000          # iteration cap
    print_every: 100           # logging stride
    min_step_norm: 1.0e-08     # minimum accepted step norm
    assert_min_step: true      # assert when steps stagnate
    rms_force: null            # explicit RMS force target
    rms_force_only: false      # rely only on RMS force convergence
    max_force_only: false      # rely only on max force convergence
    force_only: false          # skip displacement checks
    converge_to_geom_rms_thresh: 0.05   # RMS threshold when targeting geometry
    overachieve_factor: 0.0    # tighten thresholds
    check_eigval_structure: false   # validate Hessian eigenstructure
    line_search: true          # enable line search
    dump: false                # dump trajectory/restart data
    dump_restart: false        # dump restart checkpoints
    prefix: ""                 # filename prefix
    out_dir: ./result_path_search/   # output directory
    trust_radius: 0.3          # trust-region radius
    trust_update: true         # enable trust-region updates
    trust_min: 0.01            # minimum trust radius
    trust_max: 0.3             # maximum trust radius
    max_energy_incr: null      # allowed energy increase per step
    hessian_update: bfgs       # Hessian update scheme
    hessian_init: calc         # Hessian initialization source
    hessian_recalc: 100        # rebuild Hessian every N steps
    hessian_recalc_adapt: 2.0  # adaptive Hessian rebuild factor
    small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
    alpha0: 1.0                # initial micro step
    max_micro_cycles: 25       # micro-iteration limit
    rfo_overlaps: false        # enable RFO overlaps
    gediis: false              # enable GEDIIS
    gdiis: true                # enable GDIIS
    gdiis_thresh: 0.0025       # GDIIS acceptance threshold
    gediis_thresh: 0.01        # GEDIIS acceptance threshold
    gdiis_test_direction: true # test descent direction before DIIS
    adapt_step_func: false     # adaptive step scaling toggle
bond:
  device: cuda                # UMA device for bond analysis
  bond_factor: 1.2            # covalent-radius scaling
  margin_fraction: 0.05       # tolerance margin for comparisons
  delta_fraction: 0.05        # minimum relative change to flag bonds
search:
  max_depth: 10               # recursion depth limit
  stitch_rmsd_thresh: 0.0001  # RMSD threshold for stitching segments
  bridge_rmsd_thresh: 0.0001  # RMSD threshold for bridging nodes
  rmsd_align: true            # legacy alignment flag (ignored)
  max_nodes_segment: 10       # max nodes per segment
  max_nodes_bridge: 5         # max nodes per bridge
  kink_max_nodes: 3           # max nodes for kink optimizations
  max_seq_kink: 2             # max sequential kinks
```
