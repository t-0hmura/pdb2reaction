# `all` subcommand

## Overview
`pdb2reaction all` is the umbrella command that orchestrates **every pipeline stage**: pocket extraction → optional staged UMA scan → recursive GSM (`path_search`) → full-system merging → optional TS optimization + IRC (`tsopt`) → optional vibrational analysis (`freq`) → optional single-point DFT (`dft`). The command accepts multi-structure ensembles, converts single-structure scans into ordered intermediates, and can fall back to a TSOPT-only pocket workflow. All downstream tools share a single CLI surface so you can coordinate long reaction campaigns from one invocation. Format-aware XYZ/TRJ → PDB/GJF conversions across every stage are controlled by the shared `--convert-files/--no-convert-files` flag (enabled by default).

Key modes:
- **End-to-end ensemble** – Supply ≥2 PDBs/GJFs/XYZ files in reaction order plus a substrate definition; the command extracts pockets, runs GSM, merges to the parent PDB(s), and optionally runs TSOPT/freq/DFT per reactive segment.
- **Single-structure + staged scan** – Provide one PDB plus one or more `--scan-lists`; UMA scans on the extracted pocket generate intermediates that become GSM endpoints.
- **TSOPT-only pocket refinement** – Provide one input structure, omit `--scan-lists`, and enable `--tsopt True`; the command extracts the pocket (if `-c/--center` is given) and only runs TS optimization + pseudo-IRC (with optional freq/DFT) on that single system.

## Usage
```bash
pdb2reaction all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [options]
```

### Examples
```bash
# Multi-structure ensemble with explicit ligand charges and post-processing
date=$(date +%Y%m%d)
pdb2reaction all -i reactant.pdb product.pdb -c "GPP,MMT" \
    --ligand-charge "GPP:-3,MMT:-1" --multiplicity 1 --freeze-links True \
    --max-nodes 10 --max-cycles 100 --climb True --opt-mode light \
    --out-dir result_all_${date} --tsopt True --thermo True --dft True

# Single-structure staged scan followed by GSM + TSOPT/freq/DFT
pdb2reaction all -i single.pdb -c "308,309" \
    --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
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
   - The **first pocket’s total charge** is rounded to the nearest integer and propagated to scan/GSM/TSOPT (a console note appears when rounding occurs).

2. **Optional staged scan (single-input only)**
   - Each `--scan-lists` argument is a Python-like list of `(i,j,target_Å)` tuples describing a UMA scan stage. Atom indices refer to the original input PDB (1-based) and are remapped to the pocket ordering.
   - Scan inherits charge/spin, `--freeze-links`, the UMA optimizer preset (`--opt-mode`), `--dump`, `--args-yaml`, and `--preopt`. Overrides such as `--scan-out-dir`, `--scan-one-based`, `--scan-max-step-size`, `--scan-bias-k`, `--scan-relax-max-cycles`, `--scan-preopt`, and `--scan-endopt` apply per run.
   - Stage endpoints (`stage_XX/result.pdb`) become the ordered intermediates that feed the subsequent GSM step.

3. **MEP search on pockets (recursive GSM)**
   - Executes `path_search` by default using the extracted pockets (or the original structures if extraction is skipped). Relevant options: `--multiplicity`, `--freeze-links`, `--max-nodes`, `--max-cycles`, `--climb`, `--opt-mode`, `--dump`, `--preopt`, `--args-yaml`, and `--out-dir`.
   - Use `--refine-path False` to switch to a single-pass `path-opt` GSM chain without the recursive refiner.
   - For multi-input PDB runs, the full-system templates are automatically passed to `path_search` for reference merging. Single-structure scan runs reuse the original full PDB template for every stage.

4. **Merge pockets back to the full systems**
   - When reference PDB templates exist, merged `mep_w_ref*.pdb` and per-segment `mep_w_ref_seg_XX.pdb` files are emitted under `<out-dir>/path_search/`.

5. **Optional per-segment post-processing**
   - `--tsopt True`: run TS optimization on each HEI pocket, follow with EulerPC pseudo-IRC, and emit segment energy diagrams.
   - `--thermo True`: call `freq` on (R, TS, P) to obtain vibrational/thermochemistry data and a UMA Gibbs diagram.
   - `--dft True`: launch single-point DFT on (R, TS, P) and build a DFT diagram. When combined with `--thermo True`, a DFT//UMA Gibbs diagram (DFT energies + UMA thermal correction) is also produced.
   - Shared overrides include `--opt-mode`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`, and `--dft-engine` (GPU-first by default).

6. **TSOPT-only mode** (single input, `--tsopt True`, no `--scan-lists`)
   - Skips the GSM/merge stages. Runs `tsopt` on the pocket (or full input if extraction is skipped), performs EulerPC IRC, identifies the higher-energy endpoint as reactant (R), and generates the same set of energy diagrams plus optional freq/DFT outputs.

### Charge and spin precedence
- With extraction: pocket charge = rounded extractor charge; spin comes from `--multiplicity` (default 1).
- Without extraction: total system charge follows (1) numeric `--ligand-charge`, else (2) parsed from the first `.gjf`, else defaults to 0. Spin precedence becomes explicit `--multiplicity`, else `.gjf`, else 1.

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
| `--ligand-charge TEXT` | Total charge or residue-specific mapping for unknown residues (recommended). | `None` |
| `-q, --charge INT` | Force the total system charge, overriding extractor rounding / `.gjf` metadata / `--ligand-charge` (logs a warning). | _None_ |
| `--verbose BOOLEAN` | Enable INFO-level extractor logging. | `True` |
| `-m, --multiplicity INT` | Spin multiplicity forwarded to all downstream steps. | `1` |
| `--freeze-links BOOLEAN` | Freeze link parents in pocket PDBs (reused by scan/tsopt/freq). | `True` |
| `--max-nodes INT` | GSM internal nodes per segment. | `10` |
| `--max-cycles INT` | GSM maximum optimization cycles. | `300` |
| `--climb BOOLEAN` | Enable TS climbing for the first segment in each pair. | `True` |
| `--opt-mode [light\|heavy]` | Optimizer preset shared across scan, tsopt, and path_search (light → LBFGS/Dimer, heavy → RFO/RSIRFO). | `light` |
| `--dump BOOLEAN` | Dump GSM and single-structure trajectories (propagates to scan/tsopt/freq). | `False` |
| `--convert-files/--no-convert-files` | Global toggle for XYZ/TRJ → PDB/GJF companions when templates are available. | `--convert-files` |
| `--args-yaml FILE` | YAML forwarded unchanged to `path_search`, `scan`, `tsopt`, `freq`, and `dft`. | _None_ |
| `--preopt BOOLEAN` | Pre-optimise pocket endpoints before GSM (also the default for scan preopt). | `True` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | Shared UMA Hessian engine forwarded to tsopt and freq. | _None_ (uses YAML/default of `FiniteDifference`) |
| `--tsopt BOOLEAN` | Run TS optimisation + pseudo-IRC per reactive segment, or enable TSOPT-only mode (single input). | `False` |
| `--thermo BOOLEAN` | Run vibrational analysis (`freq`) on R/TS/P and build UMA Gibbs diagram. | `False` |
| `--dft BOOLEAN` | Run single-point DFT on R/TS/P and build DFT energy + optional DFT//UMA diagrams. | `False` |
| `--dft-engine [gpu\|cpu\|auto]` | Preferred backend for the DFT stage (`auto` tries GPU then CPU). | `gpu` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles` for each refinement. | _None_ |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory (resolved against `<out-dir>` when relative). | _None_ |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Override `freq --max-write`. | _None_ |
| `--freq-amplitude-ang FLOAT` | Override `freq --amplitude-ang` (Å). | _None_ |
| `--freq-n-frames INT` | Override `freq --n-frames`. | _None_ |
| `--freq-sort [value\|abs]` | Override freq mode sorting behaviour. | _None_ |
| `--freq-temperature FLOAT` | Override freq thermochemistry temperature (K). | _None_ |
| `--freq-pressure FLOAT` | Override freq thermochemistry pressure (atm). | _None_ |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Override `dft --func-basis`. | _None_ |
| `--dft-max-cycle INT` | Override `dft --max-cycle`. | _None_ |
| `--dft-conv-tol FLOAT` | Override `dft --conv-tol`. | _None_ |
| `--dft-grid-level INT` | Override `dft --grid-level`. | _None_ |
| `--scan-lists TEXT...` | One or more Python-like lists describing staged scans on the extracted pocket (single-input runs only). Each element is `(i,j,target_Å)`; indices come from the original input PDB (1-based) and are remapped internally. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory (`<out-dir>/scan`). | _None_ |
| `--scan-one-based BOOLEAN` | Force scan indexing to 1-based (`True`) or 0-based (`False`); `None` keeps the scan default (1-based). | _None_ |
| `--scan-max-step-size FLOAT` | Override scan `--max-step-size` (Å). | _None_ |
| `--scan-bias-k FLOAT` | Override the harmonic bias strength `k` (eV/Å²). | _None_ |
| `--scan-relax-max-cycles INT` | Override scan relaxation max cycles per step. | _None_ |
| `--scan-preopt BOOLEAN` | Override the scan preoptimisation toggle (otherwise `--preopt` propagates). | _None_ |
| `--scan-endopt BOOLEAN` | Override the scan end-of-stage optimisation toggle. | _None_ |

## Outputs
- `<out-dir>/summary.log`: Human-readable run digest (CLI invocation, MEP/segment stats, post-processing energies, key files);
  also stored per GSM/TSOPT branch in `<out-dir>/path_search/*/summary.log`.
- `<out-dir>/pockets/`: Per-input pocket PDBs when extraction runs.
- `<out-dir>/scan/`: Present when `--scan-lists` is used; contains staged pocket scan results (`stage_XX/result.pdb`).
- `<out-dir>/path_search/`: GSM results (trajectory, merged full-system PDBs, energy diagrams, `summary.yaml`, per-segment folders).
- `<out-dir>/path_search/tsopt_seg_XX/`: Present when post-processing is enabled; includes TS optimisation, pseudo-IRC, freq, and DFT outputs plus diagrams.
- `<out-dir>/tsopt_single/`: Present only in TSOPT-only mode; contains TS optimisation outputs, IRC endpoints, and optional freq/DFT directories.
- Console logs summarising pocket charge resolution, YAML contents, scan stages, GSM progress, and per-stage timing.

## Notes
- Always provide `--ligand-charge` (numeric or per-residue mapping) when formal charges cannot be inferred so the correct total charge propagates to scan/GSM/TSOPT/DFT.
- Reference PDB templates for merging are derived automatically from the original inputs; the explicit `--ref-pdb` option of `path_search` is intentionally hidden in this wrapper.
- Energies in diagrams are reported relative to the first state (reactant) in kcal/mol.
- Omitting `-c/--center` skips extraction and feeds the entire input structures directly to GSM/tsopt/freq/DFT; single-structure runs still require either `--scan-lists` or `--tsopt True`.
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
  coord_type: cart
calc:
  model: uma-s-1p1
  hessian_calc_mode: FiniteDifference
gs:
  max_nodes: 12
freq:
  max_write: 8
dft:
  grid_level: 6
```

```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
  model: uma-s-1p1
  task_name: omol
  device: auto
  max_neigh: null
  radius: null
  r_edges: false
  out_hess_torch: true
  freeze_atoms: null
  hessian_calc_mode: FiniteDifference
  return_partial_hessian: true
gs:
  max_nodes: 10
  perp_thresh: 0.005
  reparam_check: rms
  reparam_every: 1
  reparam_every_full: 1
  param: equi
  max_micro_cycles: 10
  reset_dlc: true
  climb: true
  climb_rms: 0.0005
  climb_lanczos: true
  climb_lanczos_rms: 0.0005
  climb_fixed: false
  scheduler: null
opt:
  type: string
  stop_in_when_full: 100
  align: false
  scale_step: global
  max_cycles: 100
  dump: false
  dump_restart: false
  reparam_thresh: 0.001
  coord_diff_thresh: 0.0
  out_dir: ./result_path_search/
  print_every: 10
sopt:
  lbfgs:
    thresh: gau
    max_cycles: 10000
    print_every: 100
    min_step_norm: 1.0e-08
    assert_min_step: true
    rms_force: null
    rms_force_only: false
    max_force_only: false
    force_only: false
    converge_to_geom_rms_thresh: 0.05
    overachieve_factor: 0.0
    check_eigval_structure: false
    line_search: true
    dump: false
    dump_restart: false
    prefix: ""
    out_dir: ./result_path_search/
    keep_last: 7
    beta: 1.0
    gamma_mult: false
    max_step: 0.3
    control_step: true
    double_damp: true
    mu_reg: null
    max_mu_reg_adaptions: 10
  rfo:
    thresh: gau
    max_cycles: 10000
    print_every: 100
    min_step_norm: 1.0e-08
    assert_min_step: true
    rms_force: null
    rms_force_only: false
    max_force_only: false
    force_only: false
    converge_to_geom_rms_thresh: 0.05
    overachieve_factor: 0.0
    check_eigval_structure: false
    line_search: true
    dump: false
    dump_restart: false
    prefix: ""
    out_dir: ./result_path_search/
    trust_radius: 0.3
    trust_update: true
    trust_min: 0.01
    trust_max: 0.3
    max_energy_incr: null
    hessian_update: bfgs
    hessian_init: calc
    hessian_recalc: 100
    hessian_recalc_adapt: 2.0
    small_eigval_thresh: 1.0e-08
    alpha0: 1.0
    max_micro_cycles: 25
    rfo_overlaps: false
    gediis: false
    gdiis: true
    gdiis_thresh: 0.0025
    gediis_thresh: 0.01
    gdiis_test_direction: true
    adapt_step_func: false
bond:
  device: cuda
  bond_factor: 1.2
  margin_fraction: 0.05
  delta_fraction: 0.05
search:
  max_depth: 10
  stitch_rmsd_thresh: 0.0001
  bridge_rmsd_thresh: 0.0001
  rmsd_align: true
  max_nodes_segment: 10
  max_nodes_bridge: 5
  kink_max_nodes: 3
  max_seq_kink: 2
```
