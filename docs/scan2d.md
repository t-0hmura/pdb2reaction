# `scan2d`

## Overview

> **Summary:** Perform a two-distance (d₁, d₂) grid scan with harmonic restraints and MLIP relaxations. Use `--scan-lists/-s` with a YAML/JSON spec file (recommended) or an inline Python literal.

### At a glance
- **Use when:** You want a 2D potential-energy map over two distances `(d₁, d₂)` — e.g. to locate a TS region or visualize the reaction landscape before MEP refinement. Input is one structure + `-s/--scan-lists scan2d.yaml` (recommended), or a single `--scan-lists/-s` inline literal containing exactly two quadruples.
- **Method:** Linear grids built with `--max-step-size`; each axis is reordered so the point closest to the (pre)optimized structure is visited first. Each grid point is relaxed with the appropriate harmonic restraints active (MLIP backend, UMA by default). Values written to `surface.csv` are always evaluated **without bias**, so grid points are directly comparable.
- **Outputs:** `surface.csv` plus `scan2d_map.png` (2D contour) and `scan2d_landscape.html` (3D surface), and per-point structures under `grid/`.
- **Defaults:** `--opt-mode grad` (LBFGS), `--no-preopt`, `--max-step-size 0.20 Å`, `--bias-k 300 eV·Å⁻²`, `--thresh baker`, `--baseline min`, `--out-dir ./result_scan2d/`. Grid size grows quickly as `(high − low) / --max-step-size` increases.
- **Next step:** Inspect `scan2d_map.png` / `scan2d_landscape.html` for a TS-region candidate, then refine with `tsopt` (or chain via `pdb2reaction all`).

`scan2d` constructs linear grids for both distances using `--max-step-size`, relaxes each grid point with the appropriate restraints active, and records unbiased MLIP energies for visualization. The default backend is UMA; select an alternative with `-b/--backend`. Use `--opt-mode hess` when you need RFOptimizer instead of LBFGS.

For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB/GJF output conversion.

## Minimal example
```bash
pdb2reaction scan2d -i input.pdb -q 0 -s scan2d.yaml -o ./result_scan2d/
```

## Output checklist
- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_i000_j000.xyz`
- `result_scan2d/scan2d_map.png` and `result_scan2d/scan2d_landscape.html`

## Common examples

1. **Run from a YAML spec file** -- see [Examples](#examples) below.
2. **Run with an inline literal** -- see [Examples](#examples) below.
3. **Enable `--dump`** to store inner trajectories by d1 step — see [Examples](#examples) below.

> **Note:** Add `--print-parsed` when you want to verify parsed pair targets from `--scan-lists/-s`.

## Usage
```bash
pdb2reaction scan2d -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan2d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### Examples
```bash
# Recommended: YAML/JSON spec file
cat > scan2d.yaml << 'YAML'
one_based: true
pairs:
 - ["TYR,285,CA", "SAM,309,C10", 1.30, 3.10]
 - ["TYR,285,CB", "SAM,309,C11", 1.20, 3.20]
YAML
pdb2reaction scan2d -i input.pdb -q 0 -s scan2d.yaml

# Alternative: inline Python literal
pdb2reaction scan2d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]'

# LBFGS, dumped inner trajectories, and Plotly outputs
pdb2reaction scan2d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]' \
 --max-step-size 0.20 --dump -o ./result_scan2d/ --opt-mode grad \
 --preopt --baseline min
```

## Scan-list spec

`scan2d` accepts exactly **two** quadruples `(i, j, low_Å, high_Å)` (under the `pairs` key for YAML/JSON, or as a single inline literal). Unlike `scan`, only **one literal** is accepted (no multi-stage support).

For the YAML/JSON file format, inline Python literal syntax, atom selectors, and quoting rules,
see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`.

## Workflow
1. Load the input geometry via `geom_loader`, resolve charge/spin, and optionally
    run an unbiased preoptimization when `--preopt`. If `-q` is omitted but
    `--ligand-charge/-l` is provided, the structure is treated as an enzyme–substrate
    complex and `extract.py`’s charge summary derives the total charge before the
    scan (for PDB inputs, or XYZ/GJF when `--ref-pdb` is supplied). The preoptimized
    structure is saved under `grid/preopt_i###_j###.*` and its unbiased energy is
    stored in `surface.csv` with indices `i = j = -1`.
2. Parse targets from `--scan-lists/-s` (YAML/JSON file or inline literal) into two quadruples, normalize indices
    (1-based by default). For PDB inputs, each atom entry can be an integer index
    or a selector string like `'TYR,285,CA'`; delimiters may be spaces, commas,
    slashes, backticks, or backslashes, and token order is flexible (fallback
    assumes resname, resseq, atom). Construct linear grids with
    `ceil(|high − low| / h) + 1` points (both endpoints included), where
    `h = --max-step-size`. Zero-length spans collapse to a single point.
    Each axis is then reordered so that the distance closest to the preoptimized
    geometry is indexed as `i = 0` / `j = 0`.
3. Iterate over every `d1[i]` (nearest-first ordering). For each value, relax
    the system with **only the d₁ restraint** active, snapshot that geometry, then
    run the inner loop over `d2[j]` with **both restraints** applied starting from
    the nearest previously converged structure.
4. At each `(i, j)` pair, store the biased-optimization result under
    `<out-dir>/grid/point_i###_j###.xyz`, record whether the bias converged, and
    evaluate the MLIP energy without bias. Optional per-outer-step inner
    trajectories are saved as `inner_path_d1_###_trj.xyz` when `--dump`.
5. After all points are visited, write `<out-dir>/surface.csv` with columns
    `i,j,d1_A,d2_A,energy_hartree,bias_converged,energy_kcal,d1_label,d2_label`, shifting the kcal
    reference via `--baseline {min|first}`. With `--baseline first`, the reference
    is the first grid entry (`i = j = 0` after reordering), not necessarily
    `(low₁, low₂)`. Generate `scan2d_map.png` (2D contour) and
    `scan2d_landscape.html` (3D surface) in `<out-dir>/`. Use `--zmin/--zmax` to
    clamp the color scale.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required |
| `-q, --charge INT` | Total charge (CLI > template/`--ligand-charge/-l`). Overrides `--ligand-charge/-l` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only; `workers_per_node` forwarded to the parallel predictor). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1. Inherits the `.gjf` template value when available; defaults to `1` when omitted. | `.gjf` template value or `1` |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (recommended) or **single** inline Python literal with two quadruples `(i,j,lowÅ,highÅ)`. `i`/`j` can be integer indices or PDB atom selectors like `'TYR,285,CA'`. | Required |
| `--one-based/--zero-based` | Interpret `(i, j)` indices as 1- or 0-based. | `True` |
| `--print-parsed/--no-print-parsed` | Print parsed pair tuples after `--scan-lists/-s` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum change allowed for either distance per increment (Å). Determines the grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². | `300` |
| `--relax-max-cycles INT` | Maximum optimizer cycles during each biased relaxation. Used unless YAML sets `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `grad` → LBFGS, `hess` → RFOptimizer. | `grad` |
| `--freeze-links/--no-freeze-links` | When the input is PDB, freeze parents of link hydrogens. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--dump/--no-dump` | Write `inner_path_d1_###_trj.xyz` for each outer step. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `-o, --out-dir TEXT` | Output directory root for grids and plots. | `./result_scan2d/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. | `False` |
| `--baseline {min,first}` | Shift kcal/mol energies so the global min or first grid point is zero. | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | Manual limits for the contour/surface color scale (kcal/mol). | Autoscaled |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical knobs to those documented for
  [YAML Reference](yaml-reference.md). `opt.dump` can be set in YAML for optimizer dumps;
  scan trajectory output is controlled by `--dump`.

### Section `bias`
- `k` (`300`): Harmonic strength in eV·Å⁻².

## Outputs
```
out_dir/ (default:./result_scan2d/)
├─ surface.csv # Structured grid table
├─ scan2d_map.png # 2D contour (requires Kaleido; the run stops if PNG export fails)
├─ scan2d_landscape.html # 3D surface visualization
├─ grid/point_i###_j###.xyz # Relaxed geometries for every (i, j) pair
├─ grid/point_i###_j###.pdb # PDB companions when conversion is enabled and templates exist
├─ grid/point_i###_j###.gjf # Gaussian companions when templates exist and conversion is enabled
├─ grid/preopt_i###_j###.xyz # Starting structure (present when --preopt is True)
├─ grid/preopt_i###_j###.pdb # PDB companion when conversion is enabled
├─ grid/preopt_i###_j###.gjf # Gaussian companion when templates exist and conversion is enabled
└─ grid/inner_path_d1_###_trj.xyz # Present only when --dump is True (mirrored to .pdb for PDB inputs with conversion)
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The MLIP backend (UMA by default) reuses the same
  `HarmonicBiasCalculator` as the 1D scan.
- Ångström limits are converted to Bohr internally to cap LBFGS steps and RFO
  trust radii; Optimizer scratch files live under temporary directories.
- The bias is always removed before final energies are recorded so you can reuse
  `surface.csv` in downstream fitting or visualization scripts.
- `--freeze-links` merges user `freeze_atoms` with detected link-H parents for
  PDB inputs, keeping extracted active site models rigid.


```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 device: auto # MLIP device selection
opt:
 thresh: baker # convergence preset (default: baker)
 max_cycles: 10000 # optimizer cycle cap
 dump: false # optimizer dumps (scan trajectories are controlled by --dump)
 out_dir: ./result_scan2d/ # output directory
lbfgs:
 max_step: 0.3 # maximum step length
 out_dir: ./result_scan2d/ # LBFGS-specific output directory
rfo:
 trust_radius: 0.10 # trust-region radius
 out_dir: ./result_scan2d/ # RFO-specific output directory
bias:
 k: 300.0 # harmonic bias strength (eV·Å⁻²)
```

More YAML options for `opt` are available in [YAML Reference](yaml-reference.md).
`--relax-max-cycles` applies only when explicitly provided **and** YAML does not set `opt.max_cycles` (default `10000`).

## See Also
- [scan](scan.md) -- 1D bond-distance scan
- [scan3d](scan3d.md) -- 3D distance-grid scan
- [opt](opt.md) -- single-structure optimization before/after scans
- [all](all.md) -- end-to-end workflow wrapper
- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
