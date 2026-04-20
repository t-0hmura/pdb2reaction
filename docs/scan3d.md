# `scan3d`

## Overview

> **Summary:** Perform a three-distance (d₁, d₂, d₃) grid scan with harmonic restraints and MLIP relaxations. Use `-s/--scan-lists` with a YAML/JSON spec file (recommended) or an inline Python literal; or plot an existing `surface.csv` via `--csv`.

### At a glance
- **Input:** One structure + `-s/--scan-lists scan3d.yaml` (recommended) or one `-s/--scan-lists` inline literal (three quadruples), unless you use `--csv` to plot only.
- **Grid ordering:** Values are reordered so points closest to the (pre)optimized structure are visited first.
- **Energies:** Recorded energies are evaluated **without bias**, so grid points are directly comparable.
- **Outputs:** `surface.csv`, per-point geometries under `grid/`, and an HTML isosurface plot (`scan3d_density.html`).
- **Caution:** 3D grids grow very quickly; consider coarser `--max-step-size` or smaller ranges first.

`scan3d` nests loops over d₁ → d₂ → d₃ and relaxes each point with the appropriate restraints active. The default optimizer is LBFGS (`--opt-mode grad`); switch to `--opt-mode hess` for RFOptimizer.

For XYZ/GJF inputs, `--ref-pdb` supplies a reference PDB topology while keeping XYZ coordinates, enabling format-aware PDB/GJF output conversion.

## Minimal example
```bash
pdb2reaction scan3d -i input.pdb -q 0 -s scan3d.yaml -o ./result_scan3d/
```

## Output checklist
- `result_scan3d/surface.csv`
- `result_scan3d/grid/point_i000_j000_k000.xyz`
- `result_scan3d/scan3d_density.html`

## Common examples
1. **Run from a YAML spec file** -- see [Examples](#examples) below.
2. **Run with an inline literal** -- see [Examples](#examples) below.
3. **Enable `--dump`** to keep inner d3 trajectories per `(d1,d2)` slice -- see [Examples](#examples) below.

> **Note:** Add `--print-parsed` when you want to verify parsed pair targets from `-s/--scan-lists`.

## Usage
```bash
pdb2reaction scan3d [-i INPUT.{pdb|xyz|trj|...}] [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan3d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] [--csv PATH]
```
Note: `-i/--input` and `-s/--scan-lists` are required unless `--csv` is provided.

### Examples
```bash
# Recommended: YAML/JSON spec file
cat > scan3d.yaml << 'YAML'
one_based: true
pairs:
 - ["TYR,285,CA", "SAM,309,C10", 1.30, 3.10]
 - ["TYR,285,CB", "SAM,309,C11", 1.20, 3.20]
 - ["TYR,285,CG", "SAM,309,C12", 1.10, 3.00]
YAML
pdb2reaction scan3d -i input.pdb -q 0 -s scan3d.yaml

# Alternative: inline Python literal
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]'

# LBFGS relaxations, dumped inner trajectories, and an HTML isosurface plot
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]' \
 --max-step-size 0.20 --dump -o ./result_scan3d/ --opt-mode grad \
 --preopt --baseline min

# Plot only from an existing surface.csv (skip new energy evaluation)
pdb2reaction scan3d --csv ./result_scan3d/surface.csv --zmin -10 --zmax 40 -o ./result_scan3d/
```

## YAML/JSON spec file format (recommended)

```yaml
one_based: true # optional; defaults to CLI --one-based
pairs:
 - [1, 5, 1.30, 3.10]
 - [2, 8, 1.20, 3.20]
 - [3, 12, 1.10, 3.00]
```

- `pairs` is required and must contain exactly 3 quadruples.
- Each quadruple is `(i, j, low_Å, high_Å)`.
- Indices may be integers or PDB selectors, same as inline literals.

## Inline Python literal format

`-s/--scan-lists` also accepts a **single inline Python literal** string. Shell quoting matters.

### Basic structure

The literal is a Python list of exactly **three** quadruples `(atom1, atom2, low_Å, high_Å)`:

```
-s '[(atom1, atom2, low_Å, high_Å), (atom3, atom4, low_Å, high_Å), (atom5, atom6, low_Å, high_Å)]'
```

- Wrap the entire literal in **single quotes** so the shell does not interpret parentheses or spaces.
- Each quadruple defines one scan axis: the distance between `atom1`–`atom2` is scanned from `low_Å` to `high_Å`.
- Unlike `scan`, only **one literal** is accepted (no multi-stage support).

### Specifying atoms

Atoms can be given as **integer indices** or **PDB selector strings**:

| Method | Example | Notes |
| --- | --- | --- |
| Integer index | `(1, 5, 1.30, 3.10)` | 1-based by default (`--one-based`) |
| PDB selector | `("TYR,285,CA", "SAM,309,C10", 1.30, 3.10)` | Residue name, residue number, atom name |

PDB selector tokens can be separated by any of: comma `,`, space, slash `/`, backtick `` ` ``, or backslash `\`. Token order is flexible.

```bash
# All of these specify the same atom:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # order is flexible
```

### Quoting rules

```bash
# Correct: single-quote the list, double-quote selector strings inside
-s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]'

# Correct: integer indices need no inner quotes
-s '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20), (3, 12, 1.10, 3.00)]'

# Avoid: double-quoting the outer literal requires escaping inner quotes
-s "[(\"TYR,285,CA\",\"SAM,309,C10\",1.30,3.10),...]"
```

## Workflow
1. Load the structure through `geom_loader`, resolve charge/spin from CLI or
    embedded Gaussian templates, and optionally run an unbiased preoptimization
    when `--preopt`. If `-q` is omitted but `--ligand-charge` is provided, the
    structure is treated as an enzyme–substrate complex and `extract.py`’s charge
    summary derives the total charge before scanning (for PDB inputs, or XYZ/GJF
    when `--ref-pdb` is supplied).
2. Parse targets from `-s/--scan-lists` (YAML/JSON file or inline literal; default 1-based indices unless
    `--zero-based` is passed) into three quadruples. For PDB inputs, each
    atom entry can be an integer index or a selector string like `'TYR,285,CA'`;
    delimiters may be spaces, commas, slashes, backticks, or backslashes, and
    token order is flexible (fallback assumes resname, resseq, atom). Build each linear grid using
    `h = --max-step-size` and reorder the values so the ones closest to the
    starting distances are visited first.
3. Outer loop over `d1[i]`: relax with only the d₁ restraint active, starting
    from the previously scanned geometry whose d₁ value is closest. Snapshot that
    structure.
4. Middle loop over `d2[j]`: relax with d₁ and d₂ restraints, starting from the
    closest (d₁, d₂) geometry. Snapshot that result.
5. Inner loop over `d3[k]`: relax with all three restraints, measure the
    unbiased energy (bias removed for evaluation), and write the constrained
    geometry and convergence flag.
6. After the scan completes, assemble `surface.csv`, apply the kcal/mol
    baseline shift (`--baseline {min|first}`), and generate a 3D RBF-interpolated
    isosurface plot (`scan3d_density.html`) honoring `--zmin/--zmax`. When
    `--csv` is provided, only this plotting step runs.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required unless `--csv` is provided |
| `-q, --charge INT` | Total charge (CLI > template/`--ligand-charge`). Overrides `--ligand-charge` when both are set. | Required unless template/derivation applies |
| `-l, --ligand-charge TEXT` | Per-residue charge mapping (e.g., `GPP:-3,SAM:1`). Automatically derives the total system charge from PDB residue charges — no manual counting needed. Used when `-q` is omitted (PDB inputs or XYZ/GJF with `--ref-pdb`). | _None_ |
| `--workers`, `--workers-per-node` | MLIP predictor parallelism (workers > 1 disables analytic Hessians; UMA backend only; `workers_per_node` forwarded to the parallel predictor). See {ref}`workers-fd-downgrade` for diagnostic notes. | `1`, `1` |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1. Inherits the `.gjf` template value when available; defaults to `1` when omitted. | `.gjf` template value or `1` |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (recommended) or **single** inline Python literal with three quadruples `(i,j,lowÅ,highÅ)`. `i`/`j` can be integer indices or PDB atom selectors like `'TYR,285,CA'`. | Required unless `--csv` is provided |
| `--one-based/--zero-based` | Interpret `(i, j)` indices as 1- or 0-based. | `True` |
| `--print-parsed/--no-print-parsed` | Print parsed pair tuples after `-s/--scan-lists` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum change allowed per distance increment (Å). Controls grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². | `300` |
| `--relax-max-cycles INT` | Maximum optimizer cycles during each biased relaxation. Used unless YAML sets `opt.max_cycles`. | `10000` |
| `--opt-mode TEXT` | `grad` → LBFGS, `hess` → RFOptimizer. | `grad` |
| `--freeze-links/--no-freeze-links` | When the input is PDB, freeze parents of link hydrogens. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| `--dump/--no-dump` | Write `inner_path_d1_###_d2_###_trj.xyz` for each (d₁, d₂). | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/GJF companions for PDB/Gaussian inputs. | `True` |
| `--ref-pdb FILE` | Reference PDB topology to use when the input is XYZ/GJF (keeps XYZ coordinates). | _None_ |
| `-o, --out-dir TEXT` | Output directory root for grids and plots. | `./result_scan3d/` |
| `--csv PATH` | Load an existing `surface.csv` and only plot it (no new scan). `-i/--input` and `-s/--scan-lists` become optional. | _None_ |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--solvent TEXT` | Implicit solvent name for xTB correction (e.g. `water`). `none` to disable. | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB solvent model. | `alpb` |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. **Scope-dependent default:** `False` standalone; flipped to `True` when invoked via `pdb2reaction all` (see [`all` → Scan Options](all.md#scan-options-single-input-runs)). | `False` |
| `--baseline {min,first}` | Shift kcal/mol energies so the global min or `(i,j,k)=(0,0,0)` is zero. | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | Manual limits for the isosurface color bands (kcal/mol). | Autoscaled |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical knobs to those documented for
  [YAML Reference](yaml-reference.md). `opt.dump` can be set in YAML for optimizer dumps;
  scan trajectory output is controlled by `--dump`.

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
 out_dir: ./result_scan3d/ # output directory
lbfgs:
 max_step: 0.3 # maximum step length
 out_dir: ./result_scan3d/ # LBFGS-specific output directory
rfo:
 trust_radius: 0.10 # trust-region radius
 out_dir: ./result_scan3d/ # RFO-specific output directory
bias:
 k: 300.0 # harmonic bias strength (eV·Å⁻²)
```

`--relax-max-cycles` applies only when explicitly provided **and** YAML does not set `opt.max_cycles` (default `10000`).

### Section `bias`
- `k` (`300`): Harmonic strength in eV·Å⁻².

## Outputs
```
out_dir/ (default:./result_scan3d/)
├─ surface.csv # Grid metadata; may include a reference row (i=j=k=-1)
├─ scan3d_density.html # 3D energy isosurface visualization
├─ grid/point_i###_j###_k###.xyz # Relaxed geometry for each grid point (Å×100 tags)
├─ grid/point_i###_j###_k###.pdb # PDB companions when conversion is enabled and templates exist
├─ grid/point_i###_j###_k###.gjf # Gaussian companions when templates exist and conversion is enabled
├─ grid/preopt_i###_j###_k###.xyz # Starting structure saved before scanning (preoptimized when --preopt is True)
└─ grid/inner_path_d1_###_d2_###_trj.xyz # Present only when --dump is True (mirrored to .pdb for PDB inputs with conversion)
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The MLIP backend (UMA by default) reuses the same
  `HarmonicBiasCalculator` as the 1D/2D scans.
- Ångström limits are converted to Bohr internally to cap LBFGS steps and RFO
  trust radii; optimizer scratch files live under temporary directories.
- `--baseline` defaults to the global minimum; `--baseline first` anchors the
  `(i,j,k)=(0,0,0)` grid point when present.
- 3D visualization uses RBF interpolation on a 50×50×50 grid with
  semi-transparent step-colored isosurfaces (no cross-sectional planes).
- `--freeze-links` merges user `freeze_atoms` with detected link-H parents for
  PDB inputs, keeping extracted active site models rigid.

## See Also
- [scan](scan.md) -- 1D bond-distance scan
- [scan2d](scan2d.md) -- 2D distance-grid scan
- [opt](opt.md) -- single-structure optimization before/after scans
- [all](all.md) -- end-to-end workflow wrapper
- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
