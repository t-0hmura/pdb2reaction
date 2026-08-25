# `scan3d`

Perform a three-distance grid scan with harmonic restraints and machine-learned interatomic potential (MLIP) relaxations, sampling the potential energy on a 3D grid over the three distances `(d₁, d₂, d₃)`. Use it when such a volume is needed, or when an existing `surface.csv` needs re-plotting.

There are two ways to drive the command. To run a new scan, supply targets with `--scan-lists/-s`, either as a YAML/JSON spec file (recommended) or as an inline Python literal. To re-plot an existing `surface.csv` without re-evaluating energies, pass it via `--csv`. During a scan, `scan3d` nests loops over d₁ → d₂ → d₃ and relaxes each point with the appropriate harmonic restraints active.

The default optimizer is L-BFGS (`--opt-mode grad`); switch to `--opt-mode hess` for RFOptimizer.

mmCIF inputs use the internal PDB bridge and emit CIF with restored IDs. For XYZ/GJF inputs, `--ref-pdb` accepts a PDB or mmCIF reference topology while keeping XYZ coordinates.

## Examples

Command form:

```bash
pdb2reaction scan3d [-i INPUT.{pdb|xyz|trj|...}] [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] \
 [-s/--scan-lists scan3d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] [--csv PATH]
```

Recommended: YAML/JSON spec file.

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
```

Alternative: inline Python literal.

```bash
# Alternative: inline Python literal
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]'
```

Plot only from an existing `surface.csv` (skip new energy evaluation).

```bash
# Plot only from an existing surface.csv (skip new energy evaluation)
pdb2reaction scan3d --csv ./result_scan3d/surface.csv --zmin -10 --zmax 40 -o ./result_scan3d/
```

## Workflow

1. Load the structure through `geom_loader`, resolve charge/spin from CLI or
    embedded Gaussian templates, and optionally run an unbiased preoptimization
    when `--preopt`. If `-q` is omitted but `--ligand-charge/-l` is provided, the
    structure is treated as an enzyme–substrate complex and `extract.py`’s charge
    summary derives the total charge before scanning (for PDB/mmCIF inputs, or XYZ/GJF
    when `--ref-pdb` is supplied).
2. Parse targets from `--scan-lists/-s` (YAML/JSON file or inline literal; default 1-based indices unless
    `--zero-based` is passed) into three quadruples. For PDB/mmCIF topology inputs, each
    atom entry can be an integer index or a selector string like `'TYR,285,CA'`;
    delimiters may be spaces, commas, slashes, backticks, or backslashes, and
    token order is flexible; use positional `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` for
    repeated names or numbering. Build each linear grid using
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
6. After the scan completes, assemble `surface.csv` (columns:
    `i,j,k,d1_A,d2_A,d3_A,energy_hartree,bias_converged,is_preopt,energy_kcal,d1_label,d2_label,d3_label`),
    apply the kcal/mol baseline shift (`--baseline {min|first}`), and generate a
    3D RBF-interpolated isosurface plot (`scan3d_density.html`) honoring
    `--zmin/--zmax`. When `--csv` is provided, only this plotting step runs.

Plot-only input requires `d1_A`, `d2_A`, `d3_A`, and either
`energy_hartree` or `energy_kcal`. Rows marked `is_preopt=true`, explicitly
unconverged rows, and non-finite rows are excluded; older CSVs without
provenance are accepted. Interpolation requires at least four
unique non-coplanar usable points spanning all three axes.

## Outputs

The key deliverables are `surface.csv`, per-point geometries under `grid/`, and the `scan3d_density.html` isosurface plot.
With `--out-json`, a fresh scan records an explicit `grid_points[]` mapping in
`result.json`. The isosurfaces are interpolated; interactive clients should use
the mapped computed grid points, rather than arbitrary isosurface vertices, when
selecting a corresponding structure.

```
out_dir/ (default:./result_scan3d/)
├─ surface.csv # Grid metadata; may include a reference row (i=j=k=-1)
├─ scan3d_density.html # 3D energy isosurface visualization (you can open it with a web browser)
├─ grid/point_i###_j###_k###.xyz # Relaxed geometry for each grid point (Å×100 tags)
├─ grid/point_i###_j###_k###.pdb # PDB companions when conversion is enabled and templates exist
├─ grid/point_i###_j###_k###.cif # Bridge-input companions with original IDs
├─ grid/point_i###_j###_k###.gjf # Gaussian companions when templates exist and conversion is enabled
├─ grid/preopt_i###_j###_k###.xyz # Starting structure saved before scanning (preoptimized when --preopt is True)
└─ grid/inner_path_d1_###_d2_###_trj.xyz # Present only when --dump is True (format companions require topology + conversion)
```

Grid-point geometries use `Å×100` tags, so `point_i130_j310_k200.xyz` corresponds to d₁=1.30, d₂=3.10, d₃=2.00 Å. If two points have the same rounded tags, later filenames append `_grid_III_JJJ_KKK` with the zero-based grid indices.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| **Input & charge** | | |
| `-i, --input PATH` | Structure file accepted by `geom_loader`. | Required unless `--csv` is provided |
| `-q, --charge INT` | Total charge for a fresh scan (CLI > template/`--ligand-charge/-l`). Overrides `--ligand-charge/-l` when both are set. Not needed in plot-only `--csv` mode. | Required for a fresh scan unless template/derivation applies |
| `-l, --ligand-charge TEXT` | For a fresh scan, either a scalar integer (e.g., `-1`) for the total ligand charge, or a per-residue mapping (e.g., `GPP:-3,SAM:1`) that derives the total from PDB/mmCIF residue metadata. Used when `-q` is omitted. Not used in plot-only `--csv` mode. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity 2S+1. Inherits the `.gjf` template value when available; defaults to `1` when omitted. | `.gjf` template value or `1` |
| **Backend & compute** | | |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend. | `uma` |
| `--workers`, `--workers-per-node` | UMA predictor parallelism; `workers_per_node` is forwarded to the parallel predictor. `workers > 1` cannot be combined with an explicit analytical Hessian request. See {ref}`workers-analytical-error`. | `1`, `1` |
| **Active-region freezing** | | |
| `--freeze-links/--no-freeze-links` | When the input is PDB, freeze parents of cap hydrogens. | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze explicitly (e.g., `'1,3,5'`). Complements `--freeze-links`; applies to any input format. | _None_ |
| **Scan targets** | | |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (recommended) or **single** inline Python literal with three quadruples `(i,j,lowÅ,highÅ)`. `i`/`j` can be integer indices, three-field selectors, or positional `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`. | Required unless `--csv` is provided |
| `--one-based/--zero-based` | Interpret `(i, j)` indices as 1- or 0-based. | `True` |
| `--print-parsed/--no-print-parsed` | Print parsed pair tuples after `--scan-lists/-s` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum change allowed per distance increment (Å). Controls grid density. | `0.20` |
| **Refinement** | | |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV·Å⁻². | `300` |
| `--opt-mode TEXT` | `grad` → L-BFGS, `hess` → RFOptimizer. | `grad` |
| `--relax-max-cycles INT` | Maximum optimizer cycles during each biased relaxation. An explicit value overrides YAML `opt.max_cycles`. | `100000` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `baker` |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. | `False` |
| **Merge & alignment** | | |
| `--ref-pdb FILE` | Reference PDB or mmCIF topology for XYZ/GJF input (keeps XYZ coordinates). | _None_ |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ → PDB/CIF/GJF companions. | `True` |
| **Output & config** | | |
| `-o, --out-dir TEXT` | Output directory root for grids and plots. | `./result_scan3d/` |
| `--csv PATH` | Load an existing `surface.csv` and only plot it (no new scan). `-i/--input` and `--scan-lists/-s` become optional. | _None_ |
| `--dump/--no-dump` | Write `inner_path_d1_###_d2_###_trj.xyz` for each (d₁, d₂). | `False` |
| `--baseline {min,first}` | Shift Hartree-backed energies so the global min or first eligible grid point is zero. With kcal-only `--csv` input, retain the supplied zero. | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | Manual limits for the isosurface color bands (kcal/mol). | Autoscaled |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. See [JSON Output Schema](json-output.md) for the schema. | `False` |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |

The full flag list is in the generated [command reference](reference/commands/index.md).

## YAML configuration

### Shared YAML sections
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: identical knobs to those documented for
  [YAML Reference](yaml-reference.md), except run-scoped `opt.dump` is ignored.
  Use `--dump` for scan trajectory output. Set the command-owned output
  directory with `-o/--out-dir`; optimizer `out_dir` YAML keys are ignored.

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/cap detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p2 # uma-s-1p2 | uma-m-1p1
 device: auto # MLIP device selection
opt:
 thresh: baker # convergence preset (default: baker)
 max_cycles: 100000 # optimizer cycle cap
 dump: false # optimizer dumps (scan trajectories are controlled by --dump)
lbfgs:
 max_step: 0.3 # maximum step length
rfo:
 trust_radius: 0.10 # trust-region radius
bias:
 k: 300.0 # harmonic bias strength (eV·Å⁻²)
```

## Notes

- `scan3d` accepts exactly **three** quadruples `(i, j, low_Å, high_Å)` (under the `pairs` key for YAML/JSON, or as a single inline literal). Unlike `scan`, only **one literal** is accepted (no multi-stage support). For the YAML/JSON file format, inline Python literal syntax, atom selectors, and quoting rules, see {ref}`CLI Conventions: Scan-list spec <scan-list-spec>`.
- 3D grids grow very quickly; consider coarser `--max-step-size` or smaller ranges first.
- The MLIP backend (UMA by default) reuses the same
  `HarmonicBiasCalculator` as the 1D/2D scans.
- Ångström limits are converted to Bohr internally to cap L-BFGS steps and RFO
  trust radii; optimizer scratch files live under temporary directories.
- `--baseline` defaults to the global minimum; `--baseline first` anchors the
  `(i,j,k)=(0,0,0)` grid point when eligible, otherwise it falls back to the
  eligible minimum. For plot-only CSV containing `energy_kcal` but no
  `energy_hartree`, the supplied zero is retained.
- 3D visualization uses RBF interpolation on a 50×50×50 grid with
  semi-transparent step-colored isosurfaces (no cross-sectional planes).
- `--freeze-links` merges user `freeze_atoms` with detected cap-H parents for
  PDB inputs, keeping extracted active site models rigid.
- Add `--print-parsed` when you want to verify parsed pair targets from `--scan-lists/-s`.
- An explicitly provided `--relax-max-cycles` overrides YAML `opt.max_cycles`; when omitted, YAML wins, then the default `100000` applies.

## See Also
- [scan](scan.md) -- Concerted multi-distance and multistage restrained scans
- [scan2d](scan2d.md) -- 2D distance-grid scan
- [opt](opt.md) -- single-structure optimization before/after scans
- [all](all.md) -- end-to-end workflow wrapper
- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
