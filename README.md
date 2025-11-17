# pdb2reaction

`pdb2reaction` is a Python toolkit for building enzymatic reaction models directly from protein databank (PDB) structures and
running an automated pipeline of pocket extraction, reaction-path exploration, transition-state refinement, vibrational
analysis, and DFT single-point post-processing. The command line interface (CLI) wraps a set of UMA-based geometry optimisers
and graph string methods so that multi-step enzymatic reaction mechanisms can be generated with minimal manual intervention.

## Installation

```bash
pip install torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128
pip install git+https://github.com/t-0hmura/pdb2reaction.git
huggingface-cli login
```

> **Note**: Several dependencies (e.g., `torch`, `fairchem-core`, `gpu4pyscf-cuda12x`) expect a CUDA-capable environment. Refer
to each project's installation guide when configuring GPUs or HPC nodes.

## Usage

The CLI is exposed via the `pdb2reaction` entry point declared in `pyproject.toml`. The Click group in
[`pdb2reaction/cli.py`](pdb2reaction/cli.py) sets `all` as its default subcommand, so running `pdb2reaction ...` or
`pdb2reaction all ...` is equivalent. All workflows require `-i/--input` (full PDBs in reaction order) and `-c/--center`
(substrate definition for pocket extraction).

### Workflow modes

#### Multi-structure GSM pipeline

```bash
pdb2reaction -i R.pdb I1.pdb I2.pdb P.pdb \
             -c "A:123,B:456" \
             --ligand-charge "GPP:-3,MMT:-1" \
             --out-dir ./result_all \
             --tsopt True --thermo True --dft True
```

Multiple full systems are processed in reaction order. The command extracts pockets, runs the recursive GSM-based
`path_search`, merges the pocket minimum-energy path back into the full structures, and (optionally) executes TS
optimisation, vibrational analysis, and DFT single points for each reactive segment.

#### Single-structure + staged scan (feeds GSM)

```bash
pdb2reaction -i SINGLE.pdb \
             -c "308,309" \
             --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
             --spin 1 --out-dir ./result_scan_all \
             --tsopt True --thermo True --dft True
```

Providing exactly one input PDB alongside `--scan-lists` performs a staged distance scan **on the extracted pocket** using
UMA. Each `stage_XX/result.pdb` becomes an ordered intermediate/product candidate that feeds the subsequent `path_search`
run before the results are merged back into the full system.

#### Single-structure TSOPT-only mode

```bash
pdb2reaction -i SINGLE.pdb \
             -c "GPP,MMT" \
             --ligand-charge "GPP:-3,MMT:-1" \
             --tsopt True --thermo True --dft True \
             --out-dir ./result_tsopt_only
```

Supplying a single input **without** `--scan-lists` while setting `--tsopt True` skips the path search entirely. The tool
optimises the pocket TS, performs a pseudo-IRC to minimise both ends, and can run `freq`/`dft` on the resulting R/TS/P trio
to build UMA, Gibbs, DFT, and DFT//UMA diagrams.

### Important flags and behaviours

- `-i/--input PATH...`: Two or more full PDBs for GSM mode, or one PDB when paired with `--scan-lists` **or** `--tsopt True`.
  A single `-i` may be followed by multiple filenames.
- `-c/--center TEXT`: Substrate definition for pocket extraction. Accepts PDB paths, residue IDs (e.g., `A:123,B:456`), or
  residue names (`GPP,MMT`).
- `--ligand-charge TEXT`: Total charge or mapping (e.g., `GPP:-3,MMT:-1`). The first pocket’s total charge is rounded to an
  integer and reused for scan/GSM/TS optimisation.
- `--scan-lists TEXT...`: One or more Python-style lists describing staged scans for single-input runs. Each `(i, j, target_Å)`
  tuple uses indices from the original full PDB (1-based) and is auto-remapped onto the extracted pocket.
- `--tsopt/--thermo/--dft BOOLEAN`: Enable TS optimisation + pseudo-IRC, vibrational analysis (UMA Gibbs diagram), and DFT
  single-point post-processing (adds a DFT//UMA Gibbs diagram when combined with `--thermo True`).
- `--sopt-mode` / `--opt-mode`: Optimiser presets. `--opt-mode` is forwarded to the staged scan and TS optimisation; when it is
  omitted, the scan inherits LBFGS or RFO based on `--sopt-mode`, and TS optimisation defaults to `light`.
- Shared knobs (`--freeze-links`, `--dump`, `--pre-opt`, `--hessian-calc-mode`, etc.) propagate to scan, GSM, TS optimisation,
  and post-processing stages.
- `--args-yaml FILE`: A YAML file forwarded unchanged to `path_search`, `scan`, `ts_opt`, `freq`, and `dft`, allowing shared
  UMA/GSM configuration blocks (see [`docs/all.md`](docs/all.md)).
- Outputs: `<out-dir>/pockets/` holds extracted pockets, `<out-dir>/scan/` stores staged scan results,
  `<out-dir>/path_search/` contains GSM outputs (plus per-segment TS/freq/DFT folders), and `<out-dir>/tsopt_single/` is used by
  TSOPT-only runs.

Single-input runs require either `--scan-lists` (staged scan feeding GSM) or `--tsopt True` (TSOPT-only mode). Refer to
[`docs/all.md`](docs/all.md) for a full option matrix, YAML schemas, and output details.

### CLI subcommands

| Subcommand | Summary | Documentation |
| --- | --- | --- |
| `all` | End-to-end workflow orchestrator that chains pocket extraction, GSM search, TS/freq/DFT post-processing, and staged scans. | [`docs/all.md`](docs/all.md) |
| `scan` | Perform staged biased scans on pocket models to create additional intermediates. | [`docs/scan.md`](docs/scan.md) |
| `opt` | Optimise a single structure with UMA (LBFGS/RFO presets). | [`docs/opt.md`](docs/opt.md) |
| `path-opt` | Run UMA optimisation on a specific path segment or snapshot. | [`docs/path_opt.md`](docs/path_opt.md) |
| `path-search` | Launch the GSM-based reaction path search and pocket/full-system merging. | [`docs/path_search.md`](docs/path_search.md) |
| `ts-opt` | Refine transition states (with optional pseudo-IRC propagation). | [`docs/ts_opt.md`](docs/ts_opt.md) |
| `freq` | Compute vibrational modes, thermochemistry, and UMA energy diagrams. | [`docs/freq.md`](docs/freq.md) |
| `irc` | Follow an intrinsic reaction coordinate starting from a TS structure. | [`docs/irc.md`](docs/irc.md) |
| `extract` | Extract catalytic pockets from full PDB structures (also used implicitly by `all`). | [`docs/extract.md`](docs/extract.md) |
| `trj2fig` | Convert trajectory data to interactive diagrams (Plotly/Kaleido). | [`docs/trj2fig.md`](docs/trj2fig.md) |
| `add-elem-info` | Augment PDB files with missing element metadata. | [`docs/add_elem_info.md`](docs/add_elem_info.md) |
| `dft` | Run single-point DFT calculations on UMA geometries using PySCF/gpu4pyscf. | [`docs/dft.md`](docs/dft.md) |
| `2d-scan` | Explore two harmonic distance restraints simultaneously and build 2D PES grids with UMA relaxations. | [`docs/scan2d.md`](docs/scan2d.md) |

Each subcommand accepts `-h/--help` for inline usage hints and can also consume `--args-yaml` files that match the schemas
documented above.
