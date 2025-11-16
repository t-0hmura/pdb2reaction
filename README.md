# pdb2reaction

`pdb2reaction` is a Python toolkit for building enzymatic reaction models directly from protein databank (PDB) structures and
running an automated pipeline of pocket extraction, reaction-path exploration, transition-state refinement, vibrational
analysis, and DFT single-point post-processing. The command line interface (CLI) wraps a set of UMA-based geometry optimisers
and graph string methods so that multi-step enzymatic reaction mechanisms can be generated with minimal manual intervention.

## Installation

```bash
# Clone this repository
git clone https://github.com/t-0hmura/pdb2reaction.git
cd pdb2reaction

# Install with Python 3.11+
pip install .
```

> **Note**: Several dependencies (e.g., `torch`, `fairchem-core`, `gpu4pyscf-cuda12x`) expect a CUDA-capable environment. Refer
to each project's installation guide when configuring GPUs or HPC nodes.

## Usage

The CLI is exposed via the `pdb2reaction` entry point. Running `pdb2reaction` without arguments is equivalent to invoking
`pdb2reaction all`, which executes the full pocket-to-product workflow:

1. Extract catalytic pockets from the supplied PDB structures.
2. (Optional) Run staged distance scans to build additional intermediates for single-structure inputs.
3. Search reaction paths with the GSM-based `path_search` engine and merge pocket structures back into the full models.
4. (Optional) Refine transition states, generate pseudo-IRC pathways, compute vibrational thermochemistry, and evaluate
   DFT single-point energies for every reactive segment.

### Default pipeline

```bash
# Example multi-structure workflow (reactant, intermediates, product)
pdb2reaction all -i R.pdb I1.pdb I2.pdb P.pdb \
                -c "A:123,B:456" \
                --ligand-charge 0 \
                --out-dir ./result_all \
                --tsopt True --thermo True --dft True
```

Key options:

- `-i, --input PATH...`: Reaction-ordered PDB files (≥2 for GSM mode, 1 when combined with `--scan-lists` or `--tsopt True`).
- `-c, --center TEXT`: Substrate specification (residue IDs or names) used to carve the binding pocket.
- `--ligand-charge TEXT`: Total or per-residue charge assignment propagated through scan/GSM/TSOPT.
- `--tsopt/--thermo/--dft BOOLEAN`: Enable TS optimisation, vibrational analysis, and DFT single-point post-processing.
- `--args-yaml FILE`: Shared YAML file that forwards UMA/GSM configuration blocks to every invoked subcommand
  (see [`docs/all.md`](docs/all.md)).

Refer to [`docs/all.md`](docs/all.md) for the exhaustive option list, YAML schema, and output description.

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

Each subcommand accepts `-h/--help` for inline usage hints and can also consume `--args-yaml` files that match the schemas
documented above.
