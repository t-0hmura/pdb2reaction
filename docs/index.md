# pdb2reaction Documentation

[GitHub](https://github.com/t-0hmura/pdb2reaction) · [ChemRxiv preprint](https://doi.org/10.26434/chemrxiv.15003538/v1) · [Open in Google Colab](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb)

*Version: v{{ release }}*

---

<img src="./overview.png" alt="pdb2reaction workflow overview" width="90%">

**pdb2reaction** is a Python CLI toolkit for exploring candidate enzymatic reaction pathways from PDB structures using machine-learning interatomic potentials (MLIPs).

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
installation
quickstart-all
quickstart-scan
quickstart-tsopt-freq
freeze-atoms
recipes-common-errors
troubleshooting
cli-conventions
cif
reproducibility
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
extract
fix-altloc
add-elem-info
opt
tsopt
path-opt
path-search
scan
scan2d
scan3d
freq
irc
dft
sp
trj2fig
energy-diagram
bond-summary
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

reference/commands/index
yaml-reference
json-output
uma-pysis
backends
architecture
output-layout
mcp_server
hpc-example
glossary
```

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

ja/index
ja/getting-started
ja/installation
ja/quickstart-all
ja/quickstart-scan
ja/quickstart-tsopt-freq
ja/freeze-atoms
ja/recipes-common-errors
ja/troubleshooting
ja/cli-conventions
ja/cif
ja/reproducibility
```

```{toctree}
:maxdepth: 2
:caption: コマンド
:hidden:

ja/all
ja/extract
ja/fix-altloc
ja/add-elem-info
ja/opt
ja/tsopt
ja/path-opt
ja/path-search
ja/scan
ja/scan2d
ja/scan3d
ja/freq
ja/irc
ja/dft
ja/sp
ja/trj2fig
ja/energy-diagram
ja/bond-summary
```

```{toctree}
:maxdepth: 2
:caption: リファレンス
:hidden:

ja/yaml-reference
ja/json-output
ja/uma-pysis
ja/backends
ja/architecture
ja/output-layout
ja/mcp_server
ja/hpc-example
ja/glossary
```

## Quick start

| Goal | Workflow |
|------|----------|
| **First end-to-end run** | [Quickstart: all](quickstart-all.md) |
| **Build an MEP from a single structure using a scan** | [Quickstart: scan](quickstart-scan.md) |
| **TS candidate available** | [Quickstart: TS-only mode](quickstart-tsopt-freq.md) |
| **Run failure / error** | [Common Error Recipes](recipes-common-errors.md) |

See [Installation](installation.md) for prerequisites.

## Subcommands

| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | Optional extraction; endpoint-MEP, scan-list, or TS-only entry mode; optional TS/IRC, thermochemistry, and DFT stages |
| [`extract`](extract.md) | Extract active site model (binding pocket) from protein–ligand complex |
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate locations |
| [`add-elem-info`](add-elem-info.md) | Repair PDB element columns (77–78) |
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS or RFO; optional flatten) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer or RS-P-RFO; optional flatten) |
| [`path-opt`](path-opt.md) | Single-step MEP optimization via GSM or DMF (from 2 structures) |
| [`path-search`](path-search.md) | Recursive multi-step MEP search with automatic refinement (2+ structures) |
| [`scan`](scan.md) | Restrained distance scan supporting concerted multi-distance and multistage scans |
| [`scan2d`](scan2d.md) | Two-dimensional energy-landscape exploration and PES mapping |
| [`scan3d`](scan3d.md) | Three-dimensional energy-landscape exploration and PES mapping |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`sp`](sp.md) | Single-point MLIP energy + forces / Hessian |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |
| [`energy-diagram`](energy-diagram.md) | Draw an energy diagram from numeric values |
| [`bond-summary`](bond-summary.md) | Detect and report covalent bond changes between consecutive structures |

## Configuration and reference

| Topic | Page |
|-------|------|
| **CLI conventions and input requirements** | [CLI Conventions](cli-conventions.md) · [mmCIF and large structures](cif.md) |
| **Cluster boundary atoms (cap hydrogens, `--freeze-atoms`)** | [Frozen Atoms](freeze-atoms.md) |
| **Common errors and fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **MLIP backend settings** | [MLIP Calculator](uma-pysis.md) |
| **Terminology** | [Glossary](glossary.md) |

## System requirements

### Hardware

- **OS:** Linux.
- **GPU (recommended):** an NVIDIA driver compatible with the backend and PyTorch wheel. CPU execution is also supported but slower.
- **VRAM / RAM:** depends on the model, system, Hessian mode, precision, and parallelism. Measure peak use on a representative calculation.

### Software

- Python >= 3.11.
- CPU or CUDA-enabled PyTorch. Prebuilt wheels include their CUDA runtime; a local toolkit is normally needed only for source builds.

See [Installation](installation.md) for setup.

## Agent skills

`skills/` contains guides for CLI commands, structure I/O, backends, workflows, output analysis, and HPC use.
See the [Skills index](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/README.md) for installation and the full list.

## Citation

```bibtex
@misc{ohmura2026pdb2reaction,
  author = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title  = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  year   = {2026}, doi = {10.26434/chemrxiv.15003538/v1}, note = {ChemRxiv preprint}
}
```

To cite the software or a specific release, use the Zenodo record:

```bibtex
@software{ohmura2026pdb2reaction_software,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  version      = {0.4.15},
  url          = {https://github.com/t-0hmura/pdb2reaction},
  license      = {GPL-3.0},
  doi          = {10.5281/zenodo.19197865}
}
```

## License

GNU General Public License v3 (GPL-3.0).

## Getting Help

```bash
# General help
pdb2reaction --help

# Command help
pdb2reaction <subcommand> --help

# Advanced options (dry-run, internal tuning, etc.)
pdb2reaction <subcommand> --help-advanced
```
