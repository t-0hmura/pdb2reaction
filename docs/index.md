# pdb2reaction Documentation

*Version: v0.3.9*

---

<img src="./overview.png" alt="pdb2reaction workflow overview" width="90%">

**pdb2reaction** is a Python CLI toolkit for automated enzymatic reaction-path elucidation directly from PDB structures using machine-learning interatomic potentials (MLIPs).

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
ja/hpc-example
ja/glossary
```


## Start here

| Goal | Workflow |
|------|----------|
| **First end-to-end run** | [Quickstart: all](quickstart-all.md) |
| **Reactant only** | [Quickstart: scan](quickstart-scan.md) |
| **TS candidate available** | [Quickstart: TS-only mode](quickstart-tsopt-freq.md) |
| **Run failure / error** | [Common Error Recipes](recipes-common-errors.md) |

Refer to [Installation](installation.md) for prerequisites.

## Subcommands

| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | End-to-end workflow: extraction → scan → MEP → TS optimization → IRC → thermochemistry → DFT |
| [`extract`](extract.md) | Extract active site model (binding pocket) from protein–ligand complex |
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate locations |
| [`add-elem-info`](add-elem-info.md) | Repair PDB element columns (77–78) |
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS or RFO. [+ optional flatten]) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer or RS-I-RFO. [+ optional flatten]) |
| [`path-opt`](path-opt.md) | Single-step MEP optimization via GSM or DMF (from 2 structures) |
| [`path-search`](path-search.md) | Recursive multi-step MEP search with automatic refinement (2+ structures) |
| [`scan`](scan.md) | 1D bond-length driven scan with restraints |
| [`scan2d`](scan2d.md) | 2D distance grid scan |
| [`scan3d`](scan3d.md) | 3D distance grid scan |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |
| [`energy-diagram`](energy-diagram.md) | Draw an energy diagram from numeric values |
| [`bond-summary`](bond-summary.md) | Detect and report covalent bond changes between consecutive structures |

## Configuration & Reference

| Topic | Page |
|-------|------|
| **CLI conventions & input requirements** | [CLI Conventions](cli-conventions.md) |
| **Cluster boundary atoms (link hydrogens, `--freeze-atoms`)** | [Frozen Atoms](freeze-atoms.md) |
| **Common errors & fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **MLIP backend settings** | [MLIP Calculator](uma-pysis.md) |
| **Terminology** | [Glossary](glossary.md) |

## System Requirements

### Hardware
- **OS**: Linux
- **GPU**: CUDA 12.x compatible
- **VRAM**: 8 GB+ recommended
- **RAM**: 16 GB+ recommended

### Software
- Python >= 3.11
- PyTorch with CUDA support
- CUDA 12.x toolkit

For setup, see [Installation](installation.md).

## Agent Skills

`pdb2reaction` ships AI-agent instructions under `skills/` covering CLI subcommands, structure I/O, backend installation, workflows, output parsing, and HPC operation. See [`skills/README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/README.md) for the full skill index and installation instructions.

## Citation

A preprint describing `pdb2reaction` is in preparation. Currently, if you find this work helpful for your research, please cite the software itself:

```bibtex
@software{ohmura2026pdb2reaction,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  month        = {5},
  version      = {0.3.9},
  url          = {https://github.com/t-0hmura/pdb2reaction},
  license      = {GPL-3.0},
  doi          = {10.5281/zenodo.19197865}
}
```

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)**.

## Getting Help

```bash
# General help
pdb2reaction --help

# Command help
pdb2reaction <subcommand> --help

# Advanced options (dry-run, internal tuning, etc.)
pdb2reaction <subcommand> --help-advanced
```

For issues and feature requests, visit the [GitHub repository](https://github.com/t-0hmura/pdb2reaction).
