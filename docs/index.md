# pdb2reaction Documentation

*Version: v{{ version }}*

**pdb2reaction** is a Python CLI toolkit for automated enzymatic reaction-path modeling directly from PDB structures using machine-learning interatomic potentials (MLIPs).

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
installation
quickstart-all
quickstart-scan
quickstart-tsopt-freq
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
:maxdepth: 1
:caption: Language
:hidden:

ja/index
```


---

## Documentation Guide

| Topic | Page |
|-------|------|
| **Getting Started** | [Getting Started](getting-started.md) |
| **Installation** | [Installation](installation.md) |
| **Symptom-first failure routing** | [Common Error Recipes](recipes-common-errors.md) |
| **Common errors & fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI conventions & input requirements** | [CLI Conventions](cli-conventions.md) |

---

## Quick Start by Goal

| Objectives | Command | Guide |
|-------------------------|---------|-------|
| First run (end-to-end) | `pdb2reaction all` | [Quickstart: all](quickstart-all.md) |
| Single-structure staged scan + MEP + TS | `pdb2reaction all --scan-lists` | [Quickstart: scan](quickstart-scan.md) |
| TS optimization and validation | `pdb2reaction tsopt` | [Quickstart: tsopt](quickstart-tsopt-freq.md) |
| Run complete reaction path search from PDB | `pdb2reaction all` | [all](all.md) |
| Extract active site as cluster model from protein-ligand complex | `pdb2reaction extract` | [extract](extract.md) |
| Optimize a single structure | `pdb2reaction opt` | [opt](opt.md) |
| Find and optimize a transition state | `pdb2reaction tsopt` | [tsopt](tsopt.md) |
| Search for minimum energy path to obtain TS candidate | `pdb2reaction path-search` | [path-search](path-search.md) |
| Run IRC from a transition state | `pdb2reaction irc` | [irc](irc.md) |
| Visualize energy profile | `pdb2reaction trj2fig` | [trj2fig](trj2fig.md) |
| Draw state energy diagram from numeric values | `pdb2reaction energy-diagram` | [energy-diagram](energy-diagram.md) |
| Diagnose failures by symptom | — | [Common Error Recipes](recipes-common-errors.md) |
| Resolve common errors | — | [Troubleshooting](troubleshooting.md) |
| Look up abbreviations and terms | — | [Glossary](glossary.md) |
| Understand CLI conventions (flags, precedence, atom/residue selectors) | — | [CLI Conventions](cli-conventions.md) |

---

## CLI Subcommands

### Main Workflow
| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | End-to-end workflow: extraction → scan → MEP → TS optimization → IRC → thermochemistry → DFT |

### Structure Preparation
| Subcommand | Description |
|------------|-------------|
| [`extract`](extract.md) | Extract active site model (binding pocket) from protein–ligand complex |
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate locations |
| [`add-elem-info`](add-elem-info.md) | Repair PDB element columns (77–78) |

### Geometry Optimization
| Subcommand | Description |
|------------|-------------|
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS or RFO. [+ optional flatten]) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer or RS-I-RFO. [+ optional flatten]) |

### Path Search & Optimization
| Subcommand | Description |
|------------|-------------|
| [`path-opt`](path-opt.md) | Single-step MEP optimization via GSM or DMF (from 2 structures) |
| [`path-search`](path-search.md) | Recursive multi-step MEP search with automatic refinement (2+ structures) |

### Scans
| Subcommand | Description |
|------------|-------------|
| [`scan`](scan.md) | 1D bond-length driven scan with restraints |
| [`scan2d`](scan2d.md) | 2D distance grid scan |
| [`scan3d`](scan3d.md) | 3D distance grid scan |

### Analysis & Post-processing
| Subcommand | Description |
|------------|-------------|
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |
| [`energy-diagram`](energy-diagram.md) | Draw an energy diagram from numeric values |
| [`bond-summary`](bond-summary.md) | Detect and report covalent bond changes between consecutive structures |

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **MLIP backend settings** | [MLIP Calculator](uma-pysis.md) |
| **Terminology** | [Glossary](glossary.md) |

---

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

For setup, see [Installation](installation.md). For runnable command examples, workflow modes, and the output-directory layout, see [Getting Started](getting-started.md).

---

## Agent Skills

`pdb2reaction` ships AI-agent instructions under `.claude/skills/`
covering the CLI subcommands, structure I/O (PDB / XYZ / GJF), backend
installation (UMA / Orb / MACE / AIMNet2 / DFT / xtb), canonical
workflows, output parsing, and HPC operation. Copy the contents into
your project repository or `~/.claude/skills/` to make them available
to your agent platform (Claude Code, Cursor, etc.).

---

## Citation

A preprint describing `pdb2reaction` is in preparation. Currently, if you find this work helpful for your research, please cite the software itself:

```bibtex
@software{ohmura2026pdb2reaction,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  month        = {4},
  version      = {0.3.6},
  url          = {https://github.com/t-0hmura/pdb2reaction},
  license      = {GPL-3.0},
  doi          = {10.5281/zenodo.19197865}
}
```

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)**.

---

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

---

*This documentation is versioned with each release.*
