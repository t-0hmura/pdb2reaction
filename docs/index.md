# pdb2reaction Documentation

**Version: {{ version }}**

**pdb2reaction** is a Python CLI toolkit for automated enzymatic reaction-path modeling directly from PDB structures using machine-learning interatomic potentials (MLIPs).

```{toctree}
:maxdepth: 2
:caption: Contents
:hidden:

getting-started
concepts
cli-conventions
troubleshooting
all
extract
add_elem_info
opt
tsopt
path_opt
path_search
scan
scan2d
scan3d
irc
freq
dft
trj2fig
yaml-reference
uma_pysis
glossary
ja/index
```

---

## Quick Start by Goal

| What do you want to do? | Command | Guide |
|-------------------------|---------|-------|
| Run complete reaction path search from PDB | `pdb2reaction all` | [all.md](all.md) |
| Extract QM region from protein-ligand complex | `pdb2reaction extract` | [extract.md](extract.md) |
| Optimize a single structure | `pdb2reaction opt` | [opt.md](opt.md) |
| Find and optimize a transition state | `pdb2reaction tsopt` | [tsopt.md](tsopt.md) |
| Search for minimum energy path | `pdb2reaction path-search` | [path_search.md](path_search.md) |
| Run IRC from a transition state | `pdb2reaction irc` | [irc.md](irc.md) |
| Visualize energy profile | `pdb2reaction trj2fig` | [trj2fig.md](trj2fig.md) |

---

## Documentation Guide

| Topic | Page |
|-------|------|
| **Installation & first run** | [Getting Started](getting-started.md) |
| **Key terms & workflow overview** | [Concepts & Workflow](concepts.md) |
| **Common errors & fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI conventions & input requirements** | [CLI Conventions](cli-conventions.md) |

---

## CLI Subcommands

### Main Workflow
| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | End-to-end workflow: extraction → MEP → TS optimization → IRC → freq → DFT |

### Structure Preparation
| Subcommand | Description |
|------------|-------------|
| [`extract`](extract.md) | Extract active-site pocket (cluster model) from protein–ligand complex |
| [`add-elem-info`](add_elem_info.md) | Repair PDB element columns (77-78) |

### Geometry Optimization
| Subcommand | Description |
|------------|-------------|
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS / RFO) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer / RS-I-RFO) |

### Path Search & Optimization
| Subcommand | Description |
|------------|-------------|
| [`path-opt`](path_opt.md) | MEP optimization via GSM or DMF (two structures) |
| [`path-search`](path_search.md) | Recursive MEP search with automatic refinement (2+ structures) |

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

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **UMA calculator settings** | [UMA Calculator](uma_pysis.md) |
| **Terminology** | [Glossary](glossary.md) |

---

## Getting Help

```bash
pdb2reaction --help
pdb2reaction <subcommand> --help
```

For issues and feature requests, visit the [GitHub repository](https://github.com/t-0hmura/pdb2reaction).
