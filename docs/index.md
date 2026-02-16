# pdb2reaction Documentation

**pdb2reaction** is a Python CLI toolkit for automated enzymatic reaction-path modeling directly from PDB structures using machine-learning interatomic potentials (MLIPs).

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
concepts
troubleshooting
cli-conventions
ja/getting-started
ja/concepts
ja/troubleshooting
ja/cli-conventions
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
extract
fix_altloc
add_elem_info
opt
tsopt
path_opt
path_search
scan
scan2d
scan3d
freq
irc
dft
trj2fig
energy-diagram
ja/all
ja/extract
ja/fix_altloc
ja/add_elem_info
ja/opt
ja/tsopt
ja/path_opt
ja/path_search
ja/scan
ja/scan2d
ja/scan3d
ja/freq
ja/irc
ja/dft
ja/trj2fig
ja/energy-diagram
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

yaml-reference
uma_pysis
glossary
ja/yaml-reference
ja/uma_pysis
ja/glossary
```

```{toctree}
:maxdepth: 1
:caption: Language
:hidden:

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
| Draw state energy diagram from numeric values | `pdb2reaction energy-diagram` | [energy-diagram.md](energy-diagram.md) |
| Understand the big picture (concepts & terms) | — | [Concepts & Workflow](concepts.md) |
| Resolve common errors | — | [Troubleshooting](troubleshooting.md) |
| Look up abbreviations and terms | — | [Glossary](glossary.md) |

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
| [`energy-diagram`](energy-diagram.md) | Build an energy diagram from numeric input values |

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **UMA calculator settings** | [UMA Calculator](uma_pysis.md) |
| **Terminology** | [Glossary](glossary.md) |

---

## System Requirements

### Hardware
- **OS**: Linux (Ubuntu 20.04+ or CentOS 8+ tested)
- **GPU**: CUDA 12.x compatible
- **VRAM**: Minimum 8 GB (16 GB+ recommended for 1000+ atoms)
- **RAM**: 16 GB+ recommended

### Software
- Python 3.11
- PyTorch with CUDA support
- CUDA 12.x toolkit

---

## Quick Examples

### Basic MEP search
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### Full workflow with TS optimization
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True --thermo True --dft True
```

### Single-structure scan mode
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### TS-only optimization
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True
```

---

## Key Concepts

### Charge and spin
- Use `--ligand-charge` to specify unknown residue charges: `'SAM:1,GPP:-3'`
- Use `-q/--charge` to override the total charge
- Spin multiplicity is set with `-m/--mult` (the `all` command) or `-m/--multiplicity` (other subcommands); default is `1`

### Boolean options
All boolean CLI options must be explicitly set to `True` or `False`:
```bash
--tsopt True --thermo True --dft False
```

### YAML configuration
Advanced settings can be provided with `--args-yaml`.
```bash
pdb2reaction all -i R.pdb P.pdb -c 'LIG' --args-yaml config.yaml
```
See the [YAML Reference](yaml-reference.md) for all options.

---

## Output Structure

Typical `pdb2reaction all` output:
```
result_all/
├── summary.log              # Human-readable summary
├── summary.yaml             # Machine-readable summary
├── pockets/                 # Extracted cluster models
├── scan/                    # (Optional) scan results
├── path_search/             # MEP trajectories and diagrams
│   ├── mep.trj              # MEP trajectory
│   ├── mep.pdb              # MEP in PDB format
│   ├── mep_w_ref.pdb        # MEP merged with full system
│   ├── mep_plot.png         # Energy profile plot
│   └── seg_*/               # Per-segment details
└── path_search/post_seg_*/  # Post-processing outputs
    ├── tsopt/               # TS optimization results
    ├── irc/                 # IRC trajectories
    ├── freq/                # Vibrational modes
    └── dft/                 # DFT results
```

---

## Citation

A preprint describing `pdb2reaction` is in preparation. Please check back later for citation details.

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)** and is derived from Pysisyphus.

---

## References

1. Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. [arXiv:2506.23971](http://arxiv.org/abs/2506.23971)
2. Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). [DOI:10.1002/qua.26390](https://doi.org/10.1002/qua.26390)

---

## Getting Help

```bash
# General help
pdb2reaction --help

# Command help
pdb2reaction <subcommand> --help
```

For issues and feature requests, visit the [GitHub repository](https://github.com/t-0hmura/pdb2reaction).
