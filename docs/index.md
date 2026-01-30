# pdb2reaction Documentation

**pdb2reaction** is a Python CLI toolkit for automated enzymatic reaction-path modeling directly from PDB structures using Machine Learning Interatomic Potentials (MLIP).

```{toctree}
:maxdepth: 2
:caption: Contents
:hidden:

getting-started
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
```

---

## Quick Navigation

### Getting Started

- [**Getting Started**](getting-started.md) - Installation, quick start, and overview
- [**System Requirements**](#system-requirements) - Hardware and software prerequisites

### Main Workflow

- [`all`](all.md) - **End-to-end workflow**: extraction → scan → MEP search → TS optimization → IRC → thermochemistry → DFT

### CLI Subcommands

#### Structure Preparation
| Command | Description |
|---------|-------------|
| [`extract`](extract.md) | Extract active-site cluster from protein-ligand complex |
| [`add-elem-info`](add_elem_info.md) | Repair PDB element columns (77-78) |

#### Geometry Optimization
| Command | Description |
|---------|-------------|
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS / RFO) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer / RS-I-RFO) |

#### Path Search & Optimization
| Command | Description |
|---------|-------------|
| [`path-opt`](path_opt.md) | MEP optimization via GSM or DMF |
| [`path-search`](path_search.md) | Recursive MEP search with automatic refinement |

#### Scans
| Command | Description |
|---------|-------------|
| [`scan`](scan.md) | 1D bond-length driven scan with restraints |
| [`scan2d`](scan2d.md) | 2D distance grid scan |
| [`scan3d`](scan3d.md) | 3D distance grid scan |

#### Analysis & Post-processing
| Command | Description |
|---------|-------------|
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |

### Configuration

- [**YAML Reference**](yaml-reference.md) - Complete YAML configuration options
- [**UMA Calculator**](uma_pysis.md) - UMA machine-learning potential settings

---

## System Requirements

### Hardware
- **OS**: Linux (Ubuntu 20.04+, CentOS 8+ tested)
- **GPU**: CUDA 12.x compatible (RTX 30xx/40xx, A100, H100 tested)
- **VRAM**: Minimum 8 GB (16 GB+ recommended for >1000 atoms)
- **RAM**: 16 GB+ recommended

### Software
- Python 3.11
- PyTorch with CUDA support
- CUDA 12.x toolkit

---

## Quick Examples

### Basic MEP Search
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### Full Workflow with TS Refinement
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True --thermo True --dft True
```

### Single-Structure Scan Mode
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### TS Optimization Only
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True
```

---

## Key Concepts

### Charge and Spin
- Use `--ligand-charge` to specify charges for unknown residues: `'SAM:1,GPP:-3'`
- Use `-q/--charge` to override the total system charge
- Use `-m/--mult` (or `-m/--multiplicity` on subcommands) to set spin multiplicity (default: 1)

### Boolean Options
All boolean CLI options require explicit `True` or `False`:
```bash
--tsopt True --thermo True --dft False
```

### YAML Configuration
Advanced settings can be provided via `--args-yaml`:
```bash
pdb2reaction all -i R.pdb P.pdb -c 'LIG' --args-yaml config.yaml
```
See [YAML Reference](yaml-reference.md) for all options.

---

## Output Structure

A typical `pdb2reaction all` run produces:
```
result_all/
├── summary.log              # Human-readable summary
├── summary.yaml             # Machine-readable summary
├── pockets/                 # Extracted cluster models
├── scan/                    # (Optional) Scan results
├── path_search/             # MEP trajectories and diagrams
│   ├── mep.trj              # MEP trajectory
│   ├── mep.pdb              # MEP as PDB
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

A preprint describing `pdb2reaction` will be released soon. Please check back for citation details.

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)** derived from Pysisyphus.

---

## References

1. Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. [arXiv:2506.23971](http://arxiv.org/abs/2506.23971)
2. Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). [DOI:10.1002/qua.26390](https://doi.org/10.1002/qua.26390)

---

## Getting Help

```bash
# General help
pdb2reaction --help

# Subcommand help
pdb2reaction <subcommand> --help
```

For issues and feature requests, please visit the [GitHub repository](https://github.com/t-0hmura/pdb2reaction).
