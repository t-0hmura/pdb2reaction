# pdb2reaction Documentation

*Version: v0.3.3*

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
ja/getting-started
ja/installation
ja/quickstart-all
ja/quickstart-scan
ja/quickstart-tsopt-freq
ja/recipes-common-errors
ja/troubleshooting
ja/cli-conventions
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
:caption: Reference
:hidden:

reference/commands/index
reference/yaml
yaml-reference
uma-pysis
glossary
ja/yaml-reference
ja/uma-pysis
ja/glossary
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
| Single-structure staged scan | `pdb2reaction scan` | [Quickstart: scan](quickstart-scan.md) |
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

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML schema** | [YAML Schema](reference/yaml.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **MLIP backend settings** | [MLIP Calculator](uma-pysis.md) |
| **Terminology** | [Glossary](glossary.md) |

---

## System Requirements

### Hardware
- **OS**: Linux
- **GPU**: CUDA 12.x compatible
- **VRAM**: Minimum 8 GB recommended
- **RAM**: 16 GB+ recommended

### Software
- Python >= 3.11
- PyTorch with CUDA support
- CUDA 12.x toolkit

---

## Quick Examples

### Basic MEP search
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

### Full workflow with TS optimization
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft
```

### Full workflow from reaction coordinates scan with single structure.
```bash
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]'
```

### Full workflow from single TS candidate structure.
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt
```

---

## Key flags in input

### Charge and spin
- Use `-l/--ligand-charge` to specify ligand charges: `'SAM:1,GPP:-3'`
- Use `-q/--charge` to override the net charge
- Spin multiplicity is set with `-m/--multiplicity` (default `1`)

### Boolean options
Boolean CLI options accept both toggle form (`--flag` / `--no-flag`) and explicit values (`--flag True/False`). Toggle form is recommended for new scripts:
```bash
--tsopt --thermo --no-dft
--tsopt True --thermo yes --dft 0
```

### YAML configuration
See the [YAML Reference](yaml-reference.md) for all options.

---

## Output Structure

Typical `pdb2reaction all` output:
```
result_all/
├── summary.log # User-friendly summary
├── summary.yaml # YAML-format summary
├── pockets/ # Extracted cluster models
├── scan/ # (Optional) scan results
├─┬ path_search/ # MEP trajectories and diagrams
│ ├── mep_trj.xyz # MEP trajectory
│ ├── mep.pdb # MEP in PDB format
│ ├── mep_w_ref.pdb # MEP merged with full system
│ ├── mep_plot.png # Energy profile plot
│ └── seg_*/ # Per-segment details
└┬── path_search/post_seg_*/ # Post-processing outputs
 ├── tsopt/ # TS optimization results
 ├── irc/ # IRC trajectories
 ├── freq/ # Vibrational analysis
 └── dft/ # DFT results
```

---

## Citation

A preprint describing `pdb2reaction` is in preparation. Please check back later for citation details.

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

*Note: This documentation is under active development. Some sections may be incomplete or subject to change.*
