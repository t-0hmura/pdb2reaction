# [pdb2reaction]{.p2r-wordmark} Documentation

:::{container} p2r-hero-meta
[Version: v{{ release }}]{.p2r-pill} [GitHub](https://github.com/t-0hmura/pdb2reaction){.p2r-meta-gh} · [ACS Omega paper](https://doi.org/10.1021/acsomega.6c08242){.p2r-meta-paper}
:::

:::{container} p2r-hero
<img src="./overview.png" alt="pdb2reaction workflow overview" class="p2r-hero-figure">

{.p2r-tagline}
**pdb2reaction** is a Python CLI toolkit for exploring candidate enzymatic reaction pathways from PDB structures using machine-learning interatomic potentials (MLIPs).

{.p2r-lead}
New to pdb2reaction? Start with [Getting Started](getting-started.md).

{.p2r-cta}
[Getting Started](getting-started.md){.p2r-btn .p2r-btn-primary} [Installation](installation.md){.p2r-btn .p2r-btn-install} [Open in Google Colab](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb){.p2r-btn .p2r-btn-colab}
:::

```{toctree}
:maxdepth: 2
:caption: Introduction
:hidden:

getting-started
installation
```

```{toctree}
:maxdepth: 2
:caption: Quickstart
:hidden:

Endpoint mode <quickstart-all>
Scan-list mode <quickstart-scan>
TS-only mode <quickstart-tsopt>
```

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

model-setup
mechanism-tips
dft-backend
troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
fix-altloc
add-elem-info
extract
opt
scan
scan2d
scan3d
path-opt
path-search
tsopt
irc
freq
dft
sp
bond-summary
trj2fig
energy-diagram
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

cli-conventions
reference/commands/index
yaml-reference
json-output
output-layout
backends
hpc-example
mcp_server
glossary
architecture
```

```{toctree}
:maxdepth: 2
:caption: 導入
:hidden:

ja/getting-started
ja/installation
```

```{toctree}
:maxdepth: 2
:caption: クイックスタート
:hidden:

Endpoint モード <ja/quickstart-all>
Scan-list モード <ja/quickstart-scan>
TS-only モード <ja/quickstart-tsopt>
```

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

ja/model-setup
ja/mechanism-tips
ja/dft-backend
ja/troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: コマンド
:hidden:

ja/all
ja/fix-altloc
ja/add-elem-info
ja/extract
ja/opt
ja/scan
ja/scan2d
ja/scan3d
ja/path-opt
ja/path-search
ja/tsopt
ja/irc
ja/freq
ja/dft
ja/sp
ja/bond-summary
ja/trj2fig
ja/energy-diagram
```

```{toctree}
:maxdepth: 2
:caption: リファレンス
:hidden:

ja/cli-conventions
ja/yaml-reference
ja/json-output
ja/output-layout
ja/backends
ja/hpc-example
ja/mcp_server
ja/glossary
ja/architecture
```

## Quick start

::::{container} p2r-cards
:::{container} p2r-card p2r-card-endpoint
**Analyze the mechanism end to end from the structures before and after the reaction**

<!-- p2r-mode-stages endpoint -->

[Quickstart: all in Endpoint mode](quickstart-all.md)
:::

:::{container} p2r-card p2r-card-scan
**Analyze the mechanism end to end from one structure**

<!-- p2r-mode-stages scan -->

[Quickstart: all in Scan-list mode](quickstart-scan.md)
:::

:::{container} p2r-card p2r-card-tsonly
**Analyze the mechanism end to end from a TS structure**

<!-- p2r-mode-stages tsonly -->

[Quickstart: TS-only mode](quickstart-tsopt.md)
:::
::::

| Goal | Page |
|------|------|
| **Build, trim, or extend the cluster model** | [Building the cluster model](model-setup.md) |
| **Study a mechanism, or the TS search fails** | [Tips for studying reaction mechanisms](mechanism-tips.md) |
| **Optimize the TS structure with DFT** | [Optimize the TS structure with DFT](dft-backend.md) |
| **A run failed** | [Troubleshooting](troubleshooting.md) |

## Subcommands

<!-- p2r-stage-strip -->

| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | Optional extraction; one of the three [input modes](getting-started.md#choosing-an-input-mode) (multi-structure MEP search, single structure + scan, TS-only mode); optional TS/IRC, thermochemistry, and DFT stages |
| [`extract`](extract.md) | Extract active site model (binding pocket) from protein–ligand complex |
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate locations |
| [`add-elem-info`](add-elem-info.md) | Repair PDB element columns (77–78) |
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS or RFO; optional `--flatten` removes leftover imaginary modes) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer or RS-P-RFO; optional `--flatten` removes extra imaginary modes) |
| [`path-opt`](path-opt.md) | Single-step MEP optimization via GSM or DMF (from 2 structures) |
| [`path-search`](path-search.md) | Recursive multi-step MEP search with automatic refinement (2+ structures) |
| [`scan`](scan.md) | Restrained distance scan supporting concerted multi-distance and multistage scans |
| [`scan2d`](scan2d.md) | Two-dimensional energy-landscape exploration and PES mapping |
| [`scan3d`](scan3d.md) | Three-dimensional energy-landscape exploration and PES mapping |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`sp`](sp.md) | Single-point calculator energy + forces / Hessian (MLIP or `-b dft`) |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |
| [`energy-diagram`](energy-diagram.md) | Draw an energy diagram from numeric values |
| [`bond-summary`](bond-summary.md) | Detect and report covalent bond changes between consecutive structures |

## Configuration and reference

| Topic | Page |
|-------|------|
| **Common options and input requirements** | [Common options and selectors](cli-conventions.md) |
| **Frozen atoms and distance restraints (`--freeze-atoms`, `--distance-restraint`)** | {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>` |
| **Common errors and fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **MLIP backend settings** | [MLIP Backends](backends.md) |
| **Files each command writes** | [Output Directory Layout](output-layout.md) |
| **Keys of `result.json` and `summary.json`** | [JSON Output Reference](json-output.md) |
| **Running on several GPU nodes (PBS + Ray)** | [HPC example](hpc-example.md) |
| **Calling pdb2reaction from an AI agent (MCP)** | [MCP server](mcp_server.md) |
| **Code structure (for developers)** | [Architecture](architecture.md) |
| **Terminology** | [Glossary](glossary.md) |

## System requirements

### Hardware

- **OS:** Linux (on Windows, install it in Linux under WSL2).
- **GPU:** an NVIDIA driver compatible with the backend and PyTorch wheel. CPU execution is also supported but slower.
- **VRAM / RAM:** depends on the model, the system size, and the Hessian mode; measure the peak on a representative run.

### Software

- Python >= 3.11.
- CPU or CUDA-enabled PyTorch. Prebuilt wheels include their CUDA runtime; a local toolkit is normally needed only for source builds.

See [Installation](installation.md) for setup.

## Agent skills

`skills/` contains guides for CLI commands, structure I/O, backends, workflows, output analysis, and HPC use.
To install them, tell your AI agent:

> Import `https://github.com/t-0hmura/pdb2reaction/tree/main/skills` as skills, and install pdb2reaction by following `pdb2reaction-install-backends`.

If you cloned the GitHub repository, you can give the local `skills/` path instead. Then you can ask, for example:

> Read *the paper*, build a model from the PDB structure *PDB ID*, and study the mechanism of *the reaction step* with the pdb2reaction skills.

## Citation

```bibtex
@article{ohmura2026pdb2reaction,
  author  = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title   = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  journal = {ACS Omega}, year = {2026}, doi = {10.1021/acsomega.6c08242}
}
```

To cite the software or a specific release, use the Zenodo record:

```bibtex
@software{ohmura2026pdb2reaction_software,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  version      = {0.5.0},
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

# Advanced options (internal tuning)
pdb2reaction <subcommand> --help-advanced
```

Report problems and feature requests on [GitHub Issues](https://github.com/t-0hmura/pdb2reaction/issues).
