# pdb2reaction Documentation

*Version: v0.4.4* — Python CLI for enzymatic reaction-path elucidation from PDB structures using machine-learning interatomic potentials (MLIPs).

<img src="./overview.png" alt="pdb2reaction workflow overview" width="90%">

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

| Goal | Page |
|---|---|
| Install + run a first end-to-end pipeline | [Installation](installation.md) · [Getting Started](getting-started.md) |
| End-to-end pipeline from a PDB | [Quickstart: all](quickstart-all.md) |
| Reactant only — staged distance scan | [Quickstart: scan](quickstart-scan.md) |
| TS candidate available — `tsopt → IRC → freq` | [Quickstart: TS-only](quickstart-tsopt-freq.md) |
| Choosing precision / TS route / imaginary-mode fix / controlled comparison | [`tsopt`](tsopt.md) |
| Staged-vs-concerted scan / barrier direction | [`scan`](scan.md) |
| Run failure / error | [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) |
| CLI conventions / YAML / Glossary | [CLI Conventions](cli-conventions.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md) |
| Bit-reproducible runs (`--deterministic`) | [Reproducibility](reproducibility.md) |
| MLIP backend settings / HPC examples | [MLIP Calculator](uma-pysis.md) · [HPC Examples](hpc-example.md) |
| Cluster boundary atoms (cap H, `--freeze-atoms`) | [Frozen Atoms](freeze-atoms.md) |

## Subcommands

| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | End-to-end workflow: extraction → scan → MEP → TS optimization → IRC → thermochemistry → DFT |
| [`extract`](extract.md) | Extract active site model (binding pocket) from protein–ligand complex |
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate locations |
| [`add-elem-info`](add-elem-info.md) | Repair PDB element columns (77–78) |
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS or RFO; optional flatten) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer or RS-I-RFO; optional flatten) |
| [`path-opt`](path-opt.md) | Single-step MEP optimization via GSM or DMF (from 2 structures) |
| [`path-search`](path-search.md) | Recursive multi-step MEP search with automatic refinement (2+ structures) |
| [`scan`](scan.md) | 1D bond-length driven scan with restraints |
| [`scan2d`](scan2d.md) | 2D distance grid scan |
| [`scan3d`](scan3d.md) | 3D distance grid scan |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`sp`](sp.md) | Single-point MLIP energy + forces / Hessian |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |
| [`energy-diagram`](energy-diagram.md) | Draw an energy diagram from numeric values |
| [`bond-summary`](bond-summary.md) | Detect and report covalent bond changes between consecutive structures |

## Getting Help

```bash
# General help
pdb2reaction --help

# Command help
pdb2reaction <subcommand> --help

# Advanced options (dry-run, internal tuning, etc.)
pdb2reaction <subcommand> --help-advanced
```

## Citation

```bibtex
@misc{ohmura2026pdb2reaction,
  author = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title  = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  year   = {2026}, doi = {10.26434/chemrxiv.15003538}, note = {ChemRxiv preprint}
}
```

A Zenodo software record is also available (DOI `10.5281/zenodo.19197865`).

## License

GNU General Public License v3 (GPL-3.0).
