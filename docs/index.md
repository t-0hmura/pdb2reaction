# pdb2reaction Documentation

*Version: v0.4.0* — Python CLI for enzymatic reaction-path elucidation from PDB structures using machine-learning interatomic potentials (MLIPs).

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
| Run failure / error | [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) |
| CLI conventions / YAML / Glossary | [CLI Conventions](cli-conventions.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md) |
| Bit-reproducible runs (`--deterministic`) | [Reproducibility](reproducibility.md) |
| MLIP backend settings / HPC examples | [MLIP Calculator](uma-pysis.md) · [HPC Examples](hpc-example.md) |
| Cluster boundary atoms (link H, `--freeze-atoms`) | [Frozen Atoms](freeze-atoms.md) |

## Subcommands

| Subcommand | Role |
|---|---|
| [`all`](all.md) | End-to-end: extract → scan → MEP → TS → IRC → freq → DFT |
| [`extract`](extract.md) · [`fix-altloc`](fix-altloc.md) · [`add-elem-info`](add-elem-info.md) | Structure preparation |
| [`opt`](opt.md) · [`tsopt`](tsopt.md) | Geometry / TS optimisation |
| [`path-opt`](path-opt.md) · [`path-search`](path-search.md) | MEP via GSM/DMF / recursive refinement |
| [`scan`](scan.md) · [`scan2d`](scan2d.md) · [`scan3d`](scan3d.md) | 1D / 2D / 3D bond-distance scans |
| [`freq`](freq.md) · [`irc`](irc.md) | Vibrational analysis + thermochemistry / IRC (EulerPC) |
| [`dft`](dft.md) · [`sp`](sp.md) | Single-point DFT / single-point MLIP energy + forces |
| [`bond-summary`](bond-summary.md) | Bond-change report between consecutive structures |
| [`trj2fig`](trj2fig.md) · [`energy-diagram`](energy-diagram.md) | Energy plot / R→TS→P diagram |

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
