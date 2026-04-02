# pdb2reaction: automated reaction-path modeling directly from PDB structures

## Overview

`pdb2reaction` is a Python CLI toolkit for turning **PDB structures** into **enzymatic reaction pathways** with machine-learning interatomic potentials (MLIPs). Each workflow step is also available as an [individual subcommand](#cli-subcommands) ([`opt`](docs/opt.md), [`scan`](docs/scan.md), [`scan2d`](docs/scan2d.md), [`path-search`](docs/path_search.md), [`tsopt`](docs/tsopt.md), [`freq`](docs/freq.md), [`irc`](docs/irc.md), [`dft`](docs/dft.md), [`energy-diagram`](docs/energy_diagram.md), [etc.](#cli-subcommands)) for fine-grained control.

A **single command** can generate a first-pass enzymatic reaction path:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","SAM,309,C10",2.20)]'
```
---

The full workflow — **MEP search → TS optimization → IRC → thermochemistry → single-point DFT** — can be run in one command:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

---

Given **(i) two or more PDB files** (R → ... → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt`**, `pdb2reaction` automatically:

- extracts an **active-site pocket** around user-defined substrates to build a **cluster model**,
- explores **minimum-energy paths (MEPs)** with GSM or DMF,
- *optionally* optimizes **transition states**, runs **vibrational analysis**, **IRC**, and **single-point DFT**,

using machine-learning interatomic potentials (MLIPs).

### Related tools

| Tool | Use case | Repository |
|------|----------|------------|
| **mlmm-toolkit** | ML/MM (ONIOM) with full protein environment — automates MM parameter generation and ML region assignment from a single PDB input | <https://github.com/t-0hmura/mlmm_toolkit> |
| **UMA–Pysisyphus Interface** | YAML-input-based reaction mechanism analysis for small molecules | <https://github.com/t-0hmura/uma_pysis> |

Both `pdb2reaction` and `mlmm-toolkit` include a custom GPU-optimized pysisyphus fork for geometry optimization, TS search, and IRC. This bundled fork is **not compatible** with the upstream pysisyphus package; do not install them side by side.

> **Important (prerequisites):**
> - Input PDB files must already contain **hydrogen atoms**.
> - When providing multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ).
> - Boolean CLI options accept both `--flag` / `--no-flag` and value style `--flag True/False` (`yes/no`, `1/0` are also accepted). Prefer toggle style in new scripts.
> - The workflow also works for small-molecule systems. If you omit `--center/-c` and `--ligand-charge`, you can use `.xyz` or `.gjf` inputs as well.

## Documentation

- [**Getting Started**](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/getting_started.md) — Quick start and workflow overview
- [**Installation**](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/installation.md) — Setup and dependency installation
- [**Concepts & Workflow**](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/concepts.md) — Key terms: pockets, templates, segments, and stages
- [**CLI Command Reference**](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/reference/commands/index.md)
- [**YAML Schema**](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/reference/yaml.md)
- [**Troubleshooting**](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/troubleshooting.md) — Common errors and fixes
- **Full command index**: [docs/index.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/index.md)

***This software is still under development. Please use it at your own risk.***

---

## Installation

`pdb2reaction` requires Linux with a CUDA-capable GPU.

### Prerequisites

- Python >= 3.11
- CUDA 12.x

### Minimal setup (CUDA 12.9, torch 2.8.0)

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
pip install pdb2reaction
plotly_get_chrome -y
huggingface-cli login
```

### For DMF method

```bash
conda create -n pdb2reaction python=3.11 -y
conda activate pdb2reaction
conda install -c conda-forge cyipopt -y
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
pip install pdb2reaction
plotly_get_chrome -y
```

### DFT single-point (`pdb2reaction dft`)

DFT dependencies are **not** installed by default. To use `pdb2reaction dft`, install the `[dft]` extra:

```bash
pip install "pdb2reaction[dft]"
```

This installs PySCF, GPU4PySCF (x86_64 only), and related CUDA libraries. Note that DFT single-point calculations are practical only for systems up to **~500 atoms**; larger systems will require prohibitive compute time and memory.

For detailed installation instructions, see [Installation](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/installation.md).

### Supported ML potentials

| Potential | Repository | Install extra |
|-----------|------------|---------------|
| **UMA** (default) | <https://github.com/facebookresearch/fairchem> | *(included)* |
| **ORB** | <https://github.com/orbital-materials/orb-models> | `pip install "pdb2reaction[orb]"` |
| **MACE** | <https://github.com/ACEsuit/mace> | `pip install "pdb2reaction[mace]"` |
| **AIMNet2** | <https://github.com/isayevlab/aimnetcentral> | `pip install "pdb2reaction[aimnet2]"` |

---

## Quick Examples

### Full workflow (multi-structure)
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

### Scan mode (single structure)
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","SAM,309,C10",2.20)]'
```

### TS optimization only
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --tsopt
```

### Step-by-step workflow

**1. Extract active-site pocket (cluster model)** — [`extract`](docs/extract.md)
```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 6.0
```

**2. Optimize geometry** — [`opt`](docs/opt.md)
```bash
pdb2reaction opt -i pocket.pdb -l 'SAM:1,GPP:-3'
```

**3. MEP search** — [`path-search`](docs/path_search.md)
```bash
pdb2reaction path-search -i R.pdb P.pdb -l 'SAM:1,GPP:-3'
```

**4. TS optimization** — [`tsopt`](docs/tsopt.md)
```bash
pdb2reaction tsopt -i hei.pdb -l 'SAM:1,GPP:-3'
```

**5. Frequency analysis** — [`freq`](docs/freq.md)
```bash
pdb2reaction freq -i ts_optimized.pdb -l 'SAM:1,GPP:-3'
```

**6. IRC** — [`irc`](docs/irc.md)
```bash
pdb2reaction irc -i ts_optimized.pdb -l 'SAM:1,GPP:-3'
```

**7. DFT single-point** — [`dft`](docs/dft.md)
```bash
pdb2reaction dft -i optimized.pdb -l 'SAM:1,GPP:-3'
```

---

## CLI Subcommands

### Workflow

| Subcommand | Role | Documentation |
|---|---|---|
| `all` | End-to-end: extraction → MEP → TS → IRC → freq → DFT | [docs/all.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/all.md) |

### Structure Preparation

| Subcommand | Role | Documentation |
|---|---|---|
| `extract` | Extract active-site pocket (cluster model) | [docs/extract.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/extract.md) |
| `fix-altloc` | Resolve alternate conformations in PDB files | [docs/fix_altloc.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/fix_altloc.md) |
| `add-elem-info` | Add/repair PDB element columns (77–78) | [docs/add_elem_info.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/add_elem_info.md) |

### Optimization & Path Search

| Subcommand | Role | Documentation |
|---|---|---|
| `opt` | Geometry optimization (L-BFGS or RFO) | [docs/opt.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/opt.md) |
| `tsopt` | TS optimization (Dimer or RS-I-RFO) | [docs/tsopt.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/tsopt.md) |
| `path-opt` | MEP optimization via GSM or DMF | [docs/path_opt.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/path_opt.md) |
| `path-search` | Recursive MEP search with refinement | [docs/path_search.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/path_search.md) |
| `scan` | 1D bond-length driven scan | [docs/scan.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/scan.md) |
| `scan2d` | 2D distance grid scan | [docs/scan2d.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/scan2d.md) |
| `scan3d` | 3D distance grid scan | [docs/scan3d.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/scan3d.md) |

### Analysis

| Subcommand | Role | Documentation |
|---|---|---|
| `freq` | Vibrational frequency analysis + thermochemistry | [docs/freq.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/freq.md) |
| `irc` | IRC calculation (EulerPC) | [docs/irc.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/irc.md) |
| `dft` | Single-point DFT (GPU4PySCF / PySCF) | [docs/dft.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/dft.md) |

### Visualization

| Subcommand | Role | Documentation |
|---|---|---|
| `trj2fig` | Energy plot from XYZ trajectory | [docs/trj2fig.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/trj2fig.md) |
| `energy-diagram` | Energy diagram from numeric values | [docs/energy_diagram.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/energy_diagram.md) |


> **Tip:** In [`tsopt`](docs/tsopt.md), [`freq`](docs/freq.md), and [`irc`](docs/irc.md), setting **`--hessian-calc-mode Analytical`** is strongly recommended when you have enough VRAM.

---

## HPC / Multi-GPU

On HPC clusters or multi-GPU workstations, `pdb2reaction` can parallelize UMA inference across nodes. Set `workers` and `workers_per_node` to enable parallel inference; see [docs/uma_pysis.md](https://github.com/t-0hmura/pdb2reaction/blob/main/docs/uma_pysis.md) for details.

---

## Getting Help

```bash
pdb2reaction --help
pdb2reaction <subcommand> --help
pdb2reaction <subcommand> --help-advanced
pdb2reaction all --help-advanced
# Shorthand alias (equivalent to pdb2reaction)
p2r --help
# Equivalent module invocation
python -m pdb2reaction --help
```

`pdb2reaction all --help` shows core options. Use `pdb2reaction all --help-advanced` for the full option list.
`scan`, `scan2d`, `scan3d`, and the calculation commands (`opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`) now follow the same progressive-help pattern (`--help` core, `--help-advanced` full). `add-elem-info`, `trj2fig`, and `energy-diagram` also use the same pattern. `extract` and `fix-altloc` also support progressive help (`--help` core, `--help-advanced` full parser options).

> If you encounter any issues, please open an issue at <https://github.com/t-0hmura/pdb2reaction/issues>.

---

## Citation

A preprint describing `pdb2reaction` is in preparation. Currently, if you find this work helpful for your research, please cite the software itself:

```bibtex
@software{ohmura2026pdb2reaction,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  month        = {3},
  version      = {0.3.2},
  url          = {https://github.com/t-0hmura/pdb2reaction},
  license      = {GPL-3.0},
  doi          = {10.5281/zenodo.19197878}
}
```

---
<!-- 
## References

[1] Wood, B. M., Dzamba, M., Fu, X., Gao, M., Shuaibi, M., Barroso-Luque, L., Abdelmaqsoud, K., Gharakhanyan, V., Kitchin, J. R., Levine, D. S., Michel, K., Sriram, A., Cohen, T., Das, A., Rizvi, A., Sahoo, S. J., Ulissi, Z. W., & Zitnick, C. L. (2025). UMA: A Family of Universal Models for Atoms. https://arxiv.org/abs/2506.23971  
[2] Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. International Journal of Quantum Chemistry, 121(3). https://doi.org/10.1002/qua.26390

--- -->

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)**.  

***This software is still under development. Please use it at your own risk.***  