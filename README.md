# pdb2reaction: automated reaction-path modeling directly from PDB structures

## Overview

`pdb2reaction` is a Python CLI toolkit for turning **PDB structures** into **enzymatic reaction pathways** with machine-learning interatomic potentials (MLIPs).

A **single command** can generate a first-pass enzymatic reaction path:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

---

The full workflow — **MEP search → TS optimization → IRC → thermochemistry → single-point DFT** — can be run in one command:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

---

Given **(i) two or more PDB files** (R → ... → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt`**, `pdb2reaction` automatically:

- extracts an **active-site pocket** around user-defined substrates to build a **cluster model**,
- explores **minimum-energy paths (MEPs)** with GSM or DMF,
- *optionally* optimizes **transition states**, runs **vibrational analysis**, **IRC**, and **single-point DFT**,

using Meta's **UMA** machine-learning interatomic potential.

> **Expectation setting for TS search**
> - Treat single-command outputs as a strong initial guess, not guaranteed final TS validation.
> - Always validate TS candidates with frequency analysis and IRC before mechanistic interpretation.

> **Important (prerequisites):**
> - Input PDB files must already contain **hydrogen atoms**.
> - When providing multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ).
> - Boolean CLI options accept both `--flag` / `--no-flag` and value style `--flag True/False` (`yes/no`, `1/0` are also accepted). Prefer toggle style in new scripts.
> - The workflow also works for small-molecule systems. If you omit `--center/-c` and `--ligand-charge`, you can use `.xyz` or `.gjf` inputs as well.

## Documentation

- [**Getting Started**](docs/getting_started.md) — Quick start and workflow overview
- [**Installation**](docs/installation.md) — Setup and dependency installation
- [**Concepts & Workflow**](docs/concepts.md) — Key terms: pockets, templates, segments, and stages
- [**CLI Command Reference**](docs/reference/commands/index.md)
- [**YAML Schema**](docs/reference/yaml.md)
- [**Troubleshooting**](docs/troubleshooting.md) — Common errors and fixes
- **Full command index**: [docs/index.md](docs/index.md)

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
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
huggingface-cli login
```

### For DMF method

```bash
conda create -n pdb2reaction python=3.11 -y
conda activate pdb2reaction
conda install -c conda-forge cyipopt -y
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

For detailed installation instructions, see [Installation](docs/installation.md).

---

## Quick Examples

### Full workflow (multi-structure)
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

### Scan mode (single structure)
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### TS optimization only
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt
```

### Generate a starter config template
```bash
pdb2reaction init --out pdb2reaction_all.config.yaml
pdb2reaction all --config pdb2reaction_all.config.yaml --dry-run
```

### Step-by-step workflow
```bash
# 1. Extract active-site pocket (cluster model)
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' -r 6.0

# 2. Optimize geometry
pdb2reaction opt -i pocket.pdb

# 3. MEP search
pdb2reaction path-search -i R.pdb P.pdb

# 4. TS optimization
pdb2reaction tsopt -i hei.pdb

# 5. Frequency analysis
pdb2reaction freq -i ts_optimized.pdb

# 6. IRC
pdb2reaction irc -i ts_optimized.pdb

# 7. DFT single-point
pdb2reaction dft -i optimized.pdb
```

---

## CLI Subcommands

### Workflow

| Subcommand | Role | Documentation |
|---|---|---|
| `all` | End-to-end: extraction → MEP → TS → IRC → freq → DFT | [docs/all.md](docs/all.md) |
| `init` | Generate a starter YAML template for `pdb2reaction all` | [docs/init.md](docs/init.md) |

### Structure Preparation

| Subcommand | Role | Documentation |
|---|---|---|
| `extract` | Extract active-site pocket (cluster model) | [docs/extract.md](docs/extract.md) |
| `fix-altloc` | Resolve alternate conformations in PDB files | [docs/fix_altloc.md](docs/fix_altloc.md) |
| `add-elem-info` | Add/repair PDB element columns (77–78) | [docs/add_elem_info.md](docs/add_elem_info.md) |

### Optimization & Path Search

| Subcommand | Role | Documentation |
|---|---|---|
| `opt` | Geometry optimization (L-BFGS or RFO) | [docs/opt.md](docs/opt.md) |
| `tsopt` | TS optimization (Dimer or RS-I-RFO) | [docs/tsopt.md](docs/tsopt.md) |
| `path-opt` | MEP optimization via GSM or DMF | [docs/path_opt.md](docs/path_opt.md) |
| `path-search` | Recursive MEP search with refinement | [docs/path_search.md](docs/path_search.md) |
| `scan` | 1D bond-length driven scan | [docs/scan.md](docs/scan.md) |
| `scan2d` | 2D distance grid scan | [docs/scan2d.md](docs/scan2d.md) |
| `scan3d` | 3D distance grid scan | [docs/scan3d.md](docs/scan3d.md) |

### Analysis

| Subcommand | Role | Documentation |
|---|---|---|
| `freq` | Vibrational frequency analysis + thermochemistry | [docs/freq.md](docs/freq.md) |
| `irc` | IRC calculation (EulerPC) | [docs/irc.md](docs/irc.md) |
| `dft` | Single-point DFT (GPU4PySCF / PySCF) | [docs/dft.md](docs/dft.md) |

### Visualization

| Subcommand | Role | Documentation |
|---|---|---|
| `trj2fig` | Energy plot from XYZ trajectory | [docs/trj2fig.md](docs/trj2fig.md) |
| `energy-diagram` | Energy diagram from numeric values | [docs/energy_diagram.md](docs/energy_diagram.md) |


> **Tip:** In `tsopt`, `freq`, and `irc`, setting **`--hessian-calc-mode Analytical`** is strongly recommended when you have enough VRAM.

---

## HPC / Multi-GPU

On HPC clusters or multi-GPU workstations, `pdb2reaction` can parallelize UMA inference across nodes. Set `workers` and `workers_per_node` to enable parallel inference; see [docs/uma_pysis.md](docs/uma_pysis.md) for details.

---

## Getting Help

```bash
pdb2reaction --help
pdb2reaction <subcommand> --help
pdb2reaction <subcommand> --help-advanced
pdb2reaction all --help-advanced
# Equivalent module invocation
python -m pdb2reaction --help
```

`pdb2reaction all --help` shows core options. Use `pdb2reaction all --help-advanced` for the full option list.
`scan`, `scan2d`, `scan3d`, and the calculation commands (`opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`) now follow the same progressive-help pattern (`--help` core, `--help-advanced` full). `add-elem-info`, `trj2fig`, and `energy-diagram` also use the same pattern. `extract` and `fix-altloc` also support progressive help (`--help` core, `--help-advanced` full parser options).

> If you encounter any issues, please open an issue at <https://github.com/t-0hmura/pdb2reaction/issues>.

---

## Docs / Smoke Checks

To keep docs and CLI behavior in sync:

```bash
python scripts/check_intro_template.py
python scripts/check_all_scan_contract.py
python scripts/check_markdown_links.py
python scripts/smoke_docs_commands.py
```

To verify trajectory dump behavior across `opt` / `tsopt` routes:

```bash
# Optional: override per-case timeout in seconds (default: 120)
export PDB2REACTION_DUMP_CASE_TIMEOUT_SEC=180
python scripts/smoke_dump_trajectories.py
```

---

## Citation

A preprint describing `pdb2reaction` is in preparation. Please check back for citation details once it is available.

---

## References

[1] Wood, B. M., Dzamba, M., Fu, X., Gao, M., Shuaibi, M., Barroso-Luque, L., Abdelmaqsoud, K., Gharakhanyan, V., Kitchin, J. R., Levine, D. S., Michel, K., Sriram, A., Cohen, T., Das, A., Rizvi, A., Sahoo, S. J., Ulissi, Z. W., & Zitnick, C. L. (2025). UMA: A Family of Universal Models for Atoms. http://arxiv.org/abs/2506.23971
[2] Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. International Journal of Quantum Chemistry, 121(3). https://doi.org/10.1002/qua.26390

---

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)**.
