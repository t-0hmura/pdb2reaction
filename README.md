# pdb2reaction: automated reaction-path modeling directly from PDB structures

## Overview

`pdb2reaction` is a Python CLI toolkit for turning **PDB structures** into **enzymatic reaction pathways** with machine-learning interatomic potentials (MLIPs).

In most cases, a **single command** like the one below is enough:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```
to model an enzymatic reaction pathway.

---
You can also run **MEP search → TS optimization → IRC → thermochemistry → single-point DFT calculations** in one command by adding `--tsopt True --thermo True --dft True`, for example:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True --thermo True --dft True
```
---

Given **(i) two or more full protein–ligand PDBs** `.pdb` (R → … → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt True`**, `pdb2reaction` automatically:

- extracts an **active-site pocket** around user‑defined substrates to build a **cluster model**,
- explores **minimum‑energy paths (MEPs)** with path optimization methods such as the Growing String Method (GSM) and Direct Max Flux (DMF),
- _optionally_ optimizes **transition states**, runs **vibrational analysis**, **IRC calculations**, and **single‑point DFT calculations**,

using Meta's UMA machine-learning interatomic potential (MLIP).

All of this is exposed through a command‑line interface (CLI) designed so that a **multi‑step enzymatic reaction mechanism** can be generated with minimal manual intervention. It can also handle small molecular systems. You can also use `.xyz` or `.gjf` inputs when you run workflows on full structures (i.e., omit `--center/-c` and `--ligand-charge`).

On **HPC clusters or multi‑GPU workstations**, `pdb2reaction` can process large cluster models (and optionally **full protein–ligand complexes**) by parallelizing UMA inference across nodes. Set `workers` and `workers_per_node` to enable parallel inference; see [`docs/uma_pysis.md`](docs/uma_pysis.md) for configuration details.

> **Important (prerequisites):**
> - Input PDB files must already contain **hydrogen atoms**.
> - When you provide multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ); otherwise an error is raised.
> - Boolean CLI options are passed explicitly as `True`/`False` (e.g., `--tsopt True`).


## Documentation

For detailed documentation, please refer to:

- [**Getting Started**](docs/getting-started.md) — Installation, quick start, and workflow overview
- [**Concepts & Workflow**](docs/concepts.md) — Mental model of pockets, templates, segments, and stages
- [**Troubleshooting**](docs/troubleshooting.md) — Common errors and fixes
- **Full command index**: [docs/index.md](docs/index.md) (English) / [docs/ja/index.md](docs/ja/index.md) (日本語)


***This software is still under development. Please use it at your own risk.***

---

## Quick Installation

`pdb2reaction` is intended for Linux environments with a CUDA‑capable GPU.

### Minimal setup (CUDA 12.9)

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

Log in to **Hugging Face Hub** to download UMA models:

```bash
huggingface-cli login
```

### For DMF method

If you want to use Direct Max Flux (DMF) for MEP search, install cyipopt first:

```bash
conda create -n pdb2reaction python=3.11 -y
conda activate pdb2reaction
conda install -c conda-forge cyipopt -y
pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

For detailed installation instructions, see [Getting Started](docs/getting-started.md).

---

## Quick Examples

### Multi-structure MEP search
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### Full workflow with TS optimization, thermochemistry, and DFT
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True --thermo True --dft True
```

### Single-structure scan mode
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### TS optimization only
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True
```

---

## CLI Subcommands

| Subcommand   | Role                                                                                  | Documentation                             |
| ------------ | ------------------------------------------------------------------------------------- | ----------------------------------------- |
| `all`        | End-to-end workflow: extraction → MEP search → TS optimization → IRC → freq → DFT              | [`docs/all.md`](docs/all.md)              |
| `extract`    | Extract active-site pocket (cluster model) from protein–ligand complex                               | [`docs/extract.md`](docs/extract.md)      |
| `opt`        | Single-structure geometry optimization (L-BFGS or RFO)                                | [`docs/opt.md`](docs/opt.md)              |
| `tsopt`      | Transition state optimization (Dimer or RS-I-RFO)                                     | [`docs/tsopt.md`](docs/tsopt.md)          |
| `path-opt`   | MEP optimization via GSM or DMF                                                       | [`docs/path_opt.md`](docs/path_opt.md)    |
| `path-search`| Recursive MEP search with automatic refinement                                        | [`docs/path_search.md`](docs/path_search.md) |
| `scan`       | 1D bond-length driven scan with restraints                                            | [`docs/scan.md`](docs/scan.md)            |
| `scan2d`     | 2D distance grid scan                                                                 | [`docs/scan2d.md`](docs/scan2d.md)        |
| `scan3d`     | 3D distance grid scan                                                                 | [`docs/scan3d.md`](docs/scan3d.md)        |
| `irc`        | IRC calculation with EulerPC                                                          | [`docs/irc.md`](docs/irc.md)              |
| `freq`       | Vibrational frequency analysis and thermochemistry                                    | [`docs/freq.md`](docs/freq.md)            |
| `dft`        | Single-point DFT using GPU4PySCF (with CPU PySCF fallback)                            | [`docs/dft.md`](docs/dft.md)              |
| `trj2fig`    | Plot ΔE or E from an XYZ trajectory                                                   | [`docs/trj2fig.md`](docs/trj2fig.md)      |
| `add-elem-info` | Add or repair PDB element columns (77–78)                                          | [`docs/add_elem_info.md`](docs/add_elem_info.md) |

> **Important:** Subcommands (except `all`) assume **cluster models** generated by `extract`. In these models, the atom closest to the Link‑H cap is automatically **frozen**. If you construct a cluster model yourself, set the Link‑H residue name to `LKH` and atom name to `HL`, or specify atoms to freeze via `--args-yaml` → `geom.freeze_atoms`.

> **Tip:** In `all`, `tsopt`, `freq`, and `irc`, setting **`--hessian-calc-mode Analytical`** is strongly recommended when you have enough VRAM.

---

## Getting Help

```bash
pdb2reaction --help
pdb2reaction <subcommand> --help
```

For detailed workflows, argument schemas, and example YAML files, consult the documentation files in `docs/`. For UMA calculator options, see [`docs/uma_pysis.md`](docs/uma_pysis.md).

> If you encounter any issues, please open an issue at <https://github.com/t-0hmura/pdb2reaction/issues>.

---

## Citation

A preprint describing `pdb2reaction` is in preparation. Please check back for citation details once it is available.

---

## References

[1] Wood, B. M., Dzamba, M., Fu, X., Gao, M., Shuaibi, M., Barroso-Luque, L., Abdelmaqsoud, K., Gharakhanyan, V., Kitchin, J. R., Levine, D. S., Michel, K., Sriram, A., Cohen, T., Das, A., Rizvi, A., Sahoo, S. J., Ulissi, Z. W., & Zitnick, C. L. (2025). UMA: A Family of Universal Models for Atoms. http://arxiv.org/abs/2506.23971
[2] Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. International Journal of Quantum Chemistry, 121(3). https://doi.org/10.1002/qua.26390

---

## License

`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)** derived from Pysisyphus.
