# Getting Started

## Overview

`pdb2reaction` is a Python CLI toolkit for turning **PDB structures** into **enzymatic reaction pathways** using machine-learning interatomic potentials (MLIPs).

In most cases, a **single command** like the one below is enough to model an enzymatic reaction pathway:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

---
You can also run **MEP search → TS refinement → IRC → thermochemistry → DFT single-point** in one go by adding `--tsopt True --thermo True --dft True`:
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True --thermo True --dft True
```
---

Given **(i) two or more full protein–ligand PDBs** `.pdb` (R → … → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt True`**, `pdb2reaction` automatically:

- extracts an **active site** around user‑defined substrates to build a **cluster model**,
- explores **minimum‑energy paths (MEPs)** with path optimization methods such as the Growing String Method (GSM) and Direct Max Flux (DMF),
- _optionally_ refines **transition states**, runs **vibrational analysis**, **IRC calculations**, and **DFT single‑point calculations**.

All calculations use Meta's UMA machine-learning interatomic potential (MLIP).

All of this is exposed through a command‑line interface (CLI) designed so that a **multi‑step enzymatic reaction mechanism** can be generated with minimal manual intervention. It can also handle small molecular systems. You can also use `.xyz` or `.gjf` inputs when you run workflows on full structures (i.e., omit `--center/-c` and `--ligand-charge`).

On **HPC clusters or multi‑GPU workstations**, `pdb2reaction` can process large cluster models (and optionally **full protein–ligand complexes**) by parallelizing UMA inference itself across nodes. Set `workers` and `workers_per_node` to enable parallel calculation; see [UMA Calculator](uma_pysis.md) for configuration details.

```{important}
- Input PDB files must already contain **hydrogen atoms**.
- When you provide multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ); otherwise an error is raised.
- Boolean CLI options are passed explicitly as `True`/`False` (e.g., `--tsopt True`).
```

### Recommended tools for hydrogen addition

If your PDB lacks hydrogen atoms, use one of the following tools before running pdb2reaction:

| Tool | Example Command | Notes |
|------|-----------------|-------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | Fast, widely used for crystallographic structures |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | Adds hydrogens and assigns partial charges |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | General-purpose cheminformatics toolkit |

To ensure identical atom ordering across multiple PDB inputs, apply the same hydrogen-addition tool with consistent settings to all structures.

```{warning}
This software is still under development. Please use it at your own risk.
```

---

## Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters) with a CUDA‑capable GPU. Several dependencies – notably **PyTorch**, **fairchem‑core (UMA)**, and **gpu4pyscf‑cuda12x** – expect a working CUDA installation.

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

### Quick start

Below is a minimal setup example that works on many CUDA 12.9 clusters. Adjust module names and versions to match your system. This example assumes the default GSM MEP mode (no DMF). For DMF, install cyipopt via conda first.

```bash
# 1) Install a CUDA-enabled PyTorch build
# 2) Install pdb2reaction from GitHub
# 3) Install a headless Chrome for Plotly figure export

pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded. Either:

```bash
# Hugging Face CLI
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

or

```bash
# Classic CLI
huggingface-cli login
```

You only need to do this once per machine / environment.

- If you want to use the Direct Max Flux (DMF) method for MEP search, create a conda environment and install cyipopt before installing pdb2reaction.
  ```bash
  # Create and activate a dedicated conda environment
  conda create -n pdb2reaction python=3.11 -y
  conda activate pdb2reaction

  # Install cyipopt (required for the DMF method in MEP search)
  conda install -c conda-forge cyipopt -y
  ```

- If you are on an HPC cluster that uses *environment modules*, load CUDA **before** installing PyTorch, like this:
  ```bash
  module load cuda/12.9
  ```


### Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load CUDA (if you use environment modules on an HPC cluster)**

   ```bash
   module load cuda/12.9
   ```

2. **Create and activate a conda environment**

   ```bash
   conda create -n pdb2reaction python=3.11 -y
   conda activate pdb2reaction
   ```

3. **Install cyipopt**
   Required if you want to use the DMF method in MEP search.

   ```bash
   conda install -c conda-forge cyipopt -y
   ```

4. **Install PyTorch with the right CUDA build**

   For CUDA 12.9:

   ```bash
   pip install torch --index-url https://download.pytorch.org/whl/cu129
   ```

   (You may use another compatible version if your cluster recommends it.)

5. **Install `pdb2reaction` itself and Chrome for visualization**

   ```bash
   pip install git+https://github.com/t-0hmura/pdb2reaction.git
   plotly_get_chrome -y
   ```

6. **Log in to Hugging Face Hub (UMA model)**

   ```bash
   huggingface-cli login
   ```

   See also:

   - <https://github.com/facebookresearch/fairchem>
   - <https://huggingface.co/facebook/UMA>
   - <https://huggingface.co/docs/hub/security-tokens>

---

## Command line basics

The main entry point is the `pdb2reaction` command, installed via `pip`. Internally it uses the **Click** library, and the default subcommand is `all`.

That means:

```bash
pdb2reaction [OPTIONS] ...
# is equivalent to
pdb2reaction all [OPTIONS] ...
```

The `all` workflow acts as an **orchestrator**: it chains cluster extraction, MEP search, TS optimization, vibrational analysis, and optional DFT single points into a single command.

All high‑level workflows share two important options when you want cluster extraction:

- `-i/--input`: one or more **full structures** (reactant, intermediate(s), product).
- `-c/--center`: how to define the **substrate / extraction center** (e.g., residue names or residue IDs).

If you omit `--center/-c`, cluster extraction is skipped and the **full input structure** is used directly.

---

## Main workflow modes

### Multi‑structure MEP pipeline (reactant → product)

Use this when you already have several full PDB structures along a putative reaction coordinate (e.g., R → I1 → I2 → P).

**Minimal example**

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

**Richer example**

```bash
pdb2reaction -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --out-dir ./result_all --tsopt True --thermo True --dft True
```

Behavior:

- takes two or more **full systems** in reaction order,
- extracts catalytic cluster models for each structure,
- performs a **recursive MEP search** via `path_search` by default,
- optionally switches to a **single‑pass** `path-opt` run with `--refine-path False`,
- when PDB templates are available, merges the cluster-model MEP back into the **full system**,
- optionally runs TS optimization, vibrational analysis, and DFT single points for each segment.

This is the recommended mode when you can generate reasonably spaced intermediates (e.g., from docking, MD, or manual modeling).

```{important}
`pdb2reaction` assumes that multiple input PDBs contain **exactly the same atoms in the same order** (only coordinates may differ). If any non-coordinate fields differ across inputs, an error is raised. Input PDB files must also contain **hydrogen atoms**.
```

---

### Single‑structure + staged scan (feeds MEP refinement)

Use this when you only have **one PDB structure**, but you know which inter‑atomic distances should change along the reaction.

Provide a single `-i` together with `--scan-lists`:

**Minimal example**

```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]'
```

**Richer example**

```bash
pdb2reaction -i SINGLE.pdb -c 'SAM,GPP' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]' --mult 1 --out-dir ./result_scan_all --tsopt True --thermo True --dft True
```

Key points:

- `--scan-lists` describes **staged distance scans** on the extracted cluster model.
- Each tuple `(i, j, target_Å)` is:
  - a PDB atom selector string like `'TYR,285,CA'` (**delimiters can be: space/comma/slash/backtick/backslash ` ` `,` `/` `` ` `` `\`**) **or** a 1‑based atom index,
  - automatically remapped to the cluster-model indices.
- Supplying one `--scan-lists` literal runs a single scan stage; multiple literals run sequential stages. Pass multiple literals after a single flag (repeated flags are not accepted).
- Each stage writes a `stage_XX/result.pdb`, which is treated as a candidate intermediate or product.
- The default `all` workflow refines the concatenated stages with recursive `path_search`.
- With `--refine-path False`, it instead performs a single-pass `path-opt` chain and skips the recursive refiner (no merged `mep_w_ref*.pdb`).

This mode is useful for building reaction paths starting from a single structure.

---

### Single‑structure TSOPT‑only mode

Use this when you already have a **transition state candidate** and only want to refine it and proceed to IRC calculations.

Provide exactly one PDB and enable `--tsopt`:

**Minimal example**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True
```

**Richer example**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt True --thermo True --dft True --out-dir ./result_tsopt_only
```

Behavior:

- skips the MEP/path search entirely,
- refines the **cluster-model TS** with TS optimization,
- runs an **IRC** in both directions and optimizes both ends to relax down to R and P minima,
- can then perform `freq` and `dft` on the R/TS/P,
- produces UMA, Gibbs, and DFT//UMA energy diagrams.

Outputs such as `energy_diagram_*_all.png` and `irc_plot_all.png` are mirrored under the top‑level `--out-dir`.

```{important}
Single‑input runs require **either** `--scan-lists` (staged scan → GSM) **or** `--tsopt True` (TSOPT‑only). Supplying only a single `-i` without one of these will not trigger a full workflow.
```

---

## Important CLI options and behaviors

Below are the most commonly used options across workflows.

| Option | Description |
|--------|-------------|
| `-i, --input PATH...` | Input structures. **≥ 2 PDBs** → MEP search; **1 PDB + `--scan-lists`** → staged scan → GSM; **1 PDB + `--tsopt True`** → TSOPT‑only mode. |
| `-c, --center TEXT` | Defines the substrate / extraction center. Supports residue names (`'SAM,GPP'`), residue IDs (`A:123,B:456`), or PDB paths. |
| `--ligand-charge TEXT` | Charge info: mapping (`'SAM:1,GPP:-3'`) or single integer. |
| `-q, --charge INT` | Hard override of total system charge. |
| `-m, --mult INT` | Spin multiplicity (e.g., `1` for singlet). Note: Use `--multiplicity` in subcommands other than `all`. |
| `--scan-lists TEXT...` | Staged distance scans for single‑input runs. |
| `--out-dir PATH` | Top‑level output directory. |
| `--tsopt {True\|False}` | Enable TS optimization and IRC. |
| `--thermo {True\|False}` | Run vibrational analysis and thermochemistry. |
| `--dft {True\|False}` | Perform DFT single‑point calculations. |
| `--refine-path {True\|False}` | Recursive MEP refinement (default) vs single‑pass. |
| `--opt-mode light\|heavy` | Optimization method: Light (LBFGS/Dimer) or Heavy (RFO/RS-I-RFO). |
| `--mep-mode gsm\|dmf` | MEP method: Growing String Method or Direct Max Flux. |
| `--hessian-calc-mode Analytical\|FiniteDifference` | Hessian calculation mode. **Analytical recommended when VRAM available.** |

For a full matrix of options and YAML schemas, see [all](all.md).

---

## Run summaries

Every `pdb2reaction all` run writes:

- `summary.log` – formatted summary for quick inspection, and
- `summary.yaml` – YAML version summary.

They typically contain:

- the exact CLI command invoked,
- global MEP statistics (e.g. maximum barrier, path length),
- per‑segment barrier heights and key bond changes,
- energies from UMA, thermochemistry, and DFT post‑processing (where enabled).

Each `path_search` segment directory also gets its own `summary.log` and `summary.yaml`, so you can inspect local refinements independently.

---

## CLI commands

While most users will primarily call `pdb2reaction all`, the CLI also exposes subcommands like `pdb2reaction opt`. Each subcommand supports `-h/--help`.

| Subcommand | Role | Documentation |
|------------|------|---------------|
| `all` | End-to-end workflow | [all](all.md) |
| `extract` | Extract active-site cluster | [extract](extract.md) |
| `opt` | Geometry optimization | [opt](opt.md) |
| `tsopt` | Transition state optimization | [tsopt](tsopt.md) |
| `path-opt` | MEP optimization (GSM/DMF) | [path_opt](path_opt.md) |
| `path-search` | Recursive MEP search | [path_search](path_search.md) |
| `scan` | 1D bond-length scan | [scan](scan.md) |
| `scan2d` | 2D distance scan | [scan2d](scan2d.md) |
| `scan3d` | 3D distance scan | [scan3d](scan3d.md) |
| `irc` | IRC calculation | [irc](irc.md) |
| `freq` | Vibrational analysis | [freq](freq.md) |
| `dft` | DFT single-point | [dft](dft.md) |
| `trj2fig` | Plot energy profiles | [trj2fig](trj2fig.md) |
| `add-elem-info` | Repair PDB element columns | [add_elem_info](add_elem_info.md) |

```{important}
The subcommands (everything **except** `all`) assume you feed them **cluster models** generated by `extract`. In these cluster models, the atom closest to the Link‑H cap is automatically **frozen**. If you construct a cluster model yourself, set the Link‑H residue name to `LKH` and the atom name to `HL`, or specify the atoms to freeze via `--args-yaml` → `geom.freeze_atoms`.
```

```{tip}
In `all`, `tsopt`, `freq` and `irc`, setting **`--hessian-calc-mode Analytical`** is strongly recommended when you have enough VRAM.
```

---

## Getting help

For any subcommand:

```bash
pdb2reaction <subcommand> --help
```

This prints the available options, defaults, and a short description. For detailed UMA calculator options, see [UMA Calculator](uma_pysis.md).

If you encounter any issues, please open an Issue on the [GitHub repository](https://github.com/t-0hmura/pdb2reaction).
