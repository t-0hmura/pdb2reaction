# pdb2reaction: automated reaction-path modeling directly from PDB structures

## Overview

`pdb2reaction` is a Python CLI toolkit for turning **PDB structures** into **enzymatic reaction pathways** with Machine Learning Interatomic Potential (MLIP).  
  
Basically, you just need a **single command** such as  
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```  
for modeling enzymatic reaction pathways.  
  
Furthermore, you can proceed MEP search --> TS refinement --> IRC --> thermochemistry analysis --> DFT single point calculation with **single command** such as  
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' **--tsopt True --thermo True --dft True**
```

---

Given **(i) two or more full protein–ligand PDBs** `.pdb` (R → … → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt True`**, `pdb2reaction` automatically:

- extracts an **active site** around user‑defined substrates to build a **cluster model**,
- explores **minimum‑energy paths (MEPs)** with path optimization methods such as the Growing String Method (GSM) and Direct Max Flux (DMF),
- _optionally_ refines **transition states**, runs **vibrational analysis**, **IRC calculations** and **DFT single‑point calculations**

with a **machine learning interatomic potential (MLIP)** using Meta’s UMA model.
   
All of this is exposed through a command‑line interface (CLI) designed so that a **multi‑step enzymatic reaction mechanism** can be generated with minimal manual intervention. Of course, this toolkit can handle small molecular systems. You can also use `.xyz` or `.gjf` format input structures when you run workflows on the full structure (i.e., omit `--center/-c` and `--ligand-charge`).

On **HPC clusters or multi‑GPU workstations**, `pdb2reaction` can process large cluster models (and optionally **full protein–ligand complexes**) by parallelizing UMA inference itself across nodes. Set `workers` and `workers_per_node` to enable parallel calculation; see [`docs/uma_pysis.md`](docs/uma_pysis.md) for configuration details.

> **Important (prerequisites):**  
> - Input PDB files must already contain **hydrogen atoms**.  
> - When you provide multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ); otherwise an error is raised.  
> - Boolean CLI options are passed explicitly as `True`/`False` (e.g., `--tsopt True`).

***This software is still under development. Please use it at your own risk.***

---

## 1. Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters) with a CUDA‑capable GPU. Several dependencies – notably **PyTorch**, **fairchem‑core (UMA)**, and **gpu4pyscf‑cuda12x** – expect a working CUDA installation.

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

### 1.1 Quick start

Below is a minimal setup example that works on many CUDA 12.9 clusters. Adjust module names and versions to match your system. Below assumes the default GSM MEP mode (no DMF). For DMF, install cyipopt via conda first.

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

- If you want to use Direct Max Flux method for MEP search, create conda environment and install cyipopt before installation.
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


### 1.2 Step‑by‑step installation

If you prefer to build the environment piece by piece:

1. **Load CUDA (when you use environment modules on HPC)**

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

## 2. Command line basics

The main entry point is the `pdb2reaction` command, installed via `pip`. Internally it uses the **Click** library, and the default subcommand is `all`.

That means:

```bash
pdb2reaction [OPTIONS] ...
# is equivalent to
pdb2reaction all [OPTIONS] ...
```

The `all` workflow is an **orchestrator**: it chains cluster extraction, MEP search, TS optimization, vibrational analysis, and optional DFT single points into a single command.

All high‑level workflows share two important options when you want cluster extraction:

- `-i/--input`: one or more **full structures** (reactant, intermediate(s), product).
- `-c/--center`: how to define the **substrate / extraction center** (e.g., residue names or residue IDs).

If you omit `--center/-c`, cluster extraction is skipped and the **full input structure** is used directly.

---

## 3. Main workflow modes

### 3.1 Multi‑structure MEP pipeline (reactant → product)

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

> **Important:** `pdb2reaction` assumes that multiple input PDBs contain **exactly the same atoms in the same order** (only coordinates may differ). If any non-coordinate fields differ across inputs, an error is raised. Input PDB files must also contain **hydrogen atoms**.  
---

### 3.2 Single‑structure + staged scan (feeds MEP refinement)

Use this when you only have **one PDB structure**, but you know which inter‑atomic distances should change along the reaction.

Provide a single `-i` together with `--scan-lists`:

**Minimal example**

```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]'
```

**Richer example**

```bash
pdb2reaction -i SINGLE.pdb -c 'SAM,GPP' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' --mult 1 --out-dir ./result_scan_all --tsopt True --thermo True --dft True
```

Key points:

- `--scan-lists` describes **staged distance scans** on the extracted cluster model.
- Each tuple `(i, j, target_Å)` is:
  - a PDB atom selector string like `'TYR,285,CA'` (**delimiters can be: space/comma/slash/backtick/backslash `space` `,` `/` `` ` `` `\`**) **or** a 1‑based atom index,  
  - automatically remapped to the cluster-model indices.
- Each stage writes a `stage_XX/result.pdb`, which is treated as a candidate intermediate or product.
- The default `all` workflow refines the concatenated stages with recursive `path_search`.
- With `--refine-path False`, it instead performs a single-pass `path-opt` chain and skips the recursive refiner (no merged `mep_w_ref*.pdb`).

This mode is useful for building approximate reaction paths starting from a single experimental structure.

---

### 3.3 Single‑structure TSOPT‑only mode

Use this when you already have a **transition state candidate** and only want to refine it, without running a full path search.

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

> **Important:** Single‑input runs require **either** `--scan-lists` (staged scan → GSM) **or** `--tsopt True` (TSOPT‑only). Supplying only a single `-i` without one of these will not trigger a full workflow.

---

## 4. Important CLI options and behaviors

Below are the most commonly used options across workflows.

- `-i, --input PATH...`  
  Input structures. Interpretation depends on how many you provide:

  - **≥ 2 PDBs** → MEP search (GSM by default, DMF with `--mep-mode dmf`) (reactant/intermediates/product).
  - **1 PDB + `--scan-lists`** → staged scan → GSM.
  - **1 PDB + `--tsopt True`** → TSOPT‑only mode.

  If `--center/-c` is omitted, cluster extraction is skipped and the **full input structure** is used directly. In this mode, `.xyz` and `.gjf` inputs are also accepted; when using them, omit `--center/-c` and `--ligand-charge`.

- `-c, --center TEXT`  
  Defines the substrate / extraction center for cluster extraction. Supports:

  - PDB paths (the atoms in the given PDB define the center structure used during cluster extraction),
  - residue IDs, e.g. `A:123,B:456`,
  - residue names, e.g. `'SAM,GPP'`.

- `--ligand-charge TEXT`  
  Total charge information, either:

  - a single integer (total cluster-model charge), or
  - a mapping, e.g. `'SAM:1,GPP:-3'`.

  The total charge of the first cluster model is summed and used for scan, GSM, and TS optimization runs.

- `-q, --charge INT`  
  Hard override of the total system charge. This bypasses:

  - extractor rounding,
  - `.gjf` metadata,
  - and `--ligand-charge` resolution.

  Use when you want full manual control of the charge.

- `--mult INT`  
  Spin multiplicity for QM regions (e.g., `--mult 1` for singlet). Used for scan and GSM runs.

- `--scan-lists TEXT...`  
  One or more Python‑style lists describing **staged scans** for single‑input runs. Example:

  ```bash
  --scan-lists '[(10,55,2.20),(23,34,1.80)]'
  ```

  Each tuple describes a harmonic distance restraint between atoms `i` and `j` driven to a target in Å. Indices are 1‑based in the original full PDB and are automatically remapped onto the cluster model.

- `--out-dir PATH`
  Top‑level output directory. All intermediate files, logs, and figures are placed here.

- `--tsopt BOOLEAN`  
  Enable TS optimization and IRC propagation. Required for TSOPT‑only mode, but also useful in multi‑structure workflows to refine TS along the path.

- `--thermo BOOLEAN`  
  Run vibrational analysis and compute thermochemistry on UMA geometries. Produces Gibbs free energies and corresponding energy diagrams.

- `--dft BOOLEAN`  
  Perform DFT single‑point calculations on UMA optimized structures via PySCF / gpu4pyscf. When combined with `--thermo True`, this adds DFT//UMA Gibbs diagrams.

- `--refine-path BOOLEAN`
  Switch between:
  - **recursive MEP refinement** with `path_search` (default),
  - **single‑pass MEP refinement** with `path-opt` (simple MEP) when set to `False`.

  When `--refine-path True` (default) and full‑system PDB templates are available, merged MEP snapshots (`mep_w_ref*.pdb`) are written under `<out-dir>/path_search/`.

- `--opt-mode light | heavy`
  Switch optimization / TS refinement methods between Light (LBFGS and Dimer) and Heavy (Hessian-using RFO and RS-I-RFO) algorithms.

- `--mep-mode gsm | dmf`
  Switch MEP refinement between the Growing String Method (GSM) and Direct Max Flux (DMF).

For a full matrix of options and YAML schemas, see `docs/all.md` in the repository.

---

## 5. Run summaries (`summary.log` and `summary.yaml`)

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

## 6. CLI subcommands

While most users will primarily call `pdb2reaction all`, the CLI also exposes subcommands like `pdb2reaction opt`. Each subcommand supports `-h/--help` (see `docs/*.md` for exact schemas).

| Subcommand   | Role                                                                                  | Documentation                             |
| ------------ | ------------------------------------------------------------------------------------- | ----------------------------------------- |
| `all`        | Extraction → (optional staged scan) → MEP search → merge to full PDBs in one shot. (Optionally, following TS opt, IRC, Thermochemistry analysis, DFT single point calculation) | [`docs/all.md`](docs/all.md)              |
| `scan`       | Bond‑length driven scan with staged harmonic restraints and relaxation.               | [`docs/scan.md`](docs/scan.md)            |
| `opt`        | Single‑structure geometry optimization using LBFGS or RFO.                            | [`docs/opt.md`](docs/opt.md)              |
| `path-opt`   | MEP optimization via the GSM or DMF.                                                  | [`docs/path_opt.md`](docs/path_opt.md)    |
| `path-search`| Multistep MEP search via recursive GSM or DMF segmentation.                           | [`docs/path_search.md`](docs/path_search.md) |
| `tsopt`      | Transition state optimization (Dimer or RS-I-RFO).                                    | [`docs/tsopt.md`](docs/tsopt.md)          |
| `freq`       | Vibrational frequency analysis and mode writer (+ thermochemistry summary).           | [`docs/freq.md`](docs/freq.md)            |
| `irc`        | IRC calculation with EulerPC.                                                         | [`docs/irc.md`](docs/irc.md)              |
| `extract`    | Extract an active site structure around substrate residues and build cluster model.   | [`docs/extract.md`](docs/extract.md)      |
| `trj2fig`    | Plot ΔE or E from an XYZ trajectory and export figure/CSV.                            | [`docs/trj2fig.md`](docs/trj2fig.md)      |
| `add-elem-info` | Add or repair PDB element columns (77–78) using Biopython.                         | [`docs/add_elem_info.md`](docs/add_elem_info.md) |
| `dft`        | Single-point DFT using GPU4PySCF (with CPU PySCF fallback).                           | [`docs/dft.md`](docs/dft.md)              |
| `scan2d`     | 2D distance scan with harmonic restraints.                                            | [`docs/scan2d.md`](docs/scan2d.md)        |
| `scan3d`     | 3D distance scan with harmonic restraints.                                            | [`docs/scan3d.md`](docs/scan3d.md)        |

> **Important:** The subcommands (everything **except** `all`) assume you feed them **cluster models** generated by `extract`. In these cluster models, the atom closest to the Link‑H cap is automatically **frozen**. If you construct a cluster model yourself, set the Link‑H residue name to `LKH` and the atom name to `HL`, or specify the atoms to freeze via `--args-yaml` → `geom.freeze_atoms`.

---

## 7. Getting help

For any subcommand:

```bash
pdb2reaction <subcommand> --help
```

This prints the available options, defaults, and a short description. For detailed workflows, argument schemas, and example YAML files, consult the `docs/*.md` files in the repository (e.g. `docs/all.md`, `docs/scan.md`). For detailed UMA calculator options and defaults, see [`docs/uma_pysis.md`](docs/uma_pysis.md).

> If you encounter any issues—such as dependency conflicts, library bugs, or other uncertainties—please open an Issue on the repository. I will address them as promptly as possible.

## Citation
A preprint describing `pdb2reaction` will be released soon. Please check back for citation details once it is available.

## References
[1] Wood, B. M., Dzamba, M., Fu, X., Gao, M., Shuaibi, M., Barroso-Luque, L., Abdelmaqsoud, K., Gharakhanyan, V., Kitchin, J. R., Levine, D. S., Michel, K., Sriram, A., Cohen, T., Das, A., Rizvi, A., Sahoo, S. J., Ulissi, Z. W., & Zitnick, C. L. (2025). UMA: A Family of Universal Models for Atoms. http://arxiv.org/abs/2506.23971   
[2] Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. International Journal of Quantum Chemistry, 121(3). https://doi.org/10.1002/qua.26390

## License
`pdb2reaction` is distributed under the **GNU General Public License version 3 (GPL-3.0)** derived from Pysisyphus.  
