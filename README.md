# pdb2reaction: automated reaction-path modelling directly from PDB structures

## Overview

`pdb2reaction` is a Python CLI toolkit for turning **PDB structures** into **enzymatic reaction pathways** with Machine Learning Interatomic Potential (MLIP). 
Basically, you just need a **single command** such as `pdb2reaction -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3"`

---

Given one or more full protein–ligand PDBs (reactant, intermediates, product), it automatically:

- extracts a **catalytic pocket** around user‑defined substrates,
- builds a **machine‑learned potential** using Meta’s UMA model,
- explores **minimum‑energy paths (MEPs)** with growing string / path optimisation methods,
- refines **transition states**, runs **vibrational analysis**, and
- optionally performs **DFT single‑point calculations** for higher‑level energetics.

All of this is exposed through a command‑line interface (CLI) designed so that a **multi‑step enzymatic mechanism** can be generated with minimal manual intervention.

Typical use cases:

- mapping out plausible reaction mechanisms from a handful of crystal structures or snapshots (with hydrogen atoms),
- refining TS candidates obtained from other tools,
- scanning key distances to generate intermediates,
- building Gibbs and DFT//ML energy diagrams for publication.

---

## 1. Installation

`pdb2reaction` is intended for Linux‑like environments (local workstations or HPC clusters) with a CUDA‑capable GPU. Several dependencies – notably **PyTorch**, **fairchem‑core (UMA)**, and **gpu4pyscf‑cuda12x** – expect a working CUDA installation.

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

### 1.1 Quick start

Below is a minimal setup that works on many CUDA 12.8 clusters. Adjust module names and versions to match your system.

```bash
# 1) Install a CUDA-enabled PyTorch build
# 2) Install pdb2reaction from GitHub
# 3) Install a headless Chrome for Plotly figure export

pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu128
pip install git+https://github.com/t-0hmura/pdb2reaction.git
plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded. Either:

```bash
# Hugging Face CLI
hf auth login --token "<YOUR_ACCESS_TOKEN>" --add-to-git-credential
```

or

```bash
# Classic CLI
huggingface-cli login
```

You only need to do this once per machine / environment.

- If you want to use Direct Max flux method for MEP search, create conda environment and install cyiopt before installation.
  ```bash
  # Create and activate a dedicated conda environment
  conda create -n pdb2reaction python=3.11 -y
  conda activate pdb2reaction

  # Install cyipopt (required for the DMF method in MEP search)
  conda install -c conda-forge cyipopt -y
  ```
  
- If you are on an HPC cluster that uses *environment modules*, load CUDA **before** installing PyTorch, like this:
  ```bash
  module load cuda/12.8
  ```


### 1.2 Step‑by‑step installation

If you prefer to build the environment piece by piece:

1. **Load CUDA (HPC modules)**

   ```bash
   module load cuda/12.8
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

   For CUDA 12.8:

   ```bash
   pip install torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128
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

The `all` workflow is an **orchestrator**: it chains pocket extraction, MEP search, TS optimisation, vibrational analysis, and optional DFT single points into a single command.

All high‑level workflows share two important options:

- `-i/--input`: one or more **full PDB structures** (reactant, intermediate(s), product).
- `-c/--center`: how to define the **substrate / pocket center** (e.g., residue names or residue IDs).

Unless you intentionally skip extraction, you must supply both.

---

## 3. Main workflow modes

### 3.1 Multi‑structure GSM pipeline (reactant → product)

Use this when you already have several full PDB structures along a putative reaction coordinate (e.g., R → I1 → I2 → P).

**Minimal example**

```bash
pdb2reaction -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3"
```

**Richer example**

```bash
pdb2reaction -i R.pdb I1.pdb I2.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --outdir ./result_all --tsopt True --thermo True --dft True
```

Behaviour:

- takes two or more **full systems** in reaction order,
- extracts catalytic pockets for each structure,
- performs a **single‑pass GSM / MEP search** via `path-opt` by default,
- optionally switches to a recursive refiner (`path_search`) with `--refine-path True`,
- when PDB templates are available, merges the pocket‑level MEP back into the **full system**,
- optionally runs TS optimisation, vibrational analysis, and DFT single points for each segment.

This is the recommended mode when you can generate reasonably spaced intermediates (e.g., from docking, MD, or manual modelling).

> You need to input structures in **same atomic order**.  
> And PDB structures must have **hydrogen atoms**.  
---

### 3.2 Single‑structure + staged scan (feeds GSM)

Use this when you only have **one PDB structure**, but you know which inter‑atomic distances should change along the reaction.

Provide a single `-i` together with `--scan-lists`:

**Minimal example**

```bash
pdb2reaction -i R.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --scan-lists "[(10,55,2.20),(23,34,1.80)]"
```

**Richer example**

```bash
pdb2reaction -i SINGLE.pdb -c "SAM,GPP" --scan-lists "[(10,55,2.20),(23,34,1.80)]" --mult 1 --outdir ./result_scan_all --tsopt True --thermo True --dft True
```

Key points:

- `--scan-lists` describes **staged distance scans** on the extracted pocket, using UMA.
- Each tuple `(i, j, target_Å)` is:
  - 1‑based indices taken from the original full PDB,
  - automatically remapped to the pocket indices.
- Each stage writes a `stage_XX/result.pdb`, which is treated as a candidate intermediate or product.
- The default `all` workflow then concatenates these stages using `path-opt`.
- With `--refine-path True`, it instead runs the recursive `path_search` refiner and (when possible) writes merged full‑system `mep_w_ref*.pdb` files under `<outdir>/path_search/`.

This mode is useful for building approximate reaction paths starting from a single experimental structure.

---

### 3.3 Single‑structure TSOPT‑only mode

Use this when you already have a **transition state candidate** and only want to refine it, without running a full path search.

Provide exactly one PDB and enable `--tsopt`:

**Minimal example**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3"
```

**Richer example**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --tsopt True --thermo True --dft True --outdir ./result_tsopt_only
```

Behaviour:

- skips the MEP/path search entirely,
- refines the **pocket TS** with UMA‑based TS optimisation,
- runs a **pseudo‑IRC** in both directions to relax down to R and P minima,
- can then perform `freq` and `dft` on the R/TS/P trio,
- produces UMA, Gibbs, and DFT//UMA energy diagrams, similar to the full `all` workflow.

Outputs such as `energy_diagram_*_all.png` and `irc_plot_all.png` are mirrored under the top‑level `--outdir`.

> **Important:** Single‑input runs require **either** `--scan-lists` (staged scan → GSM) **or** `--tsopt True` (TSOPT‑only). Supplying only a single `-i` without one of these will not trigger a full workflow.

---

## 4. Important CLI options and behaviours

Below are the most commonly used options across workflows.

- `-i, --input PATH...`  
  Input structures. Interpretation depends on how many you provide:

  - **≥ 2 PDBs** → GSM mode (reactant/intermediates/product).
  - **1 PDB + `--scan-lists`** → staged scan → GSM.
  - **1 PDB + `--tsopt True`** → TSOPT‑only mode.

  When pocket extraction is skipped, XYZ/GJF inputs are also accepted.

- `-c, --center TEXT`  
  Defines the substrate / pocket centre for extraction. Supports:

  - PDB paths,
  - residue IDs, e.g. `A:123,B:456`,
  - residue names, e.g. `"SAM,GPP"`.

- `--ligand-charge TEXT`  
  Total charge information, either:

  - a single integer (total pocket charge), or
  - a mapping, e.g. `"SAM:1,GPP:-3"`.

  The total charge of the first pocket is rounded to an integer and reused for scan, GSM, and TS optimisation runs.

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
  --scan-lists "[(10,55,2.20),(23,34,1.80)]"
  ```

  Each tuple describes a harmonic distance restraint between atoms `i` and `j` driven to a target in Å. Indices are 1‑based in the original full PDB and are automatically remapped onto the pocket.

- `--outdir PATH`  
  Top‑level output directory. All intermediate files, logs, and figures are placed here.

- `--tsopt BOOLEAN`  
  Enable TS optimisation and pseudo‑IRC propagation. Required for TSOPT‑only mode, but also useful in multi‑structure workflows to refine TS along the path.

- `--thermo BOOLEAN`  
  Run vibrational analysis and compute thermochemistry on UMA geometries. Produces Gibbs free energies and corresponding energy diagrams.

- `--dft BOOLEAN`  
  Perform DFT single‑point calculations on UMA optimised structures via PySCF / gpu4pyscf. When combined with `--thermo True`, this adds DFT//UMA Gibbs diagrams.

- `--refine-path BOOLEAN`  
  Switch between:

  - **single‑pass GSM** with `path-opt` (simple MEP),
  - **recursive GSM / MEP refinement** with `path_search`.

  When `--refine-path True` and full‑system PDB templates are available, merged MEP snapshots (`mep_w_ref*.pdb`) are written under `<outdir>/path_search/`.

For a full matrix of options and YAML schemas, see `docs/all.md` in the repository.

---

## 5. Run summaries (`summary.log` and `summary.yaml`)

Every `pdb2reaction all` run writes a human‑readable:

- `summary.log` – formatted for quick inspection, and
- `summary.yaml` – YAML version of the same information.

They typically contain:

- the exact CLI command invoked,
- global MEP statistics (e.g. maximum barrier, path length),
- per‑segment barrier heights and key bond changes,
- energies from UMA, thermochemistry, and DFT post‑processing (where enabled).

Each `path_search` segment directory also gets its own `summary.log` and `summary.yaml`, so you can inspect local refinements independently.

---

## 6. CLI subcommands

While most users will primarily call `pdb2reaction all`, the CLI also exposes lower‑level building blocks like `pdb2reaction opt`. Each subcommand supports `-h/--help` and can read arguments from YAML files via `--args-yaml` (see `docs/*.md` for exact schemas).

| Subcommand   | Role                | Documentation            |
| ------------ | ---------------------------------------------------------------------------- | ------------------------ |
| `all`        | High‑level workflow orchestrator: extraction → GSM / path search → TS/freq/DFT. | `docs/all.md`           |
| `scan`       | Staged biased scans on pocket models to generate additional intermediates.   | `docs/scan.md`          |
| `opt`        | Optimise a single structure (LBFGS / RFO presets).       | `docs/opt.md`           |
| `path-opt`   | UMA optimisation on a specific path segment or snapshot.        | `docs/path_opt.md`      |
| `path-search`| Recursive GSM‑based path search plus pocket/full‑system merging.             | `docs/path_search.md`   |
| `tsopt`      | Transition‑state refinement (Dimer+LBFGS or RS-I-RFO presets).            | `docs/tsopt.md`         |
| `freq`       | Vibrational analysis, thermochemistry.    | `docs/freq.md`          |
| `irc`        | Intrinsic reaction coordinate following from a TS structure.    | `docs/irc.md`           |
| `extract`    | Extract active sites from full PDB structures.     | `docs/extract.md`       |
| `trj2fig`    | Convert trajectory data to interactive figures (Plotly / Kaleido).           | `docs/trj2fig.md`       |
| `add-elem-info` | Add missing element metadata to PDB files.     | `docs/add_elem_info.md` |
| `dft`        | DFT single‑point calculations on UMA geometries using PySCF / gpu4pyscf.     | `docs/dft.md`           |
| `scan2d`     | 2D PES grid: two distance scan with restraints and relaxation.       | `docs/scan2d.md`        |
| `scan3d`     | 3D grid over three distances; can also visualise existing `surface.csv` data.| `docs/scan3d.md`        |

In practice, you can:

- prototype with lower‑level subcommands (`scan`, `scan2d`, `opt`),
- then wrap everything into a reproducible `pdb2reaction all` run for production.

---

## 7. Getting help

For any subcommand:

```bash
pdb2reaction <subcommand> --help
```

This prints the available options, defaults, and a short description. For detailed workflows, argument schemas, and example YAML files, consult the `docs/*.md` files in the repository (e.g. `docs/all.md`, `docs/scan.md`).
