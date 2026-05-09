# Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters) with a CUDA-capable GPU. Several dependencies – notably **PyTorch**, **fairchem-core (UMA)**, and **gpu4pyscf-cuda12x** – expect a working CUDA installation.

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

## Quick start

Below is a minimal setup example that works on many CUDA 12.9 clusters. Adjust module names and versions to match your system. This example assumes the default GSM MEP mode (`--mep-mode gsm`). For DMF (`--mep-mode dmf`), install cyipopt via conda first.

### Required

```bash
# 1) Install a CUDA-enabled PyTorch build
# 2) Install pdb2reaction
# 3) Install headless Chrome for Plotly static image export (PNG)
#    Downloads ~150 MB Chromium binary; requires internet access.

pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install pdb2reaction
plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded (requires a free HF account with read-only token; you may need to accept the UMA model license at <https://huggingface.co/facebook/UMA>):

```bash
hf auth login
# or, with an access token in scripts:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

(Requires `huggingface_hub >= 0.34`; in `huggingface_hub >= 1.0` the legacy `huggingface-cli login` entry point exits with a deprecation notice.)

You only need to do this once per machine / environment.

### Optional

- If you want to use the Direct Max Flux (DMF) method for MEP search, create a conda environment and install cyipopt before installing pdb2reaction.

  ```bash
  # Create and activate a dedicated conda environment
  conda create -n <your-env> python=3.11 -y
  conda activate <your-env>

  # Install cyipopt (required for the DMF method in MEP search)
  conda install -c conda-forge cyipopt -y
  ```

- If you are on an HPC cluster that uses *environment modules*, load CUDA **before** installing PyTorch. Run `module avail cuda` to see which CUDA versions your site provides, then load the one matching your target PyTorch wheel (e.g. `cu126` ↔ CUDA 12.6, `cu129` ↔ CUDA 12.9):

  ```bash
  module load cuda/<your-version>   # e.g. cuda/12.6 or cuda/12.9
  ```

> **Tip:** UMA is the default MLIP backend. To use ORB or AIMNet2, install the corresponding extra (e.g. `pip install "pdb2reaction[orb]"`) and pass `-b/--backend orb` to any command. See step 7 below.

```{warning}
**MACE:** `mace-torch` requires `e3nn==0.4.4`, which conflicts with `fairchem-core`'s `e3nn` pin (older `mace-torch` versions). The canonical recipe is `pip uninstall -y fairchem-core && pip install mace-torch` in a dedicated conda env — UMA and MACE cannot share the same environment unless `mace-torch` ≥ 0.3.8 (which advertises coexistence with `fairchem-core`). The `--no-deps mace-torch` variant seen in some older notes is not recommended; it leaves torch-scatter / e3nn unpinned.
```


(step-by-step-installation)=
## Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load CUDA (if you use environment modules on an HPC cluster)**

    Run `module avail cuda` to see what is provided, then load the
    version matching your target PyTorch wheel (e.g. `cu126` for CUDA
    12.6, `cu129` for CUDA 12.9):

    ```bash
    module load cuda/<your-version>
    ```

2. **Create and activate a conda environment**

    ```bash
    conda create -n <your-env> python=3.11 -y
    conda activate <your-env>
    ```

3. **Install cyipopt**
    Required if you want to use the DMF method (`--mep-mode dmf`) in MEP search. You can skip this step if you only use GSM.

    ```bash
    conda install -c conda-forge cyipopt -y
    ```

4. **Install PyTorch with the right CUDA build**

    For CUDA 12.9:

    ```bash
    pip install torch --index-url https://download.pytorch.org/whl/cu129
    ```

    PyTorch must be built for your CUDA driver version. Check compatibility at [PyTorch Get Started](https://pytorch.org/get-started/locally/). CPU-only execution is supported but significantly slower (10-100x).

5. **Install `pdb2reaction` itself and Chrome for visualization**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

6. **Log in to Hugging Face Hub (UMA model)**

    ```bash
    hf auth login
    ```

    See also:

    - <https://github.com/facebookresearch/fairchem>
    - <https://huggingface.co/facebook/UMA>
    - <https://huggingface.co/docs/hub/security-tokens>

7. **(Optional) Install additional MLIP backends**

    pdb2reaction uses UMA by default. To use alternative backends, install the corresponding optional dependency:

    ```bash
    # ORB backend
    pip install "pdb2reaction[orb]"

    # AIMNet2 backend
    pip install "pdb2reaction[aimnet]"

    # MACE backend (use a separate conda environment because mace-torch
    # pins e3nn==0.4.4 which conflicts with UMA's fairchem-core)
    pip uninstall -y fairchem-core && pip install mace-torch

    # DFT single-point post-processing (`--dft` / `pdb2reaction dft`)
    # Installs gpu4pyscf-cuda12x, PySCF, and related dependencies.
    # Note: gpu4pyscf-cuda12x publishes x86_64 wheels on PyPI; on
    # aarch64 build from source (https://github.com/pyscf/gpu4pyscf).
    pip install "pdb2reaction[dft]"
    ```

    To enable implicit solvent corrections, install [xTB](https://github.com/grimme-lab/xtb) and ensure the `xtb` command is available on your `PATH`.

    #### Installing xTB

    **For ALPB solvation model** (recommended starting point):

    ```bash
    conda install -c conda-forge xtb
    ```

    **For CPCM-X solvation model** the conda-forge xtb does not include CPCM-X; build from source. See {ref}`recipe-cpcmx-build` for the recipe.

    To use a custom xTB binary, set the `xtb_cmd` key in your YAML config or use `calc.xtb_cmd` in Python.

8. **Verify installation**

    ```bash
    pdb2reaction --version
    ```

    This should display the installed version. To verify GPU access:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    If `CUDA: False`, check that the correct CUDA module is loaded and the PyTorch build matches your CUDA driver version.

## Next steps

- [Getting Started](getting-started.md) — project overview, pipeline stages, and workflow modes
- [Quickstart: `pdb2reaction all`](quickstart-all.md) — run the end-to-end workflow from two PDBs
- [Quickstart: single-structure staged scan](quickstart-scan.md) — `--scan-lists/-s` driven MEP from one PDB
- [Quickstart: TS optimization](quickstart-tsopt-freq.md) — optimize and validate a TS candidate
- [CLI Conventions](cli-conventions.md) — flag precedence, atom/residue selectors, shared options
- [Troubleshooting](troubleshooting.md) and [Common Error Recipes](recipes-common-errors.md)

