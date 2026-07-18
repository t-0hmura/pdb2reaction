# Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters), and production runs normally use a CUDA-capable GPU. Prebuilt **PyTorch** wheels include their CUDA runtime libraries: they need a compatible NVIDIA driver, but not a local CUDA toolkit. A toolkit is needed only when building a CUDA extension or GPU package from source.

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

## Quick start

This example assumes the default GSM MEP mode (`--mep-mode gsm`). For DMF (`--mep-mode dmf`), install cyipopt via conda first. PyTorch 2.8.0 publishes `cu126`, `cu128`, and `cu129` wheels; choose the site's tested index for its driver and GPU architecture.

### Required

```bash
# 1) Install a CUDA-enabled PyTorch build
# 2) Install pdb2reaction
# 3) Install headless Chrome for Plotly static image export (PNG)
#    Downloads a Chromium binary; requires internet access.

TORCH_INDEX=cu126  # use cu128/cu129 when required by the GPU/site stack
pip install 'torch==2.8.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
pip install pdb2reaction
plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded (requires a free HF account with read-only token; you may need to accept the UMA model license at <https://huggingface.co/facebook/UMA>):

```bash
hf auth login
# or, with an access token in scripts:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

(Recent `huggingface_hub` releases ship the `hf` CLI; older versions still expose `huggingface-cli login`, which is being deprecated.)

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

- If an HPC site requires a CUDA module for locally built extensions, load the exact module documented by that site. Do not load a second CUDA runtime merely to install a prebuilt PyTorch wheel; first test the wheel with the NVIDIA driver alone.

  ```bash
  module load cuda/<your-version>   # e.g. cuda/12.6 or cuda/12.9
  ```

> **Tip:** UMA is the default MLIP backend. To use ORB or AIMNet2, install the corresponding extra (e.g. `pip install "pdb2reaction[orb]"`) and pass `-b/--backend orb` to any command. See step 7 below.

```{warning}
**MACE:** `mace-torch` requires `e3nn==0.4.4`, which conflicts with `fairchem-core`'s `e3nn>=0.5` pin (UMA). The two cannot coexist, so MACE needs a dedicated conda env; the canonical recipe is `pip uninstall -y fairchem-core && pip install mace-torch` in that env.
```


(step-by-step-installation)=
## Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load a CUDA toolkit only when the site/build requires one**

    A prebuilt PyTorch wheel does not require `nvcc`. If a dependency must be
    built from source, use `module avail cuda` and load the compiler/toolkit
    combination documented by the cluster:

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

    Conservative pre-Blackwell example:

    ```bash
    pip install 'torch==2.8.0' --index-url https://download.pytorch.org/whl/cu126
    ```

    The official 2.8.0 matrix also provides `cu128`, `cu129`, and `cpu`.
    Select by driver and GPU architecture, then verify with
    `torch.cuda.is_available()`; the `nvidia-smi` "CUDA Version" banner is not a
    local-toolkit version selector. See [PyTorch's version matrix](https://pytorch.org/get-started/previous-versions/).

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
    conda create -n <mace-env> python=3.11 -y && conda activate <mace-env> \
        && pip install pdb2reaction \
        && pip uninstall -y fairchem-core \
        && pip install mace-torch

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

    **For CPCM-X solvation model**, the conda-forge xtb does not include CPCM-X; build from source. See {ref}`recipe-cpcmx-build` for the recipe.

    To use a custom xTB binary, set the `xtb_cmd` key in your YAML config or use `calc.xtb_cmd` in Python.

8. **Verify installation**

    ```bash
    pdb2reaction --version
    ```

    This should display the installed version. To verify GPU access:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    If `CUDA: False`, inspect the installed wheel, scheduler GPU visibility,
    driver, and environment libraries before changing versions:

    ```bash
    python -m torch.utils.collect_env
    python -m pip check
    ```

## System requirements

**GPU / CUDA / VRAM.** Use one of PyTorch 2.8.0's official CUDA wheels (`cu126`, `cu128`, or `cu129`) that supports both the driver and GPU architecture. Newer GPU architectures may require a newer wheel; a local toolkit with the same numeric version is not required for prebuilt wheels. Required VRAM depends on backend/model, atom count, Hessian mode, precision, and active degrees of freedom. Pilot a representative production stage and monitor peak allocation; the smoke suite is a correctness check, not a production-memory estimate.

**RAM.** Size host memory from a representative run; dense Hessians, model loading, and concurrent worker/process stages can dominate.

**Disk.** Budget from the selected environment, backend weight caches, generated trajectories/Hessians, and optional Chromium installed by `plotly_get_chrome`. Check actual cache/environment sizes on the target filesystem before production runs.

CPU-only execution works but is usually much slower. Benchmark the selected
backend and model; there is no reliable fixed GPU/CPU ratio.

## Next steps

- [Getting Started](getting-started.md) — project overview, pipeline stages, and workflow modes
- [Quickstart: `pdb2reaction all`](quickstart-all.md) — run the end-to-end workflow from two PDBs
- [Quickstart: single-structure staged scan](quickstart-scan.md) — `--scan-lists/-s` driven MEP from one PDB
- [Quickstart: TS-only mode](quickstart-tsopt-freq.md) — validate a TS candidate end-to-end
- [CLI Conventions](cli-conventions.md) — flag precedence, atom/residue selectors, shared options
- [Troubleshooting](troubleshooting.md) and [Common Error Recipes](recipes-common-errors.md)
