# Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters), and production runs normally use a CUDA-capable GPU. Prebuilt **PyTorch** wheels include their CUDA runtime libraries: they need a compatible NVIDIA driver, but not a local CUDA toolkit.

## Quick start

For PyTorch, `nvidia-smi` shows `CUDA Version` at its top right, the newest CUDA the driver supports. Choose a wheel at or below it (`cu126`, `cu130`, or `cu132`). The commands below use the recommended `cu130`.

### Required

```bash
# 1) Create and activate a conda environment
# 2) Install a CUDA-enabled PyTorch build
# 3) Install pdb2reaction
# 4) Install headless Chrome for Plotly static image export (PNG)
#    Downloads a Chromium binary; requires internet access.

conda create -n p2r python=3.12 -y
conda activate p2r
TORCH_INDEX=cu130  # recommended; or cu126 / cu132
pip install 'torch==2.13.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
pip install pdb2reaction
plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded. This needs a free HF account with a read-only token, and you may need to accept the license on the UMA model page:

```bash
hf auth login
# or, with an access token in scripts:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

You only need to do this once per machine / environment.

Check the installation with `pdb2reaction --version`, which prints the installed version.

### Optional

For DMF, also install cyipopt right after activating the environment, as in step 3 of {ref}`Step-by-step installation <step-by-step-installation>`. ORB, AIMNet2, MACE, and DFT are installed in step 7 of the same section.

(step-by-step-installation)=
## Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load a CUDA toolkit only when the site/build requires one**

    A prebuilt PyTorch wheel does not require `nvcc`. If a dependency must be
    built from source, use `module avail cuda` and load the compiler/toolkit
    combination documented by the cluster:

    ```bash
    module load cuda/<your-version>   # e.g. cuda/12.6 or cuda/12.9
    ```

2. **Create and activate a conda environment**

    ```bash
    conda create -n <your-env> python=3.12 -y
    conda activate <your-env>
    ```

3. **Install cyipopt**
    Required if you want to use the DMF method (`--mep-mode dmf`) in MEP search. You can skip this step if you only use GSM. If `--mep-mode dmf` still stops with an import error, see {ref}`Installation / environment <installation-environment-problems>`.

    ```bash
    conda install -c conda-forge cyipopt -y
    ```

4. **Install PyTorch with the right CUDA build**

    Recommended example (`cu130`):

    ```bash
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    ```

    The official 2.13.0 matrix also provides `cu126`, `cu132`, and `cpu`.
    Choose the wheel with the `nvidia-smi` rule in the Quick start above, then check GPU access in step 8. See [PyTorch's version matrix](https://pytorch.org/get-started/previous-versions/).

5. **Install `pdb2reaction` itself and Chrome for visualization**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

6. **Log in to Hugging Face Hub (UMA model)**

    ```bash
    hf auth login
    ```

    For license requirements and non-interactive login, see the Required section above.

    Refer to the upstream projects for additional details:

    - fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
    - Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

7. **(Optional) Install additional MLIP backends**

    pdb2reaction uses UMA by default. For another backend, install its extra and select it with `-b/--backend` (for example, `-b orb`):

    ```bash
    # ORB backend (Requires Python 3.11 or 3.12; 3.12 recommended)
    pip install --only-binary=dm-tree "pdb2reaction[orb]"

    # AIMNet2 backend
    pip install "pdb2reaction[aimnet]"

    # MACE backend (use a separate conda environment because mace-torch
    # pins e3nn==0.4.4 which conflicts with UMA's fairchem-core;
    # pick the PyTorch wheel as in step 4)
    conda create -n <mace-env> python=3.11 -y && conda activate <mace-env> \
        && pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130 \
        && pip install pdb2reaction \
        && pip uninstall -y fairchem-core \
        && pip install 'mace-torch>=0.3.8'

    # DFT calculator and post-processing (`-b dft`, `--dft`, `pdb2reaction dft`)
    # [dft] installs the CUDA 13 GPU4PySCF build on Linux x86_64, for the cu130 / cu132
    # PyTorch wheels of step 4; with the cu126 wheel, install [dft-cuda12] instead.
    # On aarch64, GPU4PySCF must be built from source (https://github.com/pyscf/gpu4pyscf).
    # Without a GPU, [dft] still installs PySCF; run DFT with `--dft-engine cpu`.
    pip install "pdb2reaction[dft]"
    ```

    For when to use DFT and how to check an MLIP TS with it, see [Refine an MLIP TS with DFT](dft-backend.md).

8. **Verify installation**

    ```bash
    pdb2reaction --version
    hf auth whoami
    ```

    The first line should display the installed version, and the second your Hugging Face user name (the login that UMA downloads use). To verify GPU access:

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

**GPU / CUDA.** An NVIDIA GPU whose driver supports the chosen wheel (see Quick start); newer GPU architectures may need a newer wheel. CPU-only execution works but is usually much slower.

**VRAM, RAM, and disk.** Memory grows with the model, the atom count, and the Hessian mode, and the disk holds the environment, the model weights, and the trajectories and Hessians; run one representative calculation on the target node and watch the peak use.

## Next steps

- [Getting Started](getting-started.md) — the shortest run, and which page to read next
- [Quickstart: `pdb2reaction all`](quickstart-all.md) — build an MEP from R and P
- [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md) — build a path from one structure
- [Quickstart: TS-only mode](quickstart-tsopt.md) — optimize and check a TS candidate
- [Refine an MLIP TS with DFT](dft-backend.md) — refine and check the TS with DFT
- [Common options and selectors](cli-conventions.md) — shared options, and how to give residues and atoms
- [HPC example](hpc-example.md) — job script that spreads UMA workers over several GPU nodes
- [Troubleshooting](troubleshooting.md) — common errors and what to try
