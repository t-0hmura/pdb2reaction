# Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters), and production runs normally use a CUDA-capable GPU. Prebuilt **PyTorch** wheels include their CUDA runtime libraries: they need a compatible NVIDIA driver, but not a local CUDA toolkit.

## Quick install

`nvidia-smi` shows `CUDA Version` at its top right: the newest CUDA the driver supports. Choose a PyTorch wheel at or below it (`cu126`, `cu130`, or `cu132`). The commands below use the recommended `cu130`.

### Required

```bash
# 1) Create and activate a Python 3.12 environment (conda here; venv also works)
# 2) Install a CUDA-enabled PyTorch build
# 3) Install pdb2reaction
# 4) Install headless Chrome for Plotly static image export (PNG, JPEG, SVG, PDF)
#    Downloads a Chromium binary; requires internet access.

conda create -n p2r -c conda-forge python=3.12
conda activate p2r

pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130  # recommended; for cu126 / cu132, change the end
pip install pdb2reaction

plotly_get_chrome -y
```

UMA on a GPU and `--compile` need a compiler, and image export (PNG, JPEG, SVG, PDF) needs Chrome's libraries. Install them first if they are missing: {ref}`System prerequisites <system-prerequisites>`.

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded. This needs a free HF account with a read-only token, and you may need to accept the license on the UMA model page:

```bash
hf auth login
# or, with an access token in scripts:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

You only need to do this once per machine.

Check the installation with `pdb2reaction --version`, which prints the installed version.

### Optional

For DMF (`--mep-mode dmf`) or `--solvent` (xTB), create the environment instead with cyipopt, xTB, and PyTorch from conda-forge, then install pdb2reaction:

```bash
conda create -n p2r -c conda-forge python=3.12 \
  cyipopt xtb "pytorch=2.13.0=cuda130*"
conda activate p2r
pip install pdb2reaction
plotly_get_chrome -y
```

If the `CUDA Version` of `nvidia-smi` is below 13.0, replace `cuda130` with `cuda129`. ORB, AIMNet2, MACE, and DFT are installed in step 7 of {ref}`Step-by-step installation <step-by-step-installation>`.

(system-prerequisites)=
## System prerequisites

`pip install pdb2reaction` does not install a compiler or the libraries Chrome needs. A workstation usually has them; a minimal container or a bare cluster node often does not.

- **A compiler for UMA and `--compile`.** UMA on a GPU and `--compile` on a GPU need a C compiler (`gcc`); `--compile` on a CPU needs a C++ compiler (`g++`).
- **Libraries for Chrome.** PNG, JPEG, SVG and PDF figures (`energy-diagram`, `trj2fig`, and the diagrams of `all`) are drawn by headless Chrome, which `plotly_get_chrome -y` downloads. Chrome loads the operating-system libraries below.

Ubuntu/Debian (on Ubuntu 22.04 and Debian 12, write `libasound2` instead of `libasound2t64`; drop `sudo` if you are root):

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates fonts-liberation libasound2t64 libatk-bridge2.0-0 \
  libatk1.0-0 libatspi2.0-0 libcairo2 libcups2 libcurl4 libdbus-1-3 libexpat1 libgbm1 \
  libglib2.0-0 libgtk-3-0 libnspr4 libnss3 libpango-1.0-0 libudev1 libvulkan1 libx11-6 libxcb1 \
  libxcomposite1 libxdamage1 libxext6 libxfixes3 libxkbcommon0 libxrandr2 wget xdg-utils
```

Rocky/RHEL:

```bash
sudo dnf install -y alsa-lib at-spi2-atk at-spi2-core atk ca-certificates cairo cups-libs \
  dbus-libs expat glib2 gtk3 libX11 libXcomposite libXdamage libXext libXfixes libXrandr \
  libgcc liberation-fonts libxcb libxkbcommon mesa-libgbm nspr nss nss-util pango systemd-libs \
  vulkan-loader wget xdg-utils
```

Use the existing compiler. If gcc or CUDA is missing, load them from the system (on HPC, for example with `module load`). If the system has no C compiler or CUDA, or you do not want to use them, `conda install -c conda-forge cuda-nvcc` installs both into the environment; reactivate the environment afterwards. On HPC, ask the administrator for any missing Chrome libraries.

If a compiler or CUDA is missing, `pdb2reaction` stops and prints `Missing system prerequisite:` followed by these two options.

(step-by-step-installation)=
## Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load gcc and CUDA on HPC**

    UMA on a GPU compiles code at run time and needs gcc and CUDA. On HPC, check
    `module avail` and load the compiler and CUDA modules of the cluster:

    ```bash
    module load gcc cuda/<your-version>   # e.g. cuda/12.6 or cuda/12.9
    ```

    If the system has no C compiler or CUDA, or you do not want to use them, run
    `conda install -c conda-forge cuda-nvcc` after step 2 and reactivate the environment.

2. **Create and activate a conda environment**

    ```bash
    conda create -n <your-env> -c conda-forge python=3.12
    conda activate <your-env>
    ```

3. **Install cyipopt, xTB, and PyTorch from conda-forge**
    Needed only for DMF (`--mep-mode dmf`) or the xTB solvent correction (`--solvent`); it also installs PyTorch.

    ```bash
    conda install -c conda-forge cyipopt xtb "pytorch=2.13.0=cuda130*"   # cuda129 if nvidia-smi shows CUDA below 13.0
    ```

4. **Install PyTorch with the right CUDA build**

    Skip this step if you installed PyTorch in step 3. Recommended example (`cu130`):

    ```bash
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    ```

    The official 2.13.0 matrix also provides `cu126`, `cu132`, and `cpu`.
    See [PyTorch's version matrix](https://pytorch.org/get-started/previous-versions/).

5. **Install `pdb2reaction` itself and Chrome for visualization**

    ```bash
    pip install pdb2reaction
    plotly_get_chrome -y
    ```

    The compiler that UMA and `--compile` need and the libraries Chrome needs are system packages: see {ref}`System prerequisites <system-prerequisites>`.

    On aarch64, the `warp-lang` 1.18.0 wheel, which UMA loads through `fairchem-core`, needs glibc 2.35 and stops the import with `GLIBC_2.35' not found` on older systems (`ldd --version`); there, install `warp-lang==1.17.0`.

6. **Log in to Hugging Face Hub (UMA model)**

    ```bash
    hf auth login
    ```

    For license requirements and non-interactive login, see the Required section above.

    Refer to the upstream projects for additional details:

    - fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
    - Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

7. **(Optional) Install additional MLIP backends**

    pdb2reaction uses UMA by default. For another backend, install its package and select it with `-b/--backend` (for example, `-b orb`):

    **ORB**: `orb-models` conflicts with UMA's `fairchem-core`; install it in a separate Python 3.12 environment. Pick the PyTorch wheel as in step 4.

    ```bash
    conda create -n orb-env -c conda-forge python=3.12
    conda activate orb-env
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    pip install pdb2reaction
    pip uninstall fairchem-core && pip install orb-models
    ```

    **AIMNet2**:

    ```bash
    pip install "pdb2reaction[aimnet]"
    ```

    **MACE**: `mace-torch` pins `e3nn==0.4.4`, which conflicts with UMA's `fairchem-core`, so install it in a separate conda environment. Pick the PyTorch wheel as in step 4.

    ```bash
    conda create -n mace-env -c conda-forge python=3.12
    conda activate mace-env
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    pip install pdb2reaction
    pip uninstall fairchem-core && pip install mace-torch
    ```

    **DFT** (`-b dft`, `--dft`, `pdb2reaction dft`): `[dft]` installs the CUDA 13 GPU4PySCF build on Linux x86_64, for the cu130 / cu132 PyTorch wheels of step 4 and the `cuda130` build of step 3; with the cu126 wheel or the `cuda129` build, install `[dft-cuda12]` instead. On aarch64, build [GPU4PySCF](https://github.com/pyscf/gpu4pyscf) from source. Without a GPU, `[dft]` still installs PySCF; run DFT with `--dft-engine cpu`.

    ```bash
    pip install "pdb2reaction[dft]"
    ```

    For when to use DFT and how to check an MLIP TS with it, see [Validate the MLIP TS structure with DFT](dft-backend.md).

8. **Verify installation**

    ```bash
    pdb2reaction --version
    hf auth whoami
    ```

    The first line should display the installed version, and the second your Hugging Face user name (the login used for UMA downloads). To verify GPU access:

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

**OS.** Linux, or WSL on Windows.

**GPU / CUDA.** An NVIDIA GPU whose driver supports the chosen wheel (see Quick install); newer GPU architectures may need a newer wheel. CPU-only execution works but is usually much slower.

**VRAM, RAM, and disk.** Memory depends on the model, the atom count, and the Hessian mode, and the disk holds the environment, the model weights, and the trajectories and Hessians. As a rough guide for a Hessian with the default UMA model, the default finite-difference Hessian (`FiniteDifference`) fits about 900 atoms in 8 GB, 2,000 in 16 GB, 3,000 in 24 GB, and 10,000 in 96 GB of VRAM, and an analytical Hessian (`Analytical`) about 200, 400, 600, and 1,500 atoms. Even when the model fits in VRAM, above about 1,000 atoms each Hessian takes an hour or more, and a full mechanism study takes longer. For a very large system, run one representative calculation on the target node and watch the peak use.

## Next steps

- [Getting Started](getting-started.md) — the shortest run, and which page to read next
- [Quickstart: `pdb2reaction all`](quickstart-all.md) — build an MEP from R and P
- [Quickstart: `pdb2reaction all --scan-lists`](quickstart-scan.md) — build a path from one structure
- [Quickstart: TS-only mode](quickstart-tsopt.md) — optimize and check a TS candidate
- [Validate the MLIP TS structure with DFT](dft-backend.md) — refine and check the TS with DFT
- [Common options and selectors](cli-conventions.md) — shared options, and how to specify residues and atoms
- [HPC example](hpc-example.md) — job scripts for one GPU, and for UMA workers over several GPU nodes
- [Troubleshooting](troubleshooting.md) — common errors and what to try
