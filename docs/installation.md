# Installation

`pdb2reaction` is intended for Linux environments (local workstations or HPC clusters) with a CUDA‑capable GPU. Several dependencies – notably **PyTorch**, **fairchem‑core (UMA)**, and **gpu4pyscf‑cuda12x** – expect a working CUDA installation.

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

## Quick start

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

> **Tip:** UMA is the default MLIP backend. To use ORB, MACE, or AIMNet2 instead, install the corresponding extra (e.g. `pip install 'pdb2reaction[orb]'`) and pass `--backend orb` to any command. See [Installation](#step-by-step-installation) step 7.

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


## Step-by-step installation

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

7. **(Optional) Install additional MLIP backends**

 pdb2reaction uses UMA by default. To use alternative backends, install the corresponding optional dependency:

 ```bash
 # ORB backend
 pip install 'pdb2reaction[orb]'

 # AIMNet2 backend
 pip install 'pdb2reaction[aimnet2]'

 # MACE: pip uninstall fairchem-core && pip install mace-torch (separate env required)
 ```

 To enable implicit solvent corrections, install [xTB](https://github.com/grimme-lab/xtb) and ensure the `xtb` command is available on your `PATH`.

8. **Verify installation**

 ```bash
 pdb2reaction --version
 ```

 This should display the installed version.
