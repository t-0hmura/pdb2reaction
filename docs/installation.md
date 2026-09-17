# Installation

The standard setup uses Linux, Python 3.12, and an NVIDIA GPU. The example below uses the PyTorch CUDA 13 wheel and requires a compatible NVIDIA driver. For another driver/GPU combination, select a wheel from [PyTorch's version matrix](https://pytorch.org/get-started/previous-versions/). The wheel includes its CUDA runtime.

(step-by-step-installation)=
## Quick start

### Required

Accept the [UMA model license](https://huggingface.co/facebook/UMA), then copy this block into a terminal. `hf auth login` prompts for your Hugging Face credentials.

```bash
conda create -n pdb2reaction python=3.12 pip -y
conda activate pdb2reaction
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cu130
pip install pdb2reaction
plotly_get_chrome -y
hf auth login
pdb2reaction --version
```

This sets up the default UMA backend and Chrome for PNG export. Reactivate the environment with `conda activate pdb2reaction` in a new terminal.

### Optional

Run only the commands needed for your workflow, in the same active environment:

| Feature | Install |
|---|---|
| ORB (`-b orb`) | `pip install --only-binary=dm-tree "pdb2reaction[orb]"` |
| AIMNet2 (`-b aimnet2`) | `pip install "pdb2reaction[aimnet]"` |
| DFT (`--dft` / `pdb2reaction dft`) | `pip install "pdb2reaction[dft]"` |
| MCP server | `pip install "pdb2reaction[mcp]"` |
| DMF paths (`--mep-mode dmf`) | `conda install -c conda-forge cyipopt "numpy>=2,<2.5" -y` |

ORB requires Python 3.11 or 3.12; the quick-start environment uses 3.12. ORB and AIMNet2 do not need UMA authentication. PyDMF is already a core dependency; DMF additionally needs cyipopt.

The DFT extra installs GPU4PySCF on x86_64. On aarch64, use CPU PySCF with `--engine cpu` for the `dft` command. See [DFT](dft.md).

MACE needs a separate environment because its e3nn dependency conflicts with UMA. See the [MACE installation recipe](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/pdb2reaction-install-backends/mace.md).

## Verify GPU access

```bash
python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
python -m pip check
```

On a cluster, run the GPU check inside an allocated GPU job. If CUDA is unavailable, check the driver, selected PyTorch wheel and scheduler GPU allocation; diagnostic details are available with `python -m torch.utils.collect_env`.

## System requirements

**GPU / CUDA / VRAM.** Use one of PyTorch 2.13.0's official CUDA wheels (`cu126`, `cu130`, or `cu132`) that supports both the driver and GPU architecture. Newer GPU architectures may require a newer wheel; a local toolkit with the same numeric version is not required for prebuilt wheels. Required VRAM depends on backend/model, atom count, Hessian mode, precision, and active degrees of freedom. Pilot a representative production stage and monitor peak allocation; the smoke suite is a correctness check, not a production-memory estimate.

**RAM.** Size host memory from a representative run; dense Hessians, model loading, and concurrent worker/process stages can dominate.

**Disk.** Budget from the selected environment, backend weight caches, generated trajectories/Hessians, and optional Chromium installed by `plotly_get_chrome`. Check actual cache/environment sizes on the target filesystem before production runs.

CPU-only execution works but is usually much slower. Benchmark the selected
backend and model; there is no reliable fixed GPU/CPU ratio.

## Next steps

- [Getting Started](getting-started.md) — project overview, pipeline stages, and workflow modes
- [Quickstart: `pdb2reaction all`](quickstart-all.md) — run the end-to-end workflow from two PDBs
- [Quickstart: scan-defined single-structure workflow](quickstart-scan.md) — `--scan-lists/-s` driven MEP from one structure
- [Quickstart: TS-only mode](quickstart-tsopt-freq.md) — validate a TS candidate end-to-end
- [CLI Conventions](cli-conventions.md) — flag precedence, atom/residue selectors, shared options
- [Troubleshooting](troubleshooting.md) and [Common Error Recipes](recipes-common-errors.md)
