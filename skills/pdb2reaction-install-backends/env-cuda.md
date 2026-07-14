# CUDA + PyTorch setup (env-cuda.md)

Prerequisite: driver version, CPU architecture, and CUDA source from
`pdb2reaction-env-detect/SKILL.md`.

## Step 1. Pick an official PyTorch 2.8 wheel

`pdb2reaction` pins `torch~=2.8.0`. PyTorch's official 2.8.0 matrix
publishes Linux/Windows wheels for `cu126`, `cu128`, `cu129`, and `cpu`.
It does not publish a 2.8.0 wheel on `cu118`, `cu121`, or `cu124`.

Pick the index supported by the site's driver **and the GPU architecture**:

- use the cluster administrator's tested module/wheel combination when one is
  supplied;
- `cu126` is the conservative starting point for pre-Blackwell hardware;
- use `cu128` or `cu129` when the GPU architecture or a dependency explicitly
  requires it (for example, a wheel built for newer Blackwell devices);
- use `cpu` only when no NVIDIA GPU is assigned.

Do not convert the `nvidia-smi` "CUDA Version" banner directly into a wheel
index: it reports the newest CUDA version the driver advertises, not a locally
installed toolkit. NVIDIA documents CUDA 12.x minor-version compatibility from
Linux driver 525, but PTX and features introduced with newer toolkits can still
require a newer driver. The decisive check is the smoke test in Step 3.

## Step 2. Decide whether a local CUDA toolkit is needed

A prebuilt PyTorch wheel contains its CUDA runtime dependencies. Normal
pdb2reaction inference therefore needs an NVIDIA driver, not `nvcc`,
`CUDA_HOME`, or a separately installed toolkit.

Load/install a toolkit only if pip must build a CUDA extension or you are
building PyTorch/GPU4PySCF from source. In that case, use one of these setups:

**Setup A — HPC modulefile** (`module avail cuda` had hits):

```bash
module load <CUDA_MODULE>           # exact name from `module avail cuda`
module load gcc                     # toolchain for the CUDA build (see note below)
module load <OPENMPI_MODULE>        # only when running multi-node Ray (`--workers > 1`); single-node runs do not need it
nvcc --version                      # confirm
echo "$CUDA_HOME"                   # often set by the module
```

Add `module load <CUDA_MODULE>` and a compatible compiler to PBS/SLURM jobs
only when the runtime or a locally built extension depends on that module.
Prebuilt PyTorch wheels should first be tested without injecting a second CUDA
runtime through a module. OpenMPI is unrelated to CUDA visibility; load it only
as required by the site's external multi-node Ray procedure.

**Setup B — system install** (e.g. `/usr/local/cuda` on a workstation):

```bash
export CUDA_HOME=/usr/local/cuda
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$LD_LIBRARY_PATH"
nvcc --version
```

**Setup C — toolkit inside a conda env** (only for builds):

```bash
conda activate <YOUR_ENV>
conda install -c nvidia cuda-toolkit=<MAJOR.MINOR>   # e.g. 12.6
```

## Step 3. Install torch matching `<cu_index>`

```bash
pip install torch --index-url https://download.pytorch.org/whl/<cu_index>
```

Verify:

```bash
python -c "
import torch
print('torch    :', torch.__version__)
print('cuda     :', torch.version.cuda)
print('available:', torch.cuda.is_available())
print('device 0 :', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')
"
```

If `torch.cuda.is_available()` is `False` despite a working `nvidia-smi`, do
not immediately swap wheel indexes. Capture the evidence first:

```bash
python -m torch.utils.collect_env
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.__file__)"
python -m pip check
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES-<unset>}"
```

Common causes include a CPU wheel, an unassigned/hidden GPU, a driver/runtime
mismatch, an unsupported GPU architecture, or mixed libraries from a module or
`LD_LIBRARY_PATH`. Choose another index only after identifying which applies.

## Step 4. Diagnose library-loading collisions

Torch wheels install CUDA libraries under `site-packages/nvidia/` and contain
their own preload logic. A system/module `LD_LIBRARY_PATH` can still select an
incompatible `libcusolver`, `libcudnn`, `libnvrtc`, or `libnvJitLink` first.
Compare against a clean environment:

```bash
env -u LD_LIBRARY_PATH python -c "import torch; print(torch.cuda.is_available())"
python -m pip check
```

If the clean check works, remove the conflicting CUDA module/path entry from
the job rather than hard-coding one torch subdirectory. Do not use
`PYTORCH_NO_CUDA_PRELOAD`; it is not a documented PyTorch control variable.
If a source-built extension genuinely needs a toolkit path, make that path
match the toolkit used to build the extension and retest the complete stack.

Symptoms that you have this problem:
`OSError: libcusolver.so.11: cannot open shared object file`,
`Could not load symbol cublasLtCreate`,
`undefined symbol: cusparseLoggerSetCallback`.

## Step 5. CPU-only fallback

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

`pdb2reaction` runs MLIP backends on CPU but is usually much slower; measure a
representative structure rather than assuming a fixed slowdown. For DFT
(`pdb2reaction dft`), CPU PySCF is **not** an
automatic fallback — pass `--engine cpu` explicitly when the GPU
backend is unavailable; with the default `--engine gpu` the command
raises a `ClickException` rather than silently falling back. See
`dft.md`.

## Architecture quirk: aarch64

If `uname -m` reports `aarch64` (e.g. ARM-based servers, Apple Silicon
under Linux containers, some HPC nodes):

- torch wheels exist for aarch64 + CUDA on recent versions; check
  `https://download.pytorch.org/whl/torch/`.
- **The `gpu4pyscf-cuda12x` PyPI wheel is x86_64 only.** The packaged extra
  therefore uses CPU PySCF on aarch64. A compatible source-built GPU4PySCF
  environment is an expert, opt-in alternative and must be validated first;
  see `dft.md`.
- UMA / Orb / MACE / AIMNet2 wheels: check the backend's PyPI page.

## Loaded-state checks for an existing env

```bash
conda activate <YOUR_ENV>
python - <<'PY'
import torch, sys
print(f"python   : {sys.version.split()[0]}")
print(f"torch    : {torch.__version__}")
print(f"cuda     : {torch.version.cuda}")
print(f"cudnn    : {torch.backends.cudnn.version()}")
print(f"available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"device 0 : {torch.cuda.get_device_name(0)} ({torch.cuda.get_device_properties(0).total_memory // 1024**3} GB)")
PY
```

## See also

- `core.md` — install `pdb2reaction` itself (after torch is healthy).
- Backend mds (`uma.md`, `mace.md`, …) — extras that piggyback on the
  torch you just installed.
- `dft.md` — `gpu4pyscf-cuda12x` install + aarch64 fallback.
