# Backends and environment

Per-backend install steps, CUDA and PyTorch pairing, and environment probes.
The install order, verification, and the failure table are in
[SKILL.md](SKILL.md).

## Core package

- Python 3.11 or newer; 3.12 is recommended. ORB needs 3.11 or 3.12.
- PyPI install: [Install order](SKILL.md#install-order). PyTorch must be
  installed first ([CUDA and PyTorch](#cuda-and-pytorch)).
- DMF: `conda install -c conda-forge cyipopt -y`; `pydmf` is a core dependency.

Install from source for development:

```bash
git clone https://github.com/t-0hmura/pdb2reaction.git pdb2reaction
cd pdb2reaction
pip install -e '.[orb,aimnet,dft]'
```

Upgrade with `pip install --upgrade pdb2reaction` and confirm with
`pdb2reaction --version`. Across minor versions, also re-check
`pdb2reaction <subcommand> --help` and the `summary.json` keys
([outputs](../pdb2reaction-overview/outputs.md)). To remove it, run
`pip uninstall pdb2reaction`, or drop the whole env with
`conda env remove -n <YOUR_ENV>`.

A combined env for UMA, ORB, AIMNet2, DFT, and xTB:

```yaml
name: <YOUR_ENV>
channels: [conda-forge, nvidia]
dependencies:
  - python=3.12
  - xtb                                # only for the xTB solvent correction
  - pip
  - pip:
      - --extra-index-url https://download.pytorch.org/whl/<cu_index>
      - torch==2.13.0
      - pdb2reaction[orb,aimnet,dft]
```

`<cu_index>` is one of `cu126`, `cu130`, `cu132`, and `cpu`. For MACE, use the
same template with `python=3.11` and plain `pdb2reaction`, then do the swap in
[MACE](#mace-separate-environment) inside the new env.

## CUDA and PyTorch

PyTorch 2.13.0 publishes Linux wheels for `cu126`, `cu130`, `cu132`, and `cpu`.
`nvidia-smi` shows `CUDA Version` at its top right, the newest CUDA the driver
supports; choose a CUDA wheel at or below it. `cu130` is the recommended choice.
Use the site administrator's tested module and wheel pair when one is supplied,
and `cpu` only when no NVIDIA GPU is assigned.

A prebuilt PyTorch wheel contains its CUDA runtime dependencies. Normal runs
therefore need the NVIDIA driver, not `nvcc`, `CUDA_HOME`, or a local toolkit.
Load a toolkit only when pip must build a CUDA extension, or when you build
PyTorch or GPU4PySCF from source:

```bash
module load <CUDA_MODULE> gcc            # HPC modulefile (from `module avail cuda`)
export CUDA_HOME=/usr/local/cuda         # or a system install
export PATH="$CUDA_HOME/bin:$PATH"
conda install -c nvidia cuda-toolkit=<MAJOR.MINOR>   # or inside the conda env
nvcc --version
```

Add a CUDA module to PBS or SLURM jobs only when a locally built extension
needs it; test prebuilt wheels without a second CUDA runtime from a module.
OpenMPI is unrelated to CUDA and is needed only for the site's multi-node Ray
setup.

Install and check PyTorch:

```bash
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/<cu_index>
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

If `available` is `False` while `nvidia-smi` works, capture the evidence before
you switch wheels:

```bash
python -m torch.utils.collect_env
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.__file__)"
python -m pip check
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES-<unset>}"
```

Common causes are a CPU wheel, an unassigned or hidden GPU, a driver and wheel
mismatch, an unsupported GPU architecture, and mixed libraries from a module or
`LD_LIBRARY_PATH`.

Wheels install their CUDA libraries under `site-packages/nvidia/` and preload
them, but a module or `LD_LIBRARY_PATH` can still put an incompatible
`libcusolver`, `libcudnn`, `libnvrtc`, or `libnvJitLink` first. The symptoms
are `OSError: libcusolver.so.11: cannot open shared object file`,
`Could not load symbol cublasLtCreate`, and
`undefined symbol: cusparseLoggerSetCallback`. Compare with a clean
environment:

```bash
env -u LD_LIBRARY_PATH python -c "import torch; print(torch.cuda.is_available())"
```

If the clean check works, remove the conflicting module or path entry from the
job rather than hard-coding one torch subdirectory. Do not use `PYTORCH_NO_CUDA_PRELOAD`;
it is not a documented PyTorch control. A source-built extension that needs a
toolkit path must use the toolkit it was built with.

CPU only: install from the `cpu` index. MLIP backends run on CPU but usually
much slower; measure a representative structure. DFT does not fall back to CPU:
with the default `--dft-engine gpu`, a GPU failure stops the run, so pass
`--dft-engine cpu` explicitly.

On aarch64 (`uname -m`), PyTorch publishes CUDA wheels for recent versions, but
the `gpu4pyscf-cuda13x` wheel is x86_64 only, so the `[dft]` extra gives CPU
PySCF there ([DFT](#dft-pyscf-gpu4pyscf)). Check each MLIP backend's PyPI page
for aarch64 wheels.

## UMA

UMA is the default backend and ships through `fairchem-core`, a core
dependency, so it needs no extra. Being the default is not a claim that it is
the most accurate model for every system. It runs the OMol25 task, which covers
organic and inorganic molecules, transition-metal complexes, and electrolytes;
check that the target chemistry lies in that domain.

Every UMA checkpoint is in one gated Hugging Face repo, `facebook/UMA`. Open
<https://huggingface.co/facebook/UMA> and accept the FAIR Chemistry License v1
(approval is manual), with the same account that owns the token. Then:

```bash
hf auth login                 # paste a Read token from huggingface.co/settings/tokens
hf auth whoami                # prints the account that UMA downloads use
python -c "from importlib.metadata import version; print('fairchem-core:', version('fairchem-core'))"
```

`hf` comes from `huggingface_hub`; if the command is missing, run
`pip install 'huggingface_hub[cli]'`. The token is cached in
`~/.cache/huggingface/`, and later runs and batch jobs pick it up. On
`GatedRepoError` or `401 Client Error: Unauthorized`, check that the token's
account has been granted access on `facebook/UMA`, then log in again. The
individual models are checkpoint files inside that repo, not separate repos.

Pick a model with `--backend-model` or `calc.model` in a `--config` YAML:
`uma-s-1p2` (default) or `uma-m-1p1`. `p` replaces the decimal point
(`1p2` is 1.2). The `facebook/UMA` repo is the authoritative checkpoint list.
UMA runs in fp32 by default; `--precision fp64` also forces an fp64 Hessian and
can change TS optimization and Hessian numerics.

Several UMA workers can share a heavy MEP search, through `--uma-workers` and
`--uma-workers-per-node` or a YAML file:

```yaml
calc:
  workers: 4
  workers_per_node: 4
```

On one node, the Ray worker pool starts locally and all workers must see the
same GPUs (for example `CUDA_VISIBLE_DEVICES=0,1,2,3`). pdb2reaction does not
start a cross-node Ray cluster itself; start it under the scheduler and export
`RAY_ADDRESS` ([HPC example](../../docs/hpc-example.md)). Workers add process
and communication overhead, so benchmark before assuming they cut wall time.
With more than one worker, `--hessian-calc-mode Analytical` raises
`BackendError` (a `RuntimeError` subclass) rather than changing the explicitly
requested method; use `FiniteDifference` (the default) or one worker.

Pitfalls:

- `freq` runs out of VRAM: compare Hessian modes and model sizes on a pilot,
  reduce the Hessian target, or move Hessian assembly to CPU.
- The first call is slow: checkpoint download and cache filling
  (`~/.cache/huggingface/hub/`) or model start-up dominate.
- `Ray actor died` is not specific. Capture the actor traceback and the
  scheduler and Ray logs, check that every worker has the same env and visible
  devices, and reproduce with one worker before changing CUDA packages.

## ORB

Use Python 3.12 (orb-models 0.7 or newer) or 3.11 (orb-models 0.5.x). Weights
download on first use, with no login.

```bash
pip install 'pdb2reaction[orb]'   # or: pip install orb-models
python -c "import orb_models; print('orb_models:', orb_models.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='orb', charge=0, spin=1)"
```

The default model is `orb_v3_conservative_omol`: `conservative` means forces
are derived from the energy, and `omol` means training on OMol25. That names
the model family, not its accuracy for a given reaction. pdb2reaction runs ORB
in fp64 by default. `--precision fp32` selects ORB's `float32-high` mode, which
lowers CUDA matmul precision, so compare it with fp64 on the target system.

Pitfalls:

- Extra imaginary modes after a finite-difference Hessian: check the precision
  in use, inspect the modes, and recompute before you classify the point.
- n_imag above 1 after `tsopt` is not a first-order saddle even when the
  optimizer converged. Retry from a better MEP seed or with `--flatten`, and
  cross-check with UMA or MACE if the result stays ambiguous.

## MACE (separate environment)

`mace-torch` pins `e3nn==0.4.4`, while UMA's `fairchem-core` needs `e3nn>=0.5`,
so MACE needs its own env. UMA does not run there.

```bash
conda create -n <YOUR_MACE_ENV> python=3.11 -y
conda activate <YOUR_MACE_ENV>
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/<cu_index>
pip install pdb2reaction                 # pulls fairchem-core
pip uninstall -y fairchem-core           # remove UMA's e3nn pin
pip install 'mace-torch>=0.3.8'          # installs the e3nn MACE needs
python -c "import mace; print('mace:', mace.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='mace', charge=0, spin=1)"
```

`pip check` then reports `fairchem-core` as missing; that is expected, so do not
reinstall it here. `pip install --upgrade pdb2reaction` can pull it back; repeat
the uninstall and install and the smoke check afterwards. If UMA and MACE end up
in one env, imports fail with an `e3nn` version error; remove the env with
`conda env remove -n <YOUR_MACE_ENV>` and rebuild it.

The default model is `MACE-OMOL-0`, run in fp64. MACE has no multi-GPU sharding, and the first call can include download, model
loading, and runtime start-up. A model being available does not put a reaction
inside its reliable domain; validate stationary points with an independent
frequency calculation and IRC.

Pitfalls:

- `RuntimeError: Expected all tensors to be on the same device`: capture the
  full traceback and check that the requested device and any custom tensors
  agree; reproduce in a fresh process before changing the env.
- Hessians are slow in fp64: benchmark `--precision fp32` on the target system
  and keep the independent frequency and IRC checks.

## AIMNet2

The default `aimnet2` model covers H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br,
and I. It does not cover Mg or first-row transition metals such as Zn, Mn, or
Fe. `aimnet2-pd` adds Pd with a different training domain, not other metals.

```bash
pip install 'pdb2reaction[aimnet]'       # or: pip install 'aimnet>=0.2.0'
python -c "import aimnet; print('aimnet:', aimnet.__version__)"
python -c "from pdb2reaction.backends import create_calculator; create_calculator(backend='aimnet2', charge=0, spin=1)"
```

Other models in `aimnet>=0.2.0`, chosen with `--backend-model`, are
`aimnet2-2025` (general organic), `aimnet2-nse` (open-shell and radicals),
`aimnet2-pd` (Pd), and `aimnet2-rxn` (H/C/N/O reactions; closed-shell,
net-neutral systems only, for relative energies within one composition). The
model is never switched automatically by multiplicity. AIMNet2 runs in fp32
only; `--precision fp64` and `--deterministic` both stop the run with an error.

Use AIMNet2 when the selected model covers every element and the electronic
state, for a CPU-capable baseline, or to pre-screen candidates in its domain.
Do not use it for clusters with unsupported elements, and always confirm TS
curvature with an independent frequency calculation and IRC.

Pitfalls:

- `KeyError` on an element: the model does not cover it; choose a model or
  backend whose element list includes the whole cluster.
- Give the total cluster charge `-q` and multiplicity `-m`; they are model
  inputs, not per-atom charges
  ([charge and multiplicity](../pdb2reaction-model-setup/SKILL.md#charge-and-multiplicity)).
- Radicals: select `aimnet2-nse`, pass the real multiplicity, and validate
  independently.
- Do not use `aimnet2-rxn` outside closed-shell, net-neutral H/C/N/O systems.

## DFT (PySCF, GPU4PySCF)

DFT is optional. `pdb2reaction dft` gives single points with population
analysis, and `-b dft` runs the geometry commands with PySCF on CPU or
GPU4PySCF on CUDA x86_64.

```bash
pip install 'pdb2reaction[dft]'            # cu130 / cu132 wheels
pip install 'pdb2reaction[dft-cuda12]'     # cu126 wheel
python -c "import pyscf; print('pyscf:', pyscf.__version__)"
python -c "import gpu4pyscf; print('gpu4pyscf:', gpu4pyscf.__version__)"
python -c "import cupy; print('cupy:', cupy.__version__)"
```

`[dft]` pulls `pyscf>=2.13.0` and `basis-set-exchange>=0.11` everywhere, and
`gpu4pyscf-cuda13x>=1.8.1,<2` with `cupy-cuda13x>=13.6,<15` on Linux x86_64
only. `cutensor-cu13` is a separate optional install; add it only when the
GPU4PySCF path in use requires it. On aarch64 the install succeeds without
GPU4PySCF, leaving CPU PySCF; build GPU4PySCF from source
(<https://github.com/pyscf/gpu4pyscf>) to use `--dft-engine gpu` there.

A source build that works with Python 3.12 and a CUDA 12 toolkit module:

```bash
pip install 'pyscf>=2.13.0' pyscf-dispersion cupy-cuda12x
git clone --depth 1 --branch v1.8.1 https://github.com/pyscf/gpu4pyscf.git
cd gpu4pyscf
cmake -S gpu4pyscf/lib -B build/temp.gpu4pyscf -DCUDA_ARCHITECTURES=90-real -DBUILD_LIBXC=ON
cmake --build build/temp.gpu4pyscf -j "$(nproc)"
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
```

Set `CUDA_ARCHITECTURES` to the GPU's compute capability (`90-real` for
Hopper). Build in place and use `PYTHONPATH`: the package's `setup.py`
expects the separate libxc wheel, which has no aarch64 build. The libxc step
downloads its sources, so the build node needs network access. Before
production, run one small GPU SCF in the same environment, for example water
with `wb97m-v/def2-svp`.

`--dft-engine gpu` is the default. When the GPU path fails it stops with
`Use --dft-engine cpu to explicitly run on CPU.` and does not fall back. Use
`cpu` on aarch64, without a compatible GPU wheel, or for an explicit CPU run;
it is usually slower for large hybrid-DFT jobs, but measure rather than assume
a factor.

```bash
pdb2reaction dft -i ts.pdb -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-svp' --dft-engine gpu
```

The default `--func-basis` is `wb97m-v/def2-svp`; `pdb2reaction dft --help`
lists the other defaults.

| Symptom | Likely cause | Fix |
|---|---|---|
| `OSError: libcusolver.so.11 not found` | Missing or mixed CUDA wheel dependency, or a library collision | Clean-environment check in [CUDA and PyTorch](#cuda-and-pytorch); do not hard-code a guessed library path |
| `cupy.cuda.runtime.CUDARuntimeError: invalid device ordinal` | Device index outside the scheduler-visible set | Keep the scheduler's `CUDA_VISIBLE_DEVICES` and use a local ordinal (usually 0 in a one-GPU job) |
| `RuntimeError: CUDA out of memory` mid-SCF | Method and system exceed VRAM | Same method with `--dft-engine cpu` or a larger GPU; a smaller basis or grid is a different method and must be labeled and revalidated |
| SCF stalls or fails near start-up | The message does not identify the dependency | Capture the full log, run `pip check`, compare with GPU4PySCF's official requirements before adding CUDA libraries |
| aarch64: GPU requested but no `gpu4pyscf` | x86_64-only wheel | `--dft-engine cpu` or a source build |

Atom count alone does not set memory or wall time; basis, elements, functional,
grid, density fitting, and GPU model all matter. Run one representative single
point with memory and VRAM monitoring before a batch. There is no multi-GPU
SCF, so do not plan capacity by summing VRAM across GPUs.

## xTB solvent correction

`--solvent NAME` adds an xTB solvent-minus-vacuum delta to the base MLIP surface:

```text
ΔE = E_xTB(solvent) - E_xTB(vacuum)
E_total = E_base + ΔE
```

Forces and Hessians use the matching difference. The correction is off by
default and expensive, because each evaluation runs xTB in solvent and in
vacuum. It is meant mainly for small molecules in solution; about 200–300 atoms
is a practical upper size, not a hard limit, so benchmark the actual system.
`--solvent-model` selects `alpb` (default) or `cpcmx`. The option models a
named bulk solvent, so use it on an enzyme cluster only when that dielectric is
justified. To compare a solution barrier with a cluster barrier, keep the
reacting species, charge, multiplicity, backend and model, and energy
references the same.

pdb2reaction calls the standalone `xtb` executable, not its Python bindings,
and it must stay on `PATH` in batch jobs:

```bash
conda install -c conda-forge xtb
xtb --version
```

The YAML keys are in the [YAML reference](../../docs/yaml-reference.md); the CLI
options are listed by `--help-advanced`. When xTB SCC convergence is poor, pass
extra options through the command, for example
`--solvent-xtb-cmd 'xtb --etemp 1000'`.

The correction applies to the geometry and energy commands.
`trj2fig` also accepts it, but only uses it when `-q` or `-m` triggers frame
energy recomputation. It is not accepted by `dft`, `extract`, or
structure-only utilities; `dft` has its own `--solvent` for PySCF PCM and SMD.

## Probe the compute environment

Use these probes only when the host is unknown. Run the report once, then the
follow-up commands for the scheduler it found.

```bash
{
  echo "=== Scheduler ==="
  if [[ -n "${PBS_JOBID:-}${PBS_ENVIRONMENT:-}" ]]; then
    SCHED=pbs
  elif [[ -n "${SLURM_JOB_ID:-}${SLURM_CLUSTER_NAME:-}" ]]; then
    SCHED=slurm
  else
    has_pbs=0; has_slurm=0
    command -v qsub >/dev/null && has_pbs=1
    command -v sbatch >/dev/null && has_slurm=1
    if (( has_pbs && has_slurm )); then SCHED=ambiguous
    elif (( has_pbs )); then SCHED=pbs
    elif (( has_slurm )); then SCHED=slurm
    else SCHED=local
    fi
  fi
  echo "scheduler: $SCHED"
  [[ "$SCHED" != ambiguous ]] || echo "Both clients are visible; select from site documentation or an active allocation before submitting." >&2

  echo; echo "=== Architecture ==="
  uname -mrs
  command -v lscpu >/dev/null && lscpu | grep -E "^(Architecture|Model name|CPU\(s\)):"

  echo; echo "=== GPU ==="
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>&1 || echo "no GPU"

  echo; echo "=== CUDA toolkit ==="
  command -v module >/dev/null && module avail cuda 2>&1 | head -20
  command -v nvcc >/dev/null && nvcc --version
  ls "$(conda info --base 2>/dev/null)/envs"/*/bin/nvcc 2>/dev/null

  echo; echo "=== PBS queues (if PBS) ==="
  if [[ "$SCHED" == pbs ]]; then
    command -v qstat >/dev/null && qstat -Q 2>/dev/null
    command -v pbsnodes >/dev/null && \
      pbsnodes -a 2>/dev/null | grep -E "^[^[:space:]]|^[[:space:]]+(np|properties|gpus|resources_available\.(ncpus|ngpus|mem)|totalmem)" | head -80
  fi

  echo; echo "=== SLURM partitions (if SLURM) ==="
  if [[ "$SCHED" == slurm ]]; then
    command -v sinfo >/dev/null && sinfo -o "%P %l %c %m %N %G" 2>/dev/null
  fi

  echo; echo "=== Conda envs with pdb2reaction ==="
  while IFS= read -r env_prefix; do
    conda run -p "$env_prefix" python -c \
      'import pdb2reaction; print(pdb2reaction.__version__)' 2>/dev/null \
      && printf 'prefix: %s\n' "$env_prefix"
  done < <(conda env list --json 2>/dev/null | python -c \
    'import json,sys; print("\n".join(json.load(sys.stdin)["envs"]))')

  echo; echo "=== Loaded modules ==="
  command -v module >/dev/null && module list 2>&1
} 2>&1
```

Follow-up per scheduler:

```bash
qstat -Qf <YOUR_QUEUE>            # PBS: full queue config, including resources_max.walltime
qstat -u "$USER"                  # PBS: your jobs
scontrol show partition           # SLURM: full partition table
squeue -u "$USER"                 # SLURM: your jobs
```

Reading the report:

- aarch64: the GPU4PySCF wheel is missing, so DFT uses CPU PySCF unless a
  source-built GPU4PySCF passes an import and a representative SCF. The MLIP
  backends work where compatible wheels exist.
- No GPU: run DFT with `--dft-engine cpu`; MLIP backends run on CPU. Drop the
  GPU request from the job preamble.
- With a GPU, note the driver version and VRAM; they bound the wheel and the
  model size.
- No toolkit found does not block a prebuilt CUDA wheel
  ([CUDA and PyTorch](#cuda-and-pytorch)).
- `pbsnodes -a` and `sinfo` give node capacities, which are upper bounds, not
  requests. Size requests from the measured workload and queue policy; most
  single-backend geometry jobs request one GPU.
- The env that imports pdb2reaction is `<YOUR_ENV>`. If none does, follow
  [Install order](SKILL.md#install-order).

| Placeholder | How to fill it |
|---|---|
| `<YOUR_QUEUE>` | A queue from `qstat -Q` (PBS) whose `resources_max.walltime` covers the job |
| `<YOUR_PARTITION>` | A partition from `sinfo` (SLURM) whose time limit covers the job |
| `<N_NODES>` | Nodes for a tested multi-node worker setup, within site limits; not from task count alone |
| `<NCPU>` | CPU budget sized from measured CPU-side work and xTB or DFT threading, not the whole node by default |
| `<NGPU>` | GPUs the backend and workflow can use; usually 1 unless a tested UMA worker setup is used |
| `<MEM>` | Measured peak RAM plus headroom, within node capacity |
| `<CUDA_MODULE>` | Empty for prebuilt wheels; for a source build only, the exact module from `module avail` with its compiler |
| `<YOUR_ENV>` | The conda env that imported pdb2reaction |
| `<HH:MM:SS>` | Estimated walltime, capped by the queue limit |

Do not save the raw report in a project or repository: it can contain private
host names, paths, scheduler policy, and env names. If you need a file, write it
under `${TMPDIR:-/tmp}` with mode `0600`, redact it, and delete it after copying
the placeholder values.

## See also

- [SKILL.md](SKILL.md): install order, verification, failure table.
- [pdb2reaction-hpc](../pdb2reaction-hpc/SKILL.md): job scripts that use the placeholders above.
- [tsopt](../pdb2reaction-cli/tsopt.md) and [freq](../pdb2reaction-cli/freq.md): Hessian mode and TS checks per backend.
- [dft](../pdb2reaction-cli/dft.md): the `dft` command and DFT//MLIP single points.
- Docs: [Installation](../../docs/installation.md), [MLIP Backends](../../docs/backends.md), [Refine an MLIP TS with DFT](../../docs/dft-backend.md).
