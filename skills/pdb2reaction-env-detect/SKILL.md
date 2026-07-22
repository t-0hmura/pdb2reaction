---
name: pdb2reaction-env-detect
description: Fallback skill for detecting the current compute environment (local / PBS / SLURM, CPU architecture, GPU, CUDA, conda env). Use only when the env is unknown — other skills assume placeholders are filled by the user or context.
---

# pdb2reaction Environment Detection

## Purpose

Most `pdb2reaction` skills assume you already know your environment and
fill in placeholders like `<YOUR_QUEUE>`, `<NCPU>`, `<NGPU>`, `<CUDA_MODULE>`,
`<YOUR_ENV>` directly. **Use this skill only when those are unknown** — for
example, the first time you run on a new host, or when an automated agent
has no prior context about the cluster.

After running this discovery sequence, you should be able to populate
every placeholder used by `pdb2reaction-hpc/SKILL.md`,
`pdb2reaction-install-backends/{env-cuda,core}.md`, and the per-backend
install pages.

## Discovery sequence

Run these one by one. The goal is to populate the placeholder table at
the bottom.

### 1. Scheduler

```bash
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
```

### 2. Architecture and OS

```bash
uname -m                          # x86_64 / aarch64
uname -s                          # Linux / Darwin
lscpu | grep -E "^(Architecture|Model name|CPU\(s\)):"
```

On `aarch64` (ARM64), the packaged **`gpu4pyscf-cuda12x` PyPI wheel is not
available**, so the supported extra uses CPU PySCF. An expert may instead use
a separately source-built GPU4PySCF environment, but must validate imports and
a representative SCF before selecting `--engine gpu`. Other backends (UMA /
MACE / Orb / AIMNet2) work when compatible wheels exist for the platform.

### 3. GPU

```bash
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>/dev/null \
  || echo "no NVIDIA GPU detected"
```

If no GPU: stay on CPU — for DFT use `--engine cpu`; MLIP backends auto-fall back to CPU. Omit `gpus=N` from the PBS preamble.
If a GPU is present, note the **driver version** and **VRAM** — both
constrain which torch CUDA index and which model size you can use.

### 4. CUDA toolkit

Three places it might live:

```bash
# 4a. HPC modulefile (most common on managed clusters)
command -v module >/dev/null && module avail cuda 2>&1 | head -20

# 4b. System-wide install (apt/yum, /usr/local/cuda)
command -v nvcc && nvcc --version

# 4c. Inside a conda env
ls "$(conda info --base 2>/dev/null)/envs"/*/bin/nvcc 2>/dev/null
```

If none of the three returns anything, a local CUDA **toolkit** is not
installed. That does not prevent a prebuilt CUDA-enabled PyTorch wheel from
using the GPU: the wheel supplies its CUDA runtime libraries and needs the
NVIDIA driver. Install/load a toolkit only when a dependency must compile a
CUDA extension (or when building PyTorch/GPU4PySCF from source). See
`pdb2reaction-install-backends/env-cuda.md` before changing the environment.

### 5. PBS scheduler details (when `SCHED=pbs`)

```bash
qstat -Q                          # list queues (counts only)
qstat -Qf <queue>                 # full queue config including resources_max.walltime
pbsnodes -a 2>/dev/null | grep -E "^[a-z0-9]|^ *(np|properties|gpus)" | head
qstat -u "$USER"                  # your running / queued jobs
```

The output of `pbsnodes -a` gives node capacities. They are upper bounds, not
automatic requests: choose `<NCPU>`, `<NGPU>`, and memory from the measured
workload and queue policy (most single-backend geometry jobs request one GPU).

### 6. SLURM scheduler details (when `SCHED=slurm`)

```bash
sinfo -o "%P %l %N %G"            # partition, max time, nodes, gres
scontrol show partition           # full partition table
squeue -u "$USER"
```

### 7. Conda envs

```bash
conda env list                    # all envs
# is pdb2reaction installed somewhere?
for env in $(conda env list | awk '/^[a-zA-Z]/{print $1}'); do
  conda run -n "$env" python -c 'import pdb2reaction; print("'"$env"': ", pdb2reaction.__version__)' 2>/dev/null
done
```

The env that successfully imports `pdb2reaction` is your `<YOUR_ENV>`.
If none does, see `pdb2reaction-install-backends/SKILL.md`.

### 8. Loaded modules

```bash
command -v module >/dev/null && module list 2>&1
```

## Mapping discovery output to placeholders

| Placeholder | How to fill it from the output above |
|---|---|
| `<YOUR_QUEUE>` | A queue from `qstat -Q` (PBS); inspect with `qstat -Qf <queue>` to check `resources_max.walltime` covers your job |
| `<YOUR_PARTITION>` | A partition from `sinfo -o "%P %l %N %G"` (SLURM) whose `TIMELIMIT` covers your job |
| `<N_NODES>` | Number of nodes for a tested multi-node dispatcher/worker setup, within queue/site limits; do not infer it from task count alone |
| `<NCPU>` | Requested CPU budget, no larger than node/queue capacity; size from measured CPU-side work and xTB/DFT threading, not automatically the whole node |
| `<NGPU>` | Requested GPU count supported by the chosen backend/workflow; usually 1 unless using a tested UMA worker configuration |
| `<MEM>` | Measured peak RAM plus headroom, no larger than queue/node capacity (`pbsnodes -a` / `sinfo -o "%m"` reveal capacity, not the required request) |
| `<CUDA_MODULE>` | Leave empty for prebuilt wheels. For a source-built CUDA extension only, use the exact site module verified with `module avail` and its compatible compiler. |
| `<YOUR_ENV>` | The conda env that imported `pdb2reaction` in step 7 |
| `<HH:MM:SS>` | Your estimated walltime, capped by the queue's `resources_max.walltime` |

## Recipe: full one-shot probe

Paste this into the host once. It prints a discovery report to stdout; scheduler policy
and a representative pilot run are still needed to choose resource requests.

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

  echo
  echo "=== Architecture ==="
  uname -mrs
  command -v lscpu >/dev/null && lscpu | grep -E "^(Architecture|Model name|CPU\(s\)):"

  echo
  echo "=== GPU ==="
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>&1 || echo "no GPU"

  echo
  echo "=== CUDA modulefiles ==="
  command -v module >/dev/null && module avail cuda 2>&1 | head -20
  command -v nvcc >/dev/null && nvcc --version

  echo
  echo "=== PBS queues (if PBS) ==="
  if [[ "$SCHED" == pbs ]]; then
    command -v qstat >/dev/null && qstat -Q 2>/dev/null
    command -v pbsnodes >/dev/null && \
      pbsnodes -a 2>/dev/null | grep -E "^[^[:space:]]|^[[:space:]]+(np|properties|gpus|resources_available\.(ncpus|ngpus|mem)|totalmem)" | head -80
  fi

  echo
  echo "=== SLURM partitions (if SLURM) ==="
  if [[ "$SCHED" == slurm ]]; then
    command -v sinfo >/dev/null && sinfo -o "%P %l %c %m %N %G" 2>/dev/null
  fi

  echo
  echo "=== Conda envs with pdb2reaction ==="
  for env in $(conda env list 2>/dev/null | awk '/^[a-zA-Z]/{print $1}'); do
    conda run -n "$env" python -c \
      'import pdb2reaction; print("'"$env"':", pdb2reaction.__version__)' 2>/dev/null
  done
} 2>&1
```

Do not persist the raw report in a project/repository: it can contain private
hostnames, paths, scheduler policy, and environment names. If a file is needed
temporarily, create it under `${TMPDIR:-/tmp}` with mode `0600`, redact it, and
delete it after extracting the placeholder values.

## See also
- `pdb2reaction-hpc/SKILL.md` — uses `<YOUR_QUEUE>`, `<NCPU>`, `<NGPU>`,
  `<MEM>`, `<HH:MM:SS>`, `<CUDA_MODULE>`, `<YOUR_ENV>` placeholders.
- [`pdb2reaction-install-backends/env-cuda.md`](../pdb2reaction-install-backends/env-cuda.md) — uses the driver,
  GPU architecture, and official PyTorch matrix to choose a wheel; a toolkit module is only for source builds.
- [`pdb2reaction-install-backends/dft.md`](../pdb2reaction-install-backends/dft.md) — uses `uname -m` to decide
  between `gpu4pyscf-cuda12x` (x86_64) and CPU PySCF (aarch64).
