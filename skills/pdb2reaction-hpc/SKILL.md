---
name: pdb2reaction-hpc
description: PBS (Torque / PBSPro) and SLURM submission for pdb2reaction — preamble templates with placeholders, walltime budgeting, CPU vs GPU choice, job monitoring, and the dynamic-dispatch (flock + pbsdsh) recipe in `dynamic-dispatch.md`. TRIGGER on cluster submission / `qsub` / `sbatch` / walltime / preamble / multi-job dispatch / multiple predictor worker questions. SKIP for local single-machine runs, install setup, or output parsing.
---

# pdb2reaction HPC

## Purpose

Submit `pdb2reaction` as a PBS / SLURM job with 1 node / 1 GPU.
Placeholders filled from `pdb2reaction-env-detect/SKILL.md`.

## PBS preamble templates

Torque and PBSPro resource syntax is not interchangeable. Select the syntax
shown by working jobs or `man qsub` on the target cluster; do not put both forms
in one script.

### Torque-style resources

```bash
#!/usr/bin/env bash
#PBS -N <jobname>
#PBS -q <YOUR_QUEUE>
#PBS -l nodes=1:ppn=<NCPU>:gpus=<NGPU>,mem=<MEM>GB,walltime=<HH:MM:SS>
#PBS -j oe
# Leave -o unset so PBS uses its job-ID-qualified default output name.

set -euo pipefail
cd "${PBS_O_WORKDIR}"

# CUDA + toolchain: HPC modulefiles (env-detect outputs <CUDA_MODULE>)
# - gcc: load when the system default is too old for the CUDA toolkit or
#   when pip will compile a C/CUDA extension from source.
# - <OPENMPI_MODULE>: only for multi-node Ray (`workers > 1`); omit on
#   single-node jobs.
# Replace every placeholder. Remove modules the site/runtime does not require.
command -v module >/dev/null && module load <CUDA_MODULE>
# module load <COMPILER_MODULE>    # only for locally built extensions
# module load <OPENMPI_MODULE>     # only for the site's external multi-node Ray setup

# Conda env (env-detect outputs <YOUR_ENV>)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>

# Fail-closed GPU preflight before the long job. For an intentionally CPU-only
# template, remove these two checks together with the GPU resource request.
which pdb2reaction
pdb2reaction --version
python -c "import sys, torch; ok=torch.cuda.is_available(); print('cuda:', ok, 'device:', torch.cuda.get_device_name(0) if ok else 'unavailable'); sys.exit(0 if ok else 2)"
nvidia-smi -L

# Optional: torch CUDA tuning
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

RUN_ID=${PBS_JOBID:-manual_$(date +%Y%m%dT%H%M%S)_$$}
RUN_ID=${RUN_ID//[^A-Za-z0-9_.-]/_}
pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all > "pdb2reaction.${RUN_ID}.log" 2>&1
```

### PBSPro-style resources

Use the same body as above, but replace the Torque resource line with:

```bash
#PBS -l select=1:ncpus=<NCPU>:ngpus=<NGPU>:mem=<MEM>gb
#PBS -l walltime=<HH:MM:SS>
```

## SLURM preamble template

```bash
#!/usr/bin/env bash
#SBATCH --job-name=<jobname>
#SBATCH --partition=<YOUR_PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<NCPU>
#SBATCH --gres=gpu:<NGPU>
#SBATCH --mem=<MEM>G
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=%x.%j.out
set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
# CUDA + toolchain (see PBS template above for when gcc / OpenMPI are needed)
command -v module >/dev/null && module load <CUDA_MODULE>
# module load <COMPILER_MODULE>    # only for locally built extensions
# module load <OPENMPI_MODULE>     # only for the site's external multi-node Ray setup
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>

# Fail-closed GPU preflight. Remove only for an explicitly CPU-only template.
which pdb2reaction
pdb2reaction --version
python -c "import sys, torch; ok=torch.cuda.is_available(); print('cuda:', ok, 'device:', torch.cuda.get_device_name(0) if ok else 'unavailable'); sys.exit(0 if ok else 2)"
nvidia-smi -L

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all
```

## Walltime budgeting

Do not infer a production walltime from atom count alone. Optimizer
convergence, backend/model, Hessian mode, path segmentation, and hardware
dominate. Use the stages below to construct a measured budget:

| Stage | Pilot measurement | Scaling risk |
|---|---|---|
| `extract` | Run once locally/login node if site policy allows | Mostly file parsing; usually minor |
| `path-search` / `path-opt` | Pilot one representative segment | Nodes × cycles × recursive segment count; convergence is nonlinear |
| `tsopt` | Pilot the chosen optimizer on one candidate | Hessian builds/rebuilds and recovery branches dominate |
| `irc` | Pilot both directions | Each direction has its own cycle cap and may stop differently |
| `freq` | Time one Hessian with the production mode | Dense matrix memory grows quadratically with active degrees of freedom |
| `dft` (GPU) | Run one representative single point | Basis, functional, grid, elements, and VRAM dominate |
| `dft` (CPU) | Benchmark one representative structure | Cost depends on basis, functional, grid, elements, and CPU/GPU hardware |

Sum the measured stage costs for the expected number of segments, then add
recovery/retry headroom. Do not assume UMA walltime falls inversely with
`workers`: geometry optimization has sequential work and worker overhead.
For many segments, fan out independent work (parallel `seg_NN/` jobs or the
`dynamic-dispatch.md` recipe) within site policy.

## CPU vs GPU choice

| Workload | CPU | GPU |
|---|---|---|
| MLIP inference (any backend) | Supported but often substantially slower; benchmark the actual model/system | Recommended for production-scale workloads |
| `pdb2reaction dft` | Supported with `--engine cpu`; pilot actual method/system | GPU recommended when the GPU4PySCF stack supports the requested calculation |
| Analytical Hessian (supported built-in backend) | Can avoid GPU memory pressure but may be impractical | Autograd can use the accelerator, but speed and peak memory are backend/model/system dependent; pilot it |

Check [`pdb2reaction-install-backends/dft.md`](../pdb2reaction-install-backends/dft.md) for `--engine gpu` / `cpu`
specifics. On aarch64 the packaged extra provides CPU PySCF; only a separately
source-built and locally validated GPU4PySCF environment can enable GPU DFT.

## Monitoring and control

PBS:

```bash
qstat -u "$USER"                 # state: Q (queued), R (running), C (complete)
qstat -f <jobid>                 # full job info
qdel <jobid>                     # only after the verification below
```

SLURM:

```bash
squeue -u "$USER"
scontrol show job <jobid>
scancel <jobid>                  # only after the verification below
```

Before cancellation, inspect the exact job ID and verify its owner, job name,
cluster/server, and submitted command or working directory in `qstat -f` /
`scontrol show job`. Scheduler authorization differs by site and administrator
accounts may be able to cancel other users' jobs; ownership is not a substitute
for this check. Never cancel by a guessed ID or broad name match.

## Failed jobs / restart

If a job fails or is killed by walltime, individual stages can be
re-run on the partial output by invoking the standalone subcommands
directly:

- `tsopt`, `freq`, `irc`, `dft` — re-run on the previous output.
- `path-search` — has no resume; restart it with two or more structures in
  reaction order (`-i A.pdb -i B.pdb …`; the CLI requires at least two `-i`
  paths). Reuse the original endpoints, or split the partial multi-model
  `mep.pdb` / `mep_trj.xyz` trajectory into single-structure files and pass
  those frames as separate `-i` flags.

For walltime-truncated `all` runs, point `--out-dir` at a persistent
location and pick up where you left off by chaining the appropriate
subcommands against the artifacts that `all` already produced.

## Parallel job submission patterns

### Fan-out (one job per task)

```bash
shopt -s nullglob
tasks=(seg_*.pdb)
(( ${#tasks[@]} > 0 )) || { echo "No seg_*.pdb inputs; submitting nothing." >&2; exit 2; }
MAX_SUBMIT=${MAX_SUBMIT:-100}
(( ${#tasks[@]} <= MAX_SUBMIT )) || {
    echo "Refusing ${#tasks[@]} submissions (MAX_SUBMIT=$MAX_SUBMIT)." >&2
    exit 2
}
for ts in "${tasks[@]}"; do
    if jobid=$(qsub -v TS="$ts" generic_dft.sh); then
        echo "submitted $ts as $jobid"
    else
        echo "qsub failed for $ts; stopping fan-out." >&2
        exit 2
    fi
done
```

Each `qsub` produces an independent PBS job; the scheduler load-balances
them. For a larger collection, use a site-supported throttled job array or
the bounded dynamic-dispatch recipe below instead of an unbounded submit loop.

### Dynamic dispatch (one job, N nodes pull tasks)

When you have many short tasks and want to amortize the queue wait,
use the flock + pbsdsh pattern documented in `dynamic-dispatch.md`. One
qsub grabs N nodes, each node runs a worker that pulls tasks from a
shared list with file-lock-protected counter increment.

## Multi-node MLIP inference (`workers > 1`, UMA only)

Most subcommands that touch geometry expose `--workers` and
`--workers-per-node`, which spin up a Ray cluster of UMA predictor
workers. **Two important caveats:**

1. The `workers` / `workers_per_node` flags are filtered to the UMA
   backend (see `pdb2reaction/backends/__init__.py:_BACKEND_ACCEPTED_KEYS`).
   ORB / MACE / AIMNet2 silently drop the kwarg.
2. With UMA, `workers > 1` plus an explicit
   `hessian_calc_mode=Analytical` request raises `BackendError` (a
   `RuntimeError` subclass) before the requested method can be changed; use
   `FiniteDifference` (the default) or drop to `workers = 1`.
   See `docs/uma-pysis.md` `(workers-analytical-error)=`.

The full PBS + OpenMPI + Ray bootstrap is in
`docs/hpc-example.md`; the schematic in this skill is single-node only.

## Useful environment variables

| Variable | Purpose |
|---|---|
| `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` | Reduce torch memory fragmentation |
| `CUDA_VISIBLE_DEVICES` | Normally leave the scheduler-provided mapping unchanged. Set it manually only outside scheduler isolation or as part of a tested worker-launch scheme; device indices inside a job are local to that mapping. |
| `OMP_NUM_THREADS=<NCPU>` | Limit OpenMP threads (avoid oversubscription) |
| `MKL_NUM_THREADS=<NCPU>` | Intel MKL thread cap |
| `LD_LIBRARY_PATH` | Leave unchanged for prebuilt wheels unless a diagnosed site-specific module/build requires it; see env-cuda.md |

## ssh-based remote submission

Generally avoided (depends on per-user ssh config). If your cluster requires `ssh <login> qsub`,
add that as a wrapper around the PBS / SLURM command above; do **not**
embed it inside the skill template.

## See also

- `dynamic-dispatch.md` — flock + pbsdsh template for many short tasks.
- `pdb2reaction-env-detect/SKILL.md` — discover queue / module / env
  values for the placeholders above.
- [`pdb2reaction-install-backends/env-cuda.md`](../pdb2reaction-install-backends/env-cuda.md) — driver / torch CUDA
  pairing.
- [`pdb2reaction-cli/all.md`](../pdb2reaction-cli/all.md) — the typical workload submitted to HPC.
