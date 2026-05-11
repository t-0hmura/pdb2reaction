---
name: pdb2reaction-hpc
description: PBS (Torque / PBSPro) and SLURM submission for pdb2reaction — preamble templates with placeholders, walltime budgeting, CPU vs GPU choice, job monitoring, and the dynamic-dispatch (flock + pbsdsh) recipe in `dynamic-dispatch.md`. TRIGGER on cluster submission / `qsub` / `sbatch` / walltime / preamble / `--resume` / multi-job dispatch / `workers > 1` questions. SKIP for local single-machine runs, install setup, or output parsing.
---

# pdb2reaction HPC

## Overview

Submit `pdb2reaction` as a PBS / SLURM job with 1 node / 1 GPU.
Placeholders filled from `pdb2reaction-env-detect/SKILL.md`.

## PBS preamble template (Torque / PBSPro)

```bash
#!/usr/bin/env bash
#PBS -N <jobname>
#PBS -q <YOUR_QUEUE>
#PBS -l nodes=1:ppn=<NCPU>:gpus=<NGPU>,mem=<MEM>GB,walltime=<HH:MM:SS>

set -euo pipefail
cd "${PBS_O_WORKDIR}"

# CUDA + toolchain: HPC modulefiles (env-detect outputs <CUDA_MODULE>)
# - gcc: load when the system default is too old for the CUDA toolkit or
#   when pip will compile a C/CUDA extension from source.
# - <OPENMPI_MODULE>: only for multi-node Ray (`workers > 1`); omit on
#   single-node jobs.
command -v module >/dev/null && module load <CUDA_MODULE> gcc <OPENMPI_MODULE>

# Conda env (env-detect outputs <YOUR_ENV>)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>

# Preflight (sanity-check env before the long job; output captured to pbs.out)
which pdb2reaction
pdb2reaction --version
python -c "import torch; print('cuda:', torch.cuda.is_available(), 'device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')"
nvidia-smi -L 2>/dev/null || true

# Optional: torch CUDA tuning
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all > pdb2reaction.log 2>&1
```

PBSPro syntax differs slightly (`#PBS -l select=1:ncpus=<NCPU>:ngpus=<NGPU>:mem=<MEM>gb`).
Both are accepted by most modern Torque + PBSPro installations; check
`man qsub` on your cluster.

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
command -v module >/dev/null && module load <CUDA_MODULE> gcc <OPENMPI_MODULE>
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>

# Preflight
which pdb2reaction
pdb2reaction --version
python -c "import torch; print('cuda:', torch.cuda.is_available(), 'device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')"
nvidia-smi -L 2>/dev/null || true

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

pdb2reaction all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all
```

## Walltime budgeting

Empirical rough-cuts for clusters of ~200–700 atoms with UMA-s-1.1
on a single mid-range GPU. Adjust generously.

| Stage | Per-segment time | Notes |
|---|---|---|
| `extract` | < 1 min | Pure Python, CPU |
| `path-search` (GSM) | 5–30 min | Scales with `--max-nodes` |
| `path-search` (DMF) | 10–60 min | Slower than GSM but more robust |
| `tsopt` (RS-I-RFO, default) | 5–60 min | Full-Hessian rebuilds dominate |
| `tsopt` (Dimer, alternative) | 1–10 min | Initial Hessian + dimer rotation, with periodic Hessian refreshes when curvature drifts; cheaper per cycle than full-Hessian RS-I-RFO when full-Hessian recomputation is too expensive |
| `irc` | 5–30 min | Forward + backward; `--max-cycles` cap 125 each direction |
| `freq` | 5–30 min | Hessian once + diagonalization |
| `dft` (ωB97M-V/def2-tzvpd, GPU) | 1–10 h | Heavy; `def2-svp` shrinks this 3–10× |
| `dft` (CPU) | 10–100× GPU time | Use only for small clusters (< 100 atoms); larger systems hit the upper bound |

For an `all` run with 2 segments + DFT: budget **6–24 h walltime**.
For pure MLIP `all` (no DFT): **2–6 h** is usually enough. For
benchmark-scale runs with **≥10 segments** (e.g. the 23-step
benchmark), a single `qsub` will exceed typical 24–72 h walltime caps;
fan out across multiple jobs (parallel `seg_NN/` `qsub`s, or the
`dynamic-dispatch.md` recipe).

## CPU vs GPU choice

| Workload | CPU | GPU |
|---|---|---|
| MLIP inference (any backend) | Possible but ~50–200× slower | **Required for production** |
| `pdb2reaction dft` with ωB97M-V on > 200 atoms | Slow (10–100 h) | Recommended |
| `pdb2reaction dft` with cheap functional / small molecule | Fine | Marginal speedup |
| Hessian (analytical, UMA) | OK if VRAM-limited | Faster |

Check [`pdb2reaction-install-backends/dft.md`](../pdb2reaction-install-backends/dft.md) for `--engine gpu` / `cpu`
specifics, including the aarch64 caveat (CPU PySCF only).

## Monitoring and control

PBS:

```bash
qstat -u "$USER"                 # state: Q (queued), R (running), C (complete)
qstat -f <jobid>                 # full job info
qdel <jobid>                     # cancel YOUR job
```

SLURM:

```bash
squeue -u "$USER"
scontrol show job <jobid>
scancel <jobid>
```

**`qdel` / `scancel` only kill jobs you own.** Don't kill jobs from
other users on shared clusters.

## Failed jobs / restart

`pdb2reaction all` supports `--resume` to continue from a partially
completed `--out-dir`: pass the same `-o <out_dir> --resume` and the
pipeline skips stages whose outputs already exist on disk. Note that
when extraction is skipped by `--resume`, charges must be supplied
explicitly via `-q/--charge` (and `--ligand-charge` when applicable),
since the original CLI invocation is not persisted.

If `--resume` cannot pick up the run (e.g., the partial output is in
an unexpected layout), individual stages support manual continuation
as a fallback:

- `tsopt`, `freq`, `irc`, `dft` — re-run on the previous output.
- `path-search` — pass the partial `mep.pdb` as `-i`.

For walltime-truncated jobs, write `--out-dir` to a persistent
location and re-submit with `--resume` after a longer walltime.

## Parallel job submission patterns

### Fan-out (one job per task)

```bash
for ts in seg_*.pdb; do
    jobid=$(qsub -v TS="$ts" generic_dft.sh)
    echo "submitted $ts as $jobid"
done
```

Each `qsub` produces an independent PBS job; the scheduler load-balances
them.

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
2. When `workers > 1`, requesting `hessian_calc_mode=Analytical` raises
   a `RuntimeError` (`pdb2reaction/backends/uma.py:160`); use
   `FiniteDifference` (the default) or drop to `workers = 1`.
   See `docs/uma-pysis.md` `(workers-fd-downgrade)=`.

The full PBS + OpenMPI + Ray bootstrap is in
`docs/hpc-example.md`; the schematic in this skill is single-node only.

## Useful environment variables

| Variable | Purpose |
|---|---|
| `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` | Reduce torch memory fragmentation |
| `CUDA_VISIBLE_DEVICES=0` | Restrict to a single GPU per worker |
| `OMP_NUM_THREADS=<NCPU>` | Limit OpenMP threads (avoid oversubscription) |
| `MKL_NUM_THREADS=<NCPU>` | Intel MKL thread cap |
| `LD_LIBRARY_PATH=<torch lib>:...` | Override system CUDA libs (see env-cuda.md) |

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
