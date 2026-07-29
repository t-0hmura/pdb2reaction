# Dynamic dispatch: flock + pbsdsh (dynamic-dispatch.md)

Use one multi-node PBS allocation for many independent, short tasks. One
worker runs per allocated node and claims numbered tasks through `flock`.
Prefer a scheduler-native job array when the site supports it.

This pattern requires `pbsdsh`, a shared directory visible at the same path on
every node, and a filesystem on which `flock` is reliable. Confirm all three
site-specific assumptions first. Each task's `run.sh` must be **idempotent** or
publish outputs atomically: recovery is deliberately at-least-once, so a task
killed after publishing output but before its completion marker may run again.

Fill `<N_NODES>`, `<NCPU>`, `<NGPU>`, `<MEM>`, `<HH:MM:SS>`, and
`<YOUR_QUEUE>` from the environment detection skill. Set the shell defaults
`TASK_LIST_FILE`, `CONDA_SH`, and `P2R_CONDA_ENV` below to the verified local
values. Leave `CUDA_MODULE` empty for prebuilt wheels; set it only when a
locally built extension requires the site's toolkit module.

## State model

Every submission uses `.pdb2reaction_dispatch/<PBS_JOBID>/`; concurrent jobs
cannot overwrite one another. The dispatcher copies the task list into that
directory and records task IDs under `inflight/`, `completed/`, `failed/`, and
`interrupted/`. A retry list is derived from **every ID without a completion
marker**, not merely from commands that returned nonzero. This also recovers a
claim abandoned by node loss, a signal, or walltime termination.

If PBS kills the whole allocation before the parent epilogue runs, execute the
preserved `recover.sh` manually with that run directory. Do not resume only
from `next_task`; that counter records claims, not completions.

## `dispatcher.sh`

```bash
#!/usr/bin/env bash
#PBS -N pdb2reaction_dispatch
#PBS -q <YOUR_QUEUE>
# PBSPro form (replace with the site's verified syntax)
#PBS -l select=<N_NODES>:ncpus=<NCPU>:ngpus=<NGPU>:mem=<MEM>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -j oe

set -euo pipefail
cd "${PBS_O_WORKDIR}"

ROOT_DIR="$PWD"
TASK_LIST_FILE=${TASK_LIST_FILE:-tasks.txt}
CUDA_MODULE=${CUDA_MODULE:-}
CONDA_SH=${CONDA_SH:-${HOME}/apps/miniconda3/etc/profile.d/conda.sh}
P2R_CONDA_ENV=${P2R_CONDA_ENV:-p2r}
TASK_LIST="${ROOT_DIR}/${TASK_LIST_FILE}"
RUN_ID=${PBS_JOBID:-manual_$(date +%Y%m%dT%H%M%S)_$$}
RUN_ID=${RUN_ID//[^A-Za-z0-9_.-]/_}
STATE_DIR="${ROOT_DIR}/.pdb2reaction_dispatch/${RUN_ID}"
mkdir -p "${ROOT_DIR}/.pdb2reaction_dispatch"
if ! mkdir "${STATE_DIR}"; then
    echo "State directory already exists; refusing to collide: ${STATE_DIR}" >&2
    exit 2
fi
mkdir "${STATE_DIR}/inflight" "${STATE_DIR}/completed" \
      "${STATE_DIR}/failed" "${STATE_DIR}/interrupted"

[[ -s "${TASK_LIST}" ]] || { echo "Missing/empty task list: ${TASK_LIST}" >&2; exit 2; }
if grep -nE '^[[:space:]]*$' "${TASK_LIST}" >&2; then
    echo "Task list contains a blank line; remove it before submission." >&2
    exit 2
fi
command -v flock >/dev/null || { echo "flock not found" >&2; exit 2; }
command -v pbsdsh >/dev/null || { echo "pbsdsh not found" >&2; exit 2; }
[[ -r "${CONDA_SH}" ]] || { echo "Conda init script not readable: ${CONDA_SH}" >&2; exit 2; }
[[ -n "${PBS_NODEFILE:-}" && -s "${PBS_NODEFILE}" ]] || {
    echo "PBS_NODEFILE is missing/empty; this template requires a PBS allocation." >&2
    exit 2
}

cp -- "${TASK_LIST}" "${STATE_DIR}/tasks.txt"
TASK_SNAPSHOT="${STATE_DIR}/tasks.txt"
TOTAL_TASKS=$(awk 'END {print NR}' "${TASK_SNAPSHOT}")
printf '1\n' > "${STATE_DIR}/next_task"
: > "${STATE_DIR}/lock"

cat > "${STATE_DIR}/recover.sh" <<'RECOVER_EOF'
#!/usr/bin/env bash
set -euo pipefail
STATE_DIR=$(cd "$(dirname "$0")" && pwd)
TASK_LIST="${STATE_DIR}/tasks.txt"
TOTAL_TASKS=$(awk 'END {print NR}' "${TASK_LIST}")
retry_tmp="${STATE_DIR}/retry_tasks.txt.tmp.$$"
: > "${retry_tmp}"
for ((task_id=1; task_id<=TOTAL_TASKS; task_id++)); do
    if [[ ! -f "${STATE_DIR}/completed/${task_id}" ]]; then
        sed -n "${task_id}p" "${TASK_LIST}" >> "${retry_tmp}"
    fi
done
mv -- "${retry_tmp}" "${STATE_DIR}/retry_tasks.txt"
incomplete=$(awk 'END {print NR+0}' "${STATE_DIR}/retry_tasks.txt")
echo "Recovery: completed=$((TOTAL_TASKS-incomplete))/${TOTAL_TASKS}; retry=${incomplete}"
echo "Retry list: ${STATE_DIR}/retry_tasks.txt"
(( incomplete == 0 ))
RECOVER_EOF
chmod +x "${STATE_DIR}/recover.sh"

cat > "${STATE_DIR}/worker.sh" <<'WORKER_EOF'
#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$1"
STATE_DIR="$2"
TOTAL_TASKS="$3"
CONDA_SH="$4"
P2R_CONDA_ENV="$5"
CUDA_MODULE="$6"
TASK_LIST="${STATE_DIR}/tasks.txt"
LOCK_FILE="${STATE_DIR}/lock"

if command -v module >/dev/null && [[ -n "${CUDA_MODULE}" ]]; then
    module load "${CUDA_MODULE}"
fi
source "${CONDA_SH}"
conda activate "${P2R_CONDA_ENV}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# This template requests one or more GPUs per worker node. Fail closed instead
# of silently running a long task on CPU when a module/device is unavailable.
python -c "import sys, torch; ok=torch.cuda.is_available(); print('cuda:', ok, 'device:', torch.cuda.get_device_name(0) if ok else 'unavailable'); sys.exit(0 if ok else 2)"
nvidia-smi -L

active_task=""
record_interrupt() {
    if [[ -n "${active_task}" && ! -f "${STATE_DIR}/completed/${active_task}" ]]; then
        printf '%s\n' "$(hostname) $(date -Is)" > "${STATE_DIR}/interrupted/${active_task}"
    fi
}
trap 'record_interrupt' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM HUP

failures=0
while true; do
    task_id=""
    {
        flock -x 200
        current=$(cat "${STATE_DIR}/next_task")
        if (( current <= TOTAL_TASKS )); then
            task_id="${current}"
            printf '%s\n' "$((current + 1))" > "${STATE_DIR}/next_task"
            printf '%s %s\n' "$(hostname)" "$(date -Is)" > "${STATE_DIR}/inflight/${task_id}"
        fi
    } 200>"${LOCK_FILE}"
    [[ -n "${task_id}" ]] || break
    active_task="${task_id}"
    rel_path=$(sed -n "${task_id}p" "${TASK_LIST}")

    case "${rel_path}" in
        /*|..|../*|*/../*|*/..)
            echo "Invalid task path outside ROOT_DIR: ${rel_path}" >&2
            printf '%s\n' "invalid path" > "${STATE_DIR}/failed/${task_id}"
            rm -f "${STATE_DIR}/inflight/${task_id}"
            active_task=""
            failures=$((failures + 1))
            continue
            ;;
    esac
    task_dir="${ROOT_DIR}/${rel_path}"
    if [[ ! -d "${task_dir}" || ! -f "${task_dir}/run.sh" ]]; then
        echo "Missing task directory or run.sh: ${rel_path}" >&2
        printf '%s\n' "missing run.sh" > "${STATE_DIR}/failed/${task_id}"
        rm -f "${STATE_DIR}/inflight/${task_id}"
        active_task=""
        failures=$((failures + 1))
        continue
    fi

    echo "[$(date -Is)] $(hostname): START ${task_id} ${rel_path}"
    if (cd "${task_dir}" && bash run.sh); then
        : > "${STATE_DIR}/completed/${task_id}"
        rm -f "${STATE_DIR}/inflight/${task_id}" \
              "${STATE_DIR}/failed/${task_id}" \
              "${STATE_DIR}/interrupted/${task_id}"
        active_task=""
        echo "[$(date -Is)] $(hostname): DONE  ${task_id} ${rel_path}"
    else
        rc=$?
        printf 'exit=%s\n' "${rc}" > "${STATE_DIR}/failed/${task_id}"
        rm -f "${STATE_DIR}/inflight/${task_id}"
        active_task=""
        failures=$((failures + 1))
        echo "[$(date -Is)] $(hostname): ERROR ${task_id} ${rel_path} (exit=${rc})" >&2
    fi
done
(( failures == 0 ))
WORKER_EOF
chmod +x "${STATE_DIR}/worker.sh"

mapfile -t vnode_indices < <(awk '!seen[$1]++ {print NR-1}' "${PBS_NODEFILE}")
pids=()
for idx in "${vnode_indices[@]}"; do
    pbsdsh -n "${idx}" -- "${STATE_DIR}/worker.sh" \
        "${ROOT_DIR}" "${STATE_DIR}" "${TOTAL_TASKS}" \
        "${CONDA_SH}" "${P2R_CONDA_ENV}" "${CUDA_MODULE}" &
    pids+=("$!")
done

worker_status=0
for pid in "${pids[@]}"; do
    wait "${pid}" || worker_status=1
done

recovery_status=0
"${STATE_DIR}/recover.sh" || recovery_status=1
if (( worker_status != 0 || recovery_status != 0 )); then
    echo "Dispatcher incomplete. Preserve ${STATE_DIR} and resubmit its retry_tasks.txt." >&2
    exit 1
fi
echo "Dispatcher complete: ${TOTAL_TASKS}/${TOTAL_TASKS}; state kept at ${STATE_DIR}."
```

## Task-list format and recovery

The task list contains one non-empty path, relative to the submission
directory, per line. Each directory contains its own `run.sh`:

```text
dft_R/
dft_TS1/
dft_IM/
dft_P/
seg_03/dft_TS/
```

After a complete or interrupted run:

```bash
STATE_DIR=".pdb2reaction_dispatch/1234.server"  # replace with the exact job ID
"${STATE_DIR}/recover.sh"
cat "${STATE_DIR}/retry_tasks.txt"
```

Verify the exact job ID before selecting a state directory. Resubmit the retry
file as a new task list. Keep the old directory as a run record; a new PBS
job receives a different run directory.

## SLURM analog

Prefer an array because the scheduler owns task identity and retry state:

```bash
#SBATCH --job-name=pdb2reaction_array
#SBATCH --array=1-<N_TASKS>%<MAX_CONCURRENT>
#SBATCH ...

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
rel_path=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "<TASK_LIST_FILE>")
cd "${rel_path}" && bash run.sh
```

Use dynamic dispatch only when many independent tasks are short enough that
individual queue startup dominates. For dependent tasks, use a workflow
engine rather than this counter pattern.

## Cross-references

- `SKILL.md` — single-job PBS / SLURM templates and cancellation safety.
- `pdb2reaction-env-detect/SKILL.md` — how to fill site-specific placeholders.
