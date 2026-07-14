# HPC example: PBS + Open MPI + Ray

For large-batch or multi-node `pdb2reaction` runs, `workers` / `workers_per_node` (see {ref}`MLIP Calculator <configuration-reference>`) can be scaled across nodes by launching a Ray cluster under your scheduler.

- `workers` — total number of UMA predictor processes across all nodes (default `1`).
- `workers-per-node` — how many of those run on each node (default `1`); controls per-node GPU/memory pressure.

```{warning}
When you run the UMA backend with `workers > 1`, requesting `hessian_calc_mode="Analytical"` raises `BackendError` (a `RuntimeError` subclass) because the parallel predictor exposes no autograd model. Use `workers = 1` for an analytical Hessian, or select `FiniteDifference`. ORB / MACE / AIMNet2 do not accept `workers` / `workers_per_node` and are unaffected by this rule. See {ref}`hessian-evaluation`.
```

The following PBS script illustrates one way to build a multi-node Ray cluster on an Open MPI–equipped HPC system. **Treat it as a template**: you will need to adjust module names, conda path, ports, and resource requests to match your environment.

```bash
#!/usr/bin/env bash
# One MPI rank and one allocated GPU per node; adapt ncpus/ngpus to the site.
#PBS -l select=4:ncpus=72:mpiprocs=1:ngpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N pdb2reaction

set -euo pipefail
cd -- "${PBS_O_WORKDIR:?PBS_O_WORKDIR is unset}"

# --- Environment setting ---
source /etc/profile.d/modules.sh
module purge
CUDA_MODULE=${CUDA_MODULE:-cuda/12.6}        # change to the site's module name
P2R_CONDA_ENV=${P2R_CONDA_ENV:-p2r}         # change to the installed environment
module load gcc ompi "${CUDA_MODULE}"
source ~/apps/miniconda3/etc/profile.d/conda.sh
conda activate "${P2R_CONDA_ENV}"
# -------------------


# --- Ray setting ---
# Stable CUDA/NCCL
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET

# Fail closed: never select an unallocated device.
if [[ -z "${CUDA_VISIBLE_DEVICES:-}" || "${CUDA_VISIBLE_DEVICES}" == "NoDevFiles" ]]; then
 echo "PBS allocation exposed no GPU (CUDA_VISIBLE_DEVICES is empty)." >&2
 exit 2
fi
export GPUS_PER_NODE="$(awk -F',' '{print NF}' <<< "${CUDA_VISIBLE_DEVICES}")"

# --- Nodes ---
mapfile -t NODES < <(awk '!seen[$0]++' "$PBS_NODEFILE")
NNODES="${#NODES[@]}"
TOTAL_WORKERS=$((NNODES * GPUS_PER_NODE))

HEAD_NODE="${NODES[0]}"
HEAD_IP="$(getent ahostsv4 "${HEAD_NODE}" | awk 'NR==1{print $1}')"

# --- Ports (avoid collisions: derive from PBS_JOBID) ---
JOBTAG="${PBS_JOBID%%.*}"
JOBNUM="${JOBTAG//[^0-9]/}"; JOBNUM="${JOBNUM:-0}"
BASE_PORT=$((20000 + (JOBNUM % 20000)))

RAY_PORT="${BASE_PORT}"
RAY_OBJECT_MANAGER_PORT=$((BASE_PORT + 1))
RAY_NODE_MANAGER_PORT=$((BASE_PORT + 2))
RAY_RUNTIME_ENV_AGENT_PORT=$((BASE_PORT + 3))
RAY_METRICS_EXPORT_PORT=$((BASE_PORT + 6))
RAY_MIN_WORKER_PORT=$((BASE_PORT + 100))
RAY_MAX_WORKER_PORT=$((BASE_PORT + 999))

RAY_TEMP_DIR="/tmp/ray_${JOBTAG}"
RAY_HEAD_ADDR="${HEAD_IP}:${RAY_PORT}"

# For ray.init(address="auto") / ray status
export RAY_ADDRESS="${RAY_HEAD_ADDR}"
# (optional but handy for tmp-heavy workloads)
export TMPDIR="${RAY_TEMP_DIR}"

echo "Nodes(${NNODES}): ${NODES[*]}"
echo "Ray head: ${RAY_HEAD_ADDR}"
echo "Ray temp: ${RAY_TEMP_DIR}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES} (GPUS_PER_NODE=${GPUS_PER_NODE})"

MPI=(mpirun --bind-to none -np "${NNODES}" --map-by ppr:1:node)
BASH=(bash --noprofile --norc -c)

cleanup() {
 echo "Stopping Ray..."
 if [[ -n "${RAY_LAUNCH_PID:-}" ]]; then
  # The launcher was started in its own session below. Stop only this job's
  # process group; a bare `ray stop -f` can kill another job on a shared node.
  kill -TERM -- "-${RAY_LAUNCH_PID}" >/dev/null 2>&1 || true
  wait "${RAY_LAUNCH_PID}" 2>/dev/null || true
 fi
}
trap cleanup EXIT

# Unique ports and temp paths isolate this job; never stop unrelated Ray jobs.
"${MPI[@]}" "${BASH[@]}" "mkdir -p '${RAY_TEMP_DIR}'"
command -v setsid >/dev/null

# --- Launch Ray (rank0=head) ---
setsid "${MPI[@]}" "${BASH[@]}" "

# Keep env stable inside remote shell as well
export PYTHONPATH='${PYTHONPATH:-}'
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET
export TMPDIR='${RAY_TEMP_DIR}'

# Avoid NCCL \"duplicate GPU\" when hostid is identical across nodes
export NCCL_HOSTID=\$(hostname -s)

# Per-node GPU count. Never fabricate a device in the remote rank.
if [[ -z \"\${CUDA_VISIBLE_DEVICES:-}\" || \"\${CUDA_VISIBLE_DEVICES}\" == \"NoDevFiles\" ]]; then
 echo \"[\$(hostname -s)] no scheduler-visible GPU\" >&2
 exit 2
fi
GPUS=\$(awk -F',' '{print NF}' <<<\"\${CUDA_VISIBLE_DEVICES}\")

HOST=\$(hostname -s)
IP=\$(getent ahostsv4 \"\${HOST}\" | awk 'NR==1{print \$1}')

echo \"[\${HOST}] IP=\${IP} CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES} (GPUS=\${GPUS}) NCCL_HOSTID=\${NCCL_HOSTID}\"

if [[ \"\${OMPI_COMM_WORLD_RANK:-0}\" == \"0\" ]]; then
 echo \"[\${HOST}] ray HEAD on ${HEAD_IP}:${RAY_PORT}\"
 ray start --head --node-ip-address='${HEAD_IP}' --port='${RAY_PORT}' \
 --object-manager-port='${RAY_OBJECT_MANAGER_PORT}' --node-manager-port='${RAY_NODE_MANAGER_PORT}' \
 --runtime-env-agent-port='${RAY_RUNTIME_ENV_AGENT_PORT}' \
 --metrics-export-port='${RAY_METRICS_EXPORT_PORT}' \
 --min-worker-port='${RAY_MIN_WORKER_PORT}' --max-worker-port='${RAY_MAX_WORKER_PORT}' \
 --num-gpus=\"\${GPUS}\" \
 --temp-dir='${RAY_TEMP_DIR}' \
 --disable-usage-stats --include-dashboard=false --block
else
 connected=0
 for _attempt in \$(seq 1 120); do
  if (echo > /dev/tcp/${HEAD_IP}/${RAY_PORT}) >/dev/null 2>&1; then connected=1; break; fi
  sleep 1
 done
 (( connected == 1 )) || { echo \"Ray head did not become reachable\" >&2; exit 2; }
 echo \"[\${HOST}] ray WORKER -> ${RAY_HEAD_ADDR}\"
 ray start --address='${RAY_HEAD_ADDR}' --node-ip-address=\"\${IP}\" \
 --object-manager-port='${RAY_OBJECT_MANAGER_PORT}' --node-manager-port='${RAY_NODE_MANAGER_PORT}' \
 --runtime-env-agent-port='${RAY_RUNTIME_ENV_AGENT_PORT}' \
 --metrics-export-port='${RAY_METRICS_EXPORT_PORT}' \
 --min-worker-port='${RAY_MIN_WORKER_PORT}' --max-worker-port='${RAY_MAX_WORKER_PORT}' \
 --num-gpus=\"\${GPUS}\" \
 --temp-dir='${RAY_TEMP_DIR}' \
 --disable-usage-stats --block
fi
" &

RAY_LAUNCH_PID=$!

# Bounded readiness gate: require every allocated node and GPU resource.
export EXPECTED_RAY_NODES="${NNODES}"
export EXPECTED_RAY_GPUS="${TOTAL_WORKERS}"
ready=0
for _attempt in $(seq 1 120); do
 if ! kill -0 "${RAY_LAUNCH_PID}" 2>/dev/null; then
  echo "Ray launcher exited before the cluster became ready." >&2
  break
 fi
 if python - <<'PY'
import os
import ray

ray.init(address="auto", ignore_reinit_error=True, logging_level="ERROR")
live_nodes = sum(bool(node.get("Alive")) for node in ray.nodes())
gpu_total = float(ray.cluster_resources().get("GPU", 0.0))
expected_nodes = int(os.environ["EXPECTED_RAY_NODES"])
expected_gpus = float(os.environ["EXPECTED_RAY_GPUS"])
ray.shutdown()
if live_nodes < expected_nodes or gpu_total < expected_gpus:
    raise SystemExit(1)
PY
 then
  ready=1
  break
 fi
 sleep 2
done
(( ready == 1 )) || { echo "Ray readiness timed out." >&2; exit 2; }
ray status
# --- Ray setup end ---

pdb2reaction opt -i test.pdb -q -5 -m 1 \
 --workers "${TOTAL_WORKERS}" --workers-per-node "${GPUS_PER_NODE}"
```

## Walltime budgeting

The 24 h template above is an example ceiling, not a measured target. Pick a budget from representative pilots on the target stack:

- **Cluster-model `opt` / `tsopt`**: time the selected backend/model, Hessian mode, precision, and convergence settings on a representative structure.
- **`pdb2reaction all` end-to-end** (extract → MEP → TSOPT → IRC → freq → DFT): time a representative segment. The bundled DFT command is not a general multi-GPU SCF driver, so requesting more GPUs does not make that stage scale automatically.
- **MEP (`path-search` / `path-opt`)**: cost grows with `--max-nodes`, optimizer iterations, and recursive segment count; measure one representative segment before budgeting the full mechanism.

UMA `workers` can improve inference throughput for suitable large/batched work,
but optimizer stages contain sequential work and do not scale inversely with
worker count. Benchmark before reserving more nodes. ORB / MACE / AIMNet2 do
not use pdb2reaction's UMA worker pool.

## Precision in scheduled jobs

Choose precision by backend and purpose, then measure its cost on the allocated
GPU. Leaving it unset preserves the backend defaults (UMA/AIMNet2 fp32,
ORB/MACE fp64). Explicit fp32 is useful for screening but downgrades ORB/MACE
and their finite-difference Hessians; it is not a final-validation setting.
Precision does not imply determinism or prove a first-order saddle. See
[Reproducibility → Choosing precision by backend and purpose](reproducibility.md#choosing-precision-by-backend-and-purpose).

## See Also

- [MLIP Calculator](uma-pysis.md) — configuration reference and Hessian evaluation notes
- [opt](opt.md) / [all](all.md) — subcommands that honor `workers` / `workers_per_node`
