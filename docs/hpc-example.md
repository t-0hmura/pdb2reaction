# HPC example: PBS + Open MPI + Ray

For large-batch or multi-node `pdb2reaction` runs, `workers` / `workers_per_node` (see {ref}`MLIP Calculator <configuration-reference>`) can be scaled across nodes by launching a Ray cluster under your scheduler.

- `workers` — total number of UMA predictor processes across all nodes (default `1`).
- `workers-per-node` — how many of those run on each node (default `1`); controls per-node GPU/memory pressure.

```{warning}
When you run the UMA backend with `workers > 1`, requesting `hessian_calc_mode="Analytical"` raises a `RuntimeError` (it is not silently downgraded). Drop to `workers = 1` if you need analytical Hessians, or use `FiniteDifference` (the default). See {ref}`hessian-evaluation`. ORB / MACE / AIMNet2 do not accept `workers` / `workers_per_node` and are unaffected by this rule.
```

The following PBS script illustrates one way to build a multi-node Ray cluster on an Open MPI–equipped HPC system. **Treat it as a template**: you will need to adjust module names, conda path, ports, and resource requests to match your environment.

```bash
#!/bin/bash
#PBS -l select=4:mpiprocs=72
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N pdb2reaction

cd "$PBS_O_WORKDIR"

# --- Environment setting ---
source /etc/profile.d/modules.sh
module purge
module load gcc ompi cuda/<your-version>     # e.g. cuda/12.6 or cuda/12.9
source ~/apps/miniconda3/etc/profile.d/conda.sh
conda activate <your-env>
# -------------------


# --- Ray setting ---
# Stable CUDA/NCCL
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET

# CUDA_VISIBLE_DEVICES fallback (if scheduler doesn't set)
if [[ -z "${CUDA_VISIBLE_DEVICES:-}" || "${CUDA_VISIBLE_DEVICES}" == "NoDevFiles" ]]; then
 export CUDA_VISIBLE_DEVICES=0
fi
export GPUS_PER_NODE="$(awk -F',' '{print NF}' <<< "${CUDA_VISIBLE_DEVICES}")"

# --- Nodes ---
mapfile -t NODES < <(awk '!seen[$0]++' "$PBS_NODEFILE")
NNODES="${#NODES[@]}"

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
 [[ -n "${RAY_LAUNCH_PID:-}" ]] && kill "${RAY_LAUNCH_PID}" >/dev/null 2>&1 || true
 "${MPI[@]}" "${BASH[@]}" "ray stop -f >/dev/null 2>&1 || true" || true
}
trap cleanup EXIT

# Prepare node-local /tmp + stop any leftover ray
"${MPI[@]}" "${BASH[@]}" "mkdir -p '${RAY_TEMP_DIR}'; ray stop -f >/dev/null 2>&1 || true"

# --- Launch Ray (rank0=head) ---
"${MPI[@]}" "${BASH[@]}" "

# Keep env stable inside remote shell as well
export PYTHONPATH='${PYTHONPATH}'
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET
export TMPDIR='${RAY_TEMP_DIR}'

# Avoid NCCL \"duplicate GPU\" when hostid is identical across nodes
export NCCL_HOSTID=\$(hostname -s)

# Per-node GPU count
if [[ -z \"\${CUDA_VISIBLE_DEVICES:-}\" || \"\${CUDA_VISIBLE_DEVICES}\" == \"NoDevFiles\" ]]; then
 export CUDA_VISIBLE_DEVICES=0
fi
GPUS=\$(awk -F',' '{print NF}' <<<\"\${CUDA_VISIBLE_DEVICES}\")

HOST=\$(hostname -s)
IP=\$(getent ahostsv4 \"\${HOST}\" | awk 'NR==1{print \$1}')

echo \"[\${HOST}] IP=\${IP} CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES} (GPUS=\${GPUS}) NCCL_HOSTID=\${NCCL_HOSTID}\"

if [[ \"\${OMPI_COMM_WORLD_RANK:-0}\" == \"0\" ]]; then
 echo \"[\${HOST}] ray HEAD on \${HEAD_IP}:\${RAY_PORT}\"
 ray start --head --node-ip-address='${HEAD_IP}' --port='${RAY_PORT}' \
 --object-manager-port='${RAY_OBJECT_MANAGER_PORT}' --node-manager-port='${RAY_NODE_MANAGER_PORT}' \
 --runtime-env-agent-port='${RAY_RUNTIME_ENV_AGENT_PORT}' \
 --metrics-export-port='${RAY_METRICS_EXPORT_PORT}' \
 --min-worker-port='${RAY_MIN_WORKER_PORT}' --max-worker-port='${RAY_MAX_WORKER_PORT}' \
 --num-gpus=\"\${GPUS}\" \
 --temp-dir='${RAY_TEMP_DIR}' \
 --disable-usage-stats --include-dashboard=false --block
else
 until (echo > /dev/tcp/\${HEAD_IP}/\${RAY_PORT}) >/dev/null 2>&1; do sleep 1; done
 echo \"[\${HOST}] ray WORKER -> \${RAY_HEAD_ADDR}\"
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

sleep 10 # Wait for workers
ray status || true
# --- Ray setup end ---

pdb2reaction opt -i test.pdb -q -5 -m 1 --workers ${NNODES} --workers-per-node ${GPUS_PER_NODE}
```

## Walltime budgeting

The 24 h template above is a default ceiling, not a target. Most jobs finish well under that; pick a budget that fits your system's wall-clock pattern:

- **Cluster-model `opt` / `tsopt`** (~50–100 atoms, single GPU): minutes to a few hours.
- **`pdb2reaction all` end-to-end** (extract → MEP → TSOPT → IRC → freq → DFT) on a small substrate: typically a few hours; high-end multi-GPU nodes can shorten the DFT stage substantially.
- **MEP (`path-search` / `path-opt`)**: scales with `--max-nodes` (images per segment) and `--max-cycles` (GSM optimizer iterations) — recursive `path-search` campaigns multiply both by segment count, so multistep mechanisms can occupy a single GPU for many hours.

Walltime scales roughly inversely with effective parallelism (the total `workers` count) on the UMA backend. ORB / MACE / AIMNet2 do not parallelize across workers, so adding more nodes does not shorten their wall-clock time.

## Precision on datacenter GPUs

On the HPC datacenter cards these templates target (H100 / H200 / A100), run production TS optimizations and Hessians with `--precision fp64`: the fp64 throughput cost is small on these cards and it gives deterministic-grade, low numerical-noise results. Keep the default `--precision fp32` for screening and on consumer GPUs (RTX 50xx / 40xx), where fp64 is substantially slower. See [Reproducibility → Choosing precision by GPU class](reproducibility.md#choosing-precision-by-gpu-class) for the routing details and the `--deterministic` pairing.

## See Also

- [MLIP Calculator](uma-pysis.md) — configuration reference and Hessian evaluation notes
- [opt](opt.md) / [all](all.md) — subcommands that honor `workers` / `workers_per_node`
