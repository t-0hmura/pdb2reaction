# `uma_pysis` calculator

## Overview
`uma_pysis` exposes Meta's UMA machine-learning interatomic potentials to PySisyphus as an ASE-compatible calculator. It returns energies, forces, and Hessians (via analytical autograd or finite differences) in Hartree units while handling device placement, graph construction, and unit conversions internally. The calculator is used throughout `pdb2reaction` for optimization, path searches, thermochemistry, and trajectory post-processing.

## Quick start
```python
from pdb2reaction.uma_pysis import uma_pysis

# Build a calculator for a neutral singlet system on GPU when available
calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1", device="auto")

energy_only = calc.get_energy(["C", "O"], coords_bohr)
forces = calc.get_forces(["C", "O"], coords_bohr)
hessian = calc.get_hessian(["C", "O"], coords_bohr)
```
- Coordinates are supplied in **Bohr**; the wrapper converts to Å for UMA and converts energies/derivatives back to Hartree / Hartree·Bohr⁻¹ / Hartree·Bohr⁻².
- Attach the calculator to a `pysisyphus` geometry object or call it directly as above.

## Key features
- **UMA backend** – loads pretrained UMA checkpoints via FAIR-Chem's `pretrained_mlip` helpers and forwards charge/spin metadata in the AtomicData batch.
- **Device handling** – `device="auto"` selects CUDA when available, otherwise CPU. Graph construction happens on the chosen device; when `workers>1`, the parallel predictor manages device transfers.
- **Hessian modes** – `hessian_calc_mode="Analytical"` uses second-order autograd on the selected device; `"FiniteDifference"` (default) computes central differences of forces. Analytical mode is automatically disabled when multiple inference workers are requested.
- **Freeze atoms** – provide 0-based indices via `freeze_atoms`; frozen atoms receive zeroed forces. Hessians either drop frozen degrees of freedom (`return_partial_hessian=True`) or zero corresponding columns in the full matrix.
- **Precision control** – energies and forces are always returned as float64. Set `hessian_double=False` to obtain the Hessian in the model's native dtype (typically float32).
- **Multi-worker inference** – `workers>1` spawns FAIR-Chem's `ParallelMLIPPredictUnit` with `workers_per_node` workers per node, useful for batch throughput. Analytical Hessians are skipped in this mode.

## HPC example: PBS + Open MPI + Ray

`workers` / `workers_per_node` can be scaled across nodes by launching a Ray cluster under your scheduler. The following PBS script illustrates one way to build a multi-node Ray cluster on an Open MPI–equipped HPC system (adjust module names, ports, and resource requests to match your environment):

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
module load gcc ompi cuda/12.9 
source ~/apps/miniconda3/etc/profile.d/conda.sh
conda activate pdb2reaction
# -------------------


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
export NCCL_HOSTID=$(hostname -s)

# Per-node GPU count
if [[ -z \"${CUDA_VISIBLE_DEVICES:-}\" || \"${CUDA_VISIBLE_DEVICES}\" == \"NoDevFiles\" ]]; then
  export CUDA_VISIBLE_DEVICES=0
fi
GPUS=$(awk -F',' '{print NF}' <<<"${CUDA_VISIBLE_DEVICES}")

HOST=$(hostname -s)
IP=$(getent ahostsv4 "${HOST}" | awk 'NR==1{print $1}')

echo "[${HOST}] IP=${IP} CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES} (GPUS=${GPUS}) NCCL_HOSTID=${NCCL_HOSTID}"

if [[ \"${OMPI_COMM_WORLD_RANK:-0}\" == \"0\" ]]; then
  echo "[${HOST}] ray HEAD on ${HEAD_IP}:${RAY_PORT}"
  ray start --head --node-ip-address='${HEAD_IP}' --port='${RAY_PORT}' \
    --object-manager-port='${RAY_OBJECT_MANAGER_PORT}' --node-manager-port='${RAY_NODE_MANAGER_PORT}' \
    --runtime-env-agent-port='${RAY_RUNTIME_ENV_AGENT_PORT}' \
    --metrics-export-port='${RAY_METRICS_EXPORT_PORT}' \
    --min-worker-port='${RAY_MIN_WORKER_PORT}' --max-worker-port='${RAY_MAX_WORKER_PORT}' \
    --num-gpus="${GPUS}" \
    --temp-dir='${RAY_TEMP_DIR}' \
    --disable-usage-stats --include-dashboard=false --block
else
  until (echo > /dev/tcp/${HEAD_IP}/${RAY_PORT}) >/dev/null 2>&1; do sleep 1; done
  echo "[${HOST}] ray WORKER -> ${RAY_HEAD_ADDR}"
  ray start --address='${RAY_HEAD_ADDR}' --node-ip-address="${IP}" \
    --object-manager-port='${RAY_OBJECT_MANAGER_PORT}' --node-manager-port='${RAY_NODE_MANAGER_PORT}' \
    --runtime-env-agent-port='${RAY_RUNTIME_ENV_AGENT_PORT}' \
    --metrics-export-port='${RAY_METRICS_EXPORT_PORT}' \
    --min-worker-port='${RAY_MIN_WORKER_PORT}' --max-worker-port='${RAY_MAX_WORKER_PORT}' \
    --num-gpus="${GPUS}" \
    --temp-dir='${RAY_TEMP_DIR}' \
    --disable-usage-stats --block
fi
" &

RAY_LAUNCH_PID=$!

sleep 10 # Wait for workers
ray status || true

pdb2reaction opt -i test.pdb -q -5 -m 1
```

## Configuration reference
Common constructor keywords (defaults shown in the rightmost column):

| Option | Description | Default |
| --- | --- | --- |
| `charge` | Total system charge. | `0` |
| `spin` | Spin multiplicity (2S+1). | `1` |
| `model` | UMA pretrained model name. | `"uma-s-1p1"` |
| `task_name` | Task tag recorded in UMA batches. | `"omol"` |
| `device` | "cuda", "cpu", or automatic selection. | `"auto"` |
| `workers` / `workers_per_node` | Parallel UMA predictors; when `workers>1`, analytical Hessians are disabled. | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | Optional overrides for UMA neighborhood construction. | `None`, `None`, `False` |
| `freeze_atoms` | List of 0-based atom indices to freeze. | `None` |
| `hessian_calc_mode` | "Analytical" or "FiniteDifference" for Hessian evaluation. | `"FiniteDifference"` |
| `return_partial_hessian` | Return only the active-DOF Hessian block instead of the full matrix. | `False` |
| `hessian_double` | Assemble and return the Hessian in float64 precision. | `True` |
| `out_hess_torch` | Return Hessians as `torch.Tensor` objects. | `True` |

## CLI and YAML usage
`uma_pysis` is registered as a PySisyphus calculator entry point. With a YAML input file matching these keywords you can run:

```bash
uma_pysis input.yaml
```

Within `pdb2reaction` commands (e.g., `all`, `opt`, `path-opt`), calculator settings can be supplied via `--args-yaml` under `calc.kwargs` to reuse the same UMA configuration across stages.
