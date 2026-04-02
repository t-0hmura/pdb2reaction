# MLIP Calculator

## Overview
`pdb2reaction` supports multiple machine-learning interatomic potentials (MLIPs) as calculator backends for PySisyphus. The default backend is **UMA** (Meta's Universal Models for Atoms), but **ORB**, **MACE**, and **AIMNet2** are also available. Each backend returns energies, forces, and Hessian matrices in hartree-based atomic units while handling device placement and unit conversions internally. The calculator is used throughout `pdb2reaction` for optimization, path searches, thermochemistry, and trajectory post-processing.

## Quick start
```python
import numpy as np
from pdb2reaction.uma_pysis import uma_pysis

# Example: a neutral singlet diatomic on GPU when available
calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1", device="auto")

# uma_pysis expects coordinates in Bohr (shape: [n_atoms, 3])
coords_bohr = np.array([
 [0.0, 0.0, 0.0],
 [2.2, 0.0, 0.0], # ~1.16 Å
])

symbols = ["C", "O"]

# NOTE: These methods return dicts; extract values with the appropriate key
energy_h = calc.get_energy(symbols, coords_bohr)["energy"] # float (hartree)
forces_h_bohr = calc.get_forces(symbols, coords_bohr)["forces"] # ndarray (hartree/bohr)
hessian_h_bohr2 = calc.get_hessian(symbols, coords_bohr)["hessian"] # ndarray (hartree/bohr²)
```

- Coordinates are supplied in **bohr**; the wrapper converts to Angstrom for UMA and converts energies/derivatives back to hartree / hartree bohr⁻¹ / hartree bohr⁻².
- Attach the calculator to a `pysisyphus` geometry object or call it directly as above.

## Backend selection

Select a backend with `-b/--backend` on any command, or set `calc.backend` in YAML:

```bash
# UMA (default)
pdb2reaction opt -i input.pdb -q 0

# ORB
pdb2reaction opt -i input.pdb -q 0 -b orb

# MACE
pdb2reaction opt -i input.pdb -q 0 -b mace

# AIMNet2
pdb2reaction opt -i input.pdb -q 0 -b aimnet2
```

| Backend | Install | Analytical Hessian | Multi-worker | Notes |
|---------|---------|-------------------|-------------|-------|
| **UMA** | included | Yes | Yes | Full feature set via fairchem |
| **ORB** | `pip install "pdb2reaction[orb]"` | No (FD only) | No | orb-models |
| **MACE** | `pip install "pdb2reaction[mace]"` | No (FD only) | No | mace-torch |
| **AIMNet2** | `pip install "pdb2reaction[aimnet2]"` | No (FD only) | No | aimnet |

### Implicit solvent correction

All backends support xTB-based implicit solvent corrections via `--solvent`:

```bash
pdb2reaction opt -i input.pdb -q 0 --solvent water
pdb2reaction opt -i input.pdb -q 0 -b orb --solvent water --solvent-model cpcmx
```

The correction uses a delta approach: ΔE = E_xTB(solvent) - E_xTB(vacuum), added to the MLIP energy/forces/Hessian. Requires `xtb` to be installed and accessible on `PATH`.

## Key features
- **MLIP backends** – the default UMA backend loads pretrained UMA checkpoints via FAIR-Chem's `pretrained_mlip` helpers and forwards charge/spin metadata in the AtomicData batch. Alternative backends (ORB, MACE, AIMNet2) are available via `-b/--backend`.
- **Device handling** – `device="auto"` selects CUDA when available, otherwise CPU. Graph construction happens on the chosen device; when `workers>1`, the parallel predictor manages device transfers.
### Hessian evaluation

`hessian_calc_mode="Analytical"` uses second-order autograd on the selected device; `"FiniteDifference"` (default) computes central differences of forces. Analytical mode is automatically disabled when multiple inference workers are requested.
- **Freeze atoms** – provide 1-based indices via `freeze_atoms`; frozen atoms receive zeroed forces. The Hessian matrix either drops frozen degrees of freedom (`return_partial_hessian=True`) or zeroes the corresponding rows and columns in the full matrix.
- **Precision control** – energies and forces are always returned as float64. Set `hessian_double=False` to obtain the Hessian matrix in the model's native dtype (typically float32).
- **Multi-worker inference** – `workers>1` spawns FAIR-Chem's `ParallelMLIPPredictUnit` with `workers_per_node` workers per node, useful for batch throughput. **Warning:** when `workers>1`, analytical Hessians are silently switched to finite differences (`force_fd=True`) even if `hessian_calc_mode="Analytical"` is set. No warning is printed — check your logs if Hessian timings are unexpectedly long.

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
# --- Ray setup end ---

pdb2reaction opt -i test.pdb -q -5 -m 1
```

## Configuration reference
Common constructor keywords (defaults shown in the rightmost column):

| Option | Description | Default |
| --- | --- | --- |
| `backend` | MLIP backend engine. | `"uma"` |
| `charge` | Total system charge. | `0` |
| `spin` | Spin multiplicity (2S+1). | `1` |
| `model` | UMA pretrained model name (`uma-s-1p1`, `uma-m-1p1`). | `"uma-s-1p1"` |
| `task_name` | Task tag recorded in UMA batches. | `"omol"` |
| `device` | "cuda", "cpu", or automatic selection. | `"auto"` |
| `workers` / `workers_per_node` | Parallel UMA predictors; when `workers>1`, analytical Hessians are disabled. | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | Optional overrides for UMA neighborhood construction. | `None`, `None`, `False` |
| `freeze_atoms` | List of 1-based atom indices to freeze. | _None_ |
| `hessian_calc_mode` | "Analytical" or "FiniteDifference" for Hessian evaluation. | `"FiniteDifference"` |
| `return_partial_hessian` | Return only the active-DOF Hessian block instead of the full matrix. | `True` |
| `hessian_double` | Assemble and return the Hessian in float64 precision. | `True` |
| `out_hess_torch` | Return Hessians as `torch.Tensor` objects. | `True` |
| `print_timing` | Print Hessian computation timing breakdown. | `True` |
| `solvent` | Implicit solvent name (e.g. `"water"`) or `"none"`. | `"none"` |
| `solvent_model` | xTB solvent model: `"alpb"` or `"cpcmx"`. | `"alpb"` |


---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [opt](opt.md) -- Single-structure optimization using an MLIP backend
- [path-opt](path-opt.md) -- MEP optimization with MLIP backend
- [all](all.md) -- End-to-end workflow using MLIP across stages
