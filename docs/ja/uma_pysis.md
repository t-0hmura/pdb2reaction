# `uma_pysis` 計算機

## 概要
`uma_pysis` は、MetaのUMA機械学習ポテンシャルをPySisyphus向けのASE互換計算機として公開します。エネルギー/力/ヘシアン（解析自動微分または有限差分）を Hartree 単位で返し、デバイス配置・グラフ構築・単位変換を内部で処理します。`pdb2reaction` では最適化、経路探索、熱化学、軌跡後処理など広範に利用されます。

## クイックスタート
```python
import numpy as np
from pdb2reaction.uma_pysis import uma_pysis

# 例: 中性一重項の2原子系（GPUが利用可能ならGPU、なければCPU）
calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1", device="auto")

# uma_pysis には Bohr 単位の座標（形状: [n_atoms, 3]）を渡します
coords_bohr = np.array([
    [0.0, 0.0, 0.0],
    [2.2, 0.0, 0.0],  # 約 1.16 Å
])

symbols = ["C", "O"]

energy_h = calc.get_energy(symbols, coords_bohr)
forces_h_bohr = calc.get_forces(symbols, coords_bohr)
hessian_h_bohr2 = calc.get_hessian(symbols, coords_bohr)
```

- 座標は **Bohr** で与えます。ラッパー内部で Å に変換し、UMA計算後に Hartree / Hartree·Bohr⁻¹ / Hartree·Bohr⁻² に戻します。
- `pysisyphus` の geometry オブジェクトにアタッチするか、上記のように直接呼び出せます。

## 主な特徴
- **UMAバックエンド** – FAIR-Chem の `pretrained_mlip` ヘルパーでUMAチェックポイントを読み込み、AtomicData バッチに電荷/スピン情報を付与。
- **デバイス処理** – `device="auto"` はCUDAがあればGPU、なければCPUを選択。グラフ構築は選択デバイス上で行い、`workers>1` では並列予測器が転送を管理。
- **ヘシアンモード** – `hessian_calc_mode="Analytical"` で2階自動微分、`"FiniteDifference"`（デフォルト）は力の中心差分。`workers>1` の場合は解析ヘシアンは無効化されます。
- **凍結原子** – `freeze_atoms` に0始まりの原子インデックスを渡すと、凍結原子の力がゼロ化。`return_partial_hessian=True` で凍結自由度を除いたヘシアンを返すか、フル行列で該当行/列をゼロ化できます。
- **精度制御** – エネルギー/力は常にfloat64。`hessian_double=False` でヘシアンをモデルのネイティブdtype（通常float32）で返します。
- **マルチワーカー推論** – `workers>1` で FAIR-Chem の `ParallelMLIPPredictUnit` を起動し、`workers_per_node` をノードごとに指定可能。解析ヘシアンはこのモードでは無効です。

## HPC例: PBS + Open MPI + Ray

`workers` / `workers_per_node` は、スケジューラ環境でRayクラスタを立ち上げてスケールできます。以下はOpen MPIを使うPBSスクリプトの一例です（モジュール名、ポート、リソース要求は環境に合わせて調整してください）。

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

## 設定リファレンス
代表的なコンストラクタ引数（右端はデフォルト値）:

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `charge` | 総電荷 | `0` |
| `spin` | スピン多重度（2S+1） | `1` |
| `model` | UMAプリトレイン済みモデル名 | `"uma-s-1p1"` |
| `task_name` | UMAバッチに記録されるタスクタグ | `"omol"` |
| `device` | `"cuda"` / `"cpu"` / `"auto"` | `"auto"` |
| `workers` / `workers_per_node` | 並列UMA予測器（`workers>1` で解析ヘシアン無効） | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | 近傍構築のオプション上書き | `None`, `None`, `False` |
| `freeze_atoms` | 0始まりの凍結原子インデックス | _None_ |
| `hessian_calc_mode` | `"Analytical"` または `"FiniteDifference"` | `"FiniteDifference"` |
| `return_partial_hessian` | アクティブ自由度のみ返す | `False` |
| `hessian_double` | ヘシアンをfloat64で返す | `True` |
| `out_hess_torch` | ヘシアンを `torch.Tensor` で返す | `True` |

## CLI / YAML での利用
`uma_pysis` はPySisyphusの計算機エントリポイントとして登録されています。以下のようにYAML 入力で起動できます。

```bash
uma_pysis input.yaml
```

`pdb2reaction` の各コマンド（`all`, `opt`, `path-opt` など）では、`--args-yaml` の `calc` キー配下に同等の設定を渡すことで、同一のUMA設定を再利用できます。
