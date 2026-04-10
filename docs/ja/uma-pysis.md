# MLIP 計算機

## 概要
`pdb2reaction` は複数の機械学習原子間ポテンシャル（MLIP）を PySisyphus 向けの計算機バックエンドとしてサポートします。デフォルトバックエンドは **UMA**（Meta の Universal Models for Atoms）ですが、**ORB**、**MACE**、**AIMNet2** も利用可能です。各バックエンドはエネルギー/力/ヘシアンを Hartree 単位で返し、デバイス配置・単位変換を内部で処理します。`pdb2reaction` の最適化、経路探索、熱化学、軌跡後処理など広範に利用されます。

## クイックスタート
```python
import numpy as np
from pdb2reaction.uma_pysis import uma_pysis

# 例: 中性一重項の2原子系（GPUが利用可能ならGPU、なければCPU）
calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1", device="auto")

# uma_pysis には Bohr 単位の座標（形状: [n_atoms, 3]）を渡します
coords_bohr = np.array([
 [0.0, 0.0, 0.0],
 [2.2, 0.0, 0.0], # 約 1.16 Å
])

symbols = ["C", "O"]

# 注: これらのメソッドは dict を返すため、適切なキーで値を取り出します
energy_h = calc.get_energy(symbols, coords_bohr)["energy"] # float (Hartree)
forces_h_bohr = calc.get_forces(symbols, coords_bohr)["forces"] # ndarray (Hartree/Bohr)
hessian_h_bohr2 = calc.get_hessian(symbols, coords_bohr)["hessian"] # ndarray (Hartree/Bohr²)
```

- 座標は **Bohr** で与えます。ラッパー内部で Å に変換し、UMA計算後に Hartree / Hartree·Bohr⁻¹ / Hartree·Bohr⁻² に戻します。
- `pysisyphus` の geometry オブジェクトにアタッチするか、上記のように直接呼び出せます。

## Python API: Calculator Factory

`backends` モジュールは MLIP 計算機をプログラムから生成するファクトリを提供します:

```python
from pdb2reaction.backends import create_calculator, create_ase_calculator
```

| Function | Description |
|----------|-------------|
| `create_calculator(backend="uma", **kwargs)` | Create a PySisyphus-compatible MLIP calculator. Accepts `charge`, `spin`, `model`, `device`, `solvent`, `solvent_model`, `hessian_calc_mode`, `freeze_atoms`, and other backend-specific kwargs. Unknown keys are silently filtered per-backend. |
| `create_ase_calculator(backend="uma", **kwargs)` | Create an ASE-compatible MLIP calculator (used for DMF workflows and ASE-based tools). Same kwargs as `create_calculator`. |

### 例

```python
from pdb2reaction.backends import create_calculator

# 解析ヘシアン付き UMA 計算機
calc = create_calculator(
    backend="uma",
    charge=0,
    spin=1,
    device="auto",
    hessian_calc_mode="Analytical",
)

# 暗黙溶媒付き ORB 計算機
calc_orb = create_calculator(
    backend="orb",
    charge=-1,
    spin=1,
    solvent="water",
    solvent_model="alpb",
)
```

返される計算機は pysisyphus の計算機インターフェース（`get_energy`, `get_forces`, `get_hessian`）を実装しています。すべてのメソッドは `(atoms: List[str], coords: np.ndarray)` を受け取り、coords は **Bohr** 単位です。戻り値は `"energy"`（Hartree）、`"forces"`（Hartree/Bohr）、`"hessian"`（Hartree/Bohr²）を含む dict です。

## バックエンド選択

任意のコマンドで `-b/--backend` を指定するか、YAML で `calc.backend` を設定します:

```bash
# UMA（デフォルト）
pdb2reaction opt -i input.pdb -q 0

# ORB
pdb2reaction opt -i input.pdb -q 0 -b orb

# MACE
pdb2reaction opt -i input.pdb -q 0 -b mace

# AIMNet2
pdb2reaction opt -i input.pdb -q 0 -b aimnet2
```

| バックエンド | インストール | 解析ヘシアン | マルチワーカー | 備考 |
|---------|---------|-------------------|-------------|-------|
| **UMA** | 同梱 | あり | あり | fairchem による完全機能 |
| **ORB** | `pip install "pdb2reaction[orb]"` | なし（有限差分のみ） | なし | orb-models |
| **MACE** | `pip install --no-deps mace-torch` | なし（有限差分のみ） | なし | mace-torch |
| **AIMNet2** | `pip install "pdb2reaction[aimnet]"` | なし（有限差分のみ） | なし | aimnet |

### 暗黙溶媒補正

すべてのバックエンドで xTB ベースの暗黙溶媒補正が `--solvent` により利用可能です:

```bash
pdb2reaction opt -i input.pdb -q 0 --solvent water
pdb2reaction opt -i input.pdb -q 0 -b orb --solvent water --solvent-model cpcmx
```

デルタアプローチによる補正: ΔE = E_xTB(溶媒) - E_xTB(真空) を MLIP エネルギー/力/ヘシアンに加算します。`xtb` が `PATH` 上にインストールされている必要があります。

## 主な特徴
- **MLIPバックエンド** – デフォルトの UMA バックエンドは FAIR-Chem の `pretrained_mlip` ヘルパーでUMAチェックポイントを読み込み、AtomicData バッチに電荷/スピン情報を付与。代替バックエンド（ORB、MACE、AIMNet2）は `-b/--backend` で利用可能。
- **デバイス処理** – `device="auto"` はCUDAがあればGPU、なければCPUを選択。グラフ構築は選択デバイス上で行い、`workers>1` では並列予測器が転送を管理。
### ヘシアンモード

`hessian_calc_mode="Analytical"` で2階自動微分、`"FiniteDifference"`（デフォルト）は力の中心差分。`workers>1` の場合は解析ヘシアンは無効化されます。
- **凍結原子** – `freeze_atoms` に1始まりの原子インデックスを渡すと、凍結原子の力がゼロ化。`return_partial_hessian=True` で凍結自由度を除いたヘシアンを返すか、フル行列で該当行/列をゼロ化できます。
- **精度制御** – エネルギー/力は常にfloat64。`hessian_double=False` でヘシアンをモデルのネイティブdtype（通常float32）で返します。
- **マルチワーカー推論** – `workers>1` で FAIR-Chem の `ParallelMLIPPredictUnit` を起動し、`workers_per_node` をノードごとに指定可能。バッチ処理速度の向上に有効です。**注意:** `workers>1` の場合、`hessian_calc_mode="Analytical"` を指定していても解析ヘシアンは暗黙的に有限差分（`force_fd=True`）へ切り替わります。警告は出力されないため、ヘシアン計算時間が想定より長い場合はログを確認してください。

## HPC での使用例: PBS + Open MPI + Ray

`workers` / `workers_per_node` は、スケジューラ配下で Ray クラスタを構築することでノード間に分散できます。以下は Open MPI を使う PBS スクリプトの一例です（モジュール名、ポート、リソース要求は環境に合わせて調整してください）。

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

## 設定リファレンス
代表的なコンストラクタ引数（右端はデフォルト値）:

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `backend` | MLIP バックエンドエンジン | `"uma"` |
| `charge` | 総電荷 | `0` |
| `spin` | スピン多重度（2S+1） | `1` |
| `model` | UMAモデル名 (`uma-s-1p1`, `uma-m-1p1`) | `"uma-s-1p1"` |
| `task_name` | UMAバッチに記録されるタスクタグ | `"omol"` |
| `device` | `"cuda"` / `"cpu"` / `"auto"` | `"auto"` |
| `workers` / `workers_per_node` | 並列UMA予測器（`workers>1` で解析ヘシアン無効） | `1` / `1` |
| `max_neigh`, `radius`, `r_edges` | 近傍構築のオプション上書き | `None`, `None`, `False` |
| `freeze_atoms` | 1始まりの凍結原子インデックス | _None_ |
| `hessian_calc_mode` | `"Analytical"` または `"FiniteDifference"` | `"FiniteDifference"` |
| `return_partial_hessian` | アクティブ自由度のみ返す | `True` |
| `hessian_double` | ヘシアンをfloat64で返す | `True` |
| `out_hess_torch` | ヘシアンを `torch.Tensor` で返す | `True` |
| `print_timing` | ヘシアン計算のタイミング内訳を表示 | `True` |
| `solvent` | 暗黙溶媒名（例: `"water"`）または `"none"` | `"none"` |
| `solvent_model` | xTB 溶媒モデル: `"alpb"` または `"cpcmx"` | `"alpb"` |


---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [opt](opt.md) -- MLIP バックエンドを使う単一構造最適化
- [path-opt](path-opt.md) -- MLIP バックエンドを使う MEP 最適化
- [all](all.md) -- MLIP を複数段で使う一気通貫ワークフロー
