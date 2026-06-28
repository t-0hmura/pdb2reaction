# HPC 実行例: PBS + Open MPI + Ray

大規模バッチや複数ノードの `pdb2reaction` 実行では、`workers` / `workers_per_node`（{ref}`MLIP 計算機 <ja-configuration-reference>` 参照）をスケジューラ配下の Ray クラスタでノード間に分散できます。

- `workers` — 全ノードを通じた UMA 予測プロセスの総数（デフォルト `1`）。
- `workers-per-node` — そのうち各ノードで動作する数（デフォルト `1`）。ノードあたりの GPU / メモリ負荷を制御します。

```{warning}
UMA バックエンドを `workers > 1` で実行している状態では `hessian_calc_mode="Analytical"` を明示指定しても解析ヘシアンは警告なく有限差分へダウングレードされます。解析ヘシアンが必要なら `workers = 1` に下げるか、デフォルトの `FiniteDifference` を使用してください。{ref}`ja-hessian-evaluation` を参照してください。ORB / MACE / AIMNet2 は `workers` / `workers_per_node` を受け付けないため、この規則は適用されません。
```

以下の PBS スクリプトは Open MPI を使用して複数ノードで Ray クラスタを構築する一例です。**テンプレートとして扱ってください**: モジュール名、conda パス、ポート、PBS リソース要求は環境に合わせて調整が必要です。

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
module load gcc ompi cuda/<your-version>     # 例: cuda/12.6 または cuda/12.9
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

## ウォールタイム見積り

上の 24 時間テンプレートはデフォルトの上限であり目標値ではありません。多くのジョブはこれを大きく下回って完了します。実行環境の wall-clock パターンに合わせて選んでください。

- **クラスターモデルの `opt` / `tsopt`**（~50–100 原子、単一 GPU）: 数分〜数時間
- **`pdb2reaction all` 一気通貫**（extract → MEP → TSOPT → IRC → freq → DFT、小型基質）: 通常数時間。ハイエンド multi-GPU ノードでは DFT 段が大きく短縮可能
- **MEP（`path-search` / `path-opt`）**: `--max-nodes`（セグメントあたりのイメージ数）と `--max-cycles`（GSM 最適化サイクル数）の双方でスケール。再帰的 `path-search` キャンペーンではこれにセグメント数が掛かるため、多段階反応では単一 GPU で何時間にも及び得る

UMA backend の場合、ウォールタイムは有効並列度（`workers` の総数）に概ね反比例します。ORB / MACE / AIMNet2 は worker 並列を持たないので、ノードを増やしても wall-clock は短くなりません。

## データセンター GPU での精度

これらのテンプレートが対象とする HPC データセンターカード（H100 / H200 / A100）では、本番の TS 最適化と Hessian を `--precision fp64` で実行してください。これらのカードでは fp64 のスループットコストが小さく、判断に足る品質で数値ノイズの少ない結果が得られます。スクリーニングやコンシューマー GPU（RTX 50xx / 40xx、fp64 が大幅に遅い）ではデフォルトの `--precision fp32` のままにします。振り分けの詳細と `--deterministic` との併用は {ref}`再現性: GPU クラスによる精度の選択 <ja-precision-by-gpu-class>` を参照してください。

## 関連項目

- [MLIP 計算機](uma-pysis.md) — 設定リファレンスとヘシアン評価モード
- [opt](opt.md) / [all](all.md) — `workers` / `workers_per_node` を取るサブコマンド
