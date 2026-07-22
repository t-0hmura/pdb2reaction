# HPC 実行例: PBS + Open MPI + Ray

大規模バッチや複数ノードの `pdb2reaction` 実行では、`workers` / `workers_per_node`（{ref}`MLIP 計算機 <ja-configuration-reference>` 参照）をスケジューラ配下の Ray クラスタでノード間に分散できます。

- `workers` — 全ノードにわたる UMA 予測プロセスの総数（デフォルト `1`）。
- `workers-per-node` — そのうち各ノードで動作する数（デフォルト `1`）。ノードあたりの GPU / メモリ負荷を制御します。

```{warning}
UMA バックエンドを `workers > 1` で実行している状態で `hessian_calc_mode="Analytical"` を指定すると、並列 predictor が autograd model を持たないため `BackendError`（`RuntimeError` のサブクラス）で停止します。解析 Hessian には `workers = 1`、並列実行には `FiniteDifference` を指定してください。{ref}`ja-hessian-evaluation` を参照してください。ORB / MACE / AIMNet2 は `workers` / `workers_per_node` を受け付けないため、この規則は適用されません。
```

以下の PBS スクリプトは Open MPI を使用して複数ノードで Ray クラスタを構築する一例です。**テンプレートとして扱ってください**: モジュール名、conda パス、ポート、PBS リソース要求は環境に合わせて調整が必要です。

```bash
#!/usr/bin/env bash
# nodeごとにMPI rank 1、GPU 1を明示要求。ncpus/ngpusはsiteに合わせる。
#PBS -l select=4:ncpus=72:mpiprocs=1:ngpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N pdb2reaction

set -euo pipefail
cd -- "${PBS_O_WORKDIR:?PBS_O_WORKDIR is unset}"

# --- Environment setting ---
source /etc/profile.d/modules.sh
module purge
P2R_CONDA_ENV=${P2R_CONDA_ENV:-p2r}         # install 済み env 名に変更
module load ompi                             # この MPI/Ray テンプレートで必要な site module
# prebuilt PyTorch wheel に CUDA toolkit module は不要。source build extension の場合のみ:
# module load <CUDA_MODULE> <COMPILER_MODULE>
source ~/apps/miniconda3/etc/profile.d/conda.sh
conda activate "${P2R_CONDA_ENV}"
# -------------------


# --- Ray setting ---
# Stable CUDA/NCCL
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET

# 未割当GPUを勝手に選ばずfail closedにする
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
  # 下で独立sessionとして起動した、このjobのprocess groupだけを停止する。
  # bare `ray stop -f` は共有node上の別jobまで停止し得るため使わない。
  kill -TERM -- "-${RAY_LAUNCH_PID}" >/dev/null 2>&1 || true
  wait "${RAY_LAUNCH_PID}" 2>/dev/null || true
 fi
}
trap cleanup EXIT

# job固有のport/tmp pathで分離し、他のRay jobは停止しない
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

# remote rankでも未割当GPUを捏造しない
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

# 全node/GPU resourceを確認するbounded readiness gate
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

## ウォールタイム見積り

上の24時間templateは例示の上限であり、実測targetではありません。対象stackで代表pilotを行ってbudgetを決めます。

- **cluster modelの `opt` / `tsopt`**: 選択backend/model、Hessian mode、precision、収束設定を代表構造で実測
- **`pdb2reaction all` 一気通貫**（extract → MEP → TSOPT → IRC → freq → DFT）: 代表 segment を実測してください。同梱 DFT command は一般的な multi-GPU SCF driver ではないため、GPU request を増やしても自動的には scale しません。
- **MEP（`path-search` / `path-opt`）**: costは `--max-nodes`、optimizer iteration、再帰segment数で増えます。全mechanismのbudget前に代表segmentを実測

UMA `workers` は大規模/batched workload の inference throughput を改善し得ますが、optimizer stage には逐次処理があり、worker 数に反比例して短縮するわけではありません。node を増やす前に benchmark してください。ORB / MACE / AIMNet2 は pdb2reaction の UMA worker pool を使用しません。

## scheduler job での精度

精度は backend と用途で選び、割り当て GPU 上の cost を実測します。未指定なら tested default（UMA/AIMNet2 fp32、ORB/MACE fp64）を保ちます。明示的 fp32 は screening には使えますが ORB/MACE の既定を下げ、その finite-difference Hessian を最終検証には使えません。precision は determinism も first-order saddle も保証しません。{ref}`backend と用途による精度の選択 <ja-precision-by-gpu-class>` を参照してください。

## 関連項目

- [MLIP 計算機](uma-pysis.md) — 設定リファレンスとHessian評価モード
- [opt](opt.md) / [all](all.md) — `workers` / `workers_per_node` を取るサブコマンド
