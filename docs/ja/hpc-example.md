# HPC 実行例: PBS + Open MPI + Ray

PBS と Open MPI の上で Ray のクラスターを立て、UMA のワーカーを複数のノードに広げるジョブスクリプトの例です。1 ノードのジョブには、この Ray の設定は要りません。

ワーカーは次の 2 つで指定します。

- `--uma-workers` — 全ノードにわたる、UMA のエネルギーと力を計算するワーカープロセスの総数。
- `--uma-workers-per-node` — そのうち各ノードで動作する数。ノードあたりの GPU / メモリ負荷を制御します。

**下のスクリプトはテンプレートとして扱ってください**: モジュール名、conda パス、ポート、PBS リソース要求は環境に合わせて調整が必要です。最後のコマンドの `test.pdb` と `-q -5 -m 1` は、自分のクラスターモデルとその総電荷・スピン多重度に置き換えてください。スクリプトは `run_p2r.pbs` として入力と同じディレクトリに保存し、そのディレクトリから `qsub run_p2r.pbs` で投入します。

```bash
#!/usr/bin/env bash
# ノードごとに MPI のランク 1 つと GPU 1 つを要求する。ncpus と ngpus は計算機に合わせる。
#PBS -l select=4:ncpus=72:mpiprocs=1:ngpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N pdb2reaction

set -euo pipefail
cd -- "${PBS_O_WORKDIR:?PBS_O_WORKDIR is unset}"

# --- 環境の設定 ---
source /etc/profile.d/modules.sh
module purge
P2R_CONDA_ENV=${P2R_CONDA_ENV:-p2r}         # インストールした環境の名前に変える
module load ompi                             # この MPI/Ray のテンプレートに必要な計算機のモジュール
# ビルド済みの PyTorch の wheel なら CUDA toolkit のモジュールは不要。拡張をソースからビルドするときだけ:
# module load <CUDA_MODULE> <COMPILER_MODULE>
source ~/apps/miniconda3/etc/profile.d/conda.sh
conda activate "${P2R_CONDA_ENV}"
# -------------------


# --- Ray の設定 ---
# CUDA と NCCL の動作を安定させる
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET

# 割り当てられていない GPU を選ばず、GPU が無ければ止める
if [[ -z "${CUDA_VISIBLE_DEVICES:-}" || "${CUDA_VISIBLE_DEVICES}" == "NoDevFiles" ]]; then
 echo "PBS allocation exposed no GPU (CUDA_VISIBLE_DEVICES is empty)." >&2
 exit 2
fi
export GPUS_PER_NODE="$(awk -F',' '{print NF}' <<< "${CUDA_VISIBLE_DEVICES}")"

# --- ノード ---
mapfile -t NODES < <(awk '!seen[$0]++' "$PBS_NODEFILE")
NNODES="${#NODES[@]}"
TOTAL_WORKERS=$((NNODES * GPUS_PER_NODE))

HEAD_NODE="${NODES[0]}"
HEAD_IP="$(getent ahostsv4 "${HEAD_NODE}" | awk 'NR==1{print $1}')"

# --- ポート（衝突を避けるため PBS_JOBID から決める） ---
JOBTAG="${PBS_JOBID%%.*}"
JOBNUM="${JOBTAG//[^0-9]/}"; JOBNUM="${JOBNUM:-0}"
PORT_SLOT=$((JOBNUM % 20))
BASE_PORT=$((20000 + PORT_SLOT * 1000))

RAY_PORT="${BASE_PORT}"
RAY_OBJECT_MANAGER_PORT=$((BASE_PORT + 1))
RAY_NODE_MANAGER_PORT=$((BASE_PORT + 2))
RAY_RUNTIME_ENV_AGENT_PORT=$((BASE_PORT + 3))
RAY_METRICS_EXPORT_PORT=$((BASE_PORT + 6))
RAY_MIN_WORKER_PORT=$((BASE_PORT + 100))
RAY_MAX_WORKER_PORT=$((BASE_PORT + 999))

RAY_TEMP_DIR="/tmp/ray_${JOBTAG}"
RAY_HEAD_ADDR="${HEAD_IP}:${RAY_PORT}"

# ray.init(address="auto") と ray status で使う
export RAY_ADDRESS="${RAY_HEAD_ADDR}"
# （任意。一時ファイルを多く使う計算で便利）
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
  # 下で独立したセッションとして起動した、このジョブのプロセスグループだけを止める。
  # `ray stop -f` をそのまま使うと共有ノード上のほかのジョブまで止めることがあるので使わない。
  kill -TERM -- "-${RAY_LAUNCH_PID}" >/dev/null 2>&1 || true
  wait "${RAY_LAUNCH_PID}" 2>/dev/null || true
 fi
}
trap cleanup EXIT

# ジョブごとのポートと一時ディレクトリで分け、ほかの Ray のジョブは止めない
"${MPI[@]}" "${BASH[@]}" "mkdir -p '${RAY_TEMP_DIR}'"
command -v setsid >/dev/null

# --- Ray の起動（rank0 がヘッド） ---
setsid "${MPI[@]}" "${BASH[@]}" "

# リモートのシェルの中でも同じ環境変数にそろえる
export PYTHONPATH='${PYTHONPATH:-}'
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export NCCL_SOCKET_FAMILY=AF_INET
export TMPDIR='${RAY_TEMP_DIR}'

# ノード間で hostid が同じときに NCCL が出す \"duplicate GPU\" を避ける
export NCCL_HOSTID=\$(hostname -s)

# リモートのランクでも、割り当てられていない GPU を作らない
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

# 割り当てた全ノードと全 GPU がそろうまで待つ（回数の上限つき）
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
# --- Ray の設定ここまで ---

pdb2reaction opt -i test.pdb -q -5 -m 1 \
 --uma-workers "${TOTAL_WORKERS}" --uma-workers-per-node "${GPUS_PER_NODE}"
```

ジョブのログに `[opt] Converged!` が出て、ジョブが終了コード 0 で終われば成功です。最適化した構造は `result_opt/final_geometry.xyz` にあります。

## ウォールタイム見積り

スクリプトの 24 時間は例として置いた上限で、実測の値ではありません。対象の計算機で短い試し計算を行い、時間の見積りを決めてください。

- **クラスターモデルの `opt` / `tsopt`**: 選んだバックエンドとモデル、Hessian の計算法、精度、収束設定で、代表的な構造の計算時間を測ります。
- **`pdb2reaction all` の通し実行**: extract → MEP → TS → IRC → freq → DFT の代表的なセグメント 1 つの時間を測ります。DFT の段は複数 GPU で SCF を回す仕組みではないため、GPU を増やしても速くなりません。
- **MEP（`path-search` / `path-opt`）**: 計算量は `--max-nodes`、最適化の反復回数、再帰セグメントの数とともに増えます。機構全体を見積る前に、セグメント 1 つの時間を測ります。

pdb2reaction は、最適化の各ステップ、経路の各イメージ、有限差分 Hessian の各変位を、1 構造ずつ UMA に渡します。そのためワーカーを増やして速くなりうるのは 1 回の評価だけで、複数の構造を同時に計算することはありません。ノードを増やす前に、ベンチマークで確かめてください。

## ジョブでの精度

精度はバックエンドと用途で選び、割り当てられた GPU で計算時間を測ってください。詳しくは {ref}`MLIP バックエンド › 精度 <ja-precision-by-gpu-class>` を参照してください。

## 使用上の注意点

* **1 ノードのジョブ**: ジョブスクリプトは `pdb2reaction` のコマンドを実行するだけで、GPU 1 つならデフォルトの `--uma-workers 1` のまま、そのノードの GPU を N 個使うなら `--uma-workers N --uma-workers-per-node N` を付けます。GPU 1 つのジョブの雛形は [`skills/pdb2reaction-hpc`](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/pdb2reaction-hpc/SKILL.md) にあります。
* **Ray のクラスターが立ち上がらないとき**: `opt` を始める前に終了コード 2 で止まります。
* **解析 Hessian はワーカー 1 つで使います**: `--uma-workers` を 2 以上にすると解析 Hessian は使えず、`--hessian-calc-mode Analytical` はエラーで止まります。`--uma-workers 1` にするか、`FiniteDifference` の Hessian を使ってください。
* **ワーカーを使うのは UMA だけです**: ORB / MACE / AIMNet2 では、1 以外の `workers` / `workers_per_node` は警告を出して無視されます。

## 関連ドキュメント

- [MLIP バックエンド](backends.md) — 設定の一覧と Hessian 評価モード
- [トラブルシューティング](troubleshooting.md) — ワーカーと解析 Hessian、GPU メモリなど、実行時のエラー
- [opt](opt.md) · [tsopt](tsopt.md) · [irc](irc.md) · [freq](freq.md) · [sp](sp.md) · [all](all.md) · [path-opt](path-opt.md) · [path-search](path-search.md) · [scan](scan.md) · [scan2d](scan2d.md) · [scan3d](scan3d.md) — `--uma-workers` / `--uma-workers-per-node` を取るサブコマンド
