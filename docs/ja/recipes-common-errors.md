# 典型エラー別レシピ

症状は分かるものの、どのコマンドページから見ればよいか迷うときの入口として使ってください。
詳細は [トラブルシューティング](troubleshooting.md) を並行して参照してください。

## 早見表

各行の「詳細」列は [トラブルシューティング](troubleshooting.md) の該当セクションへ直接リンクします。

| 症状 | 最初にやること | 詳細（セクション） |
| --- | --- | --- |
| **入力 / 抽出** | | |
| 元素列欠落で抽出が止まる | 元の PDB に `add-elem-info` を適用してください | {ref}`入力 / 抽出の問題 <ts-input-extraction>` |
| `[multi] Atom count mismatch` / `[multi] Atom order mismatch` | 同じ前処理ツール・設定で全 PDB を再生成。最初に原子順序を固定したら並べ替えない | {ref}`入力 / 抽出の問題 <ts-input-extraction>` |
| **電荷 / スピン** | | |
| `-q/--charge is required` 系エラー | `-q/--charge` または `-l/--ligand-charge` を明示指定してください | {ref}`電荷 / スピンの問題 <ts-charge-spin>` |
| 計算は通るが状態/エネルギーが不自然 | [CLI 規約](cli-conventions.md) の電荷解決順序を再確認してください | {ref}`電荷 / スピンの問題 <ts-charge-spin>` |
| **計算 / 収束** | | |
| UMA の `--workers > 1` と `--hessian-calc-mode Analytical` の併用で `BackendError` | 解析 Hessian には `--workers 1`、並列実行には `FiniteDifference` を指定 | {ref}`workers > 1 のHessian制約 <ja-workers-analytical-error>` |
| 実行時に CUDA OOM | `--radius` を縮小して再抽出（extract / all のみ）、`--opt-mode grad` に切替、有限差分 Hessian のまま、または VRAM の大きい GPU へ | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| TS は収束したが小さい虚振動が複数残る | `--flatten` を追加（`tsopt`、`opt`、`pdb2reaction all` 共通） | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| TSOPT が収束しない | L-BFGS/Dimer: `max_step` を**縮小**。RFO/RS-P-RFO（`tsopt` デフォルト）/RS-I-RFO: `trust_radius`/`trust_min`/`trust_max` を**縮小**。サイクル上限を増やし、TS 品質を確認 | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| IRC が正常に終了しない | 単独: `--step-size` / `--max-cycles`。`all`: `--irc-step-size` / YAML `irc.max_cycles`（`all --max-cycles` は MEP 用）。開始構造の虚振動が 1 本か確認 | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| エネルギー平坦化後に opt/TSOPT が `stalled` で停止する | 未収束として扱い、最終構造と力・ステップ条件を確認してから、適切な閾値またはoptimizer設定で再実行 | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| MEP 探索（GSM/DMF）が失敗 | `--max-nodes` をデフォルト 20 から増やす、`--preopt` 有効化（デフォルト: `all`/`path-search`/`path-opt` で `True`、`scan*` で `False`）、別の `--mep-mode` を試す | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| **インストール / 環境** | | |
| DMF モードの import エラー（`cyipopt`）、または `No module named 'dmf'` | `conda install -c conda-forge cyipopt`（`pydmf` は `pdb2reaction` の依存として自動インストールされる） | {ref}`インストール / 環境の問題 <ts-install-env>` |
| UMA モデルで 401/403 / gated repo エラー | `hf auth login` でログインし、UMA モデルのライセンスに同意してください | {ref}`インストール / 環境の問題 <ts-install-env>` |
| `e3nn` / `fairchem-core` の import 競合（UMA env に MACE を入れた） | MACE 専用 conda env を使用（`mace-torch` は `e3nn==0.4.4` を pin し、`fairchem-core` の `e3nn>=0.5` と共存不可）。`pip uninstall -y fairchem-core && pip install mace-torch` | {ref}`インストール / 環境の問題 <ts-install-env>` |
| CUDA/GPU 実行時エラー | `torch.cuda.is_available()` と CUDA バージョンの整合を確認してください | {ref}`インストール / 環境の問題 <ts-install-env>` |
| ORB / AIMNet2 / MACE の import エラー | ORB: `pip install "pdb2reaction[orb]"`。AIMNet2: `pip install "pdb2reaction[aimnet]"`。MACE: 上の専用 env 手順 | {ref}`インストール / 環境の問題 <ts-install-env>` |
| 図の出力失敗 | `plotly_get_chrome -y` で Chrome ランタイムを導入してください | {ref}`インストール / 環境の問題 <ts-install-env>` |

## レシピ 1: MEP 前に抽出で止まる

- 兆候:
 - 元素情報不足、原子数不一致、活性部位モデル（バインディングポケット）の抽出結果が空に近い。
- 最初の確認:
 - 入力構造が同じ前処理フローで作られ、原子順が揃っているか。
 - `extract` / `all` 前に元素列が埋まっているか。
- 典型的な修正手順:
 - `pdb2reaction add-elem-info -i input.pdb -o input_fixed.pdb` で元素列を修復 — 抽出再実行 — 活性部位モデルサイズ（`--radius`）/残基選択（`--selected-resn`）を再確認。残基名・ID・chain付きselectorの指定形式は CLI 規約の {ref}`ja-selected-resn-takes-ids` を参照。

## レシピ 2: 電荷/スピンの解決で止まる

- 兆候:
 - 特に非 `.gjf` 入力で総電荷未解決エラーが出る。
- 最初の確認:
 - 対象状態に対して総電荷・多重度が妥当か。
 - `--ligand-charge/-l` の各残基キーが妥当か検証する。
 - 結果が物理的に不自然な場合は [CLI 規約](cli-conventions.md) の電荷解決順序を再確認。
- 典型的な修正手順:
 - 重要な実行では `-q/--charge` または `-l/--ligand-charge` と `-m` を明示し、scan/path/tsopt を再試行。

## レシピ 3: 環境依存で止まる

- 兆候:
 - DMF import 失敗、CUDA 不整合、図出力バックエンド不在。
- 最初の確認:
 - 実行環境にオプション依存パッケージ（cyipopt 等）が入っているか。
 - GPU 可視性と PyTorch CUDA 互換性に問題がないか。
- 典型的な修正手順:
 - 先に環境を修復し、`pdb2reaction --version` や `python -c "import torch; print(torch.cuda.is_available())"` で確認後、`--dry-run` で事前チェックしてから本実行。

## レシピ 4: 収束・後処理で止まる

- 兆候:
 - TSOPT が停滞、IRC が不安定、MEP 精密化が途中停止。
- 最初の確認:
 - Cartesian PHVA の負の振動数が1本だけで、対応するモードが反応座標方向の変位を示すか。5 cm⁻¹の閾値はモードファイル、平坦化、任意の鞍点回復に使いますが、最終的な本数は変えません。
 - 虚振動数の本数が誤っている場合（偽の 2 本目の小さいモード、または支配的な反応モードが無い）は、`--precision fp64` で精度を上げる、および／または `--coord-type dlc` に切り替え、残った小さいモードは `--flatten` を追加します。これらの手法は補完的です。{ref}`tsopt: 最適化後に虚振動数の本数が誤っている場合 <ja-wrong-imaginary-mode-count>` を参照。
 - ステップサイズ / 信頼半径（YAML キー `max_step`, `trust_radius`/`trust_min`/`trust_max`）と、最適化モード / フラット化（CLI フラグ `--opt-mode`, `--flatten`）は併用してください。YAML セクション構成は [YAML リファレンス](yaml-reference.md)、正規の修正手順は {ref}`計算 / 収束の問題 <ts-calc-conv>` を参照。
 - エネルギー平坦化により `stalled` で停止した場合も未収束です。{ref}`計算 / 収束の問題 <ts-calc-conv>` に従って確認してから再実行してください。
- 典型的な修正手順:
 - 小規模ケースで条件を詰め、安定化後に本番条件へ戻す。
