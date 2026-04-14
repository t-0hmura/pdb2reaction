# 典型エラー別レシピ

症状は分かるが、どのコマンドページから見ればよいか迷うときの入口です。
詳細は [トラブルシューティング](troubleshooting.md) を並行して参照してください。

## 早見表

各行の「詳細」列は [トラブルシューティング](troubleshooting.md) の該当セクションへ直接リンクします。

| 症状 | 最初にやること | 詳細（セクション） |
| --- | --- | --- |
| 元素カラム欠落で抽出が止まる | 元の PDB に `add-elem-info` を適用してください | {ref}`入力 / 抽出の問題 <ts-input-extraction>` |
| 「電荷が必須」系エラー | `-q/--charge` または `-l/--ligand-charge` を明示指定してください | {ref}`電荷 / スピンの問題 <ts-charge-spin>` |
| 計算は通るが状態/エネルギーが不自然 | [CLI 規約](cli-conventions.md) の電荷解決順序を再確認してください | {ref}`入力 / 抽出の問題 <ts-input-extraction>` |
| DMF モードの import エラー（`cyipopt`） | `conda install -c conda-forge cyipopt` を実行してください | {ref}`インストール / 環境の問題 <ts-install-env>` |
| TSOPT が収束しない | LBFGS/Dimer: `max_step` を**縮小**。RFO/RS-I-RFO: `trust_radius`/`trust_min`/`trust_max` を**縮小**。サイクル上限を増やし、TS 品質を確認 | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| IRC が正常に終了しない | `--step-size` を縮小、`--max-cycles` を増加、虚振動数が 1 本のみか確認 | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| opt/TSOPT が `max_cycles` で停止し、`max(force)` が閾値をわずかに超える | 通常は `opt.energy_plateau` フォールバック（v0.3.5 新機能）が自動で処理します。手動回避は `--thresh gau` または `--thresh gau_loose` | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| MEP 探索（GSM/DMF）が失敗 | `--max-nodes` をデフォルト 20 から増やす、`--preopt` 有効化、別の `--mep-mode` を試す | {ref}`計算 / 収束の問題 <ts-calc-conv>` |
| CUDA/GPU 実行時エラー | `torch.cuda.is_available()` と CUDA バージョンの整合を確認してください | {ref}`インストール / 環境の問題 <ts-install-env>` |
| 図の出力失敗 | `plotly_get_chrome -y` で Chrome ランタイムを導入してください | {ref}`インストール / 環境の問題 <ts-install-env>` |

## レシピ 1: MEP 前に抽出で止まる

- 兆候:
 - 元素情報不足、原子数不一致、活性部位モデル（バインディングポケット）抽出結果が空に近い。
- 最初の確認:
 - 入力構造が同じ前処理フローで作られ、原子順が揃っているか。
 - `extract` / `all` 前に元素カラムが埋まっているか。
- 典型的な修正手順:
 - `pdb2reaction add-elem-info -i input.pdb -o input_fixed.pdb` で元素列を修復 -> 抽出再実行 -> 活性部位モデルサイズ（`--radius`）/残基選択（`--selected-resn`）を再確認。

## レシピ 2: 電荷/スピンの解決で止まる

- 兆候:
 - 特に非 `.gjf` 入力で総電荷未解決エラーが出る。
- 最初の確認:
 - 対象状態に対して総電荷・多重度が妥当か。
 - `--ligand-charge` の残基キーが入力構造と一致しているか。
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
 - TS 候補が虚振動数 1 本（|ν| >= 100 cm⁻¹）のみを持ち、対応する虚振動モードが反応座標方向の変位を示すか。
 - LBFGS/Dimer の場合は `max_step` を**縮小**してください（YAML セクション: `lbfgs` / `hessian_dimer`）。RFO/RS-I-RFO の場合は `trust_radius`/`trust_min`/`trust_max` を**縮小**してください（YAML セクション: `rfo` / `rsirfo`）。サイクル上限の確認も。セクション構成は [YAML リファレンス](yaml-reference.md) を参照してください。
 - `max_cycles` 到達時に力のノルムが閾値をわずかに超えているだけでエネルギーがフラット化している場合、`opt.energy_plateau` フォールバック（v0.3.5 新機能）が自動で収束と判定するはずです。効かない場合は `--thresh gau` または `--thresh gau_loose` で力の閾値を緩めてください。
- 典型的な修正手順:
 - 小規模ケースで条件を詰め、安定化後に本番条件へ戻す。
