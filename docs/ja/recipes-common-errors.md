# 典型エラー別レシピ

症状は分かるが、どのコマンドページから見ればよいか迷うときの入口です。
詳細は [トラブルシューティング](troubleshooting.md) を並行して参照してください。

## 早見表

| 症状 | 最初にやること | 詳細 |
| --- | --- | --- |
| 元素カラム欠落で抽出が止まる | 元の PDB に `add-elem-info` を適用してください | [トラブルシューティング](troubleshooting.md) |
| 「電荷が必須」系エラー | `-q/--charge` または `-l/--ligand-charge` を明示指定してください | [トラブルシューティング](troubleshooting.md) |
| 計算は通るが状態/エネルギーが不自然 | [CLI 規約](cli-conventions.md) の電荷解決順序を再確認してください | [トラブルシューティング](troubleshooting.md) |
| DMF モードの import エラー（`cyipopt`） | `conda install -c conda-forge cyipopt` を実行してください | [トラブルシューティング](troubleshooting.md) |
| TSOPT/IRC が収束しない | LBFGS/Dimer: `max_step` を調整。RFO/RS-I-RFO: `trust_radius`/`trust_min`/`trust_max` を調整。サイクル上限を増やし、TS 品質を確認 | [トラブルシューティング](troubleshooting.md) |
| CUDA/GPU 実行時エラー | `torch.cuda.is_available()` と CUDA バージョンの整合を確認してください | [トラブルシューティング](troubleshooting.md) |
| 図の出力失敗 | `plotly_get_chrome -y` で Chrome ランタイムを導入してください | [トラブルシューティング](troubleshooting.md) |

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
 - LBFGS/Dimer の場合 `max_step`、RFO/RS-I-RFO の場合 `trust_radius`/`trust_min`/`trust_max` が厳し過ぎないか。サイクル上限の確認も。
- 典型的な修正手順:
 - 小規模ケースで条件を詰め、安定化後に本番条件へ戻す。
