# `irc` サブコマンド

## 概要
UMAを使用してEulerPCベースの固有反応座標（IRC）積分を実行します。CLIは意図的に狭く設計されています: 以下に記載されていないものはすべてYAMLで提供する必要があります。

## 使用法
```bash
pdb2reaction irc -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [--ligand-charge <number|'RES:Q,...'>] \
                 [--workers N] [--workers-per-node N] [-m 2S+1]
                 [--max-cycles N] [--step-size Δs] [--root k]
                 [--forward True|False] [--backward True|False]
                 [--freeze-links True|False]
                 [--out-dir DIR]
                 [--convert-files {True\|False}] [--ref-pdb FILE]
                 [--hessian-calc-mode Analytical|FiniteDifference]
                 [--args-yaml FILE]
```

### 例
```bash
# 順方向のみ、有限差分ヘシアン、大きいステップサイズ
pdb2reaction irc -i ts.xyz -q -1 -m 2 --forward True --backward False \
                --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB入力: 完成軌跡と方向別軌跡もPDBとしてエクスポート
pdb2reaction irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## ワークフロー
1. **入力準備** – `geom_loader` がサポートする任意のフォーマットを受け入れ。ソースが `.pdb` の場合、EulerPC軌跡は元のトポロジーを使用してPDBに自動変換
2. **設定マージ** – デフォルト → CLI → YAML（`geom`、`calc`、`irc`）
3. **IRC積分** – EulerPCが `irc.forward/backward`、`irc.step_length`、`irc.root`、およびUMAを通じて設定されたヘシアンワークフローに従って順方向/逆方向分岐を積分
4. **出力** – 軌跡（`finished`、`forward`、`backward`）は `.trj` として書き込まれ、PDB入力の場合は `.pdb` にミラーリング

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる遷移状態構造 | 必須 |
| `-q, --charge INT` | 総電荷; `calc.charge` をオーバーライド | テンプレート/導出が適用されない限り必須 |
| `--ligand-charge TEXT` | `-q` が省略された場合に使用される総電荷または残基名ごとのマッピング | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1） | `.gjf` テンプレート値または `1` |
| `--max-cycles INT` | 最大IRCステップ | `125` |
| `--step-size FLOAT` | 質量重み付き座標でのステップ長 | `0.10` |
| `--root INT` | 初期変位の虚数モードインデックス | `0` |
| `--forward {True\|False}` | 順方向分岐を実行 | `True` |
| `--backward {True\|False}` | 逆方向分岐を実行 | `True` |
| `--freeze-links {True\|False}` | PDB入力用、リンクH親を凍結 | `True` |
| `--out-dir TEXT` | 出力ディレクトリ | `./result_irc/` |
| `--convert-files {True\|False}` | PDB入力用のXYZ/TRJ → PDBコンパニオンをトグル | `True` |
| `--ref-pdb FILE` | 入力がXYZ/GJFの場合に使用する参照PDBトポロジー | _None_ |
| `--hessian-calc-mode CHOICE` | UMAヘシアンモード | `FiniteDifference` |
| `--args-yaml FILE` | YAMLオーバーライド | _None_ |

## 出力
```
out_dir/ (デフォルト: ./result_irc/)
├─ irc_data.h5                # dump_every > 0の場合に書き込み
├─ <prefix>finished_irc.trj   # 完全なIRC軌跡
├─ <prefix>forward_irc.trj    # 順方向分岐が実行された場合
├─ <prefix>backward_irc.trj   # 逆方向分岐が実行された場合
└─ *.pdb                      # PDB入力用の軌跡コンパニオン（変換有効時）
```

## 注意事項
- VRAMが十分な場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨
- `--freeze-links` はPDB入力にのみ適用
