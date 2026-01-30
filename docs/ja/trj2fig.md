# `trj2fig` サブコマンド

## 概要
`trj2fig` はXYZ軌跡から整形済みのエネルギープロファイルを生成します。デフォルトでは、各フレームのコメント行に含まれる Hartree エネルギーを読み取り、kcal/mol または Hartree に変換し、必要に応じて基準フレームに対する相対値にし、静的/インタラクティブ図とCSVを出力します。基準は最初のフレーム（`init`）、`--reverse-x` 使用時の最後のフレーム、または任意の明示インデックスを指定できます。`-q/--charge` と `-m/--multiplicity` を与えると、コメント行ではなく `uma_pysis` によって各フレームのエネルギーを再計算します。

## 使用法
```bash
pdb2reaction trj2fig -i TRAJECTORY.xyz [-o OUTPUTS...] [-q CHARGE] [-m MULT] [options]
```

### 例
```bash
# デフォルトPNG、最初のフレームを基準としたΔE
pdb2reaction trj2fig -i traj.xyz

# CSV + SVG、フレーム5を基準、Hartree表記
pdb2reaction trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

# 複数形式 + x軸を反転（基準は最後のフレーム）
pdb2reaction trj2fig -i traj.xyz --reverse-x True -o energy.png energy.html energy.pdf

# UMAで全フレームのエネルギーを再計算
pdb2reaction trj2fig -i traj.xyz -q 0 -m 1 -o energy.png
```

## ワークフロー
1. XYZ軌跡を解析し、各フレームのコメント行から最初の浮動小数点数を読み取ります（科学表記は未対応）。`-q/-m` がある場合は `uma_pysis` で再計算した Hartree エネルギーを使用します。エネルギーが取得できない場合は実行を中断します。
2. 参照指定を正規化:
   - `init` → フレーム `0`（`--reverse-x` が有効な場合は最後のフレーム）
   - `None`/`none`/`null` → 絶対エネルギー（参照なし）
   - 整数リテラル → 0始まりのフレーム番号
3. エネルギーを kcal/mol または Hartree に変換し、参照指定がある場合は参照値を差し引いて ΔE を作成します。
4. Plotly 図を作成（強調ティック、スプライン補間、マーカー、タイトルなし）し、指定された拡張子へ書き出します。
5. 必要に応じて `frame`, `energy_hartree` と、指定単位での ΔE/E 列を含むCSVを出力します。

## CLIオプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | コメント行にエネルギーを含むXYZ軌跡 | 必須 |
| `-o, --out PATH` | 反復可能な出力ファイル名（`.png/.jpg/.jpeg/.html/.svg/.pdf/.csv`） | `energy.png` |
| _extra arguments_ | オプション後に続く位置引数（`-o` リストにマージ） | _None_ |
| `--unit {kcal,hartree}` | 出力単位 | `kcal` |
| `-r, --reference TEXT` | 参照指定（`init`、`None`、または0始まり整数） | `init` |
| `-q, --charge INT` | 総電荷。指定時は `uma_pysis` で再計算 | _None_ |
| `-m, --multiplicity INT` | スピン多重度（2S+1）。指定時は `uma_pysis` で再計算 | _None_ |
| `--reverse-x {True\|False}` | x軸を反転し、`init` の参照を最後のフレームに変更 | `False` |

## 出力
```
<output>.[png|jpg|jpeg|html|svg|pdf]  # 指定された拡張子すべてを出力（デフォルト: energy.png）
<output>.csv                          # CSV指定時のみ
```
- `-o` や位置引数が無い場合は `energy.png` をカレントディレクトリに書き出します。
- CSV には `frame`, `energy_hartree` と、参照がある場合は `delta_kcal`/`delta_hartree`、参照がない場合は `energy_kcal`/`energy_hartree` が含まれます。
- パース失敗や未対応拡張子がある場合、コンソールに診断メッセージが出力されます。

## 注意事項
- エネルギーはコメント行内の最初の10進数から取得されます。不正なコメント行はエラーになります。
- 未対応拡張子がある場合は実行が中断されます。`.png` は Plotly の `scale=2` で高解像度出力されます。
- `--reverse-x` は軸の向きと `-r init` の解釈の両方に影響します。
- 旧 `--output-peak` オプションは削除されています。
