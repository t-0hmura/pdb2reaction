# `energy-diagram`

`pdb2reaction energy-diagram` は与えた数値エネルギーだけを入力として状態エネルギーダイアグラムを描画します。構造ファイルの読み込みや `--thermo` / `--dft` のような量子 / 熱力学 / MLIP 計算は一切行いません。状態エネルギーの数値（例: `summary.json` から取り出した値）が既にあり、整形済みダイアグラムだけが欲しいときに使います。出力は画像ファイル 1 つと、任意で機械可読なサイドカーです。

## 実行例

```bash
# コマンド形式
pdb2reaction energy-diagram -i VALUES... [-o OUTPUT] [--label-x...] [--label-y...]
```

```bash
# リスト文字列で指定（その場限りの作図に推奨）
pdb2reaction energy-diagram -i "[0, 12.5, 4.3]" -o energy.png
```

```bash
# フラグ繰り返しで指定（値ごとに -i を 1 つ）
pdb2reaction energy-diagram -i 0 -i 12.5 -i 4.3 -o energy.png
```

```bash
# X/Yラベルを指定（--label-x / --label-y 共通）
pdb2reaction energy-diagram -i "[0, 12.5, 4.3]" --label-x "['R','TS','P']" --label-y "ΔE (kcal/mol)" -o energy.png
```

## 処理の流れ
1. `-i/--input` から値を収集します（繰り返し指定、またはリスト文字列に対応）。
2. 全値を float として解釈し、2 点未満なら早期にエラーを返します。
3. 任意の `--label-x` を解釈します。未指定時は `S1`, `S2`,... を自動生成します。
4. `--label-x` の個数と値の個数の一致を検証し、図を描画します。
5. `-o/--output` に画像を保存し、保存先パスを表示します。

## 出力
```
OUTPUT.(png|jpg|jpeg|svg|pdf)
result.json   # 任意のサイドカー（--out-json 指定時）。status / n_points / files。per-point の energies・labels は含まれない
```
- `-o/--output` を省略した場合、カレントディレクトリに `energy_diagram.png` を出力します。
- 出力拡張子がない場合は `.png` が自動で補完されます。
- 必要なら親ディレクトリを自動作成します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input TEXT...` | 数値入力（複数引数またはリスト形式文字列） | 必須 |
| `-o, --output PATH` | 出力画像パス（`.png/.jpg/.jpeg/.svg/.pdf`） | `energy_diagram.png` |
| `--label-x TEXT...` | X 軸状態ラベル（入力値と同じ個数が必要） | `S1, S2,...` |
| `--label-y TEXT` | Y 軸ラベル | `ΔE (kcal/mol)` |
| `--out-json/--no-out-json` | 出力画像の隣に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

すべてのフラグ一覧は生成された [コマンドリファレンス](../reference/commands/index.md) を参照してください。

## 注記
- 入力順がそのまま描画順になります。
- 入力値は最低 2 点必要です。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [trj2fig](trj2fig.md) -- 軌跡エネルギーからプロファイルを描画
- [all](all.md) -- エネルギーダイアグラム出力を含む end-to-end 実行
