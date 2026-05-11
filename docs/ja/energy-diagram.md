# `energy-diagram`

## 概要

> **要約:** 数値エネルギーだけを入力として状態エネルギーダイアグラムを描画します（構造ファイル不要、量子/熱力学計算なし）。

### 要点
- **想定場面:** 状態エネルギーの数値（例: `summary.json` から取り出した値）が既にあり、整形済みダイアグラムだけが欲しいとき。
- **手法:** Matplotlib による描画のみ。構造ファイルの読み込み、QM / 熱力学 / MLIP 計算は一切行いません。
- **主な出力:** 画像ファイル 1 つ（`.png` / `.jpg` / `.jpeg` / `.svg` / `.pdf`）。`--out-json` 指定時は `result.json` も出力。`-o` 省略時の既定出力は `energy_diagram.png`。
- **デフォルト値:** `-o energy_diagram.png`、`--label-x` は `S1`, `S2`, ... の自動採番、`--label-y "ΔE (kcal/mol)"`、`--out-json False`。
- **次のステップ:** エネルギー軌跡も合わせて描画したい場合は [`trj2fig`](trj2fig.md) を使用。`all` の出力（`summary.json`）と組み合わせれば [`all`](all.md) パイプラインの最終可視化として機能します。

`pdb2reaction energy-diagram` は与えた数値を可視化するのみで、構造ファイルの読み込みや `--thermo` / `--dft` のような計算は行いません。

## 使用法
```bash
pdb2reaction energy-diagram -i VALUES... [-o OUTPUT] [--label-x...] [--label-y...]
```

## 例
```bash
# 数値を複数引数で指定
pdb2reaction energy-diagram -i 0 12.5 4.3 -o energy.png

# リスト文字列で指定
pdb2reaction energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png

# X/Yラベルを指定
pdb2reaction energy-diagram -i 0 12.5 4.3 --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## ワークフロー
1. `-i/--input` から値を収集します（繰り返し指定、1 フラグ後の複数値、リスト文字列に対応）。
2. 全値を float として解釈し、2 点未満なら早期にエラーを返します。
3. 任意の `--label-x` を解釈します。未指定時は `S1`, `S2`,... を自動生成します。
4. `--label-x` の個数と値の個数の一致を検証し、図を描画します。
5. `-o/--output` に画像を保存し、保存先パスを表示します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input TEXT...` | 数値入力（複数引数またはリスト形式文字列） | 必須 |
| `-o, --output PATH` | 出力画像パス（`.png/.jpg/.jpeg/.svg/.pdf`） | `energy_diagram.png` |
| `--label-x TEXT...` | X 軸状態ラベル（入力値と同じ個数が必要） | `S1, S2,...` |
| `--label-y TEXT` | Y 軸ラベル | `ΔE (kcal/mol)` |
| `--out-json/--no-out-json` | 出力画像の隣に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

## 出力
```
OUTPUT.(png|jpg|jpeg|svg|pdf)
```
- `-o/--output` を省略した場合、カレントディレクトリに `energy_diagram.png` を出力します。
- 出力拡張子がない場合は `.png` が自動で補完されます。
- 必要なら親ディレクトリを自動作成します。

## 注意事項
- 入力順がそのまま描画順になります。
- 入力値は最低 2 点必要です。
- 構造ファイルの読み込みやエネルギー計算は行いません。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [trj2fig](trj2fig.md) -- 軌跡エネルギーからプロファイルを描画
- [all](all.md) -- エネルギーダイアグラム出力を含む end-to-end 実行
