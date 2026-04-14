# クイックスタート: `pdb2reaction all`

## 目的

2 つの完全系 PDB（反応物 R と生成物 P）から、end-to-end のワークフローを 1 回実行します。

## 前提条件

- pdb2reaction がインストール済みであること（[インストール](installation.md) を参照）
- **水素原子が追加済み**の 2 つの PDB ファイル（反応物 R と生成物 P）
- すべての入力 PDB で同じ原子が同じ順序で含まれていること

> **ファイル名について:** 例の `1.R.pdb` と `3.P.pdb` は GPP C6-メチル基転移酵素 BezA のサンプルディレクトリ（[`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples)）に同梱された反応物/生成物 PDB に対応します（`1.R.pdb` = 反応物状態、`3.P.pdb` = 生成物状態、多段階向けに `2.*.pdb` の中間状態も利用可能）。ご自身の反応では、2 つ以上の全系 PDB に置き換えてください。

## 最小コマンド

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

後処理（TS 最適化、振動解析・熱化学、DFT 一点計算）まで同時に実行する場合:

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

> **VRAM 注意:** `--dft` は抽出されたクラスターモデルに対して GPU4PySCF の一点計算を走らせます。~200 原子を超える系では VRAM が 24 GB 未満の GPU で容易に OOM に至ります。`CUDA out of memory` になった場合は、`--dft` を外して基底を縮小または原子数を削減した上で別途 `pdb2reaction dft` を実行するか、より大容量の GPU ノードへ移してください。事前に `[dft]` extra のインストールも必要です（[インストール](installation.md) 手順 7 を参照）。

## 期待される出力

成功時のディレクトリ構造:

```text
result_all/
├── summary.log                    # テキストサマリ
├── summary.json                   # 機械可読な結果
├── path_search/
│   ├── mep.pdb                    # MEP 軌跡
│   ├── energy_diagram_UMA_all.png # エネルギープロファイル
│   ├── summary.json               # path-search の結果
│   └── post_seg_01/               # 後処理（--tsopt 時）
│       ├── ts/final_geometry.pdb
│       ├── irc/finished_irc_trj.xyz
│       └── freq/
└── seg_01/                        # IRC 最適化 R/TS/P 構造
    ├── reactant.pdb
    ├── ts.pdb
    └── product.pdb
```

**確認ポイント:**

1. `summary.log` — `status: success` と障壁高さ (kcal/mol)
2. `seg_01/*.pdb` — PyMOL で R/TS/P 構造を確認
3. `energy_diagram_*.png` — 明確な障壁があるエネルギープロファイル

**成功時のターミナル出力例:**

```
[all] Elapsed for Whole Pipeline: 00:05:06.123
```

`--tsopt` が有効な場合:

```
[Imaginary modes] n=1  ([-425.9])
```

|振動数| >= 100 cm⁻¹ の虚振動 1 つが有効な TS を示します。

## 補足

- `pdb2reaction all --help` は主要オプション、`pdb2reaction all --help-advanced` は全オプションを表示します。

## 次の導線

- 単一構造の段階的スキャン: [クイックスタート: `pdb2reaction scan`](quickstart-scan.md)
- TS 最適化と検証: [クイックスタート: `pdb2reaction tsopt`](quickstart-tsopt-freq.md)
- 全オプション: [all](all.md)
