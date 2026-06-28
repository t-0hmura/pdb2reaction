# クイックスタート: `pdb2reaction all`（Endpoint モード）

## 目的

2 つの完全系 PDB（反応物 R と生成物 P）から、end-to-end のワークフローを 1 回実行します。

## 前提条件

- pdb2reaction がインストール済みであること（[インストール](installation.md) を参照）
- **水素原子が追加済み**の 2 つの PDB ファイル（反応物 R と生成物 P）
- すべての入力 PDB で同じ原子が同じ順序で含まれていること

> **ファイル名について:** 例の `1.R.pdb` と `3.P.pdb` は geranyl pyrophosphate (GPP) C6-メチル基転移酵素 BezA のサンプルディレクトリ（[`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples)）に同梱された反応物/生成物 PDB に対応します（`1.R.pdb` = 反応物状態、`3.P.pdb` = 生成物状態、追加の反応物/生成物/中間体構造を含む run 向けに `2.*.pdb` の中間状態も利用可能）。ご自身の反応では、2 つ以上の全系 PDB に置き換えてください。下記コマンドをそのまま試すには、まず同梱例を取得してください: `git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction/examples`。

## 最小コマンド

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

### （オプション）同一実行で後処理まで行う

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

> **VRAM 注意:** `--dft` は抽出されたクラスターモデルに対して GPU4PySCF の一点計算を走らせます。~200 原子を超える系では VRAM が 24 GB 未満の GPU で容易に OOM を起こします。`CUDA out of memory` になった場合は、`--dft` を外して基底を縮小または原子数を削減した上で別途 `pdb2reaction dft` を実行するか、より大容量の GPU ノードへ移してください。事前に `[dft]` extra のインストールも必要です（[インストール](installation.md) 手順 7 を参照）。

## 期待される出力

成功時のディレクトリ構造:

```text
result_all/
├── summary.log                    # テキストサマリ
├── summary.json                   # 機械可読な結果
├── mep.pdb                        # マージ済み MEP 経路（ルート直下に昇格）
├── energy_diagram_MEP.png         # 全セグメントの MEP エネルギープロファイル
├── segments/
│   └── seg_01/                    # 反応セグメント別の成果物
│       ├── reactant.pdb           # 正規 IRC 最適化 R/TS/P
│       ├── ts.pdb
│       ├── product.pdb
│       ├── ts/final_geometry.pdb  # --tsopt 時
│       ├── irc/finished_irc_trj.xyz
│       └── freq/                  # --thermo 時
└── _work/                         # パイプライン作業領域（削除可）
    └── path_opt/                  # MEP エンジン生出力
        └── summary.json           # path-opt エンジンの結果
```

**確認ポイント:**

1. `summary.json` — `status` フィールド（`"success"` / `"partial"` / `"failed"`）とセグメントごとの `barrier_kcal` を確認。`summary.log` は同じ情報を人間可読形式でミラーします
2. `segments/seg_01/*.pdb` — PyMOL で R/TS/P 構造が化学的に妥当か確認
3. `energy_diagram_*.png` — 明確な障壁があるエネルギープロファイル

**成功時のターミナル出力例:**

```
[all] Elapsed for Whole Pipeline: HH:MM:SS.sss
```

（実行時間は系サイズ、GPU、有効化したステージで変動します。）

`--tsopt` が有効な場合:

```
[Imaginary modes] n=1  ([-425.9])
```

一次鞍点では反応座標方向に虚振動が **ちょうど 1 つ** 現れます。IRC 検証（`--tsopt` の一部として自動的に実行される）により、その鞍点が想定どおりの反応物・生成物を結ぶことを確認します。

## 補足

- `pdb2reaction all --help` は主要オプション、`pdb2reaction all --help-advanced` は全オプションを表示します。

## 次のステップ

- 単一構造の段階的スキャン: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- TS 候補の検証: [クイックスタート: TS のみモード](quickstart-tsopt-freq.md)
- 全オプション: [all](all.md)
