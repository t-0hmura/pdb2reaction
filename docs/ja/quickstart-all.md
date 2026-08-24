# クイックスタート: `pdb2reaction all`（Endpoint モード）

## 目的

2 つの完全系 PDB（反応物 R と生成物 P）から、end-to-end のワークフローを 1 回実行します。

## 前提条件

- pdb2reaction がインストール済みであること（[インストール](installation.md) を参照）
- **水素原子が追加済み**の 2 つの PDB ファイル（反応物 R と生成物 P）
- すべての入力 PDB で同じ原子が同じ順序で含まれていること

> **ファイル名について:** 例の `1.R.pdb` と `3.P.pdb` は geranyl pyrophosphate (GPP) C6-メチル基転移酵素 BezA のサンプルディレクトリ（[`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples)）に同梱された反応物/生成物 PDB に対応します（`1.R.pdb` = 反応物状態、`3.P.pdb` = 生成物状態、追加の反応物/生成物/中間体構造を含む実行向けに `2.*.pdb` の中間状態も利用可能）。ご自身の反応では、2 つ以上の全系 PDB に置き換えてください。下記コマンドをそのまま試すには、まず同梱例を取得してください: `git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction/examples`。

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

> **VRAM 注意:** `--dft` は抽出clusterに対してGPU4PySCF の一点計算を
> 実行します。必要memoryは構造・基底・汎関数・精度・software stackに依存するため、
> 対象nodeで代表構造をpilot実行してpeak memoryを測定してください。OOM時は、より
> 小さい基底／縮小clusterで `pdb2reaction dft` を単独実行するか、より大きいnodeを
> 使用します。`[dft]` extraのinstallも必要です
> （[インストール](installation.md) 手順7）。

## 期待される出力

成功時のディレクトリ構造:

```text
result_all/
├── summary.log                    # テキストサマリ
├── summary.json                   # 機械可読な結果
├── mep.pdb                        # マージ済み MEP 経路（ルート直下に配置）
├── energy_diagram_MEP.png         # 全セグメントの MEP エネルギープロファイル
└── _work/                         # パイプライン作業領域（削除可）
    └── path_opt/                  # MEP エンジン生出力
        ├── hei_seg_01.{xyz,pdb}   # MEP の最高エネルギー像
        └── summary.json           # path-opt エンジンの結果
```

最小コマンドは MEP ステージ終了時に停止するため、`segments/` は作成しません。
`--tsopt` を付け、反応セグメントの検証に成功すると
`segments/seg_01/{reactant.pdb,ts.pdb,product.pdb}`、`ts/`、`irc/` が追加され、
`--thermo` も付けると `freq/` が追加されます。

**確認ポイント:**

1. `summary.json` — 利用可否は `scientific_status` と `scientific_status_reasons` で判定します。path mode の `segments[].barrier_kcal` は生の MEP 電子障壁であり、要求した後処理の結果は `rate_limiting_step` と `post_segments` を確認します。`status` は互換性のために残されています
2. `_work/path_opt/hei_seg_01.pdb` — 最高エネルギー像を確認。`--tsopt` 時は正規 `segments/seg_01/*.pdb` の R/TS/P も確認
3. `energy_diagram_*.png` — 明確な障壁があるエネルギープロファイル

**成功時のターミナル出力例:**

```
[time] Elapsed Time for Whole Pipeline: HH:MM:SS.sss
```

（実行時間は系サイズ、GPU、有効化したステージで変動します。）

`--tsopt` が有効な場合:

```
[Imaginary modes] n=1 ([-425.9])
```

一次鞍点では反応座標方向に虚振動が **ちょうど 1 つ** 現れます。IRC 検証（`--tsopt` の一部として自動的に実行される）により、その鞍点が想定どおりの反応物・生成物を結ぶことを確認します。

## 補足

- `pdb2reaction all --help` は主要オプション、`pdb2reaction all --help-advanced` は全オプションを表示します。

## 次のステップ

- 単一構造の段階的スキャン: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
- TS 候補の検証: [クイックスタート: TS のみモード](quickstart-tsopt-freq.md)
- 全オプション: [all](all.md)
