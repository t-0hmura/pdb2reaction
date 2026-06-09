# クイックスタート: `pdb2reaction all --scan-lists`（Scan モード）

## 目的

`pdb2reaction all` のフルワークフローを、`-s/--scan-lists` で 1 つ以上の原子間距離を駆動して単一構造から実行します。これにより、段階的スキャン → 最小エネルギー経路（MEP）精密化 → （任意）遷移状態（TS）最適化 + 固有反応座標（IRC）が自動で連鎖します。

## 事前に必要なもの

- 入力構造: `.pdb`
- 対象状態に対応した電荷（`-q/--charge` または `-l/--ligand-charge`）・多重度（`-m`）

## 方法 A: `-s/--scan-lists` インラインリテラル（既定）

`-s/--scan-lists` はコマンドライン上で Python リテラル文字列を直接受け取ります。原子セレクタの構文（残基/原子トークン、区切り文字、順序）と外側/内側のクォートのルールについては、{ref}`CLI 規約: スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

### 基本構文

各リテラルは 3 要素タプル `(atom1, atom2, target_distance_Å)` のリストです。3 番目の要素は必ず **ångström** 単位の目標距離で、ちょうど 3 要素が必要です。1 リテラル = 1 ステージ。

```bash
# 単一ステージ、整数原子インデックス（デフォルトで1-based）
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[(1, 5, 1.35)]' -o ./result_scan

# 単一ステージ、PDBセレクタ文字列
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 \
 -s '[("TYR,285,CA", "SAM,309,C10", 1.35)]' -o ./result_scan
```

タンパク質–リガンド複合体 PDB を入力とする場合は `-c/--center` の指定が必須です。低分子 `.pdb` / `.xyz` では `-c` を省略し、`-l` の代わりに `-q` を直接渡してください。`-m/--multiplicity` のデフォルトは `1`（一重項）ですが、ここでは例として明示的に指定しています。

### 複数ステージ

複数のリテラルを渡すと、各リテラルが順番に 1 ステージとして実行されます:

```bash
# ステージ1: 1つの結合を 1.35 Å に駆動
# ステージ2: 2つの結合を同時に駆動
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 -s \
  '[("TYR,285,CA","SAM,309,C10",1.35)]' \
  '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
  -o ./result_scan
```

ステージは順番に実行され、各々は前ステージの緩和結果から開始します。

## 期待される出力

```text
result_scan/
├── summary.log
├── summary.json
├── mep.pdb                        # MEP 全体パス（ルートへ昇格）
├── energy_diagram_MEP.png         # MEP エネルギー図（ルートへ昇格）
├── segments/
│   └── seg_01/                    # 反応セグメントごとの成果物
└── _work/                         # パイプラインの作業領域（削除可）
    ├── scan/
    │   ├── preopt/                # 事前最適化構造
    │   ├── stage_01/              # スキャンステージ 1
    │   │   ├── scan_trj.xyz       # スキャン軌跡
    │   │   └── scan.pdb
    │   └── stage_02/              # スキャンステージ 2（マルチステージ時）
    └── path_search/               # MEP 探索（--refine-path False 時は path_opt/）
```

**確認ポイント:**

1. `_work/scan/stage_01/scan_trj.xyz` — 結合距離の変化を PyMOL で確認
2. `mep.pdb` — 最適化後の MEP 軌跡（出力ルートへ昇格）
3. `summary.log` — 障壁高さと結合変化

**ヒント:** `--print-parsed` を付けて（Ctrl-C で中断して）スキャン設定を事前確認:

```bash
pdb2reaction scan -i input.pdb -q 0 -s '[(1, 5, 1.35)]' --print-parsed
```

## 補足

- `-s/--scan-lists` は `all` ではインライン Python リテラルのみを受け取ります。単独の `scan` サブコマンドはこれに加えて YAML/JSON スペックファイルパスも受け取れます（[scan](scan.md) を参照）。
- スキャン step のデフォルトは両コマンドで同じで、フラグ名のみが異なります（`all` では prefix 付き、単独の `scan` では prefix なし）:

  | コマンド | 距離ごとの step 幅 | 調和バイアス強度 |
  | --- | --- | --- |
  | `pdb2reaction all` | `--scan-max-step-size 0.20 Å` | `--scan-bias-k 300 eV/Å²` |
  | 単独の `pdb2reaction scan` | `--max-step-size 0.20 Å` | `--bias-k 300 eV/Å²` |

  どちらの形式でも、または YAML の `bias` ブロックで上書き可能です（[scan](scan.md) / [yaml-reference](yaml-reference.md#bias) 参照）。
- 各スキャン stage は、最終リラックス構造に対する bond-change check（`has_bond_change`）で終了します。stage ごとの結果は scan ログに記録され、`--out-json` 指定時には scan 出力ディレクトリに書き出される集約 `result.json`（その `stages` 配列内）にも記録されます。
- 再帰的 MEP refinement（`path-search`）は scan の端点を**無条件**に消費します。実行されるかどうかのゲートは、scan stage の bond-change フラグ（`has_bond_change`）ではなく `--refine-path` です。
- `pdb2reaction all --help-advanced` で全オプション（スキャン制御を含む）を確認できます。
- 単独の `scan` サブコマンド（MEP 精密化なし）については [scan](scan.md) を参照してください。

## 次の導線

- 全オプション: [all](all.md)
- TS 候補の検証: [クイックスタート: TS のみモード](quickstart-tsopt-freq.md)
