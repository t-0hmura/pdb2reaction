# クイックスタート: 単一構造スキャンワークフロー（`--scan-lists`）

## 目的

`pdb2reaction all` のフルワークフローを、`-s/--scan-lists` で原子間距離を駆動して単一構造から実行します。段階的スキャン → MEP 精密化 → （任意）TS 最適化 + IRC が自動で連鎖します。

## 事前に必要なもの

- 入力構造: `.pdb`
- 対象状態に対応した電荷（`-q/--charge` または `-l/--ligand-charge`）・多重度（`-m`）

---

## 方法 A: `-s/--scan-lists` インラインリテラル

`-s/--scan-lists` はコマンドライン上で Python リテラル文字列を直接受け取ります。

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

タンパク質–リガンド複合体 PDB を入力とする場合は `-c/--center` の指定が必須です。低分子 `.pdb` / `.xyz` では `-c` を省略し、`-l` の代わりに `-q` を直接渡してください。`-m/--multiplicity` のデフォルトは `1`（一重項）ですが、例ではここでも明示的に示しています。

### PDBセレクタ

原子は残基名・残基番号・原子名で指定します。区切り文字はカンマ・スペース・スラッシュなどから選べます:

```bash
"TYR,285,CA"     # カンマ区切り
"TYR 285 CA"     # スペース区切り
"TYR/285/CA"     # スラッシュ区切り
"285,TYR,CA"     # 順序は自由
```

### 複数ステージ

複数のリテラルを渡すと、各リテラルが順番に1ステージとして実行されます:

```bash
# ステージ1: 1つの結合を 1.35 Å に駆動
# ステージ2: 2つの結合を同時に駆動
pdb2reaction -i input.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -m 1 -s \
  '[("TYR,285,CA","SAM,309,C10",1.35)]' \
  '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
  -o ./result_scan
```

ステージは順番に実行され、各ステージは前ステージの緩和結果から開始します。

### クォーティングルール

```bash
# 正しい: 外側をシングルクォート、内側のセレクタ文字列をダブルクォート
-s '[("TYR,285,CA","SAM,309,C10",1.35)]'

# 正しい: 整数インデックスは内側のクォート不要
-s '[(1, 5, 2.0)]'

# 避ける: 外側をダブルクォートにすると内側のエスケープが必要
-s "[(\"TYR,285,CA\",\"SAM,309,C10\",1.35)]"
```

> **Tip:** `--print-parsed` を付けると、スキャン対象が正しくパースされたか実行前に確認できます。

---

## 方法 B: YAML スペックファイル（複雑なスキャン向け）

### 1. `scan.yaml` を作成

ステージごとの目標距離（単位: Å）を順番に定義します。

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "SAM,309,C10", 1.35]]
 - [["TYR,285,CA", "SAM,309,C10", 2.20], ["TYR,285,CB", "SAM,309,C11", 1.80]]
```

### 2. 実行

```bash
pdb2reaction -i input.pdb -q 0 -m 1 -s scan.yaml -o ./result_scan
```

---

## 期待される出力

```text
result_scan/
├── summary.log
├── summary.json
├── scan/
│   ├── preopt/                    # 事前最適化構造
│   ├── stage_01/                  # スキャンステージ 1
│   │   ├── scan_trj.xyz
│   │   └── scan.pdb
│   └── stage_02/                  # ステージ 2（マルチステージ時）
└── path_search/                   # MEP 探索（デフォルト、再帰的）; --refine-path False 時は path_opt/
    ├── mep.pdb
    └── energy_diagram_UMA_all.png
```

**確認ポイント:**

1. `scan/stage_01/scan_trj.xyz` — 結合距離の変化を PyMOL で確認
2. `path_search/mep.pdb` — 最適化後の MEP 軌跡
3. `summary.log` — 障壁高さと結合変化

**ヒント:** `--print-parsed --dry-run` でスキャン設定を事前確認:

```bash
pdb2reaction scan -i input.pdb -q 0 -s '[(1, 5, 1.35)]' --print-parsed --dry-run
```

## 補足

- `-s/--scan-lists` は YAML/JSON ファイルパスまたはインライン Python リテラルを受け付けます（両方を同時に指定することはできません）。
- `pdb2reaction all --help-advanced` で全オプション（スキャン制御を含む）を確認できます。
- 単独の `scan` サブコマンド（MEP 精密化なし）については [scan](scan.md) を参照してください。

## 次の導線

- 全オプション: [all](all.md)
- TS 最適化と検証: [クイックスタート: `pdb2reaction tsopt`](quickstart-tsopt-freq.md)
