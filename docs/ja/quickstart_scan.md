# クイックスタート: `pdb2reaction scan` — 単一構造の段階的スキャン

## 目的

単一構造に対して、原子間距離を段階的にスキャン（拘束付き緩和）します。

## 事前に必要なもの

- 入力構造: `input.pdb`
- 対象状態に対応した電荷・多重度（`-q`, `-m`）

---

## 方法 A: `--spec`（YAMLファイル、複雑なスキャンに推奨）

### 1. `scan.yaml` を作成

ステージごとの目標距離（単位: Angstrom）を順番に定義します。

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "SAM,309,C10", 1.35]]
 - [["TYR,285,CA", "SAM,309,C10", 2.20], ["TYR,285,CB", "SAM,309,C11", 1.80]]
```

### 2. 実行

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --out-dir ./result_scan
```

---

## 方法 B: `--scan-lists`（CLIインライン指定）

`--scan-lists` はコマンドライン上で Python リテラル文字列を直接受け取ります。

### 基本構文

各リテラルは `(atom1, atom2, target_distance)` のタプルリストです（距離の単位は Angstrom）。1 リテラル = 1 ステージ。

```bash
# 単一ステージ、整数原子インデックス（デフォルトで1-based）
pdb2reaction scan -i input.pdb -q 0 --scan-lists '[(1, 5, 1.35)]' --out-dir ./result_scan

# 単一ステージ、PDBセレクタ文字列
pdb2reaction scan -i input.pdb -q 0 --scan-lists '[("TYR,285,CA", "SAM,309,C10", 1.35)]' --out-dir ./result_scan
```

### PDBセレクタ

原子は残基名・残基番号・原子名で指定します。区切り文字は柔軟に選べます:

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
pdb2reaction scan -i input.pdb -q 0 --scan-lists \
  '[("TYR,285,CA","SAM,309,C10",1.35)]' \
  '[("TYR,285,CA","SAM,309,C10",2.20),("TYR,285,CB","SAM,309,C11",1.80)]' \
  --out-dir ./result_scan
```

ステージは順番に実行され、各ステージは前ステージの緩和結果から開始します。

### クォーティングルール

```bash
# 正しい: 外側をシングルクォート、内側のセレクタ文字列をダブルクォート
--scan-lists '[("TYR,285,CA","SAM,309,C10",1.35)]'

# 正しい: 整数インデックスは内側のクォート不要
--scan-lists '[(1, 5, 2.0)]'

# 避ける: 外側をダブルクォートにすると内側のエスケープが必要
--scan-lists "[(\"TYR,285,CA\",\"SAM,309,C10\",1.35)]"
```

> **Tip:** `--print-parsed` を付けると、スキャン対象が正しくパースされたか実行前に確認できます。

---

## まず確認する出力

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`（複数ステージの場合）
- `--dump` 指定時は `scan_trj.xyz` / `scan.pdb`

## 補足

- `--spec` と `--scan-lists` は排他的です（どちらか一方のみ使用）。
- 詳細オプションは `pdb2reaction scan --help-advanced` で確認できます。
- 入力フォーマットの詳細は [scan](scan.md) を参照してください。

## 次の導線

- 経路精密化は [all](all.md) または [path-search](path_search.md) を参照してください。
