# クイックスタート: `pdb2reaction scan` — 段階的スキャン

## 目的

単一構造に対して、結合距離を段階的にスキャンします。
入力形式は **`--spec`**（YAML ファイル、推奨）と **`--scan-lists`**（インライン Python リテラル）の 2 種類があります。

## 事前に必要なもの

- 入力構造: `input.pdb`
- 対象状態に対応した電荷・多重度（`-q`, `-m`）

---

## 方法 A: `--spec`（YAML ファイル — 推奨）

### 1. `scan.yaml` を作成

```yaml
one_based: true
stages:
 - [["TYR,285,CA", "MMT,309,C10", 1.35]]
 - [["TYR,285,CA", "MMT,309,C10", 2.20], ["TYR,285,CB", "MMT,309,C11", 1.80]]
```

### 2. 実行

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --print-parsed --out-dir ./result_scan
```

---

## 方法 B: `--scan-lists`（インライン記法）

YAML ファイルを用意せず、同じ 2 ステージスキャンをコマンドラインだけで表現できます:

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 \
 --scan-lists '[("TYR,285,CA","MMT,309,C10",1.35)]' \
              '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]' \
 --print-parsed --out-dir ./result_scan
```

`--scan-lists` の後ろに並べた各 Python リテラル文字列が 1 ステージに対応します。
**フラグを繰り返さず**、すべてのステージを連続する引数として渡してください。

---

## まず確認する出力

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- `--dump` 指定時は `scan_trj.xyz` / `scan.pdb`

## 使い分けガイド

| | `--spec`（YAML） | `--scan-lists`（インライン） |
|---|---|---|
| **向いている用途** | 多段・複雑なスキャン | 簡単な 1 ステージ実行、スクリプト |
| **可読性** | 読みやすく、バージョン管理しやすい | シェルのクォートが煩雑になりがち |
| **再現性** | YAML ファイル自体が記録になる | コマンドライン全体を保存する必要あり |

## 補足

- 詳細オプションは `pdb2reaction scan --help-advanced` で確認できます。
- フォーマットの詳細は [scan](scan.md) を参照してください。

## 次の導線

- 経路精密化は [all](all.md) または [path-search](path_search.md) を参照してください。
