# クイックスタート: `pdb2reaction scan` + `--spec`

## 目的

単一構造に対して、YAML で定義した距離ターゲットへ段階的にスキャンします。

## 事前に必要なもの

- 入力構造: `input.pdb`
- 対象状態に対応した電荷・多重度（`-q`, `-m`）

## 1. `scan.yaml` を作成

```yaml
one_based: true
stages:
  - [["TYR,285,CA", "MMT,309,C10", 1.35]]
  - [["TYR,285,CA", "MMT,309,C10", 2.20], ["TYR,285,CB", "MMT,309,C11", 1.80]]
```

## 2. 実行

```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml --print-parsed --out-dir ./result_scan
```

## まず確認する出力

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- `--dump` 指定時は `scan.trj` / `scan.pdb`

## 補足

- 推奨は `--spec` です。legacy 互換として `--scan-lists` も残っています。
- 詳細オプションは `pdb2reaction scan --help-advanced` で確認できます。

## 次の導線

- 経路精密化は [all](all.md) または [path-search](path_search.md) を参照してください。
