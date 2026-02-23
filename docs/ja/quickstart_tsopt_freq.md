# クイックスタート: `pdb2reaction tsopt` -> `pdb2reaction freq`

## 目的

TS 候補を最適化し、振動解析で妥当性を確認します。

## 事前に必要なもの

- TS 候補構造: `ts_guess.pdb`
- 対象状態に対応した電荷・多重度（`-q`, `-m`）

## 1. TS 最適化

```bash
pdb2reaction tsopt -i ts_guess.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

## 2. 最適化結果に対して振動解析

```bash
pdb2reaction freq -i ./result_tsopt/final_geometry.pdb -q 0 -m 1 --out-dir ./result_freq
```

## まず確認する出力

- `result_tsopt/final_geometry.pdb`
- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz` と `result_freq/mode_*.pdb`

一次鞍点として妥当な TS では、虚振動（負の cm^-1）は 1 本になることが期待されます。

## 補足

- VRAM に余裕がある場合は `--hessian-calc-mode Analytical` を推奨します。
- 全オプションは `pdb2reaction tsopt --help-advanced` と `pdb2reaction freq --help-advanced` を参照してください。

## 次の導線

- 反応経路追跡は [irc](irc.md)、一括実行は [all](all.md) を参照してください。
