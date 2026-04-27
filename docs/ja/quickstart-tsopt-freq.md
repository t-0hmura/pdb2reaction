# クイックスタート: `pdb2reaction tsopt`

## 目的

TS 候補を最適化し、一次鞍点（first-order saddle point）であることを確認します。

## 事前に必要なもの

- TS 候補構造: `.pdb`
- 対象状態に対応した電荷（`-q/--charge` または `-l/--ligand-charge`）・多重度（`-m`）

## 最小コマンド

```bash
pdb2reaction tsopt -i ts_guess.pdb -q 0 -m 1 --out-dir ./result_tsopt
```

`tsopt` は最適化の最後に自動でヘシアン行列（Hessian）の計算と虚振動数の確認を実行します。ターミナル出力で以下のような行を確認してください。

```
[Imaginary modes] n=1  ([-593.1])
```

### （任意）個別の振動解析

全振動モードの一覧や熱化学補正（零点エネルギー (ZPE)、ギブズ自由エネルギーなど; `all` コマンドの `--thermo` に相当）が必要な場合は、別途 `freq` を実行してください。虚振動数の確認だけであれば、上記の `tsopt` の出力で十分です。

```bash
pdb2reaction freq -i ./result_tsopt/final_geometry.pdb -q 0 -m 1 --out-dir ./result_freq
```

## 期待される出力

```text
result_tsopt/
├── final_geometry.xyz     # 最適化済み TS
├── final_geometry.pdb     # PDB 形式（入力が PDB の場合）
└── vib/
    ├── imag_-593.10cm-1.pdb       # 虚振動モードアニメーション
    └── imag_-593.10cm-1_trj.xyz
```

**確認ポイント:**

1. ターミナル出力: **`n=1`** かつ |振動数| >= 100 cm⁻¹ で化学的に有効な一次鞍点（注: `tsopt` の内部虚振動検出は |ν| >= 5 cm⁻¹ で数えるため、100 cm⁻¹ は TS 品質の*判定基準*であり検出閾値ではありません。詳細は [典型エラー別レシピ](recipes-common-errors.md) を参照）
2. `vib/imag_*.pdb` — PyMOL でアニメーション; 期待される結合の切断/形成に対応する振動であること
3. `n=0`（虚振動なし）: 極小に収束。別の初期構造を試す
4. `n>1`（複数の虚振動）: `--flatten` で余分なモードの平坦化を試行（反復回数の上限は `--config` で YAML を渡し、`hessian_dimer.flatten_max_iter` を指定 — Dimer と RS-I-RFO の両モードともこのキーを参照します）

## 補足

- VRAM に余裕がある場合は `--hessian-calc-mode Analytical` を使用してください（デフォルトは `FiniteDifference`）。
- 全オプションは `pdb2reaction tsopt --help-advanced` と `pdb2reaction freq --help-advanced` を参照してください。

## 次の導線

- 反応経路追跡は [irc](irc.md)、一括実行は [all](all.md) を参照してください。
