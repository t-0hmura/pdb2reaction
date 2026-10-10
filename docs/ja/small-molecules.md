# 小分子の反応機構解析

`pdb2reaction` は、酵素のクラスターモデルだけでなく、小分子の反応の経路も解析できます。反応物（R）と生成物（P）の XYZ/GJF/PDB/CIF 構造を渡し、`-c` を付けずに `all` を実行します。`-c` が無いので、クラスターモデルの切り出しはスキップされ、入力構造そのままで解析が行われます。

## 最初の例

[`examples/aromatic_claisen/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples/aromatic_claisen) に、芳香族 Claisen 転位（アリルフェニルエーテル → 6-アリルシクロヘキサ-2,4-ジエン-1-オン）の反応物と生成物があります。以下のコマンドは `examples/` の中で実行します。

```bash
pdb2reaction all -i aromatic_claisen/reactant.xyz aromatic_claisen/product.xyz -q 0 --tsopt --thermo
```

最小エネルギー経路（MEP）を探索し、遷移状態（TS）の最適化・固有反応座標（IRC）・振動解析と熱化学補正まで進みます。端末の最後のほうの `====== Pipeline summary ======` の下に `Scientific status: success` と出れば成功で、`summary.json` の `scientific_status` にも同じ値が入ります。

## 入力のポイント

* **形式**: XYZ・GJF・PDB・CIF を使えます。複数の構造は、同じ原子が同じ順に並んでいる必要があります。
* **電荷とスピン**: 電荷は `-q`、スピン多重度は `-m`（デフォルトは 1）で指定します。
* **ほかの入力モード**: 1 つの構造からスキャンで経路を作る [Scan-list モード](quickstart-scan.md) と、TS の候補から始める [TS-only モード](quickstart-tsopt.md) も、同じように `-c` を付けずに使えます。
* **溶液中の反応**: `--solvent` で xTB の溶媒和の補正を足せます（[xTB 溶媒補正](backends.md#xtb-溶媒補正)）。

## 関連ドキュメント

* [`all`](all.md) — オプションの説明
* [クイックスタート: `pdb2reaction all`](quickstart-all.md) — 反応の前後の構造から始める流れ
