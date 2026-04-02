# クイックスタート: `pdb2reaction all`

## 目的

2 つの完全系 PDB（反応物 R と生成物 P）から、end-to-end のワークフローを 1 回実行します。

## 前提条件

- pdb2reaction がインストール済みであること（[インストール](installation.md) を参照）
- **水素原子が追加済み**の 2 つの PDB ファイル（反応物 R と生成物 P）
- すべての入力 PDB で同じ原子が同じ順序で含まれていること

## 最小コマンド

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

後処理（TS 最適化、振動解析・熱化学、DFT 一点計算）まで同時に実行する場合:

```bash
pdb2reaction all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## まず確認する出力

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb`（または `result_all/path_search/seg_*/`）

## 補足

- `pdb2reaction all --help` は主要オプション、`pdb2reaction all --help-advanced` は全オプションを表示します。

## 次の導線

- 単一構造の段階的スキャン: [クイックスタート: `pdb2reaction scan`](quickstart-scan.md)
- TS 最適化と検証: [クイックスタート: `pdb2reaction tsopt`](quickstart-tsopt-freq.md)
- 全オプション: [all](all.md)
