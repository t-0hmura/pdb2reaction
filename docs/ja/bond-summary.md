# `bond-summary`

連続する分子構造間の共有結合変化を検出・報告します。

## 書式

```bash
pdb2reaction bond-summary -i R.xyz P.xyz
pdb2reaction bond-summary -i R.xyz TS.xyz P.xyz
pdb2reaction bond-summary -i R.pdb IM1.pdb IM2.pdb P.pdb
```

## 説明

`bond-summary` は入力構造の連続ペアを比較し、形成または切断された結合を報告します。*N* 個の入力ファイルに対して *N − 1* 個の比較ブロック（A→B, B→C, …）を生成します。

結合認識は元素固有の共有結合半径とスケーリングファクターを使用します。距離はオングストローム（Å）で報告されます。

対応フォーマット: **XYZ**、**PDB**、**GJF**（拡張子から自動検出）。

## オプション

| オプション | 説明 | デフォルト |
|--------|-------------|---------|
| `-i, --input FILE` | 入力構造ファイル（2つ以上必須） | — |
| `--device TEXT` | 計算デバイス（`cpu`, `cuda`） | `cpu` |
| `--bond-factor FLOAT` | 共有結合半径の和に対するスケーリングファクター | `1.20` |
| `--one-based / --zero-based` | 出力の原子インデックス規約 | `--one-based` |
| `--out-json / --no-out-json` | **機械可読な JSON を標準出力へ出力**（テキストレポートの代わり）。このサブコマンドは `result.json` ファイルを書き出しません（永続化したい場合は stdout をリダイレクトしてください）。このオプションは **`pdb2reaction bond-summary --help` には表示されません**（`--help-advanced` で確認できます）。スキーマは [JSON 出力スキーマ → bond-summary](json-output.md#bond-summary) を参照。 | `False` |

## 例

### 2構造の比較

```bash
pdb2reaction bond-summary -i 1.R.xyz 3.P.xyz
```

出力:
```
============================================================
  1.R.xyz  →  3.P.xyz
============================================================
Bond formed (2):
  - O14-H106 : 1.502 Å --> 1.011 Å
  - P95-O107 : 3.477 Å --> 1.523 Å
Bond broken (2):
  - P95-O97 : 1.585 Å --> 3.270 Å
  - H106-O107 : 1.034 Å --> 1.673 Å
```

### 多構造比較（反応経路）

```bash
pdb2reaction bond-summary -i 1.R.xyz 3.IM1.xyz 5.IM2.xyz 7.P.xyz
```

R→IM1, IM1→IM2, IM2→P の3つの比較ブロックが出力されます。

## 注意事項

- すべての入力構造は**同一の原子数と元素順序**である必要があります。
- 結合検出は `all` ワークフローの IRC 端点検証で使用される内部 `bond_changes` モジュールと同じアルゴリズムを使用します。
- 境界的な結合（例: 金属配位 2.0–2.4 Å）の感度を調整するには、`--bond-factor` を大きくしてください（例: `1.30`）。

## Python API

結合変化検出関数は Python から直接利用できます:

```python
from pdb2reaction.bond_changes import compare_structures, has_bond_change, summarize_changes
```

| Function | Description |
|----------|-------------|
| `compare_structures(geom1, geom2, device="cuda", bond_factor=1.20)` | Detect covalent bonds formed or broken between two pysisyphus geometries. Returns a `BondChangeResult` with `formed_covalent` and `broken_covalent` (sets of 0-based index pairs). |
| `has_bond_change(geom_start, geom_end, bond_cfg)` | Convenience wrapper: returns `(has_changes: bool, summary_text: str)`. `bond_cfg` accepts keys: `device`, `bond_factor`, `margin_fraction`, `delta_fraction`. |
| `summarize_changes(geom, result, one_based=True)` | Format bond changes as a text report with distances in Angstrom. |

### 例

```python
from pysisyphus.helpers import geom_loader
from pdb2reaction.bond_changes import has_bond_change

geom_r = geom_loader("R.xyz")
geom_p = geom_loader("P.xyz")

changed, summary = has_bond_change(geom_r, geom_p, {"device": "cpu", "bond_factor": 1.20})
if changed:
    print(summary)
```

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [irc](irc.md) -- 結合検出で端点検証される IRC 軌跡
- [all](all.md) -- 内部で結合変化検証を使用する一気通貫ワークフロー
- [trj2fig](trj2fig.md) -- 軌跡からエネルギープロファイルを可視化
