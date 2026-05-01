# 凍結原子（Frozen Atoms）

## 概要

> **要約:** クラスターモデルでは切り出し境界の少数の原子を固定しないと、最適化器が末端のフラグメントを非物理的な構造へ引き込んでしまいます。`pdb2reaction` は **リンク水素**（`extract` が切断結合に付加）と、3 通りの `freeze_atoms` 指定でこれを扱います。

`extract`サブコマンドで、タンパク質から残基を切り出すと、境界の結合は **リンク水素**（残基 `LKH`、原子 `HL`、元の結合ベクトル方向に 1.09 Å）でキャップされます。リンク水素の親原子を自由なまま放置すると、勾配降下がキャップと親原子を一緒に動かして境界形状が変形します。関連原子を凍結すると、最適化・ MEP 探索・ IRC ・振動解析を通じて境界が固定されます。

## 凍結原子の 3 通りの指定方法

### 1. `--freeze-links/--no-freeze-links`（デフォルト `True`）

`extract` が付加したリンク水素の親原子を自動凍結。下流の全サブコマンドでデフォルト有効です。

```bash
pdb2reaction extract -i complex.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -o model.pdb
pdb2reaction opt -i model.pdb -q 0 -m 1   # --freeze-links は True、LKH 親原子を自動凍結
```

XYZ/GJF 入力には `LKH` レコードがないため `--freeze-links` は無効で、次の 2 つの方法を使ってください。`--ref-pdb FILE` を渡すと XYZ/GJF 実行が PDB トポロジーを継承し、リンク水素検出が復活します。

### 2. `--freeze-atoms 'i,j,k,...'`（CLI 明示指定）

カンマ区切りの **1 始まり** 原子インデックス。任意の入力形式に適用可。`--freeze-links` と併用すると和集合が凍結されます。

```bash
pdb2reaction tsopt -i ts_candidate.xyz -q 0 -m 1 \
  --freeze-atoms '12,15,28,29,42'
```

### 3. YAML `geom.freeze_atoms`（`--config` 経由）

```yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42]   # 1 始まりインデックス
```

長いリストや、設定ファイルに同梱したい場合に有用。CLI と YAML のリストは置換ではなくマージされます。

```bash
pdb2reaction tsopt -i ts.xyz -q 0 -m 1 --config tsopt.yaml
```

## 3 つのソースの組み合わせ方

実行時に凍結される原子集合は以下の **和集合**:

- YAML `geom.freeze_atoms`（`--config FILE`）
- CLI `--freeze-atoms`
- `--freeze-links` 有効時に検出された `LKH` 親原子

どれかを優先する仕組みはなく、いずれかのソースに登録された原子はすべて凍結されます。

## 計算への効果

- **力（Force）:** 凍結 DOF に対する力をゼロ化（最適化器は動かせません）。
- **ヘシアン:** 凍結 DOF の行・列は除去される（`calc.return_partial_hessian: true`。`freq` ではデフォルト、`irc` では強制）か、フル行列でゼロ化されます。
- **振動解析:** 凍結原子があるとき `freq` は自動で Partial Hessian Vibrational Analysis（PHVA）を実行し、active ブロックのみ対角化します。
- **MEP / IRC:** 凍結原子は全イメージ・全ステップで初期座標を保持します。

## サブコマンド対応表

| サブコマンド | `--freeze-links`（PDB） | `--freeze-atoms`（任意の入力） | YAML `geom.freeze_atoms` |
|---|:---:|:---:|:---:|
| [`extract`](extract.md) | （`LKH/HL` を挿入。フラグは下流用） | n/a | n/a |
| [`opt`](opt.md) | yes | yes | yes |
| [`tsopt`](tsopt.md) | yes | yes | yes |
| [`freq`](freq.md) | yes（PHVA を起動） | yes | yes |
| [`irc`](irc.md) | yes | yes | yes |
| [`path-opt`](path-opt.md) | yes | yes | yes |
| [`path-search`](path-search.md) | yes | yes | yes |
| [`scan`](scan.md) / [`scan2d`](scan2d.md) / [`scan3d`](scan3d.md) | yes | yes | yes |
| [`all`](all.md) | yes | yes | yes |

## よくある落とし穴

- **`LKH/HL` レコードを手動削除した場合**: `--freeze-links` が凍結対象を見つけられません。`--freeze-atoms` で明示指定するか、`extract` を再実行してください。
- **原子インデックスの再番号化**: `--freeze-atoms` と `geom.freeze_atoms` は 1 始まりで入力の原子順に依存します。再抽出で順序が変わったらインデックスを再生成してください。
- **トポロジーなしの XYZ/GJF**: `LKH` レコードが無く `--freeze-links` は no-op です。`--ref-pdb FILE` または明示的な `--freeze-atoms` を使用してください。
- **`--no-freeze-links` の使いどころ**: 自動凍結を切る診断目的のみ。本番のクラスターモデル実行では `--freeze-links` を有効のままにしてください。

## 関連項目

- [`extract`](extract.md) — リンク水素の挿入箇所（残基 `LKH`、原子 `HL`、1.09 Å）。アルゴリズム詳細は {ref}`リンク水素と凍結原子 <ja-link-hydrogen-and-frozen-atoms>` を参照。
- [YAML リファレンス](yaml-reference.md) — `geom.freeze_atoms` スキーマとマージ順序。
- [CLI 規約](cli-conventions.md) — サブコマンド共通のフラグ規約。
- [用語集](glossary.md) — 活性部位モデル、クラスターモデル、リンク水素。
