# 概念とワークフロー

このページでは、`pdb2reaction` の主要な概念——ポケット、テンプレート、セグメント、イメージ（image）——と、`all` コマンドが各サブコマンドをどのように連携させるかを説明します。

---

## ワークフローの全体像

多くのケースでは次の流れになります。

```text
全系入力 (PDB/XYZ/GJF)
 │
 ├─ (任意) ポケット抽出 [extract] ← --center/-c を使う場合は PDB が必要
 │ ↓
 │ ポケット/クラスターモデル (PDB)
 │ │
 │ ├─ (任意) 段階的スキャン [scan] ← 単一構造ワークフロー
 │ │ ↓
 │ │ 順序付けられた中間体
 │ │ ↓
 │ └─ MEP 探索 [path-search] または [path-opt]
 │ ↓
 │ MEP 経路 (mep_trj.xyz) + エネルギーダイアグラム
 │ ↓
 └─ (任意) TS 最適化 + IRC [tsopt] → [irc]
 └─ (任意) 熱化学 [freq]
 └─ (任意) DFT 一点計算 [dft]
```

各ステージはサブコマンドとして単独実行できます。また `pdb2reaction all` を使うと、複数ステージをまとめて実行できます。

```{important}
遷移状態（一次鞍点: First-Order Saddle Point）: HEI や `tsopt` の出力は **TS 候補** として扱い、`irc`（両端が意図した極小構造へ到達すること）で検証してから解釈してください。`tsopt` は内部で虚振動数チェックを実施済みです — 出力に虚振動数が 1 本（|ν| >= 100 cm⁻¹）あることを確認してください。虚振動数が複数残る場合は `--flatten` の適用を検討してください。
```

---

## 重要な概念

### 全系とポケット（クラスターモデル）
- **全系（Full System）**: 元の入力構造（酵素ならタンパク質–基質複合体など）。
- **ポケット（Pocket）**: 基質周辺の抽出範囲。`-c/--center` と `-r/--radius` で定義される。
- **クラスターモデル（Cluster Model）**: ポケットから切り出された計算対象の部分系。リンク水素でキャップされ、MEP/TS 探索に使用します。

ポケット抽出は主に以下で制御します。
- `-c/--center`: 基質の指定（残基ID、残基名、または基質のみのPDB）
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`

### イメージ（Image）とセグメント（Segment）
- **イメージ（Image）**: 経路上の 1 つの構造（1 ノード）。chain-of-states 法で離散化された各点に対応します。
- **セグメント（Segment）**: 隣接する 2 つの端点を結ぶ MEP 区間。多構造入力はセグメントに分解されます（例: R → I1, I1 → I2, ...）。

### テンプレートとファイル変換（`--convert-files`）
`pdb2reaction` は軌跡（trajectory; 例: `mep_trj.xyz`, `irc_trj.xyz`）を出力します。PDB テンプレートや Gaussian 入力がある場合、対応する付随ファイルも出力できます。

- PDB テンプレートがあるとき → `.pdb` 付随ファイル
- Gaussian テンプレートがあるとき → `.gjf` 付随ファイル

この挙動は `--convert-files/--no-convert-files`（デフォルト: `True`）で制御します。

---

## 代表的な3つの使い方

### 1) 複数構造の MEP 探索（R → … → P）
すでに反応座標に沿った **2つ以上** の構造がある場合。

例:

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### 2) 単一構造の段階的スキャン → MEP
入力が **1つ** しかないが、スキャンにより端点系列を生成できる場合。

例:

```bash
pdb2reaction -i holo.pdb -c '308,309' \
 --scan-lists '[("TYR,285,CA","SAM,309,C10",2.20)]'
```

### 3) TSOPT のみ（ポケット TS 最適化）
TS 候補がすでにある、あるいは1構造で TS 最適化だけ試したい場合。

例:

```bash
pdb2reaction -i ts_guess.pdb -c 'SAM,GPP' --tsopt
```

---

## `all` と個別サブコマンドの使い分け

### `pdb2reaction all` を推奨するケース
- extract → MEP → TSOPT/IRC → freq/DFT まで **まとめて** 実行したい
- 出力ディレクトリやログ管理を1コマンドに寄せたい

### 個別サブコマンドを推奨するケース
- あるステージだけを **切り出してデバッグ** したい（例: `extract` のみ）
- 独自の端点生成など、ワークフローを **組み替えたい**

---

## CLI の注意点

```{important}
- ブール値オプションは `--flag` / `--no-flag` と `--flag True/False`（`yes/no`, `1/0` 含む）の両方を受理します。新規スクリプトでは toggle 形式を推奨します。
- 複数 PDB を与える場合、各ファイルは **同じ原子が同じ順序** で並んでいることが重要です（座標だけが異なる）。
- 酵素の反応機構解析では、水素を含んだ入力 PDB を用意することを強く推奨します。
```

---

## 次に読むページ

### 入門
- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting_started.md) — クイックスタートとワークフロー概要
- [典型エラー別レシピ](recipes_common_errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法

### 主要サブコマンド
| サブコマンド | 用途 | ドキュメント |
|------------|------|-------------|
| `all` | end-to-endワークフロー | [all.md](all.md) |
| `extract` | ポケット抽出 | [extract.md](extract.md) |
| `path-search` | 再帰的 MEP 探索 | [path_search.md](path_search.md) |
| `tsopt` | TS 最適化 | [tsopt.md](tsopt.md) |
| `freq` | 振動解析 | [freq.md](freq.md) |
| `dft` | DFT 一点計算 | [dft.md](dft.md) |

### リファレンス
- [YAML リファレンス](yaml_reference.md) — 全オプションの YAML 設定
- [用語集](glossary.md) — 用語リファレンス
