# 概念とワークフロー

このページでは、`pdb2reaction` を使ううえでの **全体像（メンタルモデル）** を説明します。  
「ポケット」「テンプレート」「セグメント」「画像（image）」が何を指すのか、そしてトップレベルの `all` が各サブコマンドをどう組み合わせるのかを把握するためのページです。

---

## パイプラインの全体像

多くのケースでは次の流れになります。

```
全系入力 (PDB/XYZ/GJF)
   │
   ├─ (任意) ポケット抽出        [extract]     ← --center/-c を使う場合は PDB が必要
   │        ↓
   │   ポケット/クラスタモデル (PDB)
   │        │
   │        ├─ (任意) 段階的スキャン [scan]     ← 単一構造ワークフロー
   │        │        ↓
   │        │   順序付けられた中間体
   │        │
   │        └─ MEP 探索           [path-search] または [path-opt]
   │                 ↓
   │            MEP 軌跡 (mep.trj) + エネルギーダイアグラム
   │
   └─ (任意) TS 精密化 + IRC      [tsopt] → [irc]
             └─ (任意) 熱化学     [freq]
             └─ (任意) DFT SP    [dft]
```

各ステージはサブコマンドとして単独実行できます。また `pdb2reaction all` を使うと、複数ステージをまとめて実行できます。

---

## 重要な概念

### 全系とポケット（クラスタモデル）
- **全系**: 元の入力構造（酵素ならタンパク質–基質複合体など）。
- **ポケット / クラスタモデル**: 基質周辺を切り出した小系。MEP/TS 探索を軽量化します。

ポケット抽出は主に以下で制御します。
- `-c/--center`: 基質の指定（残基ID、残基名、または基質のみのPDB）
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`

### 画像（image）とセグメント
- **画像（image）**: 経路上の1つの幾何（1ノード）。
- **セグメント**: 2つの端点を結ぶ MEP。多構造入力は、隣接する端点ペアごとのセグメントに分解されます。

### テンプレートとファイル変換（`--convert-files`）
`pdb2reaction` は軌跡（例: `mep.trj`, `irc.trj`）を出力します。  
PDB テンプレートや Gaussian 入力がある場合、必要に応じて companion ファイルも出力できます。

- PDB テンプレートがあるとき: `.pdb` companion
- Gaussian テンプレートがあるとき: `.gjf` companion

この挙動は `--convert-files {True|False}`（デフォルト: `True`）で制御します。

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
  --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### 3) TSOPT のみ（ポケット精密化）
TS 候補がすでにある、あるいは1構造で TS 精密化だけ試したい場合。

例:

```bash
pdb2reaction -i ts_guess.pdb -c 'SAM,GPP' --tsopt True
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

## 知っておくと詰まりにくい CLI の癖

```{important}
- ブール値オプションは `True`/`False` を明示して渡します（例: `--tsopt True`）。
- 複数 PDB を与える場合、各ファイルは **同じ原子が同じ順序** で並んでいることが重要です（座標だけが異なる）。
- 酵素の反応機構解析では、水素を含んだ入力 PDB を用意することを強く推奨します。
```

---

## 次に読むページ

- [はじめに](getting-started.md) — インストールと初回実行
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [`all`](all.md) — エンドツーエンドのラッパー
- [`path-search`](path_search.md) — 再帰的 MEP 探索の詳細
- [YAMLリファレンス](yaml-reference.md) — 全オプションのYAML設定
