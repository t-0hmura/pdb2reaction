# はじめに

## 概要

`pdb2reaction` は、機械学習原子間ポテンシャル（MLIP: Machine Learning Interatomic Potential）を用いて **PDB 構造** から **酵素反応経路** を自動的に構築する Python 製の CLI ツールキットです。

多くのケースで、次のような **1 コマンド** から反応経路の初期案を得られます。
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

---
さらに `--tsopt --thermo --dft` を追加すると、**MEP 探索 → TS 最適化 → IRC（固有反応座標）→ 熱化学解析 → DFT 一点計算** までまとめて実行できます。
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

入力として、(i) 反応順に並べたタンパク質–リガンド複合体の PDB を 2 つ以上（R → … → P）、(ii) `--scan-lists` を指定した 1 つの PDB、または (iii) TS 候補 1 構造 + `--tsopt` を与えると、`pdb2reaction` が次を自動化します。

- ユーザーが指定した基質の周辺から **活性部位ポケット**（抽出範囲）を切り出し、計算用の **クラスターモデル**（Cluster Model）を構築
- Growing String Method (GSM) や Direct Max Flux (DMF) などの経路最適化手法で **最小エネルギー経路 (MEP: Minimum Energy Path)** を探索
- 必要に応じて **遷移状態（TS: Transition State）** を最適化し、**IRC（固有反応座標: Intrinsic Reaction Coordinate）計算**・**振動解析**・**DFT 一点計算** を実行

ポテンシャルエネルギー面（PES: Potential Energy Surface）の計算には機械学習原子間ポテンシャル（MLIP）を用います。デフォルトのバックエンドは Meta の **UMA** ですが、`--backend` により **ORB**、**MACE**、**AIMNet2** も選択できます。`--solvent` オプションで xTB ベースの暗黙溶媒補正を適用することも可能です。想定される主な用途は以下の通りです。

- DFT 等の量子化学計算では検証に時間がかかる規模の**反応機構解析の試行錯誤**
- 量子化学計算に向けた**初期構造の作成**（反応物・TS・生成物のクラスターモデル）
- 基質バリアントや酵素変異体にわたる**反応経路のハイスループット計算**

一連の処理は CLI から呼び出せるように統一されており、手作業を最小化して **多段階の酵素反応メカニズム** を組み立てられるように設計しています。抽出を行わない全系ワークフロー（`--center/-c` と `--ligand-charge` を省略）では `.xyz` / `.gjf` 入力も利用できます。小分子系にもそのまま適用可能です。

**HPC クラスターやマルチ GPU 環境**では、UMA 推論をノード間で並列化することで、大規模なクラスターモデル（必要なら **完全なタンパク質–リガンド複合体**）にもスケールできます。`workers` と `workers_per_node` で並列度を設定してください（詳細は [MLIP バックエンド](uma_pysis.md)）。`--backend` により代替バックエンド（ORB、MACE、AIMNet2）を選択することもできます。

```{important}
- 入力 PDB ファイルには**水素原子**が含まれている必要があります。
- 複数の PDB を提供する場合、**同じ原子が同じ順序**で含まれている必要があります（座標のみが異なる状態）。一致しない場合はエラーになります。
```

```{tip}
初めて使う場合は、まず [概念とワークフロー](concepts.md) を読むと全体像が掴みやすいです。
症状から切り分ける場合は、まず [典型エラー別レシピ](recipes_common_errors.md) を参照してください。
セットアップや実行でエラーに遭遇したら [トラブルシューティング](troubleshooting.md) も参照してください。
```

### CLI の規約

| 規約 | 例 | 備考 |
|-----|-----|------|
| **残基セレクタ** | `'SAM,GPP'`, `'A:123,B:456'` | 複数値はシェル展開防止のためクォート |
| **電荷マッピング** | `--ligand-charge 'SAM:1,GPP:-3'` | コロンで名前と電荷を区切り、カンマでエントリを区切る |
| **原子セレクタ** | `'TYR,285,CA'` または `'TYR 285 CA'` | 区切り文字: 空白、カンマ、スラッシュ、バッククォート、バックスラッシュ |

詳細は [CLI 規約](cli_conventions.md) を参照してください。

補足: CLI サブコマンド名は `path-search`（ハイフン区切り）ですが、ドキュメントファイル名は [`path_search.md`](path_search.md)（アンダースコア区切り）です。


### 水素原子付与の推奨ツール

PDB に水素原子がない場合は、pdb2reaction を実行する前に次のいずれかを使ってください。

| ツール | コマンド例 | 備考 |
|--------|------------|------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速、結晶構造に広く使用 |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | 水素を追加し部分電荷を割り当て |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用ケモインフォマティクスツールキット |

複数の PDB 入力で同一の原子順序を確保するには、すべての構造に同じ水素付与ツールを一貫した設定で適用してください。

```{warning}
このソフトウェアはまだ開発中です。自己責任でご使用ください。
```

---

## 推奨クイックスタート導線

セットアップと依存関係の詳細は [インストール](installation.md) を参照してください。

- [クイックスタート: `pdb2reaction all`](quickstart_all.md)
- [クイックスタート: `pdb2reaction scan` で単一構造の段階的スキャン](quickstart_scan.md)
- [クイックスタート: `pdb2reaction tsopt`（TS 最適化と検証）](quickstart_tsopt_freq.md)

---

## コマンドの基本構成

`pip` でインストールされる `pdb2reaction` コマンドが主な起点です。短縮エイリアス **`p2r`** も利用可能で、すべてのコマンドをどちらの名前でも実行できます。内部的には **Click** ライブラリを使用しており、デフォルトのサブコマンドは `all` です。

つまり:

```bash
pdb2reaction [OPTIONS]...
# は以下と同等
pdb2reaction all [OPTIONS]...
```

`all` は、クラスター抽出から MEP 探索、TS 最適化、IRC、振動解析、DFT 一点計算までを 1 コマンドで一括実行するサブコマンドです。

クラスター抽出を行う場合、ワークフロー全体で共通の重要オプションが 2 つあります:

- `-i/--input`: 1つ以上の**完全構造**（反応物、中間体、生成物）
- `-c/--center`: **基質/抽出中心**の定義方法（例: 残基名または残基ID）

`--center/-c` を省略すると、クラスター抽出はスキップされ、**完全な入力構造**が直接使用されます。

---

## 主要なワークフロー

### 複数構造MEPワークフロー（反応物 → 生成物）

推定反応座標に沿った複数の完全な PDB 構造（例: R → I1 → I2 → P）がすでにある場合に使用します。

**最小例**

```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

**詳細例**

```bash
pdb2reaction -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --out-dir ./result_all --tsopt --thermo --dft
```

処理の流れ:

- 反応順に並んだ 2 つ以上の**完全系**を受け取る
- 各構造から触媒クラスターモデルを抽出
- デフォルトで `path-search` による**再帰的 MEP 探索**を実行（出力は `path_search/`）
- `--no-refine-path` を指定すると**シングルパス** `path-opt` に切り替え
- PDB テンプレートがある場合、クラスターモデル MEP を**完全系**にマージ
- 必要に応じて各セグメントで TS 最適化、IRC、振動解析、DFT 一点計算を実行

ドッキング、MD、手動モデリングなどで中間体を準備できる場合に推奨するモードです。

```{important}
`pdb2reaction` は複数の入力 PDB が**まったく同じ原子を同じ順序**で含むことを前提としています（座標のみが異なる状態）。座標以外のフィールドが入力間で異なる場合はエラーになります。入力 PDB ファイルには**水素原子**も含まれている必要があります。
```

---

### 単一構造 + 段階的スキャン（MEP 精密化への入力）

**1 つの PDB 構造**しかないが、反応に沿ってどの原子間距離が変化するかが分かっている場合に使用します。

`--scan-lists` と一緒に単一の `-i` を指定します:

**最小例**

```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","SAM 309 C10",2.20),("TYR 285 CB","SAM 309 C11",1.80)]' --scan-lists '[("TYR 285 CB","SAM 309 C11",1.20)]'
```

**詳細例**

```bash
pdb2reaction -i SINGLE.pdb -c 'SAM,GPP' --scan-lists '[("TYR 285 CA","SAM 309 C10",2.20),("TYR 285 CB","SAM 309 C11",1.80)]' --scan-lists '[("TYR 285 CB","SAM 309 C11",1.20)]' --multiplicity 1 --out-dir ./result_scan_all --tsopt --thermo --dft
```

補足:

- `--scan-lists` は抽出されたクラスターモデルでの**段階的距離スキャン**を記述
- 各タプル `(i, j, target_Å)` は:
 - `'TYR,285,CA'` のようなPDB原子セレクター文字列（**区切り文字**: スペース/カンマ/スラッシュ/バッククォート/バックスラッシュ ` ` `,` `/` `` ` `` `\`）**または**1始まりの原子インデックス
 - クラスターモデルのインデックスに自動的に再マッピング
- 1 つの `--scan-lists` リテラルで 1 ステージを実行。複数リテラルを渡すと順次ステージとして実行されます。複数ステージは `--scan-lists` を繰り返して指定してください
- 各ステージは `stage_XX/result.pdb` を書き出し、候補中間体または生成物として扱われる
- デフォルトの `all` ワークフローは連結されたステージを再帰的 `path-search` で精密化
- `--no-refine-path` を使用すると、シングルパス `path-opt` チェーンを実行し、再帰的精密化をスキップ（マージされた `mep_w_ref*.pdb` なし）

このモードは単一構造から反応経路を構築するのに便利です。

---

### 単一構造 TSOPT のみモード

すでに**遷移状態候補**があり、それを最適化して IRC 計算を行いたい場合に使用します。

PDB を 1 つだけ指定し、`--tsopt` を有効にします:

**最小例**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt
```

**詳細例**

```bash
pdb2reaction -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt --thermo --dft --out-dir ./result_tsopt_only
```

処理の流れ:

- MEP/経路探索は実行しません。
- TS 最適化で**クラスターモデル上の TS** を収束させます。
- 両方向で **IRC** を実行し、両端点を最適化して R および P 極小に緩和します。
- R/TS/P に対して `freq` と `dft` を実行できます。
- UMA、Gibbs、DFT//UMA エネルギーダイアグラムを生成します。

`energy_diagram_*_all.png` や `irc_plot_all.png` などの出力は、トップレベルの `--out-dir` の下にもコピーされます。

```{important}
単一入力実行には **`--scan-lists`**（段階的スキャン → GSM）**または** `--tsopt`（TSOPT のみ）のいずれかが必要です。これらのいずれも指定せずに単一の `-i` のみを渡しても、ワークフローは実行されません。
```

---

## 重要な CLI オプションと動作

以下はワークフロー全体で最もよく使用されるオプションです。

| オプション | 説明 |
|----------|------|
| `-i, --input PATH...` | 入力構造。**2 つ以上の PDB** → MEP 探索; **1 つの PDB + `--scan-lists`** → 段階的スキャン → GSM; **1 つの PDB + `--tsopt`** → TSOPT のみモード |
| `-c, --center TEXT` | 基質/抽出中心を定義。残基名（`'SAM,GPP'`）、残基ID（`A:123,B:456`）、または PDB パスをサポート |
| `--ligand-charge TEXT` | 電荷情報: マッピング（`'SAM:1,GPP:-3'`）または単一整数 |
| `-q, --charge INT` | 総電荷の強制上書き |
| `-m, --multiplicity INT` | スピン多重度（例: 一重項は `1`）。 |
| `--scan-lists TEXT...` | 単一入力実行用の段階的距離スキャン |
| `--out-dir PATH` | トップレベル出力ディレクトリ |
| `--tsopt/--no-tsopt` | TS 最適化と IRC を有効化 |
| `--thermo/--no-thermo` | 振動解析と熱化学を実行 |
| `--dft/--no-dft` | DFT 一点計算を実行 |
| `--refine-path/--no-refine-path` | 再帰的 MEP 精密化（デフォルト: `True`） vs シングルパス |
| `--opt-mode grad\|hess` | `all` でのワークフロープリセット（`grad` -> LBFGS/Dimer、`hess` -> RFO/RS-I-RFO、デフォルト `grad`）。コマンド個別実行では `opt --opt-mode grad|hess`、`tsopt --opt-mode grad|hess` を推奨 |
| `--mep-mode gsm\|dmf` | MEP 手法（デフォルト: `gsm`）: Growing String Method または Direct Max Flux |
| `--hessian-calc-mode Analytical\|FiniteDifference` | ヘシアン行列の計算モード（デフォルト: `FiniteDifference`）。詳細は [MLIP 計算機](uma_pysis.md#ヘシアンモード) を参照 |

すべてのオプションと YAML スキーマについては [all](all.md) および [YAML リファレンス](yaml_reference.md) を参照してください。

---

## サマリーファイル

`pdb2reaction all` を実行すると、以下のサマリーファイルが出力されます。

- `summary.log` – 人間が読みやすい形式の結果要約
- `summary.yaml` – 機械可読な YAML 形式の結果要約

主な記載内容:

- 実行した CLI コマンド
- MEP 全体の統計（最大障壁、経路長など）
- セグメントごとの障壁高さと主要な結合変化
- UMA、熱化学、DFT 後処理で得られたエネルギー（有効な場合）

`path_search/` 配下の各セグメントディレクトリにも `summary.log` と `summary.yaml` があり、個別のセグメントの精密化結果を確認できます。

---

## CLI サブコマンド

ほとんどのユーザーは `pdb2reaction all` を主に使用しますが、`pdb2reaction opt` などの個別サブコマンドも利用できます。各サブコマンドは `-h/--help` に対応しています。
`pdb2reaction all --help` は主要オプションのみを表示し、`pdb2reaction all --help-advanced` で全オプションを表示できます。
`scan` / `scan2d` / `scan3d` と計算系サブコマンド（`opt` / `path-opt` / `path-search` / `tsopt` / `freq` / `irc` / `dft`）に加え、`add-elem-info` / `trj2fig` / `energy-diagram` も同様に `--help` は主要オプションのみ、`--help-advanced` で全オプションを表示します。`extract` と `fix-altloc` も段階的 help に対応し、`--help-advanced` で parser の全オプションを表示します。

| サブコマンド | 役割 | ドキュメント |
|------------|------|------------|
| `all` | end-to-endワークフロー | [all](all.md) |
| `extract` | 活性部位ポケットからクラスターモデルを抽出 | [extract](extract.md) |
| `opt` | 構造最適化 | [opt](opt.md) |
| `tsopt` | 遷移状態最適化 | [tsopt](tsopt.md) |
| `path-opt` | MEP最適化 (GSM/DMF) | [path_opt](path_opt.md) |
| `path-search` | 再帰的 MEP 探索 | [path_search](path_search.md) |
| `scan` | 1D結合長スキャン | [scan](scan.md) |
| `scan2d` | 2D距離スキャン | [scan2d](scan2d.md) |
| `scan3d` | 3D距離スキャン | [scan3d](scan3d.md) |
| `irc` | IRC 計算 | [irc](irc.md) |
| `freq` | 振動解析 | [freq](freq.md) |
| `dft` | DFT 一点計算 | [dft](dft.md) |
| `trj2fig` | エネルギープロファイルプロット | [trj2fig](trj2fig.md) |
| `energy-diagram` | 数値から状態エネルギーダイアグラムを描画 | [energy-diagram](energy_diagram.md) |
| `fix-altloc` | PDB代替位置指示子の解決 | [fix_altloc](fix_altloc.md) |
| `add-elem-info` | PDB元素カラム修復 | [add_elem_info](add_elem_info.md) |

```{tip}
ヘシアン評価モードの詳細は [MLIP 計算機](uma_pysis.md#ヘシアンモード) を参照してください。
```

---

## 早見表

**よく使うコマンドパターン:**

```bash
# 基本的な MEP 探索（2 構造以上）
pdb2reaction -i R.pdb P.pdb -c 'SUBSTRATE' --ligand-charge 'SUB:-1'

# 後処理付きフルワークフロー
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft

# 単一構造 + 段階的スキャン
pdb2reaction -i SINGLE.pdb -c 'LIG' --scan-lists '[("RES1,100,CA","LIG,200,C1",2.0)]'

# TS のみ最適化
pdb2reaction -i TS.pdb -c 'LIG' --tsopt --thermo
```

**主要オプション一覧:**

| オプション | 用途 |
|----------|------|
| `-i` | 入力構造 |
| `-c` | 活性部位ポケット抽出用の基質定義 |
| `--ligand-charge` | 基質電荷（例: `'SAM:1,GPP:-3'`） |
| `--tsopt` | TS 最適化 + IRC を有効化 |
| `--thermo` | 振動解析を実行 |
| `--dft` | DFT 一点計算を実行 |
| `--out-dir` | 出力ディレクトリ |

---

## ヘルプ

任意のサブコマンドについて:

```bash
pdb2reaction <subcommand> --help
pdb2reaction <subcommand> --help-advanced
pdb2reaction all --help-advanced
```

`all` では `--help` は短縮版です。全オプションを確認するときは `--help-advanced` を使用してください。UMA バックエンドの詳細オプションについては [MLIP バックエンド](uma_pysis.md) を参照してください。

問題が発生した場合は、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) でIssueを開いてください。
