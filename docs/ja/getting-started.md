# はじめに

## 概要

<img src="../overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` は、機械学習原子間ポテンシャル（MLIP: Machine Learning Interatomic Potential）を用いて **PDB 構造から酵素反応経路をモデリング** する Python 製の CLI ツールキットです。

多くのケースで、次のような **1 コマンド** から反応経路の初期案を得られます。
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

---
さらに `--tsopt --thermo --dft` を追加すると、**最小エネルギー経路（MEP: Minimum Energy Path）探索 → 遷移状態（TS: Transition State）最適化 → 固有反応座標（IRC: Intrinsic Reaction Coordinate） → 振動解析・熱化学 → DFT 一点計算** までまとめて実行できます。
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

> **実行例:** [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) ディレクトリに GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）の完全な `all` ワークフロースクリプト（MEP およびスキャンパイプライン）があります。

入力として、(i) 反応順に並べたタンパク質–リガンド複合体の PDB を 2 つ以上（R → … → P）、(ii) `--scan-lists/-s` を指定した 1 つの PDB、または (iii) TS 候補 1 構造 + `--tsopt` を与えると、`pdb2reaction` が次を自動化します。

- ユーザーが指定した基質の周辺から **活性部位モデル（バインディングポケット）** を切り出し、計算用の **クラスターモデル**（Cluster Model）を構築
- Growing String Method (GSM) や Direct Max Flux (DMF) などの経路最適化手法で **最小エネルギー経路 (MEP: Minimum Energy Path)** を探索
- 必要に応じて **遷移状態（TS: Transition State）** を最適化し、**IRC（固有反応座標: Intrinsic Reaction Coordinate）計算**・**振動解析**・**DFT 一点計算** を実行

計算には機械学習原子間ポテンシャル（MLIP）を用います。デフォルトのバックエンドは Meta の **UMA** ですが、`-b/--backend` により **ORB**、**MACE**、**AIMNet2** も選択できます。想定される主な用途は以下の通りです。

- DFT 等の量子化学計算では検証に時間がかかる規模の**反応機構解析の試行錯誤**
- 量子化学計算に向けた**初期構造の作成**（反応物・TS・生成物のクラスターモデル）
- 基質バリアントや酵素変異体にわたる**反応経路の大量計算**

本 CLI は最小限の手動設定で**多段階の酵素反応機構**を生成します。小分子系にもそのまま適用可能です。抽出を行わない全系ワークフロー（`--center/-c` と `--ligand-charge/-l` を省略）では `.xyz` / `.gjf` 入力も利用できます。

**HPC クラスターやマルチ GPU 環境**では、UMA 推論をノード間で並列化することで、大規模なクラスターモデル（必要なら **完全なタンパク質–リガンド複合体**）にも対応できます。`workers` と `workers_per_node` で並列度を設定してください（詳細は [MLIP バックエンド](uma-pysis.md)）。

### パイプライン概要

`all` サブコマンドは以下のステージを自動実行します:

```text
PDB (R, P)
  |
  v
[extract]  活性部位モデル抽出（クラスターモデル）
  |
  v
[scan]  （オプション, --scan-lists/-s）段階的距離拘束スキャン
  |
  v
[path-search]  MEP 探索（再帰的 path-search、デフォルト; `--refine-path False` で path-opt に切替）
  |
  v
[tsopt]  TS 最適化 (RS-I-RFO; 代替として Dimer)
  |
  v
[irc]  固有反応座標
  |
  v
[freq]  振動解析 + 熱化学 (R, TS, P)
  |
  v
[dft]  一点 DFT エネルギー（オプション, --dft）
```

各ステージは単独のサブコマンドとしても実行できます。`all` はこれらを統合し、`summary.json` と `summary.log` を出力します。

### 主要な出力ファイル

| ファイル | 説明 |
|---------|------|
| `summary.json` | 機械可読な結果（障壁、エネルギー、結合変化、環境情報） |
| `summary.log` | ディレクトリツリー付きテキストサマリ |
| `seg_XX/` | 反応ステップごとの IRC 最適化 R/TS/P 構造 |
| `mep.pdb` | PyMOL/VMD で表示可能な MEP 軌跡 |
| `energy_diagram_*.png` | エネルギープロファイル図（電子/Gibbs 補正） |

```{important}
- 入力 PDB ファイルには**水素原子**が含まれている必要があります。
- 複数の PDB を提供する場合、**同じ原子が同じ順序**で含まれている必要があります（座標のみが異なる状態）。一致しない場合はエラーになります。
```

```{tip}
症状から切り分ける場合は、まず [典型エラー別レシピ](recipes-common-errors.md) を参照してください。
セットアップや実行でエラーに遭遇したら [トラブルシューティング](troubleshooting.md) も参照してください。
```

### CLI の規約

| 規約 | 例 | 備考 |
|-----|-----|------|
| **残基セレクタ** | `'SAM,GPP'`, `'A:123,B:456'` | 複数値はシェル展開防止のためクォート |
| **電荷マッピング** | `-l 'SAM:1,GPP:-3'` | コロンで名前と電荷を区切り、カンマでエントリを区切る |
| **原子セレクタ** | `'TYR,285,CA'` または `'TYR 285 CA'` | 区切り文字: 空白、カンマ、スラッシュ、バッククォート、バックスラッシュ |

詳細は [CLI 規約](cli-conventions.md) を参照してください。

### 水素原子付与の推奨ツール

PDB に水素原子がない場合は、pdb2reaction を実行する前に次のいずれかを使ってください。

| ツール | コマンド例 | 備考 |
|--------|------------|------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速、結晶構造に広く使用 |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | 水素を追加し部分電荷を割り当て |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用化学情報処理ツールキット |
| **PyMOL** | `h_add`（PyMOL 内） | 分子可視化ツール（水素付加機能あり） |
| **tleap** (AmberTools) | `tleap -f leapin` | Amber 力場準備ツール |

複数の PDB 入力で同一の原子順序を確保するには、すべての構造に同じ水素付与ツールを一貫した設定で適用してください。

## 推奨クイックスタート導線

セットアップと依存関係の詳細は [インストール](installation.md) を参照してください。

- [クイックスタート: `pdb2reaction all`](quickstart-all.md)
- [クイックスタート: `pdb2reaction all --scan-lists/-s` で単一構造の段階的スキャン+MEP+TS](quickstart-scan.md)
- [クイックスタート: `pdb2reaction tsopt`（TS 最適化と検証）](quickstart-tsopt-freq.md)

## コマンドの基本構成

`pip` でインストールされる `pdb2reaction` コマンドが主な起点です。短縮エイリアス **`p2r`** も `pdb2reaction` パッケージが同じ setuptools entry point で登録しており（`pip install pdb2reaction` 直後から両方利用可能）、すべてのコマンドをどちらの名前でも実行できます。内部的には **Click** ライブラリを使用しており、デフォルトのサブコマンドは `all` です。

つまり:

```bash
pdb2reaction [OPTIONS]...
# は以下と同等
pdb2reaction all [OPTIONS]...
```

[`all`](all.md) は、クラスター抽出、MEP 探索、TS 最適化、振動解析、必要に応じた DFT までを 1 コマンドで一括実行するサブコマンドです。

クラスター抽出を行う場合、ワークフロー全体で共通の重要オプションが 2 つあります:

- `-i/--input`: 1つ以上の**完全構造**（反応物、中間体、生成物）
- `-c/--center`: **基質/抽出中心**の定義方法（例: 残基名または残基ID）

`--center/-c` を省略すると、クラスター抽出はスキップされ、**完全な入力構造**が直接使用されます。

## 主要なワークフロー

| モード | 概要 | クイックスタート |
|------|-----|--------------|
| **複数構造 MEP**（2 つ以上の PDB） | 反応座標に沿った複数の PDB（R → … → P）を受け取り、各構造のクラスターモデル抽出 → 再帰的 MEP 探索 → 必要に応じてセグメントごとに TS / IRC / freq / DFT を実行。 | [クイックスタート: `pdb2reaction all`](quickstart-all.md) |
| **単一構造 + 段階的スキャン**（1 PDB + `--scan-lists/-s`） | 1 つの PDB をクラスターモデル上で段階的距離スキャンにかけ、各ステージを再帰的 `path-search`（`--refine-path False` で単一パス `path-opt`）に渡して MEP を構築。 | [クイックスタート: 単一構造の段階的スキャン](quickstart-scan.md) |
| **単一構造 TSOPT のみ**（1 PDB + `--tsopt`） | MEP/経路探索を完全にスキップし、TS 候補を最適化 → 双方向 IRC → 端点最適化、必要なら R/TS/P に freq / DFT を実行。 | [クイックスタート: TS 最適化](quickstart-tsopt-freq.md) |

```{important}
単一入力実行には **`--scan-lists/-s`**（段階的スキャン → GSM）**または** `--tsopt`（TSOPT のみ）のいずれかが必要です。これらのいずれも指定せずに単一の `-i` のみを渡しても、ワークフローは実行されません。
```

## 重要な CLI オプションと動作

以下はワークフロー全体で最もよく使用されるオプションです。

| オプション | 説明 |
|----------|------|
| `-i, --input PATH...` | 入力構造。**2 つ以上の PDB** → MEP 探索; **1 つの PDB + `--scan-lists/-s`** → 段階的スキャン → GSM; **1 つの PDB + `--tsopt`** → TSOPT のみモード |
| `-c, --center TEXT` | 基質/抽出中心を定義。残基名（`'SAM,GPP'`）、残基ID（`A:123,B:456`）、または PDB パスをサポート |
| `-l, --ligand-charge TEXT` | 電荷情報: マッピング（`'SAM:1,GPP:-3'`）または単一整数 |
| `-q, --charge INT` | 総電荷の強制上書き |
| `-m, --multiplicity INT` | スピン多重度（例: 一重項は `1`）。 |
| `--tsopt/--no-tsopt` | TS 最適化と IRC を有効化 |
| `-b, --backend TEXT` | MLIP バックエンドの選択（`uma`, `orb`, `mace`, `aimnet2`） |

オプションの完全な一覧は [CLI 規約](cli-conventions.md) と [自動生成 CLI リファレンス](../reference/commands/index.md) を参照してください。

## サマリーファイル

`pdb2reaction all` を実行すると、以下のサマリーファイルが出力されます。

- `summary.log` – 結果要約
- `summary.json` – JSON 結果

主な記載内容:

- 実行した CLI コマンド
- MEP 全体の統計（最大障壁、経路長など）
- セグメントごとの障壁高さと主要な結合変化
- MLIP バックエンド、熱化学、DFT 後処理で得られたエネルギー（有効な場合）

`path_search/`（`--refine-path False` 使用時は `path_opt/`）配下の各セグメントディレクトリにも `summary.log` と `summary.json` があり、個別のセグメントの精密化結果を確認できます。

## CLI サブコマンド

ほとんどのユーザは `pdb2reaction all` を主に使います。CLI は個別サブコマンドも提供しており、各コマンドは `-h/--help` に対応しています（計算/スキャン/抽出/ユーティリティ系は `--help-advanced` で全オプションを表示）。サブコマンド一覧と各ドキュメントへのリンクは [ドキュメントトップ](index.md#サブコマンド) を参照してください。

## エージェントスキル

`pdb2reaction` は `.claude/skills/` に AI エージェント向けの指示書を同梱しており、CLI サブコマンド、構造 I/O（PDB / XYZ / GJF）、バックエンドインストール（UMA / Orb / MACE / AIMNet2 / DFT / xtb）、標準的なワークフロー、出力解析、HPC 運用をカバーしています。Claude Code や Cursor などのエージェントプラットフォームで使うには、`.claude/skills/` をプロジェクトリポジトリまたは `~/.claude/skills/` にコピーしてください。

## ヘルプ

任意のサブコマンドについて:

```bash
pdb2reaction <subcommand> --help
pdb2reaction <subcommand> --help-advanced
pdb2reaction all --help-advanced
```

MLIP バックエンドの詳細オプションについては [MLIP バックエンド](uma-pysis.md) を参照してください。

問題が発生した場合は、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) でIssueを開いてください。
