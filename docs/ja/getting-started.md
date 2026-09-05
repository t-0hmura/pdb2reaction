# はじめに

## 概要

<img src="../overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` は、機械学習原子間ポテンシャル（MLIP: Machine Learning Interatomic Potential）を用いて **PDB / mmCIF 構造から酵素反応経路の候補を探索する** Python 製の CLI ツールキットです。MLIP は DFT 参照データ（エネルギー・原子間力、および周期境界条件の学習データを持つ foundation model では応力テンソルも）で学習されたニューラルネットワークモデルで、DFT のポテンシャルエネルギー曲面をごくわずかな計算コストで近似します。

多くのケースでは、次のような **1 コマンド** で反応経路の初期案を得られます。
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

---
さらに `--tsopt --thermo --dft` を追加すると、**最小エネルギー経路（MEP: Minimum Energy Path）探索 → 遷移状態（TS: Transition State）最適化 → 固有反応座標（IRC: Intrinsic Reaction Coordinate） → 熱化学補正 → DFT 一点計算** までまとめて実行できます。
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

> **実行例:** [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) ディレクトリに GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）の完全な `all` ワークフロースクリプト（MEP およびスキャンパイプライン）があります。

入力として、(i) 反応順に並べた PDB / mmCIF 構造を 2 つ以上（R → … → P）、(ii) `--scan-lists/-s` を指定した 1 構造、または (iii) TS 候補 1 構造 + `--tsopt` を与えると、`pdb2reaction` が次を自動化します。

- ユーザーが指定した基質の周辺から **活性部位モデル（バインディングポケット）** を切り出し、計算用の **クラスターモデル**（Cluster Model）を構築
- Growing String Method (GSM) や Direct Max Flux (DMF) などの経路最適化手法で **最小エネルギー経路 (MEP: Minimum Energy Path)** を探索
- 必要に応じて **遷移状態（TS: Transition State）** を最適化し、**IRC（固有反応座標: Intrinsic Reaction Coordinate）計算**・**振動解析**・**DFT 一点計算** を実行

計算には機械学習原子間ポテンシャル（MLIP）を用います。デフォルトのバックエンドは Meta の **UMA** ですが、`-b/--backend` により **ORB**、**MACE**、**AIMNet2** も選択できます。クラスターモデルの TS 最適化・IRC 検証・QRRHO 熱化学を単一の GPU で実行できるかは、系の大きさ、バックエンド／モデル、Hessian の計算法、計算精度、ハードウェアに依存します。

- DFT 等の量子化学計算では検証に時間がかかる規模の**反応機構解析の試行錯誤**
- 量子化学計算に向けた**初期構造の作成**（反応物・ TS ・生成物のクラスターモデル）
- 基質バリアントや酵素変異体にわたる**反応経路の大量計算**

本 CLI は多段階反応経路の候補を探索し、TS・IRC 計算による検証を支援します。小分子系や、ユーザーが自分で構築したクラスターモデルにも適用できます。活性部位抽出を行わない場合は `--center/-c` を省略します。総電荷は `-q`、PDB/mmCIF 全体へ適用する `-l`、YAML の `calc.charge`、または `.gjf` の設定から決定します。

**HPC クラスターやマルチ GPU 環境**では、`workers` と
`workers_per_node` により UMA 推論をノード間で並列化できます。扱える系の規模は
モデル、VRAM、Hessian の計算法、通信負荷に依存し、全タンパク質を扱えるとは
限りません（[MLIP バックエンド](uma-pysis.md)）。

### パイプライン概要

`all` サブコマンドは以下のステージを自動実行します:

```text
入力構造
  |
  v
[extract]  `-c` 指定時のみ活性部位モデルを抽出
  |
  v
[scan]  `--scan-lists/-s` 指定時のみ段階的距離拘束スキャン
  |
  v
[path-opt/path-search]  TS-only 以外で MEP 探索
  |
  v
[tsopt]  `--tsopt` 指定時のみ TS 最適化
  |
  v
[irc]  `--tsopt` 指定時のみ固有反応座標
  |
  v
[freq]  `--tsopt --thermo` 指定時のみ振動解析 + 熱化学
  |
  v
[dft]  `--tsopt --dft` 指定時のみ一点 DFT エネルギー
```

名前の付いた計算ステージ（`extract`、scan 系、path 系、`tsopt`、`irc`、`freq`、`dft`）は単独のサブコマンドとしても実行できます。`all` は選択されたステージを統合し、`summary.json` と `summary.log` を出力します。

### 主要な出力ファイル

| ファイル | 説明 |
|---------|------|
| `summary.json` | 反応障壁、エネルギー、結合変化、環境情報 |
| `summary.log` | ディレクトリツリー付きテキストサマリ |
| `segments/seg_XX/` | 反応セグメントの後処理結果。R/TS/P 構造は TSOPT・IRC・端点処理の成功後に生成 |
| `mep.pdb` / `mep.cif` | MEP 軌跡。内部変換した入力の元の鎖・残基 ID は CIF に保持 |
| `energy_diagram_*.png` | エネルギープロファイル図（電子/Gibbs 補正） |

```{important}
- 入力 PDB / mmCIF には**水素原子**が含まれている必要があります。
- 複数構造には**同じ原子が同じ順序**で含まれている必要があります（座標のみが異なる状態）。
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
| **原子セレクタ** | `'TYR,285,CA'` または `'A:TYR:285:CA'` | 重複時は位置固定の `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` |

詳細は [CLI 規約](cli-conventions.md) を参照してください。

### 水素原子付与の推奨ツール

PDB に水素原子がない場合は、pdb2reaction を実行する前に次のいずれかを使ってください。

| ツール | コマンド例 | 備考 |
|--------|------------|------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速、結晶構造に広く使用 |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用化学情報処理ツールキット |
| **PyMOL** | `h_add`（PyMOL 内） | 分子可視化ツール（水素付加機能あり） |
| **tleap** (AmberTools) | `tleap -f leapin` | Amber 力場準備ツール |

複数の PDB 入力で同一の原子順序を確保するには、すべての構造に同じ水素付与ツールを一貫した設定で適用してください。

## 推奨クイックスタート導線

セットアップと依存関係の詳細は [インストール](installation.md) を参照してください。

- [クイックスタート: `pdb2reaction all`](quickstart-all.md)
- [クイックスタート: スキャンを起点とする `pdb2reaction all`](quickstart-scan.md)
- [クイックスタート: TS のみモード（`pdb2reaction all --tsopt`）](quickstart-tsopt-freq.md)
- [Colab GUI](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb) — PDB/mmCIF upload、3Dでchain付き選択、Validate、実行

## コマンドの基本構成

`pip install pdb2reaction` で `pdb2reaction` と短縮名 `p2r` が使えます。どちらも同じコマンドを実行し、サブコマンドを省略すると `all` が選ばれます。

つまり:

```bash
pdb2reaction [OPTIONS]...
# は以下と同等
pdb2reaction all [OPTIONS]...
```

[`all`](all.md) は、クラスター抽出、MEP 探索、TS 最適化、振動解析、必要に応じた DFT までを 1 コマンドで一括実行するサブコマンドです。

クラスター抽出を行う場合、ワークフロー全体で共通の重要オプションが 2 つあります:

- `-i/--input`: 1 つ以上の**完全構造**（反応物、中間体、生成物）
- `-c/--center`: **基質/抽出中心**の定義方法（例: 残基名または残基 ID）

`--center/-c` を省略すると、クラスター抽出はスキップされ、**完全な入力構造**が直接使用されます。

## 主要なワークフロー

| モード | 概要 | クイックスタート |
|------|-----|--------------|
| **複数構造 MEP**（2 つ以上の PDB） | 反応座標に沿った複数の PDB（R → … → P）を受け取り、各構造のクラスターモデル抽出 → MEP 探索（デフォルトは単一パス path-opt; `--refine-path` で再帰的 path-search）→ 必要に応じてセグメントごとに TS / IRC / freq / DFT を実行 | [クイックスタート: `pdb2reaction all`](quickstart-all.md) |
| **単一構造 + スキャン定義**（1 構造 + `--scan-lists/-s`） | 1 つの構造を距離拘束スキャンにかけ、各ステージを単一パス `path-opt`（`--refine-path` で再帰的 `path-search`）に渡して MEP を構築 | [クイックスタート: スキャン起点ワークフロー](quickstart-scan.md) |
| **単一構造 TSOPT のみ**（1 PDB + `--tsopt`） | MEP/経路探索を完全にスキップし、TS 候補を最適化 → 双方向 IRC → 端点最適化、必要なら R/TS/P に freq / DFT を実行 | [クイックスタート: TS のみモード](quickstart-tsopt-freq.md) |

```{important}
単一入力実行には **`--scan-lists/-s`**（段階的スキャン → GSM）**または** `--tsopt`（TSOPT のみ）のいずれかが必要です。これらのいずれも指定せずに単一の `-i` のみを渡しても、ワークフローは実行されません。
```

## 重要な CLI オプションと動作

以下はワークフロー全体で最もよく使用されるオプションです。

| オプション | 説明 |
|----------|------|
| `-i, --input PATH...` | 入力構造。**2 つ以上の PDB** → MEP 探索; **1 つの PDB + `--scan-lists/-s`** → 段階的スキャン → GSM; **1 つの PDB + `--tsopt`** → TSOPT のみモード |
| `-c, --center TEXT` | 基質/抽出中心を定義。残基名（`'SAM,GPP'`）、残基 ID（`A:123,B:456`）、または PDB パスをサポート |
| `-l, --ligand-charge TEXT` | 電荷情報: マッピング（`'SAM:1,GPP:-3'`）または単一整数 |
| `-q, --charge INT` | 総電荷の強制上書き |
| `-m, --multiplicity INT` | スピン多重度（例: 一重項は `1`） |
| `--tsopt` | TS 最適化と IRC を有効化 |
| `-b, --backend TEXT` | MLIP バックエンドの選択（`uma`, `orb`, `mace`, `aimnet2`） |

オプションの完全な一覧は [CLI 規約](cli-conventions.md) と [自動生成 CLI リファレンス](../reference/commands/index.md) を参照してください。

## サマリーファイル

実行全体の結果は、出力ディレクトリ直下の次のファイルで確認します。入力検証などで早期に停止すると、作成されない場合があります。

- `summary.log` – 結果要約
- `summary.json` – JSON 結果

主な記載内容:

- 実行した CLI コマンド
- MEP 全体の統計（最大障壁、経路長など）
- セグメントごとの障壁高さと主要な結合変化
- MLIP バックエンド、熱化学、DFT 後処理で得られたエネルギー（有効な場合）

`segments/seg_NN/` には、指定した後処理の結果が保存されます。各セグメントに実行全体のサマリーがあるとは限らないため、出力ディレクトリ直下の `summary.log` / `summary.json` と各処理の結果を確認してください。

## CLI サブコマンド

ほとんどのユーザーは `pdb2reaction all` を主に使います。CLI は個別サブコマンドも提供しており、各コマンドは `-h/--help` に対応しています（計算/スキャン/抽出/ユーティリティ系は `--help-advanced` で全オプションを表示）。サブコマンド一覧と各ドキュメントへのリンクは [ドキュメントトップ](index.md#サブコマンド) を参照してください。

## エージェントスキル

`pdb2reaction` は `skills/` に AI エージェント向けの指示書を同梱しており、CLI サブコマンド、構造 I/O（PDB / mmCIF / XYZ / GJF）、バックエンドインストール（UMA / Orb / MACE / AIMNet2 / DFT / xtb）、標準的なワークフロー、出力解析、HPC 運用をカバーしています。Claude Code・Codex・Cursor などのエージェントで使うには、`skills/` をプロジェクト直下に（例: Claude Code なら `.claude/skills/` に）または `~/.claude/skills/` にコピーしてください。

## ヘルプ

任意のサブコマンドについて:

```bash
pdb2reaction <subcommand> --help
pdb2reaction <subcommand> --help-advanced
pdb2reaction all --help-advanced
```

MLIP バックエンドの詳細オプションについては [MLIP バックエンド](uma-pysis.md) を参照してください。

問題が発生した場合は、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) で Issue を開いてください。
