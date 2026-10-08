# はじめに

## 概要

<img src="../overview.png" alt="pdb2reaction workflow overview" width="90%">

`pdb2reaction` は、機械学習原子間ポテンシャル（MLIP）を活用し、**PDB / mmCIF 構造から酵素の反応経路候補を自動探索する** Python 製 CLI ツールキットです。

DFT（密度汎関数法）の計算データを学習したニューラルネットワークを用いることで、DFT レベルのポテンシャルエネルギー曲面をごくわずかな計算コストで近似し、高速な経路探索を実現します。

多くのケースでは、次のような **1 コマンド** で反応経路の初期案を得られます。

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

端末の出力の最後のほうに `Scientific status: success` と出れば、求めた段はすべて収束しています。

---

さらに `--tsopt --thermo --dft` を付けると、**最小エネルギー経路（MEP）探索 → 遷移状態（TS）最適化 → 固有反応座標（IRC） → 振動解析・熱化学補正 → DFT 一点計算** までを一貫して自動実行できます。

```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```

---

> **実行例:** [`examples/`](https://github.com/t-0hmura/pdb2reaction/tree/main/examples) ディレクトリに、上のコマンドで使う構造（`1.R.pdb`、`3.P.pdb`）と、GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）を題材とした一連のワークフロースクリプト（MEP 探索とスキャン）を用意しています。[インストール](installation.md)の後、`git clone https://github.com/t-0hmura/pdb2reaction && cd pdb2reaction/examples` で取得し、その中で上のコマンドを実行してください。

### 主な用途

* DFT 等の量子化学計算では検証に時間がかかる規模の**反応機構解析の試行錯誤**
* 量子化学計算に向けた**初期構造の作成**（反応物・TS・生成物のクラスターモデル）
* 基質バリアントや酵素変異体にわたる**反応経路の大量計算**

### 主な自動化機能

入力として「(1) 反応順に並べた複数の PDB 構造（R → … → P）」「(2) 単一構造 ＋ 距離スキャン指定」「(3) 単一構造 ＋ TS 最適化指定」のいずれかを与えることで、以下を自動処理します。

1. **クラスターモデル構築**: 指定した基質周辺から活性部位（バインディングポケット）を自動切り出し
2. **最小エネルギー経路（MEP）探索**: Growing String Method (GSM) や Direct Max Flux (DMF) による経路探索
3. **高精度検証**: 遷移状態（TS）の構造最適化、IRC 計算、振動解析、DFT 一点計算

MLIP で妥当な経路が見つかったら、その TS をそのまま DFT での TS 構造最適化にもっていくことにも `pdb2reaction` は対応しています。TS 最適化 → IRC → 端点の最適化 → 振動数計算のワークフローを、GPU4PySCF を用いることで GPU で高速化された DFT 計算により実行可能です。詳しくは [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。

---

## ワークフローとパイプライン

### パイプラインの流れ

全工程を一括実行する `all` サブコマンド（デフォルト動作）は、以下のステージを順次実行します。

```text
入力構造（PDB / mmCIF）
  │
  ▼
[extract] 抽出ステージ: -c 指定時のみ活性部位モデルを切り出し
  │
  ▼
[scan] スキャンステージ: -s 指定時のみ段階的距離拘束スキャンを実施
  │
  ▼
[path-opt / path-search] 経路探索: TS-only モード以外で MEP（最小エネルギー経路）を探索
  │
  ▼
[tsopt] TS 最適化: --tsopt 指定時のみ遷移状態を精密化
  │
  ▼
[irc] IRC 計算: --tsopt 指定時のみ固有反応座標を追跡し、端点を最適化
  │
  ▼
[freq] 振動解析: --tsopt --thermo 指定時のみ熱化学補正を計算
  │
  ▼
[dft] DFT 一点計算: --tsopt --dft 指定時のみ DFT エネルギーを算出
```

各ステージは [`extract`](extract.md)、[`tsopt`](tsopt.md)、[`irc`](irc.md) などのサブコマンドとして単独でも実行できます。一覧は [サブコマンド](index.md#サブコマンド) にあります。

---

## クイックスタート導線

環境構築の詳細は [インストールガイド](installation.md) を参照してください。

* **Web ブラウザで手軽に試す**: [Colab GUI ノートブック](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb)（3D プレビューで残基を選択）
* **複数の PDB 構造から始める**: [クイックスタート: `pdb2reaction all`](quickstart-all.md)
* **1 つの PDB 構造からスキャンで探索する**: [クイックスタート: `pdb2reaction all --scan-lists`](quickstart-scan.md)
* **TS 候補構造を最適化・検証する**: [クイックスタート: TS-only モード](quickstart-tsopt.md)

---

## コマンドの基本構成

インストール後は `pdb2reaction` および短縮コマンド `p2r` が利用できます。サブコマンドを省略した場合、自動的に `all` が呼び出されます。

```bash
# 以下の 2 つは同一の処理を行います
pdb2reaction [OPTIONS]...
pdb2reaction all [OPTIONS]...
```

### 入力モードの選び方

| 実行モード | 入力条件 | 主な動作 |
| --- | --- | --- |
| **複数構造 MEP 探索** | 2 つ以上の PDB（`-i R.pdb P.pdb`） | 各構造から同じクラスターモデルを切り出し（`-c` のとき）、その間の MEP を探索 |
| **単一構造 ＋ スキャン** | 1 つの PDB ＋ `--scan-lists`（`-s`） | 指定結合の距離を段階的に変化させて経路を生成 |
| **TS-only モード** | 1 つの PDB ＋ `--tsopt` | MEP 探索をスキップし、TS 候補の最適化・IRC を直接実行 |

> **注意:** 単一構造を入力するときは、`--scan-lists/-s` か `--tsopt` が必要です。どちらも無いとエラーになります。

---

## 基本的な CLI オプション

| オプション | 引数の例 | 説明 |
| --- | --- | --- |
| `-i, --input` | `1.R.pdb 3.P.pdb` | 入力構造ファイル（PDB / mmCIF）。複数指定可能 |
| `-c, --center` | `'SAM,GPP'` / `'A:SAM:123'` | 抽出中心（基質残基名・残基 ID・入力と同じ座標の基質だけを入れた PDB ファイル）。省略時は切り出しを行わず構造全体を使用 |
| `-l, --ligand-charge` | `'SAM:1,GPP:-3'` | リガンドごとの形式電荷マッピング（標準残基とイオンの電荷は自動で数えます） |
| `-q, --charge` | `-2` | 抽出モデル全体の総電荷（自動判定を上書きする場合に指定） |
| `-m, --multiplicity` | `1` | スピン多重度（デフォルト: `1`、一重項） |
| `--tsopt` | （フラグ） | TS 最適化と IRC 計算を有効化 |
| `--thermo` | （フラグ） | 振動解析と QRRHO（準剛体ローター・調和振動子）モデルによる熱化学補正を実行（`--tsopt` と併用） |
| `--dft` | （フラグ） | 得られた構造に対して一点 DFT 計算を実行（`--tsopt` と併用） |
| `-b, --backend` | `uma` / `orb` / `mace` | 使用する MLIP バックエンドを指定（デフォルト: `uma`） |

`--dft` には、{ref}`詳細なインストール手順 <ja-step-by-step-installation>` の手順 7 で入れる DFT 用の追加パッケージが要ります。

構文ルールの詳細は [共通オプションと残基・原子の指定](cli-conventions.md)、全オプションの一覧は [`all` の CLI リファレンス](../reference/commands/all.md) を参照してください。

---

## 入力構造に関する重要事項

### 1. 水素原子の付加（必須）

入力構造には**全原子の水素が含まれている必要があります**。結晶構造など水素が欠落している構造を使用する場合は、事前に以下のツール等で付加してください。

| 推奨ツール | コマンド例 | 特徴 |
| --- | --- | --- |
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速で結晶構造の水素付加に広く使われる |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用的な化学情報処理ツール |
| **PyMOL** | コマンドラインで `h_add` | ビジュアルを確認しながら付加可能 |
| **tleap** (AmberTools) | `tleap -f leapin` | Amber 力場に基づく精密な付加 |

`all` は空の元素欄（77–78 列）を自分で埋めます。`extract` などのコマンドを単独で使う前には、[`add-elem-info`](add-elem-info.md) で埋めてください。代替位置（altLoc）は、どのコマンドも PDB を読むときに残基ごとに {ref}`平均占有率の最も高いもの <ja-mmcif-input>` を 1 つ選びます。選んだ結果をファイルに残すには [`fix-altloc`](fix-altloc.md) を使います。

### 2. 原子の並び順の一致（複数構造入力時）

反応物（R）や生成物（P）など複数の構造を入力する場合、**すべての構造で同一の原子が同じ順序で並んでいる必要があります**（座標値のみが異なる状態）。水素付加ツールを用いる際は、すべての構造に対して同一の設定で処理してください。反応で別の残基に移る原子も、R での残基名と原子名のままにします。同梱例では、GPP から Glu186 に移る水素は `3.P.pdb` でも `GPP 321` の `H11` です。

mmCIF（`.cif`・`.mmcif`）と、残基が 10,000 以上や原子が 99,999 以上の大きな構造も、PDB と同じコマンドで扱えます。詳しくは {ref}`mmCIF の入力 <ja-mmcif-input>` を参照してください。

---

## 出力ファイルの構成

実行が終わると、`-o` の出力ディレクトリに次のファイルができます。既定は `./result_all/` です。主なファイルは [出力ディレクトリのレイアウト](output-layout.md)、`summary.json` のすべての欄は {ref}`JSON 出力リファレンス <ja-summary-json-path-search-all>` にあります。

| 出力ファイル / フォルダ | 内容 |
| --- | --- |
| `summary.log` | テキスト形式のサマリー（ディレクトリ構成、各段階の進行状況） |
| `summary.json` | 機械可読形式の結果（反応障壁、各状態のエネルギー、結合変化） |
| `energy_diagram_*.png` | 生成されたエネルギープロファイル図（電子エネルギー / Gibbs 補正） |
| `mep_trj.pdb` / `mep_trj.cif` | 最小エネルギー経路（MEP）のアニメーション軌跡ファイル |
| `segments/seg_NN/` | 反応セグメントごとの詳細結果（最適化された R/TS/P 構造、IRC 軌跡など。`--tsopt` のとき） |

端末の出力の最後のほうにある `====== Pipeline summary ======` の下の `Scientific status:` の行は、`summary.json` の `scientific_status` と同じ値です。求めた段がすべて収束すると `success` になり、`--tsopt` のときは TS の n_imag が 1 であることも条件です。そうでなければ `partial` か `failed` になり、理由は `scientific_status_reasons` に出ます。虚振動のモードが狙った結合を動かしているか、IRC の両端が狙った R と P かは、`segments[].bond_changes` を見て自分で確かめてください。開くファイルは各クイックスタートにあります。

---

## AI エージェント連携（Skills）

`pdb2reaction` には、AI エージェント（Claude Code、Codex、Cursor など）向けの設定指示書が `skills/` ディレクトリに同梱されています。

CLI サブコマンドの仕様、PDB/mmCIF/XYZ/GJF の入出力ルール、バックエンドの環境構築手順、HPC 並列化のベストプラクティスが定義されています。`skills/` をエージェントに読み込ませることで、エージェントを通じた自然言語指示による計算実行・解析が可能になります。配置場所とスキルの一覧は [`skills/README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/README.md) を参照してください。MCP のクライアントからコマンドをツールとして呼ぶ方法は [MCP サーバー](mcp_server.md) にあります。

---

## トラブルシューティングとサポート

実行中にエラーが発生した場合は、以下のドキュメントを参照してください。

* {ref}`トラブルシューティング <ja-troubleshooting-quick-table>`: エラー症状別の対処法と、インストールや環境起因の不具合の解決手順
* [MLIP バックエンド](backends.md): GPU メモリや並列計算に関する詳細。複数の GPU ノードで動かすジョブスクリプトは [HPC 実行例](hpc-example.md)

コマンドの全オプションを確認したい場合は、ヘルプオプションを利用してください。

```bash
pdb2reaction <subcommand> --help
pdb2reaction all --help-advanced
```

解決しない問題やバグの報告は、[GitHub Issues](https://github.com/t-0hmura/pdb2reaction/issues) にて受け付けています。
