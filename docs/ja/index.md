---
orphan: true
---

# [pdb2reaction]{.p2r-wordmark} ドキュメント

:::{container} p2r-hero-meta
[バージョン: v{{ release }}]{.p2r-pill} [GitHub](https://github.com/t-0hmura/pdb2reaction){.p2r-meta-gh} [ACS Omega 論文](https://doi.org/10.1021/acsomega.6c08242){.p2r-meta-paper}
:::

:::{container} p2r-hero
<img src="../overview.png" alt="pdb2reaction ワークフロー概要" class="p2r-hero-figure">

{.p2r-tagline}
**pdb2reaction** は、機械学習原子間ポテンシャル（MLIP）を使用して、酵素複合体などの PDB 構造や小分子の XYZ 構造などから反応機構解析を行うための Python 製 CLI ツールキットです。

{.p2r-lead}
初めての方は [はじめに](getting-started.md) からお読みください。

{.p2r-cta}
[はじめに](getting-started.md){.p2r-btn .p2r-btn-primary} [インストール](installation.md){.p2r-btn .p2r-btn-install} [Google Colabで実行](https://colab.research.google.com/github/t-0hmura/pdb2reaction/blob/main/examples/pdb2reaction_colab.ipynb){.p2r-btn .p2r-btn-colab}
:::

## クイックスタート

::::{container} p2r-cards
:::{container} p2r-card p2r-card-endpoint
**反応の前後の構造から反応機構解析を一気通貫で行う**

<!-- p2r-mode-stages endpoint -->

[クイックスタート: all の Endpoint モード](quickstart-all.md)
:::

:::{container} p2r-card p2r-card-scan
**1 つの構造から一気通貫で反応機構解析を行う**

<!-- p2r-mode-stages scan -->

[クイックスタート: all の Scan-list モード](quickstart-scan.md)
:::

:::{container} p2r-card p2r-card-tsonly
**TS 構造から一気通貫で反応機構解析を行う**

<!-- p2r-mode-stages tsonly -->

[クイックスタート: TS-only モード](quickstart-tsopt.md)
:::
::::

| 目的 | ページ |
|------|------|
| **クラスターモデルを組む・削る・広げる** | [クラスターモデルの組み方](model-setup.md) |
| **反応機構を調べる・TS が取れない** | [反応機構を調べるコツ](mechanism-tips.md) |
| **求めた TS 構造を DFT で構造最適化する** | [求めた TS 構造を DFT で構造最適化する](dft-backend.md) |
| **計算が失敗した** | [トラブルシューティング](troubleshooting.md) |

## サブコマンド

<!-- p2r-stage-strip -->

| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | 抽出（任意）と、3 つの[入力モード](getting-started.md#入力モードの選び方)（複数構造 MEP 探索・単一構造 ＋ スキャン・TS-only モード）のいずれか、任意の TS/IRC・熱化学・DFT を統括 |
| [`extract`](extract.md) | タンパク質–リガンド複合体から活性部位モデル（バインディングポケット）を抽出 |
| [`fix-altloc`](fix-altloc.md) | PDB の代替位置指示子を解決 |
| [`add-elem-info`](add-elem-info.md) | PDB の元素列（77–78）を修復 |
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS または RFO。任意の `--flatten` で残った虚振動を除く） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer または RS-P-RFO。任意の `--flatten` で余分な虚振動を除く） |
| [`path-opt`](path-opt.md) | GSM または DMF による 1 段階の MEP 最適化（2 構造から） |
| [`path-search`](path-search.md) | 自動精密化を伴う多段階の再帰的 MEP 探索（2 構造以上） |
| [`scan`](scan.md) | 拘束付き距離スキャン（複数距離の協奏スキャン・多段階スキャンに対応） |
| [`scan2d`](scan2d.md) | 2 次元のエネルギー地形の探索・PES マッピング |
| [`scan3d`](scan3d.md) | 3 次元のエネルギー地形の探索・PES マッピング |
| [`freq`](freq.md) | 振動解析と熱化学 |
| [`irc`](irc.md) | 固有反応座標（IRC: Intrinsic Reaction Coordinate）計算 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`sp`](sp.md) | MLIP または `-b dft` による一点計算（エネルギー + 力 / Hessian） |
| [`trj2fig`](trj2fig.md) | XYZ 軌跡からエネルギープロファイルをプロット |
| [`energy-diagram`](energy-diagram.md) | 数値入力からエネルギーダイアグラムを作成 |
| [`bond-summary`](bond-summary.md) | 連続構造間の共有結合変化を検出・レポート |

## 設定・リファレンス

| トピック | ページ |
|-------|------|
| **共通オプションと入力要件** | [共通オプションと残基・原子の指定](cli-conventions.md) |
| **原子の固定と距離の拘束（`--freeze-atoms`・`--distance-restraint`）** | {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` |
| **よくあるエラーと対処** | [トラブルシューティング](troubleshooting.md) |
| **CLI コマンドリファレンス（英語のみ、自動生成）** | [コマンドリファレンス（英語のみ）](../reference/commands/index.md) |
| **YAML 設定オプション** | [YAML リファレンス](yaml-reference.md) |
| **MLIP バックエンド設定** | [MLIP バックエンド](backends.md) |
| **各コマンドが書き出すファイル** | [出力ディレクトリのレイアウト](output-layout.md) |
| **`result.json` と `summary.json` の欄** | [JSON 出力リファレンス](json-output.md) |
| **複数の GPU ノードで動かす（PBS + Ray）** | [HPC 実行例](hpc-example.md) |
| **AI エージェントから呼ぶ（MCP）** | [MCP サーバー](mcp_server.md) |
| **コードの構成（開発者向け）** | [アーキテクチャ](architecture.md) |
| **用語** | [用語集](glossary.md) |

## システム要件

### ハードウェア
- **OS**: Linux（Windows では WSL2 上の Linux に導入してください）
- **GPU**: 使用するバックエンドと PyTorch wheel に対応する NVIDIA ドライバー。CPU のみでも実行可能ですが低速です
- **VRAM / RAM**: モデル、系の大きさ、Hessian の計算方式で変わります。代表的な計算で最大使用量を測ってください

### ソフトウェア
- Python >= 3.11
- CPU 版または CUDA 対応の PyTorch。ビルド済みの wheel は CUDA ランタイムを含むので、手元の CUDA toolkit は通常いりません（ソースからビルドするときだけ必要です）

セットアップは [インストール](installation.md) を参照してください。

## エージェントスキル

`pdb2reaction` は、CLI サブコマンド・構造 I/O・バックエンドインストール・ワークフロー・出力解析・HPC 運用をカバーする AI エージェント向けの手順書を `skills/` に同梱しています。導入するときは、AI エージェントに次のように指示してください。

> `https://github.com/t-0hmura/pdb2reaction/tree/main/skills` をスキルとして取り込み、`pdb2reaction-install-backends` の手順に従って pdb2reaction をインストールして

GitHub のリポジトリを clone 済みなら、URL の代わりに手元の `skills/` の path を渡しても構いません。導入した後は、たとえば次のように頼めます。

> 〈論文〉を読んで、〈PDB ID〉の構造からモデルを作成し、〈反応段階〉の経路について、pdb2reaction のスキルを用いて反応機構解析を行ってください。

## 引用

`pdb2reaction` を研究で利用する場合は、ACS Omega の論文を引用してください:

```bibtex
@article{ohmura2026pdb2reaction,
  author       = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title        = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  journal      = {ACS Omega},
  year         = {2026},
  doi          = {10.1021/acsomega.6c08242}
}
```

ソフトウェアまたは特定のリリースを引用する場合は、Zenodo レコードを使用してください:

```bibtex
@software{ohmura2026pdb2reaction_software,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  version      = {0.5.0},
  url          = {https://github.com/t-0hmura/pdb2reaction},
  license      = {GPL-3.0},
  doi          = {10.5281/zenodo.19197865}
}
```

## ライセンス

`pdb2reaction` は **GNU General Public License version 3 (GPL-3.0)** の下で配布されています。

## ヘルプ

```bash
# 一般的なヘルプ
pdb2reaction --help

# コマンドのヘルプ
pdb2reaction <subcommand> --help

# 詳細オプション（内部チューニング用）
pdb2reaction <subcommand> --help-advanced
```

問題や機能リクエストについては、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) を参照してください。
