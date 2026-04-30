# pdb2reaction ドキュメント

*バージョン: v0.3.8*

---

<img src="../overview.png" alt="pdb2reaction ワークフロー概要" width="90%">

**pdb2reaction** は、機械学習原子間ポテンシャル（MLIP: Machine Learning Interatomic Potential）を使用して、PDB 構造から酵素反応経路を自動モデリングする Python 製 CLI ツールキットです。

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

getting-started
installation
quickstart-all
quickstart-scan
quickstart-tsopt-freq
freeze-atoms
recipes-common-errors
troubleshooting
cli-conventions
```

```{toctree}
:maxdepth: 2
:caption: コマンド
:hidden:

all
extract
fix-altloc
add-elem-info
opt
tsopt
path-opt
path-search
scan
scan2d
scan3d
freq
irc
dft
trj2fig
energy-diagram
bond-summary
```

```{toctree}
:maxdepth: 2
:caption: リファレンス
:hidden:

../reference/commands/index
yaml-reference
json-output
uma-pysis
hpc-example
glossary
```

## ここから始める

| 目的 | ワークフロー |
|------|------------|
| **End-to-end 初回実行** | [クイックスタート: all](quickstart-all.md) |
| **反応物のみ** | [クイックスタート: scan](quickstart-scan.md) |
| **TS 候補あり** | [クイックスタート: tsopt + freq](quickstart-tsopt-freq.md) |
| **実行失敗 / エラー** | [典型エラー別レシピ](recipes-common-errors.md) |

前提条件は [Installation](installation.md) を参照してください。

## サブコマンド

| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | end-to-end ワークフロー: 抽出 → スキャン → MEP 探索 → TS 最適化 → IRC → 熱化学 → DFT |
| [`extract`](extract.md) | タンパク質–リガンド複合体から活性部位モデル（バインディングポケット）を抽出 |
| [`fix-altloc`](fix-altloc.md) | PDB の代替位置指示子を解決 |
| [`add-elem-info`](add-elem-info.md) | PDB の元素カラム（77–78）を修復 |
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS or RFO。[+ 任意flatten]） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer or RS-I-RFO。[+ 任意flatten]） |
| [`path-opt`](path-opt.md) | GSM または DMF による 1段階の MEP 最適化（2構造から） |
| [`path-search`](path-search.md) | 自動精密化を伴う多段階の再帰的 MEP 探索（2構造以上） |
| [`scan`](scan.md) | 拘束付き 1D 結合長スキャン |
| [`scan2d`](scan2d.md) | 2D 距離グリッドスキャン |
| [`scan3d`](scan3d.md) | 3D 距離グリッドスキャン |
| [`freq`](freq.md) | 振動解析と熱化学 |
| [`irc`](irc.md) | 固有反応座標（IRC: Intrinsic Reaction Coordinate）計算 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`trj2fig`](trj2fig.md) | XYZ軌跡からエネルギープロファイルをプロット |
| [`energy-diagram`](energy-diagram.md) | 数値入力からエネルギーダイアグラムを作成 |
| [`bond-summary`](bond-summary.md) | 連続構造間の共有結合変化を検出・レポート |

## 設定・リファレンス

| トピック | ページ |
|-------|------|
| **CLI 規約と入力要件** | [CLI 規約](cli-conventions.md) |
| **クラスター境界の凍結原子（リンク水素・`--freeze-atoms`）** | [凍結原子](freeze-atoms.md) |
| **よくあるエラーと対処** | [トラブルシューティング](troubleshooting.md) |
| **CLI コマンドリファレンス（英語のみ、自動生成）** | [コマンドリファレンス（英語のみ）](../reference/commands/index.md) |
| **YAML 設定オプション** | [YAML リファレンス](yaml-reference.md) |
| **MLIP バックエンド設定** | [MLIP 計算機](uma-pysis.md) |
| **用語** | [用語集](glossary.md) |

## システム要件

### ハードウェア
- **OS**: Linux
- **GPU**: CUDA 12.x 互換
- **VRAM**: 8 GB 以上推奨
- **RAM**: 16 GB 以上推奨

### ソフトウェア
- Python >= 3.11
- CUDA サポート付き PyTorch
- CUDA 12.x ツールキット

セットアップは [インストール](installation.md) を参照してください。

## エージェントスキル

`pdb2reaction` は、CLI サブコマンド・構造 I/O・バックエンドインストール・ワークフロー・出力解析・HPC 運用をカバーする AI エージェント向け命令を `.claude/skills/` に同梱しています。完全なスキル索引とインストール手順は [`.claude/skills/README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/.claude/skills/README.md) を参照してください。

## 引用

`pdb2reaction` を説明するプレプリントを準備中です。現時点では、本ソフトウェア自体を引用してください:

```bibtex
@software{ohmura2026pdb2reaction,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  month        = {4},
  version      = {0.3.8},
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

# 詳細オプション（dry-run、内部チューニング等）
pdb2reaction <subcommand> --help-advanced
```

問題や機能リクエストについては、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) をご覧ください。
