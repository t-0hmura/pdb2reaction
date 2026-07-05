# pdb2reaction ドキュメント

*バージョン: v0.4.0*

---

<img src="../overview.png" alt="pdb2reaction ワークフロー概要" width="90%">

**pdb2reaction** は、機械学習原子間ポテンシャル（MLIP: Machine Learning Interatomic Potential）を使用して、PDB 構造から酵素反応経路を自動的に解明する Python 製 CLI ツールキットです。

## クイックスタート

| 目的 | ワークフロー |
|------|------------|
| **End-to-end 初回実行** | [クイックスタート: all](quickstart-all.md) |
| **反応物のみ** | [クイックスタート: scan](quickstart-scan.md) |
| **TS 候補あり** | [クイックスタート: TS のみモード](quickstart-tsopt-freq.md) |
| **実行失敗 / エラー** | [典型エラー別レシピ](recipes-common-errors.md) |

前提条件は [Installation](installation.md) を参照してください。

## サブコマンド

| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | end-to-end ワークフロー: 抽出 → スキャン → MEP 探索 → TS 最適化 → IRC → 熱化学 → DFT |
| [`extract`](extract.md) | タンパク質–リガンド複合体から活性部位モデル（バインディングポケット）を抽出 |
| [`fix-altloc`](fix-altloc.md) | PDB の代替位置指示子を解決 |
| [`add-elem-info`](add-elem-info.md) | PDB の元素カラム（77–78）を修復 |
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS または RFO。[+ 任意 flatten]） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer または RS-I-RFO。[+ 任意 flatten]） |
| [`path-opt`](path-opt.md) | GSM または DMF による 1 段階の MEP 最適化（2 構造から） |
| [`path-search`](path-search.md) | 自動精密化を伴う多段階の再帰的 MEP 探索（2 構造以上） |
| [`scan`](scan.md) | 拘束付き 1D 結合長スキャン |
| [`scan2d`](scan2d.md) | 2D 距離グリッドスキャン |
| [`scan3d`](scan3d.md) | 3D 距離グリッドスキャン |
| [`freq`](freq.md) | 振動解析と熱化学 |
| [`irc`](irc.md) | 固有反応座標（IRC: Intrinsic Reaction Coordinate）計算 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`sp`](sp.md) | MLIP による一点計算（エネルギー + 力 / Hessian） |
| [`trj2fig`](trj2fig.md) | XYZ 軌跡からエネルギープロファイルをプロット |
| [`energy-diagram`](energy-diagram.md) | 数値入力からエネルギーダイアグラムを作成 |
| [`bond-summary`](bond-summary.md) | 連続構造間の共有結合変化を検出・レポート |

## 設定・リファレンス

| トピック | ページ |
|-------|------|
| **CLI 規約と入力要件** | [CLI 規約](cli-conventions.md) |
| **クラスター境界の凍結原子（キャップ水素・ `--freeze-atoms`）** | [凍結原子](freeze-atoms.md) |
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

`pdb2reaction` は、CLI サブコマンド・構造 I/O・バックエンドインストール・ワークフロー・出力解析・HPC 運用をカバーする AI エージェント向けの手順書を `skills/` に同梱しています。完全なスキル索引とインストール手順は [`skills/README.md`](https://github.com/t-0hmura/pdb2reaction/blob/main/skills/README.md) を参照してください。

## 引用

`pdb2reaction` を研究で利用する場合は、ChemRxiv プレプリントを引用してください:

```bibtex
@misc{ohmura2026pdb2reaction,
  author       = {Ohmura, Takuto and Sato, Hajime and Terada, Tohru},
  title        = {pdb2reaction: End-to-End Reaction-Path Elucidation from PDB Structures Using Machine-Learning Interatomic Potentials},
  year         = {2026},
  doi          = {10.26434/chemrxiv.15003538},
  note         = {ChemRxiv preprint}
}
```

ソフトウェアまたは特定のリリースを引用する場合は、Zenodo レコードを使用してください:

```bibtex
@software{ohmura2026pdb2reaction_software,
  author       = {Ohmura, Takuto},
  title        = {pdb2reaction},
  year         = {2026},
  month        = {6},
  version      = {0.4.0},
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

問題や機能リクエストについては、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) を参照してください。
