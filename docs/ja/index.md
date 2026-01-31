# pdb2reaction ドキュメント

**pdb2reaction** は、機械学習原子間ポテンシャル (MLIP) を使用して、PDB構造から酵素反応経路を自動モデリングするPython CLIツールキットです。

```{toctree}
:maxdepth: 2
:caption: 目次
:hidden:

getting-started
concepts
troubleshooting
all
extract
add_elem_info
opt
tsopt
path_opt
path_search
scan
scan2d
scan3d
irc
freq
dft
trj2fig
yaml-reference
uma_pysis
glossary
```

---

## 目的別クイックスタート

| やりたいこと | 推奨コマンド | ガイド |
|--------------|--------------|--------|
| PDBから完全な反応経路探索を実行 | `pdb2reaction all` | [all.md](all.md) |
| タンパク質-リガンド複合体からQM領域を抽出 | `pdb2reaction extract` | [extract.md](extract.md) |
| 単一構造を最適化 | `pdb2reaction opt` | [opt.md](opt.md) |
| 遷移状態を探索・精緻化 | `pdb2reaction tsopt` | [tsopt.md](tsopt.md) |
| 最小エネルギー経路を探索 | `pdb2reaction path-search` | [path_search.md](path_search.md) |
| 遷移状態からIRCを実行 | `pdb2reaction irc` | [irc.md](irc.md) |
| エネルギープロファイルを可視化 | `pdb2reaction trj2fig` | [trj2fig.md](trj2fig.md) |
| 全体像（概念・用語）を把握したい | — | [概念とワークフロー](concepts.md) |
| よくあるエラーを解決したい | — | [トラブルシューティング](troubleshooting.md) |
| 略語や用語を調べる | — | [用語集](glossary.md) |

---

## クイックナビゲーション

### はじめに

- [**はじめに**](getting-started.md) - インストール、クイックスタート、概要
- [**概念とワークフロー**](concepts.md) - ポケット、テンプレート、セグメント、各ステージの全体像
- [**トラブルシューティング**](troubleshooting.md) - よくあるエラーと対処法
- [**システム要件**](#システム要件) - ハードウェアとソフトウェアの前提条件

### メインワークフロー

- [`all`](all.md) - **エンドツーエンドワークフロー**: 抽出 → スキャン → MEP探索 → TS最適化 → IRC → 熱化学 → DFT

### CLIサブコマンド

#### 構造準備
| サブコマンド | 説明 |
|---------|------|
| [`extract`](extract.md) | タンパク質-リガンド複合体から活性部位クラスターを抽出 |
| [`add-elem-info`](add_elem_info.md) | PDB元素カラム（77-78）を修復 |

#### 構造最適化
| サブコマンド | 説明 |
|---------|------|
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS / RFO） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer / RS-I-RFO） |

#### 経路探索・最適化
| サブコマンド | 説明 |
|---------|------|
| [`path-opt`](path_opt.md) | GSMまたはDMFによるMEP最適化 |
| [`path-search`](path_search.md) | 自動精密化を伴う再帰的MEP探索 |

#### スキャン
| サブコマンド | 説明 |
|---------|------|
| [`scan`](scan.md) | 拘束条件付き1D結合長スキャン |
| [`scan2d`](scan2d.md) | 2D距離グリッドスキャン |
| [`scan3d`](scan3d.md) | 3D距離グリッドスキャン |

#### 解析・後処理
| サブコマンド | 説明 |
|---------|------|
| [`irc`](irc.md) | 固有反応座標（IRC）計算 |
| [`freq`](freq.md) | 振動数解析と熱化学 |
| [`dft`](dft.md) | 一点DFT計算（GPU4PySCF / PySCF） |
| [`trj2fig`](trj2fig.md) | XYZ軌跡からエネルギープロファイルをプロット |

### 設定・リファレンス

- [**YAMLリファレンス**](yaml-reference.md) - 全サブコマンドのYAML設定オプション
- [**UMA計算機**](uma_pysis.md) - UMA機械学習ポテンシャル設定
- [**用語集**](glossary.md) - 略語と技術用語の定義

---

## システム要件

### ハードウェア
- **OS**: Linux（Ubuntu 20.04+、CentOS 8+で動作確認）
- **GPU**: CUDA 12.x互換（RTX 30xx/40xx、A100、H100で動作確認）
- **VRAM**: 最小8 GB（1000原子以上には16 GB以上推奨）
- **RAM**: 16 GB以上推奨

### ソフトウェア
- Python 3.11
- CUDAサポート付きPyTorch
- CUDA 12.xツールキット

---

## クイック例

### 基本的なMEP探索
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### TS精密化を含む完全ワークフロー
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True --thermo True --dft True
```

### 単一構造スキャンモード
```bash
pdb2reaction -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### TS最適化のみ
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True
```

---

## 重要な概念

### 電荷とスピン
- 未知残基の電荷を指定するには `--ligand-charge` を使用: `'SAM:1,GPP:-3'`
- 総電荷を上書きするには `-q/--charge` を使用
- スピン多重度は `-m/--mult`（`all` コマンド）または `-m/--multiplicity`（他のサブコマンド）で設定（デフォルト: 1）

### ブール値オプション
すべてのブール値CLIオプションは明示的に `True` または `False` を指定する必要があります:
```bash
--tsopt True --thermo True --dft False
```

### YAML設定
高度な設定は `--args-yaml` で指定できます:
```bash
pdb2reaction all -i R.pdb P.pdb -c 'LIG' --args-yaml config.yaml
```
すべてのオプションについては [YAMLリファレンス](yaml-reference.md) を参照してください。

---

## 出力構造

典型的な `pdb2reaction all` 実行の出力:
```
result_all/
├── summary.log              # 人間が読めるサマリー
├── summary.yaml             # 機械可読サマリー
├── pockets/                 # 抽出されたクラスターモデル
├── scan/                    # （オプション）スキャン結果
├── path_search/             # MEP軌跡とダイアグラム
│   ├── mep.trj              # MEP軌跡
│   ├── mep.pdb              # PDB形式のMEP
│   ├── mep_w_ref.pdb        # 全系とマージされたMEP
│   ├── mep_plot.png         # エネルギープロファイルプロット
│   └── seg_*/               # セグメントごとの詳細
└── path_search/post_seg_*/  # 後処理出力
    ├── tsopt/               # TS最適化結果
    ├── irc/                 # IRC軌跡
    ├── freq/                # 振動モード
    └── dft/                 # DFT結果
```

---

## 引用

`pdb2reaction` を説明するプレプリントを準備中です。引用情報については後日ご確認ください。

## ライセンス

`pdb2reaction` は Pysisyphus から派生した **GNU General Public License version 3 (GPL-3.0)** の下で配布されています。

---

## 参考文献

1. Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. [arXiv:2506.23971](http://arxiv.org/abs/2506.23971)
2. Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). [DOI:10.1002/qua.26390](https://doi.org/10.1002/qua.26390)

---

## ヘルプ

```bash
# 一般的なヘルプ
pdb2reaction --help

# コマンドのヘルプ
pdb2reaction <command> --help
```

問題や機能リクエストについては、[GitHubリポジトリ](https://github.com/t-0hmura/pdb2reaction) をご覧ください。
