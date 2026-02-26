# pdb2reaction ドキュメント

**pdb2reaction** は、機械学習原子間ポテンシャル (MLIP) を使用して、PDB 構造から酵素反応経路を自動モデリングする Python 製 CLI ツールキットです。

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting_started
quickstart_all
quickstart_scan
quickstart_tsopt_freq
concepts
recipes_common_errors
troubleshooting
cli_conventions
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
init
extract
fix_altloc
add_elem_info
opt
tsopt
path_opt
path_search
scan
scan2d
scan3d
freq
irc
dft
trj2fig
energy_diagram
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

yaml_reference
uma_pysis
glossary
```


---

## 目的別ガイド

| 目的 | 推奨コマンド | ガイド |
|--------------|--------------|--------|
| 最初の 1 回を実行（end-to-end） | `pdb2reaction all` | [クイックスタート: all](quickstart_all.md) |
| 単一構造の段階的スキャン | `pdb2reaction scan` | [クイックスタート: scan](quickstart_scan.md) |
| TS 検証（`tsopt` -> `freq`） | `pdb2reaction tsopt`, `pdb2reaction freq` | [クイックスタート: tsopt -> freq](quickstart_tsopt_freq.md) |
| PDB から反応経路探索を一通り実行 | `pdb2reaction all` | [all.md](all.md) |
| タンパク質-リガンド複合体からQM領域を抽出 | `pdb2reaction extract` | [extract.md](extract.md) |
| 単一構造を最適化 | `pdb2reaction opt` | [opt.md](opt.md) |
| 遷移状態を探索・最適化 | `pdb2reaction tsopt` | [tsopt.md](tsopt.md) |
| 最小エネルギー経路を探索 | `pdb2reaction path-search` | [path_search.md](path_search.md) |
| 遷移状態からIRCを実行 | `pdb2reaction irc` | [irc.md](irc.md) |
| エネルギープロファイルを可視化 | `pdb2reaction trj2fig` | [trj2fig.md](trj2fig.md) |
| 数値から状態エネルギーダイアグラムを描画 | `pdb2reaction energy-diagram` | [energy_diagram.md](energy_diagram.md) |
| 症状からエラー対処を探す | — | [典型エラー別レシピ](recipes_common_errors.md) |
| 全体像（概念・用語）を把握したい | — | [概念とワークフロー](concepts.md) |
| よくあるエラーを解決したい | — | [トラブルシューティング](troubleshooting.md) |
| 略語や用語を調べる | — | [用語集](glossary.md) |

---

## ドキュメント案内

### はじめに

- [**はじめに**](getting_started.md) - インストール、クイックスタート、概要
- [**概念とワークフロー**](concepts.md) - ポケット、テンプレート、セグメント、各ステージの全体像
- [**典型エラー別レシピ**](recipes_common_errors.md) - 症状別の最短対処ルート
- [**CLI 規約**](cli_conventions.md) - ブール値オプション、セレクタ、電荷指定などの共通規約
- [**トラブルシューティング**](troubleshooting.md) - よくあるエラーと対処法
- [**システム要件**](#システム要件) - 必要なハードウェア・ソフトウェア

### メインワークフロー

- [`all`](all.md) - **end-to-endワークフロー**: 抽出 → スキャン → MEP 探索 → TS 最適化 → IRC → 熱化学 → DFT

### CLI サブコマンド

#### 構造準備
| サブコマンド | 説明 |
|---------|------|
| [`extract`](extract.md) | タンパク質-リガンド複合体から活性部位ポケット（クラスターモデル）を抽出 |
| [`add-elem-info`](add_elem_info.md) | PDB の元素カラム（77–78）を修復 |

#### 構造最適化
| サブコマンド | 説明 |
|---------|------|
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS / RFO / hybrid + 任意flatten） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer / RS-I-RFO / hybrid、flattenは任意） |

#### 経路探索・最適化
| サブコマンド | 説明 |
|---------|------|
| [`path-opt`](path_opt.md) | GSM または DMF による 1段階の MEP 最適化 |
| [`path-search`](path_search.md) | 自動精密化を伴う多段階の再帰的 MEP 探索 |

#### スキャン
| サブコマンド | 説明 |
|---------|------|
| [`scan`](scan.md) | 拘束付き 1D 結合長スキャン |
| [`scan2d`](scan2d.md) | 2D 距離グリッドスキャン |
| [`scan3d`](scan3d.md) | 3D 距離グリッドスキャン |

#### 解析・後処理
| サブコマンド | 説明 |
|---------|------|
| [`irc`](irc.md) | 固有反応座標（IRC）計算 |
| [`freq`](freq.md) | 振動数解析と熱化学 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`trj2fig`](trj2fig.md) | XYZ軌跡からエネルギープロファイルをプロット |
| [`energy-diagram`](energy_diagram.md) | 数値入力から状態エネルギーダイアグラムを作成 |

### 設定・リファレンス

- [**CLI コマンドリファレンス（自動生成）**](../reference/commands/index.md)
- [**YAML スキーマ（自動生成）**](../reference/yaml.md)
- [**YAML リファレンス**](yaml_reference.md) - 全サブコマンドの YAML 設定オプション
- [**UMA 計算機**](uma_pysis.md) - UMA 機械学習ポテンシャル設定
- [**用語集**](glossary.md) - 略語と技術用語の定義

---

## システム要件

### ハードウェア
- **OS**: Linux（Ubuntu 20.04+、CentOS 8+で動作確認）
- **GPU**: CUDA 12.x 互換 - **VRAM**: 最小8 GB（1000原子以上には16 GB以上推奨）
- **RAM**: 16 GB以上推奨

### ソフトウェア
- Python 3.11
- CUDA サポート付き PyTorch
- CUDA 12.x ツールキット

---

## 実行例

### 基本的な MEP 探索
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### TS 最適化を含む完全ワークフロー
```bash
pdb2reaction -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft
```

### 単一構造スキャンモード
```bash
pdb2reaction scan -i input.pdb -q 0 -m 1 --spec scan.yaml
```

### TS 最適化のみ
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --tsopt
```

---

## 重要な概念

### 電荷とスピン
- 基質残基の電荷を指定するには `--ligand-charge` を使用: `'SAM:1,GPP:-3'`
- 総電荷を上書きするには `-q/--charge` を使用
- スピン多重度は `-m/--multiplicity`（デフォルト: 1）で設定

### ブール値オプション
ブール値 CLI オプションは `--flag` / `--no-flag` と `--flag True/False` の両方を受理します（新規スクリプトでは toggle 形式を推奨）:
```bash
--tsopt --thermo --no-dft
--tsopt True --thermo yes --dft 0
```

### YAML 設定
すべてのオプションについては [YAML リファレンス](yaml_reference.md) を参照してください。

---

## 出力構造

典型的な `pdb2reaction all` の出力:
```
result_all/
├── summary.log # テキスト形式の結果要約
├── summary.yaml # YAML 形式の結果要約
├── pockets/ # 抽出されたクラスターモデル
├── scan/ # （オプション）スキャン結果
├── path_search/ # MEP軌跡とダイアグラム
│ ├── mep_trj.xyz # MEP軌跡
│ ├── mep.pdb # PDB形式のMEP
│ ├── mep_w_ref.pdb # 全系とマージされたMEP
│ ├── mep_plot.png # エネルギープロファイルプロット
│ └── seg_*/ # セグメントごとの詳細
└── path_search/post_seg_*/ # 後処理出力
 ├── tsopt/ # TS最適化結果
 ├── irc/ # IRC軌跡
 ├── freq/ # 振動モード
 └── dft/ # DFT結果
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
