# pdb2reaction ドキュメント

*バージョン: v{{ version }}*

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

yaml-reference
json-output
uma-pysis
hpc-example
glossary
```

---

## ドキュメント案内

| トピック | ページ |
|-------|------|
| **はじめに** | [はじめに](getting-started.md) |
| **インストール** | [インストール](installation.md) |
| **症状起点の切り分け導線** | [典型エラー別レシピ](recipes-common-errors.md) |
| **よくあるエラーと対処** | [トラブルシューティング](troubleshooting.md) |
| **CLI 規約と入力要件** | [CLI 規約](cli-conventions.md) |

---

## 目的別ガイド

| 目的 | 推奨コマンド | ガイド |
|--------------|--------------|--------|
| 最初の 1 回を実行（end-to-end） | `pdb2reaction all` | [クイックスタート: all](quickstart-all.md) |
| 単一構造の段階的スキャン | `pdb2reaction scan` | [クイックスタート: scan](quickstart-scan.md) |
| TS 最適化と検証 | `pdb2reaction tsopt` | [クイックスタート: tsopt](quickstart-tsopt-freq.md) |
| PDB から反応経路探索を一通り実行 | `pdb2reaction all` | [all](all.md) |
| タンパク質–リガンド複合体から活性部位をクラスターモデルとして抽出 | `pdb2reaction extract` | [extract](extract.md) |
| 単一構造を最適化 | `pdb2reaction opt` | [opt](opt.md) |
| 遷移状態を探索・最適化 | `pdb2reaction tsopt` | [tsopt](tsopt.md) |
| 最小エネルギー経路を探索し TS 候補を取得 | `pdb2reaction path-search` | [path-search](path-search.md) |
| 遷移状態からIRCを実行 | `pdb2reaction irc` | [irc](irc.md) |
| エネルギープロファイルを可視化 | `pdb2reaction trj2fig` | [trj2fig](trj2fig.md) |
| 数値から状態エネルギーダイアグラムを描画 | `pdb2reaction energy-diagram` | [energy-diagram](energy-diagram.md) |
| 症状からエラー対処を探す | — | [典型エラー別レシピ](recipes-common-errors.md) |
| よくあるエラーを解決したい | — | [トラブルシューティング](troubleshooting.md) |
| 略語や用語を調べる | — | [用語集](glossary.md) |
| CLI の規約（フラグ、優先順位、原子/残基セレクタ）を調べる | — | [CLI 規約](cli-conventions.md) |

---

## CLI サブコマンド

### メインワークフロー
| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | end-to-end ワークフロー: 抽出 → スキャン → MEP 探索 → TS 最適化 → IRC → 熱化学 → DFT |

### 構造準備
| サブコマンド | 説明 |
|---------|------|
| [`extract`](extract.md) | タンパク質–リガンド複合体から活性部位モデル（バインディングポケット）を抽出 |
| [`fix-altloc`](fix-altloc.md) | PDB の代替位置指示子を解決 |
| [`add-elem-info`](add-elem-info.md) | PDB の元素カラム（77–78）を修復 |

### 構造最適化
| サブコマンド | 説明 |
|---------|------|
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS or RFO。[+ 任意flatten]） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer or RS-I-RFO。[+ 任意flatten]） |

### 経路探索・最適化
| サブコマンド | 説明 |
|---------|------|
| [`path-opt`](path-opt.md) | GSM または DMF による 1段階の MEP 最適化（2構造から） |
| [`path-search`](path-search.md) | 自動精密化を伴う多段階の再帰的 MEP 探索（2構造以上） |

### スキャン
| サブコマンド | 説明 |
|---------|------|
| [`scan`](scan.md) | 拘束付き 1D 結合長スキャン |
| [`scan2d`](scan2d.md) | 2D 距離グリッドスキャン |
| [`scan3d`](scan3d.md) | 3D 距離グリッドスキャン |

### 解析・後処理
| サブコマンド | 説明 |
|---------|------|
| [`irc`](irc.md) | 固有反応座標（IRC: Intrinsic Reaction Coordinate）計算 |
| [`freq`](freq.md) | 振動解析と熱化学 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`trj2fig`](trj2fig.md) | XYZ軌跡からエネルギープロファイルをプロット |
| [`energy-diagram`](energy-diagram.md) | 数値入力からエネルギーダイアグラムを作成 |
| [`bond-summary`](bond-summary.md) | 連続構造間の共有結合変化を検出・レポート |

---

## 設定・リファレンス

| トピック | ページ |
|-------|------|
| **CLI コマンドリファレンス（英語のみ、自動生成）** | [コマンドリファレンス（英語のみ）](../reference/commands/index.md) |
| **YAML 設定オプション** | [YAML リファレンス](yaml-reference.md) |
| **MLIP バックエンド設定** | [MLIP 計算機](uma-pysis.md) |
| **用語** | [用語集](glossary.md) |

---

## システム要件

### ハードウェア
- **OS**: Linux
- **GPU**: CUDA 12.x 互換
- **VRAM**: 8 GB 以上推奨
- **RAM**: 16 GB以上推奨

### ソフトウェア
- Python >= 3.11
- CUDA サポート付き PyTorch
- CUDA 12.x ツールキット

---

## 実行例

### 基本的な MEP 探索
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

### MEP 探索からの TS 最適化を含むワークフロー
```bash
pdb2reaction -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft
```

### 単一構造からの反応座標スキャンによるワークフロー
```bash
pdb2reaction -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("CS1 SAM 320","GPP 321 C7",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
 --tsopt --thermo --dft
```

### 単一 TS 候補構造からのワークフロー
```bash
pdb2reaction -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft
```

---

## 入力の主要フラグ

### 電荷とスピン
- 基質残基の電荷を指定するには `-l/--ligand-charge` を使用: `'SAM:1,GPP:-3'`
- 総電荷を上書きするには `-q/--charge` を使用
- スピン多重度（Spin Multiplicity）は `-m/--multiplicity`（デフォルト: 1 = 一重項）で設定

### ブール値オプション
ブール値 CLI オプションは `--flag` / `--no-flag` と `--flag True/False` の両方を受理します（新規スクリプトでは `--flag` / `--no-flag` 形式を推奨）:
```bash
--tsopt --thermo --no-dft
--tsopt True --thermo yes --dft 0
```

### YAML 設定
すべてのオプションについては [YAML リファレンス](yaml-reference.md) を参照してください。

---

## 出力構造

典型的な `pdb2reaction all` の出力（デフォルトは `--refine-path True` で `path_search/` を使用。`--refine-path False` では代わりに `path_opt/` が使われます — このフラグは `pdb2reaction all` にのみ属します、定義は {ref}`mep-search-options` を参照）:

```
result_all/
├── summary.log               # 結果要約
├── summary.json              # JSON 結果
├── models/                   # 抽出されたクラスターモデル
├── scan/                     # （オプション）スキャン結果
└── path_search/              # MEP 探索結果（デフォルト）
    ├── mep_trj.xyz           # 最終 MEP 軌跡
    ├── mep.pdb               # PDB 形式の MEP
    ├── mep_plot.png          # エネルギープロファイル図
    ├── seg_*/                # セグメントごとの MEP 詳細 (入力と MEP 出力)
    │   └── structures/       # 抽出された定常点構造
    └── post_seg_*/           # セグメントごとの後処理
        ├── ts/               # TS 最適化結果
        ├── irc/              # IRC 軌跡
        ├── freq/             # 振動解析
        └── dft/              # DFT 結果
```

`--refine-path False` を指定した場合、`path_search/` サブツリーは同じレイアウトの `path_opt/` に置き換わります。

---

## 引用

`pdb2reaction` を説明するプレプリントを準備中です。引用情報については後日ご確認ください。

## ライセンス

`pdb2reaction` は **GNU General Public License version 3 (GPL-3.0)** の下で配布されています。

---

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

---

*Note: 本ドキュメントは整備中です。一部未完成の箇所や今後変更される箇所があります。*
