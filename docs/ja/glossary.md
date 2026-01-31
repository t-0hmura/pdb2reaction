# 用語集

このページでは、pdb2reaction ドキュメント内で使われる略語・専門用語を簡潔に説明します。

---

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path | 反応物から生成物へ至る最小エネルギー経路（ポテンシャルエネルギー面上の最も低い経路）。 |
| **TS** | Transition State | 反応座標に沿ったエネルギー極大に対応する一次の鞍点。 |
| **IRC** | Intrinsic Reaction Coordinate | TS から反応物側・生成物側へ向かう、質量重み付き最急降下経路。TS の接続検証によく使われます。 |
| **GSM** | Growing String Method | 端点からストリング（画像列）を伸長・最適化して MEP を近似する手法。 |
| **DMF** | Direct Max Flux | 反応座標方向のフラックスを最大化することで MEP を最適化する chain-of-states 手法。pdb2reaction では `--mep-mode dmf` で選択します。 |
| **NEB** | Nudged Elastic Band | 画像間にばね力を導入し、画像間隔を保ちながら経路を最適化する chain-of-states 手法。 |

---

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する（多くはニューラルネットの）原子間ポテンシャル。 |
| **UMA** | Universal Machine-learning potential for Atoms | Meta が公開している事前学習 MLIP 群。pdb2reaction のデフォルト計算機バックエンドです。 |

---

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算。 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法。 |
| **Hessian** | — | エネルギーの二階微分行列。振動解析や TS 最適化に使用します。 |

---

## 構造生物学・入出力形式

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース）。 |
| **XYZ** | — | 元素記号と直交座標を並べたシンプルなテキスト形式。 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。pdb2reaction は電荷/多重度と座標の読み取りに利用します。 |

---

## 単位・定数

| 用語 | 説明 |
|------|------|
| **Hartree** | 原子単位系のエネルギー。1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV。 |
| **kcal/mol** | 反応エネルギー表現でよく使われる単位。 |
| **Bohr** | 原子単位系の長さ。1 Bohr ≈ 0.529 Å。 |
| **Å（オングストローム）** | 10⁻¹⁰ m。原子間距離の表現でよく使われる長さ単位。 |

---

## 関連ページ

- [はじめに](getting-started.md) — インストールと初回実行
- [概念とワークフロー](concepts.md) — ポケット抽出、MEP 探索、後処理の全体像
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAMLリファレンス](yaml-reference.md) — 設定ファイルの仕様
- [UMA計算機](uma_pysis.md) — UMA ポテンシャルの詳細
