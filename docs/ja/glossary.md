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
| **HEI** | Highest-Energy Image | MEP 上でエネルギーが最大の画像。TS の初期推定としてよく使われます。 |
| **画像（Image）** | — | 経路上の1つの幾何（1ノード）。 |
| **セグメント** | — | 2つの隣接する端点を結ぶ MEP（例: R → I1, I1 → I2, …）。 |

---

## 最適化アルゴリズム

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **L-BFGS** | Limited-memory BFGS | 勾配履歴からヘシアンを近似する準ニュートン法。`--opt-mode light` で使用。 |
| **RFO** | Rational Function Optimization | 明示的なヘシアン情報を使用する信頼領域最適化法。`--opt-mode heavy` で使用。 |
| **RS-I-RFO** | Restricted-Step Image-RFO | 1つの負固有値方向に沿う、鞍点（TS）最適化用の RFO 変種。 |
| **Dimer** | Dimer Method | 完全なヘシアンを計算せずに最低曲率モードを推定する TS 最適化法。`--opt-mode light` の TSOPT で使用。 |

---

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する（多くはニューラルネットの）原子間ポテンシャル。 |
| **UMA** | Universal Machine-learning potential for Atoms | Meta が公開している事前学習 MLIP 群。pdb2reaction のデフォルト計算機バックエンドです。 |
| **解析ヘシアン** | Analytical Hessian | エネルギーの正確な二階微分を計算。高速だが VRAM を多く消費。 |
| **有限差分** | Finite Difference | 微小変位による微分近似。低速だがメモリ効率が良い。 |

---

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算。 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法。 |
| **Hessian** | — | エネルギーの二階微分行列。振動解析や TS 最適化に使用します。 |
| **SP** | Single Point | 固定構造での計算（最適化なし）。高精度エネルギー補正によく使用。 |
| **スピン多重度** | Spin Multiplicity | 2S+1（S は全スピン）。シングレット = 1、ダブレット = 2、トリプレット = 3 など。 |

---

## 構造生物学・ポケット抽出

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース）。 |
| **XYZ** | — | 元素記号と直交座標を並べたシンプルなテキスト形式。 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。pdb2reaction は電荷/多重度と座標の読み取りに利用します。 |
| **ポケット** | Active-site Pocket | 基質周辺を切り出した部分構造。MEP/TS 探索の計算量削減に使用。「クラスターモデル」とも。 |
| **クラスターモデル** | Cluster Model | ポケットの別名。酵素–基質複合体から計算可能なサイズに切り出した部分系。 |
| **リンク水素** | Link Hydrogen | ポケット抽出時に切断された結合をキャップするために付加する水素原子。 |
| **主鎖** | Backbone | タンパク質の主骨格（N–Cα–C–O 原子）。`--exclude-backbone True` で除外可能。 |

---

## 熱化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ZPE** | Zero-Point Energy | 0 K での振動エネルギー。電子エネルギーへの量子補正。 |
| **ギブズエネルギー** | Free Energy (G) | G = H − TS。熱・エントロピー寄与を含む。 |
| **エンタルピー** | (H) | H = E + PV。定圧での全熱含量。 |
| **エントロピー** | (S) | 無秩序さの尺度。ギブズエネルギーに −TS として寄与。 |

---

## 単位・定数

| 用語 | 説明 |
|------|------|
| **Hartree** | 原子単位系のエネルギー。1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV。 |
| **kcal/mol** | 反応エネルギー表現でよく使われる単位。 |
| **kJ/mol** | キロジュール/モル。1 kcal/mol ≈ 4.184 kJ/mol。 |
| **eV** | 電子ボルト。1 eV ≈ 23.06 kcal/mol。 |
| **Bohr** | 原子単位系の長さ。1 Bohr ≈ 0.529 Å。 |
| **Å（オングストローム）** | 10⁻¹⁰ m。原子間距離の表現でよく使われる長さ単位。 |

---

## CLI 規則

| 用語 | 説明 |
|------|------|
| **真偽値オプション** | `True` または `False`（大文字始まり）を取る CLI フラグ。例: `--tsopt True`。 |
| **残基セレクタ** | `'SAM,GPP'`（名前）や `'A:123,B:456'`（チェーン:ID）のような指定方法。 |
| **原子セレクタ** | `'TYR,285,CA'` のように残基名・番号・原子名で特定の原子を指定する方法。 |

---

## 関連ページ

- [はじめに](getting-started.md) — インストールと初回実行
- [概念とワークフロー](concepts.md) — ポケット抽出、MEP 探索、後処理の全体像
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) — 設定ファイルの仕様
- [UMA 計算機](uma_pysis.md) — UMA ポテンシャルの詳細
