# 用語集

このページでは、pdb2reaction ドキュメント内で使われる略語・専門用語を簡潔に説明します。

---

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path | 反応物から生成物へ至る最小エネルギー経路（ポテンシャルエネルギー面上の最も低い経路）。 |
| **TS** | Transition State | ポテンシャルエネルギー面上の一次鞍点（first-order saddle point）。反応座標方向にのみ負の曲率（虚振動数）を 1 つ持つ停留点。 |
| **IRC** | Intrinsic Reaction Coordinate（固有反応座標） | TS から反応物側・生成物側へ向かう、質量重み付き最急降下経路。TS が意図した反応物と生成物を接続していることの検証に使用します。 |
| **GSM** | Growing String Method | 端点からストリング（イメージ列）を伸長・最適化して MEP を近似する手法。 |
| **DMF** | Direct Max Flux | 反応座標方向のフラックスを最大化することで MEP を最適化する chain-of-states 手法。pdb2reaction では `--mep-mode dmf` で選択します。 |
| **NEB** | Nudged Elastic Band | イメージ間にばね力を導入し、間隔を保ちながら経路を最適化する chain-of-states 手法。 |
| **HEI** | Highest-Energy Image | MEP 上でエネルギーが最大のイメージ。TS の初期推定としてよく使われます。 |
| **イメージ（Image）** | — | 経路上の 1 つの構造（1 ノード）。chain-of-states 法で離散化された各点。 |
| **セグメント** | — | 2 つの隣接する端点を結ぶ MEP 区間（例: R → I1, I1 → I2, …）。 |
| **反応セグメント** | Reactive Segment | 端点間で共有結合の変化が検出されるセグメント。反応セグメントのみが TS 最適化に進みます。 |
| **ブリッジセグメント** | Bridge Segment | 非隣接の中間体を結び、未解決の結合変化を含むセグメント。`path-search` がすべての反応領域を分離するまで再帰的に分割します。 |
| **キンク** | Kink | MEP 上で共有結合の変化は検出されないが幾何的な歪みが残る領域。`path-search` は線形補間ノードを挿入し、完全なストリング計算の代わりに個別に最適化します。 |

---

## 最適化アルゴリズム

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **L-BFGS** | Limited-memory BFGS | 勾配履歴からヘシアンを近似する準ニュートン法。`opt --opt-mode grad` で使用。 |
| **RFO** | Rational Function Optimization | 明示的なヘシアン情報を使用する信頼領域最適化法。`opt --opt-mode hess` で使用。 |
| **RS-I-RFO** | Restricted-Step Image-RFO | ヘシアン行列の 1 つの負固有値方向に沿って一次鞍点を探索する RFO 変種。`tsopt --opt-mode hess` で使用。 |
| **Dimer** | Dimer Method | 完全なヘシアンを計算せずに最低曲率モードを推定する TS 最適化法。`tsopt --opt-mode grad` で使用。 |
| **EulerPC** | Euler Predictor-Corrector | IRC 計算の積分スキーム。勾配方向への予測ステップと経路を修正する補正ステップの 2 段階で構成されます。 |
| **PHVA** | Partial Hessian Vibrational Analysis（部分ヘシアン振動解析） | 凍結されていない活性自由度のみで振動解析を行う手法。`freeze_atoms` 設定時に自動適用されます。 |
| **DLC** | Delocalized Internal Coordinates（非局在化内部座標） | 原子間距離・角度・二面角から構成される冗長内部座標系。`coord_type: dlc` で利用可能（デフォルトは `cart` = デカルト座標）。 |

---

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する（多くはニューラルネットの）原子間ポテンシャル。 |
| **UMA** | Universal Machine-learning potential for Atoms | Meta が公開している事前学習 MLIP 群。pdb2reaction のデフォルト計算バックエンドです。 |
| **解析ヘシアン** | Analytical Hessian | エネルギーの正確な二階微分を計算。高速だが VRAM を多く消費。`--hessian-calc-mode Analytical` で選択。 |
| **有限差分** | Finite Difference | 微小変位による微分近似。低速だがメモリ効率が良い。`--hessian-calc-mode FiniteDifference`（デフォルト）で選択。 |

---

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算。 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法。 |
| **Hessian（ヘシアン行列）** | — | エネルギーの二階微分行列。固有値から振動数を、固有ベクトルから振動モード（変位ベクトル）を得ます。振動解析や TS 最適化に使用します。 |
| **SP** | Single Point | 固定構造での計算（最適化なし）。高精度エネルギー補正によく使用。 |
| **スピン多重度** | Spin Multiplicity | 2S+1（S は全スピン量子数）。一重項（singlet）= 1、二重項（doublet）= 2、三重項（triplet）= 3 など。`-m/--multiplicity` で指定（デフォルト: 1）。 |

---

## 構造生物学・ポケット抽出

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース）。 |
| **XYZ** | — | 元素記号と直交座標を並べたシンプルなテキスト形式。 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。pdb2reaction は電荷/多重度と座標の読み取りに利用します。 |
| **ポケット** | Active-site Pocket | `-c/--center` と `-r/--radius` で定義される基質周辺の抽出範囲。 |
| **クラスターモデル** | Cluster Model | ポケットから切り出され、リンク水素でキャップされた計算対象の部分系。MEP/TS 探索の計算量削減に使用します。 |
| **リンク水素** | Link Hydrogen | ポケット抽出時に切断された結合をキャップするために付加する水素原子。 |
| **主鎖** | Backbone | タンパク質の主骨格（N–Cα–C–O 原子）。`--exclude-backbone` で除外可能。 |

---

## 熱化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ZPE** | Zero-Point Energy（零点エネルギー） | 0 K での振動エネルギー。電子エネルギーへの量子補正。 |
| **ギブズエネルギー** | Gibbs Free Energy (G) | G = H - TS。熱・エントロピー寄与を含む自由エネルギー。 |
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
| **Å（オングストローム）** | 10⁻¹⁰ m。原子間距離の標準単位。 |
| **cm⁻¹** | 波数（逆センチメートル）。振動数の標準単位。虚振動数は負の値で表されます。 |

---

## CLI 規則

| 用語 | 説明 |
|------|------|
| **ブール値オプション** | トグル形式（`--flag` / `--no-flag`）または値形式（`True`/`False`、`yes`/`no`、`1`/`0`）を受け付ける CLI フラグ。例: `--tsopt`。 |
| **残基セレクタ** | `'SAM,GPP'`（名前）や `'A:123,B:456'`（チェーン:ID）のような指定方法。 |
| **原子セレクタ** | `'TYR,285,CA'` のように残基名・番号・原子名で特定の原子を指定する方法。 |

---

## 関連ページ

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting_started.md) — クイックスタートと初回実行
- [概念とワークフロー](concepts.md) — ポケット抽出、MEP 探索、後処理の全体像
- [典型エラー別レシピ](recipes_common_errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml_reference.md) — 設定ファイルの仕様
- [MLIP 計算機](uma_pysis.md) — MLIP バックエンドの詳細
