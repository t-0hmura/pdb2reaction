# 用語集

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path | 反応物から生成物へ至る最小エネルギー経路（ポテンシャルエネルギー面上の最も低い経路） |
| **TS** | Transition State | ポテンシャルエネルギー面上の一次鞍点（first-order saddle point）。反応座標方向にのみ負の曲率（虚振動数）を 1 つ持つ停留点 |
| **IRC** | Intrinsic Reaction Coordinate（固有反応座標） | 古典的には TS から反応物側・生成物側へ向かう**質量重み付き**最急降下経路として定義され、TS が意図した端点を接続することの検証に使用します。pdb2reaction の EulerPC 積分は**質量重みなしのデカルト座標**で動作します（`irc.md`: `--step-size` は質量重みなしデカルトの Bohr）。`geom.coord_type` は `cart` に強制されます |
| **GSM** | Growing String Method | 端点からストリング（イメージ列）を伸長・最適化して MEP を近似する手法 |
| **DMF** | Direct Max Flux | 反応座標方向のフラックスを最大化することで MEP を最適化する chain-of-states 手法。pdb2reaction では `--mep-mode dmf` で選択します |
| **HEI** | Highest-Energy Image | MEP 上でエネルギーが最大のイメージ。TS の初期推定としてよく使われます |
| **イメージ（Image）** | — | 経路上の 1 つの構造（1 ノード）。chain-of-states 法で離散化された各点 |
| **セグメント** | — | 2 つの隣接する端点を結ぶ MEP 区間（例: R → I1, I1 → I2, …） |
| **反応セグメント** | Reactive Segment | 端点間で共有結合の変化が検出されるセグメント。反応セグメントのみが TS 最適化に進みます |
| **ブリッジセグメント** | Bridge Segment | 非隣接の中間体を結び、未解決の結合変化を含むセグメント。`path-search` がすべての反応領域を分離するまで再帰的に分割します |
| **ねじれ** | Kink | MEP 上で共有結合の変化は検出されないが幾何的な歪みが残る領域。`path-search` は線形補間ノードを挿入し、完全なストリング計算の代わりに個別に最適化します |
| **PES** | Potential Energy Surface（ポテンシャルエネルギー面） | 原子配置に対するエネルギーの超曲面。MEP は PES 上の最低エネルギー経路 |

## 最適化アルゴリズム

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **BFGS** | Broyden-Fletcher-Goldfarb-Shanno | 準ニュートン型のHessian更新スキーム（`hessian_update: bfgs`） |
| **L-BFGS** | Limited-memory BFGS | 勾配履歴からHessianを近似する準ニュートン法。`opt --opt-mode grad` で使用 |
| **RFO** | Rational Function Optimization | 明示的なHessian情報を使用する信頼領域最適化法。`opt --opt-mode hess` で使用 |
| **RS-I-RFO** | Restricted-Step Image-RFO | Hessian行列の 1 つの負固有値方向に沿って一次鞍点を探索する RFO 変種。`tsopt --opt-mode rsirfo` で選択（`hess` のデフォルトは RS-I-RFO） |
| **Dimer** | Dimer Method | 完全なHessianを計算せずに最低曲率モードを推定する TS 最適化法。`tsopt --opt-mode grad` で使用。pdb2reaction は Hessian Guided Dimer 変種を使い、活性部分空間の正確なHessianを周期的に評価してダイマー方向を更新します |
| **Bofill** | Bofill Update | SR1（対称ランク 1）と PSB（Powell-symmetric-Broyden）を混合したHessian更新スキーム。鞍点探索に適します。`hessian_update: bofill` で選択され、RS-I-RFO と Dimer のフラット化ループで使用されます |
| **SR1** | Symmetric Rank-One | ランク 1 のHessian更新スキーム。Bofill の 2 要素のうち 1 つ |
| **PSB** | Powell-Symmetric-Broyden | 対称Hessian更新スキーム。Bofill の 2 要素のうちもう 1 つ |
| **EulerPC** | Euler Predictor-Corrector | IRC 計算の積分スキーム。勾配方向への予測ステップと経路を修正する補正ステップの 2 段階で構成されます |
| **PHVA** | Partial Hessian Vibrational Analysis（部分Hessian振動解析） | 凍結されていない活性自由度のみで振動解析を行う手法。`freeze_atoms` 設定時に自動適用されます |
| **Active DOF** | Active Degrees of Freedom（活性自由度） | `freeze_atoms` に列挙されていない原子が持つ 3N 個のデカルト座標。PHVA、部分Hessian TS 最適化、解析的Hessian経路はすべてこの活性部分空間のみで動作し、凍結原子は縮約Hessianの行・列を寄与しません |
| **DLC** | Delocalized Internal Coordinates（非局在化内部座標） | 原子間距離・角度・二面角から構成される冗長内部座標系。`coord_type: dlc` で利用可能（デフォルトは `cart` = デカルト座標） |

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する（多くはニューラルネットの）原子間ポテンシャル |
| **UMA** | Universal Models for Atoms | Meta が公開している事前学習 MLIP 群。pdb2reaction のデフォルト計算バックエンドです |
| **ORB** | ORB Models | Orbital Materials の MLIP バックエンド。`-b orb` で選択 |
| **MACE** | MACE | Equivariant message-passing MLIP。`-b mace` で選択 |
| **AIMNet2** | Atoms-in-Molecules Neural Network Potential, 2nd generation | 電荷対応ニューラルネットワークポテンシャル（Anstine et al., *Chem. Sci.* 2025）。`-b aimnet2` で選択 |
| **xTB** | Extended Tight Binding | 半経験的量子化学手法。pdb2reaction では `--solvent` による暗黙溶媒補正に使用 |
| **fairchem** | — | Meta がオープンソースで公開している基盤モデルツールキット。UMA 系のチェックポイントを提供します。pdb2reaction は UMA 予測器のロードに `fairchem-core` へ依存します |
| **ASE** | Atomic Simulation Environment | pdb2reaction の MLIP バックエンド全てが利用する Calculator API を提供する Python フレームワーク（Larsen et al., *J. Phys. Condens. Matter* 2017）。 |
| **task_name** | — | UMA の推論バッチに記録されるタスクタグ（YAML: `calc.task_name`、デフォルト `omol`）。チェックポイントが学習したタスク/プリセットを選択します |
| **解析Hessian** | Analytical Hessian | エネルギーの正確な二階微分を計算。高速だが VRAM を多く消費。`--hessian-calc-mode Analytical` で選択 |
| **有限差分** | Finite Difference | 微小変位による微分近似。低速だがメモリ効率が良い。`--hessian-calc-mode FiniteDifference`（デフォルト）で選択 |

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法 |
| **DFT//MLIP** | — | 複合計算法の表記: MLIP で最適化した構造の上で DFT による単点エネルギーを評価する手法。MLIP の構造最適化・分子動力学と DFT のエネルギー精密化を組み合わせます。`//` 区切りは量子化学の標準慣習「エネルギー計算レベル // 構造最適化レベル」に従います |
| **Hessian（Hessian行列）** | — | エネルギーの二階微分行列。固有値から振動数を、固有ベクトルから振動モード（変位ベクトル）を得ます。振動解析や TS 最適化に使用します |
| **SP** | Single Point | 固定構造での計算（最適化なし）。高精度エネルギー補正によく使用 |
| **スピン多重度** | Spin Multiplicity | 2S+1（S は全スピン量子数）。一重項（singlet）= 1、二重項（doublet）= 2、三重項（triplet）= 3 など。`-m/--multiplicity` で指定（デフォルト: 1） |
| **ALPB** | Analytical Linearized Poisson-Boltzmann | xTB で利用可能な暗黙溶媒モデル（`--solvent-model alpb`、デフォルト） |
| **CPCM-X** | 拡張型 Conductor-like Polarizable Continuum Solvation Model | xTB で利用可能な暗黙溶媒モデル（`--solvent-model cpcmx`）。"X" は "eXtended" を意味し、CPCM に COSMO-RS 由来の σ-profile と SMD 流の非静電項を結合することで任意溶媒に対応する溶媒和自由エネルギーを与える（Stahn, Ehlert, Grimme, *J. Phys. Chem. A* 2023）。 |
| **cyipopt** | — | IPOPT 内点法ソルバの Python バインディング。DMF（`--mep-mode dmf`）経路精密化パイプラインが依存します |
| **IPOPT** | Interior Point OPTimizer | 非線形制約付き最適化のオープンソース solver（Wächter & Biegler 2006）。DMF 経路精密化で `cyipopt` 経由で使用されます。 |
| **SCF** | Self-Consistent Field | DFT/HF で電子波動関数を反復収束させる手続き。`pdb2reaction dft` では `--max-cycle` / `--conv-tol` で制御されます。 |

## 構造生物学・活性部位モデル抽出

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース） |
| **XYZ** | — | 元素記号と直交座標を並べたシンプルなテキスト形式 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。pdb2reaction は電荷/多重度と座標の読み取りに利用します |
| **活性部位モデル** | Active Site Model（バインディングポケット） | `-c/--center` と `-r/--radius` で定義される基質周辺の抽出範囲。pdb2reaction ドキュメントでは「活性部位モデル」と「クラスターモデル」を下流の計算入力の意味でほぼ同義に使います。厳密には「活性部位モデル」は幾何学的な選択領域、「クラスターモデル」はそれをキャップ水素でキャップして QM/MLIP 計算入力として完成させた段階を指します |
| **クラスターモデル** | Cluster Model | 抽出した活性部位モデルの切断された共有結合をキャップ水素でキャップし、QM/MLIP 計算入力として整えた部分系。実質的には「活性部位モデル」（キャップ前の幾何選択）と「クラスターモデル」（キャップ後、計算実行可能）は同じ対象の 2 段階です |
| **キャップ水素** | Cap Hydrogen | 活性部位モデル抽出時に切断された結合をキャップするために付加する水素原子 |
| **主鎖** | Backbone | タンパク質の主骨格（N–Cα–C–O 原子）。活性部位モデル抽出時に `--exclude-backbone` で除外可能 |

## 熱化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ZPE** | Zero-Point Energy（零点エネルギー） | 0 K での振動エネルギー。電子エネルギーへの量子補正 |
| **ギブズ自由エネルギー** | Gibbs Free Energy (G) | G = H - TS。熱・エントロピー寄与を含む自由エネルギー |
| **エンタルピー** | (H) | H = E + PV。定圧での全熱含量 |
| **エントロピー** | (S) | 無秩序さの尺度。ギブズ自由エネルギーに −TS として寄与 |
| **QRRHO** | Quasi-Rigid-Rotor Harmonic Oscillator | Grimme の低振動数補正を含む熱化学近似。`freq` で自動適用 |

## 単位・定数

| 用語 | 説明 |
|------|------|
| **Hartree** | 原子単位系のエネルギー。1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV |
| **RMSD** | Root-Mean-Square Deviation。`path-search` のセグメント stitch/bridge 類似性指標（`stitch_rmsd_thresh`, `bridge_rmsd_thresh`）として使用されます。 |
| **MAE** | Mean Absolute Error。ベンチマークや回帰のレポートで使用されます。 |
| **CI** | Confidence Interval（統計）。量子化学の "configuration interaction" とは別で、pdb2reaction ではベンチマーク統計の文脈でのみ使用します。 |
| **kcal/mol** | 反応エネルギー表現でよく使われる単位 |
| **kJ/mol** | キロジュール/モル。1 kcal/mol ≈ 4.184 kJ/mol |
| **eV** | 電子ボルト。1 eV ≈ 23.06 kcal/mol |
| **Bohr** | 原子単位系の長さ。1 Bohr ≈ 0.529 Å |
| **Å（オングストローム）** | 10⁻¹⁰ m。原子間距離の標準単位 |
| **cm⁻¹** | 波数（逆センチメートル）。振動数の標準単位。虚振動数は負の値で表されます |
| **虚振動数** | Hessian行列の負の固有値に対応する振動数。TS では 1 本のみ存在（一次鞍点）。負の cm⁻¹ 値で報告されます。検出閾値は `hessian_dimer.neg_freq_thresh_cm`（デフォルト 5 cm⁻¹） |

(ja-frequency-thresholds)=
### 振動数閾値: 5 cm⁻¹（虚振動検出）と 100 cm⁻¹（QRRHO rotor cutoff）

`pdb2reaction` 内には独立した 2 つの cm⁻¹ 閾値があり、作用するモード集合と用途が全く異なります:

| 閾値 | 役割 | 定義場所 |
|------|------|----------|
| **5 cm⁻¹** | *虚振動モード検出カットオフ*。絶対値がこれ未満の負の固有値は虚振動と見なされない（剛体運動または数値ノイズとして扱われる） | `pdb2reaction/core/defaults.py` の `hessian_dimer.neg_freq_thresh_cm = 5.0`。YAML で調整可 |
| **100 cm⁻¹** | *QRRHO rotor cutoff*（Grimme）。`freq` の熱化学計算において、これ未満の **正の** 低振動モードは harmonic-oscillator から自由回転子の entropy へ滑らかに移行する。entropy / Gibbs 自由エネルギーのみに影響 | `thermoanalysis/config.py` の `ROTOR_CUT_DEFAULT = 100.0` |

## CLI 規則

| 用語 | 説明 |
|------|------|
| **ブール値オプション** | 切替形式（`--flag` / `--no-flag`）または値形式（`True`/`False`、`yes`/`no`、`1`/`0`）を受け付ける CLI フラグ。例: `--tsopt` |
| **残基セレクタ** | `'SAM,GPP'`（名前）や `'A:123,B:456'`（チェーン:ID）のような指定方法 |
| **原子セレクタ** | `'TYR,285,CA'` のように残基名・番号・原子名で特定の原子を指定する方法 |

## 関連ページ

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting-started.md) — クイックスタートと初回実行
- [all](all.md) — 活性部位モデル抽出、MEP 探索、後処理の全体像
- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) — 設定ファイルの仕様
- [MLIP 計算機](uma-pysis.md) — MLIP バックエンドの詳細
