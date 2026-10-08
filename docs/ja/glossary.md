# 用語集

ドキュメントに出てくる略語・手法名・単位を、分野ごとに 1 行で引くページです。オプション、出力の欄、ステータスの値は、各コマンドのページと [JSON 出力の一覧](json-output.md) にあります。

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path | 反応物から生成物へ至る最小エネルギー経路（ポテンシャルエネルギー面上の最も低い経路） |
| **TS** | Transition State | ポテンシャルエネルギー面上の一次鞍点（first-order saddle point）。反応座標方向にのみ負の曲率（虚振動数）を 1 つ持つ停留点 |
| **n_imag** | Number of imaginary modes（虚振動の数） | 虚振動の分類基準（既定では ν < −5.00 cm⁻¹）より下の振動モードの数。TS では n_imag = 1 で、`result.json` には `tsopt` では `n_imaginary_modes`、`freq` では `n_imaginary` として記録されます |
| **IRC** | Intrinsic Reaction Coordinate（固有反応座標） | 古典的には TS から反応物側・生成物側へ向かう**質量重み付き**最急降下経路として定義され、TS が意図した端点を接続することの検証に使用します。pdb2reaction の EulerPC 積分も質量重み付き座標を進めます。`--step-size` は質量加重しないデカルト座標で測る長さ（Bohr）です |
| **GSM** | Growing String Method | 端点からストリング（イメージ列）を伸長・最適化して MEP を近似する手法 |
| **DMF** | Direct Max Flux | 反応座標方向のフラックスを最大化することで MEP を最適化する chain-of-states 手法。pdb2reaction では `--mep-mode dmf` で選択します |
| **FB-ENM** | Flat-Bottom Elastic Network Model | DMF の初期経路を作るモデル（Koda & Saito, *J. Chem. Theory Comput.* 2024）。CFB-ENM はその相関を入れた変種です |
| **HEI** | Highest-Energy Image | MEP 上でエネルギーが最大のイメージ。TS の初期推定としてよく使われます。HEI±1 はその両隣のイメージです |
| **イメージ（Image）** | — | 経路上の 1 つの構造（1 ノード）。chain-of-states 法で離散化された各点 |
| **COS** | Chain-of-States | イメージの鎖をまとめて最適化する経路の手法。GSM、ストリングオプティマイザ（`stopt`）、DMF など |
| **セグメント** | — | 2 つの隣接する端点を結ぶ MEP 区間（例: R → I1, I1 → I2, …） |
| **反応セグメント** | Reactive Segment | TS 候補を持つセグメント。ブリッジとねじれを除くすべての MEP のセグメント（ふつうは両端で共有結合が変わる区間）と、[TS-only モード](quickstart-tsopt.md)（1 構造に `all --tsopt` を付けた実行）で入力した TS がこれに当たり、`all` は要求した TS 最適化・IRC・熱化学・DFT をこれにだけ行います |
| **ブリッジセグメント** | Bridge Segment | 隣り合うセグメントの間をつなぐ短い MEP。`path-search` がセグメントを 1 本の経路につなぐとき、前のセグメントの終わりと次のセグメントの始まりが一致せず、その間で結合が変わらなければ、そのすき間をブリッジで埋めます（[path-search の処理の仕組み](path-search.md#処理の仕組みと計算仕様) の 5） |
| **ねじれ** | Kink | 配座だけが変わる経路の区間（セグメント）。`path-search` が HEI の両側で最適化した 2 つの構造（End1 と End2。[path-search の処理の仕組み](path-search.md#処理の仕組みと計算仕様) の 2）の間で、共有結合が変わらない区間を指します。`path-search` は新しい GSM・DMF の経路の代わりに、線形補間のノードを数個（`search.kink_max_nodes`、既定 3）入れて 1 つずつ最適化します |
| **PES** | Potential Energy Surface（ポテンシャルエネルギー面） | 原子配置に対するエネルギーの超曲面。MEP は PES 上の最低エネルギー経路 |

## 最適化アルゴリズム

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **BFGS** | Broyden-Fletcher-Goldfarb-Shanno | 準ニュートン型の Hessian 更新スキーム（`hessian_update: bfgs`） |
| **TS-BFGS** | Transition-State BFGS | ステップ方向の曲率が正であることを要求しない BFGS 型の Hessian 更新スキーム。モデル Hessian が不定値のままでも更新できます。RFO による極小化の既定（`hessian_update: ts_bfgs`） |
| **L-BFGS** | Limited-memory BFGS | 勾配履歴から Hessian を近似する準ニュートン法。`opt --opt-mode grad` で使用 |
| **RFO** | Rational Function Optimization | 明示的な Hessian 情報を使用する信頼領域最適化法。`opt --opt-mode hess` で使用 |
| **RS-I-RFO** | Restricted-Step Image-RFO | Hessian 行列の 1 つの負固有値方向に沿って一次鞍点を探索する RFO 変種。`tsopt --opt-mode rsirfo` で選択（`hess` のデフォルトは RS-P-RFO） |
| **RS-P-RFO** | Restricted-Step Partitioned RFO（制限ステップ分割 RFO） | `tsopt` の既定の TS 最適化法（`--opt-mode hess` または `rsprfo`）。TS モードに沿ってエネルギーを最大化し、ほかのすべてのモードに沿って最小化する RFO 変種です（Banerjee et al. 1985） |
| **TRIM** | Trust-Region Image Minimization（信頼領域イメージ最小化） | TS モードに沿った Hessian の固有値と勾配の符号を反転させてから最小化する、信頼領域の TS 最適化法（Helgaker 1991）。`tsopt --opt-mode trim` で選択 |
| **Flatten** | — | `--flatten`（`opt`・`tsopt`・`all`、既定で無効）を付けると、最適化の後に余分な虚振動のモードに沿って構造をずらして最適化し直し、`opt` では虚振動が無くなるまで、`tsopt` では 1 つになるまで、または回数の上限まで繰り返します |
| **Dimer** | Dimer Method | 勾配を使うステップで最低曲率モードを追う TS 最適化法。pdb2reaction は Hessian-guided Dimer 変種を使い、動ける原子の厳密な Hessian で最初にダイマー方向を決め、`hessian_dimer.update_interval_hessian` ステップ（既定 500）ごとに更新します。`tsopt --opt-mode dimer` で選択（`grad` は別名） |
| **Bofill** | Bofill Update | SR1（対称ランク 1）と PSB（Powell-symmetric-Broyden）を混合した Hessian 更新スキーム。鞍点探索に適します。`rsirfo` セクション（RS-P-RFO・RS-I-RFO・TRIM）と `irc` の `hessian_update` の既定値です |
| **SR1** | Symmetric Rank-One | ランク 1 の Hessian 更新スキーム。Bofill の 2 要素のうち 1 つ |
| **PSB** | Powell-Symmetric-Broyden | 対称 Hessian 更新スキーム。Bofill の 2 要素のうちもう 1 つ |
| **EulerPC** | Euler Predictor-Corrector | IRC 計算の積分スキーム。勾配方向への予測ステップと経路を修正する補正ステップの 2 段階で構成されます |
| **PHVA** | Partial Hessian Vibrational Analysis（部分 Hessian 振動解析） | 凍結されていない活性自由度のみで振動解析を行う手法。`freeze_atoms` 設定時に自動適用されます |
| **Active DOF** | Active Degrees of Freedom（活性自由度） | `freeze_atoms` に列挙されていない原子が持つ 3N 個のデカルト座標。PHVA、部分 Hessian TS 最適化、解析的 Hessian 経路はすべてこの活性部分空間のみで動作し、凍結原子は縮約 Hessian の行・列を寄与しません |
| **DLC** | Delocalized Internal Coordinates（非局在化内部座標） | 原子間距離・角度・二面角（基本内部座標）の線形結合で作る、冗長性のない内部座標系。`--coord-type dlc`（YAML の `geom.coord_type: dlc`）で選択。デフォルトは `cart`（デカルト座標） |

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する（多くはニューラルネットの）原子間ポテンシャル |
| **UMA** | Universal Models for Atoms | Meta が公開している事前学習 MLIP 群。pdb2reaction のデフォルト計算バックエンドです |
| **ORB** | ORB Models | Orbital Materials の MLIP バックエンド。`-b orb` で選択 |
| **MACE** | MACE | 同変メッセージパッシング MLIP。`-b mace` で選択 |
| **AIMNet2** | Atoms-in-Molecules Neural Network Potential, 2nd generation | 電荷対応ニューラルネットワークポテンシャル（Anstine et al., *Chem. Sci.* 2025）。`-b aimnet2` で選択 |
| **fairchem** | — | Meta がオープンソースで公開している基盤モデルツールキット。UMA 系のチェックポイントを提供します。pdb2reaction は UMA 予測器のロードに `fairchem-core` へ依存します |
| **ASE** | Atomic Simulation Environment | pdb2reaction の MLIP バックエンド全てが利用する Calculator API を提供する Python フレームワーク（Larsen et al., *J. Phys. Condens. Matter* 2017）。 |
| **task_name** | — | UMA の推論バッチに記録されるタスクタグ（YAML: `calc.task_name`、デフォルト `omol`）。チェックポイントが学習したタスク/プリセットを選択します |
| **解析 Hessian** | Analytical Hessian | 選択した MLIP のエネルギーを自動微分で二階微分し、有限変位の打切り誤差を避ける（浮動小数点と自動微分の誤差は残る）。計算時間と GPU などのメモリ使用量は、バックエンド・モデル・系によって変わります。`--hessian-calc-mode Analytical` で選択 |
| **有限差分** | Finite Difference | 原子を有限変位させて Hessian を近似する、どの環境でも使えるデフォルト。通常は GPU などのメモリ使用量の最大値が小さいが、計算時間と変位による誤差は設定によって変わります。`--hessian-calc-mode FiniteDifference` で選択 |

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法 |
| **DFT//MLIP** | — | 複合計算法の表記: MLIP で最適化した構造の上で DFT による一点エネルギーを評価する手法。MLIP の構造最適化・分子動力学と、より高い理論レベルでの DFT エネルギー評価を組み合わせます。`//` 区切りは量子化学の標準慣習「エネルギー計算レベル // 構造最適化レベル」に従います |
| **MLIP_Gibbs**、**DFT//MLIP_Gibbs** | — | `summary.json` の手法の名前で、`MEP`・`MLIP`・`DFT` と並びます。`MLIP_Gibbs` は MLIP のギブズ自由エネルギー、`DFT//MLIP_Gibbs` は DFT のエネルギーに MLIP の熱補正を足した値です（`-b dft` のときは、`MLIP` と `MLIP_Gibbs` の代わりに `DFT` と `DFT_Gibbs` になります） |
| **Hessian（Hessian 行列）** | — | 原子座標に関するエネルギーの二階微分行列。固有値から振動数、固有ベクトルから振動モード（変位ベクトル）が得られます。振動解析や TS 最適化に使用します |
| **SP** | Single Point | 固定構造での計算（最適化なし）。より高い理論レベルでのエネルギー評価によく使用 |
| **スピン多重度** | Spin Multiplicity | 2S+1（S は全スピン量子数）。一重項（singlet）= 1、二重項（doublet）= 2、三重項（triplet）= 3 など。`-m/--multiplicity` で指定（デフォルト: 1） |
| **cyipopt** | — | IPOPT 内点法ソルバの Python バインディング。DMF（`--mep-mode dmf`）経路精密化パイプラインが依存します |
| **IPOPT** | Interior Point OPTimizer | 非線形制約付き最適化のオープンソースのソルバ（Wächter & Biegler 2006）。DMF 経路精密化で `cyipopt` 経由で使用されます。 |
| **SCF** | Self-Consistent Field | DFT/HF で電子波動関数を反復収束させる手続き。`pdb2reaction dft` では `--scf-max-cycles` / `--scf-tol` で制御されます。 |

## 構造生物学・活性部位モデル抽出

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース） |
| **XYZ** | — | 元素記号と直交座標を並べたシンプルなテキスト形式 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。pdb2reaction は電荷/多重度と座標の読み取りに利用します |
| **活性部位モデル** | Active Site Model（バインディングポケット） | `-c/--center` と `-r/--radius` で選ぶ基質周辺の領域。このドキュメントでは、活性部位モデルとクラスターモデルはどちらも、`extract` が書き出す、切断した結合にキャップ水素を付けたモデルを指します |
| **クラスターモデル** | Cluster Model | 抽出した領域の切断された共有結合をキャップ水素でキャップし、QM/MLIP 計算入力として整えた部分系 |
| **キャップ水素** | Cap Hydrogen | 活性部位モデル抽出時に切断された結合をキャップするために付加する水素原子。link 水素とも呼ぶ（`--add-linkh`、`--freeze-links`、`n_link_hydrogens`） |
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
| **RMSD** | Root-Mean-Square Deviation。`path-search` でセグメントを接合（stitch）・橋渡し（bridge）するときの類似性指標（`stitch_rmsd_thresh`, `bridge_rmsd_thresh`）として使用されます。 |
| **kcal/mol** | 反応エネルギー表現でよく使われる単位 |
| **kJ/mol** | キロジュール/モル。1 kcal/mol ≈ 4.184 kJ/mol |
| **eV** | 電子ボルト。1 eV ≈ 23.06 kcal/mol |
| **Bohr** | 原子単位系の長さ。1 Bohr ≈ 0.529 Å |
| **Å（オングストローム）** | 10⁻¹⁰ m。原子間距離の標準単位 |
| **cm⁻¹** | 波数（逆センチメートル）。振動数の標準単位。虚振動数は負の値で表されます |
| **虚振動数** | Hessian 行列の負の固有値に対応する振動数。TS では 1 本のみ存在（一次鞍点）。負の cm⁻¹ 値で報告されます。 |

(ja-frequency-thresholds)=
### 虚振動の分類基準と QRRHO のローター閾値

| 閾値 | 役割 | 定義場所 |
|------|------|----------|
| **ν < −5.00 cm⁻¹** | 既定の虚振動分類基準。 | 設定で変えられます（`freq.zero_cutoff_cm`） |
| **100 cm⁻¹** | *QRRHO のローター閾値*（Grimme）。`freq` の熱化学計算では、これ未満の **正の** 低振動モードのエントロピーを、調和振動子の値から自由回転子の値へ滑らかに切り替えます。変わるのはエントロピーとギブズ自由エネルギーだけです | 固定（pdb2reaction の設定では変えられません） |

## CLI 規則

ブール値オプション、残基セレクタ、原子セレクタの書き方は [共通オプションと残基・原子の指定](cli-conventions.md) にまとめています。

## 使用上の注意点

* **虚振動と負の符号**: 振動数はすべて符号付きで出力されます。負の値の総数 `n_negative_modes`（`result.json`）は n_imag とは別の診断値で、収束の判定を変えません。

## 関連ドキュメント

- [はじめに](getting-started.md) — 最短の実行と次に読むページ
- [インストール](installation.md) — セットアップと依存関係
- [all](all.md) — ワークフロー全体
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML 設定の一覧](yaml-reference.md) — 設定ファイルの仕様
- [MLIP バックエンド](backends.md) — MLIP バックエンドの詳細
