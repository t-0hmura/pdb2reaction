# 用語集

pdb2reaction ドキュメントで使用される略語と専門用語の定義です。

---

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path（最小エネルギー経路） | 反応物から生成物へ遷移状態を経由する最低エネルギー経路。 |
| **TS** | Transition State（遷移状態） | ポテンシャルエネルギー曲面上の鞍点で、反応座標に沿った最高エネルギー点。 |
| **IRC** | Intrinsic Reaction Coordinate（固有反応座標） | TSから反応物・生成物方向への質量重み付き最急降下経路。 |
| **GSM** | Growing String Method | 両端から中央に向かってイメージを反復的に成長させることでMEPを探索するアルゴリズム。 |
| **NEB** | Nudged Elastic Band | バネ力を使用して反応経路に沿ったイメージ間隔を維持するchain-of-states法。 |

---

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential（機械学習原子間ポテンシャル） | 量子力学データで訓練されたニューラルネットワークベースのエネルギー/力予測器。 |
| **UMA** | Universal Machine-learning potential for Atoms | pdb2reactionでデフォルト計算機として使用されるMetaの事前学習済みMLIPファミリー。 |
| **DMF** | Distance Matrix Fingerprints（距離行列フィンガープリント） | 原子間距離に基づく構造記述子。GSMでの類似度チェックに使用。 |

---

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics（量子力学） | DFT、HF、post-HF等の第一原理電子構造計算。 |
| **DFT** | Density Functional Theory（密度汎関数理論） | 電子密度汎関数を介して電子構造をモデル化する量子力学的手法。 |
| **Hessian**（ヘシアン） | — | 原子座標に対するエネルギーの二次微分行列。振動解析とTS最適化に使用。 |

---

## 構造生物学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | 高分子3D構造のファイル形式およびデータベース。 |
| **XYZ** | — | 原子記号とデカルト座標を列挙するシンプルなテキスト形式。 |
| **GJF** | Gaussian Job File | Gaussianの入力形式。pdb2reactionはこれらのファイルから電荷/多重度と座標を読み取る。 |

---

## 単位・定数

| 用語 | 説明 |
|------|------|
| **Hartree** | エネルギーの原子単位。1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV。 |
| **kcal/mol** | キロカロリー毎モル。反応エネルギー論の一般的な単位。 |
| **Bohr** | 長さの原子単位。1 Bohr ≈ 0.529 Å。 |
| **Angstrom (Å)** | 10⁻¹⁰ メートル。原子間距離の標準単位。 |

---

## 関連項目

- [はじめに](getting-started.md) — インストールと初回実行
- [YAML リファレンス](../yaml-reference.md) — 設定ファイル形式
- [UMA 計算機](uma_pysis.md) — 機械学習ポテンシャルの詳細
