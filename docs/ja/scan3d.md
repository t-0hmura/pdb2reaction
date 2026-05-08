# `scan3d`

## 概要

> **要約:** 調和拘束と MLIP 緩和により、3 距離（d₁, d₂, d₃）のグリッドスキャンを行います。`-s/--scan-lists`（YAML/JSON ファイルパス、推奨）またはインライン Python リテラルでターゲットを指定します。`--csv` では既存 `surface.csv` の可視化のみ実行します。

### 要点
- **想定場面:** 3 つの距離 `(d₁, d₂, d₃)` 上で 3D ポテンシャル体積を得たいとき、または既存の `surface.csv` を再プロットしたいとき。入力は 1 つの構造 + `-s scan3d.yaml`（推奨）または `-s/--scan-lists` の **単一** インラインリテラル（四つ組は 3 つ）。`--csv` 指定時はプロットのみで実行可能です。
- **計算手法:** d₁ → d₂ → d₃ のネストループで `--max-step-size` から線形グリッドを構成し、事前最適化構造に近い値が先に走査されるよう各軸を並べ替えます。各点は対応する調和拘束を有効にして MLIP バックエンド（デフォルト: UMA）で緩和し、記録されるエネルギーは **バイアスを除去して**評価されるため、格子点間で直接比較できます。
- **主な出力:** `surface.csv`、`grid/` 配下の各点の構造、HTML の等値面図（`scan3d_density.html`）。
- **デフォルト値:** `--opt-mode grad`（LBFGS）、`--no-preopt`、`--max-step-size 0.20 Å`、`--bias-k 300 eV·Å⁻²`、`--thresh baker`、`--baseline min`、`--out-dir ./result_scan3d/`。3D グリッドは点数が急激に増加するため、まず `--max-step-size` を大きくするか範囲を狭めることを検討してください。
- **次のステップ:** `scan3d_density.html` で低エネルギーチャネルを確認し、2D の `scan2d` スライスで絞り込むか、`tsopt` で TS 候補を精密化します。

`scan3d` は d₁ → d₂ → d₃ の順にループをネストし、対応する拘束をかけて各格子点を緩和します。デフォルトは LBFGS（`--opt-mode grad`）で、RFOptimizer が必要な場合は `--opt-mode hess` を指定してください。

XYZ/GJF 入力では、`--ref-pdb` で参照 PDB トポロジーを指定すると、XYZ 座標を保持したまま PDB/GJF へのフォーマット対応変換が可能になります。

## 最小例
```bash
pdb2reaction scan3d -i input.pdb -q 0 -s scan3d.yaml -o ./result_scan3d/
```

## 出力の見方
- `result_scan3d/surface.csv`
- `result_scan3d/grid/point_iDDD_jDDD_kDDD.xyz`（`DDD = round(d × 100)`、例: d₁=1.30, d₂=3.10, d₃=2.00 Å → `point_i130_j310_k200.xyz`）
- `result_scan3d/scan3d_density.html`

> **Note:** `-s/--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。

## よくある例
```bash
# 推奨: YAML/JSON spec
cat > scan3d.yaml << 'YAML'
one_based: true
pairs:
 - ["TYR,285,CA", "SAM,309,C10", 1.30, 3.10]
 - ["TYR,285,CB", "SAM,309,C11", 1.20, 3.20]
 - ["TYR,285,CG", "SAM,309,C12", 1.10, 3.00]
YAML
pdb2reaction scan3d -i input.pdb -q 0 -s scan3d.yaml

# 代替: Python リテラル
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]'

# LBFGS 緩和、内側軌跡ダンプ、HTML 等値面プロット
pdb2reaction scan3d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20),("TYR,285,CG","SAM,309,C12",1.10,3.00)]' \
 --max-step-size 0.20 --dump -o ./result_scan3d/ --opt-mode grad \
 --preopt --baseline min

# 既存 surface.csv からのプロットのみ（スキャンしない）
pdb2reaction scan3d --csv ./result_scan3d/surface.csv --zmin -10 --zmax 40 -o ./result_scan3d/
```

## 使用法
```bash
pdb2reaction scan3d [-i INPUT.{pdb|xyz|trj|...}] [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan3d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE] [--csv PATH]
```
注: `-i/--input` と `-s/--scan-lists` は `--csv` が指定されていない限り必須です。

## スキャンリスト仕様

`scan3d` はちょうど **3 つ**の四つ組 `(i, j, low_Å, high_Å)` を受け付けます（YAML/JSON では `pairs` キー、インラインでは単一リテラル）。`scan` と異なり、リテラルは **1 つだけ**を受け付けます（複数ステージは非対応）。

YAML/JSON ファイル書式、インライン Python リテラル構文、原子セレクタ、クォート規則については
{ref}`CLI 規約: スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

## ワークフロー
1. `geom_loader` で構造を読み込み、CLI または Gaussian テンプレートから電荷とスピンを解決します。`--preopt` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` が与えられている場合、構造は酵素--基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付き XYZ/GJF）で `extract.py` の電荷サマリーから総電荷を導出します。
2. `-s/--scan-lists`（YAML/JSON ファイルパスまたはインライン Python リテラル、デフォルト 1 始まり、`--zero-based` で 0 始まり）を 3 つの四つ組に解析します。PDB 入力では、各原子指定に整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を使用できます。区切りは空白・カンマ・スラッシュ・バッククォート・バックスラッシュのいずれも可で、トークン順序は任意です（不明な場合は resname, resseq, atom の順と仮定）。`h = --max-step-size` で各距離の線形グリッドを生成し、開始距離に近い値が先に走査されるよう並べ替えます。
3. 外側ループで `d1[i]` を走査し、**d₁ 拘束のみ**を適用して緩和します。近い d₁ 値の既存構造から開始します。
4. 中間ループで `d2[j]` を走査し、**d₁ + d₂ 拘束**を適用して緩和します。近い (d₁, d₂) の構造から開始します。
5. 内側ループで `d3[k]` を走査し、**3 つの拘束すべて**を適用して緩和します。バイアスを除去したエネルギーを測定し、構造と収束フラグを書き出します。
6. 完了後に `surface.csv`（カラム: `i,j,k,d1_A,d2_A,d3_A,energy_hartree,bias_converged,energy_kcal,d1_label,d2_label,d3_label`）を組み立て、`--baseline {min|first}` で kcal/mol の基準を設定し、`--zmin/--zmax` に従った 3D RBF 補間等値面図 `scan3d_density.html` を生成します。`--csv` が指定された場合、この可視化ステップのみを実行します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | `--csv` 未指定時は必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート/`--ligand-charge`）。両方指定時は `-q` が優先 | テンプレート/導出がない場合は必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送）。診断上の注意は {ref}`ja-workers-fd-downgrade` を参照 | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（推奨）または **単一**のインライン Python リテラルで 3 つの四つ組 `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ | `--csv` 未指定時に必須 |
| `--one-based/--zero-based` | `(i, j)` のインデックスを 1 始まり/0 始まりとして解釈 | `True` |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のペア情報を表示 | `False` |
| `--max-step-size FLOAT` | 各距離の 1 増分あたりの最大変化量（Å）。グリッド密度を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル数。YAML で `opt.max_cycles` が指定されていない場合に使用 | `10000` |
| `--opt-mode TEXT` | `grad` → LBFGS、`hess` → RFOptimizer | `grad` |
| `--freeze-links/--no-freeze-links` | PDB 入力時にリンク水素の親原子を凍結 | `True` |
| `--freeze-atoms TEXT` | 凍結する原子の 1 始まりインデックスをカンマ区切りで明示的に指定（例: `'1,3,5'`）。`--freeze-links` と併用可、任意の入力形式に適用 | _None_ |
| `--dump/--no-dump` | 各 (d₁, d₂) ペアの `inner_path_d1_###_d2_###_trj.xyz` を保存 | `False` |
| `--convert-files/--no-convert-files` | PDB/Gaussian 入力で XYZ/TRJ → PDB/GJF 変換を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB トポロジー（XYZ 座標を保持） | _None_ |
| `-o, --out-dir TEXT` | グリッドとプロットの出力ディレクトリ | `./result_scan3d/` |
| `--csv PATH` | 既存の `surface.csv` を読み込みプロットのみ実行（新規スキャンなし）。指定時は `-i/--input` と `-s/--scan-lists` が任意になります | _None_ |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用） | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--preopt/--no-preopt` | スキャン前に無バイアス最適化を実行 | `False` |
| `--baseline {min,first}` | kcal/mol の基準をグローバル最小値または `(i,j,k)=(0,0,0)` に設定 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | 等値面の色範囲（kcal/mol） | 自動 |
| `--out-json/--no-out-json` | `out_dir` に機械可読な `result.json` を書き出す。スキーマは [JSON 出力スキーマ](json-output.md) を参照 | `False` |

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml-reference.md) と同じキーを使用します。`opt.dump` は YAML で設定可能ですが、軌跡出力は `--dump` で制御します。

```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 1-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p1 # uma-s-1p1 | uma-m-1p1
 device: auto # MLIP device selection
opt:
 thresh: baker # convergence preset (default: baker)
 max_cycles: 10000 # optimizer cycle cap
 dump: false # optimizer dumps (scan trajectories are controlled by --dump)
 out_dir: ./result_scan3d/ # output directory
lbfgs:
 max_step: 0.3 # maximum step length
 out_dir: ./result_scan3d/ # LBFGS-specific output directory
rfo:
 trust_radius: 0.10 # trust-region radius
 out_dir: ./result_scan3d/ # RFO-specific output directory
bias:
 k: 300.0 # harmonic bias strength (eV·Å⁻²)
```

`--relax-max-cycles` は**明示的に指定され**、かつ YAML で `opt.max_cycles` が設定されていない場合にのみ適用されます（デフォルト `10000`）。

### セクション `bias`
- `k`（`300`）: 調和バイアス強度（eV·Å⁻²）。

## 出力
```
out_dir/ (デフォルト:./result_scan3d/)
├─ surface.csv # グリッドメタデータ（i=j=k=-1 の参照行を含む場合あり）
├─ scan3d_density.html # 3D エネルギー等値面の可視化
├─ grid/point_i###_j###_k###.xyz # 各グリッド点の緩和構造（Å×100 タグ）
├─ grid/point_i###_j###_k###.pdb # 変換有効時の PDB コンパニオン
├─ grid/point_i###_j###_k###.gjf # テンプレートがある場合の Gaussian コンパニオン
├─ grid/preopt_i###_j###_k###.xyz # スキャン開始前の構造（--preopt の場合は最適化済み）
└─ grid/inner_path_d1_###_d2_###_trj.xyz # --dump の場合のみ（PDB 入力で変換有効時は .pdb も生成）
```

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 計算エンジンは MLIP バックエンド（デフォルト: UMA）で、1D/2D スキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å 単位の制限値は内部で Bohr に変換され、LBFGS ステップや RFO 信頼半径の制御に使われます。最適化の一時ファイルはテンポラリディレクトリに配置されます。
- `--baseline` はデフォルトでグローバル最小値を基準としてゼロにします。`--baseline first` は `(i,j,k)=(0,0,0)` の格子点を基準にします。
- 3D 可視化は 50×50×50 グリッドでの RBF 補間と、半透明の段階的等値面を使用します（断面表示はありません）。
- `--freeze-links` はユーザー指定の `freeze_atoms` にリンク水素親原子をマージし、抽出された活性部位モデルの境界を固定します。

## 関連項目
- [scan](scan.md) -- 1D 結合距離スキャン
- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [opt](opt.md) -- スキャン前後の単一構造最適化
- [all](all.md) -- 一気通貫ワークフロー
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
