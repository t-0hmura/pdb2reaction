# `scan2d`

## 概要

> **要約:** 調和拘束と MLIP 緩和（デフォルト: UMA、`-b/--backend` で ORB・MACE・AIMNet2 も選択可能）により、2 距離（d₁, d₂）のグリッドスキャンを行います。`-s/--scan-lists`（YAML/JSON ファイルパス、推奨）またはインライン Python リテラル を使用します。

### 要点
- **入力:** 1 つの構造 + `-s scan2d.yaml`（推奨）または `-s/--scan-lists` の **単一** インラインリテラル（四つ組はちょうど 2 つ）。
- **訪問順:** 各軸は（事前最適化された）構造に最も近い点を先に訪れるよう並べ替えられます。
- **エネルギー:** `surface.csv` の値は常に **バイアスなし**で評価されるため、格子点間で直接比較できます。
- **主な出力:** `surface.csv`、`scan2d_map.png`、`scan2d_landscape.html`、および `grid/` 配下の各点の構造。
- **注意:** `(high − low) / --max-step-size` が大きいと格子点数が急増します。

`scan2d` は `--max-step-size` に基づいて両軸の線形グリッドを作成し、各格子点を拘束付きで緩和してバイアスなしの UMA エネルギーを記録します。可視化用の出力も生成します。RFOptimizer を使用する場合は `--opt-mode hess` を指定してください。

XYZ/GJF 入力では、`--ref-pdb` で参照 PDB トポロジーを指定すると、XYZ 座標を保持したまま PDB/GJF へのフォーマット対応変換が可能になります。

## 最小例
```bash
pdb2reaction scan2d -i input.pdb -q 0 -s scan2d.yaml --out-dir ./result_scan2d/
```

## 出力の見方
- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_i000_j000.xyz`
- `result_scan2d/scan2d_map.png` と `result_scan2d/scan2d_landscape.html`

## よくある例

1. **YAML spec から実行する** — 下記の[例](#例)を参照。
2. **インラインリテラルを使う** — 下記の[例](#例)を参照。
3. **`--dump` を有効にして d1 ごとの内側軌跡を保存する** — 下記の[例](#例)を参照。

> **Note:** `-s/--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。

## 使用法
```bash
pdb2reaction scan2d -i INPUT.{pdb|xyz|trj|...} [-q CHARGE] [-l, --ligand-charge <number|'RES:Q,...'>] [-m MULT] \
 [-b/--backend uma|orb|mace|aimnet2] [--solvent SOLVENT] [--solvent-model alpb|cpcmx] \
 [-s/--scan-lists scan2d.yaml | '[(i,j,lowÅ,highÅ), (i,j,lowÅ,highÅ)]'] [options] \
 [--convert-files/--no-convert-files] [--ref-pdb FILE]
```

### 例
```bash
# 推奨: YAML/JSON spec
cat > scan2d.yaml << 'YAML'
one_based: true
pairs:
 - ["TYR,285,CA", "SAM,309,C10", 1.30, 3.10]
 - ["TYR,285,CB", "SAM,309,C11", 1.20, 3.20]
YAML
pdb2reaction scan2d -i input.pdb -q 0 -s scan2d.yaml

# 代替: Python リテラル
pdb2reaction scan2d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]'

# LBFGS、内側軌跡ダンプ、Plotly 出力
pdb2reaction scan2d -i input.pdb -q 0 \
 -s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]' \
 --max-step-size 0.20 --dump --out-dir ./result_scan2d/ --opt-mode grad \
 --preopt --baseline min
```

## YAML/JSON スペックファイルの書式（推奨）

```yaml
one_based: true # 任意。未指定時は CLI の --one-based を使用
pairs:
 - [1, 5, 1.30, 3.10]
 - [2, 8, 1.20, 3.20]
```

- `pairs` は必須で、ちょうど 2 つの四つ組を含む必要があります。
- 各四つ組は `(i, j, low_Å, high_Å)` です。
- インデックスは整数または PDB セレクタのどちらでも指定できます（インラインリテラルと同じ）。

## インライン Python リテラルの書式

`-s/--scan-lists` に **1 つのインライン Python リテラル文字列**を渡すこともできます。シェルのクォート処理に注意が必要です。

### 基本構造

リテラルは、ちょうど **2 つ**の四つ組 `(原子1, 原子2, 下限Å, 上限Å)` を含む Python リストです。

```
-s '[(原子1, 原子2, 下限Å, 上限Å), (原子3, 原子4, 下限Å, 上限Å)]'
```

- リスト全体を **シングルクォート** `'...'` で囲みます（シェルがカッコや空白を解釈しないようにするため）。
- 各四つ組は 1 つのスキャン軸を定義し、`原子1`–`原子2` 間の距離を `下限Å` から `上限Å` まで走査します。
- `scan` と異なり、リテラルは **1 つだけ**を受け付けます（複数ステージは非対応）。

### 原子の指定方法

原子は**整数インデックス**または **PDB セレクタ文字列**で指定します。

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 1.30, 3.10)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクタ | `("TYR,285,CA", "SAM,309,C10", 1.30, 3.10)` | 残基名、残基番号、原子名の三つ組 |

PDB セレクタのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれでも区切れます。トークンの順序も任意です。

```bash
# 以下はすべて同じ原子を指定:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # 順序は自由
```

### クォートの規則

```bash
# 正しい: シングルクォートでリストを囲み、セレクタはダブルクォート
-s '[("TYR,285,CA","SAM,309,C10",1.30,3.10),("TYR,285,CB","SAM,309,C11",1.20,3.20)]'

# 正しい: 整数インデックスなら内側のダブルクォートは不要
-s '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20)]'

# 非推奨: ダブルクォートで外側を囲むとエスケープが必要
-s "[(\"TYR,285,CA\",\"SAM,309,C10\",1.30,3.10),...]"
```

## ワークフロー
1. `geom_loader` で入力構造をロードし、電荷とスピンを解決します。`--preopt` の場合は無バイアスの事前最適化を実行します。`-q` が省略され `--ligand-charge` がある場合、構造は酵素--基質複合体として扱われ、PDB 入力（または `--ref-pdb` 付き XYZ/GJF）では `extract.py` の電荷サマリーから総電荷を導出します。事前最適化構造は `grid/preopt_i###_j###.*` に保存され、`surface.csv` には `i = j = -1` のエントリとしてバイアスなしエネルギーが記録されます。
2. `-s/--scan-lists`（YAML/JSON ファイルパスまたはインライン Python リテラル）を 2 つの四つ組に解析し、インデックスを正規化します（デフォルトは 1 始まり）。PDB 入力では、各エントリに整数インデックスまたは `'TYR,285,CA'` のようなセレクタ文字列を指定できます。区切りは空白・カンマ・スラッシュ・バッククォート・バックスラッシュのいずれも可で、トークン順序は任意です（フォールバックは resname, resseq, atom を想定）。線形グリッドは `ceil(|high − low| / h) + 1` 点（両端を含む）で構成します（`h = --max-step-size`）。長さ 0 の範囲は 1 点に縮退します。その後、各軸は事前最適化構造に最も近い距離が `i = 0` / `j = 0` になるよう並べ替えられます。
3. 外側ループで `d1[i]`（近い順）を走査します。各値で **d₁ 拘束のみ**を適用して緩和し、その構造をスナップショットとして保存します。次に内側ループで `d2[j]` を走査し、**d₁ と d₂ の両拘束**を適用して、最も近い既収束構造から緩和を開始します。
4. 各 `(i, j)` について、`<out-dir>/grid/point_i###_j###.xyz` に構造を保存し、バイアス収束の可否を記録し、バイアスを除去した UMA エネルギーを評価します。`--dump` の場合、外側ループごとの内側軌跡が `inner_path_d1_###_trj.xyz` として保存されます。
5. すべての点を走査したら、`<out-dir>/surface.csv` を作成します。`--baseline {min|first}` で kcal/mol の基準を設定します。`--baseline first` では基準点が `(low₁, low₂)` ではなく再並べ替え後の最初の格子点（`i = j = 0`）になります。`scan2d_map.png`（2D コンター）と `scan2d_landscape.html`（3D サーフェス）を `<out-dir>/` に生成します。`--zmin/--zmax` でカラースケールを固定できます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け入れる構造ファイル | 必須 |
| `-q, --charge INT` | 総電荷（CLI > テンプレート/`--ligand-charge`）。両方指定時は `-q` が優先 | テンプレート/導出がない場合は必須 |
| `-l, --ligand-charge TEXT` | 残基別電荷マッピング（例: `GPP:-3,SAM:1`）。PDB の残基電荷から全系の電荷を自動導出します（手動計算不要）。`-q` 省略時に使用（PDB 入力、または `--ref-pdb` 付き XYZ/GJF） | _None_ |
| `--workers`, `--workers-per-node` | UMA 予測器の並列度（workers > 1 で解析ヘシアン無効; `workers_per_node` は並列予測器へ転送） | `1`, `1` |
| `-m, --multiplicity INT` | スピン多重度 2S+1。`.gjf` テンプレートがあれば継承し、未指定時は `1` | `.gjf` テンプレート値または `1` |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（推奨）または **単一**のインライン Python リテラルで 2 つの四つ組 `(i,j,lowÅ,highÅ)` を指定。`i`/`j` は整数インデックスまたは PDB セレクタ（`'TYR,285,CA'`） | 必須 |
| `--one-based/--zero-based` | `(i, j)` のインデックスを 1 始まり/0 始まりとして解釈 | `True` |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のペア情報を表示。 | `False` |
| `--max-step-size FLOAT` | 各距離の 1 増分あたりの最大変化量（Å）。グリッド密度を決定 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV·Å⁻²） | `300` |
| `--relax-max-cycles INT` | 各バイアス緩和の最大最適化サイクル数。YAML で `opt.max_cycles` が指定されていない場合に使用 | `10000` |
| `--opt-mode TEXT` | `grad` → LBFGS、`hess` → RFOptimizer | `grad` |
| `--freeze-links/--no-freeze-links` | PDB 入力時にリンク水素の親原子を凍結 | `True` |
| `--dump/--no-dump` | 外側ループごとの `inner_path_d1_###_trj.xyz` を保存 | `False` |
| `--convert-files/--no-convert-files` | PDB/Gaussian 入力で XYZ/TRJ → PDB/GJF 変換を切り替え | `True` |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB トポロジー（XYZ 座標を保持） | _None_ |
| `-o, --out-dir TEXT` | 出力ディレクトリ | `./result_scan2d/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） | `baker` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用） | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP バックエンド | `uma` |
| `--solvent TEXT` | xTB 暗黙溶媒（例: `water`）。`none` で無効化 | `none` |
| `--solvent-model {alpb,cpcmx}` | xTB 溶媒モデル | `alpb` |
| `--preopt/--no-preopt` | スキャン前に無バイアス最適化を実行 | `True` |
| `--baseline {min,first}` | kcal/mol の基準をグローバル最小値または最初の格子点に設定 | `min` |
| `--zmin FLOAT`, `--zmax FLOAT` | カラースケールの下限/上限（kcal/mol） | 自動 |

### 共有 YAML セクション
- `geom`, `calc`, `opt`, `lbfgs`, `rfo`: [YAML リファレンス](yaml_reference.md) と同じキーを使用します。`opt.dump` は YAML で設定可能ですが、スキャン軌跡の出力は `--dump` で制御します。

### セクション `bias`
- `k`（`300`）: 調和バイアス強度（eV·Å⁻²）。

## 出力
```
out_dir/ (デフォルト:./result_scan2d/)
├─ surface.csv # 構造化グリッド表
├─ scan2d_map.png # 2D コンター（Kaleido 必須; PNG 出力に失敗すると実行が停止）
├─ scan2d_landscape.html # 3D サーフェス可視化
├─ grid/point_i###_j###.xyz # 各 (i, j) の緩和構造
├─ grid/point_i###_j###.pdb # 変換有効時の PDB コンパニオン
├─ grid/point_i###_j###.gjf # テンプレートがある場合の Gaussian コンパニオン
└─ grid/inner_path_d1_###_trj.xyz # --dump の場合のみ（PDB 入力時は.pdb にも変換）
```

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- MLIP バックエンド（デフォルト: UMA、`-b/--backend` で切替可能）が計算エンジンで、1D スキャンと同じ `HarmonicBiasCalculator` を再利用します。
- Å 単位の制限値は内部で Bohr に変換され、LBFGS ステップや RFO 信頼半径の制御に使われます。最適化の一時ファイルはテンポラリディレクトリに配置されます。
- バイアスはエネルギー記録前に除去されるため、`surface.csv` を下流のフィッティングや可視化スクリプトにそのまま利用できます。
- `--freeze-links` はユーザー指定の `freeze_atoms` にリンク水素親原子をマージし、抽出ポケットの境界を固定します。


```yaml
geom:
 coord_type: cart # coordinate type: cartesian vs dlc internals
 freeze_atoms: [] # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0 # total charge (CLI/template override)
 spin: 1 # spin multiplicity 2S+1
 model: uma-s-1p2 # uma-s-1p1 | uma-s-1p2 | uma-m-1p1
 device: auto # UMA device selection
opt:
 thresh: baker # convergence preset (default: baker)
 max_cycles: 10000 # optimizer cycle cap
 dump: false # optimizer dumps (scan trajectories are controlled by --dump)
 out_dir: ./result_scan2d/ # output directory
lbfgs:
 max_step: 0.3 # maximum step length
 out_dir: ./result_scan2d/ # LBFGS-specific output directory
rfo:
 trust_radius: 0.1 # trust-region radius
 out_dir: ./result_scan2d/ # RFO-specific output directory
bias:
 k: 300.0 # harmonic bias strength (eV·Å⁻²)
```

`opt` の詳細は [YAML リファレンス](yaml_reference.md) を参照してください。
`--relax-max-cycles` は**明示的に指定され**、かつ YAML で `opt.max_cycles` が設定されていない場合にのみ適用されます（デフォルト `10000`）。

## 関連項目
- [scan](scan.md) -- 1D 結合距離スキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- スキャン前後の単一構造最適化
- [all](all.md) -- end-to-endワークフロー
- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
